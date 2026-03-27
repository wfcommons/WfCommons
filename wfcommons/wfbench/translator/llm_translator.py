"""
LLM-based forward translator for WFCommons WFBench.

Translates WfFormat JSON INTO WMS-specific workflow code using an LLM.
Extends the abstract Translator class, matching the same interface as
the hand-coded translators (Nextflow, Pegasus, CWL, etc.).
"""

import json
import logging
import pathlib
import re
from typing import Optional, Union, List

from dotenv import load_dotenv

from .abstract_translator import Translator
from .utils.llm_client import LLMClient, client_from_yaml
from .utils.wfinstances import retrieve_instances
from .skills.loader import SkillLoader
from ...common import Workflow

load_dotenv()

FORWARD_SKILLS_DIR = pathlib.Path(__file__).resolve().parent / "skills" / "forward"


class LLMTranslator(Translator):
    """
    An LLM-based WfFormat-to-WMS translator.

    Uses an LLM to translate a WfFormat workflow specification into
    executable code for a target workflow management system (e.g., Nextflow,
    CWL, Airflow). Supports the same ``translate(output_folder)`` interface
    as all other concrete translators.

    :param workflow: Workflow benchmark object or path to the WfFormat JSON instance.
    :type workflow: Union[Workflow, pathlib.Path]
    :param target_system: Target WMS name (e.g. "nextflow", "cwl", "airflow").
    :type target_system: str
    :param llm_client: A pre-configured LLMClient. Provide either this or ``model_name``.
    :type llm_client: LLMClient, optional
    :param model_name: Key from models.yaml to auto-build an LLMClient.
    :type model_name: str, optional
    :param models_file: Path to a custom models YAML file.
    :type models_file: str or pathlib.Path, optional
    :param context_examples: List of example WMS scripts/files to include as style reference.
    :type context_examples: List[str], optional
    :param examples_instances: WfInstances queries to fetch WfFormat examples as LLM context.
        Accepts application names (``"blast"``), specific instance filenames
        (``"blast-chameleon-small-001"``), or full repo paths (``"makeflow/blast"``).
        Can be a single string or a list.
    :type examples_instances: Union[str, List[str]], optional
    :param num_examples: Max number of WfInstances examples to include. Defaults to 3.
    :type num_examples: int
    :param system_prompt: Override the default system instructions (bypasses skills).
    :type system_prompt: str, optional
    :param logger: Logger instance.
    :type logger: logging.Logger, optional
    """

    def __init__(self,
                 workflow: Union[Workflow, pathlib.Path],
                 target_system: str,
                 llm_client: Optional[LLMClient] = None,
                 model_name: Optional[str] = None,
                 models_file: Optional[Union[str, pathlib.Path]] = None,
                 context_examples: Optional[List[str]] = None,
                 examples_instances: Optional[Union[str, List[str]]] = None,
                 num_examples: int = 3,
                 system_prompt: Optional[str] = None,
                 logger: Optional[logging.Logger] = None) -> None:
        """Create an object of the LLM forward translator."""
        super().__init__(workflow, logger)

        if llm_client is None and model_name is None:
            raise ValueError("Provide either llm_client or model_name.")
        if llm_client is not None and model_name is not None:
            raise ValueError("Provide only one of llm_client or model_name, not both.")

        if model_name is not None:
            llm_client = client_from_yaml(model_name, models_file)

        self.llm = llm_client
        self.target_system = target_system.lower()
        self.context_examples = context_examples or []
        if isinstance(examples_instances, str):
            examples_instances = [examples_instances]
        self.examples_instances = examples_instances
        self.num_examples = num_examples
        self._system_prompt_override = system_prompt
        self._skill_loader = SkillLoader(skills_dir=FORWARD_SKILLS_DIR)

    def translate(self, output_folder: pathlib.Path) -> None:
        """
        Translate a workflow benchmark description (WfFormat) into a WMS workflow
        application using an LLM.

        :param output_folder: The path to the folder in which the workflow benchmark will be generated.
        :type output_folder: pathlib.Path
        """
        output_folder.mkdir(parents=True, exist_ok=True)

        wfformat_json = json.dumps(self._trim_workflow_json(), indent=1)

        # Fetch WfInstances examples if configured
        wf_examples = []
        if self.examples_instances:
            wf_examples = retrieve_instances(
                self.examples_instances,
                num_examples=self.num_examples,
                score_against=wfformat_json,
            )

        prompt = self._build_prompt(wfformat_json, wf_examples=wf_examples)

        self.logger.info(f"Sending WfFormat to LLM for {self.target_system} translation...")
        raw_output = self.llm.complete(prompt)

        code = self._extract_code(raw_output)

        output_file = output_folder / self._default_output_filename()
        self._write_output_file(code, output_file)

        # Skip _copy_binary_files and _generate_input_files — the LLM
        # translator produces standalone WMS code, not a wfbench harness.

    def _build_prompt(self,
                      wfformat_json: str,
                      wf_examples: Optional[list] = None) -> str:
        """
        Build the full LLM prompt from skills, WfFormat input, examples, and context.

        :param wfformat_json: The serialized WfFormat JSON.
        :type wfformat_json: str
        :param wf_examples: WfInstances examples fetched from GitHub.
        :type wf_examples: list, optional
        :return: The composed prompt.
        :rtype: str
        """
        if self._system_prompt_override is not None:
            system_prompt = self._system_prompt_override
        else:
            system_prompt = self._skill_loader.compose_prompt(
                trace_text=self.target_system,
                skill_name=self.target_system,
            )

        prompt = system_prompt.strip() + "\n\n"

        prompt += f"=== TARGET SYSTEM: {self.target_system.upper()} ===\n\n"

        # Include WfInstances examples so the LLM sees real WfFormat structures
        if wf_examples:
            prompt += "=== REFERENCE WFFORMAT INSTANCES (from WfInstances) ===\n"
            prompt += ("These are real WfFormat workflow instances for reference. "
                       "They show the structure and conventions used in practice.\n")
            for i, ex in enumerate(wf_examples, 1):
                prompt += f"\n--- Instance {i}: {ex['filename']} ---\n"
                prompt += ex["content"][:8000]  # safety truncation
                prompt += "\n"
            prompt += "\n"

        prompt += "=== WFFORMAT INPUT (translate this) ===\n"
        prompt += wfformat_json + "\n"

        if self.context_examples:
            prompt += "\n=== CONTEXT EXAMPLES (target system code for reference) ===\n"
            for i, example in enumerate(self.context_examples, 1):
                prompt += f"\n--- Example {i} ---\n"
                prompt += example[:10000]  # safety truncation
                prompt += "\n"

        prompt += (
            f"\n=== OUTPUT REQUIREMENTS ===\n"
            f"Produce ONLY valid, executable {self.target_system} workflow code.\n"
            f"The code must faithfully represent the workflow from the WfFormat input.\n"
            f"Respect all task dependencies, input/output files, and task ordering.\n"
            f"Do not include explanations or markdown formatting.\n"
        )

        return prompt

    def _trim_workflow_json(self) -> dict:
        """
        Produce a trimmed copy of the workflow JSON for the LLM prompt.

        Strips fields that are not needed for code generation:
        - Execution tasks: keep only id, command, coreCount
        - Specification tasks: simplify inputFiles/outputFiles to just IDs
        - Drop the top-level files array and machines list
        """
        import copy
        wf = copy.deepcopy(self.workflow.workflow_json)

        spec = wf.get("workflow", {}).get("specification", {})

        # Simplify spec tasks: file references -> just IDs
        if "tasks" in spec:
            for task in spec["tasks"]:
                for key in ("inputFiles", "outputFiles"):
                    if key in task and task[key]:
                        task[key] = [f["id"] if isinstance(f, dict) else f
                                     for f in task[key]]
            # Drop the redundant top-level files list
            spec.pop("files", None)

        # Strip execution tasks to essentials
        exe = wf.get("workflow", {}).get("execution", {})
        if "tasks" in exe:
            keep_keys = {"id", "command", "coreCount"}
            exe["tasks"] = [
                {k: v for k, v in t.items() if k in keep_keys}
                for t in exe["tasks"]
            ]
        exe.pop("machines", None)

        return wf

    @staticmethod
    def _extract_code(output: str) -> str:
        """
        Extract code from LLM output, stripping markdown code blocks if present.

        :param output: Raw LLM output.
        :type output: str
        :return: Extracted code.
        :rtype: str
        """
        code_block_match = re.search(r'```(?:\w+)?\s*([\s\S]*?)\s*```', output)
        if code_block_match:
            return code_block_match.group(1).strip()
        return output.strip()

    def _default_output_filename(self) -> str:
        """
        Return the default output filename for the target system.

        :return: Filename string.
        :rtype: str
        """
        extensions = {
            "nextflow": "workflow.nf",
            "cwl": "workflow.cwl",
            "airflow": "workflow_dag.py",
            "bash": "run_workflow.sh",
            "pegasus": "run_workflow.py",
            "parsl": "run_workflow.py",
            "dask": "run_workflow.py",
            "swift_t": "workflow.swift",
            "pycompss": "run_workflow.py",
            "taskvine": "run_workflow.py",
            "radical_pilot": "run_workflow.py",
        }
        return extensions.get(self.target_system, f"workflow.{self.target_system}")

    @staticmethod
    def available_skills() -> list[str]:
        """Return the list of available forward skill names."""
        return SkillLoader(skills_dir=FORWARD_SKILLS_DIR).available_skills()

"""
LLM-based backward translator for WFCommons WFBench.

Translates WMS traces/logs INTO WfFormat JSON using an LLM.
This is the reverse direction of the standard translators.
"""

from pathlib import Path

from dotenv import load_dotenv
from typing import Optional, Dict, Any, List
from wfcommons.wfbench.translator.utils.llm_client import (
    LLMClient, client_from_yaml,
)
from wfcommons.wfbench.translator.utils.wfinstances import (
    load_instances, retrieve_instances, _word_overlap,
)
from wfcommons.wfbench.translator.skills.loader import SkillLoader

load_dotenv()  # loads .env from cwd (project root)

BACKWARD_SKILLS_DIR = Path(__file__).resolve().parent / "skills" / "backward"


class LLMBackwardTranslator():
    """
    • Uses existing WfFormat as examples.
    • Accepts a trace from a NEW workflow system.
    • Sends all of this as grounding context to an LLM.
    • Produces a new trace in WfFormat automatically.

    The user does not implement translation logic.
    """

    def __init__(self,
                llm_client: LLMClient | None = None,
                model_name: str | None = None,
                models_file: str | Path | None = None,
                examples_instances: Optional[List[str]] = None,
                num_examples: int = 3,
                system_prompt: Optional[str] = None,
                skill_name: Optional[str] = None,
                **kwargs
            ):
        """
        Parameters
        ----------
        llm_client : LLMClient, optional
            A pre-configured LLMClient instance. Either this or
            ``model_name`` must be provided.
        model_name : str, optional
            Key from models.yaml (e.g. "qwen3", "ollama/llama3").
            The matching config is used to build an LLMClient automatically.
        models_file : str or Path, optional
            Path to a custom models YAML file.  Defaults to the
            ``models.yaml`` shipped alongside this module.
        example_instances : List[str]
             URLs pointing to translator examples or benchmarks:
              - raw GitHub links
              - JSON traces
              - scripts
        num_examples : int, optional
            Number of example instances to include in the prompt.
        system_prompt : str, optional
            Override the default system instructions for the LLM.
            When provided, the skills system is bypassed entirely.
        skill_name : str, optional
            Explicit skill to use (e.g. "nextflow", "cwl").
            Auto-detected from trace content if not specified.
        kwargs : dict
            Additional parameters passed to the parent Translator if needed.
        """
        super().__init__(**kwargs)

        if llm_client is None and model_name is None:
            raise ValueError("Provide either llm_client or model_name.")
        if llm_client is not None and model_name is not None:
            raise ValueError("Provide only one of llm_client or model_name, not both.")

        if model_name is not None:
            llm_client = client_from_yaml(model_name, models_file)

        self.llm = llm_client
        if isinstance(examples_instances, str):
            examples_instances = [examples_instances]
        self.examples_instances = examples_instances
        self.num_examples = num_examples
        self.skill_name = skill_name
        self._skill_loader = SkillLoader(skills_dir=BACKWARD_SKILLS_DIR)
        self._system_prompt_override = system_prompt

    @staticmethod
    def available_skills() -> list[str]:
        """Return the list of available backward skill names."""
        return SkillLoader(skills_dir=BACKWARD_SKILLS_DIR).available_skills()
    
    def translate(self, trace, metadata=None, json_schema: dict | None = None, **kwargs):
        # --- Normalize the trace into a string ---
        if isinstance(trace, dict):
            import json
            trace_text = json.dumps(trace, indent=2)
        elif isinstance(trace, (list, tuple)):
            trace_text = "\n".join(map(str, trace))
        else:
            # assume Python code or any raw text
            trace_text = str(trace)


        examples = self._retrieve_examples(trace_text)

        prompt = self._build_prompt(
            trace=trace_text,
            examples=examples,
            metadata=metadata,
        )

        response_format = None
        if json_schema is not None:
            response_format = {
                "type": "json_schema",
                "json_schema": {
                    "name": "WfFormat",
                    "schema": json_schema
                }
            }

        output = self.llm.complete(prompt, response_format=response_format)
        return output
    
    def _retrieve_examples(self, trace_text: str) -> List[Dict[str, Any]]:
        """
        Fetch and rank WfInstances examples by relevance to the trace.

        Returns an empty list if no examples_instances were provided or
        none of the paths exist in the WfInstances repository.
        """
        return retrieve_instances(
            self.examples_instances or [],
            num_examples=self.num_examples,
            score_against=trace_text,
        )
    
    def _build_prompt(
        self,
        trace: str,
        examples: List[Dict[str, Any]] = [],
        metadata: Optional[Dict[str, Any]] = None,
    ) -> str:

        if self._system_prompt_override is not None:
            # Explicit override: use as-is (backward compat)
            system_prompt = self._system_prompt_override
        else:
            # Compose from skill files
            skill_hint = self.skill_name
            if not skill_hint and metadata and "source_system" in metadata:
                skill_hint = metadata["source_system"].lower()
            system_prompt = self._skill_loader.compose_prompt(
                trace_text=trace,
                skill_name=skill_hint,
            )

        prompt = system_prompt.strip() + "\n\n"

        prompt += "=== EXAMPLE WORKFLOW INSTANCES (WFFORMAT) ===\n"
        for i, ex in enumerate(examples, 1):
            prompt += f"\n--- Example {i} ---\n"
            prompt += f"Source URL: {ex['url']}\n"
            prompt += "Content:\n"
            prompt += ex["content"][:5000]  # safety truncation
            prompt += "\n"

        prompt += "\n=== NEW WORKFLOW TRACE TO TRANSLATE (e.g., dispel4py) ===\n"
        prompt += trace + "\n"

        if metadata:
            prompt += "\n=== ADDITIONAL METADATA ===\n"
            for k, v in metadata.items():
                prompt += f"- {k}: {v}\n"

        prompt += (
            "\n=== OUTPUT REQUIREMENTS ===\n"
            "Produce ONLY a JSON object compatible with WorkflowBenchmark.from_dict().\n"
            "Infer tasks, dependencies, runtimes, and workflow structure.\n"
            "Do not include explanations.\n"
        )

        return prompt
    
    def _parse_llm_output(self, output: str) -> Dict[str, Any]:
        import json
        import re

        # Try direct parse first
        try:
            return json.loads(output)
        except json.JSONDecodeError:
            pass

        # Extract JSON from markdown code blocks
        code_block_match = re.search(r'```(?:json)?\s*([\s\S]*?)\s*```', output)
        if code_block_match:
            try:
                return json.loads(code_block_match.group(1))
            except json.JSONDecodeError:
                pass

        # Extract JSON object by finding first { and last }
        first_brace = output.find('{')
        last_brace = output.rfind('}')
        if first_brace != -1 and last_brace != -1 and last_brace > first_brace:
            try:
                return json.loads(output[first_brace:last_brace + 1])
            except json.JSONDecodeError:
                pass

        raise ValueError("Could not extract valid JSON from LLM output.")
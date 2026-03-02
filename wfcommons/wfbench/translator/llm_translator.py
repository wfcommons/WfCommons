"""
LLM-based translator scaffolding for WFCommons WFBench.

This translator inherits from the existing abstract Translator class
from wfcommons.wfbench.translator and implements the required interface.
"""

import os
import re
from pathlib import Path

import requests
import yaml
from dotenv import load_dotenv
from typing import Optional, Dict, Any, List
from wfcommons.wfbench.translator.utils.llm_client import LLMClient

load_dotenv()  # loads .env from cwd (project root)

MODELS_YAML = Path("models.yaml")


class LLMTranslator():
    """
    • Uses existing WfFormat as examples.
    • Accepts a trace from a NEW workflow system.
    • Sends all of this as grounding context to an LLM.
    • Produces a new recipe automatically.

    The user does not implement translation logic.
    """

    def __init__(self,
                llm_client: LLMClient | None = None,
                model_name: str | None = None,
                models_file: str | Path | None = None,
                examples_instances: Optional[List[str]] = None,
                num_examples: int = 3,
                system_prompt: Optional[str] = None,
                **kwargs
            ):
        """
        Parameters
        ----------
        llm_client : LLMClient, optional
            A pre-configured LLMClient instance.  Either this or
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
        kwargs : dict
            Additional parameters passed to the parent Translator if needed.
        """
        super().__init__(**kwargs)

        if llm_client is None and model_name is None:
            raise ValueError("Provide either llm_client or model_name.")
        if llm_client is not None and model_name is not None:
            raise ValueError("Provide only one of llm_client or model_name, not both.")

        if model_name is not None:
            llm_client = self._client_from_yaml(model_name, models_file)

        self.llm = llm_client
        self.examples_instances = examples_instances
        self.num_examples = num_examples
        self.system_prompt = system_prompt or DEFAULT_SYSTEM_PROMPT

    # ------------------------------------------------------------------ #
    #  YAML helpers                                                       #
    # ------------------------------------------------------------------ #

    @staticmethod
    def available_models(models_file: str | Path | None = None) -> list[str]:
        """Return the list of model keys defined in models.yaml."""
        cfg = LLMTranslator._load_models_yaml(models_file)
        return list(cfg.keys())

    @staticmethod
    def _load_models_yaml(models_file: str | Path | None = None) -> dict:
        path = Path(models_file) if models_file else MODELS_YAML
        with open(path) as f:
            return yaml.safe_load(f)

    @staticmethod
    def _resolve_env(value: str) -> str:
        """Replace ``${VAR}`` placeholders with environment variables."""
        def _replace(m):
            var = m.group(1)
            val = os.environ.get(var)
            if val is None:
                raise EnvironmentError(
                    f"Environment variable '{var}' is not set "
                    f"(required by models.yaml)."
                )
            return val
        return re.sub(r"\$\{(\w+)\}", _replace, value)

    @staticmethod
    def _client_from_yaml(model_name: str,
                          models_file: str | Path | None = None) -> LLMClient:
        cfg = LLMTranslator._load_models_yaml(models_file)
        if model_name not in cfg:
            raise KeyError(
                f"Model '{model_name}' not found in models.yaml. "
                f"Available: {list(cfg.keys())}"
            )
        entry = cfg[model_name]
        api_key = LLMTranslator._resolve_env(str(entry["api_key"]))
        base_url = entry.get("base_url")
        if base_url and base_url != "null":
            base_url = LLMTranslator._resolve_env(str(base_url))
        else:
            base_url = None
        return LLMClient(
            model=entry["model"],
            api_key=api_key,
            base_url=base_url,
        )
    
    def _load_examples(self,
                       path_list: List[str] , 
                       ref: str="main") -> List[str]:
        """
        Load and return the content from the provided URLs.
        
        Parameters
        ----------
        path_list : List[str]
            List of paths within the repository to fetch files from.
        
        Returns
        -------
        List[Dict[str, Any]]
            List of dictionaries with 'url', 'filename', and 'content' keys.
        """
        all_examples = {}
        for path in path_list:
            examples = self._fetch_examples_from_path(
                path=path,
                ref=ref
            )
            all_examples[path] = examples
        return all_examples

     
    def _fetch_examples_from_path(self,
                                  path: str,
                                  ref: str="main") -> List[Dict[str, Any]]:
        """
        Fetch Python files from a specific path in a GitHub repository.
        
        Parameters
        ----------
        path : str
            Path within the repository to fetch files from.
        ref : str, optional
            Git reference (branch, tag, or commit SHA). Defaults to "main".

        Returns         
        -------
        List[Dict[str, Any]]
            List of dictionaries with 'url', 'filename', and 'content' keys.
        """
        print(f"Fetching examples from path: {path} at ref: {ref}")
        url = f"https://api.github.com/repos/wfcommons/WfInstances/contents/{path}?ref={ref}"
        listing = requests.get(url).json()

        examples = []
        for item in listing:
            if item["type"] == "file" and item["name"].endswith(".json") and not item["name"].endswith(".md"):
                raw = requests.get(item["download_url"]).text
                examples.append({
                    "url": item["download_url"],
                    "filename": item["name"],
                    "content": raw
                })
        return examples
    
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


        prompt = self._build_prompt(
            trace=trace_text,
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
    
    def _retrieve_examples(self, trace_text: str):
        """
        Simple scoring method to choose top-k (num_examples) examples
        Replace with embeddings if desired.
        """
        examples = self._load_examples(self.examples_instances)
        flat_examples = []
        if isinstance(examples, dict):
            for v in examples.values():
                flat_examples.extend(v)
        else:
            flat_examples = list(examples)

        results = []
        for example in flat_examples:
            score = self._similarity(trace_text, example["content"])
            results.append((score, example))

        results.sort(reverse=True, key=lambda x: x[0])
        return [ex for _, ex in results[: self.num_examples]]
    
    @staticmethod
    def _similarity(a: str, b: str) -> float:
        """
        Very naive similarity based on word overlap. Replace with embeddings if desired.
        """
        return len(set(a.split()) & set(b.split()))
    
    def _build_prompt(
        self,
        trace: str,
        examples: List[Dict[str, Any]] = [],
        metadata: Optional[Dict[str, Any]] = None,
    ) -> str:

        prompt = self.system_prompt.strip() + "\n\n"

        prompt += "=== EXAMPLE TRANSLATORS (FROM URLS) ===\n"
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
    
DEFAULT_SYSTEM_PROMPT = """
You are an expert software engineer specializing in workflow systems.
Translate workflow definitions/traces into WfCommons WfFormat 1.5 JSON.

OUTPUT THIS EXACT STRUCTURE:
{
  "name": "<workflow_name - REQUIRED>",
  "schemaVersion": "1.5",
  "workflow": {
    "specification": {
      "tasks": [
        {
          "name": "<task_name>",
          "id": "<task_id>",
          "parents": [],
          "children": [],
          "inputFiles": [],
          "outputFiles": []
        }
      ],
      "files": []
    },
    "execution": {
      "makespanInSeconds": <number or 0 if unknown>,
      "executedAt": "<timestamp or "1970-01-01T00:00:00Z" if unknown>",
      "tasks": [
        {
          "id": "<task_id matching specification>",
          "runtimeInSeconds": <number or 0 if unknown>,
          "executedAt": "<timestamp or "1970-01-01T00:00:00Z" if unknown>",
          "command": {
            "program": "<program name>",
            "arguments": []
          },
          "coreCount": <number or 1>,
          "avgCPU": <percentage or 0>,
          "readBytes": <number or 0>,
          "writtenBytes": <number or 0>,
          "memoryInBytes": <number or 0>,
          "machines": ["unknown"]
        }
      ],
      "machines": [
        {
          "nodeName": "unknown"
        }
      ]
    }
  }
}

RULES:
1. Use EXACTLY this structure - do not add or rename fields
2. "name" is REQUIRED - use the workflow name from the file, or the filename from metadata
3. "schemaVersion" is always "1.5"
4. Do NOT include optional top-level fields like "description", "createdAt", "author", "runtimeSystem" unless explicitly provided in source
5. For arrays not found, use empty array []
6. For numbers not found, use 0
7. For timestamp strings not found, use "1970-01-01T00:00:00Z"
8. Each task in specification MUST have: name, id, parents, children
9. Each task in execution MUST have: id (matching specification), runtimeInSeconds
10. Infer task dependencies from data flow (channels, inputs/outputs)
11. Only populate execution fields if runtime data exists in the source - otherwise use 0 or placeholder values
12. Output ONLY valid JSON - no explanations or markdown
"""
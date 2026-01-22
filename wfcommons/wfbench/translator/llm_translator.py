"""
LLM-based translator scaffolding for WFCommons WFBench.

This translator inherits from the existing abstract Translator class
from wfcommons.wfbench.translator and implements the required interface.
"""

from os import path
import requests
from wfcommons.wfbench.bench import WorkflowBenchmark
from typing import Optional, Dict, Any, List
from wfcommons.wfbench.translator.utils.llm_client import LLMClient


class LLMTranslator():
    """
    • Uses existing WfFormat as examples.
    • Accepts a trace from a NEW workflow system.
    • Sends all of this as grounding context to an LLM.
    • Produces a new recipe automatically.

    The user does not implement translation logic.
    """

    def __init__(self,
                llm_client: LLMClient,
                examples_instances: Optional[List[str]],
                num_examples: int = 3,
                system_prompt: Optional[str] = None,
                **kwargs
            ):
        """
        Parameters
        ----------
        llm_client : Any
            An object with `.complete(prompt: str) -> str`.
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

        self.llm = llm_client
        self.examples_instances = examples_instances
        self.num_examples = num_examples
        self.system_prompt = system_prompt or DEFAULT_SYSTEM_PROMPT
    
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

        # grounding_examples = self._retrieve_examples(trace_text)

        prompt = self._build_prompt(
            trace=trace_text,
            # examples=grounding_examples,
            metadata=metadata,
        )

        output = self.llm.complete(
            prompt,
            response_format={
                "type": "json_schema",
                "json_schema": {
                    "name": "WfFormat",
                    "schema": json_schema
                }
            }
        )
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
        Very naive similarity; replace with embeddings for real use.
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
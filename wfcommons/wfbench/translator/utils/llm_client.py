"""
LLM Client and model configuration utilities for WfCommons LLM Translators.

Uses the OpenAI SDK which is compatible with many LLM providers:
- OpenAI (GPT models)
- Ollama (local models like Llama3) via OpenAI-compatible endpoint
- vLLM, LM Studio, LocalAI, and other OpenAI-compatible servers
"""

import os
import re
from pathlib import Path

import yaml
from openai import OpenAI
from pydantic import BaseModel

MODELS_YAML = Path(__file__).resolve().parent.parent / "models.yaml"


class LLMClient:
    """
    Universal LLM client using OpenAI SDK.

    Examples
    --------
    # OpenAI
    client = LLMClient(model="gpt-4o-mini", api_key="sk-...")

    # Ollama (local Llama3)
    client = LLMClient(
        model="llama3",
        base_url="http://localhost:11434/v1",
        api_key="ollama"  # Ollama requires a dummy key
    )

    """

    def __init__(self,
                 model: str = "gpt-4o-mini",
                 api_key: str | None = None,
                 base_url: str | None = None):
        """
        Initialize the LLM client.

        Parameters
        ----------
        model : str, optional
            Model name to use. Defaults to "gpt-4o-mini".
        api_key : str, optional
            API key. For OpenAI, uses OPENAI_API_KEY env var if not provided.
            For local servers like Ollama, use any non-empty string (e.g., "ollama").
        base_url : str, optional
            Base URL for the API. Defaults to OpenAI's API.
            For Ollama: "http://localhost:11434/v1"
            For vLLM: "http://localhost:8000/v1"
        """
        self.model = model
        self.client = OpenAI(api_key=api_key, base_url=base_url)

    def complete(self,
                 prompt: str,
                 response_format: dict | type[BaseModel] | None = None,
                 temperature: float = 0.0) -> str:
        """
        Generate a completion for the given prompt.

        Parameters
        ----------
        prompt : str
            The input prompt to send to the LLM.
        response_format : dict or BaseModel, optional
            Format specification for structured output (e.g., JSON schema).
            Note: Not all providers support this parameter.
        temperature : float, optional
            Sampling temperature (default: 0.0 for deterministic output).
        Returns
        -------
        str
            The generated completion text.
        """
        kwargs = {
            "model": self.model,
            "messages": [{"role": "user", "content": prompt}],
            "temperature": temperature
        }

        if response_format is not None:
            kwargs["response_format"] = response_format

        response = self.client.chat.completions.create(**kwargs)
        return response.choices[0].message.content


def load_models_yaml(models_file: str | Path | None = None) -> dict:
    """Load and parse a models YAML configuration file."""
    path = Path(models_file) if models_file else MODELS_YAML
    with open(path) as f:
        return yaml.safe_load(f)


def resolve_env(value: str) -> str:
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


def client_from_yaml(model_name: str,
                     models_file: str | Path | None = None) -> LLMClient:
    """Build an LLMClient from a models.yaml entry."""
    cfg = load_models_yaml(models_file)
    if model_name not in cfg:
        raise KeyError(
            f"Model '{model_name}' not found in models.yaml. "
            f"Available: {list(cfg.keys())}"
        )
    entry = cfg[model_name]
    api_key = resolve_env(str(entry["api_key"]))
    base_url = entry.get("base_url")
    if base_url and base_url != "null":
        base_url = resolve_env(str(base_url))
    else:
        base_url = None
    return LLMClient(
        model=entry["model"],
        api_key=api_key,
        base_url=base_url,
    )


def available_models(models_file: str | Path | None = None) -> list[str]:
    """Return the list of model keys defined in models.yaml."""
    cfg = load_models_yaml(models_file)
    return list(cfg.keys())

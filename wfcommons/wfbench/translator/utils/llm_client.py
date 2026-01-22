"""
LLM Client for WfCommons LLM Translator.

Uses the OpenAI SDK which is compatible with many LLM providers:
- OpenAI (GPT models)
- Ollama (local models like Llama3) via OpenAI-compatible endpoint
- vLLM, LM Studio, LocalAI, and other OpenAI-compatible servers
"""

from openai import OpenAI
from pydantic import BaseModel


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
            Sampling temperature (0.0 = deterministic). Defaults to 0.0.

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

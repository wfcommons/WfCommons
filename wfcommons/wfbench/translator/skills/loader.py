"""
Skill file loader for the LLM translator.

Discovers, parses, and composes .md skill files from the skills/ directory
into prompts for LLM-based workflow translation.
"""

import re
from pathlib import Path
from typing import Optional


SKILLS_DIR = Path(__file__).resolve().parent


class Skill:
    """Parsed representation of a single .md skill file."""

    def __init__(self, name: str, filepath: Path, triggers: list[str],
                 description: str, domain_knowledge: str, examples: str,
                 validation: dict | None = None):
        self.name = name
        self.filepath = filepath
        self.triggers = triggers
        self.description = description
        self.domain_knowledge = domain_knowledge
        self.examples = examples
        self.validation = validation or {}

    def matches(self, text: str) -> bool:
        """Check if any trigger appears in the text (case-insensitive)."""
        text_lower = text.lower()
        return any(t.strip().lower() in text_lower
                   for t in self.triggers if t.strip())

    @property
    def required_elements(self) -> list[str]:
        """Strings that must appear in generated code."""
        return self.validation.get("required_elements", [])

    @property
    def syntax_check_command(self) -> str | None:
        """Shell command for syntax checking. '{file}' is substituted with the output file."""
        return self.validation.get("syntax_check_command")


class SkillLoader:
    """Loads and composes skill files for prompt building."""

    def __init__(self, skills_dir: Path | str | None = None):
        self.skills_dir = Path(skills_dir) if skills_dir else SKILLS_DIR
        self._skills: dict[str, Skill] = {}
        self._base_skill: Skill | None = None
        self._load_all()

    def _load_all(self):
        """Discover and parse all .md files in the skills directory."""
        for md_file in sorted(self.skills_dir.glob("*.md")):
            skill = self._parse_skill_file(md_file)
            if md_file.stem == "base":
                self._base_skill = skill
            else:
                self._skills[md_file.stem] = skill

    @staticmethod
    def _parse_skill_file(filepath: Path) -> Skill:
        """Parse a markdown skill file into a Skill object."""
        content = filepath.read_text()

        def extract_section(header: str) -> str:
            pattern = rf"^## {re.escape(header)}\s*\n(.*?)(?=^## |\Z)"
            match = re.search(pattern, content, re.MULTILINE | re.DOTALL)
            return match.group(1).strip() if match else ""

        # Extract name from H1
        name_match = re.search(r"^# Skill:\s*(.+)", content, re.MULTILINE)
        name = name_match.group(1).strip() if name_match else filepath.stem

        triggers_raw = extract_section("Triggers")
        triggers = [t.strip() for t in triggers_raw.split(",") if t.strip()]

        validation_raw = extract_section("Validation")
        validation = SkillLoader._parse_validation(validation_raw) if validation_raw else {}

        return Skill(
            name=name,
            filepath=filepath,
            triggers=triggers,
            description=extract_section("Description"),
            domain_knowledge=extract_section("Domain Knowledge"),
            examples=extract_section("Examples"),
            validation=validation,
        )

    @staticmethod
    def _parse_validation(text: str) -> dict:
        """Parse the ## Validation section into a dict of rules."""
        result = {}

        req_match = re.search(
            r"###\s*Required elements\s*\n(.*?)(?=^###|\Z)",
            text, re.MULTILINE | re.DOTALL)
        if req_match:
            result["required_elements"] = [
                line.lstrip("- ").strip()
                for line in req_match.group(1).strip().splitlines()
                if line.strip().startswith("-")
            ]

        cmd_match = re.search(
            r"###\s*Syntax check\s*\n.*?command:\s*(.+?)(?:\n|$)",
            text, re.MULTILINE | re.DOTALL)
        if cmd_match:
            result["syntax_check_command"] = cmd_match.group(1).strip()

        return result

    def detect_skills(self, trace_text: str) -> list[Skill]:
        """Auto-detect which system-specific skills match the trace."""
        return [skill for skill in self._skills.values()
                if skill.matches(trace_text)]

    def get_skill(self, name: str) -> Optional[Skill]:
        """Get a specific skill by filename stem (e.g. 'nextflow')."""
        return self._skills.get(name)

    def compose_prompt(self, trace_text: str,
                       skill_name: Optional[str] = None) -> str:
        """
        Compose the system prompt from base + detected/specified skills.

        Parameters
        ----------
        trace_text : str
            The workflow trace text, used for auto-detection.
        skill_name : str, optional
            Explicit skill name to use (bypasses auto-detection).

        Returns
        -------
        str
            Combined skill content to use as the system prompt.
        """
        sections = []

        # Always include base
        if self._base_skill:
            sections.append(self._base_skill.domain_knowledge)

        # Determine which system skill(s) to use
        system_skills = []
        if skill_name:
            skill = self.get_skill(skill_name)
            if skill:
                system_skills = [skill]
        else:
            system_skills = self.detect_skills(trace_text)

        for skill in system_skills:
            sections.append(
                f"=== SYSTEM-SPECIFIC KNOWLEDGE: {skill.name} ===\n"
                f"{skill.domain_knowledge}"
            )
            if skill.examples:
                sections.append(
                    f"=== SYSTEM-SPECIFIC EXAMPLES: {skill.name} ===\n"
                    f"{skill.examples}"
                )

        return "\n\n".join(sections)

    def available_skills(self) -> list[str]:
        """Return names of all available system-specific skills."""
        return sorted(self._skills.keys())

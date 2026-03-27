"""
Utilities for fetching WfFormat example instances from the WfInstances GitHub repository.

Shared by both the forward and backward LLM translators.

The WfInstances repo is organized as::

    <wms>/<app>/<app>-<platform>-<variant>-<run>.json

For example::

    makeflow/blast/blast-chameleon-large-001.json
    pegasus/montage/montage-chameleon-2.0-001.json
    nextflow/rnaseq-dirt02-001.json          (flat, no app subdir)

Users can refer to instances by:
- Application name: ``"blast"`` — searches across all WMS directories
- Specific filename: ``"blast-chameleon-small-001"`` — finds that exact instance
- Full repo path: ``"makeflow/blast"`` — used directly
"""

import logging
from typing import Dict, Any, List, Optional

import requests

logger = logging.getLogger(__name__)

WFINSTANCES_API = "https://api.github.com/repos/wfcommons/WfInstances/contents"

# Known WMS top-level directories in the WfInstances repo
_WMS_DIRS = ["pegasus", "makeflow", "nextflow", "helloworld"]


def resolve_instance_query(query: str,
                           ref: str = "main") -> List[Dict[str, Any]]:
    """
    Resolve a user-friendly instance query into fetched WfFormat instances.

    Accepts any of:
    - An application name (e.g. ``"blast"``, ``"montage"``)
    - A specific instance filename with or without ``.json`` (e.g. ``"blast-chameleon-small-001"``)
    - A full repo path (e.g. ``"makeflow/blast"``)

    :param query: The instance query string.
    :type query: str
    :param ref: Git reference. Defaults to "main".
    :type ref: str
    :return: List of instance dicts with 'url', 'filename', 'content'.
    :rtype: List[Dict[str, Any]]
    """
    query = query.strip()

    # If it contains a slash, treat as a direct repo path
    if "/" in query:
        return fetch_instances_from_path(query, ref=ref)

    # Try as a specific filename first (e.g. "blast-chameleon-small-001")
    filename = query if query.endswith(".json") else query + ".json"
    result = _find_file_by_name(filename, ref=ref)
    if result:
        return [result]

    # Otherwise treat as an application name — search across WMS dirs
    return _find_app(query, ref=ref)


def _find_file_by_name(filename: str,
                       ref: str = "main") -> Optional[Dict[str, Any]]:
    """
    Search for a specific JSON file across all WMS directories.

    Uses the GitHub search-like approach: list each WMS dir, look for
    subdirs matching the app prefix, then check for the exact file.
    """
    # Extract app name from filename (first segment before '-')
    app_name = filename.split("-")[0]

    for wms in _WMS_DIRS:
        # Try <wms>/<app>/<filename> (nested layout like makeflow/blast/)
        path = f"{wms}/{app_name}/{filename}"
        url = f"{WFINSTANCES_API}/{path}?ref={ref}"
        resp = requests.get(url)
        if resp.status_code == 200:
            data = resp.json()
            if isinstance(data, dict) and data.get("download_url"):
                raw = requests.get(data["download_url"]).text
                return {
                    "url": data["download_url"],
                    "filename": data["name"],
                    "content": raw,
                }

        # Try <wms>/<filename> (flat layout like nextflow/)
        path = f"{wms}/{filename}"
        url = f"{WFINSTANCES_API}/{path}?ref={ref}"
        resp = requests.get(url)
        if resp.status_code == 200:
            data = resp.json()
            if isinstance(data, dict) and data.get("download_url"):
                raw = requests.get(data["download_url"]).text
                return {
                    "url": data["download_url"],
                    "filename": data["name"],
                    "content": raw,
                }

    return None


def _find_app(app_name: str,
              ref: str = "main") -> List[Dict[str, Any]]:
    """
    Search for an application directory across all WMS top-level dirs.

    Tries ``<wms>/<app_name>`` for each known WMS. Also checks flat dirs
    (like nextflow/) for files whose name starts with ``app_name-``.
    """
    for wms in _WMS_DIRS:
        # Try nested: <wms>/<app_name>/
        instances = fetch_instances_from_path(f"{wms}/{app_name}", ref=ref)
        if instances:
            return instances

    # Try flat dirs — files named <app_name>-*.json directly under a WMS dir
    for wms in _WMS_DIRS:
        all_in_dir = fetch_instances_from_path(wms, ref=ref)
        matched = [ex for ex in all_in_dir
                   if ex["filename"].startswith(f"{app_name}-")]
        if matched:
            return matched

    logger.warning(f"Application '{app_name}' not found in WfInstances.")
    return []


def fetch_instances_from_path(path: str,
                              ref: str = "main") -> List[Dict[str, Any]]:
    """
    Fetch JSON files from a specific path in the WfInstances GitHub repository.

    :param path: Path within the repository (e.g. "makeflow/blast").
    :type path: str
    :param ref: Git reference (branch, tag, or commit SHA). Defaults to "main".
    :type ref: str
    :return: List of dicts with 'url', 'filename', and 'content' keys.
    :rtype: List[Dict[str, Any]]
    """
    logger.info(f"Fetching instances from WfInstances: {path} (ref={ref})")
    url = f"{WFINSTANCES_API}/{path}?ref={ref}"
    resp = requests.get(url)

    if resp.status_code == 404:
        return []

    if resp.status_code != 200:
        logger.warning(f"GitHub API returned status {resp.status_code} for path '{path}'.")
        return []

    listing = resp.json()

    if isinstance(listing, dict) and "message" in listing:
        logger.warning(f"GitHub API error for path '{path}': {listing['message']}.")
        return []

    if not isinstance(listing, list):
        return []

    examples = []
    for item in listing:
        if (item["type"] == "file"
                and item["name"].endswith(".json")
                and not item["name"].endswith(".md")):
            raw = requests.get(item["download_url"]).text
            examples.append({
                "url": item["download_url"],
                "filename": item["name"],
                "content": raw,
            })
    return examples


def load_instances(queries: List[str],
                   ref: str = "main") -> Dict[str, List[Dict[str, Any]]]:
    """
    Resolve multiple instance queries and fetch instances.

    Each query can be an application name, specific filename, or repo path.

    :param queries: List of instance queries (e.g. ["blast", "montage"]).
    :type queries: List[str]
    :param ref: Git reference. Defaults to "main".
    :type ref: str
    :return: Dict mapping each query to its list of fetched instances.
    :rtype: Dict[str, List[Dict[str, Any]]]
    """
    all_examples = {}
    for query in queries:
        all_examples[query] = resolve_instance_query(query, ref=ref)
    return all_examples


def retrieve_instances(queries: List[str],
                       num_examples: int = 3,
                       score_against: Optional[str] = None,
                       ref: str = "main") -> List[Dict[str, Any]]:
    """
    Fetch instances from WfInstances and return the top-k most relevant.

    Each query in ``queries`` can be:
    - An application name: ``"blast"``
    - A specific instance: ``"blast-chameleon-small-001"``
    - A full repo path: ``"makeflow/blast"``

    If ``score_against`` is provided, instances are ranked by word overlap
    with that text. Otherwise they are returned in fetch order, capped at
    ``num_examples``.

    :param queries: Instance queries to resolve and fetch.
    :type queries: List[str]
    :param num_examples: Max number of instances to return.
    :type num_examples: int
    :param score_against: Text to score similarity against (optional).
    :type score_against: str, optional
    :param ref: Git reference. Defaults to "main".
    :type ref: str
    :return: List of instance dicts with 'url', 'filename', 'content'.
    :rtype: List[Dict[str, Any]]
    """
    if not queries:
        return []

    by_query = load_instances(queries, ref=ref)
    flat = []
    for instances in by_query.values():
        flat.extend(instances)

    if not flat:
        logger.warning("No valid instances found from any of the provided queries.")
        return []

    if score_against:
        scored = [(_word_overlap(score_against, ex["content"]), ex) for ex in flat]
        scored.sort(reverse=True, key=lambda x: x[0])
        return [ex for _, ex in scored[:num_examples]]

    return flat[:num_examples]


def _word_overlap(a: str, b: str) -> float:
    """Naive similarity based on word overlap."""
    return len(set(a.split()) & set(b.split()))

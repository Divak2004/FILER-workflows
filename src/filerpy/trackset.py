"""
filerpy.trackset
~~~~~~~~~~~~~~~~
Read and write FILER Track Set manifest files (.trackset.json / .trackset.tsv).

A Track Set is a pinned, versioned list of FILER track identifiers plus the
metadata filters and genome build used to produce it.  It is the canonical
reproducibility artifact for any FILER-based analysis.

Minimal trackset.json schema::

    {
      "name": "immune_atac_hg38",
      "genome_build": "hg38",
      "created_at": "2025-02-24T00:00:00Z",
      "filters": {
        "assayType": "ATAC-seq",
        "tissueCategory": "blood"
      },
      "track_ids": ["ENCFF001ABC", "ENCFF002DEF"],
      "filer_release": null          // optional; fill if known
    }
"""
from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import pandas as pd

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Schema helpers
# ---------------------------------------------------------------------------

REQUIRED_FIELDS = {"genome_build", "track_ids"}

def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def make_trackset(
    df: pd.DataFrame,
    genome_build: str,
    filters: dict,
    name: str = "trackset",
    filer_release: Optional[str] = None,
) -> dict:
    """
    Build a Track Set manifest dict from a search_tracks() DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Output of ``search_tracks()``. Must contain an ``identifier`` column.
    genome_build : str
        Genome build used in the query.
    filters : dict
        The filter dict passed to ``search_tracks()``, for provenance.
    name : str
        Human-readable name for this track set.
    filer_release : str, optional
        FILER release / snapshot date, if known.

    Returns
    -------
    dict
        Manifest ready for ``save_trackset()``.

    Example
    -------
    >>> ts = make_trackset(df, "hg38", {"assayType": "ATAC-seq"}, name="atac_blood")
    """
    if "identifier" not in df.columns:
        raise ValueError("DataFrame must have an 'identifier' column.")

    return {
        "name":          name,
        "genome_build":  genome_build,
        "created_at":    _now_utc(),
        "filters":       filters,
        "track_count":   len(df),
        "track_ids":     df["identifier"].tolist(),
        "filer_release": filer_release,
    }


def save_trackset(
    trackset: dict,
    out_dir: str | Path = ".",
    stem: Optional[str] = None,
    write_tsv: bool = True,
    df: Optional[pd.DataFrame] = None,
) -> tuple[Path, Optional[Path]]:
    """
    Write a Track Set to disk as JSON (and optionally TSV).

    Parameters
    ----------
    trackset : dict
        Manifest produced by ``make_trackset()``.
    out_dir : str or Path
        Directory to write into (created if absent).
    stem : str, optional
        Base filename stem. Defaults to ``trackset["name"]``.
    write_tsv : bool
        If True and ``df`` is supplied, also write a TSV of the full
        track metadata alongside the JSON.
    df : pd.DataFrame, optional
        Full metadata DataFrame (from ``search_tracks()``). Required when
        ``write_tsv=True``.

    Returns
    -------
    (json_path, tsv_path)
        Paths to the written files. ``tsv_path`` is None if not written.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = stem or trackset.get("name", "trackset")

    json_path = out_dir / f"{stem}.trackset.json"
    json_path.write_text(json.dumps(trackset, indent=2))
    log.info("Wrote %s", json_path)

    tsv_path: Optional[Path] = None
    if write_tsv and df is not None:
        tsv_path = out_dir / f"{stem}.trackset.tsv"
        df.to_csv(tsv_path, sep="\t", index=False)
        log.info("Wrote %s", tsv_path)

    return json_path, tsv_path


def load_trackset(path: str | Path) -> dict:
    """
    Load a Track Set manifest from a JSON file and validate required fields.

    Parameters
    ----------
    path : str or Path
        Path to a ``*.trackset.json`` file.

    Returns
    -------
    dict
        Parsed manifest.

    Raises
    ------
    ValueError
        If required fields are missing.
    """
    path = Path(path)
    trackset = json.loads(path.read_text())
    missing = REQUIRED_FIELDS - set(trackset)
    if missing:
        raise ValueError(f"Track set {path} is missing required fields: {missing}")
    return trackset


def inspect_trackset(trackset: dict) -> pd.DataFrame:
    """
    Print a human-readable summary of a Track Set manifest.

    Returns a one-row summary DataFrame (useful in notebooks).
    """
    summary = {
        "name":          trackset.get("name", "—"),
        "genome_build":  trackset.get("genome_build", "—"),
        "created_at":    trackset.get("created_at", "—"),
        "track_count":   trackset.get("track_count", len(trackset.get("track_ids", []))),
        "filer_release": trackset.get("filer_release", "—"),
        "filters":       json.dumps(trackset.get("filters", {})),
    }
    df = pd.DataFrame([summary])
    print(df.to_string(index=False))
    return df
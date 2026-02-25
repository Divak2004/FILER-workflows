"""
filerpy.client
~~~~~~~~~~~~~~
Thin HTTP client for the FILER web API endpoints.

Typical usage::

    from filerpy.client import search_tracks, get_metadata_values

    df = search_tracks("hg38", assayType="ATAC-seq", cellType="CD14+ monocyte")
    print(df[["identifier", "assay", "cell_type", "tabix_file_url"]])
"""
from __future__ import annotations

import time
import logging
from typing import Optional

import requests
import pandas as pd

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BASE_URL  = "https://tf.lisanwanglab.org/FILER2"
ENDPOINTS = {
    "metadata":        f"{BASE_URL}/get_metadata.php",
    "region_overlaps": f"{BASE_URL}/get_overlapping_tracks_by_coord.php",
    "overlaps":        f"{BASE_URL}/get_overlaps.php",
    "data_region":     f"{BASE_URL}/get_data_region.php",
}

# Columns returned by get_metadata.php that we surface prominently
METADATA_STANDARD_COLS = [
    "identifier",
    "genome_build",
    "assay",
    "cell_type",
    "biosample_type",
    "tissue_category",
    "system_category",
    "life_stage",
    "data_source",
    "data_category",
    "classification",
    "output_type",
    "track_name",
    "processed_file_download_url",
    "tabix_file_url",
]

# Map friendly Python names → PHP GET parameter names for get_metadata.php
_METADATA_PARAM_MAP: dict[str, str] = {
    "assay":           "assayType",
    "cell_type":       "cellType",
    "tissue_category": "tissueCategory",
    "data_source":     "dataSource",
    "track_id":        "trackID",
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get(url: str, params: dict, timeout: int = 60, retries: int = 3) -> requests.Response:
    """GET with simple exponential-backoff retry on transient errors."""
    delay = 2
    last_exc: Optional[Exception] = None
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            r.raise_for_status()
            return r
        except requests.exceptions.RequestException as exc:
            last_exc = exc
            if attempt < retries - 1:
                log.warning("Request failed (%s), retrying in %ds…", exc, delay)
                time.sleep(delay)
                delay *= 2
    raise last_exc  # type: ignore[misc]


def _reorder(df: pd.DataFrame, standard_cols: list[str]) -> pd.DataFrame:
    """Put standard columns first; append any extras."""
    present = [c for c in standard_cols if c in df.columns]
    extras  = [c for c in df.columns   if c not in standard_cols]
    return df[present + extras]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def search_tracks(
    genome_build: str = "hg38",
    filter_string: Optional[str] = None,
    **kwargs: str,
) -> pd.DataFrame:
    """
    Query the FILER metadata endpoint and return matching tracks.

    Parameters
    ----------
    genome_build : str
        ``"hg38"`` (default) or ``"hg19"``.
    filter_string : str, optional
        Raw jq filter expression, e.g.::

            '.data_source == "ENCODE" and .assay == "ATAC-seq"'

        When supplied, named keyword arguments are ignored.
    **kwargs
        Named metadata filters using **PHP parameter names**:

        =================  ============================
        ``assayType``      e.g. ``"ATAC-seq"``
        ``cellType``       e.g. ``"CD14+ monocyte"``
        ``tissueCategory`` e.g. ``"blood"``
        ``dataSource``     e.g. ``"ENCODE"``
        ``trackID``        specific track identifier
        =================  ============================

        Alternatively you may use the friendly Python names
        (``assay``, ``cell_type``, ``tissue_category``, ``data_source``,
        ``track_id``) and they will be translated automatically.

    Returns
    -------
    pd.DataFrame
        One row per matching track. Standard columns appear first:
        ``identifier``, ``genome_build``, ``assay``, ``cell_type``,
        ``biosample_type``, ``tissue_category``, ``system_category``,
        ``life_stage``, ``data_source``, ``data_category``,
        ``classification``, ``output_type``, ``track_name``,
        ``processed_file_download_url``, ``tabix_file_url``.

    Examples
    --------
    >>> df = search_tracks("hg38", assayType="ATAC-seq", cellType="CD14+ monocyte")
    >>> df = search_tracks("hg38", filter_string='.data_source == "DASHR2"')
    """
    params: dict[str, str] = {
        "genomeBuild":  genome_build,
        "outputFormat": "json",
    }

    if filter_string:
        params["filterString"] = filter_string
    else:
        for key, val in kwargs.items():
            # Accept both friendly names and PHP names
            php_key = _METADATA_PARAM_MAP.get(key, key)
            params[php_key] = val

    log.debug("search_tracks params: %s", params)
    r = _get(ENDPOINTS["metadata"], params)

    data = r.json()
    if not data:
        log.info("search_tracks: no results for params=%s", params)
        return pd.DataFrame(columns=METADATA_STANDARD_COLS)

    df = pd.DataFrame(data)
    return _reorder(df, METADATA_STANDARD_COLS)


def get_metadata_values(
    field: str,
    query: Optional[str] = None,
    filters: Optional[str] = None,
    sort: str = "count",
    limit: int = 50,
    include_counts: bool = True,
) -> list[dict]:
    """
    List valid values for a metadata field (for building filter UIs / autocomplete).

    Parameters
    ----------
    field : str
        Metadata field name, e.g. ``"cell_type"``, ``"assay"``.
    query : str, optional
        Prefix/contains search within values, e.g. ``"mono"``.
    filters : str, optional
        Comma-separated context filters, e.g. ``"assay:ATAC-seq,tissue:blood"``.
    sort : str
        ``"count"`` (default) or ``"alpha"``.
    limit : int
        Maximum number of values to return (default 50).
    include_counts : bool
        If True, each entry includes a ``count`` field (default True).

    Returns
    -------
    list[dict]
        Each dict has at least ``value``; optionally ``count``.

    Examples
    --------
    >>> vals = get_metadata_values("cell_type", query="mono", filters="assay:ATAC-seq")
    >>> vals = get_metadata_values("assay", include_counts=True)
    """
    params: dict = {
        "field":          field,
        "sort":           sort,
        "limit":          limit,
        "include_counts": "true" if include_counts else "false",
    }
    if query:
        params["q"] = query
    if filters:
        params["filters"] = filters

    # NOTE: endpoint path TBC — update once confirmed
    url = f"{BASE_URL}/api/metadata/values"
    log.debug("get_metadata_values params: %s", params)
    r = _get(url, params)
    return r.json().get("values", [])
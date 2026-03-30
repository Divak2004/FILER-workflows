"""
filerpy.client
~~~~~~~~~~~~~~
Thin HTTP client for the FILER web API endpoints.

This module is responsible only for making HTTP requests and returning raw
data. All business logic (batching, flattening, filtering, error handling)
lives in the individual recipe scripts.

Typical usage::

    from filerpy.client import search_tracks, get_metadata, get_overlaps,
                               get_overlapping_tracks, get_metadata_values

    df  = search_tracks("hg38", assayType="ATAC-seq", cellType="CD14+ monocyte")
    raw = get_overlaps("chr19:44905791-44909393", ["NGBLPL2W2SM2WC"])
"""
from __future__ import annotations

import logging
import time
from typing import Optional

import pandas as pd
import requests

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
# Public API — metadata
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
        One row per matching track, with standard columns first.

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
            php_key = _METADATA_PARAM_MAP.get(key, key)
            params[php_key] = val

    log.debug("search_tracks params: %s", params)
    r = _get(ENDPOINTS["metadata"], params)

    data = r.json()
    if not data:
        log.info("search_tracks: no results for params=%s", params)
        return pd.DataFrame(columns=METADATA_STANDARD_COLS)

    return _reorder(pd.DataFrame(data), METADATA_STANDARD_COLS)


def get_metadata(track_id: str, genome_build: str = "hg38") -> dict:
    """
    Fetch metadata for a single track by ID.

    Parameters
    ----------
    track_id : str
        FILER track identifier, e.g. ``"NGBLPL2W2SM2WC"``.
    genome_build : str
        ``"hg38"`` (default) or ``"hg19"``.

    Returns
    -------
    dict
        Raw metadata record for the track.

    Examples
    --------
    >>> record = get_metadata("NGBLPL2W2SM2WC")
    """
    r = _get(ENDPOINTS["metadata"], {"genomeBuild": genome_build, "trackID": track_id})
    return r.json()[0]


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


# ---------------------------------------------------------------------------
# Public API — overlaps
# ---------------------------------------------------------------------------

def get_overlaps(region: str, track_ids: list[str]) -> list[dict]:
    """
    Fetch raw overlap records for a batch of track IDs and a region.

    This is a single GET request with no batching or flattening — that
    logic lives in ``filer_find_overlaps.py``.

    Parameters
    ----------
    region : str
        Genomic region in ``chrN:start-end`` format,
        e.g. ``"chr19:44905791-44909393"``.
    track_ids : list[str]
        Track identifiers for this batch.

    Returns
    -------
    list[dict]
        Raw response from the server. Each dict has ``Identifier``,
        ``queryRegion``, and ``features`` (a list of interval records).

    Examples
    --------
    >>> records = get_overlaps("chr19:44905791-44909393", ["NGBLPL2W2SM2WC"])
    """
    r = _get(
        ENDPOINTS["overlaps"],
        params={"region": region, "trackIDs": ",".join(track_ids)},
        timeout=300,
    )
    return r.json()


def get_overlapping_tracks(
    region: str,
    genome_build: str,
    **params,
) -> list[dict]:
    """
    Find all tracks in a genome build whose intervals overlap a region.

    This is a single GET request with no filtering or DataFrame logic —
    that lives in ``filer_coordinate_search.py``.

    Parameters
    ----------
    region : str
        Genomic region in ``chrN:start-end`` format,
        e.g. ``"chr1:100000-200000"``.
    genome_build : str
        ``"hg38"`` or ``"hg19"``.
    **params
        Additional GET parameters passed directly to the endpoint, e.g.
        ``filterString``, ``fullMetadata``, ``countOnly``.

    Returns
    -------
    list[dict]
        Raw response from the server, one dict per overlapping track.

    Examples
    --------
    >>> records = get_overlapping_tracks("chr1:100000-200000", "hg38", countOnly=1)
    """
    r = _get(
        ENDPOINTS["region_overlaps"],
        params={"region": region, "genomeBuild": genome_build, "outputFormat": "json", **params},
        timeout=300,
    )
    return r.json()


# ---------------------------------------------------------------------------
# Public convenience helpers — overlaps
# ---------------------------------------------------------------------------

HIT_SEP = "@@@"
DEFAULT_CHUNK_SIZE = 250


def load_ids_from_file(path: str, id_col: str) -> list[str]:
    """
    Load track IDs from a TSV or JSON file.

    This mirrors the CLI script's `--file/--id-col` behavior so notebooks can
    share the same inputs.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"File '{path}' not found.")

    if path.endswith(".json"):
        with open(path, "r") as f:
            data = json.load(f)
        df = pd.DataFrame(data)
    else:
        df = pd.read_csv(path, sep="\t")

    if id_col not in df.columns:
        available = ", ".join(df.columns.tolist())
        raise KeyError(f"Column '{id_col}' not found in '{path}'. Available columns: {available}")

    return df[id_col].dropna().unique().tolist()


def _feature_to_hit_string(feature: dict) -> str:
    """
    Serialize a feature dict into a @@@-delimited hit string.

    We rely on the API's field ordering by joining dict values in iteration
    order. Field sets/order can vary by track type, so we do not impose any
    fixed schema here.
    """
    return HIT_SEP.join("" if v is None else str(v) for v in feature.values())


def _normalize_overlaps_records(records: list[dict]) -> pd.DataFrame:
    """Flatten API overlap records into Identifier/queryRegion/hitString rows."""
    rows: list[dict] = []
    for rec in records:
        identifier = rec.get("Identifier", rec.get("identifier", ""))
        query_region = rec.get("queryRegion", rec.get("query_region", ""))

        features = rec.get("features")
        if features is None:
            # Flat record — treat everything except Identifier/queryRegion as a feature.
            feature = {
                k: v
                for k, v in rec.items()
                if k not in ("Identifier", "identifier", "queryRegion", "query_region")
            }
            features = [feature]

        for feat in features:
            rows.append(
                {
                    "Identifier": identifier,
                    "queryRegion": query_region,
                    "hitString": _feature_to_hit_string(feat),
                }
            )

    return pd.DataFrame(rows, columns=["Identifier", "queryRegion", "hitString"])


def _normalize_overlaps_response(data) -> pd.DataFrame:
    """Handle list vs wrapped payload shapes returned by the server."""
    if isinstance(data, list):
        return _normalize_overlaps_records(data)
    if isinstance(data, dict):
        for key in ("data", "results", "records", "overlaps"):
            if key in data and isinstance(data[key], list):
                return _normalize_overlaps_records(data[key])
        return _normalize_overlaps_records([data])
    raise TypeError(f"Unexpected overlaps response type: {type(data)}")


def fetch_overlaps(
    region: str,
    ids: list[str],
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> pd.DataFrame:
    """
    Fetch overlaps and return a flat DataFrame with a `hitString` column.

    Each feature dict in the API response is converted into `hitString` by
    joining its values with "@@@", preserving API-provided ordering.
    """
    if not ids:
        return pd.DataFrame(columns=["Identifier", "queryRegion", "hitString"])

    chunks = [ids[i : i + chunk_size] for i in range(0, len(ids), chunk_size)]
    frames: list[pd.DataFrame] = []

    for chunk in chunks:
        # Use the standard endpoint registry (like other APIs in this module)
        r = _get(
            ENDPOINTS["overlaps"],
            params={"region": region, "trackIDs": ",".join(chunk)},
            timeout=300,
        )
        data = r.json()
        df = _normalize_overlaps_response(data)
        if not df.empty:
            frames.append(df)

    return (
        pd.concat(frames, ignore_index=True)
        if frames
        else pd.DataFrame(columns=["Identifier", "queryRegion", "hitString"])
    )
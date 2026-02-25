"""
filerpy — thin Python client for the FILER functional genomics database.

Quick start::

    from filerpy import search_tracks

    df = search_tracks("hg38", assayType="ATAC-seq", cellType="CD14+ monocyte")
    print(df[["identifier", "assay", "cell_type", "tabix_file_url"]])
"""
from filerpy.client import search_tracks, get_metadata_values  # noqa: F401
from filerpy.trackset import make_trackset, save_trackset, load_trackset, inspect_trackset  # noqa: F401

__version__ = "0.1.0"
"""
Integration tests for Recipe 1 — Track Discovery via Metadata Filters.

These tests hit the real FILER endpoint at:
    https://tf.lisanwanglab.org/FILER2/get_metadata.php

They are skipped automatically when the endpoint is unreachable (e.g. in CI
without external network access).

Run manually:
    pytest tests/integration/test_recipe01_integration.py -v
"""
import pytest
import requests
import pandas as pd

from filerpy.client import search_tracks, ENDPOINTS


# ---------------------------------------------------------------------------
# Skip entire module if FILER endpoint is unreachable
# ---------------------------------------------------------------------------

def _filer_reachable() -> bool:
    try:
        r = requests.get(
            ENDPOINTS["metadata"],
            params={"genomeBuild": "hg38", "trackID": "NGBLPL2W2SM2WC"},
            timeout=10,
        )
        return r.status_code == 200 and len(r.json()) > 0
    except Exception:
        return False


pytestmark = pytest.mark.skipif(
    not _filer_reachable(),
    reason="FILER endpoint not reachable — skipping integration tests",
)


# ---------------------------------------------------------------------------
# 1. Single known track lookup
# ---------------------------------------------------------------------------

class TestKnownTrackLookup:
    """Verify a specific known track returns the expected metadata."""

    KNOWN_ID = "NGBLPL2W2SM2WC"

    def test_known_track_returns_one_result(self):
        df = search_tracks("hg38", trackID=self.KNOWN_ID)
        assert len(df) == 1, f"Expected 1 result for trackID={self.KNOWN_ID}"

    def test_known_track_identifier_correct(self):
        df = search_tracks("hg38", trackID=self.KNOWN_ID)
        assert df.iloc[0]["identifier"] == self.KNOWN_ID

    def test_known_track_genome_build_correct(self):
        df = search_tracks("hg38", trackID=self.KNOWN_ID)
        assert df.iloc[0]["genome_build"] == "hg38"

    def test_known_track_has_valid_download_url(self):
        df = search_tracks("hg38", trackID=self.KNOWN_ID)
        url = df.iloc[0]["processed_file_download_url"]
        assert url.startswith("https://tf.lisanwanglab.org"), \
            f"Unexpected download URL: {url}"

    def test_known_track_has_valid_tabix_url(self):
        df = search_tracks("hg38", trackID=self.KNOWN_ID)
        url = df.iloc[0]["tabix_file_url"]
        assert url.endswith(".tbi"), f"tabix_file_url should end with .tbi: {url}"


# ---------------------------------------------------------------------------
# 2. Metadata filter queries
# ---------------------------------------------------------------------------

class TestMetadataFilters:
    """Verify that filter parameters return correctly scoped results."""

    def test_assay_filter_returns_results(self):
        df = search_tracks("hg38", assayType="ATAC-seq")
        assert len(df) > 0

    def test_assay_filter_respected(self):
        df = search_tracks("hg38", assayType="ATAC-seq")
        assert (df["assay"] == "ATAC-seq").all(), \
            "All returned tracks should have assay == ATAC-seq"

    def test_cell_type_filter_respected(self):
        df = search_tracks("hg38", assayType="ATAC-seq", cellType="CD14+ monocyte")
        if not df.empty:
            assert (df["cell_type"] == "CD14+ monocyte").all(), \
                "All returned tracks should have cell_type == CD14+ monocyte"

    def test_genome_build_hg19_works(self):
        df = search_tracks("hg19", assayType="ATAC-seq")
        assert isinstance(df, pd.DataFrame)
        if not df.empty:
            assert (df["genome_build"] == "hg19").all()

    def test_tissue_category_filter_respected(self):
        df = search_tracks("hg38", tissueCategory="Blood")
        if not df.empty:
            assert (df["tissue_category"] == "Blood").all()

    def test_data_source_filter_respected(self):
        df = search_tracks("hg38", dataSource="ENCODE")
        if not df.empty:
            assert (df["data_source"] == "ENCODE").all()


# ---------------------------------------------------------------------------
# 3. Response schema
# ---------------------------------------------------------------------------

class TestResponseSchema:
    """Verify the API returns the expected columns."""

    EXPECTED_COLS = [
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

    def test_all_expected_columns_present(self):
        df = search_tracks("hg38", trackID="NGBLPL2W2SM2WC")
        for col in self.EXPECTED_COLS:
            assert col in df.columns, f"Missing expected column: {col}"

    def test_identifiers_are_unique(self):
        df = search_tracks("hg38", assayType="ATAC-seq", tissueCategory="Blood")
        assert df["identifier"].nunique() == len(df), \
            "Duplicate identifiers found in results"


# ---------------------------------------------------------------------------
# 4. filterString equivalence
# ---------------------------------------------------------------------------

class TestFilterStringEquivalence:
    """Verify that named params and raw filterString produce the same results."""

    def test_assay_filter_equivalent(self):
        df_named = search_tracks("hg38", assayType="ATAC-seq", cellType="CD14+ monocyte")
        df_filter = search_tracks(
            "hg38",
            filter_string='.assay == "ATAC-seq" and .cell_type == "CD14+ monocyte"',
        )
        ids_named  = sorted(df_named["identifier"].tolist())
        ids_filter = sorted(df_filter["identifier"].tolist())
        assert ids_named == ids_filter, \
            "Named params and filterString returned different track sets"
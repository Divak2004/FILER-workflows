"""
Unit tests for Recipe 1 — Track Discovery via Metadata Filters.

Tests cover:
- search_tracks() response parsing and column ordering
- Parameter translation (friendly names → PHP param names)
- filterString override behaviour
- Empty response handling
- make_trackset() manifest structure
- save_trackset() / load_trackset() round-trip
- inspect_trackset() output

No network calls are made — all HTTP responses are mocked.

Run:
    pytest tests/unit/test_recipe01_unit.py -v
"""
import json
import pytest
from unittest.mock import patch, MagicMock

import pandas as pd

from filerpy.client import (
    search_tracks,
    ENDPOINTS,
    METADATA_STANDARD_COLS,
)
from filerpy.trackset import (
    make_trackset,
    save_trackset,
    load_trackset,
    inspect_trackset,
)


# ---------------------------------------------------------------------------
# Shared fixture — one real-schema track record + one synthetic
# ---------------------------------------------------------------------------

MOCK_TRACKS = [
    {
        "identifier":                  "NGBLPL2W2SM2WC",
        "data_source":                 "Blueprint",
        "file_name":                   "57418.Blueprint.ERS792033.WGB-Seq.peak_calls.bed.gz",
        "number_of_intervals":         738551,
        "bp_covered":                  1857857024,
        "output_type":                 "peaks",
        "genome_build":                "hg38",
        "cell_type":                   "B cell",
        "biosample_type":              "Primary cell",
        "biosamples_term_id":          "CL_2000006",
        "tissue_category":             "Blood",
        "system_category":             "Cardiovascular",
        "encode_experiment_id":        "Not applicable",
        "biological_replicate":        "Not applicable",
        "technical_replicate":         "Not applicable",
        "antibody":                    "Not applicable",
        "assay":                       "WGB-Seq",
        "file_format":                 "bed bed4",
        "file_size":                   8541765,
        "downloaded_date":             "03/12/2021",
        "release_date":                "08/31/2016",
        "date_added_to_filer":         "04/10/2021",
        "processed_file_download_url": "https://tf.lisanwanglab.org/GADB/FILER2/Annotationtracks/IHEC/IHEC_Blueprint/WGB-Seq/bed4/hg38/57418.Blueprint.ERS792033.WGB-Seq.peak_calls.bed.gz",
        "processed_file_md5":          "30bf3adf841fd5c17a1832fef34f6c76",
        "link_out_url":                "https://ihec-epigenomes.org/",
        "data_category":               "Methylation",
        "classification":              "WGB-Seq peaks",
        "track_description":           "Original_cell_type=Germinal center B cell;ID=ERS792033",
        "life_stage":                  "Child",
        "track_name":                  "Blueprint B cell WGB-Seq peaks (bed4) [Life stage: Child]",
        "tabix_file_url":              "https://tf.lisanwanglab.org/GADB/FILER2/Annotationtracks/IHEC/IHEC_Blueprint/WGB-Seq/bed4/hg38/57418.Blueprint.ERS792033.WGB-Seq.peak_calls.bed.gz.tbi",
    },
    {
        "identifier":                  "ENCFF002DEF",
        "data_source":                 "ENCODE",
        "file_name":                   "ENCFF002DEF.bed.gz",
        "number_of_intervals":         120000,
        "bp_covered":                  500000000,
        "output_type":                 "peaks",
        "genome_build":                "hg38",
        "cell_type":                   "CD14+ monocyte",
        "biosample_type":              "Primary cell",
        "biosamples_term_id":          "CL_0002057",
        "tissue_category":             "Blood",
        "system_category":             "Immune",
        "encode_experiment_id":        "ENCSR000ABC",
        "biological_replicate":        "1",
        "technical_replicate":         "1",
        "antibody":                    "Not applicable",
        "assay":                       "ATAC-seq",
        "file_format":                 "bed narrowPeak",
        "file_size":                   2000000,
        "downloaded_date":             "01/01/2024",
        "release_date":                "01/01/2023",
        "date_added_to_filer":         "02/01/2024",
        "processed_file_download_url": "https://tf.lisanwanglab.org/GADB/FILER2/ENCODE/ATAC-seq/hg38/ENCFF002DEF.bed.gz",
        "processed_file_md5":          "abc123",
        "link_out_url":                "https://www.encodeproject.org/",
        "data_category":               "Chromatin accessibility",
        "classification":              "ATAC-seq peaks",
        "track_description":           "CD14+ monocyte ATAC-seq peaks",
        "life_stage":                  "Adult",
        "track_name":                  "ENCODE CD14+ monocyte ATAC-seq peaks [Life stage: Adult]",
        "tabix_file_url":              "https://tf.lisanwanglab.org/GADB/FILER2/ENCODE/ATAC-seq/hg38/ENCFF002DEF.bed.gz.tbi",
    },
]


def _mock_response(data):
    mock = MagicMock()
    mock.json.return_value = data
    mock.raise_for_status.return_value = None
    return mock


# ---------------------------------------------------------------------------
# search_tracks()
# ---------------------------------------------------------------------------

class TestSearchTracks:

    @patch("filerpy.client.requests.get")
    def test_returns_dataframe(self, mock_get):
        mock_get.return_value = _mock_response(MOCK_TRACKS)
        df = search_tracks("hg38", assayType="ATAC-seq")
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 2

    @patch("filerpy.client.requests.get")
    def test_standard_columns_present(self, mock_get):
        mock_get.return_value = _mock_response(MOCK_TRACKS)
        df = search_tracks("hg38", assayType="ATAC-seq")
        for col in [
            "identifier", "assay", "cell_type", "biosample_type",
            "tissue_category", "data_source",
            "processed_file_download_url", "tabix_file_url",
        ]:
            assert col in df.columns, f"Missing column: {col}"

    @patch("filerpy.client.requests.get")
    def test_standard_columns_come_first(self, mock_get):
        mock_get.return_value = _mock_response(MOCK_TRACKS)
        df = search_tracks("hg38", assayType="ATAC-seq")
        for i, col in enumerate(["identifier", "genome_build", "assay"]):
            assert df.columns[i] == col, \
                f"Expected column {i} to be '{col}', got '{df.columns[i]}'"

    @patch("filerpy.client.requests.get")
    def test_correct_endpoint_called(self, mock_get):
        mock_get.return_value = _mock_response(MOCK_TRACKS)
        search_tracks("hg38", assayType="ATAC-seq")
        called_url = mock_get.call_args[0][0]
        assert called_url == ENDPOINTS["metadata"]

    @patch("filerpy.client.requests.get")
    def test_genome_build_passed_as_genomeBuild(self, mock_get):
        mock_get.return_value = _mock_response(MOCK_TRACKS)
        search_tracks("hg19", assayType="ATAC-seq")
        params = mock_get.call_args[1]["params"]
        assert params["genomeBuild"] == "hg19"

    @patch("filerpy.client.requests.get")
    def test_output_format_always_json(self, mock_get):
        mock_get.return_value = _mock_response(MOCK_TRACKS)
        search_tracks("hg38", assayType="ATAC-seq")
        params = mock_get.call_args[1]["params"]
        assert params["outputFormat"] == "json"

    @patch("filerpy.client.requests.get")
    def test_friendly_name_assay_translated(self, mock_get):
        mock_get.return_value = _mock_response(MOCK_TRACKS)
        search_tracks("hg38", assay="ATAC-seq")
        params = mock_get.call_args[1]["params"]
        assert "assayType" in params
        assert params["assayType"] == "ATAC-seq"

    @patch("filerpy.client.requests.get")
    def test_friendly_name_cell_type_translated(self, mock_get):
        mock_get.return_value = _mock_response(MOCK_TRACKS)
        search_tracks("hg38", cell_type="CD14+ monocyte")
        params = mock_get.call_args[1]["params"]
        assert "cellType" in params
        assert params["cellType"] == "CD14+ monocyte"

    @patch("filerpy.client.requests.get")
    def test_friendly_name_data_source_translated(self, mock_get):
        mock_get.return_value = _mock_response(MOCK_TRACKS)
        search_tracks("hg38", data_source="ENCODE")
        params = mock_get.call_args[1]["params"]
        assert "dataSource" in params

    @patch("filerpy.client.requests.get")
    def test_filter_string_sent_as_filterString(self, mock_get):
        mock_get.return_value = _mock_response(MOCK_TRACKS)
        search_tracks("hg38", filter_string='.data_source == "ENCODE"')
        params = mock_get.call_args[1]["params"]
        assert "filterString" in params
        assert params["filterString"] == '.data_source == "ENCODE"'

    @patch("filerpy.client.requests.get")
    def test_filter_string_overrides_named_kwargs(self, mock_get):
        """When filterString is provided, named kwargs should not appear in params."""
        mock_get.return_value = _mock_response(MOCK_TRACKS)
        search_tracks("hg38", filter_string='.data_source == "ENCODE"', assayType="ATAC-seq")
        params = mock_get.call_args[1]["params"]
        assert "filterString" in params
        assert "assayType" not in params

    @patch("filerpy.client.requests.get")
    def test_empty_response_returns_empty_dataframe(self, mock_get):
        mock_get.return_value = _mock_response([])
        df = search_tracks("hg38", assayType="ATAC-seq")
        assert isinstance(df, pd.DataFrame)
        assert df.empty
        assert "identifier" in df.columns

    @patch("filerpy.client.requests.get")
    def test_extra_columns_appended_after_standard(self, mock_get):
        """Non-standard columns like file_name should appear after standard cols."""
        mock_get.return_value = _mock_response(MOCK_TRACKS)
        df = search_tracks("hg38", assayType="ATAC-seq")
        std_indices = [df.columns.get_loc(c) for c in METADATA_STANDARD_COLS if c in df.columns]
        extra_indices = [df.columns.get_loc(c) for c in ["file_name", "file_size"] if c in df.columns]
        if std_indices and extra_indices:
            assert max(std_indices) < min(extra_indices)


# ---------------------------------------------------------------------------
# make_trackset()
# ---------------------------------------------------------------------------

class TestMakeTrackset:

    def test_required_fields_present(self):
        df = pd.DataFrame(MOCK_TRACKS)
        ts = make_trackset(df, "hg38", {"assayType": "ATAC-seq"}, name="test_ts")
        for field in ["name", "genome_build", "created_at", "filters", "track_count", "track_ids"]:
            assert field in ts, f"Missing field: {field}"

    def test_genome_build_correct(self):
        df = pd.DataFrame(MOCK_TRACKS)
        ts = make_trackset(df, "hg38", {"assayType": "ATAC-seq"})
        assert ts["genome_build"] == "hg38"

    def test_track_count_matches_dataframe(self):
        df = pd.DataFrame(MOCK_TRACKS)
        ts = make_trackset(df, "hg38", {})
        assert ts["track_count"] == len(df)

    def test_track_ids_match_identifiers(self):
        df = pd.DataFrame(MOCK_TRACKS)
        ts = make_trackset(df, "hg38", {})
        assert ts["track_ids"] == ["NGBLPL2W2SM2WC", "ENCFF002DEF"]

    def test_filters_stored(self):
        df = pd.DataFrame(MOCK_TRACKS)
        filters = {"assayType": "ATAC-seq", "tissueCategory": "Blood"}
        ts = make_trackset(df, "hg38", filters)
        assert ts["filters"] == filters

    def test_created_at_is_iso_string(self):
        df = pd.DataFrame(MOCK_TRACKS)
        ts = make_trackset(df, "hg38", {})
        assert "T" in ts["created_at"]  # basic ISO 8601 check
        assert ts["created_at"].endswith("Z")

    def test_missing_identifier_column_raises(self):
        df = pd.DataFrame([{"not_id": "x"}])
        with pytest.raises(ValueError, match="identifier"):
            make_trackset(df, "hg38", {})


# ---------------------------------------------------------------------------
# save_trackset() / load_trackset()
# ---------------------------------------------------------------------------

class TestSaveLoadTrackset:

    def test_json_file_created(self, tmp_path):
        df = pd.DataFrame(MOCK_TRACKS)
        ts = make_trackset(df, "hg38", {}, name="test_ts")
        json_path, _ = save_trackset(ts, out_dir=tmp_path, write_tsv=False)
        assert json_path.exists()
        assert json_path.suffix == ".json"

    def test_tsv_file_created_when_requested(self, tmp_path):
        df = pd.DataFrame(MOCK_TRACKS)
        ts = make_trackset(df, "hg38", {}, name="test_ts")
        _, tsv_path = save_trackset(ts, out_dir=tmp_path, write_tsv=True, df=df)
        assert tsv_path is not None
        assert tsv_path.exists()

    def test_round_trip_preserves_genome_build(self, tmp_path):
        df = pd.DataFrame(MOCK_TRACKS)
        ts = make_trackset(df, "hg38", {}, name="test_ts")
        json_path, _ = save_trackset(ts, out_dir=tmp_path, write_tsv=False)
        loaded = load_trackset(json_path)
        assert loaded["genome_build"] == "hg38"

    def test_round_trip_preserves_track_ids(self, tmp_path):
        df = pd.DataFrame(MOCK_TRACKS)
        ts = make_trackset(df, "hg38", {}, name="test_ts")
        json_path, _ = save_trackset(ts, out_dir=tmp_path, write_tsv=False)
        loaded = load_trackset(json_path)
        assert loaded["track_ids"] == ["NGBLPL2W2SM2WC", "ENCFF002DEF"]

    def test_tsv_contains_correct_identifiers(self, tmp_path):
        df = pd.DataFrame(MOCK_TRACKS)
        ts = make_trackset(df, "hg38", {}, name="test_ts")
        _, tsv_path = save_trackset(ts, out_dir=tmp_path, write_tsv=True, df=df)
        loaded_df = pd.read_csv(tsv_path, sep="\t")
        assert list(loaded_df["identifier"]) == ["NGBLPL2W2SM2WC", "ENCFF002DEF"]

    def test_load_missing_required_field_raises(self, tmp_path):
        bad_ts = {"name": "bad"}  # missing genome_build and track_ids
        path = tmp_path / "bad.trackset.json"
        path.write_text(json.dumps(bad_ts))
        with pytest.raises(ValueError, match="missing required fields"):
            load_trackset(path)


# ---------------------------------------------------------------------------
# inspect_trackset()
# ---------------------------------------------------------------------------

class TestInspectTrackset:

    def test_returns_one_row_dataframe(self, capsys):
        df = pd.DataFrame(MOCK_TRACKS)
        ts = make_trackset(df, "hg38", {"assayType": "ATAC-seq"}, name="test_ts")
        result = inspect_trackset(ts)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 1

    def test_output_contains_genome_build(self, capsys):
        df = pd.DataFrame(MOCK_TRACKS)
        ts = make_trackset(df, "hg38", {}, name="test_ts")
        inspect_trackset(ts)
        captured = capsys.readouterr()
        assert "hg38" in captured.out

    def test_output_contains_track_count(self, capsys):
        df = pd.DataFrame(MOCK_TRACKS)
        ts = make_trackset(df, "hg38", {}, name="test_ts")
        inspect_trackset(ts)
        captured = capsys.readouterr()
        assert "2" in captured.out
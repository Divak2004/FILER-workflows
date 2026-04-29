#!/usr/bin/env python3
"""
filer_install.py — Unified CLI for the FILER genomic annotation repository.

Subcommands
-----------
  search    Query the FILER API for tracks matching metadata filters.
  install   Download and index annotation tracks.

              Two modes:
                (a) Template mode  — wraps install_filer.sh with a metadata template URL/file.
                (b) From-tracks mode — reads a Recipe 3 or Recipe 10 TSV produced by
                    filer_filter_then_overlaps.py and downloads + indexes only those tracks directly
                    in Python (no shell script required).

              In from-tracks mode a separate track-metadata TSV is written
              alongside the downloads (default: <target-dir>/track_metadata.tsv).

Typical workflows
-----------------
  # Template mode (bulk install everything in a metadata template)
  python filer_install.py install FILER_data metadata.hg38.template filer.ini

  # From-tracks mode (install a biologically filtered subset)
  python filer_filter_then_overlaps.py --genome-build hg38 --region "chr1:100000-200000" \\
      --assay ATAC-seq --tissue-category Blood --out results.tsv
  python filer_install.py install --from-tracks results.tsv --target-dir FILER_data \\
      --giggle /path/to/giggle --tabix /path/to/tabix

Run `python filer_install.py <subcommand> --help` for full option details.
"""

import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path


# ──────────────────────────────────────────────────────────────────────────────
# ANSI helpers
# ──────────────────────────────────────────────────────────────────────────────

USE_COLOR = sys.stdout.isatty()

def _c(code, text):  return f"\033[{code}m{text}\033[0m" if USE_COLOR else text
def green(t):        return _c("32", t)
def yellow(t):       return _c("33", t)
def red(t):          return _c("31", t)
def cyan(t):         return _c("36", t)
def bold(t):         return _c("1",  t)
def dim(t):          return _c("2",  t)
def magenta(t):      return _c("35", t)


# ──────────────────────────────────────────────────────────────────────────────
# Shared utilities
# ──────────────────────────────────────────────────────────────────────────────

def print_header(title: str) -> None:
    width = 52
    print(bold(f"\n╔{'═' * width}╗"))
    print(bold(f"║  {title:<{width - 2}}║"))
    print(bold(f"╚{'═' * width}╝\n"))


def print_kv(label: str, value: str) -> None:
    print(f"  {dim(f'{label:<16}')} {cyan(value)}")

def format_bytes(n):
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} PB"

def estimate_and_confirm(tracks: list[dict], target_dir: str) -> bool:
    """
    Deduplicate on 'identifier', sum file_size, check free disk space,
    then prompt the user to confirm before downloading.
    Returns True if the user confirms, False otherwise.
    """
    seen = {}
    for row in tracks:
        ident = row.get("identifier")
        if ident and ident not in seen:
            seen[ident] = row

    total_bytes = 0
    missing_size = 0
    for row in seen.values():
        raw = row.get("file_size", "").strip()
        if raw.isdigit():
            total_bytes += int(raw)
        else:
            missing_size += 1

    # Free space on the target volume
    Path(target_dir).mkdir(parents=True, exist_ok=True)
    free_bytes = shutil.disk_usage(target_dir).free

    print("\n── Storage estimate ─────────────────────────────")
    print(f"  Unique tracks to download : {len(seen)}")
    if missing_size:
        print(f"  Tracks with unknown size  : {missing_size} (not counted)")
    print(f"  Estimated download size   : {format_bytes(total_bytes)}")
    print(f"  Free space on target      : {format_bytes(free_bytes)}")

    if total_bytes > free_bytes:
        print("\n  ⚠️  WARNING: estimated size exceeds available disk space.")

    print("─────────────────────────────────────────────────\n")

    answer = input("Proceed with installation? [y/N] ").strip().lower()
    return answer in ("y", "yes")

# ──────────────────────────────────────────────────────────────────────────────
# `filer search` — query FILER metadata API
# ──────────────────────────────────────────────────────────────────────────────

FILER_BASE_URL = "https://tf.lisanwanglab.org/FILER2"
METADATA_ENDPOINT = f"{FILER_BASE_URL}/get_metadata.php"

# Columns included in the default TSV output (subset of full API response)
DEFAULT_TSV_COLUMNS = [
    "identifier",
    "genome_build",
    "assay",
    "cell_type",
    "biosample_type",
    "tissue_category",
    "life_stage",
    "data_source",
    "data_category",
    "classification",
    "track_name",
    "number_of_intervals",
    "file_size",
    "processed_file_download_url",
    "tabix_file_url",
]


def cmd_search(args: argparse.Namespace) -> int:
    """Execute the `filer search` subcommand."""

    print_header("FILER Search")
    print_kv("Genome build", args.genome_build)
    if args.assay:          print_kv("Assay",          args.assay)
    if args.cell_type:      print_kv("Cell type",      args.cell_type)
    if args.tissue_category:print_kv("Tissue category",args.tissue_category)
    if args.data_source:    print_kv("Data source",    args.data_source)
    if args.track_id:       print_kv("Track ID",       args.track_id)
    if args.filter_string:  print_kv("Filter string",  args.filter_string)
    print()

    # Build query parameters
    params: dict[str, str] = {
        "genomeBuild":  args.genome_build,
        "outputFormat": "json",
    }
    if args.assay:            params["assayType"]       = args.assay
    if args.cell_type:        params["cellType"]        = args.cell_type
    if args.tissue_category:  params["tissueCategory"]  = args.tissue_category
    if args.data_source:      params["dataSource"]      = args.data_source
    if args.track_id:         params["trackID"]         = args.track_id
    if args.filter_string:    params["filterString"]    = args.filter_string

    query_string = urllib.parse.urlencode(params)
    url = f"{METADATA_ENDPOINT}?{query_string}"

    # Fetch
    print(f"  {cyan('→')} Querying FILER API …", flush=True)
    t0 = time.time()
    try:
        with urllib.request.urlopen(url, timeout=60) as resp:
            raw = resp.read().decode("utf-8")
    except Exception as exc:
        print(red(f"\n  ✗  API request failed: {exc}"))
        return 1
    elapsed = time.time() - t0

    # Parse response
    try:
        data = json.loads(raw)
    except json.JSONDecodeError as exc:
        print(red(f"\n  ✗  Could not parse API response as JSON: {exc}"))
        print(dim(f"     First 200 chars: {raw[:200]}"))
        return 1

    if not isinstance(data, list):
        # API may return an error object
        print(red(f"\n  ✗  Unexpected API response: {data}"))
        return 1

    n = len(data)
    print(f"  {green('✓')} {bold(str(n))} track(s) found  {dim(f'({elapsed:.1f}s)')}\n")

    if n == 0:
        print(yellow("  No tracks matched your filters. Try broadening the search."))
        return 0

    # Determine output columns
    columns = DEFAULT_TSV_COLUMNS if not args.columns else [c.strip() for c in args.columns.split(",")]

    # Resolve output path
    out_path = Path(args.out) if args.out else None
    if out_path:
        out_path.parent.mkdir(parents=True, exist_ok=True)

    # Write TSV
    tsv_path = out_path or Path("tracks.tsv")
    _write_tsv(data, columns, tsv_path)
    print(f"  {green('✓')} TSV saved → {cyan(str(tsv_path))}")

    # Optionally write JSON
    if args.json:
        json_path = tsv_path.with_suffix(".json")
        json_path.write_text(json.dumps(data, indent=2))
        print(f"  {green('✓')} JSON saved → {cyan(str(json_path))}")

    # Print a preview table
    if not args.quiet:
        _print_preview(data, columns, max_rows=args.preview_rows)

    # Summary stats
    print()
    _print_summary(data)

    return 0


def _write_tsv(records: list[dict], columns: list[str], path: Path) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def _print_preview(records: list[dict], columns: list[str], max_rows: int = 10) -> None:
    """Print a compact aligned preview of the first N records."""
    preview_cols = [c for c in ["identifier", "assay", "cell_type", "tissue_category", "data_source", "track_name"]
                    if c in columns]
    if not preview_cols:
        preview_cols = columns[:5]

    # Compute column widths
    widths = {c: len(c) for c in preview_cols}
    rows_to_show = records[:max_rows]
    for rec in rows_to_show:
        for c in preview_cols:
            widths[c] = max(widths[c], len(str(rec.get(c, ""))))

    max_width = 28
    widths = {c: min(w, max_width) for c, w in widths.items()}

    def fmt(val, w):
        s = str(val) if val is not None else ""
        return s[:w].ljust(w) if len(s) <= w else s[:w - 1] + "…"

    sep = "  " + "  ".join("─" * widths[c] for c in preview_cols)
    header = "  " + "  ".join(bold(fmt(c, widths[c])) for c in preview_cols)

    print(f"\n  {dim(f'Preview (first {min(max_rows, len(records))} of {len(records)} tracks):')}")
    print(sep)
    print(header)
    print(sep)
    for rec in rows_to_show:
        row = "  " + "  ".join(fmt(rec.get(c, ""), widths[c]) for c in preview_cols)
        print(row)
    if len(records) > max_rows:
        print(f"  {dim(f'… and {len(records) - max_rows} more rows in the output file.')}")
    print(sep)


def _print_summary(records: list[dict]) -> None:
    """Print aggregate breakdowns by assay and data_source."""
    def count_by(field):
        counts: dict[str, int] = {}
        for r in records:
            v = r.get(field) or "unknown"
            counts[v] = counts.get(v, 0) + 1
        return sorted(counts.items(), key=lambda x: -x[1])

    print(f"  {bold('Breakdown by assay:')}")
    for val, cnt in count_by("assay")[:8]:
        bar = "█" * min(cnt, 30)
        print(f"    {fmt_left(val, 28)} {cyan(bar)} {cnt}")

    print(f"\n  {bold('Breakdown by data_source:')}")
    for val, cnt in count_by("data_source")[:8]:
        bar = "█" * min(cnt, 30)
        print(f"    {fmt_left(val, 28)} {magenta(bar)} {cnt}")
    print()


def fmt_left(s, w):
    s = str(s)
    return s[:w].ljust(w) if len(s) <= w else s[:w - 1] + "…"


# ──────────────────────────────────────────────────────────────────────────────
# `filer install` — shell script wrapper (from filer_install.py)
# ──────────────────────────────────────────────────────────────────────────────

class DownloadTracker:
    _RE_DOWNLOADING = re.compile(r"Dowloading (.+?) \[(\d+)/(\d+)\]")
    _RE_WARNING_SIZE = re.compile(r"\*\*\*WARNING: file size mismatch")
    _RE_WARNING_RETRY = re.compile(r"removing existing file and re-downloading")
    _RE_ERROR_MD5    = re.compile(r"\*\*\*ERROR: Downloading")
    _RE_DONE         = re.compile(r"Dowloading completed")

    def __init__(self):
        self.current = self.total = self.warnings = self.errors = 0
        self.current_file = ""

    def feed(self, line: str) -> bool:
        m = self._RE_DOWNLOADING.search(line)
        if m:
            self.current_file = m.group(1).split("/")[-1]
            self.current, self.total = int(m.group(2)), int(m.group(3))
            pct = int(100 * self.current / self.total) if self.total else 0
            filled = int(28 * pct / 100)
            bar_str = "█" * filled + "░" * (28 - filled)
            bar = green(bar_str) if pct == 100 else cyan(bar_str)
            print(f"\r  {cyan('↓')} [{bar}] {pct:3d}%  {dim(self.current_file[:50]):<52}", end="", flush=True)
            return True
        if self._RE_WARNING_SIZE.search(line) or self._RE_WARNING_RETRY.search(line):
            self.warnings += 1
            print(f"\n  {yellow('⚠')}  {line.strip()}")
            return True
        if self._RE_ERROR_MD5.search(line):
            self.errors += 1
            print(f"\n  {red('✗')}  {line.strip()}")
            return True
        if self._RE_DONE.search(line):
            print(f"\r  {green('✓')} Download complete — "
                  f"{self.total} track(s)  "
                  f"warnings={yellow(str(self.warnings))}  "
                  f"errors={red(str(self.errors)) if self.errors else green('0')}")
            return True
        return False


class IndexTracker:
    _RE_GIGGLE_START = re.compile(r"Starting Giggle indexing")
    _RE_GIGGLE_DIR   = re.compile(r"Indexing (.+?) \[(\d+)/(\d+)\]")
    _RE_GIGGLE_SKIP  = re.compile(r"\*\*\*WARNING: SKIPPING directory (.+?)\.")
    _RE_TABIX_START  = re.compile(r"Creating tabix index")

    def __init__(self):
        self.phase = None
        self.warnings = 0

    def feed(self, line: str) -> bool:
        if self._RE_GIGGLE_START.search(line):
            self.phase = "giggle"
            print(f"\n  {bold('◆ Giggle indexing')}")
            return True
        m = self._RE_GIGGLE_DIR.search(line)
        if m:
            d, cur, tot = m.group(1).split("/")[-1], int(m.group(2)), int(m.group(3))
            print(f"    [{cyan(str(cur).rjust(len(str(tot))))}/{tot}]  {dim(d)}")
            return True
        m2 = self._RE_GIGGLE_SKIP.search(line)
        if m2:
            self.warnings += 1
            print(f"    {yellow('⚠')}  Skipping (no .bed.gz): {dim(m2.group(1).split('/')[-1])}")
            return True
        if self._RE_TABIX_START.search(line):
            self.phase = "tabix"
            print(f"\n  {bold('◆ Tabix indexing individual tracks …')}")
            return True
        return False


_PHASE_PATTERNS = {
    "metadata": re.compile(r"Found \d+ track records"),
    "space":    re.compile(r"Required space="),
    "download": re.compile(r"Starting dowloading"),
    "giggle":   re.compile(r"Starting Giggle indexing"),
    "tabix":    re.compile(r"Creating tabix index"),
    "done":     re.compile(r"FILER root directory="),
}


def cmd_install(args: argparse.Namespace) -> int:
    """Route to from-tracks or template install mode."""
    if args.from_tracks:
        return _install_from_tracks(args)
    return _install_from_template(args)


# ── From-tracks install (Recipe 3 / Recipe 10 TSV input) ─────────────────────

# Columns required for a validated install. wget_command is preferred but
# processed_file_download_url is the fallback for Recipe 10 output.
_REQUIRED_COLS    = {"identifier", "file_name", "file_size", "processed_file_md5"}
_PREFERRED_COLS   = _REQUIRED_COLS | {"wget_command"}


def _detect_tsv_source(columns: set[str]) -> str:
    """Return 'r3', 'r10', or raise if the TSV lacks required columns."""
    missing = _REQUIRED_COLS - columns
    if missing:
        raise ValueError(
            f"Input TSV is missing required columns: {', '.join(sorted(missing))}.\n"
        )
    if "wget_command" in columns:
        return "r3"   # wget_command present verbatim
    if "processed_file_download_url" in columns:
        return "r10"  # reconstruct target dir from URL
    raise ValueError("TSV has file_name/file_size/md5 but no wget_command or download URL.")


def _derive_rel_dir(url: str, base_url: str = "https://tf.lisanwanglab.org/GADB/") -> str:
    """
    Derive the relative target directory from a processed_file_download_url.
    Strips the base URL prefix and the filename, leaving the relative path
    that mirrors the server's directory structure locally.

    Example:
      https://tf.lisanwanglab.org/GADB/Annotationtracks/ENCODE/data/ATAC-seq/narrowpeak/hg38/ENCFF124SPE.bed.gz
      → Annotationtracks/ENCODE/data/ATAC-seq/narrowpeak/hg38
    """
    if url.startswith(base_url):
        rel = url[len(base_url):]
    else:
        # Fallback: strip scheme + host
        from urllib.parse import urlparse
        rel = urlparse(url).path.lstrip("/")
    return str(Path(rel).parent)


def _md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _download_file(url: str, dest: Path, expected_size: int,
                   expected_md5: str, track_num: int, total: int,
                   verbose: bool) -> bool:
    """
    Download url → dest, verifying size and MD5.
    Returns True on success, False on failure (non-fatal; caller decides).
    """
    dest.parent.mkdir(parents=True, exist_ok=True)

    # Skip if already present with correct size + MD5
    if dest.exists():
        if dest.stat().st_size == expected_size:
            actual_md5 = _md5_file(dest)
            if actual_md5 == expected_md5:
                if verbose:
                    print(f"\n  {dim(f'  [{track_num}/{total}] Already present, skipping: {dest.name}')}")
                return True
            else:
                print(f"\n  {yellow('⚠')}  MD5 mismatch on existing file, re-downloading: {dest.name}")
        else:
            print(f"\n  {yellow('⚠')}  Size mismatch on existing file, re-downloading: {dest.name}")

    # Download
    fname = dest.name[:50]
    pct_label = f"[{track_num}/{total}]"
    try:
        with urllib.request.urlopen(url, timeout=120) as resp, dest.open("wb") as out:
            downloaded = 0
            while True:
                chunk = resp.read(1 << 16)
                if not chunk:
                    break
                out.write(chunk)
                downloaded += len(chunk)
                if expected_size:
                    pct = int(100 * downloaded / expected_size)
                    filled = int(28 * pct / 100)
                    bar = cyan("█" * filled + "░" * (28 - filled))
                    print(f"\r  {cyan('↓')} {pct_label} [{bar}] {pct:3d}%  {dim(fname):<52}",
                          end="", flush=True)
    except Exception as exc:
        print(f"\n  {red('✗')}  Download failed for {dest.name}: {exc}")
        if dest.exists():
            dest.unlink()
        return False

    # Verify size
    actual_size = dest.stat().st_size
    if actual_size != expected_size:
        print(f"\n  {red('✗')}  Size mismatch: expected {expected_size}, got {actual_size} — {dest.name}")
        dest.unlink()
        return False

    # Verify MD5
    actual_md5 = _md5_file(dest)
    if actual_md5 != expected_md5:
        print(f"\n  {red('✗')}  MD5 mismatch: expected {expected_md5}, got {actual_md5} — {dest.name}")
        dest.unlink()
        return False

    return True


def _run_giggle_index(directories: list[Path], giggle_bin: str) -> None:
    """Run giggle index on each directory that contains .bed.gz files."""
    print(f"\n  {bold('◆ Giggle indexing')}")
    dirs_with_beds = [d for d in directories if any(d.glob("*.bed.gz"))]
    if not dirs_with_beds:
        print(f"  {yellow('⚠')}  No directories with .bed.gz files found; skipping Giggle indexing.")
        return

    for i, d in enumerate(dirs_with_beds, 1):
        index_dir = d / "giggle_index"
        if index_dir.exists():
            shutil.rmtree(index_dir)
        print(f"    [{cyan(str(i).rjust(len(str(len(dirs_with_beds)))))}/{len(dirs_with_beds)}]  {dim(d.name)}")
        result = subprocess.run(
            [giggle_bin, "index", "-s", "-i", str(d / "*.bed.gz"), "-o", str(index_dir)],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            print(f"    {red('✗')}  giggle failed on {d.name}: {result.stderr.strip()}")


def _run_tabix_index(files: list[Path], tabix_bin: str) -> None:
    """Run tabix -f on each .bed.gz or .vcf.gz file."""
    print(f"\n  {bold('◆ Tabix indexing individual tracks …')}")
    for f in files:
        ftype = "vcf" if f.name.endswith(".vcf.gz") else "bed"
        result = subprocess.run(
            [tabix_bin, "-f", "-p", ftype, str(f)],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            print(f"  {yellow('⚠')}  tabix failed on {f.name}: {result.stderr.strip()}")
    print(f"  {green('✓')} Tabix indexing complete.")


def _write_track_metadata(tracks: list[dict], dest_paths: dict[str, Path],
                          meta_out: Path) -> None:
    """
    Write a single TSV containing one row per downloaded track with all
    available metadata columns plus the resolved local path.

    Parameters
    ----------
    tracks : list[dict]
        The deduplicated track rows (original TSV dicts).
    dest_paths : dict[str, Path]
        Mapping of identifier → local file path for successfully downloaded
        tracks.  Tracks not present here are still included but with an
        empty local_path (they may have failed to download).
    meta_out : Path
        Output TSV path.
    """
    if not tracks:
        return

    # Collect the union of all keys across rows, preserving first-seen order,
    # and prepend local_path right after identifier.
    seen_cols: dict[str, None] = {}
    for row in tracks:
        for k in row:
            seen_cols.setdefault(k, None)
    all_cols = list(seen_cols)

    # Ensure identifier is first, then local_path, then everything else.
    ordered_cols = ["identifier", "local_path"]
    ordered_cols += [c for c in all_cols if c not in ("identifier", "local_path")]

    meta_out.parent.mkdir(parents=True, exist_ok=True)
    with meta_out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=ordered_cols, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for row in tracks:
            out_row = dict(row)
            ident = row.get("identifier", "")
            out_row["local_path"] = str(dest_paths.get(ident, ""))
            writer.writerow(out_row)

    print(f"  {green('✓')} Track metadata ({len(tracks)} tracks) saved → {cyan(str(meta_out))}")


def _install_from_tracks(args: argparse.Namespace) -> int:
    """Download and index tracks from a Recipe 3 / Recipe 10 TSV."""

    tsv_path = Path(args.from_tracks)
    if not tsv_path.exists():
        print(red(f"  ✗  Tracks file not found: {tsv_path}"))
        return 2

    # Read TSV
    with tsv_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        columns = set(reader.fieldnames or [])

    # Detect source and validate columns
    try:
        source = _detect_tsv_source(columns)
    except ValueError as exc:
        print(red(f"  ✗  {exc}"))
        return 2

    # Deduplicate on identifier (Recipe 10 has one row per interval, not per track)
    seen_ids: set[str] = set()
    tracks = []
    for row in rows:
        tid = row.get("identifier", "")
        if tid and tid not in seen_ids:
            seen_ids.add(tid)
            tracks.append(row)

    total_unique = len(tracks)

    # Rank by num_overlaps and take top N if requested
    if args.top is not None and args.top > 0:
        if "num_overlaps" in columns:
            def _num_overlaps_key(row: dict) -> int:
                val = row.get("num_overlaps", "")
                try:
                    return int(val) if val else 0
                except ValueError:
                    return 0
            tracks.sort(key=_num_overlaps_key, reverse=True)
        else:
            print(yellow(
                f"  ⚠  --top {args.top} set but num_overlaps column missing; "
                f"taking first {args.top} rows in input order."
            ))
        tracks = tracks[: args.top]

    target_dir = Path(args.target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)

    # Resolve metadata output path
    meta_out = Path(args.meta_out) if args.meta_out else target_dir / "track_metadata.tsv"

    print_header("FILER Install — From Tracks")
    print_kv("Source TSV",   str(tsv_path))
    print_kv("Source type",  f"Recipe {'3' if source == 'r3' else '10'} output")
    if args.top is not None and args.top > 0:
        print_kv("Tracks",   f"{len(tracks)} (top {args.top} by num_overlaps of {total_unique} unique)")
    else:
        print_kv("Tracks",   str(len(tracks)))
    print_kv("Target dir",   str(target_dir))
    print_kv("Metadata out", str(meta_out))
    print_kv("Giggle",       args.giggle)
    print_kv("Tabix",        args.tabix)
    print()

    # Check available disk space vs required
    required_bytes = sum(int(r.get("file_size", 0) or 0) for r in tracks)
    stat = shutil.disk_usage(target_dir)
    print(f"  {dim('Disk space required:')}  {required_bytes:,} bytes")
    print(f"  {dim('Disk space available:')} {stat.free:,} bytes")
    if stat.free < required_bytes:
        print(red(f"\n  ✗  Not enough disk space. Need {required_bytes:,}, have {stat.free:,}."))
        return 1
    print()

    if not estimate_and_confirm(tracks, str(target_dir)):
        print(yellow("  Installation cancelled."))
        return 0

    # Download phase
    print(f"  {bold('◆ Downloading tracks')}")
    downloaded_files: list[Path] = []
    download_dirs:    set[Path]  = set()
    dest_paths:       dict[str, Path] = {}   # identifier → local file path
    warnings = errors = 0

    for i, row in enumerate(tracks, 1):
        fname         = row["file_name"]
        expected_size = int(row.get("file_size", 0) or 0)
        expected_md5  = row.get("processed_file_md5", "").strip()

        # Determine URL and relative target directory
        if source == "r3":
            wget_cmd = row.get("wget_command", "")
            # wget_command format: "wget <url> -P <rel_dir>"
            parts   = wget_cmd.split()
            url     = parts[1] if len(parts) >= 2 else row.get("processed_file_download_url", "")
            rel_dir = parts[3] if len(parts) >= 4 else _derive_rel_dir(
                url if url else row.get("processed_file_download_url", "")
            )
        else:
            url     = row.get("processed_file_download_url", "")
            rel_dir = _derive_rel_dir(url)

        dest = target_dir / rel_dir / fname
        download_dirs.add(dest.parent)

        ok = _download_file(url, dest, expected_size, expected_md5, i, len(tracks), args.verbose)
        if ok:
            downloaded_files.append(dest)
            dest_paths[row.get("identifier", "")] = dest
        else:
            errors += 1
            if not args.keep_going:
                print(red(f"\n  ✗  Stopping after first error. Use --keep-going to continue."))
                # Write metadata for whatever we managed to download so far
                _write_track_metadata(tracks, dest_paths, meta_out)
                return 1

    print(f"\r  {green('✓')} Download complete — "
          f"{len(downloaded_files)}/{len(tracks)} track(s)  "
          f"errors={red(str(errors)) if errors else green('0')}\n")

    # Write track-metadata TSV (always, even if some downloads failed)
    _write_track_metadata(tracks, dest_paths, meta_out)

    if args.skip_index:
        print(dim("  Skipping indexing (--skip-index set)."))
        return 0 if errors == 0 else 1

    # Giggle indexing
    if not shutil.which(args.giggle):
        print(yellow(f"  ⚠  giggle not found at '{args.giggle}'; skipping Giggle indexing."))
    else:
        _run_giggle_index(sorted(download_dirs), args.giggle)

    # Tabix indexing
    if not shutil.which(args.tabix):
        print(yellow(f"  ⚠  tabix not found at '{args.tabix}'; skipping tabix indexing."))
    else:
        indexable = [f for f in downloaded_files if f.name.endswith((".bed.gz", ".vcf.gz"))]
        _run_tabix_index(indexable, args.tabix)

    return 0 if errors == 0 else 1


# ── Template-mode install (wraps install_filer.sh) ────────────────────────────

def _install_from_template(args: argparse.Namespace) -> int:
    """Wrap install_filer.sh for bulk template-based installs."""

    script_path = Path(args.script).resolve()
    if not script_path.exists():
        print(red(f"  ✗  Shell script not found: {script_path}"))
        return 2

    cmd = [
        "bash", str(script_path),
        args.target_dir_pos,
        args.metadata_url,
        args.config,
        "1" if args.force_overwrite else "0",
        "1" if args.force_restart  else "0",
        "1" if args.skip_download  else "0",
    ]

    print_header("FILER Install — Template Mode")
    print_kv("Target dir",  args.target_dir_pos)
    print_kv("Metadata",    args.metadata_url)
    print_kv("Config",      args.config)
    flags = [f for f, v in [("force-overwrite", args.force_overwrite),
                              ("force-restart",   args.force_restart),
                              ("skip-download",   args.skip_download)] if v]
    if flags:
        print_kv("Flags", ", ".join(flags))
    print()

    dl_tracker  = DownloadTracker()
    idx_tracker = IndexTracker()
    current_phase = None

    log_fh = open(args.log, "w") if args.log else None

    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, bufsize=1)
        for raw_line in proc.stdout:
            line = raw_line.rstrip()
            if log_fh:
                log_fh.write(raw_line)
                log_fh.flush()

            for phase_key, pattern in _PHASE_PATTERNS.items():
                if pattern.search(line) and phase_key != current_phase:
                    current_phase = phase_key
                    if phase_key not in ("download", "giggle", "tabix", "done"):
                        label = {"metadata": "Preparing metadata", "space": "Checking disk space"}.get(phase_key, phase_key)
                        print(f"\n  {bold('◆')} {label}")

            if current_phase == "download" and dl_tracker.feed(line):
                continue
            if current_phase in ("giggle", "tabix", "done") and idx_tracker.feed(line):
                continue

            if any(kw in line for kw in ("ERROR", "WARNING", "Log file=", "FILER root", "FILER metadata")):
                prefix = red("✗  ") if "ERROR" in line else yellow("⚠  ") if "WARNING" in line else dim("   ")
                print(f"{prefix}{line.strip()}")
            elif args.verbose:
                print(f"  {dim(line)}")

        proc.wait()
        return proc.returncode
    finally:
        if log_fh:
            log_fh.close()


# ──────────────────────────────────────────────────────────────────────────────
# Argument parser
# ──────────────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="filer_install",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="subcommand", metavar="<subcommand>")
    sub.required = True

    # ── search ────────────────────────────────────────────────────────────────
    sp = sub.add_parser(
        "search",
        help="Query the FILER API for tracks matching metadata filters.",
        description="""
Query the FILER metadata API and save matching tracks to a TSV (and optionally JSON).

Examples:
  # ATAC-seq tracks from ENCODE on hg38
  python filer_install.py search --genome-build hg38 --assay ATAC-seq --data-source ENCODE --out tracks.tsv

  # Hi-C tracks in brain tissue
  python filer_install.py search --genome-build hg38 --assay Hi-C --tissue-category Brain --out tracks.tsv

  # Look up a specific track by ID
  python filer_install.py search --genome-build hg38 --track-id NGBLPL2W2SM2WC --out tracks.tsv

  # Power-user raw jq filter
  python filer_install.py search --genome-build hg38 \\
      --filter-string '.data_source == "ENCODE" and .assay == "ATAC-seq"' --out tracks.tsv
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sp.set_defaults(func=cmd_search)

    required = sp.add_argument_group("Required")
    required.add_argument("--genome-build", required=True, metavar="BUILD",
                          choices=["hg19", "hg38"],
                          help="Reference genome build. Choices: hg19, hg38.")

    filters = sp.add_argument_group("Metadata filters (all optional; combine freely)")
    filters.add_argument("--assay",            metavar="ASSAY",
                         help="Assay type, e.g. 'ATAC-seq', 'ChIP-seq', 'DNase-seq', 'Hi-C', 'WGB-Seq'.")
    filters.add_argument("--cell-type",        metavar="CELL_TYPE",
                         help="Cell type, e.g. 'CD14+ monocyte', 'K562', 'T cell'.")
    filters.add_argument("--tissue-category",  metavar="TISSUE",
                         help="Tissue category, e.g. 'Blood', 'Brain'.")
    filters.add_argument("--data-source",      metavar="SOURCE",
                         help="Data source, e.g. 'ENCODE', 'Roadmap', 'DASHR2', 'Blueprint'.")
    filters.add_argument("--track-id",         metavar="ID",
                         help="Exact FILER track identifier, e.g. 'NGBLPL2W2SM2WC'.")
    filters.add_argument("--filter-string",    metavar="JQ_FILTER",
                         help="Raw jq-style filter string for advanced queries. "
                              "Overrides named filters if both are supplied.")

    out_grp = sp.add_argument_group("Output")
    out_grp.add_argument("--out", metavar="FILE", default="tracks.tsv",
                         help="Path for the output TSV file (default: tracks.tsv). "
                              "Parent directories are created automatically.")
    out_grp.add_argument("--json", action="store_true",
                         help="Also save a JSON file alongside the TSV.")
    out_grp.add_argument("--columns", metavar="COL1,COL2,...",
                         help="Comma-separated list of API columns to include in the TSV. "
                              "Defaults to a curated set of the most useful columns.")
    out_grp.add_argument("--preview-rows", type=int, default=10, metavar="N",
                         help="Number of rows to preview in the terminal (default: 10).")
    out_grp.add_argument("--quiet", "-q", action="store_true",
                         help="Suppress the preview table and summary stats.")

    # ── install ───────────────────────────────────────────────────────────────
    ip = sub.add_parser(
        "install",
        help="Download and index annotation tracks (template mode or from-tracks mode).",
        description="""
Download and index FILER annotation tracks. Two modes are available:

  FROM-TRACKS MODE (recommended for biologically filtered installs)
  ----------------------------------------------------------------
  Reads a Recipe 3 or Recipe 10 TSV produced by filer_filter_then_overlaps.py and downloads
  only those tracks in pure Python — no shell script required.

  python filer_install.py install --from-tracks results.tsv \\
      --target-dir FILER_data --giggle /path/to/giggle --tabix /path/to/tabix

  # Skip indexing (download only)
  python filer_install.py install --from-tracks results.tsv --target-dir FILER_data \\
      --giggle giggle --tabix tabix --skip-index

  # Continue past individual download errors instead of stopping
  python filer_install.py install --from-tracks results.tsv --target-dir FILER_data \\
      --giggle giggle --tabix tabix --keep-going

  TEMPLATE MODE (bulk install everything in a metadata template)
  --------------------------------------------------------------
  Wraps install_filer.sh with a metadata template URL or local file.

  python filer_install.py install FILER_test \\
      https://tf.lisanwanglab.org/FILER/test_metadata.hg19.template filer.ini

  # Resume after interruption
  python filer_install.py install FILER_test metadata.hg19.template filer.ini \\
      --force-restart --skip-download
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ip.set_defaults(func=cmd_install)

    # ── From-tracks mode args ─────────────────────────────────────────────────
    ft = ip.add_argument_group("From-tracks mode (Recipe 3 / Recipe 10 TSV input)")
    ft.add_argument("--from-tracks", metavar="TSV_FILE",
                    help="Path to a Recipe 3 or Recipe 10 output TSV. "
                         "When set, all template-mode positional arguments are ignored.")
    ft.add_argument("--target-dir",  metavar="DIR", default="FILER_data",
                    help="Directory to install tracks into (default: FILER_data).")
    ft.add_argument("--meta-out",    metavar="FILE", default=None,
                    help="Path for the track-metadata TSV (default: <target-dir>/track_metadata.tsv). "
                         "Contains one row per track with all input metadata columns plus local_path.")
    ft.add_argument("--giggle",      metavar="PATH", default="giggle",
                    help="Path or name of the giggle binary (default: giggle).")
    ft.add_argument("--tabix",       metavar="PATH", default="tabix",
                    help="Path or name of the tabix binary (default: tabix).")
    ft.add_argument("--skip-index",  action="store_true",
                    help="Download files but skip Giggle and tabix indexing.")
    ft.add_argument("--keep-going",  action="store_true",
                    help="Continue past individual download failures rather than stopping.")
    ft.add_argument("--top",         metavar="N", type=int, default=None,
                    help="Install only the top N tracks ranked by num_overlaps "
                         "(descending). Applied after deduplication. If num_overlaps "
                         "is absent from the input, takes the first N rows instead.")

    # ── Template mode args (positional, all optional when --from-tracks is set) ─
    tm = ip.add_argument_group("Template mode (wraps install_filer.sh)")
    ip.add_argument("target_dir_pos",   metavar="TARGET_DIR",          nargs="?", default=None,
                    help="Directory where FILER data will be installed.")
    ip.add_argument("metadata_url",     metavar="METADATA_URL_OR_FILE", nargs="?", default=None,
                    help="URL or local path to a FILER metadata template file.")
    ip.add_argument("config",           metavar="CONFIG_FILE",          nargs="?", default=None,
                    help="Path to filer.ini (defines GIGGLE, TABIX paths, etc.).")

    tm.add_argument("--force-overwrite", action="store_true",
                    help="Delete and completely overwrite an existing TARGET_DIR.")
    tm.add_argument("--force-restart",   action="store_true",
                    help="Continue within an existing TARGET_DIR.")
    tm.add_argument("--skip-download",   action="store_true",
                    help="Skip downloading and go straight to indexing.")
    tm.add_argument("--script", metavar="PATH", default="install_filer.sh",
                    help="Path to install_filer.sh (default: ./install_filer.sh).")

    shared = ip.add_argument_group("Shared options")
    shared.add_argument("--verbose", "-v", action="store_true",
                        help="Print detailed output.")
    shared.add_argument("--log", metavar="FILE",
                        help="Write full output to FILE in addition to the console "
                             "(template mode only).")

    return parser


# ──────────────────────────────────────────────────────────────────────────────
# Entrypoint
# ──────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = build_parser()
    args   = parser.parse_args()

    # Guards
    if args.subcommand == "install":
        if not args.from_tracks:
            # Template mode: positional args are required
            if not args.target_dir_pos or not args.metadata_url or not args.config:
                parser.error(
                    "Template mode requires TARGET_DIR, METADATA_URL_OR_FILE, and CONFIG_FILE.\n"
                    "  Use --from-tracks TSV_FILE for from-tracks mode."
                )
        if args.subcommand == "install" and not args.from_tracks \
                and args.force_overwrite and args.force_restart:
            parser.error("--force-overwrite and --force-restart are mutually exclusive.")

    t0 = time.time()
    rc = args.func(args)
    elapsed = time.time() - t0

    mins, secs = divmod(int(elapsed), 60)
    hrs,  mins = divmod(mins, 60)
    elapsed_str = (f"{hrs}h {mins}m {secs}s" if hrs else
                   f"{mins}m {secs}s"         if mins else
                   f"{secs}s")

    print()
    if rc == 0:
        print(green(bold(f"  ✓  Done in {elapsed_str}.")))
    else:
        print(red(bold(f"  ✗  Failed (exit code {rc}) after {elapsed_str}.")))
        sys.exit(rc)


if __name__ == "__main__":
    main()
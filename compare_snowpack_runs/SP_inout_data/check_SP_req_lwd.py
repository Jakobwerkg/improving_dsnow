#!/usr/bin/env python3
"""
End-to-end SNOWPACK data preparation for LWD Tirol .smet(.gz) files.
────────────────────────────────────────────────────────────────────
STEP 1  Strip spurious .gz suffix
         Files are plain ASCII, not compressed; only the extension is wrong.

STEP 2  Check SNOWPACK forcing requirements per station / station-pair
         Files whose first 4 characters match are treated as a pair.
         Required variables:
           TA, RH, VW                         (always required)
           ISWR or RSWR or NET_SW             (shortwave radiation)
           ILWR or TSS                        (longwave / surface temp)
           PSUM or HS                         (precip or snow height)
         A pair passes if the two files together cover all requirements.

STEP 3  Delete files that cannot satisfy requirements even as a pair.

STEP 4  Merge station pairs into a single .smet file.
         - "Windstation" in name  → VW, VW_MAX, DW from that file;
           everything else from the partner (Schneestation / lower alt).
         - No explicit label      → lower-altitude station is primary;
           wind fields pulled from the higher station.
         - Overlapping fields (RH, TA, TD): primary wins; gaps filled
           from secondary.
         Output file carries the primary station's metadata.

STEP 5  Write station_summary.csv
         One row per output station.  For every variable a cell shows
         where the value comes from:
           snow    = primary station (Schneestation / lower alt)
           wind    = secondary station (Windstation / higher alt)
           single  = file was never part of a pair
           (blank) = variable not present

Usage:
    python check_SP_req_lwd.py [path/to/smet/directory]
    Default directory: LWD_all/  next to this script.
"""

import gzip
import os
import re
import shutil
import sys
import pandas as pd
from collections import defaultdict

# ── Configuration ─────────────────────────────────────────────────────────────

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SMET_DIR   = sys.argv[1] if len(sys.argv) > 1 else os.path.join(SCRIPT_DIR, "LWD_all")
CSV_OUT    = os.path.join(SMET_DIR, "station_summary.csv")

WIND_FIELDS   = {"VW", "VW_MAX", "DW"}
NODATA        = -777

# Variables displayed as columns in the CSV (in this order)
TRACK_VARS = [
    "TA", "RH", "VW", "VW_MAX", "DW",
    "ISWR", "RSWR", "NET_SW", "ILWR",
    "TSS", "HS", "PSUM", "TSG", "TD", "P",
]
SNOW_TEMP_RE = re.compile(r"^TS\.\d+$")


# ── Helpers ───────────────────────────────────────────────────────────────────

def read_header_only(filepath):
    """Return (header_dict, fields_set) without reading DATA."""
    header, fields = {}, set()
    with open(filepath, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if line == "[DATA]":
                break
            if line.startswith("#") or line.startswith("SMET") or line.startswith("["):
                continue
            if "=" in line:
                k, _, v = line.partition("=")
                k, v = k.strip(), v.strip()
                if k == "fields":
                    fields = set(v.split()) - {"timestamp"}
                else:
                    header[k] = v
    return header, fields


def read_smet_full(filepath):
    """Return (header_dict, ordered_field_list, DataFrame indexed by timestamp)."""
    header, comment_lines, fields, data_rows = {}, [], [], []
    in_data = False

    with open(filepath, encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if line == "[DATA]":
                in_data = True
                continue
            if in_data:
                if line:
                    data_rows.append(line.split())
            else:
                if line.startswith("#") or line.startswith("SMET") or line.startswith("["):
                    continue
                if "=" in line:
                    k, _, v = line.partition("=")
                    k, v = k.strip(), v.strip()
                    if k == "fields":
                        fields = v.split()
                    else:
                        header[k] = v

    if not fields or not data_rows:
        return header, [], pd.DataFrame()

    # Deduplicate fields (rare but possible, e.g. SSON1 has DW twice)
    seen, unique_fields = set(), []
    for f in fields:
        if f not in seen:
            seen.add(f)
            unique_fields.append(f)

    df = pd.DataFrame(
        [row[:len(unique_fields)] for row in data_rows],
        columns=unique_fields,
    )
    df["timestamp"] = pd.to_datetime(df["timestamp"], utc=True, errors="coerce")
    df = df.dropna(subset=["timestamp"]).set_index("timestamp")

    nodata_val = float(header.get("nodata", NODATA))
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce").replace(nodata_val, float("nan"))

    non_ts = [c for c in unique_fields if c != "timestamp"]
    return header, non_ts, df


def write_smet(filepath, header, fields, df):
    """Write a SMET 1.2 ASCII file."""
    nodata_val = int(header.get("nodata", NODATA))
    key_order  = ["station_id", "station_name", "latitude", "longitude", "altitude", "source"]
    ordered    = {k: header[k] for k in key_order if k in header}
    for k, v in header.items():
        if k not in ordered and k != "nodata":
            ordered[k] = v

    with open(filepath, "w", encoding="utf-8") as f:
        f.write("SMET 1.2 ASCII\n[HEADER]\n")
        for k, v in ordered.items():
            f.write(f"{k} = {v}\n")
        f.write(f"nodata = {nodata_val}\n")
        f.write(f"fields = timestamp {' '.join(fields)}\n")
        f.write("[DATA]\n")
        for ts, row in df.iterrows():
            vals = []
            for col in fields:
                v = row.get(col, float("nan"))
                vals.append(str(nodata_val) if pd.isna(v) else f"{v:.6g}")
            f.write(f"{ts.strftime('%Y-%m-%dT%H:%M:%SZ')} {' '.join(vals)}\n")


def check_requirements(fields):
    """Return (ok, missing_list) for a set of field names."""
    missing = []
    for v in ["TA", "RH", "VW"]:
        if v not in fields:
            missing.append(v)
    if not (fields & {"ISWR", "RSWR", "NET_SW"}):
        missing.append("ISWR/RSWR/NET_SW")
    if not (fields & {"ILWR", "TSS"}):
        missing.append("ILWR/TSS")
    if not (fields & {"PSUM", "HS"}):
        missing.append("PSUM/HS")
    return len(missing) == 0, missing


def is_wind_station(name):
    return "windstation" in name.lower()


def section(title):
    print(f"\n{'═'*70}")
    print(f"  {title}")
    print(f"{'═'*70}")


# ── STEP 1 – Strip .gz suffix (decompress or rename) ─────────────────────────

def step1_strip_gz():
    section("STEP 1 · Remove .gz from filenames (decompress if real gzip, else rename)")
    processed = 0
    for fname in sorted(os.listdir(SMET_DIR)):
        if not fname.endswith(".gz"):
            continue
        src = os.path.join(SMET_DIR, fname)
        dst = os.path.join(SMET_DIR, fname[:-3])
        # Try real gzip decompression first
        try:
            with gzip.open(src, "rb") as f_in, open(dst, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(src)
            print(f"  decompressed  {fname}  →  {fname[:-3]}")
        except (gzip.BadGzipFile, OSError):
            # Not actually gzip – just strip the extension
            os.rename(src, dst)
            print(f"  renamed       {fname}  →  {fname[:-3]}  (not gzip, extension only)")
        processed += 1
    if processed == 0:
        print("  No .gz files found – already done or not needed.")
    else:
        print(f"  {processed} files processed.")


# ── STEP 2 – Check requirements ───────────────────────────────────────────────

def step2_check_requirements():
    section("STEP 2 · Check SNOWPACK forcing requirements")

    smet_files = sorted(f for f in os.listdir(SMET_DIR) if f.endswith(".smet"))
    file_fields = {}
    for fname in smet_files:
        _, fields = read_header_only(os.path.join(SMET_DIR, fname))
        file_fields[fname] = fields

    groups = defaultdict(list)
    for fname in smet_files:
        groups[fname[:4]].append(fname)

    passing_single, passing_pair, failing = [], [], []

    for prefix, members in sorted(groups.items()):
        if len(members) == 1:
            fname  = members[0]
            ok, missing = check_requirements(file_fields[fname])
            entry = {"files": [fname], "missing": missing}
            (passing_single if ok else failing).append(entry)

        else:
            # First try each file individually
            indiv_ok = []
            for fname in members:
                ok, missing = check_requirements(file_fields[fname])
                if ok:
                    indiv_ok.append({"files": [fname], "missing": []})

            passing_single.extend(indiv_ok)
            indiv_ok_names = {e["files"][0] for e in indiv_ok}
            remaining = [f for f in members if f not in indiv_ok_names]

            if remaining:
                combined = set()
                for f in members:
                    combined |= file_fields[f]
                ok, missing = check_requirements(combined)
                if ok:
                    passing_pair.append({"files": members, "missing": []})
                else:
                    # Neither alone nor together → all fail
                    failing.append({"files": members, "missing": missing})

    # Print
    print(f"\n  Passing (single file): {len(passing_single)}")
    for e in passing_single:
        print(f"    ✓  {e['files'][0]}")

    print(f"\n  Passing (pair, combined): {len(passing_pair)}")
    for e in passing_pair:
        print(f"    ✓  {e['files'][0]}  +  {e['files'][1]}")

    print(f"\n  Failing: {len(failing)}")
    for e in failing:
        names = " + ".join(e["files"])
        print(f"    ✗  {names}")
        print(f"       missing: {', '.join(e['missing'])}")

    # Build the set of files to keep
    keep = set()
    for e in passing_single + passing_pair:
        keep.update(e["files"])

    return keep, file_fields


# ── STEP 3 – Delete failing files ─────────────────────────────────────────────

def step3_delete_failing(keep):
    section("STEP 3 · Delete stations that cannot satisfy requirements")
    smet_files = sorted(f for f in os.listdir(SMET_DIR) if f.endswith(".smet"))
    deleted = 0
    for fname in smet_files:
        if fname not in keep:
            os.remove(os.path.join(SMET_DIR, fname))
            print(f"  deleted  {fname}")
            deleted += 1
    if deleted == 0:
        print("  Nothing to delete.")
    else:
        print(f"  {deleted} files deleted.")


# ── STEP 4 – Merge station pairs ──────────────────────────────────────────────

def step4_merge_pairs():
    """
    Returns merge_log: dict  output_fname → {
        "primary_file":  original filename of primary station,
        "primary_name":  station_name of primary,
        "wind_file":     original filename of wind/secondary station (or None),
        "wind_name":     station_name of wind/secondary (or None),
        "was_merged":    bool,
        "fields":        list of output field names,
    }
    """
    section("STEP 4 · Merge station pairs")

    smet_files = sorted(f for f in os.listdir(SMET_DIR) if f.endswith(".smet"))
    groups = defaultdict(list)
    for fname in smet_files:
        groups[fname[:4]].append(fname)

    merge_log = {}

    for prefix, members in sorted(groups.items()):
        if len(members) == 1:
            fname = members[0]
            h, _, _ = read_smet_full(os.path.join(SMET_DIR, fname))
            merge_log[fname] = {
                "primary_file": fname,
                "primary_name": h.get("station_name", ""),
                "primary_alt":  float(h.get("altitude", 0)),
                "wind_file":    None,
                "wind_name":    None,
                "was_merged":   False,
            }
            print(f"  [single]  {fname}  ({h.get('station_name','?')})")
            continue

        if len(members) > 2:
            print(f"  [WARN]  prefix={prefix!r}: {len(members)} files – skipped (manual review needed)")
            for m in members:
                print(f"            {m}")
            continue

        f1, f2 = members
        h1, flds1, df1 = read_smet_full(os.path.join(SMET_DIR, f1))
        h2, flds2, df2 = read_smet_full(os.path.join(SMET_DIR, f2))

        name1 = h1.get("station_name", "")
        name2 = h2.get("station_name", "")
        alt1  = float(h1.get("altitude", 9999))
        alt2  = float(h2.get("altitude", 9999))

        # Decide primary vs secondary (wind) station
        if is_wind_station(name1) and not is_wind_station(name2):
            ph, pflds, pdf, pfile = h2, flds2, df2, f2
            wh, wflds, wdf, wfile = h1, flds1, df1, f1
        elif is_wind_station(name2) and not is_wind_station(name1):
            ph, pflds, pdf, pfile = h1, flds1, df1, f1
            wh, wflds, wdf, wfile = h2, flds2, df2, f2
        elif alt1 <= alt2:
            ph, pflds, pdf, pfile = h1, flds1, df1, f1
            wh, wflds, wdf, wfile = h2, flds2, df2, f2
        else:
            ph, pflds, pdf, pfile = h2, flds2, df2, f2
            wh, wflds, wdf, wfile = h1, flds1, df1, f1

        # Wind fields available in the secondary
        avail_wind = [c for c in wdf.columns if c in WIND_FIELDS]
        extra_wind = [c for c in avail_wind if c not in pdf.columns]

        out_fields = pflds + extra_wind
        merged     = pdf.copy()

        # Add wind fields not in primary
        for col in extra_wind:
            merged = merged.join(wdf[[col]], how="outer")

        # Fill gaps in shared wind fields from secondary
        for col in avail_wind:
            if col in merged.columns and col in wdf.columns:
                merged[col] = merged[col].combine_first(wdf[col])

        merged = merged.sort_index()

        # Output filename = primary station_id
        primary_id = ph.get("station_id",
                             os.path.splitext(os.path.basename(pfile))[0])
        out_fname  = f"{primary_id}.smet"
        out_path   = os.path.join(SMET_DIR, out_fname)

        # Write via tmp to avoid clobbering a source file in place
        tmp = out_path + ".tmp"
        write_smet(tmp, ph, out_fields, merged)
        os.remove(os.path.join(SMET_DIR, f1))
        os.remove(os.path.join(SMET_DIR, f2))
        os.rename(tmp, out_path)

        merge_log[out_fname] = {
            "primary_file": pfile,
            "primary_name": ph.get("station_name", ""),
            "primary_alt":  float(ph.get("altitude", 0)),
            "wind_file":    wfile,
            "wind_name":    wh.get("station_name", ""),
            "was_merged":   True,
            "fields":       out_fields,
        }

        print(f"  [merged]  {f1}  +  {f2}")
        print(f"            primary  = {ph.get('station_name','?')}  "
              f"(alt {ph.get('altitude','?')} m)  [{pfile}]")
        print(f"            wind src = {wh.get('station_name','?')}  "
              f"(alt {wh.get('altitude','?')} m)  [{wfile}]")
        print(f"            output   = {out_fname}  "
              f"rows={len(merged)}  fields={' '.join(out_fields)}")

    return merge_log


# ── STEP 5 – Write summary CSV ────────────────────────────────────────────────

def step5_write_csv(merge_log):
    section("STEP 5 · Write station_summary.csv")

    rows = []
    for out_fname, info in sorted(merge_log.items()):
        fpath = os.path.join(SMET_DIR, out_fname)
        header, field_set = read_header_only(fpath)

        was_merged   = info["was_merged"]
        primary_file = info["primary_file"]
        wind_file    = info.get("wind_file") or ""
        primary_name = info["primary_name"]
        wind_name    = info.get("wind_name") or ""

        def source_label(var):
            if not was_merged:
                return "single" if var in field_set else ""
            if var in field_set:
                return "wind" if var in WIND_FIELDS else "snow"
            return ""

        snow_temp_fields = sorted(f for f in field_set if SNOW_TEMP_RE.match(f))

        row = {
            "output_file":        out_fname,
            "station_id":         header.get("station_id", ""),
            "station_name":       header.get("station_name", ""),
            "altitude_m":         header.get("altitude", ""),
            "latitude":           header.get("latitude", ""),
            "longitude":          header.get("longitude", ""),
            "was_merged":         "yes" if was_merged else "no",
            "primary_source_file": primary_file,
            "primary_station_name": primary_name,
            "wind_source_file":   wind_file,
            "wind_station_name":  wind_name,
        }

        for var in TRACK_VARS:
            row[var] = source_label(var)

        row["snow_temp_layers"] = ", ".join(snow_temp_fields) if snow_temp_fields else ""

        # Requirement satisfaction flags
        row["SW_ok"]  = "yes" if field_set & {"ISWR", "RSWR", "NET_SW"} else "MISSING"
        row["LW_ok"]  = "yes" if field_set & {"ILWR", "TSS"}            else "MISSING"
        row["SH_ok"]  = "yes" if field_set & {"PSUM", "HS"}             else "MISSING"

        rows.append(row)

    col_order = (
        ["output_file", "station_id", "station_name", "altitude_m",
         "latitude", "longitude", "was_merged",
         "primary_source_file", "primary_station_name",
         "wind_source_file",    "wind_station_name"]
        + TRACK_VARS
        + ["snow_temp_layers", "SW_ok", "LW_ok", "SH_ok"]
    )

    df = pd.DataFrame(rows, columns=col_order).sort_values("station_id").reset_index(drop=True)
    df.to_csv(CSV_OUT, index=False)
    print(f"  Saved {CSV_OUT}")
    print(f"  {len(df)} stations  ·  {len(col_order)} columns\n")
    print(df.to_string(index=False))


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    print(f"\n{'━'*70}")
    print(f"  LWD SNOWPACK data preparation")
    print(f"  Directory: {SMET_DIR}")
    print(f"{'━'*70}")

    step1_strip_gz()
    keep, _       = step2_check_requirements()
    step3_delete_failing(keep)
    merge_log     = step4_merge_pairs()
    step5_write_csv(merge_log)

    section("DONE")
    n = sum(1 for f in os.listdir(SMET_DIR) if f.endswith(".smet"))
    print(f"  {n} ready-to-use .smet files in {SMET_DIR}")
    print(f"  Summary CSV: {CSV_OUT}\n")


if __name__ == "__main__":
    main()

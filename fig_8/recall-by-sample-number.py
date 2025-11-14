#!/usr/bin/env python3
import argparse, os
from collections import defaultdict, OrderedDict
import pysam
import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(
        description=("Sensitivity by SN for DEL/INS/DUP for a target sample. "
                     "Writes FN lists, a standard table, an adjusted table (removing FNs shared by all callers), "
                     "and two sets of line plots (one per table; one figure per SVTYPE)."))
    p.add_argument("benchmark_vcf")
    p.add_argument("target_sample")
    p.add_argument("compare_files", nargs="+", help="caller-compare txt files (IDs correctly called)")
    p.add_argument("--outdir", required=True, help="Output directory. Will contain fn/, tables, and plots/.")
    p.add_argument("--id-delim", default=";", help="Delimiter when VCF ID has multiple tokens (default ';').")
    p.add_argument("--dup-ids", default=None,
                   help="Optional file with IDs (one per line) to force SVTYPE=DUP, overriding VCF.")
    return p.parse_args()


colors = ['#ff7f0e', '#1f77b4', '#2ca02c']

def capitalize(s):
    return s[0].upper() + s[1:]

# ---------------- Core helpers ----------------

def sample_has_alt(smp):
    gt = smp.get("GT")
    if not gt:
        return False
    return any(a is not None and a > 0 for a in gt)

def normalize_svtype(rec):
    svtype = rec.info.get("SVTYPE")
    if isinstance(svtype, (list, tuple)):
        svtype = svtype[0] if svtype else None
    if isinstance(svtype, bytes):
        svtype = svtype.decode()
    if isinstance(svtype, str):
        st = svtype.upper()
        if st in ("DEL", "INS", "DUP"):
            return st
    alts = rec.alts or ()
    if alts:
        alt0 = str(alts[0]).strip("<>").upper()
        if alt0 in ("DEL", "INS", "DUP"):
            return alt0
    return None

def read_id_set(path):
    ids = set()
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            ids.add(s.split()[0])
    return ids

def index_benchmark(vcf_path, target_sample, id_delim, dup_override_ids):
    vcf = pysam.VariantFile(vcf_path)
    if target_sample not in vcf.header.samples:
        raise SystemExit(f"Target sample '{target_sample}' not in VCF samples.")
    samples = list(vcf.header.samples)

    svtypes = ("DEL", "INS", "DUP")
    data = {svt: {"sn_to_ids_target": defaultdict(set)} for svt in svtypes}

    for rec in vcf.fetch():
        vid = rec.id
        if not vid or vid == ".":
            continue
        if not sample_has_alt(rec.samples[target_sample]):
            continue

        sn = sum(1 for s in samples if sample_has_alt(rec.samples[s]))
        if sn == 0:
            continue

        rec_ids = [i for i in str(vid).split(id_delim) if i]
        if dup_override_ids and any(i in dup_override_ids for i in rec_ids):
            svt = "DUP"
        else:
            svt = normalize_svtype(rec)
        if svt not in data:
            continue

        for iid in rec_ids:
            data[svt]["sn_to_ids_target"][sn].add(iid)

    return data

def caller_label(path):
    base = os.path.basename(path)
    return base[:-len("-compare.txt")] if base.endswith("-compare.txt") else os.path.splitext(base)[0]


# ---------------- Table building ----------------

def compute_rows(data, callers, dump_fn_dir=None, adjusted=False):
    """
    Build table rows.
    If adjusted=True, exclude FNs shared by all callers from denominators (treat as benchmark FPs).
    Optionally write per-caller FN lists to dump_fn_dir (only for the non-adjusted run).
    Returns (header, rows).
    """
    all_sns = sorted({sn for svt in data for sn in data[svt]["sn_to_ids_target"].keys()})

    overall_ids = {svt: set() for svt in data}
    overall_tp = {svt: {lbl: 0 for lbl in callers} for svt in data}
    overall_fn_by_caller = {svt: {lbl: set() for lbl in callers} for svt in data}

    # Header: TOTAL_<SVTYPE>, then for each caller: <caller>_<SVTYPE> [and <caller>_SHARED_FN_%_<SVTYPE> if standard]
    header = ["SN"]
    for svt in ("DEL", "INS", "DUP"):
        header.append(f"TOTAL_{svt}")
        for lbl in callers:
            header.append(f"{lbl}_{svt}")
            if not adjusted:
                header.append(f"{lbl}_SHARED_FN_%_{svt}")

    rows = []
    for sn in all_sns:
        row = [str(sn)]
        for svt in ("DEL", "INS", "DUP"):
            ids_raw = set(data[svt]["sn_to_ids_target"].get(sn, set()))

            # FN sets (raw)
            fn_sets = {lbl: (ids_raw - idset) for lbl, idset in callers.items()}

            # shared by all callers (raw)
            if callers and ids_raw:
                iter_sets = list(fn_sets.values())
                shared_all_raw = set(iter_sets[0]).intersection(*iter_sets[1:]) if iter_sets else set()
            else:
                shared_all_raw = set()

            # Adjust denominator if requested
            ids = set(ids_raw)
            if adjusted and shared_all_raw:
                ids -= shared_all_raw

            denom = len(ids)
            row.append(str(denom))
            overall_ids[svt].update(ids)

            # Append per-caller *in header order*: sensitivity, then (if standard) that caller's shared-FN%
            for lbl, idset in callers.items():
                tp = len(ids & idset)
                sens = tp / denom if denom else 0.0
                row.append(f"{sens:.6f}")
                overall_tp[svt][lbl] += tp

                if not adjusted:
                    caller_fn = ids_raw - idset
                    overall_fn_by_caller[svt][lbl].update(caller_fn)
                    if dump_fn_dir:
                        os.makedirs(dump_fn_dir, exist_ok=True)
                        out_path = os.path.join(dump_fn_dir, f"{lbl}-{svt}-{sn}.fns")
                        with open(out_path, "w") as fh:
                            for iid in sorted(caller_fn):
                                fh.write(iid + "\n")
                    pct = (len(shared_all_raw & caller_fn) / len(caller_fn) * 100.0) if caller_fn else 0.0
                    row.append(f"{pct:.2f}%")

        rows.append(row)

    # OVERALL row (same header order)
    over_row = ["OVERALL"]
    for svt in ("DEL", "INS", "DUP"):
        total = len(overall_ids[svt])
        over_row.append(str(total))

        # overall shared-all (for % in standard table)
        if not adjusted:
            if callers:
                iter_sets = list(overall_fn_by_caller[svt].values())
                shared_all_overall = set(iter_sets[0]).intersection(*iter_sets[1:]) if iter_sets else set()
            else:
                shared_all_overall = set()

        for lbl in callers:
            sens = overall_tp[svt][lbl] / total if total else 0.0
            over_row.append(f"{sens:.6f}")
            if not adjusted:
                caller_fn_overall = overall_fn_by_caller[svt][lbl]
                pct = (len(shared_all_overall & caller_fn_overall) / len(caller_fn_overall) * 100.0) if caller_fn_overall else 0.0
                over_row.append(f"{pct:.2f}%")

    rows.append(over_row)
    return header, rows


def write_tsv(path, header, rows):
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(r) + "\n")


# ---------------- Plotting ----------------

def _to_float(cell):
    """Parse floats that may have a trailing '%'."""
    s = str(cell).strip()
    if s.endswith("%"):
        s = s[:-1]
    return float(s)

def plot_tables(out_png_prefix, header, rows, svtypes=("DEL", "INS", "DUP")):
    """
    Make line plots (SN vs sensitivity) per SVTYPE.
    One PNG per SVTYPE, one line per caller.
    Auto-discovers all '<caller>_<SVTYPE>' columns; ignores 'TOTAL_*' and '*_SHARED_FN_%_*'.
    """
    os.makedirs(os.path.dirname(out_png_prefix), exist_ok=True)

    sn_rows = [r for r in rows if r and r[0] != "OVERALL"]
    name_to_idx = {h: i for i, h in enumerate(header)}

    for svt in svtypes:
        # Sensitivity columns for this SVTYPE
        caller_cols = {}
        for i, col in enumerate(header):
            if not col or col == "SN":         continue
            if col.startswith("TOTAL_"):       continue
            if "_SHARED_FN_%_" in col:         continue
            if col.endswith(f"_{svt}"):
                caller = col[:-(len(svt) + 1)]
                if caller:
                    caller_cols[caller] = i

        # Optional: skip rows with zero denom if TOTAL_<SVTYPE> exists
        total_idx = name_to_idx.get(f"TOTAL_{svt}", None)

        # Debug (appears on stdout)
        print(f"[plot:{svt}] columns -> " + ", ".join(f"{c}@{i}" for c, i in sorted(caller_cols.items())))

        if not caller_cols:
            continue

        series = {caller: {"x": [], "y": []} for caller in sorted(caller_cols)}
        for r in sn_rows:
            try:
                sn_val = int(r[0])
            except (ValueError, TypeError, IndexError):
                continue

            if total_idx is not None:
                try:
                    denom = int(r[total_idx])
                    if denom == 0:
                        continue
                except (ValueError, TypeError, IndexError):
                    pass

            for caller, idx in caller_cols.items():
                if idx >= len(r):
                    continue
                try:
                    y = _to_float(r[idx])
                except (ValueError, TypeError):
                    continue
                series[caller]["x"].append(sn_val)
                series[caller]["y"].append(y)

        if all(len(v["x"]) == 0 for v in series.values()):
            continue

        plt.figure(figsize=(8, 6))
        color_i = 0
        for caller, data in series.items():
            if data["x"]:
                xs, ys = zip(*sorted(zip(data["x"], data["y"])))
                plt.plot(xs, ys, marker="o", label=capitalize(caller), color=colors[color_i])
                color_i += 1
        plt.xlabel("Sample count")
        plt.ylabel("Sensitivity")
        # plt.title(f"Sensitivity vs SC ({svt})")
        plt.grid(True, linestyle=":", linewidth=0.8)
        plt.ylim(0, 1)
        plt.xlim(left=0)
        # plt.legend(loc="best", frameon=True)
        out_path = f"{out_png_prefix}_{svt}.png"
        plt.tight_layout()
        plt.savefig(out_path, dpi=150)
        plt.close()


# ---------------- Main ----------------

def main():
    a = parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    fn_dir = os.path.join(a.outdir, "fn")
    plots_dir = os.path.join(a.outdir, "plots")

    dup_override_ids = read_id_set(a.dup_ids) if a.dup_ids else None
    callers = OrderedDict((caller_label(p), read_id_set(p)) for p in a.compare_files)
    data = index_benchmark(a.benchmark_vcf, a.target_sample, a.id_delim, dup_override_ids)

    # Standard (raw) table + FN dumps
    header_std, rows_std = compute_rows(data, callers, dump_fn_dir=fn_dir, adjusted=False)
    std_path = os.path.join(a.outdir, "sensitivity.tsv")
    write_tsv(std_path, header_std, rows_std)

    # Plots (one figure per SVTYPE) for each table
    plot_tables(os.path.join(plots_dir, "standard"), header_std, rows_std)

    print("Wrote:")
    print(f"  {std_path}")
    print(f"  FN lists in {fn_dir}/")
    print(f"  Plots in {plots_dir}/")


if __name__ == "__main__":
    main()

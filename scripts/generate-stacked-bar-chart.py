#!/usr/bin/env python3
import os
import argparse
import matplotlib.pyplot as plt
from collections import defaultdict, OrderedDict
import numpy as np
from matplotlib.patches import Patch

# --------------------------
# Mapping & data readers
# --------------------------

def read_mapping_file(mapping_file):
    """
    Reads lines of:
        CID  COLOR  LEGEND  METHOD_LABEL
    - Blank line => visual gap between caller groups.
    Returns:
      cid_mapping[cid] = {
          "legend": str,      # caller display name (groups multiple CIDs)
          "color": str,       # hex color to use for that caller
          "method": str,      # arbitrary label (e.g., 'strict', 'loose', etc.)
          "is_gap": bool
      }
      legend_sequence = list with legends in display order and None for gaps
      method_order    = list of method labels in order of first appearance
    """
    cid_mapping = OrderedDict()
    legend_sequence = []
    seen_legends = set()
    method_order = []
    seen_methods = set()

    with open(mapping_file, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line:
                # visual gap marker
                legend_sequence.append(None)
                continue

            parts = line.split(maxsplit=3)
            if len(parts) < 4:
                raise ValueError("Mapping file must have 4 columns: CID COLOR LEGEND METHOD_LABEL")

            cid, color, legend_name, method_label = parts
            method_label = method_label.strip()

            cid_mapping[cid] = {
                "legend": legend_name,
                "color": color,
                "method": method_label,
                "is_gap": False
            }

            # Preserve legend order (first time only)
            if legend_name not in seen_legends:
                legend_sequence.append(legend_name)
                seen_legends.add(legend_name)

            # Preserve global method order (first appearance wins)
            if method_label not in seen_methods:
                method_order.append(method_label)
                seen_methods.add(method_label)

    return cid_mapping, legend_sequence, method_order


def read_sid_mapping(slabel_sid_path):
    with open(slabel_sid_path, 'r') as f:
        return OrderedDict((line.split()[0], ' '.join(line.split()[1:])) for line in f)


# --------------------------
# Legend helpers
# --------------------------

from matplotlib.patches import Patch

def build_caller_handles(cid_mapping, legend_sequence):
    """
    One handle per caller (legend), in the visual order from the mapping.
    Single column by design.
    """
    from collections import OrderedDict
    legend_to_color = OrderedDict()
    for cid, info in cid_mapping.items():
        lg = info["legend"]
        if lg not in legend_to_color:
            legend_to_color[lg] = info["color"]

    handles = []
    seen = set()
    for item in legend_sequence:
        if item is None:
            continue
        lg = item
        if lg in seen:
            continue
        if lg in legend_to_color:
            handles.append(Patch(facecolor=legend_to_color[lg],
                                 edgecolor='black', linewidth=0.5, label=lg))
            seen.add(lg)
    return handles


def build_method_handles(method_order):
    """
    Methods use neutral, opaque swatches; hatch indicates higher layers.
    Single column by design.
    """
    neutral = '#BBBBBB'
    hatches = ['', '//', '\\\\', 'xx', '++', '..', 'oo']  # base has no hatch
    handles = []
    for i, m in enumerate(method_order):
        handles.append(Patch(facecolor=neutral,
                             edgecolor='black',
                             linewidth=0.5,
                             hatch=hatches[i % len(hatches)],
                             label=m))
    return handles


from matplotlib.patches import Patch
from matplotlib.transforms import Bbox

def _save_single_column_legend(out_folder, handles, filename):
    """
    Save a legend-only PDF with:
      - ONE vertical column (ncol=1)
      - minimal whitespace (figure sized to legend bbox)
      - robust handling so color patches don't disappear
    """
    import os
    import matplotlib.pyplot as plt

    os.makedirs(out_folder, exist_ok=True)

    # Fallback: tiny stub if nothing to render
    if not handles:
        fig = plt.figure(figsize=(1, 1))
        fig.savefig(os.path.join(out_folder, filename),
                    bbox_inches="tight", pad_inches=0.01)
        plt.close(fig)
        return

    # Start tiny; we’ll resize to the legend’s true bbox
    fig = plt.figure(figsize=(1, 1), dpi=200)
    # Single, figure-level legend (no axes at all)
    leg = fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.0, 1.0),
        frameon=False,
        ncol=1,                 # ONE COLUMN, always
        fontsize=10,
        handlelength=1.8,
        handletextpad=0.6,
        borderpad=0.0,
        labelspacing=0.25,
        columnspacing=0.8
    )

    # Draw once and measure legend bbox in pixels
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bbox_px = leg.get_window_extent(renderer=renderer)

    # Slight padding factor (2%) to avoid clipping outlines/hatches
    bbox_px = bbox_px.expanded(1.02, 1.02)

    # Convert to inches and resize the whole figure to exactly fit
    width_in = bbox_px.width / fig.dpi
    height_in = bbox_px.height / fig.dpi
    fig.set_size_inches(width_in, height_in)

    # Re-anchor legend to the top-left of the (now resized) figure
    leg.set_bbox_to_anchor((0.0, 1.0), transform=fig.transFigure)

    # Save tightly; no extra padding
    out_path = os.path.join(out_folder, filename)
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.0)
    plt.close(fig)

def save_methods_legend(out_folder, method_handles, filename="fig.legend.methods.pdf"):
    _save_single_column_legend(out_folder, method_handles, filename)

def save_callers_legend(out_folder, caller_handles, filename="fig.legend.callers.pdf"):
    _save_single_column_legend(out_folder, caller_handles, filename)


# --------------------------
# Parsing result files & F1
# --------------------------

def parse_files(slabel_sid_path, clabel_cid_path, out_folder,
                legend_outside=False,        # kept for completeness (ignored when separate_legend=True)
                separate_legend=False):

    sid_mapping = read_sid_mapping(slabel_sid_path)
    cid_mapping, legend_sequence, method_order = read_mapping_file(clabel_cid_path)

    # Build legend handles once
    caller_handles = build_caller_handles(cid_mapping, legend_sequence)
    method_handles = build_method_handles(method_order)

    # data[SID][CID] = {metric: value}
    from collections import defaultdict
    data = defaultdict(lambda: defaultdict(dict))
    metrics_to_collect = ["DEL", "DUP", "INS", "INV"]

    # Collect values from <CID>/<SID>.txt (unchanged)
    for cid, info in cid_mapping.items():
        for sid in sid_mapping:
            fp = os.path.join(cid, f"{sid}.txt")
            if not os.path.exists(fp):
                continue
            with open(fp, 'r') as f:
                for line in f:
                    for metric in metrics_to_collect:
                        if f"{metric} SENSITIVITY" in line:
                            try:
                                data[sid][cid][f"{metric}_SENSITIVITY"] = float(line.split()[4].strip("()"))
                            except Exception:
                                pass
                        elif f"{metric} PRECISION" in line:
                            try:
                                data[sid][cid][f"{metric}_PRECISION"] = float(line.split()[4].strip("()"))
                            except Exception:
                                pass
            # F1
            for metric in metrics_to_collect:
                s_key, p_key = f"{metric}_SENSITIVITY", f"{metric}_PRECISION"
                if s_key in data[sid][cid] and p_key in data[sid][cid]:
                    s, p = data[sid][cid][s_key], data[sid][cid][p_key]
                    data[sid][cid][f"{metric}_F1"] = 2*(s*p)/(s+p) if (s+p) > 0 else 0.0

    # Render plots (no legends on panels if separate_legend=True)
    metrics = [
        "DEL_SENSITIVITY", "DEL_PRECISION", "DEL_F1",
        "DUP_SENSITIVITY", "DUP_PRECISION", "DUP_F1",
        "INS_SENSITIVITY", "INS_PRECISION", "INS_F1",
        "INV_SENSITIVITY", "INV_PRECISION", "INV_F1"
    ]
    for metric in metrics:
        plot_stacked_barchart(
            data=data,
            metric=metric,
            out_folder=out_folder,
            cid_mapping=cid_mapping,
            sid_mapping=sid_mapping,
            legend_sequence=legend_sequence,
            method_order=method_order,
            legend_outside=legend_outside,
            show_legends=not separate_legend
        )

    # Write separate PDFs when requested
    if separate_legend:
        save_methods_legend(out_folder, method_handles, filename="fig.legend.methods.pdf")
        save_callers_legend(out_folder, caller_handles, filename="fig.legend.callers.pdf")


# --------------------------
# Plotting: stacked bars per legend (caller), N methods
# --------------------------

def plot_stacked_barchart(
    data,
    metric,
    out_folder,
    cid_mapping,
    sid_mapping,
    legend_sequence,
    method_order,
    gap_frac=0.4,
    legend_outside=False,
    show_legends=True
):
    """
    For each caller (legend), stack methods in the order they first appear in the mapping file.
    Base = first method; each subsequent method stacks the positive delta over the previous method.
    Total height equals the last method's value. Supports 2+ methods.

    Legends:
      - If show_legends=True:
          - "Methods" (neutral gray swatches) = top-left INSIDE the axes
          - "Callers" (actual caller colors)   = top-right INSIDE the axes (default),
                                                 or outside on the right if legend_outside=True.
      - If show_legends=False: no legends are drawn on the panel.
    """
    # ---- X-axis = samples in the order from slabel_sid_path (stable) ----
    sids_in_order = list(sid_mapping.keys())
    labels = [sid_mapping[sid] for sid in sids_in_order]

    # ---- Group CIDs by legend and method ----
    legend_to_cids = defaultdict(dict)  # legend -> {method_label: cid}
    for cid, info in cid_mapping.items():
        legend_to_cids[info["legend"]][info["method"]] = cid

    # ---- Visual legend order with gaps preserved (None = gap) ----
    visual_legend_seq = []
    for item in legend_sequence:
        if item is None:
            visual_legend_seq.append(None)
        else:
            if item in legend_to_cids:  # keep only legends that actually exist
                visual_legend_seq.append(item)

    # ---- Count groups and gaps ----
    real_groups = [lg for lg in visual_legend_seq if lg is not None]
    N = len(real_groups)
    G = sum(1 for lg in visual_legend_seq if lg is None)

    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(labels))

    # ---- Width and offsets per legend group ----
    width = 0.8 / max(N, 1)
    gap_width = gap_frac * width
    total_span = (N - 1) * width + G * gap_width
    start_offset = - total_span / 2.0

    gaps_seen = 0
    groups_seen = 0

    def series_for_cid(cid):
        if cid is None:
            return np.zeros(len(labels), dtype=float)
        return np.array(
            [data.get(sid, {}).get(cid, {}).get(metric, 0.0) for sid in sids_in_order],
            dtype=float
        )

    # Hatches for successive stacked layers (beyond the base)
    hatch_cycle = ['//', '\\\\', 'xx', '..', '++', '**', 'oo']
    alpha_base = 1.0
    alpha_top = 0.35

    # ---- Draw each legend group as a single stacked bar per sample ----
    for lg in visual_legend_seq:
        if lg is None:
            gaps_seen += 1
            continue

        cids_by_method = legend_to_cids[lg]
        # Pick a consistent color per caller (from any method under that legend)
        any_cid = next(iter(cids_by_method.values()))
        color = cid_mapping[any_cid]["color"]

        offset = start_offset + groups_seen * width + gaps_seen * gap_width

        # Stack in global method order (only methods present for this legend)
        present_methods = [m for m in method_order if m in cids_by_method]

        prev_vals = None
        for i, m in enumerate(present_methods):
            vals = series_for_cid(cids_by_method[m])
            if i == 0:
                # Base layer
                ax.bar(
                    x + offset,
                    vals,
                    width,
                    label=lg,  # caller legend text (handles built separately if needed)
                    color=color,
                    linewidth=0
                )
                prev_vals = vals
            else:
                # Stack only the positive delta
                delta = np.maximum(vals - prev_vals, 0.0)
                ax.bar(
                    x + offset,
                    delta,
                    width,
                    bottom=prev_vals,
                    color=color,
                    alpha=alpha_top,
                    hatch=hatch_cycle[(i - 1) % len(hatch_cycle)],
                    linewidth=0
                )
                prev_vals = np.maximum(prev_vals, vals)

        groups_seen += 1

    # ---- Axes / labels ----
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)

    ylabel = metric.split('_')[1]
    ylabel = ylabel[0] + ylabel[1:].lower()
    ax.set_ylabel(ylabel, fontsize=16)

    # ---- Legends on the panel (optional) ----
    extra_artists = []
    if show_legends:
        # Methods legend (neutral)
        method_handles = build_method_handles(method_order)
        leg_methods = ax.legend(handles=method_handles, loc='upper left',
                                fontsize=10, frameon=False, title="")

        # Callers legend (colors)
        caller_handles = []
        caller_seen = set()
        for lg in real_groups:
            if lg in caller_seen:
                continue
            any_cid = next(iter(legend_to_cids[lg].values()))
            color = cid_mapping[any_cid]["color"]
            caller_handles.append(Patch(facecolor=color, edgecolor='none', label=lg))
            caller_seen.add(lg)

        if legend_outside:
            leg_callers = ax.legend(handles=caller_handles, title="",
                                    ncol=1, fontsize=10, frameon=False,
                                    loc='center left', bbox_to_anchor=(1.02, 0.5),
                                    borderaxespad=0.0)
            ax.add_artist(leg_methods)
            ax.add_artist(leg_callers)
            extra_artists = [leg_callers]
        else:
            leg_callers = ax.legend(handles=caller_handles, title="",
                                    ncol=1, fontsize=10, frameon=False, loc='upper right')
            ax.add_artist(leg_methods)
            ax.add_artist(leg_callers)

    os.makedirs(out_folder, exist_ok=True)
    save_kwargs = dict(bbox_inches='tight')
    if extra_artists:
        save_kwargs['bbox_extra_artists'] = extra_artists

    plt.savefig(os.path.join(out_folder, f"fig.{metric.lower()}.pdf"), **save_kwargs)
    plt.close(fig)


# --------------------------
# CLI
# --------------------------

def main():
    parser = argparse.ArgumentParser(
        description=("Generate stacked charts per caller. Methods are arbitrary strings from the mapping; "
                     "stacking order is the order of first appearance. Base = first method; higher methods show "
                     "only the positive delta over the previous. Optionally write a separate legend PDF.")
    )
    parser.add_argument("slabel_sid_path", help="Path to the file containing SLabel and SID mappings")
    parser.add_argument("clabel_cid_path", help="Path to the file containing CLabel and CID mappings (with METHOD label)")
    parser.add_argument("out_folder", help="Output folder for generated charts")
    parser.add_argument("--legend-outside", action="store_true",
                        help="Place the Callers legend outside the axes (right side) on each panel (ignored if --separate-legend).")
    parser.add_argument("--separate-legend", action="store_true",
                        help="Write a standalone legend PDF (fig.legend.pdf) and suppress legends on the metric panels.")
    args = parser.parse_args()

    parse_files(args.slabel_sid_path, args.clabel_cid_path, args.out_folder,
                legend_outside=args.legend_outside,
                separate_legend=args.separate_legend)


if __name__ == "__main__":
    main()

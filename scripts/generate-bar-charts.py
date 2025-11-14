import os
import argparse
import matplotlib.pyplot as plt
from collections import defaultdict, OrderedDict
import numpy as np

def read_mapping_file(mapping_file):
    folder_to_info = OrderedDict()
    with open(mapping_file, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:  # blank line = gap
                gap_id = f"GAP_{len(folder_to_info)}"
                folder_to_info[gap_id] = {"legend": None, "color": None, "is_gap": True}
                continue
            folder, color, legend_name = line.split(maxsplit=2)
            folder_to_info[folder] = {"legend": legend_name, "color": color, "is_gap": False}
    return folder_to_info

def parse_files(slabel_sid_path, clabel_cid_path, out_folder):
    with open(slabel_sid_path, 'r') as f:
        sid_mapping = {line.split()[0]: ' '.join(line.split()[1:]) for line in f}

    cid_mapping = read_mapping_file(clabel_cid_path)

    # Initialize data structure to hold results
    data = defaultdict(lambda: defaultdict(dict))  # data[SID][CID] = {metric: value}

    # Collect data from each CID/SID.txt file
    metrics_to_collect = ["DEL", "DUP", "INS", "INV"]
    for cid, info in cid_mapping.items():
        for sid in sid_mapping:
            file_path = f"{cid}/{sid}.txt"
            if os.path.exists(file_path):
                with open(file_path, 'r') as f:
                    for line in f:
                        for metric in metrics_to_collect:
                            if f"{metric} SENSITIVITY" in line:
                                value = float(line.split()[4].strip("()"))
                                data[sid][cid][f"{metric}_SENSITIVITY"] = value
                            elif f"{metric} PRECISION" in line:
                                value = float(line.split()[4].strip("()"))
                                data[sid][cid][f"{metric}_PRECISION"] = value

                    # Calculate F1 score for each metric
                    for metric in metrics_to_collect:
                        sensitivity_key = f"{metric}_SENSITIVITY"
                        precision_key = f"{metric}_PRECISION"
                        if sensitivity_key in data[sid][cid] and precision_key in data[sid][cid]:
                            sensitivity = data[sid][cid][sensitivity_key]
                            precision = data[sid][cid][precision_key]
                            if sensitivity + precision > 0:
                                f1_score = 2 * (sensitivity * precision) / (sensitivity + precision)
                            else:
                                f1_score = 0
                            data[sid][cid][f"{metric}_F1"] = f1_score

    # Generate plots for each metric
    metrics = ["DEL_SENSITIVITY", "DEL_PRECISION", "DEL_F1", "DUP_SENSITIVITY", "DUP_PRECISION", "DUP_F1", \
    "INS_SENSITIVITY", "INS_PRECISION", "INS_F1", "INV_SENSITIVITY", "INV_PRECISION", "INV_F1"]
    for metric in metrics:
        plot_barchart(data, metric, out_folder, cid_mapping, sid_mapping)

def plot_barchart(data, metric, out_folder, cid_mapping, sid_mapping, gap_frac=0.4):
    # Prepare data
    sids_in_order = list(data.keys())
    labels = [sid_mapping[sid] for sid in sids_in_order]
    cids = list(cid_mapping.keys())
    metric_data = {
        cid: [data[sid].get(cid, {}).get(metric, 0.0) for sid in sids_in_order]
        for cid in cids
    }

    # Count bars and gaps
    real_cids = [cid for cid in cids if not cid_mapping[cid].get("is_gap", False)]
    G = sum(1 for cid in cids if cid_mapping[cid].get("is_gap", False))
    N = len(real_cids)

    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(labels))

    # Base bar width (per real CID)
    width = 0.8 / max(N, 1)
    gap_width = gap_frac * width

    # Total span of the group (distance between first and last bar centers)
    total_span = (N - 1) * width + G * gap_width
    start_offset = - total_span / 2.0  # center group around x

    # Place each CID: advance by width for real bars, gap_width for gaps
    gaps_seen = 0
    real_seen = 0
    for cid in cids:
        info = cid_mapping[cid]
        if info.get("is_gap", False):
            gaps_seen += 1
            continue

        offset = start_offset + real_seen * width + gaps_seen * gap_width
        ax.bar(
            x + offset,
            metric_data[cid],
            width,
            label=info["legend"],
            color=info["color"]
        )
        real_seen += 1

    # Axes/labels
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    ylabel = metric.split('_')[1]
    ylabel = ylabel[0] + ylabel[1:].lower()
    ax.set_ylabel(ylabel, fontsize=16)

    # Legend policy (as before)
    if 'DEL' in metric and 'PRECISION' not in metric:
        ax.legend()

    os.makedirs(out_folder, exist_ok=True)
    plt.savefig(os.path.join(out_folder, f"fig.{metric.lower()}.pdf"), bbox_inches='tight')
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description="Generate sensitivity and precision charts from data files.")
    parser.add_argument("slabel_sid_path", help="Path to the file containing SLabel and SID mappings")
    parser.add_argument("clabel_cid_path", help="Path to the file containing CLabel and CID mappings")
    parser.add_argument("out_folder", help="Output folder for generated charts")
    
    args = parser.parse_args()
    parse_files(args.slabel_sid_path, args.clabel_cid_path, args.out_folder)

if __name__ == "__main__":
    main()

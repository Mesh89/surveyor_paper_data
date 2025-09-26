import os, re, argparse
import matplotlib.pyplot as plt
from collections import OrderedDict

def extract_metrics(file_path):
    metrics = {}
    with open(file_path, 'r') as file:
        content = file.read()
        
        def grab(pattern, default=0.0):
            m = re.search(pattern, content)
            if m and m.group(1).strip() != "":
                return float(m.group(1))
            return default

        # DEL
        ds = grab(r'DEL SENSITIVITY: \d+/\d+ = ([0-9.eE+-]*)')
        dp = grab(r'DEL PRECISION: \d+/\d+ = ([0-9.eE+-]*)')
        metrics['DEL'] = (ds, dp)

        # DUP
        us = grab(r'DUP SENSITIVITY: \d+/\d+ = ([0-9.eE+-]*)')
        up = grab(r'DUP PRECISION: \d+/\d+ = ([0-9.eE+-]*)')
        metrics['DUP'] = (us, up)

        # INS
        is_ = grab(r'INS SENSITIVITY: \d+/\d+ = ([0-9.eE+-]*)')
        ip = grab(r'INS PRECISION: \d+/\d+ = ([0-9.eE+-]*)')
        metrics['INS'] = (is_, ip)
    
    return metrics

def plot_metric(metric, data, output_prefix, legend_order):
    plt.figure(figsize=(8, 7))  # Width: 8, Height: 7 for more space

    # Loop through the sorted data and plot each point
    for legend_name, color in legend_order:
        print(f'Plotting {metric} for {legend_name} with color {color}')
        x = [p[0] for p in data[legend_name]]
        y = [p[1] for p in data[legend_name]]
        plt.scatter(x, y, label=legend_name, color=color)

    plt.xlabel('Sensitivity', fontsize=16)
    plt.ylabel('Precision', fontsize=16)
    plt.grid(True)

    # Ensure axes start at (0,0)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.tick_params(axis='both', which='major', labelsize=12)

    # Place the legend inside the plot
    if metric == 'DEL':
        plt.legend(loc="upper left", fancybox=True, shadow=True, fontsize=14)

    # Adjust the layout to make space for the legend and ensure nothing is cut off
    plt.tight_layout()

    output_file = f'{output_prefix}.{metric}.pdf'
    plt.savefig(output_file, bbox_inches='tight')  # Save the plot with the legend inside
    plt.close()
    print(f'Plot saved as {output_file}')


def read_mapping_file(mapping_file):
    folder_to_info = OrderedDict()
    with open(mapping_file, 'r') as file:
        for line in file:
            parts = line.strip().split(maxsplit=2)
            folder, color, legend_name = parts
            folder_to_info[folder] = {"legend": legend_name, "color": color}
    return folder_to_info

def process_directory(directory, output_prefix, callers_mapping):
    del_data = {}
    dup_data = {}
    ins_data = {}

    for caller, info in callers_mapping.items():
        subdir = os.path.join(directory, caller)  # Construct the full path
        if not os.path.isdir(subdir):
            raise FileNotFoundError(f"Subdirectory {subdir} does not exist.")

        legend_name = info["legend"]  # Get the custom legend name
        
        for file in os.listdir(subdir):
            if file.endswith(".txt"):
                file_path = os.path.join(subdir, file)
                metrics = extract_metrics(file_path)
                del_data.setdefault(legend_name, []).append(metrics['DEL'])
                dup_data.setdefault(legend_name, []).append(metrics['DUP'])
                ins_data.setdefault(legend_name, []).append(metrics['INS'])

    legend_order = [(info["legend"], info["color"]) for info in callers_mapping.values()]

    # Plot the data for DEL, DUP, INS, with a legend for each software
    if del_data:
        plot_metric('DEL', del_data, output_prefix, legend_order)
    if dup_data:
        plot_metric('DUP', dup_data, output_prefix, legend_order)
    if ins_data:
        plot_metric('INS', ins_data, output_prefix, legend_order)

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description='Plot Sensitivity vs Precision for DEL, DUP, and INS from .txt files.')
    parser.add_argument('directory', type=str, help='The directory containing subdirectories with .txt files.')
    parser.add_argument('output_prefix', type=str, help='The output file prefix for the generated plots.')
    parser.add_argument('mapping_file', type=str, help='A text file that maps folder names to legend names.')

    args = parser.parse_args()

    # Read the mapping file
    callers_mapping = read_mapping_file(args.mapping_file)

    # Process the directory and plot the data
    process_directory(args.directory, args.output_prefix, callers_mapping)

if __name__ == "__main__":
    main()

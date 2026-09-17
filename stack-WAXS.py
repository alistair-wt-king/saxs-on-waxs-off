# Usage - python3 stack-WAXS.py -range 0.5 3 -xy "Q (Å-1)" "Intensity (cm-1)" -rev -files LiBr.3H2O-Dry.xy blue LiBr.3H2O-Wet.xy Magenta NaOH-Urea-Dry.xy Green NaOH-Urea-Wet.xy Orange Pulp-Dry.xy Red
# bt Alistair WT King
import sys
import matplotlib.pyplot as plt
import os

def read_xy_file(filename):
    x_vals = []
    y_vals = []
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    x = float(parts[0])
                    y = float(parts[1])
                    x_vals.append(x)
                    y_vals.append(y)
                except ValueError:
                    continue
    return x_vals, y_vals

def plot_stacked_scatter(file_color_pairs, xlabel='X', ylabel='Y (stacked)', x_range=None):
    plt.figure(figsize=(10, 6))
    offset = 0
    max_y_total = 0
    min_y_total = 0

    for file, color in file_color_pairs:
        x, y = read_xy_file(file)
        y_offset = [val + offset for val in y]
        label = os.path.splitext(os.path.basename(file))[0]  # Remove .xy extension
        plt.scatter(x, y_offset, label=label, s=10, color=color)
        if y:
            offset += max(y) * 1.1
            max_y_total = max(max_y_total, max(y_offset))
            min_y_total = min(min_y_total, min(y_offset))

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if x_range:
        plt.xlim(x_range)
    else:
        plt.xlim(min_x_total, max_x_total * 1.2)
    plt.legend(loc='upper right')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

if __name__ == "__main__":
    args = sys.argv[1:]

    reverse = False
    xlabel = 'X'
    ylabel = 'Y (stacked)'
    x_range = None
    file_color_pairs = []

    # Parse -range flag
    if '-range' in args:
        range_index = args.index('-range')
        try:
            xmin = float(args[range_index + 1])
            xmax = float(args[range_index + 2])
            x_range = (xmin, xmax)
            del args[range_index:range_index + 3]
        except (IndexError, ValueError):
            print("Error: -range flag must be followed by two numeric values (xmin xmax).")
            sys.exit(1)

    # Parse -xy flag
    if '-xy' in args:
        xy_index = args.index('-xy')
        try:
            xlabel = args[xy_index + 1]
            ylabel = args[xy_index + 2]
            del args[xy_index:xy_index + 3]
        except IndexError:
            print("Error: -xy flag must be followed by two quoted axis labels.")
            sys.exit(1)

    # Parse -rev flag
    if '-rev' in args:
        reverse = True
        args.remove('-rev')

    # Parse -files flag
    if '-files' in args:
        files_index = args.index('-files')
        file_args = args[files_index + 1:]
        args = args[:files_index]
        if len(file_args) % 2 != 0:
            print("Error: Each file must be paired with a color.")
            sys.exit(1)
        file_color_pairs = [(file_args[i], file_args[i+1]) for i in range(0, len(file_args), 2)]
    else:
        print("Error: Missing -files flag with file-color pairs.")
        sys.exit(1)

    if reverse:
        file_color_pairs.reverse()

    plot_stacked_scatter(file_color_pairs, xlabel, ylabel, x_range)
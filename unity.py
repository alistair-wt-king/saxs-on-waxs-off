#!/usr/bin/env python3
# Written by Alistair King & Co-Pilot

# Script to combine SAXS & WAXS data regions, before or after correction
# Usage: python3 unity.py -is SAXS_data.xy -iw WAXS_data.xy -c 0.2 -log
# The -c flag and argument defines the q value where the SAXS data stops and WAXS starts
# log or normal plotting is possible


import argparse
import numpy as np
import matplotlib.pyplot as plt
import os


def load_xy_file(path):
    """Load (q, I) data from a 2-column .xy file."""
    q_vals, i_vals = [], []
    with open(path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            try:
                q_vals.append(float(parts[0]))
                i_vals.append(float(parts[1]))
            except ValueError:
                pass
    return np.array(q_vals), np.array(i_vals)


def write_xy_file(path, q, I):
    """Write combined q, I data to output file."""
    with open(path, "w") as f:
        for qq, ii in zip(q, I):
            f.write(f"{qq:.6e}  {ii:.6e}\n")


def coeff_str(v):
    return f"{v:.3f}".replace(".", "p")


def main():

    parser = argparse.ArgumentParser(description="Combine SAXS + WAXS into one .xy file")

    parser.add_argument("-is", dest="saxs_file", required=True,
                        help="Input SAXS .xy file")
    parser.add_argument("-iw", dest="waxs_file", required=True,
                        help="Input WAXS .xy file")
    parser.add_argument("-c", dest="cutoff", required=True, type=float,
                        help="q-cutoff where SAXS ends and WAXS begins")
    parser.add_argument("-log", action="store_true",
                        help="Plot log10(Intensity) instead of linear")

    args = parser.parse_args()

    # -------------------------------------------------------
    # Load both datasets
    # -------------------------------------------------------
    q_saxs, I_saxs = load_xy_file(args.saxs_file)
    q_waxs, I_waxs = load_xy_file(args.waxs_file)

    # Sort data
    idx = np.argsort(q_saxs); q_saxs, I_saxs = q_saxs[idx], I_saxs[idx]
    idx = np.argsort(q_waxs); q_waxs, I_waxs = q_waxs[idx], I_waxs[idx]

    # -------------------------------------------------------
    # Apply cutoff
    # -------------------------------------------------------
    cutoff = args.cutoff

    mask_saxs = q_saxs <= cutoff
    mask_waxs = q_waxs >= cutoff

    Q_combined = np.concatenate([q_saxs[mask_saxs], q_waxs[mask_waxs]])
    I_combined = np.concatenate([I_saxs[mask_saxs], I_waxs[mask_waxs]])

    # Sort final combined data
    idx = np.argsort(Q_combined)
    Q_combined = Q_combined[idx]
    I_combined = I_combined[idx]

    print(f"SAXS points used: {np.sum(mask_saxs)}")
    print(f"WAXS points used: {np.sum(mask_waxs)}")
    print(f"Total combined points: {len(Q_combined)}")

    # -------------------------------------------------------
    # Save combined data
    # -------------------------------------------------------
    out_file = (
        os.path.splitext(args.saxs_file)[0]
        + "_combined_cut"
        + coeff_str(cutoff)
        + ".xy"
    )

    write_xy_file(out_file, Q_combined, I_combined)
    print("Saved combined file:", out_file)

    # -------------------------------------------------------
    # Plot result
    # -------------------------------------------------------
    plt.close("all")
    plt.figure(figsize=(7,5))

    if args.log:
        y = np.log10(np.clip(I_combined, 1e-30, None))
        plt.plot(Q_combined, y, label="Combined (log10)")
        plt.ylabel("log10(I)")
    else:
        plt.plot(Q_combined, I_combined, label="Combined")
        plt.ylabel("Intensity")

    plt.xlabel("q")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()

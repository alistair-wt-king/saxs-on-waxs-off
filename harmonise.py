#!/usr/bin/env python3
# Written by Alistair King & Co-Pilot

# Script to smooth oscillation errors caused by poor background fitting in the WAXS region
# Usage: python3 harmonise.py -i Combined_SAXS-WAXS_data.xy -qrange-smooth 0.3 0.45 -log
# Boundary-based polynomial reconstruction is used, i.e., fitting uses an extended
# Range of forward and backward points, outside of anchor points, to avoid fitting errors

import argparse
import numpy as np
import matplotlib.pyplot as plt
import os


def load_xy_file(path):
    """Load q, I from 2-column ASCII file."""
    q_vals, i_vals = [], []
    with open(path, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                q_vals.append(float(parts[0]))
                i_vals.append(float(parts[1]))
            except Exception:
                pass
    return np.array(q_vals), np.array(i_vals)


def write_xy_file(path, q, I):
    """Save q, I to output file."""
    with open(path, "w") as f:
        for qq, ii in zip(q, I):
            f.write(f"{qq:.6e}  {ii:.6e}\n")


def smooth_region_boundary_aware(q, I, qmin, qmax, degree=3, n_anchor=10):
    """
    Strong polynomial smoothing between qmin–qmax using anchor points
    outside the region for stable fitting.
    """

    print("\n=== Boundary-Based Reconstruction ===")
    print(f"Range: {qmin} – {qmax}")
    print(f"Degree: {degree}, Anchors each side: {n_anchor}")

    mask_smooth = (q >= qmin) & (q <= qmax)
    mask_left = q < qmin
    mask_right = q > qmax

    if np.sum(mask_smooth) < degree + 2:
        print("WARNING: Too few points in smoothing region. Returning original data.")
        return I.copy()

    # LEFT ANCHORS
    if np.sum(mask_left) >= n_anchor:
        left_q = q[mask_left][-n_anchor:]
        left_I = I[mask_left][-n_anchor:]
    else:
        left_q = q[mask_left]
        left_I = I[mask_left]

    # RIGHT ANCHORS
    if np.sum(mask_right) >= n_anchor:
        right_q = q[mask_right][:n_anchor]
        right_I = I[mask_right][:n_anchor]
    else:
        right_q = q[mask_right]
        right_I = I[mask_right]

    # Combine anchor points
    q_fit = np.concatenate([left_q, right_q])
    I_fit = np.concatenate([left_I, right_I])

    print(f"Anchor points used in fit: {len(q_fit)}")

    # Polynomial fit
    coeff = np.polyfit(q_fit, I_fit, degree)

    # Apply polynomial inside smoothing region
    q_in = q[mask_smooth]
    I_poly = np.polyval(coeff, q_in)

    # Replace smoothed region
    I_new = I.copy()
    I_new[mask_smooth] = I_poly

    return I_new


def main():
    parser = argparse.ArgumentParser(description="Boundary-based SAXS/WAXS reconstruction")
    parser.add_argument("-i", required=True, help="Input I(q) file")
    parser.add_argument("-qrange-smooth", nargs=2, type=float, required=True,
                        metavar=("QMIN", "QMAX"),
                        help="q-range to strongly smooth")
    parser.add_argument("-degree", type=int, default=3,
                        help="Polynomial degree for smoothing (default=3)")
    parser.add_argument("-anchors", type=int, default=10,
                        help="Anchor points on each side (default=10)")
    parser.add_argument("-log", action="store_true",
                        help="Plot log10(I)")
    args = parser.parse_args()

    qmin, qmax = args.qrange_smooth

    # Load data
    q, I = load_xy_file(args.i)
    idx = np.argsort(q)
    q, I = q[idx], I[idx]

    # Smooth
    I_new = smooth_region_boundary_aware(
        q, I, qmin, qmax,
        degree=args.degree,
        n_anchor=args.anchors
    )

    # Save
    outname = (
        os.path.splitext(args.i)[0]
        + f"_harmonised_{qmin:.3f}_{qmax:.3f}.xy"
    )
    write_xy_file(outname, q, I_new)
    print("Saved:", outname)

    # Plot
    plt.close("all")
    plt.figure(figsize=(7,5))

    if args.log:
        I_safe = np.where(I_new > 0, I_new, 1e-30)
        plt.plot(q, np.log10(I_safe), label="Harmonised (log10)")
        plt.ylabel("log10(I)")
    else:
        plt.plot(q, I_new, label="Harmonised")
        plt.ylabel("Intensity")

    plt.xlabel("q")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()

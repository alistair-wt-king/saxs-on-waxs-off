
# Written by Alistair King & Co-Pilot
#!/usr/bin/env python3

# Script to subtract sample holder and solvent diffraction intensity from raw SAXS & WAXS data
# Usage: python3 saxs-waxs-sub.py -i Sample_WAXS.dat 1.0 -H Kapton_WAXS.dat 0.4 -s H2O_WAXS.dat 0.7 -log
# Accepts co-efficients for multiplication of the subtraction intensities
# Plotting can be using -log10(I) but this suffers from loss of negative points if fitting is poor
# Alternatively, the plotting region can also be defined using a I range, e.g., -Irange -0.01 0.02
# Smoothing of the noisy high-q values (SAXS) can be performed using Adaptive LOESS or Savitzky-Golay
# E.g., -savitzky-golay 20  OR  -loess-adaptive 0.5

import argparse
import numpy as np
import matplotlib.pyplot as plt
import os

try:
    from scipy.signal import savgol_filter
    HAVE_SG = True
except ImportError:
    HAVE_SG = False


# =====================================================================
# FAST + SAFE LOESS (NumPy-only)
# =====================================================================
def loess_fast(x, y, frac):
    n = len(x)
    if n < 5:
        return y.copy()

    span = max(5, int(frac * n))
    span = min(span, n - 1)
    half = span // 2

    y_out = np.zeros_like(y)

    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half)
        xw = x[lo:hi]
        yw = y[lo:hi]

        if len(xw) < 3:
            y_out[i] = y[i]
            continue

        d = np.abs(xw - x[i])
        dmax = d.max() if len(d) > 0 else 1
        if dmax == 0:
            y_out[i] = yw[0]
            continue

        w = (1 - (d/dmax)**3)**3
        w = np.clip(w, 1e-6, None)

        X = np.column_stack((np.ones_like(xw), xw))
        W = np.diag(w)
        XtW = X.T @ W

        try:
            beta, *_ = np.linalg.lstsq(XtW @ X, XtW @ yw, rcond=None)
            y_out[i] = beta[0] + beta[1] * x[i]
        except Exception:
            y_out[i] = y[i]

    return y_out


# =====================================================================
# ADAPTIVE LOESS
# =====================================================================
def loess_adaptive(Q, I, strength=1.0):

    print("\n--- LOESS ADAPTIVE ---")
    print(f"strength = {strength}")
    print(f"Input pts: {len(Q)}")

    strength = max(0.2, min(strength, 3.0))
    low_frac  = min(0.8, max(0.02, 0.05 * strength))
    mid_frac  = min(0.8, max(0.02, 0.15 * strength))
    high_frac = min(0.8, max(0.02, 0.35 * strength))

    low_limit = 0.10
    mid_limit = 0.20

    mask_low  = Q < low_limit
    mask_mid  = (Q >= low_limit) & (Q < mid_limit)
    mask_high = Q >= mid_limit

    print(f"low region:  {np.sum(mask_low)}")
    print(f"mid region:  {np.sum(mask_mid)}")
    print(f"high region: {np.sum(mask_high)}")

    I_out = np.zeros_like(I)

    if np.sum(mask_low) > 5:
        I_out[mask_low] = loess_fast(Q[mask_low], I[mask_low], low_frac)
    else:
        I_out[mask_low] = I[mask_low]

    if np.sum(mask_mid) > 5:
        I_out[mask_mid] = loess_fast(Q[mask_mid], I[mask_mid], mid_frac)
    else:
        I_out[mask_mid] = I[mask_mid]

    if np.sum(mask_high) > 5:
        I_out[mask_high] = loess_fast(Q[mask_high], I[mask_high], high_frac)
    else:
        I_out[mask_high] = I[mask_high]

    print(f"Output pts: {len(I_out)}")
    print("--- LOESS DONE ---\n")

    return I_out


# =====================================================================
# FILE I/O
# =====================================================================
def load_xy_file(path):
    q_vals, i_vals = [], []
    with open(path, "r") as f:
        for line in f:
            p = line.split()
            if len(p) < 2:
                continue
            try:
                q_vals.append(float(p[0]))
                i_vals.append(float(p[1]))
            except:
                continue
    return np.array(q_vals), np.array(i_vals)


def write_xy_file(path, q, I):
    with open(path, "w") as f:
        for qq, ii in zip(q, I):
            f.write(f"{qq:.6e} {ii:.6e}\n")


def coeff_str(c):
    return f"{c:.3f}".replace(".", "p")


# =====================================================================
# MAIN
# =====================================================================
def main():

    parser = argparse.ArgumentParser(description="SAXS subtraction + smoothing DIAGNOSTIC")
    parser.add_argument("-i", nargs=2, required=True)
    parser.add_argument("-H", nargs=2, required=True)
    parser.add_argument("-s", nargs=2, required=True)
    parser.add_argument("-log", action="store_true")
    parser.add_argument("-Irange", nargs=2, type=float)
    parser.add_argument("-savitzky-golay", nargs="?", const="11")
    parser.add_argument("-loess-adaptive", nargs="?", const="1.0")
    args = parser.parse_args()

    # --------------------------------------------------------
    # Load
    # --------------------------------------------------------
    q_i, I_i = load_xy_file(args.i[0])
    q_h, I_h = load_xy_file(args.H[0])
    q_s, I_s = load_xy_file(args.s[0])

    print("\n--- LOAD ---")
    print("sample:", len(q_i))
    print("holder:", len(q_h))
    print("solvent:", len(q_s))

    # Sort
    q_i, I_i = q_i[np.argsort(q_i)], I_i[np.argsort(q_i)]
    q_h, I_h = q_h[np.argsort(q_h)], I_h[np.argsort(q_h)]
    q_s, I_s = q_s[np.argsort(q_s)], I_s[np.argsort(q_s)]

    # --------------------------------------------------------
    # Subtraction
    # --------------------------------------------------------
    print("\n--- SUBTRACTION ---")
    I_corr = (
        I_i * float(args.i[1])
        - np.interp(q_i, q_h, I_h) * float(args.H[1])
        - np.interp(q_i, q_s, I_s) * float(args.s[1])
    )

    Q = q_i.copy()
    I = I_corr.copy()

    print("After subtraction:", len(Q))
    print("Negative values:", np.sum(I < 0))
    print("NaNs:", np.sum(~np.isfinite(I)))

    # --------------------------------------------------------
    # Remove low-q
    # --------------------------------------------------------
    q_threshold = 0.0071
    mask = Q >= q_threshold
    Q, I = Q[mask], I[mask]

    print("\nAfter q-threshold:", len(Q))

    # --------------------------------------------------------
    # Only remove NaN / Inf
    # --------------------------------------------------------
    mask = np.isfinite(I)
    Q, I = Q[mask], I[mask]

    print("After finite filter:", len(Q))

    # --------------------------------------------------------
    # Smoothing
    # --------------------------------------------------------
    tag = ""

    if args.savitzky_golay:
        win = int(args.savitzky_golay)
        print("\n--- SG smoothing ---")
        print("window =", win)
        if HAVE_SG:
            I = savgol_filter(I, window_length=win, polyorder=3)
            tag = f"_SG{win}"
        else:
            print("SciPy missing")

    elif args.loess_adaptive:
        strength = float(args.loess_adaptive)
        print("\n--- LOESS smoothing ---")
        I = loess_adaptive(Q, I, strength)
        tag = f"_LOESS{strength}"

    print("After smoothing:", len(Q))

    # --------------------------------------------------------
    # Save file
    # --------------------------------------------------------
    out = (
        os.path.splitext(args.i[0])[0]
        + f"_i{coeff_str(float(args.i[1]))}"
        + f"_H{coeff_str(float(args.H[1]))}"
        + f"_s{coeff_str(float(args.s[1]))}"
        + tag
        + ".xy"
    )

    write_xy_file(out, Q, I)
    print("\nSaved:", out)

    # --------------------------------------------------------
    # PLOT
    # --------------------------------------------------------
    print("\n--- PLOT BLOCK ENTERED ---")
    print("Q length:", len(Q))
    print("I length:", len(I))

    if len(Q) == 0:
        print("\nERROR: No data to plot!")
        return

    plt.close("all")
    plt.figure(figsize=(7, 5))

    if args.log:
        y = np.log10(np.clip(I, 1e-30, None))
        plt.plot(Q, y, label="log10(I)")
        plt.ylabel("log10(I)")
    else:
        plt.plot(Q, I, label="Intensity")
        plt.ylabel("Intensity")
        if args.Irange:
            plt.ylim(args.Irange[0], args.Irange[1])

    plt.xlabel("q")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()

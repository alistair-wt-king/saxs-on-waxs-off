#!/usr/bin/env python3

# please use the appropriate flags for conversion between the values
# -2t2q for conversion from 2theta to q
# -q22t for conversion of q to 2theta

import sys
from pathlib import Path
import numpy as np

# Cu Kα wavelength (Å)
LAMBDA_XRAY = 1.5406


def load_xy_file(filename):
    """
    Load a two-column XY file while automatically
    skipping any non-numeric header lines.
    """

    skiprows = 0

    with open(filename, "r") as f:

        for line in f:

            line = line.strip()

            # Skip blank lines
            if not line:
                skiprows += 1
                continue

            fields = line.split()

            try:
                float(fields[0])
                float(fields[1])
                break

            except (ValueError, IndexError):
                skiprows += 1

    return np.loadtxt(filename, skiprows=skiprows)


def twotheta_to_q(two_theta):
    """
    Convert 2θ (degrees) to q (Å⁻¹)
    """

    theta_rad = np.deg2rad(two_theta / 2.0)

    q = (4.0 * np.pi / LAMBDA_XRAY) * np.sin(theta_rad)

    return q


def q_to_twotheta(q):
    """
    Convert q (Å⁻¹) to 2θ (degrees)
    """

    arg = q * LAMBDA_XRAY / (4.0 * np.pi)

    if np.any(np.abs(arg) > 1.0):
        raise ValueError(
            "One or more q values are too large for the selected wavelength."
        )

    theta_rad = np.arcsin(arg)

    two_theta = np.rad2deg(2.0 * theta_rad)

    return two_theta


def main():

    if len(sys.argv) != 3:

        print(
            "\nUsage:\n"
            "  python3 2theta2q.py -2t2q input.xy\n"
            "  python3 2theta2q.py -q22t input.xy\n"
        )

        sys.exit(1)

    mode = sys.argv[1]
    infile = Path(sys.argv[2])

    data = load_xy_file(infile)

    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(
            f"{infile} does not appear to contain two-column XY data."
        )

    x = data[:, 0]
    intensity = data[:, 1]

    if mode == "-2t2q":

        x_new = twotheta_to_q(x)

        outfile = infile.with_name(
            infile.stem + "_q.xy"
        )

    elif mode == "-q22t":

        x_new = q_to_twotheta(x)

        outfile = infile.with_name(
            infile.stem + "_2theta.xy"
        )

    else:

        print(f"Unknown option: {mode}")

        sys.exit(1)

    output = np.column_stack((x_new, intensity))

    np.savetxt(
        outfile,
        output,
        fmt="%.6f"
    )

    print(f"Conversion complete: {outfile}")


if __name__ == "__main__":
    main()
# saxs-on-waxs-off
Small scripting toolkit for processing SAXS &amp; WAXS data, in particular related to 'disordered' cellulosic samples

**2theta2q.py** is a script for conversion of xy data from q range into 2theta (with the -q22t flag) or from 2theta to q (with the -2t2q flag).

usage: python3 2theta2q.py -q22t LiBr-3H2O_FD_norm.dat

**stack-WAXS.py** is a script for stacking xy patterns for publication images using matplotlib. The range is set with '-range', xy axes titles are set with '-xy', stack order can be reversed with '-rev' and the colours are set with the matplotlib CSS colour palette names (https://matplotlib.org/stable/gallery/color/named_colors.html).

usage: python3 stack-WAXS.py -range 5 50 -xy "2Theta (degrees)" "Area Normalised Intensity" -rev -files NaOH-Urea_FD_norm_bs.xy rosybrown NaOH-Urea_ND_norm_bs-scaled.xy lightcoral NaOH-Urea_1Gaussian.xy deeppink NaOH-Urea_2Gaussian.xy magenta NaOH-Urea_3Gaussian.xy blueviolet NaOH-Urea_4Gaussian.xy blue

**unity.py** is a script for merging SAXS & WAXS data. The value at which the data is merged at is set with the '-c' flag. SAXS data is set with '-is' and WAXS data is set with '-iw'. Data is output as xy data but for visualisation the '-log' flag can be used for ease of viewing and outputting an image.

usage: python3 unity.py -is NaOH-Urea_FD_SAXS.dat -iw NaOH-Urea_FD_WAXS.dat -c 0.15 -log

**saxs-waxs-sub.py** is a script for subtraction of background data from experimental data, _e.g_., glass, kapton or water subtraction. Scaling factors can be used for extraction of the data can be used. There are flags for input data (-i), holder data (-h) and solvent data (-s). Scaling factors are put after the file input names. '-log' can be used for ease of visualisation although negative intensity points will not plot. This can be visually problematic in regions with high noise, _e.g_., at higher q in SAXS. In this case plot without '-log', the I range can be set with '-Irange' to focus on an I range of interest, or, smoothing can be performed: Adaptive LOESS, _e.g_., '-loess-adaptive 0.5', OR, Savitzky-Golay with _e.g_., '-savitzky-golay 20' (vary the numbers to change the intensity of smoothing). If you don't want these options for visualisation as log, then strip the negative points from the data. Smoothing is is only really needed if you need to smooth the data for visualisation in the SAXS range. It is not really of value for SAXS model fitting. The scripts are not ideal for pre-processing for SAXS model fitting as the error data is stripped, although this could be changed in the future, if needed.

usage: python3 saxs-waxs-sub.py -i NaOH-Urea_ND_SAXS_combined_cut0p150.xy 1.0 -H Kapton_SAXS_combined_cut0p150.xy 1.0 -s H2O_SAXS_combined_cut0p150_i1p000_H1p000_s0p000_harmonised_0.250_0.450.xy 0.9 -Irange -0.005 0.02

**harmonise.py** is a script for removal of noise or artefacts, over short q-ranges. This uses boundary-based polynomial reconstruction between anchor points, defined with the '-qrange-smooth' flag. This is not recommended to be used purely for cosmetic reasons.

usage: python3 harmonise.py -i H2O_SAXS_combined_cut0p150_i1p000_H1p000_s0p000.xy -qrange-smooth 0.25 0.45 -log

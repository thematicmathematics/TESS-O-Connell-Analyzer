"""
A fully automated, terminal-based pipeline for analyzing TESS light curves 
for the O'Connell effect (starspot evolution) using TIC IDs. #O’Connell, D. J. K. (1951). Publications of the Riverview College Observatory, 2(6), 85.
Features:
- Dynamically downloads data from MAST via Lightkurve using TIC ID.
- Interactive terminal prompts if command-line arguments are not provided.
- Automated light curve cleaning, normalization, and phase folding.
- Computes both global and local O'Connell effect metrics.
- Generates publication-ready analysis panels for each sector and an asymmetry evolution trend.
"""
from __future__ import annotations
import argparse
import math
import os
import sys
import warnings
from dataclasses import dataclass
from typing import Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import lightkurve as lk

warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=UserWarning)

@dataclass
class TargetParams:
    """Data class to store dynamic parameters for the target.
    Attributes:
        tic (str): The TIC ID of the target (e.g., "12345678").
        period (float): The orbital period of the system in days.
        epoch (float): The reference epoch (T0) in BJD.
        author (str): The data author/pipeline (e.g., "QLP").
    """
    tic: str
    period: float
    epoch: float
    author: str = "AUTO"

def parse_arguments() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Automated TESS O'Connell Effect Pipeline.")
    ap.add_argument("--tic", help="TIC ID of the target star (without 'TIC' prefix, e.g., '12345678').")
    ap.add_argument("--period", type=float, help="Orbital period in days.")
    ap.add_argument("--epoch", type=float, help="Reference epoch (T0) in BJD.")
    ap.add_argument("--author", default="AUTO", help="Data author/pipeline (default: AUTO).")
    ap.add_argument("--outdir", help="Output directory name.")
    return ap.parse_args()

def get_interactive_inputs(args: argparse.Namespace) -> TargetParams:
    if args.tic and args.period and args.epoch:
        return TargetParams(tic=args.tic, period=args.period, epoch=args.epoch, author=args.author)
    print("\n" + "="*50)
    print(" TESS O'CONNELL ASYMMETRY ANALYSIS SYSTEM ")
    print("="*50)
    print("Terminal arguments are missing. Please enter the information manually:\n")

    try:
        tic = input(" TIC Number (Enter only the digits, e.g., 12345678): ").strip()
        period_str = input("Period (in days, e.g., 3.1294323): ").strip()
        epoch_str = input("T0 Reference Time (in BJD format, e.g., 2452964.6091): ").strip()
        author = input("Data Source / Author [If you leave it blank, an automatic search will be performed (SPOC -> TESS-SPOC -> QLP)]: ").strip()
        if not author:
            author = "AUTO"
        return TargetParams(
            tic=tic,
            period=float(period_str),
            epoch=float(epoch_str),
            author=author
        )
    
    except ValueError:
        print("\nERROR: Period and T0 values must be numeric. The program is terminating.", file=sys.stderr)
        sys.exit(1)

    except KeyboardInterrupt:
        print("\n\nThe transaction was canceled by the user.", file=sys.stderr)
        sys.exit(0)

def calculate_metrics(
    lc: lk.LightCurve,
    folded_lc: lk.FoldedLightCurve,
    params: TargetParams
) -> Tuple[float, float, float, float, np.ndarray, np.ndarray,
           float, float, float, float, float, float]:
    
    t, f = lc.time.value, lc.flux.value
    t0_btjd = params.epoch - 2457000.0

    # --- ÖNCE ph_wrapped ve peak1/peak2 hesapla ---
    ph = folded_lc.phase.value   
    fl = folded_lc.flux.value    
    ph_norm    = ph / params.period            
    ph_wrapped = np.where(ph_norm < 0, ph_norm + 1.0, ph_norm)

    def find_true_maximum(p_arr, f_arr, search_lo, search_hi,
                          n_bins=100, smooth_hw=5, meas_w=0.05):
        edges = np.linspace(search_lo, search_hi, n_bins + 1)
        bin_centers = 0.5 * (edges[:-1] + edges[1:])
        bin_flux = np.full(n_bins, np.nan)
        for i in range(n_bins):
            mask_b = (p_arr >= edges[i]) & (p_arr < edges[i + 1])
            if np.sum(mask_b) >= 2:
                bin_flux[i] = np.nanmedian(f_arr[mask_b])
        valid = np.isfinite(bin_flux)
        if valid.sum() < 5:
            mid = (search_lo + search_hi) / 2.0
            raw_mask = (p_arr > mid - meas_w) & (p_arr < mid + meas_w)
            if np.any(raw_mask):
                return mid, np.nanmedian(f_arr[raw_mask]), \
                       np.nanstd(f_arr[raw_mask]) / np.sqrt(np.sum(raw_mask)), True
            return mid, 0.0, 0.0, True
        bin_flux[~valid] = np.interp(bin_centers[~valid], bin_centers[valid], bin_flux[valid])
        kernel_size = 2 * smooth_hw + 1
        smoothed = np.convolve(bin_flux, np.ones(kernel_size) / kernel_size, mode='same')
        peak_idx  = np.argmax(smoothed)
        true_peak = bin_centers[peak_idx]
        meas_mask = (p_arr > true_peak - meas_w) & (p_arr < true_peak + meas_w)
        if np.sum(meas_mask) < 3:
            meas_mask = (p_arr > true_peak - meas_w * 1.5) & (p_arr < true_peak + meas_w * 1.5)
        if np.sum(meas_mask) == 0:
            return true_peak, 0.0, 0.0, True
        fv = np.nanmedian(f_arr[meas_mask])
        fe = np.nanstd(f_arr[meas_mask]) / np.sqrt(np.sum(meas_mask))
        return true_peak, fv, fe, False

    def fit_minima(p_arr, f_arr, expected, window=0.05):
        if expected == 0.0:
            mask = (p_arr > 1.0 - window) | (p_arr < window)
            p_fit = np.where(p_arr[mask] > 0.9, p_arr[mask] - 1.0, p_arr[mask])
        else:
            mask = (p_arr > expected - window) & (p_arr < expected + window)
            p_fit = p_arr[mask]
        f_fit = f_arr[mask]
        if len(p_fit) > 10:
            try:
                z, cov = np.polyfit(p_fit, f_fit, 2, cov=True)
                a, b = z[0], z[1]
                if a > 0:
                    min_p = -b / (2 * a)
                    sigma_a = np.sqrt(cov[0, 0])
                    sigma_b = np.sqrt(cov[1, 1])
                    err_p = abs(min_p) * np.sqrt((sigma_b / b)**2 + (sigma_a / a)**2)
                    if abs(min_p - expected) < window:
                        res_p = min_p if min_p >= 0 else min_p + 1.0
                        return res_p, err_p
            except Exception:
                pass
        if np.any(mask):
            return p_arr[mask][np.argmin(f_fit)], 0.0
        return expected, 0.0

    min1_phase, min1_err = fit_minima(ph_wrapped, fl, 0.0)
    min2_phase, min2_err = fit_minima(ph_wrapped, fl, 0.5)

    shift1 = min1_phase if min1_phase < 0.5 else min1_phase - 1.0
    shift1_err = min1_err
    shift2 = min2_phase - 0.5
    shift2_err = min2_err

    eclipse_hw  = 0.15
    edge_margin = 0.03
    min1_ph = (1.0 + shift1) % 1.0
    min2_ph = (0.5 + shift2) % 1.0

    max1_lo = (min2_ph + eclipse_hw + edge_margin) % 1.0
    max1_hi = (min1_ph - eclipse_hw - edge_margin) % 1.0
    max2_lo = (min1_ph + eclipse_hw + edge_margin) % 1.0
    max2_hi = (min2_ph - eclipse_hw - edge_margin) % 1.0

    def _valid_range(lo, hi, fallback_lo, fallback_hi):
        span = (hi - lo) % 1.0
        if span < 0.05 or span > 0.60:
            return fallback_lo, fallback_hi
        return lo, hi

    max1_lo, max1_hi = _valid_range(max1_lo, max1_hi, 0.15, 0.40)
    max2_lo, max2_hi = _valid_range(max2_lo, max2_hi, 0.60, 0.85)

    peak1, g_m1, err_g1, fb1 = find_true_maximum(ph_wrapped, fl, max1_lo, max1_hi)
    peak2, g_m2, err_g2, fb2 = find_true_maximum(ph_wrapped, fl, max2_lo, max2_hi)

    if g_m1 != 0 and g_m2 != 0 and not (fb1 and fb2):
        global_diff = ((g_m2 - g_m1) / g_m1) * 100
        if g_m2 != 0:
            global_diff_err = abs(
                (g_m2 / g_m1) * np.sqrt((err_g2 / g_m2)**2 + (err_g1 / g_m1)**2)
            ) * 100
        else:
            global_diff_err = 0.0
    else:
        global_diff, global_diff_err = 0.0, 0.0
        peak1, peak2 = 0.25, 0.75

    n_range = range(
        math.floor((t.min() - t0_btjd) / params.period),
        math.ceil((t.max() - t0_btjd) / params.period)
    )

    m1_vals, m2_vals = [], []
    for n in n_range:
        tm1 = t0_btjd + n * params.period + peak1 * params.period
        tm2 = t0_btjd + n * params.period + peak2 * params.period
        w = 0.05 * params.period

        mask1 = (t >= tm1 - w) & (t <= tm1 + w)
        mask2 = (t >= tm2 - w) & (t <= tm2 + w)

        if np.any(mask1):
            m1_vals.append([np.nanmean(t[mask1]), np.nanmedian(f[mask1]), np.nanstd(f[mask1]) / np.sqrt(np.sum(mask1))])
        if np.any(mask2):
            m2_vals.append([np.nanmean(t[mask2]), np.nanmedian(f[mask2]), np.nanstd(f[mask2]) / np.sqrt(np.sum(mask2))])

    m1_pts, m2_pts = np.array(m1_vals), np.array(m2_vals)

    if m1_pts.size > 0 and m2_pts.size > 0:
        mean_m1, mean_m2 = np.mean(m1_pts[:, 1]), np.mean(m2_pts[:, 1])
        err_m1 = np.sqrt(np.sum(m1_pts[:, 2]**2)) / len(m1_pts)
        err_m2 = np.sqrt(np.sum(m2_pts[:, 2]**2)) / len(m2_pts)
        local_diff = ((mean_m2 - mean_m1) / mean_m1) * 100 if mean_m1 != 0 else 0.0
        if mean_m1 != 0 and mean_m2 != 0:
            local_diff_err = abs(
                (mean_m2 / mean_m1) * np.sqrt((err_m2 / mean_m2)**2 + (err_m1 / mean_m1)**2)
            ) * 100
        else:
            local_diff_err = 0.0
    else:
        local_diff, local_diff_err = 0.0, 0.0

    return (global_diff, global_diff_err, local_diff, local_diff_err,
            m1_pts, m2_pts, shift1, shift1_err, shift2, shift2_err,
            peak1, peak2)

def create_analysis_panel(
    lc: lk.LightCurve,
    folded: lk.FoldedLightCurve,
    metrics: Tuple,
    out_dir: str,
    tic_id: str
) -> None:
    
    g_diff, g_err, l_diff, l_diff_err, m1_pts, m2_pts, \
        s1, s1_err, s2, s2_err, peak1, peak2 = metrics
    
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2)

    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(lc.time.value, lc.flux.value, '.', color='black', alpha=0.6, markersize=2, zorder=1, label='Data')

    if m1_pts.size > 0:
        ax1.scatter(m1_pts[:, 0], m1_pts[:, 1], c='dodgerblue', s=90, 
                    edgecolors='white', linewidth=1.2, zorder=5, 
                    label=f'Max I  (ph={peak1:.3f})')
    if m2_pts.size > 0:
        ax1.scatter(m2_pts[:, 0], m2_pts[:, 1], c='limegreen', s=90, 
                    edgecolors='white', linewidth=1.2, zorder=5, 
                    label=f'Max II (ph={peak2:.3f})')
                    
    ax1.set_title(f"TIC {tic_id} - Sector {lc.sector} Time Series (Local Diff: {l_diff:.4f}%)", 
                  fontsize=18, fontweight='bold')
    ax1.set_ylabel("Flux", fontsize=18, fontweight='bold')
    ax1.set_xlabel("Time [BTJD]",fontsize=18, fontweight='bold')
    ax1.tick_params(axis='both', which='major', labelsize=16)
    ax1.legend(loc='lower right', framealpha=0.9, edgecolor='black')

    ax2 = fig.add_subplot(gs[1, 0])
    try:
        _period = float(folded.period.value)
    except Exception:
        _period = 1.0

    ph_n = folded.phase.value
    fl_n = folded.flux.value
    ph_norm_plot = ph_n / _period
    ph_wrap_plot = np.where(ph_norm_plot < 0, ph_norm_plot + 1.0, ph_norm_plot)

    w_show = 0.04
    f1_show = fl_n[(ph_wrap_plot > peak1 - w_show) & (ph_wrap_plot < peak1 + w_show)]
    f2_show = fl_n[(ph_wrap_plot > peak2 - w_show) & (ph_wrap_plot < peak2 + w_show)]

    v1 = np.nanmedian(f1_show) if len(f1_show) > 0 else np.nan
    v2 = np.nanmedian(f2_show) if len(f2_show) > 0 else np.nan
    norm_factor = max(v1, v2) if np.isfinite(v1) and np.isfinite(v2) else 1.0

    try:
        binned = folded.bin(time_bin_size=0.01)
        ax2.scatter(binned.phase.value, binned.flux.value / norm_factor,
                    c='black', s=8, alpha=0.9, zorder=2, label='Binned Data')
    except Exception:
        pass 

    p1_jd = (peak1 if peak1 <= 0.5 else peak1 - 1.0) * _period
    p2_jd = (peak2 if peak2 <= 0.5 else peak2 - 1.0) * _period

    if np.isfinite(v1):
        ax2.axvline(p1_jd, color='blue',  ls='--', lw=1.2, alpha=0.7, label=f'Max I  ph={peak1:.3f}')
        ax2.scatter([p1_jd], [v1 / norm_factor], c='blue',  s=80, zorder=5)
    if np.isfinite(v2):
        ax2.axvline(p2_jd, color='green', ls='--', lw=1.2, alpha=0.7, label=f'Max II ph={peak2:.3f}')
        ax2.scatter([p2_jd], [v2 / norm_factor], c='green', s=80, zorder=5)

    ax2.legend(fontsize=10)
    ax2.set_title(f"Phase Folded - (Global Diff: {g_diff:.4f}%)", fontsize=20, fontweight='bold')
    ax2.set_xlabel("Phase", fontsize=18, fontweight='bold')
    ax2.set_ylabel("Normalized Flux", fontsize=18, fontweight='bold')
    ax2.tick_params(axis='both', which='major', labelsize=16)

    ax3 = fig.add_subplot(gs[1, 1])
    ax3.axis('off')

    snr = abs(g_diff) / g_err if g_err > 0 else 0
    abs_asym = abs(g_diff)
    if abs_asym < 0.5:                               #Ricker, G. R., et al. (2015). Journal of Astronomical Telescopes, Instruments, and Systems, 1(1), 014003. 
        verdict = "No Significant Spot (<0.5% Asym)"
    else:
        if snr < 3:
            verdict = "Insignificant Signal (<3σ)" #Bevington, P. R., & Robinson, D. K. (2003). Data Reduction and Error Analysis for the Physical Sciences.
        elif snr < 5:
            verdict = "Marginal Spot Activity (3σ-5σ)"
        else:
            verdict = "Confirmed Spot Activity (>5σ)"

    text = (
        f"==== SECTOR {lc.sector:02d} ANALYSIS REPORT ====\n\n"
        f"[ O'Connell Effect ]\n"
        f"  Global Diff: {g_diff:+.4f} % (± {g_err:.4f})\n"
        f"  Local Diff : {l_diff:+.4f} % (± {l_diff_err:.4f})\n"
        f"  SNR        : {snr:.2f}\n\n"
        f"[ Peak Phases ]\n"
        f"  Maximum I  : {peak1:.4f}\n"
        f"  Maximum II : {peak2:.4f}\n\n"
        f"[ Minima Shifts ]\n"
        f"  Primary     : {s1:+.4f} ph\n"
        f"  Secondary   : {s2:+.4f} ph\n\n"
        f"Status: {verdict}"
    )

    ax3.text(0.10, 0.99, text, fontsize=20, fontweight='bold',family='monospace',
             bbox=dict(facecolor='#f4f4f4', edgecolor='black', alpha=0.9, boxstyle='round,pad=0.2'),
             va='top', transform=ax3.transAxes)
    
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"S{lc.sector:02d}_Evidence.png"), dpi=330)
    plt.savefig(os.path.join(out_dir, f"S{lc.sector:02d}_Evidence.eps"))
    plt.close(fig)

def run_pipeline() -> None:
    args = parse_arguments()
    params = get_interactive_inputs(args)
    target_search = f"TIC {params.tic}"
    out_dir = args.outdir if args.outdir else f"TIC_{params.tic}_Analysis"
    os.makedirs(out_dir, exist_ok=True)

    print(f"\n Searching for data from {target_search} on MAST servers...")
    authors_to_try = ["SPOC", "TESS-SPOC", "QLP"] if params.author == "AUTO" else [params.author]
    search_result = None

    for auth in authors_to_try:
        print(f"  Trying author: {auth}...")
        try:
            temp_result = lk.search_lightcurve(target_search, author=auth)
            if len(temp_result) > 0:
                search_result = temp_result
                params.author = auth
                print(f" Data found in {auth}!")
                break
        except Exception as e:
            print(f"  ⚠️ Warning while searching {auth}: {e}")

    if search_result is None or len(search_result) == 0:
        print(f"❌ ERROR: No light curve data could be found for {target_search} in any common pipeline.", file=sys.stderr)
        sys.exit(1)
        
    print(f" {len(search_result)} sector(s) found. Downloading (This may take some time)...")
    
    try:
        lcs = search_result.download_all()
    except Exception as e:
        print(f"❌ ERROR: Data could not be downloaded. Check your internet connection or the MAST servers. Details: {e}", file=sys.stderr)
        sys.exit(1)

    if lcs is None:
        print("ERROR: Data could not be downloaded. Check your internet connection or the MAST servers.",
              file=sys.stderr)
        sys.exit(1)

    summary_list = []
    print("\n Data is being processed and panels are being created...")

    for lc in lcs:
        try:
            lc_clean = lc.remove_nans().remove_outliers(sigma_lower=100, sigma_upper=3).normalize()
            try:
                quality_mask = lc_clean.quality == 0
                lc_quality = lc_clean[quality_mask]
                if len(lc_quality) < 50:
                    print(f"    ⚠️  Sector {lc_clean.sector}: quality==0 filter "
                          f"left too few points ({len(lc_quality)}); filter skipped.")
                    lc_quality = lc_clean
            except Exception:
                lc_quality = lc_clean

            epoch_btjd = params.epoch - 2457000.0
            folded = lc_quality.fold(period=params.period, epoch_time=epoch_btjd)
            metrics = calculate_metrics(lc_clean, folded, params)

            g_diff, g_err, l_diff, l_diff_err, m1_pts, m2_pts, \
                s1, s1_err, s2, s2_err, peak1, peak2 = metrics
            
            create_analysis_panel(lc_clean, folded, metrics, out_dir, params.tic)

            summary_list.append({
                'Sector': lc_clean.sector,
                'Global_Diff': round(float(g_diff), 4),  
                'Global_Err': round(float(g_err), 4),    
                'Local_Diff': round(float(l_diff), 4),
                'Error_Pct': round(float(l_diff_err), 4),
                'Min1_Shift': round(float(s1), 4),
                'Min1_Err': round(float(s1_err), 4),
                'Min2_Shift': round(float(s2), 4),
                'Min2_Err': round(float(s2_err), 4),
                'Time_Mean': round(float(np.mean(lc_clean.time.value)), 4)
            })

            print(f"### Sector {lc_clean.sector} completed : \n"
                  f"(Global: {g_diff:.3f}% | Shifts: {s1:+.3f}, {s2:+.3f})")
        except Exception as e:
            sector_label = getattr(lc, 'sector', '?')
            print(f"    ⚠️  Sector {sector_label} skipped due to an error: {e}", file=sys.stderr)
            continue

    if len(summary_list) > 0:
        df = pd.DataFrame(summary_list).sort_values('Sector')
        csv_file_path = os.path.join(out_dir, f"TIC_{params.tic}_OConnell_Report.csv")
        df.to_csv(csv_file_path, index=False)
        print(f" The results were saved as a CSV file: {csv_file_path}")
    if len(df) > 1:
        plt.figure(figsize=(12, 6))
        x_indices = np.arange(len(df)) 
        plt.axhline(0, color='red', ls='--')
        plt.errorbar(x_indices, df['Global_Diff'], yerr=df['Global_Err'],
             fmt='D-', ms=6, lw=2, color='purple',
             ecolor='purple', elinewidth=1.5, capsize=5,
             label='Asymmetry Trend')
        
        for i, y_val in enumerate(df['Global_Diff']):
            plt.annotate(f"{y_val:.2f}%", 
                         (x_indices[i], y_val), 
                         textcoords="offset points", 
                         xytext=(0, 15), 
                         ha='center', 
                         fontweight='bold',
                         fontsize=12)
            
        plt.title(f"TIC {params.tic} - Asymmetry Graph", fontsize=14, fontweight='bold')
        plt.xlabel("TESS Sector Number")
        plt.ylabel("The O'Connell Difference (Max II - Max I) %")
        plt.grid(True, alpha=0.3)
        plt.xticks(x_indices, df['Sector'], rotation=0)   
        y_min, y_max = plt.ylim()
        plt.ylim(y_min - 0.25, y_max + 0.35)
        plt.tight_layout()

        trend_file = os.path.join(out_dir, f"TIC_{params.tic}_Asymmetry_Evolution.png")
        plt.savefig(trend_file, dpi=330)
        trend_file_eps = os.path.join(out_dir, f"TIC_{params.tic}_Asymmetry_Evolution.eps")
        plt.savefig(trend_file_eps)
        plt.close()

    print(f"\n Operation completed successfully! "
          f"All evidence panels have been saved to the '{out_dir}' folder.")
    
if __name__ == "__main__":
    run_pipeline()

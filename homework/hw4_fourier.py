"""
HW4 – Fourier Analysis of an Eclipsing Binary Light Curve

Loads a TESS light curve (FITS file) for an eclipsing binary star, selects a
single well-sampled observing epoch, and computes its Fourier transform and
power spectrum. The light curve is then reconstructed with the inverse
transform using as few Fourier coefficients as possible to find the minimum
needed to capture the eclipsing behavior. The analysis is repeated on a second
epoch to check that the results are consistent.

Because TESS cadence is nearly even but has gaps from missing observations,
the gaps are filled by linear interpolation onto a uniform time grid and the
analysis is redone to see how much the gaps affect the result.

Data source: https://tessebs.villanova.edu/

Usage examples:
    python hw4.py                                  # full light curve, default epoch
    python hw4.py --epoch 2750 2765 --plot ft      # FT magnitude of a different epoch
    python hw4.py --plot power                     # power spectrum of default epoch
    python hw4.py --plot all                       # all three figures

DISCLAIMER: Used more AI (Claude Code) than normal on this one to help with plot design
"""

# Star File Name: 0007849727.fits

import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from numpy.fft import rfft, rfftfreq, irfft

DEFAULT_EPOCH = (2720.0, 2735.0)   # inside the valid 2710-2800 range, gap at 2732-2735 


# func to read the data file FITS via hdul
def load_data(filename="0007849727.fits"):
    
    here = Path(__file__).parent
    with fits.open(here / filename) as hdul:
        times = np.array(hdul[1].data['times'], dtype=float)
        fluxes = np.array(hdul[1].data['fluxes'], dtype=float)
    good = np.isfinite(times) & np.isfinite(fluxes)
    return times[good], fluxes[good]

# plot of the full spectrum, used to find data points (epochs) of star data
def plot_full(times, fluxes, epoch):
    
    plt.figure(figsize=(12, 5))
    plt.plot(times, fluxes, '.', ms=2, color='tab:blue', label='TESS flux')
    plt.axvspan(epoch[0], epoch[1], color='tab:orange', alpha=0.3,
                label=f'chosen epoch ({epoch[0]:g}-{epoch[1]:g})')
    plt.xlabel("Time (days)")
    plt.ylabel("Flux")
    plt.title("Full Light Curve: 0007849727")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()


# func to do the FT of the data, using library in np 
def do_fft(times, fluxes, epoch):
    
    mask = (times >= epoch[0]) & (times <= epoch[1])
    t = times[mask]
    f = fluxes[mask]
    if len(t) < 10:
        raise ValueError(f"Only {len(t)} points in {epoch}; pick a denser range.")

    c = rfft(f - f.mean())      # subtract mean so the k=0 spike doesn't dominate
    freq = rfftfreq(len(f), d=np.median(np.diff(t)))   # cycles/day
    return t, f, freq, c

# plot of the epoch data to see holes or gaps in chosen range of epoch
def plot_epoch(t, f, epoch):
    
    plt.figure(figsize=(10, 4))
    plt.plot(t, f, '.-', ms=3, lw=0.5, color='tab:blue', label='selected epoch')
    plt.xlabel("Time (days)")
    plt.ylabel("Flux")
    plt.title(f"Selected Epoch: {epoch[0]:g} to {epoch[1]:g}")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

# plot of the transformed data, shows frequency peaks
def plot_ft(freq, c):
    
    plt.figure(figsize=(10, 4))
    plt.plot(freq, np.abs(c), '-', color='tab:green', label=r'$|c_k|$')
    plt.xlabel("Frequency (cycles/day)")
    plt.ylabel("Amplitude |c_k|")
    plt.title("Fourier Transform (magnitude)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

# plot of power spectrum of the FT
def plot_power(freq, c, n_harmonics=15):
    
    power = np.abs(c) ** 2
    k_peak = np.argmax(power[1:]) + 1            # skip k=0
    f0 = freq[k_peak]                            # strongest frequency

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.vlines(freq, power.min(), power, color='tab:red', lw=1)           
    ax.semilogy(freq, power, 'o', ms=3, color='tab:red', label=r'$|c_k|^2$')
    ax.axvline(f0, color='k', ls='--', lw=1,
               label=f'peak: {f0:.3f} /day (P = {1/f0:.3f} d)')

    # mark integer multiples of the peak frequency (harmonics)
    for n in range(2, n_harmonics + 1):
        ax.axvline(n * f0, color='gray', ls=':', lw=0.8,
                   label='harmonics (n x peak)' if n == 2 else None)

    ax.set_xlim(0, (n_harmonics + 1) * f0)       # zoom on the interesting region
    ax.set_yscale('log')
    ax.set_xlabel("Frequency (cycles/day)")
    ax.set_ylabel("Power |c_k|^2")
    ax.set_title("Power Spectrum (zoomed to harmonics)")
    ax.grid(True, alpha=0.3, which='both')
    ax.legend()
    fig.tight_layout()


# plot of FT data back to constructing original signal
def plot_recon(t, f, c, ks):
    
    fig, axes = plt.subplots(len(ks), 1, figsize=(10, 2.4 * len(ks)), sharex=True)
    axes = np.atleast_1d(axes)

    print("\n  K    rms residual   (smaller = better fit)")
    for ax, K in zip(axes, ks):
        c_cut = np.zeros_like(c)
        c_cut[:K + 1] = c[:K + 1]                 # keep k = 0..K, zero the rest
        recon = irfft(c_cut, n=len(f)) + f.mean() # back to time domain, add mean back
        rms = np.sqrt(np.mean((f - recon) ** 2))
        print(f"{K:3d}    {rms:.6f}")

        ax.plot(t, f, '.', ms=3, color='gray', label='data')
        ax.plot(t, recon, '-', lw=1.8, color='tab:red',
                label=f'inverse FFT, K = {K}  (rms = {rms:.4f})')
        ax.set_ylabel("Flux")
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right')

    axes[-1].set_xlabel("Time (days)")
    axes[0].set_title("Reconstruction from the first K Fourier coefficients")
    fig.tight_layout()


def parse_args():
    
    parser = argparse.ArgumentParser(
        description="Fourier analysis of a TESS eclipsing binary light curve.")
    parser.add_argument(
        "--epoch", nargs=2, type=float, metavar=("START", "END"),
        default=DEFAULT_EPOCH,
        help="Time range to analyze. Valid data spans 2710 to 2800. "
             "Recommended ranges: 2720 2735 (default), 2750 2765, 2775 2790. "
             "Pick ~10-20 days that are densely sampled.")
    parser.add_argument(
        "--plot", choices=["full", "ft", "power", "recon", "all"], default="full",
        help="full = whole TESS light curve; ft = Fourier transform magnitude; "
             "power = power spectrum of the FT; recon = inverse FFT using few "
             "coefficients; all = every figure (default: full)")
    parser.add_argument(
        "--ks", nargs="+", type=int, default=[1, 3, 5, 10, 20],
        help="Numbers of Fourier coefficients K to keep for --plot recon "
             "(default: 1 3 5 10 20). Try smaller values to find the minimum.")
    return parser.parse_args()


def main():
    args = parse_args()
    epoch = tuple(args.epoch)
    if not (2710 <= epoch[0] < epoch[1] <= 2800):
        raise SystemExit("Epoch must satisfy 2710 <= START < END <= 2800.")

    times, fluxes = load_data()

    if args.plot in ("full", "all"):
        plot_full(times, fluxes, epoch)

    if args.plot in ("ft", "power", "recon", "all"):
        t, f, freq, c = do_fft(times, fluxes, epoch)
        if args.plot == "all":
            plot_epoch(t, f, epoch)
        if args.plot in ("ft", "all"):
            plot_ft(freq, c)
        if args.plot in ("power", "all"):
            plot_power(freq, c)
        if args.plot in ("recon", "all"):
            plot_recon(t, f, c, args.ks)

    plt.show()


if __name__ == "__main__":
    main()
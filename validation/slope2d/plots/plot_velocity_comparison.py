#!/usr/bin/env python3
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt


def load_csv_flexible(path):
    # Read header
    with open(path, 'r') as f:
        header = f.readline().strip()
    headers = [h.strip() for h in header.split(',')]
    name_to_idx = {name: i for i, name in enumerate(headers)}

    # Determine indices for required fields if present
    req_names = ['time', 'v_tangent', 'v_theory', 'error_percent']
    idx = {n: name_to_idx.get(n, None) for n in req_names}

    t_list, vtan_list, vthe_list, err_list = [], [], [], []
    with open(path, 'r') as f:
        next(f)  # skip header
        for line in f:
            if not line.strip():
                continue
            parts = [p.strip() for p in line.split(',')]
            # safely fetch by index
            def getf(i):
                try:
                    return float(parts[i]) if i is not None and i < len(parts) and parts[i] != '' else np.nan
                except ValueError:
                    return np.nan
            t = getf(idx['time'])
            vtan = getf(idx['v_tangent'])
            vthe = getf(idx['v_theory'])
            err = getf(idx['error_percent'])
            t_list.append(t)
            vtan_list.append(vtan)
            vthe_list.append(vthe)
            err_list.append(err)

    t = np.asarray(t_list, dtype=float) # time [s]
    vtan = np.asarray(vtan_list, dtype=float) # velocity [m/s] (numeric)
    vthe = np.asarray(vthe_list, dtype=float) # velocity [m/s] (theory)
    err = np.asarray(err_list, dtype=float) # relative error [%]

    # Compute error if not present or invalid
    if np.all(~np.isfinite(err)):
        with np.errstate(divide='ignore', invalid='ignore'):
            err = np.abs(vtan - vthe) / np.where(np.abs(vthe) > 1e-12, np.abs(vthe), 1.0) * 100.0
    else:
        bad = ~np.isfinite(err)
        if np.any(bad):
            fix = np.abs(vtan - vthe) / np.where(np.abs(vthe) > 1e-12, np.abs(vthe), 1.0) * 100.0
            err[bad] = fix[bad]

    return t, vtan, vthe, err


def plot_comparison(csv_path, out_png):
    t, vtan, vthe, err = load_csv_flexible(csv_path)

    # mask non-finite
    m = np.isfinite(t) & np.isfinite(vtan) & np.isfinite(vthe)
    t, vtan, vthe = t[m], vtan[m], vthe[m]
    err = err[m]

    fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

    axs[0].plot(t, vtan, label='v_tangent (numeric)', color='tab:blue', lw=1.8)
    axs[0].plot(t, vthe, label='v_theory', color='tab:red', lw=1.5, ls='--')
    axs[0].set_ylabel('Velocity [m/s]')
    axs[0].grid(True, ls=':')
    axs[0].legend()
    # axs[0].set_xlim(0, 0.0002)
    # axs[0].set_ylim(-0.01, 0.0)

    axs[1].plot(t, err, label='relative error [%]', color='tab:green', lw=1.5)
    axs[1].set_xlabel('Time [s]')
    axs[1].set_ylabel('Error [%]')
    # axs[1].set_ylim(-0.1, 0.1)
    axs[1].grid(True, ls=':')
    # axs[1].set_xlim(0, 0.002)
    
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    fig.savefig(out_png, dpi=150)
    print(f'Saved: {out_png}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', default='./data/pem_slope_trace.csv', help='input CSV path')
    ap.add_argument('--out', default='./plots/velocity_comparison.png', help='output PNG path')
    args = ap.parse_args()
    plot_comparison(args.csv, args.out)


if __name__ == '__main__':
    main()

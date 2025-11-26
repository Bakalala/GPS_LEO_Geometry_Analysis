# pip install gnss-lib-py matplotlib numpy
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timezone, timedelta

import gnss_lib_py as glp

# =========================
# USER INPUTS
# =========================
lat_deg  = 43.6532       # Toronto
lon_deg  = -79.3832
alt_m    = 100.0         # WGS-84 height (m)

# Time window (UTC). These should overlap the SP3 file's coverage.
# The sample SP3 is for 2021-04-28; adjust as you like.
start_utc = datetime(2021, 4, 28, 18, 0, 0, tzinfo=timezone.utc)
end_utc   = datetime(2021, 4, 28, 23, 0, 0, tzinfo=timezone.utc)
dt_sec    = 60

elev_mask_deg = 10.0
use_constellation = "gps"  # gps | glonass | galileo | beidou | qzss | all

# Local SP3 path (sample CODE final MGEX)
sp3_path = "data/COD0MGXFIN_20211180000_01D_05M_ORB.SP3"

# =========================
# Load SP3 & prepare epoch list
# =========================
sp3 = glp.Sp3(sp3_path)   # uses gnss_lib_py's parser
unique_ms = np.unique(sp3["gps_millis"])

# If a time window is defined, clip the SP3 epochs to it
if start_utc and end_utc:
    start_ms = glp.datetime_to_gps_millis(start_utc)
    end_ms   = glp.datetime_to_gps_millis(end_utc)
    sel = (unique_ms >= start_ms) & (unique_ms <= end_ms)
    epoch_ms = unique_ms[sel]
    # Optional resampling: pick closest epochs to a regular grid
    if dt_sec:
        # build a target grid
        grid = np.arange(start_ms, end_ms + 1e-6, dt_sec * 1000.0)
        # for each grid time, pick the nearest SP3 epoch
        idx = np.searchsorted(unique_ms, grid, side="left")
        idx = np.clip(idx, 0, len(unique_ms) - 1)
        # choose nearer of idx and idx-1
        idx_minus = np.clip(idx - 1, 0, len(unique_ms) - 1)
        use_idx = np.where(
            np.abs(unique_ms[idx] - grid) < np.abs(unique_ms[idx_minus] - grid),
            idx, idx_minus
        )
        epoch_ms = np.unique(unique_ms[use_idx])
else:
    epoch_ms = unique_ms

if epoch_ms.size == 0:
    raise RuntimeError("No SP3 epochs in the requested time window. Adjust start_utc/end_utc.")

# =========================
# Receiver ECEF (meters)
# =========================
rx_ecef = glp.geodetic_to_ecef(np.array([[lat_deg, lon_deg, alt_m]]))[0].reshape(3, 1)

# =========================
# Constellation prefix map
# =========================
const_prefix = {
    "gps": "G", "glonass": "R", "galileo": "E",
    "beidou": "C", "qzss": "J"
}

# =========================
# Accumulators
# =========================
dop_times = []
visible_counts = []
pdop_vals, hdop_vals, vdop_vals = [], [], []

# =========================
# Loop over epochs
# =========================
for t_ms in epoch_ms:
    epoch_slice = sp3.where("gps_millis", float(t_ms), "eq")

    # Pull arrays (one per SV at this epoch)
    try:
        x = np.asarray(epoch_slice["x_sv_m"]).reshape(-1)
        y = np.asarray(epoch_slice["y_sv_m"]).reshape(-1)
        z = np.asarray(epoch_slice["z_sv_m"]).reshape(-1)
        sv_ids = epoch_slice["sv_id"]  # list-like
    except KeyError:
        continue

    if len(x) == 0:
        continue

    # Coerce sv_ids to strings
    sv_ids = [s.decode() if hasattr(s, "decode") else str(s) for s in sv_ids]

    # Constellation filter
    if use_constellation.lower() != "all":
        pref = const_prefix.get(use_constellation.lower(), None)
        if pref is None:
            pref = "G"
        keep_idx = [i for i, sv in enumerate(sv_ids) if sv.startswith(pref)]
    else:
        keep_idx = list(range(len(sv_ids)))

    if not keep_idx:
        visible_counts.append(0)
        dop_times.append(t_ms)
        pdop_vals.append(np.nan); hdop_vals.append(np.nan); vdop_vals.append(np.nan)
        continue

    # Build 3×N and convert km→m if needed (SP3 files are usually in km)
    pos_sv = np.vstack([x[keep_idx], y[keep_idx], z[keep_idx]])
    if np.nanmean(np.abs(pos_sv)) < 1e6:  # heuristic: km-scale values ~1e4
        pos_sv = pos_sv * 1000.0

    # Elevation/azimuth
    el_az = glp.ecef_to_el_az(rx_ecef, pos_sv)  # rows: [el(deg); az(deg)]
    el_deg = el_az[0, :]
    az_deg = el_az[1, :]

    # Elevation mask
    above = el_deg >= elev_mask_deg
    n_vis = int(np.sum(above))
    visible_counts.append(n_vis)
    dop_times.append(t_ms)

    if n_vis < 4:
        pdop_vals.append(np.nan); hdop_vals.append(np.nan); vdop_vals.append(np.nan)
        continue

    # Build NavData for DOP
    nav_epoch = glp.NavData()
    nav_epoch["gps_millis"] = np.full(n_vis, float(t_ms))
    nav_epoch["el_sv_deg"]  = el_deg[above]
    nav_epoch["az_sv_deg"]  = az_deg[above]

    dop_nav = glp.get_dop(nav_epoch, PDOP=True, HDOP=True, VDOP=True)
    pdop_vals.append(float(dop_nav["PDOP", 0]))
    hdop_vals.append(float(dop_nav["HDOP", 0]))
    vdop_vals.append(float(dop_nav["VDOP", 0]))

# =========================
# Post-process & availability
# =========================
dop_times = np.array(dop_times, dtype=float)
utc_times = glp.gps_millis_to_datetime(dop_times)

visible_counts = np.array(visible_counts, dtype=float)
pdop_vals = np.array(pdop_vals, dtype=float)

pdop_thresh = 2.0
avail_ge4 = 100.0 * np.mean(visible_counts >= 4) if visible_counts.size else 0.0
avail_pdop = 100.0 * np.mean((visible_counts >= 4) & (pdop_vals <= pdop_thresh)) if visible_counts.size else 0.0

print(f"Availability (≥4 SV): {avail_ge4:.1f}%")
print(f"Availability (PDOP ≤ {pdop_thresh}): {avail_pdop:.1f}%")

# =========================
# Plots
# =========================
fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

ax[0].plot(utc_times, visible_counts, marker="o")
ax[0].set_ylabel("# visible SVs")
ax[0].grid(True, linestyle="--", alpha=0.4)

ax[1].plot(utc_times, pdop_vals, marker="o")
ax[1].axhline(pdop_thresh, linestyle="--")
ax[1].set_ylabel("PDOP")
ax[1].set_xlabel("UTC time")
ax[1].grid(True, linestyle="--", alpha=0.4)

plt.tight_layout()
plt.show()

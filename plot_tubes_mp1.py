"""
Plot all flux-tube points for mp=1: geographic latitude vs altitude.

Each lp shell forms an arch from its southern footpoint (90 km) through
the apex (at the equator) to its northern footpoint (90 km).

Uses dipole-field approximation for geographic latitude:
    r = L * R_earth * cos^2(lambda)
=> cos^2(lambda) = (R_earth + h) / (L * R_earth)
=> lambda = +/- arccos( sqrt( (R_earth+h) / (L*R_earth) ) )
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

R_EARTH = 6371.2  # km

# ---------------------------------------------------------------------------
# 1.  Reconstruct L-values for all 67 lp shells
# ---------------------------------------------------------------------------
data = []
with open('tiegcm_defined_apex_heights') as f:
    for line in f:
        row = line.split()
        if len(row) >= 4 and row[0].lstrip('-').isdigit():
            data.append((int(row[0]), float(row[1]), float(row[2]), float(row[3])))

rows48 = data[:48]
Largest_L = 0.0
Smallest_L = 1000.0
i_largest = i_smallest = -1

for i, row in enumerate(rows48):
    L = row[3] / R_EARTH
    if 0 < L <= 4.0:
        if L > Largest_L:
            Largest_L = L
            i_largest = i
        if L < Smallest_L:
            Smallest_L = L
            i_smallest = i

N_tubes = i_smallest - i_largest + 1
L_selected = [rows48[i_largest + i][3] / R_EARTH for i in range(N_tubes)]
NLP = 2 * N_tubes - 1
L_all = [0.0] * NLP
for i in range(N_tubes):
    L_all[2 * i] = L_selected[i]
for i in range(N_tubes - 1):
    L_all[2 * i + 1] = (L_selected[i] + L_selected[i + 1]) / 2.0

# ---------------------------------------------------------------------------
# 2.  Height formula (mirrors generate_apex_coordinates.f)
# ---------------------------------------------------------------------------
def tube_heights(HA_km, lp_idx):
    """Return array of altitudes (km) from footpoint to just below apex."""
    factor = 0.2 if lp_idx < 11 else 0.4
    HA_m = HA_km * 1000.0
    height = 90000.0
    heights = []
    for icount in range(1, 1001):
        step = (float(icount - 1) * np.sqrt(max(HA_m - height, 0.0))) * factor
        height += step
        if height > HA_m:
            break
        heights.append(height / 1000.0)
    return np.array(heights)

# ---------------------------------------------------------------------------
# 3.  Build all tube arches: (lat, height) for both hemispheres + apex
# ---------------------------------------------------------------------------
cmap = plt.colormaps['plasma_r']
norm = mcolors.Normalize(vmin=1, vmax=NLP)

tubes = []
for lp_idx, L in enumerate(L_all, start=1):
    HA_km = max((L - 1.0) * R_EARTH, 0.0)
    color = cmap(norm(lp_idx))
    heights = tube_heights(HA_km, lp_idx)

    ratio = np.clip((R_EARTH + heights) / (L * R_EARTH), 0.0, 1.0)
    lat_deg = np.degrees(np.arccos(np.sqrt(ratio)))

    all_lats    = np.concatenate([-lat_deg, [0.0],  lat_deg[::-1]])
    all_heights = np.concatenate([heights,  [HA_km], heights[::-1]])
    tubes.append((lp_idx, HA_km, all_lats, all_heights, color))

# ---------------------------------------------------------------------------
# 4.  Three-panel plot
# ---------------------------------------------------------------------------
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(22, 8))
fig.subplots_adjust(left=0.05, right=0.88, wspace=0.28)
fig.suptitle('GIP flux-tube grid, mp=1  (dipole approximation)', fontsize=14, y=0.98)

panels = [
    (ax1, None,  20000, 'Full range (90 – apex)'),
    (ax2, None,   2000, 'Zoomed: 90 – 2000 km'),
    (ax3, None,    500, 'Zoomed: 90 – 500 km'),
]

for ax, _, ymax, title in panels:
    ax.set_xlabel('Magnetic latitude (°)  [dipole approx.]', fontsize=11)
    ax.set_ylabel('Altitude (km)', fontsize=11)
    ax.set_title(title, fontsize=11)
    ax.set_xlim([-65, 65])
    ax.axvline(0, color='gray', linewidth=0.6, linestyle='--', alpha=0.5)
    ax.axhline(90,  color='gray', linewidth=0.8, linestyle='-',  alpha=0.4, label='90 km footpoint')
    ax.axhline(200, color='red',  linewidth=0.9, linestyle=':',  alpha=0.7, label='200 km (E-F top)')
    ax.set_ylim([80, ymax])
    ax.grid(True, alpha=0.15)
    ax.legend(fontsize=8, loc='upper right')

for lp_idx, HA_km, lats, heights, color in tubes:
    lw = 0.6 if lp_idx <= 10 else 0.5
    ms = 1.8 if lp_idx <= 10 else 1.2

    for ax, _, ymax, _ in panels:
        mask = heights <= ymax
        if not mask.any():
            continue
        ax.plot(lats[mask], heights[mask], '-', color=color, linewidth=lw, alpha=0.8)
        ax.scatter(lats[mask], heights[mask], s=ms, color=color, alpha=0.7, zorder=3)

# Annotate selected lp shells on the full-range panel (right side, above apex)
for lp_target, offset_x in [(1, 2), (5, 2), (10, 2), (20, 2), (40, 3), (55, 4), (67, 4)]:
    L = L_all[lp_target - 1]
    HA_km = max((L - 1.0) * R_EARTH, 0.0)
    color = cmap(norm(lp_target))
    if HA_km < 20000:
        ax1.annotate(f' lp={lp_target}', xy=(0, HA_km),
                     fontsize=7, color=color,
                     va='center', ha='left')

# Colourbar on the right
sm = plt.cm.ScalarMappable(cmap='plasma_r', norm=norm)
sm.set_array([])
cbar_ax = fig.add_axes([0.90, 0.15, 0.018, 0.70])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('lp shell  (1 = outermost,  67 = innermost)', fontsize=10)
cbar.set_ticks([1, 10, 20, 30, 40, 50, 60, 67])
plt.savefig('tube_plot_mp1.png', dpi=180, bbox_inches='tight')
print("Saved: tube_plot_mp1.png")
plt.show()

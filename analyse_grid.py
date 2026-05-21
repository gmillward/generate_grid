"""
Analyse current GIP apex coordinate grid spacing.

Reads tiegcm_defined_apex_heights, reconstructs L-values for all 67 lp
shells (same logic as generate_apex_coordinates.f), then simulates the
current sqrt height formula to show point spacing for each shell.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

R_EARTH = 6371.2  # km

# ---------------------------------------------------------------------------
# 1.  Read tiegcm_defined_apex_heights (97 rows)
# ---------------------------------------------------------------------------
data = []
with open('tiegcm_defined_apex_heights') as f:
    for line in f:
        row = line.split()
        if len(row) >= 4 and row[0].lstrip('-').isdigit():
            idx = int(row[0])
            lat = float(row[1])
            apex_h = float(row[2])   # km
            apex_r = float(row[3])   # km
            data.append((idx, lat, apex_h, apex_r))

# ---------------------------------------------------------------------------
# 2.  Reproduce tube selection from generate_apex_coordinates.f
#     Read first 48 rows (southern hemisphere).  Select L <= 4 and L > 0.
#     Mirrors the Fortran exactly: i_largest = row index with max L (outermost),
#     i_smallest = row index with min L (innermost).
# ---------------------------------------------------------------------------
rows48 = data[:48]   # 1-indexed rows 1..48 in Fortran; 0-indexed 0..47 here

# In Fortran: Apex_L_value(i) = apex_radius / R_EARTH  (radius in km)
# i starts at 1; we use 0-based i here but replicate the loop logic.
Largest_L  = 0.0
Smallest_L = 1000.0
i_largest  = -1   # Fortran i_largest: index where L is maximum
i_smallest = -1   # Fortran i_smallest: index where L is minimum
i_first    = -1   # i_first_index_tube_less_than_L4

for i, row in enumerate(rows48):   # i is 0-based
    apex_r = row[3]                 # km
    L = apex_r / R_EARTH
    if 0 < L <= 4.0:
        if i_first == -1:
            i_first = i
        if L < Smallest_L:
            Smallest_L = L
            i_smallest = i
        if L > Largest_L:
            Largest_L = L
            i_largest  = i

# tiegcm_L_value(i) = Apex_L_value(i + i_largest)  for i = 1..iNumber_of_tubes
N_tubes = i_smallest - i_largest + 1   # should be 34
L_selected = [rows48[i_largest + i][3] / R_EARTH for i in range(N_tubes)]

# Double: tiegcm_L_value2 has 2*N_tubes-1 = 67 entries
# Odd 1-based positions (0-based 0,2,4,...) = L_selected[0], [1], [2], ...
# Even 1-based positions (0-based 1,3,5,...) = midpoints
NLP = 2 * N_tubes - 1
L_doubled2 = [0.0] * NLP
for i in range(N_tubes):
    L_doubled2[2 * i] = L_selected[i]
for i in range(N_tubes - 1):
    L_doubled2[2 * i + 1] = (L_selected[i] + L_selected[i + 1]) / 2.0

print(f"N_tubes = {N_tubes},  NLP = {NLP}")
print(f"L range: {L_doubled2[0]:.4f} (lp=1, outermost) to {L_doubled2[-1]:.4f} (lp={NLP}, innermost)")
print(f"Apex height range: {(L_doubled2[0]-1)*R_EARTH:.0f} km to {(L_doubled2[-1]-1)*R_EARTH:.0f} km")

# ---------------------------------------------------------------------------
# 3.  Simulate the height formula for each lp shell
# ---------------------------------------------------------------------------
def compute_heights(HA_km, lp_index):
    """
    Returns list of altitudes (km) from footpoint (90 km) to apex.
    lp_index is 1-based (lp=1 outermost, lp=67 innermost).
    """
    factor = 0.2 if lp_index < 11 else 0.4
    HA_m = HA_km * 1000.0
    height = 90000.0  # metres
    heights = [90.0]  # km
    for iht in range(1, 1001):
        step = (float(iht - 1) * np.sqrt(HA_m - height)) * factor
        height += step
        if height > HA_m:
            break
        heights.append(height / 1000.0)
    return heights

# ---------------------------------------------------------------------------
# 4.  Summary table
# ---------------------------------------------------------------------------
print(f"\n{'lp':>4} {'apex_km':>9} {'L':>6} {'pts/hemi':>9} {'avg_dh_km':>10} {'min_dh_km':>11}")
print("-" * 55)

results = []
for lp, L in enumerate(L_doubled2, start=1):
    HA_km = (L - 1.0) * R_EARTH
    if HA_km < 0.1:
        HA_km = 0.0
    heights = compute_heights(HA_km, lp)
    diffs = np.diff(heights) if len(heights) > 1 else np.array([0.0])
    n_pts = len(heights)
    avg_dh = diffs.mean() if len(diffs) else 0.0
    min_dh = diffs.min() if len(diffs) else 0.0
    results.append((lp, HA_km, L, heights, n_pts, avg_dh, min_dh))
    print(f"{lp:>4} {HA_km:>9.1f} {L:>6.4f} {n_pts:>9} {avg_dh:>10.2f} {min_dh:>11.2f}")

total_pts = sum(r[4] * 2 - 1 for r in results)  # both hemis + apex
print(f"\nTotal packed points (npts2 estimate): {total_pts}")

# ---------------------------------------------------------------------------
# 5.  Plots
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(18, 12))
fig.suptitle("Current GIP Grid: Along-tube Point Spacing (sqrt formula)", fontsize=14)

gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.35)

# 5a. Points per hemisphere vs lp
ax1 = fig.add_subplot(gs[0, 0])
lps = [r[0] for r in results]
pts = [r[4] for r in results]
ax1.bar(lps, pts, color='steelblue', width=0.8)
ax1.set_xlabel("lp shell")
ax1.set_ylabel("Points per hemisphere")
ax1.set_title("Points per hemisphere vs lp")
ax1.axvline(x=11, color='red', linestyle='--', alpha=0.6, label='factor change (lp=11)')
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3)

# 5b. Average step size vs apex height
ax2 = fig.add_subplot(gs[0, 1])
apexes = [r[1] for r in results]
avg_dhs = [r[5] for r in results]
ax2.scatter(apexes, avg_dhs, s=20, color='darkorange')
ax2.set_xlabel("Apex height (km)")
ax2.set_ylabel("Average dh (km)")
ax2.set_title("Average step size vs apex height")
ax2.set_xscale('log')
ax2.axvline(x=200, color='red', linestyle='--', alpha=0.5, label='200 km (E-F region top)')
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)

# 5c. Height profiles for selected lp shells — height spacing
ax3 = fig.add_subplot(gs[1, :])
highlight_lps = [1, 5, 10, 20, 30, 40, 47, 49, 55, 60, 65, 67]
cmap = plt.cm.viridis
for lp_target in highlight_lps:
    r = results[lp_target - 1]
    lp, HA_km, L, heights, n_pts, avg_dh, min_dh = r
    if len(heights) < 2:
        continue
    diffs = np.diff(heights)
    h_mid = [(heights[i] + heights[i+1]) / 2.0 for i in range(len(diffs))]
    color = cmap(lp / 67.0)
    ax3.semilogy(h_mid, diffs, '-o', ms=3, color=color,
                 label=f"lp={lp} (apex={HA_km:.0f}km)")

ax3.set_xlabel("Height (km)")
ax3.set_ylabel("Step size dh (km) [log scale]")
ax3.set_title("Along-tube step size vs height for selected lp shells")
ax3.axvline(x=90, color='gray', linestyle=':', alpha=0.5)
ax3.axvline(x=200, color='red', linestyle='--', alpha=0.5, label='200 km (E-F region)')
ax3.set_xlim([85, 3000])
ax3.legend(fontsize=7, ncol=2, loc='upper left')
ax3.grid(True, alpha=0.3)

plt.savefig('grid_analysis_current.png', dpi=150, bbox_inches='tight')
print("\nSaved: grid_analysis_current.png")
plt.show()

import numpy as np
import matplotlib.pyplot as plt

# --- Assume you have your data like this ---
# Example Data (replace with your actual data)
l_values = np.geomspace(10, 2000, 50) # Example x-values (multipoles)
# Make sure these are numpy arrays for easier calculation
bispectrum_ex = np.array((l_values**-2.5) * (1 + np.random.normal(0, 0.05, size=l_values.shape))) # Your calculated values
bispectrum_tosh = np.array(l_values**-2.5) # Toshiya's reference values
# --- End Example Data ---

# Labels and Title
x_label = "Multipole moment $l$"
y_label_main = "$l_1 l_2 l_3 B_{l_1 l_2 l_3} / (2\pi)^2$" # Adjust as needed
y_label_ratio = "Ratio (My Calc / Toshiya)"
plot_title = "Comparison of Weak Lensing Bispectra"

# --- Create the Figure and Axes objects ---
# Create a figure with 2 subplots stacked vertically (2 rows, 1 column)
# sharex=True links the x-axes of the two plots
# gridspec_kw allows setting relative heights (e.g., main plot 3 times taller than ratio)
fig, axs = plt.subplots(
    nrows=2,
    ncols=1,
    sharex=True,
    figsize=(8, 6), # Adjust figure size as needed
    gridspec_kw={'height_ratios': [3, 1]} # Main plot is 3x taller than ratio plot
)

# --- Plot on the top Axes (axs[0]) ---
ax_main = axs[0]
ax_main.plot(l_values, bispectrum_ex, label='My Calculation (ex)', marker='o', linestyle='-', markersize=4)
ax_main.plot(l_values, bispectrum_tosh, label='Toshiya Ref (tosh)', marker='', linestyle='--', color='k')

# Customize the top plot (use 'ax_main.' instead of 'plt.')
ax_main.set_ylabel(y_label_main)
ax_main.set_yscale('log')  # Often bispectra are plotted on log scales
ax_main.set_xscale('log')
ax_main.set_title(plot_title)
ax_main.legend()
ax_main.grid(True, which='both', linestyle=':', linewidth=0.5)

# --- Calculate and Plot the Ratio on the bottom Axes (axs[1]) ---
ax_ratio = axs[1]

# Calculate ratio, handle potential division by zero if bispectrum_tosh can be 0
ratio = np.divide(bispectrum_ex, bispectrum_tosh,
                  out=np.full_like(bispectrum_ex, np.nan), # Fill with NaN where division is invalid
                  where=bispectrum_tosh!=0)

ax_ratio.plot(l_values, ratio, marker='.', linestyle='-', color='r', markersize=5)

# Add a horizontal line at y=1 for reference
ax_ratio.axhline(1, linestyle='--', color='grey', linewidth=1)

# Customize the bottom plot
ax_ratio.set_xlabel(x_label)
ax_ratio.set_ylabel(y_label_ratio)
# ax_ratio.set_xscale('log') # Already set by sharex=True
# Set y-limits for the ratio plot to focus the view (optional)
# For example, if you expect ratios close to 1:
ratio_ylim_padding = 0.1 # Adjust padding as needed
min_ratio = np.nanmin(ratio) if not np.all(np.isnan(ratio)) else 0.9
max_ratio = np.nanmax(ratio) if not np.all(np.isnan(ratio)) else 1.1
ax_ratio.set_ylim(max(0, min_ratio - ratio_ylim_padding), max_ratio + ratio_ylim_padding) # Ensure lower limit isn't negative

ax_ratio.grid(True, linestyle=':', linewidth=0.5)


# --- Final Adjustments ---
# Remove space between subplots
# fig.subplots_adjust(hspace=0) # Alternative: plt.tight_layout() often works well

# Improve layout to prevent labels overlapping
plt.tight_layout()

# Show the plot
plt.show()
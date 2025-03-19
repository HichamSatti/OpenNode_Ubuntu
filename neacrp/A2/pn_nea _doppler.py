import numpy as np
import matplotlib.pyplot as plt

# File path to your data file
file_path = 'dop1x1x1.txt'
file_path2 ='dop1x1x2.txt'
file_path3 = 'dop2x2.txt'
file_path4 = 'dop4x4.txt'

# Read the data file using numpy (skipping the first two header lines)
data = np.genfromtxt(file_path, skip_header=2)
data2 = np.genfromtxt(file_path2, skip_header=2)
data3 = np.genfromtxt(file_path3, skip_header=2)
data4 = np.genfromtxt(file_path4, skip_header=2)

# Extract columns for time and relative power
time = data[:, 0]              # First column: time
rel_power = data[:, 1]         # Second column: relative power
rel_power2 = data2[:, 1]       # Second column: relative power
rel_power3 = data3[:, 1]       # Second column: relative power
rel_power4 = data4[:, 1]       # Second column: relative power

# Second set of data for comparison
time2 = [0.5, 1, 1.5, 2, 2.5]
rel_power_ref = [0.9997826086956519, 0.2541784302653869, 0.2012718802936193, 
0.19472896668548834, 0.18545454545454554]

# Generate the full time list from 0 to 60.0 in increments of 0.25
result = []  # List for new time points
new_rel_power = []  # List for new relative power values
rel_power_index = 0  # Index to track the corresponding rel_power2 values


time1 = [
    0.015600624024960999,
    0.39781591263650545,
    0.7722308892355695,
    1.154446177847114,
    1.606864274570983,
    1.903276131045242,
    2.5663026521060845,
    3.0967238689547583,
    3.5959438377535102,
    3.9859594383775354,
    4.492979719188767,
    4.742589703588144,
    4.976599063962559
]
time2 = [
    0.078003120124805,
    0.33541341653666146,
    0.6786271450858035,
    0.9750390015600624,
    1.4196567862714509,
    2.145085803432137,
    2.714508580343214,
    3.2761310452418098,
    3.8377535101404057,
    4.251170046801872,
    4.609984399375975,
    4.84399375975039,
    4.976599063962559
]

rel_power_ref1 = [
    546.4949704142012,
    548.7040433925049,
    550.3608481262328,
    551.4653846153847,
    552.8460552268245,
    553.1221893491124,
    554.2267258382643,
    554.7789940828403,
    555.0551282051282,
    555.3312623274162,
    555.8835305719922,
    555.8835305719922,
    555.8835305719922
]
rel_power_ref2 = [
    545.3904339250494,
    547.0472386587771,
    548.151775147929,
    549.2563116370809,
    550.6369822485208,
    552.0176528599605,
    552.8460552268245,
    553.3983234714004,
    553.6744575936884,
    553.9505917159763,
    553.9505917159763,
    553.9505917159763,
    553.9505917159763
]


#for t in np.arange(0, 2.5, 0.5):
#    if t in time2:
#        result.append(t)
#        new_rel_power.append(rel_power_ref[rel_power_index])
#        rel_power_index += 1
#    else:
#        result.append(None)
#        new_rel_power.append(None)

# Plot for Relative Power
plt.figure(figsize=(12, 7))

# Plot Time vs. rel_power1
plt.plot(time[:len(rel_power3)], rel_power3, color='red', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Radial Node 2x2')
plt.plot(time[:len(rel_power4)], rel_power4, color='green', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Radial Node 4x4')
plt.plot(time[:len(rel_power)], rel_power, color='blue', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Axial Node 1x1')
plt.plot(time[:len(rel_power2)], rel_power2, color='magenta', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Axial Node 1x2')
# Plot Time vs. rel_power2 (second dataset)
# plt.plot(result[:len(new_rel_power)], new_rel_power, color='cyan', marker='*', linestyle='None', linewidth=1.0, markersize=8, label='SPANDEX')

plt.plot(time1, rel_power_ref1, color='cyan', marker='*', linestyle='--', linewidth=0.6, markersize=8, label='PARCS')
plt.plot(time2, rel_power_ref2, color='crimson', marker='*', linestyle='--', linewidth=0.6, markersize=8, label='CORCA-K')

# Enhancing the plot
plt.xlabel('Time (s)', fontsize=14, fontweight='bold')
plt.ylabel('Avg. Doppler Temperature', fontsize=14, fontweight='bold')

# Adding grid, legend, and fine-tuning axis
plt.grid(True, linestyle='--', linewidth=0.6, alpha=0.7)
plt.xticks(np.arange(0, 5.5, 0.5), fontsize=12)  # Time axis ticks every 5 seconds
plt.yticks(np.arange(485, 635, 10), fontsize=12)  # Power axis ticks every 0.2
plt.ylim(485, 625)  # Limit y-axis range
plt.legend(fontsize=12, loc='lower right')

# Set x-axis limits to start from 0 without extra space
plt.xlim(0, 5)

# Saving the enhanced plot
plt.savefig('enhanced_relative_power_plot.png', dpi=300)
plt.show()

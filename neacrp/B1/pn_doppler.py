import numpy as np
import matplotlib.pyplot as plt

# File path to your data file
file_path = 'dop4x4.txt'
file_path2 = 'dop2x2.txt'
file_path3 = 'dop1x1.txt'
#file_path4 = 'nea1x1x2.txt'
#file_path5 = 'nea1x1x4.txt'
#file_path3 = 'nea1x2.txt'

# Read the data file using numpy (skipping the first two header lines)
data = np.genfromtxt(file_path, skip_header=2)
data2 = np.genfromtxt(file_path2, skip_header=2)
data3 = np.genfromtxt(file_path3, skip_header=2)
#data4 = np.genfromtxt(file_path4, skip_header=2)
#data5 = np.genfromtxt(file_path5, skip_header=2)

# Extract columns for time and relative power
time = data[:, 0]              # First column: time
rel_power = data[:, 1]         # Second column: relative power
rel_power2 = data2[:, 1]         # Second column: relative power
rel_power3 = data3[:, 1]         # Second column: relative power
#rel_power4 = data4[:, 1]         # Second column: relative power
#rel_power5 = data5[:, 1]         # Second column: relative power

# Second set of data for comparison
time2 = [0.5, 1, 1.5, 2, 2.5]
#rel_power_ref = [0.9997826086956519, 0.2541784302653869, 0.2012718802936193, 
#0.19472896668548834,0.19420783645655876]

time_ref = [
    0.05124450951683748,
    0.424597364568082,
    0.5929721815519766,
    0.8931185944363105,
    1.6691068814055638,
    2.225475841874085,
    3.2137628111273795,
    3.909224011713031,
    4.333821376281113,
    4.904831625183016
]
rel_power_ref = [
    285.7388888888889,
    286.47962962962964,
    300.18333333333334,
    307.5907407407407,
    314.9981481481482,
    317.9611111111111,
    321.2944444444445,
    323.51666666666665,
    324.2574074074074,
    324.9981481481482
]

time_ref1 = [
    0.10980966325036604,
    0.424525452544082,
    0.5051244509516838,
    0.5710102489019033,
    1.0248901903367496,
    1.7642752562225477,
    2.320644216691069,
    3.28696925329429,
    3.843338213762811,
    4.428989751098097,
    4.79502196193265
]

rel_power_ref1 = [
    285.7388888888889,
    286.4796225466964,
    292.40555555555557,
    299.0722222222222,
    309.0722222222222,
    314.9981481481482,
    317.5907407407407,
    320.92407407407404,
    322.40555555555557,
    323.51666666666665,
    324.2574074074074
]

time_ref2 = [
    0.7061503416856493,
    0.8086560364464693,
    0.9453302961275627,
    1.2300683371298406,
    1.7084282460136675,
    2.5056947608200457,
    3.0068337129840548,
    3.507972665148064,
    3.997722095671982,
    4.510250569476082,
    4.9886104783599095
]

rel_power_ref2 = [
    288.97815912636503,
    294.12636505460216,
    299.3915756630265,
    304.53978159126365,
    309.33697347893917,
    314.3681747269891,
    316.591263650546,
    318.22932917316695,
    319.63338533541344,
    320.9204368174727,
    321.73946957878314
]

# Generate the full time list from 0 to 60.0 in increments of 0.25
result = []  # List for new time points
new_rel_power = []  # List for new relative power values
rel_power_index = 0  # Index to track the corresponding rel_power2 values

#for t in np.arange(0, 3.0, 0.5):
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
plt.plot(time[:len(rel_power3)], rel_power3, color='green', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Axial Node 1x1')
plt.plot(time[:len(rel_power2)], rel_power2, color='red', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Radial Node 2x2')
plt.plot(time[:len(rel_power)], rel_power, color='blue', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Radial Node 4x4')
#plt.plot(time[:len(rel_power4)], rel_power4, color='yellow', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Axial Node 1x1x2')
#plt.plot(time[:len(rel_power5)], rel_power5, color='magenta', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Axial Node 1x1x4')
plt.xlim(0, max(time))  # Start x-axis at 0

# Plot Time vs. rel_power2 (second dataset)
#plt.plot(result[:len(new_rel_power)], new_rel_power, color='cyan', marker='*', linestyle='None', linewidth=1.0, markersize=8, label='PANTHER')

# Plot reference results
plt.plot(time_ref, rel_power_ref, color='cyan', marker='*', linestyle='--', linewidth=1.0, markersize=8, label='PARCS')
plt.plot(time_ref1, rel_power_ref1, color='crimson', marker='*', linestyle='--', linewidth=1.0, markersize=8, label='CORCA-K')
plt.plot(time_ref2, rel_power_ref2, color='teal', marker='*', linestyle='--', linewidth=1.0, markersize=8, label='PANTHER')

# Enhancing the plot
plt.title('', fontsize=18, fontweight='bold', color='navy')
plt.xlabel('Time (s)', fontsize=14, fontweight='bold')
plt.ylabel('Fuel Doppler Temperature (°C)', fontsize=14, fontweight='bold')

# Adding grid, legend, and fine-tuning axis
plt.grid(True, linestyle='--', linewidth=0.6, alpha=0.7)
plt.xticks(np.arange(0, 5.5, 0.5), fontsize=12)  # Time axis ticks every 5 seconds
plt.yticks(np.arange(285, 340, 5), fontsize=12)  # Power axis ticks every 0.2
plt.ylim(285, 330)  # Limit y-axis range
plt.legend(fontsize=12, loc='lower right')
plt.tight_layout()

# Saving the enhanced plot
plt.savefig('enhanced_relative_power_plot.png', dpi=500)
plt.show()



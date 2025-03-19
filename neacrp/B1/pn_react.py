import numpy as np
import matplotlib.pyplot as plt

# File path to your data file
file_path_ = 'react4x4.txt'
file_path_2 = 'react2x2.txt'
file_path_3 = 'react1x4.txt'
#file_path4 = 'nea1x1x2.txt'
#file_path5 = 'nea1x1x4.txt'
#file_path3 = 'nea1x2.txt'

# Read the data file using numpy (skipping the first two header lines)
data_ = np.genfromtxt(file_path_, skip_header=2)
data_2 = np.genfromtxt(file_path_2, skip_header=2)
data_3 = np.genfromtxt(file_path_3, skip_header=2)
#data4 = np.genfromtxt(file_path4, skip_header=2)
#data5 = np.genfromtxt(file_path5, skip_header=2)

# Extract columns for time and relative power
time_ = data_[:, 0]              # First column: time
rel_power_ = data_[:, 1]         # Second column: relative power
rel_power_2 = data_2[:, 1]         # Second column: relative power
rel_power_3 = data_3[:, 1]         # Second column: relative power
#rel_power4 = data4[:, 1]         # Second column: relative power
#rel_power5 = data5[:, 1]         # Second column: relative power

# Second set of data for comparison
time2 = [0.5, 1, 1.5, 2, 2.5]
#rel_power_ref = [0.9997826086956519, 0.2541784302653869, 0.2012718802936193, 
#0.19472896668548834,0.19420783645655876]

time_ref_ = [
    0.0182648401826484, 0.0365296803652968, 0.0578386605783866, 
    0.076103500761035, 0.1035007610350076, 0.2343987823439878, 
    0.3835616438356164, 0.4687975646879756, 0.5814307458143074, 
    0.6605783866057838, 0.7640791476407914, 0.8706240487062404, 
    1.0076103500761033, 1.1933028919330289, 1.4033485540334856, 
    1.5981735159817352, 1.7990867579908674, 1.9969558599695585
]
rec_power_ref = [
    0.060115606936416176, 0.3907514450867052, 0.776878612716763, 
    1.005780346820809, 1.07514450867052, 1.0797687861271674, 
    1.0774566473988438, 1.052023121387283, 0.9132947976878611, 
    0.8254335260115606, 0.7445086705202312, 0.6728323699421964, 
    0.6173410404624277, 0.5757225433526011, 0.5433526011560693, 
    0.5248554913294797, 0.5063583815028901, 0.4878612716763005
]

time_ref_1 = [
    0.02735562310030395,
    0.03343465045592705,
    0.0425531914893617,
    0.05167173252279635,
    0.057750759878419454,
    0.0790273556231003,
    0.14285714285714285,
    0.4376899696048632,
    0.5379939209726444,
    0.6109422492401215,
    0.7477203647416414,
    0.8693009118541033,
    0.9939209726443768,
    1.2948328267477203,
    1.5075987841945289,
    1.7325227963525835,
    1.966565349544073
]
rec_power_ref1 = [
    0.18728323699421964,
    0.3005780346820809,
    0.5294797687861271,
    0.7052023121387282,
    0.8208092485549132,
    0.9919075144508669,
    1.0774566473988438,
    1.0635838150289016,
    0.9710982658959536,
    0.8716763005780346,
    0.7421965317919075,
    0.6658959537572253,
    0.6127167630057803,
    0.5549132947976878,
    0.5317919075144508,
    0.5086705202312138,
    0.4878612716763005
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
#plt.plot(time[:len(rel_power3)], rel_power3, color='green', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Axial Node 1x4')
plt.plot(time_[:len(rel_power_2)], rel_power_2, color='red', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Radial Node 2x2')
plt.plot(time_[:len(rel_power_)], rel_power_, color='blue', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Radial Node 4x4')
#plt.plot(time[:len(rel_power4)], rel_power4, color='yellow', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Axial Node 1x1x2')
#plt.plot(time[:len(rel_power5)], rel_power5, color='magenta', marker='o', linestyle='-', linewidth=1.0, markersize=1, label='Axial Node 1x1x4')
plt.xlim(0, max(time_))  # Start x-axis at 0

# Plot Time vs. rel_power2 (second dataset)
#plt.plot(result[:len(new_rel_power)], new_rel_power, color='cyan', marker='*', linestyle='None', linewidth=1.0, markersize=8, label='PANTHER')

# Plot reference results
plt.plot(time_ref_, rec_power_ref, color='cyan', marker='*', linestyle='--', linewidth=1.0, markersize=8, label='PARCS')
plt.plot(time_ref_1, rec_power_ref1, color='crimson', marker='*', linestyle='--', linewidth=1.0, markersize=8, label='CORCA-K')

# Enhancing the plot
plt.title('', fontsize=18, fontweight='bold', color='navy')
plt.xlabel('Time (s)', fontsize=14, fontweight='bold')
plt.ylabel('Reactivity ($)', fontsize=14, fontweight='bold')

# Adding grid, legend, and fine-tuning axis
plt.grid(True, linestyle='--', linewidth=0.6, alpha=0.7)
plt.xticks(np.arange(0, 5.5, 0.5), fontsize=12)  # Time axis ticks every 5 seconds
plt.yticks(np.arange(0, 1.4, 0.2), fontsize=12)  # Power axis ticks every 0.2
plt.ylim(0, 1.2)  # Limit y-axis range
plt.legend(fontsize=12, loc='upper right')
plt.tight_layout()

# Saving the enhanced plot
plt.savefig('enhanced_relative_power_plot.png', dpi=500)
plt.show()



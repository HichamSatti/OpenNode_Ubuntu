import numpy as np
import matplotlib.pyplot as plt

# Define file paths
#power_files1 = ['maxf1x1.txt']
power_files2 = ['maxf4x4.txt', 'maxf2x2.txt']
#reactivity_files1 = ['dop1x1.txt']
reactivity_files2 = ['dop2x2.txt', 'dop4x4.txt']

#power_names1 = ['Radial Node 1x1']
#power_names_1 = ['Axial Node 1x2', 'Axial Node 1x4']
power_names2 = ['Radial Node 2x2', 'Radial Node 4x4']
#react_names1 = ['Radial Node 1x4']
react_names2 = ['Radial Node 4x4', 'Radial Node 2x2']

# Initialize dictionaries to store data
reactivity_data1 = {}
reactivity_data2 = {}
power_data1 = {}
power_data_1 = {}
power_data2 = {}


colors1 = ['black']
colors_1 = ['dimgray','black']
colors__1 = ['dimgray','black']
colors2 = ['blue']
colors_2 = ['blue','dodgerblue']




time_ref_ = [
    0.005,
    0.19101123595505565,
    0.3707865168539321,
    0.5955056179775279,
    0.7640449438202245,
    0.9101123595505615,
    0.9987640449438195
    ]
rel_power_ref_ = [
    542.9725609756098,
    543.6585365853658,
    544.420731707317,
    545.2591463414634,
    545.7164634146342,
    546.0975609756098,
    546.3262195121952
    ]

time_ref_1=[
    0.022471910112359637,
    0.16853932584269665,
    0.292134831460674,
    0.43820224719101103,
    0.584269662921348,
    0.7191011235955052,
    0.8764044943820214,
    0.9887640449438195
    ]
rel_power_ref_1=[
    543.9634146341464,
    544.2682926829268,
    544.4969512195122,
    544.8018292682926,
    545.030487804878,
    545.1067073170732,
    545.3353658536586,
    545.4878048780488
    ]


time3 = [
    0.001,
    0.1125,
    0.2,
    0.28750000000000003,
    0.3875,
    0.5,
    0.6125,
    0.7375,
    0.8500000000000001,
    0.925,
    1,
    1
    ]
rel_power_ref3 = [
    1574,
    1574.400684931507,
    1574.8458904109589,
    1575.335616438356,
    1575.914383561644,
    1576.4931506849316,
    1577.1164383561643,
    1577.695205479452,
    1578.2294520547946,
    1578.541095890411,
    1578.9417808219177,
    1578.9417808219177
    ]
time1 = [
    0,
    0.1125,
    0.25,
    0.41250000000000003,
    0.5875,
    0.7125,
    0.875,
    1
    ]
rel_power_ref1 =[
    1578.6301369863013,
    1578.8527397260275,
    1579.1198630136987,
    1579.4315068493152,
    1579.6986301369864,
    1579.8321917808219,
    1580.054794520548,
    1580.2328767123288
    ]




time_ref = [
    0.04038772213247173,
    0.33117932148626816,
    0.4361873990306947,
    0.4765751211631664,
    0.5169628432956381,
    0.5492730210016155,
    0.6946688206785138,
    0.8885298869143781,
    1.308562197092084,
    1.8416801292407108,
    2.4232633279483036,
    2.8109854604200324,
    3.263327948303716,
    3.618739903069467,
    3.9176090468497575,
    4.361873990306947,
    4.975767366720517
]
rel_power_ref = [
    285.75539568345323,
    286.1870503597122,
    288.3453237410072,
    294.24460431654677,
    303.45323741007195,
    312.23021582733816,
    319.28057553956836,
    324.0287769784173,
    330.7913669064748,
    336.1151079136691,
    340.431654676259,
    342.58992805755395,
    344.46043165467626,
    345.8992805755396,
    346.76258992805754,
    348.3453237410072,
    350.0719424460432
]






time_ref_2 = [
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
rel_power_ref_2 = [
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


# Read power files
#for idx, file_path in enumerate(power_files1):
#    try:
#        name = power_names_1[idx]
#        data = np.genfromtxt(file_path, skip_header=2)
#        time = data[:, 0]  # Time column
#        relative_power_1 = data[:, 1]  # Relative power column
#        power_data_1[name] = (time, relative_power_1)
#    except Exception as e:
#        print(f"Error reading file {file_path}: {e}")
for idx, file_path in enumerate(power_files2):
    try:
        name = power_names2[idx]
        data = np.genfromtxt(file_path, skip_header=2)
        time = data[:, 0]  # Time column
        relative_power2 = data[:, 1]  # Relative power column
        power_data2[name] = (time, relative_power2)
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")

# Read reactivity files
for idx_, file_path_ in enumerate(reactivity_files2):
    try:
        name = react_names2[idx_]
        data = np.genfromtxt(file_path_, skip_header=2)
        time = data[:, 0]  # Time column
        reactivity2 = data[:, 1]  # Reactivity column
        reactivity_data2[name] = (time, reactivity2)
    except Exception as e:
        print(f"Error reading file {file_path_}: {e}")

# Plotting with dual y-axes
plt.figure(figsize=(12, 7))


# Plot relative power on the left y-axis
ax1 = plt.gca()  # Get the current axes
color_index_ = 0  # Initialize color index
for key, (time, rel_power1) in power_data1.items():
    ax1.plot(time, rel_power1, label=key, linestyle='--', marker='o',linewidth=1.2, markersize=1,
    color=colors1[color_index_ % len(colors1)])
    color_index_ += 1  # Update color index
#color_index_3 = 0  # Initialize color index
#for key, (time, rel_power_1) in power_data_1.items():
#    ax1.plot(time, rel_power_1, label=key, linestyle='-.', marker='o',linewidth=1.2, markersize=1,
#    color=colors__1[color_index_3 % len(colors__1)])
#    color_index_3 += 1  # Update color index
color_index_2 = 0  # Initialize color index
for key, (time, rel_power2) in power_data2.items():
    ax1.plot(time, rel_power2, label=key, linestyle='-', marker='o',linewidth=1.2, markersize=1,
    color=colors_1[color_index_2 % len(colors_1)])
    color_index_2 += 1  # Update color index

ax1.set_xlabel('Time (s)', fontsize=14, fontweight='bold')
ax1.set_ylabel('Maximum Fuel Temperature (°C)', fontsize=14, fontweight='bold', color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.grid(True, linestyle='--', alpha=0.7)
# Set the limits of the x and y axes to start from 0
ax1.set_xlim(left=0)  # x-axis starts from 0
#ax1.set_ylim(bottom=200)  # Left y-axis starts from 0
#ax1.set_ylim(top=800)  # Left y-axis starts from 0
ax1.set_ylim(1300, 1700)
# Set y-axis ticks (200 to 500 by step of 50)
yticks = range(1300, 1701, 50)  # 501 is used to include 500
ax1.set_yticks(yticks)


ax1.set_xlim(right=5)  # x-axis starts from 0

# Create a second y-axis for reactivity
ax2 = ax1.twinx()
color_index = 0  # Initialize color index
for key, (time, reactivity1) in reactivity_data1.items():
    ax2.plot(time, reactivity1, label=key, linestyle='--', marker='x',linewidth=1.2, markersize=1, 
    color=colors2[color_index % len(colors2)])
    color_index += 1  # Update color index

color_index2 = 0  # Initialize color index
for key, (time, reactivity2) in reactivity_data2.items():
    ax2.plot(time, reactivity2, label=key, linestyle='-', marker='x',linewidth=1.2, markersize=1, 
    color=colors_2[color_index2 % len(colors_2)])
    color_index2 += 1  # Update color index


ax2.set_ylabel('Fuel Doppler Temperature (°C)', fontsize=14, fontweight='bold', color='blue')
ax2.tick_params(axis='y', labelcolor='blue')
# Set the right y-axis to start from 0
ax2.set_xlim(left=0)  # x-axis starts from 0

ax2.set_ylim(500, 800)
# Set y-axis ticks (200 to 500 by step of 50)
yticks = range(500, 801, 50)  # 501 is used to include 500
ax2.set_yticks(yticks)

ax2.set_xlim(right=1)  # x-axis starts from 0

# Set fontweight='bold' for both x and y ticks on ax1
for label in ax1.get_xticklabels() + ax1.get_yticklabels():
    label.set_fontweight('bold')

# Set fontweight='bold' for both x and y ticks on ax2
for label in ax2.get_yticklabels():
    label.set_fontweight('bold')

# Plot reference results of relative power
ax1.plot(time3, rel_power_ref3, color='black', marker='+', linestyle='--', linewidth=0.6, markersize=8, label='PANTHER')
ax1.plot(time1, rel_power_ref1, color='black', marker='*', linestyle='--', linewidth=0.6, markersize=8, label='OKAPI')
# Plot reference results of reactivity
ax2.plot(time_ref_, rel_power_ref_, color='blue', marker='+', linestyle='--', linewidth=0.6, markersize=8, label='PANTHER')
ax2.plot(time_ref_1, rel_power_ref_1, color='blue', marker='*', linestyle='--', linewidth=0.6, markersize=8, label='OKAPI')
#ax2.plot(-20,-20,color='chartreuse',  marker='^', linestyle='', linewidth=0.6, label='PANTHER', markersize=6)


# Add a combined legend
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc='center right', fontsize=9)


# Title and save the plot
#plt.title('Relative Power and Reactivity Over Time', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('dual_y_axes_plot.png', dpi=500)
plt.show()

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




time1 = [
    0.7311320754716981,
    0.8490566037735849,
    0.9669811320754718,
    1.0731132075471699,
    1.214622641509434,
    1.3443396226415096,
    1.5212264150943398,
    1.7452830188679247,
    1.8867924528301887,
    2.099056603773585,
    2.535377358490566,
    3.018867924528302,
    3.785377358490566,
    4.410377358490567,
    4.988207547169812
]
rel_power_ref1 = [
    400.64516129032256,
    412.9032258064516,
    427.4193548387097,
    439.6774193548387,
    452.5806451612903,
    464.51612903225805,
    475.48387096774195,
    486.1290322580645,
    494.83870967741933,
    505.80645161290323,
    518.7096774193549,
    531.9354838709678,
    547.4193548387096,
    557.4193548387096,
    564.8387096774194
]
time3 = [
    0.7783018867924529,
    0.8372641509433962,
    0.9433962264150944,
    1.1202830188679247,
    1.25,
    1.3915094339622642,
    1.509433962264151,
    1.7099056603773586,
    1.8867924528301887,
    2.0872641509433962,
    2.4646226415094343,
    2.9009433962264155,
    3.6320754716981134,
    4.268867924528302,
    4.988207547169812
]
rel_power_ref3 = [
    401.61290322580646,
    411.61290322580646,
    423.8709677419355,
    439.6774193548387,
    450.9677419354839,
    462.5806451612903,
    472.9032258064516,
    483.5483870967742,
    492.9032258064516,
    503.2258064516129,
    519.0322580645161,
    532.258064516129,
    548.3870967741935,
    558.7096774193549,
    568.7096774193549
]



time_ref_ = [
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
rel_power_ref_ = [
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



time_ref_1 = [
    0.7412060301507538,
    0.8165829145728644,
    0.9296482412060302,
    1.0678391959798996,
    1.1809045226130654,
    1.4824120603015076,
    1.8341708542713568,
    2.3366834170854274,
    2.864321608040201,
    3.341708542713568,
    3.7311557788944727,
    4.208542713567839,
    4.937185929648241
]
rel_power_ref_1 = [
    321.27334465195247,
    323.6502546689304,
    326.2818336162988,
    329.3378607809847,
    331.1205432937182,
    335.36502546689303,
    339.52461799660443,
    343.1748726655348,
    346.1460101867572,
    348.1833616298812,
    349.6264855687606,
    351.15449915110355,
    353.0220713073005
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
ax1.set_ylim(280, 760)
# Set y-axis ticks (200 to 500 by step of 50)
yticks = range(280, 761, 80)  # 501 is used to include 500
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

ax2.set_ylim(200, 500)
# Set y-axis ticks (200 to 500 by step of 50)
yticks = range(200, 501, 50)  # 501 is used to include 500
ax2.set_yticks(yticks)

ax2.set_xlim(right=5)  # x-axis starts from 0

# Set fontweight='bold' for both x and y ticks on ax1
for label in ax1.get_xticklabels() + ax1.get_yticklabels():
    label.set_fontweight('bold')

# Set fontweight='bold' for both x and y ticks on ax2
for label in ax2.get_yticklabels():
    label.set_fontweight('bold')

# Plot reference results of relative power
ax1.plot(time1, rel_power_ref1, color='black', marker='+', linestyle='--', linewidth=0.6, markersize=8, label='PANTHER')
ax1.plot(time3, rel_power_ref3, color='black', marker='*', linestyle='--', linewidth=0.6, markersize=8, label='QUANDRY')
# Plot reference results of reactivity
ax2.plot(time_ref_, rel_power_ref_, color='blue', marker='+', linestyle='--', linewidth=0.6, markersize=8, label='PANTHER')
ax2.plot(time_ref_1, rel_power_ref_1, color='blue', marker='*', linestyle='--', linewidth=0.6, markersize=8, label='QUANDRY')
#ax2.plot(-20,-20,color='chartreuse',  marker='^', linestyle='', linewidth=0.6, label='PANTHER', markersize=6)


# Add a combined legend
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc='upper left', fontsize=9)


# Title and save the plot
#plt.title('Relative Power and Reactivity Over Time', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('dual_y_axes_plot.png', dpi=500)
plt.show()

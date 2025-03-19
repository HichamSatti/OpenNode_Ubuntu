import numpy as np
import matplotlib.pyplot as plt

# Define file paths
#power_files1 = ['nea1x1.txt']
#power_files_1 = ['nea1x1x2.txt', 'nea1x1x4.txt']
power_files2 = ['nea2x2.txt', 'nea4x4.txt']
#reactivity_files1 = ['react1x4.txt']
reactivity_files2 = ['react2x2.txt', 'react4x4.txt']

#power_names1 = ['Radial Node 1x1']
#power_names_1 = ['Axial Node 1x2', 'Axial Node 1x4']
power_names2 = ['Radial Node 2x2', 'Radial Node 4x4']
#react_names1 = ['Axial Node 1x4']
react_names2 = ['Radial Node 2x2', 'Radial Node 4x4']

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




time2 = [
    0.01845924843085204,
    0.020593258136874633,
    0.022727267842897227,
    0.02806232792670738,
    0.03233034733875257,
    0.03766533578505538,
    0.04833545595267569,
    0.05473748507074347,
    0.06220655486057622,
    0.07287660339068919,
    0.09848486313797498,
    0.15503622780383472,
    0.215855683519247,
    0.26920614108233387,
    0.32469053671393594,
    0.3876440021353708,
    0.44739648881652544,
    0.4804737467161367,
    0.49967983407034
    ]
rel_power_ref2 = [
    1.0202051507509673,
    1.0227401109842154,
    1.0268485047267033,
    1.0315687847137451,
    1.0369883662003738,
    1.0443310216756967,
    1.0495757755866417,
    1.0534219294327956,
    1.0575303217080925,
    1.0605897632012,
    1.063736615547767,
    1.063736615547767,
    1.0632121401566725,
    1.062687664765578,
    1.06181353862469,
    1.0610268248044528,
    1.0602401124514067,
    1.0596282247396616,
    1.0592785739898682
    ]
time4=[
    0.5077262087888228,
    0.6622518666551853,
    0.783664671823055,
    0.91611460052712,
    0.9933774294603014
    ]
rel_power_ref4=[
    1.0591743560694824,
    1.0565610851681149,
    1.0549019598411595,
    1.0532428345142042,
    1.0521870256467787
    ]

time1 = [
    0.0093,
    0.02164,
    0.0354,
    0.0577,
    0.07726301419441967,
    0.17439325832871538,
    0.33222990504694594,
    0.4779252712483895,
    0.5871967958994722,
    0.6843270400337679,
    0.8178811257184245,
    0.975717772436655
    ]
rel_power_ref1= [
    1.0112297248140165,
    1.03,
    1.0427531060261679,
    1.0546988083802462,
    1.0618330472861541,
    1.0623307848842407,
    1.0600080094265032,
    1.0575193214360703,
    1.055528371043724,
    1.0543669833148552,
    1.0525419454552043,
    1.0510487326609446
    ]









# Read power files
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
color_index_3 = 0  # Initialize color index
for key, (time, rel_power_1) in power_data_1.items():
    ax1.plot(time, rel_power_1, label=key, linestyle='-.', marker='o',linewidth=1.2, markersize=1,
    color=colors__1[color_index_3 % len(colors__1)])
    color_index_3 += 1  # Update color index
color_index_2 = 0  # Initialize color index
for key, (time, rel_power2) in power_data2.items():
    ax1.plot(time, rel_power2, label=key, linestyle='-', marker='o',linewidth=1.2, markersize=1,
    color=colors_1[color_index_2 % len(colors_1)])
    color_index_2 += 1  # Update color index

ax1.set_xlabel('Time (s)', fontsize=14, fontweight='bold')
ax1.set_ylabel('Relative Power', fontsize=14, fontweight='bold', color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.grid(True, linestyle='--', alpha=0.7)
# Set the limits of the x and y axes to start from 0
ax1.set_xlim(left=0)  # x-axis starts from 0
ax1.set_ylim(bottom=1)  # Left y-axis starts from 0
ax1.set_ylim(top=1.07)  # Left y-axis starts from 0
ax1.set_xlim(right=1)  # x-axis starts from 0

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


ax2.set_ylabel('Reactivity ($)', fontsize=14, fontweight='bold', color='blue')
ax2.tick_params(axis='y', labelcolor='blue')
# Set the right y-axis to start from 0
ax2.set_ylim(bottom=0)
ax2.set_ylim(top=0.050)  # Left y-axis starts from 0

# Set fontweight='bold' for both x and y ticks on ax1
for label in ax1.get_xticklabels() + ax1.get_yticklabels():
    label.set_fontweight('bold')

# Set fontweight='bold' for both x and y ticks on ax2
for label in ax2.get_yticklabels():
    label.set_fontweight('bold')


# Plot reference results of relative power
ax1.plot(time2, rel_power_ref2, color='black', marker='+', linestyle=':', linewidth=0.6, label='PANTHER', markersize=8)
ax1.plot(time4, rel_power_ref4, color='black', marker='+', linestyle=':', linewidth=0.6, markersize=8)
#ax1.plot(time3, rel_power_ref3, color='black', marker='+', linestyle=':', linewidth=0.6, markersize=8)
ax1.plot(time1, rel_power_ref1, color='black', marker='*', linestyle=':', linewidth=0.6, label='Reference', markersize=6)
# Plot reference results of reactivity

# Add a combined legend
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc='upper right', fontsize=9)


# Title and save the plot
#plt.title('Relative Power and Reactivity Over Time', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('dual_y_axes_plot.png', dpi=500)
plt.show()

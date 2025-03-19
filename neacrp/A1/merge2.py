import numpy as np
import matplotlib.pyplot as plt

# Define file paths
power_files1 = ['maxf1x1.txt']
power_files2 = ['maxf4x4.txt', 'maxf2x2.txt']
reactivity_files1 = ['dop1x1.txt']
reactivity_files2 = ['dop2x2.txt', 'dop4x4.txt']

power_names1 = ['Radial Node 1x1']
#power_names_1 = ['Axial Node 1x2', 'Axial Node 1x4']
power_names2 = ['Radial Node 2x2', 'Radial Node 4x4']
react_names1 = ['Radial Node 1x4']
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


time2 = [
    0.22693997071742314,0.3879941434846267,
    0.4758418740849195,
    0.5051244509516838,
    0.5710102489019033,
    0.746705710102489,
    1.0322108345534406,
    1.4494875549048316,
    2.4158125915080526,
    2.957540263543192,
    3.49194729136164,
    4.062957540263543,
    4.773060029282577,4.992679355783309
]
rel_power_ref2 = [
    284.47081784386614,287.2589219330855,
    303.9875464684015,
    330.939219330855,
    381.12509293680296,
    446.18085501858735,
    498.2254646840148,
    549.3407063197026,
    628.3369888475836,
    657.1473977695167,
    678.5228624535316,
    698.0395910780669,
    717.5563197026022,723.1325278810409
]
time1 = [
    0.1390922401171303,
    0.3733528550512445,
    0.5124450951683748,
    0.5710102489019033,
    0.7320644216691069,
    1.2884333821376281,
    1.896046852122987,
    2.8184480234260616,
    3.7042459736456808,
    4.136163982430454,
    4.882869692532943,4.985358711566618
]
rel_power_ref1 = [
    284.47081784386614,
    285.40018587360595,
    323.5042750929368,
    367.1845724907063,
    429.45223048327136,
    513.0953531598512,
    569.7868029739777,
    622.760780669145,
    653.4299256505576,
    666.4410780669145,
    683.1697026022305,685.0284386617101
]
time3 = [0.022779043280182234,0.38724373576309795,
    0.5694760820045558,
    0.6719817767653758,
    0.7858769931662871,
    1.0820045558086562,
    1.5148063781321186,
    2.5170842824601367,
    3.029612756264237,
    4.009111617312073,
    4.4874715261959,
    4.9886104783599095
]
rel_power_ref3 = [285.8569051580699,287.52079866888516,
    300.83194675540767,
    350.7487520798669,
    391.5141430948419,
    453.0782029950083,
    503.82695507487523,
    580.3660565723794,
    606.1564059900167,
    640.2662229617305,
    653.5773710482529,
    664.3926788685524
]

time_ref_ = [
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
rel_power_ref_ = [
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
time_ref_1 = [
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
rel_power_ref_1 = [
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
for idx, file_path in enumerate(power_files1):
    try:
        name = power_names1[idx]
        data = np.genfromtxt(file_path, skip_header=2)
        time = data[:, 0]  # Time column
        relative_power1 = data[:, 1]  # Relative power column
        power_data1[name] = (time, relative_power1)
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
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
for idx_, file_path_ in enumerate(reactivity_files1):
    try:
        name = react_names1[idx_]
        data = np.genfromtxt(file_path_, skip_header=2)
        time = data[:, 0]  # Time column
        reactivity1 = data[:, 1]  # Reactivity column
        reactivity_data1[name] = (time, reactivity1)
    except Exception as e:
        print(f"Error reading file {file_path_}: {e}")
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
ax1.plot(time1, rel_power_ref1, color='black', marker='+', linestyle='--', linewidth=0.6, markersize=8, label='PARCS')
ax1.plot(time3, rel_power_ref3, color='black', marker='*', linestyle='--', linewidth=0.6, markersize=8, label='CORCA-K')
ax1.plot(time2, rel_power_ref2, color='black', marker='o', linestyle='--', linewidth=0.6, markersize=6, label='PANTHER')
# Plot reference results of reactivity
ax2.plot(time_ref_, rel_power_ref_, color='blue', marker='+', linestyle='--', linewidth=0.6, markersize=8, label='PARCS')
ax2.plot(time_ref_1, rel_power_ref_1, color='blue', marker='*', linestyle='--', linewidth=0.6, markersize=8, label='CORCA-K')
ax2.plot(time_ref_2, rel_power_ref_2, color='blue', marker='o', linestyle='--', linewidth=0.6, markersize=6, label='PANTHER')
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

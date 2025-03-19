import numpy as np
import matplotlib.pyplot as plt

# Define file paths
reactivity_files1 = ['dop1x1x1.txt', 'dop1x1x2.txt']
reactivity_files2 = ['dop2x2.txt', 'dop4x4.txt']
power_files1 = ['nea1x1.txt']
power_files2 = ['nea2x2.txt']

react_names1 = ['Radial Node 1x1', 'Axial Node 1x2']
#power_names_1 = ['Axial Node 1x2', 'Axial Node 1x4']
react_names2 = ['Radial Node 2x2', 'Radial Node 4x4']
power_names1 = ['Radial Node 1x1']
power_names2 = ['Radial Node 2x2']

# Initialize dictionaries to store data
reactivity_data1 = {}
reactivity_data2 = {}
power_data1 = {}
power_data_1 = {}
power_data2 = {}


colors_1 = ['black']
colors1 = ['black']
colors__1 = ['black']
colors2 = ['dodgerblue', 'blue']
colors_2 = ['dodgerblue','blue']

#bleu
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
time3 = [
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
rel_power_ref3 = [
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
time2 = [
    0.12779552715654952,
    0.3194888178913738,
    0.4792332268370607,
    0.8945686900958466,
    1.3498402555910542,
    1.7651757188498403,
    2.356230031948882,
    2.8354632587859423,
    3.2587859424920125,
    3.642172523961661,
    3.9536741214057507,
    4.464856230031948,
    4.760383386581469,
    4.976038338658147
]
rel_power_ref2 = [
    546.5576208178438,
    547.5613382899628,
    548.3420074349442,
    549.903345724907,
    551.1301115241636,
    552,
    552.9368029739777,
    553.449814126394,
    553.8289962825279,
    554.0743494423792,
    554.2750929368029,
    554.4981412639405,
    554.631970260223,
    554.6765799256506
]

#black
time_ref_ = [
    0.0390015600624025,
    0.4914196567862715,
    0.9048361934477379,
    1.2792511700468019,
    1.513260530421217,
    1.903276131045242,
    2.3400936037441498,
    2.784711388455538,
    3.291731669266771,
    3.7597503900156006,
    4.227769110764431,
    4.6801872074882995,
    4.968798751950078
]
rel_power_ref_ = [
    1672.4121301775147,
    1677.1458579881655,
    1680.6961538461537,
    1683.0630177514793,
    1685.4298816568046,
    1686.6133136094672,
    1688.9801775147928,
    1691.3470414201183,
    1691.3470414201183,
    1693.7139053254436,
    1693.7139053254436,
    1693.7139053254436,
    1693.7139053254436
]
time_ref_1 = [
    0.046801872074883,
    0.3666146645865835,
    0.842433697347894,
    1.4118564742589703,
    2.215288611544462,
    2.7535101404056164,
    3.3151326053042123,
    3.884555382215289,
    4.282371294851794,
    4.648985959438377,
    4.968798751950078
]
rel_power_ref_1 = [
    1871.228698224852,
    1875.9624260355029,
    1879.512721893491,
    1883.0630177514793,
    1888.9801775147928,
    1891.3470414201183,
    1893.7139053254436,
    1893.7139053254436,
    1893.7139053254436,
    1894.8973372781065,
    1894.8973372781065
]
time_ref_2 = [
    0.078003120124805,
    0.31201248049922,
    0.5694227769110765,
    0.7566302652106085,
    0.982839313572543,
    1.2090483619344774,
    1.5054602184087365,
    1.8018720748829953,
    2.1138845553822154,
    2.3946957878315134,
    2.6989079563182528,
    2.995319812792512,
    3.299531981279251,
    3.6037441497659906,
    3.9079563182527304,
    4.204368174726989,
    4.508580343213729,
    4.804992199687987
]
rel_power_ref_2 = [
    1694.8973372781065,
    1696.0807692307692,
    1699.6310650887574,
    1700.81449704142,
    1701.9979289940827,
    1704.3647928994083,
    1706.7316568047336,
    1707.9150887573965,
    1710.2819526627218,
    1712.6488165680473,
    1713.83224852071,
    1713.83224852071,
    1713.83224852071,
    1716.1991124260353,
    1716.1991124260353,
    1716.1991124260353,
    1716.1991124260353,
    1716.1991124260353
]
time_ref_3 = [
    0.0463821892393321,
    0.25974025974025977,
    0.49165120593692024,
    0.7792207792207793,
    1.103896103896104,
    1.4471243042671615,
    1.9944341372912802,
    2.4118738404452693,
    3.0241187384044528,
    3.4972170686456403,
    4.01669758812616,
    4.471243042671614,
    4.962894248608534
]
rel_power_ref_3 = [
    1671.9782608695652,
    1674.1847826086957,
    1676.6195652173913,
    1679.054347826087,
    1681.5652173913043,
    1683.8478260869565,
    1686.7391304347825,
    1688.2608695652175,
    1690.0869565217392,
    1690.9239130434783,
    1691.5326086956522,
    1691.8369565217392,
    1691.9130434782608
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
ax1.set_ylim(1640, 2000)
# Set y-axis ticks (200 to 500 by step of 50)
yticks = range(1600, 2001, 50)  # 501 is used to include 500
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

ax2.set_ylim(480, 620)
# Set y-axis ticks (200 to 500 by step of 50)
yticks = range(480, 621, 20)  # 501 is used to include 500
ax2.set_yticks(yticks)

ax2.set_xlim(right=5)  # x-axis starts from 0

# Set fontweight='bold' for both x and y ticks on ax1
for label in ax1.get_xticklabels() + ax1.get_yticklabels():
    label.set_fontweight('bold')

# Set fontweight='bold' for both x and y ticks on ax2
for label in ax2.get_yticklabels():
    label.set_fontweight('bold')

# Plot reference results of relative power
ax1.plot(time_ref_, rel_power_ref_, color='black', marker='+', linestyle='--', linewidth=0.6, markersize=8, label='PARCS')
#ax1.plot(time_ref_1, rel_power_ref_1, color='black', marker='*', linestyle='--', linewidth=0.6, markersize=8, label='CORCA-K(Pin)')
#ax1.plot(time_ref_2, rel_power_ref_2, color='black', marker='^', linestyle='--', linewidth=0.6, markersize=6, label='PANTHER(Nodal)')
ax1.plot(time_ref_3, rel_power_ref_3, color='black', marker='o', linestyle='--', linewidth=0.6, markersize=6, label='PANTHER')
# Plot reference results of reactivity
ax2.plot(time1, rel_power_ref1, color='blue', marker='+', linestyle='--', linewidth=0.6, markersize=8, label='PARCS')
ax2.plot(time3, rel_power_ref3, color='blue', marker='*', linestyle='--', linewidth=0.6, markersize=8, label='CORCA-K')
ax2.plot(time2, rel_power_ref2, color='blue', marker='o', linestyle='--', linewidth=0.6, markersize=6, label='PANTHER')
#ax2.plot(-20,-20,color='chartreuse',  marker='^', linestyle='', linewidth=0.6, label='PANTHER', markersize=6)


# Add a combined legend
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc='upper right', fontsize=9)


# Title and save the plot
#plt.title('Relative Power and Reactivity Over Time', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('dual_y_axes_plot.png', dpi=500)
plt.show()

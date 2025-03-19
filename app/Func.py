#!usr/bin/python3
# -*- coding: utf-8 -*-  
import json
import numpy 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patches as mpatch
from matplotlib.patches import Ellipse, Arc
import mpl_toolkits.mplot3d.art3d as art3d
import math 
import matplotlib.pyplot
import matplotlib
matplotlib.use('Qt5Agg')
from pathlib import Path
import matplotlib.lines as mlines

# Initialize empty lists for sizes
xsize = []
ysize = []
zsize = []

# Read directory and method files
try:
    with open('app/link/script.dir', "r") as f:
        Drct = f.read().strip()
    with open('app/link/script00.py', "r") as f:
        Methods = f.read().strip()
except FileNotFoundError as e:
    print(f"Error: {e}")
    exit(1)

FileName = Path(Drct).stem

# Load JSON data if method is 'NEM'
if Methods == 'NEM':
    try:
        with open(Drct) as json_data:
            data = json.load(json_data)
            params = data.get('Data', {}).get('Parameters', {})
            xsize = params.get('X-Assembly Size', [])
            ysize = params.get('Y-Assembly Size', [])
            zsize = params.get('Z-Assembly Size', [])

            if not xsize or not ysize or not zsize:
                raise ValueError("One or more size lists are empty in JSON data.")
    except FileNotFoundError:
        print(f"Error: File not found at {Drct}")
        exit(1)
    except (json.JSONDecodeError, KeyError, ValueError) as e:
        print(f"Error reading JSON: {e}")
        exit(1)

# Reverse ysize if it has elements
if ysize:
    ysize.reverse()

def add_one_by_one_gen(l):
    cumsum = 0
    for i in l:
        cumsum += i
        yield cumsum

# Ensure xsize and ysize have elements before accessing indices
if xsize:
    xsize[0] *= 2
else:
    raise ValueError("xsize is empty, cannot modify its first element.")

# Mirror elements around the pivot in width_x
pivot_x = xsize[0]
mirrored_width_x = xsize[1:] if len(xsize) > 1 else []
xsize = [0.0] + mirrored_width_x[::-1] + [pivot_x] + mirrored_width_x

if ysize:
    if FileName in ('IAEA_2D', 'Ad_IAEA_2D', 'LWM'):
        ysize[-1] *= 2
        pivot_y = ysize[-1]
        mirrored_width_y = ysize[:-1] if len(ysize) > 1 else []
    else:
        ysize[0] *= 2
        pivot_y = ysize[0]
        mirrored_width_y = ysize[1:] if len(ysize) > 1 else []
    
    ysize = [0.0] + mirrored_width_y[::-1] + [pivot_y] + mirrored_width_y
else:
    raise ValueError("ysize is empty, cannot modify its elements.")

# Ensure zsize is not empty before modifying
if zsize:
    zsize = [0.0] + zsize
else:
    raise ValueError("zsize is empty, cannot modify it.")

# Generate cumulative sums
xs = list(add_one_by_one_gen(xsize))
zs = list(add_one_by_one_gen(zsize))

if ysize != xsize:
    ysize.reverse()
ys = list(add_one_by_one_gen(ysize))





def AxPower1(self):
    data = np.loadtxt('app/Output/AxialPower.out')
    max_columns = len(data[0]) - 1
    max_rows = len(data) - 1
    x = [data[0][colnum + 1] for colnum in range(max_columns)]
    y = [data[rownum+1][0] for rownum in range(max_rows)] 
    z = [[data[rownum+1][colnum + 1] for rownum in range(max_rows)] for colnum in range(max_columns)]
    x = np.array(x).T
    y = np.array(y).T
    z = np.array(z).T
    x,y = np.meshgrid(x, y)
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)
    ax.set_ylabel('Core Height(cm)') 
    im = ax.imshow(z, interpolation='bilinear', cmap='jet', 
             origin='lower', extent=[0, abs(x).max(), 0,abs(y).max()])
    clb = fig.colorbar(im)
    clb.set_label('Power Density')
    plt.title("Axial Power Distribution")
#    plt.savefig('2Dimage.png',dpi=500)

    # Show the plot.
    plt.show()
    
import numpy as np
import matplotlib.pyplot as plt

def AxPower2(self):
    import numpy as np
    import matplotlib.pyplot as plt

    # Load data from the file
    try:
        data = np.loadtxt('app/Output/AxialPower.out')
    except Exception as e:
        print(f"Error loading data: {e}")
        return

    # Process the data
    if data.ndim != 2 or data.shape[1] < 2:
        print("Error: Data file must contain at least two columns.")
        return

    x = data[1:, 0]  # Core Height
    y = data[1:, 1]  # Power

    # Create the plot
    plt.figure(figsize=(8, 6))

    # Plot the data points with lines connecting them
    plt.plot(x, y, marker='o', linestyle='-', color='b', linewidth=1.5, markersize=6)

    # Set axis limits correctly
    plt.xlim(0, 400)
    plt.ylim(-0.1, 2)
    plt.yticks(np.arange(0, 2.1, 1))

    # Add labels and title
    plt.xlabel('Core Height (cm)', fontsize=12)
    plt.ylabel('Power Density', fontsize=12)
    plt.title('Axial Power Distribution', fontsize=14, fontweight='bold')

    # Add grid and adjust layout
    plt.grid(True)
    plt.tight_layout()

    # Save the plot with high resolution
    plt.savefig('2Dimage.png', dpi=500)

    # Display the plot
    plt.show()


#def AxPower2(self):
#    import numpy as np
#    import matplotlib.pyplot as plt
#
#    # Load OpenNode data
#    try:
#        data_opennode = np.loadtxt('app/Output/AxialPower.out')
#    except Exception as e:
#        print(f"Error loading OpenNode data: {e}")
#        return
#
#    # Load KOMODO data
#    try:
#        data_komodo = np.loadtxt('app/Output/KOMODO_AxialPower.out')
#    except Exception as e:
#        print(f"Error loading KOMODO data: {e}")
#        return
#
#    # Process the data for both OpenNode and KOMODO
#    if data_opennode.ndim != 2 or data_opennode.shape[1] < 2:
#        print("Error: OpenNode data file must contain at least two columns.")
#        return
#
#    if data_komodo.ndim != 2 or data_komodo.shape[1] < 2:
#        print("Error: KOMODO data file must contain at least two columns.")
#        return
#
#    # Extract data columns for both sets
#    x_opennode = data_opennode[1:, 0]  
#    y_opennode = data_opennode[1:, 1]  
#
#    x_komodo = data_komodo[1:, 0]  
#    y_komodo = data_komodo[1:, 1]  
#
#    # Create the plot for comparison
#    plt.figure(figsize=(8, 6))
#
#    # Plot OpenNode data
#    plt.plot(x_opennode, y_opennode, marker='o', linestyle='-', color='b', linewidth=1.5, markersize=6, label='OpenNode')
#
#    # Plot KOMODO data
#    plt.plot(x_komodo, y_komodo, marker='s', linestyle='--', color='r', linewidth=1.5, markersize=6, label='KOMODO')
#
#    # Set axis limits and ticks
#    plt.xlim(0, 400)
#    plt.ylim(-0.1, 10)
#    plt.yticks(np.arange(0, 10.1, 1))
#
#    # Add labels, title, and legend
#    plt.xlabel('Core Height (cm)', fontsize=12)
#    plt.ylabel('Power Density', fontsize=12)
#    plt.title('Axial Power Distribution Comparison', fontsize=14, fontweight='bold')
#    plt.legend(fontsize=12)
#
#    # Add grid and adjust layout
#    plt.grid(True)
#    plt.tight_layout()
#
#    # Save the plot with high resolution
#    plt.savefig('AxialPower_Comparison.png', dpi=500)
#
#    # Display the plot
#    plt.show()
    
def extend_and_mirror_layer(layer):
    # Reverse each row
    reversed_rows_matrix = [row[::-1] for row in layer]
    
    # Extend each row by appending its mirror image, excluding the first value
    mirrored_matrix = [row[:-1] + row[::-1] for row in reversed_rows_matrix]
    
    # Combine the top half and the bottom half to form the full matrix
    full_matrix = mirrored_matrix + mirrored_matrix[::-1][1:]

    return full_matrix
def process_plannar(plannar):
    extended_plannar_layers = [extend_and_mirror_layer(layer) for layer in plannar]
    return extended_plannar_layers


def RdPower(self):
    # Load the data from the file
    data = np.loadtxt('app/Output/RadialPower.out')
    
    # Extract the coordinates
    x_coords = data[0, 1:]  # X-coordinates from the first row (excluding the first element)
    y_coords = data[1:, 0]  # Y-coordinates from the first column (excluding the first element)
    
    # Extract the power values
    power_values = data[1:, 1:]  # Power values excluding the first row and column
    
    # Step 1: Horizontal mirroring (mirror along the y-axis)
    half_distribution = np.hstack((power_values[:, ::-1], power_values))
    
    # Step 2: Vertical mirroring (mirror along the x-axis)
    full_distribution = np.vstack((half_distribution[::-1, :], half_distribution))
    
    # Extend coordinates for the full distribution
    extended_x_coords = np.hstack((-x_coords[::-1], x_coords))
    extended_y_coords = np.hstack((-y_coords[::-1], y_coords))
    
    # Plotting the full radial power distribution
    fig, ax = plt.subplots()
    cax = ax.imshow(full_distribution, interpolation='bilinear', cmap='jet', origin='lower',
                    extent=[extended_x_coords.min(), extended_x_coords.max(), extended_y_coords.min(), extended_y_coords.max()])
    
    # Add colorbar
    clb = fig.colorbar(cax)
    clb.set_label('Power Density')

    # Set labels and title
    ax.set_xlabel('X [cm]')
    ax.set_ylabel('Y [cm]')
    ax.set_title('Full Radial Power Distribution')

    # Set custom ticks for the X and Y axes
    ax.set_xticks(extended_x_coords)
    ax.set_xticklabels([f'{x:.1f}' for x in extended_x_coords], rotation=45, ha='right')
    ax.set_yticks(extended_y_coords)
    ax.set_yticklabels([f'{y:.1f}' for y in extended_y_coords])
    plt.savefig('2Dimage.png',dpi=500)

    plt.show()

def Flux(self):
    with open(Drct) as json_data:
         datas = json.load(json_data)
         ng = datas['Data']['Parameters']['Number of Energy Groups']         
    data = []
    # M00 = open('app/link/script00.py', "r" ).read()
    # if M00 == 'NEM':
        # data = np.loadtxt('app/Output/FluxGr1.out')
    # else:
        # QMessageBox.warning(self, "Warning", "select the calculation method.")

    for i in range(ng):   
        # Load the data from the file
        data = np.loadtxt('app/Output/FluxGr '+str(i+1))
        # Extract the coordinates
        x_coords = data[0, 1:]  # X-coordinates from the first row (excluding the first element)
        y_coords = data[1:, 0]  # Y-coordinates from the first column (excluding the first element)
        # Extract the power values
        power_values = data[1:, 1:]  # Power values excluding the first row and column
        # Step 1: Horizontal mirroring (mirror along the y-axis)
        half_distribution = np.hstack((power_values[:, ::-1], power_values))
        
        # Step 2: Vertical mirroring (mirror along the x-axis)
        full_distribution = np.vstack((half_distribution[::-1, :], half_distribution))
        
        # Extend coordinates for the full distribution
        extended_x_coords = np.hstack((-x_coords[::-1], x_coords))
        extended_y_coords = np.hstack((-y_coords[::-1], y_coords))
        # Plotting the full radial power distribution
        fig, ax = plt.subplots()
        cax = ax.imshow(full_distribution, interpolation='bilinear', cmap='jet', origin='lower',
                        extent=[extended_x_coords.min(), extended_x_coords.max(), extended_y_coords.min(), extended_y_coords.max()])
        # Add colorbar
        clb = fig.colorbar(cax)
        clb.set_label('Flux Density')

        # Set labels and title
        ax.set_xlabel('X [cm]')
        ax.set_ylabel('Y [cm]')
        ax.set_title('Energy Group '+str(i+1),  loc='right', fontsize=10)
        
        # Set custom ticks for the X and Y axes
        ax.set_xticks(extended_x_coords)
        ax.set_xticklabels([f'{x:.1f}' for x in extended_x_coords], rotation=45, ha='right')
        ax.set_yticks(extended_y_coords)
        ax.set_yticklabels([f'{y:.1f}' for y in extended_y_coords])

        plt.title('Radial Flux Distribution', loc='left')
        plt.savefig(f'2Dimage_Group{i+1}.png', dpi=500)  # Save with high resolution (500 dpi)
        plt.show()




def Visualization2D(self):
    assm =[[]]
    plannar = []
    pin = []
    regmat = []
    nom    = []
    width_x = []
    width_y = []
    with open(Drct) as json_data:
        data = json.load(json_data)
        nmat = data['Data']['Parameters']['Number of Materials']
        nx1 = data['Data']['Parameters']['Number of X-Assembly']
        ny1 = data['Data']['Parameters']['Number of Y-Assembly']
        np  = data['Data']['Parameters']['Number of Planar'] 
        # xsize = data['Data']['Parameters']['X-Assembly Size']
        # ysize = data['Data']['Parameters']['Y-Assembly Size']
        plannar = data['Data']['XY_Assembly']
        
        if FileName in ('IAEA_2D', 'Ad_IAEA_2D', 'LWM'):
            plannar = [array[::-1] for array in plannar]

        nx = nx1*2 - 1
        ny = ny1*2 - 1
        # Process the plannar list
        extended_plannar = process_plannar(plannar)
     #   width_x = data['Data']['Parameters']['X-Assembly Size']
     #   width_y = data['Data']['Parameters']['Y-Assembly Size']
        
        # if FileName == ('IAEA_2D' or 'Ad_IAEA_2D'):
            # plannar[0].reverse()

        # for i in range(np):
            # plannar.append(data['Data']['Assemblies'][i]['Assembly'])
        # for p in range(np):
        # for i in range(npc):
            # pin.append(data['data']['PinCells'][i]['mat_fill'])




        width_x1 = data['Data']['Parameters']['X-Assembly Size']
        width_y1 = data['Data']['Parameters']['Y-Assembly Size']
        
        # Multiply the first element in width_x and the last element in width_y by 2
        width_x1[0] *= 2
        width_y1[-1] *= 2
        
        # Mirror the elements around the pivot in width_x
        pivot_x = width_x1[0]
        mirrored_width_x = width_x1[1:]  # Elements to the right of the pivot
        width_x1 = mirrored_width_x[::-1] + [pivot_x] + mirrored_width_x  # Mirroring
        
        # Mirror the elements around the pivot in width_y
        pivot_y = width_y1[-1]
        mirrored_width_y = width_y1[:-1]  # Elements to the left of the pivot
        width_y1 = mirrored_width_y[::-1] + [pivot_y] + mirrored_width_y  # Mirroring

        for p in range(np) :
            Height = len(extended_plannar[p])
            width = len(extended_plannar[p])
            assm =numpy.zeros(((Height),(width)),dtype='i')

            i1=0
            # for j in range(len(core[0])):
            for k in range(len(extended_plannar[p])):
                    # for m in range(len(pin[0])):
                i2=0
                        # for i in range(len(core[0])):
                for l in range(len(extended_plannar[p])):
                                # for n in range(len(pin[0])):
                    assm[i1][i2] = extended_plannar[p][k][l]
                    i2+=1
                i1+=1
            # nx = len(np.amax(data['Data']['Core'],axis=0))
            # ny = len(np.amax(data['Data']['Core'],axis=1))
                
            nxx = len(extended_plannar[p])
            nyy = len(extended_plannar[p])
            NX = nxx
            NY = nyy
        #    width_x.append([21.606,21.606,21.606,21.606,21.606,21.606,21.606,21.606,21.606,21.606,21.606,21.606,21.606,21.606,21.606,21.606,21.606])
        #    width_y.append([21.606,21.606,21.606,21.606,21.606,21.606,21.606,21.606,21.606,21.606,21.606,21.606,21.606,21.606,21.606,21.606,21.606])
            width_x.append(width_x1)
            width_y.append(width_y1)

            # nx =  len(width_x[0])*NX
            # ny =  len(width_y[0])*NY
            xcm  = width_x[0]*NX
            ycm  = width_x[0]*NY
            for j in range(nmat):
                nom.append(data['Data']['Materials'][j]['Name'])
            # for i in range(npc):
                # regmat.append(data['data']['PinCells'][i]['mat_fill'])
            fig, ax = plt.subplots()
            widthx = [0]
            widthy = [0]
            somx = [0]
            somy = [0]
            for i in range(nx):
                widthx.append(xcm[i])
                somx.append(sum(widthx))
            for j in range(ny):
                widthy.append(ycm[j])
                somy.append(sum(widthy))
            rectangles = []
            red_patch = []
            mx = somx[0]
            my = somy[0]

            # Draw lines between assemblies
            for i in range(nx + 1):
                ax.axvline(x=somx[i], color='k', linestyle='--', linewidth=0.5)
            for j in range(ny + 1):
                ax.axhline(y=somy[j], color='k', linestyle='--', linewidth=0.5)
            
            
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.72, box.height])
            for i in range(nx):
                for j in range(ny):
                    rectangles.append(mpatch.Rectangle((somx[i],somy[j]), xcm[i], ycm[j],linewidth=0.0,edgecolor='k'))
          #  colr = ['#33A1C9','#DC143C','#00C78C','#FFEC8B','#00FFFF','#FA8072', 

            if FileName in ('A1', 'A2', 'B1', 'B2', 'A1t', 'A2t', 'B1t', 'B2t', 'LWM'):
                colr = ['royalblue', 'gold', 'darkOrange', 'Red', 'Pink', 'deepskyblue', 
                        'Limegreen', 'Brown', 'dodgerblue', 'SeaGreen', 'SkyBlue',
                        '#ffb3e6', '#c2c2f0', '#ffb3e6', '#c2c2f0', '#ffb3e6']
            else:
                colr = ['Blue', 'Red', 'Green', 'Yellow', 'SkyBlue', 'Pink', 
                        '#FFC1C1', '#ACFFC7', 'Black', '#c2c2f0', '#ffb3e6', 
                        '#c2c2f0', '#ffb3e6', '#c2c2f0', '#ffb3e6', 'White']
                        
            #colr = ['#00C78C','#DC143C','#33A1C9','#FFEC8B','#00FFFF','#FA8072', 
            #        '#FFC1C1','#ACFFC7','Black', '#c2c2f0', '#ffb3e6', 
            #        '#c2c2f0','#ffb3e6', '#c2c2f0','#ffb3e6','White']
            

            for i in range(nmat):
                red_patch.append(mpatches.Patch(color=colr[i], label=nom[i]))
            n=nx-1
            m=0
            for r in rectangles:
                ax.add_artist(r) 
                if assm[n][m]== 0 :
                    r.set_facecolor(color='white')
                else :
                    r.set_facecolor(color=colr[int(assm[n][m])-1])

                rx, ry = r.get_xy()
                cx = rx + r.get_width()/2.0
                cy = ry + r.get_height()/2.0
                n-=1
                if (n<0):
                    n=nx-1
                    m+=1
            ax.set_ylim((min(somy), max(somy)))
            ax.set_xlim((min(somx), max(somx)))
            #ax.set_xticklabels([])
            #ax.set_yticklabels([]) 
            ax.set_xlabel('X [cm]')
            ax.set_ylabel('Y [cm]')
            ax.set_title('Planar '+str(p+1),  loc='left', fontsize=10)
            ax.set_title('Radial Geometry',  loc='center')
            clb = plt.legend(handles=red_patch, loc='center left', title="Materials", fontsize=6.8, bbox_to_anchor=(1, 0.5))
            plt.xticks(xs, rotation=45)
            plt.yticks(ys)
            plt.savefig(f'2Dimage_planar{p+1}.png', dpi=500)  # Save with high resolution (500 dpi)
            plt.show()

def AxiGeo2D(self):
    assm =[[]]
    plannar = []
    pin = []
    regmat = []
    nom    = []
    width_x = []
    width_y = []
    with open(Drct) as json_data:
        data = json.load(json_data)
        nmat = data['Data']['Parameters']['Number of Materials']
        nx1 = data['Data']['Parameters']['Number of X-Assembly']  # Number of Pin Cell
        nz = data['Data']['Parameters']['Number of Z-Assembly']  # Number of Pin Cell
        np  = data['Data']['Parameters']['Number of Planar'] 
        # xsize = data['Data']['Parameters']['X-Assembly Size']
        # ysize = data['Data']['Parameters']['Y-Assembly Size']
        # plannar = data['Data']['Core2']

        plannar.append(data['Data']['Z_Assembly'])
        extended_plannar = process_plannar(plannar)

        nx = nx1*2 - 1
        width_x1 = data['Data']['Parameters']['X-Assembly Size']
        
        # Multiply the first element in width_x and the last element in width_y by 2
        width_x1[0] *= 2
        
        # Mirror the elements around the pivot in width_x
        pivot_x = width_x1[0]
        mirrored_width_x = width_x1[1:]  # Elements to the right of the pivot
        width_x1 = mirrored_width_x[::-1] + [pivot_x] + mirrored_width_x  # Mirroring
        
        width  = nx
        Height = nz
        assm =numpy.zeros(((Height),(width)),dtype='i')

        i1=0
        for k in range(nz):
            i2=0
            for l in range(nx):
                assm[i1][i2] = extended_plannar[0][k][l]
                i2+=1
            i1+=1                

        width_x.append(width_x1)
        width_y.append(data['Data']['Parameters']['Z-Assembly Size'])
        xcm  = width_x[0]*nx
        ycm  = width_y[0]*nz
        for j in range(nmat):
            nom.append(data['Data']['Materials'][j]['Name'])
        fig, ax = plt.subplots(figsize=(10, 8))
        widthx = [0]
        widthy = [0]
        somx = [0]
        somy = [0]
        for i in range(nx):
            widthx.append(xcm[i])
            somx.append(sum(widthx))
        for j in range(nz):
            widthy.append(ycm[j])
            somy.append(sum(widthy))
        rectangles = []
        red_patch = []
        mx = somx[0]
        my = somy[0]

        # Draw lines between assemblies
        for i in range(nx + 1):
            ax.axvline(x=somx[i], color='k', linestyle='--', linewidth=0.5)
        for j in range(nz + 1):
            ax.axhline(y=somy[j], color='k', linestyle='--', linewidth=0.5)
            
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
        for i in range(nx):    
            for j in range(nz):

                rectangles.append(mpatch.Rectangle((somx[i],somy[j]), xcm[i], ycm[j],linewidth=0.0,edgecolor='k'))
        #colr = ['#33A1C9','#DC143C','#00C78C','#FFEC8B','#00FFFF','#FA8072', 
        #            '#FFC1C1','#ACFFC7','Black', '#c2c2f0', '#ffb3e6',
        #                '#ffb3e6', '#c2c2f0','#ffb3e6', '#c2c2f0','#ffb3e6']
       # colr = ['Blue','Red','Green','Yellow','SkyBlue','Pink', 
       #         'White', 'Black','#c2c2f0','#ffb3e6', '#c2c2f0',
       #         '#ffb3e6', '#c2c2f0','#ffb3e6', '#c2c2f0','#ffb3e6']
       # colr = ['#33A1C9','#DC143C','#00C78C','#FFEC8B','#00FFFF','#FA8072', 
       #             '#FFC1C1','#ACFFC7','Black', '#c2c2f0', '#ffb3e6',
       #                 '#ffb3e6', '#c2c2f0','#ffb3e6', '#c2c2f0','#ffb3e6']
   #+     colr = ['royalblue','gold','darkOrange','Red','Pink','deepskyblue', 
           # colr = ['Blue','Red','Green','Yellow','SkyBlue','Pink', 
           #     'Limegreen', 'Brown','dodgerblue','SeaGreen', 'SkyBlue',
           #     '#ffb3e6', '#c2c2f0','#ffb3e6', '#c2c2f0','#ffb3e6']

        if FileName in ('A1', 'A2', 'B1', 'B2', 'A1t', 'A2t', 'B1t', 'B2t', 'LWM'):
            colr = ['royalblue', 'gold', 'darkOrange', 'Red', 'Pink', 'deepskyblue', 
                    'Limegreen', 'Brown', 'dodgerblue', 'SeaGreen', 'SkyBlue',
                    '#ffb3e6', '#c2c2f0', '#ffb3e6', '#c2c2f0', '#ffb3e6']
        else:
            colr = ['Blue', 'Red', 'Green', 'Yellow', 'SkyBlue', 'Pink', 
                    '#FFC1C1', '#ACFFC7', 'Black', '#c2c2f0', '#ffb3e6', 
                    '#c2c2f0', '#ffb3e6', '#c2c2f0', '#ffb3e6', 'White']

        for i in range(nmat):
            red_patch.append(mpatches.Patch(color=colr[i], label=nom[i]))
        n=nz-1
        m=0
        for r in rectangles:
            ax.add_artist(r) 
            if assm[n][m]== 0 :
                r.set_facecolor(color='white')
            else :
                r.set_facecolor(color=colr[int(assm[n][m])-1])

            rx, ry = r.get_xy()
            cx = rx + r.get_width()/2.0
            cy = ry + r.get_height()/2.0
            n-=1
            if (n<0):
                n=nz-1
                m+=1
        ax.set_ylim((min(somy), max(somy)))
        ax.set_xlim((min(somx), max(somx)))
                #ax.set_xticklabels([])
                #ax.set_yticklabels([]) 
        ax.set_xlabel('X [cm]')
        ax.set_ylabel('Z [cm]')
        ax.set_title('Axial Geometry ',  loc='center')
        clb = plt.legend(handles=red_patch, loc='center left', title="Materials", fontsize='small', bbox_to_anchor=(1, 0.5))
        plt.xticks(xs, rotation=45)
        plt.yticks(zs)
        plt.savefig('2Dimage.png',dpi=500)
        plt.show()






def VisualizationCR(self):
    assm = [[]]
    plannar = []
    bmap = []
    bpos = []
    pin = []
    regmat = []
    nom = []
    width_x = []
    width_y = []

    # Define unique colors for each distinct bpos value
    colors = ['Lime', 'cyan', 'Fuchsia', 'lightcoral', 'orange', 'brown', 'cyan']  # Add more colors if needed
    bpos_to_color = {}

    # Load data
    with open(Drct) as json_data:
        data = json.load(json_data)
        nmat = data['Data']['Parameters']['Number of Materials']
        nx = data['Data']['Parameters']['Number of X-Assembly']
        ny = data['Data']['Parameters']['Number of Y-Assembly']
        npp = data['Data']['Parameters']['Number of Planar']
        plannar = data['Data']['XY_Assembly']
        bmap = data['Data']['CBCSearch']['Radial CR bank map']
        bpos = data['Data']['CBCSearch']['CR bank position']
        fbpos = data['Data']['REject']['Final CR bank position']

        if FileName in ('IAEA_2D', 'Ad_IAEA_2D', 'LWM'):
            plannar = [array[::-1] for array in plannar]
        if FileName in ('IAEA_2D', 'Ad_IAEA_2D', 'LWM'):
            bmap = [array[::-1] for array in bmap]

        # Create mapping from bpos values to colors
        unique_bpos = sorted(set(bpos))
        color_index = 0
        for value in unique_bpos:
            bpos_to_color[str(value)] = colors[color_index % len(colors)]
            color_index += 1

        # Detect differences between bpos and fbpos
        differences = [i for i in range(len(bpos)) if bpos[i] != fbpos[i]]

        bmap_array = np.array(bmap)
        reversed_horizontal = np.fliplr(bmap_array)
        reversed_vertical = np.flipud(bmap_array)
        reversed_both = np.flipud(np.fliplr(bmap_array))

        for p in range(npp):
            Height = len(plannar[p])
            width = len(plannar[p])
            assm = np.zeros((Height, width), dtype='i')

            i1 = 0
            for k in range(len(plannar[p])):
                i2 = 0
                for l in range(len(plannar[p])):
                    assm[i1][i2] = plannar[p][k][l]
                    i2 += 1
                i1 += 1

            nxx = len(plannar[p])
            nyy = len(plannar[p])
            NX = nxx
            NY = nyy
            width_x.append(data['Data']['Parameters']['X-Assembly Size'])
            width_y.append(data['Data']['Parameters']['Y-Assembly Size'])
            xcm = width_x[0] * NX
            ycm = width_x[0] * NY

            for j in range(nmat):
                nom.append(data['Data']['Materials'][j]['Name'])

            fig, ax = plt.subplots()
            widthx = [0]
            widthy = [0]
            somx = [0]
            somy = [0]
            for i in range(nx):
                widthx.append(xcm[i])
                somx.append(sum(widthx))
            for j in range(ny):
                widthy.append(ycm[j])
                somy.append(sum(widthy))

            rectangles = []
            mx = somx[0]
            my = somy[0]

            # Draw lines between assemblies
            for i in range(nx + 1):
                ax.axvline(x=somx[i], color='k', linestyle='--', linewidth=0.5)
            for j in range(ny + 1):
                ax.axhline(y=somy[j], color='k', linestyle='--', linewidth=0.5)

            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.72, box.height])

            for i in range(nx):
                for j in range(ny):
                    rectangles.append(
                        mpatch.Rectangle(
                            (somx[i], somy[j]),
                            xcm[i],
                            ycm[j],
                            linewidth=0.0,
                            edgecolor='k',
                        )
                    )

            colr = [
                'royalblue',
                'gold',
                'darkOrange',
                'Red',
                'Pink',
                'deepskyblue',
                'Limegreen',
                'Brown',
                'dodgerblue',
                'SeaGreen',
                'SkyBlue',
            ]

            legend_entries = []
            color_to_values = {color: set() for color in bpos_to_color.values()}
            ejected_entries = []  # To store entries for "CA to be ejected" labeled with 'x'

            n = nx - 1
            m = 0
            for r in rectangles:
                ax.add_artist(r)
                if assm[n][m] == 0:

                    r.set_facecolor(color='white')
                else:
                    r.set_facecolor(color=colr[int(assm[n][m]) - 1])

                # Get the center of the rectangle
                rx, ry = r.get_xy()
                cx = rx + r.get_width() / 2.0
                cy = ry + r.get_height() / 2.0

                # Add colored points for CR bank types
                if 0 <= n < len(reversed_both) and 0 <= m < len(reversed_both[n]):
                    bmap_value = reversed_both[n][m]
                    if 0 < bmap_value <= len(bpos):  # Note: 0 is excluded for bpos as index
                        bpos_value = bpos[bmap_value - 1]
                        color = bpos_to_color.get(str(bpos_value), 'grey')
                        ax.plot(cx, cy, 'o', color=color, markersize=10)
                        color_to_values[color].add(bpos_value)

                        # Highlight differences with a marker and label for ejected
                        if bmap_value - 1 in differences:
                            ax.plot(cx, cy, marker='x', color='black', markersize=12)
                            ejected_entries.append("CA to be ejected")

                n -= 1
                if n < 0:
                    n = nx - 1
                    m += 1

            # Create legend entries for colors and their corresponding values
            legend_entries = []
            for color, values in color_to_values.items():
                # Adding the label with both positions and ejected label (only once)
                positions_label = f'Position in Steps: {", ".join(map(str, sorted(values)))}'
                legend_entries.append(
                    mlines.Line2D(
                        [], [], color=color, marker='o', linestyle='None', markersize=10, label=positions_label
                    )
                )

            # Add a single entry for "CA to be ejected" if it exists
            if ejected_entries:
                ejected_entry = mlines.Line2D([], [], color='black', marker='x', linestyle='None', markersize=12, label="CA to be ejected")
                clb = plt.legend(
                    handles=legend_entries + [ejected_entry],
                    loc='center left',
                    title="CR Bank Types",
                    fontsize=8.7,
                    bbox_to_anchor=(1, 0.5),
                )
            else:
                clb = plt.legend(
                    handles=legend_entries,
                    loc='center left',
                    title="CR Bank Types",
                    fontsize=8.7,
                    bbox_to_anchor=(1, 0.5),
                )

            ax.set_ylim((min(somy), max(somy)))
            ax.set_xlim((min(somx), max(somx)))
            ax.set_xlabel('X [cm]')
            ax.set_ylabel('Y [cm]')
            ax.set_title('Planar ' + str(p + 1), loc='left', fontsize=10)
            ax.set_title('Radial Geometry', loc='center')
            plt.savefig(f'2Dimage_planar{p+1}.png', dpi=500)
            plt.show()

def AxiCR(self):
    # Initialize variables
    assm = np.zeros((10, 6), dtype='i')  # Adjust dimensions as needed
    plannar = []
    nom = []
    width_x = []
    width_y = []

    # Define colors and labels for CR banks
    colors = ['Lime', 'cyan', 'Fuchsia', 'lightcoral', 'orange', 'brown', 'cyan']
    color_labels = ['Bank 1', 'Bank 2', 'Bank 3', 'Bank 4', 'Bank 5', 'Bank 6', 'Bank 7']

    # Read data from file
    try:
        with open(Drct) as json_data:
            data = json.load(json_data)
            nmat = data['Data']['Parameters']['Number of Materials']
            nx = data['Data']['Parameters']['Number of X-Assembly']
            nz = data['Data']['Parameters']['Number of Z-Assembly']
            bpos = data['Data']['CBCSearch']['CR bank position']
            fbpos = data['Data']['REject']['Final CR bank position']
            
            plannar.append(data['Data']['Z_Assembly'])
            assm = np.array(plannar[0])

            width_x.append(data['Data']['Parameters']['X-Assembly Size'])
            width_y.append(data['Data']['Parameters']['Z-Assembly Size'])
            xcm = width_x[0]
            ycm = width_y[0]
            
            for j in range(nmat):
                nom.append(data['Data']['Materials'][j]['Name'])
            
            bmap = data['Data']['CBCSearch']['Radial CR bank map']
            if FileName not in ('IAEA_2D', 'Ad_IAEA_2D', 'LWM'):
                bmap = [array[::-1] for array in bmap]

    except (FileNotFoundError, KeyError) as e:
        print(f"Error reading data: {e}")
        return

    # Create a single figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))  # 1 row, 2 columns

    def plot_axial_geometry(ax, time, bpos, bmap):
        widthx = [0]
        widthy = [0]
        somx = [0]
        somy = [0]

        for i in range(nx):
            widthx.append(xcm[i])
            somx.append(sum(widthx))

        for j in range(nz):
            widthy.append(ycm[j])
            somy.append(sum(widthy))

        rectangles = []
        colr = ['royalblue', 'gold', 'darkOrange', 'Red', 'Pink', 'deepskyblue', 
                'Limegreen', 'Brown', 'dodgerblue', 'SeaGreen', 'SkyBlue',
                '#ffb3e6', '#c2c2f0', '#ffb3e6', '#c2c2f0', '#ffb3e6']

        for i in range(nx):
            for j in range(nz):
                rectangles.append(mpatches.Rectangle((somx[i], somy[j]), xcm[i], ycm[j], linewidth=0.0, edgecolor='k'))

        n = nz - 1
        m = 0
        for r in rectangles:
            ax.add_artist(r)
            r.set_facecolor('white' if assm[n][m] == 0 else colr[int(assm[n][m]) - 1])

            n -= 1
            if n < 0:
                n = nz - 1
                m += 1

        reversed_somy = somy[::-1]

        # Create a mapping from values to colors
        unique_values = sorted(set(bpos))
        value_to_color = {val: colors[i % len(colors)] for i, val in enumerate(unique_values)}
        
        # Track added colors and labels
        added_colors = set()
        legend_handles = []

        for z in reversed(range(nz)):
            for i, pos in enumerate(bmap[0]):
                if pos > 0:
                    bank_type = int(pos)
                    if bank_type <= len(colors):
                        condition_value = bpos[bank_type - 1]
                        color = value_to_color.get(condition_value, 'black')  # Default to 'black' if value is not found
                        if reversed_somy[z] > condition_value:
                            ax.plot([somx[i] + xcm[i] / 2] * 2, [somy[nz - z - 1], somy[nz - z - 1] + ycm[z]], color=color, linewidth=6)
                            if color not in added_colors:
                                added_colors.add(color)
                                label = f'Position in Steps={condition_value}'
                                # Add marker to the legend
                                legend_handles.append(mlines.Line2D([], [], color=color, marker='o', markersize=10, linestyle='None', label=label))

        for i in range(nx + 1):
            ax.axvline(x=somx[i], color='k', linestyle='--', linewidth=0.5)
        for j in range(nz + 1):
            ax.axhline(y=somy[j], color='k', linestyle='--', linewidth=0.5)

        ax.set_xlabel('X [cm]')
        ax.set_ylabel('Z [cm]')
        ax.set_title(f'Axial Geometry with CR Bank Types (t={time})', loc='center')

        ax.set_ylim((min(somy), max(somy)))
        ax.set_xlim((min(somx), max(somx)))
        ax.set_xticks(somx)
        ax.set_yticks(somy)
        
        # Rotate x-ticks
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

        # Add legend outside the lower right of the plot area
        ax.legend(handles=legend_handles, loc='upper right', bbox_to_anchor=(1, 1), title='CR Bank Types')

    plot_axial_geometry(ax1, '0s', bpos, bmap)
    plot_axial_geometry(ax2, '5s', fbpos, bmap)

    # Adjust layout to make space for the legend
    plt.tight_layout(rect=[0, 0, 0.85, 1])  # Increase the right margin to fit the legend
    plt.savefig('2Dimage.png',dpi=500)

    plt.show()

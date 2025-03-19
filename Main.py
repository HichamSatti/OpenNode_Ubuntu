#! /usr/bin/python3
#! -*- coding:utf-8 -*-
import NEM
import os
import json
import numpy as np
from datetime import datetime

start  = datetime.now()

xdiv = []
ydiv = []
zdiv = []
xsize = []
ysize = []
zsize = []
zpln = []
assm = []
x_siga = []
x_sigtr = []
x_sigf = []
x_nu_sigf = []
chi = []
x_sigs = []
ds = []
sps = []
xpos = []
ypos = []
zpos = []
nk = 0
csiga = []
csigtr = []
csigf = []
cnuf = []
csigs = []
fsiga = []
fsigtr = []
fsigf = []
fnuf = []
fsigs = []
msiga = []
msigtr = []
msigf = []
mnuf = []
msigs = []
lsiga = []
lsigtr = []
lsigf = []
lnuf = []
lsigs = []
dsiga = []
dsigtr = []
dsigf = []
dnuf = []
dsigs = []
fbpos = []
tmove = []
bspeed = []
ibeta = []
lamb = []
velo = []
bmap =[]
stab= []
bpos= []


PATH    = open(os.getcwd()+'/app/link/script.dir', "r" ).read()
Methods = open(os.getcwd()+'/app/link/script00.py', "r" ).read()

if Methods == 'NEM':
    with open(PATH) as json_data:
        data = json.load(json_data)
        mode = data['Data']['Parameters']['Calculation Mode']
        ng = data['Data']['Parameters']['Number of Energy Groups']
        nmat = data['Data']['Parameters']['Number of Materials']
        order = data['Data']['Parameters']['Polynomial Order']
        nx = data['Data']['Parameters']['Number of X-Assembly']
        ny = data['Data']['Parameters']['Number of Y-Assembly']
        nz = data['Data']['Parameters']['Number of Z-Assembly']
        npp = data['Data']['Parameters']['Number of Planar']
        x_east = data['Data']['Boundary Condition']['X_East']
        x_west = data['Data']['Boundary Condition']['X_West']
        y_north = data['Data']['Boundary Condition']['Y_North']
        y_south = data['Data']['Boundary Condition']['Y_South']
        z_top = data['Data']['Boundary Condition']['Z_Top']
        z_bott = data['Data']['Boundary Condition']['Z_Bottom']
        xdiv = data['Data']['Parameters']['X-Assembly Division']
        ydiv = data['Data']['Parameters']['Y-Assembly Division']
        zdiv = data['Data']['Parameters']['Z-Assembly Division']
        xsize = data['Data']['Parameters']['X-Assembly Size']
        ysize = data['Data']['Parameters']['Y-Assembly Size']
        zsize = data['Data']['Parameters']['Z-Assembly Size']
        zpln = data['Data']['Parameters']['Planar Assignement to Z']
        assm = data['Data']['XY_Assembly']
        # assmx = data['Data']['Z_Assembly']
        for i in range(nmat):
            x_siga.append(data['Data']['Materials'][i]['Absorption_XS'])
            x_sigtr.append(data['Data']['Materials'][i]['Transport_XS'])
            x_sigf.append(data['Data']['Materials'][i]['Fission_XS'])
            x_nu_sigf.append(data['Data']['Materials'][i]['Nu*Fission_XS'])
            chi.append(data['Data']['Materials'][i]['Chi'])
            x_sigs.append(data['Data']['Materials'][i]['Scattering_XS'])


        if mode == 'Fixed Source':
            ns = data['Data']['Fixed Source Parameters']['Source Number']
            npor = data['Data']['Fixed Source Parameters']['Radial FS Number']
            npox = data['Data']['Fixed Source Parameters']['Axial FS Number']
            for i in range(ns):
                ds.append(data['Data']['Fixed Source Parameters']['Source Parameters'][i]['Source Density'])
                sps.append(data['Data']['Fixed Source Parameters']['Source Parameters'][i]['Source Energy Spectrum'])
                xpos.append(data['Data']['Fixed Source Parameters']['Source Parameters'][i]['Radial FS Position (x)'])
                ypos.append(data['Data']['Fixed Source Parameters']['Source Parameters'][i]['Radial FS Position (y)'])
                zpos.append(data['Data']['Fixed Source Parameters']['Source Parameters'][i]['Axial FS Position'])
        elif mode == 'CBC':
            nb = data['Data']['CBCSearch']['Number of CR banks']
            rbcon = data['Data']['CBCSearch']['Boron Concentration Reference']
            pos0 = data['Data']['CBCSearch']['Zero step position']
            ssize = data['Data']['CBCSearch']['step size']
            bpos = data['Data']['CBCSearch']['CR bank position']
            bmap = data['Data']['CBCSearch']['Radial CR bank map']
            nstep = data['Data']['CBCSearch']['Number of steps']
            cftem = data['Data']['CBCSearch']['average fuel temperature']
            rftem = data['Data']['CBCSearch']['fuel temperature reference']
            cmtem = data['Data']['CBCSearch']['average Moderator temperature']
            rmtem = data['Data']['CBCSearch']['Moderator temperature reference']
            ccden = data['Data']['CBCSearch']['average Coolant Density']
            rcden = data['Data']['CBCSearch']['Coolant Density reference']
            for i in range(nmat):
                csiga.append(data['Data']['CBCSearch']['B_Materials'][i]['B_Absorption_XS'])
                csigtr.append(data['Data']['CBCSearch']['B_Materials'][i]['B_Transport_XS'])
                csigf.append(data['Data']['CBCSearch']['B_Materials'][i]['B_Fission_XS'])
                cnuf.append(data['Data']['CBCSearch']['B_Materials'][i]['B_Nu*Fission_XS'])
                csigs.append(data['Data']['CBCSearch']['B_Materials'][i]['B_Scattering_XS'])
                fsiga.append(data['Data']['CBCSearch']['F_Materials'][i]['F_Absorption_XS'])
                fsigtr.append(data['Data']['CBCSearch']['F_Materials'][i]['F_Transport_XS'])
                fsigf.append(data['Data']['CBCSearch']['F_Materials'][i]['F_Fission_XS'])
                fnuf.append(data['Data']['CBCSearch']['F_Materials'][i]['F_Nu*Fission_XS'])
                fsigs.append(data['Data']['CBCSearch']['F_Materials'][i]['F_Scattering_XS'])
                msiga.append(data['Data']['CBCSearch']['M_Materials'][i]['M_Absorption_XS'])
                msigtr.append(data['Data']['CBCSearch']['M_Materials'][i]['M_Transport_XS'])
                msigf.append(data['Data']['CBCSearch']['M_Materials'][i]['M_Fission_XS'])
                mnuf.append(data['Data']['CBCSearch']['M_Materials'][i]['M_Nu*Fission_XS'])
                msigs.append(data['Data']['CBCSearch']['M_Materials'][i]['M_Scattering_XS'])
                lsiga.append(data['Data']['CBCSearch']['C_Materials'][i]['C_Absorption_XS'])
                lsigtr.append(data['Data']['CBCSearch']['C_Materials'][i]['C_Transport_XS'])
                lsigf.append(data['Data']['CBCSearch']['C_Materials'][i]['C_Fission_XS'])
                lnuf.append(data['Data']['CBCSearch']['C_Materials'][i]['C_Nu*Fission_XS'])
                lsigs.append(data['Data']['CBCSearch']['C_Materials'][i]['C_Scattering_XS'])
            #--------------------------------------------------------
                dsiga.append(data['Data']['REject']['J_Materials'][i]['J_Absorption_XS'])
                dsigtr.append(data['Data']['REject']['J_Materials'][i]['J_Transport_XS'])
                dsigf.append(data['Data']['REject']['J_Materials'][i]['J_Fission_XS'])
                dnuf.append(data['Data']['REject']['J_Materials'][i]['J_Nu*Fission_XS'])
                dsigs.append(data['Data']['REject']['J_Materials'][i]['J_Scattering_XS'])
        elif mode == 'CBCTH':
            nb = data['Data']['CBCSearch']['Number of CR banks']
            ntem = data['Data']['CBCSearch']['temperature steam table']
            pra = data['Data']['CBCSearch']['Number of Parameters']
            rbcon = data['Data']['CBCSearch']['Boron Concentration Reference']
            pos0 = data['Data']['CBCSearch']['Zero step position']
            ssize = data['Data']['CBCSearch']['step size']
            bpos = data['Data']['CBCSearch']['CR bank position']
            bmap = data['Data']['CBCSearch']['Radial CR bank map']
            nstep = data['Data']['CBCSearch']['Number of steps']
            cftem = data['Data']['CBCSearch']['average fuel temperature']
            rftem = data['Data']['CBCSearch']['fuel temperature reference']
            cmtem = data['Data']['CBCSearch']['average Moderator temperature']
            rmtem = data['Data']['CBCSearch']['Moderator temperature reference']
            ccden = data['Data']['CBCSearch']['average Coolant Density']
            rcden = data['Data']['CBCSearch']['Coolant Density reference']
            #--------------------------------------------------------
            cf = data['Data']['ThHydraulic']['heat fraction deposited into coolant']
            tin = data['Data']['ThHydraulic']['coolant inlet temperature']
            nfpin = data['Data']['ThHydraulic']['Number of fuel pin']
            cmflow = data['Data']['ThHydraulic']['Fuel Assembly mass flow rate']
            rf = data['Data']['ThHydraulic']['Fuel meat radius']
            tg = data['Data']['ThHydraulic']['gap thickness']
            tc = data['Data']['ThHydraulic']['clad thickness']
            ppitch = data['Data']['ThHydraulic']['pin picth']
            powr = data['Data']['ThHydraulic']['reactor full thermal power']
            ppow = data['Data']['ThHydraulic']['percent power']
            stab = data['Data']['ThHydraulic']['Steam table matrix']
            
            
            
            #--------------------------------------------------------
            for i in range(nmat):
                csiga.append(data['Data']['CBCSearch']['B_Materials'][i]['B_Absorption_XS'])
                csigtr.append(data['Data']['CBCSearch']['B_Materials'][i]['B_Transport_XS'])
                csigf.append(data['Data']['CBCSearch']['B_Materials'][i]['B_Fission_XS'])
                cnuf.append(data['Data']['CBCSearch']['B_Materials'][i]['B_Nu*Fission_XS'])
                csigs.append(data['Data']['CBCSearch']['B_Materials'][i]['B_Scattering_XS'])
                fsiga.append(data['Data']['CBCSearch']['F_Materials'][i]['F_Absorption_XS'])
                fsigtr.append(data['Data']['CBCSearch']['F_Materials'][i]['F_Transport_XS'])
                fsigf.append(data['Data']['CBCSearch']['F_Materials'][i]['F_Fission_XS'])
                fnuf.append(data['Data']['CBCSearch']['F_Materials'][i]['F_Nu*Fission_XS'])
                fsigs.append(data['Data']['CBCSearch']['F_Materials'][i]['F_Scattering_XS'])
                msiga.append(data['Data']['CBCSearch']['M_Materials'][i]['M_Absorption_XS'])
                msigtr.append(data['Data']['CBCSearch']['M_Materials'][i]['M_Transport_XS'])
                msigf.append(data['Data']['CBCSearch']['M_Materials'][i]['M_Fission_XS'])
                mnuf.append(data['Data']['CBCSearch']['M_Materials'][i]['M_Nu*Fission_XS'])
                msigs.append(data['Data']['CBCSearch']['M_Materials'][i]['M_Scattering_XS'])
                lsiga.append(data['Data']['CBCSearch']['C_Materials'][i]['C_Absorption_XS'])
                lsigtr.append(data['Data']['CBCSearch']['C_Materials'][i]['C_Transport_XS'])
                lsigf.append(data['Data']['CBCSearch']['C_Materials'][i]['C_Fission_XS'])
                lnuf.append(data['Data']['CBCSearch']['C_Materials'][i]['C_Nu*Fission_XS'])
                lsigs.append(data['Data']['CBCSearch']['C_Materials'][i]['C_Scattering_XS'])
            #--------------------------------------------------------
                dsiga.append(data['Data']['REject']['J_Materials'][i]['J_Absorption_XS'])
                dsigtr.append(data['Data']['REject']['J_Materials'][i]['J_Transport_XS'])
                dsigf.append(data['Data']['REject']['J_Materials'][i]['J_Fission_XS'])
                dnuf.append(data['Data']['REject']['J_Materials'][i]['J_Nu*Fission_XS'])
                dsigs.append(data['Data']['REject']['J_Materials'][i]['J_Scattering_XS'])
                
                
        elif mode == 'RODEJCT':
            nb = data['Data']['CBCSearch']['Number of CR banks']
            pos0 = data['Data']['CBCSearch']['Zero step position']
            ssize = data['Data']['CBCSearch']['step size']
            bpos = data['Data']['CBCSearch']['CR bank position']
            bmap = data['Data']['CBCSearch']['Radial CR bank map']
            nstep = data['Data']['CBCSearch']['Number of steps']
            #--------------------------------------------------------
            fbpos = data['Data']['REject']['Final CR bank position']
            tmove = data['Data']['REject']['Time when CR bank starts moving']
            bspeed = data['Data']['REject']['CR bank movement speed']
            ibeta = data['Data']['REject']['delayed neutron fraction']
            lamb = data['Data']['REject']['precusor decay constant']
            velo = data['Data']['REject']['Neutron velocity']
            ttot = data['Data']['REject']['Total Simulation Time']
            tstep1 = data['Data']['REject']['First Time Step']
            tstep2 = data['Data']['REject']['Second Time Step']
            tdiv = data['Data']['REject']['When Second Time Step Apply']
            for i in range(nmat):
                dsiga.append(data['Data']['REject']['J_Materials'][i]['J_Absorption_XS'])
                dsigtr.append(data['Data']['REject']['J_Materials'][i]['J_Transport_XS'])
                dsigf.append(data['Data']['REject']['J_Materials'][i]['J_Fission_XS'])
                dnuf.append(data['Data']['REject']['J_Materials'][i]['J_Nu*Fission_XS'])
                dsigs.append(data['Data']['REject']['J_Materials'][i]['J_Scattering_XS'])

            # Convert lists to NumPy arrays and ensure Fortran-contiguous layout
            x_sigs = np.asfortranarray(np.array(x_sigs, dtype=np.float64).reshape(nmat, ng, ng))
            x_siga = np.asfortranarray(np.array(x_siga, dtype=np.float64).reshape(nmat, ng))
            x_sigtr = np.asfortranarray(np.array(x_sigtr, dtype=np.float64).reshape(nmat, ng))
            x_sigf = np.asfortranarray(np.array(x_sigf, dtype=np.float64).reshape(nmat, ng))
            x_nu_sigf = np.asfortranarray(np.array(x_nu_sigf, dtype=np.float64).reshape(nmat, ng))
            chi = np.asfortranarray(np.array(chi, dtype=np.float64).reshape(nmat, ng))
            dsiga = np.asfortranarray(np.array(dsiga, dtype=np.float64).reshape(nmat, ng))
            dsigtr = np.asfortranarray(np.array(dsigtr, dtype=np.float64).reshape(nmat, ng))
            dsigf = np.asfortranarray(np.array(dsigf, dtype=np.float64).reshape(nmat, ng))
            dnuf = np.asfortranarray(np.array(dnuf, dtype=np.float64).reshape(nmat, ng))
            dsigs = np.asfortranarray(np.array(dsigs, dtype=np.float64).reshape(nmat, ng, ng))

            
        elif mode == 'THRODEJCT':
            nb = data['Data']['CBCSearch']['Number of CR banks']
            ntem = data['Data']['CBCSearch']['temperature steam table']
            pra = data['Data']['CBCSearch']['Number of Parameters']
            bcon = data['Data']['CBCSearch']['Boron Concentration']
            rbcon = data['Data']['CBCSearch']['Boron Concentration Reference']
            pos0 = data['Data']['CBCSearch']['Zero step position']
            ssize = data['Data']['CBCSearch']['step size']
            bpos = data['Data']['CBCSearch']['CR bank position']
            bmap = data['Data']['CBCSearch']['Radial CR bank map']
            nstep = data['Data']['CBCSearch']['Number of steps']
            cftem = data['Data']['CBCSearch']['average fuel temperature']
            rftem = data['Data']['CBCSearch']['fuel temperature reference']
            cmtem = data['Data']['CBCSearch']['average Moderator temperature']
            rmtem = data['Data']['CBCSearch']['Moderator temperature reference']
            ccden = data['Data']['CBCSearch']['average Coolant Density']
            rcden = data['Data']['CBCSearch']['Coolant Density reference']
            #--------------------------------------------------------
            fbpos = data['Data']['REject']['Final CR bank position']
            tmove = data['Data']['REject']['Time when CR bank starts moving']
            bspeed = data['Data']['REject']['CR bank movement speed']
            ibeta = data['Data']['REject']['delayed neutron fraction']
            lamb = data['Data']['REject']['precusor decay constant']
            velo = data['Data']['REject']['Neutron velocity']
            ttot = data['Data']['REject']['Total Simulation Time']
            tstep1 = data['Data']['REject']['First Time Step']
            tstep2 = data['Data']['REject']['Second Time Step']
            tdiv = data['Data']['REject']['When Second Time Step Apply']
            #--------------------------------------------------------
            cf = data['Data']['ThHydraulic']['heat fraction deposited into coolant']
            tin = data['Data']['ThHydraulic']['coolant inlet temperature']
            nfpin = data['Data']['ThHydraulic']['Number of fuel pin']
            cmflow = data['Data']['ThHydraulic']['Fuel Assembly mass flow rate']
            rf = data['Data']['ThHydraulic']['Fuel meat radius']
            tg = data['Data']['ThHydraulic']['gap thickness']
            tc = data['Data']['ThHydraulic']['clad thickness']
            ppitch = data['Data']['ThHydraulic']['pin picth']
            powr = data['Data']['ThHydraulic']['reactor full thermal power']
            ppow = data['Data']['ThHydraulic']['percent power']
            stab = data['Data']['ThHydraulic']['Steam table matrix']
            #--------------------------------------------------------
            for i in range(nmat):
                csiga.append(data['Data']['CBCSearch']['B_Materials'][i]['B_Absorption_XS'])
                csigtr.append(data['Data']['CBCSearch']['B_Materials'][i]['B_Transport_XS'])
                csigf.append(data['Data']['CBCSearch']['B_Materials'][i]['B_Fission_XS'])
                cnuf.append(data['Data']['CBCSearch']['B_Materials'][i]['B_Nu*Fission_XS'])
                csigs.append(data['Data']['CBCSearch']['B_Materials'][i]['B_Scattering_XS'])
                fsiga.append(data['Data']['CBCSearch']['F_Materials'][i]['F_Absorption_XS'])
                fsigtr.append(data['Data']['CBCSearch']['F_Materials'][i]['F_Transport_XS'])
                fsigf.append(data['Data']['CBCSearch']['F_Materials'][i]['F_Fission_XS'])
                fnuf.append(data['Data']['CBCSearch']['F_Materials'][i]['F_Nu*Fission_XS'])
                fsigs.append(data['Data']['CBCSearch']['F_Materials'][i]['F_Scattering_XS'])
                msiga.append(data['Data']['CBCSearch']['M_Materials'][i]['M_Absorption_XS'])
                msigtr.append(data['Data']['CBCSearch']['M_Materials'][i]['M_Transport_XS'])
                msigf.append(data['Data']['CBCSearch']['M_Materials'][i]['M_Fission_XS'])
                mnuf.append(data['Data']['CBCSearch']['M_Materials'][i]['M_Nu*Fission_XS'])
                msigs.append(data['Data']['CBCSearch']['M_Materials'][i]['M_Scattering_XS'])
                lsiga.append(data['Data']['CBCSearch']['C_Materials'][i]['C_Absorption_XS'])
                lsigtr.append(data['Data']['CBCSearch']['C_Materials'][i]['C_Transport_XS'])
                lsigf.append(data['Data']['CBCSearch']['C_Materials'][i]['C_Fission_XS'])
                lnuf.append(data['Data']['CBCSearch']['C_Materials'][i]['C_Nu*Fission_XS'])
                lsigs.append(data['Data']['CBCSearch']['C_Materials'][i]['C_Scattering_XS'])
                dsiga.append(data['Data']['REject']['J_Materials'][i]['J_Absorption_XS'])
                dsigtr.append(data['Data']['REject']['J_Materials'][i]['J_Transport_XS'])
                dsigf.append(data['Data']['REject']['J_Materials'][i]['J_Fission_XS'])
                dnuf.append(data['Data']['REject']['J_Materials'][i]['J_Nu*Fission_XS'])
                dsigs.append(data['Data']['REject']['J_Materials'][i]['J_Scattering_XS'])


            # Convert lists to NumPy arrays and ensure Fortran-contiguous layout
            x_sigs = np.asfortranarray(np.array(x_sigs, dtype=np.float64).reshape(nmat, ng, ng))
            x_siga = np.asfortranarray(np.array(x_siga, dtype=np.float64).reshape(nmat, ng))
            x_sigtr = np.asfortranarray(np.array(x_sigtr, dtype=np.float64).reshape(nmat, ng))
            x_sigf = np.asfortranarray(np.array(x_sigf, dtype=np.float64).reshape(nmat, ng))
            x_nu_sigf = np.asfortranarray(np.array(x_nu_sigf, dtype=np.float64).reshape(nmat, ng))
            chi = np.asfortranarray(np.array(chi, dtype=np.float64).reshape(nmat, ng))
            dsiga = np.asfortranarray(np.array(dsiga, dtype=np.float64).reshape(nmat, ng))
            dsigtr = np.asfortranarray(np.array(dsigtr, dtype=np.float64).reshape(nmat, ng))
            dsigf = np.asfortranarray(np.array(dsigf, dtype=np.float64).reshape(nmat, ng))
            dnuf = np.asfortranarray(np.array(dnuf, dtype=np.float64).reshape(nmat, ng))
            dsigs = np.asfortranarray(np.array(dsigs, dtype=np.float64).reshape(nmat, ng, ng))
            
            csiga = np.asfortranarray(np.array(csiga, dtype=np.float64))
            csigtr = np.asfortranarray(np.array(csigtr, dtype=np.float64))
            csigf = np.asfortranarray(np.array(csigf, dtype=np.float64))
            cnuf = np.asfortranarray(np.array(cnuf, dtype=np.float64))
            csigs = np.asfortranarray(np.array(csigs, dtype=np.float64))
            fsiga = np.asfortranarray(np.array(fsiga, dtype=np.float64))
            fsigtr = np.asfortranarray(np.array(fsigtr, dtype=np.float64))
            fsigf = np.asfortranarray(np.array(fsigf, dtype=np.float64))
            fnuf = np.asfortranarray(np.array(fnuf, dtype=np.float64))
            fsigs = np.asfortranarray(np.array(fsigs, dtype=np.float64))
            msiga = np.asfortranarray(np.array(msiga, dtype=np.float64))
            msigtr = np.asfortranarray(np.array(msigtr, dtype=np.float64))
            msigf = np.asfortranarray(np.array(msigf, dtype=np.float64))
            mnuf = np.asfortranarray(np.array(mnuf, dtype=np.float64))
            msigs = np.asfortranarray(np.array(msigs, dtype=np.float64))
            lsiga = np.asfortranarray(np.array(lsiga, dtype=np.float64))
            lsigtr = np.asfortranarray(np.array(lsigtr, dtype=np.float64))
            lsigf = np.asfortranarray(np.array(lsigf, dtype=np.float64))
            lnuf = np.asfortranarray(np.array(lnuf, dtype=np.float64))
            lsigs = np.asfortranarray(np.array(lsigs, dtype=np.float64))
            dsiga = np.asfortranarray(np.array(dsiga, dtype=np.float64))
            dsigtr = np.asfortranarray(np.array(dsigtr, dtype=np.float64))
            dsigf = np.asfortranarray(np.array(dsigf, dtype=np.float64))
            dnuf = np.asfortranarray(np.array(dnuf, dtype=np.float64))
            dsigs = np.asfortranarray(np.array(dsigs, dtype=np.float64))


# ns=2
# ds=[10.0,20.0]
# sps = [[1.0,0.0],[1.0,0.0]]
# xpos = [[2, 2, 1, 0],[9, 8, 8, 0]]
# ypos = [[1, 0, 0, 0],[9, 0, 0, 0]]
# zpos = [[10, 1],[10, 1]]
""" Modes """
# mode = 'Forward'
# mode = 'Adjoint'
# mode = 'Fixed Source'

""" Forward """
# nx = 9
# ny = 9
# nz = 19
# ng = 2
# nmat = 5
# np=4
# xdiv = [1]*9
# ydiv = [1]*9
# zdiv = [1]*19
# xsize=[20]*8
# xsize.insert(0,10)
# ysize=[20]*8
# ysize.insert(8,10)
# zsize = [20]*19
# zpln=[1]
# zpln.extend([2]*13)
# zpln.extend([3]*4)
# zpln.extend([4])
# assm=[[[4,4,4,4,0,0,0,0,0],[4,4,4,4,4,4,0,0,0],[4,4,4,4,4,4,4,0,0],[4,4,4,4,4,4,4,4,0],[4,4,4,4,4,4,4,4,0],[4,4,4,4,4,4,4,4,4],[4,4,4,4,4,4,4,4,4],[4,4,4,4,4,4,4,4,4],[4,4,4,4,4,4,4,4,4]],[[4,4,4,4,0,0,0,0,0],[1,1,1,4,4,4,0,0,0],[2,2,1,1,1,4,4,0,0],[2,2,2,2,1,1,4,4,0],[3,2,2,2,3,1,1,4,0],[2,2,2,2,2,2,1,4,4],[2,2,2,2,2,2,1,1,4],[2,2,2,2,2,2,2,1,4],[3,2,2,2,3,2,2,1,4]],[[4,4,4,4,0,0,0,0,0],[1,1,1,4,4,4,0,0,0],[2,2,1,1,1,4,4,0,0],[2,2,2,2,1,1,4,4,0],[3,2,2,2,3,1,1,4,0],[2,2,2,2,2,2,1,4,4],[2,2,3,2,2,2,1,1,4],[2,2,2,2,2,2,2,1,4],[3,2,2,2,3,2,2,1,4]],[[4,4,4,4,0,0,0,0,0],[4,4,4,4,4,4,0,0,0],[4,4,4,4,4,4,4,0,0],[4,4,4,4,4,4,4,4,0],[5,4,4,4,5,4,4,4,0],[4,4,4,4,4,4,4,4,4],[4,4,5,4,4,4,4,4,4],[4,4,4,4,4,4,4,4,4],[5,4,4,4,5,4,4,4,4]]]
# """# assm=[[[4,4,4,4,4,4,4,4,4],[4,4,4,4,4,4,4,4,4],[4,4,4,4,4,4,4,4,4],[4,4,4,4,4,4,4,4,4],[4,4,4,4,4,4,4,4,0],[4,4,4,4,4,4,4,4,0],[4,4,4,4,4,4,4,0,0],[4,4,4,4,4,4,0,0,0],[4,4,4,4,0,0,0,0,0]],[[3,2,2,2,3,2,2,1,4],[2,2,2,2,2,2,2,1,4],[2,2,2,2,2,2,1,1,4],[2,2,2,2,2,2,1,4,4],[3,2,2,2,3,1,1,4,0],[2,2,2,2,1,1,4,4,0],[2,2,1,1,1,4,4,0,0],[1,1,1,4,4,4,0,0,0],[4,4,4,4,0,0,0,0,0]],[[3,2,2,2,3,2,2,1,4],[2,2,2,2,2,2,2,1,4],[2,2,2,2,2,2,1,1,4],[2,2,2,2,2,2,1,4,4],[3,2,2,2,3,1,1,4,0],[2,2,2,2,1,1,4,4,0],[2,2,1,1,1,4,4,0,0],[1,1,1,4,4,4,0,0,0],[4,4,4,4,0,0,0,0,0]],[[5,4,4,4,5,4,4,4,4],[4,4,4,4,4,4,4,4,4],[4,4,5,4,4,4,4,4,4],[4,4,4,4,4,4,4,4,4],[5,4,4,4,5,4,4,4,0],[4,4,4,4,4,4,4,4,0],[4,4,4,4,4,4,4,0,0],[4,4,4,4,4,4,0,0,0],[4,4,4,4,0,0,0,0,0]]]"""
# """# assm=[[[4,4,4,4,4,4,4,4,4],[4,4,4,4,4,4,4,4,4],[4,4,4,4,4,4,4,4,4],[4,4,4,4,4,4,4,4,4],[0,4,4,4,4,4,4,4,4],[0,4,4,4,4,4,4,4,4],[0,0,4,4,4,4,4,4,4],[0,0,0,4,4,4,4,4,4],[0,0,0,0,0,4,4,4,4]],[[4,1,2,2,3,2,2,2,3],[4,1,2,2,2,2,2,2,2],[4,1,1,2,2,2,2,2,2],[4,4,1,2,2,2,2,2,2],[0,4,1,1,3,2,2,2,3],[0,4,4,1,1,2,2,2,2],[0,0,4,4,1,1,1,2,2],[0,0,0,4,4,4,1,1,1],[0,0,0,0,0,4,4,4,4]],[[4,1,2,2,3,2,2,2,3],[4,1,2,2,2,2,2,2,2],[4,1,1,2,2,2,2,2,2],[4,4,1,2,2,2,2,2,2],[0,4,1,1,3,2,2,2,3],[0,4,4,1,1,2,2,2,2],[0,0,4,4,1,1,1,2,2],[0,0,0,4,4,4,1,1,1],[0,0,0,0,0,4,4,4,4]],[[4,4,4,4,5,4,4,4,5],[4,4,4,4,4,4,4,4,4],[4,4,4,4,4,4,5,4,4],[4,4,4,4,4,4,4,4,4],[0,4,4,4,5,4,4,4,5],[0,4,4,4,4,4,4,4,4],[0,0,4,4,4,4,4,4,4],[0,0,0,4,4,4,4,4,4],[0,0,0,0,0,4,4,4,4]]]"""
# x_sigtr=[[0.222222,0.833333],[0.222222,0.833333],[0.222222,0.833333],[0.166667,1.111111],[0.166667,1.111111]]
# x_siga=[[0.010,0.080],[0.010,0.085],[0.0100,0.1300],[0.000,0.010],[0.000,0.055]]
# x_nu_sigf=[[0.000,0.135],[0.000,0.135],[0.000,0.135],[0.000,0.000],[0.000,0.000]]
# x_sigf =[[0.000,0.135],[0.000,0.135],[0.000,0.135],[0.000,0.000],[0.000,0.000]]
# chi =[[1.0,0.0],[1.0,0.0],[1.0,0.0],[0.0,0.0],[0.0,0.0]]
# x_sigs =[[[0.1922,0.020],[0.000,0.7533]],[[0.1922,0.020],[0.000,0.7483]],[[0.1922,0.020],[0.000,0.7033]],[[0.1267,0.040],[0.000,1.1011]],[[0.000,0.040],[0.000,0.000]]]
# x_east = 1
# x_west = 2
# y_north = 2
# y_south = 1
# z_top = 1
# z_bott = 1
# order = 4
#nf=6



"""Constructing wrapper function "nodxyz"..."""
nxx,nyy,nzz = NEM.nodxyz(xdiv,ydiv,zdiv,[nx,ny,nz])

"""Constructing wrapper function "asmblg_delta"..."""
node,mnum,delx,dely,delz = NEM.asmblg_delta(nmat,nxx,nyy,nzz,xsize,ysize,zsize,xdiv,ydiv,zdiv,zpln,assm,[npp,nx,ny,nz])

"""Constructing wrapper function "stagg"..."""
y_smax,y_smin,x_smax,x_smin = NEM.stagg(mnum,[nxx,nyy,nzz])

"""Constructing wrapper function "nod"..."""
nk = NEM.nod(nzz,y_smax,y_smin,[nyy])
#print('nk=',nk)

"""Constructing wrapper function "posinod"..."""
ix,iy,iz,xyz,mat = NEM.posinod(y_smax,y_smin,nk,mnum,[nxx,nyy,nzz])

"""Constructing wrapper function "deltv"..."""
delv = NEM.deltv(delx,dely,delz,ix,iy,iz,[nk,nxx,nyy,nzz])

""" Constructing wrapper function "xd_xsigr"..."""
xd,x_sigr = NEM.xd_xsigr(x_siga,x_sigs,x_sigtr,[ng,nk])
# for g in range(ng):
    # print ('g=', g, ' x_sigr=', x_sigr[0][g])

"""Constructing wrapper function "title1"..."""
NEM.title1()

"""Constructing wrapper function "timestamp"..."""
NEM.timestamp()

"""Constructing wrapper function "output"..."""
NEM.output(nx,ny,delx,dely,delz,xd,x_sigr,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,zdiv,node,zpln,x_east,x_west,y_north,y_south,z_top,z_bott,mode,[ng,np,nmat,nz,nxx,nyy,nzz])

"""Constructing wrapper function "init"..."""
# keff,jo,ji,l0,al,f0,fx1,fy1,fz1,fx2,fy2,fz2 = NEM.init(ng,nk)
                  
"""Constructing wrapper function "xs_updt"..."""
# sigs,siga,sigtr,sigf,nu_sigf,d,sigr = NEM.xs_updt(mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,[ng,nk,nmat])

"""Constructing wrapper function "d_sigr"..."""
# d,sigr = NEM.d_sigr(siga,sigs,sigtr,[ng,nk])

"""Constructing wrapper function "nodal_coup4"..."""
# p4,r4 = NEM.nodal_coup4(delx,dely,delz,d,sigr,ix,iy,iz,[ng,nk,nxx,nyy,nzz])

"""Constructing wrapper function "nodal_coup2"..."""
# p2,r2 = NEM.nodal_coup2(delx,dely,delz,d,sigr,ix,iy,iz,[ng,nk,nxx,nyy,nzz])

"""Constructing wrapper function "outer"..."""
# NEM.outer(keff,f0,fx1,fy1,fz1,fx2,fy2,fz2,order,r2,p2,r4,p4,mat,al,jo,ji,l0,d,sigr,chi,sigs,nu_sigf,ix,iy,iz,delx,dely,delz,delv,xyz,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott,[ng,nk,nmat,nxx,nyy,nzz])

# if mode == 'Forward':
"""Constructing wrapper function "forward"..."""
    # app.NEM.forward(order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,ix,iy,iz,delx,dely,delz,delv,x_smax,x_smin,xdiv,ydiv,zdiv,y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott,mode,[ng,nk,nx,ny,nz,nmat,nxx,nyy,nzz])

# if mode == 'Adjoint':
"""Constructing wrapper function "adjoint"..."""
    # app.NEM.adjoint(order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,ix,iy,iz,delx,dely,delz,delv,x_smax,x_smin,xdiv,ydiv,zdiv,y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott,mode,[ng,nk,nx,ny,nz,nmat,nxx,nyy,nzz])

# if mode == 'Fixed Source':
"""Constructing wrapper function "fixedsrc"..."""
    # app.NEM.fixedsrc(order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,ix,iy,iz,delx,dely,delz,delv,x_smax,x_smin,y_smax,y_smin,xdiv,ydiv,zdiv,xyz,x_east,x_west,y_north,y_south,z_top,z_bott,ds,sps,zpos,xpos,ypos,mode,[ng,nk,nx,ny,nz,nmat,nxx,nyy,nzz,ns,npor,npox])

#######################   Print Extra Sources   #######################
    #Constructing wrapper function "exsrc"...
if mode == 'Fixed Source':
    esrc = NEM.exsrc(nk,ds,sps,xyz,xdiv,ydiv,zdiv,xpos,ypos,zpos,[ng,nx,ny,nz,nxx,nyy,nzz,ns,npox,npor])
    NEM.fixedsrc(order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,ix,iy,iz,delx,dely,delz,delv,
             x_smax,x_smin,y_smax,y_smin,xdiv,ydiv,zdiv,xyz,x_east,x_west,y_north,y_south,z_top,z_bott,
             mode,esrc,[ng,nk,nx,ny,nz,nmat,nxx,nyy,nzz])

    #Constructing wrapper function "cbsearch"...
elif mode == 'CBC':
       NEM.cbsearch(order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,ix,iy,iz,delx,dely,delz,delv,
                   x_smax,x_smin,xdiv,ydiv,zdiv,y_smax,y_smin,xyz,x_east,x_west,y_north,
                   y_south,z_top,z_bott,mode,rbcon,csiga,csigtr,csigf,cnuf,csigs,cftem,rftem,fsiga,
                   fsigtr,fsigf,fnuf,fsigs,msiga,msigtr,msigf,mnuf,msigs,rmtem,cmtem,lsiga,lsigtr,
                   lsigf,lnuf,lsigs,rcden,ccden,dsiga,dsigtr,dsigf,dnuf,dsigs,pos0,ssize,bpos,bmap,
                   nstep,[ng,nk,nb,nx,ny,nz,nmat,nxx,nyy,nzz])

    #Constructing wrapper function "cbsearchth"...
elif mode == 'CBCTH':
       NEM.cbsearchth(order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,ix,iy,iz,delx,dely,delz,delv,
                      x_smax,x_smin,xdiv,ydiv,zdiv,xsize,ysize,y_smax,y_smin,xyz,x_east,x_west,y_north,
                      y_south,z_top,z_bott,mode,rbcon,csiga,csigtr,csigf,cnuf,csigs,cftem,rftem,fsiga,
                      fsigtr,fsigf,fnuf,fsigs,msiga,msigtr,msigf,mnuf,msigs,rmtem,cmtem,lsiga,lsigtr,
                      lsigf,lnuf,lsigs,rcden,ccden,dsiga,dsigtr,dsigf,dnuf,dsigs,pos0,ssize,bpos,bmap,
                      nstep,stab,cf,tin,nfpin,cmflow,rf,tg,tc,ppitch,powr,ppow,
                      [ng,nk,nb,nx,ny,nz,nmat,nxx,nyy,nzz,ntem,pra])
 
    #Constructing wrapper function "rod_eject"...
elif mode == 'RODEJCT':
       NEM.rod_eject(order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,ix,iy,iz,delx,dely,delz,delv,
                     x_smax,x_smin,xdiv,ydiv,zdiv,xsize,ysize,y_smax,y_smin,xyz,x_east,x_west,y_north,
                     y_south,z_top,z_bott,mode,dsiga,dsigtr,dsigf,dnuf,dsigs,pos0,ssize,bpos,bmap,
                     nstep,fbpos,tmove,bspeed,ibeta,lamb,velo,ttot,tstep1,tstep2,tdiv,
                     [ng,nk,nb,nx,ny,nz,nmat,nxx,nyy,nzz])
elif mode == 'THRODEJCT':
    #Constructing wrapper function "throd_eject"...
       NEM.throd_eject(order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,ix,iy,iz,delx,dely,delz,delv,
                       x_smax,x_smin,xdiv,ydiv,zdiv,xsize,ysize,y_smax,y_smin,xyz,x_east,x_west,y_north,
                       y_south,z_top,z_bott,mode,dsiga,dsigtr,dsigf,dnuf,pos0,ssize,bpos,bmap,nstep,fbpos,
                       tmove,bspeed,ibeta,lamb,velo,ttot,tstep1,tstep2,tdiv,stab,cf,tin,cftem,cmtem,ccden,
                       nfpin,cmflow,rf,tg,tc,ppitch,csiga,csigtr,csigf,cnuf,fsiga,fsigtr,fsigf,fnuf,msiga,
                       msigtr,msigf,mnuf,lsiga,lsigtr,lsigf,lnuf,bcon,rbcon,rftem,rmtem,rcden,csigs,fsigs,
                       msigs,lsigs,dsigs,powr,ppow,[ng,nk,nb,nx,ny,nz,nmat,nxx,nyy,nzz,ntem,pra])
else: 
    #Constructing wrapper function "modes"...
 #   app.NEM.modes(order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,ix,iy,iz,delx,dely,delz,delv,x_smax,x_smin,xdiv,ydiv,zdiv,y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott,mode,[ng,nk,nx,ny,nz,nmat,nxx,nyy,nzz])

    #Constructing wrapper function "forward"...
    NEM.forward(order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,ix,iy,iz,delx,dely,delz,delv,
                x_smax,x_smin,xdiv,ydiv,zdiv,y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,
                z_bott,mode,[ng,nk,nx,ny,nz,nmat,nxx,nyy,nzz])
"""Constructing wrapper function "outer4fx"..."""
# f0,fx1,fy1,fz1,fx2,fy2,fz2 = NEM.outerfx(order,r2,p2,r4,p4,chi,mat,d,sigr,sigs,siga,nu_sigf,ix,iy,iz,delx,dely,delz,delv,exsrc,x_smax,x_smin,y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott,[ng,nk,nmat,nxx,nyy,nzz])

"""Constructing wrapper function "outeradj"..."""
# f0,fx1,fy1,fz1,fx2,fy2,fz2 = NEM.outeradj(order,r2,p2,r4,p4,chi,mat,d,sigr,sigs,nu_sigf,ix,iy,iz,delx,dely,delz,delv,x_smax,x_smin,y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott,[ng,nk,nmat,nxx,nyy,nzz])

"""Constructing wrapper function "title2"..."""
NEM.title2()

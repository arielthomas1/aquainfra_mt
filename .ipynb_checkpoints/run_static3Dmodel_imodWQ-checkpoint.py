""" -------------------------------------------------------------------------------------------------------------
    The main script that creates and runs the 3D model.
-------------------------------------------------------------------------------------------------------------     """

#   Import python libraries
import os
import pandas as pd
#import pyevtk.vtk
import numpy as np
import math
import xarray as xr
import flopy
import flopy.utils.binaryfile as bf
import sys
#sys.path.append(r'g:/_NZ_Canterbury/Scripts')  # i needed to add this otherwise it wouldnt work somehow
sys.path.append('C:/Users/athom04/WorkPC-Master/canterbury_3d/Scripts')
#sys.path.append(r'/home/dzamrsky/_scripts_NZ')  # i needed to add this otherwise it wouldnt work somehow
import cbury_tools as cb
# import _tools_inputData as tlz_input
# import _tools_seawatModel as tlz_model
#import _tool_plotting as tlz_tp

#   the indexes are defined from the job script
#a = int((sys.argv[1]))
#time_start = int((sys.argv[2]))
#time_end = int((sys.argv[2]))

#a = 0
#time_start = 0
#time_end = 10000

#   define the model name
model_name = '_Canterbury_testModel'
#   define grid dimensions for row, col (dx and dy) and lay thickness (dz)
dx = 500
dy = 250
dz = 20
#%%
#   define the directories
# input_dir = r'/projects/0/qt16165/_dzamrsky/_NZ_Canterbury/_input'
# out_dir = r'/projects/0/qt16165/_dzamrsky/_NZ_Canterbury'
# seawat_exe_dir = r'/home/dzamrsky/swtv4'

input_dir = r'C:/Users\athom04\WorkPC-Master\canterbury_3d\Models\Final'
out_dir = r'C:/Users\athom04\WorkPC-Master\canterbury_3d\Outputs'
seawat_exe_dir = r'g:\Water_Nexus\Modelling\swt_v4_00_05\exe\swt_v4x64.exe'

main_dir = os.path.join(out_dir, model_name)
#   create the output direcotry (if it doesnt exist yet)
os.makedirs(main_dir, exist_ok=True)

#   create the IBOUND array from the Petrel output, using the cmod_full
#fields = ['I', 'J', 'K', 'X', 'Y', 'Z', 'Vsh']
#df_in = pd.read_csv(os.path.join(input_dir, 'cmod_full'), header = 8, sep = " ", index_col = False, names = fields)
#df_in['Z'] = df_in['Z'].round(1)
#df_in['Vsh'] = df_in['Vsh'].round(2)

# loading model exported from petrel (gslib format) to dataframe
fields=['I','J','K','X','Y','Z','Vsh','Code','Por'] #now including compacted porosity
df_in=pd.read_csv(os.path.join(input_dir, 'cmod_full'), header=10, sep=" ",index_col=False,names=fields)

# clean-up/rounding excess decimal points
df_in['Z']=df_in['Z'].round(1)
df_in['Vsh']=df_in['Vsh'].round(2)
df_in['Code']=round(df_in['Code'])
df_in['Por']=df_in['Por'].round(2)



#%%
#   make a 3D numpy array with the dimensions from the input file
"""----------   Checl that the IJK matches the lay row column!  """
lay_arr = np.unique(df_in.K.values)  # unique layer numbers
row_arr = np.unique(df_in.I.values)  # unique row numbers
col_arr = np.unique(df_in.J.values)  # unique column numbers
full_ibound_arr = np.zeros((lay_arr.shape[0], row_arr.shape[0], col_arr.shape[0]))
#%%
#   get unique z values (i guess those are the top of the layer right)
z_vals = np.unique(df_in.Z.values)
"""----------   When you check the z_vals there is a -99 value i guess thats an error right?  """
top_elev = np.nanmax(z_vals)  #  the max value will be the top elevation of the model
bot_elev = np.arange(top_elev, np.nanmin(z_vals) - 2 * dz, -dz)[1:]  # have to do 2 * dz so it keeps the right end value

#   define the model stress periods, names etc.
sp_names = ['cmod_1', 'cmod_2', 'cmod_3', 'cmod_4', 'cmod_full']  # just set it to match the csv names
sp_time_dur = [1e05, 1e05, 1e05, 1e05]  # duriation of each stress period - in years
sp_sea_level = [-130., -80., -20., 0.]  # sea level for each stress period

#   create the model directory for the stress period model
model_name_sp = sp_names[a]
sp_dir = os.path.join(main_dir, model_name_sp)
os.makedirs(sp_dir, exist_ok=True)
sea_level = sp_sea_level[a]  # get the sea level for the stress period

"""--------- TODO : check if part of the simulation already ran - so if output .nc file already exists 
                    create the full_ibound_arr, top_arr etc. only once and save as .npy and load it
"""

#   read in the csv file with active geology for the stress period
df_sp_in = pd.read_csv(os.path.join(input_dir, model_name_sp), header = 8, sep = " ", index_col = False, names = fields)
#   add the lithology codes to the dataframe
df = cb.get_lith_v2(df_sp_in, 0.3, 0.6)
df['hyd_k'] = pd.NaT
df['por'] = pd.NaT
df.loc[df.Code == 1, ['hyd_k', 'por']] = [float(0.1), float(0.3)]
df.loc[df.Code == 2, ['hyd_k', 'por']] = [float(0.01), float(0.45)]
df.loc[df.Code == 3, ['hyd_k', 'por']] = [float(0.0001), float(0.60)]

#   now we can use this dataframe to fill in the ibound, poro and hk arrays as copy of the full ibound array
ibound_arr = np.copy(full_ibound_arr)
hk_arr = np.copy(full_ibound_arr)
poro_arr = np.copy(full_ibound_arr)
#   now gor row by row and if the hk_val is not nan or 0 then adapt the ibound (and other) array
for row in df.iterrows():
    if row[1]['Vsh'] > 0.:
        hk_arr[row[1]['K'] - 1, row[1]['I'] - 1, row[1]['J'] - 1] = row[1]['hyd_k']
        poro_arr[row[1]['K'] - 1, row[1]['I'] - 1, row[1]['J'] - 1] = row[1]['por']
        ibound_arr[row[1]['K'] - 1, row[1]['I'] - 1, row[1]['J'] - 1] = 1

"""

    ----- Fix the vtk file export later, for some reason it only exports the outside cells

#   create a vtk file for the HK array
z_arr = np.insert(bot_elev, 0, top_elev)
z_arr = np.repeat(z_arr[:, np.newaxis], row_arr.shape[0], axis=1)
z_arr = np.repeat(z_arr[:, :, np.newaxis], col_arr.shape[0], axis=2)
z_arr = np.insert(z_arr, -1, np.array([z_arr[:, -1, :]]), axis=1)
z_arr = np.insert(z_arr, -1, np.array([z_arr[:, :, -1]]), axis=2)

#   create the geology plots and vtk files
main_res_dir = os.path.join(out_dir, model_name + 'results')
os.makedirs(main_res_dir, exist_ok = True)
vtk_dir = os.path.join(main_res_dir, '_vtk')
os.makedirs(vtk_dir, exist_ok=True)
x_arr = np.arange(0, hk_arr.shape[2] * dx + dx, dx)
y_arr = np.arange(0, hk_arr.shape[1] * dy + dy, dy)
#   change the 0 values in the hk_arr to -99
hk_arr[hk_arr == 0] = -99
tlz_tp.create_VtkFile(z_arr, x_arr, y_arr, hk_arr, os.path.join(out_dir, 'hk_arr')

plot_hk_dir = os.path.join(out_dir, '_plots_HK')
os.makedirs(plot_hk_dir, exist_ok=True)
for lay_hk in range(hk_arr.shape[0]):
    tlz_tp.plot_2Darray_HK(hk_arr[lay_hk, :, :], 'Hk for layer ' + str(lay_hk),
                           os.path.join(plot_hk_dir, 'Hk_lay_' + str(lay_hk) + '.png'))

"""

"""--------- TODO : Change the values that go into the new cells based on new hk_arr
"""

#   define starting concentrations and starting heads (those will be changed later on when we start the model run)
if time_start == 0 and a == 0:
    sconc_arr = full_ibound_arr * 35.
    strt_arr = full_ibound_arr * 0
#   if the time_start is 0 and a is non zero then read the sconc and strt arrays from the final mini stress period
#   of the previous stress period
elif time_start == 0 and a != 0:
    sp_old = sp_names[a - 1]
    mini_modelname_old = 'SP_00' + str(int(sp_time_dur[a] - 10000)) + '_to_0' + str(int(sp_time_dur[a]))
    temp_dir_old = os.path.join(main_dir, sp_old, mini_modelname_old)
    #   load the sconc and strt arrays, if there are any nan values set them to 0
    sconc_dir = os.path.join(temp_dir_old, '_sconc_arr.npy')
    sconc_arr = np.load(sconc_dir, allow_pickle=True)
    #   set all the nan values to 0, all new cells should then also be 0. automatically
    sconc_arr[np.isnan(sconc_arr)] = 0
    strt_dir = os.path.join(temp_dir_old, '_strt_arr.npy')
    strt_arr = np.load(strt_dir, allow_pickle=True)
    #   set all the nan values to 0, all new cells should then also be 0. automatically
    strt_arr[np.isnan(strt_arr)] = 0
#   otherwise go to the previous mini stress period of the main stress period
else:
    sp_old = sp_names[a]
    #   make sure the strings are correct
    if time_start - 10000 == 0:
        time_start_str_old = '0000000'
    else:
        time_start_str_old = '00' + str(int(time_start - 10000))
    #   also the end string
    time_end_str_old = '00' + str(time_start)
    #   load the sconc and strt arrays, if there are any nan values set them to 0
    mini_modelname_old = 'SP_' + time_start_str_old + '_to_' + time_end_str_old
    temp_dir_old = os.path.join(main_dir, sp_old, mini_modelname_old)
    #   load the sconc and strt arrays, if there are any nan values set them to 0
    sconc_dir = os.path.join(temp_dir_old, '_sconc_arr.npy')
    sconc_arr = np.load(sconc_dir, allow_pickle=True)
    sconc_arr[np.isnan(sconc_arr)] = 0
    strt_dir = os.path.join(temp_dir_old, '_strt_arr.npy')
    strt_arr = np.load(strt_dir, allow_pickle=True)
    strt_arr[np.isnan(strt_arr)] = 0

#   define the name of the mini stress period
if time_start == 0:
    time_start_str = '0000000'
elif 10000 <= time_start < 100000:
    time_start_str = '00' + str(time_start)
elif 100000 <= time_start < 1000000:
    time_start_str = '0' + str(time_start)
#   do the same for the end time
time_end = time_start + 10000
if 10000 <= time_end < 100000:
    time_end_str = '00' + str(time_end)
elif 100000 <= time_end < 1000000:
    time_end_str = '0' + str(time_end)

#   define the modelname and folder for the mini stress period
mini_modelname = 'SP_' + time_start_str + '_to_' + time_end_str
mini_sp_dir = os.path.join(sp_dir, mini_modelname)

#   create the SEAWAT model object and start creating individual packages
mswt = flopy.seawat.Seawat(mini_modelname, 'nam_swt', model_ws = mini_sp_dir, exe_name = seawat_exe_dir)

#   make the DIS package, define the input first
nlay = lay_arr.shape[0]
nrow = row_arr.shape[0]
ncol = col_arr.shape[0]
delr = [dx] * ncol
delc = [dy] * nrow
perlen = [365.25 * 9999]
nstp = [10]
nper = len(perlen)
dis = flopy.modflow.ModflowDis(mswt, nlay, nrow, ncol, nper = 1, delr = delr, delc = delc, top = top_elev,
                               botm = bot_elev, perlen = perlen, nstp = nstp)

#   next create the BAS package
bas = flopy.modflow.ModflowBas(mswt, ibound = ibound_arr, strt = strt_arr)

#   in this setup we will use the BCF package
vk_arr = hk_arr * 0.1
lpf = flopy.modflow.ModflowLpf(mswt, laytyp = 0, hk = hk_arr, vka = vk_arr, ipakcb = 1)

#   create the icbund array
icbund_arr = ibound_arr

#   create the GHB package input here, also start creating the SSM package
itype = flopy.mt3d.Mt3dSsm.itype_dict()
ghb_input_lst = []
chb_input_lst = []
ssmdata = []
cond_const = 1000000.
#   the inland part on the edges of the active model domain will be assigned the topographical head
#   for each row check the first active cell - and if it is above sea level then assign fresh head
col_idx = [0, ncol - 1]
for i in range(nrow):
    #   select the active cells only
    for col in col_idx:
        row_cells = [t for t in ibound_arr[:, i, col].tolist() if t == 1]
        if len(row_cells) > 0:
            lay_idx = ibound_arr[:, i, col].tolist().index(1)
            if ibound_arr[lay_idx, i, col] == 1 and bot_elev[lay_idx] + dz >= sea_level:
                #   check if the top elevation of the cell is above sea level (thats why we do bot_elev + dz)
                for k in range(nlay):
                    if ibound_arr[k, i, col] == 1:
                        cond_val = hk_arr[k, i, col] * dz * 1000
                        ghb_input_lst.append([k, i, col, bot_elev[k] + dz, cond_val])
                        ssmdata.append([k, i, 0, 0.0, itype['GHB']])
                        icbund_arr[k, i, col] = -1
#   now do the same but for columns
row_idx = [0, nrow - 1]
for i in range(ncol):
    #   select the active cells only
    for row in row_idx:
        col_cells = [t for t in ibound_arr[:, row, i].tolist() if t == 1]
        if len(col_cells) > 0:
            lay_idx = ibound_arr[:, row, i].tolist().index(1)
            if ibound_arr[lay_idx, row, i] == 1 and bot_elev[lay_idx] + dz >= sea_level:
                for k in range(nlay):
                    if bot_elev[k] + dz >= sea_level:
                        cond_val = hk_arr[k, row_idx, i] * dz * 1000
                        ghb_input_lst.append([k, row_idx, i, bot_elev[k] + dz, cond_val])
                        ssmdata.append([k, row_idx, i, 0.0, itype['GHB']])
                        icbund_arr[k, row_idx, i] = -1
#   now check for the offshore domain and set all the cells below sea level to saltwater concentration and
#   head equal to sea level. Only in the top layer
for i in range(nrow):
    for j in range(ncol):
        lay_cells = [t for t in ibound_arr[:, i, j].tolist() if t == 1]
        if len(lay_cells) > 0:
            lay_idx = ibound_arr[:, i, j].tolist().index(1)
            if ibound_arr[lay_idx, i, j] == 1 and bot_elev[lay_idx] + dz < sea_level:
                cond_val = (vk_arr[lay_idx, i, j] * dx * dy) / dz
                ghb_input_lst.append([lay_idx, i, j, sea_level, cond_val])
                ssmdata.append([lay_idx, i, j, 35.0, itype['GHB']])
#   write the final output dictionary, inlcude each stress period
ghb_arr_in = {}
for d in range(len(perlen)):
    ghb_arr_in[d] = ghb_input_lst
ghb = flopy.modflow.ModflowGhb(mswt, ipakcb = 1, stress_period_data = ghb_arr_in)

#   the RCH package
rch_val = 0.00025
rch_arr = np.zeros((1, ibound_arr.shape[1], ibound_arr.shape[2]))
#   only apply recharge to the cells above sea level, for each column find the first active layer
for i in range(ibound_arr.shape[1]):
    for j in range(ibound_arr.shape[2]):
        lay_lst = [t for t in ibound_arr[:, i, j].tolist() if t == 1]
        if len(lay_lst) > 0:
            top_act_lay = ibound_arr[:, i, j].tolist().index(1)
            top_act_elev = top_elev - top_act_lay * dz
            #   if the top elevation of the cell is above sea level then assign recharge to it
            if top_act_elev >= sea_level:
                rch_arr[top_act_lay, i, j] = rch_val
rch = flopy.modflow.ModflowRch(mswt, nrchop = 3, ipakcb = 1, rech = rch_arr)

#   the DRN package is assigned only to cells that receive recharge - cells with elev above sea level
drn_input_lst = []
for i in range(ibound_arr.shape[1]):
    for j in range(ibound_arr.shape[2]):
        lay_lst = [t for t in ibound_arr[:, i, j].tolist() if t == 1]
        if len(lay_lst) > 0:
            top_act_lay = ibound_arr[:, i, j].tolist().index(1)
            top_act_elev = top_elev - top_act_lay * dz
            #   if the top elevation of the cell is above sea level then assign recharge to it
            if top_act_elev >= sea_level:
                cond_cell = (vk_arr[top_act_lay, i, j] * dx * dy) / dz
                drn_input_lst.append([int(top_act_lay), i, j, top_act_elev, cond_cell])
#   write the final output dictionary, inlcude each stress period
if len(drn_input_lst) > 0:
    drn_arr_in = {}
    for c in range(len(perlen)):
        drn_arr_in[c] = drn_input_lst
    drn = flopy.modflow.ModflowDrn(mswt, ipakcb=1, stress_period_data=drn_arr_in)

#   write the OC package
ihedfm = 1  # a code for the format in which heads will be printed.
iddnfm = 0  # a code for the format in which drawdowns will be printed.
extension = ['oc', 'hds', 'ddn', 'cbc']
unitnumber = [14, 30, 52, 51]
#   create the dictionary that defines how to write the output file
spd = {(0, 0): ['SAVE HEAD', 'SAVE BUDGET', 'PRINT HEAD', 'PRINT BUDGET', 'SAVE HEADTEC', 'SAVE CONCTEC',
                'SAVE VXTEC', 'SAVE VYTEC', 'SAVE VZTEC']}
for t in range(0, nper):
    per = t  # + 1
    #   to save space on disk, every 10th timestep is saved
    for g in range(1, nstp[t] + 1):
        spd[(per, int(g))] = ['SAVE HEAD', 'SAVE BUDGET', 'PRINT HEAD', 'PRINT BUDGET', 'SAVE HEADTEC', 'SAVE CONCTEC',
                              'SAVE VXTEC', 'SAVE VYTEC', 'SAVE VZTEC']
oc = flopy.modflow.ModflowOc(mswt, ihedfm=ihedfm, stress_period_data=spd, unitnumber=unitnumber, compact=True)

#   the BTN package
porosity = 0.3
dt0 = 365.25
nprs = 1
ifmtcn = 0
chkmas = False
nprmas = 10
nprobs = 10
timprs_lst = list(np.linspace(1., perlen[0], 10, endpoint=True, dtype=int))
btn = flopy.mt3d.Mt3dBtn(mswt, nprs=nprs, timprs=timprs_lst, prsity=porosity, sconc=sconc_arr,
                         ifmtcn=ifmtcn, chkmas=chkmas, nprobs=nprobs, nprmas=nprmas, dt0=0)

#   write the ADV package
adv = flopy.mt3d.Mt3dAdv(mswt, mixelm=0, mxpart=2000000)

#   write the ADV package
dmcoef = 0.0000864  # effective molecular diffusion coefficient [M2/D]
al = 1.
trpt = 0.1
trpv = 0.1
dsp = flopy.mt3d.Mt3dDsp(mswt, al=al, trpt=trpt, trpv=trpv, dmcoef=dmcoef)

#   write the VDF package
iwtable = 0
densemin = 1000.
densemax = 1025.
denseref = 1000.
denseslp = 0.7143
firstdt = 0.001
vdf = flopy.seawat.SeawatVdf(mswt, iwtable=iwtable, densemin=densemin, densemax=densemax,
                             denseref=denseref, denseslp=denseslp, firstdt=firstdt)

#   write the SSM package
ssm_rch_in = np.copy(rch_arr) * 0.0
ssmdata_dict = {0: ssmdata, 1: ssmdata}
ssm = flopy.mt3d.Mt3dSsm(mswt, crch=ssm_rch_in, stress_period_data=ssmdata_dict)

#   write packages and run model
mswt.write_input()

#   write the ascii file with vertical sum of active cells in IBOUND
ibound_arr_sum = np.sum(ibound_arr, axis=0, dtype=np.int32)
ibound_arr_sum = ibound_arr_sum.astype(str)
with open(os.path.join(mini_sp_dir, 'LOAD.ASC'), 'wb') as f:
    f.write(ibound_arr_sum)

#   create the pksf and pkst files - change it in case the grid discretization changes
pksf_lines = ['ISOLVER 1', 'NPC 2', 'MXITER 200', 'RELAX .98', 'HCLOSEPKS 0.01', 'RCLOSEPKS 10000.0', 'PARTOPT 0',
              'PARTDATA', 'external 40 1. (free) -1', 'GNCOL 274', 'GNROW 223', 'GDELR', '1000', 'GDELC', '1000',
              'NOVLAPADV 2', 'END']
pkst_lines = ['ISOLVER 2', 'NPC 2', 'MXITER 1000', 'INNERIT 50', 'RELAX .98', 'RCLOSEPKS 1.0E-05',
              'HCLOSEPKS 1.0E+12', 'RELATIVE-L2NORM', 'END']
# 'CCLOSEPKS=0.00001'

with open(os.path.join(mini_sp_dir, mini_modelname + '.pksf'), 'w') as f:
    for line in pksf_lines:
        f.write(line)
        f.write('\n')

with open(os.path.join(mini_sp_dir, mini_modelname + '.pkst'), 'w') as f:
    for line in pkst_lines:
        f.write(line)
        f.write('\n')

#   open the nam_swt file and append these three lines
nam_lines = ['PKSF          	  27 ' + mini_modelname + '.pksf', 'PKST              35 ' +
             mini_modelname + '.pkst', 'DATA 40 LOAD.ASC']

with open(os.path.join(mini_sp_dir, mini_modelname + '.nam_swt'), 'a') as f:
    for line in nam_lines:
        f.write(line)
        f.write('\n')

#   save the ibound_arr (if it doesnt exist)
ibound_arr_dir = os.path.join(sp_dir, 'ibound_arr.npy')
if not os.path.isfile(ibound_arr_dir):
    np.save(ibound_arr_dir, ibound_arr)

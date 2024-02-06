# %% CONFIG
from mat4py import loadmat
from parflow import Run, ParflowBinaryReader

# Set your project path here :
proj_path = "/home/rappl/PROJETS/Wank1D_ex/"

scripts_path = proj_path + "scripts/"
exec(open(scripts_path+'fonctions_comp.py').read())

# %% OBS

obs_path = proj_path + "obs/"            
mat_obs = loadmat(obs_path+'Eval_mil_2005_2012.mat')
obs = read_obs(mat_obs)

# %% SIMU SETTINGS

pfidb = Run.from_definition(proj_path + 'tcl_pfidb/'+'Wank1D.pfidb')
simus_path = proj_path + 'simus/'

#Set your simulation directories
simu_dirs = ['LRH/', 'MT/']

# %% LOAD DATA

#simu1
with ParflowBinaryReader(glob.glob(simus_path+simu_dirs[0]+'*.press.00000.pfb')[0]) as s: h = s.header
simu1 = SIMU(simus_path+simu_dirs[0], simu_dirs[0][:-1], pfidb, h,'snow')

#simu2
with ParflowBinaryReader(glob.glob(simus_path+simu_dirs[1]+'*.press.00000.pfb')[0]) as s: h = s.header
simu2 = SIMU(simus_path+simu_dirs[1], simu_dirs[0][:-1], pfidb, h)

# %% NETCDF 
simu1.clm.to_netcdf(simus_path+'clm1.nc')
simu1.pf.to_netcdf(simus_path+'pf1.nc')

simu2.clm.to_netcdf(simus_path+'clm2.nc')
simu2.pf.to_netcdf(simus_path+'pf2.nc')



"""
Created on 2024_02_05
@author: rappl
"""

import numpy as np
import glob
import xarray as xr
import pandas as pd
#import datetime as dt
from parflow import read_pfb, read_pfb_sequence
import datetime

def process_alb (ds):
    Alb= ds.surfalb
    Time = ds.time.to_series()
    for i in range(len(Time)) :
        if (Time[i].hour < 10) or (Time[i].hour>14):
            Alb[i] = float('nan')
    ds['Alb_process'] = Alb
    return(ds)

######################## Lecture OBSERVATIONS ##############################

def attrs_unit_obs(ds) :
    for var in ds.data_vars :
        if ('flux' in var) or ('W' in var) or ('Rn' in var) or ('LE' in var) or ('H' in var) :
            ds[var].attrs['units']= 'W/m^2'
        elif 'TS' in var :
            ds[var].attrs['units']= '°C'
        else :
            ds[var].attrs['units']= '-'
    return ds

def pre_process(ds):
    for var in ds.data_vars :
        ds[var] = xr.where(ds[var]>-1000, ds[var], float('nan'))
        if ('alb' in var) or ('SMC' in var) : 
            ds[var] = xr.where(ds[var]<=1, ds[var], float('nan'))
            ds[var] = xr.where(ds[var]>=0, ds[var], float('nan'))
    return(ds)


def dates_mat(submat,pos_DOY) :
    start = submat[0]
    end = submat[-1]
    date_start = datetime.datetime(int(start[0]),1,1) + datetime.timedelta(days=(start[pos_DOY]-1))
    date_start = date_start.replace(microsecond=0)
    date_end = datetime.datetime(int(end[0]),1,1) + datetime.timedelta(days = (end[pos_DOY]-1))
    date_end = date_end.replace(microsecond=0)
    return(date_start, date_end)


def read_obs (mat) :
    date_start, date_end  = dates_mat(mat['FluxNRJ_mil'],4)
    mat['FluxEC_mil']=np.array(mat['FluxEC_mil'])
    mat['FluxNRJ_mil']=np.array(mat['FluxNRJ_mil'])
    mat['SMC_mil']=np.array(mat['SMC_mil'])
    mat['Tempsol_mil']=np.array(mat['Tempsol_mil'])
    ds = xr.Dataset({ "H": (("time"), mat['FluxEC_mil'][1823:,-2]) ,
                                "LE": (("time"), mat['FluxEC_mil'][1823:,-1]) ,
                                "Rn": (("time"), mat['FluxNRJ_mil'][:,5]) ,
                                "surfalb": (("time"), mat['FluxNRJ_mil'][:,6]) ,
                                "fluxG1": (("time"), mat['FluxNRJ_mil'][:,7]) ,
                                "fluxG2": (("time"), mat['FluxNRJ_mil'][:,8]) ,
                                "fluxG3": (("time"), mat['FluxNRJ_mil'][:,9]) ,
                                "SWin": (("time"), mat['FluxNRJ_mil'][:,10]) ,
                                "SWout": (("time"), mat['FluxNRJ_mil'][:,11]) ,
                                "LWin": (("time"), mat['FluxNRJ_mil'][:,12]) ,
                                "LWout": (("time"), mat['FluxNRJ_mil'][:,13]) ,
                                "SMC10": (("time"), mat['SMC_mil'][:,5]) ,
                                "SMC50": (("time"), mat['SMC_mil'][:,6]) ,
                                "SMC100": (("time"), mat['SMC_mil'][:,7]) ,
                                "SMC150": (("time"), mat['SMC_mil'][:,8]) ,
                                "SMC200": (("time"), mat['SMC_mil'][:,9]) ,
                                "SMC250": (("time"), mat['SMC_mil'][:,10]) ,
                                "TS10": (("time"), mat['Tempsol_mil'][:,5]) ,
                                "TS50": (("time"), mat['Tempsol_mil'][:,6]) ,
                                "TS100": (("time"), mat['Tempsol_mil'][:,7]) ,
                                "TS150": (("time"), mat['Tempsol_mil'][:,8]) #,
                                #"TS200": (("time"), mat['Tempsol_mil'][:,9]) ,
                                #"TS250": (("time"), mat['Tempsol_mil'][:,10]) 
                                } 
                                ,
            coords={"time": pd.date_range(date_start, periods=len(mat['FluxNRJ_mil']), end = date_end),
                "reference_time": pd.Timestamp(date_start)})
    ds = pre_process(ds)
    ds = ds.assign(evap_tot = lambda x: x.LE*1800/(2.504*10**6))
    ds = attrs_unit_obs(ds)
    ds = process_alb(ds) 
    return(ds)


########################### Lecture SIMULATION ##############################


############ Mini fonctions clm

def calculate_OVFLOW(width,manning,slope,h):
    #simple function which return the overland flux from one cell to another across one cell
    return ((width/manning)*(slope**0.5)*h**(5./3.))

def calculate_runoff_from_press(htop,slopex,slopey,dx,dy,m):
    #compute runoff at single location, htop dim is time
    runoff=[]
    Sy = np.abs(slopey)
    Sx = np.abs(slopex)
    for h in htop :
        if h > 0 :
            runoff.append(np.abs(calculate_OVFLOW(dy,m,Sx,h))+np.abs(calculate_OVFLOW(dx,m,Sy,h)))
        else:
            runoff.append(0.)
    return (np.asarray(runoff))


############### Mini fonctions forçages 
def process_LAI(wdir) :
    lai = []
    with open(wdir+'lai.dat', "r") as f:
        for line in f:
            data = line.split()
            lai.append( float(data[-1]) )
    return(lai)

def process_SAI(wdir) :
    sai = []
    with open(wdir+'sai.dat', "r") as f:
        for line in f:
            data = line.split()
            sai.append( float(data[-1]) )
    return(sai[:105168])

################# CLASS SIMU

class SIMU :
    "définition d'une simulation"
    def __init__(self, simu_path, name, pfidb, h, src ='', date_start="2006-01-01"):
        self.path = simu_path
        self.name = name
        self.pfidb = pfidb
        self.src = src
        self.manning = pfidb.Mannings.Geom.domain.Value
        self.h = h
        self.date_start = date_start
        self.get_z()
        self.read_forcings()
        self.read_clm_outputs()
        self.date_end = str(self.clm.time[-1])[36:46]
        self.read_and_process_pf_outputs()

    def get_z(self):
        self.var_dz = read_pfb(glob.glob(self.path+'*mult*.pfb')[0])
        var_dz_vec = self.var_dz[:,0,0]
        var_dz_vec = var_dz_vec[::-1]
        dz = var_dz_vec * self.h['dz']
        self.z = np.cumsum(dz) - var_dz_vec/2 
        return()

############################################# Lecture CLM

    def attrs_unit_clm(self) :
        ds = self.clm
        for var in ds.data_vars :
            if ('lh' in var) or ('lw' in var) or ('sh' in var) or ('soil' in var) :
                ds[var].attrs['units']='W/m^2'
            elif ('evap'in var) or ('tran' in var) or ('infl' in var):
                ds[var].attrs['units']='mm/h'
            elif ('TS' in var) or ('t_' in var) :
                ds[var] = ds[var]-272.15
                ds[var].attrs['units']='°C'
            else :
                ds[var].attrs['units']='-'
        ds.swe_out.attrs['units']='mm'
        ds.surfalb.attrs['units']='-'
        ds.htop.attrs['units']='m'
        ds.y.attrs['units']='m'
        ds.Q.attrs['units']='m^3/CLMtimestep'
        self.clm = ds
        return

    def read_clm_outputs(self):
        h = self.h
        files = np.sort(glob.glob(self.path+'*.clm_output.*.pfb'))
        clms = read_pfb_sequence(files)
        # create dataset
        if self.src=='snow' :
            ds = xr.Dataset({"lh_tot": (("time","y"), clms[:,0,:,0]), # latent heat flux
                            "lwrad_out": (("time","y"), clms[:,1,:,0]), 
                            "sh_tot": (("time","y"), clms[:,2,:,0]), #sensible heat flux
                            "soil_grnd": (("time","y"), clms[:,3,:,0]),
                            "evap_tot": (("time","y"), clms[:,4,:,0]*60*60),
                            "evap_grnd": (("time","y"), clms[:,5,:,0]*60*60),
                            "evap_soi": (("time","y"), clms[:,6,:,0]*60*60),
                            "evap_veg": (("time","y"), clms[:,7,:,0]*60*60),
                            "tran_veg": (("time","y"), clms[:,8,:,0]*60*60),
                            "infl": (("time","y"), clms[:,9,:,0]*60*60),
                            "swe_out": (("time","y"), clms[:,10,:,0]),
                            "surfalb":(("time","y"), clms[:,13,:,0]),
                            #"ndvi":(("time","y"), clms[:,15,:,0]),
                            "t_grnd": (("time","y"), clms[:,17,:,0]), 
                            "htop": (("time","y"), clms[:,18,:,0]),
                            #"TS10": (("time","y"), clms[:,37,:,0]),
                            #"TS50": (("time","y"), clms[:,30,:,0]),
                            #"TS100": (("time","y"), clms[:,24,:,0]),
                            #"TS150": (("time","y"), clms[:,20,:,0]),
                            },
                coords={
                    "y":np.arange(start = h['y'],stop = h['y']+ h['ny']*h['dy'],step=h['dy']),
                    "z":-self.z,
                    "time": pd.date_range(self.date_start, periods=len(files),freq='1H'),
                    "reference_time": pd.Timestamp(self.date_start)})
        else : 
            print('MTreading')
            ds = xr.Dataset({"lh_tot": (("time","y"), clms[:,0,:,0]), # latent heat flux
                            "lwrad_out": (("time","y"), clms[:,1,:,0]), 
                            "sh_tot": (("time","y"), clms[:,2,:,0]), #sensible heat flux
                            "soil_grnd": (("time","y"), clms[:,3,:,0]),
                            "evap_tot": (("time","y"), clms[:,4,:,0]*60*60),
                            "evap_grnd": (("time","y"), clms[:,5,:,0]*60*60),
                            "evap_soi": (("time","y"), clms[:,6,:,0]*60*60),
                            "evap_veg": (("time","y"), clms[:,7,:,0]*60*60),
                            "tran_veg": (("time","y"), clms[:,8,:,0]*60*60),
                            "infl": (("time","y"), clms[:,9,:,0]*60*60),
                            "swe_out": (("time","y"), clms[:,10,:,0]),
                            "surfalb":(("time","y"), clms[:,11,:,0]),
                            "t_grnd": (("time","y"), clms[:,12,:,0]), #13
                            "htop": (("time","y"), clms[:,13,:,0]),#14
                            #"TS10": (("time","y"), clms[:,37,:,0]),
                            #"TS50": (("time","y"), clms[:,30,:,0]),
                            #"TS100": (("time","y"), clms[:,24,:,0]),
                            #"TS150": (("time","y"), clms[:,20,:,0]),
                            },
                coords={
                    "y":np.arange(start = h['y'],stop = h['y']+ h['ny']*h['dy'],step=h['dy']),
                    "z":-self.z,
                    "time": pd.date_range(self.date_start, periods=len(files),freq='1H'),
                    "reference_time": pd.Timestamp(self.date_start)})
        ds = ds.assign(slopex=(("y"),read_pfb(glob.glob(self.path+'*slope_x*.pfb')[0])[0,:,0]))
        ds = ds.assign(slopey=(("y"),read_pfb(glob.glob(self.path+'*slope_y*.pfb')[0])[0,:,0]))
        ds = ds.assign(mask=(("z","y","x"),read_pfb(glob.glob(self.path+'*mask*.pfb')[0])))
        ds = ds.assign(Q=lambda x: ('time',calculate_runoff_from_press(x.htop.data[:,0],
                                                                 x.slopex.data[0],
                                                                 x.slopey.data[0],
                                                                 10,10,self.manning)))
        ds = ds.assign(Rn = lambda x: x.lh_tot+ x.sh_tot +x.soil_grnd)
        ds = ds.assign(swrad_out = lambda x: x.Rn - x.lwrad_out)
        ds = pre_process(ds)
        ds = process_alb(ds)
        self.clm=ds
        self.attrs_unit_clm()
        return 

############################################# Lecture FORCINGS

    def attrs_unit_forc(self):
        ds = self.forc
        for var in ds.data_vars :
            if ('SW'in var) or ('LW' in var) or ('H' in var) :
                ds[var].attrs['units'] = 'W/m^2'
            elif 'Wind' in var :
                ds[var].attrs['units'] = 'm/s'
            else :
                ds[var].attrs['units'] = '-'
        ds.P.attrs['units'] = 'mm'
        ds.T.attrs['units'] = '°C'
        ds.Press.attrs['units'] = 'm'
        self.forc = ds
        return 

    def read_forcings(self) :
        with open(self.path+'forcagePF.txt.0') as f:
            forc = np.array([[float(s) for s in ln.split()] for ln in f])
        forc = np.array(forc)
        LAI = process_LAI(self.path)
        SAI = process_SAI(self.path)
        ds = xr.Dataset ({"SWin": (("time"), forc[:,0]), 
                        "LWin": (("time"), forc[:,1]), 
                        "P": (("time"), forc[:,2]),
                        "T": (("time"), forc[:,3]),
                        "WindX": (("time"), forc[:,4]),
                        "WindY": (("time"), forc[:,5]),
                        "Press": (("time"), forc[:,6]),
                        "Hs": (("time"), forc[:,7]),
                        "LAI": (("time"), LAI),
                        "SAI": (("time"), SAI)},
            coords = {"time": pd.date_range(self.date_start, periods=len(forc), freq='30min'),
                "reference_time": pd.Timestamp(self.date_start)})
        self.forc= ds
        self.attrs_unit_forc()
        return ()

######################################### Lecture PARFLOW

    
    def process_ETR(self, dur) : 
        # calcul de la transpiration par horizon , La transp tot, l'ETP tot sur la colonne
        P= self.forc.P.sel(time=dur)*30*60 # conversion 30min en sec
        h = self.h
        Transpsum =  xr.where(np.abs(self.pf.z)<=np.abs(self.pf.z[self.pfidb.Solver.CLM.RootZoneNZ-1]), self.pf.evaptranssum*self.pf.vdz,0).sum(dim="z").data
        self.pf = self.pf.assign(Transpsum = (("time","y"),Transpsum))
        self.pf = self.pf.assign(ETRsum=(("time","y"),Transpsum))
        # On somme les valeur d'ETP en les multipliant par la hauteur de maille 
        for k in range(1,6) :
            #récupère l'objet horizon pour chaque horizon et stocke dans H, on somme sur les mailles de l'horizon
            H = self.pfidb.Geom['H'+str(k)]
            Transp = xr.where((np.abs(self.pf.z[30-H.Upper.Z]) <= np.abs(self.pf.z)) & (np.abs(self.pf.z) <= np.abs(self.pf.z[30-H.Lower.Z-1])),self.pf.evaptranssum*self.pf.vdz,0).sum(dim="z").data
            #on stocke dans le dataset pour récupérer facilement les index (a ce stade ce sont des ETP pas de Transp)
            self.pf['Transp'+str(k)] = (("time","y"),Transp)
        #on convertit pour avoir l'etp en mm
        for var in self.pf.data_vars :
            if 'Transp' in var :
                self.pf[var] = self.pf[var].mean(dim='y').sel(time=dur)
                self.pf[var] = - self.pf[var]*1000/h['dx']/h['dy'] #m3 -> mm
                #on retire P et evap_grnd pour avoir accès uniquement à la transpiration (uniquement dans la première maille)
                if ('1' in var) or ('sum'in var) :
                    self.pf[var] = self.pf[var] + P.resample(time ='D').sum() - self.clm.evap_grnd.resample(time ='D').sum()
        self.pf['ETRsum']=  self.pf['Transpsum'] + self.clm.evap_grnd.resample(time ='D').sum()
        self.attrs_unit_pf()
        return

    def process_WB(self):
        dur = slice(self.date_start,self.date_end)
        #calcul Transp et ETRsum
        self.process_ETR(dur)
        # calcul storage
        storage_sum = self.pf.storage.sel(time=dur).mean(dim='y').to_dataframe()['storage']*1000/self.h['dx']/self.h['dy']/self.h['ny']
        storage = storage_sum - storage_sum[0]
        P = self.forc.P.sel(time=dur).resample(time ='Y').sum()*30*60# conv 30min en s
        ETRsum = self.pf.ETRsum.sel(time=dur).isel(y=0).resample(time ='Y').sum()
        #print(P)
        Q = self.pf.Q.sel(time=dur).resample(time ='Y').sum()
        years = np.sort(storage.index.year.unique())
        WB = pd.DataFrame({'P':P,'ETR':ETRsum,'RunOff':Q,'WSC':np.nan, 'Inf':np.nan, 'CR':np.nan}, index=[str(year) for year in years])
        for countyear, year in enumerate(years):
            if countyear ==0:
                WB.WSC[countyear] = storage[storage.index.year==year][-1] - storage[storage.index.year==year][3]
            else:
                WB.WSC[countyear] = storage[storage.index.year==year][-1] - storage[storage.index.year==year-1][-1] 
            WB.CR[countyear] = Q[countyear]/P[countyear]*100
        WB['Inf'] = WB['P']- WB['ETR']-WB['RunOff']-WB['WSC']
        WB.dropna(inplace=True)
        #print (WB.CR)
        self.WB = WB
        return

    def attrs_unit_pf(self) :
        ds = self.pf
        for var in ds.data_vars :
            if ('press' in var) or ('y' in var) or ('z' in var) :
                ds[var].attrs['units']= 'm'
            elif 'storage' in var :
                ds[var].attrs['units']= 'm^3'
            elif ('satur' in var) or ('soil_moisture' in var) :
                ds[var].attrs['units']= '$m^3/m^3$'
            elif ('evap' in var) or ('overland' in var) :
                ds[var].attrs['units']= 'm^3/PFtimestep' 
            elif 'Transp' in var :
                ds[var].attrs['units']='mm'
            else :
                ds[var].attrs['units']='-'
        self.pf = ds
        return()

    def read_and_process_pf_outputs(self):
        h = self.h
        sdir = self.path
        # get pressure (pop(0) skips initial condition)
        files = np.sort(glob.glob(sdir+'*.press.*.pfb'))[1::]
        press = read_pfb_sequence(files)
        # get saturation
        files = np.sort(glob.glob(sdir+'*.satur.*.pfb'))[1::]
        satur = read_pfb_sequence(files)
        # get evaptranssum
        files = np.sort(glob.glob(sdir+'*.evaptranssum.*.pfb'))
        evaptranssum = read_pfb_sequence(files)
        # get overlandsum
        files = np.sort(glob.glob(sdir+'*.overlandsum.*.pfb'))
        overlandsum = read_pfb_sequence(files)
            
        # create dataset
        ds = xr.Dataset({"press": (("time","z","y"), press[:,::-1,:,0]),
                        "satur": (("time","z","y"), satur[:,::-1,:,0]),
                        "evaptranssum": (("time","z","y"), evaptranssum[:,::-1,:,0]),
                        "overlandsum":(("time","y"), overlandsum[:,0,:,0]),
                        'vdz': (("z","y"),self.var_dz[::-1,:,0])},
            coords={"x":np.arange(start = h['x'],stop = h['x']+h['nx']*h['dx'],step=h['dx']),
                "y":np.arange(start = h['y'],stop = h['y']+h['ny']*h['dy'],step=h['dy']),
                "z":-self.z, 
                "time": pd.date_range(self.date_start, periods=len(files)),
                "reference_time": pd.Timestamp(self.date_start)})
    
        # add WTD:
        ds = ds.assign(WTD = self.z[-1] - ds.press.isel(z=-1))
        # add auxiliary variables:
        ds = ds.assign(poro=(("z","y"),read_pfb(glob.glob(sdir+'*poro*.pfb')[0])[::-1,:,0]))
        ds = ds.assign(specstor=(("z","y"),read_pfb(glob.glob(sdir+'*specific*.pfb')[0])[::-1,:,0]))    
        # compute storage:
        ds = ds.assign(storage_cbyc=lambda x: h['dx']*h['dy']*h['dz']*x.satur*x.vdz*x.poro + \
                                            h['dx']*h['dy']*h['dz']*x.vdz*x.specstor*x.satur*x.press + \
                                            xr.where((x.z==x.z[0]) & (x.press>0),h['dx']*h['dy']*x.press,0))
        ds = ds.assign(storage=lambda x: x.storage_cbyc.sum(dim=('z')))
        ds = ds.assign(Q = lambda x : x.overlandsum.isel(y=0)*1000/h['dy']/h['ny']/h['dx'])
        ds = ds.assign(soil_moisture = lambda x: x.satur*x.poro)
        self.pf = ds
        self.process_WB()
        self.attrs_unit_pf()
        return ds


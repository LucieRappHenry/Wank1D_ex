#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 16:40:30 2024

@author: rappl
"""
import os,shutil
import time
from parflow import Run

class RUN :
    def __init__(self, proj_path, forc_path, pfidb, simu_name):
        self.proj_path = proj_path
        self.forc_path = forc_path
        self.pfidb = pfidb
        self.simu_name = simu_name

def run_funct(run_obj) :
    forc_path = run_obj.forc_path
    simu_path = run_obj.proj_path + 'simus/' + run_obj.simu_name + '/'
    if os.path.exists(simu_path):
        shutil.rmtree(simu_path) 
    os.mkdir(simu_path)
    wank = Run.from_definition(pfidb)
    # Copie les fichiers dans le répertoire de simulation
    shutil.copy(forc_path+'lai.dat',simu_path+'lai.dat')
    shutil.copy(forc_path+'sai.dat',simu_path+'sai.dat')
    shutil.copy(forc_path+'z0m.dat',simu_path+'z0m.dat')
    shutil.copy(forc_path+'displa.dat',simu_path+'displa.dat')
    shutil.copy(forc_path+'forc.txt',simu_path+'forcagePF.txt.0')
    shutil.copy(forc_path+'veg_map.pfb',simu_path+'veg_map.pfb')
    shutil.copy(forc_path+'drv_vegm.dat',simu_path+'drv_vegm.dat')
    shutil.copy(forc_path+'drv_vegp.dat',simu_path+'drv_vegp.dat')
    shutil.copy(forc_path+'drv_clmin.dat',simu_path+'drv_clmin.dat')
    shutil.copy(run_obj.pfidb,simu_path+'wank.pfidb')
    # distribute files that need be distributed
    wank.dist(simu_path +'veg_map.pfb')
    #OverwriteDrvClmin
    #lance la simu et sort le temps de calcul
    time1 = time.process_time()
    # lance la simu
    wank.run(working_directory=simu_path,skip_validation=True)
    time2 = time.process_time()
    print('temps de calcul pour ' + simu_name + ': ' , time2 - time1)
    with open(simu_path + "simu_time.out.txt", "a") as file:
        file.write("Temps de calcul via python : " + str(time2 - time1))



###################################################################
#Set wanted parflow version 
vers = 'LRH'
os.environ['PARFLOW_DIR'] = '/home/rappl/PFTree/Lucie/parflow-dev_'+vers+'/'

proj_path = "/home/rappl/PROJETS/Wank1D_ex/"
forc_path = proj_path + "/forcings/forc_"+vers+"/"
pfidb = proj_path +'tcl_pfidb/' + 'Wank1D.pfidb'


simu_name =vers
run_obj = RUN(proj_path, forc_path, pfidb, simu_name)
run_funct(run_obj)

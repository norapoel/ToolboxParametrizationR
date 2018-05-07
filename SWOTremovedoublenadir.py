# SWOTdatapreparation.py
"""The SWOTdatapreparation module is a toolbox developed specifically in preparation of the SWOT mission. It provides a toolbox to prepare SWOT like data for the data assimilation process. 

The main function is SWOTdatapreparation (same name as the module itself), and for standard applications, the user should not need to call other module functions.
# AUTHOR:
Nora Poel (1,2)
(1) CNRS/UGA/IRD/G-INP, IGE, Grenoble, France
(2) Institute for Computer Science, Potsdam, Germany

# HISTORY:
- April 2018: version 1
""" 

####################
# import libraries #
####################

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import xarray as xr
import re
from netCDF4 import Dataset
import os

#import SWOTaugmentation as swotaug
import SWOTdenoise as swotd


def read_ens(ensemble):
    """Read arrays from an ensemble of netcdf files.
    
    Parameters:
    ----------
    ensemble: input ensemble name
    
    Returns:
    -------
    ssh_model, ssh_obs, lat, lon, x_ac, time_sec: arrays
    """
# function reads an ensemble and returns a matrix containing the ssh at all grid points of all members #

### verify name of ensemble ###
    #regexEnsname = re.compile(r'[\w]\d{4,4}.nc.bas') #<name><nnnn>.nc.bas
    ensname = ensemble.split('/')[-1]
    #ensdir = ensemble.split(ensname)[0]
    #match = regexEnsname.search(ensname)
    #if match == None:
        #raise NameError("Ensemble not correctly named - see http://pp.ige-grenoble.fr/pageperso/brankarj/SESAM/ for naming input files")

### extract number of members of ensemble ###
    enssizelist = re.findall(r'\d{4,4}', ensname)
    enssize = int(enssizelist[0])

### get grid dimension ###
    membername1 = ensemble+'/vctgridSWOT0001.nc' 
    """eventually remove _denoised"""
    with xr.open_dataset(membername1,mask_and_scale=False) as dsmember1:
        time = dsmember1.time.size
        nc = dsmember1.nC.size
        lat = dsmember1.lat[:,:].values
        lon = dsmember1.lon[:,:].values
        x_ac = dsmember1.x_ac[:].values
        x_al = dsmember1.x_al[:].values
        time_sec = dsmember1.time_sec[:].values
        # first member of ensemble is later used to mask all members
        fill_value = dsmember1.ssh_model._FillValue
### load data ###
    ssh_model = np.zeros([enssize,time,nc]) # [layer, time, nC]
    ssh_obs = np.zeros([enssize,time,nc])
    for k in range (1,enssize+1):
        membername = (ensemble+'/vctgridSWOT{:04d}.nc').format(k) 
        with xr.open_dataset(membername,mask_and_scale=False) as ds:
            buf_model = ds.ssh_model[:,:]
            buf_obs = ds.ssh_obs[:,:]
            ssh_model[k-1,:,:] = buf_model
            ssh_obs[k-1,:,:] = buf_obs
    
    ssh_obs = ma.masked_where(ssh_obs==fill_value, ssh_obs)
    ssh_model = ma.masked_where(ssh_model==fill_value, ssh_model)
    ma.set_fill_value(ssh_obs, fill_value)
    ma.set_fill_value(ssh_model, fill_value)

### return data ###
    return ssh_model, ssh_obs, lat, lon, x_ac, x_al, time_sec


def write_ens(ensemble, ssh_model, ssh_obs, lat, lon, x_ac, x_al, time_sec):
    """Write arrays to netcdf files.
    
    Parameters:
    ----------
    ensemble: Name of ensemble
    ssh_model, ssh_obs, lat, lon, x_ac, time_sec: arrays
    
    Returns:
    -------
    outdir: Name of output directory"""
    
    _FillValue = ssh_model.fill_value
    
### create dir ###
    ensname = ensemble.split('/')[-1]
    ensdir = ensemble.split(ensname)[0]
    outdir = ensdir + 'rmdc_' + ensname #removedcenter_
    try:
        os.stat(outdir)
    except:
        os.mkdir(outdir)
          
### extract dimensions from input data ###
    enssize = ssh_model.shape[0]
    dimtime = ssh_model.shape[1] # ssh_model:[layer,time,nC]
    dimnc = ssh_model.shape[2]
    
### write data to netcdf ###
    for k in range (1,enssize+1):
        # create dataset
        ds = Dataset((outdir+'/vctgridSWOT{:04d}.nc').format(k), 'w') 
        # create dimensions #
        time = ds.createDimension("time", dimtime)
        nc = ds.createDimension("nC", dimnc)
        # create variables #
        varlat = ds.createVariable("lat","f8",("time", "nC"),fill_value=_FillValue)
        varlon = ds.createVariable("lon","f8",("time", "nC"),fill_value=_FillValue)
        varsshmodel = ds.createVariable("ssh_model","f8",("time", "nC"),fill_value=_FillValue)
        varsshobs = ds.createVariable("ssh_obs","f8",("time", "nC"),fill_value=_FillValue)
        varxac = ds.createVariable("x_ac","f8",("nC"),fill_value=_FillValue)
        varxal = ds.createVariable("x_al","f8",("time"),fill_value=_FillValue)
        vartimesec = ds.createVariable("time_sec","f8",("time"),fill_value=_FillValue)
        # fill variables #
        varlat[:] = lat[:,:]
        varlon[:] = lon[:,:]
        varsshmodel[:] = ssh_model[k-1,:,:]
        varsshobs[:] = ssh_obs[k-1,:,:]
        varxac[:] = x_ac[:]
        varxal[:] = x_al[:]
        vartimesec[:] = time_sec[:]
        ds.close()    
### return ###    
    return outdir

    
####################
# rm_centre_colomn #
####################

def rm_centre_colomn(ssh_model, ssh_obs, lon, lat, x_ac, x_al, time_sec):
    
    fill_value = ssh_model.fill_value
    
    enssize = ssh_model.shape[0]
    dimtime = ssh_model.shape[1]
    dimnc = ssh_model.shape[2]-1 
    
    lon = np.delete(lon, (60), axis=1)
    lat = np.delete(lat, (60), axis=1)
    x_ac = np.delete(x_ac, (60))
    
    ssh_model_rmd = np.zeros([enssize,dimtime,dimnc]) # [layer, time, nC]
    ssh_obs_rmd = np.zeros([enssize,dimtime,dimnc]) 
    
    for k in range (0,enssize):
        ssh_model_rmd[k,:,:] = np.delete(ssh_model[k,:,:], (60), axis=1)
        ssh_obs_rmd[k,:,:] = np.delete(ssh_obs[k,:,:], (60), axis=1)
        
    ssh_model_rmd = ma.masked_where(ssh_model_rmd==fill_value,ssh_model_rmd)
    ma.set_fill_value(ssh_model_rmd, fill_value)
    ssh_obs_rmd = ma.masked_where(ssh_obs_rmd==fill_value,ssh_obs_rmd)
    ma.set_fill_value(ssh_obs_rmd, fill_value)
    
    return ssh_model_rmd, ssh_obs_rmd, lon, lat, x_ac, x_al, time_sec


#######################
# SWOTdatapreparation #
#######################

def SWOTremovedoublenadir(ensemble):
    ssh_model, ssh_obs, lat, lon, x_ac, x_al, time_sec = read_ens(ensemble)
    
    ssh_model_rmd, ssh_obs_rmd, lon, lat, x_ac, x_al, time_sec = rm_centre_colomn(ssh_model, ssh_obs, lat, lon, x_ac, x_al, time_sec) #mask
    
    outdir = write_ens(ensemble, ssh_model_rmd, ssh_obs_rmd, lon, lat, x_ac, x_al, time_sec)
    
    print 'Ensemble without double nadir in ', outdir
    
    return 
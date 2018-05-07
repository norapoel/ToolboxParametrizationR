# SWOTaugmentation.py
"""The SWOTaugmentation module is a toolbox developed specifically in preparation of the SWOT mission. It provides a toolbox to add gradient observations to SWOT data. The main function is SWOTaugmentation (same name as the module itself), and for standard applications, the user should not need to call other module functions.
# AUTHOR:
Nora Poel (1,2)
(1) CNRS/UGA/IRD/G-INP, IGE, Grenoble, France
(2) Institute for Computer Science, Potsdam, Germany

# HISTORY:
- March 2018: version 1
""" 

### ATTENTION: Module might not be up to date -- compare to SWOTaugmentation.py!!


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

###########
# readens #
###########

def read_ens(ensemble,varname='ssh_model'):
    """Read arrays from an ensemble of netcdf files.
    
    Parameters:
    ----------
    ensemble: input ensemble name
    
    Returns:
    -------
    ssh_all, lat, lon: arrays
    """
# function reads an ensemble and returns a matrix containing the ssh at all grid points of all members #

### verify name of ensemble ###
    regexEnsname = re.compile(r'[\w]\d{4,4}.nc.bas') #<name><nnnn>.nc.bas
    ensname = ensemble.split('/')[-1]
    #ensdir = ensemble.split(ensname)[0]
    match = regexEnsname.search(ensname)
    if match == None:
        raise NameError("Ensemble not correctly named - see http://pp.ige-grenoble.fr/pageperso/brankarj/SESAM/ for naming input files")

### extract number of members of ensemble ###
    enssizelist = re.findall(r'\d{4,4}', ensname)
    enssize = int(enssizelist[0])
    
### verify number of members (?)###
    'TODO: prove if <nnnn> is equal number of members (+1 for mean in ...0000) in ensemble'

### verify backslashs in input (?)###
    'TODO: verify backslashs in input names / add+remove backslashs'

### get grid dimension ###
    membername1 = ensemble+'/vctgridSWOT0001.nc'
    with xr.open_dataset(membername1) as dsmember1:
        time = dsmember1.time.size
        nc = dsmember1.nC.size
        lat = dsmember1.lat[:,:].values
        lon = dsmember1.lon[:,:].values
    
### load data ###
    ssh_all = np.zeros([enssize,time,nc]) # [layer, time, nC]
    for k in range (1,enssize):
        membername = (ensemble+'/vctgridSWOT{:04d}.nc').format(k)
        with xr.open_dataset(membername) as ds:
            if (varname == 'ssh_obs'):
                buf = ds.ssh_obs[:,:]
            else:
                buf = ds.ssh_model[:,:]
            ssh_all[k-1,:,:] = buf
        
### return data ###
    return ssh_all, lat, lon


###########
# readobs #
###########

def read_obs(observation):
    """Read arrays from a netcdf file.
    
    Parameters:
    ----------
    observation: input observation name
    
    Returns:
    -------
    ssh_obs, lat, lon: arrays
    """
# function reads a observation and returns a matrix containing the ssh at all gridpoints of the observation #
    with xr.open_dataset(observation) as ds:
        time = ds.time.size
        nc = ds.nC.size
        lat = ds.lat[:,:].values
        lon = ds.lon[:,:].values
        ssh_obs = np.zeros([time,nc])
        ssh_obs[:,:] = ds.ssh_obs[:,:]
    
### return data ###
    return ssh_obs, lat, lon


############
# cut_edge #
############

def cut_edge(ssh_all):
    """Temporary function. Reduces the size of the data.
    
    Parameters:
    ----------
    ssh_all: array.
    
    Returns:
    -------
    reduced_ssh_all: array
    """
    # function to create a subdomain of the original area
    if (ssh_all.ndim == 3):
        enssize = ssh_all.shape[0]
        dimtime = ssh_all.shape[1]
        dimnc = ssh_all.shape[2]
        reduced_ssh_all = np.ones([enssize,dimtime,dimnc])
        reduced_ssh_all[:,:,:] = np.nan
        reduced_ssh_all[:,90:1210,6:95] = ssh_all[:,90:1210,6:95]
    
    if (ssh_all.ndim == 2):
        dimtime = ssh_all.shape[0]
        dimnc = ssh_all.shape[1]
        reduced_ssh_all = np.ones([dimtime,dimnc])
        reduced_ssh_all[:,:] = np.nan
        reduced_ssh_all[90:1210,6:95] = ssh_all[90:1210,6:95]
    
    ma.masked_where(np.isnan(reduced_ssh_all),reduced_ssh_all)
    
    return reduced_ssh_all

###############
# generate_d1 #
###############

def generate_d1(ssh_all):
    """Calculate the first derivations of the ssh across and along track (right difference)
    
    Parameters:
    ----------
    ssh_all: array
    
    Returns:
    -------
    d1c,d1a: arrays
    
    PROBLEM: SWOT simulator does not fill grid points that are not within the swaths. (ncdump -v ssh_obs xyz.nc shows a _ at these points). ncview fills these points with the _FillValue=2147483647.0 whereas python would like to fill them with missing_value which is not declared in the SWOT simulator output files. Python therefore takes an own value=9.96920996839e+36. SOLUTION: i) mask arrays with _PythonFillValue, ii) add a missing_value=2147483647.0 to variable ssh_obs in netcdf files and adjust this script to mask nan"""

    _FillValue = 2147483647.
    _PythonFillValue = 9969209968386869046778552952102584320
    
    # input: ensemble
    if (ssh_all.ndim == 3):
        isobs = False
        enssize = ssh_all.shape[0]
        dimtime = ssh_all.shape[1]
        dimnc = ssh_all.shape[2]
    # input: single observation
    if (ssh_all.ndim == 2):
        isobs = True
        enssize = 1
        dimtime = ssh_all.shape[0]
        dimnc = ssh_all.shape[1]
        tmp = np.zeros([enssize,dimtime,dimnc])
        tmp[0,:,:] = ssh_all
        ssh_all = tmp

    ### create first derivations across track ###
    
    ## calculation seperatly on left (a) and right (b) swath ##
    #operand1a = ssh_all[:,:,1:dimnc/2]
    #operand1b = ssh_all[:,:,dimnc/2+1:]
    #operand2a = ssh_all[:,:,:dimnc/2-1]
    #operand2b = ssh_all[:,:,dimnc/2:dimnc-1]
    # mask values
    #operand1a = ma.masked_where(operand1a==_PythonFillValue, operand1a)
    #operand1b = ma.masked_where(operand1b==_PythonFillValue, operand1b)
    #operand2a = ma.masked_where(operand2a==_PythonFillValue, operand2a)
    #operand2b = ma.masked_where(operand2b==_PythonFillValue, operand2b)
    # calculate derivation
    #d1ctmpa = operand1a-operand2a
    #d1ctmpb = operand1b-operand2b
    # return on grid with dimensions (dimtime, dimnc-1)
    #d1c = np.zeros([enssize,dimtime,dimnc])
    #d1c[:,:,:dimnc/2-1] = d1ctmpa[:,:,:] # fill grid from index [0...49]
    #d1c[:,:,dimnc/2:dimnc-1] = d1ctmpb[:,:,:] # fill grid from [51...100]
    # fill middle [51] column with _FillValue
    #d1c[:,:,dimnc/2-1] = _FillValue
    # fill last [101] column with _FillValue
    #d1c[:,:,dimnc-1] = _FillValue
    # mask the fill value
    #d1c = ma.masked_where((d1c==_FillValue), d1c)
    
    ## calculation on both swaths in once ##
    operand1 = ssh_all[:,:,1:]
    operand2 = ssh_all[:,:,:dimnc-1]
    # mask values
    operand1 = ma.masked_where(operand1==_PythonFillValue, operand1)
    operand2 = ma.masked_where(operand2==_PythonFillValue, operand2)
    # calculate derivation
    d1c = operand1-operand2
    # devide middle column (left of the 20km gap between swaths) by 20
    d1c[:,:,dimnc/2-1] = d1c[:,:,dimnc/2-1]/20.
    # mask the fill value
    d1c = ma.masked_where((d1c==_FillValue), d1c)

    ### create first derivations along track ###
    operand1 = ssh_all[:,1:,:]
    operand2 = ssh_all[:,:dimtime-1,:]
    # mask value
    operand1 = ma.masked_where(operand1==_PythonFillValue, operand1)
    operand2 = ma.masked_where(operand2==_PythonFillValue, operand2)
    # calculate derivation
    #d1atmp = operand1 - operand2
    d1a = operand1 - operand2
    # return on original grid
    #d1a = np.zeros([enssize,dimtime,dimnc])
    #d1a[:,:dimtime-1,:] = d1atmp[:,:,:]
    # fill last row with _FillValue
    #d1a[:,dimtime-1,:] = _FillValue
    # mask the fill value
    d1a = ma.masked_where(d1a==_FillValue,d1a)
    
    ### return data ###
    # reduce dimensions of d1a, d1c if input is a single observation
    if (isobs):
        return d1a[0,:,:], d1c[0,:,:]
    return d1a,d1c


###############
# generate_d2 #
###############

def generate_d2(ssh_all):
    """Calculate the second derivations of the ssh across and along track (central difference)
    
    Parameters:
    ----------
    ssh_all: array
    
    Returns:
    -------
    d2c,d2a: arrays"""

    _FillValue = 2147483647.
    _PythonFillValue = 9969209968386869046778552952102584320
    
    # input: ensemble
    if (ssh_all.ndim == 3):
        isobs = False
        enssize = ssh_all.shape[0]
        dimtime = ssh_all.shape[1]
        dimnc = ssh_all.shape[2]
    
    # input: single observation
    if (ssh_all.ndim == 2):
        isobs = True
        enssize = 1
        dimtime = ssh_all.shape[0]
        dimnc = ssh_all.shape[1]
        tmp = np.zeros([enssize,dimtime,dimnc])
        tmp[0,:,:] = ssh_all
        ssh_all = tmp
    
### create second derivation across track ###
    ## calculation seperatly on left (a) and right (b) swath ##
    #operand1a = ssh_all[:,:,2:dimnc/2] 
    #operand1b = ssh_all[:,:,dimnc/2+2:]
    #operand2a = ssh_all[:,:,1:dimnc/2-1]*2
    #operand2b = ssh_all[:,:,dimnc/2+1:dimnc-1]*2
    #operand3a = ssh_all[:,:,:dimnc/2-2]
    #operand3b = ssh_all[:,:,dimnc/2:dimnc-2]
    ## mask value
    #operand1a = ma.masked_where(operand1a==_PythonFillValue, operand1a)
    #operand1b = ma.masked_where(operand1b==_PythonFillValue, operand1b)
    #operand2a = ma.masked_where(operand2a==_PythonFillValue, operand2a)
    #operand2b = ma.masked_where(operand2b==_PythonFillValue, operand2b)
    #operand3a = ma.masked_where(operand3a==_PythonFillValue, operand3a)
    #operand3b = ma.masked_where(operand3b==_PythonFillValue, operand3b)
    ## calculate derivation
    #d2ctmpa = operand1a-operand2a+operand3a
    #d2ctmpb = operand1b-operand2b+operand3b
    ## return on original grid
    #d2c = np.zeros([enssize,dimtime,dimnc])
    #d2c[:,:,1:dimnc/2-1] = d2ctmpa[:,:,:] # fill grid from index [1...49]
    #d2c[:,:,dimnc/2+1:dimnc-1] = d2ctmpb[:,:,:] # fill grid from index [51...100]
    ## fill first, two middle, last column with _FillValue
    #d2c[:,:,0] = _FillValue # 0
    #d2c[:,:,dimnc/2-1] = _FillValue # 49
    #d2c[:,:,dimnc/2] = _FillValue # 50 
    #d2c[:,:,dimnc-1] = _FillValue # 101
    ## mask the fill value
    #d2c = ma.masked_where(d2c==_FillValue, d2c)
    
    ## calculation on both swaths in once ##
    operand1 = ssh_all[:,:,2:] 
    operand2 = ssh_all[:,:,1:dimnc-1]*2
    operand3 = ssh_all[:,:,:dimnc-2]
    # mask value
    operand1 = ma.masked_where(operand1==_PythonFillValue, operand1)
    operand2 = ma.masked_where(operand2==_PythonFillValue, operand2)
    operand3 = ma.masked_where(operand3==_PythonFillValue, operand3)
    # calculate derivation
    d2c = operand1-operand2+operand3
    # devide two middle columns (left and right of the 20km gap between swaths) by 20
    d2c[:,:,dimnc/2-2] = d2c[:,:,dimnc/2-2]/201. # /2-1 = 49 # /2 = 50
    d2c[:,:,dimnc/2-1] = d2c[:,:,dimnc/2-1]/201. # /2 = 50 # /2+1 = 51
    
    # mask the fill value
    d2c = ma.masked_where(d2c==_FillValue, d2c)
    
### create second derivation along track ###
    operand1 = ssh_all[:,2:,:] 
    operand2 = ssh_all[:,1:dimtime-1,:]*2
    operand3 = ssh_all[:,:dimtime-2,:]
    # mask value
    operand1 = ma.masked_where(operand1==_PythonFillValue, operand1)
    operand2 = ma.masked_where(operand2==_PythonFillValue, operand2)
    operand3 = ma.masked_where(operand1==_PythonFillValue, operand3)
    # calculate derivation
    #d2atmp = operand1-operand2+operand3
    d2a = operand1-operand2+operand3
    # return on original grid
    #d2a = np.zeros([enssize,dimtime,dimnc])
    #d2a[:,1:dimtime-1,:] = d2atmp[:,:,:]
    # fill first and last row with _FillValue
    #d2a[:,0,:] = _FillValue
    #d2a[:,dimtime-1,:] = _FillValue
    # mask the fill value
    d2a = ma.masked_where(d2a==_FillValue,d2a) 
    
### return data ###
    # reduce dimensions of d1a, d1c if input is a single observation
    if (isobs):
        return d2a[0,:,:],d2c[0,:,:]
    if (ssh_all.ndim == 3):
        return d2a,d2c
    
############
# writeens #
############

def write_ens(ensemble,lat,lon,d0=None,d1a=None,d1c=None,d2a=None,d2c=None):
    """Write arrays to netcdf files.
    
    Parameters:
    ----------
    ensemble: Name of ensemble
    lat,lon,d0=None,d1a=None,d1c=None,d2a=None,d2c=None: arrays
    
    Returns:
    -------
    outdir: Name of output directory"""
    ### create dir ###
    ensname = ensemble.split('/')[-1]
    ensdir = ensemble.split(ensname)[0]
    outdir = ensdir + 'ext_' + ensname
    try:
        os.stat(outdir)
    except:
        os.mkdir(outdir)
        
### extract number of members in ensemble
    enssizelist = re.findall(r'\d{4,4}', ensname)
    enssize = int(enssizelist[0])
    
### extract dimensions from input data ###
    dimtime = d0.shape[1] # d0:[layer,time,nC]
    dimnc = d0.shape[2]
    
### write data to netcdf ###
    for k in range (1,enssize+1):
        # create dataset
        ds = Dataset((outdir+'/vctgridSWOT{:04d}.nc').format(k), 'w') 
        # create dimensions #
        time = ds.createDimension("time", dimtime)
        time1 = ds.createDimension("time1", dimtime-1)
        time2 = ds.createDimension("time2", dimtime-2)
        nc = ds.createDimension("nC", dimnc)
        nc1 = ds.createDimension("nC1", dimnc-1)
        nc2 = ds.createDimension("nC2", dimnc-2)
        # create variables #
        'TODO: check if d0,d1a,... != None'
        varlat = ds.createVariable("lat","f8",("time", "nC"))
        varlon = ds.createVariable("lon","f8",("time", "nC"))
        varssh = ds.createVariable("ssh","f8",("time", "nC"))
        vard1ssha = ds.createVariable("d1ssha","f8",("time1", "nC"))
        vard1sshc = ds.createVariable("d1sshc","f8",("time", "nC1"))
        vard2ssha = ds.createVariable("d2ssha","f8",("time2", "nC"))
        vard2sshc = ds.createVariable("d2sshc","f8",("time", "nC2"))
        # fill variables #
        varlat[:] = lat[:,:]
        varlon[:] = lon[:,:]
        varssh[:] = d0[k-1,:,:]
        vard1ssha[:] = d1a[k-1,:,:]
        vard1sshc[:] = d1c[k-1,:,:]
        vard2ssha[:] = d2a[k-1,:,:]
        vard2sshc[:] = d2c[k-1,:,:]
    ds.close()    
### return ###    
    return outdir

############
# writeobs #
############

def write_obs(observation,lat,lon,d0=None,d1a=None,d1c=None,d2a=None,d2c=None):
    """Write array to netcdf file.
    
    Parameters:
    ----------
    observation: Name of observation
    lat,lon,d0=None,d1a=None,d1c=None,d2a=None,d2c=None: arrays
    
    Returns:
    -------
    outdir: Name of output directory"""
    ### extract dimensions from input data
    dimtime = d0.shape[0] # d0:[time,nC]
    dimnc = d0.shape[1]
    ### write data to netcdf
    # create dataset
    obsname = observation.split('/')[-1]
    obsdir = observation.split(obsname)[0]
    obsname = obsname.split('.nc')[0]+'_extended.nc'
    ds = Dataset(obsdir+obsname, 'w') 
    # create dimensions #
    time = ds.createDimension("time", dimtime)
    time1 = ds.createDimension("time1", dimtime-1)
    time2 = ds.createDimension("time2", dimtime-2)
    nc = ds.createDimension("nC", dimnc)
    nc1 = ds.createDimension("nC1", dimnc-1)
    nc2 = ds.createDimension("nC2", dimnc-2)
    # create variables #
    'TODO: check if d0,d1a,... != None'
    varlat = ds.createVariable("lat","f8",("time", "nC"))
    varlon = ds.createVariable("lon","f8",("time", "nC"))
    varssh = ds.createVariable("ssh","f8",("time", "nC"))
    vard1ssha = ds.createVariable("d1ssha","f8",("time1", "nC"))
    vard1sshc = ds.createVariable("d1sshc","f8",("time", "nC1"))
    vard2ssha = ds.createVariable("d2ssha","f8",("time2", "nC"))
    vard2sshc = ds.createVariable("d2sshc","f8",("time", "nC2"))
    # fill variables #
    varlat[:] = lat[:,:]
    varlon[:] = lon[:,:]
    varssh[:] = d0[:,:]
    vard1ssha[:] = d1a[:,:]
    vard1sshc[:] = d1c[:,:]
    vard2ssha[:] = d2a[:,:]
    vard2sshc[:] = d2c[:,:]
    ds.close()
### return ###    
    return obsdir+obsname

####################
# SWOTaugmentation #
####################

def SWOTaugmentation(filename,obs=False,ens=False):
    """Main function: reads the specified netcdf file(s), 
    calculates and annex first and second derivatives to inout file(s),
    write arrays to netcdf file(s).
    
    Parameters:
    ----------
    filename: Name of ensemble or observation
    obs,ens: bool
    
    Returns:
    -------"""
    if ((obs == False) & (ens == False)) | ((obs == True) & (ens == True)):
        print "Please specify the input type. Set obs=True OR ens=True"
        return
    if obs == True:
        ssh,lat,lon = read_obs(filename)
        d1a,d1c = generate_d1(ssh)
        d2a,d2c = generate_d2(ssh)
        outfilename = write_obs(filename,lat,lon,ssh,d1a,d1c,d2a,d2c)
        print 'Extended observation in ', outfilename
    else:
        ssh,lat,lon = read_ens(filename)
        d1a,d1c = generate_d1(ssh)
        d2a,d2c = generate_d2(ssh)
        outfilename = write_ens(filename,lat,lon,ssh,d1a,d1c,d2a,d2c)
        print 'Extended ensemble in ', outfilename
    return
    
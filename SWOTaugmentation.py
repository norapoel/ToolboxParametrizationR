# uncompyle6 version 3.1.2
# Python bytecode 2.7 (62211)
# Decompiled from: Python 2.7.14 |Anaconda custom (64-bit)| (default, Oct 16 2017, 17:29:19) 
# [GCC 7.2.0]
# Embedded file name: SWOTaugmentation.py
# Compiled at: 2018-05-03 16:16:50
"""The SWOTaugmentation module is a toolbox developed specifically in preparation of the SWOT mission. It provides a toolbox to add gradient observations to SWOT data. The main function is SWOTaugmentation (same name as the module itself), and for standard applications, the user should not need to call other module functions.
# AUTHOR:
Nora Poel (1,2)
(1) CNRS/UGA/IRD/G-INP, IGE, Grenoble, France
(2) Institute for Computer Science, Potsdam, Germany

# HISTORY:
- March 2018: version 1
"""
import numpy as np, numpy.ma as ma, matplotlib.pyplot as plt, xarray as xr, re
from netCDF4 import Dataset
import os

def read_ens(ensemble, varname='ssh_model'):
    """Read arrays from an ensemble of netcdf files.
    
    Parameters:
    ----------
    ensemble: input ensemble name
    
    Returns:
    -------
    ssh_all, lat, lon: arrays
    """
    regexEnsname = re.compile('[\\w]\\d{4,4}.nc.bas')
    ensname = ensemble.split('/')[-1]
    match = regexEnsname.search(ensname)
    if match == None:
        raise NameError('Ensemble not correctly named - see http://pp.ige-grenoble.fr/pageperso/brankarj/SESAM/ for naming input files')
    enssizelist = re.findall('\\d{4,4}', ensname)
    enssize = int(enssizelist[0])
    membername1 = ensemble + '/vctgridSWOT0001.nc'
    with xr.open_dataset(membername1, mask_and_scale=False) as (dsmember1):
        time = dsmember1.time.size
        nc = dsmember1.nC.size
        lat = dsmember1.lat[:, :].values
        lon = dsmember1.lon[:, :].values
        if varname == 'ssh_obs':
            fill_value = dsmember1.ssh_obs._FillValue
        else:
            fill_value = dsmember1.ssh_model._FillValue
    ssh_all = np.zeros([enssize, time, nc])
    for k in range(1, enssize + 1):
        membername = (ensemble + '/vctgridSWOT{:04d}.nc').format(k)
        with xr.open_dataset(membername, mask_and_scale=False) as (ds):
            if varname == 'ssh_obs':
                buf = ds.ssh_obs[:, :]
            else:
                buf = ds.ssh_model[:, :]
            ssh_all[k - 1, :, :] = buf

    ssh_all = ma.masked_where(ssh_all == fill_value, ssh_all)
    ma.set_fill_value(ssh_all, fill_value)
    return (
     ssh_all, lat, lon)


def read_obs(observation):
    """Read arrays from a netcdf file.
    
    Parameters:
    ----------
    observation: input observation name
    
    Returns:
    -------
    ssh_obs, lat, lon: arrays
    """
    with xr.open_dataset(observation, mask_and_scale=False) as (ds):
        time = ds.time.size
        nc = ds.nC.size
        lat = ds.lat[:, :].values
        lon = ds.lon[:, :].values
        ssh_obs = np.zeros([time, nc])
        ssh_obs[:, :] = ds.ssh_obs[:, :]
        fill_value = ds.ssh_obs._FillValue
        ssh_obs = ma.masked_where(ssh_obs == fill_value, ssh_obs)
        ma.set_fill_value(ssh_obs, fill_value)
    return (
     ssh_obs, lat, lon)


def write_ens(ensemble, lat, lon, d0, d1a, d1c, d2a, d2c):
    """Write arrays to netcdf files.
    
    Parameters:
    ----------
    ensemble: Name of ensemble
    lat,lon,d0,d1a,d1c,d2a,d2c: arrays
    
    Returns:
    -------
    outdir: Name of output directory"""
    _FillValue = d0.fill_value
    ensname = ensemble.split('/')[-1]
    ensdir = ensemble.split(ensname)[0]
    outdir = ensdir + 'ext_' + ensname
    try:
        os.stat(outdir)
    except:
        os.mkdir(outdir)

    enssizelist = re.findall('\\d{4,4}', ensname)
    enssize = int(enssizelist[0])
    dimtime = d0.shape[1]
    dimnc = d0.shape[2]
    for k in range(1, enssize + 1):
        ds0 = Dataset((outdir + '/vctgridSWOT{:04d}.nc').format(k), 'w')
        ds1a = Dataset((outdir + '/vctgridD1a{:04d}.nc').format(k), 'w')
        ds1c = Dataset((outdir + '/vctgridD1c{:04d}.nc').format(k), 'w')
        ds2a = Dataset((outdir + '/vctgridD2a{:04d}.nc').format(k), 'w')
        ds2c = Dataset((outdir + '/vctgridD2c{:04d}.nc').format(k), 'w')
        ds0.createDimension('time', dimtime)
        ds0.createDimension('nC', dimnc)
        ds1a.createDimension('time', dimtime - 1)
        ds1a.createDimension('nC', dimnc)
        ds1c.createDimension('time', dimtime)
        ds1c.createDimension('nC', dimnc - 1)
        ds2a.createDimension('time', dimtime - 2)
        ds2a.createDimension('nC', dimnc)
        ds2c.createDimension('time', dimtime)
        ds2c.createDimension('nC', dimnc - 2)
        varlat0 = ds0.createVariable('lat', 'f8', ('time', 'nC'), fill_value=_FillValue)
        varlon0 = ds0.createVariable('lon', 'f8', ('time', 'nC'), fill_value=_FillValue)
        varlat1a = ds1a.createVariable('lat', 'f8', ('time', 'nC'), fill_value=_FillValue)
        varlon1a = ds1a.createVariable('lon', 'f8', ('time', 'nC'), fill_value=_FillValue)
        varlat1c = ds1c.createVariable('lat', 'f8', ('time', 'nC'), fill_value=_FillValue)
        varlon1c = ds1c.createVariable('lon', 'f8', ('time', 'nC'), fill_value=_FillValue)
        varlat2a = ds2a.createVariable('lat', 'f8', ('time', 'nC'), fill_value=_FillValue)
        varlon2a = ds2a.createVariable('lon', 'f8', ('time', 'nC'), fill_value=_FillValue)
        varlat2c = ds2c.createVariable('lat', 'f8', ('time', 'nC'), fill_value=_FillValue)
        varlon2c = ds2c.createVariable('lon', 'f8', ('time', 'nC'), fill_value=_FillValue)
        varssh = ds0.createVariable('ssh', 'f8', ('time', 'nC'), fill_value=_FillValue)
        vard1ssha = ds1a.createVariable('d1a', 'f8', ('time', 'nC'), fill_value=_FillValue)
        vard1sshc = ds1c.createVariable('d1c', 'f8', ('time', 'nC'), fill_value=_FillValue)
        vard2ssha = ds2a.createVariable('d2a', 'f8', ('time', 'nC'), fill_value=_FillValue)
        vard2sshc = ds2c.createVariable('d2c', 'f8', ('time', 'nC'), fill_value=_FillValue)
        varlat0[:] = lat[:, :]
        varlon0[:] = lon[:, :]
        varlat1a[:] = np.delete(lat, dimtime - 1, axis=0)[:, :]
        varlon1a[:] = np.delete(lon, dimtime - 1, axis=0)[:, :]
        varlat1c[:] = np.delete(lat, dimnc - 1, axis=1)[:, :]
        varlon1c[:] = np.delete(lon, dimnc - 1, axis=1)[:, :]
        varlat2a[:] = np.delete(lat, [0, dimtime - 1], axis=0)[:, :]
        varlon2a[:] = np.delete(lon, [0, dimtime - 1], axis=0)[:, :]
        varlat2c[:] = np.delete(lat, [0, dimnc - 1], axis=1)[:, :]
        varlon2c[:] = np.delete(lon, [0, dimnc - 1], axis=1)[:, :]
        varssh[:] = d0[k - 1, :, :]
        vard1ssha[:] = d1a[k - 1, :, :]
        vard1sshc[:] = d1c[k - 1, :, :]
        vard2ssha[:] = d2a[k - 1, :, :]
        vard2sshc[:] = d2c[k - 1, :, :]
        ds0.close()
        ds1a.close()
        ds1c.close()
        ds2a.close()
        ds2c.close()

    return outdir


def write_obs(observation, lat, lon, d0=None, d1a=None, d1c=None, d2a=None, d2c=None):
    """Write array to netcdf file.
    
    Parameters:
    ----------
    observation: Name of observation
    lat,lon,d0=None,d1a=None,d1c=None,d2a=None,d2c=None: arrays
    
    Returns:
    -------
    outdir: Name of output directory"""
    _FillValue = d0.fill_value
    dimtime = d0.shape[0]
    dimnc = d0.shape[1]
    obsname = observation.split('/')[-1]
    obsdir = observation.split(obsname)[0]
    obsname = obsname.split('.nc')[0]
    ds0 = Dataset(obsdir + obsname + '_SWOT.nc', 'w')
    ds1a = Dataset(obsdir + obsname + '_D1a.nc', 'w')
    ds1c = Dataset(obsdir + obsname + '_D1c.nc', 'w')
    ds2a = Dataset(obsdir + obsname + '_D2a.nc', 'w')
    ds2c = Dataset(obsdir + obsname + '_D2c.nc', 'w')
    ds0.createDimension('time', dimtime)
    ds0.createDimension('nC', dimnc)
    ds1a.createDimension('time', dimtime - 1)
    ds1a.createDimension('nC', dimnc)
    ds1c.createDimension('time', dimtime)
    ds1c.createDimension('nC', dimnc - 1)
    ds2a.createDimension('time', dimtime - 2)
    ds2a.createDimension('nC', dimnc)
    ds2c.createDimension('time', dimtime)
    ds2c.createDimension('nC', dimnc - 2)
    varlat0 = ds0.createVariable('lat', 'f8', ('time', 'nC'), fill_value=_FillValue)
    varlon0 = ds0.createVariable('lon', 'f8', ('time', 'nC'), fill_value=_FillValue)
    varlat1a = ds1a.createVariable('lat', 'f8', ('time', 'nC'), fill_value=_FillValue)
    varlon1a = ds1a.createVariable('lon', 'f8', ('time', 'nC'), fill_value=_FillValue)
    varlat1c = ds1c.createVariable('lat', 'f8', ('time', 'nC'), fill_value=_FillValue)
    varlon1c = ds1c.createVariable('lon', 'f8', ('time', 'nC'), fill_value=_FillValue)
    varlat2a = ds2a.createVariable('lat', 'f8', ('time', 'nC'), fill_value=_FillValue)
    varlon2a = ds2a.createVariable('lon', 'f8', ('time', 'nC'), fill_value=_FillValue)
    varlat2c = ds2c.createVariable('lat', 'f8', ('time', 'nC'), fill_value=_FillValue)
    varlon2c = ds2c.createVariable('lon', 'f8', ('time', 'nC'), fill_value=_FillValue)
    varssh = ds0.createVariable('ssh', 'f8', ('time', 'nC'), fill_value=_FillValue)
    vard1ssha = ds1a.createVariable('d1a', 'f8', ('time', 'nC'), fill_value=_FillValue)
    vard1sshc = ds1c.createVariable('d1c', 'f8', ('time', 'nC'), fill_value=_FillValue)
    vard2ssha = ds2a.createVariable('d2a', 'f8', ('time', 'nC'), fill_value=_FillValue)
    vard2sshc = ds2c.createVariable('d2c', 'f8', ('time', 'nC'), fill_value=_FillValue)
    varlat0[:] = lat[:, :]
    varlon0[:] = lon[:, :]
    varlat1a[:] = np.delete(lat, dimtime - 1, axis=0)[:, :]
    varlon1a[:] = np.delete(lon, dimtime - 1, axis=0)[:, :]
    varlat1c[:] = np.delete(lat, dimnc - 1, axis=1)[:, :]
    varlon1c[:] = np.delete(lon, dimnc - 1, axis=1)[:, :]
    varlat2a[:] = np.delete(lat, [0, dimtime - 1], axis=0)[:, :]
    varlon2a[:] = np.delete(lon, [0, dimtime - 1], axis=0)[:, :]
    varlat2c[:] = np.delete(lat, [0, dimnc - 1], axis=1)[:, :]
    varlon2c[:] = np.delete(lon, [0, dimnc - 1], axis=1)[:, :]
    varssh[:] = d0[:, :]
    vard1ssha[:] = d1a[:, :]
    vard1sshc[:] = d1c[:, :]
    vard2ssha[:] = d2a[:, :]
    vard2sshc[:] = d2c[:, :]
    ds0.close()
    ds1a.close()
    ds1c.close()
    ds2a.close()
    ds2c.close()
    return obsdir


def generate_d1(ssh_all):
    """Calculate the first derivations of the ssh across and along track (right difference)
    
    Parameters:
    ----------
    ssh_all: array
    
    Returns:
    -------
    d1c,d1a: arrays
    
    PROBLEM: SWOT simulator does not fill grid points that are not within the swaths. (ncdump -v ssh_obs xyz.nc shows a _ at these points). ncview fills these points with the _FillValue=2147483647.0 whereas python would like to fill them with missing_value which is not declared in the SWOT simulator output files. Python therefore takes an own value=9.96920996839e+36. SOLUTION: i) mask arrays with _PythonFillValue, ii) add a missing_value=2147483647.0 to variable ssh_obs in netcdf files and adjust this script to mask nan"""
    fill_value = ssh_all.fill_value
    if ssh_all.ndim == 3:
        isobs = False
        enssize = ssh_all.shape[0]
        dimtime = ssh_all.shape[1]
        dimnc = ssh_all.shape[2]
    if ssh_all.ndim == 2:
        isobs = True
        enssize = 1
        dimtime = ssh_all.shape[0]
        dimnc = ssh_all.shape[1]
        tmp = np.zeros([enssize, dimtime, dimnc])
        tmp[0, :, :] = ssh_all
        ssh_all = tmp
    ssh_all = ma.masked_where(ssh_all == fill_value, ssh_all)
    operand1 = ssh_all[:, :, 1:]
    operand2 = ssh_all[:, :, :dimnc - 1]
    d1c = operand1 - operand2
    d1c = ma.masked_where(d1c == fill_value, d1c)
    ma.set_fill_value(d1c, fill_value)
    operand1 = ssh_all[:, 1:, :]
    operand2 = ssh_all[:, :dimtime - 1, :]
    d1a = operand1 - operand2
    d1a = ma.masked_where(d1a == fill_value, d1a)
    ma.set_fill_value(d1a, fill_value)
    if isobs:
        return (d1a[0, :, :], d1c[0, :, :])
    return (d1a, d1c)


def generate_d2(ssh_all):
    """Calculate the second derivations of the ssh across and along track (central difference)
    
    Parameters:
    ----------
    ssh_all: array
    
    Returns:
    -------
    d2c,d2a: arrays"""
    fill_value = ssh_all.fill_value
    if ssh_all.ndim == 3:
        isobs = False
        enssize = ssh_all.shape[0]
        dimtime = ssh_all.shape[1]
        dimnc = ssh_all.shape[2]
    if ssh_all.ndim == 2:
        isobs = True
        enssize = 1
        dimtime = ssh_all.shape[0]
        dimnc = ssh_all.shape[1]
        tmp = np.zeros([enssize, dimtime, dimnc])
        tmp[0, :, :] = ssh_all
        ssh_all = tmp
    ssh_all = ma.masked_where(ssh_all == fill_value, ssh_all)
    operand1 = ssh_all[:, :, 2:]
    operand2 = ssh_all[:, :, 1:dimnc - 1] * 2
    operand3 = ssh_all[:, :, :dimnc - 2]
    d2c = operand1 - operand2 + operand3
    d2c = ma.masked_where(d2c == fill_value, d2c)
    ma.set_fill_value(d2c, fill_value)
    operand1 = ssh_all[:, 2:, :]
    operand2 = ssh_all[:, 1:dimtime - 1, :] * 2
    operand3 = ssh_all[:, :dimtime - 2, :]
    d2a = operand1 - operand2 + operand3
    d2a = ma.masked_where(d2a == fill_value, d2a)
    ma.set_fill_value(d2a, fill_value)
    if isobs:
        return (d2a[0, :, :], d2c[0, :, :])
    if ssh_all.ndim == 3:
        return (d2a, d2c)


def write_Rplus(observation, varC, alpha0a, alpha0c, alpha1a, alpha1c, alpha2a, alpha2c):
    """ Function writes 5 netcdf files, each one containing the error covariances of an other (gradient) observation.
    
    Parameters:
    ----------
    observation: filename of an observation for which Rplus will be created
    varC: variance across track, array
    a0a,a0c,a1a,a1c,a2a,a2c: parameter for R-matrix, double
    
    Returns:
    """
    with xr.open_dataset(observation) as (ds):
        dimtime = ds.time.size
        dimnc = ds.nC.size
    row0 = varC * alpha0a
    R0 = np.tile(row0, (dimtime, 1))
    row1c = varC * alpha1c
    row1c = np.delete(row1c, dimnc - 1)
    R1c = np.tile(row1c, (dimtime, 1))
    row2c = varC * alpha2c
    row2c = np.delete(row2c, [0, dimnc - 1])
    R2c = np.tile(row2c, (dimtime, 1))
    ratio1 = alpha1a / alpha0a
    row1a = row0 * ratio1
    R1a = np.tile(row1a, (dimtime - 1, 1))
    ratio2 = alpha2a / alpha0a
    row2a = row0 * ratio2
    R2a = np.tile(row2a, (dimtime - 2, 1))
    obsname = observation.split('/')[-1]
    outdir = observation.split(obsname)[0]
    obsname = 'R_gridSWOT.nc'
    ds0 = Dataset(outdir + obsname, 'w')
    time = ds0.createDimension('time', dimtime)
    nc = ds0.createDimension('nC', dimnc)
    Ra0 = ds0.createVariable('ssh', 'f8', ('time', 'nC'))
    Ra0[:] = R0[:, :]
    ds0.close()
    obsname = 'R_gridD1a.nc'
    ds1a = Dataset(outdir + obsname, 'w')
    time = ds1a.createDimension('time', dimtime - 1)
    nc = ds1a.createDimension('nC', dimnc)
    Ra1a = ds1a.createVariable('d1a', 'f8', ('time', 'nC'))
    Ra1a[:] = R1a[:, :]
    ds1a.close()
    obsname = 'R_gridD1c.nc'
    ds1c = Dataset(outdir + obsname, 'w')
    time = ds1c.createDimension('time', dimtime)
    nc = ds1c.createDimension('nC', dimnc - 1)
    Ra1c = ds1c.createVariable('d1c', 'f8', ('time', 'nC'))
    Ra1c[:] = R1c[:, :]
    ds1c.close()
    obsname = 'R_gridD2a.nc'
    ds2a = Dataset(outdir + obsname, 'w')
    time = ds2a.createDimension('time', dimtime - 2)
    nc = ds2a.createDimension('nC', dimnc)
    Ra2a = ds2a.createVariable('d2a', 'f8', ('time', 'nC'))
    Ra2a[:] = R2a[:, :]
    ds2a.close()
    obsname = 'R_gridD2c.nc'
    ds2c = Dataset(outdir + obsname, 'w')
    time = ds2c.createDimension('time', dimtime)
    nc = ds2c.createDimension('nC', dimnc - 2)
    Ra2c = ds2c.createVariable('d2c', 'f8', ('time', 'nC'))
    Ra2c[:] = R2c[:, :]
    ds2c.close()
    print 'R in ', outdir


def SWOTaugmentation(filename, obs=False, ens=False):
    """Main function: reads the specified netcdf file(s), 
    calculates and annex first and second derivatives to inout file(s),
    write arrays to netcdf file(s).
    
    Parameters:
    ----------
    filename: Name of ensemble or observation
    obs,ens: boolean
    
    Returns:
    -------"""
    if (obs == False) & (ens == False) | (obs == True) & (ens == True):
        print 'Please specify the input type. Set obs=True OR ens=True'
        return
    if obs == True:
        ssh, lat, lon = read_obs(filename)
        d1a, d1c = generate_d1(ssh)
        d2a, d2c = generate_d2(ssh)
        outfilename = write_obs(filename, lat, lon, ssh, d1a, d1c, d2a, d2c)
        print 'Extended observation in ', outfilename
    else:
        ssh, lat, lon = read_ens(filename)
        d1a, d1c = generate_d1(ssh)
        d2a, d2c = generate_d2(ssh)
        outfilename = write_ens(filename, lat, lon, ssh, d1a, d1c, d2a, d2c)
        print 'Extended ensemble in ', outfilename
# okay decompiling SWOTaugmentation.pyc

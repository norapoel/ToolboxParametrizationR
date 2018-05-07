# uncompyle6 version 3.1.2
# Python bytecode 2.7 (62211)
# Decompiled from: Python 2.7.14 |Anaconda custom (64-bit)| (default, Oct 16 2017, 17:29:19) 
# [GCC 7.2.0]
# Embedded file name: SWOTdenoise.py
# Compiled at: 2018-04-27 12:31:41
"""The SWOTdenoise module is a toolbox developed specifically in preparation of the SWOT mission. It provides a toolbox to remove small-scale noise from SWOT data. The main function is SWOTdenoise (same name as the module itself), and for standard applications, the user should not need to call other module functions. Optionally, other functions that can be directly useful are read_data (to read data from a Netcdf file) and fill_nadir_gap: this function fills the lon and lat arrays in the SWOT nadir gap, and introduces fill values in the SSH array. Look at dedicated helps.

# AUTHORS:
Laura Gomez Navarro (1,2), Emmanuel Cosme (1), Nicolas Papadakis (3)
(1) CNRS/UGA/IRD/G-INP, IGE, Grenoble, France
(2) IMEDEA (CSIC-UIB), Esporles, Spain
(3) CNRS/Univ. Bordeaux/B-INP, IMB, Bordeaux, France

# HISTORY:
- March 2018: version 1
- March 2018: version 2 

# CHANGES in version 2: netcdf dimensions and variables now adapted from SWOT simulator version 2.3 to version 3.0
"""
import numpy as np
from netCDF4 import Dataset
from scipy import ndimage as nd
from scipy.interpolate import RectBivariateSpline
from types import *
import sys

def read_data(filename, *args):
    """Read arrays from netcdf file.
    
    Parameters:
    ----------
    filename: input file name
    *args: strings, variables to be read as named in the netcdf file.
    
    Returns:
    -------
    arrays. The number of output arrays must be identical to the number of variables.
    """
    fid = Dataset(filename)
    output = []
    for entry in args:
        output.append(fid.variables[entry][:])

    fid.close()
    return tuple(output)


def write_data(filename, ssh_d, lon_d, lat_d, x_ac_d, time_d):
    """
    Write SSH in output file.
    
    Parameters:
    ----------
    filename: output filename
    ssh_d, lon_d, lat_d, x_ac_d, time_d: standard SWOT data arrays. See SWOTdenoise function.
    
    Returns:
    -------
    Outpur file name.
    """
    _FillValue = ssh_d.fill_value
    rootname = filename.split('.nc')[0]
    filenameout = rootname + '_denoised.nc'
    x_al_r = read_data(filename, 'x_al')
    fid = Dataset(filenameout, 'w', format='NETCDF4')
    fid.description = 'Filtered SWOT data'
    fid.creator_name = 'SWOTdenoise module'
    time = fid.createDimension('time', len(time_d))
    x_ac = fid.createDimension('nC', len(x_ac_d))
    lat = fid.createVariable('lat', 'f8', ('time', 'nC'))
    lat.long_name = 'Latitude'
    lat.units = 'degrees_north'
    lat[:] = lat_d
    lon = fid.createVariable('lon', 'f8', ('time', 'nC'))
    lon.long_name = 'Longitude'
    lon.units = 'degrees_east'
    lon[:] = lon_d
    vtime = fid.createVariable('time_sec', 'f8', 'time')
    vtime.long_name = 'time from beginning of simulation (in s)'
    vtime.units = 's'
    vtime[:] = time_d * 86400
    x_al = fid.createVariable('x_al', 'f8', 'time')
    x_al.long_name = 'Along track distance from the beginning of the pass'
    x_al.units = 'km'
    x_al[:] = x_al_r
    vx_ac = fid.createVariable('x_ac', 'f8', 'nC')
    vx_ac.long_name = 'Across track distance from nadir'
    vx_ac.units = 'km'
    vx_ac[:] = x_ac_d
    ssh = fid.createVariable('ssh_obs', 'f8', ('time', 'nC'), fill_value=_FillValue)
    ssh.long_name = 'SSH denoised'
    ssh.units = 'm'
    ssh[:] = ssh_d
    fid.close()
    return filenameout


def copy_arrays(*args):
    """numpy-copy arrays.
    
    Parameters:
    ----------
    *args: arrays to copy.
    
    Returns:
    -------
    arrays. The number of output arrays must be identical to the number of inputs.
    """
    output = []
    for entry in args:
        output.append(entry.copy())

    return tuple(output)


def fill_nadir_gap(ssh, lon, lat, x_ac, time, method='fill_value'):
    """
    Fill the nadir gap in the middle of SWOT swath.
    Longitude and latitude are interpolated linearly. For SSH, there are two options:
        If the gap is already filled in the input arrays, it returns the input arrays.
        Parameters:
    ----------
    ssh, lon, lat, x_ac, time: input masked arrays from SWOT data. See SWOTdenoise function.
    method: method used to fill SSH array in the gap. Two options:
        - 'fill_value': the gap is filled with the fill value of SSH masked array;
        - 'interp': the gap is filled with a 2D, linear interpolation.
        Returns:
    -------
    ssh_f, lon_f, lat_f, x_ac_f: Filled SSH (masked), lon, lat 2D arrays, and across-track coordinates.
    """
    nhsw = len(x_ac) / 2
    step = abs(x_ac[nhsw + 1] - x_ac[nhsw])
    ins = np.arange(x_ac[nhsw - 1], x_ac[nhsw], step)[1:]
    nins = len(ins)
    if nins == 0:
        lon_f = lon
        lat_f = lat
        x_ac_f = x_ac
        ssh_f = ssh
    else:
        x_ac_f = np.insert(x_ac, nhsw, ins)
        lon_f = RectBivariateSpline(time, x_ac, lon)(time, x_ac_f)
        lat_f = RectBivariateSpline(time, x_ac, lat)(time, x_ac_f)
        if np.ma.isMaskedArray(ssh) == False:
            ssh = np.ma.asarray(ssh)
            print 'ssh had to be masked'
        if method == 'interp':
            ssh_f = np.ma.masked_values(RectBivariateSpline(time, x_ac, ssh)(time, x_ac_f), ssh.fill_value)
        else:
            ins_ssh = np.full((nins, len(time)), ssh.fill_value, dtype='float32')
            ssh_f = np.ma.masked_values(np.insert(ssh, nhsw, ins_ssh, axis=1), ssh.fill_value)
    return (ssh_f, lon_f, lat_f, x_ac_f)


def empty_nadir_gap(ssh_f, x_ac_f, ssh, x_ac):
    """
    Remove entries of the nadir gap from ssh array.
        Parameters:
    ----------
    ssh_f: input 2D masked array of SSH data with filled gap
    x_ac_f: across-track coordinates of ssh_f
    ssh: 2D masked array of original SWOT SSH, with the gap
    x_ac: across-track coordinates of ssh
        Returns:
    -------
    2D masked array is of the same shape as the initial SWOT array.
    """
    ninter = len(x_ac_f) - len(x_ac)
    if ninter != 0:
        nx = (np.shape(ssh_f)[1] - ninter) / 2
        ssh_out = np.concatenate([ssh_f[:, 0:nx], ssh_f[:, -nx:]], axis=1)
        ssh_out = np.ma.array(ssh_out, mask=ssh.mask, fill_value=ssh.fill_value)
    else:
        ssh_out = ssh_f
        return ssh_out


def convolution_filter(ssh, param, method):
    """
    Filter an image with a convolution of a generic function (Gaussian or boxcar).
    The input image can contain gaps (masked values).
    Gaps are filled with 0. An array of 1 is created with gaps set to 0. Both are filtered and divided. Inspired from
    http://stackoverflow.com/questions/18697532/gaussian-filtering-a-image-with-nan-in-python
    This function calls scipy.ndimage.
    
    Parameters:
    ----------
    ssh: 2D masked array to filter
    param: parameter for the method:
        - standard deviation for the Gaussian
        - box size for boxcar
    method: Gaussian or boxcar.
    
    Returns:
    -------
    2D ndarray (not a masked array).
    """
    assert np.ma.any(ssh.mask), 'u must be a masked array'
    mask = np.flatnonzero(ssh.mask)
    v = ssh.data.copy()
    v.flat[mask] = 0
    w = np.ones_like(ssh.data)
    w.flat[mask] = 0
    if method == 'boxcar':
        param = int(param)
        v[:] = nd.generic_filter(v, function=np.nanmean, size=param)
        w[:] = nd.generic_filter(w, function=np.nanmean, size=param)
    else:
        if method == 'gaussian':
            v[:] = nd.gaussian_filter(v, sigma=param)
            w[:] = nd.gaussian_filter(w, sigma=param)
        else:
            write_error_and_exit(2)
    w = np.clip(w, 1e-08, 1.0)
    return v / w


def gradx(I):
    """
    Calculates the gradient in the x-direction of an image I and gives as output M.
    In order to keep the size of the initial image the last row is left as 0s.
    """
    m, n = I.shape
    M = np.zeros([m, n])
    M[0:-1, :] = np.subtract(I[1::, :], I[0:-1, :])
    return M


def grady(I):
    """
    Calculates the gradient in the y-direction of an image I and gives as output M.
    In order to keep the size of the initial image the last column is left as 0s.
    """
    m, n = I.shape
    M = np.zeros([m, n])
    M[:, 0:-1] = np.subtract(I[:, 1::], I[:, 0:-1])
    return M


def div(px, py):
    """
    Calculates the divergence of a 2D field. 
    For the specific application of image denoising, the calculation follows Chambolle (REF)
    ## BELOW, TO BE CLARIFIED
    The x component of M (Mx) first row is = to the first row of px.
    The x component of M (Mx) last row is = to - the before last row of px. (last one = 0)
    The y component of M (My) first column is = to the first column of py.
    The y component of M (My) last column is = to - the before last column of py. (last one = 0)
    ??#(de sorte que div=-(grad)^*)
    Parameters: two 2D ndarray
    Returns: 2D ndarray
    """
    m, n = px.shape
    M = np.zeros([m, n])
    Mx = np.zeros([m, n])
    My = np.zeros([m, n])
    Mx[1:m - 1, :] = px[1:m - 1, :] - px[0:m - 2, :]
    Mx[0, :] = px[0, :]
    Mx[m - 1, :] = -px[m - 2, :]
    My[:, 1:n - 1] = py[:, 1:n - 1] - py[:, 0:n - 2]
    My[:, 0] = py[:, 0]
    My[:, n - 1] = -py[:, n - 2]
    M = Mx + My
    return M


def laplacian(u):
    """
    Calculates laplacian using the divergence and gradient functions defined in the module.
    Parameter: 2D ndarray
    Returns: 2D ndarray
    """
    Ml = div(gradx(u), grady(u))
    return Ml


def variational_regularization_filter(ssh, param, itermax=10000, epsilon=1e-09):
    """
    Apply variational regularization filter. 
    
    
    Parameters:
    ----------
    ssh: masked array with nadir gap filled.
    param: 3-entry tuple for first, second, and third order terms of the cost function, respectively.
    itermax: maximum number of iterations in the gradient descent method.
    epsilon: for convergence criterium.
    
    Returns:
    -------
    ssh_d: 2D ndarray containing denoised ssh data (ssh_d is not a masked array!)
    """
    ssh_d = convolution_filter(ssh, 10.0, method='gaussian')
    tau = np.min((1.0 / (1 + 8 * param[0]), 1.0 / (1 + 64 * param[1]), 1.0 / (1 + 512 * param[2])))
    print tau
    mask = 1 - ssh.mask
    iteration = 1
    while iteration < itermax:
        iteration += 1
        ssh_tmp = np.copy(ssh_d)
        lap_tmp = laplacian(ssh_tmp)
        bilap_tmp = laplacian(lap_tmp)
        incr = mask * (ssh.data - ssh_tmp) + param[0] * lap_tmp - param[1] * bilap_tmp + param[2] * laplacian(bilap_tmp)
        ssh_d = ssh_tmp + tau * incr
        norm = np.ma.sum(mask * incr * incr) / np.sum(mask)
        if norm < epsilon:
            break

    print iteration, norm / epsilon
    return ssh_d


def write_error_and_exit(nb):
    """Function called in case of error, to guide the user towards appropriate adjustment."""
    if nb == 1:
        print 'You must provide a SWOT file name OR SSH, lon, lat, x_ac and time arrays. SSH must be a masked array.'
    if nb == 2:
        print 'The filtering method is not correctly set.'
    if nb == 3:
        print 'For the variational regularization filter, lambd must be a 3-entry tuple.'
    sys.exit()


def SWOTdenoise(*args, **kwargs):
    """
    Perform denoising of SWOT data.
    
    Parameters:
    ----------
    *args: name of file containing the SWOT SSH field to denoise (optional). Example of use:
        SWOTdenoise(filename)
        denoise data in file 'filename' and write an output file in the same directory. 
    
        The output file is named 'foo_denoised.nc' if the input file name is 'foo.nc'.
        
    **kwargs include:
    - ssh : input ssh array (2D)
    - lon : input longitude array (2D)
    - lat : input latitude array (2D)
    - x_ac : input across-track coordinates (1D)
    - time : input along-track coordinates (1D)
    The above-mentioned arguments are mandatory if no file name is given. They are exactly in the format provided by the SWOT simulator for ocean science version 3.0.
    
    Other keywords arguments are:
    - method: gaussian, boxcar, or var_reg (default);
    - param: number for gaussian and boxcar; 3-entry tuple for var_reg (default: (1.5, 0, 0); underinvestigation) ;
    - itermax: only for var_reg: maximum number of iterations in the gradient descent algortihm (default: 10000);
    - epsilon: only for var_reg: convergence criterium for the gradient descent algortihm (default: 1e-9);
    - inpainting: if True, the nadir gap is inpainted. If False, it is not and the returned SSH array is of the same shape as the original one. If the SWOTdenoise function is called using arrays (see above description) with inpainting=True, then it returns SSH, lon, and lat arrays. If it is called using arrays with inpainting=False, it returns only SSH, since lon and lat arrays are the same as for the input field. Default is False.
    
    The algorithms are detailed in the scientific documentation.
    
    """
    file_input = len(args) == 1
    if file_input:
        if type(args[0]) is not str:
            write_error_and_exit(1)
        filename = args[0]
        swotfile = filename.split('/')[-1]
        swotdir = filename.split(swotfile)[0]
        listvar = ('ssh_obs', 'lon', 'lat', 'x_ac', 'time_sec')
        ssh, lon, lat, x_ac, time = read_data(filename, *listvar)
    else:
        ssh = kwargs.get('ssh', None)
        lon = kwargs.get('lon', None)
        lat = kwargs.get('lat', None)
        x_ac = kwargs.get('x_ac', None)
        time = kwargs.get('time', None)
        if any((isinstance(ssh, NoneType), isinstance(lon, NoneType), isinstance(lat, NoneType),
         isinstance(x_ac, NoneType), isinstance(time, NoneType))):
            write_error_and_exit(1)
        method = kwargs.get('method', 'var_reg')
        param = kwargs.get('param', (1.5, 0, 0))
        itermax = kwargs.get('itermax', 10000)
        epsilon = kwargs.get('epsilon', 1e-09)
        inpainting = kwargs.get('inpainting', False)
        ssh_f, lon_f, lat_f, x_ac_f = fill_nadir_gap(ssh, lon, lat, x_ac, time)
        if method == 'do_nothing':
            ssh_d, lon_d, lat_d, x_ac_d, time_d = copy_arrays(ssh, lon, lat, x_ac, time)
        if method == 'boxcar':
            ssh_d = convolution_filter(ssh_f, param, method='boxcar')
        if method == 'gaussian':
            ssh_d = convolution_filter(ssh_f, param, method='gaussian')
        if method == 'var_reg':
            if len(param) is not 3:
                write_error_and_exit(3)
            ssh_d = variational_regularization_filter(ssh_f, param, itermax=itermax, epsilon=epsilon)
        if inpainting is True:
            ssh_tmp, _, _, _ = fill_nadir_gap(ssh, lon, lat, x_ac, time, method='interp')
            ssh_d = np.ma.array(ssh_d, mask=ssh_tmp.mask, fill_value=ssh.fill_value)
            lon_d, lat_d, x_ac_d = copy_arrays(lon_f, lat_f, x_ac_f)
        else:
            ssh_d = np.ma.array(ssh_d, mask=ssh_f.mask, fill_value=ssh.fill_value)
            ssh_d = empty_nadir_gap(ssh_d, x_ac_f, ssh, x_ac)
            lon_d, lat_d, x_ac_d = copy_arrays(lon, lat, x_ac)
        if file_input:
            filenameout = write_data(filename, ssh_d, lon_d, lat_d, x_ac_d, time)
            print 'Filtered field in ', filenameout
        else:
            if inpainting is True:
                return (ssh_d, lon_d, lat_d)
        return ssh_d
    return
# okay decompiling SWOTdenoise.pyc

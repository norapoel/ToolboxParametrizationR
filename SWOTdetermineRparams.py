# SWOTdetermineRparams.py
"""The SWOTparameterR module is a toolbox developed specifically in preparation of the SWOT mission. It provides functions to find optimal parameters for the R+ matrix. The main function is SWOTdetermineRparams (same name as the module itself), and for standard applications, the user should not need to call other module functions.
# AUTHOR:
Nora Poel (1,2)
(1) CNRS/UGA/IRD/G-INP, IGE, Grenoble, France
(2) Institute for Computer Science, Potsdam, Germany

# HISTORY:
- March 2018: version 1
""" 

####################
# import libraries #
####################

import numpy as np
import numpy.ma as ma
import numpy.linalg as linalg
import matplotlib.pyplot as plt
import xarray as xr
import re
from netCDF4 import Dataset
import os

# for cost function
from matplotlib import cm
import scipy.optimize as so

import SWOTaugmentation as swota
import SWOTdenoise as swotd

'''The value of fill_value that is provided by the SWOT simulator. To be adapted when it is changed'''
fill_value=2147483647.0
##############
# create_TTc #
##############

def create_TTc(dimnc):
    t = np.diag(np.ones(dimnc))
    t = ma.masked_where(t==fill_value,t) # necessary for SWOTaugmentation
    invalid,T1cT = swota.generate_d1(t)
    T1c = np.transpose(T1cT)
    dotT1c = np.dot(T1cT,T1c)
    
    invalid,T2cT = swota.generate_d2(t)
    T2c = np.transpose(T2cT)
    dotT2c = np.dot(T2cT,T2c)
    
    return dotT1c, dotT2c

##############
# create_TTa #
##############

def create_TTa(dimtime):
    t = np.diag(np.ones(dimtime))
    t = ma.masked_where(t==fill_value,t) # necessary for SWOTaugmentation
    T1a,invalid = swota.generate_d1(t)
    T1aT = np.transpose(T1a)
    dotT1a = np.dot(T1aT,T1a)
    
    T2a,invalid = swota.generate_d2(t)
    T2aT = np.transpose(T2a)
    dotT2a = np.dot(T2aT,T2a)
    
    return dotT1a, dotT2a

##############
# corrAcross #
##############

def corrAcross(ssh,time):
    '''Function to calculate the correlation along one line (across track)
    
    Parameters:
    -----------
    ssh: 3d-array
    time: float. To chose a point
    
    Output:
    -------
    corr,cov: 2d-arrays (dimension: width of trace x width of trace)'''
    # number of members
    enssize = ssh.shape[0]
    # width of trace
    dimnc = ssh.shape[2]
    
    # set dimensions of final correlation and covariance matrix
    corr = np.zeros([dimnc,dimnc])
    cov = np.zeros([dimnc,dimnc])
    # calculate correlation between each point along the line 
    for i in range (0,dimnc):
        # xi - E(xi)
        pointmean = np.nanmean(ssh[:,time,i])
        pointvar = np.var(ssh[:,time,i])
        operand1 = ssh[:,time,i]-pointmean
        for j in range (0,dimnc):
            # xj - E(xj)
            operand2 = ssh[:,time,j]-np.nanmean(ssh[:,time,j])
            # Cov(xi,xj) = E[operand1*operand2^T]
            cov_tmp = (np.sum(operand1*np.transpose(operand2)))/(enssize)
            cov[i,j] = cov_tmp
            # Corr(xi,xj) = Cov(xi,xj)/sqrt(Var(xi)*Var(xj))
            var2 = np.var(ssh[:,time,j])
            corr_tmp = cov_tmp/np.sqrt(pointvar*var2)
            corr[i,j] = corr_tmp
    
    return corr,cov

#############
# corrAlong #
#############

def corrAlong(ssh):
    '''Function to calculate the correlation along one line (along track)
    
    Parameters:
    -----------
    ssh: 3d-array
    
    Output:
    -------
    corr,cov: 2d-arrays (dimension: width of trace x width of trace)'''
    # number of members
    enssize = ssh.shape[0]
    # length of trace
    dimtime = ssh.shape[1]
    ref_line = 0 # reference line
    # set dimensions of final correlation and covariance matrix
    corr = np.zeros([dimtime,dimtime])
    cov = np.zeros([dimtime,dimtime])
    # calculate correlation between each point along the reference line 
    for i in range (0,dimtime):
        # xi - E(xi)
        pointmean = np.nanmean(ssh[:,i,ref_line])
        pointvar = np.var(ssh[:,i,ref_line])
        operand1 = ssh[:,i,ref_line]-pointmean
        for j in range (0,dimtime):
            # xj - E(xj)
            operand2 = ssh[:,j,ref_line]-np.nanmean(ssh[:,j,ref_line])
            # Cov(xi,xj) = E[operand1*operand2^T]
            cov_tmp = (np.sum(operand1*np.transpose(operand2)))/(enssize)
            cov[i,j] = cov_tmp
            # Corr(xi,xj) = Cov(xi,xj)/sqrt(Var(xi)*Var(xj))
            var2 = np.var(ssh[:,j,ref_line])
            corr_tmp = cov_tmp/np.sqrt(pointvar*var2)
            corr[i,j] = corr_tmp
    
    return corr,cov


#############
# plot_corr #
#############

def plot_corr(data,cblabel):
    
    size = data.shape[0]
    step = np.arange(0,size,1,dtype=np.int)
    X, Y = np.meshgrid(step, step)
    Z = np.array(data)
    plt.contourf(X, Y, Z, 500,vmin=0.9,vmax=1., cmap='YlOrRd')
    plt.colorbar(label=cblabel)
    plt.gca().invert_yaxis() # y axis starts in upper left corner
    plt.axhline(color='grey', linestyle='dashed', linewidth=1.)
    plt.axvline(color='grey', linestyle='dashed', linewidth=1.)
    

############
# plot_cov #
############

def plot_cov_across(data,cblabel):
    
    level = np.arange(-0.01,0.011,0.00025)
    #level = np.arange(-0.0075,0.00775,0.00025)
    
    size = data.shape[0]
    step = np.arange(0,size,1,dtype=np.int)
    X, Y = np.meshgrid(step, step)
    Z = np.array(data)
    plt.contourf(X, Y, Z, level, cmap='seismic')
    plt.colorbar(label=cblabel)
    plt.gca().invert_yaxis() # y axis starts in upper left corner
    plt.axhline(color='grey', linestyle='dashed', linewidth=1.)
    plt.axvline(color='grey', linestyle='dashed', linewidth=1.)
    
    
##################
# plot_cov_along #
##################    

def plot_cov_along(data,cblabel):
    
    #level = np.arange(-0.008,0.00825,0.00025)
    
    size = data.shape[0]
    step = np.arange(0,size,1,dtype=np.int)
    X, Y = np.meshgrid(step, step)
    Z = np.array(data)
    plt.contourf(X, Y, Z, 50, cmap='YlOrRd')
    plt.colorbar(label=cblabel)
    plt.gca().invert_yaxis() # y axis starts in upper left corner
    plt.axhline(color='grey', linestyle='dashed', linewidth=1.)
    plt.axvline(color='grey', linestyle='dashed', linewidth=1.)
    
    
################
# gaspari_cohn #
################
    
def gaspari_cohn(r):
    """Gaspari-Cohn function. Used as weight in the cost function."""
    if type(r) is float:
        ra = np.array([r])
    else:
        ra = r
    ra = np.abs(ra)
    gp = np.zeros_like(ra)
    i=np.where(ra<=1.)[0]
    gp[i]=-0.25*ra[i]**5+0.5*ra[i]**4+0.625*ra[i]**3-5./3.*ra[i]**2+1.
    i=np.where((ra>1.)*(ra<=2.))[0]
    gp[i]=1./12.*ra[i]**5-0.5*ra[i]**4+0.625*ra[i]**3+5./3.*ra[i]**2-5.*ra[i]+4.-2./3./ra[i]
    if type(r) is float:
        gp = gp[0]

    return gp


############################
# gaspari_cohn_diag_matrix #
############################

def gaspari_cohn_diag_matrix(dim,x=np.arange(-30.,30.,1.)):
    """Generate a matrix with the Gaspari-Cohn function centered along its diagonal."""
    gc = gaspari_cohn(x/x[-1]*2)
    dimGC = gc.shape[0]
    weight = np.diag(np.ones(dim))
    for i in range(1, dimGC/2-1):#1 to 29 filled with values, 30 to end equal zero
            diag_pos = np.diag(np.ones(dim-i)*gc[dimGC/2+(i+1)],k=i) #starts at 31
            diag_neg = np.diag(np.ones(dim-i)*gc[dimGC/2-(i+1)],k=-i) #starts at 29
            weight += diag_pos+diag_neg
    return weight

#!/bin/ksh
#
set -x
#---------------------------------------------------------#
# Run first analysis to update ssh on swot trace
#---------------------------------------------------------#
# Prepare localization (modify localization.cfg)
rm -f L10zone.czon
rm -f L10partition_*
sesam -mode zone -incfg L10localization.cfg -outpartvar L10partition_#.nc -outzon L10zone.czon

# Run analysis with localization
rm -fr L10UPENS0060.nc.bas
mkdir L10UPENS0060.nc.bas
rm xa_ana1_grid*.nc

sesam -mode lroa -outvar xa_ana1_#.nc -invarref zero_#.nc -inxbas ENS0060.nc.bas -outxbas L10UPENS0060.nc.bas -indta rmdc_obs_p182_SWOT0004_denoised_#.ncdta -oestd R_#.nc -inpartvar L10partition_#.nc -inzon L10zone.czon -disable TTTFT


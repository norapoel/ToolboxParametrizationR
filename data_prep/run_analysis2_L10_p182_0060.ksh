#!/bin/ksh
#

sesam -mode oper -invar mask_#.nc -outvar stdR_0.0001_#.nc -typeoper cst_0.0001

mv xa_ana1_gridSWOT.nc xa_ana1_gridSWOT.ncdbs

mv mask_gridT.nc mask.nc

ncap2 -s "lon=lon-360.0" xa_ana1_gridSWOT.ncdbs xa_ana1_gridSWOT.ncdbs

sesam -mode obsv -indbs xa_ana1_gridSWOT.ncdbs -outobs obs_ana1_#.cobs -affectobs SWOT

rm L10zone.czon L10partition*
sesam -mode zone -incfg L10localization.cfg -outpartvar L10partition_#.nc -outzon L10zone.czon

rm -fr DERUPENS0060.nc.bas ; mkdir DERUPENS0060.nc.bas
rm -f xa_ana2_gridT.nc

sesam -mode lroa -inxbas ENS0060.nc.bas -invarref cst_0.00_#.nc -inobs obs_ana1_#.cobs -configobs obs_ana1_#.cobs -outxbas DERUPENS0060.nc.bas -outvar xa_ana2_#.nc -oestd stdR_0.0001_#.nc -inpartvar L10partition_#.nc -inzon L10zone.czon -disable TTTFT

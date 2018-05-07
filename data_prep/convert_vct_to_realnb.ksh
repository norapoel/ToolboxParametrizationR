#!/bin/ksh
#

dir="/home/poeln/Work/obs_pass182/OBS0060"

for k in $(seq 1 60); do
  
  tag=$(printf "%04d" $k)

  # conversion to real numbers
  ncap -O -s "ssh_obs=(ssh_obs*1.0);ssh_model=(ssh_model*1.0);lon=(lon*1.0);lat=(lat*1.0)" $dir/vctgridSWOT$tag.nc $dir/vctgridSWOT$tag.nc
  
done


#!/bin/ksh
#
set -x
outdir="/home/poeln/Work/swotsimulator/output/p182_from_xt_0060"

let imem=1
while [ $imem -le 8 ] ; do
  swotsimulator params_pass182_withgap.py

  tag=`echo $imem | awk '{printf("%04d", $1)}'`
  mv $outdir/swot292-OSMOSIS_c01_p182.nc $outdir/obs_p182_SWOT$tag.nc
  let imem=$imem+1
done


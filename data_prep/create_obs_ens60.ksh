#!/bin/ksh
#
set -x
indir="/home/poeln/Work/swotsimulator/input/ENS0060.nc.bas"
outdir="/home/poeln/Work/swotsimulator/output/OBS0060"


for k in $(seq 2 60); do
  ### change list_of_file.txt
  tag=$(printf "%04d" $k)
  
  cp $indir/list_tmp.txt $indir/list_of_file.txt

  sed -i "s/NB/$tag/g" $indir/list_of_file.txt
  
  ### run swot simulator
  swotsimulator params_pass182_nogap.py

  ### change name of output
  cp $outdir/swot292-OSMOSIS_c01_p182.nc $outdir/vctgridSWOT$tag.nc
done


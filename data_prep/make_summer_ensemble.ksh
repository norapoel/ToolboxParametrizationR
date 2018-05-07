#!/bin/bash
#

outdir="ENS0060.nc.bas"
let day1=1
let imem=1
while [ $imem -le 10 ] ; do 
  tag=`echo $imem | awk '{printf("%04d", $1)}'`
  echo $tag $day1
  ncks -d time_counter,$day1,$day1 alldata.nc $outdir/vctgridT$tag.nc
  let days=$day1+2
  let imem=$imem+1
done 

let day2=375
while [ $imem -le 60 ] ; do
  tag=`echo $imem | awk '{printf("%04d", $1)}'`
  echo $tag $day2
  ncks -d time_counter,$day2,$day2 alldata.nc $outdir/vctgridT$tag.nc
  let imem=$imem+1
  let day2=$day2+2
done

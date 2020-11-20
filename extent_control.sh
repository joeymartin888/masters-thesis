#!/bin/bash

cd ~/Data/CanCM4_control
for year in `seq 2300 2449`; do
    filein=sc_dhfp1e_e001_${year}_m01_${year}_m12_sicn.nc
    fileout_sie=sc_dhfp1e_e001_${year}_SIE.nc
#	    fileout_sia=SIA_monthly_${model}_i${year}${mon}.nc
    echo $fileout_sie
#	    echo $fileout_sia
    ncks -v sicn $filein tmp.nc  #extracts SIC
    cdo -gtc,0.15 tmp.nc tmp2.nc 
    cdo  -fldmean -sellonlatbox,0,360,0,90 tmp2.nc SIE/$fileout_sie 
#	    cdo  -fldmean -sellonlatbox,0,360,0,90 tmp.nc  SIA/$fileout_sia >& /dev/null

    /bin/rm tmp.nc tmp2.nc
done
       
    

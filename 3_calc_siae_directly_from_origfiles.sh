#!/bin/bash

#modellist='CanCM3 CanCM4'
modellist='CanCM3 CanCM4'

for model in $modellist; do
    cd /HOME/rms/DATA/CanSIPS/monthly/
    for year in `seq 1979 2010`; do
	for mon in `seq 1 12`; do
	    mon=`echo $mon|awk -F':' '{printf "%2.2d",$1}' -`
            filein=sic/orig/sic_monthly_CCCma-${model}_CHFP_original_grid_i${year}${mon}01.nc
	    fileout_sie=SIE_monthly_${model}_i${year}${mon}.nc
#	    fileout_sia=SIA_monthly_${model}_i${year}${mon}.nc
	    echo $fileout_sie
#	    echo $fileout_sia

	    ncks -v sic $filein tmp.nc
	    cdo -gtc,0.15 tmp.nc tmp2.nc
	    cdo  -fldmean -sellonlatbox,0,360,0,90 tmp2.nc SIE/$fileout_sie 
#	    cdo  -fldmean -sellonlatbox,0,360,0,90 tmp.nc  SIA/$fileout_sia >& /dev/null

	    /bin/rm tmp.nc tmp2.nc
	done
    done
done        
    

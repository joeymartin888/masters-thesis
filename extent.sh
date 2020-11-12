#!/bin/bash

#modellist='CanCM3 CanCM4'
modellist='CanCM4'
#modellist='GEM-NEMO'

for model in $modellist; do
    echo $model
    cd Data/GEM-NEMO
    for year in `seq 1979 2019`; do
	for mon in `seq 1 12`; do
	    mon=`echo $mon|awk -F':' '{printf "%2.2d",$1}' -`
            if [[ $model == 'GEM-NEMO' ]]; then
	         echo $model
		 filein=sic_monthly_CCCma-GEM_NEMO_1x1_grid_i${year}${mon}01.nc
	    else
		 if [ $year -gt 2012 ]; then
			filein=sic_monthly_CCCma-${model}_NEW_OPER_original_grid_i${year}${mon}01.nc
	   	 else
	    		filein=sic_monthly_CCCma-${model}_NEW_original_grid_i${year}${mon}01.nc
	    	fi
	    fi
	    fileout_sie=SIE_monthly_${model}_i${year}${mon}.nc
#	    fileout_sia=SIA_monthly_${model}_i${year}${mon}.nc
	    echo $fileout_sie
#	    echo $fileout_sia

	    ncks -v sic $filein tmp.nc  #extracts SIC
	    cdo -gtc,0.15 tmp.nc tmp2.nc 
	    cdo  -fldmean -sellonlatbox,0,360,0,90 tmp2.nc SIE/$fileout_sie 
#	    cdo  -fldmean -sellonlatbox,0,360,0,90 tmp.nc  SIA/$fileout_sia >& /dev/null

	    /bin/rm tmp.nc tmp2.nc
	done
    done
done        
    

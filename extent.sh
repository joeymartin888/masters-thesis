#!/bin/bash

#modellist='CanCM3 CanCM4'
#modellist='CanCM3 CanCM4'
#modellist='GEM-NEMO'
modellist='PPrun'
years=(2312 2435 2413 2356 2387 2332)

for model in $modellist; do
    cd /home/josmarti/Data/$model
    for year in ${years[@]}; do
	if [ $model='PPrun' ]; then
		for e in `seq -w 001 012`; do
			for i in `seq 0 2`; do
				filein=sc_dhfp1e_e${e}_i${year}_m01_$((year+i))_m01_$((year+i))_m12_sicn.nc
				fileout_sie=SIE_${model}_e${e}_i${year}_$((year+i)).nc
				echo $fileout_sie

				ncks -v sicn $filein tmp.nc  #extracts SIC
				cdo -gtc,0.15 tmp.nc tmp2.nc 
				cdo  -fldmean -sellonlatbox,0,360,0,90 tmp2.nc SIE/$fileout_sie 
#				cdo  -fldmean -sellonlatbox,0,360,0,90 tmp.nc  SIA/$fileout_sia >& /dev/null
				
				/bin/rm tmp.nc tmp2.nc
			done
		done
	else
		for mon in `seq 1 12`; do
		    mon=`echo $mon|awk -F':' '{printf "%2.2d",$1}' -`
		    if [ $model='GEM-NEMO' ]; then
			 filein=sic_monthly_CCCma-GEM_NEMO_1x1_grid_i${year}${mon}01.nc
		    else
			 if [ $model='CanCM3' ] && [ $year -gt 2012 ]; then
				filein=sic_monthly_CCCma-${model}_NEW_OPER_original_grid_i${year}${mon}01.nc
		   	 else
		    		filein=sic_monthly_CCCma-${model}_NEW_original_grid_i${year}${mon}01.nc
		    	fi
		    fi
	#	    fileout_sie=SIE_monthly_${model}_i${year}${mon}.nc
		    fileout_sia=SIA_monthly_${model}_i${year}${mon}.nc
	#	    echo $fileout_sie
		    echo $fileout_sia
	
		    ncks -v sic $filein tmp.nc  #extracts SIC
		    cdo -gtc,0.15 tmp.nc tmp2.nc 
	#	    cdo  -fldmean -sellonlatbox,0,360,0,90 tmp2.nc SIE/$fileout_sie 
		    cdo  -fldmean -sellonlatbox,0,360,0,90 tmp.nc  SIA/$fileout_sia >& /dev/null
	
	    /bin/rm tmp.nc tmp2.nc	
		done
	fi    
      done
done        
    

#!/bin/bash

#modellist='CanCM3 CanCM4'
#modellist='CanCM3 CanCM4'
#modellist='GEM-NEMO'
modellist='PP_run_full'
years=(2312 2435 2413 2356 2387 2332)


for model in $modellist; do
    cd /home/josmarti/Data/$model
    for year in ${years[@]}; do
	for month in `seq -w 01 2 12`; do
		if [ $model = 'PP_run_full' ]; then
			for e in `seq -w 001 012`; do
				if [ $month = 01 ]; then
					end=12
					for i in `seq 0 2`; do
						filein=sc_dhfp1e_e${e}_i${year}_m${month}_$((year+i))_m${month}_$((year+i))_m${end}_sicn.nc
					
						ncks -v sicn $filein tmp3.nc  #extracts SIC
						cdo -gtc,0.15 tmp3.nc tmp4.nc 
						cdo  -fldmean -sellonlatbox,0,360,0,90 tmp4.nc tmpy${i}.nc 
		#				cdo  -fldmean -sellonlatbox,0,360,0,90 tmp.nc  SIA/$fileout_sia >& /dev/null
						
						/bin/rm tmp3.nc tmp4.nc 
					done
				
				fileout_sie=SIE_${model}_e${e}_i${year}${month}_$((year+2))${end}.nc
						
				echo $fileout_sie
				
				cdo cat tmpy*.nc SIE/$fileout_sie

				/bin/rm tmpy*.nc	
				else
					end=$(printf "%02d" $((10#$month-1)))	
					filein=sc_dhfp1e_e${e}_i${year}_m${month}_${year}_m${month}_$((year+3))_m${end}_sicn.nc
					fileout_sie=SIE_${model}_e${e}_i${year}${month}_$((year+3))${end}.nc
									
					echo $fileout_sie

					ncks -v sicn $filein tmp.nc  #extracts SIC
					cdo -gtc,0.15 tmp.nc tmp2.nc 
					cdo  -fldmean -sellonlatbox,0,360,0,90 tmp2.nc SIE/$fileout_sie 
	#				cdo  -fldmean -sellonlatbox,0,360,0,90 tmp.nc  SIA/$fileout_sia >& /dev/null
					
					/bin/rm tmp.nc tmp2.nc
				fi
			done
		else
			for mon in `seq 1 12`; do
			    mon=`echo $mon|awk -F':' '{printf "%2.2d",$1}' -`
			    if [ $model = 'GEM-NEMO' ]; then
				 filein=sic_monthly_CCCma-GEM_NEMO_1x1_grid_i${year}${mon}01.nc
			    else
				 if [ $model = 'CanCM3' ] && [ $year -gt 2012 ]; then
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
done    
    

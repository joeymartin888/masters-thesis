#!/bin/bash

regions=`seq 1 14`

#Select OLD or NEW
#echo "Which version? "
#read versions
#versions='OLD NEW'
#versions="GEM-NEMO"

for reg in $regions; do 
cdo -setmisstoc,0 -eqc,${reg} ~/Data/1x1_reg_mask.nc ~/Data/Observations/tmpmask.nc

filein=~/Data/Observations/had2cis_1x1_198001_202004_sicn.nc

echo $reg
            
	    fileout_sie=~/Data/Observations/Observed_SIE_${reg}.nc
#	    fileout_sia=SIA_monthly_${model}_i${year}${mon}.nc
	    echo $fileout_sie
#	    echo $fileout_sia

	    cdo remapbil,~/Data/1x1_reg_mask.nc $filein tmp2.nc
	    cdo -setmisstoc,1 tmp2.nc tmp2fixed.nc
	    ncks -v LSMASK ~/Data/lsmask_cansipsv2_sea.nc lsmask_cansipsv2_maskonly.nc
	    cdo mul tmp2fixed.nc lsmask_cansipsv2_maskonly.nc tmpmulmask.nc
	    cdo -gtc,0.15 tmpmulmask.nc tmp3.nc  	    
	    cdo mul tmp3.nc ~/Data/Observations/tmpmask.nc tmpmul.nc
	    #ncks -v SICN tmpmul.nc tmp.nc  #extracts SIC 
	    cdo  -fldmean -sellonlatbox,0,360,0,90 tmpmul.nc $fileout_sie 
#	    cdo  -fldmean -sellonlatbox,0,360,0,90 tmp.nc  SIA/$fileout_sia >& /dev/null

	    /bin/rm tmp.nc tmp2.nc tmpmul.nc tmp3.nc tmpmulmask.nc lsmask_cansipsv2_maskonly.nc tmp2fixed.nc
  
rm ~/Data/Observations/tmpmask.nc

done   
    

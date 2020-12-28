#!/bin/bash

echo "Which region?"
read region
echo "r?"
read reg

cdo -setmisstoc,0 -eqc,${reg} ~/Data/1x1_reg_mask.nc ~/Data/Observations/tmpmask.nc

cd ~/Data/CanCM4_control
for year in `seq 2300 2449`; do
    filein=sc_dhfp1e_e001_${year}_m01_${year}_m12_sicn.nc
    fileout_sie=sc_dhfp1e_e001_${region}_${year}_SIE.nc
#	    fileout_sia=SIA_monthly_${model}_i${year}${mon}.nc
    echo $fileout_sie
#	    echo $fileout_sia

    cdo remapbil,~/Data/1x1_reg_mask.nc $filein tmp3.nc
    ncks -v sicn tmp3.nc tmp.nc  #extracts SIC
    cdo -gtc,0.15 tmp.nc tmp2.nc 
    cdo mul tmp2.nc ~/Data/Observations/tmpmask.nc tmpmul.nc
    cdo  -fldmean -sellonlatbox,0,360,0,90 tmpmul.nc SIE/${region}/$fileout_sie 
#	    cdo  -fldmean -sellonlatbox,0,360,0,90 tmp.nc  SIA/$fileout_sia >& /dev/null

    /bin/rm tmp.nc tmp2.nc tmp3.nc tmpmul.nc
done

rm ~/Data/Observations/tmpmask.nc  
       
    

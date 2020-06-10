
#!/bin/bash

region='5GRE'
echo $region
r=5
echo $r

#Select OLD or NEW
echo "Which version? "
read version
ncap2 -s 'where(REG!=5) REG=0;' Data/iceregions_128_64.nc -O Data/${version}/tmpmask2div.nc
cdo divc,5 Data/${version}/tmpmask2div.nc Data/${version}/tmpmask.nc
rm Data/${version}/tmpmask2div.nc

#modellist='CanCM3 CanCM4'
if [ $version = "GEM-NEMO" ]
then
modellist='GEM_NEMO'
years=`seq 1980 2010`
else
modellist='CanCM3 CanCM4'
years=`seq 1979 2010`
fi


for model in $modellist; do
    cd Data/${version}
    if [ ! -d "SIA/$region" ]
    then
    	mkdir SIA/$region
    fi
    for year in $years; do
	for mon in `seq 1 12`; do
	    mon=`echo $mon|awk -F':' '{printf "%2.2d",$1}' -`
            if [ $version = "OLD" ]
	    then
	    	filein=sic_monthly_CCCma-${model}_CHFP_original_grid_i${year}${mon}01.nc
            elif [ $version = "GEM-NEMO" ]
	    then
	    	filein=sic_monthly_CCCma-${model}_1x1_grid_i${year}${mon}01.nc
	    elif [ $version = "NEW" ]
	    then
	    filein=sic_monthly_CCCma-${model}_NEW_original_grid_i${year}${mon}01.nc
	    fi
#	    fileout_sie=SIE_monthly_${model}_i${year}${mon}.nc
	    fileout_sia=SIA_monthly_${region}_${model}_i${year}${mon}.nc
#	    echo $fileout_sie
	    echo $fileout_sia

	    ncks -v sic $filein tmp.nc  #extracts SIC
	    tcat=()
	    ecat=()
  	    for e in $(seq 0 9)
	    do
	    	for t in $(seq 0 11)
        	do
			ncks -v sic -d time,$t -d ensemble,$e tmp.nc tmp1.nc
               		cdo -invertlat tmp1.nc tmplat.nc
                	cdo mul tmplat.nc tmpmask.nc tmpt${t}.nc
                	tcat+="tmpt${t}.nc "
                	rm tmp1.nc tmplat.nc
		done
		cdo cat $tcat tmptoute${e}.nc
		rm $tcat
		tcat=()
		ecat+="tmptoute${e}.nc "
	    done
	    cdo merge $ecat tmpout.nc
	    rm $ecat
	    cdo -gtc,0.15 tmpout.nc tmp2.nc 
#	    cdo  -fldmean -sellonlatbox,0,360,0,90 tmp2.nc SIE/$fileout_sie 
	    cdo  -fldmean -sellonlatbox,0,360,0,90 tmpout.nc  SIA/$region/$fileout_sia >& /dev/null

	    /bin/rm tmp.nc tmp2.nc tmpout.nc
	done
    done
done        
    

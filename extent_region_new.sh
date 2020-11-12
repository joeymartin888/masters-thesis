
#!/bin/bash

echo "Which region? "
read region
echo "r? " 
read r

#Select OLD or NEW
#echo "Which version? "
#read versions
versions='OLD NEW GEM-NEMO'


for version in $versions; do 
cdo -setctomiss,0 -eqc,${r} ~/1x1_reg_mask.nc ~/Data/${version}/tmpmask.nc

#modellist='CanCM3 CanCM4'
if [ $version = "GEM-NEMO" ]
then
modellist='GEM_NEMO'
years=`seq 1980 2018`
else
modellist='CanCM3 CanCM4'
years=`seq 1979 2018`
fi

for model in $modellist; do
    if [ $version='GEM-NEMO' ]
    then
	cd ~/Data/${version}
    else
	cd ~/Data/${version}/1x1
    fi
    if [ ! -d "/home/josmarti/Data/${version}/1x1/SIE/${region}" ]
    then
    	mkdir ~/Data/${version}/1x1/SIE/${region}
    fi
    for year in $years; do
	for mon in `seq 1 12`; do
	    mon=`echo $mon|awk -F':' '{printf "%2.2d",$1}' -`
            if [ $version = "OLD" ]
	    then
	    	filein=~/Data/OLD/1x1/sic_monthly_CCCma-${model}_NEW_1x1_grid_i${year}${mon}01.nc
            elif [ $version = "GEM-NEMO" ]
	    then
	    	filein=~/Data/GEM-NEMO/1x1/sic_monthly_CCCma-${model}_1x1_grid_i${year}${mon}01.nc
	    elif [ $version = "NEW" ]
	    then
	        filein=~/Data/NEW/1x1/sic_monthly_CCCma-${model}_NEW_1x1_grid_i${year}${mon}01.nc
	    fi
	    fileout_sie=~/Data/${version}/1x1/SIE/${region}/SIE_monthly_${region}_${model}_i${year}${mon}.nc
#	    fileout_sia=SIA_monthly_${region}_${model}_i${year}${mon}.nc
	    echo $fileout_sie
#	    echo $fileout_sia

	    ncks -v sic $filein tmp.nc  #extracts SIC
	    tcat=()
	    ecat=()
  	    for e in $(seq 0 9)
	    do
	    	for t in $(seq 0 11)
        	do
			if [ $version = "GEM-NEMO" ]
			then
				ncks -v sic -d time,$t -d ensemble,$e tmp.nc tmplat.nc
				cdo -invertlat tmplat.nc tmp1.nc
				rm tmplat.nc
			else
				ncks -v sic -d time_counter,$t -d ensemble,$e tmp.nc tmp1.nc
               		fi
			cdo -gtc,0.15 tmp1.nc tmp2.nc 
			cdo mul tmp2.nc tmpmask.nc tmpout.nc
			cdo  -fldmean -sellonlatbox,0,360,0,90 tmpout.nc tmpt${t}.nc
                	tcat+="tmpt${t}.nc "
                	rm tmp1.nc tmp2.nc tmpout.nc
		done
		cdo cat $tcat tmptoute${e}.nc
		rm $tcat
		tcat=()
		ecat+="tmptoute${e}.nc "
	    done
	    cdo merge $ecat $fileout_sie
	    rm $ecat  
#	    cdo  -fldmean -sellonlatbox,0,360,0,90 tmpout.nc  SIA/$region/$fileout_sia >& /dev/null

	    /bin/rm tmp.nc 
	done
    done
done
cd         
rm ~/Data/${version}/tmpmask.nc
done   

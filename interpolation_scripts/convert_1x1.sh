years="2020"
ms="08"

for m in ${ms} ; do
for y in ${years} ; do

cp sic_daily_CCCma-CanCM4_NEW_original_grid_i${y}${m}01.nc input.nc

for e in 1 2 3 4 5 6 7 8 9 10 ; do

ncks -h -F -d ensemble,${e} input.nc tmp.nc
ncwa -h -a ensemble tmp.nc tmp1.nc
ncwa -h -a time_bnd tmp1.nc tmp2.nc

# mask out
cdo ifthen maskn.nc tmp2.nc tmp3.nc
ncatted -a coordinates,sic,c,c,"nav_lon nav_lat" tmp3.nc

cdo fillmiss2,4 tmp3.nc tmp4.nc
cdo remapbil,mygrid_1deg tmp4.nc tmp5.nc

cdo ifthen mask1x1n.nc tmp5.nc out.nc
ncrename -d reftime,time_counter -v reftime,time_counter out.nc

mv out.nc sic_daily_CCCma-CanCM4_NEW_1x1_grid_i${y}${m}01_${e}.nc
rm tmp*.nc

done
rm input.nc

done
done

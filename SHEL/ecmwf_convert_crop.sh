#!/bin/bash

export HOMEPATH=/data1/2020/KMA/SFR-07
export SHELPATH=${HOMEPATH}/SHEL
export DAINPATH=${HOMEPATH}/DAIN
export DAIOPATH=${HOMEPATH}/DAIO
export DASVPATH=${HOMEPATH}/DASV

export GRIB_TO_NETCDF=/home/shpark/Downloads/grib_api/grib_api-1.28.0-Source/tools/grib_to_netcdf 
export GRIB_COPY=/home/shpark/Downloads/grib_api/grib_api-1.28.0-Source/tools/grib_copy
#srtdate=`date +"%Y%m%d%H"`
srtdate=2020040800
SYY=${srtdate:0:4}
SMM=${srtdate:4:2}
SDD=${srtdate:6:2}
SHH=${srtdate:8:2}

if [ $SHH -lt 12 ]; then
   SHH=`printf "%02d" 00`
else
   SHH=12
fi
#ECMWF Enseble Model data merge
cat ${HOMEPATH}/202004/0800/ecmw_ewam_acd2_*.${srtdate}00.gb1 > ${DAIOPATH}/ecmw_ewam_total_tmp.${srtdate}00.gb1
#ECMWF Enseble Model variable extraction
${GRIB_COPY} -w shortName=swh ${DAIOPATH}/ecmw_ewam_total_tmp.${srtdate}00.gb1 ${DAIOPATH}/swh_tmp.gb1
${GRIB_COPY} -w shortName=mwd ${DAIOPATH}/ecmw_ewam_total_tmp.${srtdate}00.gb1 ${DAIOPATH}/mwd_tmp.gb1
${GRIB_COPY} -w shortName=mwp ${DAIOPATH}/ecmw_ewam_total_tmp.${srtdate}00.gb1 ${DAIOPATH}/mwp_tmp.gb1
${GRIB_COPY} -w shortName=wind ${DAIOPATH}/ecmw_ewam_total_tmp.${srtdate}00.gb1 ${DAIOPATH}/wind_tmp.gb1
${GRIB_COPY} -w shortName=dwi ${DAIOPATH}/ecmw_ewam_total_tmp.${srtdate}00.gb1 ${DAIOPATH}/dwi_tmp.gb1
rm -f ${DAIOPATH}/ecmw_ewam_total_tmp.${srtdate}00.gb1
cat ${DAIOPATH}/*_tmp.gb1 > ${DAIOPATH}/ecmw_ewam_total.${srtdate}00.gb1
rm -f ${DAIOPATH}/*_tmp.gb1
#ECMWF Enseble netcdf convert
$GRIB_TO_NETCDF -o ${DAIOPATH}/ecmwf_ense_tmp.nc ${DAIOPATH}/ecmw_ewam_total.202004080000.gb1
#ECMWF Enseble cropping
cdo -sellonlatbox,120,142,20,50 ${DAIOPATH}/ecmwf_ense_tmp.nc ${DAIOPATH}/ecmwf_ense.nc
exit
#ECMWF olny Model data merge
#cat ${HOMEPATH}/${SYY}${SMM}/${SDD}${SHH}/ecmw_hwam_glob_h*.${srtdate}00.gb1 > ${DAIOPATH}/ecmw_hwam_total_tmp.${srtdate}00.gb1
#${GRIB_COPY} -w shortName=swh ${DAIOPATH}/ecmw_hwam_total_tmp.${srtdate}00.gb1 ${DAIOPATH}/swh_tmp.gb1
#${GRIB_COPY} -w shortName=mwd ${DAIOPATH}/ecmw_hwam_total_tmp.${srtdate}00.gb1 ${DAIOPATH}/mwd_tmp.gb1
#${GRIB_COPY} -w shortName=mwp ${DAIOPATH}/ecmw_hwam_total_tmp.${srtdate}00.gb1 ${DAIOPATH}/mwp_tmp.gb1
#${GRIB_COPY} -w shortName=wind ${DAIOPATH}/ecmw_hwam_total_tmp.${srtdate}00.gb1 ${DAIOPATH}/wind_tmp.gb1
#${GRIB_COPY} -w shortName=dwi ${DAIOPATH}/ecmw_hwam_total_tmp.${srtdate}00.gb1 ${DAIOPATH}/dwi_tmp.gb1
#rm -f ${DAIOPATH}/ecmw_hwam_total_tmp.${srtdate}00.gb1
#cat ${DAIOPATH}/*_tmp.gb1 > ${DAIOPATH}/ecmw_hwam_total.${srtdate}00.gb1

#ECMWF only netcdf convert
#$GRIB_TO_NETCDF -o ${DAIOPATH}/ecmwf_only_tmp.nc ${DAIOPATH}/ecmw_hwam_total.${srtdate}00.gb1
$GRIB_TO_NETCDF -o ${DAIOPATH}/ecmwf_only_tmp.nc ${HOMEPATH}/20200309/ecmw_hwam_glob_h000.202003090000.gb1

#ECMWF only cropping
cdo -sellonlatbox,120,142,20,50 ${DAIOPATH}/ecmwf_only_tmp.nc ${DAIOPATH}/ecmwf_only_total.nc

#done

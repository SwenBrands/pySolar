#!/bin/bash

source ${HOME}/.bashrc

#check python version
echo "Your Python version is:"
python --version
#launch to queue on lustre
RUNDIR=${HOME}/datos/tareas/proyectos/pticlima/pyPTIclima/radiation
STOREDIR=${HOME}/datos/tareas/proyectos/pticlima/radiation/netcdf
FILENAME=rsds_day_aemet_20171201_20221231.nc
FIGDIR=${HOME}/datos/tareas/proyectos/pticlima/radiation/figs
FIGNAME=rsds_day_aemet_coverage_20171201_20221231.pdf
REMOTEDIR=/lustre/gmeteo/WORK/DATA/PTI-Clima/OBSERVATIONS/AEMET-Stations_PTI
cd ${RUNDIR}

#convert csv files from AEMET to netCDF file named ${FILENAME} located in ${STOREDIR}
echo "Invoking csv2nc.py with bash..."
python csv2nc.py > logfile_csv2nc.log

#sent local file to remote directory
echo "Sending ${STOREDIR}/${FILENAME} to swen@ui.sci.unican.es:${REMOTEDIR}..."
rsync -av ${STOREDIR}/${FILENAME} swen@ui.sci.unican.es:${REMOTEDIR}

#crop pdf file
echo "Cropping ${FIGDIR}/${FIGNAME}..."
pdfcrop ${FIGDIR}/${FIGNAME} ${FIGDIR}/${FIGNAME}
sleep 2

#sent to ui
echo "Sending ${FIGDIR}/${FIGNAME} to swen@ui.sci.unican.es:${REMOTEDIR}..."
rsync -av ${FIGDIR}/${FIGNAME} swen@ui.sci.unican.es:${REMOTEDIR}

#and exit
echo "launchme.sh has run successfully, exiting now..."
exit 0

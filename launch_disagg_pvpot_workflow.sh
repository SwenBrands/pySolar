#!/bin/bash

if [ "$1" == "-h" ]; then
  echo "Description: `basename $0` This script launches the entire workflow beginning with the disaggregation of ERA5-Land flux data from CDS, followed by the PVpot index caculation and upload of the result files to lustre via the rsync command."
  exit 0
fi

source ${HOME}/.bashrc
#check python version
echo "Your Python version is:"
python --version

#set path to the run directory as well as to the local and remote directories containing the netCDF files to be transferred
RUNDIR=${HOME}/datos/tareas/proyectos/pticlima/pyPTIclima/pySolar
LOGDIR=${RUNDIR}/LOG
STOREDIR=${HOME}/datos/OBSData/era5_land_disagg
REMOTEDIR1=/lustre/gmeteo/PTICLIMA/DATA/REANALYSIS/ERA5-Land/data
REMOTEDIR2=/lustre/gmeteo/PTICLIMA/DATA/REANALYSIS/ERA5-Land/data_derived

#user options
cleandirs="no" #delete the variable output directories
runpython="no" #run the 3 distinct Python scripts
syncronize="yes" #upload results to lustre

## EXECUTE #############################################################
cd ${RUNDIR}

#optionally delete previously created result folders
if [ ${cleandirs} == "yes" ]
then
    echo "Upon user request, the variable directories located at ${STOREDIR} are deleted"
    for domain in Iberia Canarias
    do
        for aggreg in hour day
        do
            for variable in tp ssrd pvpot
            do
                #define source directory
                sourcedir=${STOREDIR}/${domain}/${aggreg}/${variable}
                echo "Deleting directory ${sourcedir} ...."
                rm -r ${sourcedir}
            done
        done
    done
else
    echo "Upon user request, the variable directories located at ${STOREDIR} are not deleted"
fi

##launch the relevant Python scripts sequentially

if [ ${runpython} == "yes" ]
then
    echo "Upon user request, the Python scripts will be launched."
    
    echo "launching disagg.py..."
    python disagg.py > ${LOGDIR}/logfile_disagg.log
    echo "disagg.py has terminated !"

    echo "launching pvpot_calculator.py..."
    python pvpot_calculator.py > ${LOGDIR}/logfile_pvpot_calculator.log
    echo "pvpot_calculator.py has terminated !"

    echo "launching plot_timeseries_aspects.py..."
    python plot_timeseries_aspects.py > ${LOGDIR}/logfile_plot_timeseries_aspects.log
    echo "plot_timeseries_aspects.py has terminated as well!"
else
    echo "Upon user request, the Python scripts will not be launched."
fi

#optionally loop through domain, temporal aggregation and variable and upload to lustre
if [ ${syncronize} == "yes" ]
then
    for domain in Canarias Iberia
    do
        for aggreg in day hour
        do
            for variable in ssrd tp pvpot
            do
                #define source directory
                sourcedir=${STOREDIR}/${domain}/${aggreg}/${variable}
                ##set directory on remote server depending on the considered variable
                if [ ${variable} == ssrd ] || [ ${variable} == tp ]
                then
                    echo loop1
                    echo variable is ${variable}
                    tardir=swen@ui.sci.unican.es:${REMOTEDIR1}/${domain}/${aggreg}/${variable}_pySolar
                elif [ ${variable} == pvpot ]
                then
                    echo loop2
                    echo variable is ${variable}
                    tardir=swen@ui.sci.unican.es:${REMOTEDIR2}/${domain}/${aggreg}/${variable}
                else
                    echo "ERROR: Unknown entry for <variable>, exiting now !"
                    exit 1
                fi
                echo "Source directory: ${sourcedir}"
                echo "Remote destination directory: ${tardir}"
                chmod 777 -R ${sourcedir}
                rsync -av --delete ${sourcedir}/ ${tardir}
                #rsync -av ${sourcedir}/ ${tardir}
                sleep 2
            done
        done
    done
else
    echo "Upon user request, the output variable directories located at ${STOREDIR} are not uploaded to the remote server."
fi

echo "nc_to_lustre.sh has run successfully, exiting now..."
exit 0

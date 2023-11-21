#!/bin/bash

source ${HOME}/.bashrc

#check python version
echo "Your Python version is:"
python --version
#launch to queue on lustre
RUNDIR=${HOME}/datos/tareas/proyectos/pticlima/pyPTIclima/pySolar
STOREDIR=${HOME}/datos/OBSData/era5_land_disagg
REMOTEDIR1=/lustre/gmeteo/PTICLIMA/DATA/REANALYSIS/ERA5-Land/data
REMOTEDIR2=/lustre/gmeteo/PTICLIMA/DATA/REANALYSIS/ERA5-Land/data_derived

## EXECUTE #############################################################
cd ${RUNDIR}
#loop through domain, temporal aggregation and variable
for domain in Iberia Canarias
do
    for aggreg in hour
    do
        for variable in ssrd pvpot
        do
            #copy source directory
            sourcedir=${STOREDIR}/${domain}/${aggreg}/${variable}
            #copydir=${STOREDIR}/${domain}/${aggreg}/${variable}_pySolar
            #echo "copying ${sourcedir} to ${copydir}"
            #cp -r ${sourcedir} ${copydir}
            
            ##set directory on remote server depending on the considered variable
            if [ ${variable}==ssrd ]
            then
                tardir=swen@ui.sci.unican.es:${REMOTEDIR1}/${domain}/${aggreg}/${variable}_pySolar
            elif [ ${variable}==pvpot ]
            then
                tardir=swen@ui.sci.unican.es:${REMOTEDIR2}/${domain}/${aggreg}/${variable}
            else
                echo "ERROR: Unknown entry for <variable>, exiting now !"
                exit 1
            fi
            echo ${copydir}
            echo ${sourcedir}
            echo ${tardir}
            rsync -I -av --delete ${sourcedir} ${tardir}
            sleep 2
        done
    done
done

echo "nc_to_lustre.sh has run successfully, exiting now..."
exit 0

#!/bin/bash  --
exec 2> Errors_cleanExp.txt
###############################################

#  Function:  Run a set of experiments from a pre-defined directory tree

###############################################

## Current and run directory
cdir=$(pwd)
rundir=$1

#Check for, and remove, trailing slashes
inpL=${#rundir}
if test "${rundir:$inpL-1:1}" == "/" ; then
    rundir=${rundir:0:$(($inpL-1))}
fi

#Directory containing runs
workdir=$cdir/$rundir
if ! cd $workdir; then
    echo 'Invalid Run Directory:' $rundir 
    exit
elif [ "$rundir" == '' ]; then
    echo 'Run Directory must be specified'
    exit
fi

echo 'Cleaning the' $rundir 'directory'

cd $workdir
rm *.log
for dir in `ls -d ./*/`; do
    echo 'Moving in' $dir
    cd $dir
    rm OAR.MF1DS_*.std*
    rm MF1D.out
    # Check for the output directories
    if [ -d /data/failles/$USER/$rundir/${dir:2:-1} ]; then
        rm -r /data/failles/$USER/$rundir/${dir:2:-1}
    fi
    if [ `ls -ld Figs | wc -l` == 1 ]; then
        rm -r Figs
    fi
    if [ -d database ]; then
        #ls database
        rm -r database
    fi
    if ! [ -f MF1Dsub.sh ]; then rm MF1Dsub.sh; fi
    cd ../
done

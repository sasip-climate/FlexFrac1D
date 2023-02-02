#!/bin/bash  --
exec 2> $1/Errors_subTest.txt
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

## Set numbers of jobs to run at the same time
maxjobs=2
if ! [ "$2" ==  '' ]; then
    maxjobs=$2
fi

#Directory containing runs
workdir=$cdir/$rundir
if ! cd $workdir; then
    echo 'Invalid Run Directory:' $workdir
    exit
elif [ "$rundir" == '' ]; then
    echo 'Run Directory must be specified'
    exit
fi

#Path to the code directory
compdir='/home/auclaije/Objs1D'
if ! cd $compdir; then
    echo 'Invalid Code Directory'
    exit
fi

echo 'Submitting' $maxjobs 'jobs from the' $rundir 'directory'
echo 'using code from' $compdir

cd $workdir
logfile=$workdir/$(date +%y%m%d%H%M%S)'.log'
echo $(date) > $logfile

for dir in `ls -d ./*/`; do
    echo 'Moving in' $dir
    cd $dir
    # Check for the proper directories
    if ! [ -d Figs ]; then
        if ! [ -d /data/failles/$USER/$rundir ]; then
            mkdir /data/failles/$USER/$rundir
        fi
        if ! [ -d /data/failles/$USER/$rundir/${dir:2:-1} ]; then
	        mkdir /data/failles/$USER/$rundir/${dir:2:-1}
	        mkdir /data/failles/$USER/$rundir/${dir:2:-1}/Figs
        fi
        ln -s /data/failles/$USER/$rundir/${dir:2:-1}/Figs Figs

        echo 'Creating figures directory'
        mkdir Figs/Floes
        mkdir Figs/Summary
        mkdir Figs/Spec
    fi
    if ! [ -d database ]; then
        echo 'Creating database directory'
        mkdir database
        mkdir database/temp
    fi

    # Setup run submission
    if ! [ -f config.py ]; then	ln -s $compdir/config.py .; fi
    if ! [ -f MF1DSpecExp.py ]; then ln -s $compdir/MF1DSpecExp.py .; fi
    if ! [ -f MF1Dtest.sh ]; then
        cp $cdir/MF1Dtest.sh .
        echo "tar -chf /data/failles/$USER/data/test${dir:2:-1}.tar ./*" >> MF1Dtest.sh
    fi

    jobname=MF1DS_test_${dir:2:-1}
    # Check if run is needed
    if ! [ -f MF1D.out ]; then
        NeedRun=true
        msg="Submitting $jobname"
    elif [ $(tail MF1D.out | grep -ic 'Simulation completed') == 1 ]; then
        NeedRun=false
        msg="Skipping already completed $jobname"
    elif [ $(($(oarstat -u $USER | grep -c $jobname))) -ge 1 ]; then
        NeedRun=false
        msg="Skipping already running $jobname"
    elif [ $(grep -ic aborted MF1D.out) -ge 1 ]; then
        NeedRun=true
        msg="Resubmitting $jobname because of singular matrices"
    else
        jobsdone=`grep -ic "Time taken" MF1D.out`
        ((jobsdone+=`grep -ic "Reading" MF1D.out`))
        jobswanted=`grep Launching MF1D.out | cut -d ' ' -f 2`
        if [ "$jobswanted" == '' ]; then jobswanted=0; fi
        if [ $jobsdone -lt $jobswanted ] || [ $jobswanted -eq 0 ]; then
            NeedRun=true
            msg="Submitting $jobname. Only $jobsdone out of $jobswanted completed"
        else
            NeedRun=false
            msg="Skipping $jobname. $jobsdone out of $jobswanted completed"
        fi
    fi

    echo $(date  +'%F %T'): $msg

    # Submit run if needed
    if $NeedRun; then
        jobnum=$(($(oarstat -u $USER | grep -c MF1DS_)))
        while [ $jobnum -ge $maxjobs ]
        do
            echo -n $jobnum
	    sleep 600s
            jobnum=$(($(oarstat -u $USER | grep -c MF1DS_)))
        done

        rm OAR*std*
        oarsub -S -n $jobname --project iste-equ-failles ./MF1Dtest.sh -p "dedicated='none'" -t devel
    fi

    echo $(date  +'%F %T'): $msg >> $logfile

    cd ../
done

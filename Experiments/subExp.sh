#!/bin/bash  --
exec 2> $1/Errors_subExp.txt
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
maxjobs=5
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

#Check for another script running on the same directory
#-gt 2, 1 for script, 1 for.. grep?
subs=`ps -ef | grep $USER | grep subExp | grep -c $rundir`
if [ $subs -gt 2 ]; then
    echo "Another submission script is already running in $rundir"
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
    if ! [ -d /data/failles/$USER/data/$rundir ]; then
        mkdir /data/failles/$USER/data/$rundir
    fi

    # Setup run submission
    if ! [ -f config.py ]; then	ln -s $compdir/config.py .; fi
    if ! [ -f MF1DSpecExp.py ]; then ln -s $compdir/MF1DSpecExp.py .; fi
    if ! [ -f MF1Dsub.sh ]; then
        cp $cdir/MF1Dsub.sh .
        echo "tar -chf /data/failles/$USER/data/$rundir/${dir:2:-1}.tar ./*" >> MF1Dsub.sh
    fi

    jobname=MF1DS_${dir:2:-1}
    # Check if run is needed
    if [ $(($(oarstat -u $USER | grep -c $jobname))) -ge 1 ]; then
        NeedRun=false
        msg="Skipping already running $jobname"
    elif ! [ -f MF1D.out ]; then
	    NeedRun=true
        msg="Submitting $jobname"
    elif [ $(tail MF1D.out | grep -ic 'Simulation completed') == 1 ]; then
        NeedRun=false
        msg="Skipping already completed $jobname"
    elif [ $(grep -ic aborted MF1D.out) -ge 1 ]; then
        NeedRun=true
        msg="Resubmitting $jobname because of singular matrices"
    else 
        jobsdone=`grep -ic "Time taken" MF1D.out`
        ((jobsdone+=`grep -ic "Reading" MF1D.out`))
        jobswanted=`grep Launching MF1D.out | cut -d ' ' -f 2`
        if [ "$jobswanted" == '' ]; then jobswanted=0; fi
        if [ $jobsdone -lt $jobswanted ] || [ $jobswanted -eq 0 ] ; then         
            NeedRun=true
            if [ $(grep -ic KILLED OAR*err) -ge 1 ]; then
                t0L=`grep walltime MF1Dsub.sh`
                t0=${t0L: -8:2}
                if [ $jobsdone -lt 1 ]; then
                    progM="`tail -n 1 MF1D.out` \n"
                    prog=`grep -o '#' <<< $progM | wc -l`
                    if [ $prog -lt 1]; then
                        jdf=0.1
                    else
                        jdf=$(bc -l <<< "scale=1; $prog/10")
                    fi
                else
                    jdf=$jobsdone
                fi
                tN=$(bc -l <<< "scale=0; $t0*1.2*$jobswanted/$jdf")
                msg="Resubmitting $jobname. Only $jobsdone out of $jobswanted completed \n $progM Asking for $tN hours instead of $t0"
                sed -i -e "s/#OAR.*/#OAR -l \/cpu=1\/core=2,walltime=$tN:00:00/g" ./MF1Dsub.sh
            else
                msg="Resubmitting $jobname. Only $jobsdone out of $jobswanted completed"
            fi
        else
            NeedRun=false
            msg="Skipping $jobname. $jobsdone out of $jobswanted completed"
        fi
    fi

    echo -e $(date  +'%F %T'): $msg

    # Submit run if needed
    if $NeedRun; then
        jobnum=$(($(oarstat -u $USER | grep -c MF1DS_)))
        while [ $jobnum -ge $maxjobs ]
        do
            echo -n $jobnum
	    sleep 600s
            jobnum=$(($(oarstat -u USER | grep -c MF1DS_)))
        done

        rm OAR*std*
        oarsub -S -n $jobname --project iste-equ-failles ./MF1Dsub.sh
    fi

    echo $(date  +'%F %T'): $msg >> $logfile

    cd ../
done

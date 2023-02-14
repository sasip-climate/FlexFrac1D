#!/bin/bash  --
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
    echo 'Invalid Run Directory:' $workdir
    exit
elif [ "$rundir" == '' ]; then
    echo 'Run Directory must be specified'
    exit
fi

echo 'Checking status of the' $rundir 'directory'

exp_done=0
exp_run=0
exp_int=0
exp_left=0
exp_kill=0
exp_bug=0

cd $workdir

for dir in `ls -d ./*/`; do
    cd $dir

    jobname=MF1DS_${dir:2:-1}
    # Check if run is needed
    if ! [ -f MF1D.out ]; then
        msg="$jobname not run"
        ((exp_left+=1))
    elif [ $(tail MF1D.out | grep -ic 'Simulation completed') == 1 ]; then
        msg="$jobname completed"
        ((exp_done+=1))
    else
        jobsdone=`grep -ic "Time taken" MF1D.out`
        ((jobsdone+=`grep -ic "Reading" MF1D.out`))
        jobswanted=`grep Launching MF1D.out | cut -d ' ' -f 2`
        if [ "$jobswanted" == '' ]; then jobswanted=0; fi
        jobsmsg="$jobsdone out of $jobswanted completed"

        ((exp_int+=1))
        if [ $(grep -ic aborted MF1D.out) -ge 1 ]; then
            msg="$jobname failed because of singular matrices: $jobsmsg"
        elif [ $jobsdone -lt $jobswanted ] || [ $jobswanted -eq 0 ]; then         
            if [ $(($(oarstat -u $USER | grep -c $jobname))) -ge 1 ]; then
                msg="$jobname currently running: $jobsmsg \n `tail -n 1 MF1D.out`"
                ((exp_int-=1))
                ((exp_run+=1))
            else
                msg="$jobname incomplete: Only $jobsmsg"
                if [ $jobsdone -lt 1 ]; then
                    msg="$msg \n `tail -n 1 MF1D.out`"
                fi
                if [ -f OAR*err ]; then
                    msg="$msg \n `tail -n 1 OAR*err`"
                    if [ $(grep -ic KILLED OAR*err) -ge 1 ]; then
                        ((exp_int-=1))
                        ((exp_kill+=1))
                    elif [ $(grep -ic Error OAR*err) -ge 1 ]; then
                        ((exp_int-=1))
                        ((exp_bug+=1))
                    fi
                fi
            fi
        else
            msg="$jobname complete: $jobsmsg"
            ((exp_int-=1))
            ((exp_done+=1))
        fi
    fi

    echo -e "$msg"

    cd ../
done

echo Completed: $exp_done
echo Running: $exp_run
echo Incomplete: $exp_int
echo Bugged: $exp_bug
echo Killed: $exp_kill
echo Remaining: $exp_left
echo Total: $(($exp_done+$exp_run+$exp_int+$exp_bug+$exp_kill+$exp_left))

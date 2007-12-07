 #!/bin/bash

# redirect output to 'startup-pe' log files

SO=startup-pe.o$JOB_ID 
SE=startup-pe.e$JOB_ID

WHICH_MPIRUN="`which mpirun`"

if [ $MPIRUN_RANK = 0 ]; then 
echo ""  >>$SO 2>>$SE
echo "+ ============= begin of mpiwrap.sh" >>$SO 2>>$SE
echo "displaying info for MPIRUN_RANK=0 only" >>$SO 2>>$SE
echo "+ --------- checking environment" >>$SO 2>>$SE
cat >>$SO 2>>$SE <<EOF
HOST= $HOSTNAME
CWD= $PWD
TMPDIR= $TMPDIR
PATH= $PATH
LD_LIBRARY_PATH= $LD_LIBRARY_PATH
JOB_ID=$JOB_ID
JOB_NAME=$JOB_NAME
WHICH_MPIRUN= $MPIRUN_CMD 
MPIRUN_MPD= $MPIRUN_MPD
MPIRUN_HOST= $MPIRUN_HOST
MPIRUN_PORT= $MPIRUN_PORT
MPIRUN_PROCESSES= $MPIRUN_PROCESSES
MPIRUN_RANK= $MPIRUN_RANK
MPIRUN_NPROCS= $MPIRUN_NPROCS
MPIRUN_ID= $MPIRUN_NPROCS
EOF
echo "+ --------- executable $1 depends on dynamics libraries:"  >>$SO 2>>$SE
ldd $1 >>$SO 2>>$SE
echo "+ --------- executing: $@" >>$SO 2>>$SE
echo "delaying execution by 1 sec ..." >>$SO 2>>$SE
fi
sleep 1
echo "running (rank $MPIRUN_RANK/$MPIRUN_NPROCS)" >>$SO 2>>$SE
#numactl --cpubind 1 --membind 1 $@
numactl --membind 1 $@
ERR=$?
sleep 1

if [ $ERR != 0 ]; then
echo "(rank $MPIRUN_RANK/$MPIRUN_NPROCS) failed with error code: $ERR" >>$SO 2>>$SE
echo "+ --------- diagnostic output: $@" >>$SO 2>>$SE
vstat -v >>$SO 2>>$SE
fi

if [ $MPIRUN_RANK = 0 ]; then 
echo "+ ============= end of mpiwrap.sh" >>$SO 2>>$SE
fi

exit $ERR


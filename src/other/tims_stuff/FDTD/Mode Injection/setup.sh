#!/bin/sh

# Set HOME location.
HOMEDIR=/data/users/twp11/ITO_MIM/Passive

echo 'Home directory: ' $HOMEDIR
# Create Queued directory
WORKDIR=/data/users/twp11/Queued
mkdir $WORKDIR

# Create WORK directory.
DATE=$(date '+%H:%M:%S')
WORKDIR=/data/users/twp11/Queued/Queued_$DATE
mkdir $WORKDIR
echo 'Work directory: ' $WORKDIR

# Copy files from Input to WORK.
cp *.lua *.qscript *.sh $WORKDIR

# Change to Tools.
TDIR=$HOMEDIR/Input/Tools
cd $TDIR

# Copy files from Tools to WORK.
cp meta luacfg *.m *.txt *.gpl *.in $WORKDIR


# Change to WORK.
cd $WORKDIR
# Run Script - passing through HOME and WORK directories.
echo -n 'ID Directory: '
qsub runcode.qscript -v HOMEDIR=$HOMEDIR,WORKDIR=$WORKDIR



echo 'Script Submitted.'

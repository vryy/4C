#!/usr/local/bin/bash
#
# File: mpi_hello_job.sh
#
# This is a demo submit script for MPI jobs in the
# Computer Center of the University of Ulm
# It expects GridEngine as the underlaying scheduler.
#
# Start of this script after compiling hello.c or hello.f with:
# 
#         qsub mpi_hello_job.sh 
# 
# Check the manpage of qsub(1)
# ----------------------------------------------------------------------
#
# Preparation, done e.g. interactively:
# 
# - Compilation: 
#
#     mpcc  -o hello hello.c         -lmpi
#     mpCC  -o hello hello.cc  -mt   -lmpi   
#     mpf77 -o hello -dalign hello.f -lmpi
#     mpf90 -o hello -dalign hello.f -lmpi
#
#     Note   For the Fortran interface, the  dalign option is necessary to avoid the 
#            possibility of bus errors. (The underlying C or C++ routines in Sun MPI 
#            internals assume that parameters and buffer types passed as REALs are 
#            double-aligned.)
#     
# - Parallel interactive execution (here on 2 processors):
#     
#     mprun -np 2 ./hello
# 
# ----------------------------------------------------------------------
# 
# Statusinfo, after the batchjob is queued:
# 
#     qstat -u "username"
#     qstat -j "number of job as returned by qsub" 
# 
# ----------------------------------------------------------------------
# add resource requests; '#$' is interpreted by qsub(1) so it's NOT a comment

# give the job a name

#$ -N cca

# setup resource requests like number of CPUs, memory, ...
# eg 512M of virtual memory and 1 minute of CPU time, and 2-8 parallel CPUs 

#$ -l h_vmem=512M
#$ -l h_cpu=00:05:00
#$ -pe mpi 1-2

# the following is necessary for parallel jobs 
#$ -l cre=true

# we assigned special resource like availability of software and the like to some nodes
#	request a queue with gaussian support:		gaussian=true
#	request a queue CADENCE... support:		etechnik=true
#	request a queue with prefetching enabled:	prefetch=true
#	request a queue not running beta OS:		nobeta=true
#
# the defaults are
#	no prefetching
#	ok to submit to Solaris 9 box
#	no special software requirements
#
# the software dependencies just leave more flexibility to the host
# admins but the software is very likely to be available at every node

#$ -l prefetch=true
# --- uncommented since not neeeded here --- #$ -l gaussian=true


# to change the default location and name of
#	stderr:<jobname>.e<job_id>
#	stdout:<jobname>.o<job_id>
# eg discard stdout

# --- uncommented since output is wanted --- #$ -o /dev/null


# by default we get mail when a job is started/finished/aborted
# but if you don't want it

#$ -m n 


# ----------------------------------------------------------------------
# the real work starts here

echo '...the real work starts here' `/users/stgt/<name>/ccarat/cca.exe`

# make special 'computer center extensions' available and setup environment

. /etc/bash_profile
echo '...after ". /etc/bash_profile"' 


# make the script abort on any error
# set -e


# if you need special packages use the same 'option' commands as
# interactively; you can get a list by 'options'
# eg for Gaussian 98
# --- uncommented since not neeeded here --- option g98


# you may want to create a working diretory for this job
# in the local scratch area
#
WORKBASE=/work/$USER
WORKDIR=$JOB_NAME.$JOB_ID
cd $WORKBASE
mkdir $WORKDIR
cd $WORKDIR


# if you need files from the submit host to be copied to the
# local scratch area on the execution host add something like
#
# rcp$COD_O_HOST:<original_path>/<original_name> $WORKBASE/$WORKDIR
rcp $COD_O_HOST:/users/stgt/<name>/ccarat/dyn_shell_lind.dat $WORKBASE/$WORKDIR
rcp $COD_O_HOST:/users/stgt/<name>/ccarat/cca.exe $WORKBASE/$WORKDIR
#
# be aware that by default the batch system changes the current
# directory to the one from which the job was submitted
#
# also check rsync(1) if it's more than a single file


# ok, start the real work (my app here is just 'date')
# (output gets discarded because we redirected stdout to /dev/null)

echo '...preparing the parallel application start' `/users/stgt/<name>/ccarat/cca.exe`
echo '...this job is running with' NSLOTS=$NSLOTS , $TMPDIR/machines:
cat $TMPDIR/machines
echo 'EOF' 
 
echo '...starting the application' `/users/stgt/<name>/ccarat/cca.exe`
 
# the following command starts your MPI-executable in parallel 
#                       --------------------------------------       
 
$CODINE_ROOT/mpi/MPRUN -np $NSLOTS -Mf $TMPDIR/machines ./cca.exe dyn_shell_lind.dat testout
#                                                       ------- 
 
# you may substitute "hello" by your executable and you may add additional
# arguments, but you must not change the other parts of this command. 
 
echo '...end of the application' `/users/stgt/<name>/ccarat/cca.exe`


# do cleanup if required
#
# cd $WORKBASE/$WORKDIR
# ### add code to save important files ###
# cd $WORKBASE
# rm -rf $WORKBASE/$WORKDIR

echo '...end of clean up' `/users/stgt/gee/ccarat/cca.exe`

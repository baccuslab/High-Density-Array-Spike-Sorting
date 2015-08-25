#!/bin/bash

while true 
do 
    echo 'Waiting for the job.sh file: ' $1/$2/$2_job.sh '...'
    [ -f  $1/$2/$2_job.sh ] && break
    sleep 10
done

# Move old job log files to a directory
# cd log
# mkdir $2
cd ~/


# lastdir=$(ls -d [0-9][0-9][0-9][0-9] | tail -1)
# newdir=$((++lastdir))
# dir=$(printf "%04u" $newdir)
# #mkdir $(printf "%04u" $newdir) 
# #mv job* $(printf "%04u" $newdir)
# mkdir $dir
# mv job* $dir

ls -l $1/$2/$2_job.sh
qsub $1/$2/$2_job.sh
echo 'Job  submitted!'

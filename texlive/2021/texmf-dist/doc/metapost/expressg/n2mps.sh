#! /bin/sh

#####################################################################
# Shell script n2mps.sh
# Call as: n2mps.sh basename
# List each file basename.* in the directory and run the Perl script
# to copy each basename.N to basenameN.mps, where N is an integer
#
# Copyright 2000, Mauro S. Costa and Peter R. Wilson
#####################################################################

basename=${1:?"A file basename is required."}
extname=N.mps
echo Files $basename.N, where N is a number, will be copied to $basename$extname
for file in `ls $basename.*`
do
 n2mpsprl.prl $file
done

####################### end shell script ###########################


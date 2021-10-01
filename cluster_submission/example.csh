#!/bin/csh
#
# Script:  R_loops.csh
# Usage: For submitting multiple batch R jobs to LSF,
#     an example of looping over years and models
# Author:  NCSU HPC
#
# 
## To run, type:
#     ./R_loops.csh [# years] [# models]
#  Script must have execute permissions, i.e.,
#     chmod u+x R_loops.csh

module load R

if ($#argv != 2) then
	echo "Usage: You need to feed two arguments to this program which is"
 	echo "the number of years and the number of models. For example,"
    	echo "./R_loops.csh 2 2"
    	exit 1
endif 
 
# Specify number of jobs to submit
set numYears = $1

set numModels = $2

# Initialize year loop counter 
set year = 1

while ($year <= $numYears)


    # Initialize model loop counter
    set model = 1
    while ($model <= $numModels)

       echo "Submit job year = $year and model = $model"

       bsub -n 1 -W 30 -oo out.year=$year.model=$model -eo err.year=$year.model=$model "Rscript codehpc.R $year $model"

       @ model++

    end 
 
    @ year++
    
end

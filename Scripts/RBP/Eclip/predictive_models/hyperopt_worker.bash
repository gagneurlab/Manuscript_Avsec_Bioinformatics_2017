#!/bin/bash
#
# Run the hyperopt mongo worker
#
# author: avsec
############################################

#SBATCH --mem=16G
#SBATCH -c 4

# TODO - edit
MONGO_IP=localhost
THIS_DIR=$PWD
#

DB=RBP__Eclip
export PYTHONPATH=${THIS_DIR}


# run the worker
hyperopt-mongo-worker \
    --mongo=${MONGO_IP}:27017/$DB \
    --poll-interval=1 \
    --reserve-timeout=3600

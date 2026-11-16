#!/bin/bash

# check the genome installation and prepare it if necessary
runWorkflowStep 1 check_genome $MODULES_DIR/genome/check_genome.sh

# set derivative environment variables and file paths
export PIPELINE_SHARED_DIR=${PIPELINE_DIR}/shared
source ${PIPELINE_SHARED_DIR}/workflow.sh

# run sniffles
runWorkflowStep 2 sniffles sniffles.sh

# run sniffles
runWorkflowStep 3 cutesv cutesv.sh

# clean up
rm -fr $TMP_DIR_WRK_SMALL
rm -fr $TMP_DIR_WRK_SHM

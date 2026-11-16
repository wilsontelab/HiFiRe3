#!/bin/bash

# check the genome installation and prepare it if necessary
runWorkflowStep 1 check_genome $MODULES_DIR/genome/check_genome.sh

# set derivative environment variables and file paths
export PIPELINE_SHARED_DIR=${PIPELINE_DIR}/shared
source ${PIPELINE_SHARED_DIR}/workflow.sh
source $MODULES_DIR/align/set_read_file_vars.sh
mkdir -p $PLOTS_DIR
export ALIGN_DIR=${ACTION_DIR}/align

# set working directory to READ_FILE_DIR to avoid too-long argument list with multiple read files
cd ${READ_FILE_DIR}

# align read sequences to reference genome
runWorkflowStep 2 align align/align.sh

# align read sequences to reference genome
runWorkflowStep 3 align align/index.sh

# clean up
rm -fr $TMP_DIR_WRK_SMALL
rm -fr $TMP_DIR_WRK_SHM

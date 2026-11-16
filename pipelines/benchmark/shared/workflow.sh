# initialize actions common to benchmark actions

# set derivative environment variables and file paths
source $MODULES_DIR/rust/set_rust_vars.sh
source $MODULES_DIR/genome/set_genome_vars.sh
source $MODULES_DIR/align/set_alignment_vars.sh
source $MODULES_DIR/REs/set_re_vars.sh
source $MODULES_DIR/library/set_library_vars.sh
source $MODULES_DIR/analyze/set_analysis_vars.sh
source $MODULES_DIR/benchmark/set_benchmark_vars.sh

# set temporary directories
source $MODULES_DIR/utilities/shell/create_temp_dir_small.sh
source $MODULES_DIR/utilities/shell/create_temp_dir_shm.sh

# set other working parameters
export SORT_RAM_PER_CPU_INT=$(((TOTAL_RAM_INT - 50000000000) / N_CPU))
export SAMTOOLS_CPU=$((N_CPU - 1))
export ALN_CPU=$(( N_CPU - 1 )) # minimap2 uses 1 core for IO

#!/bin/bash

# ===================並列化の設定===========================                                         
#SBATCH -p F1cpu                                                                                     
#SBATCH -N 1                                                                                         
#SBATCH -n 1                                                                                         
#SBATCH -c 128                                                                                       
#SBATCH --job-name="MHdata2(distorted)"                                                              
#SBATCH --mail-type=BEGIN                                                                            
#SBATCH --mail-type=END                                                                              
#SBATCH --mail-user=6222530@ed.tus.ac.jp 

SYS_NUM="$1"
SYS_SITE_A="$2"
SYS_SITE_B="$3"
MIN_UPS_SPIN="$4"
MAX_UP_SPIN="$5"
J_RED="$6"
J_GREEN="$7"
J_BLUE="$8"
DIR_JSET0="$9"
DIR_JSET1="$10"
DIR_JSET2="$11"
DIR_OUTPUT_EVAL="$12"
DIR_OUTPUT_TIME="$13"
DIR_OUTPUT_SXX_REL="$14"
DIR_OUTPUT_SZZ_REL="$15"
START_UP_SPIN="$16"
END_UP_SPIN="$17"
LANCZOS_TYPE="$18"

srun $PLOT_MH_EXE_FILE "$SYS_NUM" "$SYS_SITE_A" "$SYS_SITE_B" "$MIN_UP_SPIN" "$MAX_UP_SPIN" "$J_RED" "$J_GREEN" "$J_BLUE" "$DIR_JSET0" "$DIR_JSET1" "$DIR_JSET2" "$DIR_OUTPUT_EVAL" "$DIR_OUTPUT_TIME" "$DIR_OUTPUT_SXX_REL" "$DIR_OUTPUT_SZZ_REL" "$START_UP_SPIN" "$END_UP_SPIN" "$LANCZOS_TYPE"

#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# 1. Check inputs
if [ -z "$1" ]; then
    echo "Error: Please specify the mode."
    echo "Usage:   ./run.sh [MODE] [OPTIONAL_RESTART_DIR]"
    echo "Example Fresh:   ./run.sh I"
    echo "Example Restart: ./run.sh I results/ModeI_12-05-2024_10-00-00"
    exit 1
fi

MODE=$1
RESTART_DIR=$2

# CONFIGURATION
SOLVER_NAME="pffm_solver"
SOURCE_FILE="apps/main.cpp"

# 2. Automated Build
echo "---------------------------------------"
echo "Checking Build System..."

mkdir -p "$SCRIPT_DIR/build"
cd "$SCRIPT_DIR/build"

if [ ! -f "Makefile" ]; then
    echo "Makefile not found. Configuring with CMake..."
    cmake ..
fi

echo "Compiling..."
make

# Verify executable exists
if [ ! -f "$SOLVER_NAME" ]; then
    echo "Error: Compilation failed. $SOLVER_NAME not found."
    exit 1
fi

# 3. Handle Directory Logic
cd "$SCRIPT_DIR" # Go back to project root

if [ -z "$RESTART_DIR" ]; then
    # --- CASE A: FRESH RUN ---
    TIMESTAMP=$(date "+%d-%m-%Y_%H-%M-%S")
    OUTPUT_DIR="results/mode-${MODE}_${TIMESTAMP}"
    
    echo "---------------------------------------"
    echo "Initialising FRESH run for mode-$MODE"
    echo "Output directory: $OUTPUT_DIR"
    echo "---------------------------------------"
    
    mkdir -p "$OUTPUT_DIR"
else
    # --- CASE B: RESTART ---
    OUTPUT_DIR="$RESTART_DIR"
    
    if [ ! -d "$OUTPUT_DIR" ]; then
        echo "Error: Restart directory '$OUTPUT_DIR' does not exist!"
        exit 1
    fi
    
    echo "---------------------------------------"
    echo "Resuming run in EXISTING directory"
    echo "Target Directory: $OUTPUT_DIR"
    echo "---------------------------------------"
fi

# 4. Snapshot Parameters
echo "Snapshotting parameters from $SOURCE_FILE..."
cp "$SCRIPT_DIR/$SOURCE_FILE" "$OUTPUT_DIR/source_backup.cpp"

# 5. Prepare Execution
cd "$OUTPUT_DIR"
FULL_PATH=$(pwd)

if [ -z "$RESTART_DIR" ]; then
    # FRESH: Pass Mode + Output Directory
    EXEC_ARGS="$MODE $FULL_PATH"
else
    # RESTART: Find last increment
    RESTART_INC=$(ls inc_${MODE}=*.csv 2>/dev/null | sed "s/inc_${MODE}=\(.*\)\.csv/\1/" | sort -n | tail -1)
    if [ -z "$RESTART_INC" ]; then 
        echo "Warning: No output files found. Defaulting restart to 0."
        RESTART_INC=0
    else
        echo ">> Auto-detected restart increment: $RESTART_INC"
    fi

    # RESTART: Pass Mode + Directory + Increment
    EXEC_ARGS="$MODE $FULL_PATH $RESTART_INC"
fi

echo "Running $SOLVER_NAME in mode-$MODE..."

# Thread controls (OpenBLAS/Eigen optimisation)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# 6. Run
"$SCRIPT_DIR/build/$SOLVER_NAME" $EXEC_ARGS

echo "---------------------------------------"
echo "Done."
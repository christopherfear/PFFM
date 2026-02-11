#!/bin/bash

# Get the absolute path of the directory containing this script
# (This ensures the script works no matter where you call it from)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# 1. Check inputs
if [ -z "$1" ]; then
    echo "Error: Please specify the mode."
    echo "Usage:   ./run.sh [MODE] [OPTIONAL_RESTART_DIR]"
    echo "Example Fresh:   ./run.sh I"
    echo "Example Restart: ./run.sh I results/modeI_12-05-2024_10-00-00"
    exit 1
fi

MODE=$1
RESTART_DIR=$2

# 2. Set Compiler Variables
if [ "$MODE" == "I" ]; then
    SOLVER_NAME="solver_modeI"
    SOURCE_SNAPSHOT="apps/main_modeI.cpp"
elif [ "$MODE" == "II" ]; then
    SOLVER_NAME="solver_modeII"
    SOURCE_SNAPSHOT="apps/main_modeI.cpp"
else
    echo "Error: Invalid mode '$MODE'. Please use 'I' or 'II'."
    exit 1
fi

# 3. Automated Build (Using Absolute Paths)
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

# 4. Handle Directory Logic
cd "$SCRIPT_DIR" # Go back to project root

if [ -z "$RESTART_DIR" ]; then
    # --- CASE A: FRESH RUN (No directory provided) ---
    TIMESTAMP=$(date "+%d-%m-%Y_%H-%M-%S")
    OUTPUT_DIR="results/Mode${MODE}_${TIMESTAMP}"
    
    echo "---------------------------------------"
    echo "Initializing FRESH run for Mode $MODE"
    echo "Output Directory: $OUTPUT_DIR"
    echo "---------------------------------------"
    
    mkdir -p "$OUTPUT_DIR"
else
    # --- CASE B: RESTART (Directory provided) ---
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

# 5. Snapshot & Run
echo "Snapshotting parameters from $SOURCE_SNAPSHOT..."
cp "$SCRIPT_DIR/$SOURCE_SNAPSHOT" "$OUTPUT_DIR/source_backup.cpp"

cd "$OUTPUT_DIR"

# --- CRITICAL FIX: Get the Absolute Path ---
# Since we just cd'd into the folder, $(pwd) gives us the full, real path.
FULL_PATH=$(pwd)

if [ -z "$RESTART_DIR" ]; then
    # --- FRESH RUN ---
    # Pass the full path so C++ writes to the correct place
    EXEC_ARGS="$FULL_PATH"
    
else
    # --- RESTART RUN ---
    RESTART_INC=$(ls inc=*.csv 2>/dev/null | sed 's/inc=\(.*\)\.csv/\1/' | sort -n | tail -1)
    
    if [ -z "$RESTART_INC" ]; then 
        echo "Warning: No output files found. Defaulting restart to 0."
        RESTART_INC=0
    else
        echo ">> Auto-detected restart increment: $RESTART_INC"
    fi

    # Pass the full path + increment
    EXEC_ARGS="$FULL_PATH $RESTART_INC"
fi

echo "Running $SOLVER_NAME..."

# Thread controls
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Run the solver
"$SCRIPT_DIR/build/$SOLVER_NAME" $EXEC_ARGS

echo "---------------------------------------"
echo "Done."
#!/bin/bash

# 1. Check if the user provided an input (I or II)
if [ -z "$1" ]; then
    echo "Error: Please specify the mode."
    echo "Usage: ./run.sh I   OR   ./run.sh II"
    exit 1
fi

# 2. Set variables based on the input
MODE=$1

if [ "$MODE" == "I" ]; then
    SOLVER_NAME="solver_modeI"
    SOURCE_FILE="apps/main_modeI.cpp"
elif [ "$MODE" == "II" ]; then
    SOLVER_NAME="solver_modeII"
    SOURCE_FILE="apps/main_modeII.cpp"
else
    echo "Error: Invalid mode '$MODE'. Please use 'I' or 'II'."
    exit 1
fi

# 3. Automated build
echo "---------------------------------------"
echo "Checking Build System..."

mkdir -p build
cd build

# Check if Makefile exists. If NOT, we must run CMake.
if [ ! -f "Makefile" ]; then
    echo "Makefile not found. Configuring with CMake..."
    cmake ..
fi

# Run Make. 
# NOTE: If you changed CMakeLists.txt, Make will detect it 
# and re-run CMake automatically here.
echo "Compiling..."
make -j4

cd ..
echo "Build Complete."

# 4. Setup Timestamp and Directories
TIMESTAMP=$(date "+%d-%m-%Y_%H-%M-%S")
OUTPUT_DIR="results/Mode${MODE}_${TIMESTAMP}"

echo "---------------------------------------"
echo "Initializing run for Mode $MODE"
echo "Output Directory: $OUTPUT_DIR"
echo "---------------------------------------"

mkdir -p "$OUTPUT_DIR"

# 5. Snapshot the specific source code for this mode
echo "Snapshotting parameters from $SOURCE_FILE..."
cp "$SOURCE_FILE" "$OUTPUT_DIR/parameters_snapshot.cpp"

# 6. Move into dir and run
cd "$OUTPUT_DIR"

# Check if the executable exists first
if [ ! -f "../../build/$SOLVER_NAME" ]; then
    echo "Error: Executable ../../build/$SOLVER_NAME not found!"
    echo "Did you run 'cmake ..' and 'make' in the build folder?"
    exit 1
fi

echo "Running $SOLVER_NAME..."
../../build/$SOLVER_NAME

echo "---------------------------------------"
echo "Done."

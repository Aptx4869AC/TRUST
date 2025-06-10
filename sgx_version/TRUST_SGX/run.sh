#!/bin/bash
# Check if ./app exists and then clean up
if [ -f "./app" ]; then
    make clean
fi

# Build the project
make

# Run the application
./app

# Check if ./app exists and then clean up
if [ -f "./app" ]; then
    make clean
fi

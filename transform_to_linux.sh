#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

# Assign the arguments to variables
filename=$1 # File with bacterial species
# Check and convert Windows line endings (CRLF) to Linux line endings (LF)
if file "$filename" | grep -q "CRLF"; then
    echo "Converting Windows line endings to Linux line endings..."
    sed -i 's/\r$//' "$filename"
fi
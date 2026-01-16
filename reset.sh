#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: $0 <test_name>"
    echo "Example: $0 test1"
    exit 1
fi

TEST=$1

python transplant.py ${TEST}/base ${TEST}

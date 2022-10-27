#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PROJECT_DIR=$(dirname $SCRIPT_DIR)

for i in $(find $PROJECT_DIR -name "*.nim");
do
    if [[ -e "${i%.nim}" ]]; then 
        echo "removing ${i%.nim}"
        rm "${i%.nim}"
    fi
done
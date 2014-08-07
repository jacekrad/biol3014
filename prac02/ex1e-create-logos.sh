#!/bin/bash

python exercise1e.py 1>/dev/null 2>alignment_filenames.txt

for alignment_filename in `cat alignment_filenames.txt`; do
    export logo_filename=`echo ${alignment_filename} | cut -f 1 -d "."`.svg
    echo "generating ${logo_filename}"
    weblogo --format=svg < ${alignment_filename} > ${logo_filename}
done


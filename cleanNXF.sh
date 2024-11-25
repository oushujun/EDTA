#!/usr/bin/env bash

rm -rf .nextflow*
echo "Cleaned .nextflow..."
rm -rf .nextflow.pid
echo "Cleaned .nextflow.pid..."
for i in $(ls work | grep -v "conda");
do
    rm -rf "work/$i"
done
echo "Cleaned work..."

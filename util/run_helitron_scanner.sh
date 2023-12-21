#!/bin/bash -login
### This script was modified from https://github.com/mcstitzer/maize_v4_TE_annotation/blob/master/helitron/run_helitron_scanner.sh
### Original author: Michelle Stitzer, Apr 11, 2018
### Modifier: Shujun Ou (shujun.ou.1@gmail.com), May 1, 2019

### specify the genome file
GENOME=$1

### the base path of this script
path=$(dirname "$0")

## where to find HelitronScanner.jar
HSDIR=$path/../bin/HelitronScanner

### preset CPU and max memory
CPU=4
MEMGB=150 #Gb

### allow user to specify CPU number to run HelitronScanner
if [ ! -z "$2" ];
	then CPU=$2
fi

###########################
##   DIRECT ORIENTATION  ##
###########################

##find helitron heads
### will load each chromosome into memory, without splitting into 1Mb batches (-buffer_size option ==0) 
java -Xmx${MEMGB}g -jar ${HSDIR}/HelitronScanner.jar scanHead -lcv_filepath ${HSDIR}/TrainingSet/head.lcvs -g $GENOME -buffer_size 0 -output ${GENOME}.HelitronScanner.head

## helitron tails
java -Xmx${MEMGB}g -jar ${HSDIR}/HelitronScanner.jar scanTail -lcv_filepath ${HSDIR}/TrainingSet/tail.lcvs -g $GENOME -buffer_size 0 -output ${GENOME}.HelitronScanner.tail

## pair the ends to generate possible helitrons
java -Xmx${MEMGB}g -jar ${HSDIR}/HelitronScanner.jar pairends -head_score ${GENOME}.HelitronScanner.head -tail_score ${GENOME}.HelitronScanner.tail -output ${GENOME}.HelitronScanner.pairends

## draw the helitrons into fastas
java -Xmx${MEMGB}g -jar ${HSDIR}/HelitronScanner.jar draw -pscore ${GENOME}.HelitronScanner.pairends -g $GENOME -output ${GENOME}.HelitronScanner.draw -pure_helitron
 
############################
##    REVERSE COMPLEMENT  ##
############################

##find helitron heads
### will load each chromosome into memory, without splitting into 1Mb batches (-buffer_size option ==0) 
java -Xmx${MEMGB}g -jar ${HSDIR}/HelitronScanner.jar scanHead -lcv_filepath ${HSDIR}/TrainingSet/head.lcvs -g $GENOME -buffer_size 0 --rc -output ${GENOME}.HelitronScanner.rc.head

## helitron tails
java -Xmx${MEMGB}g -jar ${HSDIR}/HelitronScanner.jar scanTail -lcv_filepath ${HSDIR}/TrainingSet/tail.lcvs -g $GENOME -buffer_size 0 --rc -output ${GENOME}.HelitronScanner.rc.tail

## pair the ends to generate possible helitrons
java -Xmx${MEMGB}g -jar ${HSDIR}/HelitronScanner.jar pairends -head_score ${GENOME}.HelitronScanner.rc.head -tail_score ${GENOME}.HelitronScanner.rc.tail --rc -output ${GENOME}.HelitronScanner.rc.pairends

## draw the helitrons
java -Xmx${MEMGB}g -jar ${HSDIR}/HelitronScanner.jar draw -pscore ${GENOME}.HelitronScanner.rc.pairends -g $GENOME -output ${GENOME}.HelitronScanner.draw.rc -pure_helitron
 

#########################
##   tab format output ##
######################### 

### will read in both $GENOME.HelitronScanner.draw.hel.fa $GENOME.HelitronScanner.draw.rc.hel.fa and filter out candidates based on prediction scores (min = 12) and target site (AT or TT).
perl $path/format_helitronscanner_out.pl -genome $GENOME -sitefilter 1 -minscore 12 -keepshorter 1 -extout 0




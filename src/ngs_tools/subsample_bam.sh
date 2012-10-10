#!/bin/bash

bamfile=$1
subs=$2

samtools view $bamfile | ./subsample $subs

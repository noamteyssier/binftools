#!/usr/bin/env bash

fn=$1

fgrep 'GENOME_TERRITORY' -A 1 $fn | fgrep 'GENOME' -v | sed "s/$/	$fn/" 

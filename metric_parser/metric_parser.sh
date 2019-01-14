#!/usr/bin/env bash

# input first arg or stdin
[ $# -ge 1 -a -f "$1" ] && input="$1" || input="-"

fgrep 'GENOME_TERRITORY' -A 1 | fgrep 'GENOME' -v 

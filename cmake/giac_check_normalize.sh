#!/bin/sh
# Normalize giac check output for stable diff comparison.
# Usage: giac_check_normalize.sh <mode> <infile> <outfile>
# Modes: geo, cas_floats

set -e

mode=$1
infile=$2
outfile=$3

if [ -z "$mode" ] || [ -z "$infile" ] || [ -z "$outfile" ]; then
  echo "Usage: $0 <geo|cas_floats> <infile> <outfile>" >&2
  exit 1
fi

case "$mode" in
  geo)
    perl -pe '
      if (/pnt\(pnt/) {
        s/\[536870[0-9]+\]/[ID]/g;
        s/,536870[0-9]+,/,ID,/g;
        s/,[0-9]{1,10}\]/,L]/g;
        s/,[0-9]{1,10},/,L,/g;
      }
    ' "$infile" > "$outfile"
    ;;
  cas_floats)
    perl -pe '
      s/(?<![\w.])(-?\d+\.\d+)(?![\w.])/
        my $n = $1;
        sprintf("%.10g", $n + 0.0)
      /ge
    ' "$infile" > "$outfile"
    ;;
  *)
    echo "Unknown mode: $mode" >&2
    exit 1
    ;;
esac

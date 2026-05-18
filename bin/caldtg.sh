#!/bin/bash

# 1. Print help message if arguments are missing
if [ $# -lt 2 ]; then
    echo "Usage: $(basename "$0") DATE_STRING [+|-]TT"
    echo "DATE_STRING format: YYMMDDHH or YYYYMMDDHH"
    exit 0
fi

dtg=$1
offset=$2

if [[ ${#dtg} -ne 8 ]] && [[ ${#dtg} -ne 10 ]]; then
  echo "DATE_STRING length incorrect; use 8-digit or 10-digit format"
  exit 1
fi

#For YYMMDDHH input, the YY was converted to YYYY with epoch setting:
#YY between 69~99: 19YY
#YY between 00~68: 20YY
if [ ${#dtg} -eq 8 ]; then
  short=1
  yy_short=${dtg:0:2}
  remain=${dtg:2:6}
  if [ $yy_short -ge 70 ]; then
    dtg="19$yy_short$remain"
  else
    dtg="20$yy_short$remain"
  fi
fi

# Extract parts for formatting (YYYY-MM-DD HH:00)
yy=${dtg:0:4}
mm=${dtg:4:2}
dd=${dtg:6:2}
hh=${dtg:8:2}

# Parse the offset and decide on the syntax
# Strip the leading '+' if it exists, leave the '-' for math
clean_offset=${offset#+}

# Detect if it's a subtraction to use "ago" for maximum compatibility
if [[ "$clean_offset" == -* ]]; then
    # Remove the minus sign and use 'ago'
    abs_offset=${clean_offset#-}
    calc_string="$yy-$mm-$dd $hh:00 $abs_offset hours ago"
else
    # Use standard addition
    calc_string="$yy-$mm-$dd $hh:00 $clean_offset hours"
fi

output=`date -d "$calc_string" +"%Y%m%d%H" 2> /dev/null`

# Check if the date command actually succeeded
if [ $? -ne 0 ]; then
    echo "Error: Invalid date or offset provided."
    exit 1
fi

# Execute the calculation
# Use %y for 2-digit year output to match original script's YYMMDDHH
if [[ "$short" -eq 1 ]]; then
  echo ${output:2:8}
else
  echo $output
fi

#!/bin/bash

# 1. Print help message if arguments are missing
if [ $# -lt 2 ]; then
    echo "Usage: $(basename "$0") yymmddhh [+|-]TT"
    exit 0
fi

dtg=$1
offset=$2

# Extract parts for formatting (YY-MM-DD HH:00)
# We assume 20xx for the century
yy=${dtg:0:2}
mm=${dtg:2:2}
dd=${dtg:4:2}
hh=${dtg:6:2}

# 2. Parse the offset and decide on the syntax
# Strip the leading '+' if it exists, leave the '-' for math
clean_offset=${offset#+}

# 3. Detect if it's a subtraction to use "ago" for maximum compatibility
if [[ "$clean_offset" == -* ]]; then
    # Remove the minus sign and use 'ago'
    abs_offset=${clean_offset#-}
    calc_string="20$yy-$mm-$dd $hh:00 $abs_offset hours ago"
else
    # Use standard addition
    calc_string="20$yy-$mm-$dd $hh:00 $clean_offset hours"
fi

# Execute the calculation
# Use %y for 2-digit year output to match original script's YYMMDDHH
date -d "$calc_string" +"%y%m%d%H" 2>/dev/null

# Check if the date command actually succeeded
if [ $? -ne 0 ]; then
    echo "Error: Invalid date or offset provided."
    exit 1
fi

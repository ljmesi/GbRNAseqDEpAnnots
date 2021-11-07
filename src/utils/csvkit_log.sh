#!/bin/bash
# https://jtrocinski.com/posts/Bash-Error_handling_in_scripts.html
set -euo pipefail

echo -e "csvstat $1 $2:\n" &&
csvstat "$1" "$2" &&
echo "-----------------------------------------" &&
echo -e "csvlook $1 $2:\n" &&
csvlook "$1" "$2" | sed -n 1,10p
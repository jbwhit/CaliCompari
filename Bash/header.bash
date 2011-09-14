#!/bin/bash
# pre-parses header files so later python code can read it in

# Directories to be set
HEAD_DIR="Headers/"

for file in $(ls 2004*)
do
  echo $file
  TEMP="$(echo "$file" | awk -F . '{print $2}')"
  sed 's/=//g; s/[ \t]*$//; /^$/d' $file | awk -F \' '{print $1, $2}' | awk '{if (NF > 1) print $1, $2}' > $TEMP.header

done

echo "Done"

# sed 's/=//g; s/[ \t]*$//; /^$/d' 2004.0065 | awk -F \' '{print $1, $2}' | awk '{if (NF > 1) print $1, $2}' > tryout

# Parsing practice
# sed 's/=//g' 2004.0065 > try
# sed 's/[ \t]*$//' try > try2
# sed '/^$/d' try2 > try3
# awk -F \' '{print $1, $2}' try3 > try4
# awk '{if (NF > 1) print $1, $2}' try4 > tryout
# 


# read the necessary defines from the input file and 
# write them to $definefile
awk '/necessary defines:/ { for (i=4; i<=NF; i++) print $i }' "$inputfile" > "$definefile"
awk '/#define / {print $3 "=" $4}' "$inputfile" | sed -e 's/(//g' -e 's/)//g' >> "$definefile"

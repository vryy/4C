
# read the necessary defines from the input file and 
# write them to $definefile
awk '/necessary defines:/ { for (i=4; i<=NF; i++) print $i }' "$inputfile" > "$definefile"
awk '/#define / {print $3 "=" $4}' "$inputfile" | sed -e 's/(//g' -e 's/)//g' >> "$definefile"


# get those that are defined in the PROBLEM SIZE section
MAXNOD=`awk '
  /BEGIN/ { look=0 } 
  /^--*--/ { look=0 }
  /^--*--PROBLEM SIZE/ { look=1 }
  /^MAXNOD[ \t]*[0-9]*/ { if (look==1) { print $2 } }
' "$inputfile"`

if [ "x$MAXNOD" = x ]; then
  echo $0: Warning: MAXNOD undefined
else
  if grep MAXNOD "$definefile" 2>&1 > /dev/null ; then
    sed -e "s/MAXNOD=[0-9]*/MAXNOD=$MAXNOD/" "$definefile" > "$definefile".new
    mv "$definefile".new "$definefile"
  else
    echo MAXNOD=$MAXNOD >> "$definefile"
  fi
fi


MAXELE=`awk '
  /BEGIN/ { look=0 } 
  /^--*--/ { look=0 }
  /^--*--PROBLEM SIZE/ { look=1 }
  /^MAXELE[ \t]*[0-9]*/ { if (look==1) { print $2 } }
' "$inputfile"`

if [ "x$MAXELE" = x ]; then
  echo $0: Warning: MAXELE undefined
else
  if grep MAXELE "$definefile" 2>&1 > /dev/null ; then
    sed -e "s/MAXELE=[0-9]*/MAXELE=$MAXELE/" "$definefile" > "$definefile".new
    mv "$definefile".new "$definefile"
  else
    echo MAXELE=$MAXELE >> "$definefile"
  fi
fi


MAXDOFPERNODE=`awk '
  /BEGIN/ { look=0 } 
  /^--*--/ { look=0 }
  /^--*--PROBLEM SIZE/ { look=1 }
  /^MAXDOFPERNODE[ \t]*[0-9]*/ { if (look==1) { print $2 } }
' "$inputfile"`

if [ "x$MAXDOFPERNODE" = x ]; then
  echo $0: Warning: MAXDOFPERNODE undefined
else
  if grep MAXDOFPERNODE "$definefile" 2>&1 > /dev/null ; then
    sed -e "s/MAXDOFPERNODE=[0-9]*/MAXDOFPERNODE=$MAXDOFPERNODE/" "$definefile" > "$definefile".new
    mv "$definefile".new "$definefile"
  else
    echo MAXDOFPERNODE=$MAXDOFPERNODE >> "$definefile"
  fi
fi


MAXGAUSS=`awk '
  /BEGIN/ { look=0 } 
  /^--*--/ { look=0 }
  /^--*--PROBLEM SIZE/ { look=1 }
  /^MAXGAUSS[ \t]*[0-9]*/ { if (look==1) { print $2 } }
' "$inputfile"`

if [ "x$MAXGAUSS" = x ]; then
  echo $0: Warning: MAXGAUSS undefined
else
  if grep MAXGAUSS "$definefile" 2>&1 > /dev/null ; then
    sed -e "s/MAXGAUSS=[0-9]*/MAXGAUSS=$MAXGAUSS/" "$definefile" > "$definefile".new
    mv "$definefile".new "$definefile"
  else
    echo MAXGAUSS=$MAXGAUSS >> "$definefile"
  fi
fi


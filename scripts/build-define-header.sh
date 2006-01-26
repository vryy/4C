
# get the define flags

mkdir -p $DEST/src/headers
cat > $DEST/src/headers/compile_settings.h <<EOF
/*----------------------------------------------------------------------*/
/*
 Created by $0 on `date`. Do not edit.
*/
/*----------------------------------------------------------------------*/

#ifndef COMPILE_SETTINGS_H
#define COMPILE_SETTINGS_H

#define DEFINE_STRING "`sed -e 's/#.*//' "$definefile" | awk '{ print $1 }' | sort | awk 'BEGIN { defs = "" }
END { print defs }
  { if (length($1) > 0) { defs = defs"\\\\n\\\\t"$1 } }' | sort`"

#define CREATOR "`whoami`@`hostname`"
#define CREATION_DATE "`date`"
#define CONFIGURATION "$configfile"

#endif

EOF

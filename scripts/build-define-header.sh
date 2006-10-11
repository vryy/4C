
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

#define DEFINE_STRING "`echo $DEFINES | sed -e 's/ /\n/g' -e 's/-D//g' | awk '{ print $1 }' | sort | awk 'BEGIN { defs = "" }
END { print defs }
  { if ( $NF ) { if ( $1 != "" ) { defs = defs"\\\\n\\\\t"$1 } } }' | sort`"

`echo $DEFINES | sed -e 's/ /\n/g' -e 's/=/ /g' -e 's/-D//g' | awk '{ print "/*#define " $0 " */" }'`

#define CREATOR "`whoami`@`hostname`"
#define CREATION_DATE "`date`"
#define CONFIGURATION "$configfile"

#endif

EOF

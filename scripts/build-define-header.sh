
# get the define flags

cat > src/headers/compile_settings.h <<EOF
/*----------------------------------------------------------------------*/
/*
 Created by $0 on `date`. Do not edit.
*/
/*----------------------------------------------------------------------*/

#define DEFINE_STRING "`sed -e 's/#.*//' config/defines.my | awk '{ print $1 }' | sort | awk 'BEGIN { defs = "" }
END { print defs }
  { if (length($1) > 0) { defs = defs"\\\\n\\\\t"$1 } }' | sort`"

EOF

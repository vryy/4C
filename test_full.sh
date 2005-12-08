#!/bin/sh
#==========================================
#
#                           (C) by mn 01/05
#==========================================


# check that there is at least one argument
if [ "x$1" = "x" ]; then
  echo "usage: $0  config-file-seq [config-file-par] "
  exit 1
fi


# check for existance and readability of first config file
if [ ! -s "$1" ]; then
  echo $0: file not found: $1
  exit 1
fi

if [ ! -r "$1" ]; then
  echo $0: file not readable: $1
  exit 1
fi


# check for existance and readability of second config file
if [ $# -eq 2 ]; then
  if [ ! -s "$2" ]; then
    echo $0: file not found: $2
    exit 1
  fi

  if [ ! -r "$2" ]; then
    echo $0: file not readable: $2
    exit 1
  fi
fi


out="test_result.txt"
echo > $out




echo "============================================================" >> $out
echo "Testing 'all' with '$1'" >> $out
echo "============================================================" >> $out
echo "" >> $out
echo "" >> $out
./run_test all $1 >> $out

echo "" >> $out
echo "" >> $out
echo "============================================================" >> $out
echo "Testing 'rel' with '$1'" >> $out
echo "============================================================" >> $out
echo "" >> $out
echo "" >> $out
./run_test rel $1 >> $out

if [ $# -eq 2 ]; then

  echo "" >> $out
  echo "" >> $out
  echo "============================================================" >> $out
  echo "Testing 'all' with '$2'" >> $out
  echo "============================================================" >> $out
  echo "" >> $out
  echo "" >> $out
  ./run_test all $2 >> $out

  echo "" >> $out
  echo "" >> $out
  echo "============================================================" >> $out
  echo "Testing 'rel' with '$2'" >> $out
  echo "============================================================" >> $out
  echo "" >> $out
  echo "" >> $out
  ./run_test rel $2 >> $out

fi



#!/bin/bash

bibinputs=$1
bufolder=$2
jobname=$3
for auxfile in `ls $bufolder/$jobname.*.aux`
do
    echo $bibinputs && bibtex ${auxfile%*.aux}
    $bibinputs && bibtex ${auxfile%*.aux}
done
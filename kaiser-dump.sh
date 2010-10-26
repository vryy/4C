#!/bin/sh

echo "bt" > gdb-cmds.txt
echo "q" >> gdb-cmds.txt

ps ax|grep [/]baci|awk '{print "gdb -x gdb-cmds.txt " $5 " " $1}'|bash

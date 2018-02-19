#!/bin/sh
cat $1 | sed 's/\r/\n/g' > ~/.tmp_file.txt
mv ~/.tmp_file.txt $1

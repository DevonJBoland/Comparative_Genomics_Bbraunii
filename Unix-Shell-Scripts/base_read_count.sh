#!/bin/bash

reads="insert read file name here"

zgrep . ${reads} |
     awk 'NR%4==2{c++; l+=length($0)}
          END{
                print "Number of reads: "c; 
                print "Number of bases in reads: "l
              }'
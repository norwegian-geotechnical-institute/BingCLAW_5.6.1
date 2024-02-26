#!/bin/sh
for fortfile in ../fort.q*
do
  sh ./gmt6_one_fortfile_2_png.sh  ${fortfile}
done

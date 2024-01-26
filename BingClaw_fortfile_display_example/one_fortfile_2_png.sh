#!/bin/sh
#
# Steven J Gibbons
# 2023-11-10
# Rough script to make a png output from a BingCLAW fort output file
# Requires GMT installed.
#
gmtbin=/opt/gmt/gmt4/bin
gmt5bin=/opt/gmt/gmt5/bin
gmt6bin=/opt/gmt/gmt6/bin
latmin=62.00
latmax=66.00
lonmin=-2.00
lonmax=7.00
proj=-JM12
reg=-R${lonmin}/${lonmax}/${latmin}/${latmax}
#
${gmtbin}/gmtset BASEMAP_TYPE FANCY ANOT_FONT_SIZE 12 LABEL_FONT_SIZE 12 \
       LABEL_FONT 1 ANOT_FONT 1 D_FORMAT %.12lg \
       PLOT_DEGREE_FORMAT +D
#
scriptname=./one_fortfile_2_png.sh
if [ $# != 1 ]
then
  echo
  echo "USAGE: "
  echo "$scriptname   ./fort.q0179  "
  echo
  exit 1
fi
#
fortfile=$1
if test ! -r ${fortfile}
then
  echo No file ${fortfile} found ...
  exit 1
fi
claw2file3=./claw2file3.py
if test ! -r ${claw2file3}
then
  echo No file ${claw2file3} found ...
  exit 1
fi
convertpy=./format3.py
if test ! -r ${convertpy}
then
  echo No file ${convertpy} found ...
  exit 1
fi
shortfile=`basename $fortfile`
stem=`echo $shortfile | sed 's/fort.//g'`
number=`echo $stem | sed 's/q//g'`
# minutes=`multiply $number 2.0 f6.1 `
minutes=`echo $number 2.0 | awk '{printf("%6.1f\n",$1*$2)}' `
timestring="Time: $minutes minutes"
ncfile=${stem}.nc
python ${claw2file3} -c1 -f ${fortfile} -s ${ncfile} 
mindepth=0.0
maxdepth=450.0
incdepth=25.0
elevinca=100.0
elevincf=50.0
if test ! -r depth.cpt
then
  ${gmtbin}/makecpt -Chot -T${mindepth}/${maxdepth}/${incdepth} -I > depth.cpt
fi
echo $stem
fil=${stem}.ps
${gmt6bin}/grdimage $ncfile -Cdepth.cpt $reg $proj -P -K -V              > $fil
${gmt6bin}/pscoast $reg $proj -Df -N1/2/0/0/0 -Ba2g2/a1g1/neWS  -V   \
              -G120/120/120               -O  -K   >> $fil
echo " -1.8   65.8 18  0 1 ML $timestring" | ${gmt6bin}/pstext $reg $proj -O -K >> $fil
${gmt6bin}/psscale -D0.65/2.5/4.0/0.6 -Cdepth.cpt  \
         -Ba${elevinca}f${elevincf}::/:"depth(m)":  -O >> $fil
/opt/gmt/gmt6/bin/psconvert $fil -TG -A
pngfile=${stem}.png
if test ! -d png_files
then
  mkdir png_files
fi
if test -r ${pngfile}
then
  mv ${pngfile} png_files
fi

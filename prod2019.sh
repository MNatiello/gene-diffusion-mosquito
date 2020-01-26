#!/bin/bash
DIR=/home/mario/Oxitec/code
#COU=1
#for rel in  52  260
#   do
#     for m in 0.50 0.75
#        do
#	  OXINAME=X-"$rel"-"$m"
#	  TOT=$((208 + $rel))
#          time ../v2/Oxi-FellerV2 $OXINAME $TOT 10 52 $rel $m 956287$COU
#	  COU=$(($COU+1))
#        done
#   done    

cd $DIR/prod

#
cp ../oldtests/v4.0-1/IC.dat.0 IC.dat
OXINAME=X-260-0.65-0
time ../oldtests/v4.0-1/Oxi-FellerV4 $OXINAME 10 156 260 156 0.0 0.65 5962871
../mymerge $OXINAME

cp ../oldtests/v4.0-1/IC.dat.0 IC.dat
OXINAME=X-260-0.65-12
time ../oldtests/v4.0-1/Oxi-FellerV4 $OXINAME 10 156 260 156 0.12 0.65 5962827
../mymerge $OXINAME

cp ../oldtests/v4.0-1/IC.dat.0 IC.dat
OXINAME=X-260-0.65-50
time ../oldtests/v4.0-1/Oxi-FellerV4 $OXINAME 10 156 260 156 0.50 0.65 5962838
../mymerge $OXINAME



gnuplot ../plots2019.gnuplot

../convert-to-pdf

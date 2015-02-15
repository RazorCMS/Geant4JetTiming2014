#!/bin/bash

fname=list_dustin_NoDuplicates.txt
fname2=list_dustin_BACON_NoDuplicates.txt

for i in `seq 1 3000`;
do
  nlines=`grep "g4simhitsEcal_"$i"_" $fname | wc -l`
  if [ $nlines -eq 1 ]; then
      aux=`grep "g4simhitsEcal_"$i"_" $fname`
      aux2=${aux#*Ecal}
      nlines2=`grep $aux2 $fname2 | wc -l`
      #match=`grep $aux2 $fname2` 
      if [ $nlines2 -gt 0 ]; then
          #echo $match
	  match=`grep $aux2 $fname2`
	  echo $match >> list_BACON_FINAL.txt
	  echo $aux >> list_g4Hits_FINAL.txt
      else
	  echo "no match"
      fi
  fi
done


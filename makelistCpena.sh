#!/bin/bash

ctr=0
ctr1=0
fname=list_dustin.txt
for i in `seq 1 5540`;
do
  #aux=`grep "default_data_HTMHT_Parked_RunB_v1_"$i"_" $fname | wc -l`
  #aux2=`grep -m 1 "default_data_HTMHT_Parked_RunB_v1_"$i"_" $fname`
  aux=`grep "g4simhitsEcal_"$i"_" $fname | wc -l`
  aux2=`grep -m 1 "g4simhitsEcal_"$i"_" $fname`
  #echo $aux2
  if [ $aux -gt 1 ];then
      echo $aux "  " $i
      ctr=$((ctr+1))
  fi
  
  if [ $aux -gt 0 ];then
      echo "$aux2" >> list_dustin_NoDuplicates.txt
      ctr1=$((ctr1+1))
  fi

done

echo "Duplicated Files:  "$ctr
echo "Unduplicated Files:  "$ctr1
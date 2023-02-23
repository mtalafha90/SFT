#!/bin/bash
#mkdir res_case1
#mkdir res_case2
#mkdir res_case3
# 7 x 15 x 11 = 1105 cases/flow type
for flow in 2
do
for u0 in 10
do
for eta in 250 
do
for tau in 8
do
for blat in 2.4
do
for bjoy in 0.15
do
python2.7 transp.py $u0 $eta $tau $flow $blat $bjoy
done
done
done
done
done
done
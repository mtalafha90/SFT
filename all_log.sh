#! /bin/bash
for u0 in 11
do
for eta in 250 450 600
do
for tau in 8 1000
do
python3 -W ignore all_log.py $u0 $eta $tau
done
done
done

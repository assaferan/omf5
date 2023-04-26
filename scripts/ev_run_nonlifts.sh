export LD_LIBRARY_PATH=/home/assaferan/lib

for i in {1..999}
do
/usr/bin/time -o logs/hecke_ev_3_0_nl_200_$i.time ./src/omf5 -genus=data/qf5db.sage -isom -nonlifts -disc=$i -prec=200 > data/hecke_ev_3_0_nl_200_$i.dat 2> logs/hecke_ev_3_0_nl_200_$i.log &
done

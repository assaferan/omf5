export LD_LIBRARY_PATH=/home/assaferan/lib

for i in {1..9}
do
./src/omf5 -genus=data/qf5db.sage -isom -disc=$i -prec=100 > data/hecke_ev_3_0_$i.dat 2> logs/hecke_ev_3_0_$i.log &
done

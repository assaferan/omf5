export LD_LIBRARY_PATH=/home/assaferan/lib

for i in {1..9}
do
./src/omf5 -genus=data/qf5.db -disc=$i -hecke > data/hecke_$i.dat 2> logs/hecke_$i.log &
done

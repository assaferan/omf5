export LD_LIBRARY_PATH=/home/assaferan/lib

./src/omf5 -format=GG -quad=1,1,1,0,0,1,1,1,0,1,0,1,0,0,1,0,0,8 -p=41 -hecke -cond=1 >& ./logs/omf5_61_41.log &
./src/omf5 -format=GG -quad=1,1,1,0,0,1,1,1,0,1,0,1,0,0,1,0,0,8 -p=101 -hecke -cond=1 >& ./logs/omf5_61_101.log &
./src/omf5 -format=GG -quad=1,1,1,0,0,1,1,1,0,1,0,1,0,0,1,0,0,8 -p=199 -hecke -cond=1 >& ./logs/omf5_61_199.log &
./src/omf5 -format=GG -quad=1,1,1,0,0,1,1,1,0,1,0,1,0,0,1,0,0,8 -p=293 -hecke -cond=1 >& ./logs/omf5_61_293.log &
./src/omf5 -format=GG -quad=1,1,1,0,0,1,1,1,0,1,0,1,0,0,1,0,0,8 -p=397 -hecke -cond=1 >& ./logs/omf5_61_397.log &

./src/omf5 -quad=1,1,0,1,0,1,0,0,0,1,1,0,1,1,34 -format=GG -p=41 -hecke -cond=167 >& ./logs/omf5_167_41.log &
./src/omf5 -quad=1,1,0,1,0,1,0,0,0,1,1,0,1,1,34 -format=GG -p=101 -hecke -cond=167 >& ./logs/omf5_167_101.log &
./src/omf5 -quad=1,1,0,1,0,1,0,0,0,1,1,0,1,1,34 -format=GG -p=199 -hecke -cond=167 >& ./logs/omf5_167_199.log &

./src/omf5 -quad=1,0,0,1,0,1,0,1,0,1,0,1,1,0,33 -format=GG -p=41 -hecke -cond=2 >& ./logs/omf5_262_41.log &
./src/omf5 -quad=1,0,0,1,0,1,0,1,0,1,0,1,1,0,33 -format=GG -p=101 -hecke -cond=2 >& ./logs/omf5_262_101.log &
./src/omf5 -quad=1,0,0,1,0,1,0,1,0,1,0,1,1,0,33 -format=GG -p=131 -hecke -cond=2 >& ./logs/omf5_262_131.log &

./src/omf5 -quad=1,0,0,1,0,1,0,1,0,1,1,1,2,0,17 -format=GG -p=41 -hecke -cond=334 >& ./logs/omf5_334_41.log &
./src/omf5 -quad=1,0,0,1,0,1,0,1,0,1,1,1,2,0,17 -format=GG -p=101 -hecke -cond=334 >& ./logs/omf5_334_101.log &
./src/omf5 -quad=1,0,0,1,0,1,0,1,0,1,1,1,2,0,17 -format=GG -p=199 -hecke -cond=334 >& ./logs/omf5_334_199.log &

./src/omf5 -quad=1,0,0,1,0,1,0,0,1,1,0,1,1,0,42 -format=GG -p=41 -hecke -cond=166 >& ./logs/omf5_498_41.log &
./src/omf5 -quad=1,0,0,1,0,1,0,0,1,1,0,1,1,0,42 -format=GG -p=101 -hecke -cond=166 >& ./logs/omf5_498_101.log &
./src/omf5 -quad=1,0,0,1,0,1,0,0,1,1,0,1,1,0,42 -format=GG -p=131 -hecke -cond=166 >& ./logs/omf5_498_131.log &
./src/omf5 -quad=1,0,0,1,0,1,0,0,1,1,0,1,1,0,42 -format=GG -p=199 -hecke -cond=166 >& ./logs/omf5_498_199.log &

./src/omf5 -quad=1,1,0,1,0,1,0,0,0,1,1,1,1,0,101 -format=GG -p=31 -hecke -cond=503 >& ./logs/omf5_503_31.log &
./src/omf5 -quad=1,1,0,1,0,1,0,0,0,1,1,1,1,0,101 -format=GG -p=41 -hecke -cond=503 >& ./logs/omf5_503_41.log &
./src/omf5 -quad=1,1,0,1,0,1,0,0,0,1,1,1,1,0,101 -format=GG -p=101 -hecke -cond=503 >& ./logs/omf5_503_101.log &

./src/omf5 -format=GG -quad=1,-1,1,-1,-1,2,-2,1,1,2,-2,-2,2,1,18 -p=31 -hecke -cond=1 >& ./logs/omf5_993_31.log &
./src/omf5 -format=GG -quad=1,-1,1,-1,-1,2,-2,1,1,2,-2,-2,2,1,18 -p=41 -hecke -cond=1 >& ./logs/omf5_993_41.log &
./src/omf5 -format=GG -quad=1,-1,1,-1,-1,2,-2,1,1,2,-2,-2,2,1,18 -p=101 -hecke -cond=1 >& ./logs/omf5_993_101.log &
./src/omf5 -format=GG -quad=1,-1,1,-1,-1,2,-2,1,1,2,-2,-2,2,1,18 -p=131 -hecke -cond=1 >& ./logs/omf5_993_131.log &
./src/omf5 -format=GG -quad=1,-1,1,-1,-1,2,-2,1,1,2,-2,-2,2,1,18 -p=199 -hecke -cond=1 >& ./logs/omf5_993_199.log &

./src/omf5 -quad=1,1,1,1,1,1,1,1,1,2,-1,0,11,7,11 -format=GG -hecke -p=31 -cond=1993 >& ./logs/omf5_1993_31.log &
./src/omf5 -quad=1,1,1,1,1,1,1,1,1,2,-1,0,11,7,11 -format=GG -hecke -p=41 -cond=1993 >& ./logs/omf5_1993_41.log &
./src/omf5 -quad=1,1,1,1,1,1,1,1,1,2,-1,0,11,7,11 -format=GG -hecke -p=101 -cond=1993 >& ./logs/omf5_1993_101.log &

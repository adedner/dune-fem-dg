#/bin/bash

rm -rf data

~/bin/replace "200" "50" -- ../grids/unitcube1.dgf
./euler2 fixedTimeStep:0.0024 > run2_50.out
grep "residual" run2_50.out > run2_50.gnu
mv data res2_50
mv run2_50.* res2_50

~/bin/replace "50" "100" -- ../grids/unitcube1.dgf
./euler2 fixedTimeStep:0.0012 > run2_100.out
grep "residual" run2_100.out > run2_100.gnu
mv data res2_100
mv run2_100.* res2_100

~/bin/replace "100" "200" -- ../grids/unitcube1.dgf
./euler2 fixedTimeStep:0.0006 > run2_200.out
grep "residual" run2_200.out > run2_200.gnu
mv data res2_200
mv run2_200.* res2_200

~/bin/replace "200" "400" -- ../grids/unitcube1.dgf
./euler2 fixedTimeStep:0.0003 > run2_400.out
grep "residual" run2_400.out > run2_400.gnu
mv data res2_400
mv run2_400.* res2_400

~/bin/replace "400" "800" -- ../grids/unitcube1.dgf
./euler2 fixedTimeStep:0.00015 > run2_800.out
grep "residual" run2_800.out > run2_800.gnu
mv data res2_800
mv run2_800.* res2_800

~/bin/replace "800" "200" -- ../grids/unitcube1.dgf

####

~/bin/replace "200" "50" -- ../grids/unitcube1.dgf
./euler3 fixedTimeStep:0.0024 > run3_50.out
grep "residual" run3_50.out > run3_50.gnu
mv data res3_50
mv run3_50.* res3_50

~/bin/replace "50" "100" -- ../grids/unitcube1.dgf
./euler3 fixedTimeStep:0.0012 > run3_100.out
grep "residual" run3_100.out > run3_100.gnu
mv data res3_100
mv run3_100.* res3_100

~/bin/replace "100" "200" -- ../grids/unitcube1.dgf
./euler3 fixedTimeStep:0.0006 > run3_200.out
grep "residual" run3_200.out > run3_200.gnu
mv data res3_200
mv run3_200.* res3_200

~/bin/replace "200" "400" -- ../grids/unitcube1.dgf
./euler3 fixedTimeStep:0.0003 > run3_400.out
grep "residual" run3_400.out > run3_400.gnu
mv data res3_400
mv run3_400.* res3_400

~/bin/replace "400" "800" -- ../grids/unitcube1.dgf
./euler3 fixedTimeStep:0.00015 > run3_800.out
grep "residual" run3_800.out > run3_800.gnu
mv data res3_800
mv run3_800.* res3_800

~/bin/replace "800" "200" -- ../grids/unitcube1.dgf

mv res?_* pressure_small

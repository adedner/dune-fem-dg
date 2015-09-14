#/bin/bash

./euler2 fixedTimeStep:0.005 > run2_200.out
grep "reisdual" run2_200.out > run2_200.gnu

~/bin/replace "200" "400" -- ../grids/unitcube1.dgf
./euler2 fixedTimeStep:0.0025 > run2_400.out
grep "reisdual" run2_400.out > run2_400.gnu

~/bin/replace "400" "800" -- ../grids/unitcube1.dgf
./euler2 fixedTimeStep:0.00125 > run2_800.out
grep "reisdual" run2_800.out > run2_800.gnu

~/bin/replace "800" "200" -- ../grids/unitcube1.dgf

####

cat ../grids/unitcube1.dgf > run3_200.out
./euler3 fixedTimeStep:0.005 >> run3_200.out
grep "reisdual" run3_200.out > run3_200.gnu

~/bin/replace "200" "400" -- ../grids/unitcube1.dgf
cat ../grids/unitcube1.dgf > run3_400.out
./euler3 fixedTimeStep:0.0025 >> run3_400.out
grep "reisdual" run3_400.out > run3_400.gnu

~/bin/replace "400" "800" -- ../grids/unitcube1.dgf
cat ../grids/unitcube1.dgf > run3_800.out
./euler3 fixedTimeStep:0.00125 >> run3_800.out
grep "reisdual" run3_800.out > run3_800.gnu

~/bin/replace "800" "200" -- ../grids/unitcube1.dgf



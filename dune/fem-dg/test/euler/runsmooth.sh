#/bin/bash
~/bin/replace "200" "25" -- ../grids/unitcube1.dgf
./euler0 fixedTimeStep:0.0024 > run0_25.out
grep "residual" run0_25.out > run0_25.gnu

~/bin/replace "25" "50" -- ../grids/unitcube1.dgf
./euler0 fixedTimeStep:0.0012 > run0_50.out
grep "residual" run0_50.out > run0_50.gnu

~/bin/replace "50" "100" -- ../grids/unitcube1.dgf
./euler0 fixedTimeStep:0.0006 > run0_100.out
grep "residual" run0_100.out > run0_100.gnu

~/bin/replace "100" "200" -- ../grids/unitcube1.dgf
./euler0 fixedTimeStep:0.0003 > run0_200.out
grep "residual" run0_200.out > run0_200.gnu

~/bin/replace "200" "400" -- ../grids/unitcube1.dgf
./euler0 fixedTimeStep:0.00015 > run0_400.out
grep "residual" run0_400.out > run0_400.gnu

~/bin/replace "400" "800" -- ../grids/unitcube1.dgf
./euler0 fixedTimeStep:0.000075 > run0_800.out
grep "residual" run0_800.out > run0_800.gnu

~/bin/replace "800" "1600" -- ../grids/unitcube1.dgf
./euler0 fixedTimeStep:0.0000375 > run0_1600.out
grep "residual" run0_1600.out > run0_1600.gnu

~/bin/replace "1600" "200" -- ../grids/unitcube1.dgf

tail -qn 1 run0_25.gnu run0_50.gnu run0_100.gnu run0_200.gnu run0_400.gnu run0_800.gnu run0_1600.gnu > run0.gnu

##

~/bin/replace "200" "25" -- ../grids/unitcube1.dgf
./euler1 fixedTimeStep:0.0024 > run1_25.out
grep "residual" run1_25.out > run1_25.gnu

~/bin/replace "25" "50" -- ../grids/unitcube1.dgf
./euler1 fixedTimeStep:0.0012 > run1_50.out
grep "residual" run1_50.out > run1_50.gnu

~/bin/replace "50" "100" -- ../grids/unitcube1.dgf
./euler1 fixedTimeStep:0.0006 > run1_100.out
grep "residual" run1_100.out > run1_100.gnu

~/bin/replace "100" "200" -- ../grids/unitcube1.dgf
./euler1 fixedTimeStep:0.0003 > run1_200.out
grep "residual" run1_200.out > run1_200.gnu

~/bin/replace "200" "400" -- ../grids/unitcube1.dgf
./euler1 fixedTimeStep:0.00015 > run1_400.out
grep "residual" run1_400.out > run1_400.gnu

~/bin/replace "400" "800" -- ../grids/unitcube1.dgf
./euler1 fixedTimeStep:0.000075 > run1_800.out
grep "residual" run1_800.out > run1_800.gnu

~/bin/replace "800" "200" -- ../grids/unitcube1.dgf

tail -qn 1 run1_25.gnu run1_50.gnu run1_100.gnu run1_200.gnu run1_400.gnu run1_800.gnu > run1.gnu

####

~/bin/replace "200" "25" -- ../grids/unitcube1.dgf
./euler2 fixedTimeStep:0.0024 > run2_25.out
grep "residual" run2_25.out > run2_25.gnu

~/bin/replace "25" "50" -- ../grids/unitcube1.dgf
./euler2 fixedTimeStep:0.0012 > run2_50.out
grep "residual" run2_50.out > run2_50.gnu

~/bin/replace "50" "100" -- ../grids/unitcube1.dgf
./euler2 fixedTimeStep:0.0006 > run2_100.out
grep "residual" run2_100.out > run2_100.gnu

~/bin/replace "100" "200" -- ../grids/unitcube1.dgf
./euler2 fixedTimeStep:0.0003 > run2_200.out
grep "residual" run2_200.out > run2_200.gnu

~/bin/replace "200" "400" -- ../grids/unitcube1.dgf
./euler2 fixedTimeStep:0.00015 > run2_400.out
grep "residual" run2_400.out > run2_400.gnu

~/bin/replace "400" "800" -- ../grids/unitcube1.dgf
./euler2 fixedTimeStep:0.000075 > run2_800.out
grep "residual" run2_800.out > run2_800.gnu

~/bin/replace "800" "200" -- ../grids/unitcube1.dgf

tail -qn 1 run2_25.gnu run2_50.gnu run2_100.gnu run2_200.gnu run2_400.gnu run2_800.gnu > run2.gnu

####

~/bin/replace "200" "25" -- ../grids/unitcube1.dgf
./euler3 fixedTimeStep:0.0024 > run3_25.out
grep "residual" run3_25.out > run3_25.gnu

~/bin/replace "25" "50" -- ../grids/unitcube1.dgf
./euler3 fixedTimeStep:0.0012 > run3_50.out
grep "residual" run3_50.out > run3_50.gnu

~/bin/replace "50" "100" -- ../grids/unitcube1.dgf
./euler3 fixedTimeStep:0.0006 > run3_100.out
grep "residual" run3_100.out > run3_100.gnu

~/bin/replace "100" "200" -- ../grids/unitcube1.dgf
./euler3 fixedTimeStep:0.0003 > run3_200.out
grep "residual" run3_200.out > run3_200.gnu

~/bin/replace "200" "400" -- ../grids/unitcube1.dgf
./euler3 fixedTimeStep:0.00015 > run3_400.out
grep "residual" run3_400.out > run3_400.gnu

~/bin/replace "400" "200" -- ../grids/unitcube1.dgf

tail -qn 1 run3_25.gnu run3_50.gnu run3_100.gnu run3_200.gnu run3_400.gnu > run3.gnu

####

~/bin/replace "200" "25" -- ../grids/unitcube1.dgf
./euler4 fixedTimeStep:0.0024 > run4_25.out
grep "residual" run4_25.out > run4_25.gnu

~/bin/replace "25" "50" -- ../grids/unitcube1.dgf
./euler4 fixedTimeStep:0.0012 > run4_50.out
grep "residual" run4_50.out > run4_50.gnu

~/bin/replace "50" "100" -- ../grids/unitcube1.dgf
./euler4 fixedTimeStep:0.0006 > run4_100.out
grep "residual" run4_100.out > run4_100.gnu

~/bin/replace "100" "200" -- ../grids/unitcube1.dgf
./euler4 fixedTimeStep:0.0003 > run4_200.out
grep "residual" run4_200.out > run4_200.gnu

~/bin/replace "200" "400" -- ../grids/unitcube1.dgf
./euler4 fixedTimeStep:0.00015 > run4_400.out
grep "residual" run4_400.out > run4_400.gnu

~/bin/replace "400" "200" -- ../grids/unitcube1.dgf

tail -qn 1 run4_25.gnu run4_50.gnu run4_100.gnu run4_200.gnu run4_400.gnu > run4.gnu


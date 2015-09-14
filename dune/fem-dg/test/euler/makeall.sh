make clean
replace "POLORDER=2" "POLORDER=0" -- Makefile.am
make euleronep ; mv euleronep euler0
make clean
replace "POLORDER=0" "POLORDER=1" -- Makefile.am
make euleronep ; mv euleronep euler1
make clean
replace "POLORDER=1" "POLORDER=2" -- Makefile.am
make euleronep ; mv euleronep euler2
make clean
replace "POLORDER=2" "POLORDER=3" -- Makefile.am
make euleronep ; mv euleronep euler3
make clean
replace "POLORDER=3" "POLORDER=4" -- Makefile.am
make euleronep ; mv euleronep euler4
make clean

replace "POLORDER=4" "POLORDER=2" -- Makefile.am


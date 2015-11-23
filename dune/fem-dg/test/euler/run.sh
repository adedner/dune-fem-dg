#/bin/bash

# echo $1
name="sodOrig"
echo $name

# rm -rf $name
rm -rf data

mkdir $name

mkdir data
./euler0 femhowto.eocSteps:5 fixedTimeStep:0.0024 >& data/run0.out
mv data $name/p0
mkdir data
./euler1 femhowto.eocSteps:5 fixedTimeStep:0.0024 >& data/run0.out
mv data $name/p1
mkdir data
./euler2 femhowto.eocSteps:5 fixedTimeStep:0.0024 >& data/run0.out
mv data $name/p2
mkdir data
./euler3 femhowto.eocSteps:5 fixedTimeStep:0.0024 >& data/run0.out
mv data $name/p3
mkdir data
./euler4 femhowto.eocSteps:4 fixedTimeStep:0.0024 >& data/run0.out
mv data $name/p4
mkdir data

cp euler? parameter $name

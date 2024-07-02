#!/usr/bin/bash

## experiment to calculate the energies for  
# make clean > terminal.txt
# make >> terminal.txt

# cd run/
# for i in {1..5}
# do
#     sed -i "s/Posfile = second\/POSCAR$((i-1))/Posfile = second\/POSCAR$i/g" input.in
#     ./coulomb.x
# done
# echo "File,Reciprocal Energy,Reciprocal Energy Integral,Error,Reciprocal Energy FFT,Error" > second/second.csv
# for i in {1..5}
# do
#     sed -i "s/Posfile = second\/POSCAR$((i-1))/Posfile = second\/POSCAR$i/g" input.in
#     echo -n "POSCAR"$i >> second/second.csv
#     ./coulomb.x >> second/second.csv
#     echo " " >> second/second.csv
# done
# cd ..
# rm terminal.txt

###################################################################################

## code for multiple cells and same threads
# make clean > terminal.txt
# make >> terminal.txt

# cd run/
# echo "File,Reciprocal,Time,Reciprocal(Integral),Time,error,%Reduction in Time,Real,Time" > exp1/r1.csv

# for i in {1..14} 
# do
#     sed -i "s/Posfile = exp1\/POSCAR\.$((i-1))/Posfile = exp1\/POSCAR\.$i/g" input.in
#     echo -n "POSCAR"$i >> exp1/r1.csv
#     for j in {1..5}
#     do
#         ./coulomb.x >> exp1/r1.csv
#         echo " " >> exp1/r1.csv
#     done
#     echo " " >> exp1/r1.csv
#     echo " " >> exp1/r1.csv
# done
# cd ..
# rm terminal.txt

###################################################################################

## Same system and varying threads
cd run/
echo "# of Threads,Reciprocal PHL,Reciprocal Integral" > exp2/r4.csv
cd ..
for i in {1..24}
do
    sed -i "s/\#define NUM_THREADS $((i-1))/\#define NUM_THREADS $i/g" inc/const.h
    make clean > output.txt
    make >> output.txt
    rm output.txt
    cd run/
    echo -n $i >> exp2/r4.csv
    for j in {1..5}
    do
        ./coulomb.x >> exp2/r4.csv
        echo " " >> exp2/r4.csv
    done
    echo " " >> exp2/r4.csv
    echo " " >> exp2/r4.csv
    cd ..
done
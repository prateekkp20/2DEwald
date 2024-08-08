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
# echo "File,Reciprocal,Time,Reciprocal(Integral),Time,error,%Reduction in Time(For Integral),Reciprocal(3DFFT),Time,error,Real,Time,Total" > exp1box45/R1.csv

# for i in {1..9} 
# do
#     sed -i "s/Posfile = exp1box45\/POSCAR\.$((i-1))/Posfile = exp1box45\/POSCAR\.$i/g" input.in
#     echo -n "POSCAR"$i >> exp1box45/R1.csv
#     for j in {1..5}
#     do
#         ./coulomb.x >> exp1box45/R1.csv
#         echo " " >> exp1box45/R1.csv
#     done
#     echo " " >> exp1box45/R1.csv
#     echo " " >> exp1box45/R1.csv
# done
# cd ..
# rm terminal.txt

###################################################################################

# ## Same system and varying threads
# cd run/
# echo "# of Threads,Reciprocal PHL,Reciprocal Integral" > exp2/r4.csv
# cd ..
# for i in {1..24}
# do
#     sed -i "s/\#define NUM_THREADS $((i-1))/\#define NUM_THREADS $i/g" inc/const.h
#     make clean > output.txt
#     make >> output.txt
#     rm output.txt
#     cd run/
#     echo -n $i >> exp2/r4.csv
#     for j in {1..5}
#     do
#         ./coulomb.x >> exp2/r4.csv
#         echo " " >> exp2/r4.csv
#     done
#     echo " " >> exp2/r4.csv
#     echo " " >> exp2/r4.csv
#     cd ..
# done

######################################################

## code for multiple cells and same threads only of the real energy calculations for the folder exp1

# make clean > terminal.txt
# make >> terminal.txt

# cd run/
# echo "File,Real, Time," > exp1/RealV2D.csv

# for i in {1..14} 
# do
#     sed -i "s/Posfile = exp1\/POSCAR\.$((i-1))/Posfile = exp1\/POSCAR\.$i/g" input.in
#     echo -n "POSCAR"$i >> exp1/RealV2D.csv
#     for j in {1..5}
#     do
#         ./coulomb.x >> exp1/RealV2D.csv
#         echo " " >> exp1/RealV2D.csv
#     done
#     echo " " >> exp1/RealV2D.csv
#     echo " " >> exp1/RealV2D.csv
# done
# cd ..
# rm terminal.txt

####################################################
# New bench marking, 06 Aug task 1

make clean > terminal.txt
make >> terminal.txt
locationcsv="separate_p_m/z/out.csv"
locationinput="separate_p_m\/z"

cd run/
echo "File,E_long,E_coul,Total" > $locationcsv
lastfile=4
for i in $(seq 1 $lastfile)
do
    sed -i "s/Posfile = ${locationinput}\/POSCAR\.$((i-1))/Posfile = ${locationinput}\/POSCAR\.$i/g" input.in
    echo -n "POSCAR"$i >> $locationcsv
    ./coulomb.x >> $locationcsv
    echo " " >> $locationcsv
    if [ $i == $lastfile ];then
        sed -i "s/Posfile = ${locationinput}\/POSCAR\.$i/Posfile = ${locationinput}\/POSCAR\.0/g" input.in
    fi
done
cd ..
rm terminal.txt
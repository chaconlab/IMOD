#!/bin/bash
echo "COMPILING & LINKING ALL iMOD v1.04 TOOLS"
export IMOD_FLAGS=""
echo
echo "Entering: src"
cd src

# COMPILING LIBRARIES AND EXECUTABLES
 for i in libtools libpdb libvolume libnma libnmafit libgausscorr imodview imove imc imode nmacm nmafit
 do
  echo "  Entering: $i"
  cd $i
   for j in *
   do
    if test -d $j
    then
     echo "    Entering: $j"
     cd $j
      echo Cleaning $i $j
      make clean
      echo Making $i $j
      make all
     cd ..
     echo "    Exiting: $j"
    fi
   done
  cd ..
  echo "  Exiting: $i"
 done

cd ..
echo "Exiting: src"

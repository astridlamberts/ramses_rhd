#! /bin/tcsh
echo $WORKDIR
cd $WORKDIR
mkdir tmp_distrib
cd tmp_distrib

cp -r $HOME/svnramses/trunk ramses
cd ramses
rm -rf .svn
rm -rf */.svn
rm -rf */*/.svn
rm -rf */*/*/.svn
rm -rf */*/*/*/.svn

rm -rf *~
rm -rf */*~
rm -rf */*/*~
rm -rf */*/*/*~
rm -rf */*/*/*~/*~

rm -rf diffusion
rm -rf multimat
rm -rf stiff
rm -rf grafic2

cd patch
rm -rf diffusion
rm -rf induction

cd ../doc
cp src/ramses_ug.pdf .
rm -rf src

cd ../namelist
rm -rf *_mmat.nml
rm -rf *_stiff.nml
rm -rf *_diff.nml
rm -rf alloy.nml beltrami.nml ponomarenko.nml implosion.nml spitzer.nml

cd ../bin
rm -rf *.o
rm -rf *.mod
rm -rf *1d
rm -rf *2d
rm -rf *3d
rm -rf *_grafic
rm -rf histo
rm -rf map_gas*
rm -rf map_mass*
rm -rf amrdir
rm -rf sod
rm -rf amr2map
rm -rf amr2cube
rm -rf amr2gmsh
rm -rf part2map
rm -rf part2cube
rm -rf header
rm -rf log2col
mv Makefile.ramses Makefile
rm -rf Makefile.*

cd ../utils/f90
rm -rf *_grafic
rm -rf histo
rm -rf map_gas*
rm -rf map_mass*
rm -rf amr2map
rm -rf amrdir
rm -rf sod
rm -rf amr2cube
rm -rf amr2gmsh
rm -rf part2map
rm -rf part2cube
rm -rf header
rm -rf log2col
rm -rf hop_ramses/hop
rm -rf hop_ramses/regroup
rm -rf hop_ramses/poshalo
rm -rf hop_ramses/*.o

cd ../idl
rm -rf smooth

cd ../../..

tar cvf ramses.tar ramses
gzip ramses.tar

echo 'ramses package done' 

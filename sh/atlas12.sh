#!/bin/sh

date
#rm -f fort.*
ln -s refmodel_input_pyatlas.py fort.3
ln -s ../../atlas12/molecules.dat fort.2
ln -s ../../atlas12/fclowlines.bin fort.11
ln -s ../../atlas12/nltelines.bin fort.19
ln -s ../../atlas12/fchighlines.bin fort.21
ln -s ../../atlas12/diatomicsiwl.bin fort.31
ln -s ../../atlas12/schwenke.bin fort.41 
ln -s ../../atlas12/h2ofastfix.bin fort.51
# ln -s ../../bin/atlas12.exe


echo "Running Atlas12 - 1..."
./../../bin/atlas12.exe<<EOF >run1.log
MOLECULES ON
READ MOLECULES
READ PUNCH
READ LINES
CONVECTION OVER 1.25 0 36
ITERATIONS 1 PRINT 1 PUNCH 0
ABUNDANCES INPUT FROM FILE
BEGIN
END
EOF
# date
echo "Running Atlas12 - 2..."
./../../bin/atlas12.exe<<EOF >run2.log
MOLECULES ON
READ MOLECULES
READ PUNCH
TITLE ATLAS12 
OPACITY ON LINES
OPACITY ON XLINES
CONVECTION OVER 1.25 0 36
ITERATIONS 45
PRINT 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
PUNCH 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
SCALE MODEL 72 -6.875 0.125 5777. 4.43770
VTURB 1.0E+05
ABUNDANCES INPUT FROM FILE
BEGIN 
END
EOF

mv fort.7  model.dat 
#rm fort.*
date
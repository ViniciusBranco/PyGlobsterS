#!/bin/sh
date
echo "Retrieving files"
ln -s ../../bin/* .
ln -s ../../gfs/* .
ln -s ../../models/stars/* .   
ln -s ../../molecules/* .
ln -s ../../abunds/NGC0104/* .

echo "Setting the model"
cat flux.header > synthe.mod
head -n 3 t3800_g1.0_m07p04_asp9.mod >> synthe.mod
cat conv.txt abundances_47tuc_asp9_m07p04_2020.txt >> synthe.mod
tail -n 75 t3800_g1.0_m07p04_asp9.mod >> synthe.mod

echo "Reading Molecules and Continua"
ln -s molecules.dat fort.2
ln -s continua.dat fort.17
ln -s he1tables.dat fort.18

./xnfpelsyn.exe < synthe.mod > /dev/null

echo "Synbeg input"
./synbeg.exe << "EOF" > /dev/null
AIR       500.0     510.0     300000.   0.     0         -10 .001         0   00
AIRorVAC  WLBEG     WLEND     RESOLU    TURBV  IFNLTE LINOUT CUTOFF        NREAD
EOF

echo "Reading atomic lines"
ln -s gfallkurucz17.dat fort.11
./rgfalllinesnew.exe > /dev/null
rm fort.11

echo "Running SYNTHE"
./synthe.exe > /dev/null

cat <<"EOF" > fort.25
0.0       0.        1.        0.        0.        0.        0.        0.
0.
RHOXJ     R1        R101      PH1       PC1       PSI1      PRDDOP    PRDPOW
EOF

echo "Running SPECTRV"
./spectrv.exe< synthe.mod > /dev/null
mv fort.7 spec.dat

echo "Rotating Spectrum"
ln -s spec.dat fort.1
./rotate.exe<< "EOF" > /dev/null
    1
4.
EOF
mv ROT1 spec.dat

echo "Broadening Spectrum"
ln -s spec.dat fort.21
./broaden.exe<< "EOF" > /dev/null
GAUSSIAN  300000.    RESOLUTION
EOF
mv fort.21 spec.bin
rm fort.*

echo "Converting Spectrum to ASC"
ln -s spec.bin fort.1
./converfsynnmtoa.exe > /dev/null
mv fort.2 spec.asc

mv spec.asc ../libraries/synthetic/stars
date
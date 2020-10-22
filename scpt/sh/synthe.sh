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

echo "Reading molecule lines"

## AlH
    # echo "alhax..."
    ln -s AlH/alhax.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

    # echo "alhxx..."
    ln -s AlH/alhxx.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

## AlO
    # echo "alopatrascu..."
    ln -s AlO/alopatrascu.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

## C2
    # echo "c2ax..."
    ln -s C2/c2ax.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

    # echo "c2ba..."
    ln -s C2/c2ba.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

    # echo "c2dabrookek..."
    ln -s C2/c2dabrookek.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

    # echo "c2ea..."
    ln -s C2/c2ea.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

## CaH
    # echo "cah..."
    ln -s CaH/cah.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

## CaO
    # echo "caoyurchenko..."
    ln -s CaO/caoyurchenko.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

## CH
    # echo "ch..."
    ln -s CH/chmasseron.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

## CN
    # echo "cnaxbrookek..."
    ln -s CN/cnaxbrookek.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

    # echo "cnbxbrookek..."
    ln -s CN/cnbxbrookek.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

    # echo "cnxx12brooke..."
    ln -s CN/cnxx12brooke.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

## CO
    # echo "coax..."
    ln -s CO/coax.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

    # echo "coxx..."
    ln -s CO/coxx.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

## CrH
    # echo "crhaxbernath..."
    ln -s CrH/crhaxbernath.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

## FeH
    # echo "fehfx..."
    ln -s FeH/fehfx.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

## H2
    # echo "h2..."
    ln -s H2/h2.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

## MgH
    # echo "mgh..."
    ln -s MgH/mgh.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

## MgO
    # echo "mgodaily..."
    ln -s MgO/mgodaily.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

## NaH
    # echo "nahrivlin..."
    ln -s NaH/nahrivlin.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

## NH
    # echo "nh..."
    ln -s NH/nh.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

## OH
    # echo "ohupdate..."
    ln -s OH/ohupdate.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

## SiH
    # echo "sihnew..."
    ln -s SiH/sihnew.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

## SiO
    # echo "sioax..."
    ln -s SiO/sioax.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

    # echo "sioex..."
    ln -s SiO/sioex.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

    # echo "sioxx..."
    ln -s SiO/sioxx.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

## TiH
    # echo "tih..."
    ln -s TiH/tih.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

## VO
    # echo "vo..."
    ln -s VO/vo.asc fort.11
     ./rmolecasc.exe >  /dev/null
    rm -f fort.11

# echo "Reading tioschwenke.bin molecule"
ln -s TiO/tioschwenke.bin fort.11
# echo "Reading eschwenke.bin"
ln -s TiO/eschwenke.bin fort.48
./rschwenk.exe > /dev/null
rm fort.11
rm fort.48

# echo "Reading h2ofastfix.bin"
# ln -s H2O/h2ofastfix.bin fort.11
# ./rh2ofast.exe > /dev/null
# rm fort.11

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
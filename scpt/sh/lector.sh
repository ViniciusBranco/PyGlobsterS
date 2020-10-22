#!/bin/sh
date
ln -s ../LECTOR/* .
ln -s ../libraries/synthetic/gc/NGC1904/convwag_synthe_NGC1904_AlphaCNONaHe.asc

echo "LECTOR"
./LECTOR << "EOF" > output.txt
s
convwag_synthe_NGC1904_AlphaCNONaHe.asc
0
n
EOF

mv output.txt ../libraries/synthetic/gc/NGC1904/
mv convwag_synthe_NGC1904_AlphaCNONaHe.asc_INDICES ../libraries/synthetic/gc/NGC1904/LECTOR
mv convwag_synthe_NGC1904_AlphaCNONaHe.asc_LINE ../libraries/synthetic/gc/NGC1904/LECTOR
mv convwag_synthe_NGC1904_AlphaCNONaHe.asc_ROSE ../libraries/synthetic/gc/NGC1904/LECTOR

find -type l -delete
date

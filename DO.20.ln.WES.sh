#!/bin/sh
export RAW=/DATA/CUK/Somatic_Copy_Number_Alteration/CHET_58
rm -rf CRC4
mkdir CRC4
cd CRC4
ln -s ${RAW}/CHET_58_1D.ba[im] .
ln -s ${RAW}/CHET_58_ND.ba[im] .
cd ..

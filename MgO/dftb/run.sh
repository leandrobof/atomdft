cp bulk2 dftb_in.hsd
./dftb+
cp bulk dftb_in.hsd
./dftb+
dp_bands band.out mgo_band
efermi=`grep "Fermi energy" detailed.out | awk '{print $5}' `
python plot2.py $efermi

rm -rf build
mkdir build
cd build 
cmake ..
make -j8
cd bin
cp ../../data/network/*.yaml .
cp  ../../tests/*.inp .
cp ../../tests/*txt .
cp ../../data/planet/*txt .
cp ../../data/stellar/sun_spec.inp .
cp -f ../../data/network/*yaml /data4/ananyo/opt/share/cantera/data/.
export CANTERA_DATA=/data4/ananyo/opt/share/cantera/data/.

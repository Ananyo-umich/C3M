rm -rf build
mkdir build
cd build 
cmake ..
make -j8
cd bin
cp ../../tests/*.yaml .
cp  ../../tests/*.inp .
cp  ../../tests/*.csv .
cp -f ../../tests/*yaml /data4/ananyo/opt/share/cantera/data/.
export CANTERA_DATA=/data4/ananyo/opt/share/cantera/data/.

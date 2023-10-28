rm -rf build
mkdir build
cd build 
cmake ..
make -j8
cd bin
cp ../../tests/*.yaml .
cp  ../../tests/*.inp .
cp  ../../tests/*.csv .

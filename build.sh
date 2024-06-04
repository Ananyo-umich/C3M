rm -rf build
mkdir build
cd build
cmake ..
make -j4
cd tests
ln -s ../../data/KINETICS7 KINETICS7
ln -s ../../data/VULCAN VULCAN
ln -s ../../data/planet planet

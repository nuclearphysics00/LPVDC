#g++ fast_ion_signal.C mdt_gamma.C mdt_mt.C mdt.C \
#g++ mdt_gain.C \
g++ -std=c++17 -O2 -g alpha.cc -DGEOM_WIRE \
  -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
  -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
  -lGarfield $(root-config --glibs) \
  -o run

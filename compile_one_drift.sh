
# wire 版
g++ -std=c++17 -O2 -g one_drift.cc -DGEOM_WIRE \
  -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
  -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
  -lGarfield $(root-config --glibs) \
  -o run_one_drift

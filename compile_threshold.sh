# pap
# plate 版
g++ -std=c++17 -O2 -g main_fieldview_pap_threshold.cc -DGEOM_PLATE \
  -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
  -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
  -lGarfield $(root-config --glibs) \
  -o fieldview_pap_plate_threshold

# wire 版
g++ -std=c++17 -O2 -g main_fieldview_pap_threshold.cc -DGEOM_WIRE \
  -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
  -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
  -lGarfield $(root-config --glibs) \
  -o fieldview_pap_wirecath_threshold

# ap
# plate 版
g++ -std=c++17 -O2 -g main_fieldview_ap_threshold.cc -DGEOM_PLATE \
  -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
  -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
  -lGarfield $(root-config --glibs) \
  -o fieldview_ap_plate_threshold

# wire 版
g++ -std=c++17 -O2 -g main_fieldview_ap_threshold.cc -DGEOM_WIRE \
  -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
  -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
  -lGarfield $(root-config --glibs) \
  -o fieldview_ap_wirecath_threshold
# pap
# plate 版
g++ -std=c++17 -O2 -g main_fieldview_pap.cc -DGEOM_PLATE \
  -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
  -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
  -lGarfield $(root-config --glibs) \
  -o fieldview_pap_plate

# wire 版
g++ -std=c++17 -O2 -g main_fieldview_pap.cc -DGEOM_WIRE \
  -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
  -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
  -lGarfield $(root-config --glibs) \
  -o fieldview_pap_wirecath

# ap
# plate 版
g++ -std=c++17 -O2 -g main_fieldview_ap.cc -DGEOM_PLATE \
  -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
  -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
  -lGarfield $(root-config --glibs) \
  -o fieldview_ap_plate

# wire 版
g++ -std=c++17 -O2 -g main_fieldview_ap.cc -DGEOM_WIRE \
  -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
  -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
  -lGarfield $(root-config --glibs) \
  -o fieldview_ap_wirecath

# aa
  # plate 版
g++ -std=c++17 -O2 -g main_fieldview_a.cc -DGEOM_PLATE \
  -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
  -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
  -lGarfield $(root-config --glibs) \
  -o fieldview_a_plate

# wire 版
g++ -std=c++17 -O2 -g main_fieldview_a.cc -DGEOM_WIRE \
  -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
  -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
  -lGarfield $(root-config --glibs) \
  -o fieldview_a_wirecath
  ## plate 版
  g++ -std=c++17 -O2 -g optimize_efield_pap.cc -DGEOM_PLATE \
     -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
     -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
     -lGarfield $(root-config --glibs) \
     -o optimize_efield_pap_plate

  ## wire 版
  g++ -std=c++17 -O2 -g optimize_efield_pap.cc -DGEOM_WIRE \
     -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
     -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
     -lGarfield $(root-config --glibs) \
     -o optimize_efield_pap_wirecath

  ## plate 版
  g++ -std=c++17 -O2 -g optimize_efield_ap.cc -DGEOM_PLATE \
     -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
     -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
     -lGarfield $(root-config --glibs) \
     -o optimize_efield_ap_plate

  ## wire 版
  g++ -std=c++17 -O2 -g optimize_efield_ap.cc -DGEOM_WIRE \
     -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
     -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
     -lGarfield $(root-config --glibs) \
     -o optimize_efield_ap_wirecath

  ## plate 版
  g++ -std=c++17 -O2 -g optimize_efield_a.cc -DGEOM_PLATE \
     -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
     -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
     -lGarfield $(root-config --glibs) \
     -o optimize_efield_a_plate

  ## wire 版
  g++ -std=c++17 -O2 -g optimize_efield_a.cc -DGEOM_WIRE \
     -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
     -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
     -lGarfield $(root-config --glibs) \
     -o optimize_efield_a_wirecath
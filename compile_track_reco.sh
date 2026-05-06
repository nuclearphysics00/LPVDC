   ## nomal
   ## plate 版
   g++ -std=c++17 -O2 -g track_reco_avalanche.cc -DGEOM_PLATE \
      -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
      -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
      -lGarfield $(root-config --glibs) \
      -o track_reco_avalanche_plate
 
   ## wire 版
   g++ -std=c++17 -O2 -g track_reco_avalanche.cc -DGEOM_WIRE \
      -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
      -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
      -lGarfield $(root-config --glibs) \
      -o track_reco_avalanche_wirecath
 
   ## pap
   ## plate 版
   g++ -std=c++17 -O2 -g track_reco_avalanche_pap.cc -DGEOM_PLATE \
      -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
      -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
      -lGarfield $(root-config --glibs) \
      -o track_reco_avalanche_pap_plate
 
   ## wire 版
   g++ -std=c++17 -O2 -g track_reco_avalanche_pap.cc -DGEOM_WIRE \
      -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
      -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
      -lGarfield $(root-config --glibs) \
      -o track_reco_avalanche_pap_wirecath
 
   ## ap
   ## plate 版
   g++ -std=c++17 -O2 -g track_reco_avalanche_ap.cc -DGEOM_PLATE \
      -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
      -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
      -lGarfield $(root-config --glibs) \
      -o track_reco_avalanche_ap_plate
 
   ## wire 版
   g++ -std=c++17 -O2 -g track_reco_avalanche_ap.cc -DGEOM_WIRE \
      -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
      -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
      -lGarfield $(root-config --glibs) \
      -o track_reco_avalanche_ap_wirecath
 
   ## a
   ## plate 版
   g++ -std=c++17 -O2 -g track_reco_avalanche_a.cc -DGEOM_PLATE \
      -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
      -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
      -lGarfield $(root-config --glibs) \
      -o track_reco_avalanche_a_plate
 
   ## wire 版
   g++ -std=c++17 -O2 -g track_reco_avalanche_a.cc -DGEOM_WIRE \
      -I"$GARFIELD_HOME/include" -I. -Igeometry $(root-config --cflags) \
      -L"$GARFIELD_HOME/install/lib" -Wl,-rpath,"$GARFIELD_HOME/install/lib" \
      -lGarfield $(root-config --glibs) \
      -o track_reco_avalanche_a_wirecath
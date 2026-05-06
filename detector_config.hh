// detector_config.hh
#pragma once
#include "Garfield/MediumMagboltz.hh"
#include "detector_geometry.hh"
#include <string>
#include <cstdlib>
#include <iostream>

namespace Config {

// ガスをロード
inline bool LoadGas(Garfield::MediumMagboltz& gas,
                    const std::string& gasfile="ic4H10_100_0.1Torr.gas",
                    double torr=0.0, double atm=0.0) {
  using namespace Garfield;
  if (!gas.LoadGasFile(gasfile)) {
    std::cerr << "Gas load failed: " << gasfile << std::endl;
    return false;
  }
  if (torr > 0) gas.SetPressure(torr / 750.061683);
  if (atm  > 0) gas.SetPressure(atm * 1.01325);
  return true;
}

// 電圧一括上書き
inline void OverrideVoltages(Detector::Geometry& geo,
                             const std::vector<std::pair<std::string,double>>& list) {
  for (auto& e : geo.electrodes)
    for (auto& kv : list)
      if (e.name == kv.first) e.V = kv.second;
}

} // namespace Config

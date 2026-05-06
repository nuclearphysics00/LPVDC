// detector_geometry.hh
#pragma once
#include "Garfield/ComponentAnalyticField.hh"
#include <vector>
#include <memory>
#include <string>

namespace Detector {

enum class ElectrodeKind { PlaneY, WireRow };

struct Electrode {
  ElectrodeKind kind = ElectrodeKind::PlaneY;
  double x0 = 0.0;      // [cm]
  double y  = 0.0;      // [cm]
  double V  = 0.0;      // [V]
  double radius = 1e-3; // [cm] (Wireのみ)
  double pitch  = 0.2;  // [cm] x方向の周期目安
  std::string name;
};

struct Geometry {
  double xmin = -0.5, xmax = +0.5;
  double ymin = -1.0, ymax = +1.0;
  bool periodicX = true;
  double pitchX  = 0.2;
  std::vector<Electrode> electrodes;
};

// Geometry -> Garfield::ComponentAnalyticField
inline std::unique_ptr<Garfield::ComponentAnalyticField>
BuildField(const Geometry& geo) {
  using namespace Garfield;
  auto comp = std::make_unique<ComponentAnalyticField>();
  if (geo.periodicX) comp->SetPeriodicityX(geo.pitchX);
  for (const auto& e : geo.electrodes) {
    if (e.kind == ElectrodeKind::PlaneY)  comp->AddPlaneY(e.y, e.V, e.name);
    else if (e.kind == ElectrodeKind::WireRow) comp->AddWire(e.x0, e.y, 2.0 * e.radius, e.V, e.name);
  }
  return comp;
}

} // namespace Detector

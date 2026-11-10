// SMVac — physical and mathematical constants.
// Author: Naman Goyal (IISc Bangalore). See docs/SCIENTIFIC_METHODOLOGY.md for conventions.
#pragma once

namespace smvac {

constexpr double kPi = 3.14159265358979323846;
constexpr double kPi2 = kPi * kPi;
constexpr double kLoop = 16.0 * kPi2;           // (4 pi)^2: loop factor in the beta functions
constexpr double kEightPi2Over3 = 8.0 * kPi2 / 3.0;

// Electroweak vev used in the tree-level relations of the matching (GeV).
constexpr double kVev = 246.22;
// Upper end of the RG running: non-reduced Planck mass (GeV), as in Buttazzo et al. (2013) Fig. 2.
constexpr double kPlanckMass = 1.22e19;
// Tree-level inputs for the small Yukawa couplings (GeV); their effect on lambda(M_Pl) is < 1e-5.
constexpr double kMbTree = 4.0;
constexpr double kMtauTree = 1.777;

// Unit conversions.
constexpr double kHbarGeVs = 6.582119569e-25;     // hbar in GeV s (CODATA 2018, exact SI definition)
constexpr double kMpcKm = 3.0856775814913673e19;  // 1 Mpc in km

// Past-light-cone space-time volume in units of H0^-4 (Buttazzo et al. 2013, eq. 82).
constexpr double kLightConeVolume = 0.15;
// Default Hubble rate used by Buttazzo et al. (2013), km/s/Mpc.
constexpr double kH0Default = 67.4;

}  // namespace smvac

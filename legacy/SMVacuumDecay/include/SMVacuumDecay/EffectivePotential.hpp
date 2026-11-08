#pragma once
#include <SMVacuumDecay/RGE.hpp>

namespace SMVacuumDecay {

double V1_loop(double phi, double yt);

double get_lambda_eff(const StandardModelParameters& p);

double V_eff(double phi, RGEHelper& rge);

}

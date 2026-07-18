#include <SMVacuumDecay/FubiniLipatov.hpp>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

int main() {
    // Standard Model Point
    double SM_Mh = 125.1;
    double SM_Mt = 173.1;
    double SM_ref_FL = 2297.667484139407;
    int SM_ref_Status = 2;
    
    // Low-Mass Point
    double Low_Mh = 5.0;
    double Low_Mt = 105.0;
    double Low_ref_FL = 11327.13339117053;
    int Low_ref_Status_FL = 2; // Analytical solver gives 2, Numerical gives 3
    
    auto sm_res = SMVacuumDecay::classify_buttazzo(SM_Mh, SM_Mt);
    int sm_status = std::get<0>(sm_res);
    double sm_action = std::get<1>(sm_res);
    double sm_mh = std::get<2>(sm_res);

    auto low_res = SMVacuumDecay::classify_buttazzo(Low_Mh, Low_Mt);
    int low_status = std::get<0>(low_res);
    double low_action = std::get<1>(low_res);
    double low_mh = std::get<2>(low_res);
    
    double err_sm = std::abs((sm_action - SM_ref_FL)/SM_ref_FL);
    double err_low = std::abs((low_action - Low_ref_FL)/Low_ref_FL);
    
    std::cout << std::scientific << std::setprecision(12);
    std::cout << "SM Error: " << err_sm << " Status: " << sm_status << " (Ref: " << SM_ref_Status << ")\n";
    std::cout << "Low Error: " << err_low << " Status: " << low_status << " (Ref: " << Low_ref_Status_FL << ")\n";
    
    double max_err = std::max(err_sm, err_low);
    std::cout << "FL Validation Max Error: " << max_err << std::endl;
    
    if (err_sm < 1e-12 && err_low < 1e-12 && sm_status == SM_ref_Status && low_status == Low_ref_Status_FL) {
        std::cout << "PASS" << std::endl;
        return 0;
    } else {
        std::cout << "FAIL" << std::endl;
        return 1;
    }
}

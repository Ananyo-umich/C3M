#include <cantera/base/Solution.h>
#include <cantera/thermo.h>
#include <iostream>

int main(int argc, char** argv)
{
    // Create a new Solution object
    auto sol = Cantera::newSolution("Vulcan_CHNOS.yaml");
    auto gas = sol->thermo();

    //std::cout << gas->temperature() << std::endl;
    return 0;
}

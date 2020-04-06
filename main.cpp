#include<iostream>
#include<fstream>
#include<cmath>
#include<chrono>
#include"simulation.h"

using namespace std;
using namespace Eigen;

int main()
{
    //record running time
    auto t_start = chrono::high_resolution_clock::now();
    //
    Simulation vlasov;
    vlasov.Run();
    auto t_end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(t_end - t_start).count();
    cout << "Running Time: " << duration << "s" << endl;
    cout << "**********************************" << endl;

    return 0;
}

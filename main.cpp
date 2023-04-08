/*--------------------------------------------------------------------------
Heisenberg模型においてM-Hカーブを調べるために模型がとる各磁化ごとのエネルギー
固有値を計算する.
計算したデータは./outputにcsv形式で吐き出す
コードの実行はmain.shを用いて行う。このbashファイル内で変数の受け渡しも行う
--------------------------------------------------------------------------*/

#include <string>
#include <vector>
#include <chrono>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "Subsystem_Sz.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    // min_up_spinからmax_up_spinまでの各部分空間での最小エネルギーを調べる際の雛形
    vector<string> file = {"./settings/jset0.txt", "./settings/jset1.txt", "./settings/jset2.txt"};
    int sys_num, sys_site_A, sys_site_B, min_up_spin, max_up_spin;
    std::string dir_output;

    cout << "argc = " << argc << endl;
    if (argc == 7)
    {
        sys_num = stoi(argv[1]);
        sys_site_A = stoi(argv[2]);
        sys_site_B = stoi(argv[3]);
        min_up_spin = stoi(argv[4]);
        max_up_spin = stoi(argv[5]);
        dir_output = argv[6];

        for (int i = 0; i < 6; i++)
        {
            cout << "argv[" << i << "] = " << argv[i] << endl;
        }
    }
    else
    {
        cout << "Usage: " << argv[0] << endl;
        for (int i = 0; i < argc; i++)
        {
            cout << "argv[" << i << "] = " << argv[i] << endl;
        }
    }

    double start_make_hamiltonian, end_make_hamiltonian, total_time_make_hamiltonian;
    double total_lanczos_time;

    start_make_hamiltonian = omp_get_wtime();
    Subsystem_Sz H(sys_num, sys_site_A, sys_site_B, file, min_up_spin);
    H.sub_space_check();
    H.set_system_info();
    H.sub_hamiltonian();
    end_make_hamiltonian = omp_get_wtime();

    total_lanczos_time = H.MP_sub_lanczos_timetest(1000, dir_output);
    cout << H << endl;

    ofstream ofs(dir_output, std::ios::app);
    ofs << setw(15) << left << "Construct Hamiltonian"
        << ","
        << "Total_lanczos_time" << endl;
    ofs.close();
}
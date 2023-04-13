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

    /*--------lanczos法の各処理に要する時間を調べる場合にはこちらを使用する-------*/
    // double start = omp_get_wtime();
    // Subsystem_Sz H(sys_num, sys_site_A, sys_site_B, file, min_up_spin); // Subsystem_Szオブジェクトのupスピン本数の再設定

    // // upスピンの本数が上記の場合における部分空間とHamiltonian行列の用意
    // H.sub_space_check();
    // H.set_system_info();
    // H.sub_hamiltonian();

    // H.MP_sub_lanczos_timetest(1000, dir_output); // lanczos法による固有値計算
    // double end = omp_get_wtime();

    // cout << "totla time[sec] : " << end - start << endl;
    // cout << H << endl;

    /*--------磁化曲線をplotするためのdataを得る場合にはこちらを使用する-------*/
    MP_plot_MHcurve(sys_num, sys_site_A, sys_site_B, max_up_spin, min_up_spin, file, dir_output);
}

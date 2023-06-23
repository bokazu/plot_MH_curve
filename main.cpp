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
    vector<string> file;
    int sys_num, sys_site_A, sys_site_B, min_up_spin, max_up_spin, start_up_spin, end_up_spin;
    double J_red, J_green, J_blue;
    std::string file_output_eval, dir_output_time, dir_output_spin_sxx_rel, dir_output_spin_szz_rel, dir_output_spin_sz_rel;
    char *lanczos_type;
    cout << "argc = " << argc << endl;
    if (argc == 20)
    {
        sys_num = stoi(argv[1]);
        sys_site_A = stoi(argv[2]);
        sys_site_B = stoi(argv[3]);
        min_up_spin = stoi(argv[4]);
        max_up_spin = stoi(argv[5]);
        J_red = atof(argv[6]);
        J_green = atof(argv[7]);
        J_blue = atof(argv[8]);
        file.push_back(argv[9]);
        file.push_back(argv[10]);
        file.push_back(argv[11]);
        file_output_eval = argv[12];
        dir_output_time = argv[13];
        dir_output_spin_sxx_rel = argv[14];
        dir_output_spin_sz_rel = argv[15];
        dir_output_spin_szz_rel = argv[16];
        start_up_spin = stoi(argv[17]);
        end_up_spin = stoi(argv[18]);
        lanczos_type = argv[19];


        cout << "input : scuccess\n";
        cout << "J_red                 : " << J_red << endl;
        cout << "J_green               : " << J_green << endl;
        cout << "J_blue                : " << J_blue << endl;
        cout << "- File of eigen value : " << file_output_eval << endl;
        cout << "- File of time        : " << dir_output_time << endl;
        cout << "- Dir of <SxSx>       : " << dir_output_spin_sxx_rel << endl;
        cout << "- Dir of <Sz>         : " << dir_output_spin_sz_rel << endl;
        cout << "- Dir of <SzSz>       : " << dir_output_spin_szz_rel << endl;
        cout << "start up spin         : " << start_up_spin << endl;
        cout << "end up spin         : " << end_up_spin << endl;
    }
    else
    {
        cout << "Usage: " << argv[0] << endl;
        for (int i = 0; i < argc; i++)
        {
            cout << "argv[" << i << "] = " << argv[i] << endl;
        }
    }

    ranged_calc_gs_energy(sys_num, sys_site_A, sys_site_B, max_up_spin, start_up_spin, end_up_spin, J_red, J_green, J_blue, file, file_output_eval, dir_output_time, dir_output_spin_sxx_rel, dir_output_spin_sz_rel, dir_output_spin_szz_rel, *lanczos_type);
}

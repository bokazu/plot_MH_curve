#include <iostream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <string>

using namespace std;

void calc_MHcurve_data(double J_red, double J_green, double J_blue, std::string input_file, std::string output_file)
{
    vector<double> M;     // 磁化の情報を格納するための配列
    vector<double> e_min; // 各部分空間における基底状態のエネルギー固有値を代入するための配列

    vector<double> plot_h; // 交点のh座標を格納する(横軸)
    vector<double> plot_m; // 交点のh座標を格納する(縦軸)

    vector<double> plateau_width; // h_{i+1} - h_iの値を格納する

    // ファイルから磁化と対応するエネルギー固有値の情報を読み取る
    // inputファイルの磁化の値は飽和磁化で規格化済みであるとする
    ifstream input(input_file);
    string str_tmp;
    double mag_tmp, energy_tmp;

    while (getline(input, str_tmp))
    {
        stringstream ss;
        ss << str_tmp;
        ss >> mag_tmp >> energy_tmp;
        M.push_back(mag_tmp);
        e_min.push_back(energy_tmp);
    }

    input.close();

    // M0についてはここで格納する
    // plot_h.push_back(0);
    // plot_m.push_back(M[0]);

    for (int i = 0; i < M.size() - 1; i++)
    {
        vector<double> mi_intersec(M.size(), 50.0); // M_iとM_i+1,M_i+2...との交点のx座標を格納した配列

        for (int k = i + 1; k < M.size(); k++)
        {
            mi_intersec[k - (i + 1)] = double((e_min[k] - e_min[i]) / (k - i));
        }

        // 列挙した交点のうち、最もx座標の値が小さいものを選択して記録する
        auto min_index = std::distance(mi_intersec.begin(), std::min_element(mi_intersec.begin(), mi_intersec.end())); // min_indexの方はsize_t
        // M[i]とM[i+1]のhの交点座標 > M[i]とM[i+2]のhの交点座標の場合にM[i+1]とM[i+2],M[i+3]...のhの交点を調べるのはskipする
        int skip_itr = int(min_index) + 1;

        // skipされるぶんについても値を出力する
        for (int l = 0; l < skip_itr - 1; l++)
        {
            plot_h.push_back(mi_intersec[min_index]);
            plot_m.push_back(M[i + l + 1]);
        }

        plot_h.push_back(mi_intersec[min_index]); // 変更箇所2
        plot_m.push_back(M[i + skip_itr]);
        // 上記に該当する磁化のindexを取得する
        i += int(min_index); // M_iとM_{i+α}の交点が最小だった場合には、M_{i+α-1}に着目したloopはskipする。
    }

    //====================== プラトー幅を調べる ===============================
    double width;
    for (int i = 0; i < M.size() - 1; i++)
    {
        width = plot_h[i + 1] - plot_h[i];
        plateau_width.push_back(width);
    }

    plateau_width.push_back(0.);

    // 交点の情報の表示
    cout
        << "================================================\n";
    cout << internal << setw(15) << "h"
         << "  " << setw(15) << "M"
         << " "
         << "J_{red}"
         << " "
         << "J_{green}"
         << " "
         << "J_{blue}"
         << " "
         << "plateau width" << endl;
    cout << "------------------------------------------------\n";
    for (int i = 0; i < plot_h.size(); i++)
    {
        cout << setprecision(15) << left << setw(15) << plot_h[i]
             << " " << setw(15) << plot_m[i]
             << " " << J_red
             << " " << J_green
             << " " << J_blue
             << " " << plateau_width[i] << endl;
    }
    cout << "================================================\n";

    // 交点の情報をファイルへ書き出す
    ofstream plateau_data(output_file);
    for (int i = 0; i < plot_h.size(); i++)
    {
        plateau_data << setprecision(15) << left << plot_h[i]
                     << " " << setw(15) << plot_m[i]
                     << " " << J_red
                     << " " << J_green
                     << " " << J_blue
                     << " " << plateau_width[i] << endl;
    }
    plateau_data.close();
};

int main(int argc, char *argv[])
{
    int start_up_spin, end_up_spin;
    double J_red, J_green, J_blue;
    std::string input_file, output_file;

    if (argc == 8)
    {
        start_up_spin = atoi(argv[1]);
        end_up_spin = atoi(argv[2]);
        J_red = atof(argv[3]);
        J_green = atof(argv[4]);
        J_blue = atof(argv[5]);
        input_file = argv[6];
        output_file = argv[7];

        calc_MHcurve_data(J_red, J_green, J_blue, input_file, output_file);
    }
    else
    {
        return EXIT_FAILURE;
    }
}
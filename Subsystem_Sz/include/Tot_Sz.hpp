/*
特定の(totS_A^z, totS_B^z)での部分状態空間の情報やHamiltonian行列、Lanczosベクトルをメンバーとする
*/

#ifndef ___Class_Tot_Sz
#define ___Class_Tot_Sz

#include <iomanip>
#include <string>
#include <vector>
#include <boost/dynamic_bitset.hpp>

#include "MEIGEN.hpp"
#include "Iso_Hamiltonian.hpp"
#include "Int_Hamiltonian.hpp"
#include "Jset.hpp"

class Tot_Sz
{
public:
    static int system_num;
    static int tot_site_A;
    static int tot_site_B;
    int up_spin_pair[2]; // up_spin_pair[0] : 部分系Aのupスピンの本数
    int *bm_A;           // 部分系Aの状態空間を張るスピン状態番号を記録する
    int bm_A_size;       // 配列bm_Aのサイズ
    int *gbm_A;          // スピン状態番号がブロック行列の何列目かを指定する
    int gbm_A_size;
    int *bm_B;
    int bm_B_size;
    int *gbm_B;
    int gbm_B_size;

    /*-----Lanczosベクトル-----*/
    double **V0;
    double **V1;

    /*Hamiltonian行列用のオブジェクト*/
    Iso_Hamiltonian H_iso[2];
    Int_Hamiltonian *H_int;

    MEIGEN Eig;

    // コンストラクタ
    Tot_Sz();

    // デストラクタ
    ~Tot_Sz()
    {
        delete[] bm_A;
        delete[] gbm_A;
        delete[] bm_B;
        delete[] gbm_B;
        //[ToDo] **V0と**V1のメモリの開放はsub_lanczos内で行っている
        delete[] H_int;

        delete[] V0;
        delete[] V1;
        // std::cout << "Tot_Sz::destructed.\n";
    }

    // Hamiltonian行列の非ゼロ要素数のカウントを行う
    void sub_count_isonnz();
    void sub_count_intnnz(const int No, const int pair_num, const int system_num, const int siteA_num, const int siteB_num);
    // siteA_numはAの総サイト数

    void sub_iso_hamiltonian(const int siteA_num, const int siteB_num); //[Note] documentにはsub_Iso_Hamilotnian()と書いている
};

#endif

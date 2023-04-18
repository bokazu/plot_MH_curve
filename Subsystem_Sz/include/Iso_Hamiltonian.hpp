#ifndef ___Class_Iso_Hamiltonian
#define ___Class_Iso_Hamiltonian

#include <fstream>
#include <iomanip>
#include <iostream>

#include "Jset.hpp"
class Iso_Hamiltonian
{
public:
    int Iso_id;
    static int Iso_id_counter;

    int tot_site_num;
    int mat_dim;
    int nnz;
    int *row;
    int *col_ptr;
    double *val;

    Jset J;
    // コンストラクタ
    Iso_Hamiltonian();
    // デストラクタ
    ~Iso_Hamiltonian()
    {
        delete[] row;
        delete[] col_ptr;
        delete[] val;

        // std::cout << "Iso_Hamiltonian::destructed.\n";
    }

    // 部分空間でのS^+|m> S^-|m>、S^z|m>に相当する演算を行う
    void iso_sub_spin(const int bm_itr, int *bm, int *gbm, int gbm_size, const int tot_site_num, const int site_i, int &row_index, int &col_ptr_val, double &szz);

    /*---------------------------mm系関数-----------------------------*/
    //[Note]あとでこれはSubsystem_Szのメンバ関数にするかもしれない
    // V+= HVを計算する
    void isoA_mmprod(int dim_A, int dim_B, double **u_i, double **u_j);

    // += ΨHを計算する
    void isoB_mmprod(int dim_A, int dim_B, double **u_i, double **u_j);
    /*-----------------------------------------------------------------*/

    /*-------------------上記mm系関数のOpenMP利用version-----------------------*/
    void MP_isoA_mmprod(int dim_A, int dim_B, double **u_i, double **u_j);
    void MP_isoB_mmprod(int dim_A, int dim_B, double **u_i, double **u_j);

    /*-------------------上記mm系関数のOpenMP with scheduling利用version-----------------------*/
    void MP_schedule_isoA_mmprod(int dim_A, int dim_B, double **u_i, double **u_j);
    void MP_schedule_isoB_mmprod(int dim_A, int dim_B, double **u_i, double **u_j);

    // idの値を返す
    int iso_id() const { return Iso_id; }

    // row, col, valの初期化をおこなう
    void init();

    // jset_filenameをセットする
    void set_jset(std::string filename) { J.set_filename(filename); }

    // Iso_Hamiltonianオブジェクトの文字列表現を返却する
    std::string to_string() const;
};

// 出力ストリームにhを挿入する
std::ostream &operator<<(std::ostream &s, const Iso_Hamiltonian &h);

#endif

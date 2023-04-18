#ifndef ___Class_Int_Hamiltonian
#define ___Class_Int_Hamiltonian

#include <fstream>
#include <iomanip>
#include <iostream>

#include "Jset.hpp"

class Int_Hamiltonian
{
public:
    int Int_id;
    static int Int_id_counter;

    /*行列S^+の非ゼロ要素の行、列番号*/
    int *prow_ind_A;
    int *prow_ind_B;
    int *pcol_ind_A;
    int *pcol_ind_B;
    int *mrow_ind_A;
    int *mrow_ind_B;
    int *mcol_ind_A;
    int *mcol_ind_B;
    int *sz_A;
    int *sz_B;

    static int tot_site_A;
    static int tot_site_B;
    int dim_A;
    int dim_B;
    int nnz_pA;
    int nnz_pB;
    int nnz_mA;
    int nnz_mB;
    Jset J;
    // コンストラクタ
    Int_Hamiltonian();
    // デストラクタ
    ~Int_Hamiltonian()
    {
        delete[] prow_ind_A;
        delete[] prow_ind_B;
        delete[] pcol_ind_A;
        delete[] pcol_ind_B;
        delete[] mrow_ind_A;
        delete[] mrow_ind_B;
        delete[] mcol_ind_A;
        delete[] mcol_ind_B;
        delete[] sz_A;
        delete[] sz_B;

        // std::cout << "Int_Hamiltonian::destructed.\n";
    };

    void sub_count_intnnz();

    int int_id() const { return Int_id; }

    // 配列の初期化
    void init();

    // jset_filenameをセットする
    void set_jset(std::string filename) { J.set_filename(filename); }

    void int_dsmn_mmprod(const int nnz_mA, const int nnz_mB, int *mrow_ind_A, int *mcol_ind_A, int *mrow_ind_B, int *mcol_ind_B, double **V0, double **V1, double **V1_inc1); // S_A^-VS_B^-
    void int_rise_mmprod(const int nnz_pA, const int nnz_pB, int *prow_ind_A, int *pcol_ind_A, int *prow_ind_B, int *pcol_ind_B, double **V0, double **V1, double **V1_dic1); // S_A^+VS_B^+
    void int_zz_mmprod(const int dim_A, const int dim_B, int *sz_A, int *sz_B, double **V0, double **V1);                                                                     // S_A^zVS_B^z

    /*------上記関数のopneMP利用version-------*/
    void MP_int_dsmn_mmprod(const int nnz_mA, const int nnz_mB, int *mrow_ind_A, int *mcol_ind_A, int *mrow_ind_B, int *mcol_ind_B, double **V0, double **V1, double **V1_inc1); // S_A^-VS_B^-
    void MP_int_rise_mmprod(const int nnz_pA, const int nnz_pB, int *prow_ind_A, int *pcol_ind_A, int *prow_ind_B, int *pcol_ind_B, double **V0, double **V1, double **V1_dic1); // S_A^+VS_B^+
    void MP_int_zz_mmprod(const int dim_A, const int dim_B, int *sz_A, int *sz_B, double **V0, double **V1);                                                                     // S_A^zVS_B^z

    /*------上記関数のopneMP with schedule利用version-------*/
    void MP_schedule_int_dsmn_mmprod(const int nnz_mA, const int nnz_mB, int *mrow_ind_A, int *mcol_ind_A, int *mrow_ind_B, int *mcol_ind_B, double **V0, double **V1, double **V1_inc1); // S_A^-VS_B^-
    void MP_schedule_int_rise_mmprod(const int nnz_pA, const int nnz_pB, int *prow_ind_A, int *pcol_ind_A, int *prow_ind_B, int *pcol_ind_B, double **V0, double **V1, double **V1_dic1); // S_A^+VS_B^+
    void MP_schedule_int_zz_mmprod(const int dim_A, const int dim_B, int *sz_A, int *sz_B, double **V0, double **V1);                                                                     // S_A^zVS_B^z

    // Int_Hamiltonianオブジェクトの文字列表現を返却する
    std::string to_string() const;
};

// 出力ストリームにhを挿入する
std::ostream &
operator<<(std::ostream &s, const Int_Hamiltonian &h);

#endif

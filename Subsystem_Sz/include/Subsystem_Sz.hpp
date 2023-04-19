#ifndef ___Class_Subsystem_Sz
#define ___Class_Subsystem_Sz

#include <mkl.h>

#include <iomanip>
#include <string>
#include <vector>
#include <math.h>
#include <random>
#include <omp.h>

#include "MEIGEN.hpp"
#include "Iso_Hamiltonian.hpp"
#include "Int_Hamiltonian.hpp"
#include "Tot_Sz.hpp"
#include "Jset.hpp"

class Subsystem_Sz
{
public:
    /*系についての情報*/
    int system_num; // 部分系の個数　
    int tot_site_A; // 部分系Aのサイト数　
    int tot_site_B; // 部分系Bのサイト数
    std::vector<std::string> filename;

    /*部分空間の情報*/
    int up_spin; // 複合系のupスピンのサイト数
    double mag;
    int block_row_num; // 磁化の値がmagの部分空間の張る状態行列の行数
    int block_col_num; // 磁化の値がmagの部分空間の張る状態行列の列数
    int pair_num;      // up_spin = up_spin_A + up_spin_Bを満たす(up_spin_A, up_spin_B)の組み合わせの総数

    Tot_Sz *tot_Sz;

    /*lanczos法関連*/
    int ls_count;  // 収束までに要したstep数を表示する
    bool ls_check; // 指定した反復回数内で計算が収束したかを表すための変数
    double eigen_value;
    // MEIGEN Eig; //[Note]これもtot_Szオブジェクトごとに用意しておいたほうがいいかもしれない

    // コンストラクタ
    Subsystem_Sz(int sys_num, int site_A, int site_B, std::vector<std::string> &file, int up);

    // デストラクタ
    ~Subsystem_Sz()
    {
        delete[] tot_Sz;
        std::cout << "Subsystem is destructed\n";
    };

    /*-------------メンバ関数---------------*/
    //(totS_A^zm, tot_S_B^z)のペア数のカウントを行うための関数
    int count_No();

    void set_system_info();

    // 生成されたTot_Szオブジェクトの個数を出力する
    void print_pairnum() { std::cout << "There is " << pair_num << " Tot_Sz Objects.\n"; }
    void print_up_pair();
    void print_subspace();
    void pirnt_id();
    void print_member(); // 全メンバの値を出力する

    // 各Tot_Szオブジェクトのupスピンの本数を計算し記録する
    void check_up_pair();

    // 部分空間を張るスピン状態番号を確認する
    void sub_space_check();
    // 上記の関数で使用
    void sub_space_check_A(const int No);
    void sub_space_check_B(const int No);

    // Hamiltonian行列の非ゼロ要素数のカウントをおこなう
    void sub_count_nnz(const int No);

    // 配列のメモリの再確保を行う
    void sub_resize(const int No);

    /*Subsystem_Sz::sub_resize()中で以下の関数を使用する*/
    void calloc_splus(const int No, const int id); // S^+の配列を確保して0初期化
    void calloc_smins(const int No, const int id); // S^-の配列を確保して0初期化
    void calloc_szz(const int No, const int id);   // S^zの配列を確保して0初期化

    /*----------------Int systemのHamiltonian行列要素の計算を行う------------------*/
    void sub_int_hamiltonian(); //[Note] documentにはsub_Int_Hamilotnian()と書いている
    // Hamiltonian行列の行列要素の計算と配列への格納を行う
    void calc_mat_splus(const int No);
    void calc_mat_smins(const int No);
    void calc_mat_szz(const int No);
    /*------------------------------------------------------------------------*/

    void sub_hamiltonian();

    // lanczos法
    void sub_lanczos(const int tri_mat_dim, char c = 'N', char info_ls = 'n');
    double sub_lanczos_timetest(const int tri_mat_dim, std::string output_filename, char c = 'N', char info_ls = 'n');
    void MP_sub_lanczos(const int tri_mat_dim, char c = 'N', char info_ls = 'n');
    double MP_sub_lanczos_timetest(const int tri_mat_dim, std::string output_filename, char c = 'N', char info_ls = 'n'); // lanczos法で収束に要した時間を返す
    void MP_schedule_sub_lanczos(const int tri_mat_dim, char c = 'N', char info_ls = 'n');
    double MP_schedule_sub_lanczos_timetest(const int tri_mat_dim, std::string output_filename, char c = 'N', char info_ls = 'n'); // lanczos法で収束に要した時間を変えす

    /*------------------------mm系関数------------------------*/
    void iso_mmprod(const int No, double **V0, double **V1);                                     // Iso systemの行列-行列積
    void int_mmprod(const int No, double **V0, double **V1, double **V1_dic1, double **V1_inc1); // Int systemの行列-行列-行列積
    void int_mmzzord(const int No, double **V0, double **V1);
    void calc_alpha_evenstep(const int ls, double *alpha);
    void calc_alpha_oddstep(const int ls, double *alpha);
    void calc_beta_evenstep(const int tri_mat_dim, const int ls, double *alpha, double *beta);
    void calc_beta_oddstep(const int tri_mat_dim, const int ls, double *alpha, double *beta);

    double mm_ddot(const int row_dim, const int col_dim, double **V0, double **V1);
    void mm_dscal(double alpha, const int row_dim, const int col_dim, double **V);
    void mm_dcopy(const int row_dim, const int col_dim, double **V0, double **V1);
    void mm_daxpy(double alpha, const int row_dim, const int col_dim, double **V0, double **V1);
    double mm_dnrm2(const int row_dim, const int col_dim, double **V);
    void mm_sdz(const int row_dim, const int col_dim, double **V);
    void mm_sdz_V0(); // V[0]~V[pari_num-1]までトータルで考えて規格化する。 NG)V[0],V[1]...ごとに規格化する
    void mm_sdz_V1();
    void mm_init(const int row_dim, const int col_dim, double **V);
    /*-------------------------------------------------------*/

    // MP_sub_lanczos内のコードでは以下の関数を利用する
    /*-------------------mm関数 OpenMP利用version------------------*/
    void MP_iso_mmprod(const int No, double **V0, double **V1);                                     // Iso systemの行列-行列積
    void MP_int_mmprod(const int No, double **V0, double **V1, double **V1_dic1, double **V1_inc1); // Int systemの行列-行列-行列積
    void MP_int_mmzzord(const int No, double **V0, double **V1);

    void MP_calc_alpha_evenstep(const int ls, double *alpha);
    void MP_calc_alpha_oddstep(const int ls, double *alpha);
    void MP_calc_beta_evenstep(const int tri_mat_dim, const int ls, double *alpha, double *beta);
    void MP_calc_beta_oddstep(const int tri_mat_dim, const int ls, double *alpha, double *beta);

    double MP_mm_ddot(const int row_dim, const int col_dim, double **V0, double **V1);
    void MP_mm_dscal(double alpha, const int row_dim, const int col_dim, double **V);
    void MP_mm_dcopy(const int row_dim, const int col_dim, double **V0, double **V1);
    void MP_mm_daxpy(double alpha, const int row_dim, const int col_dim, double **V0, double **V1);
    double MP_mm_dnrm2(const int row_dim, const int col_dim, double **V);
    void MP_mm_sdz(const int row_dim, const int col_dim, double **V);
    void MP_mm_sdz_V0(); // V[0]~V[pari_num-1]までトータルで考えて規格化する。 NG)V[0],V[1]...ごとに規格化する
    void MP_mm_sdz_V1();
    void MP_mm_init(const int row_dim, const int col_dim, double **V);

    /*-------------------mm関数 OpenMP with scheduling利用version------------------*/
    void MP_schedule_iso_mmprod(const int No, double **V0, double **V1);                                     // Iso systemの行列-行列積
    void MP_schedule_int_mmprod(const int No, double **V0, double **V1, double **V1_dic1, double **V1_inc1); // Int systemの行列-行列-行列積
    void MP_schedule_int_mmzzord(const int No, double **V0, double **V1);

    void MP_schedule_calc_alpha_evenstep(const int ls, double *alpha);
    void MP_schedule_calc_alpha_oddstep(const int ls, double *alpha);
    void MP_schedule_calc_beta_evenstep(const int tri_mat_dim, const int ls, double *alpha, double *beta);
    void MP_schedule_calc_beta_oddstep(const int tri_mat_dim, const int ls, double *alpha, double *beta);

    double MP_schedule_mm_ddot(const int row_dim, const int col_dim, double **V0, double **V1);
    void MP_schedule_mm_dscal(double alpha, const int row_dim, const int col_dim, double **V);
    void MP_schedule_mm_dcopy(const int row_dim, const int col_dim, double **V0, double **V1);
    void MP_schedule_mm_daxpy(double alpha, const int row_dim, const int col_dim, double **V0, double **V1);
    double MP_schedule_mm_dnrm2(const int row_dim, const int col_dim, double **V);
    void MP_schedule_mm_sdz(const int row_dim, const int col_dim, double **V);
    void MP_schedule_mm_sdz_V0(); // V[0]~V[pari_num-1]までトータルで考えて規格化する。 NG)V[0],V[1]...ごとに規格化する
    void MP_schedule_mm_sdz_V1();
    void MP_schedule_mm_init(const int row_dim, const int col_dim, double **V);

    // spin-spin相関の計算<Ψ|S_i^zS_j^z|Ψ>
    void clac_spin_rel(const int site_num, std::string dir_output);

    // 特定のplateauについて相互作用のパラメータを変化させたときのwidthの変化を調べる
    // M-H curveについてのデータファイルが揃っていることが前提
    void measure_plateau_width(std::string dir_input, std::string dir_output);
    //  Subsystem_Szオブジェクトの文字列表現を返却する
    std::string to_string() const;
};

std::ostream &operator<<(std::ostream &s, const Subsystem_Sz &H);

template <typename T>
void vec_init(int dim, T *vec)
{
    for (int i = 0; i < dim; i++)
    {
        vec[i] = 0;
    }
}

template <typename T>
void MP_vec_init(int dim, T *vec)
{
#pragma omp parallel for
    for (int i = 0; i < dim; i++)
    {
        vec[i] = 0;
    }
}

template <typename T>
void MP_schedule_vec_init(int dim, T *vec)
{
#pragma omp parallel for schedule(runtime)
    for (int i = 0; i < dim; i++)
    {
        vec[i] = 0;
    }
}

void plot_MHcurve(int sys_num, int sys_site_A, int sys_site_B, int max_up_spin, int min_up_spin, double J_red, double J_green, double J_blue, std::vector<std::string> &file, std::string GNUPLOT_DATA_DIR);
void MP_plot_MHcurve(int sys_num, int sys_site_A, int sys_site_B, int max_up_spin, int min_up_spin, double J_red, double J_green, double J_blue, std::vector<std::string> &file, std::string GNUPLOT_DATA_DIR);
void MP_schedule_plot_MHcurve(int sys_num, int sys_site_A, int sys_site_B, int max_up_spin, int min_up_spin, double J_red, double J_green, double J_blue, std::vector<std::string> &file, std::string GNUPLOT_DATA_DIR);

// コンビネーションnCrの計算を行う
int comb(int n, int r);

#endif

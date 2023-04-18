#include "../include/Iso_Hamiltonian.hpp"

#include <mkl.h>

#include <boost/dynamic_bitset.hpp>
#include <fstream>
#include <iomanip>
#include <random>

using namespace std;

int Iso_Hamiltonian::Iso_id_counter = 0;

// コンストラクタ
Iso_Hamiltonian::Iso_Hamiltonian()
    : Iso_id(Iso_id_counter++),
      tot_site_num(0),
      mat_dim(1), // tot_site_num = 0で初期化してるのでバグる
      nnz(0),     // → tot_Szオブジェクトのコンストラクタで初期化
      row(new int[1]),
      col_ptr(new int[1]),
      val(new double[1]),
      J(" ")
{
    // std::cout << "Iso_Hamiltonian::constructed.\n";
}

void Iso_Hamiltonian::iso_sub_spin(const int bm_itr, int *bm, int *gbm, int gbm_size, const int tot_site_num, const int site_i, int &row_index, int &col_ptr_val, double &szz)
{
    int m = bm[bm_itr];
    int bm_ctr;
    boost::dynamic_bitset<> ket_m(tot_site_num, m);
    bool is_up_site_i, is_up_site_j;

    int j_line = 0;
    int Jset_line = J.get_line();

    for (int l = j_line; l < Jset_line; l++)
    {
        if (J.index(0, l) == site_i)
        {
            int site_j = J.index(1, l);
            is_up_site_i = ket_m.test(site_i);
            is_up_site_j = ket_m.test(site_j);

            // 注目する2siteのスピン状態を場合分けする　
            if (is_up_site_i == is_up_site_j)
            {
                szz += 0.25 * J.val(l);
            }
            else
            {
                boost::dynamic_bitset<> ket_m1(tot_site_num, m);

                ket_m1.flip(site_i);
                ket_m1.flip(site_j);

                int n = (int)(ket_m1.to_ulong());

                if (n < gbm_size)
                {
                    bm_ctr = gbm[n]; // spin状態が対応する部分空間の行列において何列目かを取得
                    row[row_index] = bm_ctr;
                    val[row_index] = 0.5 * J.val(l);
                    row_index++;
                    col_ptr_val++;
                } // nが対応する部分空間に存在するかを確認する

                szz -= 0.25 * J.val(l);
            }
            j_line++;
        }
    }
};

/// Iso_HamiltonianはCSC形式で保存しているが、エルミート行列であることを利用してCSR形式での行列-行列積の計算方法を使っている
void Iso_Hamiltonian::isoA_mmprod(int dim_A, int dim_B, double **u_i, double **u_j)
{
    //[memo] memory連続アクセスに配慮したloopにした
    for (int row_num = 0; row_num < dim_A; row_num++)
    {
        for (int k = col_ptr[row_num]; k < col_ptr[row_num + 1]; k++)
        {
            for (int col_num = 0; col_num < dim_B; col_num++)
            {
                u_j[row_num][col_num] += val[k] * u_i[row[k]][col_num];
            }
        }
    }
}

void Iso_Hamiltonian::isoB_mmprod(int dim_A, int dim_B, double **u_i, double **u_j)
{
    //[memo] memory連続アクセスに配慮したloopにした
    for (int row_num = 0; row_num < dim_A; row_num++)
    {
        for (int col_num = 0; col_num < dim_B; col_num++)
        {
            for (int k = col_ptr[col_num]; k < col_ptr[col_num + 1]; k++)
            {
                u_j[row_num][col_num] += u_i[row_num][row[k]] * val[k];
            }
        }
    }
}

/*----------上記mm関数のOpenMP利用version----------*/
void Iso_Hamiltonian::MP_isoA_mmprod(int dim_A, int dim_B, double **u_i, double **u_j)
{
    int k;
    int col_num;
#pragma omp parallel for private(k, col_num)
    for (int row_num = 0; row_num < dim_A; row_num++)
    {
        for (k = col_ptr[row_num]; k < col_ptr[row_num + 1]; k++)
        {
            for (col_num = 0; col_num < dim_B; col_num++)
            {
                u_j[row_num][col_num] += val[k] * u_i[row[k]][col_num];
            }
        }
    }
}

void Iso_Hamiltonian::MP_schedule_isoA_mmprod(int dim_A, int dim_B, double **u_i, double **u_j)
{
    int k;
    int col_num;
#pragma omp parallel for private(k, col_num) schedule(runtime)
    for (int row_num = 0; row_num < dim_A; row_num++)
    {
        for (k = col_ptr[row_num]; k < col_ptr[row_num + 1]; k++)
        {
            for (col_num = 0; col_num < dim_B; col_num++)
            {
                u_j[row_num][col_num] += val[k] * u_i[row[k]][col_num];
            }
        }
    }
}

void Iso_Hamiltonian::MP_isoB_mmprod(int dim_A, int dim_B, double **u_i, double **u_j)
{
    int k, col_num;

#pragma omp parallel for private(col_num, k)
    for (int row_num = 0; row_num < dim_A; row_num++)
    {
        for (col_num = 0; col_num < dim_B; col_num++)
        {
            for (k = col_ptr[col_num]; k < col_ptr[col_num + 1]; k++)
            {
                u_j[row_num][col_num] += u_i[row_num][row[k]] * val[k];
            }
        }
    }
}

void Iso_Hamiltonian::MP_schedule_isoB_mmprod(int dim_A, int dim_B, double **u_i, double **u_j)
{
    int k, col_num;

#pragma omp parallel for private(col_num, k) schedule(runtime)
    for (int row_num = 0; row_num < dim_A; row_num++)
    {
        for (col_num = 0; col_num < dim_B; col_num++)
        {
            for (k = col_ptr[col_num]; k < col_ptr[col_num + 1]; k++)
            {
                u_j[row_num][col_num] += u_i[row_num][row[k]] * val[k];
            }
        }
    }
}
// row, col_ptr,valの初期化を行う
void Iso_Hamiltonian::init()
{
    for (int i = 0; i < nnz; i++)
    {
        row[i] = 0;
        val[i] = 0.0;
    }

    for (int i = 0; i < mat_dim + 1; i++)
    {
        col_ptr[i] = 0;
    }
}

// 文字列表現を返却する
std::string Iso_Hamiltonian::to_string() const
{
    std::ostringstream s;

    s << "-------------------------------------------------------\n";
    s << "@Number of site   : " << tot_site_num << std::endl;
    s << "@Matrix dimension : " << mat_dim << std::endl;
    s << "@non zero element : " << nnz << std::endl;
    s << J.to_string() << std::endl;
    s << "-------------------------------------------------------\n";

    return s.str();
}

ostream &operator<<(ostream &s, const Iso_Hamiltonian &h)
{
    return s << h.to_string();
}

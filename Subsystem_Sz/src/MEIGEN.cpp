#include "../include/MEIGEN.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

// コピーコンストラクタ
MEIGEN::MEIGEN(const MEIGEN &x) : dim_A(x.dim_A), dim_B(x.dim_B), eigen_val(x.eigen_val)
{
    eigen_mat = new double *[dim_A];
    for (int i = 0; i < dim_A; i++)
    {
        eigen_mat[i] = new double[dim_B];
    }

    for (int i = 0; i < dim_A; i++)
    {
        for (int j = 0; j < dim_B; j++)
        {
            eigen_mat[i][j] = x.eigen_mat[i][j];
        }
    }
    std::cout << "copy has done\n";
}

// 代入演算子
MEIGEN &
MEIGEN::operator=(const MEIGEN &x)
{
    for (int i = 0; i < dim_A; i++)
    {
        for (int j = 0; j < dim_B; j++)
        {
            eigen_mat[i][j] = x.eigen_mat[i][j];
        }
    }
    return *this;
}

// 固有ベクトルの初期化(0で初期化する)
void MEIGEN::vec_init()
{
    for (int i = 0; i < dim_A; i++)
    {
        for (int j = 0; j < dim_B; j++)
        {
            eigen_mat[i][j] = 0.0;
        }
    }
}

void MEIGEN::MP_scheduled_evec_init()
{
    int i, j;
#pragma omp parallel for private(j) schedule(runtime)
    for (i = 0; i < dim_A; i++)
    {
        for (j = 0; j < dim_B; j++)
        {
            eigen_mat[i][j] = 0.0;
        }
    }
}

// 固有ベクトルの要素数を変更し、初期化する
void MEIGEN::evec_elem(int m, int n)
{
    for (int i = 0; i < dim_A; i++)
    {
        delete[] eigen_mat[i];
    }
    delete[] eigen_mat;

    dim_A = m;
    dim_B = n;
    eigen_mat = new double *[dim_A];
    for (int i = 0; i < dim_A; i++)
    {
        eigen_mat[i] = new double[dim_B];
    }

    for (int i = 0; i < dim_A; i++)
    {
        for (int j = 0; j < dim_B; j++)
        {
            eigen_mat[i][j] = 0.0;
        }
    }
}

// 文字列表現を返却する
std::string MEIGEN::to_string() const
{
    std::ostringstream s;

    s << "@eigen_value = " << eigen_val << std::endl;
    s << "@eigen_vector  \n";
    for (int i = 0; i < dim_A; i++)
    {
        for (int j = 0; j < dim_B; j++)
        {
            s << "eigen_vec[" << i * dim_B + j << "] = " << eigen_mat[i][j] << std::endl;
        }
    }
    return s.str();
}

void MEIGEN::to_file(std::string filename) const
{
    std::ofstream ofs(filename);

    ofs << "@eigen_value = " << eigen_val << std::endl;
    ofs << "@eigen_vector  \n";
    for (int i = 0; i < dim_A; i++)
    {
        for (int j = 0; j < dim_B; j++)
        {
            ofs << "eigen_vec[" << i * dim_B + j << "] = " << eigen_mat[i][j] << std::endl;
        }
    }
}

ostream &operator<<(ostream &s, const MEIGEN &x) { return s << x.to_string(); }

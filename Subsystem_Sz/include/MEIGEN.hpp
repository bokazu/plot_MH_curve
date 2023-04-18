#ifndef ___Class_MEIGEN
#define ___Class_MEIGEN

#include <mkl.h>

#include <iomanip>
#include <iostream>

class MEIGEN
{
public:
    int dim_A; // 行数
    int dim_B; // 列数

public:
    double eigen_val;
    double **eigen_mat; // M x N行列(M≧N)
    // コンストラクタ
    MEIGEN(int dimA, int dimB) : dim_A(dimA), dim_B(dimB), eigen_val(0)
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
                eigen_mat[i][j] = 0.0;
            }
        }
        // std::cout << "EIGEN::constructed" << std::endl;
    }
    // コピーコンストラクタ
    MEIGEN(const MEIGEN &eigen);

    // デストラクタ
    ~MEIGEN()
    {
        for (int i = 0; i < dim_A; i++)
        {
            delete[] eigen_mat[i];
        }
        delete[] eigen_mat;
        // std::cout << "EIGEN::destructed" << std::endl;
    }

    // 代入演算子=
    MEIGEN &operator=(const MEIGEN &x);

    // 添字演算子[]
    //  double& operator[][](int i,int j) { return eigen_mat[i][j]; }

    // // const版添字演算子[]
    // const double& operator[][](int i,int j) const { return eigen_mat[i][j]; }

    // 固有値をセットする
    void set_eval(double val) { eigen_val = val; }

    // 固有値の値を取得する
    double eval() const { return eigen_val; }

    // 固有ベクトルの第[i]成分を返却する
    double emat(int i, int j) const { return eigen_mat[i][j]; }

    double **data() { return eigen_mat; }
    // 固有ベクトルの値を初期化する
    void vec_init();

    // 固有ベクトルの要素数を変更し、初期化する
    void evec_elem(int m, int n);

    // 固有値、固有ベクトルの文字列表現を返却する
    std::string to_string() const;
    void to_file(std::string filename) const;
};

// 出力ストリームにxを挿入する
std::ostream &operator<<(std::ostream &s, const MEIGEN &x);

#endif

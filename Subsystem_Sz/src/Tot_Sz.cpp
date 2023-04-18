#include "../include/Tot_Sz.hpp"
#include <vector>
#include <iomanip>
#include <string>

using namespace std;

int Tot_Sz::system_num;
int Tot_Sz::tot_site_A;
int Tot_Sz::tot_site_B;

// コンストラクタ
Tot_Sz::Tot_Sz()
    : bm_A(new int[2]), gbm_A(new int[2]), bm_B(new int[2]), gbm_B(new int[2]), H_int(new Int_Hamiltonian[1]), Eig(2, 2)
{
    Iso_Hamiltonian::Iso_id_counter = 0;

    up_spin_pair[0] = 0;
    up_spin_pair[1] = 0;

    V0 = new double *[2];
    V1 = new double *[2];
    // cout << "Tot_Sz::constructed.\n";
}

void Tot_Sz::sub_count_isonnz()
{
    for (int id = 0; id < 2; id++)
    {
        int nnz = 0;
        int m;
        double szz;

        int Jset_line = H_iso[id].J.get_line();
        int *bm;
        if (id == 0)
            bm = bm_A;
        else
            bm = bm_B;

        for (int bm_itr = 0; bm_itr < H_iso[id].mat_dim; bm_itr++)
        {
            m = bm[bm_itr];
            szz = 0.;
            for (int site_i = 0; site_i < H_iso[id].tot_site_num; site_i++)
            {
                boost::dynamic_bitset<> ket_m(H_iso[id].tot_site_num, m);
                bool is_up_i, is_up_j; // up -> 1(true)、down -> 0

                int j_line = 0;
                for (int l = j_line; l < Jset_line; l++)
                {
                    if (H_iso[id].J.index(0, l) == site_i)
                    {
                        int site_j = H_iso[id].J.index(1, l);
                        is_up_i = ket_m.test(site_i);
                        is_up_j = ket_m.test(site_j);

                        // 注目する2site(site_i, site_j)のスピン状態で場合分け
                        if (is_up_i == is_up_j)
                        { // スピンの向きが揃っている場合
                            szz += 0.25 * H_iso[id].J.val(l);
                        }
                        else
                        {
                            // 昇降演算子を作用させたことによる状態の遷移先が部分空間内に存在するかをチェックする
                            boost::dynamic_bitset<> ket_m1(H_iso[id].tot_site_num, m);
                            ket_m1.flip(site_i);
                            ket_m1.flip(site_j);

                            int n = (int)ket_m1.to_ulong();

                            if (id == 0 && n <= gbm_A_size)
                            {
                                nnz++; // S^+S^- + S^-S^+のぶん
                            }
                            else if (id == 1 && n <= gbm_B_size)
                            {
                                nnz++; // S^+S^- + S^-S^+のぶん
                            }
                            szz -= 0.25 * H_iso[id].J.val(l);
                        }
                        j_line++;
                    }
                } // l
            }     // site_i
            if (szz != 0)
                nnz++;
        } // bm_itr
        H_iso[id].nnz = nnz;
    }
};

// S_A^+、S_A^-、S_B^+、S^B-の配列の要素数をカウントする(S^zは不要)
void Tot_Sz::sub_count_intnnz(const int No, const int pair_num, const int system_num, const int siteA_num, const int siteB_num)
{
    for (int id = 0; id < system_num - 2; id++)
    {
        int nnz = 0;
        int m_A, m_B;

        int siteA_i = H_int[id].J.index(0, 0);
        int siteB_i = H_int[id].J.index(1, 0);

        bool is_exist_sminsA, is_exist_splusA;
        bool is_exist_sminsB, is_exist_splusB;

        /*--------------------system A--------------------*/
        // Tot_Sz[No].bm_AはS^+　、S^-の列番号に対応する
        for (int bm_itr = 0; bm_itr < bm_A_size; bm_itr++)
        {
            m_A = bm_A[bm_itr];
            boost::dynamic_bitset<> ket_col_A(siteA_num, m_A); // mは列番号に対応するスピン状態

            /*site_iAがup状態 -> nnz_mA++*/
            /*site_iAがdown状態 ^> nnz_pA++*/
            // 系Aのsiteについて、S^-行列を定義でき、かつ要素が非ゼロである判定する(非ゼロならtrue)
            is_exist_sminsA = ket_col_A.test(siteA_i) && (No < pair_num - 1);
            // 系AのsiteについてS^+行列が定義でき、かつ要素が非ゼロではないか判定する
            is_exist_splusA = !(ket_col_A.test(siteA_i)) && (No > 0);

            // 条件を満たす場合にのみnnz_pA、nnz_mAをそれぞれカウント
            if (is_exist_sminsA)
                H_int[id].nnz_mA++;
            if (is_exist_splusA)
                H_int[id].nnz_pA++;
        }

        /*--------------------system B--------------------*/
        // Tot_Sz[No].bm_BはS^+、S^-の行番号に対応する　
        for (int bm_itr = 0; bm_itr < bm_B_size; bm_itr++)
        {
            m_B = bm_B[bm_itr];
            boost::dynamic_bitset<> bra_row_B(siteB_num, m_B);

            //<m|S^- = (S^+|m>)^{\dag} -> site_iがdown状態ならnnz_mB++
            //<m|S^+ = (S^-|m>)^(\dag) -> site_iがupスピン状態ならnnz_pB++
            // 系BのsiteについてS^-を定義でき、かつ要素が非ゼロであるか判定する(非ゼロならtrue)
            is_exist_sminsB = !(bra_row_B.test(siteB_i)) && (No < pair_num - 1);
            // 系BのsiteについてS^+を定義でき、かつ要素が非ゼロであるか判定する(非ゼロならtrue)
            is_exist_splusB = bra_row_B.test(siteB_i) && (No > 0);

            // 条件を満たす場合にnnz_pB、nnz_mBをカウント
            if (is_exist_sminsB)
                H_int[id].nnz_mB++;
            if (is_exist_splusB)
                H_int[id].nnz_pB++;
        }
    }
}
// tot_Sz[No]での閉じた系A,BのHamiltonian行列要素を計算する
void Tot_Sz::sub_iso_hamiltonian(const int siteA_num, const int siteB_num)
{
    int row_index = 0;
    int col_ptr_index = 0;
    int col_ptr_val = 0;
    double szz;

    /*--------------------system A--------------------*/
    for (int bm_itr = 0; bm_itr < bm_A_size; bm_itr++)
    {
        szz = 0.;
        H_iso[0].col_ptr[col_ptr_index] = col_ptr_val;
        col_ptr_index++;

        for (int site_i = 0; site_i < siteA_num; site_i++)
        {
            H_iso[0].iso_sub_spin(bm_itr, bm_A, gbm_A, gbm_A_size, siteA_num, site_i, row_index, col_ptr_val, szz);
        }
        if (szz != 0.0)
        {
            H_iso[0].row[row_index] = bm_itr;
            H_iso[0].val[row_index] = szz;
            row_index++;
            col_ptr_val++;
        }
    }
    H_iso[0].col_ptr[H_iso[0].mat_dim] = H_iso[0].nnz;

    /*--------------------system B--------------------*/
    row_index = 0;
    col_ptr_index = 0;
    col_ptr_val = 0;

    for (int bm_itr = 0; bm_itr < bm_B_size; bm_itr++)
    {
        szz = 0.;
        H_iso[1].col_ptr[col_ptr_index] = col_ptr_val;
        col_ptr_index++;

        for (int site_i = 0; site_i < siteB_num; site_i++)
        {
            H_iso[1].iso_sub_spin(bm_itr, bm_B, gbm_B, gbm_B_size, siteB_num, site_i, row_index, col_ptr_val, szz);
        }
        if (szz != 0.0)
        {
            H_iso[1].row[row_index] = bm_itr;
            H_iso[1].val[row_index] = szz;
            row_index++;
            col_ptr_val++;
        }
    }
    H_iso[1].col_ptr[H_iso[1].mat_dim] = H_iso[1].nnz;
};

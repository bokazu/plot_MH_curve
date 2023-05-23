#include "../include/Subsystem_Sz.hpp"
#include <stdio.h>
#include <chrono>
#include <map>
#include <boost/dynamic_bitset.hpp>
using namespace std;

// コンストラクタ
Subsystem_Sz::Subsystem_Sz(int sys_num, int site_A, int site_B, vector<string> &file, int up)
    : system_num(sys_num), tot_site_A(site_A), tot_site_B(site_B), up_spin(up),
      block_row_num(1 << tot_site_A), block_col_num(1 << tot_site_B), pair_num(count_No()), tot_Sz(new Tot_Sz[pair_num]), ls_count(0), ls_check(false), eigen_value(0.), filename(file)
{
  mag = up_spin - (tot_site_A + tot_site_B) / 2.;
  Tot_Sz::system_num = system_num;
  for (int No = 0; No < pair_num; No++)
  {
    delete[] tot_Sz[No].H_int;
    tot_Sz[No].H_int = new Int_Hamiltonian[system_num - 2];
    Int_Hamiltonian::Int_id_counter = 0; // H_intオブジェクトのid番号を正しくカウントするために必要
  }

  Tot_Sz::tot_site_A = tot_site_A;
  Tot_Sz::tot_site_B = tot_site_B;

  /*---------------------------Jsetファイル　のset---------------------------*/

  //@iso_system
  for (int No = 0; No < pair_num; No++)
  {
    for (int id = 0; id < 2; id++)
    {
      tot_Sz[No].H_iso[id].J.set_filename(filename[id]);
      tot_Sz[No].H_iso[id].J.set_line(tot_Sz[No].H_iso[id].J.count_lines());
      tot_Sz[No].H_iso[id].J.resize(tot_Sz[No].H_iso[id].J.get_line());
      tot_Sz[No].H_iso[id].J.set();
    }
  }
  //@int system
  // 必要な情報をset
  for (int No = 0; No < pair_num; No++)
  {
    for (int id = 0; id < system_num - 2; id++)
    {
      tot_Sz[No].H_int[id].J.set_filename(filename[2]);
      tot_Sz[No].H_int[id].J.set_line(1); // 対応するbondの本数
      tot_Sz[No].H_int[id].J.resize(1);   // bond1本文のデータしか読み込まないので1つぶんだけメモリ確保
    }
  }

  // fileから情報を読み込む
  ifstream ifs(filename[2]);
  for (int No = 0; No < pair_num; No++)
  {
    // ファイルポインタを先頭に戻す
    ifs.seekg(ios::beg);
    for (int id = 0; id < system_num - 2; id++)
    {
      ifs >> tot_Sz[No].H_int[id].J.J_index[0][0] >> tot_Sz[No].H_int[id].J.J_index[1][0] >> tot_Sz[No].H_int[id].J.J_val[0];
    }
    // EOFフラグをクリアし
    ifs.clear();
  }

  ifs.close();
  /*--------------------------------------------------------------------*/
}

//(totS_A^zm, tot_S_B^z)のペア数のカウントを行うための関数
int Subsystem_Sz::count_No()
{
  int count = 0;
  for (int up_A = 0; up_A <= tot_site_A; up_A++)
  {
    for (int up_B = 0; up_B <= tot_site_B; up_B++)
    {
      if (up_A + up_B == up_spin)
        count++;
    }
  }
  return count;
};

void Subsystem_Sz::set_system_info()
{
  Tot_Sz::tot_site_A = tot_site_A;
  for (int No = 0; No < pair_num; No++)
  {
    tot_Sz[No].H_iso[0].tot_site_num = tot_site_A;
    tot_Sz[No].H_iso[1].tot_site_num = tot_site_B;
    tot_Sz[No].H_iso[0].mat_dim = tot_Sz[No].bm_A_size;
    tot_Sz[No].H_iso[1].mat_dim = tot_Sz[No].bm_B_size;
  }
  Int_Hamiltonian::tot_site_A = tot_site_A;
  Int_Hamiltonian::tot_site_B = tot_site_B;

  for (int No = 0; No < pair_num; No++)
  {
    for (int id = 0; id < system_num - 2; id++)
    {
      tot_Sz[No].H_int[id].dim_A = tot_Sz[No].bm_A_size;
      tot_Sz[No].H_int[id].dim_B = tot_Sz[No].bm_B_size;
    }
  }
}

void Subsystem_Sz::print_up_pair()
{
  cout << "--------------------------------------\n";
  cout << "No."
       << " "
       << "up_A"
       << " "
       << "up_B" << endl;
  cout << "--------------------------------------\n";
  for (int No = 0; No < pair_num; No++)
  {
    cout << No << " " << tot_Sz[No].up_spin_pair[0] << " " << tot_Sz[No].up_spin_pair[1] << endl;
  }
  cout << "--------------------------------------\n";
}

void Subsystem_Sz::print_subspace()
{
  for (int No = 0; No < pair_num; No++)
  {
    cout << "-----------------------------\n";
    cout << "@ No = " << No << endl;
    cout << "-----------------------------\n";
    cout << "Spin state of A = \n";
    cout << "[";
    for (int i = 0; i < tot_Sz[No].bm_A_size; i++)
    {
      cout << tot_Sz[No].bm_A[i] << " , ";
    }
    cout << "]\n";
    cout << "Spin state of B = \n";
    cout << "[";
    for (int i = 0; i < tot_Sz[No].bm_B_size; i++)
    {
      cout << tot_Sz[No].bm_B[i] << " , ";
    }
    cout << "]\n";
    cout << "-----------------------------\n";
  }
}

void Subsystem_Sz::pirnt_id()
{
  for (int No = 0; No < pair_num; No++)
  {
    cout << "-----------------------------\n";
    cout << "@ No = " << No << endl;
    cout << "-----------------------------\n";
    cout << "#Iso system\n";
    for (int id = 0; id < 2; id++)
    {
      cout << "H_iso[" << tot_Sz[No].H_iso[id].iso_id() << "]" << endl;
    }
    cout << "#Int system\n";
    for (int id = 0; id < system_num - 2; id++)
    {
      cout << "H_int[" << tot_Sz[No].H_int[id].int_id() << "]" << endl;
    }
  }
}

void Subsystem_Sz::print_member()
{
  cout << "======================================================================================================================" << endl;
  cout << "              Member of Subsystem_Sz                       " << endl;
  cout << "======================================================================================================================" << endl;
  cout << "- system_num : " << system_num << endl;
  cout << "- tot_site_A : " << tot_site_A << endl;
  cout << "- tot_site_B : " << tot_site_B << endl;
  for (int i = 0; i < system_num; i++)
  {
    cout << "filename[" << i << "] : " << filename[i] << endl;
  }
  cout << "- up spin   : " << up_spin << endl;
  cout << "- mag       : " << mag << endl;
  cout << "- block_row_num : " << block_row_num << endl;
  cout << "- block_col_num : " << block_col_num << endl;
  cout << "- pair_num  : " << pair_num << endl;
  cout << "- ls_count  : " << ls_count << endl;
  cout << "- ls_check  : " << ls_check << endl;
  cout << "- eigen value : " << eigen_value << endl;
  for (int No = 0; No < pair_num; No++)
  {
    cout << "---------------------------------------------------------------------------------------------------------------------- " << endl;
    cout << "@tot_Szz[" << No << "]" << endl;
    cout << "---------------------------------------------------------------------------------------------------------------------- " << endl;
    cout << "- system_num : " << Tot_Sz::system_num << endl;
    cout << "- tot_site_A : " << Tot_Sz::tot_site_A << endl;
    cout << "- tot_site_B : " << Tot_Sz::tot_site_B << endl;
    cout << "- up_spin_pair[0] : " << tot_Sz[No].up_spin_pair[0] << endl;
    cout << "- up_spin_pair[1] : " << tot_Sz[No].up_spin_pair[1] << endl;
    cout << "- bm_A =" << endl;
    cout << "[";
    for (int i = 0; i < tot_Sz[No].bm_A_size; i++)
      cout << tot_Sz[No].bm_A[i] << ",";
    cout << "]" << endl;
    cout << "- bm_A_size : " << tot_Sz[No].bm_A_size << endl;
    cout << "- gbm_A =" << endl;
    cout << "[";
    for (int i = 0; i < tot_Sz[No].gbm_A_size; i++)
      cout << tot_Sz[No].gbm_A[i] << ",";
    cout << "]" << endl;
    cout << "- gbm_A_size : " << tot_Sz[No].gbm_A_size << endl;
    cout << "- bm_B =" << endl;
    cout << "[";
    for (int i = 0; i < tot_Sz[No].bm_B_size; i++)
      cout << tot_Sz[No].bm_B[i] << ",";
    cout << "]" << endl;
    cout << "- bm_B_size : " << tot_Sz[No].bm_B_size << endl;
    cout << "- gbm_B =" << endl;
    cout << "[";
    for (int i = 0; i < tot_Sz[No].gbm_B_size; i++)
      cout << tot_Sz[No].gbm_B[i] << ",";
    cout << "]" << endl;
    cout << "- gbm_B_size : " << tot_Sz[No].gbm_B_size << endl;
    cout << "---------------------------------------------------------------------------------------------------------------------- " << endl;
    cout << "#Iso Hamiltonian" << endl;
    cout << "---------------------------------------------------------------------------------------------------------------------- " << endl;
    for (int id = 0; id < 2; id++)
    {
      cout << "- ID : " << tot_Sz[No].H_iso[id].Iso_id << endl;
      cout << "- tot_site_num : " << tot_Sz[No].H_iso[id].tot_site_num << endl;
      cout << "- mad_dim : " << tot_Sz[No].H_iso[id].mat_dim << endl;
      cout << "- nnz : " << tot_Sz[No].H_iso[id].nnz << endl;
      tot_Sz[No].H_iso[id].J.print();
    }
    cout << "---------------------------------------------------------------------------------------------------------------------- " << endl;
    cout << "---------------------------------------------------------------------------------------------------------------------- " << endl;
    cout << "#Int Hamiltonian" << endl;
    cout << "---------------------------------------------------------------------------------------------------------------------- " << endl;
    for (int id = 0; id < system_num - 2; id++)
    {
      cout << "- ID : " << tot_Sz[No].H_int[id].Int_id << endl;
      cout << "- tot_site_A : " << tot_Sz[No].H_int[id].tot_site_A << endl;
      cout << "- tot_site_B : " << tot_Sz[No].H_int[id].tot_site_B << endl;
      cout << "- dim_A : " << tot_Sz[No].H_int[id].dim_A << endl;
      cout << "- dim_B : " << tot_Sz[No].H_int[id].dim_B << endl;
      cout << "- nnz_pA : " << tot_Sz[No].H_int[id].nnz_pA << endl;
      cout << "- nnz_pB : " << tot_Sz[No].H_int[id].nnz_pB << endl;
      cout << "- nnz_mA : " << tot_Sz[No].H_int[id].nnz_mA << endl;
      cout << "- nnz_mB : " << tot_Sz[No].H_int[id].nnz_mB << endl;
      tot_Sz[No].H_int[id].J.print();
    }
    cout << "---------------------------------------------------------------------------------------------------------------------- " << endl;
  }
  cout << "======================================================================================================================" << endl;
}

void Subsystem_Sz::check_up_pair()
{
  int No = 0;
  for (int up_A = tot_site_A; up_A >= 0; up_A--)
  {
    for (int up_B = 0; up_B <= tot_site_B; up_B++)
    {
      if (up_A + up_B == up_spin)
      {
        tot_Sz[No].up_spin_pair[0] = up_A;
        tot_Sz[No].up_spin_pair[1] = up_B;
        No++;
      }
    }
  }
}

// 部分空間を張るスピン状態番号を確認する
void Subsystem_Sz::sub_space_check()
{
  check_up_pair();
  // bm_A、bm_Bのメモリ再確保
  for (int No = 0; No < pair_num; No++)
  {
    delete[] tot_Sz[No].bm_A;
    delete[] tot_Sz[No].bm_B;

    tot_Sz[No].bm_A_size = comb(tot_site_A, tot_Sz[No].up_spin_pair[0]);
    tot_Sz[No].bm_B_size = comb(tot_site_B, tot_Sz[No].up_spin_pair[1]);

    tot_Sz[No].bm_A = new int[tot_Sz[No].bm_A_size];
    tot_Sz[No].bm_B = new int[tot_Sz[No].bm_B_size];
  }

  for (int No = 0; No < pair_num; No++)
  {
    delete[] tot_Sz[No].gbm_A;
    delete[] tot_Sz[No].gbm_B;

    sub_space_check_A(No);
    sub_space_check_B(No);
  }
};

void Subsystem_Sz::sub_space_check_A(const int No)
{
  int itr = 0;
  bool flag = false;
  // gbm_Aのサイズ確認と配列のメモリ確保
  for (int m = block_row_num - 1; m >= 0; m--)
  {
    boost::dynamic_bitset<> ket_m(Subsystem_Sz::tot_site_A, m);

    if (ket_m.count() == tot_Sz[No].up_spin_pair[0])
    { //|m>のup spinの本数が条件と一致するか
      if (itr == 0)
      { // 一番最初に条件と一致するものが現れたときの処理
        tot_Sz[No].gbm_A_size = m + 1;
        tot_Sz[No].gbm_A = new int[tot_Sz[No].gbm_A_size];
        vec_init(tot_Sz[No].gbm_A_size, tot_Sz[No].gbm_A);
        flag = true;
        itr++;
      }
    }
    if (flag == true)
      break;
  } // m

  // bm_A, gbm_Aへの要素の格納
  itr = 0;
  for (int m = 0; m < tot_Sz[No].gbm_A_size; m++)
  {
    boost::dynamic_bitset<> ket_m(Subsystem_Sz::tot_site_A, m);
    if (ket_m.count() == tot_Sz[No].up_spin_pair[0])
    {
      tot_Sz[No].bm_A[itr] = m;
      tot_Sz[No].gbm_A[m] = itr;
      itr++;
    }
  }

  if (tot_Sz[No].bm_A[tot_Sz[No].bm_A_size - 1] != tot_Sz[No].gbm_A_size - 1)
  {
    cout << "error::bm & gbm bad arrays @subspace_check_A(). " << endl;
  }
};

void Subsystem_Sz::sub_space_check_B(const int No)
{
  int itr = 0;
  bool flag = false;

  // gbm_Bのサイズ確認と配列のメモリ確保
  for (int m = block_col_num - 1; m >= 0; m--)
  {
    boost::dynamic_bitset<> ket_m(Subsystem_Sz::tot_site_B, m);

    if (ket_m.count() == tot_Sz[No].up_spin_pair[1])
    { //|m>のup spinの本数が条件と一致するか
      if (itr == 0)
      { // 一番最初に条件と一致するものが現れたときの処理
        tot_Sz[No].gbm_B_size = m + 1;
        tot_Sz[No].gbm_B = new int[tot_Sz[No].gbm_B_size];
        vec_init(tot_Sz[No].gbm_B_size, tot_Sz[No].gbm_B);
        flag = true;
        itr++;
      }
    }
    if (flag == true)
      break;
  } // m

  // bm_B, gbm_Bへの要素の格納
  itr = 0;
  for (int m = 0; m < tot_Sz[No].gbm_B_size; m++)
  {
    boost::dynamic_bitset<> ket_m(Subsystem_Sz::tot_site_B, m);
    if (ket_m.count() == tot_Sz[No].up_spin_pair[1])
    {
      tot_Sz[No].bm_B[itr] = m;
      tot_Sz[No].gbm_B[m] = itr;
      itr++;
    }
  }

  if (tot_Sz[No].bm_B[tot_Sz[No].bm_B_size - 1] != tot_Sz[No].gbm_B_size - 1)
  {
    cout << "error::bm & gbm bad arrays @subspace_check() B. " << endl;
  }
};

// Hamiltonian行列の非ゼロ要素数のカウントを行う
void Subsystem_Sz::sub_count_nnz(const int No)
{
  tot_Sz[No].sub_count_isonnz();
  tot_Sz[No].sub_count_intnnz(No, pair_num, system_num, tot_site_A, tot_site_B);
};

// S_A^+とS_B^+の配列の確保&0初期化
void Subsystem_Sz::calloc_splus(const int No, const int id)
{
  // prow_ind_A
  tot_Sz[No].H_int[id].prow_ind_A = new int[tot_Sz[No].H_int[id].nnz_pA];
  for (int i = 0; i < tot_Sz[No].H_int[id].nnz_pA; i++)
    tot_Sz[No].H_int[id].prow_ind_A[i] = 0;
  // prow_ind_B
  tot_Sz[No].H_int[id].prow_ind_B = new int[tot_Sz[No].H_int[id].nnz_pB];
  for (int i = 0; i < tot_Sz[No].H_int[id].nnz_pB; i++)
    tot_Sz[No].H_int[id].prow_ind_B[i] = 0;

  // pcol_ind_A
  tot_Sz[No].H_int[id].pcol_ind_A = new int[tot_Sz[No].H_int[id].nnz_pA];
  for (int i = 0; i < tot_Sz[No].H_int[id].nnz_pA; i++)
    tot_Sz[No].H_int[id].pcol_ind_A[i] = 0;
  // pcol_ind_B
  tot_Sz[No].H_int[id].pcol_ind_B = new int[tot_Sz[No].H_int[id].nnz_pB];
  for (int i = 0; i < tot_Sz[No].H_int[id].nnz_pB; i++)
    tot_Sz[No].H_int[id].pcol_ind_B[i] = 0;
};

// S_A^zとS_B^zの配列の確保&0初期化
void Subsystem_Sz::calloc_szz(const int No, const int id)
{
  // sz_A
  tot_Sz[No].H_int[id].sz_A = new int[tot_Sz[No].bm_A_size];
  for (int i = 0; i < tot_Sz[No].bm_A_size; i++)
    tot_Sz[No].H_int[id].sz_A[i] = 0;
  // sz_B
  tot_Sz[No].H_int[id].sz_B = new int[tot_Sz[No].bm_B_size];
  for (int i = 0; i < tot_Sz[No].bm_B_size; i++)
    tot_Sz[No].H_int[id].sz_B[i] = 0;
};

// S_A^-とS_B^-の配列の確保&0初期化
void Subsystem_Sz::calloc_smins(const int No, const int id)
{
  // mrow_ind_A
  tot_Sz[No].H_int[id].mrow_ind_A = new int[tot_Sz[No].H_int[id].nnz_mA];
  for (int i = 0; i < tot_Sz[No].H_int[id].nnz_mA; i++)
    tot_Sz[No].H_int[id].mrow_ind_A[i] = 0;
  // mrow_ind_B(S_B^-の行はNoの部分空間のスピン状態)
  tot_Sz[No].H_int[id].mrow_ind_B = new int[tot_Sz[No].H_int[id].nnz_mB];
  for (int i = 0; i < tot_Sz[No].H_int[id].nnz_mB; i++)
    tot_Sz[No].H_int[id].mrow_ind_B[i] = 0;

  // mcol_ind_A
  tot_Sz[No].H_int[id].mcol_ind_A = new int[tot_Sz[No].H_int[id].nnz_mA];
  for (int i = 0; i < tot_Sz[No].H_int[id].nnz_mA; i++)
    tot_Sz[No].H_int[id].mcol_ind_A[i] = 0;
  // mcol_ind_B
  tot_Sz[No].H_int[id].mcol_ind_B = new int[tot_Sz[No].H_int[id].nnz_mB];
  for (int i = 0; i < tot_Sz[No].H_int[id].nnz_mB; i++)
    tot_Sz[No].H_int[id].mcol_ind_B[i] = 0;
};

// 配列のメモリへの再確保を行う
void Subsystem_Sz::sub_resize(const int No)
{
  /*----------------------@Iso system----------------------*/
  for (int id = 0; id < 2; id++)
  {
    // メモリの開放
    delete[] tot_Sz[No].H_iso[id].row;
    delete[] tot_Sz[No].H_iso[id].col_ptr;
    delete[] tot_Sz[No].H_iso[id].val;

    // メモリの再確保&初期化
    //- About row
    tot_Sz[No].H_iso[id].row = new int[tot_Sz[No].H_iso[id].nnz];
    for (int i = 0; i < tot_Sz[No].H_iso[id].nnz; i++)
      tot_Sz[No].H_iso[id].row[i] = 0;

    //- About val
    tot_Sz[No].H_iso[id].val = new double[tot_Sz[No].H_iso[id].nnz];
    for (int i = 0; i < tot_Sz[No].H_iso[id].nnz; i++)
      tot_Sz[No].H_iso[id].val[i] = 0.;

    //- About col_ptr
    tot_Sz[No].H_iso[id].col_ptr = new int[tot_Sz[No].H_iso[id].mat_dim + 1];
    for (int i = 0; i < tot_Sz[No].H_iso[id].mat_dim + 1; i++)
      tot_Sz[No].H_iso[id].col_ptr[i] = 0;
  }

  /*----------------------@Int system----------------------*/
  for (int id = 0; id < system_num - 2; id++)
  {
    // メモリの開放

    delete[] tot_Sz[No].H_int[id].sz_A;
    delete[] tot_Sz[No].H_int[id].sz_B;

    // メモリの再確保&初期化
    if (No == 0 && pair_num != 1)
    { // S^-とS^zを求める(S^zは後で計算)  (No=0ではS^+を定義できない)
      delete[] tot_Sz[No].H_int[id].mrow_ind_A;
      delete[] tot_Sz[No].H_int[id].mrow_ind_B;
      delete[] tot_Sz[No].H_int[id].mcol_ind_A;
      delete[] tot_Sz[No].H_int[id].mcol_ind_B;
      calloc_smins(No, id);
    }
    else if (No == pair_num - 1 && pair_num != 1)
    { // S^+とS^zを求める(S^zは後で計算)
      delete[] tot_Sz[No].H_int[id].prow_ind_A;
      delete[] tot_Sz[No].H_int[id].prow_ind_B;
      delete[] tot_Sz[No].H_int[id].pcol_ind_A;
      delete[] tot_Sz[No].H_int[id].pcol_ind_B;
      calloc_splus(No, id);
    }
    else
    { // 上記以外であればS^+とS^-をと主に定義できる
      delete[] tot_Sz[No].H_int[id].prow_ind_A;
      delete[] tot_Sz[No].H_int[id].prow_ind_B;
      delete[] tot_Sz[No].H_int[id].pcol_ind_A;
      delete[] tot_Sz[No].H_int[id].pcol_ind_B;
      delete[] tot_Sz[No].H_int[id].mrow_ind_A;
      delete[] tot_Sz[No].H_int[id].mrow_ind_B;
      delete[] tot_Sz[No].H_int[id].mcol_ind_A;
      delete[] tot_Sz[No].H_int[id].mcol_ind_B;
      calloc_splus(No, id);
      calloc_smins(No, id);
    }
    // S^zはNoによらず定義することができる
    calloc_szz(No, id);
  }
};

void Subsystem_Sz::calc_mat_splus(const int No)
{
  int siteA_i, siteB_i;
  int pA_ind, pB_ind;
  int r_A, c_A, r_B, c_B;
  int bm_rtr_A, bm_rtr_B;
  bool is_up_spin_A, is_up_spin_B;

  for (int id = 0; id < system_num - 2; id++)
  {
    siteA_i = tot_Sz[No].H_int[id].J.index(0, 0); // 系A,BをつなぐbondのA側のサイトを指定
    siteB_i = tot_Sz[No].H_int[id].J.index(1, 0); // 系A,BをつなぐbondのB側のサイトを指定

    pA_ind = 0, pB_ind = 0;
    /*--------------------system A側--------------------*/
    for (int bm_itr = 0; bm_itr < tot_Sz[No].bm_A_size; bm_itr++)
    { // 列番号の走査
      c_A = tot_Sz[No].bm_A[bm_itr];
      boost::dynamic_bitset<> ket_m(tot_site_A, c_A); // c_Aは列番号に対応したスピン状態番号

      // site_iがupスピンならS_i^+|c_A> = 0 -> upスピンの場合を除外する
      is_up_spin_A = ket_m.test(siteA_i);
      if (!is_up_spin_A)
      {
        tot_Sz[No].H_int[id].pcol_ind_A[pA_ind] = bm_itr; // 非ゼロ要素の列番号を記録
        ket_m.flip(siteA_i);
        r_A = (int)(ket_m.to_ulong());        // スピン状態を取得する
        bm_rtr_A = tot_Sz[No - 1].gbm_A[r_A]; // 対応するスピン状態を行番号に変換

        tot_Sz[No].H_int[id].prow_ind_A[pA_ind] = bm_rtr_A; // 非ゼロ要素の行番号を記録
        pA_ind++;
      }
    }

    /*--------------------system B側--------------------*/
    // 注)列番号はtot_Sz[No-1]のスピン状態番号となっている
    for (int bm_itr = 0; bm_itr < tot_Sz[No - 1].bm_B_size; bm_itr++)
    { // 列番号を走査
      c_B = tot_Sz[No - 1].bm_B[bm_itr];
      boost::dynamic_bitset<> ket_m(tot_site_B, c_B);

      // siteB_iがup状態ならS^+|c_B> = 0 -> upスピン状態の場合を除外する
      is_up_spin_B = ket_m.test(siteB_i);
      if (!is_up_spin_B)
      {
        tot_Sz[No].H_int[id].pcol_ind_B[pB_ind] = bm_itr;
        ket_m.flip(siteB_i);
        r_B = (int)(ket_m.to_ulong());
        bm_rtr_B = tot_Sz[No].gbm_B[r_B];

        tot_Sz[No].H_int[id].prow_ind_B[pB_ind] = bm_rtr_B;
        pB_ind++;
      }
    }
  }
}

void Subsystem_Sz::calc_mat_smins(const int No)
{
  int siteA_i, siteB_i;
  int mA_ind, mB_ind;
  int r_A, c_A, r_B, c_B;
  int bm_rtr_A, bm_rtr_B;
  bool is_up_spin_A, is_up_spin_B;

  // idのループは関数の外で回してもいいかもしれない
  for (int id = 0; id < system_num - 2; id++)
  {
    siteA_i = tot_Sz[No].H_int[id].J.index(0, 0);
    siteB_i = tot_Sz[No].H_int[id].J.index(1, 0);

    mA_ind = 0, mB_ind = 0;
    /*--------------------system A側--------------------*/
    for (int bm_itr = 0; bm_itr < tot_Sz[No].bm_A_size; bm_itr++)
    { // 列番号の走査
      c_A = tot_Sz[No].bm_A[bm_itr];
      boost::dynamic_bitset<> ket_m(tot_site_A, c_A);

      // site_iがupスピンならS_i^-|c_A>は非ゼロ -> upスピンの場合のみを考える
      is_up_spin_A = ket_m.test(siteA_i);
      if (is_up_spin_A)
      {
        tot_Sz[No].H_int[id].mcol_ind_A[mA_ind] = bm_itr; // 非ゼロ要素の列番号を記録
        ket_m.flip(siteA_i);
        r_A = (int)(ket_m.to_ulong());
        bm_rtr_A = tot_Sz[No + 1].gbm_A[r_A]; // 対応するスピン状態番号に変換

        tot_Sz[No].H_int[id].mrow_ind_A[mA_ind] = bm_rtr_A; // 非ゼロ要素の行番号を記録
        mA_ind++;
      }
    }

    /*--------------------system B側--------------------*/
    for (int bm_itr = 0; bm_itr < tot_Sz[No + 1].bm_B_size; bm_itr++)
    { // 列番号の走査
      c_B = tot_Sz[No + 1].bm_B[bm_itr];
      boost::dynamic_bitset<> ket_m(tot_site_B, c_B);

      // siteB_iがupスピン状態ならS_i^-|c_B>は非ゼロ -> upスピンの場合のみを考える
      is_up_spin_B = ket_m.test(siteB_i);
      if (is_up_spin_B)
      {
        tot_Sz[No].H_int[id].mcol_ind_B[mB_ind] = bm_itr;
        ket_m.flip(siteB_i);
        r_B = (int)(ket_m.to_ulong());
        bm_rtr_B = tot_Sz[No].gbm_B[r_B];

        tot_Sz[No].H_int[id].mrow_ind_B[mB_ind] = bm_rtr_B;
        mB_ind++;
      }
    }
  }
}

void Subsystem_Sz::calc_mat_szz(const int No)
{
  int siteA_i, siteB_i;
  int r_A, r_B;
  bool is_up_spin_A, is_up_spin_B;
  for (int id = 0; id < system_num - 2; id++)
  {
    siteA_i = tot_Sz[No].H_int[id].J.index(0, 0);
    siteB_i = tot_Sz[No].H_int[id].J.index(1, 0);

    /*--------------------system A側--------------------*/
    for (int bm_itr = 0; bm_itr < tot_Sz[No].bm_A_size; bm_itr++)
    { // 列番号の走査
      r_A = tot_Sz[No].bm_A[bm_itr];
      boost::dynamic_bitset<> ket_m(tot_site_A, r_A);

      is_up_spin_A = ket_m.test(siteA_i);
      // site_iがupスピンならS_i^z|r_A>は正の値
      if (is_up_spin_A)
        tot_Sz[No].H_int[id].sz_A[bm_itr] = 1;
      else
        tot_Sz[No].H_int[id].sz_A[bm_itr] = -1;
    }

    /*--------------------system B側--------------------*/
    for (int bm_itr = 0; bm_itr < tot_Sz[No].bm_B_size; bm_itr++)
    {
      r_B = tot_Sz[No].bm_B[bm_itr];
      boost::dynamic_bitset<> ket_m(tot_site_B, r_B);

      is_up_spin_B = ket_m.test(siteB_i);
      if (is_up_spin_B)
        tot_Sz[No].H_int[id].sz_B[bm_itr] = 1;
      else
        tot_Sz[No].H_int[id].sz_B[bm_itr] = -1;
    }
  }
};

void Subsystem_Sz::sub_int_hamiltonian()
{
  if (pair_num != 1)
  {
    for (int No = 1; No < pair_num; No++) // 動作問題なし
      calc_mat_splus(No);
    for (int No = 0; No < pair_num - 1; No++) // 動作問題なし
      calc_mat_smins(No);
  }

  for (int No = 0; No < pair_num; No++) // 動作問題なし
    calc_mat_szz(No);
};

// Hamiltonian行列の行列要素の計算と配列への格納を行う
void Subsystem_Sz::sub_hamiltonian()
{
  for (int No = 0; No < pair_num; No++)
  {
    sub_count_nnz(No); // 各Tot_SzにおけるHamiltonian行列の非ゼロ要素数をカウントする
    sub_resize(No);    // sub_count_nnz()の結果をもとに、配列サイズを変更し、0初期化を行う

    /*--------------------配列の行列要素を計算する--------------------*/
    tot_Sz[No].sub_iso_hamiltonian(tot_site_A, tot_site_B);
    //[memo]sub_iso_hamiltonianもSubsytem_Szクラスのメンバ関数に後で直す
  }
  sub_int_hamiltonian();
};

// lanczos法
void Subsystem_Sz::sub_lanczos(const int tri_mat_dim, char c, char info_ls)
{
  ls_count = 0;
  double eps = 1.0;
  double err = 1.0e-15;
  bool err_checker = true;

  // 固有ベクトルの用意
  if (c == 'V')
  {
    /*MEIGENの設定*/
    /*コンストラクタで適当に確保してたメモリの開放*/
    for (int No = 0; No < pair_num; No++)
    {
      for (int row = 0; row < tot_Sz[No].Eig.dim_A; row++)
      {
        delete[] tot_Sz[No].Eig.eigen_mat[row];
      }
      delete[] tot_Sz[No].Eig.eigen_mat;
    }
    /*メモリの再確保*/
    // 行列サイズの設定
    for (int No = 0; No < pair_num; No++)
    {
      tot_Sz[No].Eig.dim_A = tot_Sz[No].bm_A_size;
      tot_Sz[No].Eig.dim_B = tot_Sz[No].bm_B_size;
    }
    // メモリの確保
    for (int No = 0; No < pair_num; No++)
    {
      tot_Sz[No].Eig.eigen_mat = new double *[tot_Sz[No].Eig.dim_A];
      for (int row = 0; row < tot_Sz[No].Eig.dim_A; row++)
      {
        tot_Sz[No].Eig.eigen_mat[row] = new double[tot_Sz[No].Eig.dim_B];
      }

      // 0初期化
      tot_Sz[No].Eig.vec_init();
    }
  }

  /*--------------------初期状態行列の用意(あとで関数として用意する)--------------------*/
  for (int No = 0; No < pair_num; No++)
  {
    /*-----コンストラクタで適当に確保していたメモリの開放-----*/
    delete[] tot_Sz[No].V0;
    delete[] tot_Sz[No].V1;

    /*-----メモリの再確保-----*/
    int row_num = tot_Sz[No].bm_A_size;
    int col_num = tot_Sz[No].bm_B_size;
    // V0
    tot_Sz[No].V0 = new double *[row_num];
    for (int i = 0; i < row_num; i++)
      tot_Sz[No].V0[i] = new double[col_num];
    for (int j = 0; j < col_num; j++)
    {
      for (int i = 0; i < row_num; i++)
      {
        tot_Sz[No].V0[i][j] = 0.;
      }
    }
    // V1
    tot_Sz[No].V1 = new double *[row_num];
    for (int i = 0; i < row_num; i++)
      tot_Sz[No].V1[i] = new double[col_num];
    for (int j = 0; j < col_num; j++)
    {
      for (int i = 0; i < row_num; i++)
      {
        tot_Sz[No].V1[i][j] = 0.;
      }
    }
  }

  random_device rand;
  mt19937 mt(rand());
  // mt19937 mt(0x00111111);
  uniform_real_distribution<> rand1(0, 1);
  for (int No = 0; No < pair_num; No++)
  {
    int row_num = tot_Sz[No].bm_A_size;
    int col_num = tot_Sz[No].bm_B_size;
    for (int j = 0; j < col_num; j++)
    {
      for (int i = 0; i < row_num; i++)
      {
        tot_Sz[No].V0[i][j] = rand1(mt);
        tot_Sz[No].V1[i][j] = 0.0;
      }
    }
  }
  mm_sdz_V0(); // 初期状態行列の規格化

  // 初期状態行列の要素を固有ベクトル用の配列にコピーしておく
  if (c == 'V')
  {
    for (int No = 0; No < pair_num; No++)
    {
      mm_dcopy(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].Eig.eigen_mat);
    }
  }

  // 三重対角行列の主対角成分
  double *alpha = new double[tri_mat_dim];
  vec_init(tri_mat_dim, alpha);

  // 三重対角行列の次対角成分
  double *beta = new double[tri_mat_dim - 1];
  vec_init(tri_mat_dim - 1, beta);

  // ls = 偶数stepでの近似固有値
  double *eval_even = new double[tri_mat_dim];
  vec_init(tri_mat_dim, eval_even);

  // ls = 奇数stepでの近似固有値
  double *eval_odd = new double[tri_mat_dim];
  vec_init(tri_mat_dim, eval_odd);

  // LAPACKに三重対角行列の主対角成分を渡す用の配列
  double *diag = new double[tri_mat_dim];
  vec_init(tri_mat_dim, diag);

  // LAPACKに三重対角行列の主対角成分を渡す用の配列
  double *sub_diag = new double[tri_mat_dim - 1];
  vec_init(tri_mat_dim - 1, sub_diag);

  // LAPACKに渡し、c = 'N'なら参照されず、'V'なら固有ベクトルが格納される
  double *tri_diag_evec;

  // 固有ベクトルを計算する場合は配列を確保する
  if (c == 'V')
  {
    tri_diag_evec = new double[tri_mat_dim * tri_mat_dim];
    vec_init(tri_mat_dim * tri_mat_dim, tri_diag_evec);
  }

  /*計算速度チェックのためのセット*/
  ofstream ofs_mmprod("time_mmprod_v0.csv");
  ofs_mmprod << "ls"
             << ","
             << "Time of Iso_mmprod"
             << ","
             << "Time of Int_mmprod"
             << ","
             << "Time of others" << endl;
  chrono::system_clock::time_point start_Iso_mmprod, end_Iso_mmprod, start_Int_mmprod, end_Int_mmprod;
  chrono::system_clock::time_point start_others, end_others;
  auto time_Iso_prod = end_Iso_mmprod - start_Iso_mmprod;
  auto time_Int_prod = end_Int_mmprod - start_Int_mmprod;
  auto time_others = end_others - start_others;

  bool is_odd;
  /*----------------lanczos Algorithm---------------*/
  for (int ls = 0; ls < tri_mat_dim; ls++)
  {
    is_odd = ls % 2; // even -> ls % 2 = 0 -> false, odd -> ls % 2 = 1 -> true
    if (err_checker)
    {
      ls_count = ls;

      // 省メモリのためのlanczosベクトル更新
      if (ls > 0)
      {
        if (is_odd)
        {
          for (int No = 0; No < pair_num; No++)
          {
            mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
          }
        }
        else
        {
          for (int No = 0; No < pair_num; No++)
          {
            mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
          }
        }
      }
      /*========================行列積計算 ＆ α、βの計算==========================*/
      // @odd step
      if (is_odd)
      {
        start_Iso_mmprod = chrono::system_clock::now();
        for (int No = 0; No < pair_num; No++)
        {
          iso_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0);
        }
        end_Iso_mmprod = chrono::system_clock::now();
        time_Iso_prod = end_Iso_mmprod - start_Iso_mmprod;

        start_Int_mmprod = chrono::system_clock::now();
        for (int No = 0; No < pair_num; No++)
        {
          int_mmzzord(No, tot_Sz[No].V1, tot_Sz[No].V0);
          if (pair_num != 1)
          {
            int_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0, tot_Sz[No - 1].V0, tot_Sz[No + 1].V0); // 動作ok
          }
        }
        end_Int_mmprod = chrono::system_clock::now();
        time_Int_prod = end_Int_mmprod - start_Int_mmprod;

        start_others = chrono::system_clock::now();
        for (int No = 0; No < pair_num; No++)
          alpha[ls] += mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].V0);
        calc_beta_oddstep(tri_mat_dim, ls, alpha, beta);

        // lanczosベクトルの更新
        for (int No = 0; No < pair_num; No++)
          mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);

        end_others = chrono::system_clock::now();
        time_others = end_others - start_others;

        ofs_mmprod << ls << "," << chrono::duration_cast<chrono::milliseconds>(time_Iso_prod).count() << ","
                   << chrono::duration_cast<chrono::milliseconds>(time_Int_prod).count() << ","
                   << chrono::duration_cast<chrono::milliseconds>(time_others).count() << endl;
      } // end odd step
      else
      {
        start_Iso_mmprod = chrono::system_clock::now();
        for (int No = 0; No < pair_num; No++)
        {
          iso_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1);
        }
        end_Iso_mmprod = chrono::system_clock::now();
        time_Iso_prod = end_Iso_mmprod - start_Iso_mmprod;

        start_Int_mmprod = chrono::system_clock::now();
        for (int No = 0; No < pair_num; No++)
        {
          int_mmzzord(No, tot_Sz[No].V0, tot_Sz[No].V1);
          if (pair_num != 1)
          {
            int_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1, tot_Sz[No - 1].V1, tot_Sz[No + 1].V1);
          }
        }
        end_Int_mmprod = chrono::system_clock::now();
        time_Int_prod = end_Int_mmprod - start_Int_mmprod;

        start_others = chrono::system_clock::now();
        for (int No = 0; No < pair_num; No++)
          alpha[ls] += mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V1);

        calc_beta_evenstep(tri_mat_dim, ls, alpha, beta);

        // lanczosベクトルの更新
        for (int No = 0; No < pair_num; No++)
          mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);

        end_others = chrono::system_clock::now();
        time_others = end_others - start_others;

        ofs_mmprod << ls << "," << chrono::duration_cast<chrono::milliseconds>(time_Iso_prod).count() << ","
                   << chrono::duration_cast<chrono::milliseconds>(time_Int_prod).count() << ","
                   << chrono::duration_cast<chrono::milliseconds>(time_others).count() << endl;
      } // even step

      /*===========================三重対角行列の数値対角化===========================*/
      vec_init(tri_mat_dim, diag);
      vec_init(tri_mat_dim - 1, sub_diag);
      int info = 0;

      if (is_odd) // odd step
      {
        cblas_dcopy(tri_mat_dim, alpha, 1, diag, 1);
        cblas_dcopy(tri_mat_dim - 1, beta, 1, sub_diag, 1);
        if (ls < tri_mat_dim - 1)
        {
          sub_diag[ls] = 0.;
          if (c == 'N')
          {
            // 固有値のみを計算する場合
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
          }
          else
          {
            // 固有ベクトルのみを計算する場合
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
          }

          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev's error." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        else
        {
          int info =
              LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                            sub_diag, tri_diag_evec, ls + 1);
          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev's error." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        cblas_dcopy(tri_mat_dim, diag, 1, eval_odd, 1);
        if (info_ls == 'y')
        {
          std::cout << "@ls = " << ls
                    << " : eigen value = " << eval_odd[0]
                    << std::endl;
        }
        else if (info_ls == 's')
        {
          cout << "@ls = " << ls << endl;
        }
      } // end of odd step
      else
      {
        // 偶数step
        cblas_dcopy(tri_mat_dim, alpha, 1, diag, 1);
        cblas_dcopy(tri_mat_dim - 1, beta, 1, sub_diag, 1);

        if (ls < tri_mat_dim - 1)
        {
          sub_diag[ls] = 0.;
          if (c == 'N')
          {
            // 固有値のみを計算する場合
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
          }
          else
          { // 固有ベクトルも計算する場合
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
          }

          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev errored." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        else
        {
          if (c == 'N')
          {
            // 固有値のみを計算する場合
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
          }
          else
          {
            // 固有ベクトルを計算する場合
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
          }

          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev errored." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        cblas_dcopy(tri_mat_dim, diag, 1, eval_even, 1);
        if (info_ls == 'y')
        {
          std::cout << "@ls = " << ls
                    << " : eigen value = " << eval_even[0]
                    << std::endl;
        }
        else if (info_ls == 's')
        {
          std::cout << "@ls = " << ls << std::endl;
        }
      } // end of even step
      /*======================================================================*/

      /*============================収束状況の確認==============================*/
      if (ls > 0)
      {
        eps = abs(eval_even[0] - eval_odd[0]);
        if (info_ls == 'y')
        {
          cout << "eps = " << std::setprecision(17) << eps << endl;
        }

        if (eps > err)
          err_checker = true;
        else
        {
          err_checker = false;
          ls_check = true;
        }
      }
      /*=====================================================================*/
    } // err_checer
    else
    {
      cout << "eps = " << eps << endl;
      --ls_count;
      break;
    }
  } // ls
  /*========================基底状態の固有値===========================*/
  if (ls_count % 2 == 0)
    eigen_value = eval_even[0];
  else
    eigen_value = eval_odd[0];

  /*========================配列リソースのリリース part1===================*/
  delete[] eval_even;
  delete[] eval_odd;

  /*======================基底状態の固有ベクトルの計算---------------------*/
  if (c == 'V')
  {
    for (int No = 0; No < pair_num; No++)
    {
      mm_init(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
      mm_init(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
      mm_dcopy(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat, tot_Sz[No].V0);
    }

    for (int ls = 0; ls < ls_count + 2; ls++)
    {
      is_odd = ls % 2;
      if (is_odd)
      {
        for (int No = 0; No < pair_num; No++)
        {
          mm_daxpy(tri_diag_evec[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].Eig.eigen_mat);
        }
      }
      else
      {
        if (ls == 0)
        {
          for (int No = 0; No < pair_num; No++)
          {
            mm_dscal(tri_diag_evec[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat);
          }
        }
        else
        {
          for (int No = 0; No < pair_num; No++)
          {
            mm_daxpy(tri_diag_evec[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].Eig.eigen_mat);
          }
        }
      }

      if (ls > 0)
      {
        if (is_odd)
        {
          for (int No = 0; No < pair_num; No++)
          {
            mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
          }
        }
        else
        {
          for (int No = 0; No < pair_num; No++)
          {
            mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
          }
        }
      }

      if (is_odd)
      {
        for (int No = 0; No < pair_num; No++)
        {
          iso_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0);
          int_mmzzord(No, tot_Sz[No].V1, tot_Sz[No].V0);
          if (pair_num != 1)
          {
            int_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0, tot_Sz[No - 1].V0, tot_Sz[No + 1].V0); // 動作ok
          }
        }
        for (int No = 0; No < pair_num; No++)
        {
          mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].V0);
        }
        // lanczosベクトルの更新
        if (ls != tri_mat_dim - 1)
        {
          for (int No = 0; No < pair_num; No++)
            mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
        }
      }
      else
      {
        for (int No = 0; No < pair_num; No++)
        {
          iso_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1);
          int_mmzzord(No, tot_Sz[No].V0, tot_Sz[No].V1);
          if (pair_num != 1)
          {
            int_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1, tot_Sz[No - 1].V1, tot_Sz[No + 1].V1); // 動作ok
          }                                                                                     // 動作OK
        }
        for (int No = 0; No < pair_num; No++)
        {
          mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V1);
        }
        // lanczosベクトルの更新
        if (ls != tri_mat_dim - 1)
        {
          for (int No = 0; No < pair_num; No++)
            mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
        }
      }
    }
    // lanczosベクトルの規格化
    double tot_norm = 0;
    double dnrm2;
    for (int No = 0; No < pair_num; No++)
    {
      tot_norm += mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat, tot_Sz[No].Eig.eigen_mat);
    }
    dnrm2 = sqrt(tot_norm);
    for (int No = 0; No < pair_num; No++)
    {
      mm_dscal(1. / dnrm2, tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat);
    }

    for (int No = 0; No < pair_num; No++)
    {
      for (int i = 0; i < tot_Sz[No].bm_A_size; i++)
      {
        for (int j = 0; j < tot_Sz[No].bm_B_size; j++)
        {
          double val = tot_Sz[No].Eig.eigen_mat[i][j];
          if (abs(val) > 1e-03)
          {
            cout << setw(18) << left << "tot_Sz[" << No << "].Eig.eigen_mat[" << i << "][" << j << "]" << setprecision(15) << val << endl;
          }
        }
      }
    }
  }
  /*========================配列リソースのリリース part2===================*/
  delete[] alpha;
  delete[] beta;
  delete[] diag;
  delete[] sub_diag;
  if (c == 'V')
    delete[] tri_diag_evec;
  for (int No = 0; No < pair_num; No++)
  {
    for (int i = 0; i < tot_Sz[No].bm_A_size; i++)
    {
      delete[] tot_Sz[No].V0[i];
      delete[] tot_Sz[No].V1[i];
    }
    delete[] tot_Sz[No].V0;
    delete[] tot_Sz[No].V1;
  }

  // Tot_Szのデストラクタでもメモリの開放を行うが、その際にdouble freeエラーがでないように
  // するためにもう一度V0,V1のメモリの確保を行っておく
  for (int No = 0; No < pair_num; No++)
  {
    tot_Sz[No].V0 = new double *[2];
    tot_Sz[No].V1 = new double *[2];
  }
};

// lanczos法 OpenMP利用version
void Subsystem_Sz::MP_sub_lanczos(const int tri_mat_dim, char c, char info_ls)
{
  ls_count = 0;
  double eps = 1.0;
  double err = 1.0e-15;
  bool err_checker = true;

  // 固有ベクトルの用意
  if (c == 'V')
  {
    /*MEIGENの設定*/
    /*コンストラクタで適当に確保してたメモリの開放*/
    for (int No = 0; No < pair_num; No++)
    {
      for (int row = 0; row < tot_Sz[No].Eig.dim_A; row++)
      {
        delete[] tot_Sz[No].Eig.eigen_mat[row];
      }
      delete[] tot_Sz[No].Eig.eigen_mat;
    }
    /*メモリの再確保*/
    // 行列サイズの設定
    for (int No = 0; No < pair_num; No++)
    {
      tot_Sz[No].Eig.dim_A = tot_Sz[No].bm_A_size;
      tot_Sz[No].Eig.dim_B = tot_Sz[No].bm_B_size;
    }
    // メモリの確保
    for (int No = 0; No < pair_num; No++)
    {
      tot_Sz[No].Eig.eigen_mat = new double *[tot_Sz[No].Eig.dim_A];
      for (int row = 0; row < tot_Sz[No].Eig.dim_A; row++)
      {
        tot_Sz[No].Eig.eigen_mat[row] = new double[tot_Sz[No].Eig.dim_B];
      }

      // 0初期化
      tot_Sz[No].Eig.vec_init();
    }
  }

  /*--------------------初期状態行列の用意(あとで関数として用意する)--------------------*/
  int j_init;
  for (int No = 0; No < pair_num; No++)
  {
    /*-----コンストラクタで適当に確保していたメモリの開放-----*/
    delete[] tot_Sz[No].V0;
    delete[] tot_Sz[No].V1;

    /*-----メモリの再確保-----*/
    int row_num = tot_Sz[No].bm_A_size;
    int col_num = tot_Sz[No].bm_B_size;
    // V0
    tot_Sz[No].V0 = new double *[row_num];
    for (int i = 0; i < row_num; i++)
      tot_Sz[No].V0[i] = new double[col_num];

#pragma omp parallel for private(j_init)
    for (int i = 0; i < row_num; i++)
    {
      for (j_init = 0; j_init < col_num; j_init++)
      {
        tot_Sz[No].V0[i][j_init] = 0.;
      }
    }
    // V1
    tot_Sz[No].V1 = new double *[row_num];
    for (int i = 0; i < row_num; i++)
      tot_Sz[No].V1[i] = new double[col_num];

#pragma omp parallel for private(j_init)
    for (int i = 0; i < row_num; i++)
    {
      for (j_init = 0; j_init < col_num; j_init++)
      {
        tot_Sz[No].V1[i][j_init] = 0.;
      }
    }
  }

  random_device rand;
  mt19937 mt(rand());
  uniform_real_distribution<> rand1(0, 1);
  for (int No = 0; No < pair_num; No++)
  {
    int row_num = tot_Sz[No].bm_A_size;
    int col_num = tot_Sz[No].bm_B_size;

    // #pragma omp parallel for private(j_init)
    for (int i = 0; i < row_num; i++)
    {
      for (j_init = 0; j_init < col_num; j_init++)
      {
        tot_Sz[No].V0[i][j_init] = rand1(mt);
        tot_Sz[No].V1[i][j_init] = 0.0;
      }
    }
  }
  MP_mm_sdz_V0(); // 初期状態行列の規格化

  // 初期状態行列の要素を固有ベクトル用の配列にコピーしておく
  if (c == 'V')
  {
    // #pragma omp parallel for
    for (int No = 0; No < pair_num; No++)
    {
      MP_mm_dcopy(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].Eig.eigen_mat);
    }
  }

  // 三重対角行列の主対角成分
  double *alpha = new double[tri_mat_dim];
  MP_vec_init(tri_mat_dim, alpha);

  // 三重対角行列の次対角成分
  double *beta = new double[tri_mat_dim - 1];
  MP_vec_init(tri_mat_dim - 1, beta);

  // ls = 偶数stepでの近似固有値
  double *eval_even = new double[tri_mat_dim];
  MP_vec_init(tri_mat_dim, eval_even);

  // ls = 奇数stepでの近似固有値
  double *eval_odd = new double[tri_mat_dim];
  MP_vec_init(tri_mat_dim, eval_odd);

  // LAPACKに三重対角行列の主対角成分を渡す用の配列
  double *diag = new double[tri_mat_dim];
  MP_vec_init(tri_mat_dim, diag);

  // LAPACKに三重対角行列の主対角成分を渡す用の配列
  double *sub_diag = new double[tri_mat_dim - 1];
  MP_vec_init(tri_mat_dim - 1, sub_diag);

  // LAPACKに渡し、c = 'N'なら参照されず、'V'なら固有ベクトルが格納される
  double *tri_diag_evec;

  // 固有ベクトルを計算する場合は配列を確保する
  if (c == 'V')
  {
    tri_diag_evec = new double[tri_mat_dim * tri_mat_dim];
    MP_vec_init(tri_mat_dim * tri_mat_dim, tri_diag_evec);
  }

  bool is_odd;
  /*----------------lanczos Algorithm---------------*/
  for (int ls = 0; ls < tri_mat_dim; ls++)
  {
    is_odd = ls % 2; // even -> ls % 2 = 0 -> false, odd -> ls % 2 = 1 -> true
    if (err_checker)
    {
      ls_count = ls;

      // 省メモリのためのlanczosベクトル更新
      if (ls > 0)
      {
        if (is_odd)
        {
          for (int No = 0; No < pair_num; No++)
          {
            MP_mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
          }
        }
        else
        {
          for (int No = 0; No < pair_num; No++)
          {
            MP_mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
          }
        }
      }
      /*========================行列積計算 ＆ α、βの計算==========================*/
      // @odd step
      if (is_odd)
      {
        for (int No = 0; No < pair_num; No++)
        {
          MP_iso_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0);
        }

        // #pragma omp parallel for
        for (int No = 0; No < pair_num; No++)
        {
          MP_int_mmzzord(No, tot_Sz[No].V1, tot_Sz[No].V0);
          if (pair_num != 1)
          {
            MP_int_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0, tot_Sz[No - 1].V0, tot_Sz[No + 1].V0); // 動作ok
          }
        }

        MP_calc_alpha_oddstep(ls, alpha);
        MP_calc_beta_oddstep(tri_mat_dim, ls, alpha, beta);

        // lanczosベクトルの更新
        for (int No = 0; No < pair_num; No++)
          MP_mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);

      } // end odd step
      else
      {
        for (int No = 0; No < pair_num; No++)
        {
          MP_iso_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1);
        }

        for (int No = 0; No < pair_num; No++)
        {
          MP_int_mmzzord(No, tot_Sz[No].V0, tot_Sz[No].V1);
          if (pair_num != 1)
          {
            MP_int_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1, tot_Sz[No - 1].V1, tot_Sz[No + 1].V1);
          }
        }

        MP_calc_alpha_evenstep(ls, alpha);
        MP_calc_beta_evenstep(tri_mat_dim, ls, alpha, beta);

        // lanczosベクトルの更新
        for (int No = 0; No < pair_num; No++)
          MP_mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
      } // even step

      /*===========================三重対角行列の数値対角化===========================*/
      MP_vec_init(tri_mat_dim, diag);
      MP_vec_init(tri_mat_dim - 1, sub_diag);
      int info = 0;

      if (is_odd) // odd step
      {
        cblas_dcopy(tri_mat_dim, alpha, 1, diag, 1);
        cblas_dcopy(tri_mat_dim - 1, beta, 1, sub_diag, 1);
        if (ls < tri_mat_dim - 1)
        {
          sub_diag[ls] = 0.;
          if (c == 'N')
          {
            // 固有値のみを計算する場合
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
          }
          else
          {
            // 固有ベクトルを計算する場合
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
          }

          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev's error." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        else
        {
          int info =
              LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                            sub_diag, tri_diag_evec, ls + 1);

          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev's error." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        cblas_dcopy(tri_mat_dim, diag, 1, eval_odd, 1);
        if (info_ls == 'y')
        {
          std::cout << "@ls = " << ls
                    << " : eigen value = " << eval_odd[0]
                    << std::endl;
        }
        else if (info_ls == 's')
        {
          cout << "@ls = " << ls << endl;
        }
      } // end of odd step
      else
      {
        // 偶数step
        cblas_dcopy(tri_mat_dim, alpha, 1, diag, 1);
        cblas_dcopy(tri_mat_dim - 1, beta, 1, sub_diag, 1);

        if (ls < tri_mat_dim - 1)
        {
          sub_diag[ls] = 0.;
          if (c == 'N')
          {
            // 固有値のみを計算する場合
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
          }
          else
          { // 固有ベクトルも計算する場合
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
          }

          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev errored." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        else
        {
          if (c == 'N')
          {
            // 固有値のみを計算する場合
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
          }
          else
          {
            // 固有ベクトルを計算する場合
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
          }

          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev errored." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        cblas_dcopy(tri_mat_dim, diag, 1, eval_even, 1);
        if (info_ls == 'y')
        {
          std::cout << "@ls = " << ls
                    << " : eigen value = " << eval_even[0]
                    << std::endl;
        }
        else if (info_ls == 's')
        {
          std::cout << "@ls = " << ls << std::endl;
        }
      } // end of even step
      /*======================================================================*/

      /*============================収束状況の確認==============================*/
      if (ls > 0)
      {
        eps = abs(eval_even[0] - eval_odd[0]);
        if (info_ls == 'y')
        {
          cout << "eps = " << std::setprecision(17) << eps << endl;
        }

        if (eps > err)
          err_checker = true;
        else
        {
          err_checker = false;
          ls_check = true;
        }
      }
      /*=====================================================================*/
    } // err_checer
    else
    {
      cout << "eps = " << eps << endl;
      --ls_count;
      break;
    }
  } // ls
  /*========================基底状態の固有値===========================*/
  if (ls_count % 2 == 0)
    eigen_value = eval_even[0];
  else
    eigen_value = eval_odd[0];

  /*========================配列リソースのリリース part1===================*/
  delete[] eval_even;
  delete[] eval_odd;

  /*======================基底状態の固有ベクトルの計算---------------------*/
  if (c == 'V')
  {
    for (int No = 0; No < pair_num; No++)
    {
      MP_mm_init(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
      MP_mm_init(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
      MP_mm_dcopy(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat, tot_Sz[No].V0);
    }

    for (int ls = 0; ls < ls_count + 2; ls++) // ls_count+2->ls_count
    {
      is_odd = ls % 2;
      if (is_odd)
      {
        for (int No = 0; No < pair_num; No++)
        {
          MP_mm_daxpy(tri_diag_evec[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].Eig.eigen_mat);
        }
      }
      else
      {
        if (ls == 0)
        {
          for (int No = 0; No < pair_num; No++)
          {
            MP_mm_dscal(tri_diag_evec[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat);
          }
        }
        else
        {
          for (int No = 0; No < pair_num; No++)
          {
            MP_mm_daxpy(tri_diag_evec[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].Eig.eigen_mat);
          }
        }
      }

      if (ls > 0)
      {
        if (is_odd)
        {
          for (int No = 0; No < pair_num; No++)
          {
            MP_mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
          }
        }
        else
        {
          for (int No = 0; No < pair_num; No++)
          {
            MP_mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
          }
        }
      }

      if (is_odd)
      {
        for (int No = 0; No < pair_num; No++)
        {
          MP_iso_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0);
          MP_int_mmzzord(No, tot_Sz[No].V1, tot_Sz[No].V0);
          if (pair_num != 1)
          {
            MP_int_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0, tot_Sz[No - 1].V0, tot_Sz[No + 1].V0); // 動作ok
          }
        }
        for (int No = 0; No < pair_num; No++)
        {
          MP_mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].V0);
        }
        // lanczosベクトルの更新
        if (ls != tri_mat_dim - 1)
        {
          for (int No = 0; No < pair_num; No++)
            MP_mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
        }
      }
      else
      {
        for (int No = 0; No < pair_num; No++)
        {
          MP_iso_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1);
          MP_int_mmzzord(No, tot_Sz[No].V0, tot_Sz[No].V1);
          if (pair_num != 1)
          {
            MP_int_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1, tot_Sz[No - 1].V1, tot_Sz[No + 1].V1); // 動作ok
          }                                                                                        // 動作OK
        }
        for (int No = 0; No < pair_num; No++)
        {
          MP_mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V1);
        }
        // lanczosベクトルの更新
        if (ls != tri_mat_dim - 1)
        {
          for (int No = 0; No < pair_num; No++)
            MP_mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
        }
      }
    }
    // lanczosベクトルの規格化
    double tot_norm = 0;
    double inv_dnrm2;
    for (int No = 0; No < pair_num; No++)
    {
      tot_norm += MP_mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat, tot_Sz[No].Eig.eigen_mat);
    }
    inv_dnrm2 = 1. / sqrt(tot_norm);

    for (int No = 0; No < pair_num; No++)
    {
      MP_mm_dscal(inv_dnrm2, tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat);
    }
  }

  /*========================配列リソースのリリース part2===================*/
  delete[] alpha;
  delete[] beta;
  delete[] diag;
  delete[] sub_diag;
  if (c == 'V')
    delete[] tri_diag_evec;
  for (int No = 0; No < pair_num; No++)
  {
    for (int i = 0; i < tot_Sz[No].bm_A_size; i++)
    {
      delete[] tot_Sz[No].V0[i];
      delete[] tot_Sz[No].V1[i];
    }
    delete[] tot_Sz[No].V0;
    delete[] tot_Sz[No].V1;
  }

  // Tot_Szのデストラクタでもメモリの開放を行うが、その際にdouble freeエラーがでないように
  // するためにもう一度V0,V1のメモリの確保を行っておく
  for (int No = 0; No < pair_num; No++)
  {
    tot_Sz[No].V0 = new double *[2];
    tot_Sz[No].V1 = new double *[2];
  }
};

// lanczos法 OpenMP with schedule利用version
void Subsystem_Sz::MP_schedule_sub_lanczos(const int tri_mat_dim, char c, char info_ls)
{
  ls_count = 0;
  double eps = 1.0;
  double err = 1.0e-15;
  bool err_checker = true;

  // 固有ベクトルの用意
  if (c == 'V')
  {
    /*MEIGENの設定*/
    /*コンストラクタで適当に確保してたメモリの開放*/
    for (int No = 0; No < pair_num; No++)
    {
      for (int row = 0; row < tot_Sz[No].Eig.dim_A; row++)
      {
        delete[] tot_Sz[No].Eig.eigen_mat[row];
      }
      delete[] tot_Sz[No].Eig.eigen_mat;
    }
    /*メモリの再確保*/
    // 行列サイズの設定
    for (int No = 0; No < pair_num; No++)
    {
      tot_Sz[No].Eig.dim_A = tot_Sz[No].bm_A_size;
      tot_Sz[No].Eig.dim_B = tot_Sz[No].bm_B_size;
    }
    // メモリの確保
    for (int No = 0; No < pair_num; No++)
    {
      tot_Sz[No].Eig.eigen_mat = new double *[tot_Sz[No].Eig.dim_A];
      for (int row = 0; row < tot_Sz[No].Eig.dim_A; row++)
      {
        tot_Sz[No].Eig.eigen_mat[row] = new double[tot_Sz[No].Eig.dim_B];
      }

      // 0初期化
      tot_Sz[No].Eig.vec_init();
    }
  }

  /*--------------------初期状態行列の用意(あとで関数として用意する)--------------------*/
  int j_init;
  for (int No = 0; No < pair_num; No++)
  {
    /*-----コンストラクタで適当に確保していたメモリの開放-----*/
    delete[] tot_Sz[No].V0;
    delete[] tot_Sz[No].V1;

    /*-----メモリの再確保-----*/
    int row_num = tot_Sz[No].bm_A_size;
    int col_num = tot_Sz[No].bm_B_size;
    // V0
    tot_Sz[No].V0 = new double *[row_num];
    for (int i = 0; i < row_num; i++)
      tot_Sz[No].V0[i] = new double[col_num];

#pragma omp parallel for private(j_init) schedule(runtime)
    for (int i = 0; i < row_num; i++)
    {
      for (j_init = 0; j_init < col_num; j_init++)
      {
        tot_Sz[No].V0[i][j_init] = 0.;
      }
    }
    // V1
    tot_Sz[No].V1 = new double *[row_num];
    for (int i = 0; i < row_num; i++)
      tot_Sz[No].V1[i] = new double[col_num];

#pragma omp parallel for private(j_init) schedule(runtime)
    for (int i = 0; i < row_num; i++)
    {
      for (j_init = 0; j_init < col_num; j_init++)
      {
        tot_Sz[No].V1[i][j_init] = 0.;
      }
    }
  }

  random_device rand;
  mt19937 mt(rand());
  // mt19937 mt(0x00111111);
  uniform_real_distribution<> rand1(0, 1);
  for (int No = 0; No < pair_num; No++)
  {
    int row_num = tot_Sz[No].bm_A_size;
    int col_num = tot_Sz[No].bm_B_size;

#pragma omp parallel for private(j_init)
    for (int i = 0; i < row_num; i++)
    {
      for (j_init = 0; j_init < col_num; j_init++)
      {
        tot_Sz[No].V0[i][j_init] = rand1(mt);
        tot_Sz[No].V1[i][j_init] = 0.0;
      }
    }
  }
  MP_schedule_mm_sdz_V0(); // 初期状態行列の規格化

  // 初期状態行列の要素を固有ベクトル用の配列にコピーしておく
  if (c == 'V')
  {
    // #pragma omp parallel for
    for (int No = 0; No < pair_num; No++)
    {
      MP_schedule_mm_dcopy(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].Eig.eigen_mat);
    }
  }

  // 三重対角行列の主対角成分
  double *alpha = new double[tri_mat_dim];
  MP_schedule_vec_init(tri_mat_dim, alpha);

  // 三重対角行列の次対角成分
  double *beta = new double[tri_mat_dim - 1];
  MP_schedule_vec_init(tri_mat_dim - 1, beta);

  // ls = 偶数stepでの近似固有値
  double *eval_even = new double[tri_mat_dim];
  MP_schedule_vec_init(tri_mat_dim, eval_even);

  // ls = 奇数stepでの近似固有値
  double *eval_odd = new double[tri_mat_dim];
  MP_schedule_vec_init(tri_mat_dim, eval_odd);

  // LAPACKに三重対角行列の主対角成分を渡す用の配列
  double *diag = new double[tri_mat_dim];
  MP_schedule_vec_init(tri_mat_dim, diag);

  // LAPACKに三重対角行列の主対角成分を渡す用の配列
  double *sub_diag = new double[tri_mat_dim - 1];
  MP_schedule_vec_init(tri_mat_dim - 1, sub_diag);

  // LAPACKに渡し、c = 'N'なら参照されず、'V'なら固有ベクトルが格納される
  double *tri_diag_evec;

  // 固有ベクトルを計算する場合は配列を確保する
  if (c == 'V')
  {
    tri_diag_evec = new double[tri_mat_dim * tri_mat_dim];
    MP_schedule_vec_init(tri_mat_dim * tri_mat_dim, tri_diag_evec);
  }

  bool is_odd;
  /*----------------lanczos Algorithm---------------*/
  for (int ls = 0; ls < tri_mat_dim; ls++)
  {
    is_odd = ls % 2; // even -> ls % 2 = 0 -> false, odd -> ls % 2 = 1 -> true
    if (err_checker)
    {
      ls_count = ls;

      // 省メモリのためのlanczosベクトル更新
      if (ls > 0)
      {
        if (is_odd)
        {
          for (int No = 0; No < pair_num; No++)
          {
            MP_schedule_mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
          }
        }
        else
        {
          for (int No = 0; No < pair_num; No++)
          {
            MP_schedule_mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
          }
        }
      }
      /*========================行列積計算 ＆ α、βの計算==========================*/
      // @odd step
      if (is_odd)
      {
        for (int No = 0; No < pair_num; No++)
        {
          MP_schedule_iso_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0);
        }

        for (int No = 0; No < pair_num; No++)
        {
          MP_schedule_int_mmzzord(No, tot_Sz[No].V1, tot_Sz[No].V0);
          if (pair_num != 1)
          {
            MP_schedule_int_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0, tot_Sz[No - 1].V0, tot_Sz[No + 1].V0); // 動作ok
          }
        }

        MP_schedule_calc_alpha_oddstep(ls, alpha);
        MP_schedule_calc_beta_oddstep(tri_mat_dim, ls, alpha, beta);

        // lanczosベクトルの更新
        for (int No = 0; No < pair_num; No++)
          MP_schedule_mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
      } // end odd step
      else
      {
        // #pragma omp parallel for
        for (int No = 0; No < pair_num; No++)
        {
          MP_schedule_iso_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1);
        }
        for (int No = 0; No < pair_num; No++)
        {
          MP_schedule_int_mmzzord(No, tot_Sz[No].V0, tot_Sz[No].V1);
          if (pair_num != 1)
          {
            MP_schedule_int_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1, tot_Sz[No - 1].V1, tot_Sz[No + 1].V1);
          }
        }

        MP_schedule_calc_alpha_evenstep(ls, alpha);
        MP_schedule_calc_beta_evenstep(tri_mat_dim, ls, alpha, beta);

        // lanczosベクトルの更新
        for (int No = 0; No < pair_num; No++)
          MP_schedule_mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
      } // even step

      /*===========================三重対角行列の数値対角化===========================*/
      MP_schedule_vec_init(tri_mat_dim, diag);
      MP_schedule_vec_init(tri_mat_dim - 1, sub_diag);
      int info = 0;

      if (is_odd) // odd step
      {
        cblas_dcopy(tri_mat_dim, alpha, 1, diag, 1);
        cblas_dcopy(tri_mat_dim - 1, beta, 1, sub_diag, 1);
        if (ls < tri_mat_dim - 1)
        {
          sub_diag[ls] = 0.;
          if (c == 'N')
          {
            // 固有値のみを計算する場合
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
          }
          else
          {
            // 固有ベクトルのみを計算する場合
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
          }

          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev's error." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        else
        {
          int info =
              LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                            sub_diag, tri_diag_evec, ls + 1);
          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev's error." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        cblas_dcopy(tri_mat_dim, diag, 1, eval_odd, 1);
        if (info_ls == 'y')
        {
          std::cout << "@ls = " << ls
                    << " : eigen value = " << eval_odd[0]
                    << std::endl;
        }
        else if (info_ls == 's')
        {
          cout << "@ls = " << ls << endl;
        }
      } // end of odd step
      else
      {
        // 偶数step
        cblas_dcopy(tri_mat_dim, alpha, 1, diag, 1);
        cblas_dcopy(tri_mat_dim - 1, beta, 1, sub_diag, 1);

        if (ls < tri_mat_dim - 1)
        {
          sub_diag[ls] = 0.;
          if (c == 'N')
          {
            // 固有値のみを計算する場合
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
          }
          else
          { // 固有ベクトルも計算する場合
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
          }

          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev errored." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        else
        {
          if (c == 'N')
          {
            // 固有値のみを計算する場合
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
          }
          else
          {
            // 固有ベクトルを計算する場合
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
          }

          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev errored." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        cblas_dcopy(tri_mat_dim, diag, 1, eval_even, 1);
        if (info_ls == 'y')
        {
          std::cout << "@ls = " << ls
                    << " : eigen value = " << eval_even[0]
                    << std::endl;
        }
        else if (info_ls == 's')
        {
          std::cout << "@ls = " << ls << std::endl;
        }
      } // end of even step
      /*======================================================================*/

      /*============================収束状況の確認==============================*/
      if (ls > 0)
      {
        eps = abs(eval_even[0] - eval_odd[0]);
        if (info_ls == 'y')
        {
          cout << "eps = " << std::setprecision(17) << eps << endl;
        }

        if (eps > err)
          err_checker = true;
        else
        {
          err_checker = false;
          ls_check = true;
        }
      }
      /*=====================================================================*/
    } // err_checer
    else
    {
      cout << "eps = " << eps << endl;
      --ls_count;
      break;
    }
  } // ls
  /*========================基底状態の固有値===========================*/
  if (ls_count % 2 == 0)
    eigen_value = eval_even[0];
  else
    eigen_value = eval_odd[0];

  /*========================配列リソースのリリース part1===================*/
  delete[] eval_even;
  delete[] eval_odd;

  /*======================基底状態の固有ベクトルの計算---------------------*/
  if (c == 'V')
  {
    for (int No = 0; No < pair_num; No++)
    {
      MP_schedule_mm_init(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
      MP_schedule_mm_init(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
      MP_schedule_mm_dcopy(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat, tot_Sz[No].V0);
    }

    for (int ls = 0; ls < ls_count + 2; ls++)
    {
      is_odd = ls % 2;
      if (is_odd)
      {
        for (int No = 0; No < pair_num; No++)
        {
          MP_schedule_mm_daxpy(tri_diag_evec[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].Eig.eigen_mat);
        }
      }
      else
      {
        if (ls == 0)
        {
          for (int No = 0; No < pair_num; No++)
          {
            MP_schedule_mm_dscal(tri_diag_evec[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat);
          }
        }
        else
        {
          for (int No = 0; No < pair_num; No++)
          {
            MP_schedule_mm_daxpy(tri_diag_evec[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].Eig.eigen_mat);
          }
        }
      }

      if (ls > 0)
      {
        if (is_odd)
        {
          for (int No = 0; No < pair_num; No++)
          {
            MP_schedule_mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
          }
        }
        else
        {
          for (int No = 0; No < pair_num; No++)
          {
            MP_schedule_mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
          }
        }
      }

      if (is_odd)
      {
        for (int No = 0; No < pair_num; No++)
        {
          MP_schedule_iso_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0);
          MP_schedule_int_mmzzord(No, tot_Sz[No].V1, tot_Sz[No].V0);
          if (pair_num != 1)
          {
            MP_schedule_int_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0, tot_Sz[No - 1].V0, tot_Sz[No + 1].V0); // 動作ok
          }
        }
        for (int No = 0; No < pair_num; No++)
        {
          MP_schedule_mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].V0);
        }
        // lanczosベクトルの更新
        if (ls != tri_mat_dim - 1)
        {
          for (int No = 0; No < pair_num; No++)
            MP_schedule_mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
        }
      }
      else
      {
        for (int No = 0; No < pair_num; No++)
        {
          MP_schedule_iso_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1);
          MP_schedule_int_mmzzord(No, tot_Sz[No].V0, tot_Sz[No].V1);
          if (pair_num != 1)
          {
            MP_schedule_int_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1, tot_Sz[No - 1].V1, tot_Sz[No + 1].V1); // 動作ok
          }                                                                                                 // 動作OK
        }
        for (int No = 0; No < pair_num; No++)
        {
          MP_schedule_mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V1);
        }
        // lanczosベクトルの更新
        if (ls != tri_mat_dim - 1)
        {
          for (int No = 0; No < pair_num; No++)
            MP_schedule_mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
        }
      }
    }
    // lanczosベクトルの規格化
    double tot_norm = 0;
    double inv_dnrm2;
    for (int No = 0; No < pair_num; No++)
    {
      tot_norm += MP_schedule_mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat, tot_Sz[No].Eig.eigen_mat);
    }
    inv_dnrm2 = 1. / sqrt(tot_norm);

    for (int No = 0; No < pair_num; No++)
    {
      MP_schedule_mm_dscal(inv_dnrm2, tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat);
    }
  }

  /*========================配列リソースのリリース part2===================*/
  delete[] alpha;
  delete[] beta;
  delete[] diag;
  delete[] sub_diag;
  if (c == 'V')
    delete[] tri_diag_evec;
  for (int No = 0; No < pair_num; No++)
  {
    for (int i = 0; i < tot_Sz[No].bm_A_size; i++)
    {
      delete[] tot_Sz[No].V0[i];
      delete[] tot_Sz[No].V1[i];
    }
    delete[] tot_Sz[No].V0;
    delete[] tot_Sz[No].V1;
  }

  // Tot_Szのデストラクタでもメモリの開放を行うが、その際にdouble freeエラーがでないように
  // するためにもう一度V0,V1のメモリの確保を行っておく
  for (int No = 0; No < pair_num; No++)
  {
    tot_Sz[No].V0 = new double *[2];
    tot_Sz[No].V1 = new double *[2];
  }
};

// lanczos法 OpenMPversionの各関数の処理に要する時間をテストする
// 時間計測にはOpenMPのライブラリ内の関数を使用している
// 乱数のシードは固定している
double Subsystem_Sz::sub_lanczos_timetest(const int tri_mat_dim, std::string output_filename, char c, char info_ls)
{
  double start_lanczos = omp_get_wtime(); // lanczos法全体の処理に要する時間の計測開始地点

  ls_count = 0;
  double eps = 1.0;
  double err = 1.0e-15;
  bool err_checker = true;

  // 各関数の処理に要する時間を各stepごとに書き出すファイルの設定
  ofstream ofs(output_filename);

  double start_iso, end_iso, time_isoprod;                     // iso_mmprodに要する時間[sec]
  double start_int, end_int, time_intprod;                     // int_mmprodに要する時間[sec]
  double start_alpha, end_alpha, time_alpha;                   // alphaの計算に要する時間[sec]
  double start_beta, end_beta, time_beta;                      // betaの計算に要する時間[sec]
  double start_renew, end_renew, time_renew;                   // lanczosベクトルの更新に要する時間[sec]
  double start_LowMemory, end_LowMemory, time_LowMemory;       // メモリ節約のためにLanczosベクトルに施す処理に要する時間[sec]
  double start_diagonalize, end_diagonalize, time_diagonalize; // LAPACKを用いた数値対角化に要する時間[sec]
  double start_1step, end_1step, time_1step_lanczos;           // lanczos法において1step当たりに要する時間

  double total_lanczos; ////lanczos法全体に要する時間

  // 固有ベクトルの用意
  if (c == 'V')
  {
    /*MEIGENの設定*/
    /*コンストラクタで適当に確保してたメモリの開放*/
    for (int No = 0; No < pair_num; No++)
    {
      for (int row = 0; row < tot_Sz[No].Eig.dim_A; row++)
      {
        delete[] tot_Sz[No].Eig.eigen_mat[row];
      }
      delete[] tot_Sz[No].Eig.eigen_mat;
    }
    /*メモリの再確保*/
    // 行列サイズの設定
    for (int No = 0; No < pair_num; No++)
    {
      tot_Sz[No].Eig.dim_A = tot_Sz[No].bm_A_size;
      tot_Sz[No].Eig.dim_B = tot_Sz[No].bm_B_size;
    }
    // メモリの確保
    for (int No = 0; No < pair_num; No++)
    {
      tot_Sz[No].Eig.eigen_mat = new double *[tot_Sz[No].Eig.dim_A];
      for (int row = 0; row < tot_Sz[No].Eig.dim_A; row++)
      {
        tot_Sz[No].Eig.eigen_mat[row] = new double[tot_Sz[No].Eig.dim_B];
      }

      // 0初期化
      tot_Sz[No].Eig.vec_init();
    }
  }

  /*--------------------初期状態行列の用意(あとで関数として用意する)--------------------*/
  int j_init;
  for (int No = 0; No < pair_num; No++)
  {
    /*-----コンストラクタで適当に確保していたメモリの開放-----*/
    delete[] tot_Sz[No].V0;
    delete[] tot_Sz[No].V1;

    /*-----メモリの再確保-----*/
    int row_num = tot_Sz[No].bm_A_size;
    int col_num = tot_Sz[No].bm_B_size;
    // V0
    tot_Sz[No].V0 = new double *[row_num];
    for (int i = 0; i < row_num; i++)
      tot_Sz[No].V0[i] = new double[col_num];

    // #pragma omp parallel for private(j_init) schedule(runtime)
    for (int i = 0; i < row_num; i++)
    {
      for (j_init = 0; j_init < col_num; j_init++)
      {
        tot_Sz[No].V0[i][j_init] = 0.;
      }
    }
    // V1
    tot_Sz[No].V1 = new double *[row_num];
    for (int i = 0; i < row_num; i++)
      tot_Sz[No].V1[i] = new double[col_num];

#pragma omp parallel for private(j_init)
    for (int i = 0; i < row_num; i++)
    {
      for (j_init = 0; j_init < col_num; j_init++)
      {
        tot_Sz[No].V1[i][j_init] = 0.;
      }
    }
  }

  random_device rand;
  // mt19937 mt(rand());  //通常はこちらを使用する
  mt19937 mt(0x00111111); // 乱数のシードを固定する場合はこちらを使用する
  uniform_real_distribution<> rand1(0, 1);
  for (int No = 0; No < pair_num; No++)
  {
    int row_num = tot_Sz[No].bm_A_size;
    int col_num = tot_Sz[No].bm_B_size;

    for (int i = 0; i < row_num; i++)
    {
      for (j_init = 0; j_init < col_num; j_init++)
      {
        tot_Sz[No].V0[i][j_init] = rand1(mt);
        tot_Sz[No].V1[i][j_init] = 0.0;
      }
    }
  }
  MP_mm_sdz_V0(); // 初期状態行列の規格化

  // 初期状態行列の要素を固有ベクトル用の配列にコピーしておく
  if (c == 'V')
  {
    for (int No = 0; No < pair_num; No++)
    {
      MP_mm_dcopy(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].Eig.eigen_mat);
    }
  }

  // 三重対角行列の主対角成分
  double *alpha = new double[tri_mat_dim];
  MP_vec_init(tri_mat_dim, alpha);

  // 三重対角行列の次対角成分
  double *beta = new double[tri_mat_dim - 1];
  MP_vec_init(tri_mat_dim - 1, beta);

  // ls = 偶数stepでの近似固有値
  double *eval_even = new double[tri_mat_dim];
  MP_vec_init(tri_mat_dim, eval_even);

  // ls = 奇数stepでの近似固有値
  double *eval_odd = new double[tri_mat_dim];
  MP_vec_init(tri_mat_dim, eval_odd);

  // LAPACKに三重対角行列の主対角成分を渡す用の配列
  double *diag = new double[tri_mat_dim];
  MP_vec_init(tri_mat_dim, diag);

  // LAPACKに三重対角行列の主対角成分を渡す用の配列
  double *sub_diag = new double[tri_mat_dim - 1];
  MP_vec_init(tri_mat_dim - 1, sub_diag);

  // LAPACKに渡し、c = 'N'なら参照されず、'V'なら固有ベクトルが格納される
  double *tri_diag_evec;

  // 固有ベクトルを計算する場合は配列を確保する
  if (c == 'V')
  {
    tri_diag_evec = new double[tri_mat_dim * tri_mat_dim];
    MP_vec_init(tri_mat_dim * tri_mat_dim, tri_diag_evec);
  }

  ofs << setw(5) << left << "ls"
      << ","
      << setw(15) << "iso"
      << ","
      << setw(15) << "int"
      << ","
      << setw(15) << "alpha"
      << ","
      << setw(15) << "beta"
      << ","
      << setw(15) << "LowMemory"
      << ","
      << setw(15) << "Renew Matrix"
      << ","
      << setw(15) << "Diagonalization"
      << ","
      << setw(15) << "1step total time" << endl;

  bool is_odd;
  /*----------------lanczos Algorithm---------------*/
  for (int ls = 0; ls < tri_mat_dim; ls++)
  {
    start_1step = omp_get_wtime();
    is_odd = ls % 2; // even -> ls % 2 = 0 -> false, odd -> ls % 2 = 1 -> true
    if (err_checker)
    {
      ls_count = ls;

      // 省メモリのためのlanczosベクトル更新
      if (ls > 0)
      {
        if (is_odd)
        {
          start_LowMemory = omp_get_wtime();
          for (int No = 0; No < pair_num; No++)
          {
            mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
          }
          end_LowMemory = omp_get_wtime();
          time_LowMemory = end_LowMemory - start_LowMemory; //[SECONDS]
        }
        else
        {
          start_LowMemory = omp_get_wtime();
          for (int No = 0; No < pair_num; No++)
          {
            mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
          }
          end_LowMemory = omp_get_wtime();
          time_LowMemory = end_LowMemory - start_LowMemory;
        }
      }
      /*========================行列積計算 ＆ α、βの計算==========================*/
      // @odd step
      if (is_odd)
      {
        start_iso = omp_get_wtime();
        for (int No = 0; No < pair_num; No++)
        {
          iso_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0);
        }
        end_iso = omp_get_wtime();
        time_isoprod = end_iso - start_iso;

        start_int = omp_get_wtime();
        for (int No = 0; No < pair_num; No++)
        {
          int_mmzzord(No, tot_Sz[No].V1, tot_Sz[No].V0);
          if (pair_num != 1)
          {
            int_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0, tot_Sz[No - 1].V0, tot_Sz[No + 1].V0); // 動作ok
          }
        }
        end_int = omp_get_wtime();
        time_intprod = end_int - start_int;

        start_alpha = omp_get_wtime();
        calc_alpha_oddstep(ls, alpha);
        end_alpha = omp_get_wtime();
        time_alpha = end_alpha - start_alpha;

        start_beta = omp_get_wtime();
        calc_beta_oddstep(tri_mat_dim, ls, alpha, beta);
        end_beta = omp_get_wtime();
        time_beta = end_beta - start_beta;

        // lanczosベクトルの更新
        start_renew = omp_get_wtime();
        for (int No = 0; No < pair_num; No++)
          mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
        end_renew = omp_get_wtime();
        time_renew = end_renew - start_renew;
      } // end odd step
      else
      {
        start_iso = omp_get_wtime();
        for (int No = 0; No < pair_num; No++)
        {
          iso_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1);
        }
        end_iso = omp_get_wtime();
        time_isoprod = end_iso - start_iso;

        start_int = omp_get_wtime();
        for (int No = 0; No < pair_num; No++)
        {
          int_mmzzord(No, tot_Sz[No].V0, tot_Sz[No].V1);
          if (pair_num != 1)
          {
            int_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1, tot_Sz[No - 1].V1, tot_Sz[No + 1].V1);
          }
        }

        end_int = omp_get_wtime();
        time_intprod = end_int - start_int;

        start_alpha = omp_get_wtime();
        calc_alpha_evenstep(ls, alpha);
        end_alpha = omp_get_wtime();
        time_alpha = end_alpha - start_alpha;

        start_beta = omp_get_wtime();
        calc_beta_evenstep(tri_mat_dim, ls, alpha, beta);
        end_beta = omp_get_wtime();
        time_beta = end_beta - start_beta;

        // lanczosベクトルの更新
        start_renew = omp_get_wtime();
        for (int No = 0; No < pair_num; No++)
          mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
        end_renew = omp_get_wtime();
        time_renew = end_renew - start_renew;
      } // even step

      /*===========================三重対角行列の数値対角化===========================*/
      vec_init(tri_mat_dim, diag);
      vec_init(tri_mat_dim - 1, sub_diag);
      int info = 0;

      if (is_odd) // odd step
      {
        cblas_dcopy(tri_mat_dim, alpha, 1, diag, 1);
        cblas_dcopy(tri_mat_dim - 1, beta, 1, sub_diag, 1);
        if (ls < tri_mat_dim - 1)
        {
          sub_diag[ls] = 0.;
          if (c == 'N')
          {
            // 固有値のみを計算する場合
            start_diagonalize = omp_get_wtime();
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
            end_diagonalize = omp_get_wtime();
            time_diagonalize = end_diagonalize - start_diagonalize;
          }
          else
          {
            // 固有ベクトルのみを計算する場合
            start_diagonalize = omp_get_wtime();
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
            end_diagonalize = omp_get_wtime();
            time_diagonalize = end_diagonalize - start_diagonalize;
          }

          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev's error." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        else
        {
          start_diagonalize = omp_get_wtime();
          int info =
              LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                            sub_diag, tri_diag_evec, ls + 1);
          end_diagonalize = omp_get_wtime();
          time_diagonalize = end_diagonalize - start_diagonalize;
          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev's error." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        cblas_dcopy(tri_mat_dim, diag, 1, eval_odd, 1);
        if (info_ls == 'y')
        {
          std::cout << "@ls = " << ls
                    << " : eigen value = " << eval_odd[0]
                    << std::endl;
        }
        else if (info_ls == 's')
        {
          cout << "@ls = " << ls << endl;
        }
      } // end of odd step
      else
      {
        // 偶数step
        cblas_dcopy(tri_mat_dim, alpha, 1, diag, 1);
        cblas_dcopy(tri_mat_dim - 1, beta, 1, sub_diag, 1);

        if (ls < tri_mat_dim - 1)
        {
          sub_diag[ls] = 0.;
          if (c == 'N')
          {
            // 固有値のみを計算する場合
            start_diagonalize = omp_get_wtime();
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
            end_diagonalize = omp_get_wtime();
            time_diagonalize = end_diagonalize - start_diagonalize;
          }
          else
          { // 固有ベクトルも計算する場合
            start_diagonalize = omp_get_wtime();
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
            end_diagonalize = omp_get_wtime();
            time_diagonalize = end_diagonalize - start_diagonalize;
          }

          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev errored." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        else
        {
          if (c == 'N')
          {
            // 固有値のみを計算する場合
            start_diagonalize = omp_get_wtime();
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
            end_diagonalize = omp_get_wtime();
            time_diagonalize = end_diagonalize - start_diagonalize;
          }
          else
          {
            // 固有ベクトルを計算する場合
            start_diagonalize = omp_get_wtime();
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
            end_diagonalize = omp_get_wtime();
            time_diagonalize = end_diagonalize - start_diagonalize;
          }

          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev errored." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        cblas_dcopy(tri_mat_dim, diag, 1, eval_even, 1);
        if (info_ls == 'y')
        {
          std::cout << "@ls = " << ls
                    << " : eigen value = " << eval_even[0]
                    << std::endl;
        }
        else if (info_ls == 's')
        {
          std::cout << "@ls = " << ls << std::endl;
        }
      } // end of even step

      end_1step = omp_get_wtime();
      time_1step_lanczos = end_1step - start_1step;

      // 1e03倍してmillisecに変換する
      ofs << setw(5) << ls << "," << setw(15) << time_isoprod * 1e03 << ","
          << setw(15) << time_intprod * 1e03 << ","
          << setw(15) << time_alpha * 1e03 << ","
          << setw(15) << time_beta * 1e03 << ","
          << setw(15) << time_LowMemory * 1e03 << ","
          << setw(15) << time_renew * 1e03 << ","
          << setw(15) << time_diagonalize * 1e03 << ","
          << setw(15) << time_1step_lanczos * 1e03 << endl;
      /*======================================================================*/

      /*============================収束状況の確認==============================*/
      if (ls > 0)
      {
        eps = abs(eval_even[0] - eval_odd[0]);
        if (info_ls == 'y')
        {
          cout << "eps = " << std::setprecision(17) << eps << endl;
        }

        if (eps > err)
          err_checker = true;
        else
        {
          err_checker = false;
          ls_check = true;
        }
      }
      /*=====================================================================*/
    } // err_checer
    else
    {
      cout << "eps = " << eps << endl;
      --ls_count;
      break;
    }
  }            // ls
  ofs.close(); // 各処理に要する時間を計測するための時間を計測
  /*========================基底状態の固有値===========================*/
  if (ls_count % 2 == 0)
    eigen_value = eval_even[0];
  else
    eigen_value = eval_odd[0];

  /*========================配列リソースのリリース part1===================*/
  delete[] eval_even;
  delete[] eval_odd;
  double end_lanczos = omp_get_wtime();
  total_lanczos = end_lanczos - start_lanczos;

  /*======================基底状態の固有ベクトルの計算---------------------*/
  if (c == 'V')
  {
    for (int No = 0; No < pair_num; No++)
    {
      mm_init(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
      mm_init(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
      mm_dcopy(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat, tot_Sz[No].V0);
    }

    for (int ls = 0; ls < ls_count + 2; ls++)
    {
      is_odd = ls % 2;
      if (is_odd)
      {
        for (int No = 0; No < pair_num; No++)
        {
          mm_daxpy(tri_diag_evec[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].Eig.eigen_mat);
        }
      }
      else
      {
        if (ls == 0)
        {
          for (int No = 0; No < pair_num; No++)
          {
            mm_dscal(tri_diag_evec[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat);
          }
        }
        else
        {
          for (int No = 0; No < pair_num; No++)
          {
            mm_daxpy(tri_diag_evec[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].Eig.eigen_mat);
          }
        }
      }

      if (ls > 0)
      {
        if (is_odd)
        {
          for (int No = 0; No < pair_num; No++)
          {
            mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
          }
        }
        else
        {
          for (int No = 0; No < pair_num; No++)
          {
            mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
          }
        }
      }

      if (is_odd)
      {
        for (int No = 0; No < pair_num; No++)
        {
          iso_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0);
          int_mmzzord(No, tot_Sz[No].V1, tot_Sz[No].V0);
          if (pair_num != 1)
          {
            int_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0, tot_Sz[No - 1].V0, tot_Sz[No + 1].V0); // 動作ok
          }
        }
        for (int No = 0; No < pair_num; No++)
        {
          mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].V0);
        }
        // lanczosベクトルの更新
        if (ls != tri_mat_dim - 1)
        {
          for (int No = 0; No < pair_num; No++)
            mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
        }
      }
      else
      {
        for (int No = 0; No < pair_num; No++)
        {
          iso_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1);
          int_mmzzord(No, tot_Sz[No].V0, tot_Sz[No].V1);
          if (pair_num != 1)
          {
            int_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1, tot_Sz[No - 1].V1, tot_Sz[No + 1].V1); // 動作ok
          }                                                                                     // 動作OK
        }
        for (int No = 0; No < pair_num; No++)
        {
          mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V1);
        }
        // lanczosベクトルの更新
        if (ls != tri_mat_dim - 1)
        {
          for (int No = 0; No < pair_num; No++)
            mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
        }
      }
    }
    // lanczosベクトルの規格化
    double tot_norm = 0;
    double inv_dnrm2;
    for (int No = 0; No < pair_num; No++)
    {
      tot_norm += mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat, tot_Sz[No].Eig.eigen_mat);
    }
    inv_dnrm2 = 1. / sqrt(tot_norm);

    for (int No = 0; No < pair_num; No++)
    {
      mm_dscal(inv_dnrm2, tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat);
    }
  }

  /*========================配列リソースのリリース part2===================*/
  delete[] alpha;
  delete[] beta;
  delete[] diag;
  delete[] sub_diag;
  if (c == 'V')
    delete[] tri_diag_evec;
  for (int No = 0; No < pair_num; No++)
  {
    for (int i = 0; i < tot_Sz[No].bm_A_size; i++)
    {
      delete[] tot_Sz[No].V0[i];
      delete[] tot_Sz[No].V1[i];
    }
    delete[] tot_Sz[No].V0;
    delete[] tot_Sz[No].V1;
  }

  // Tot_Szのデストラクタでもメモリの開放を行うが、その際にdouble freeエラーがでないように
  // するためにもう一度V0,V1のメモリの確保を行っておく
  for (int No = 0; No < pair_num; No++)
  {
    tot_Sz[No].V0 = new double *[2];
    tot_Sz[No].V1 = new double *[2];
  }
  return total_lanczos;
};

// lanczos法 OpenMPversionの各関数の処理に要する時間をテストする
// 時間計測にはOpenMPのライブラリ内の関数を使用している
// 乱数のシードは固定している
double Subsystem_Sz::MP_sub_lanczos_timetest(const int tri_mat_dim, std::string output_filename, char c, char info_ls)
{
  double start_lanczos = omp_get_wtime(); // lanczos法全体の処理に要する時間の計測開始地点

  ls_count = 0;
  double eps = 1.0;
  double err = 1.0e-15;
  bool err_checker = true;

  // 各関数の処理に要する時間を各stepごとに書き出すファイルの設定
  ofstream ofs(output_filename);

  double start_iso, end_iso, time_isoprod;                     // iso_mmprodに要する時間[sec]
  double start_int, end_int, time_intprod;                     // int_mmprodに要する時間[sec]
  double start_alpha, end_alpha, time_alpha;                   // alphaの計算に要する時間[sec]
  double start_beta, end_beta, time_beta;                      // betaの計算に要する時間[sec]
  double start_renew, end_renew, time_renew;                   // lanczosベクトルの更新に要する時間[sec]
  double start_LowMemory, end_LowMemory, time_LowMemory;       // メモリ節約のためにLanczosベクトルに施す処理に要する時間[sec]
  double start_diagonalize, end_diagonalize, time_diagonalize; // LAPACKを用いた数値対角化に要する時間[sec]
  double start_1step, end_1step, time_1step_lanczos;           // lanczos法において1step当たりに要する時間

  double total_lanczos; ////lanczos法全体に要する時間

  // 固有ベクトルの用意
  if (c == 'V')
  {
    /*MEIGENの設定*/
    /*コンストラクタで適当に確保してたメモリの開放*/
    for (int No = 0; No < pair_num; No++)
    {
      for (int row = 0; row < tot_Sz[No].Eig.dim_A; row++)
      {
        delete[] tot_Sz[No].Eig.eigen_mat[row];
      }
      delete[] tot_Sz[No].Eig.eigen_mat;
    }
    /*メモリの再確保*/
    // 行列サイズの設定
    for (int No = 0; No < pair_num; No++)
    {
      tot_Sz[No].Eig.dim_A = tot_Sz[No].bm_A_size;
      tot_Sz[No].Eig.dim_B = tot_Sz[No].bm_B_size;
    }
    // メモリの確保
    for (int No = 0; No < pair_num; No++)
    {
      tot_Sz[No].Eig.eigen_mat = new double *[tot_Sz[No].Eig.dim_A];
      for (int row = 0; row < tot_Sz[No].Eig.dim_A; row++)
      {
        tot_Sz[No].Eig.eigen_mat[row] = new double[tot_Sz[No].Eig.dim_B];
      }

      // 0初期化
      tot_Sz[No].Eig.vec_init();
    }
  }

  /*--------------------初期状態行列の用意(あとで関数として用意する)--------------------*/
  int j_init;
  for (int No = 0; No < pair_num; No++)
  {
    /*-----コンストラクタで適当に確保していたメモリの開放-----*/
    delete[] tot_Sz[No].V0;
    delete[] tot_Sz[No].V1;

    /*-----メモリの再確保-----*/
    int row_num = tot_Sz[No].bm_A_size;
    int col_num = tot_Sz[No].bm_B_size;
    // V0
    tot_Sz[No].V0 = new double *[row_num];
    for (int i = 0; i < row_num; i++)
      tot_Sz[No].V0[i] = new double[col_num];

    // #pragma omp parallel for private(j_init) schedule(runtime)
    for (int i = 0; i < row_num; i++)
    {
      for (j_init = 0; j_init < col_num; j_init++)
      {
        tot_Sz[No].V0[i][j_init] = 0.;
      }
    }
    // V1
    tot_Sz[No].V1 = new double *[row_num];
    for (int i = 0; i < row_num; i++)
      tot_Sz[No].V1[i] = new double[col_num];

#pragma omp parallel for private(j_init)
    for (int i = 0; i < row_num; i++)
    {
      for (j_init = 0; j_init < col_num; j_init++)
      {
        tot_Sz[No].V1[i][j_init] = 0.;
      }
    }
  }

  random_device rand;
  // mt19937 mt(rand());  //通常はこちらを使用する
  mt19937 mt(0x00111111); // 乱数のシードを固定する場合はこちらを使用する
  uniform_real_distribution<> rand1(0, 1);
  for (int No = 0; No < pair_num; No++)
  {
    int row_num = tot_Sz[No].bm_A_size;
    int col_num = tot_Sz[No].bm_B_size;

    for (int i = 0; i < row_num; i++)
    {
      for (j_init = 0; j_init < col_num; j_init++)
      {
        tot_Sz[No].V0[i][j_init] = rand1(mt);
        tot_Sz[No].V1[i][j_init] = 0.0;
      }
    }
  }
  MP_mm_sdz_V0(); // 初期状態行列の規格化

  // 初期状態行列の要素を固有ベクトル用の配列にコピーしておく
  if (c == 'V')
  {
    for (int No = 0; No < pair_num; No++)
    {
      MP_mm_dcopy(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].Eig.eigen_mat);
    }
  }

  // 三重対角行列の主対角成分
  double *alpha = new double[tri_mat_dim];
  MP_vec_init(tri_mat_dim, alpha);

  // 三重対角行列の次対角成分
  double *beta = new double[tri_mat_dim - 1];
  MP_vec_init(tri_mat_dim - 1, beta);

  // ls = 偶数stepでの近似固有値
  double *eval_even = new double[tri_mat_dim];
  MP_vec_init(tri_mat_dim, eval_even);

  // ls = 奇数stepでの近似固有値
  double *eval_odd = new double[tri_mat_dim];
  MP_vec_init(tri_mat_dim, eval_odd);

  // LAPACKに三重対角行列の主対角成分を渡す用の配列
  double *diag = new double[tri_mat_dim];
  MP_vec_init(tri_mat_dim, diag);

  // LAPACKに三重対角行列の主対角成分を渡す用の配列
  double *sub_diag = new double[tri_mat_dim - 1];
  MP_vec_init(tri_mat_dim - 1, sub_diag);

  // LAPACKに渡し、c = 'N'なら参照されず、'V'なら固有ベクトルが格納される
  double *tri_diag_evec;

  // 固有ベクトルを計算する場合は配列を確保する
  if (c == 'V')
  {
    tri_diag_evec = new double[tri_mat_dim * tri_mat_dim];
    MP_vec_init(tri_mat_dim * tri_mat_dim, tri_diag_evec);
  }

  ofs << setw(5) << left << "ls"
      << ","
      << setw(15) << "iso"
      << ","
      << setw(15) << "int"
      << ","
      << setw(15) << "alpha"
      << ","
      << setw(15) << "beta"
      << ","
      << setw(15) << "LowMemory"
      << ","
      << setw(15) << "Renew Matrix"
      << ","
      << setw(15) << "Diagonalization"
      << ","
      << setw(15) << "1step total time" << endl;

  bool is_odd;
  /*----------------lanczos Algorithm---------------*/
  for (int ls = 0; ls < tri_mat_dim; ls++)
  {
    start_1step = omp_get_wtime();
    is_odd = ls % 2; // even -> ls % 2 = 0 -> false, odd -> ls % 2 = 1 -> true
    if (err_checker)
    {
      ls_count = ls;

      // 省メモリのためのlanczosベクトル更新
      if (ls > 0)
      {
        if (is_odd)
        {
          start_LowMemory = omp_get_wtime();
          for (int No = 0; No < pair_num; No++)
          {
            MP_mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
          }
          end_LowMemory = omp_get_wtime();
          time_LowMemory = end_LowMemory - start_LowMemory; //[SECONDS]
        }
        else
        {
          start_LowMemory = omp_get_wtime();
          for (int No = 0; No < pair_num; No++)
          {
            MP_mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
          }
          end_LowMemory = omp_get_wtime();
          time_LowMemory = end_LowMemory - start_LowMemory;
        }
      }
      /*========================行列積計算 ＆ α、βの計算==========================*/
      // @odd step
      if (is_odd)
      {
        start_iso = omp_get_wtime();
        for (int No = 0; No < pair_num; No++)
        {
          MP_iso_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0);
        }
        end_iso = omp_get_wtime();
        time_isoprod = end_iso - start_iso;

        start_int = omp_get_wtime();
        for (int No = 0; No < pair_num; No++)
        {
          MP_int_mmzzord(No, tot_Sz[No].V1, tot_Sz[No].V0);
          if (pair_num != 1)
          {
            MP_int_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0, tot_Sz[No - 1].V0, tot_Sz[No + 1].V0); // 動作ok
          }
        }
        end_int = omp_get_wtime();
        time_intprod = end_int - start_int;

        start_alpha = omp_get_wtime();
        MP_calc_alpha_oddstep(ls, alpha);
        end_alpha = omp_get_wtime();
        time_alpha = end_alpha - start_alpha;

        start_beta = omp_get_wtime();
        MP_calc_beta_oddstep(tri_mat_dim, ls, alpha, beta);
        end_beta = omp_get_wtime();
        time_beta = end_beta - start_beta;

        // lanczosベクトルの更新
        start_renew = omp_get_wtime();
        for (int No = 0; No < pair_num; No++)
          MP_mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
        end_renew = omp_get_wtime();
        time_renew = end_renew - start_renew;
      } // end odd step
      else
      {
        start_iso = omp_get_wtime();
        for (int No = 0; No < pair_num; No++)
        {
          MP_iso_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1);
        }
        end_iso = omp_get_wtime();
        time_isoprod = end_iso - start_iso;

        start_int = omp_get_wtime();
        for (int No = 0; No < pair_num; No++)
        {
          MP_int_mmzzord(No, tot_Sz[No].V0, tot_Sz[No].V1);
          if (pair_num != 1)
          {
            MP_int_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1, tot_Sz[No - 1].V1, tot_Sz[No + 1].V1);
          }
        }

        end_int = omp_get_wtime();
        time_intprod = end_int - start_int;

        start_alpha = omp_get_wtime();
        MP_calc_alpha_evenstep(ls, alpha);
        end_alpha = omp_get_wtime();
        time_alpha = end_alpha - start_alpha;

        start_beta = omp_get_wtime();
        MP_calc_beta_evenstep(tri_mat_dim, ls, alpha, beta);
        end_beta = omp_get_wtime();
        time_beta = end_beta - start_beta;

        // lanczosベクトルの更新
        start_renew = omp_get_wtime();
        for (int No = 0; No < pair_num; No++)
          MP_mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
        end_renew = omp_get_wtime();
        time_renew = end_renew - start_renew;
      } // even step

      /*===========================三重対角行列の数値対角化===========================*/
      MP_vec_init(tri_mat_dim, diag);
      MP_vec_init(tri_mat_dim - 1, sub_diag);
      int info = 0;

      if (is_odd) // odd step
      {
        cblas_dcopy(tri_mat_dim, alpha, 1, diag, 1);
        cblas_dcopy(tri_mat_dim - 1, beta, 1, sub_diag, 1);
        if (ls < tri_mat_dim - 1)
        {
          sub_diag[ls] = 0.;
          if (c == 'N')
          {
            // 固有値のみを計算する場合
            start_diagonalize = omp_get_wtime();
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
            end_diagonalize = omp_get_wtime();
            time_diagonalize = end_diagonalize - start_diagonalize;
          }
          else
          {
            // 固有ベクトルのみを計算する場合
            start_diagonalize = omp_get_wtime();
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
            end_diagonalize = omp_get_wtime();
            time_diagonalize = end_diagonalize - start_diagonalize;
          }

          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev's error." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        else
        {
          start_diagonalize = omp_get_wtime();
          int info =
              LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                            sub_diag, tri_diag_evec, ls + 1);
          end_diagonalize = omp_get_wtime();
          time_diagonalize = end_diagonalize - start_diagonalize;
          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev's error." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        cblas_dcopy(tri_mat_dim, diag, 1, eval_odd, 1);
        if (info_ls == 'y')
        {
          std::cout << "@ls = " << ls
                    << " : eigen value = " << eval_odd[0]
                    << std::endl;
        }
        else if (info_ls == 's')
        {
          cout << "@ls = " << ls << endl;
        }
      } // end of odd step
      else
      {
        // 偶数step
        cblas_dcopy(tri_mat_dim, alpha, 1, diag, 1);
        cblas_dcopy(tri_mat_dim - 1, beta, 1, sub_diag, 1);

        if (ls < tri_mat_dim - 1)
        {
          sub_diag[ls] = 0.;
          if (c == 'N')
          {
            // 固有値のみを計算する場合
            start_diagonalize = omp_get_wtime();
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
            end_diagonalize = omp_get_wtime();
            time_diagonalize = end_diagonalize - start_diagonalize;
          }
          else
          { // 固有ベクトルも計算する場合
            start_diagonalize = omp_get_wtime();
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
            end_diagonalize = omp_get_wtime();
            time_diagonalize = end_diagonalize - start_diagonalize;
          }

          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev errored." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        else
        {
          if (c == 'N')
          {
            // 固有値のみを計算する場合
            start_diagonalize = omp_get_wtime();
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
            end_diagonalize = omp_get_wtime();
            time_diagonalize = end_diagonalize - start_diagonalize;
          }
          else
          {
            // 固有ベクトルを計算する場合
            start_diagonalize = omp_get_wtime();
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
            end_diagonalize = omp_get_wtime();
            time_diagonalize = end_diagonalize - start_diagonalize;
          }

          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev errored." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        cblas_dcopy(tri_mat_dim, diag, 1, eval_even, 1);
        if (info_ls == 'y')
        {
          std::cout << "@ls = " << ls
                    << " : eigen value = " << eval_even[0]
                    << std::endl;
        }
        else if (info_ls == 's')
        {
          std::cout << "@ls = " << ls << std::endl;
        }
      } // end of even step

      end_1step = omp_get_wtime();
      time_1step_lanczos = end_1step - start_1step;

      // 1e03倍してmillisecに変換する
      ofs << setw(5) << ls << "," << setw(15) << time_isoprod * 1e03 << ","
          << setw(15) << time_intprod * 1e03 << ","
          << setw(15) << time_alpha * 1e03 << ","
          << setw(15) << time_beta * 1e03 << ","
          << setw(15) << time_LowMemory * 1e03 << ","
          << setw(15) << time_renew * 1e03 << ","
          << setw(15) << time_diagonalize * 1e03 << ","
          << setw(15) << time_1step_lanczos * 1e03 << endl;
      /*======================================================================*/

      /*============================収束状況の確認==============================*/
      if (ls > 0)
      {
        eps = abs(eval_even[0] - eval_odd[0]);
        if (info_ls == 'y')
        {
          cout << "eps = " << std::setprecision(17) << eps << endl;
        }

        if (eps > err)
          err_checker = true;
        else
        {
          err_checker = false;
          ls_check = true;
        }
      }
      /*=====================================================================*/
    } // err_checer
    else
    {
      cout << "eps = " << eps << endl;
      --ls_count;
      break;
    }
  }            // ls
  ofs.close(); // 各処理に要する時間を計測するための時間を計測
  /*========================基底状態の固有値===========================*/
  if (ls_count % 2 == 0)
    eigen_value = eval_even[0];
  else
    eigen_value = eval_odd[0];

  /*========================配列リソースのリリース part1===================*/
  delete[] eval_even;
  delete[] eval_odd;
  double end_lanczos = omp_get_wtime();
  total_lanczos = end_lanczos - start_lanczos;

  /*======================基底状態の固有ベクトルの計算---------------------*/
  if (c == 'V')
  {
    for (int No = 0; No < pair_num; No++)
    {
      MP_mm_init(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
      MP_mm_init(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
      MP_mm_dcopy(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat, tot_Sz[No].V0);
    }

    for (int ls = 0; ls < ls_count + 2; ls++)
    {
      is_odd = ls % 2;
      if (is_odd)
      {
        for (int No = 0; No < pair_num; No++)
        {
          MP_mm_daxpy(tri_diag_evec[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].Eig.eigen_mat);
        }
      }
      else
      {
        if (ls == 0)
        {
          for (int No = 0; No < pair_num; No++)
          {
            MP_mm_dscal(tri_diag_evec[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat);
          }
        }
        else
        {
          for (int No = 0; No < pair_num; No++)
          {
            MP_mm_daxpy(tri_diag_evec[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].Eig.eigen_mat);
          }
        }
      }

      if (ls > 0)
      {
        if (is_odd)
        {
          for (int No = 0; No < pair_num; No++)
          {
            MP_mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
          }
        }
        else
        {
          for (int No = 0; No < pair_num; No++)
          {
            MP_mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
          }
        }
      }

      if (is_odd)
      {
        for (int No = 0; No < pair_num; No++)
        {
          MP_iso_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0);
          MP_int_mmzzord(No, tot_Sz[No].V1, tot_Sz[No].V0);
          if (pair_num != 1)
          {
            MP_int_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0, tot_Sz[No - 1].V0, tot_Sz[No + 1].V0); // 動作ok
          }
        }
        for (int No = 0; No < pair_num; No++)
        {
          MP_mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].V0);
        }
        // lanczosベクトルの更新
        if (ls != tri_mat_dim - 1)
        {
          for (int No = 0; No < pair_num; No++)
            MP_mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
        }
      }
      else
      {
        for (int No = 0; No < pair_num; No++)
        {
          MP_iso_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1);
          MP_int_mmzzord(No, tot_Sz[No].V0, tot_Sz[No].V1);
          if (pair_num != 1)
          {
            MP_int_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1, tot_Sz[No - 1].V1, tot_Sz[No + 1].V1); // 動作ok
          }                                                                                        // 動作OK
        }
        for (int No = 0; No < pair_num; No++)
        {
          MP_mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V1);
        }
        // lanczosベクトルの更新
        if (ls != tri_mat_dim - 1)
        {
          for (int No = 0; No < pair_num; No++)
            MP_mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
        }
      }
    }
    // lanczosベクトルの規格化
    double tot_norm = 0;
    double inv_dnrm2;
    for (int No = 0; No < pair_num; No++)
    {
      tot_norm += MP_mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat, tot_Sz[No].Eig.eigen_mat);
    }
    inv_dnrm2 = 1. / sqrt(tot_norm);

    for (int No = 0; No < pair_num; No++)
    {
      MP_mm_dscal(inv_dnrm2, tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat);
    }
  }

  /*========================配列リソースのリリース part2===================*/
  delete[] alpha;
  delete[] beta;
  delete[] diag;
  delete[] sub_diag;
  if (c == 'V')
    delete[] tri_diag_evec;
  for (int No = 0; No < pair_num; No++)
  {
    for (int i = 0; i < tot_Sz[No].bm_A_size; i++)
    {
      delete[] tot_Sz[No].V0[i];
      delete[] tot_Sz[No].V1[i];
    }
    delete[] tot_Sz[No].V0;
    delete[] tot_Sz[No].V1;
  }

  // Tot_Szのデストラクタでもメモリの開放を行うが、その際にdouble freeエラーがでないように
  // するためにもう一度V0,V1のメモリの確保を行っておく
  for (int No = 0; No < pair_num; No++)
  {
    tot_Sz[No].V0 = new double *[2];
    tot_Sz[No].V1 = new double *[2];
  }
  return total_lanczos;
};

// lanczos法 OpenMP with schedule利用versionの各関数の処理に要する時間をテストする
double Subsystem_Sz::MP_schedule_sub_lanczos_timetest(const int tri_mat_dim, std::string output_filename, char c, char info_ls)
{
  double start_lanczos = omp_get_wtime();

  ls_count = 0;
  double eps = 1.0;
  double err = 1.0e-15;
  bool err_checker = true;

  // 各関数の処理に要する時間を各stepごとに書き出すファイルの設定
  ofstream ofs(output_filename);

  double start_iso, end_iso, time_isoprod;                     // iso_mmprodに要する時間[sec]
  double start_int, end_int, time_intprod;                     // int_mmprodに要する時間[sec]
  double start_alpha, end_alpha, time_alpha;                   // alphaの計算に要する時間[sec]
  double start_beta, end_beta, time_beta;                      // betaの計算に要する時間[sec]
  double start_renew, end_renew, time_renew;                   // lanczosベクトルの更新に要する時間[sec]
  double start_LowMemory, end_LowMemory, time_LowMemory;       // メモリ節約のためにLanczosベクトルに施す処理に要する時間[sec]
  double start_diagonalize, end_diagonalize, time_diagonalize; // LAPACKを用いた数値対角化に要する時間[sec]
  double start_1step, end_1step, time_1step_lanczos;           // lanczos法において1step当たりに要する時間

  double total_lanczos; // lanczos法全体に要する時間
  // 固有ベクトルの用意
  if (c == 'V')
  {
    /*MEIGENの設定*/
    /*コンストラクタで適当に確保してたメモリの開放*/
    for (int No = 0; No < pair_num; No++)
    {
      for (int row = 0; row < tot_Sz[No].Eig.dim_A; row++)
      {
        delete[] tot_Sz[No].Eig.eigen_mat[row];
      }
      delete[] tot_Sz[No].Eig.eigen_mat;
    }
    /*メモリの再確保*/
    // 行列サイズの設定
    for (int No = 0; No < pair_num; No++)
    {
      tot_Sz[No].Eig.dim_A = tot_Sz[No].bm_A_size;
      tot_Sz[No].Eig.dim_B = tot_Sz[No].bm_B_size;
    }
    // メモリの確保
    for (int No = 0; No < pair_num; No++)
    {
      tot_Sz[No].Eig.eigen_mat = new double *[tot_Sz[No].Eig.dim_A];
      for (int row = 0; row < tot_Sz[No].Eig.dim_A; row++)
      {
        tot_Sz[No].Eig.eigen_mat[row] = new double[tot_Sz[No].Eig.dim_B];
      }

      // 0初期化
      tot_Sz[No].Eig.vec_init();
    }
  }

  /*--------------------初期状態行列の用意(あとで関数として用意する)--------------------*/
  int j_init;
  for (int No = 0; No < pair_num; No++)
  {
    /*-----コンストラクタで適当に確保していたメモリの開放-----*/
    delete[] tot_Sz[No].V0;
    delete[] tot_Sz[No].V1;

    /*-----メモリの再確保-----*/
    int row_num = tot_Sz[No].bm_A_size;
    int col_num = tot_Sz[No].bm_B_size;
    // V0
    tot_Sz[No].V0 = new double *[row_num];
    for (int i = 0; i < row_num; i++)
      tot_Sz[No].V0[i] = new double[col_num];

    // #pragma omp parallel for private(j_init) schedule(runtime)
    for (int i = 0; i < row_num; i++)
    {
      for (j_init = 0; j_init < col_num; j_init++)
      {
        tot_Sz[No].V0[i][j_init] = 0.;
      }
    }
    // V1
    tot_Sz[No].V1 = new double *[row_num];
    for (int i = 0; i < row_num; i++)
      tot_Sz[No].V1[i] = new double[col_num];

#pragma omp parallel for private(j_init) schedule(runtime)
    for (int i = 0; i < row_num; i++)
    {
      for (j_init = 0; j_init < col_num; j_init++)
      {
        tot_Sz[No].V1[i][j_init] = 0.;
      }
    }
  }

  random_device rand;
  // mt19937 mt(rand());
  mt19937 mt(0x00111111);
  uniform_real_distribution<> rand1(0, 1);
  for (int No = 0; No < pair_num; No++)
  {
    int row_num = tot_Sz[No].bm_A_size;
    int col_num = tot_Sz[No].bm_B_size;

    for (int i = 0; i < row_num; i++)
    {
      for (j_init = 0; j_init < col_num; j_init++)
      {
        tot_Sz[No].V0[i][j_init] = rand1(mt);
        tot_Sz[No].V1[i][j_init] = 0.0;
      }
    }
  }
  MP_schedule_mm_sdz_V0(); // 初期状態行列の規格化

  // 初期状態行列の要素を固有ベクトル用の配列にコピーしておく
  if (c == 'V')
  {
    // #pragma omp parallel for
    for (int No = 0; No < pair_num; No++)
    {
      MP_schedule_mm_dcopy(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].Eig.eigen_mat);
    }
  }

  // 三重対角行列の主対角成分
  double *alpha = new double[tri_mat_dim];
  MP_schedule_vec_init(tri_mat_dim, alpha);

  // 三重対角行列の次対角成分
  double *beta = new double[tri_mat_dim - 1];
  MP_schedule_vec_init(tri_mat_dim - 1, beta);

  // ls = 偶数stepでの近似固有値
  double *eval_even = new double[tri_mat_dim];
  MP_schedule_vec_init(tri_mat_dim, eval_even);

  // ls = 奇数stepでの近似固有値
  double *eval_odd = new double[tri_mat_dim];
  MP_schedule_vec_init(tri_mat_dim, eval_odd);

  // LAPACKに三重対角行列の主対角成分を渡す用の配列
  double *diag = new double[tri_mat_dim];
  MP_schedule_vec_init(tri_mat_dim, diag);

  // LAPACKに三重対角行列の主対角成分を渡す用の配列
  double *sub_diag = new double[tri_mat_dim - 1];
  MP_schedule_vec_init(tri_mat_dim - 1, sub_diag);

  // LAPACKに渡し、c = 'N'なら参照されず、'V'なら固有ベクトルが格納される
  double *tri_diag_evec;

  // 固有ベクトルを計算する場合は配列を確保する
  if (c == 'V')
  {
    tri_diag_evec = new double[tri_mat_dim * tri_mat_dim];
    MP_schedule_vec_init(tri_mat_dim * tri_mat_dim, tri_diag_evec);
  }
  else
  {
    tri_diag_evec = new double[2];
    MP_schedule_vec_init(2, tri_diag_evec);
  }
  // tri_diag_evec = new double[tri_mat_dim * tri_mat_dim];
  // MP_schedule_vec_init(tri_mat_dim * tri_mat_dim, tri_diag_evec);
  ofs << setw(5) << left << "ls"
      << ","
      << setw(15) << "iso"
      << ","
      << setw(15) << "int"
      << ","
      << setw(15) << "alpha"
      << ","
      << setw(15) << "beta"
      << ","
      << setw(15) << "LowMemory"
      << ","
      << setw(15) << "Renew Matrix"
      << ","
      << setw(15) << "Diagonalization"
      << ","
      << setw(15) << "1step total time" << endl;

  bool is_odd;
  /*----------------lanczos Algorithm---------------*/
  for (int ls = 0; ls < tri_mat_dim; ls++)
  {
    start_1step = omp_get_wtime();
    is_odd = ls % 2; // even -> ls % 2 = 0 -> false, odd -> ls % 2 = 1 -> true
    if (err_checker)
    {
      ls_count = ls;

      // 省メモリのためのlanczosベクトル更新
      if (ls > 0)
      {
        if (is_odd)
        {
          start_LowMemory = omp_get_wtime();
          for (int No = 0; No < pair_num; No++)
          {
            MP_schedule_mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
          }
          end_LowMemory = omp_get_wtime();
          time_LowMemory = end_LowMemory - start_LowMemory;
        }
        else
        {
          start_LowMemory = omp_get_wtime();
          for (int No = 0; No < pair_num; No++)
          {
            MP_schedule_mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
          }
          end_LowMemory = omp_get_wtime();
          time_LowMemory = end_LowMemory - start_LowMemory;
        }
      }
      /*========================行列積計算 ＆ α、βの計算==========================*/
      // @odd step
      if (is_odd)
      {
        start_iso = omp_get_wtime();
        for (int No = 0; No < pair_num; No++)
        {
          MP_schedule_iso_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0);
        }
        end_iso = omp_get_wtime();
        time_isoprod = end_iso - start_iso;

        start_int = omp_get_wtime();
        for (int No = 0; No < pair_num; No++)
        {
          MP_schedule_int_mmzzord(No, tot_Sz[No].V1, tot_Sz[No].V0);
          if (pair_num != 1)
          {
            MP_schedule_int_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0, tot_Sz[No - 1].V0, tot_Sz[No + 1].V0); // 動作ok
          }
        }
        end_int = omp_get_wtime();
        time_intprod = end_int - start_int;

        start_alpha = omp_get_wtime();
        MP_schedule_calc_alpha_oddstep(ls, alpha);
        end_alpha = omp_get_wtime();
        time_alpha = end_alpha - start_alpha;

        start_beta = omp_get_wtime();
        MP_schedule_calc_beta_oddstep(tri_mat_dim, ls, alpha, beta);
        end_beta = omp_get_wtime();
        time_beta = end_beta - start_beta;

        // lanczosベクトルの更新
        start_renew = omp_get_wtime();
        for (int No = 0; No < pair_num; No++)
          MP_schedule_mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
        end_renew = omp_get_wtime();
        time_renew = end_renew - start_renew;
      } // end odd step
      else
      {
        start_iso = omp_get_wtime();
        for (int No = 0; No < pair_num; No++)
        {
          MP_schedule_iso_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1);
        }
        end_iso = omp_get_wtime();
        time_isoprod = end_iso - start_iso;

        start_int = omp_get_wtime();
        for (int No = 0; No < pair_num; No++)
        {
          MP_schedule_int_mmzzord(No, tot_Sz[No].V0, tot_Sz[No].V1);
          if (pair_num != 1)
          {
            MP_schedule_int_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1, tot_Sz[No - 1].V1, tot_Sz[No + 1].V1);
          }
        }
        end_int = omp_get_wtime();
        time_intprod = end_int - start_int;

        start_alpha = omp_get_wtime();
        MP_schedule_calc_alpha_evenstep(ls, alpha);
        end_alpha = omp_get_wtime();
        time_alpha = end_alpha - start_alpha;

        start_beta = omp_get_wtime();
        MP_schedule_calc_beta_evenstep(tri_mat_dim, ls, alpha, beta);
        end_beta = omp_get_wtime();
        time_beta = end_beta - start_beta;

        // lanczosベクトルの更新
        start_renew = omp_get_wtime();
        for (int No = 0; No < pair_num; No++)
          MP_schedule_mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
        end_renew = omp_get_wtime();
        time_renew = end_renew - start_renew;
      } // even step

      /*===========================三重対角行列の数値対角化===========================*/
      MP_schedule_vec_init(tri_mat_dim, diag);
      MP_schedule_vec_init(tri_mat_dim - 1, sub_diag);
      int info = 0;

      if (is_odd) // odd step
      {
        cblas_dcopy(tri_mat_dim, alpha, 1, diag, 1);
        cblas_dcopy(tri_mat_dim - 1, beta, 1, sub_diag, 1);
        if (ls < tri_mat_dim - 1)
        {
          sub_diag[ls] = 0.;
          if (c == 'N')
          {
            // 固有値のみを計算する場合
            start_diagonalize = omp_get_wtime();
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);

            end_diagonalize = omp_get_wtime();
            time_diagonalize = end_diagonalize - start_diagonalize;
          }
          else
          {
            // 固有ベクトルのみを計算する場合
            start_diagonalize = omp_get_wtime();
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
            end_diagonalize = omp_get_wtime();
            time_diagonalize = end_diagonalize - start_diagonalize;
          }

          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev's error." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        else
        {
          start_diagonalize = omp_get_wtime();
          int info =
              LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                            sub_diag, tri_diag_evec, ls + 1);
          end_diagonalize = omp_get_wtime();
          time_diagonalize = end_diagonalize - start_diagonalize;
          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev's error." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        cblas_dcopy(tri_mat_dim, diag, 1, eval_odd, 1);
        if (info_ls == 'y')
        {
          std::cout << "@ls = " << ls
                    << " : eigen value = " << eval_odd[0]
                    << std::endl;
        }
        else if (info_ls == 's')
        {
          cout << "@ls = " << ls << endl;
        }
      } // end of odd step
      else
      {
        // 偶数step
        cblas_dcopy(tri_mat_dim, alpha, 1, diag, 1);
        cblas_dcopy(tri_mat_dim - 1, beta, 1, sub_diag, 1);

        if (ls < tri_mat_dim - 1)
        {
          sub_diag[ls] = 0.;
          if (c == 'N')
          {
            // 固有値のみを計算する場合
            start_diagonalize = omp_get_wtime();
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);

            // info =
            //     LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', tri_mat_dim, diag,
            //                   sub_diag, tri_diag_evec, tri_mat_dim);
            end_diagonalize = omp_get_wtime();
            time_diagonalize = end_diagonalize - start_diagonalize;
          }
          else
          { // 固有ベクトルも計算する場合
            start_diagonalize = omp_get_wtime();
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
            end_diagonalize = omp_get_wtime();
            time_diagonalize = end_diagonalize - start_diagonalize;
          }

          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev errored." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        else
        {
          if (c == 'N')
          {
            // 固有値のみを計算する場合
            start_diagonalize = omp_get_wtime();
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'N', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
            end_diagonalize = omp_get_wtime();
            time_diagonalize = end_diagonalize - start_diagonalize;
          }
          else
          {
            // 固有ベクトルを計算する場合
            start_diagonalize = omp_get_wtime();
            info =
                LAPACKE_dstev(LAPACK_COL_MAJOR, 'V', ls + 1, diag,
                              sub_diag, tri_diag_evec, ls + 1);
            end_diagonalize = omp_get_wtime();
            time_diagonalize = end_diagonalize - start_diagonalize;
          }

          if (info != 0)
          {
            std::cout << "@ls = " << ls
                      << " , LAPACKE_detev errored." << std::endl;
            cout << "info = " << info << endl;
          }
        }
        cblas_dcopy(tri_mat_dim, diag, 1, eval_even, 1);
        if (info_ls == 'y')
        {
          std::cout << "@ls = " << ls
                    << " : eigen value = " << eval_even[0]
                    << std::endl;
        }
        else if (info_ls == 's')
        {
          std::cout << "@ls = " << ls << std::endl;
        }
      } // end of even step

      end_1step = omp_get_wtime();
      time_1step_lanczos = end_1step - start_1step;

      // 1e03倍してmillisecに変換する
      ofs << setw(5) << ls << "," << setw(15) << time_isoprod * 1e03 << ","
          << setw(15) << time_intprod * 1e03 << ","
          << setw(15) << time_alpha * 1e03 << ","
          << setw(15) << time_beta * 1e03 << ","
          << setw(15) << time_LowMemory * 1e03 << ","
          << setw(15) << time_renew * 1e03 << ","
          << setw(15) << time_diagonalize * 1e03 << ","
          << setw(15) << time_1step_lanczos * 1e03 << endl;
      /*======================================================================*/

      /*============================収束状況の確認==============================*/
      if (ls > 0)
      {
        eps = abs(eval_even[0] - eval_odd[0]);
        if (info_ls == 'y')
        {
          cout << "eps = " << std::setprecision(17) << eps << endl;
        }

        if (eps > err)
          err_checker = true;
        else
        {
          err_checker = false;
          ls_check = true;
        }
      }
      /*=====================================================================*/
    } // err_checer
    else
    {
      cout << "eps = " << eps << endl;
      --ls_count;
      break;
    }
  } // ls

  ofs.close(); // 各処理に要する時間を出力するためのファイルをclose
  /*========================基底状態の固有値===========================*/
  if (ls_count % 2 == 0)
    eigen_value = eval_even[0];
  else
    eigen_value = eval_odd[0];

  /*========================配列リソースのリリース part1===================*/
  delete[] eval_even;
  delete[] eval_odd;

  double end_lanczos = omp_get_wtime();
  total_lanczos = end_lanczos - start_lanczos;

  /*======================基底状態の固有ベクトルの計算---------------------*/
  if (c == 'V')
  {
    for (int No = 0; No < pair_num; No++)
    {
      MP_schedule_mm_init(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
      MP_schedule_mm_init(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
      MP_schedule_mm_dcopy(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat, tot_Sz[No].V0);
    }

    for (int ls = 0; ls < ls_count + 2; ls++)
    {
      is_odd = ls % 2;
      if (is_odd)
      {
        for (int No = 0; No < pair_num; No++)
        {
          MP_schedule_mm_daxpy(tri_diag_evec[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].Eig.eigen_mat);
        }
      }
      else
      {
        if (ls == 0)
        {
          for (int No = 0; No < pair_num; No++)
          {
            MP_schedule_mm_dscal(tri_diag_evec[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat);
          }
        }
        else
        {
          for (int No = 0; No < pair_num; No++)
          {
            MP_schedule_mm_daxpy(tri_diag_evec[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].Eig.eigen_mat);
          }
        }
      }

      if (ls > 0)
      {
        if (is_odd)
        {
          for (int No = 0; No < pair_num; No++)
          {
            MP_schedule_mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
          }
        }
        else
        {
          for (int No = 0; No < pair_num; No++)
          {
            MP_schedule_mm_dscal(-beta[ls - 1], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
          }
        }
      }

      if (is_odd)
      {
        for (int No = 0; No < pair_num; No++)
        {
          MP_schedule_iso_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0);
          MP_schedule_int_mmzzord(No, tot_Sz[No].V1, tot_Sz[No].V0);
          if (pair_num != 1)
          {
            MP_schedule_int_mmprod(No, tot_Sz[No].V1, tot_Sz[No].V0, tot_Sz[No - 1].V0, tot_Sz[No + 1].V0); // 動作ok
          }
        }
        for (int No = 0; No < pair_num; No++)
        {
          MP_schedule_mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].V0);
        }
        // lanczosベクトルの更新
        if (ls != tri_mat_dim - 1)
        {
          for (int No = 0; No < pair_num; No++)
            MP_schedule_mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
        }
      }
      else
      {
        for (int No = 0; No < pair_num; No++)
        {
          MP_schedule_iso_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1);
          MP_schedule_int_mmzzord(No, tot_Sz[No].V0, tot_Sz[No].V1);
          if (pair_num != 1)
          {
            MP_schedule_int_mmprod(No, tot_Sz[No].V0, tot_Sz[No].V1, tot_Sz[No - 1].V1, tot_Sz[No + 1].V1); // 動作ok
          }                                                                                                 // 動作OK
        }
        for (int No = 0; No < pair_num; No++)
        {
          MP_schedule_mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V1);
        }
        // lanczosベクトルの更新
        if (ls != tri_mat_dim - 1)
        {
          for (int No = 0; No < pair_num; No++)
            MP_schedule_mm_dscal(1. / beta[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
        }
      }
    }
    // lanczosベクトルの規格化
    double tot_norm = 0;
    double dnrm2;
    for (int No = 0; No < pair_num; No++)
    {
      tot_norm += MP_schedule_mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat, tot_Sz[No].Eig.eigen_mat);
    }
    dnrm2 = sqrt(tot_norm);

    for (int No = 0; No < pair_num; No++)
    {
      MP_schedule_mm_dscal(1. / dnrm2, tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].Eig.eigen_mat);
    }
  }
  /*========================配列リソースのリリース part2===================*/
  delete[] alpha;
  delete[] beta;
  delete[] diag;
  delete[] sub_diag;

  delete[] tri_diag_evec;
  for (int No = 0; No < pair_num; No++)
  {
    for (int i = 0; i < tot_Sz[No].bm_A_size; i++)
    {
      delete[] tot_Sz[No].V0[i];
      delete[] tot_Sz[No].V1[i];
    }
    delete[] tot_Sz[No].V0;
    delete[] tot_Sz[No].V1;
  }

  // Tot_Szのデストラクタでもメモリの開放を行うが、その際にdouble freeエラーがでないように
  // するためにもう一度V0,V1のメモリの確保を行っておく
  for (int No = 0; No < pair_num; No++)
  {
    tot_Sz[No].V0 = new double *[2];
    tot_Sz[No].V1 = new double *[2];
  }
  return total_lanczos;
};

void Subsystem_Sz::calc_beta_evenstep(const int tri_mat_dim, const int ls, double *alpha, double *beta)
{
  for (int No = 0; No < pair_num; No++)
  {
    mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V1);
    if (ls != tri_mat_dim - 1)
      beta[ls] += mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].V1);
  }
  if (ls != tri_mat_dim - 1)
  {
    beta[ls] = sqrt(beta[ls]);
  }
};
void Subsystem_Sz::calc_beta_oddstep(const int tri_mat_dim, const int ls, double *alpha, double *beta)
{
  for (int No = 0; No < pair_num; No++)
  {
    mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].V0);
    if (ls != tri_mat_dim - 1)
      beta[ls] += mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V0);
  }
  if (ls != tri_mat_dim - 1)
  {
    beta[ls] = sqrt(beta[ls]);
  }
};

// mm系関数
void Subsystem_Sz::iso_mmprod(const int No, double **V0, double **V1)
{
  tot_Sz[No].H_iso[0].isoA_mmprod(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, V0, V1);
  tot_Sz[No].H_iso[1].isoB_mmprod(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, V0, V1);
}; // Iso systemの行列-行列積

void Subsystem_Sz::MP_iso_mmprod(const int No, double **V0, double **V1)
{
  tot_Sz[No].H_iso[0].MP_isoA_mmprod(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, V0, V1);
  tot_Sz[No].H_iso[1].MP_isoB_mmprod(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, V0, V1);
};

void Subsystem_Sz::MP_schedule_iso_mmprod(const int No, double **V0, double **V1)
{
  tot_Sz[No].H_iso[0].MP_schedule_isoA_mmprod(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, V0, V1);
  tot_Sz[No].H_iso[1].MP_schedule_isoB_mmprod(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, V0, V1);
}

// pair_num != 1のif文中で使用する
void Subsystem_Sz::int_mmprod(const int No, double **V0, double **V1, double **V1_dic1, double **V1_inc1)
{
  bool is_No_0;
  (No == 0) ? is_No_0 = true : is_No_0 = false;
  bool is_No_pairnum_1;
  (No == pair_num - 1) ? is_No_pairnum_1 = true : is_No_pairnum_1 = false;

  for (int id = 0; id < system_num - 2; id++)
  {
    if (is_No_pairnum_1) // No == pair_num-1ではS^+とS^zだけ定義可能
    {
      tot_Sz[No].H_int[id].int_rise_mmprod(tot_Sz[No].H_int[id].nnz_pA, tot_Sz[No].H_int[id].nnz_pB, tot_Sz[No].H_int[id].prow_ind_A, tot_Sz[No].H_int[id].pcol_ind_A, tot_Sz[No].H_int[id].prow_ind_B, tot_Sz[No].H_int[id].pcol_ind_B, V0, V1, V1_dic1);
    }
    else if (is_No_0) // No == 0ではS^-とS^zだけ定義可能
    {
      tot_Sz[No].H_int[id].int_dsmn_mmprod(tot_Sz[No].H_int[id].nnz_mA, tot_Sz[No].H_int[id].nnz_mB, tot_Sz[No].H_int[id].mrow_ind_A, tot_Sz[No].H_int[id].mcol_ind_A, tot_Sz[No].H_int[id].mrow_ind_B, tot_Sz[No].H_int[id].mcol_ind_B, V0, V1, V1_inc1);
    }
    else if (No > 0 && No < pair_num - 1) // 0 < No < pair_num - 1ではS^+,S^-,S^zすべて定義可能
    {
      tot_Sz[No].H_int[id].int_rise_mmprod(tot_Sz[No].H_int[id].nnz_pA, tot_Sz[No].H_int[id].nnz_pB, tot_Sz[No].H_int[id].prow_ind_A, tot_Sz[No].H_int[id].pcol_ind_A, tot_Sz[No].H_int[id].prow_ind_B, tot_Sz[No].H_int[id].pcol_ind_B, V0, V1, V1_dic1);
      // No = 3, id = 0以降でエラーが出ている
      tot_Sz[No].H_int[id].int_dsmn_mmprod(tot_Sz[No].H_int[id].nnz_mA, tot_Sz[No].H_int[id].nnz_mB, tot_Sz[No].H_int[id].mrow_ind_A, tot_Sz[No].H_int[id].mcol_ind_A, tot_Sz[No].H_int[id].mrow_ind_B, tot_Sz[No].H_int[id].mcol_ind_B, V0, V1, V1_inc1); // 動作問題あり
    }
  }
}; /// Int systemの行列-行列-行列積

void Subsystem_Sz::MP_int_mmprod(const int No, double **V0, double **V1, double **V1_dic1, double **V1_inc1)
{
  bool is_No_0;
  (No == 0) ? is_No_0 = true : is_No_0 = false;
  bool is_No_pairnum_1;
  (No == pair_num - 1) ? is_No_pairnum_1 = true : is_No_pairnum_1 = false;

  for (int id = 0; id < system_num - 2; id++)
  {
    if (is_No_pairnum_1) // No == pair_num-1ではS^+とS^zだけ定義可能
    {
      tot_Sz[No].H_int[id].MP_int_rise_mmprod(tot_Sz[No].H_int[id].nnz_pA, tot_Sz[No].H_int[id].nnz_pB, tot_Sz[No].H_int[id].prow_ind_A, tot_Sz[No].H_int[id].pcol_ind_A, tot_Sz[No].H_int[id].prow_ind_B, tot_Sz[No].H_int[id].pcol_ind_B, V0, V1, V1_dic1);
    }
    else if (is_No_0) // No == 0ではS^-とS^zだけ定義可能
    {
      tot_Sz[No].H_int[id].MP_int_dsmn_mmprod(tot_Sz[No].H_int[id].nnz_mA, tot_Sz[No].H_int[id].nnz_mB, tot_Sz[No].H_int[id].mrow_ind_A, tot_Sz[No].H_int[id].mcol_ind_A, tot_Sz[No].H_int[id].mrow_ind_B, tot_Sz[No].H_int[id].mcol_ind_B, V0, V1, V1_inc1);
    }
    else if (No > 0 && No < pair_num - 1) // 0 < No < pair_num - 1ではS^+,S^-,S^zすべて定義可能
    {
      tot_Sz[No].H_int[id].MP_int_rise_mmprod(tot_Sz[No].H_int[id].nnz_pA, tot_Sz[No].H_int[id].nnz_pB, tot_Sz[No].H_int[id].prow_ind_A, tot_Sz[No].H_int[id].pcol_ind_A, tot_Sz[No].H_int[id].prow_ind_B, tot_Sz[No].H_int[id].pcol_ind_B, V0, V1, V1_dic1);
      // No = 3, id = 0以降でエラーが出ている
      tot_Sz[No].H_int[id].MP_int_dsmn_mmprod(tot_Sz[No].H_int[id].nnz_mA, tot_Sz[No].H_int[id].nnz_mB, tot_Sz[No].H_int[id].mrow_ind_A, tot_Sz[No].H_int[id].mcol_ind_A, tot_Sz[No].H_int[id].mrow_ind_B, tot_Sz[No].H_int[id].mcol_ind_B, V0, V1, V1_inc1); // 動作問題あり
    }
  }
};

void Subsystem_Sz::MP_schedule_int_mmprod(const int No, double **V0, double **V1, double **V1_dic1, double **V1_inc1)
{
  bool is_No_0;
  (No == 0) ? is_No_0 = true : is_No_0 = false;
  bool is_No_pairnum_1;
  (No == pair_num - 1) ? is_No_pairnum_1 = true : is_No_pairnum_1 = false;

  for (int id = 0; id < system_num - 2; id++)
  {
    if (is_No_pairnum_1) // No == pair_num-1ではS^+とS^zだけ定義可能
    {
      tot_Sz[No].H_int[id].MP_schedule_int_rise_mmprod(tot_Sz[No].H_int[id].nnz_pA, tot_Sz[No].H_int[id].nnz_pB, tot_Sz[No].H_int[id].prow_ind_A, tot_Sz[No].H_int[id].pcol_ind_A, tot_Sz[No].H_int[id].prow_ind_B, tot_Sz[No].H_int[id].pcol_ind_B, V0, V1, V1_dic1);
    }
    else if (is_No_0) // No == 0ではS^-とS^zだけ定義可能
    {
      tot_Sz[No].H_int[id].MP_schedule_int_dsmn_mmprod(tot_Sz[No].H_int[id].nnz_mA, tot_Sz[No].H_int[id].nnz_mB, tot_Sz[No].H_int[id].mrow_ind_A, tot_Sz[No].H_int[id].mcol_ind_A, tot_Sz[No].H_int[id].mrow_ind_B, tot_Sz[No].H_int[id].mcol_ind_B, V0, V1, V1_inc1);
    }
    else if (No > 0 && No < pair_num - 1) // 0 < No < pair_num - 1ではS^+,S^-,S^zすべて定義可能
    {
      tot_Sz[No].H_int[id].MP_schedule_int_rise_mmprod(tot_Sz[No].H_int[id].nnz_pA, tot_Sz[No].H_int[id].nnz_pB, tot_Sz[No].H_int[id].prow_ind_A, tot_Sz[No].H_int[id].pcol_ind_A, tot_Sz[No].H_int[id].prow_ind_B, tot_Sz[No].H_int[id].pcol_ind_B, V0, V1, V1_dic1);
      // No = 3, id = 0以降でエラーが出ている
      tot_Sz[No].H_int[id].MP_schedule_int_dsmn_mmprod(tot_Sz[No].H_int[id].nnz_mA, tot_Sz[No].H_int[id].nnz_mB, tot_Sz[No].H_int[id].mrow_ind_A, tot_Sz[No].H_int[id].mcol_ind_A, tot_Sz[No].H_int[id].mrow_ind_B, tot_Sz[No].H_int[id].mcol_ind_B, V0, V1, V1_inc1); // 動作問題あり
    }
  }
};

void Subsystem_Sz::int_mmzzord(const int No, double **V0, double **V1)
{
  for (int id = 0; id < system_num - 2; id++)
  {
    tot_Sz[No].H_int[id].int_zz_mmprod(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].H_int[id].sz_A, tot_Sz[No].H_int[id].sz_B, V0, V1); // 動作OK
  }
};

void Subsystem_Sz::MP_int_mmzzord(const int No, double **V0, double **V1)
{
  // #pragma omp parallel for
  for (int id = 0; id < system_num - 2; id++)
  {
    tot_Sz[No].H_int[id].MP_int_zz_mmprod(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].H_int[id].sz_A, tot_Sz[No].H_int[id].sz_B, V0, V1); // 動作OK
  }
}

void Subsystem_Sz::MP_schedule_int_mmzzord(const int No, double **V0, double **V1)
{
  // #pragma omp parallel for
  for (int id = 0; id < system_num - 2; id++)
  {
    tot_Sz[No].H_int[id].MP_schedule_int_zz_mmprod(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].H_int[id].sz_A, tot_Sz[No].H_int[id].sz_B, V0, V1); // 動作OK
  }
}

double Subsystem_Sz::mm_ddot(const int row_dim, const int col_dim, double **V0, double **V1)
{
  double val = 0.;

  for (int j = 0; j < col_dim; j++)
  {
    for (int i = 0; i < row_dim; i++)
    {
      val += V0[i][j] * V1[i][j];
    }
  }
  return val;
};

void Subsystem_Sz::mm_dscal(double alpha, const int row_dim, const int col_dim, double **V)
{
  for (int j = 0; j < col_dim; j++)
  {
    for (int i = 0; i < row_dim; i++)
    {
      V[i][j] *= alpha;
    }
  }
};

void Subsystem_Sz::mm_dcopy(const int row_dim, const int col_dim, double **V0, double **V1)
{
  for (int j = 0; j < col_dim; j++)
  {
    for (int i = 0; i < row_dim; i++)
    {
      V1[i][j] = V0[i][j];
    }
  }
}

void Subsystem_Sz::mm_daxpy(double alpha, const int row_dim, const int col_dim, double **V0, double **V1)
{
  for (int j = 0; j < col_dim; j++)
  {
    for (int i = 0; i < row_dim; i++)
    {
      V1[i][j] += alpha * V0[i][j];
    }
  }
}

double Subsystem_Sz::mm_dnrm2(const int row_dim, const int col_dim, double **V)
{
  double val = 0.;
  for (int j = 0; j < col_dim; j++)
  {
    for (int i = 0; i < row_dim; i++)
    {
      val += V[i][j] * V[i][j];
    }
  }
  return sqrt(val);
}

void Subsystem_Sz::mm_init(const int row_dim, const int col_dim, double **V)
{
  for (int col = 0; col < col_dim; col++)
  {
    for (int row = 0; row < row_dim; row++)
    {
      V[row][col] = 0.0;
    }
  }
}

void Subsystem_Sz::MP_mm_init(const int row_dim, const int col_dim, double **V)
{
  int col;
#pragma omp parallel for private(col)
  for (int row = 0; row < row_dim; row++)
  {
    for (col = 0; col < col_dim; col++)
    {
      V[row][col] = 0.0;
    }
  }
}

void Subsystem_Sz::MP_schedule_mm_init(const int row_dim, const int col_dim, double **V)
{
  int col;
#pragma omp parallel for private(col) schedule(runtime)
  for (int row = 0; row < row_dim; row++)
  {
    for (col = 0; col < col_dim; col++)
    {
      V[row][col] = 0.0;
    }
  }
}

void Subsystem_Sz::mm_sdz(const int row_dim, const int col_dim, double **V)
{
  double a = 1. / mm_dnrm2(row_dim, col_dim, V);
  mm_dscal(a, row_dim, col_dim, V);
}

void Subsystem_Sz::mm_sdz_V0()
{ // lanczos行列V0の規格化を行う
  double tot_norm = 0;
  double dnrm2;
  for (int No = 0; No < pair_num; No++)
  {
    tot_norm += mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V0);
  }
  dnrm2 = sqrt(tot_norm);
  for (int No = 0; No < pair_num; No++)
  {
    mm_dscal(1. / dnrm2, tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
  }
}

void Subsystem_Sz::MP_mm_sdz_V0()
{
  // lanczos行列V0の規格化を行う
  double tot_norm = 0;
  double dnrm2;

  for (int No = 0; No < pair_num; No++)
  {
    tot_norm += MP_mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V0);
  }
  dnrm2 = sqrt(tot_norm);

  for (int No = 0; No < pair_num; No++)
  {
    MP_mm_dscal(1. / dnrm2, tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
  }
}

void Subsystem_Sz::MP_schedule_mm_sdz_V0()
{
  // lanczos行列V0の規格化を行う
  double tot_norm = 0;
  double dnrm2;

  for (int No = 0; No < pair_num; No++)
  {
    tot_norm += MP_schedule_mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V0);
  }
  dnrm2 = sqrt(tot_norm);

  for (int No = 0; No < pair_num; No++)
  {
    MP_schedule_mm_dscal(1. / dnrm2, tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0);
  }
}

void Subsystem_Sz::mm_sdz_V1()
{
  double tot_norm = 0;
  double dnrm2;
  for (int No = 0; No < pair_num; No++)
  {
    tot_norm += mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V1);
  }
  dnrm2 = sqrt(tot_norm);
  for (int No = 0; No < pair_num; No++)
  {
    mm_dscal(1. / dnrm2, tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
  }
}

void Subsystem_Sz::MP_mm_sdz_V1()
{
  double tot_norm = 0;
  double dnrm2;

  for (int No = 0; No < pair_num; No++)
  {
    tot_norm += MP_mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V1);
  }
  dnrm2 = sqrt(tot_norm);

  for (int No = 0; No < pair_num; No++)
  {
    MP_mm_dscal(1. / dnrm2, tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
  }
}

void Subsystem_Sz::MP_schedule_mm_sdz_V1()
{
  double tot_norm = 0;
  double dnrm2;

  for (int No = 0; No < pair_num; No++)
  {
    tot_norm += MP_schedule_mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V1);
  }
  dnrm2 = sqrt(tot_norm);

  for (int No = 0; No < pair_num; No++)
  {
    MP_schedule_mm_dscal(1. / dnrm2, tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1);
  }
}

/*----------上記関数のOpenMP利用version----------*/
double Subsystem_Sz::MP_mm_ddot(const int row_dim, const int col_dim, double **V0, double **V1)
{
  int j;
  double val = 0.;
#pragma omp parallel for private(j) reduction(+ : val)
  for (int i = 0; i < row_dim; i++)
  {
    for (j = 0; j < col_dim; j++)
    {
      val += V0[i][j] * V1[i][j];
    }
  }
  return val;
};

double Subsystem_Sz::MP_schedule_mm_ddot(const int row_dim, const int col_dim, double **V0, double **V1)
{
  int j;
  double val = 0.;
#pragma omp parallel for private(j) reduction(+ : val) schedule(runtime)
  for (int i = 0; i < row_dim; i++)
  {
    for (j = 0; j < col_dim; j++)
    {
      val += V0[i][j] * V1[i][j];
    }
  }
  return val;
};

void Subsystem_Sz::MP_mm_dscal(double alpha, const int row_dim, const int col_dim, double **V)
{
  int j;
#pragma omp parallel for private(j)
  for (int i = 0; i < row_dim; i++)
  {
    for (j = 0; j < col_dim; j++)
    {
      V[i][j] *= alpha;
    }
  }
};

void Subsystem_Sz::MP_schedule_mm_dscal(double alpha, const int row_dim, const int col_dim, double **V)
{
  int j;
#pragma omp parallel for private(j) schedule(runtime)
  for (int i = 0; i < row_dim; i++)
  {
    for (j = 0; j < col_dim; j++)
    {
      V[i][j] *= alpha;
    }
  }
};

void Subsystem_Sz::MP_mm_dcopy(const int row_dim, const int col_dim, double **V0, double **V1)
{
  int j;
#pragma omp parallel for private(j)
  for (int i = 0; i < row_dim; i++)
  {
    for (j = 0; j < col_dim; j++)
    {
      V1[i][j] = V0[i][j];
    }
  }
};

void Subsystem_Sz::MP_schedule_mm_dcopy(const int row_dim, const int col_dim, double **V0, double **V1)
{
  int j;
#pragma omp parallel for private(j) schedule(runtime)
  for (int i = 0; i < row_dim; i++)
  {
    for (j = 0; j < col_dim; j++)
    {
      V1[i][j] = V0[i][j];
    }
  }
};

void Subsystem_Sz::MP_mm_daxpy(double alpha, const int row_dim, const int col_dim, double **V0, double **V1)
{
  int j;
#pragma omp parallel for private(j)
  for (int i = 0; i < row_dim; i++)
  {
    for (j = 0; j < col_dim; j++)
    {
      V1[i][j] += alpha * V0[i][j];
    }
  }
};

void Subsystem_Sz::MP_schedule_mm_daxpy(double alpha, const int row_dim, const int col_dim, double **V0, double **V1)
{
  int j;
#pragma omp parallel for private(j) schedule(runtime)
  for (int i = 0; i < row_dim; i++)
  {
    for (j = 0; j < col_dim; j++)
    {
      V1[i][j] += alpha * V0[i][j];
    }
  }
};

double Subsystem_Sz::MP_mm_dnrm2(const int row_dim, const int col_dim, double **V)
{
  int j;
  double val = 0.;
#pragma omp parallel for private(j) reduction(+ : val)
  for (int i = 0; i < row_dim; i++)
  {
    for (j = 0; j < col_dim; j++)
    {
      val += V[i][j] * V[i][j];
    }
  }
  return sqrt(val);
};

double Subsystem_Sz::MP_schedule_mm_dnrm2(const int row_dim, const int col_dim, double **V)
{
  int j;
  double val = 0.;
#pragma omp parallel for private(j) reduction(+ : val) schedule(runtime)
  for (int i = 0; i < row_dim; i++)
  {
    for (j = 0; j < col_dim; j++)
    {
      val += V[i][j] * V[i][j];
    }
  }
  return sqrt(val);
};

void Subsystem_Sz::MP_mm_sdz(const int row_dim, const int col_dim, double **V)
{
  double a = 1. / MP_mm_dnrm2(row_dim, col_dim, V);
  MP_mm_dscal(a, row_dim, col_dim, V);
};

void Subsystem_Sz::MP_schedule_mm_sdz(const int row_dim, const int col_dim, double **V)
{
  double a = 1. / MP_schedule_mm_dnrm2(row_dim, col_dim, V);
  MP_schedule_mm_dscal(a, row_dim, col_dim, V);
};

void Subsystem_Sz::calc_alpha_evenstep(const int ls, double *alpha)
{
  double tmp = 0.;
  for (int No = 0; No < pair_num; No++)
  {
    tmp += mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V1);
  }
  alpha[ls] = tmp;
}

void Subsystem_Sz::MP_calc_alpha_evenstep(const int ls, double *alpha)
{
  double tmp = 0.;
  for (int No = 0; No < pair_num; No++)
  {
    tmp += MP_mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V1);
  }
  alpha[ls] = tmp;
}

void Subsystem_Sz::MP_schedule_calc_alpha_evenstep(const int ls, double *alpha)
{
  double tmp = 0.;
  for (int No = 0; No < pair_num; No++)
  {
    tmp += MP_schedule_mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V1);
  }
  alpha[ls] = tmp;
}

void Subsystem_Sz::calc_alpha_oddstep(const int ls, double *alpha)
{
  double tmp = 0.;
  for (int No = 0; No < pair_num; No++)
  {
    tmp += mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].V0);
  }
  alpha[ls] = tmp;
};

void Subsystem_Sz::MP_calc_alpha_oddstep(const int ls, double *alpha)
{
  double tmp = 0.;
  for (int No = 0; No < pair_num; No++)
  {
    tmp += MP_mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].V0);
  }
  alpha[ls] = tmp;
};

void Subsystem_Sz::MP_schedule_calc_alpha_oddstep(const int ls, double *alpha)
{
  double tmp = 0.;
  for (int No = 0; No < pair_num; No++)
  {
    tmp += MP_schedule_mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].V0);
  }
  alpha[ls] = tmp;
};

void Subsystem_Sz::MP_calc_beta_evenstep(const int tri_mat_dim, const int ls, double *alpha, double *beta)
{
  double tmp = 0.;
  if (ls != tri_mat_dim - 1)
  {
    // #pragma omp parallel for reduction(+ \
//                                    : tmp)
    for (int No = 0; No < pair_num; No++)
    {
      MP_mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V1);
      tmp += MP_mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].V1);
    }
    beta[ls] = sqrt(tmp);
  }
  else
  {
    // #pragma omp parallel for
    for (int No = 0; No < pair_num; No++)
    {
      MP_mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V1);
    }
  }
};

void Subsystem_Sz::MP_schedule_calc_beta_evenstep(const int tri_mat_dim, const int ls, double *alpha, double *beta)
{
  double tmp = 0.;
  if (ls != tri_mat_dim - 1)
  {
    for (int No = 0; No < pair_num; No++)
    {
      MP_schedule_mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V1);
      tmp += MP_schedule_mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].V1);
    }
    beta[ls] = sqrt(tmp);
  }
  else
  {
    // #pragma omp parallel for
    for (int No = 0; No < pair_num; No++)
    {
      MP_schedule_mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V1);
    }
  }
};

void Subsystem_Sz::MP_calc_beta_oddstep(const int tri_mat_dim, const int ls, double *alpha, double *beta)
{
  double tmp = 0;
  // 並列化を意識してloop中のif文を除外したコード
  if (ls != tri_mat_dim - 1)
  {
    for (int No = 0; No < pair_num; No++)
    {
      MP_mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].V0);
      tmp += MP_mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V0);
    }
    beta[ls] = sqrt(tmp);
  }
  else
  {
    for (int No = 0; No < pair_num; No++)
    {
      MP_mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].V0);
    }
  }
};

void Subsystem_Sz::MP_schedule_calc_beta_oddstep(const int tri_mat_dim, const int ls, double *alpha, double *beta)
{
  double tmp = 0;
  // 並列化を意識してloop中のif文を除外したコード
  if (ls != tri_mat_dim - 1)
  {
    for (int No = 0; No < pair_num; No++)
    {
      MP_schedule_mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].V0);
      tmp += MP_schedule_mm_ddot(tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V0, tot_Sz[No].V0);
    }
    beta[ls] = sqrt(tmp);
  }
  else
  {
    for (int No = 0; No < pair_num; No++)
    {
      MP_schedule_mm_daxpy(-alpha[ls], tot_Sz[No].bm_A_size, tot_Sz[No].bm_B_size, tot_Sz[No].V1, tot_Sz[No].V0);
    }
  }
};

/*磁化曲線のplotを行う*/
void plot_MHcurve(int sys_num, int sys_site_A, int sys_site_B, int max_up_spin, int min_up_spin, double J_red, double J_green, double J_blue, vector<string> &file, string GNUPLOT_DATA_DIR)
{
  vector<double> M;     // 磁化の情報を格納するための配列
  vector<double> e_min; // 各部分空間における基底状態のエネルギー固有値を代入するための配列

  vector<double> plot_h; // 交点のh座標を格納する(横軸)
  vector<double> plot_m; // 交点のM座標を格納する(縦軸)

  double abs_of_dif_Jr_Jg = abs(J_red - J_green);  //|J_red - J_green|の値
  double abs_of_dif_Jr_Jb = abs(J_red - J_blue);   //|J_red - J_blue|の値
  double abs_of_dif_Jg_Jb = abs(J_green - J_blue); //|J_green - J_blue|の値
  vector<double> plateau_width;                    // h_{i+1} - h_iの値を格納する

  for (int up = min_up_spin; up < max_up_spin + 1; up++)
  {
    Subsystem_Sz H(sys_num, sys_site_A, sys_site_B, file, up); // Subsystem_Szオブジェクトのupスピン本数の再設定

    // upスピンの本数が上記の場合における部分空間とHamiltonian行列の用意
    H.sub_space_check();
    H.set_system_info();
    H.sub_hamiltonian();

    H.sub_lanczos(500); // lanczos法による固有値計算

    // 基底状態での固有エネルギーと磁化の値の格納
    M.push_back(H.mag);
    e_min.push_back(H.eigen_value);
  }

  // 磁場の規格化を行う
  auto max_index = distance(M.begin(), max_element(M.begin(), M.end())); // min_indexの方はsize_t
  for (int i = 0; i < M.size(); i++)
  {
    M[i] /= M[max_index];
  }

  plot_h.push_back(0);
  plot_m.push_back(M[0]);

  // E-hグラフにおける交点を求める処理を以下で行う
  for (int i = 0; i < M.size() - 1; i++)
  {
    vector<double> mi_intersec(M.size(), 20.0); // M_iとM_i+1...との交点のx座標を格納した配列

    for (int k = i + 1; k < M.size(); k++) // M_iとM_i+1の交点の座標を計算して記録
    {
      mi_intersec[k - (i + 1)] = double((e_min[k] - e_min[i]) / (k - i));
    }

    // 列挙した交点のうち、最もx座標の値が小さいものを選択して記録する
    auto min_index = distance(mi_intersec.begin(), min_element(mi_intersec.begin(), mi_intersec.end())); // min_indexの方はsize_t
    // M[i]とM[i+1]のhの交点座標 > M[i]とM[i+2]のhの交点座標の場合にM[i+1]とM[i+2],M[i+3]...のhの交点を調べるのはskipする
    int skip_itr = int(min_index) + 1;

    // skipされるぶんについても値を出力する
    for (int l = 0; l < skip_itr - 1; l++)
    {
      plot_h.push_back(mi_intersec[min_index]);
      plot_m.push_back(M[i + l + 1]);
    }

    plot_h.push_back(mi_intersec[min_index]);
    plot_m.push_back(M[i + 1]);
    // 上記に該当する磁化のindexを取得する
    i += int(min_index); // M_iとM_{i+α}の交点が最小だった場合には、M_{i+α-1}に着目したloopはskipする。
  }

  //====================== プラトー幅を調べる ===============================
  double width;
  for (int i = 0; i < M.size() - 2; i++)
  {
    width = plot_h[i + 1] - plot_h[i];
    plateau_width.push_back(width);
  }

  plateau_width.push_back(0);

  // 交点の情報の表示
  cout
      << "================================================\n";
  cout << setw(5) << "h"
       << "  " << setw(16) << "M"
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
    cout << left << setw(5) << plot_h[i]
         << " " << setw(15) << plot_m[i]
         << " " << J_red
         << " " << J_green
         << " " << J_blue
         << " " << plateau_width[i] << endl;
  }
  cout << "================================================\n";

  // 交点の情報をファイルへ書き出す

  int site_num = sys_site_A + sys_site_B;
  ofstream plateau_data(GNUPLOT_DATA_DIR);
  for (int i = 0; i < plot_h.size(); i++)
  {
    plateau_data << left << plot_h[i]
                 << " " << setw(15) << plot_m[i]
                 << " " << J_red
                 << " " << J_green
                 << " " << J_blue
                 << " " << plateau_width[i] << endl;
  }

  plateau_data.close();
  // FILE *gnuplot = popen("gnuplot -persist", "w");
  // fprintf(gnuplot, "set key left top\n");
  // fprintf(gnuplot, "set xrange[0:%f]\n", plot_h[plot_h.size() - 1] + 0.5);
  // fprintf(gnuplot, "plot \"/home/kazu/Documents/numcal/My_Research_Code/class_Subsystem_Sz/plateau/kagome_18site_uniform/data.txt\" using 1:2 with steps title \"18site\"\n");
  // fprintf(gnuplot, "set title \"kagome lattice 18site uniform\"\n");
  // fprintf(gnuplot, "set term png\n");
  // fprintf(gnuplot, "set output \"/home/kazu/Documents/numcal/My_Research_Code/class_Subsystem_Sz/plateau/kagome_18site_uniform/kagome_18site_uniform_plateau.png\";replot\n");
  // fprintf(gnuplot, "set term qt\n");
  // fflush(gnuplot);
  // pclose(gnuplot);
}

/*磁化曲線のplotを行う*/
void MP_plot_MHcurve(int sys_num, int sys_site_A, int sys_site_B, int max_up_spin, int min_up_spin, double J_red, double J_green, double J_blue, vector<string> &file, string GNUPLOT_DATA_DIR)
{
  vector<double> M;     // 磁化の情報を格納するための配列
  vector<double> e_min; // 各部分空間における基底状態のエネルギー固有値を代入するための配列

  vector<double> plot_h; // 交点のh座標を格納する(横軸)
  vector<double> plot_m; // 交点のM座標を格納する(縦軸)

  double abs_of_dif_Jr_Jg = abs(J_red - J_green);  //|J_red - J_green|の値
  double abs_of_dif_Jr_Jb = abs(J_red - J_blue);   //|J_red - J_blue|の値
  double abs_of_dif_Jg_Jb = abs(J_green - J_blue); //|J_green - J_blue|の値
  vector<double> plateau_width;                    // h_{i+1} - h_iの値を格納する

  for (int up = min_up_spin; up < max_up_spin + 1; up++)
  {
    Subsystem_Sz H(sys_num, sys_site_A, sys_site_B, file, up); // Subsystem_Szオブジェクトのupスピン本数の再設定

    // upスピンの本数が上記の場合における部分空間とHamiltonian行列の用意
    H.sub_space_check();
    H.set_system_info();
    H.sub_hamiltonian();

    H.MP_sub_lanczos(1000); // lanczos法による固有値計算

    // 基底状態での固有エネルギーと磁化の値の格納
    M.push_back(H.mag);
    e_min.push_back(H.eigen_value);
  }

  // 磁場の規格化を行う
  auto max_index = distance(M.begin(), max_element(M.begin(), M.end())); // min_indexの方はsize_t
  for (int i = 0; i < M.size(); i++)
  {
    M[i] /= M[max_index];
  }

  plot_h.push_back(0);
  plot_m.push_back(M[0]);

  // E-hグラフにおける交点を求める処理を以下で行う
  for (int i = 0; i < M.size() - 1; i++)
  {
    vector<double> mi_intersec(M.size(), 20.0); // M_iとM_i+1,M_i+2...との交点のx座標を格納した配列

    for (int k = i + 1; k < M.size(); k++) // M_iとM_i+1の交点の座標を計算して記録 k < M.size() - 1 -> k < M.size()
    {
      mi_intersec[k - (i + 1)] = double((e_min[k] - e_min[i]) / (k - i));
    }

    // 列挙した交点のうち、最もx座標の値が小さいものを選択して記録する
    auto min_index = distance(mi_intersec.begin(), min_element(mi_intersec.begin(), mi_intersec.end())); // min_indexの方はsize_t
    int skip_itr = int(min_index) + 1;

    // skipされるぶんについても値を出力する
    for (int l = 0; l < skip_itr - 1; l++)
    {
      plot_h.push_back(mi_intersec[min_index]);
      plot_m.push_back(M[i + l + 1]);
    }

    plot_h.push_back(mi_intersec[min_index]);
    plot_m.push_back(M[i + skip_itr]);
    // 上記に該当する磁化のindexを取得する
    i += int(min_index); // M_iとM_{i+α}の交点が最小だった場合には、M_{i+α-1}に着目したloopはskipする。
  }

  //====================== プラトー幅を調べる ===============================
  double width;
  for (int i = 0; i < M.size() - 2; i++)
  {
    width = plot_h[i + 1] - plot_h[i];
    plateau_width.push_back(width);
  }

  plateau_width.push_back(0);

  // 交点の情報の表示
  cout
      << "================================================\n";
  cout << setw(5) << "h"
       << "  " << setw(16) << "M"
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
    cout << left << setw(5) << plot_h[i]
         << " " << setw(15) << plot_m[i]
         << " " << J_red
         << " " << J_green
         << " " << J_blue
         << " " << plateau_width[i] << endl;
  }
  cout << "================================================\n";

  // 交点の情報をファイルへ書き出す

  int site_num = sys_site_A + sys_site_B;
  ofstream plateau_data(GNUPLOT_DATA_DIR);
  for (int i = 0; i < plot_h.size(); i++)
  {
    plateau_data << left << plot_h[i]
                 << " " << setw(15) << plot_m[i]
                 << " " << J_red
                 << " " << J_green
                 << " " << J_blue
                 << " " << plateau_width[i] << endl;
  }

  plateau_data.close();

  // FILE *gnuplot = popen("gnuplot -persist", "w");
  // fprintf(gnuplot, "set key left top\n");
  // fprintf(gnuplot, "set xrange[0:%f]\n", plot_h[plot_h.size() - 1] + 0.5);
  // fprintf(gnuplot, "plot \"%s\" using 1:2 with steps title \"%s\"\n", GNUPLOT_DATA_DIR.c_str(), GNUPLOT_LEGEND.c_str());
  // fprintf(gnuplot, "set title \"%s\"\n", GNUPLOT_TITLE.c_str());
  // fprintf(gnuplot, "set term png\n");
  // fprintf(gnuplot, "set output \"%s\";replot\n", GNUPLOT_GRAPH_DIR.c_str());
  // fprintf(gnuplot, "set term qt\n");
  // fflush(gnuplot);
  // pclose(gnuplot);
}

/*磁化曲線のplotを行う*/
void MP_schedule_plot_MHcurve(int sys_num, int sys_site_A, int sys_site_B, int max_up_spin, int min_up_spin, double J_red, double J_green, double J_blue, vector<string> &file, string GNUPLOT_DATA_DIR)
{
  vector<double> M;     // 磁化の情報を格納するための配列
  vector<double> e_min; // 各部分空間における基底状態のエネルギー固有値を代入するための配列

  vector<double> plot_h; // 交点のh座標を格納する(横軸)
  vector<double> plot_m; // 交点のM座標を格納する(縦軸)

  double abs_of_dif_Jr_Jg = abs(J_red - J_green);  //|J_red - J_green|の値
  double abs_of_dif_Jr_Jb = abs(J_red - J_blue);   //|J_red - J_blue|の値
  double abs_of_dif_Jg_Jb = abs(J_green - J_blue); //|J_green - J_blue|の値
  vector<double> plateau_width;                    // h_{i+1} - h_iの値を格納する

  for (int up = min_up_spin; up < max_up_spin + 1; up++)
  {
    double start_subsystem = omp_get_wtime();
    Subsystem_Sz H(sys_num, sys_site_A, sys_site_B, file, up); // Subsystem_Szオブジェクトのupスピン本数の再設定

    // upスピンの本数が上記の場合における部分空間とHamiltonian行列の用意
    H.sub_space_check();
    H.set_system_info();
    H.sub_hamiltonian();

    H.MP_schedule_sub_lanczos(1000); // lanczos法による固有値計算
    double end_subsystem = omp_get_wtime();

    cout << "run_time = " << end_subsystem - start_subsystem << endl;
    cout << H << endl;
    // 基底状態での固有エネルギーと磁化の値の格納
    M.push_back(H.mag);
    e_min.push_back(H.eigen_value);
  }

  // 磁場の規格化を行う
  auto max_index = distance(M.begin(), max_element(M.begin(), M.end())); // min_indexの方はsize_t
  for (int i = 0; i < M.size(); i++)
  {
    M[i] /= M[max_index];
  }

  // M0についてはここで格納する
  plot_h.push_back(0);
  plot_m.push_back(M[0]);

  // E-hグラフにおける交点を求める処理を以下で行う
  // M[i]とM[i+1]の交点を求める
  for (int i = 0; i < M.size() - 1; i++)
  {
    vector<double> mi_intersec(M.size(), 20.0); // M_iとM_i+1,M_i+2...との交点のx座標を格納した配列

    for (int k = i + 1; k < M.size(); k++) // M_iとM_i+1の交点の座標を計算して記録
    {
      mi_intersec[k - (i + 1)] = double((e_min[k] - e_min[i]) / (k - i));
    }

    // 列挙した交点のうち、最もx座標の値が小さいものを選択して記録する
    auto min_index = distance(mi_intersec.begin(), min_element(mi_intersec.begin(), mi_intersec.end())); // min_indexの方はsize_t
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

  int site_num = sys_site_A + sys_site_B;
  ofstream plateau_data(GNUPLOT_DATA_DIR);
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

  // FILE *gnuplot = popen("gnuplot -persist", "w");
  // fprintf(gnuplot, "set key left top\n");
  // fprintf(gnuplot, "set xrange[0:%f]\n", plot_h[plot_h.size() - 1] + 0.5);
  // fprintf(gnuplot, "plot \"%s\" using 1:2 with steps title \"%s\"\n", GNUPLOT_DATA_DIR.c_str(), GNUPLOT_LEGEND.c_str());
  // fprintf(gnuplot, "set title \"%s\"\n", GNUPLOT_TITLE.c_str());
  // fprintf(gnuplot, "set term png\n");
  // fprintf(gnuplot, "set output \"%s\";replot\n", GNUPLOT_GRAPH_DIR.c_str());
  // fprintf(gnuplot, "set term qt\n");
  // fflush(gnuplot);
  // pclose(gnuplot);
}

// spin-spin相関の計算<Ψ|S_i^zS_j^z|Ψ>
void Subsystem_Sz::clac_szz_rel(const int site_num, std::string dir_output)
{
  double rel_ij; // サイトi,jの相関係数
  bool is_up_spin_i, is_up_spin_j;
  int sign; // S_i^zS_j^z|Ψ>の符号を代入する sign = ±1
  double evec_val;
  ofstream ofs(dir_output);

  // siteについてのloop(はじめにどのサイト同士の相関係数を計算するかを指定する)
  //[1]i,jがともに部分系Aに属するサイトの場合
  ofs << setw(3) << "i∈A"
      << ", "
      << "j∈A"
      << ","
      << setw(18) << "rel_ij"
      << ","
      << "rel_ij x 5.0"
      << "\n ";
  ofs
      << "---------------------------------------------\n";
  for (int i = 0; i < tot_site_A; i++)
  {
    for (int j = i; j < tot_site_A; j++)
    {
      rel_ij = 0.;
      for (int No = 0; No < pair_num; No++)
      {

        // 状態についてのloop
        for (int n = 0; n < tot_Sz[No].bm_A_size; n++)
        {
          for (int m = 0; m < tot_Sz[No].bm_B_size; m++)
          {
            // 状態の用意
            // スピン演算子は部分系AのサイトについてのものなのでbitはAのものだけを用意すれば良い
            int state_num_of_A = tot_Sz[No].bm_A[n];
            boost::dynamic_bitset<> ket_A(tot_site_A, state_num_of_A);

            // i,j番目のスピンの向きを確認
            is_up_spin_i = ket_A.test(i);
            is_up_spin_j = ket_A.test(j);

            if (is_up_spin_i == is_up_spin_j)
              sign = 1;
            else
              sign = -1;
            evec_val = tot_Sz[No].Eig.eigen_mat[n][m];
            rel_ij += sign * 0.25 * evec_val * evec_val;
          }
        }
      }
      ofs << setw(3) << i << ", " << j << ", " << setw(18) << setprecision(15) << rel_ij << " , " << setprecision(3) << rel_ij * 5.0 << endl;
    }
  }
  ofs << "---------------------------------------------\n";
  //[2]i,jがともに部分系Bに属するサイトの場合
  ofs << setw(3) << "i∈B"
      << ", "
      << "j∈B"
      << ","
      << setw(18) << "rel_ij"
      << ","
      << "rel_ij x 5.0"
      << "\n";
  ofs << "---------------------------------------------\n";
  for (int i = 0; i < tot_site_B; i++)
  {
    for (int j = i; j < tot_site_B; j++)
    {
      rel_ij = 0.;
      for (int No = 0; No < pair_num; No++)
      {
        // 状態についてのloop
        for (int n = 0; n < tot_Sz[No].bm_A_size; n++)
        {
          for (int m = 0; m < tot_Sz[No].bm_B_size; m++)
          {
            // 状態の用意
            // スピン演算子は部分系BのサイトについてのものなのでbitはBのものだけを用意すれば良い
            int state_num_of_B = tot_Sz[No].bm_B[m];
            boost::dynamic_bitset<> ket_B(tot_site_B, state_num_of_B);

            // i,j番目のスピンの向きを確認
            is_up_spin_i = ket_B.test(i);
            is_up_spin_j = ket_B.test(j);

            if (is_up_spin_i == is_up_spin_j)
              sign = 1;
            else
              sign = -1;

            evec_val = tot_Sz[No].Eig.eigen_mat[n][m];
            rel_ij += sign * 0.25 * evec_val * evec_val;
          }
        }
      }
      ofs << setw(3) << i << ", " << j << "," << setw(18) << setprecision(15) << rel_ij << " , " << setprecision(3) << rel_ij * 5.0 << endl;
    }
  }
  //[3]2サイトがそれぞれ部分系A,Bに属する場合
  ofs << "---------------------------------------------\n";
  ofs << setw(3) << "i∈A"
      << ", "
      << "j∈B"
      << ","
      << setw(18) << "rel_ij\n";
  ofs << "---------------------------------------------\n";
  for (int i = 0; i < tot_site_A; i++)
  {
    for (int j = 0; j < tot_site_B; j++)
    {
      rel_ij = 0.;
      // 部分空間についてのloop
      for (int No = 0; No < pair_num; No++)
      {
        // 状態についてのloop
        for (int n = 0; n < tot_Sz[No].bm_A_size; n++)
        {
          for (int m = 0; m < tot_Sz[No].bm_B_size; m++)
          {
            // 状態の用意
            int state_num_of_A = tot_Sz[No].bm_A[n];
            int state_num_of_B = tot_Sz[No].bm_B[m];
            boost::dynamic_bitset<> ket_A(tot_site_A, state_num_of_A);
            boost::dynamic_bitset<> ket_B(tot_site_B, state_num_of_B);

            // i,j番目のスピンの向きを確認
            is_up_spin_i = ket_A.test(i);
            is_up_spin_j = ket_B.test(j);

            if (is_up_spin_i == is_up_spin_j)
              sign = 1;
            else
              sign = -1;
            evec_val = tot_Sz[No].Eig.eigen_mat[n][m];
            rel_ij += sign * 0.25 * evec_val * evec_val;
          }
        }
      }
      ofs << setw(3) << i << ", " << j << ", " << setw(18) << setprecision(15) << rel_ij << " , " << setprecision(3) << rel_ij * 5.0 << endl;
    }
  }
  ofs.close();
};

void Subsystem_Sz::calc_sxx_rel(const int site_num, string dir_output)
{
  double rel_ij; // サイトi,jの相関係数
  bool is_up_spin_i, is_up_spin_j;
  double evec_val;
  ofstream ofs(dir_output);

  //[1]i,jがともに部分系Aに属するサイトの場合
  ofs << setw(3) << "i∈A"
      << ", "
      << "j∈A"
      << ","
      << setw(18) << "rel_ij"
      << ","
      << "rel_ij x 5.0"
      << "\n ";
  ofs
      << "---------------------------------------------\n";
  for (int i = 0; i < tot_site_A; i++)
  {
    for (int j = 0; j < tot_site_A; j++) // j=0 startにすることでS_i^+Sj-とS_i^-とS_j^+を一度に計算することができる
    {
      rel_ij = 0.;
      for (int No = 0; No < pair_num; No++)
      {
        // 状態についてのloop(|スピン状態> = |n>|m>)
        for (int n = 0; n < tot_Sz[No].bm_A_size; n++)
        {
          for (int m = 0; m < tot_Sz[No].bm_B_size; m++)
          {
            // 状態の用意
            // スピン演算子は部分系AのサイトについてのものなのでbitはAのものだけを用意すれば良い
            int state_num_of_A = tot_Sz[No].bm_A[n];
            boost::dynamic_bitset<> ket_A(tot_site_A, state_num_of_A);

            // i,j番目のスピンの向きを確認
            is_up_spin_i = ket_A.test(i); // bitが1(up)ならtrue、そうでないならfalse
            is_up_spin_j = ket_A.test(j);

            // i番目のスピンがdownかつj番目のスピンがupの場合
            if (is_up_spin_i == false && is_up_spin_j == true)
            {
              boost::dynamic_bitset<> ket_A1(tot_site_A, state_num_of_A);
              ket_A1.flip(i);
              ket_A1.flip(j);

              int n_trans = (int)(ket_A1.to_ulong()); //|n>の遷移先を調べる
              int bm_ctr = tot_Sz[No].gbm_A[n_trans]; // |n>の遷移先が部分空間で何番目のスピン状態かを調べる
              rel_ij += 0.25 * tot_Sz[No].Eig.eigen_mat[bm_ctr][m] * tot_Sz[No].Eig.eigen_mat[n][m];
            }

            // i番目のスピンがupかつj番目のスピンがdownの場合
            if (is_up_spin_i == true && is_up_spin_j == false)
            {
              boost::dynamic_bitset<> ket_A1(tot_site_A, state_num_of_A);
              ket_A1.flip(i);
              ket_A1.flip(j);

              int n_trans = (int)(ket_A1.to_ulong()); //|n>の遷移先を調べる
              int bm_ctr = tot_Sz[No].gbm_A[n_trans]; // |n>の遷移先が部分空間で何番目のスピン状態かを調べる
              rel_ij += 0.25 * tot_Sz[No].Eig.eigen_mat[bm_ctr][m] * tot_Sz[No].Eig.eigen_mat[n][m];
            }
          }
        }
      }
      ofs << setw(3) << i << ", " << j << ", " << setw(18) << setprecision(15) << rel_ij << endl;
    }
  }
  ofs << "---------------------------------------------\n";

  //[2]i,jがともに部分系Bに属するサイトの場合
  ofs << setw(3) << "i∈B"
      << ", "
      << "j∈B"
      << ","
      << setw(18) << "rel_ij"
      << ","
      << "rel_ij x 5.0"
      << "\n";
  ofs << "---------------------------------------------\n";
  for (int i = 0; i < tot_site_B; i++)
  {
    for (int j = 0; j < tot_site_B; j++)
    {
      rel_ij = 0.;
      for (int No = 0; No < pair_num; No++)
      {
        // 状態についてのloop
        for (int n = 0; n < tot_Sz[No].bm_A_size; n++)
        {
          for (int m = 0; m < tot_Sz[No].bm_B_size; m++)
          {
            // 状態の用意
            // スピン演算子は部分系BのサイトについてのものなのでbitはBのものだけを用意すれば良い
            int state_num_of_B = tot_Sz[No].bm_B[m];
            boost::dynamic_bitset<> ket_B(tot_site_B, state_num_of_B);

            // i,j番目のスピンの向きを確認
            is_up_spin_i = ket_B.test(i);
            is_up_spin_j = ket_B.test(j);

            // i番目のスピンがdownかつj番目のスピンがupの場合
            if (is_up_spin_i == false && is_up_spin_j == true)
            {
              boost::dynamic_bitset<> ket_B1(tot_site_B, state_num_of_B);
              ket_B1.flip(i);
              ket_B1.flip(j);

              int m_trans = (int)(ket_B1.to_ulong()); //|m>の遷移先を調べる
              int bm_ctr = tot_Sz[No].gbm_B[m_trans]; //|m>の遷移先が部分空間で何番目のスピン状態かを調べる
              rel_ij += 0.25 * tot_Sz[No].Eig.eigen_mat[n][bm_ctr] * tot_Sz[No].Eig.eigen_mat[n][m];
            }

            // i番目のスピンがupかつj番目のスピンがdownの場合
            if (is_up_spin_i == true && is_up_spin_j == false)
            {
              boost::dynamic_bitset<> ket_B1(tot_site_B, state_num_of_B);
              ket_B1.flip(i);
              ket_B1.flip(j);

              int m_trans = (int)(ket_B1.to_ulong()); //|m>の遷移先を調べる
              int bm_ctr = tot_Sz[No].gbm_B[m_trans]; //|m>の遷移先が部分空間で何番目のスピン状態かを調べる
              rel_ij += 0.25 * tot_Sz[No].Eig.eigen_mat[n][bm_ctr] * tot_Sz[No].Eig.eigen_mat[n][m];
            }
          }
        }
      }
      ofs << setw(3) << i << ", " << j << "," << setw(18) << setprecision(15) << rel_ij << " , " << setprecision(3) << rel_ij * 5.0 << endl;
    }
  }

  //[3]2サイトがそれぞれ部分系A,Bに属する場合
  ofs << "---------------------------------------------\n";
  ofs << setw(3) << "i∈A"
      << ", "
      << "j∈B"
      << ","
      << setw(18) << "rel_ij\n";
  ofs << "---------------------------------------------\n";
  for (int i = 0; i < tot_site_A; i++)
  {
    for (int j = 0; j < tot_site_B; j++)
    {
      rel_ij = 0.;
      // 部分空間についてのloop
      for (int No = 0; No < pair_num; No++)
      {
        // 状態についてのloop
        for (int n = 0; n < tot_Sz[No].bm_A_size; n++)
        {
          for (int m = 0; m < tot_Sz[No].bm_B_size; m++)
          {
            // 状態の用意
            int state_num_of_A = tot_Sz[No].bm_A[n];
            int state_num_of_B = tot_Sz[No].bm_B[m];
            boost::dynamic_bitset<> ket_A(tot_site_A, state_num_of_A);
            boost::dynamic_bitset<> ket_B(tot_site_B, state_num_of_B);

            // i,j番目のスピンの向きを確認
            is_up_spin_i = ket_A.test(i);
            is_up_spin_j = ket_B.test(j);

            // i,j番目のスピンがdownの場合
            if (No > 0)
            {
              if (is_up_spin_i == false && is_up_spin_j == true)
              {
                boost::dynamic_bitset<> ket_A1(tot_site_A, state_num_of_A);
                boost::dynamic_bitset<> ket_B1(tot_site_B, state_num_of_B);
                ket_A1.flip(i);
                ket_B1.flip(j);

                int n_trans = (int)(ket_A1.to_ulong());
                int m_trans = (int)(ket_B1.to_ulong());

                int bn_ctr = tot_Sz[No - 1].gbm_A[n_trans];
                int bm_ctr = tot_Sz[No - 1].gbm_B[m_trans];

                rel_ij += 0.25 * tot_Sz[No - 1].Eig.eigen_mat[bn_ctr][bm_ctr] * tot_Sz[No].Eig.eigen_mat[n][m];
              }
            }

            if (No < pair_num - 1)
            {
              if (is_up_spin_i == true && is_up_spin_j == false)
              {

                boost::dynamic_bitset<> ket_A1(tot_site_A, state_num_of_A);
                boost::dynamic_bitset<> ket_B1(tot_site_B, state_num_of_B);
                ket_A1.flip(i);
                ket_B1.flip(j);

                int n_trans = (int)(ket_A1.to_ulong());
                int m_trans = (int)(ket_B1.to_ulong());

                int bn_ctr = tot_Sz[No + 1].gbm_A[n_trans];
                int bm_ctr = tot_Sz[No + 1].gbm_B[m_trans];

                rel_ij += 0.25 * tot_Sz[No + 1].Eig.eigen_mat[bn_ctr][bm_ctr] * tot_Sz[No].Eig.eigen_mat[n][m];
              }
            }
          }
        }
      }
      ofs << setw(3) << i << ", " << j << "," << setw(18) << setprecision(15) << rel_ij << " , " << setprecision(3) << rel_ij * 5.0 << endl;
    }
  }
  ofs.close();
};

/*------------------------------標準出力関係------------------------------*/
string Subsystem_Sz::to_string() const
{
  ostringstream s;
  s << "\n\n";
  s << "==================================\n";
  s << "[Info of Subsystem_Sz]\n";
  s << "----------------------------------\n";
  s << "- system_num       : " << system_num << endl;
  s << "- site num of A    : " << tot_site_A << endl;
  s << "- site num of B    : " << tot_site_B << endl;
  s << "- magnetization    : " << mag << endl;
  s << "- No of subspace   : " << pair_num << endl;
  s << "Total lanczos step : " << ls_count << endl;
  s << "- lanczos method   : ";
  if (ls_check == true)
    s << "success" << endl;
  else
    s << "failure" << endl;
  s << "Eigen value of Groundstate : " << setprecision(16) << eigen_value << endl;
  s << "==================================\n";

  return s.str();
}

ostream &operator<<(ostream &s, const Subsystem_Sz &H)
{
  return s << H.to_string();
}

int comb(int n, int r)
{
  vector<vector<int>> v(n + 1, vector<int>(n + 1, 0));
  for (int i = 0; i < v.size(); i++)
  {
    v[i][0] = 1;
    v[i][i] = 1;
  }

  for (int k = 2; k < v.size(); k++)
  {
    for (int j = 1; j < k; j++)
    {
      v[k][j] = (v[k - 1][j - 1] + v[k - 1][j]);
    }
  }
  return v[n][r];
}

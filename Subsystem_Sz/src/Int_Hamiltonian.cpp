#include "../include/Int_Hamiltonian.hpp"

#include <mkl.h>

#include <boost/dynamic_bitset.hpp>
#include <fstream>
#include <iomanip>
#include <random>

int Int_Hamiltonian::Int_id_counter;
int Int_Hamiltonian::tot_site_A;
int Int_Hamiltonian::tot_site_B;

using namespace std;

// コンストラクタ
Int_Hamiltonian::Int_Hamiltonian()
    : Int_id(Int_id_counter++),
      prow_ind_A(new int[1]),
      prow_ind_B(new int[1]),
      pcol_ind_A(new int[1]),
      pcol_ind_B(new int[1]),
      mrow_ind_A(new int[1]),
      mrow_ind_B(new int[1]),
      mcol_ind_A(new int[1]),
      mcol_ind_B(new int[1]),
      sz_A(new int[1]),
      sz_B(new int[1]),
      nnz_pA(0),
      nnz_pB(0),
      nnz_mA(0),
      nnz_mB(0),
      J(" ")
{
    // std::cout << "Int_Hamiltonian::constructed.\n";
}

void Int_Hamiltonian::int_rise_mmprod(const int nnz_pA, const int nnz_pB, int *prow_ind_A, int *pcol_ind_A, int *prow_ind_B, int *pcol_ind_B, double **V0, double **V1, double **V1_dic1)
{
    for (int i = 0; i < nnz_pA; i++)
    {
        for (int j = 0; j < nnz_pB; j++)
        {
            V1_dic1[prow_ind_A[i]][pcol_ind_B[j]] += 0.5 * J.val(0) * V0[pcol_ind_A[i]][prow_ind_B[j]];
        }
    }
}

/*--------------------OpenMP利用version----------------------------*/
void Int_Hamiltonian::MP_int_rise_mmprod(const int nnz_pA, const int nnz_pB, int *prow_ind_A, int *pcol_ind_A, int *prow_ind_B, int *pcol_ind_B, double **V0, double **V1, double **V1_dic1)
{
    int j;
    double bond = J.val(0);
#pragma omp parallel for private(j)
    for (int i = 0; i < nnz_pA; i++)
    {
        for (j = 0; j < nnz_pB; j++)
        {
            V1_dic1[prow_ind_A[i]][pcol_ind_B[j]] += 0.5 * bond * V0[pcol_ind_A[i]][prow_ind_B[j]];
        }
    }
}

void Int_Hamiltonian::int_dsmn_mmprod(const int nnz_mA, const int nnz_mB, int *mrow_ind_A, int *mcol_ind_A, int *mrow_ind_B, int *mcol_ind_B, double **V0, double **V1, double **V1_inc1)
{
    for (int i = 0; i < nnz_mA; i++)
    {
        for (int j = 0; j < nnz_mB; j++)
        {
            V1_inc1[mrow_ind_A[i]][mcol_ind_B[j]] += 0.5 * J.val(0) * V0[mcol_ind_A[i]][mrow_ind_B[j]];
        }
    }
}

void Int_Hamiltonian::MP_int_dsmn_mmprod(const int nnz_mA, const int nnz_mB, int *mrow_ind_A, int *mcol_ind_A, int *mrow_ind_B, int *mcol_ind_B, double **V0, double **V1, double **V1_inc1)
{
    int j;
    double bond = J.val(0);
#pragma omp parallel for private(j)
    for (int i = 0; i < nnz_mA; i++)
    {
        for (j = 0; j < nnz_mB; j++)
        {
            V1_inc1[mrow_ind_A[i]][mcol_ind_B[j]] += 0.5 * bond * V0[mcol_ind_A[i]][mrow_ind_B[j]];
        }
    }
}

void Int_Hamiltonian::int_zz_mmprod(const int dim_A, const int dim_B, int *sz_A, int *sz_B, double **V0, double **V1)
{
    for (int i = 0; i < dim_A; i++)
    {
        for (int j = 0; j < dim_B; j++)
        {
            V1[i][j] += 0.25 * J.val(0) * sz_A[i] * sz_B[j] * V0[i][j];
        }
    }
}

void Int_Hamiltonian::MP_int_zz_mmprod(const int dim_A, const int dim_B, int *sz_A, int *sz_B, double **V0, double **V1)
{
    int j;
    double bond = J.val(0);
#pragma omp parallel for private(j)
    for (int i = 0; i < dim_A; i++)
    {
        for (j = 0; j < dim_B; j++)
        {
            V1[i][j] += 0.25 * bond * sz_A[i] * sz_B[j] * V0[i][j];
        }
    }
}

void Int_Hamiltonian::MP_schedule_int_rise_mmprod(const int nnz_pA, const int nnz_pB, int *prow_ind_A, int *pcol_ind_A, int *prow_ind_B, int *pcol_ind_B, double **V0, double **V1, double **V1_dic1)
{
    int j;
    double bond = J.val(0);
#pragma omp parallel for private(j) schedule(runtime)
    for (int i = 0; i < nnz_pA; i++)
    {
        for (j = 0; j < nnz_pB; j++)
        {
            V1_dic1[prow_ind_A[i]][pcol_ind_B[j]] += 0.5 * bond * V0[pcol_ind_A[i]][prow_ind_B[j]];
        }
    }
}

void Int_Hamiltonian::MP_schedule_int_dsmn_mmprod(const int nnz_mA, const int nnz_mB, int *mrow_ind_A, int *mcol_ind_A, int *mrow_ind_B, int *mcol_ind_B, double **V0, double **V1, double **V1_inc1)
{
    int j;
    double bond = J.val(0);
#pragma omp parallel for private(j) schedule(runtime)
    for (int i = 0; i < nnz_mA; i++)
    {
        for (j = 0; j < nnz_mB; j++)
        {
            V1_inc1[mrow_ind_A[i]][mcol_ind_B[j]] += 0.5 * bond * V0[mcol_ind_A[i]][mrow_ind_B[j]];
        }
    }
}

void Int_Hamiltonian::MP_schedule_int_zz_mmprod(const int dim_A, const int dim_B, int *sz_A, int *sz_B, double **V0, double **V1)
{
    int j;
    double bond = J.val(0);
#pragma omp parallel for private(j) schedule(runtime)
    for (int i = 0; i < dim_A; i++)
    {
        for (j = 0; j < dim_B; j++)
        {
            V1[i][j] += 0.25 * bond * sz_A[i] * sz_B[j] * V0[i][j];
        }
    }
}

//[ToDo] to_stringの実装

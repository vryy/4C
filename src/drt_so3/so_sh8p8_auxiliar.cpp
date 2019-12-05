/*----------------------------------------------------------------------*/
/*! \file
\brief 8-node solid shell element
\level 2
\maintainer Christoph Meier
*/

/*----------------------------------------------------------------------*/
/* definitions */

/*----------------------------------------------------------------------*/
/* headers */
#include "so_sh8p8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/voigt_notation.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_globalproblem.H"
#include "Epetra_Time.h"
#include "Teuchos_TimeMonitor.hpp"

using VoigtMapping = ::UTILS::VOIGT::IndexMappings;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::AxialMetricsAtOrigin(const LINALG::Matrix<NUMNOD_, NUMDIM_>& xrefe,
    LINALG::Matrix<NUMDIM_, NUMDIM_>& jac0, LINALG::Matrix<NUMDIM_, 1>& metr0)
{
  // vector of df(origin)
  static double df0_vector[NUMDIM_ * NUMNOD_] = {-0.125, -0.125, -0.125, +0.125, -0.125, -0.125,
      +0.125, +0.125, -0.125, -0.125, +0.125, -0.125, -0.125, -0.125, +0.125, +0.125, -0.125,
      +0.125, +0.125, +0.125, +0.125, -0.125, +0.125, +0.125};

  // shape function derivatives, evaluated at origin (r=s=t=0.0)
  LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> df0(df0_vector, true);  // view

  // compute Jacobian, evaluated at element origin (r=s=t=0.0)
  jac0.Multiply(df0, xrefe);

  // line metrics of r-, s- and t-axis in reference space
  for (int i = 0; i < NUMDIM_; ++i)
  {
    double metr = 0.0;
    for (int j = 0; j < NUMDIM_; ++j)
    {
      metr += jac0(i, j) * jac0(i, j);
    }
    metr0(i) = std::sqrt(metr);
  }

  // Kette rechts
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::LocalMetrics(
    const LINALG::Matrix<NUMDIM_, NUMDIM_>& jac, LINALG::Matrix<NUMDIM_, NUMDIM_>& metr)
{
  // metrics of r-, s- and t-axis in reference space
  for (int i = 0; i < NUMDIM_; ++i)
  {
    for (int j = 0; j < NUMDIM_; ++j)
    {
      double metrij = 0.0;
      for (int k = 0; k < NUMDIM_; ++k)
      {
        metrij += jac(i, k) * jac(j, k);
      }
      metr(i, j) = metrij;
    }
  }

  // Kette links
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::AnsSetup2(
    const LINALG::Matrix<NUMNOD_, NUMDIM_>& xrefe,  // material element coords
    const LINALG::Matrix<NUMNOD_, NUMDIM_>& xcurr,  // current element coords
    std::vector<LINALG::Matrix<NUMDIM_, NUMNOD_>>**
        deriv_sp,                                            // derivs eval. at all sampling points
    std::vector<LINALG::Matrix<NUMDIM_, NUMDIM_>>& jac_sps,  // jac at all sampling points
    std::vector<LINALG::Matrix<NUMDIM_, NUMDIM_>>&
        jac_cur_sps,                                          // current jac at all sampling points
    LINALG::Matrix<NUMANS_ * NUMSP2ND_, NUMDISP_>& B_ans_loc  // modified B
)
{
  const int num_sp = NUMSP2ND_;

  // static matrix object of derivs at sampling points, kept in memory
  static std::vector<LINALG::Matrix<NUMDIM_, NUMNOD_>> df_sp(num_sp);
  static bool dfsp_eval;  // flag for re-evaluate everything

  if (dfsp_eval != 0)
  {                      // if true f,df already evaluated
    *deriv_sp = &df_sp;  // return adress of static object to target of pointer
  }
  else
  {
    /*====================================================================*/
    /* 8-node hexhedra Solid-Shell node topology
     * and location of sampling points A to H                             */
    /*--------------------------------------------------------------------*/
    /*                      t
     *                      |
     *             4========|================7
     *          // |        |   DU         //||
     *        //   |        |   |        //  ||
     *      // AU  |        |   D      //    ||
     *     5=======E=================6   CU  H
     *    ||   |   |        |   DL   ||   |  ||
     *    ||   A   |     BU o--------||-- C -------s
     *    ||   |   |      |/         ||   |  ||
     *    F    AL  0----- B ---------G --CL--3
     *    ||     //     / |          ||    //
     *    ||   //     /   BL         ||  //
     *    || //     r                ||//
     *     1=========================2
     *
     */
    /*====================================================================*/
    // gp sampling point value for linear fct
    const double gpl = 1.0 / sqrt(3.0);
    // (r,s,t) gp-locations of sampling points AL,AU,BL,BU,CL,CU,DL,DU
    // numsp = 8 here set explicitly to allow direct initializing
    //                           AL,  BL,  CL,  DL,  AU,  BU,  CU,  DU
    const double r[NUMSP2ND_] = {0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0};
    const double s[NUMSP2ND_] = {-1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0};
    const double t[NUMSP2ND_] = {-gpl, -gpl, -gpl, -gpl, gpl, gpl, gpl, gpl};

    // fill up df_sp w.r.t. rst directions (NUMDIM) at each sp
    for (int i = 0; i < num_sp; ++i)
    {
      // df wrt to r "+0" for each node(0..7) at each sp [i]
      df_sp[i](0, 0) = -(1.0 - s[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](0, 1) = (1.0 - s[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](0, 2) = (1.0 + s[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](0, 3) = -(1.0 + s[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](0, 4) = -(1.0 - s[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](0, 5) = (1.0 - s[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](0, 6) = (1.0 + s[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](0, 7) = -(1.0 + s[i]) * (1.0 + t[i]) * 0.125;

      // df wrt to s "+1" for each node(0..7) at each sp [i]
      df_sp[i](1, 0) = -(1.0 - r[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](1, 1) = -(1.0 + r[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](1, 2) = (1.0 + r[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](1, 3) = (1.0 - r[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](1, 4) = -(1.0 - r[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](1, 5) = -(1.0 + r[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](1, 6) = (1.0 + r[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](1, 7) = (1.0 - r[i]) * (1.0 + t[i]) * 0.125;

      // df wrt to t "+2" for each node(0..7) at each sp [i]
      df_sp[i](2, 0) = -(1.0 - r[i]) * (1.0 - s[i]) * 0.125;
      df_sp[i](2, 1) = -(1.0 + r[i]) * (1.0 - s[i]) * 0.125;
      df_sp[i](2, 2) = -(1.0 + r[i]) * (1.0 + s[i]) * 0.125;
      df_sp[i](2, 3) = -(1.0 - r[i]) * (1.0 + s[i]) * 0.125;
      df_sp[i](2, 4) = (1.0 - r[i]) * (1.0 - s[i]) * 0.125;
      df_sp[i](2, 5) = (1.0 + r[i]) * (1.0 - s[i]) * 0.125;
      df_sp[i](2, 6) = (1.0 + r[i]) * (1.0 + s[i]) * 0.125;
      df_sp[i](2, 7) = (1.0 - r[i]) * (1.0 + s[i]) * 0.125;
    }

    // return adresses of just evaluated matrices
    *deriv_sp = &df_sp;  // return adress of static object to target of pointer
    dfsp_eval = 1;       // now all arrays are filled statically
  }

  for (int sp = 0; sp < num_sp; ++sp)
  {
    // compute (REFERENCE) Jacobian matrix at all sampling points
    jac_sps[sp].Multiply(df_sp[sp], xrefe);
    // compute CURRENT Jacobian matrix at all sampling points
    jac_cur_sps[sp].Multiply(df_sp[sp], xcurr);
  }

  /*
  ** Compute modified B-operator in local(parametric) space,
  ** evaluated at all sampling points
  */
  // loop over each sampling point
  LINALG::Matrix<NUMDIM_, NUMDIM_> jac_cur;
  for (int sp = 0; sp < num_sp; ++sp)
  {
    /* compute the CURRENT Jacobian matrix at the sampling point:
    **         [ xcurr_,r  ycurr_,r  zcurr_,r ]
    **  Jcur = [ xcurr_,s  ycurr_,s  zcurr_,s ]
    **         [ xcurr_,t  ycurr_,t  zcurr_,t ]
    ** Used to transform the global displacements into parametric space
    */
    jac_cur.Multiply(df_sp[sp], xcurr);

    // fill up B-operator
    for (int inode = 0; inode < NUMNOD_; ++inode)
    {
      for (int dim = 0; dim < NUMDIM_; ++dim)
      {
        // modify B_loc_tt = N_t.X_t
        B_ans_loc(sp * num_ans + 0, inode * 3 + dim) = df_sp[sp](2, inode) * jac_cur(2, dim);
        // modify B_loc_st = N_s.X_t + N_t.X_s
        B_ans_loc(sp * num_ans + 1, inode * 3 + dim) =
            df_sp[sp](1, inode) * jac_cur(2, dim) + df_sp[sp](2, inode) * jac_cur(1, dim);
        // modify B_loc_rt = N_r.X_t + N_t.X_r
        B_ans_loc(sp * num_ans + 2, inode * 3 + dim) =
            df_sp[sp](0, inode) * jac_cur(2, dim) + df_sp[sp](2, inode) * jac_cur(0, dim);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::AnsSetup3(
    const LINALG::Matrix<NUMNOD_, NUMDIM_>& xrefe,  // material element coords
    const LINALG::Matrix<NUMNOD_, NUMDIM_>& xcurr,  // current element coords
    std::vector<LINALG::Matrix<NUMDIM_, NUMNOD_>>**
        deriv_sp,                                            // derivs eval. at all sampling points
    std::vector<LINALG::Matrix<NUMDIM_, NUMDIM_>>& jac_sps,  // jac at all sampling points
    std::vector<LINALG::Matrix<NUMDIM_, NUMDIM_>>&
        jac_cur_sps,                                          // current jac at all sampling points
    LINALG::Matrix<NUMANS_ * NUMSP3RD_, NUMDISP_>& B_ans_loc  // modified B
)
{
  const int num_sp = NUMSP3RD_;

  // static matrix object of derivs at sampling points, kept in memory
  static std::vector<LINALG::Matrix<NUMDIM_, NUMNOD_>> df_sp(num_sp);
  static bool dfsp_eval;  // flag for re-evaluate everything

  if (dfsp_eval != 0)
  {                      // if true f,df already evaluated
    *deriv_sp = &df_sp;  // return adress of static object to target of pointer
  }
  else
  {
    /*====================================================================*/
    /* 8-node hexhedra Solid-Shell node topology
     * and location of sampling points A to H                             */
    /*--------------------------------------------------------------------*/
    /*                      t
     *                      |
     *             4========|================7
     *          // |        |   DU         //||
     *        //   |        |   |        //  ||
     *      // AU  |        |   D      //    ||
     *     5=======E=================6   CU  H
     *    ||   |   |        |   DL   ||   |  ||
     *    ||   A   |     BU o--------||-- C -------s
     *    ||   |   |      |/         ||   |  ||
     *    F    AL  0----- B ---------G --CL--3
     *    ||     //     / |          ||    //
     *    ||   //     /   BL         ||  //
     *    || //     r                ||//
     *     1=========================2
     *
     */
    /*====================================================================*/
    // gp sampling point value for linear fct
    const double gpl = 1.0 / sqrt(3.0);
    // (r,s,t) gp-locations of sampling points
    // numsp = 8 here set explicitly to allow direct initializing
    //                             EL,  FL,  GL,  HL,  EU,  FU,  GU,  HU
    const double r[NUMSP3RD_] = {-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0};
    const double s[NUMSP3RD_] = {-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0};
    const double t[NUMSP3RD_] = {-gpl, -gpl, -gpl, -gpl, gpl, gpl, gpl, gpl};

    // fill up df_sp w.r.t. rst directions (NUMDIM) at each sp
    for (int i = 0; i < num_sp; ++i)
    {
      // df wrt to r "+0" for each node(0..7) at each sp [i]
      df_sp[i](0, 0) = -(1.0 - s[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](0, 1) = (1.0 - s[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](0, 2) = (1.0 + s[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](0, 3) = -(1.0 + s[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](0, 4) = -(1.0 - s[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](0, 5) = (1.0 - s[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](0, 6) = (1.0 + s[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](0, 7) = -(1.0 + s[i]) * (1.0 + t[i]) * 0.125;

      // df wrt to s "+1" for each node(0..7) at each sp [i]
      df_sp[i](1, 0) = -(1.0 - r[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](1, 1) = -(1.0 + r[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](1, 2) = (1.0 + r[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](1, 3) = (1.0 - r[i]) * (1.0 - t[i]) * 0.125;
      df_sp[i](1, 4) = -(1.0 - r[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](1, 5) = -(1.0 + r[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](1, 6) = (1.0 + r[i]) * (1.0 + t[i]) * 0.125;
      df_sp[i](1, 7) = (1.0 - r[i]) * (1.0 + t[i]) * 0.125;

      // df wrt to t "+2" for each node(0..7) at each sp [i]
      df_sp[i](2, 0) = -(1.0 - r[i]) * (1.0 - s[i]) * 0.125;
      df_sp[i](2, 1) = -(1.0 + r[i]) * (1.0 - s[i]) * 0.125;
      df_sp[i](2, 2) = -(1.0 + r[i]) * (1.0 + s[i]) * 0.125;
      df_sp[i](2, 3) = -(1.0 - r[i]) * (1.0 + s[i]) * 0.125;
      df_sp[i](2, 4) = (1.0 - r[i]) * (1.0 - s[i]) * 0.125;
      df_sp[i](2, 5) = (1.0 + r[i]) * (1.0 - s[i]) * 0.125;
      df_sp[i](2, 6) = (1.0 + r[i]) * (1.0 + s[i]) * 0.125;
      df_sp[i](2, 7) = (1.0 - r[i]) * (1.0 + s[i]) * 0.125;
    }

    // return adresses of just evaluated matrices
    *deriv_sp = &df_sp;  // return adress of static object to target of pointer
    dfsp_eval = 1;       // now all arrays are filled statically
  }

  for (int sp = 0; sp < num_sp; ++sp)
  {
    // compute (REFERENCE) Jacobian matrix at all sampling points
    jac_sps[sp].Multiply(df_sp[sp], xrefe);
    // compute CURRENT Jacobian matrix at all sampling points
    jac_cur_sps[sp].Multiply(df_sp[sp], xcurr);
  }

  /*
  ** Compute modified B-operator in local(parametric) space,
  ** evaluated at all sampling points
  */
  // loop over each sampling point
  LINALG::Matrix<NUMDIM_, NUMDIM_> jac_cur;
  for (int sp = 0; sp < num_sp; ++sp)
  {
    /* compute the CURRENT Jacobian matrix at the sampling point:
    **         [ xcurr_,r  ycurr_,r  zcurr_,r ]
    **  Jcur = [ xcurr_,s  ycurr_,s  zcurr_,s ]
    **         [ xcurr_,t  ycurr_,t  zcurr_,t ]
    ** Used to transform the global displacements into parametric space
    */
    jac_cur.Multiply(df_sp[sp], xcurr);

    // fill up B-operator
    for (int inode = 0; inode < NUMNOD_; ++inode)
    {
      for (int dim = 0; dim < NUMDIM_; ++dim)
      {
        // modify B_loc_tt = N_t.X_t
        B_ans_loc(sp * num_ans + 0, inode * 3 + dim) = df_sp[sp](2, inode) * jac_cur(2, dim);
        // modify B_loc_st = N_s.X_t + N_t.X_s
        B_ans_loc(sp * num_ans + 1, inode * 3 + dim) =
            df_sp[sp](1, inode) * jac_cur(2, dim) + df_sp[sp](2, inode) * jac_cur(1, dim);
        // modify B_loc_rt = N_r.X_t + N_t.X_r
        B_ans_loc(sp * num_ans + 2, inode * 3 + dim) =
            df_sp[sp](0, inode) * jac_cur(2, dim) + df_sp[sp](2, inode) * jac_cur(0, dim);
      }
    }
  }

  return;
}


void DRT::ELEMENTS::So_sh8p8::Matrix2TensorToVector9Voigt_Inconsistent(
    LINALG::Matrix<NUMDFGR_, 1>& fvct, const LINALG::Matrix<NUMDIM_, NUMDIM_>& fmat,
    const bool transpose)
{
  const int* voigt9row = NULL;
  const int* voigt9col = NULL;
  if (transpose)
  {
    voigt9row = &(VOIGT9COL_INCONSISTENT_[0]);
    voigt9col = &(VOIGT9ROW_INCONSISTENT_[0]);
  }
  else
  {
    voigt9row = &(VOIGT9ROW_INCONSISTENT_[0]);
    voigt9col = &(VOIGT9COL_INCONSISTENT_[0]);
  }

  for (int ij = 0; ij < NUMDFGR_; ++ij)
  {
    const int i = voigt9row[ij];
    const int j = voigt9col[ij];
    fvct(ij, 0) = fmat(i, j);  // F_ij
  }

  return;
}

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::InvVector9VoigtDiffByItself(
    LINALG::Matrix<NUMDFGR_, NUMDFGR_>& invfderf, const LINALG::Matrix<NUMDIM_, NUMDIM_>& invfmat,
    const bool transpose)
{
  // VERIFIED

  const int* voigt9row = NULL;
  const int* voigt9col = NULL;
  if (transpose)
  {
    voigt9row = &(VOIGT9COL_INCONSISTENT_[0]);
    voigt9col = &(VOIGT9ROW_INCONSISTENT_[0]);
  }
  else
  {
    voigt9row = &(VOIGT9ROW_INCONSISTENT_[0]);
    voigt9col = &(VOIGT9COL_INCONSISTENT_[0]);
  }

  for (int kl = 0; kl < NUMDFGR_; ++kl)
  {
    const int k = VOIGT9ROW_INCONSISTENT_[kl];
    const int l = VOIGT9COL_INCONSISTENT_[kl];
    for (int ij = 0; ij < NUMDFGR_; ++ij)
    {
      const int i = voigt9row[ij];
      const int j = voigt9col[ij];
      invfderf(ij, kl) = -invfmat(i, k) * invfmat(l, j);
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::InvVector6VoigtDiffByItself(
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& invfderf,
    const LINALG::Matrix<NUMDIM_, NUMDIM_>& invfmat)
{
#if 0
  // VERIFIED

//  std::cout << std::endl;
//  std::cout << std::endl;
  for (int ij=0; ij<MAT::NUM_STRESS_3D; ++ij)
  {
//    std::cout << "[";
    const int i = VOIGT6ROW[ij];
    const int j = VOIGT6COL[ij];
    for (int kl=0; kl<MAT::NUM_STRESS_3D; ++kl)
    {
      const int k = VOIGT6ROW[kl];
      const int l = VOIGT6COL[kl];
      invfderf(ij,kl) = -0.5*(invfmat(i,k)*invfmat(l,j) + invfmat(i,l)*invfmat(k,j));
//      std::cout << "invfderf("<<ij<<","<<kl<<") = ";
//      std::cout << "-0.5*(invfmat("<<i<<","<<k<<")*invfmat("<<l<<","<<j<<")+invfmat("<<i<<","<<l<<")*invfmat("<<k<<","<<j<<"));";
//      std::cout << "ct["<<i+1<<","<<k+1<<"]*ct["<<l+1<<","<<j+1<<"]+ct["<<i+1<<","<<l+1<<"]*ct["<<k+1<<","<<j+1<<"]";
//      std::cout << std::endl;
      if (ij >= NUMDIM_)
      {
#if 0
        invfderf(ij,kl) += -0.5*(invfmat(j,k)*invfmat(l,i) + invfmat(j,l)*invfmat(k,i));
//        std::cout << "+ct["<<j+1<<","<<k+1<<"]*ct["<<l+1<<","<<i+1<<"]+ct["<<j+1<<","<<l+1<<"]*ct["<<k+1<<","<<i+1<<"]";
#else
        invfderf(ij,kl) *= 2.0;
//        std::cout << "invfderf("<<ij<<","<<kl<<") *= 2.0;";
//        std::cout << std::endl;
#endif
      }
//      std::cout << ", ";
    }
//    std::cout << "]," << std::endl;
  }
#else
  invfderf(0, 0) = -0.5 * (invfmat(0, 0) * invfmat(0, 0) + invfmat(0, 0) * invfmat(0, 0));
  invfderf(1, 0) = -0.5 * (invfmat(1, 0) * invfmat(0, 1) + invfmat(1, 0) * invfmat(0, 1));
  invfderf(2, 0) = -0.5 * (invfmat(2, 0) * invfmat(0, 2) + invfmat(2, 0) * invfmat(0, 2));
  invfderf(3, 0) = -1.0 * (invfmat(0, 0) * invfmat(0, 1) + invfmat(0, 0) * invfmat(0, 1));
  invfderf(4, 0) = -1.0 * (invfmat(1, 0) * invfmat(0, 2) + invfmat(1, 0) * invfmat(0, 2));
  invfderf(5, 0) = -1.0 * (invfmat(2, 0) * invfmat(0, 0) + invfmat(2, 0) * invfmat(0, 0));

  invfderf(0, 1) = -0.5 * (invfmat(0, 1) * invfmat(1, 0) + invfmat(0, 1) * invfmat(1, 0));
  invfderf(1, 1) = -0.5 * (invfmat(1, 1) * invfmat(1, 1) + invfmat(1, 1) * invfmat(1, 1));
  invfderf(2, 1) = -0.5 * (invfmat(2, 1) * invfmat(1, 2) + invfmat(2, 1) * invfmat(1, 2));
  invfderf(3, 1) = -1.0 * (invfmat(0, 1) * invfmat(1, 1) + invfmat(0, 1) * invfmat(1, 1));
  invfderf(4, 1) = -1.0 * (invfmat(1, 1) * invfmat(1, 2) + invfmat(1, 1) * invfmat(1, 2));
  invfderf(5, 1) = -1.0 * (invfmat(2, 1) * invfmat(1, 0) + invfmat(2, 1) * invfmat(1, 0));

  invfderf(0, 2) = -0.5 * (invfmat(0, 2) * invfmat(2, 0) + invfmat(0, 2) * invfmat(2, 0));
  invfderf(1, 2) = -0.5 * (invfmat(1, 2) * invfmat(2, 1) + invfmat(1, 2) * invfmat(2, 1));
  invfderf(2, 2) = -0.5 * (invfmat(2, 2) * invfmat(2, 2) + invfmat(2, 2) * invfmat(2, 2));
  invfderf(3, 2) = -1.0 * (invfmat(0, 2) * invfmat(2, 1) + invfmat(0, 2) * invfmat(2, 1));
  invfderf(4, 2) = -1.0 * (invfmat(1, 2) * invfmat(2, 2) + invfmat(1, 2) * invfmat(2, 2));
  invfderf(5, 2) = -1.0 * (invfmat(2, 2) * invfmat(2, 0) + invfmat(2, 2) * invfmat(2, 0));

  invfderf(0, 3) = -0.5 * (invfmat(0, 0) * invfmat(1, 0) + invfmat(0, 1) * invfmat(0, 0));
  invfderf(1, 3) = -0.5 * (invfmat(1, 0) * invfmat(1, 1) + invfmat(1, 1) * invfmat(0, 1));
  invfderf(2, 3) = -0.5 * (invfmat(2, 0) * invfmat(1, 2) + invfmat(2, 1) * invfmat(0, 2));
  invfderf(3, 3) = -1.0 * (invfmat(0, 0) * invfmat(1, 1) + invfmat(0, 1) * invfmat(0, 1));
  invfderf(4, 3) = -1.0 * (invfmat(1, 0) * invfmat(1, 2) + invfmat(1, 1) * invfmat(0, 2));
  invfderf(5, 3) = -1.0 * (invfmat(2, 0) * invfmat(1, 0) + invfmat(2, 1) * invfmat(0, 0));

  invfderf(0, 4) = -0.5 * (invfmat(0, 1) * invfmat(2, 0) + invfmat(0, 2) * invfmat(1, 0));
  invfderf(1, 4) = -0.5 * (invfmat(1, 1) * invfmat(2, 1) + invfmat(1, 2) * invfmat(1, 1));
  invfderf(2, 4) = -0.5 * (invfmat(2, 1) * invfmat(2, 2) + invfmat(2, 2) * invfmat(1, 2));
  invfderf(3, 4) = -1.0 * (invfmat(0, 1) * invfmat(2, 1) + invfmat(0, 2) * invfmat(1, 1));
  invfderf(4, 4) = -1.0 * (invfmat(1, 1) * invfmat(2, 2) + invfmat(1, 2) * invfmat(1, 2));
  invfderf(5, 4) = -1.0 * (invfmat(2, 1) * invfmat(2, 0) + invfmat(2, 2) * invfmat(1, 0));

  invfderf(0, 5) = -0.5 * (invfmat(0, 2) * invfmat(0, 0) + invfmat(0, 0) * invfmat(2, 0));
  invfderf(1, 5) = -0.5 * (invfmat(1, 2) * invfmat(0, 1) + invfmat(1, 0) * invfmat(2, 1));
  invfderf(2, 5) = -0.5 * (invfmat(2, 2) * invfmat(0, 2) + invfmat(2, 0) * invfmat(2, 2));
  invfderf(3, 5) = -1.0 * (invfmat(0, 2) * invfmat(0, 1) + invfmat(0, 0) * invfmat(2, 1));
  invfderf(4, 5) = -1.0 * (invfmat(1, 2) * invfmat(0, 2) + invfmat(1, 0) * invfmat(2, 2));
  invfderf(5, 5) = -1.0 * (invfmat(2, 2) * invfmat(0, 0) + invfmat(2, 0) * invfmat(2, 0));
#endif

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::InvVector6VoigtTwiceDiffByItself(
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D * MAT::NUM_STRESS_3D>& invbvdderb,
    const LINALG::Matrix<NUMDIM_, NUMDIM_>& ibt)
{
  // VERIFIED

  //  std::cout << std::endl;
  //  std::cout << std::endl;

  for (int kl = 0; kl < MAT::NUM_STRESS_3D; ++kl)
  {
    const int k = VoigtMapping::Voigt6ToRow(kl);
    const int l = VoigtMapping::Voigt6ToCol(kl);
    for (int mn = kl; mn < MAT::NUM_STRESS_3D; ++mn)
    {
      const int m = VoigtMapping::Voigt6ToRow(mn);
      const int n = VoigtMapping::Voigt6ToCol(mn);
      const int klmn = MAT::NUM_STRESS_3D * kl + mn;
      const int mnkl = MAT::NUM_STRESS_3D * mn + kl;
      for (int ij = 0; ij < MAT::NUM_STRESS_3D; ++ij)
      {
        //    std::cout << "[\n";
        const int i = VoigtMapping::Voigt6ToRow(ij);
        const int j = VoigtMapping::Voigt6ToCol(ij);
        double invbvdderb_ijklmn;
        invbvdderb_ijklmn = 0.25 * ((ibt(i, m) * ibt(n, k) + ibt(i, n) * ibt(m, k)) * ibt(l, j) +
                                       ibt(i, k) * (ibt(l, m) * ibt(n, j) + ibt(l, n) * ibt(m, j)) +
                                       (ibt(i, m) * ibt(n, l) + ibt(i, n) * ibt(m, l)) * ibt(k, j) +
                                       ibt(i, l) * (ibt(k, m) * ibt(n, j) + ibt(k, n) * ibt(m, j)));
        //        std::cout << ""
        //             <<
        //             "(ct["<<i+1<<","<<m+1<<"]*ct["<<n+1<<","<<k+1<<"]+ct["<<i+1<<","<<n+1<<"]*ct["<<m+1<<","<<k+1<<"])*ct["<<l+1<<","<<j+1<<"]"
        //             <<
        //             "+ct["<<i+1<<","<<k+1<<"]*(ct["<<l+1<<","<<m+1<<"]*ct["<<n+1<<","<<j+1<<"]+ct["<<l+1<<","<<n+1<<"]*ct["<<m+1<<","<<j+1<<"])"
        //             <<
        //             "+(ct["<<i+1<<","<<m+1<<"]*ct["<<n+1<<","<<l+1<<"]+ct["<<i+1<<","<<n+1<<"]*ct["<<m+1<<","<<l+1<<"])*ct["<<k+1<<","<<j+1<<"]"
        //             <<
        //             "+ct["<<i+1<<","<<l+1<<"]*(ct["<<k+1<<","<<m+1<<"]*ct["<<n+1<<","<<j+1<<"]+ct["<<k+1<<","<<n+1<<"]*ct["<<m+1<<","<<j+1<<"])"
        //             << "";
        if (ij >= NUMDIM_)  // swap 'i' and 'j'
        {
#if 0
          invbvdderb_ijklmn += 0.25*(
              ( ibt(j,m)*ibt(n,k) + ibt(j,n)*ibt(m,k) )*ibt(l,i)
              + ibt(j,k)*( ibt(l,m)*ibt(n,i) + ibt(l,n)*ibt(m,i) )
              + ( ibt(j,m)*ibt(n,l) + ibt(j,n)*ibt(m,l) )*ibt(k,i)
              + ibt(j,l)*( ibt(k,m)*ibt(n,i) + ibt(k,n)*ibt(m,i) )
            );
//          std::cout << ""
//               << "+(ct["<<j+1<<","<<m+1<<"]*ct["<<n+1<<","<<k+1<<"]+ct["<<j+1<<","<<n+1<<"]*ct["<<m+1<<","<<k+1<<"])*ct["<<l+1<<","<<i+1<<"]"
//               << "+ct["<<j+1<<","<<k+1<<"]*(ct["<<l+1<<","<<m+1<<"]*ct["<<n+1<<","<<i+1<<"]+ct["<<l+1<<","<<n+1<<"]*ct["<<m+1<<","<<i+1<<"])"
//               << "+(ct["<<j+1<<","<<m+1<<"]*ct["<<n+1<<","<<l+1<<"]+ct["<<j+1<<","<<n+1<<"]*ct["<<m+1<<","<<l+1<<"])*ct["<<k+1<<","<<i+1<<"]"
//               << "+ct["<<j+1<<","<<l+1<<"]*(ct["<<k+1<<","<<m+1<<"]*ct["<<n+1<<","<<i+1<<"]+ct["<<k+1<<","<<n+1<<"]*ct["<<m+1<<","<<i+1<<"])"
//               << "";
#else
          invbvdderb_ijklmn *= 2.0;
#endif
        }
        invbvdderb(ij, klmn) = invbvdderb_ijklmn;
        if (mn != kl) invbvdderb(ij, mnkl) = invbvdderb_ijklmn;
        //        std::cout << ",\n";
      }
      //      std::cout << "";
    }
    //    std::cout << "],\n";
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::SqVector6VoigtDiffByItself(
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& sqfderf,
    const LINALG::Matrix<NUMDIM_, NUMDIM_>& fmat, ::UTILS::VOIGT::NotationType outvoigt6)
{
  // VERIFIED

#if 0
  // identity 2-tensor
  LINALG::Matrix<NUMDIM_,NUMDIM_> id(true);
  for (int i=0; i<NUMDIM_; ++i) id(i,i) = 1.0;

  // (F.F)_{,F} with F^T=F
//  std::cout << std::endl;
  for (int ij=0; ij<MAT::NUM_STRESS_3D; ++ij)
  {
//    std::cout << "[";
    const int i = VOIGT6ROW[ij];
    const int j = VOIGT6COL[ij];
    for (int kl=0; kl<MAT::NUM_STRESS_3D; ++kl)
    {
      const int k = VOIGT6ROW[kl];
      const int l = VOIGT6COL[kl];
      sqfderf(ij,kl) = id(i,k)*fmat(l,j) + id(j,l)*fmat(i,k);
//      std::cout << "id["<<i+1<<","<<k+1<<"]*St["<<l+1<<","<<j+1<<"]+id["<<j+1<<","<<l+1<<"]*St["<<i+1<<","<<k+1<<"]";
      if ( (outvoigt6 == voigt6_strain) and (ij >= NUMDIM_) )
      {
        sqfderf(ij,kl) += id(j,k)*fmat(l,i) + id(i,l)*fmat(j,k);
//        std::cout << "+id["<<j+1<<","<<k+1<<"]*St["<<l+1<<","<<i+1<<"]+id["<<i+1<<","<<l+1<<"]*St["<<j+1<<","<<k+1<<"]";
      }
//      std::cout << ", ";
    }
//    std::cout << "]," << std::endl;
  }
#else
  if (outvoigt6 != ::UTILS::VOIGT::NotationType::strain)
    dserror("Can only produce row of strain-like type");

  sqfderf(0, 0) = 2.0 * fmat(0, 0);
  sqfderf(1, 0) = 0.0;
  sqfderf(2, 0) = 0.0;
  sqfderf(3, 0) = fmat(1, 0) + fmat(0, 1);
  sqfderf(4, 0) = 0.0;
  sqfderf(5, 0) = fmat(2, 0) + fmat(0, 2);

  sqfderf(0, 1) = 0.0;
  sqfderf(1, 1) = 2.0 * fmat(1, 1);
  sqfderf(2, 1) = 0.0;
  sqfderf(3, 1) = fmat(1, 0) + fmat(0, 1);
  sqfderf(4, 1) = fmat(2, 1) + fmat(1, 2);
  sqfderf(5, 1) = 0.0;

  sqfderf(0, 2) = 0.0;
  sqfderf(1, 2) = 0.0;
  sqfderf(2, 2) = 2.0 * fmat(2, 2);
  sqfderf(3, 2) = 0.0;
  sqfderf(4, 2) = fmat(2, 1) + fmat(1, 2);
  sqfderf(5, 2) = fmat(2, 0) + fmat(0, 2);

  sqfderf(0, 3) = fmat(0, 1);
  sqfderf(1, 3) = fmat(0, 1);
  sqfderf(2, 3) = 0.0;
  sqfderf(3, 3) = fmat(1, 1) + fmat(0, 0);
  sqfderf(4, 3) = fmat(0, 2);
  sqfderf(5, 3) = fmat(2, 1);

  sqfderf(0, 4) = 0.0;
  sqfderf(1, 4) = fmat(1, 2);
  sqfderf(2, 4) = fmat(1, 2);
  sqfderf(3, 4) = fmat(0, 2);
  sqfderf(4, 4) = fmat(2, 2) + fmat(1, 1);
  sqfderf(5, 4) = fmat(1, 0);

  sqfderf(0, 5) = fmat(2, 0);
  sqfderf(1, 5) = 0.0;
  sqfderf(2, 5) = fmat(2, 0);
  sqfderf(3, 5) = fmat(2, 1);
  sqfderf(4, 5) = fmat(1, 0);
  sqfderf(5, 5) = fmat(2, 2) + fmat(0, 0);
#endif

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
void DRT::ELEMENTS::So_sh8p8::SqVector9VoigtDiffByItself(
  LINALG::Matrix<NUMDFGR_,NUMDFGR_>& sqfderf,
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& fmat,
  const bool transpose
  )
{
  // identity 2-tensor
  LINALG::Matrix<NUMDIM_,NUMDIM_> id(true);
  for (int i=0; i<NUMDIM_; ++i) id(i,i) = 1.0;

  // (F^T.F)_{,F}
//  std::cout << std::endl;
  for (int ij=0; ij<NUMDFGR_; ++ij)
  {
//    std::cout << "[";
    const int i = VOIGT9ROW_INCONSISTENT_[ij];
    const int j = VOIGT9COL_INCONSISTENT_[ij];
    for (int kl=0; kl<NUMDFGR_; ++kl)
    {
      const int k = VOIGT9ROW_INCONSISTENT_[kl];
      const int l = VOIGT9COL_INCONSISTENT_[kl];
//      std::cout << "i=" << i << ", j=" << j << ", k=" << k << ", l=" << l << std::endl;
      if (transpose)  // swap indices of fmat
        sqfderf(ij,kl) = id(i,k)*fmat(j,l) + id(j,l)*fmat(k,i);
      else
        sqfderf(ij,kl) = id(i,k)*fmat(l,j) + id(j,l)*fmat(i,k);
//      std::cout <<
"id["<<i+1<<","<<k+1<<"]*St["<<l+1<<","<<j+1<<"]+id["<<j+1<<","<<l+1<<"]*St["<<i+1<<","<<k+1<<"]";
//      std::cout << ", ";
    }
//    std::cout << "]," << std::endl;
  }

  return;
}
*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::SqVector6VoigtTwiceDiffByItself(
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D * MAT::NUM_STRESS_3D>& sqfdderf,
    const LINALG::Matrix<NUMDIM_, NUMDIM_>& fmat)
{
#if 0
  // identity 2-tensor
  LINALG::Matrix<NUMDIM_,NUMDIM_> id(true);
  for (int i=0; i<NUMDIM_; ++i) id(i,i) = 1.0;

  // VERIFIED

  // (F^T.F)_{,FF} with F^T=F
//  std::cout << std::endl;
  for (int ij=0; ij<MAT::NUM_STRESS_3D; ++ij)
  {
//    std::cout << "[";
    const int i = VOIGT6ROW[ij];
    const int j = VOIGT6COL[ij];
    for (int kl=0; kl<MAT::NUM_STRESS_3D; ++kl)
    {
      const int k = VOIGT6ROW[kl];
      const int l = VOIGT6COL[kl];
      for (int mn=0; mn<MAT::NUM_STRESS_3D; ++mn)
      {
        const int m = VOIGT6ROW[mn];
        const int n = VOIGT6COL[mn];
        const int klmn = MAT::NUM_STRESS_3D*kl + mn;
        double sqfdderf_ijklmn = 0.25*(id(i,k)*id(l,m)*id(j,n)+id(j,l)*id(i,m)*id(k,n)
                                       +id(i,k)*id(l,n)*id(j,m)+id(j,l)*id(i,n)*id(k,m)  // swap 'm' and 'n'
                                       +id(i,l)*id(k,m)*id(j,n)+id(j,k)*id(i,m)*id(l,n)  // swap 'k' and 'l'
                                       +id(i,l)*id(k,n)*id(j,m)+id(j,k)*id(i,n)*id(l,m));  // swap 'm' and 'n' as well as 'k' and 'l'
        if (ij >= NUMDIM_)  // swap 'i' and 'j'
        {
          sqfdderf_ijklmn += 0.25*(id(j,k)*id(l,m)*id(i,n)+id(i,l)*id(j,m)*id(k,n)  // swap 'i' and 'j'
                                   +id(j,k)*id(l,n)*id(i,m)+id(i,l)*id(j,n)*id(k,m)  // swap 'i' and 'j' as well as 'm' and 'n'
                                   +id(j,l)*id(k,m)*id(i,n)+id(i,k)*id(j,m)*id(l,n)  // swap 'i' and 'j' as well as 'k' and 'l'
                                   +id(j,l)*id(k,n)*id(i,m)+id(i,k)*id(j,n)*id(l,m) );  // swap 'i' and 'j' as well as 'm' and 'n' as well as 'k' and 'l'
        }
        sqfdderf(ij,klmn) = sqfdderf_ijklmn;
//        std::cout << sqfdderf_ijklmn;
//        std::cout << ", ";
      }
//      std::cout << "\n";
    }
//    std::cout << "],\n";
  }
#else
  sqfdderf.Clear();

  sqfdderf(0, 0) = 2.0;
  sqfdderf(0, 21) = 0.5;
  sqfdderf(0, 35) = 0.5;

  sqfdderf(1, 7) = 2.0;
  sqfdderf(1, 21) = 0.5;
  sqfdderf(1, 28) = 0.5;

  sqfdderf(2, 14) = 2.0;
  sqfdderf(2, 28) = 0.5;
  sqfdderf(2, 35) = 0.5;

  sqfdderf(3, 3) = 1.0;
  sqfdderf(3, 9) = 1.0;
  sqfdderf(3, 18) = 1.0;
  sqfdderf(3, 19) = 1.0;
  sqfdderf(3, 29) = 0.5;
  sqfdderf(3, 34) = 0.5;

  sqfdderf(4, 9) = 1.0;
  sqfdderf(4, 16) = 1.0;
  sqfdderf(4, 23) = 0.5;
  sqfdderf(4, 25) = 1.0;
  sqfdderf(4, 26) = 1.0;
  sqfdderf(4, 33) = 0.5;

  sqfdderf(5, 5) = 1.0;
  sqfdderf(5, 17) = 1.0;
  sqfdderf(5, 22) = 0.5;
  sqfdderf(5, 27) = 0.5;
  sqfdderf(5, 30) = 1.0;
  sqfdderf(5, 32) = 1.0;
#endif

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::SqVector6VoigtTwiceDiffByItself(
    int* isqfdderf,  //[MAT::NUM_STRESS_3D*6];
    LINALG::Matrix<MAT::NUM_STRESS_3D, 6>& sqfdderf)
{
  isqfdderf[MAT::NUM_STRESS_3D * 0 + 0] = 0;
  sqfdderf(0, 0) = 2.0;
  isqfdderf[MAT::NUM_STRESS_3D * 0 + 1] = 21;
  sqfdderf(0, 1) = 0.5;
  isqfdderf[MAT::NUM_STRESS_3D * 0 + 2] = 35;
  sqfdderf(0, 2) = 0.5;
  isqfdderf[MAT::NUM_STRESS_3D * 0 + 3] = -1;
  sqfdderf(0, 3) = 0.0;  // dummy
  isqfdderf[MAT::NUM_STRESS_3D * 0 + 4] = -1;
  sqfdderf(0, 4) = 0.0;  // dummy
  isqfdderf[MAT::NUM_STRESS_3D * 0 + 5] = -1;
  sqfdderf(0, 5) = 0.0;  // dummy

  isqfdderf[MAT::NUM_STRESS_3D * 1 + 0] = 7;
  sqfdderf(1, 0) = 2.0;
  isqfdderf[MAT::NUM_STRESS_3D * 1 + 1] = 21;
  sqfdderf(1, 1) = 0.5;
  isqfdderf[MAT::NUM_STRESS_3D * 1 + 2] = 28;
  sqfdderf(1, 2) = 0.5;
  isqfdderf[MAT::NUM_STRESS_3D * 1 + 3] = -1;
  sqfdderf(1, 3) = 0.0;  // dummy
  isqfdderf[MAT::NUM_STRESS_3D * 1 + 4] = -1;
  sqfdderf(1, 4) = 0.0;  // dummy
  isqfdderf[MAT::NUM_STRESS_3D * 1 + 5] = -1;
  sqfdderf(1, 5) = 0.0;  // dummy

  isqfdderf[MAT::NUM_STRESS_3D * 2 + 0] = 14;
  sqfdderf(2, 0) = 2.0;
  isqfdderf[MAT::NUM_STRESS_3D * 2 + 1] = 28;
  sqfdderf(2, 1) = 0.5;
  isqfdderf[MAT::NUM_STRESS_3D * 2 + 2] = 35;
  sqfdderf(2, 2) = 0.5;
  isqfdderf[MAT::NUM_STRESS_3D * 2 + 3] = -1;
  sqfdderf(2, 3) = 0.0;  // dummy
  isqfdderf[MAT::NUM_STRESS_3D * 2 + 4] = -1;
  sqfdderf(2, 4) = 0.0;  // dummy
  isqfdderf[MAT::NUM_STRESS_3D * 2 + 5] = -1;
  sqfdderf(2, 5) = 0.0;  // dummy

  isqfdderf[MAT::NUM_STRESS_3D * 3 + 0] = 3;
  sqfdderf(3, 0) = 1.0;
  isqfdderf[MAT::NUM_STRESS_3D * 3 + 1] = 9;
  sqfdderf(3, 1) = 1.0;
  isqfdderf[MAT::NUM_STRESS_3D * 3 + 2] = 18;
  sqfdderf(3, 2) = 1.0;
  isqfdderf[MAT::NUM_STRESS_3D * 3 + 3] = 19;
  sqfdderf(3, 3) = 1.0;
  isqfdderf[MAT::NUM_STRESS_3D * 3 + 4] = 29;
  sqfdderf(3, 4) = 0.5;
  isqfdderf[MAT::NUM_STRESS_3D * 3 + 5] = 34;
  sqfdderf(3, 5) = 0.5;

  isqfdderf[MAT::NUM_STRESS_3D * 4 + 0] = 9;
  sqfdderf(4, 0) = 1.0;
  isqfdderf[MAT::NUM_STRESS_3D * 4 + 1] = 16;
  sqfdderf(4, 1) = 1.0;
  isqfdderf[MAT::NUM_STRESS_3D * 4 + 2] = 23;
  sqfdderf(4, 2) = 0.5;
  isqfdderf[MAT::NUM_STRESS_3D * 4 + 3] = 25;
  sqfdderf(4, 3) = 1.0;
  isqfdderf[MAT::NUM_STRESS_3D * 4 + 4] = 26;
  sqfdderf(4, 4) = 1.0;
  isqfdderf[MAT::NUM_STRESS_3D * 4 + 5] = 33;
  sqfdderf(4, 5) = 0.5;

  isqfdderf[MAT::NUM_STRESS_3D * 5 + 0] = 5;
  sqfdderf(5, 0) = 1.0;
  isqfdderf[MAT::NUM_STRESS_3D * 5 + 1] = 17;
  sqfdderf(5, 1) = 1.0;
  isqfdderf[MAT::NUM_STRESS_3D * 5 + 2] = 22;
  sqfdderf(5, 2) = 0.5;
  isqfdderf[MAT::NUM_STRESS_3D * 5 + 3] = 27;
  sqfdderf(5, 3) = 0.5;
  isqfdderf[MAT::NUM_STRESS_3D * 5 + 4] = 30;
  sqfdderf(5, 4) = 1.0;
  isqfdderf[MAT::NUM_STRESS_3D * 5 + 5] = 32;
  sqfdderf(5, 5) = 1.0;

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Matrix2TensorToMatrix6x9Voigt(
    LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDFGR_>& bm, const LINALG::Matrix<NUMDIM_, NUMDIM_>& bt,
    const bool transpose)
{
  // VERIFIED

  //  std::cout << std::endl;
  for (int kl = 0; kl < NUMDFGR_; ++kl)
  {
    const int k = VOIGT9ROW_INCONSISTENT_[kl];
    const int l = VOIGT9COL_INCONSISTENT_[kl];
    for (int ij = 0; ij < MAT::NUM_STRESS_3D; ++ij)
    {
      //    std::cout << "[";
      const int i = VoigtMapping::Voigt6ToRow(ij);
      const int j = VoigtMapping::Voigt6ToCol(ij);
      if (j == l)
        if (transpose)
        {
          bm(ij, kl) = bt(k, i);
          //      std::cout << "bt["<<k+1<<","<<i+1<<"]";
        }
        else
          bm(ij, kl) = bt(i, k);
      else if ((ij >= NUMDIM_) and (i == l))
        if (transpose)
        {
          bm(ij, kl) = bt(k, j);
          //      std::cout << "bt["<<k+1<<","<<j+1<<"]";
        }
        else
          bm(ij, kl) = bt(j, k);
      else
      {
        bm(ij, kl) = 0.0;
        //      std::cout << "0";
      }
      //      std::cout << ", ";
    }
    //    std::cout << "]," << std::endl;
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Matrix2TensorToLeftRightProductMatrix6x6Voigt(
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& bm,  ///< (out) 6x6 Voigt matrix
    const LINALG::Matrix<NUMDIM_, NUMDIM_>& bt,                  ///< (in) 3x3 matrix of 2-tensor
    const bool transpose,                                        ///< 3x3 input matrix is transposed
    ::UTILS::VOIGT::NotationType rowvoigt6,  ///< 6-Voigt vector layout on rows of 6x6 matrix
    ::UTILS::VOIGT::NotationType colvoigt6   ///< 6-Voigt vector layout on columns of 6x6 matrix
)
{
  for (int ab = 0; ab < MAT::NUM_STRESS_3D; ++ab)
  {
    const int a = VoigtMapping::Voigt6ToRow(ab);
    const int b = VoigtMapping::Voigt6ToCol(ab);
    for (int AB = 0; AB < MAT::NUM_STRESS_3D; ++AB)
    {
      const int A = VoigtMapping::Voigt6ToRow(AB);
      const int B = VoigtMapping::Voigt6ToCol(AB);
      if (transpose)
      {
        bm(AB, ab) = bt(A, a) * bt(B, b);
        if (ab >= NUMDIM_) bm(AB, ab) += bt(A, b) * bt(B, a);
      }
      else
      {
        bm(AB, ab) = bt(a, A) * bt(b, B);
        if (ab >= NUMDIM_) bm(AB, ab) += bt(b, A) * bt(a, B);
      }
      if ((colvoigt6 == ::UTILS::VOIGT::NotationType::stress) and (ab >= NUMDIM_))
        bm(AB, ab) *= 0.5;
      if ((rowvoigt6 == ::UTILS::VOIGT::NotationType::strain) and (AB >= NUMDIM_))
        bm(AB, ab) *= 2.0;
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::StretchTensor(
    double* detut,                              // determinant of material stretch tensor
    LINALG::Matrix<NUMDIM_, NUMDIM_>* ut,       // material stretch tensor
    LINALG::Matrix<NUMDIM_, NUMDIM_>* invut,    // inverse material stretch tensor
    const LINALG::Matrix<NUMDIM_, NUMDIM_>& ct  // right Cauchy-Green tensor
)
{
  if ((ut == NULL) and (invut == NULL))
    dserror("Senseless call: You do not want to compute anything");

  // set identity tensor
  LINALG::Matrix<NUMDIM_, NUMDIM_> it(true);
  for (int i = 0; i < NUMDIM_; ++i) it(i, i) = 1.0;

  // squared right Cauchy-Green deformation tensor
  // C^2 = C . C
  LINALG::Matrix<NUMDIM_, NUMDIM_> c2t;
  c2t.MultiplyTN(ct, ct);

  // invariants of right Cauchy-Green tensor
  // 1st principal invariant: I_C = tr(C)
  const double ci = ct(0, 0) + ct(1, 1) + ct(2, 2);
  // 2nd principal invariant: II_C = 1/2 ( tr(C)^2 - tr(C^2) )
  const double c2i = c2t(0, 0) + c2t(1, 1) + c2t(2, 2);
  const double cii = 0.5 * (ci * ci - c2i);
  // 3rd principal invariant: III_C = det(C)
  const double ciii = ct.Determinant();

  // determination of I_U acc. to [1]
  double ui = 0.0;
  {
    // auxiliar variables to get trace of material stretch tensor U
    const double xi = (2.0 * ci * ci * ci - 9.0 * ci * cii + 27.0 * ciii) / 27.0;

    double eta = (4.0 * cii * cii * cii - ci * ci * cii * cii + 4.0 * ci * ci * ci * ciii -
                     18.0 * ci * cii * ciii + 27.0 * ciii * ciii) /
                 27.0;
    if (eta < 0.0)
    {
      if (fabs(eta) < EPS6)
        eta = 0.0;
      else
        dserror("Trouble with negative eta=%g", eta);
    }

    // const double zeta = -2.0*ci/3.0
    //                   + pow(xi+sqrt(eta), 1.0/3.0)
    //                   + pow(xi-sqrt(eta), 1.0/3.0);
    double zeta = -2.0 * ci / 3.0;
    const double xiplussqrteta = 32.0 * (xi + sqrt(eta));
    if (xiplussqrteta < 0.0)
      zeta -= std::pow(fabs(xiplussqrteta), (1.0 / 3.0));
    else
      zeta += std::pow(xiplussqrteta, (1.0 / 3.0));
    const double ximinussqrteta = 32.0 * (xi - sqrt(eta));
    if (ximinussqrteta < 0.0)
      zeta -= std::pow(fabs(ximinussqrteta), (1.0 / 3.0));
    else
      zeta += std::pow(ximinussqrteta, (1.0 / 3.0));

    // invariants of material stretch tensor U
    // 1st invariant: I_U = tr(U)
    if (fabs(zeta + 2.0 * ci) < EPS12)
    {
      ui = sqrt(ci + 2.0 * sqrt(cii));
    }
    else
    {
      const double aux = sqrt(2.0 * ci + zeta);
      ui = 0.5 * (aux + sqrt(2.0 * ci - zeta + 16.0 * sqrt(ciii) / aux));
    }
  }

  // 2nd and 3rd invariant of material stretch tensor U
  // 2nd invariant: II_U = 1/2 * (I_U^2 - I_C)
  const double uii = 0.5 * (ui * ui - ci);  // OR // sqrt(cii + 2.0*sqrt(ciii)*ui);
  // 3rd invariant: III_U = det(U) = sqrt(III_C)
  const double uiii = sqrt(ciii);

  // inverse of material stretch tensor U^{-1}
  // Hoger & Carlson [1] identified
  //     U^{-1} = [ III_U^2*(III_U+I_U*I_C)
  //                + I_U^2*(I_U*III_C + III_U*II_C) ]^{-1}
  //            * [ I_U*(I_U*II_U - III_U)*C^2
  //                - (I_U*II_U - III_U)*(III_U +I_U*I_C)*C
  //                + { II_U*III_U*(III_U+I_U*I_C)
  //                    + I_U^2*(II_U*II_C+III_C) }*I ]
  //            = 1/denom * [ pc2*C^2 + pc*C + pi*I ]
  if (invut != NULL)
  {
    const double denom = uiii * uiii * (uiii + ui * ci) + ui * ui * (ui * ciii + uiii * cii);
    const double pc2 = ui * (ui * uii - uiii);
    const double pc = -(ui * uii - uiii) * (uiii + ui * ci);
    const double pi = uii * uiii * (uiii + ui * ci) + ui * ui * (uii * cii + ciii);
    for (int j = 0; j < NUMDIM_; j++)
    {
      for (int i = 0; i < NUMDIM_; i++)
      {
        const double invut_ij = (pc2 * c2t(i, j) + pc * ct(i, j) + pi * it(i, j)) / denom;
        (*invut)(i, j) = invut_ij;
      }
    }
  }

  // material stretch tensor U
  // Hoger & Carlson [1] wrote
  //     U = [ II_U * { II_U * (II_U + I_C) + II_C } + III_C ]^{-1}
  //       * [ -(I_U*II_U - III_U)*C^2
  //           + (I_U*II_U - III_U)*(II_U + I_C)*C
  //           + { I_U*III_U + III_U * ( II_U*(II_U+I_C)+II_C ) }*I ]
  //       = 1/denom * [ pc2*C^2 + pc*C + pi*I ]
  //
  // alternative:
  // U could be calculated based on R later: U = R^T . F
  if (ut != NULL)
  {
    const double denom = uii * (uii * (uii + ci) + cii) + ciii;
    const double pc2 = -(ui * uii - uiii);
    const double pc = (ui * uii - uiii) * (uii + ci);
    const double pi = ui * ciii + uiii * (uii * (uii + ci) + cii);
    for (int j = 0; j < NUMDIM_; j++)
    {
      for (int i = 0; i < NUMDIM_; i++)
      {
        const double ut_ij = (pc2 * c2t(i, j) + pc * ct(i, j) + pi * it(i, j)) / denom;
        (*ut)(i, j) = ut_ij;
      }
    }
  }

  // determinat of right stretch tensor
  if (detut != NULL)
  {
    *detut = uiii;
  }

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_sh8p8::SymSpectralDecompJacIter(LINALG::Matrix<NUMDIM_, NUMDIM_>& ew,
    LINALG::Matrix<NUMDIM_, NUMDIM_>& ev, const LINALG::Matrix<NUMDIM_, NUMDIM_>& at,
    const double itertol, const int itermax)
{
  // sum of all entries (moduli) in #at
  double asum = 0.0;

  // initialise eigenvalue tensor and eigenvector tensor
#if 0
  asum = at.Norm1();
  ew.Update(at);
  ev.Clear();
  for (int idim=0; idim<NUMDIM_; idim++) ev(idim,idim) = 1.0;
#else
  {
    asum = 0.0;
    for (int jdim = 0; jdim < NUMDIM_; jdim++)
    {
      for (int idim = 0; idim < NUMDIM_; idim++)
      {
        asum += fabs(at(idim, jdim));
        ew(idim, jdim) = at(idim, jdim);
        ev(idim, jdim) = 0.0;
      }
      ev(jdim, jdim) = 1.0;
    }
  }
#endif

  // check for trivial problem
  if (asum < EPS12)
  {
    ew.Clear();
    return 0;
  }

  // scale sum of at compenents to achieve relative convergence check
  asum /= (double)(NUMDIM_ * NUMDIM_);

  // reduce ew to diagonal (the eigenvalues)
  double itercnt = 0;  // initialise iteration index <i>
  while (itercnt < itermax)
  {
    double vsum = 0.0;  // sum of all subtriangluar entries
    // loop lower triangle
    for (int jdim = 1; jdim < NUMDIM_; jdim++)
    {
      for (int idim = 0; idim < jdim; idim++)
      {
        // sum of all triag entries
        vsum += ew(idim, jdim);

        // rotation angle th
        // 2*th = atan(2*evt(idim,jdim)/(ew[idim,idim]-ew[jdim,jdim])
        const double th = 0.5 * atan2(2.0 * ew(idim, jdim), ew(idim, idim) - ew(jdim, jdim));
        const double sith = sin(th);  // sine of rotation angle
        const double coth = cos(th);  // cosine of rotation angle

        // this defines the rotation matrix,
        // e.g.
        //
        //             [ T_{idim,idim}  0  T_{idim,jdim} ]
        //   T^<i+1> = [             0  1              0 ]
        //             [ T_{jdim,idim}  0  T_{jdim,jdim} ]
        //
        //       [ cos(th)  0  -sin(th) ]
        //     = [       0  1         0 ]
        //       [ sin(th)  0   cos(th) ]

        // update eigenvector matrix by right-multiplying with T
        // T is mostly 0 thus it is more efficient to do explicitly
        //    ev^<i+1> = ev^<i> . T^<i+1>
        for (int kdim = 0; kdim < NUMDIM_; kdim++)
        {
          const double evki = ev(kdim, idim);
          ev(kdim, idim) = coth * evki + sith * ev(kdim, jdim);
          ev(kdim, jdim) = -sith * evki + coth * ev(kdim, jdim);
        }

        // update eigenvalue tensor by right-multiplying with T and
        // left-multiplying with transposed T
        //    ew^<i+1> = transposed(T^<i+1>) . ew^<i> . T^<i+1>
        // modify "idim" and "jdim" columns
        for (int kdim = 0; kdim < NUMDIM_; kdim++)
        {
          const double ewki = ew(kdim, idim);
          ew(kdim, idim) = coth * ewki + sith * ew(kdim, jdim);
          ew(kdim, jdim) = -sith * ewki + coth * ew(kdim, jdim);
        }
        // modify diagonal terms
        ew(idim, idim) = coth * ew(idim, idim) + sith * ew(jdim, idim);
        ew(jdim, jdim) = -sith * ew(idim, jdim) + coth * ew(jdim, jdim);
        ew(idim, jdim) = 0.0;
        // make symmetric
        for (int kdim = 0; kdim < NUMDIM_; kdim++)
        {
          ew(idim, kdim) = ew(kdim, idim);
          ew(jdim, kdim) = ew(kdim, jdim);
        }
      }
    }

    // check convergence
    if (fabs(vsum) / asum < itertol)
    {
      break;
    }
    // increment iteration index
    itercnt += 1;
  }

  // check if iteration loop diverged
  int err = 1;
  if (itercnt == itermax)
  {
    err = 1;  // failed
    dserror("Divergent spectral decomposition (Jacobi's iterative method)!");
  }
  else
  {
    err = 0;  // passed
  }

  return err;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::ExtractDispAndPres(std::vector<double>& mystat,
    LINALG::Matrix<NUMDISP_, 1>& mydisp, LINALG::Matrix<NUMPRES_, 1>& mypres)
{
  for (int inod = 0; inod < NUMNOD_; ++inod)
  {
    for (int idis = 0; idis < NODDISP_; ++idis)
      mydisp(idis + (inod * NODDISP_), 0) = mystat[idis + (inod * NODDOF_)];
    mypres(inod, 0) = mystat[NODDISP_ + 0 + (inod * NODDOF_)];
  }
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::BuildElementMatrix(LINALG::Matrix<NUMDOF_, NUMDOF_>* mat,
    const LINALG::Matrix<NUMDISP_, NUMDISP_>* matdd,
    const LINALG::Matrix<NUMDISP_, NUMPRES_>* matdp,
    const LINALG::Matrix<NUMPRES_, NUMDISP_>* matpd,
    const LINALG::Matrix<NUMPRES_, NUMPRES_>* matpp)
{
  const int* d2dp = &(DISPTODISPPRES_[0]);
  const int* p2dp = &(PRESTODISPPRES_[0]);

  // k_dd
  for (int j = 0; j < NUMDISP_; ++j)
  {
    const int J = d2dp[j];
    for (int i = 0; i < NUMDISP_; ++i)
    {
      const int I = d2dp[i];
      if (matdd != NULL)
        (*mat)(I, J) = (*matdd)(i, j);
      else
        (*mat)(I, J) = 0.0;
    }
  }

  // k_dp
  if (matdp != NULL)
  {
    for (int l = 0; l < NUMPRES_; ++l)
    {
      const int L = p2dp[l];
      for (int i = 0; i < NUMDISP_; ++i)
      {
        const int I = d2dp[i];
        (*mat)(I, L) = (*matdp)(i, l);
      }
    }
  }
  else
  {
    for (int l = 0; l < NUMPRES_; ++l)
    {
      const int L = p2dp[l];
      for (int i = 0; i < NUMDISP_; ++i)
      {
        const int I = d2dp[i];
        (*mat)(I, L) = 0.0;
      }
    }
  }

  // k_pd
  if (matpd != NULL)
  {
    for (int j = 0; j < NUMDISP_; ++j)
    {
      const int J = d2dp[j];
      for (int k = 0; k < NUMPRES_; ++k)
      {
        const int K = p2dp[k];
        (*mat)(K, J) = (*matpd)(k, j);
      }
    }
  }
  else if (matdp != NULL)
  {
    for (int j = 0; j < NUMDISP_; ++j)
    {
      const int J = d2dp[j];
      for (int k = 0; k < NUMPRES_; ++k)
      {
        const int K = p2dp[k];
        (*mat)(K, J) = (*matdp)(j, k);
      }
    }
  }
  else
  {
    for (int j = 0; j < NUMDISP_; ++j)
    {
      const int J = d2dp[j];
      for (int k = 0; k < NUMPRES_; ++k)
      {
        const int K = p2dp[k];
        (*mat)(K, J) = 0.0;
      }
    }
  }

  // k_pp
  if (matpp != NULL)
  {
    for (int l = 0; l < NUMPRES_; ++l)
    {
      const int L = p2dp[l];
      for (int k = 0; k < NUMPRES_; ++k)
      {
        const int K = p2dp[k];
        (*mat)(K, L) = (*matpp)(k, l);
      }
    }
  }
  else
  {
    for (int l = 0; l < NUMPRES_; ++l)
    {
      const int L = p2dp[l];
      for (int k = 0; k < NUMPRES_; ++k)
      {
        const int K = p2dp[k];
        (*mat)(K, L) = 0.0;
      }
    }
  }

  // done
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::BuildElementVector(LINALG::Matrix<NUMDOF_, 1>* vct,
    const LINALG::Matrix<NUMDISP_, 1>* vctd, const LINALG::Matrix<NUMPRES_, 1>* vctp)
{
  const int* d2dp = &(DISPTODISPPRES_[0]);
  const int* p2dp = &(PRESTODISPPRES_[0]);

  vct->Clear();

  // r_d
  if (vctd != NULL)
  {
    for (int i = 0; i < NUMDISP_; ++i)
    {
      const int I = d2dp[i];
      (*vct)(I, 0) = (*vctd)(i, 0);
    }
  }

  // r_p
  if (vctp != NULL)
  {
    for (int k = 0; k < NUMPRES_; ++k)
    {
      const int K = p2dp[k];
      (*vct)(K, 0) = (*vctp)(k, 0);
    }
  }

  // What shall we do with a drunken sailor?
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::AssembleVolume(
    Teuchos::ParameterList& params,  ///< parameter list for in 'n' out
    const double& elevol             ///< current element volume
)
{
  const double totvol = params.get<double>("volume");
  params.set("volume", totvol + elevol);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/* DEBUG ONLY */
void DRT::ELEMENTS::So_sh8p8::GnuplotOut(
    Teuchos::ParameterList& params,  ///< parameter list for in 'n' out
    std::vector<double>& state,      ///< current state vector, i.e. displacements and pressure DOFs
    LINALG::Matrix<NUMDOF_, 1>& resid,  ///< current internal force / incompressibility residual
    LINALG::Matrix<NUMDOF_, NUMDOF_>& tangent  ///< current tangent of inter force WRT state
)
{
  // last call time: we iterate to find solution, but we only want to write final state
  static double lasttime;
  static double oldtime = 0.0;
  // store current output line
  static std::ostringstream txtline;

  // store
  static std::vector<double> laststate(NUMDOF_, 0);
  static std::vector<double> oldstate(NUMDOF_, 0);
  static LINALG::Matrix<NUMDOF_, NUMDOF_> lasttangent;

  const double time = params.get<double>("total time");
  const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
  const std::string txtname = filebase + ".sosh8p8.txt";

  // store last converged state
  if (not(time == lasttime))
  {
    oldtime = lasttime;
    oldstate = laststate;
  }

  // control file
  if (time == 0.0)
  {
    const std::string gpltname = filebase + ".sosh8p8.gplt";
    std::ofstream gpltfile;
    gpltfile.open(gpltname.c_str(), std::ios_base::out);

    // header
    gpltfile << "# Gnuplot script" << std::endl;

    // print residual over time
    gpltfile << "plot \\" << std::endl;
    // for (int d=0; d<NUMDOF_; d+=4) {  // x-dir
    for (int d = 3; d < NUMDOF_; d += 4)
    {  // pres
       // for (int d=0; d<NUMDOF_; d+=1) {  // all
      gpltfile << " \"" << txtname << "\""
               << " using " << 1 << ":" << 2 + 4 * 8 + d + 1 << " title \"res" << d + 1 << "\""
               << " with l linetype " << d / NODDOF_ + 1 << " linewidth 2";
      gpltfile << ", \"" << txtname << "\""
               << " every 5 using " << 1 << ":" << 2 + 4 * 8 + d + 1 << " title \"res" << d + 1
               << "\""
               << " with p linetype " << d / NODDOF_ + 1 << " linewidth 2";
      gpltfile << ", \"" << txtname << "\""
               << " every 5 using " << 1 << ":" << 2 + 4 * 8 + d + 1 << ":(3*($" << 1 << "-$" << 2
               << ")):(3*$" << 2 + 2 * 4 * 8 + d + 1 << ")"
               << " title \"incres" << d + 1 << "\""
               << " with vector linetype " << d / NODDOF_ + 1 << " linewidth 1";
      gpltfile << ", \"" << txtname << "\""
               << " every 5 using " << 1 << ":" << 2 + 4 * 8 + d + 1 << ":(-3*($" << 1 << "-$" << 2
               << ")):(-3*$" << 2 + 2 * 4 * 8 + d + 1 << ")"
               << " title \"incres" << d + 1 << "\""
               << " with vector nohead linetype " << d / NODDOF_ + 1 << " linewidth 1";
      if (d + 1 < NUMDOF_) gpltfile << ", \\";
      gpltfile << std::endl;
    }

    // wait for user's activity
    gpltfile << "pause -1" << std::endl;

    gpltfile.close();
  }


  // data file
  std::ofstream txtfile;
  if (time == 0.0)
  {
    txtfile.open(txtname.c_str(), std::ios_base::out);
    txtfile << "# time"
            << " disX^1 disY^1 disZ^1 pres^1 ... disX^8 disY^8 disZ^8 pres^8"
            << " fintX^1 fintY^1 fintZ^1 incom^1 ... fintX^8 fintY^8 fintZ^8 incom^8"
            << " incfintX^1 incfintY^1 incfintZ^1 incincom^1 ... incfintX^8 incfintY^8 incfintZ^8 "
               "incincom^8"
            << std::endl;
    txtfile.close();
  }
  else
  {
    if (time != lasttime)
    {
      txtfile.open(txtname.c_str(), std::ios_base::app);
      txtfile << txtline.str();
      txtfile.close();
    }
  }

  // remove stored line of last call
  txtline.str("");

  // time
  txtline << std::scientific << time << " " << std::scientific << oldtime;

  // state
  for (std::vector<double>::iterator is = state.begin(); is != state.end(); ++is)
    txtline << std::scientific << " " << *is;

  txtline << "  ";

  // residual
  for (int ir = 0; ir < NUMDOF_; ++ir) txtline << std::scientific << " " << resid(ir);

  txtline << "  ";

  // linearised residual
  {
    LINALG::Matrix<NUMDOF_, 1> os = LINALG::Matrix<NUMDOF_, 1>(&(oldstate[0]));
    LINALG::Matrix<NUMDOF_, 1> ls = LINALG::Matrix<NUMDOF_, 1>(&(laststate[0]));
    LINALG::Matrix<NUMDOF_, 1> is;
    is = ls;
    is -= os;
    LINALG::Matrix<NUMDOF_, 1> incresid;
    incresid.MultiplyNN(lasttangent, is);

    for (int ir = 0; ir < NUMDOF_; ++ir) txtline << std::scientific << " " << incresid(ir);
  }

  // finish of line
  txtline << std::endl;

  // store last time
  lasttime = time;
  laststate = state;
  lasttangent = tangent;

  //
  return;
}

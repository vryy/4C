/*----------------------------------------------------------------------*/
/*! \file
\brief evaluation routines for sosh18 elements
\level 3


*----------------------------------------------------------------------*/

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_densematrix_eigen.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_sh18.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"

FOUR_C_NAMESPACE_OPEN


int Discret::ELEMENTS::SoSh18::init_jacobian_mapping()
{
  Core::LinAlg::Matrix<NUMNOD_SOH18, NUMDIM_SOH18> xrefe;
  for (int i = 0; i < NUMNOD_SOH18; ++i)
  {
    Core::Nodes::Node** nodes = Nodes();
    if (!nodes) FOUR_C_THROW("Nodes() returned null pointer");
    xrefe(i, 0) = Nodes()[i]->X()[0];
    xrefe(i, 1) = Nodes()[i]->X()[1];
    xrefe(i, 2) = Nodes()[i]->X()[2];
  }
  //  std::cout << "ele: " << Id() << " xrefe: " << xrefe ;
  invJ_.resize(NUMGPT_SOH18);
  detJ_.resize(NUMGPT_SOH18);


  for (int gp = 0; gp < NUMGPT_SOH18; ++gp)
  {
    // reset
    invJ_[gp].clear();
    detJ_[gp] = 0.;

    // in-plane shape functions and derivatives
    Core::LinAlg::Matrix<9, 1> shapefunct_q9;
    Core::FE::shape_function<Core::FE::CellType::quad9>(xsi_[gp], shapefunct_q9);
    Core::LinAlg::Matrix<2, 9> deriv_q9;
    Core::FE::shape_function_deriv1<Core::FE::CellType::quad9>(xsi_[gp], deriv_q9);

    for (int dim = 0; dim < 3; ++dim)
      for (int k = 0; k < 9; ++k)
      {
        invJ_[gp](0, dim) +=
            .5 * deriv_q9(0, k) * (xrefe(k + 9, dim) + xrefe(k, dim)) +
            .5 * xsi_[gp](2) * deriv_q9(0, k) * (xrefe(k + 9, dim) - xrefe(k, dim));

        invJ_[gp](1, dim) +=
            .5 * deriv_q9(1, k) * (xrefe(k + 9, dim) + xrefe(k, dim)) +
            .5 * xsi_[gp](2) * deriv_q9(1, k) * (xrefe(k + 9, dim) - xrefe(k, dim));

        invJ_[gp](2, dim) += .5 * shapefunct_q9(k) * (xrefe(k + 9, dim) - xrefe(k, dim));
      }
    detJ_[gp] = invJ_[gp].invert();
    if (detJ_[gp] < 0.) return 1;
  }

  return 0;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                          seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoSh18::nlnstiffmass(std::vector<int>& lm,  ///< location matrix
    std::vector<double>& disp,                                      ///< current displacements
    std::vector<double>& residual,                                  ///< current residual displ
    Core::LinAlg::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* stiffmatrix,  ///< element stiffness matrix
    Core::LinAlg::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* massmatrix,   ///< element mass matrix
    Core::LinAlg::Matrix<NUMDOF_SOH18, 1>* force,  ///< element internal force vector
    Core::LinAlg::Matrix<NUMGPT_SOH18, Mat::NUM_STRESS_3D>* elestress,  ///< stresses at GP
    Core::LinAlg::Matrix<NUMGPT_SOH18, Mat::NUM_STRESS_3D>* elestrain,  ///< strains at GP
    Teuchos::ParameterList& params,         ///< algorithmic parameters e.g. time
    const Inpar::STR::StressType iostress,  ///< stress output option
    const Inpar::STR::StrainType iostrain   ///< strain output option
)
{
  // get parameter interface
  set_params_interface_ptr(params);

  // update element geometry
  Core::LinAlg::Matrix<NUMNOD_SOH18, NUMDIM_SOH18> xrefe;  // reference coord. of element
  Core::LinAlg::Matrix<NUMNOD_SOH18, NUMDIM_SOH18> xcurr;  // current  coord. of element

  Core::Nodes::Node** nodes = Nodes();
  for (int i = 0; i < NUMNOD_SOH18; ++i)
  {
    const auto& x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];

    xcurr(i, 0) = xrefe(i, 0) + disp[i * NODDOF_SOH18 + 0];
    xcurr(i, 1) = xrefe(i, 1) + disp[i * NODDOF_SOH18 + 1];
    xcurr(i, 2) = xrefe(i, 2) + disp[i * NODDOF_SOH18 + 2];
  }

  // we need the (residual) displacement at the previous step
  Core::LinAlg::Matrix<NUMDOF_SOH18, 1> res_d;
  for (int i = 0; i < NUMDOF_SOH18; ++i) res_d(i) = residual[i];

  // EAS stuff
  std::vector<Core::LinAlg::Matrix<6, num_eas>> M_gp(num_eas);
  Core::LinAlg::Matrix<3, 1> G3_0_contra;
  Core::LinAlg::Matrix<6, num_eas> M;
  if (eas_)
  {
    // recover EAS **************************************
    if (not IsParamsInterface())
      if (stiffmatrix)
      {
        feas_.multiply(1., Kad_, res_d, 1.);
        alpha_eas_inc_.multiply(-1., KaaInv_, feas_, 0.);
        alpha_eas_.update(1., alpha_eas_inc_, 1.);
      }
    // recover EAS **************************************

    // prepare EAS***************************************
    eas_setup(M_gp, G3_0_contra, xrefe);
    feas_.clear();
    KaaInv_.clear();
    Kad_.clear();
    // prepare EAS***************************************
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp = 0; gp < NUMGPT_SOH18; ++gp)
  {
    // in-plane shape functions and derivatives
    Core::LinAlg::Matrix<9, 1> shapefunct_q9;
    Core::FE::shape_function<Core::FE::CellType::quad9>(xsi_[gp], shapefunct_q9);
    Core::LinAlg::Matrix<2, 9> deriv_q9;
    Core::FE::shape_function_deriv1<Core::FE::CellType::quad9>(xsi_[gp], deriv_q9);

    /* get the inverse of the Jacobian matrix which looks like:
    **         [ x_,r  y_,r  z_,r ]
    **     J = [ x_,s  y_,s  z_,s ]
    **         [ x_,t  y_,t  z_,t ]
    */
    // compute the Jacobian shell-style (G^T)
    Core::LinAlg::Matrix<NUMDIM_SOH18, NUMDIM_SOH18> jac;
    for (int dim = 0; dim < 3; ++dim)
      for (int k = 0; k < 9; ++k)
      {
        jac(0, dim) += .5 * deriv_q9(0, k) * (xrefe(k + 9, dim) + xrefe(k, dim)) +
                       .5 * xsi_[gp](2) * deriv_q9(0, k) * (xrefe(k + 9, dim) - xrefe(k, dim));

        jac(1, dim) += .5 * deriv_q9(1, k) * (xrefe(k + 9, dim) + xrefe(k, dim)) +
                       .5 * xsi_[gp](2) * deriv_q9(1, k) * (xrefe(k + 9, dim) - xrefe(k, dim));

        jac(2, dim) += .5 * shapefunct_q9(k) * (xrefe(k + 9, dim) - xrefe(k, dim));
      }
    double detJ = jac.determinant();

    // transformation from local (parameter) element space to global(material) space
    // with famous 'T'-matrix already used for EAS but now evaluated at each gp
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, Mat::NUM_STRESS_3D> TinvT;
    evaluate_t(jac, TinvT);

    // **********************************************************************
    // set up B-Operator in local(parameter) element space including ANS
    // **********************************************************************
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, NUMDOF_SOH18> bop_loc(true);
    calculate_bop_loc(xcurr, xrefe, shapefunct_q9, deriv_q9, gp, bop_loc);
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, NUMDOF_SOH18> bop;
    bop.multiply(TinvT, bop_loc);

    // **************************************************************************
    // shell-like calculation of strains
    // see Diss. Koschnik page 41
    // **************************************************************************
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1> lstrain(true);
    calculate_loc_strain(xcurr, xrefe, shapefunct_q9, deriv_q9, gp, lstrain);
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1> glstrain;
    glstrain.multiply(TinvT, lstrain);
    // **************************************************************************
    // shell-like calculation of strains
    // **************************************************************************

    // EAS: enhance the strains ***********************************************
    if (eas_)
    {
      double t33 = 0.;
      for (int dim = 0; dim < 3; ++dim) t33 += jac(2, dim) * G3_0_contra(dim);

      M.multiply(t33 * t33 / detJ, TinvT, M_gp[gp], 0.);
      glstrain.multiply(1., M, alpha_eas_, 1.);
    }
    // end EAS: enhance the strains *******************************************

    // calculate the deformation gradient consistent to the modified strains
    // but only if the material needs a deformation gradient (e.g. plasticity)
    Core::LinAlg::Matrix<NUMDIM_SOH18, NUMDIM_SOH18> defgrd;
    if (Teuchos::rcp_static_cast<Mat::So3Material>(Material())->needs_defgrd() ||
        iostrain == Inpar::STR::strain_ea || iostress == Inpar::STR::stress_cauchy)
    {
      // compute the deformation gradient - shell-style
      // deformation gradient with derivatives w.r.t. local basis
      Core::LinAlg::Matrix<NUMDIM_SOH18, NUMDIM_SOH18> defgrd_loc(true);
      for (int k = 0; k < 9; ++k)
        for (int dim = 0; dim < NUMDIM_SOH18; ++dim)
        {
          defgrd_loc(dim, 0) += .5 * deriv_q9(0, k) *
                                ((xcurr(k + 9, dim) + xcurr(k, dim)) +
                                    xsi_[gp](2) * (xcurr(k + 9, dim) - xcurr(k, dim)));
          defgrd_loc(dim, 1) += .5 * deriv_q9(1, k) *
                                ((xcurr(k + 9, dim) + xcurr(k, dim)) +
                                    xsi_[gp](2) * (xcurr(k + 9, dim) - xcurr(k, dim)));
          defgrd_loc(dim, 2) += .5 * shapefunct_q9(k) * (xcurr(k + 9, dim) - xcurr(k, dim));
        }

      // displacement-based deformation gradient
      Core::LinAlg::Matrix<NUMDIM_SOH18, NUMDIM_SOH18> defgrd_disp;
      defgrd_disp.multiply_nt(defgrd_loc, invJ_[gp]);
      calc_consistent_defgrd(defgrd_disp, glstrain, defgrd);
    }

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix must be retrieved,
    ** all necessary data must be passed.
    */
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, Mat::NUM_STRESS_3D> cmat(true);
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1> stress(true);

    Teuchos::RCP<Mat::So3Material> so3mat = Teuchos::rcp_static_cast<Mat::So3Material>(Material());
    so3mat->evaluate(&defgrd, &glstrain, params, &stress, &cmat, gp, Id());
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // strain output **********************************************************
    if (elestrain)
    {
      // return gp strains if necessary
      switch (iostrain)
      {
        case Inpar::STR::strain_gl:
        {
          if (elestrain == nullptr) FOUR_C_THROW("strain data not available");
          for (int i = 0; i < 3; ++i)
          {
            (*elestrain)(gp, i) = glstrain(i);
          }
          for (int i = 3; i < 6; ++i)
          {
            (*elestrain)(gp, i) = 0.5 * glstrain(i);
          }
        }
        break;
        case Inpar::STR::strain_ea:
        {
          Core::LinAlg::Matrix<3, 3> bi;
          bi.multiply_nt(defgrd, defgrd);
          bi.invert();
          for (int i = 0; i < 3; i++) (*elestrain)(gp, i) = .5 * (1. - bi(i, i));
          (*elestrain)(gp, 3) = -bi(0, 1);
          (*elestrain)(gp, 4) = -bi(2, 1);
          (*elestrain)(gp, 5) = -bi(0, 2);
          break;
        }
        case Inpar::STR::strain_none:
          break;
        default:
          FOUR_C_THROW("requested strain option not available");
          break;
      }
    }
    // end of strain output ***************************************************

    // stress output **********************************************************
    if (elestress)
    {
      // return gp strains if necessary
      switch (iostress)
      {
        case Inpar::STR::stress_2pk:
        {
          if (elestress == nullptr) FOUR_C_THROW("stress data not available");
          for (int i = 0; i < Mat::NUM_STRESS_3D; ++i)
          {
            (*elestress)(gp, i) = stress(i);
          }
        }
        break;
        case Inpar::STR::stress_cauchy:
        {
          if (elestress == nullptr) FOUR_C_THROW("stress data not available");
          Core::LinAlg::Matrix<3, 3> pkstress;
          pkstress(0, 0) = stress(0);
          pkstress(0, 1) = stress(3);
          pkstress(0, 2) = stress(5);
          pkstress(1, 0) = pkstress(0, 1);
          pkstress(1, 1) = stress(1);
          pkstress(1, 2) = stress(4);
          pkstress(2, 0) = pkstress(0, 2);
          pkstress(2, 1) = pkstress(1, 2);
          pkstress(2, 2) = stress(2);

          Core::LinAlg::Matrix<3, 3> cauchystress;
          Core::LinAlg::Matrix<3, 3> temp;
          temp.multiply(1.0 / defgrd.determinant(), defgrd, pkstress);
          cauchystress.multiply_nt(temp, defgrd);

          (*elestress)(gp, 0) = cauchystress(0, 0);
          (*elestress)(gp, 1) = cauchystress(1, 1);
          (*elestress)(gp, 2) = cauchystress(2, 2);
          (*elestress)(gp, 3) = cauchystress(0, 1);
          (*elestress)(gp, 4) = cauchystress(1, 2);
          (*elestress)(gp, 5) = cauchystress(0, 2);
        }
        break;
        case Inpar::STR::stress_none:
          break;
        default:
          FOUR_C_THROW("requested stress option not available");
          break;
      }
    }
    // end of stress output ***************************************************

    double detJ_w = detJ * wgt_[gp];
    // update internal force vector
    if (force != nullptr) force->multiply_tn(detJ_w, bop, stress, 1.0);

    // update stiffness matrix
    if (stiffmatrix != nullptr)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, NUMDOF_SOH18> cb;
      cb.multiply(cmat, bop);
      stiffmatrix->multiply_tn(detJ_w, bop, cb, 1.0);  // standard hex8 evaluation
      // intergrate `geometric' stiffness matrix and add to keu *****************
      calculate_geo_stiff(shapefunct_q9, deriv_q9, TinvT, gp, detJ_w, stress, stiffmatrix);

      // EAS technology: integrate matrices --------------------------------- EAS
      if (eas_)
      {
        Core::LinAlg::Matrix<6, num_eas> cM;
        cM.multiply(cmat, M);
        KaaInv_.multiply_tn(detJ_w, M, cM, 1.);
        Kad_.multiply_tn(detJ_w, M, cb, 1.);
        feas_.multiply_tn(detJ_w, M, stress, 1.);
      }
      // EAS technology: integrate matrices --------------------------------- EAS
    }

    if (massmatrix != nullptr)  // evaluate mass matrix +++++++++++++++++++++++++
    {
      // shape function and derivatives
      Core::LinAlg::Matrix<NUMNOD_SOH18, 1> shapefunct;
      Core::FE::shape_function<Core::FE::CellType::hex18>(xsi_[gp], shapefunct);

      double density = Material()->Density(gp);

      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      double ifactor, massfactor;
      for (int inod = 0; inod < NUMNOD_SOH18; ++inod)
      {
        ifactor = shapefunct(inod) * factor;
        for (int jnod = 0; jnod < NUMNOD_SOH18; ++jnod)
        {
          massfactor = shapefunct(jnod) * ifactor;  // intermediate factor
          (*massmatrix)(NUMDIM_SOH18 * inod + 0, NUMDIM_SOH18 * jnod + 0) += massfactor;
          (*massmatrix)(NUMDIM_SOH18 * inod + 1, NUMDIM_SOH18 * jnod + 1) += massfactor;
          (*massmatrix)(NUMDIM_SOH18 * inod + 2, NUMDIM_SOH18 * jnod + 2) += massfactor;
        }
      }
    }  // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++

    /* =========================================================================*/
  } /* ==================================================== end of Loop over GP */
  /* =========================================================================*/

  if (stiffmatrix && eas_)
  {
    Core::LinAlg::FixedSizeSerialDenseSolver<num_eas, num_eas, 1> solve_for_KaaInv;
    solve_for_KaaInv.SetMatrix(KaaInv_);
    int err2 = solve_for_KaaInv.Factor();
    int err = solve_for_KaaInv.invert();
    if ((err != 0) || (err2 != 0)) FOUR_C_THROW("Inversion of Kaa failed");

    Core::LinAlg::Matrix<NUMDOF_SOH18, num_eas> KdaKaa;
    KdaKaa.multiply_tn(Kad_, KaaInv_);
    stiffmatrix->multiply(-1., KdaKaa, Kad_, 1.);
    force->multiply(-1., KdaKaa, feas_, 1.);
  }

  return;
}  // Discret::ELEMENTS::So_hex8::nlnstiffmass


/*----------------------------------------------------------------------*
 |  lump mass matrix (private)                              seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoSh18::soh18_lumpmass(
    Core::LinAlg::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* emass)
{
  // lump mass matrix
  if (emass != nullptr)
  {
    // we assume #elemat2 is a square matrix
    for (unsigned int c = 0; c < (*emass).numCols(); ++c)  // parse columns
    {
      double d = 0.0;
      for (unsigned int r = 0; r < (*emass).numRows(); ++r)  // parse rows
      {
        d += (*emass)(r, c);  // accumulate row entries
        (*emass)(r, c) = 0.0;
      }
      (*emass)(c, c) = d;  // apply sum of row entries on diagonal
    }
  }
}

/*----------------------------------------------------------------------*
 |  init the element (public)                               seitz 11/14 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoSh18Type::initialize(Core::FE::Discretization& dis)
{
  // here we order the nodes such that we have a positive definite jacobian
  //       maybe the python script generating the hex18 elements would be a better place for this.
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<Discret::ELEMENTS::SoSh18*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex18* failed");
    if (actele->init_jacobian_mapping() == 1) actele->flip_t();
  }
  dis.fill_complete(false, false, false);

  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<Discret::ELEMENTS::SoSh18*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_hex18* failed");
    if (actele->init_jacobian_mapping() == 1) FOUR_C_THROW("why");
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  revert the 3rd parameter direction                      seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoSh18::flip_t()
{
  if (NodeIds() == nullptr) FOUR_C_THROW("couldn't get node ids");
  // reorder nodes
  int new_nodeids[NUMNOD_SOH18];
  new_nodeids[0] = NodeIds()[9];
  new_nodeids[1] = NodeIds()[10];
  new_nodeids[2] = NodeIds()[11];
  new_nodeids[3] = NodeIds()[12];
  new_nodeids[4] = NodeIds()[13];
  new_nodeids[5] = NodeIds()[14];
  new_nodeids[6] = NodeIds()[15];
  new_nodeids[7] = NodeIds()[16];
  new_nodeids[8] = NodeIds()[17];

  new_nodeids[9] = NodeIds()[0];
  new_nodeids[10] = NodeIds()[1];
  new_nodeids[11] = NodeIds()[2];
  new_nodeids[12] = NodeIds()[3];
  new_nodeids[13] = NodeIds()[4];
  new_nodeids[14] = NodeIds()[5];
  new_nodeids[15] = NodeIds()[6];
  new_nodeids[16] = NodeIds()[7];
  new_nodeids[17] = NodeIds()[8];

  SetNodeIds(NUMNOD_SOH18, new_nodeids);
  return;
}


void Discret::ELEMENTS::SoSh18::evaluate_t(
    const Core::LinAlg::Matrix<NUMDIM_SOH18, NUMDIM_SOH18>& jac,
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, Mat::NUM_STRESS_3D>& TinvT)
{
  // build T^T transformation matrix which maps
  // between global (r,s,t)-coordinates and local (x,y,z)-coords
  // later, invert the transposed to map from local to global
  // see literature for details (e.g. Andelfinger)
  // it is based on the voigt notation for strains: xx,yy,zz,xy,yz,xz
  TinvT(0, 0) = jac(0, 0) * jac(0, 0);
  TinvT(1, 0) = jac(1, 0) * jac(1, 0);
  TinvT(2, 0) = jac(2, 0) * jac(2, 0);
  TinvT(3, 0) = 2 * jac(0, 0) * jac(1, 0);
  TinvT(4, 0) = 2 * jac(1, 0) * jac(2, 0);
  TinvT(5, 0) = 2 * jac(0, 0) * jac(2, 0);

  TinvT(0, 1) = jac(0, 1) * jac(0, 1);
  TinvT(1, 1) = jac(1, 1) * jac(1, 1);
  TinvT(2, 1) = jac(2, 1) * jac(2, 1);
  TinvT(3, 1) = 2 * jac(0, 1) * jac(1, 1);
  TinvT(4, 1) = 2 * jac(1, 1) * jac(2, 1);
  TinvT(5, 1) = 2 * jac(0, 1) * jac(2, 1);

  TinvT(0, 2) = jac(0, 2) * jac(0, 2);
  TinvT(1, 2) = jac(1, 2) * jac(1, 2);
  TinvT(2, 2) = jac(2, 2) * jac(2, 2);
  TinvT(3, 2) = 2 * jac(0, 2) * jac(1, 2);
  TinvT(4, 2) = 2 * jac(1, 2) * jac(2, 2);
  TinvT(5, 2) = 2 * jac(0, 2) * jac(2, 2);

  TinvT(0, 3) = jac(0, 0) * jac(0, 1);
  TinvT(1, 3) = jac(1, 0) * jac(1, 1);
  TinvT(2, 3) = jac(2, 0) * jac(2, 1);
  TinvT(3, 3) = jac(0, 0) * jac(1, 1) + jac(1, 0) * jac(0, 1);
  TinvT(4, 3) = jac(1, 0) * jac(2, 1) + jac(2, 0) * jac(1, 1);
  TinvT(5, 3) = jac(0, 0) * jac(2, 1) + jac(2, 0) * jac(0, 1);


  TinvT(0, 4) = jac(0, 1) * jac(0, 2);
  TinvT(1, 4) = jac(1, 1) * jac(1, 2);
  TinvT(2, 4) = jac(2, 1) * jac(2, 2);
  TinvT(3, 4) = jac(0, 1) * jac(1, 2) + jac(1, 1) * jac(0, 2);
  TinvT(4, 4) = jac(1, 1) * jac(2, 2) + jac(2, 1) * jac(1, 2);
  TinvT(5, 4) = jac(0, 1) * jac(2, 2) + jac(2, 1) * jac(0, 2);

  TinvT(0, 5) = jac(0, 0) * jac(0, 2);
  TinvT(1, 5) = jac(1, 0) * jac(1, 2);
  TinvT(2, 5) = jac(2, 0) * jac(2, 2);
  TinvT(3, 5) = jac(0, 0) * jac(1, 2) + jac(1, 0) * jac(0, 2);
  TinvT(4, 5) = jac(1, 0) * jac(2, 2) + jac(2, 0) * jac(1, 2);
  TinvT(5, 5) = jac(0, 0) * jac(2, 2) + jac(2, 0) * jac(0, 2);

  // now evaluate T^{-T} with solver
  Core::LinAlg::FixedSizeSerialDenseSolver<Mat::NUM_STRESS_3D, Mat::NUM_STRESS_3D, 1>
      solve_for_inverseT;
  solve_for_inverseT.SetMatrix(TinvT);
  int err2 = solve_for_inverseT.Factor();
  int err = solve_for_inverseT.invert();
  if ((err != 0) && (err2 != 0)) FOUR_C_THROW("Inversion of Tinv (Jacobian) failed");
  return;
}

void Discret::ELEMENTS::SoSh18::calculate_loc_strain(
    const Core::LinAlg::Matrix<NUMNOD_SOH18, NUMDIM_SOH18>& xcurr,
    const Core::LinAlg::Matrix<NUMNOD_SOH18, NUMDIM_SOH18>& xrefe,
    const Core::LinAlg::Matrix<9, 1>& shape_q9, const Core::LinAlg::Matrix<2, 9>& deriv_q9,
    const int gp, Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1>& lstrain)
{
  for (int dim = 0; dim < 3; ++dim)
    for (int k = 0; k < 9; ++k)
      for (int l = 0; l < 9; ++l)
      {
        if (dsg_membrane_)
        {
          // constant normal strain
          lstrain(0) += .125 * dsg_membrane_r_[gp % 9](k, l) *
                        ((xcurr(k + 9, dim) + xcurr(k, dim)) *
                                (xcurr(l + 9, dim) + xcurr(l, dim))  // a_alphaalpha
                            - (xrefe(k + 9, dim) + xrefe(k, dim)) *
                                  (xrefe(l + 9, dim) + xrefe(l, dim))  // A_alphaalpha
                        );
          lstrain(1) += .125 * dsg_membrane_s_[gp % 9](k, l) *
                        ((xcurr(k + 9, dim) + xcurr(k, dim)) *
                                (xcurr(l + 9, dim) + xcurr(l, dim))  // a_alphaalpha
                            - (xrefe(k + 9, dim) + xrefe(k, dim)) *
                                  (xrefe(l + 9, dim) + xrefe(l, dim))  // A_alphaalpha
                        );
          // constant in-plane shear strain

          lstrain(3) +=
              .25 * dsg_membrane_rs_[gp % 9](k, l) *
              ((xcurr(k + 9, dim) + xcurr(k, dim)) * (xcurr(l + 9, dim) + xcurr(l, dim))  // a_01
                  -
                  (xrefe(k + 9, dim) + xrefe(k, dim)) * (xrefe(l + 9, dim) + xrefe(l, dim))  // A_01
              );
        }
        else
        {
          // constant normal strain
          for (int alpha = 0; alpha < 2; ++alpha)
            lstrain(alpha) += .125 * deriv_q9(alpha, k) * deriv_q9(alpha, l) *
                              ((xcurr(k + 9, dim) + xcurr(k, dim)) *
                                      (xcurr(l + 9, dim) + xcurr(l, dim))  // a_alphaalpha
                                  - (xrefe(k + 9, dim) + xrefe(k, dim)) *
                                        (xrefe(l + 9, dim) + xrefe(l, dim))  // A_alphaalpha
                              );

          // constant in-plane shear strain
          lstrain(3) +=
              .25 * deriv_q9(0, k) * deriv_q9(1, l) *
              ((xcurr(k + 9, dim) + xcurr(k, dim)) * (xcurr(l + 9, dim) + xcurr(l, dim))  // a_01
                  -
                  (xrefe(k + 9, dim) + xrefe(k, dim)) * (xrefe(l + 9, dim) + xrefe(l, dim))  // A_01
              );
        }

        // linear normal strain
        for (int alpha = 0; alpha < 2; ++alpha)
          lstrain(alpha) +=
              xsi_[gp](2) * .125 * deriv_q9(alpha, k) * deriv_q9(alpha, l) *
              ((xcurr(k + 9, dim) + xcurr(k, dim)) * (xcurr(l + 9, dim) - xcurr(l, dim)) +
                  (xcurr(k + 9, dim) - xcurr(k, dim)) *
                      (xcurr(l + 9, dim) + xcurr(l, dim))  // b_alphaalpha
                  - (xrefe(k + 9, dim) + xrefe(k, dim)) * (xrefe(l + 9, dim) - xrefe(l, dim)) -
                  (xrefe(k + 9, dim) - xrefe(k, dim)) *
                      (xrefe(l + 9, dim) + xrefe(l, dim))  // B_alphaalpha
              );

        // linear in-plane shear strain
        lstrain(3) +=
            xsi_[gp](2) * .25 * deriv_q9(0, k) * deriv_q9(1, l) *
            ((xcurr(k + 9, dim) + xcurr(k, dim)) * (xcurr(l + 9, dim) - xcurr(l, dim)) +
                (xcurr(k + 9, dim) - xcurr(k, dim)) * (xcurr(l + 9, dim) + xcurr(l, dim))  // b_01
                - (xrefe(k + 9, dim) + xrefe(k, dim)) * (xrefe(l + 9, dim) - xrefe(l, dim)) -
                (xrefe(k + 9, dim) - xrefe(k, dim)) * (xrefe(l + 9, dim) + xrefe(l, dim))  // B_01
            );

        if (dsg_shear_)
        {
          // constant transverse shear strain
          lstrain(4) +=
              .25 * dsg_shear_s_[gp % 9](k, l) *
              ((xcurr(k + 9, dim) + xcurr(k, dim)) * (xcurr(l + 9, dim) - xcurr(l, dim)) -
                  (xrefe(k + 9, dim) + xrefe(k, dim)) * (xrefe(l + 9, dim) - xrefe(l, dim)));
          lstrain(5) +=
              .25 * dsg_shear_r_[gp % 9](k, l) *
              ((xcurr(k + 9, dim) + xcurr(k, dim)) * (xcurr(l + 9, dim) - xcurr(l, dim)) -
                  (xrefe(k + 9, dim) + xrefe(k, dim)) * (xrefe(l + 9, dim) - xrefe(l, dim)));
        }
        else
        {
          // constant transverse shear strain
          lstrain(4) +=
              .25 * deriv_q9(1, k) * shape_q9(l) *
              ((xcurr(k + 9, dim) + xcurr(k, dim)) * (xcurr(l + 9, dim) - xcurr(l, dim))  // a_12
                  -
                  (xrefe(k + 9, dim) + xrefe(k, dim)) * (xrefe(l + 9, dim) - xrefe(l, dim))  // A_12
              );
          lstrain(5) +=
              .25 * deriv_q9(0, k) * shape_q9(l) *
              ((xcurr(k + 9, dim) + xcurr(k, dim)) * (xcurr(l + 9, dim) - xcurr(l, dim))  // a_02
                  -
                  (xrefe(k + 9, dim) + xrefe(k, dim)) * (xrefe(l + 9, dim) - xrefe(l, dim))  // A_02
              );
        }

        // linear transverse shear strain
        lstrain(4) += xsi_[gp](2) * .25 * deriv_q9(1, k) * shape_q9(l) *
                      ((xcurr(k + 9, dim) - xcurr(k, dim)) *
                              (xcurr(l + 9, dim) - xcurr(l, dim))  // a_2 dot a_2,1
                          - (xrefe(k + 9, dim) - xrefe(k, dim)) *
                                (xrefe(l + 9, dim) - xrefe(l, dim))  // A_2 dot A_2,1
                      );
        lstrain(5) += xsi_[gp](2) * .25 * deriv_q9(0, k) * shape_q9(l) *
                      ((xcurr(k + 9, dim) - xcurr(k, dim)) *
                              (xcurr(l + 9, dim) - xcurr(l, dim))  // a_2 dot a_2,0
                          - (xrefe(k + 9, dim) - xrefe(k, dim)) *
                                (xrefe(l + 9, dim) - xrefe(l, dim))  // A_2 dot A_2,0
                      );

        // transverse normal strain
        if (dsg_ctl_)
        {
          lstrain(2) += .125 * dsg_transverse_t_[gp % 9](k, l) *
                        ((xcurr(k + 9, dim) - xcurr(k, dim)) *
                                (xcurr(l + 9, dim) - xcurr(l, dim))  // a_2 dot a_2
                            - (xrefe(k + 9, dim) - xrefe(k, dim)) *
                                  (xrefe(l + 9, dim) - xrefe(l, dim))  // A_2 dot A_2
                        );
        }
        else
        {
          lstrain(2) += .125 * shape_q9(k) * shape_q9(l) *
                        ((xcurr(k + 9, dim) - xcurr(k, dim)) *
                                (xcurr(l + 9, dim) - xcurr(l, dim))  // a_2 dot a_2
                            - (xrefe(k + 9, dim) - xrefe(k, dim)) *
                                  (xrefe(l + 9, dim) - xrefe(l, dim))  // A_2 dot A_2
                        );
        }
      }  // k=0..8; l=0..8; dim=0..2

  return;
}

void Discret::ELEMENTS::SoSh18::calculate_bop_loc(
    const Core::LinAlg::Matrix<NUMNOD_SOH18, NUMDIM_SOH18>& xcurr,
    const Core::LinAlg::Matrix<NUMNOD_SOH18, NUMDIM_SOH18>& xrefe,
    const Core::LinAlg::Matrix<9, 1>& shape_q9, const Core::LinAlg::Matrix<2, 9>& deriv_q9,
    const int gp, Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, NUMDOF_SOH18>& bop_loc)
{
  for (int dim = 0; dim < NUMDIM_SOH18; ++dim)
    for (int k = 0; k < 9; ++k)
      for (int l = 0; l < 9; ++l)
      {
        if (dsg_membrane_)
        {
          // constant normal strain
          bop_loc(0, (k + 9) * 3 + dim) +=
              .125 * dsg_membrane_r_[gp % 9](k, l) * (xcurr(l + 9, dim) + xcurr(l, dim));
          bop_loc(0, k * 3 + dim) +=
              .125 * dsg_membrane_r_[gp % 9](k, l) * (xcurr(l + 9, dim) + xcurr(l, dim));
          bop_loc(0, (l + 9) * 3 + dim) +=
              .125 * dsg_membrane_r_[gp % 9](k, l) * (xcurr(k + 9, dim) + xcurr(k, dim));
          bop_loc(0, l * 3 + dim) +=
              .125 * dsg_membrane_r_[gp % 9](k, l) * (xcurr(k + 9, dim) + xcurr(k, dim));

          bop_loc(1, (k + 9) * 3 + dim) +=
              .125 * dsg_membrane_s_[gp % 9](k, l) * (xcurr(l + 9, dim) + xcurr(l, dim));
          bop_loc(1, k * 3 + dim) +=
              .125 * dsg_membrane_s_[gp % 9](k, l) * (xcurr(l + 9, dim) + xcurr(l, dim));
          bop_loc(1, (l + 9) * 3 + dim) +=
              .125 * dsg_membrane_s_[gp % 9](k, l) * (xcurr(k + 9, dim) + xcurr(k, dim));
          bop_loc(1, l * 3 + dim) +=
              .125 * dsg_membrane_s_[gp % 9](k, l) * (xcurr(k + 9, dim) + xcurr(k, dim));

          // constant in-plane shear strain
          bop_loc(3, (k + 9) * 3 + dim) +=
              .25 * dsg_membrane_rs_[gp % 9](k, l) * (xcurr(l + 9, dim) + xcurr(l, dim));
          bop_loc(3, k * 3 + dim) +=
              .25 * dsg_membrane_rs_[gp % 9](k, l) * (xcurr(l + 9, dim) + xcurr(l, dim));
          bop_loc(3, (l + 9) * 3 + dim) +=
              .25 * dsg_membrane_rs_[gp % 9](k, l) * (xcurr(k + 9, dim) + xcurr(k, dim));
          bop_loc(3, l * 3 + dim) +=
              .25 * dsg_membrane_rs_[gp % 9](k, l) * (xcurr(k + 9, dim) + xcurr(k, dim));
        }
        else
        {
          // constant normal strain
          for (int alpha = 0; alpha < 2; ++alpha)
          {
            bop_loc(alpha, (k + 9) * 3 + dim) += .125 * deriv_q9(alpha, k) * deriv_q9(alpha, l) *
                                                 (xcurr(l + 9, dim) + xcurr(l, dim));
            bop_loc(alpha, k * 3 + dim) += .125 * deriv_q9(alpha, k) * deriv_q9(alpha, l) *
                                           (xcurr(l + 9, dim) + xcurr(l, dim));
            bop_loc(alpha, (l + 9) * 3 + dim) += .125 * deriv_q9(alpha, k) * deriv_q9(alpha, l) *
                                                 (xcurr(k + 9, dim) + xcurr(k, dim));
            bop_loc(alpha, l * 3 + dim) += .125 * deriv_q9(alpha, k) * deriv_q9(alpha, l) *
                                           (xcurr(k + 9, dim) + xcurr(k, dim));
          }

          // constant in-plane shear strain
          bop_loc(3, (k + 9) * 3 + dim) +=
              .25 * deriv_q9(0, k) * deriv_q9(1, l) * (xcurr(l + 9, dim) + xcurr(l, dim));
          bop_loc(3, k * 3 + dim) +=
              .25 * deriv_q9(0, k) * deriv_q9(1, l) * (xcurr(l + 9, dim) + xcurr(l, dim));
          bop_loc(3, (l + 9) * 3 + dim) +=
              .25 * deriv_q9(0, k) * deriv_q9(1, l) * (xcurr(k + 9, dim) + xcurr(k, dim));
          bop_loc(3, l * 3 + dim) +=
              .25 * deriv_q9(0, k) * deriv_q9(1, l) * (xcurr(k + 9, dim) + xcurr(k, dim));
        }

        // linear normal strain
        for (int alpha = 0; alpha < 2; ++alpha)
        {
          bop_loc(alpha, (k + 9) * 3 + dim) +=
              xsi_[gp](2) * .25 * deriv_q9(alpha, k) * deriv_q9(alpha, l) * xcurr(l + 9, dim);
          bop_loc(alpha, k * 3 + dim) +=
              -xsi_[gp](2) * .25 * deriv_q9(alpha, k) * deriv_q9(alpha, l) * xcurr(l, dim);
          bop_loc(alpha, (l + 9) * 3 + dim) +=
              xsi_[gp](2) * .25 * deriv_q9(alpha, k) * deriv_q9(alpha, l) * xcurr(k + 9, dim);
          bop_loc(alpha, l * 3 + dim) +=
              -xsi_[gp](2) * .25 * deriv_q9(alpha, k) * deriv_q9(alpha, l) * xcurr(k, dim);
        }

        // linear in-plane shear strain
        bop_loc(3, (k + 9) * 3 + dim) +=
            xsi_[gp](2) * .5 * deriv_q9(0, k) * deriv_q9(1, l) * xcurr(l + 9, dim);
        bop_loc(3, k * 3 + dim) +=
            -xsi_[gp](2) * .5 * deriv_q9(0, k) * deriv_q9(1, l) * xcurr(l, dim);
        bop_loc(3, (l + 9) * 3 + dim) +=
            xsi_[gp](2) * .5 * deriv_q9(0, k) * deriv_q9(1, l) * xcurr(k + 9, dim);
        bop_loc(3, l * 3 + dim) +=
            -xsi_[gp](2) * .5 * deriv_q9(0, k) * deriv_q9(1, l) * xcurr(k, dim);

        if (dsg_shear_)
        {
          // constant transverse shear strain
          bop_loc(4, (k + 9) * 3 + dim) +=
              .25 * dsg_shear_s_[gp % 9](k, l) * (xcurr(l + 9, dim) - xcurr(l, dim));
          bop_loc(4, k * 3 + dim) +=
              .25 * dsg_shear_s_[gp % 9](k, l) * (xcurr(l + 9, dim) - xcurr(l, dim));
          bop_loc(4, (l + 9) * 3 + dim) +=
              .25 * dsg_shear_s_[gp % 9](k, l) * (xcurr(k + 9, dim) + xcurr(k, dim));
          bop_loc(4, l * 3 + dim) +=
              -.25 * dsg_shear_s_[gp % 9](k, l) * (xcurr(k + 9, dim) + xcurr(k, dim));

          bop_loc(5, (k + 9) * 3 + dim) +=
              .25 * dsg_shear_r_[gp % 9](k, l) * (xcurr(l + 9, dim) - xcurr(l, dim));
          bop_loc(5, k * 3 + dim) +=
              .25 * dsg_shear_r_[gp % 9](k, l) * (xcurr(l + 9, dim) - xcurr(l, dim));
          bop_loc(5, (l + 9) * 3 + dim) +=
              .25 * dsg_shear_r_[gp % 9](k, l) * (xcurr(k + 9, dim) + xcurr(k, dim));
          bop_loc(5, l * 3 + dim) +=
              -.25 * dsg_shear_r_[gp % 9](k, l) * (xcurr(k + 9, dim) + xcurr(k, dim));
        }
        else
        {
          // constant transverse shear strain
          bop_loc(4, (k + 9) * 3 + dim) +=
              .25 * deriv_q9(1, k) * shape_q9(l) * (xcurr(l + 9, dim) - xcurr(l, dim));
          bop_loc(4, k * 3 + dim) +=
              .25 * deriv_q9(1, k) * shape_q9(l) * (xcurr(l + 9, dim) - xcurr(l, dim));
          bop_loc(4, (l + 9) * 3 + dim) +=
              .25 * deriv_q9(1, k) * shape_q9(l) * (xcurr(k + 9, dim) + xcurr(k, dim));
          bop_loc(4, l * 3 + dim) +=
              -.25 * deriv_q9(1, k) * shape_q9(l) * (xcurr(k + 9, dim) + xcurr(k, dim));
          bop_loc(5, (k + 9) * 3 + dim) +=
              .25 * deriv_q9(0, k) * shape_q9(l) * (xcurr(l + 9, dim) - xcurr(l, dim));
          bop_loc(5, k * 3 + dim) +=
              .25 * deriv_q9(0, k) * shape_q9(l) * (xcurr(l + 9, dim) - xcurr(l, dim));
          bop_loc(5, (l + 9) * 3 + dim) +=
              .25 * deriv_q9(0, k) * shape_q9(l) * (xcurr(k + 9, dim) + xcurr(k, dim));
          bop_loc(5, l * 3 + dim) +=
              -.25 * deriv_q9(0, k) * shape_q9(l) * (xcurr(k + 9, dim) + xcurr(k, dim));
        }

        // linear transverse shear strain
        bop_loc(4, (k + 9) * 3 + dim) +=
            xsi_[gp](2) * .25 * deriv_q9(1, k) * shape_q9(l) * (xcurr(l + 9, dim) - xcurr(l, dim));
        bop_loc(4, k * 3 + dim) +=
            -xsi_[gp](2) * .25 * deriv_q9(1, k) * shape_q9(l) * (xcurr(l + 9, dim) - xcurr(l, dim));
        bop_loc(4, (l + 9) * 3 + dim) +=
            xsi_[gp](2) * .25 * deriv_q9(1, k) * shape_q9(l) * (xcurr(k + 9, dim) - xcurr(k, dim));
        bop_loc(4, l * 3 + dim) +=
            -xsi_[gp](2) * .25 * deriv_q9(1, k) * shape_q9(l) * (xcurr(k + 9, dim) - xcurr(k, dim));
        bop_loc(5, (k + 9) * 3 + dim) +=
            xsi_[gp](2) * .25 * deriv_q9(0, k) * shape_q9(l) * (xcurr(l + 9, dim) - xcurr(l, dim));
        bop_loc(5, k * 3 + dim) +=
            -xsi_[gp](2) * .25 * deriv_q9(0, k) * shape_q9(l) * (xcurr(l + 9, dim) - xcurr(l, dim));
        bop_loc(5, (l + 9) * 3 + dim) +=
            xsi_[gp](2) * .25 * deriv_q9(0, k) * shape_q9(l) * (xcurr(k + 9, dim) - xcurr(k, dim));
        bop_loc(5, l * 3 + dim) +=
            -xsi_[gp](2) * .25 * deriv_q9(0, k) * shape_q9(l) * (xcurr(k + 9, dim) - xcurr(k, dim));

        // transverse normal strain
        if (dsg_ctl_)
        {
          bop_loc(2, (k + 9) * 3 + dim) +=
              .125 * dsg_transverse_t_[gp % 9](k, l) * (xcurr(l + 9, dim) - xcurr(l, dim));
          bop_loc(2, k * 3 + dim) -=
              .125 * dsg_transverse_t_[gp % 9](k, l) * (xcurr(l + 9, dim) - xcurr(l, dim));
          bop_loc(2, (l + 9) * 3 + dim) +=
              .125 * dsg_transverse_t_[gp % 9](k, l) * (xcurr(k + 9, dim) - xcurr(k, dim));
          bop_loc(2, l * 3 + dim) -=
              .125 * dsg_transverse_t_[gp % 9](k, l) * (xcurr(k + 9, dim) - xcurr(k, dim));
        }
        else
        {
          bop_loc(2, (k + 9) * 3 + dim) +=
              .125 * shape_q9(k) * shape_q9(l) * (xcurr(l + 9, dim) - xcurr(l, dim));
          bop_loc(2, k * 3 + dim) +=
              -.125 * shape_q9(k) * shape_q9(l) * (xcurr(l + 9, dim) - xcurr(l, dim));
          bop_loc(2, (l + 9) * 3 + dim) +=
              .125 * shape_q9(k) * shape_q9(l) * (xcurr(k + 9, dim) - xcurr(k, dim));
          bop_loc(2, l * 3 + dim) +=
              -.125 * shape_q9(k) * shape_q9(l) * (xcurr(k + 9, dim) - xcurr(k, dim));
        }
      }  // k=0..8; l=0..8; dim=0..2
  return;
}

void Discret::ELEMENTS::SoSh18::calculate_geo_stiff(const Core::LinAlg::Matrix<9, 1>& shape_q9,
    const Core::LinAlg::Matrix<2, 9>& deriv_q9,
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, Mat::NUM_STRESS_3D>& TinvT, const int gp,
    const double detJ_w, const Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1>& stress,
    Core::LinAlg::Matrix<NUMDOF_SOH18, NUMDOF_SOH18>* stiffmatrix)
{
  // intergrate `geometric' stiffness matrix and add to keu *****************
  for (int k = 0; k < 9; ++k)
    for (int l = 0; l < 9; ++l)
    {
      Core::LinAlg::Matrix<6, 1> G_kl, G_klp, G_kpl, G_kplp, G_lk, G_lkp, G_lpk, G_lpkp;

      // Normalverzerrungen
      if (dsg_membrane_)
      {
        // constant normal strain
        G_kl(0) += .125 * dsg_membrane_r_[gp % 9](k, l);
        G_klp(0) += .125 * dsg_membrane_r_[gp % 9](k, l);
        G_kpl(0) += .125 * dsg_membrane_r_[gp % 9](k, l);
        G_kplp(0) += .125 * dsg_membrane_r_[gp % 9](k, l);
        G_lk(0) += .125 * dsg_membrane_r_[gp % 9](k, l);
        G_lkp(0) += .125 * dsg_membrane_r_[gp % 9](k, l);
        G_lpk(0) += .125 * dsg_membrane_r_[gp % 9](k, l);
        G_lpkp(0) += .125 * dsg_membrane_r_[gp % 9](k, l);

        G_kl(1) += .125 * dsg_membrane_s_[gp % 9](k, l);
        G_klp(1) += .125 * dsg_membrane_s_[gp % 9](k, l);
        G_kpl(1) += .125 * dsg_membrane_s_[gp % 9](k, l);
        G_kplp(1) += .125 * dsg_membrane_s_[gp % 9](k, l);
        G_lk(1) += .125 * dsg_membrane_s_[gp % 9](k, l);
        G_lkp(1) += .125 * dsg_membrane_s_[gp % 9](k, l);
        G_lpk(1) += .125 * dsg_membrane_s_[gp % 9](k, l);
        G_lpkp(1) += .125 * dsg_membrane_s_[gp % 9](k, l);

        // constant in-plane shear strain
        G_kl(3) += .25 * dsg_membrane_rs_[gp % 9](k, l);
        G_klp(3) += .25 * dsg_membrane_rs_[gp % 9](k, l);
        G_kpl(3) += .25 * dsg_membrane_rs_[gp % 9](k, l);
        G_kplp(3) += .25 * dsg_membrane_rs_[gp % 9](k, l);
        G_lk(3) += .25 * dsg_membrane_rs_[gp % 9](k, l);
        G_lkp(3) += .25 * dsg_membrane_rs_[gp % 9](k, l);
        G_lpk(3) += .25 * dsg_membrane_rs_[gp % 9](k, l);
        G_lpkp(3) += .25 * dsg_membrane_rs_[gp % 9](k, l);
      }
      else
      {
        // constant normal strain
        for (int alpha = 0; alpha < 2; ++alpha)
        {
          G_kl(alpha) += .125 * deriv_q9(alpha, k) * deriv_q9(alpha, l);
          G_klp(alpha) += .125 * deriv_q9(alpha, k) * deriv_q9(alpha, l);
          G_kpl(alpha) += .125 * deriv_q9(alpha, k) * deriv_q9(alpha, l);
          G_kplp(alpha) += .125 * deriv_q9(alpha, k) * deriv_q9(alpha, l);
          G_lk(alpha) += .125 * deriv_q9(alpha, k) * deriv_q9(alpha, l);
          G_lkp(alpha) += .125 * deriv_q9(alpha, k) * deriv_q9(alpha, l);
          G_lpk(alpha) += .125 * deriv_q9(alpha, k) * deriv_q9(alpha, l);
          G_lpkp(alpha) += .125 * deriv_q9(alpha, k) * deriv_q9(alpha, l);
        }

        // constant in-plane shear strain
        G_kl(3) += .25 * deriv_q9(0, k) * deriv_q9(1, l);
        G_klp(3) += .25 * deriv_q9(0, k) * deriv_q9(1, l);
        G_kpl(3) += .25 * deriv_q9(0, k) * deriv_q9(1, l);
        G_kplp(3) += .25 * deriv_q9(0, k) * deriv_q9(1, l);
        G_lk(3) += .25 * deriv_q9(0, k) * deriv_q9(1, l);
        G_lkp(3) += .25 * deriv_q9(0, k) * deriv_q9(1, l);
        G_lpk(3) += .25 * deriv_q9(0, k) * deriv_q9(1, l);
        G_lpkp(3) += .25 * deriv_q9(0, k) * deriv_q9(1, l);
      }

      // linear normal strain
      for (int alpha = 0; alpha < 2; ++alpha)
      {
        G_kl(alpha) -= xsi_[gp](2) * .25 * deriv_q9(alpha, k) * deriv_q9(alpha, l);
        G_klp(alpha) += 0.;
        G_kpl(alpha) += 0.;
        G_kplp(alpha) += xsi_[gp](2) * .25 * deriv_q9(alpha, k) * deriv_q9(alpha, l);
        G_lk(alpha) -= xsi_[gp](2) * .25 * deriv_q9(alpha, k) * deriv_q9(alpha, l);
        G_lkp(alpha) += 0.;
        G_lpk(alpha) += 0.;
        G_lpkp(alpha) += xsi_[gp](2) * .25 * deriv_q9(alpha, k) * deriv_q9(alpha, l);
      }

      // linear in-plane shear strain
      G_kl(3) -= xsi_[gp](2) * .5 * deriv_q9(0, k) * deriv_q9(1, l);
      G_klp(3) += 0.;
      G_kpl(3) += 0.;
      G_kplp(3) += xsi_[gp](2) * .5 * deriv_q9(0, k) * deriv_q9(1, l);
      G_lk(3) -= xsi_[gp](2) * .5 * deriv_q9(0, k) * deriv_q9(1, l);
      G_lkp(3) += 0.;
      G_lpk(3) += 0.;
      G_lpkp(3) += xsi_[gp](2) * .5 * deriv_q9(0, k) * deriv_q9(1, l);

      if (dsg_shear_)
      {
        // constant transverse shear strain
        G_kl(4) -= .25 * dsg_shear_s_[gp % 9](k, l);
        G_klp(4) += .25 * dsg_shear_s_[gp % 9](k, l);
        G_kpl(4) -= .25 * dsg_shear_s_[gp % 9](k, l);
        G_kplp(4) += .25 * dsg_shear_s_[gp % 9](k, l);
        G_lk(4) -= .25 * dsg_shear_s_[gp % 9](k, l);
        G_lkp(4) -= .25 * dsg_shear_s_[gp % 9](k, l);
        G_lpk(4) += .25 * dsg_shear_s_[gp % 9](k, l);
        G_lpkp(4) += .25 * dsg_shear_s_[gp % 9](k, l);

        G_kl(5) -= .25 * dsg_shear_r_[gp % 9](k, l);
        G_klp(5) += .25 * dsg_shear_r_[gp % 9](k, l);
        G_kpl(5) -= .25 * dsg_shear_r_[gp % 9](k, l);
        G_kplp(5) += .25 * dsg_shear_r_[gp % 9](k, l);
        G_lk(5) -= .25 * dsg_shear_r_[gp % 9](k, l);
        G_lkp(5) -= .25 * dsg_shear_r_[gp % 9](k, l);
        G_lpk(5) += .25 * dsg_shear_r_[gp % 9](k, l);
        G_lpkp(5) += .25 * dsg_shear_r_[gp % 9](k, l);
      }
      else
      {
        // constant transverse shear strain
        G_kl(4) -= .25 * deriv_q9(1, k) * shape_q9(l);
        G_klp(4) += .25 * deriv_q9(1, k) * shape_q9(l);
        G_kpl(4) -= .25 * deriv_q9(1, k) * shape_q9(l);
        G_kplp(4) += .25 * deriv_q9(1, k) * shape_q9(l);
        G_lk(4) -= .25 * deriv_q9(1, k) * shape_q9(l);
        G_lkp(4) -= .25 * deriv_q9(1, k) * shape_q9(l);
        G_lpk(4) += .25 * deriv_q9(1, k) * shape_q9(l);
        G_lpkp(4) += .25 * deriv_q9(1, k) * shape_q9(l);

        G_kl(5) -= .25 * deriv_q9(0, k) * shape_q9(l);
        G_klp(5) += .25 * deriv_q9(0, k) * shape_q9(l);
        G_kpl(5) -= .25 * deriv_q9(0, k) * shape_q9(l);
        G_kplp(5) += .25 * deriv_q9(0, k) * shape_q9(l);
        G_lk(5) -= .25 * deriv_q9(0, k) * shape_q9(l);
        G_lkp(5) -= .25 * deriv_q9(0, k) * shape_q9(l);
        G_lpk(5) += .25 * deriv_q9(0, k) * shape_q9(l);
        G_lpkp(5) += .25 * deriv_q9(0, k) * shape_q9(l);
      }

      // linear transverse shear strain
      G_kl(4) += xsi_[gp](2) * .25 * deriv_q9(1, k) * shape_q9(l);
      G_klp(4) -= xsi_[gp](2) * .25 * deriv_q9(1, k) * shape_q9(l);
      G_kpl(4) -= xsi_[gp](2) * .25 * deriv_q9(1, k) * shape_q9(l);
      G_kplp(4) += xsi_[gp](2) * .25 * deriv_q9(1, k) * shape_q9(l);
      G_lk(4) += xsi_[gp](2) * .25 * deriv_q9(1, k) * shape_q9(l);
      G_lkp(4) -= xsi_[gp](2) * .25 * deriv_q9(1, k) * shape_q9(l);
      G_lpk(4) -= xsi_[gp](2) * .25 * deriv_q9(1, k) * shape_q9(l);
      G_lpkp(4) += xsi_[gp](2) * .25 * deriv_q9(1, k) * shape_q9(l);


      G_kl(5) += xsi_[gp](2) * .25 * deriv_q9(0, k) * shape_q9(l);
      G_klp(5) -= xsi_[gp](2) * .25 * deriv_q9(0, k) * shape_q9(l);
      G_kpl(5) -= xsi_[gp](2) * .25 * deriv_q9(0, k) * shape_q9(l);
      G_kplp(5) += xsi_[gp](2) * .25 * deriv_q9(0, k) * shape_q9(l);
      G_lk(5) += xsi_[gp](2) * .25 * deriv_q9(0, k) * shape_q9(l);
      G_lkp(5) -= xsi_[gp](2) * .25 * deriv_q9(0, k) * shape_q9(l);
      G_lpk(5) -= xsi_[gp](2) * .25 * deriv_q9(0, k) * shape_q9(l);
      G_lpkp(5) += xsi_[gp](2) * .25 * deriv_q9(0, k) * shape_q9(l);

      // transverse normal strain
      if (dsg_ctl_)
      {
        G_kl(2) += .125 * dsg_transverse_t_[gp % 9](k, l);
        G_klp(2) -= .125 * dsg_transverse_t_[gp % 9](k, l);
        G_kpl(2) -= .125 * dsg_transverse_t_[gp % 9](k, l);
        G_kplp(2) += .125 * dsg_transverse_t_[gp % 9](k, l);
        G_lk(2) += .125 * dsg_transverse_t_[gp % 9](k, l);
        G_lkp(2) -= .125 * dsg_transverse_t_[gp % 9](k, l);
        G_lpk(2) -= .125 * dsg_transverse_t_[gp % 9](k, l);
        G_lpkp(2) += .125 * dsg_transverse_t_[gp % 9](k, l);
      }
      else
      {
        G_kl(2) += .125 * shape_q9(k) * shape_q9(l);
        G_klp(2) -= .125 * shape_q9(k) * shape_q9(l);
        G_kpl(2) -= .125 * shape_q9(k) * shape_q9(l);
        G_kplp(2) += .125 * shape_q9(k) * shape_q9(l);
        G_lk(2) += .125 * shape_q9(k) * shape_q9(l);
        G_lkp(2) -= .125 * shape_q9(k) * shape_q9(l);
        G_lpk(2) -= .125 * shape_q9(k) * shape_q9(l);
        G_lpkp(2) += .125 * shape_q9(k) * shape_q9(l);
      }

      Core::LinAlg::Matrix<6, 1> G_kl_g;
      G_kl_g.multiply(TinvT, G_kl);
      const double Gkl = detJ_w * stress.dot(G_kl_g);
      Core::LinAlg::Matrix<6, 1> G_klp_g;
      G_klp_g.multiply(TinvT, G_klp);
      const double Gklp = detJ_w * stress.dot(G_klp_g);
      Core::LinAlg::Matrix<6, 1> G_kpl_g;
      G_kpl_g.multiply(TinvT, G_kpl);
      const double Gkpl = detJ_w * stress.dot(G_kpl_g);
      Core::LinAlg::Matrix<6, 1> G_kplp_g;
      G_kplp_g.multiply(TinvT, G_kplp);
      const double Gkplp = detJ_w * stress.dot(G_kplp_g);
      Core::LinAlg::Matrix<6, 1> G_lk_g;
      G_lk_g.multiply(TinvT, G_lk);
      const double Glk = detJ_w * stress.dot(G_lk_g);
      Core::LinAlg::Matrix<6, 1> G_lkp_g;
      G_lkp_g.multiply(TinvT, G_lkp);
      const double Glkp = detJ_w * stress.dot(G_lkp_g);
      Core::LinAlg::Matrix<6, 1> G_lpk_g;
      G_lpk_g.multiply(TinvT, G_lpk);
      const double Glpk = detJ_w * stress.dot(G_lpk_g);
      Core::LinAlg::Matrix<6, 1> G_lpkp_g;
      G_lpkp_g.multiply(TinvT, G_lpkp);
      const double Glpkp = detJ_w * stress.dot(G_lpkp_g);
      for (int dim = 0; dim < 3; ++dim)
      {
        (*stiffmatrix)(k * 3 + dim, l * 3 + dim) += Gkl;
        (*stiffmatrix)(k * 3 + dim, (l + 9) * 3 + dim) += Gklp;
        (*stiffmatrix)((k + 9) * 3 + dim, l * 3 + dim) += Gkpl;
        (*stiffmatrix)((k + 9) * 3 + dim, (l + 9) * 3 + dim) += Gkplp;

        (*stiffmatrix)(l * 3 + dim, k * 3 + dim) += Glk;
        (*stiffmatrix)(l * 3 + dim, (k + 9) * 3 + dim) += Glkp;
        (*stiffmatrix)((l + 9) * 3 + dim, k * 3 + dim) += Glpk;
        (*stiffmatrix)((l + 9) * 3 + dim, (k + 9) * 3 + dim) += Glpkp;
      }
    }  // k=0..8 l=0..8
  // end of integrate `geometric' stiffness******************************
}

void Discret::ELEMENTS::SoSh18::calc_consistent_defgrd(
    const Core::LinAlg::Matrix<3, 3>& defgrd_disp, Core::LinAlg::Matrix<6, 1> glstrain_mod,
    Core::LinAlg::Matrix<3, 3>& defgrd_mod)
{
  Core::LinAlg::Matrix<3, 3> R;       // rotation tensor
  Core::LinAlg::Matrix<3, 3> U_mod;   // modified right stretch tensor
  Core::LinAlg::Matrix<3, 3> U_disp;  // displacement-based right stretch tensor
  Core::LinAlg::Matrix<3, 3> EW;      // temporarily store eigenvalues
  Core::LinAlg::Matrix<3, 3> tmp;     // temporary matrix for matrix matrix matrix products
  Core::LinAlg::Matrix<3, 3> tmp2;    // temporary matrix for matrix matrix matrix products

  // ******************************************************************
  // calculate modified right stretch tensor
  // ******************************************************************
  for (int i = 0; i < 3; i++) U_mod(i, i) = 2. * glstrain_mod(i) + 1.;
  U_mod(0, 1) = glstrain_mod(3);
  U_mod(1, 0) = glstrain_mod(3);
  U_mod(1, 2) = glstrain_mod(4);
  U_mod(2, 1) = glstrain_mod(4);
  U_mod(0, 2) = glstrain_mod(5);
  U_mod(2, 0) = glstrain_mod(5);

  Core::LinAlg::SYEV(U_mod, EW, U_mod);
  for (int i = 0; i < 3; ++i) EW(i, i) = sqrt(EW(i, i));
  tmp.multiply(U_mod, EW);
  tmp2.multiply_nt(tmp, U_mod);
  U_mod.update(tmp2);

  // ******************************************************************
  // calculate displacement-based right stretch tensor
  // ******************************************************************
  U_disp.multiply_tn(defgrd_disp, defgrd_disp);

  Core::LinAlg::SYEV(U_disp, EW, U_disp);
  for (int i = 0; i < 3; ++i) EW(i, i) = sqrt(EW(i, i));
  tmp.multiply(U_disp, EW);
  tmp2.multiply_nt(tmp, U_disp);
  U_disp.update(tmp2);

  // ******************************************************************
  // compose consistent deformation gradient
  // ******************************************************************
  U_disp.invert();
  R.multiply(defgrd_disp, U_disp);
  defgrd_mod.multiply(R, U_mod);

  // you're done here
  return;
}

/*----------------------------------------------------------------------*
 |  integrate DSG integral                                  seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoSh18::integrate_dsg_shear_r(
    const int gp, Core::LinAlg::Matrix<9, 9>& dsg_shear_r)
{
  dsg_shear_r.clear();
  const std::array<double, 2> coord_refNode = {0., 0.};

  Core::LinAlg::Matrix<2, 9> deriv;
  Core::FE::shape_function_deriv1<Core::FE::CellType::quad9>(xsi_[gp], deriv);

  Core::FE::IntPointsAndWeights<1> ip(Core::FE::GaussRule1D::line_2point);
  for (int i = 0; i < 9; ++i)
  {
    const Core::LinAlg::Matrix<3, 1> coord_i = node_param_coord(i);

    // integrations with non-empty integration domain
    if (coord_i(0) != coord_refNode[0])
      // perform integration
      for (int gp = 0; gp < ip.IP().nquad; ++gp)
      {
        const double jac = .5 * (coord_i(0) - coord_refNode[0]);

        // gauss point coordinates in element parameter space
        Core::LinAlg::Matrix<2, 1> xi_gp;
        xi_gp(0) = .5 * (coord_i(0) + coord_refNode[0]) + jac * ip.IP().qxg[gp][0];
        xi_gp(1) = coord_i(1);

        // shape function
        Core::LinAlg::Matrix<9, 1> shape_gp;
        Core::FE::shape_function<Core::FE::CellType::quad9>(xi_gp, shape_gp);
        // derivative
        Core::LinAlg::Matrix<2, 9> deriv_gp;
        Core::FE::shape_function_deriv1<Core::FE::CellType::quad9>(xi_gp, deriv_gp);

        for (int k = 0; k < 9; ++k)
          for (int l = 0; l < 9; ++l)
            dsg_shear_r(k, l) +=
                deriv(0, i) * deriv_gp(0, k) * shape_gp(l) * jac * ip.IP().qwgt[gp];
      }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  integrate DSG integral                                  seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoSh18::integrate_dsg_shear_s(
    const int gp, Core::LinAlg::Matrix<9, 9>& dsg_shear_s)
{
  dsg_shear_s.clear();
  const std::array<double, 2> coord_refNode = {0., 0.};

  Core::LinAlg::Matrix<2, 9> deriv;
  Core::FE::shape_function_deriv1<Core::FE::CellType::quad9>(xsi_[gp], deriv);

  Core::FE::IntPointsAndWeights<1> ip(Core::FE::GaussRule1D::line_2point);
  for (int i = 0; i < 9; ++i)
  {
    const Core::LinAlg::Matrix<3, 1> coord_i = node_param_coord(i);

    // integrations with non-empty integration domain
    if (coord_i(1) != coord_refNode[1])
      // perform integration
      for (int gp = 0; gp < ip.IP().nquad; ++gp)
      {
        // integration jacobian
        const double jac = .5 * (coord_i(1) - coord_refNode[1]);

        // gauss point coordinates in element parameter space
        Core::LinAlg::Matrix<2, 1> xi_gp;
        xi_gp(0) = coord_i(0);
        xi_gp(1) = .5 * (coord_i(1) + coord_refNode[1]) + jac * ip.IP().qxg[gp][0];

        // shape function
        Core::LinAlg::Matrix<9, 1> shape_gp;
        Core::FE::shape_function<Core::FE::CellType::quad9>(xi_gp, shape_gp);
        // derivative
        Core::LinAlg::Matrix<2, 9> deriv_gp;
        Core::FE::shape_function_deriv1<Core::FE::CellType::quad9>(xi_gp, deriv_gp);

        for (int k = 0; k < 9; ++k)
          for (int l = 0; l < 9; ++l)
            dsg_shear_s(k, l) +=
                deriv(1, i) * deriv_gp(1, k) * shape_gp(l) * jac * ip.IP().qwgt[gp];
      }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  integrate DSG integral                                  seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoSh18::integrate_dsg_membrane_rs(
    const int gp, Core::LinAlg::Matrix<9, 9>& dsg_membrane_rs)
{
  dsg_membrane_rs.clear();

  // integration
  const std::array<double, 2> coord_refNode = {0., 0.};
  const Core::LinAlg::Matrix<18, 3> coords = node_param_coord();
  Core::FE::IntPointsAndWeights<1> ip(Core::FE::GaussRule1D::line_2point);
  Core::LinAlg::Matrix<2, 9> deriv_xieta;
  Core::FE::shape_function_deriv1<Core::FE::CellType::quad9>(xsi_[gp], deriv_xieta);

  for (int k = 0; k < 9; ++k)
    for (int l = 0; l < 9; ++l)
    {
      for (int r = 0; r < 9; ++r)
        for (int s = 0; s < 9; ++s)

          for (int g = 0; g < ip.IP().nquad; ++g)
            for (int h = 0; h < ip.IP().nquad; ++h)
            {
              const double jac_g = .5 * (coords(r, 0) - coord_refNode[0]);
              const double jac_h = .5 * (coords(s, 1) - coord_refNode[1]);
              const double g_loc =
                  .5 * (coords(r, 0) + coord_refNode[0]) + jac_g * ip.IP().qxg[g][0];
              const double h_loc =
                  .5 * (coords(s, 1) + coord_refNode[1]) + jac_h * ip.IP().qxg[h][0];

              Core::LinAlg::Matrix<2, 1> xsi_g_eta;
              xsi_g_eta(0) = g_loc;
              xsi_g_eta(1) = xsi_[gp](1);
              Core::LinAlg::Matrix<2, 1> xsi_g_h;
              xsi_g_h(0) = g_loc;
              xsi_g_h(1) = h_loc;
              Core::LinAlg::Matrix<2, 9> deriv_g_eta;
              Core::FE::shape_function_deriv1<Core::FE::CellType::quad9>(xsi_g_eta, deriv_g_eta);
              Core::LinAlg::Matrix<2, 9> deriv_g_h;
              Core::FE::shape_function_deriv1<Core::FE::CellType::quad9>(xsi_g_h, deriv_g_h);

              if (coords(r, 0) != coord_refNode[0] && coords(s, 1) != coord_refNode[1])
                dsg_membrane_rs(k, l) += deriv_xieta(0, r) * deriv_g_eta(1, s) * deriv_g_h(0, k) *
                                         deriv_g_h(1, l) * jac_g * ip.IP().qwgt[g] * jac_h *
                                         ip.IP().qwgt[h];
            }
    }

  return;
}

/*----------------------------------------------------------------------*
 |  integrate DSG integral                                  seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoSh18::integrate_dsg_membrane_r(
    const int gp, Core::LinAlg::Matrix<9, 9>& dsg_membrane_r)
{
  dsg_membrane_r.clear();
  const std::array<double, 2> coord_refNode = {0., 0.};

  Core::LinAlg::Matrix<2, 9> deriv;
  Core::FE::shape_function_deriv1<Core::FE::CellType::quad9>(xsi_[gp], deriv);

  Core::FE::IntPointsAndWeights<1> ip(Core::FE::GaussRule1D::line_2point);
  for (int i = 0; i < 9; ++i)
  {
    const Core::LinAlg::Matrix<3, 1> coord_i = node_param_coord(i);

    // integrations with non-empty integration domain
    if (coord_i(0) != coord_refNode[0])
      // perform integration
      for (int gp = 0; gp < ip.IP().nquad; ++gp)
      {
        // integration jacobian
        const double jac = .5 * (coord_i(0) - coord_refNode[0]);

        // gauss point coordinates in element parameter space
        Core::LinAlg::Matrix<2, 1> xi_gp;
        xi_gp(0) = .5 * (coord_i(0) + coord_refNode[0]) + jac * ip.IP().qxg[gp][0];
        xi_gp(1) = coord_i(1);

        // derivative
        Core::LinAlg::Matrix<2, 9> deriv_gp;
        Core::FE::shape_function_deriv1<Core::FE::CellType::quad9>(xi_gp, deriv_gp);

        // fill up array
        for (int k = 0; k < 9; ++k)
          for (int l = 0; l < 9; ++l)
            dsg_membrane_r(k, l) +=
                deriv(0, i) * deriv_gp(0, k) * deriv_gp(0, l) * jac * ip.IP().qwgt[gp];
      }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  integrate DSG integral                                  seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoSh18::integrate_dsg_membrane_s(
    const int gp, Core::LinAlg::Matrix<9, 9>& dsg_membrane_s)
{
  dsg_membrane_s.clear();
  const std::array<double, 2> coord_refNode = {0., 0.};

  Core::LinAlg::Matrix<2, 9> deriv;
  Core::FE::shape_function_deriv1<Core::FE::CellType::quad9>(xsi_[gp], deriv);

  Core::FE::IntPointsAndWeights<1> ip(Core::FE::GaussRule1D::line_2point);
  for (int i = 0; i < 9; ++i)
  {
    const Core::LinAlg::Matrix<3, 1> coord_i = node_param_coord(i);

    // integrations with non-empty integration domain
    if (coord_i(1) != coord_refNode[1])
      // perform integration
      for (int gp = 0; gp < ip.IP().nquad; ++gp)
      {
        // integration jacobian
        const double jac = .5 * (coord_i(1) - coord_refNode[1]);

        // gauss point coordinates in element parameter space
        Core::LinAlg::Matrix<2, 1> xi_gp;
        xi_gp(0) = coord_i(0);
        xi_gp(1) = .5 * (coord_i(1) + coord_refNode[1]) + jac * ip.IP().qxg[gp][0];

        // derivative
        Core::LinAlg::Matrix<2, 9> deriv_gp;
        Core::FE::shape_function_deriv1<Core::FE::CellType::quad9>(xi_gp, deriv_gp);

        // fill up array
        for (int k = 0; k < 9; ++k)
          for (int l = 0; l < 9; ++l)
            dsg_membrane_s(k, l) +=
                deriv(1, i) * deriv_gp(1, k) * deriv_gp(1, l) * jac * ip.IP().qwgt[gp];
      }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  integrate DSG integral                                  seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoSh18::integrate_dsg_transverse_t(
    const int gp, Core::LinAlg::Matrix<9, 9>& dsg_transverse_t)
{
  // reset
  dsg_transverse_t.clear();
  Core::LinAlg::Matrix<9, 1> shape;
  Core::FE::shape_function<Core::FE::CellType::quad9>(xsi_[gp], shape);

  for (int i = 0; i < 9; ++i)
  {
    const Core::LinAlg::Matrix<3, 1> coord_i = node_param_coord(i);
    Core::LinAlg::Matrix<9, 1> shape_gp;
    Core::FE::shape_function<Core::FE::CellType::quad9>(coord_i, shape_gp);
    for (int k = 0; k < 9; ++k)
      for (int l = 0; l < 9; ++l) dsg_transverse_t(k, l) += shape(i) * shape_gp(k) * shape_gp(l);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  setup all necessary DSG terms                           seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoSh18::setup_dsg()
{
  if (dsg_shear_)
  {
    dsg_shear_r_.resize(9);
    dsg_shear_s_.resize(9);
  }
  else
  {
    dsg_shear_r_.resize(0);
    dsg_shear_s_.resize(0);
  }
  if (dsg_membrane_)
  {
    dsg_membrane_r_.resize(9);
    dsg_membrane_s_.resize(9);
    dsg_membrane_rs_.resize(9);
  }
  else
  {
    dsg_membrane_r_.resize(0);
    dsg_membrane_s_.resize(0);
    dsg_membrane_rs_.resize(0);
  }
  if (dsg_ctl_)
    dsg_transverse_t_.resize(9);
  else
    dsg_transverse_t_.resize(0);

  // NOTE: For Lagrange polynomials the following arrays
  // are very sparse. To increase performance of Langrange
  // FE one might rather just code the non-zeros instead
  // of all those loops.
  // However, they represent what actually happens
  // for DSG elements. Moreover, other ansatz functions
  // e.g. Nurbs should hopefully be "unlocked" efficiently
  // without any further modification.
  for (int gp = 0; gp < 9; ++gp)
  {
    if (dsg_shear_)
    {
      integrate_dsg_shear_r(gp, dsg_shear_r_[gp]);
      integrate_dsg_shear_s(gp, dsg_shear_s_[gp]);
    }
    if (dsg_membrane_)
    {
      integrate_dsg_membrane_r(gp, dsg_membrane_r_[gp]);
      integrate_dsg_membrane_s(gp, dsg_membrane_s_[gp]);
      integrate_dsg_membrane_rs(gp, dsg_membrane_rs_[gp]);
    }
    if (dsg_ctl_)
    {
      integrate_dsg_transverse_t(gp, dsg_transverse_t_[gp]);
    }
  }

  return;
}

Core::LinAlg::Matrix<18, 3> Discret::ELEMENTS::SoSh18::node_param_coord()
{
  Core::LinAlg::Matrix<18, 3> coord;
  for (int node = 0; node < NUMNOD_SOH18; ++node)
  {
    Core::LinAlg::Matrix<3, 1> nodeCoord = node_param_coord(node);
    for (int i = 0; i < 3; ++i) coord(node, i) = nodeCoord(i);
  }
  return coord;
}

Core::LinAlg::Matrix<3, 1> Discret::ELEMENTS::SoSh18::node_param_coord(const int node)
{
  Core::LinAlg::Matrix<3, 1> coord;

  switch (node % 9)
  {
    case 0:
    case 3:
    case 7:
      coord(0) = -1.;
      break;
    case 4:
    case 6:
    case 8:
      coord(0) = +0.;
      break;
    case 1:
    case 2:
    case 5:
      coord(0) = +1.;
      break;
  }
  switch (node % 9)
  {
    case 0:
    case 1:
    case 4:
      coord(1) = -1.;
      break;
    case 5:
    case 7:
    case 8:
      coord(1) = +0.;
      break;
    case 2:
    case 3:
    case 6:
      coord(1) = +1.;
      break;
  }

  if (node < 9)
    coord(2) = -1.;
  else
    coord(2) = +1.;

  return coord;
}

/*----------------------------------------------------------------------*
 |  setup EAS terms                                         seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoSh18::eas_setup(
    std::vector<Core::LinAlg::Matrix<6, num_eas>>& M_gp,  // M-matrix evaluated at GPs
    Core::LinAlg::Matrix<3, 1>& G3_contra,
    const Core::LinAlg::Matrix<NUMNOD_SOH18, 3>& xrefe)  // material element coords
{
  // build EAS interpolation matrix M, evaluated at the GPs
  static std::vector<Core::LinAlg::Matrix<6, num_eas>> M(NUMGPT_SOH18);
  static Core::LinAlg::Matrix<3, 1> G3_c;
  static bool eval = false;
  if (!eval)
  {
    for (int gp = 0; gp < NUMGPT_SOH18; ++gp)
    {
      double r = xsi_[gp](0);
      double s = xsi_[gp](1);
      double t = xsi_[gp](2);
      M[gp](2, 0) = 1.;
      M[gp](2, 1) = r;
      M[gp](2, 2) = s;
      M[gp](2, 3) = r * s;
      M[gp](2, 4) = 1. - 3. * r * r;
      M[gp](2, 5) = 1. - 3. * s * s;
      M[gp](2, 6) = r * r * s;
      M[gp](2, 7) = r * s * s;
      M[gp](2, 8) = 1. - 9. * r * r * s * s;

      M[gp].scale(t);
    }
    eval = true;
  }
  // compute Jacobian, evaluated at element origin (r=s=t=0.0)
  Core::LinAlg::Matrix<NUMDIM_SOH18, NUMDIM_SOH18> jac0inv;
  Core::LinAlg::Matrix<2, 1> xsi_center(true);
  Core::LinAlg::Matrix<2, 9> deriv_q9;
  Core::FE::shape_function_deriv1<Core::FE::CellType::quad9>(xsi_center, deriv_q9);
  Core::LinAlg::Matrix<9, 1> shapefunct_q9;
  Core::FE::shape_function<Core::FE::CellType::quad9>(xsi_center, shapefunct_q9);
  for (int dim = 0; dim < 3; ++dim)
    for (int k = 0; k < 9; ++k)
    {
      jac0inv(0, dim) += .5 * deriv_q9(0, k) * (xrefe(k + 9, dim) + xrefe(k, dim));
      jac0inv(1, dim) += .5 * deriv_q9(1, k) * (xrefe(k + 9, dim) + xrefe(k, dim));
      jac0inv(2, dim) += .5 * shapefunct_q9(k) * (xrefe(k + 9, dim) - xrefe(k, dim));
    }
  jac0inv.invert();

  for (int dim = 0; dim < 3; ++dim) G3_c(dim) = jac0inv(2, dim);

  G3_contra = G3_c;
  M_gp = M;

  return;
}

/*----------------------------------------------------------------------*
 |  update EAS terms at the end of time step                seitz 11/14 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoSh18::update()
{
  if (eas_)
  {
    alpha_eas_delta_over_last_timestep_.update(alpha_eas_last_timestep_);
    alpha_eas_delta_over_last_timestep_.update(1., alpha_eas_, -1.);
    alpha_eas_last_timestep_.update(alpha_eas_);
    Kad_.clear();
    KaaInv_.clear();
    feas_.clear();
  }
}

void Discret::ELEMENTS::SoSh18::recover(const std::vector<double>& residual)
{
  if (not eas_) return;

  const double step_length = str_params_interface().get_step_length();
  Core::LinAlg::Matrix<NUMDOF_SOH18, 1> res_d;
  for (int i = 0; i < NUMDOF_SOH18; ++i) res_d(i) = residual[i];

  if (str_params_interface().is_default_step())
  {
    // first, store the eas state of the previous accepted Newton step
    str_params_interface().sum_into_my_previous_sol_norm(
        NOX::Nln::StatusTest::quantity_eas, num_eas, alpha_eas_.data(), Owner());

    feas_.multiply(1., Kad_, res_d, 1.);
    alpha_eas_inc_.multiply(-1., KaaInv_, feas_, 0.);
    alpha_eas_.update(step_length, alpha_eas_inc_, 1.);
  }
  else
  {
    alpha_eas_.update(-old_step_length_, alpha_eas_inc_, 1.);
    alpha_eas_inc_.scale(step_length / old_step_length_);
    alpha_eas_.update(1., alpha_eas_inc_, 1.);
  }
  old_step_length_ = step_length;

  str_params_interface().sum_into_my_update_norm(NOX::Nln::StatusTest::quantity_eas, num_eas,
      alpha_eas_inc_.data(), alpha_eas_.data(), step_length, Owner());
}

FOUR_C_NAMESPACE_CLOSE

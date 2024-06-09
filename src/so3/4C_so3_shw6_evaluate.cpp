/*----------------------------------------------------------------------*/
/*! \file
\brief
\level 1


*----------------------------------------------------------------------*/
#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_discretization_fem_general_utils_integration.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_densematrix_multiply.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_mat_micromaterial.hpp"
#include "4C_mat_viscoanisotropic.hpp"
#include "4C_so3_shw6.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_SerialDenseSolver.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              maf 04/07|
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoShw6::Evaluate(Teuchos::ParameterList& params,
    Discret::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
    Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
    Core::LinAlg::SerialDenseVector& elevec1_epetra,
    Core::LinAlg::SerialDenseVector& elevec2_epetra,
    Core::LinAlg::SerialDenseVector& elevec3_epetra)
{
  // Check whether the solid material post_setup() routine has already been called and call it if
  // not
  ensure_material_post_setup(params);

  // get parameter interface
  set_params_interface_ptr(params);

  Core::LinAlg::Matrix<NUMDOF_WEG6, NUMDOF_WEG6> elemat1(elemat1_epetra.values(), true);
  Core::LinAlg::Matrix<NUMDOF_WEG6, NUMDOF_WEG6> elemat2(elemat2_epetra.values(), true);
  Core::LinAlg::Matrix<NUMDOF_WEG6, 1> elevec1(elevec1_epetra.values(), true);
  Core::LinAlg::Matrix<NUMDOF_WEG6, 1> elevec2(elevec2_epetra.values(), true);
  Core::LinAlg::Matrix<NUMDOF_WEG6, 1> elevec3(elevec3_epetra.values(), true);

  // start with "none"
  Discret::ELEMENTS::SoWeg6::ActionType act = SoWeg6::none;

  // get the required action
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    FOUR_C_THROW("No action supplied");
  else if (action == "calc_struct_linstiff")
    act = SoWeg6::calc_struct_linstiff;
  else if (action == "calc_struct_nlnstiff")
    act = SoWeg6::calc_struct_nlnstiff;
  else if (action == "calc_struct_internalforce")
    act = SoWeg6::calc_struct_internalforce;
  else if (action == "calc_struct_linstiffmass")
    act = SoWeg6::calc_struct_linstiffmass;
  else if (action == "calc_struct_nlnstiffmass")
    act = SoWeg6::calc_struct_nlnstiffmass;
  else if (action == "calc_struct_nlnstifflmass")
    act = SoWeg6::calc_struct_nlnstifflmass;
  else if (action == "calc_struct_stress")
    act = SoWeg6::calc_struct_stress;
  else if (action == "calc_struct_eleload")
    act = SoWeg6::calc_struct_eleload;
  else if (action == "calc_struct_fsiload")
    act = SoWeg6::calc_struct_fsiload;
  else if (action == "calc_struct_update_istep")
    act = SoWeg6::calc_struct_update_istep;
  else if (action == "calc_struct_reset_istep")
    act = SoWeg6::calc_struct_reset_istep;
  else if (action == "calc_struct_recover")
    act = SoWeg6::calc_recover;
  else if (action == "calc_struct_predict")
    return 0;
  else
    FOUR_C_THROW("Unknown type of action for So_weg6");

  // what should the element do
  switch (act)
  {
    // linear stiffness
    case calc_struct_linstiff:
    {
      // need current displacement and residual forces
      std::vector<double> mydisp(lm.size());
      for (double& i : mydisp) i = 0.0;
      std::vector<double> myres(lm.size());
      for (double& myre : myres) myre = 0.0;
      soshw6_nlnstiffmass(lm, mydisp, myres, &elemat1, nullptr, &elevec1, nullptr, nullptr, nullptr,
          params, Inpar::STR::stress_none, Inpar::STR::strain_none);
    }
    break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (disp == Teuchos::null || res == Teuchos::null)
        FOUR_C_THROW("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      Core::FE::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      Core::FE::ExtractMyValues(*res, myres, lm);
      soshw6_nlnstiffmass(lm, mydisp, myres, &elemat1, nullptr, &elevec1, &elevec3, nullptr,
          nullptr, params, Inpar::STR::stress_none, Inpar::STR::strain_none);
    }
    break;

    // internal force vector only
    case calc_struct_internalforce:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (disp == Teuchos::null || res == Teuchos::null)
        FOUR_C_THROW("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      Core::FE::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      Core::FE::ExtractMyValues(*res, myres, lm);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      Core::LinAlg::Matrix<NUMDOF_WEG6, NUMDOF_WEG6> myemat(true);
      soshw6_nlnstiffmass(lm, mydisp, myres, &myemat, nullptr, &elevec1, nullptr, nullptr, nullptr,
          params, Inpar::STR::stress_none, Inpar::STR::strain_none);
    }
    break;

    // linear stiffness and consistent mass matrix
    case calc_struct_linstiffmass:
      FOUR_C_THROW("Case 'calc_struct_linstiffmass' not yet implemented");
      break;

    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (disp == Teuchos::null || res == Teuchos::null)
        FOUR_C_THROW("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      Core::FE::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      Core::FE::ExtractMyValues(*res, myres, lm);
      soshw6_nlnstiffmass(lm, mydisp, myres, &elemat1, &elemat2, &elevec1, &elevec3, nullptr,
          nullptr, params, Inpar::STR::stress_none, Inpar::STR::strain_none);
      if (act == calc_struct_nlnstifflmass) sow6_lumpmass(&elemat2);
    }
    break;

    // evaluate stresses at gauss points
    case calc_struct_stress:
    {
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      Teuchos::RCP<std::vector<char>> stressdata =
          params.get<Teuchos::RCP<std::vector<char>>>("stress", Teuchos::null);
      Teuchos::RCP<std::vector<char>> straindata =
          params.get<Teuchos::RCP<std::vector<char>>>("strain", Teuchos::null);
      if (disp == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'displacement'");
      if (stressdata == Teuchos::null) FOUR_C_THROW("Cannot get stress 'data'");
      if (straindata == Teuchos::null) FOUR_C_THROW("Cannot get strain 'data'");
      std::vector<double> mydisp(lm.size());
      Core::FE::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      Core::FE::ExtractMyValues(*res, myres, lm);
      Core::LinAlg::Matrix<NUMGPT_WEG6, Mat::NUM_STRESS_3D> stress;
      Core::LinAlg::Matrix<NUMGPT_WEG6, Mat::NUM_STRESS_3D> strain;
      auto iostress = Core::UTILS::GetAsEnum<Inpar::STR::StressType>(
          params, "iostress", Inpar::STR::stress_none);
      auto iostrain = Core::UTILS::GetAsEnum<Inpar::STR::StrainType>(
          params, "iostrain", Inpar::STR::strain_none);
      soshw6_nlnstiffmass(lm, mydisp, myres, nullptr, nullptr, nullptr, nullptr, &stress, &strain,
          params, iostress, iostrain);
      {
        Core::Communication::PackBuffer data;
        add_to_pack(data, stress);
        data.StartPacking();
        add_to_pack(data, stress);
        std::copy(data().begin(), data().end(), std::back_inserter(*stressdata));
      }
      {
        Core::Communication::PackBuffer data;
        add_to_pack(data, strain);
        data.StartPacking();
        add_to_pack(data, strain);
        std::copy(data().begin(), data().end(), std::back_inserter(*straindata));
      }
    }
    break;

    case calc_struct_eleload:
      FOUR_C_THROW("this method is not supposed to evaluate a load, use evaluate_neumann(...)");
      break;

    case calc_struct_fsiload:
      FOUR_C_THROW("Case not yet implemented");
      break;

    case calc_struct_update_istep:
    {
      // do something with internal EAS, etc parameters
      if (eastype_ == soshw6_easpoisthick)
      {
        auto* alpha = &easdata_.alpha;    // Alpha_{n+1}
        auto* alphao = &easdata_.alphao;  // Alpha_n
        // alphao := alpha
        Core::LinAlg::DenseFunctions::update<double, soshw6_easpoisthick, 1>(*alphao, *alpha);
      }
      SolidMaterial()->Update();
    }
    break;

    case calc_struct_reset_istep:
    {
      // do something with internal EAS, etc parameters
      if (eastype_ == soshw6_easpoisthick)
      {
        auto* alpha = &easdata_.alpha;    // Alpha_{n+1}
        auto* alphao = &easdata_.alphao;  // Alpha_n
        // alpha := alphao
        Core::LinAlg::DenseFunctions::update<double, soshw6_easpoisthick, 1>(*alpha, *alphao);
      }
      // Reset of history (if needed)
      SolidMaterial()->reset_step();
    }
    break;

    case calc_recover:
    {
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      std::vector<double> myres(lm.size());
      Core::FE::ExtractMyValues(*res, myres, lm);
      soshw6_recover(myres);
    }
    break;

    default:
      FOUR_C_THROW("Unknown type of action for Solid3");
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             maf 04/07|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoShw6::soshw6_nlnstiffmass(std::vector<int>& lm,  // location matrix
    std::vector<double>& disp,                                             // current displacements
    std::vector<double>& residual,                                         // current residual displ
    Core::LinAlg::Matrix<NUMDOF_WEG6, NUMDOF_WEG6>* stiffmatrix,  // element stiffness matrix
    Core::LinAlg::Matrix<NUMDOF_WEG6, NUMDOF_WEG6>* massmatrix,   // element mass matrix
    Core::LinAlg::Matrix<NUMDOF_WEG6, 1>* force,                  // element internal force vector
    Core::LinAlg::Matrix<NUMDOF_WEG6, 1>* force_str,              // structure force
    Core::LinAlg::Matrix<NUMGPT_WEG6, Mat::NUM_STRESS_3D>* elestress,  // stresses at GP
    Core::LinAlg::Matrix<NUMGPT_WEG6, Mat::NUM_STRESS_3D>* elestrain,  // strains at GP
    Teuchos::ParameterList& params,         // algorithmic parameters e.g. time
    const Inpar::STR::StressType iostress,  // stress output option
    const Inpar::STR::StrainType iostrain)  // strain output option
{
  /* ============================================================================*
  ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for Wedge_6 with 6 GAUSS POINTS*
  ** ============================================================================*/
  const static std::vector<Core::LinAlg::Matrix<NUMNOD_WEG6, 1>> shapefcts = sow6_shapefcts();
  const static std::vector<Core::LinAlg::Matrix<NUMDIM_WEG6, NUMNOD_WEG6>> derivs = sow6_derivs();
  const static std::vector<double> gpweights = sow6_weights();
  /* ============================================================================*/

  // update element geometry
  Core::LinAlg::Matrix<NUMNOD_WEG6, NUMDIM_WEG6> xrefe;  // material coord. of element
  Core::LinAlg::Matrix<NUMNOD_WEG6, NUMDIM_WEG6> xcurr;  // current  coord. of element
  Core::Nodes::Node** nodes = Nodes();
  for (int i = 0; i < NUMNOD_WEG6; ++i)
  {
    const auto& x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];

    xcurr(i, 0) = xrefe(i, 0) + disp[i * NODDOF_WEG6 + 0];
    xcurr(i, 1) = xrefe(i, 1) + disp[i * NODDOF_WEG6 + 1];
    xcurr(i, 2) = xrefe(i, 2) + disp[i * NODDOF_WEG6 + 2];
  }

  /*
  ** EAS Technology: declare, intialize, set up, and alpha history -------- EAS
  */
  // in any case declare variables, sizes etc. only in eascase
  Core::LinAlg::SerialDenseMatrix* alpha = nullptr;              // EAS alphas
  std::vector<Core::LinAlg::SerialDenseMatrix>* M_GP = nullptr;  // EAS matrix M at all GPs
  Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, soshw6_easpoisthick>
      M;                                 // EAS matrix M at current GP, fixed for sosh8
  Core::LinAlg::SerialDenseVector feas;  // EAS portion of internal forces
  Core::LinAlg::SerialDenseMatrix Kaa;   // EAS matrix Kaa
  Core::LinAlg::SerialDenseMatrix Kda;   // EAS matrix Kda
  double detJ0;                          // detJ(origin)
  Core::LinAlg::SerialDenseMatrix* oldfeas = nullptr;    // EAS history
  Core::LinAlg::SerialDenseMatrix* oldKaainv = nullptr;  // EAS history
  Core::LinAlg::SerialDenseMatrix* oldKda = nullptr;     // EAS history
  Core::LinAlg::SerialDenseMatrix* eas_inc = nullptr;    // EAS increment

  // transformation matrix T0, maps M-matrix evaluated at origin
  // between local element coords and global coords
  // here we already get the inverse transposed T0
  Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, Mat::NUM_STRESS_3D> T0invT;  // trafo matrix
  if (eastype_ == soshw6_easpoisthick)
  {
    /*
    ** EAS Update of alphas:
    ** the current alphas are (re-)evaluated out of
    ** Kaa and Kda of previous step to avoid additional element call.
    ** This corresponds to the (innermost) element update loop
    ** in the nonlinear FE-Skript page 120 (load-control alg. with EAS)
    */
    alpha = &easdata_.alpha;  // get old alpha
    // evaluate current (updated) EAS alphas (from history variables)
    // get stored EAS history
    oldfeas = &easdata_.feas;
    oldKaainv = &easdata_.invKaa;
    oldKda = &easdata_.Kda;
    eas_inc = &easdata_.eas_inc;
    if (!alpha || !oldKaainv || !oldKda || !oldfeas || !eas_inc)
      FOUR_C_THROW("Missing EAS history-data");

    // we need the (residual) displacement at the previous step
    Core::LinAlg::SerialDenseVector res_d(NUMDOF_WEG6);
    for (int i = 0; i < NUMDOF_WEG6; ++i)
    {
      res_d(i) = residual[i];
    }

    // this is a line search step, i.e. the direction of the eas increments
    // has been calculated by a Newton step and now it is only scaled
    if (not IsParamsInterface())
    {
      if (params.isParameter("alpha_ls"))
      {
        double alpha_ls = params.get<double>("alpha_ls");
        // undo step
        eas_inc->scale(-1.);
        alpha->operator+=(*eas_inc);
        // scale increment
        eas_inc->scale(-1. * alpha_ls);
        // add reduced increment
        alpha->operator+=(*eas_inc);
      }
      else
      {
        // add Kda . res_d to feas
        Core::LinAlg::DenseFunctions::multiply<double, soshw6_easpoisthick, NUMDOF_WEG6, 1>(
            1.0, *oldfeas, 1.0, *oldKda, res_d);
        // "new" alpha is: - Kaa^-1 . (feas + Kda . old_d), here: - Kaa^-1 . feas
        Core::LinAlg::DenseFunctions::multiply<double, soshw6_easpoisthick, soshw6_easpoisthick, 1>(
            0.0, *eas_inc, -1.0, *oldKaainv, *oldfeas);
        Core::LinAlg::DenseFunctions::update<double, soshw6_easpoisthick, 1>(
            1., *alpha, 1., *eas_inc);
      }
    }
    /* end of EAS Update ******************/

    // EAS portion of internal forces, also called enhacement vector s or Rtilde
    feas.size(neas_);

    // EAS matrix K_{alpha alpha}, also called Dtilde
    Kaa.shape(neas_, neas_);

    // EAS matrix K_{d alpha}
    Kda.shape(neas_, NUMDOF_WEG6);

    /* evaluation of EAS variables (which are constant for the following):
    ** -> M defining interpolation of enhanced strains alpha, evaluated at GPs
    ** -> determinant of Jacobi matrix at element origin (r=s=t=0.0)
    ** -> T0^{-T}
    */
    soshw6_eassetup(&M_GP, detJ0, T0invT, xrefe);
  }
  else if (eastype_ == soshw6_easnone)
  {
    // std::cout << "Warning: Solid-Shell Wegde6 without EAS" << std::endl;
  }
  else
    FOUR_C_THROW("Unknown EAS-type for solid wedge6");  // ------------------- EAS

  /*
  ** ANS Element technology to remedy
  *  - transverse-shear locking E_rt and E_st
  *  - trapezoidal (curvature-thickness) locking E_tt
  */
  // modified B-operator in local(parameter) element space
  // ANS modified rows of bop in local(parameter) coords
  Core::LinAlg::Matrix<num_ans * num_sp, NUMDOF_WEG6> B_ans_loc;
  // Jacobian evaluated at all ANS sampling points
  std::vector<Core::LinAlg::Matrix<NUMDIM_WEG6, NUMDIM_WEG6>> jac_sps(num_sp);
  // CURRENT Jacobian evaluated at all ANS sampling points
  std::vector<Core::LinAlg::Matrix<NUMDIM_WEG6, NUMDIM_WEG6>> jac_cur_sps(num_sp);
  // pointer to derivs evaluated at all sampling points
  std::vector<Core::LinAlg::Matrix<NUMDIM_WEG6, NUMNOD_WEG6>>* deriv_sp =
      nullptr;  // derivs eval. at all sampling points
  // evaluate all necessary variables for ANS
  soshw6_anssetup(xrefe, xcurr, &deriv_sp, jac_sps, jac_cur_sps, B_ans_loc);
  // (r,s) gp-locations of fully integrated linear 6-node wedge
  // necessary for ANS interpolation
  const Core::FE::GaussRule3D gaussrule_ = Core::FE::GaussRule3D::wedge_6point;
  const Core::FE::IntegrationPoints3D intpoints(gaussrule_);

  // check if we need to split the residuals (for Newton line search)
  // if true an additional global vector is assembled containing
  // the internal forces without the condensed EAS entries and the norm
  // of the EAS residual is calculated
  bool split_res = params.isParameter("cond_rhs_norm");

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp = 0; gp < NUMGPT_WEG6; ++gp)
  {
    /* compute the Jacobian matrix which looks like:
    **         [ x_,r  y_,r  z_,r ]
    **     J = [ x_,s  y_,s  z_,s ]
    **         [ x_,t  y_,t  z_,t ]
    */
    Core::LinAlg::Matrix<NUMDIM_WEG6, NUMDIM_WEG6> jac;
    jac.Multiply(derivs[gp], xrefe);

    // compute determinant of Jacobian by Sarrus' rule
    double detJ = jac.Determinant();
    if (abs(detJ) < 1E-16)
      FOUR_C_THROW("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0)
      FOUR_C_THROW("NEGATIVE JACOBIAN DETERMINANT");

    /* compute the CURRENT Jacobian matrix which looks like:
    **         [ xcurr_,r  ycurr_,r  zcurr_,r ]
    **  Jcur = [ xcurr_,s  ycurr_,s  zcurr_,s ]
    **         [ xcurr_,t  ycurr_,t  zcurr_,t ]
    ** Used to transform the global displacements into parametric space
    */
    Core::LinAlg::Matrix<NUMDIM_WEG6, NUMDIM_WEG6> jac_cur;
    jac_cur.Multiply(derivs[gp], xcurr);

    // need gp-locations for ANS interpolation
    const double r = intpoints.qxg[gp][0];
    const double s = intpoints.qxg[gp][1];
    // const double t = intpoints.qxg[gp][2]; // not needed

    // set up B-Operator in local(parameter) element space including ANS
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, NUMDOF_WEG6> bop_loc;
    for (int inode = 0; inode < NUMNOD_WEG6; ++inode)
    {
      for (int dim = 0; dim < NUMDIM_WEG6; ++dim)
      {
        // B_loc_rr = N_r.X_r
        bop_loc(0, inode * 3 + dim) = derivs[gp](0, inode) * jac_cur(0, dim);
        // B_loc_ss = N_s.X_s
        bop_loc(1, inode * 3 + dim) = derivs[gp](1, inode) * jac_cur(1, dim);
        // B_loc_tt = interpolation along (r x s) of ANS B_loc_tt
        //          = (1-r-s) * B_ans(SP C) + r * B_ans(SP D) + s * B_ans(SP E)
        bop_loc(2, inode * 3 + dim) = (1 - r - s) * B_ans_loc(0 + 2 * num_ans, inode * 3 + dim) +
                                      r * B_ans_loc(0 + 3 * num_ans, inode * 3 + dim) +
                                      s * B_ans_loc(0 + 4 * num_ans, inode * 3 + dim);
        // B_loc_rs = N_r.X_s + N_s.X_r
        bop_loc(3, inode * 3 + dim) =
            derivs[gp](0, inode) * jac_cur(1, dim) + derivs[gp](1, inode) * jac_cur(0, dim);
        // B_loc_st = interpolation along r of ANS B_loc_st
        //          = r * B_ans(SP B)
        bop_loc(4, inode * 3 + dim) = r * B_ans_loc(1 + 1 * num_ans, inode * 3 + dim);
        // B_loc_rt = interpolation along s of ANS B_loc_rt
        //          = s * B_ans(SP A)
        bop_loc(5, inode * 3 + dim) = s * B_ans_loc(2 + 0 * num_ans, inode * 3 + dim);

        //        // testing without ans:
        //        bop_loc(2,inode*3+dim) = deriv_gp(2,inode) * jac_cur(2,dim);
        //        bop_loc(4,inode*3+dim) = deriv_gp(1,inode) * jac_cur(2,dim)
        //                                +deriv_gp(2,inode) * jac_cur(1,dim);
        //        bop_loc(5,inode*3+dim) = deriv_gp(0,inode) * jac_cur(2,dim)
        //                                +deriv_gp(2,inode) * jac_cur(0,dim);
      }
    }

    // transformation from local (parameter) element space to global(material) space
    // with famous 'T'-matrix already used for EAS but now evaluated at each gp
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, Mat::NUM_STRESS_3D> TinvT;
    soshw6_evaluate_t(jac, TinvT);
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, NUMDOF_WEG6> bop;
    bop.Multiply(TinvT, bop_loc);

    // local GL strain vector lstrain={Err,Ess,Ett,2*Ers,2*Est,2*Ert}
    // but with modified ANS strains Ett, Est and Ert
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1> lstrain;
    // evaluate glstrains in local(parameter) coords
    // Err = 0.5 * (dx/dr * dx/dr^T - dX/dr * dX/dr^T)
    lstrain(0) = 0.5 * (+(jac_cur(0, 0) * jac_cur(0, 0) + jac_cur(0, 1) * jac_cur(0, 1) +
                            jac_cur(0, 2) * jac_cur(0, 2)) -
                           (jac(0, 0) * jac(0, 0) + jac(0, 1) * jac(0, 1) + jac(0, 2) * jac(0, 2)));
    // Ess = 0.5 * (dy/ds * dy/ds^T - dY/ds * dY/ds^T)
    lstrain(1) = 0.5 * (+(jac_cur(1, 0) * jac_cur(1, 0) + jac_cur(1, 1) * jac_cur(1, 1) +
                            jac_cur(1, 2) * jac_cur(1, 2)) -
                           (jac(1, 0) * jac(1, 0) + jac(1, 1) * jac(1, 1) + jac(1, 2) * jac(1, 2)));
    // Ers = (dx/ds * dy/dr^T - dX/ds * dY/dr^T)
    lstrain(3) = (+(jac_cur(0, 0) * jac_cur(1, 0) + jac_cur(0, 1) * jac_cur(1, 1) +
                      jac_cur(0, 2) * jac_cur(1, 2)) -
                  (jac(0, 0) * jac(1, 0) + jac(0, 1) * jac(1, 1) + jac(0, 2) * jac(1, 2)));

    //    // testing without ans:
    //    lstrain(2)= 0.5 * (
    //       +(jac_cur(2,0)*jac_cur(2,0) + jac_cur(2,1)*jac_cur(2,1) + jac_cur(2,2)*jac_cur(2,2))
    //       -(jac(2,0)*jac(2,0)         + jac(2,1)*jac(2,1)         + jac(2,2)*jac(2,2)));
    //    lstrain(4)= (
    //       +(jac_cur(1,0)*jac_cur(2,0) + jac_cur(1,1)*jac_cur(2,1) + jac_cur(1,2)*jac_cur(2,2))
    //       -(jac(1,0)*jac(2,0)         + jac(1,1)*jac(2,1)         + jac(1,2)*jac(2,2)));
    //    lstrain(5)= (
    //       +(jac_cur(0,0)*jac_cur(2,0) + jac_cur(0,1)*jac_cur(2,1) + jac_cur(0,2)*jac_cur(2,2))
    //       -(jac(0,0)*jac(2,0)         + jac(0,1)*jac(2,1)         + jac(0,2)*jac(2,2)));

    // ANS modification of strains ************************************** ANS
    double dydt_A = 0.0;
    double dYdt_A = 0.0;
    const int spA = 0;
    double dxdt_B = 0.0;
    double dXdt_B = 0.0;
    const int spB = 1;
    double dzdt_C = 0.0;
    double dZdt_C = 0.0;
    const int spC = 2;
    double dzdt_D = 0.0;
    double dZdt_D = 0.0;
    const int spD = 3;
    double dzdt_E = 0.0;
    double dZdt_E = 0.0;
    const int spE = 4;

    const int xdir = 0;  // index to matrix x-row, r-row respectively
    const int ydir = 1;  // index to matrix y-row, s-row respectively
    const int zdir = 2;  // index to matrix z-row, t-row respectively

    // vector product of rows of jacobians at corresponding sampling point
    for (int dim = 0; dim < NUMDIM_WEG6; ++dim)
    {
      dydt_A += jac_cur_sps[spA](xdir, dim) * jac_cur_sps[spA](zdir, dim);
      dYdt_A += jac_sps[spA](xdir, dim) * jac_sps[spA](zdir, dim);
      dxdt_B += jac_cur_sps[spB](ydir, dim) * jac_cur_sps[spB](zdir, dim);
      dXdt_B += jac_sps[spB](ydir, dim) * jac_sps[spB](zdir, dim);

      dzdt_C += jac_cur_sps[spC](zdir, dim) * jac_cur_sps[spC](zdir, dim);
      dZdt_C += jac_sps[spC](zdir, dim) * jac_sps[spC](zdir, dim);
      dzdt_D += jac_cur_sps[spD](zdir, dim) * jac_cur_sps[spD](zdir, dim);
      dZdt_D += jac_sps[spD](zdir, dim) * jac_sps[spD](zdir, dim);
      dzdt_E += jac_cur_sps[spE](zdir, dim) * jac_cur_sps[spE](zdir, dim);
      dZdt_E += jac_sps[spE](zdir, dim) * jac_sps[spE](zdir, dim);
    }
    // E33: remedy of curvature thickness locking
    // Ett = 0.5* ( (1-r-s) * Ett(SP C) + r * Ett(SP D) + s * Ett(SP E) )
    lstrain(2) =
        0.5 * ((1 - r - s) * (dzdt_C - dZdt_C) + r * (dzdt_D - dZdt_D) + s * (dzdt_E - dZdt_E));
    // E23: remedy of transverse shear locking
    // Est = r * Est(SP B)
    lstrain(4) = r * (dxdt_B - dXdt_B);
    // E13: remedy of transverse shear locking
    // Ert = s * Est(SP A)
    lstrain(5) = s * (dydt_A - dYdt_A);
    // ANS modification of strains ************************************** ANS

    // transformation of local glstrains 'back' to global(material) space
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1> glstrain(true);
    glstrain.Multiply(TinvT, lstrain);

    // EAS technology: "enhance the strains"  ----------------------------- EAS
    if (eastype_ == soshw6_easpoisthick)
    {
      // map local M to global, also enhancement is refered to element origin
      // M = detJ0/detJ T0^{-T} . M
      Core::LinAlg::DenseFunctions::multiply<double, Mat::NUM_STRESS_3D, Mat::NUM_STRESS_3D,
          soshw6_easpoisthick>(M.A(), detJ0 / detJ, T0invT.A(), M_GP->at(gp).values());
      // add enhanced strains = M . alpha to GL strains to "unlock" element
      Core::LinAlg::DenseFunctions::multiply<double, Mat::NUM_STRESS_3D, soshw6_easpoisthick, 1>(
          1.0, glstrain.A(), 1.0, M.A(), (*alpha).values());
    }  // ------------------------------------------------------------------ EAS

    // return gp GL strains (only possible option) if necessary
    switch (iostrain)
    {
      case Inpar::STR::strain_gl:
      {
        if (elestress == nullptr) FOUR_C_THROW("no strain data available");
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
        FOUR_C_THROW("no Euler-Almansi strains available for soshw6");
        break;
      case Inpar::STR::strain_none:
        break;
      default:
        FOUR_C_THROW("requested strain type not available");
    }

    // call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    /* Caution!! the defgrd can not be modified with ANS to remedy locking
       To get the consistent F a spectral decomposition would be necessary, see sosh8_cauchy.
       However if one only maps e.g. stresses from current to material configuration,
       I have never noticed any difference to applying just the disp_based F
       which is therefore computed and passed here (no significant add. computation time).  */
    Core::LinAlg::Matrix<NUMDIM_WEG6, NUMDIM_WEG6> defgrd;  // Caution!! disp_based!
    Core::LinAlg::Matrix<NUMDIM_WEG6, NUMNOD_WEG6> N_XYZ;
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJ_[gp], derivs[gp]);
    // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    defgrd.MultiplyTT(xcurr, N_XYZ);
    //
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, Mat::NUM_STRESS_3D> cmat(true);
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1> stress(true);

    SolidMaterial()->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, gp, Id());
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp stresses
    switch (iostress)
    {
      case Inpar::STR::stress_2pk:
      {
        if (elestress == nullptr) FOUR_C_THROW("no stress data available");
        for (int i = 0; i < Mat::NUM_STRESS_3D; ++i)
        {
          (*elestress)(gp, i) = stress(i);
        }
      }
      break;
      case Inpar::STR::stress_cauchy:
      {
        if (elestress == nullptr) FOUR_C_THROW("stress data not available");
        soshw6_cauchy(elestress, gp, defgrd, glstrain, stress);
      }
      break;
      case Inpar::STR::stress_none:
        break;
      default:
        FOUR_C_THROW("requested stress type not available");
    }

    double detJ_w = detJ * gpweights[gp];
    // update internal force vector
    if (force != nullptr)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w, bop, stress, 1.0);
    }
    if (split_res) force_str->MultiplyTN(detJ_w, bop, stress, 1.0);

    // update stiffness matrix
    if (stiffmatrix != nullptr)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, NUMDOF_WEG6> cb;
      cb.Multiply(cmat, bop);  // temporary C . B
      stiffmatrix->MultiplyTN(detJ_w, bop, cb, 1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      // here also the ANS interpolation comes into play
      for (int inod = 0; inod < NUMNOD_WEG6; ++inod)
      {
        for (int jnod = 0; jnod < NUMNOD_WEG6; ++jnod)
        {
          Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1> G_ij;
          G_ij(0) = derivs[gp](0, inod) * derivs[gp](0, jnod);  // rr-dir
          G_ij(1) = derivs[gp](1, inod) * derivs[gp](1, jnod);  // ss-dir
          G_ij(3) = derivs[gp](0, inod) * derivs[gp](1, jnod) +
                    derivs[gp](1, inod) * derivs[gp](0, jnod);  // rs-dir

          //          // testing without ANS:
          //          G_ij(2) = derivs[gp](2, inod) * derivs[gp](2, jnod); // tt-dir
          //          G_ij(4) = derivs[gp](1, inod) * derivs[gp](2, jnod)
          //                  + derivs[gp](2, inod) * derivs[gp](1, jnod); // st-dir
          //          G_ij(5) = derivs[gp](0, inod) * derivs[gp](2, jnod)
          //                  + derivs[gp](2, inod) * derivs[gp](0, jnod); // rt-dir

          // ANS modification in tt-dir
          G_ij(2) = (1 - r - s) * (*deriv_sp)[spC](zdir, inod) * (*deriv_sp)[spC](zdir, jnod) +
                    r * (*deriv_sp)[spD](zdir, inod) * (*deriv_sp)[spD](zdir, jnod) +
                    s * (*deriv_sp)[spE](zdir, inod) * (*deriv_sp)[spE](zdir, jnod);
          // ANS modification in st-dir
          G_ij(4) = r * ((*deriv_sp)[spB](ydir, inod) * (*deriv_sp)[spB](zdir, jnod) +
                            (*deriv_sp)[spB](zdir, inod) * (*deriv_sp)[spB](ydir, jnod));
          // ANS modification in rt-dir
          G_ij(5) = s * ((*deriv_sp)[spA](xdir, inod) * (*deriv_sp)[spA](zdir, jnod) +
                            (*deriv_sp)[spA](zdir, inod) * (*deriv_sp)[spA](xdir, jnod));
          // transformation of local(parameter) space 'back' to global(material) space
          Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1> G_ij_glob;
          G_ij_glob.Multiply(TinvT, G_ij);

          // Scalar Gij results from product of G_ij with stress, scaled with detJ*weights
          double Gij = detJ_w * stress.Dot(G_ij_glob);

          // add "geometric part" Gij times detJ*weights to stiffness matrix
          (*stiffmatrix)(NUMDIM_WEG6 * inod + 0, NUMDIM_WEG6 * jnod + 0) += Gij;
          (*stiffmatrix)(NUMDIM_WEG6 * inod + 1, NUMDIM_WEG6 * jnod + 1) += Gij;
          (*stiffmatrix)(NUMDIM_WEG6 * inod + 2, NUMDIM_WEG6 * jnod + 2) += Gij;
        }
      }  // end of intergrate `geometric' stiffness ******************************

      // EAS technology: integrate matrices --------------------------------- EAS
      if (eastype_ == soshw6_easpoisthick)
      {
        // integrate Kaa: Kaa += (M^T . cmat . M) * detJ * w(gp)
        Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, soshw6_easpoisthick> cM;  // temporary c . M
        cM.Multiply(cmat, M);
        Core::LinAlg::DenseFunctions::multiplyTN<double, soshw6_easpoisthick, Mat::NUM_STRESS_3D,
            soshw6_easpoisthick>(1.0, Kaa.values(), detJ_w, M.A(), cM.A());

        // integrate Kda: Kda += (M^T . cmat . B) * detJ * w(gp)
        Core::LinAlg::DenseFunctions::multiplyTN<double, soshw6_easpoisthick, Mat::NUM_STRESS_3D,
            NUMDOF_WEG6>(1.0, Kda.values(), detJ_w, M.A(), cb.A());

        // integrate feas: feas += (M^T . sigma) * detJ *wp(gp)
        Core::LinAlg::DenseFunctions::multiplyTN<double, soshw6_easpoisthick, Mat::NUM_STRESS_3D,
            1>(1.0, feas.values(), detJ_w, M.A(), stress.A());
      }  // ------------------------------------------------------------------ EAS
    }

    if (massmatrix != nullptr)
    {  // evaluate mass matrix +++++++++++++++++++++++++
      double density = Material()->Density(gp);
      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      double ifactor, massfactor;
      for (int inod = 0; inod < NUMNOD_WEG6; ++inod)
      {
        ifactor = shapefcts[gp](inod) * factor;
        for (int jnod = 0; jnod < NUMNOD_WEG6; ++jnod)
        {
          massfactor = shapefcts[gp](jnod) * ifactor;  // intermediate factor
          (*massmatrix)(NUMDIM_WEG6 * inod + 0, NUMDIM_WEG6 * jnod + 0) += massfactor;
          (*massmatrix)(NUMDIM_WEG6 * inod + 1, NUMDIM_WEG6 * jnod + 1) += massfactor;
          (*massmatrix)(NUMDIM_WEG6 * inod + 2, NUMDIM_WEG6 * jnod + 2) += massfactor;
        }
      }
    }  // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
    /* =========================================================================*/
  } /* ==================================================== end of Loop over GP */
  /* =========================================================================*/

  // rhs norm of eas equations
  if (eastype_ != soshw6_easnone && split_res)
    // only add for row-map elements
    if (params.get<int>("MyPID") == Owner())
      params.get<double>("cond_rhs_norm") += pow(Core::LinAlg::Norm2(feas), 2.);

  if (force != nullptr && stiffmatrix != nullptr)
  {
    // EAS technology: ------------------------------------------------------ EAS
    // subtract EAS matrices from disp-based Kdd to "soften" element
    if (eastype_ == soshw6_easpoisthick)
    {
      // we need the inverse of Kaa
      using ordinalType = Core::LinAlg::SerialDenseMatrix::ordinalType;
      using scalarType = Core::LinAlg::SerialDenseMatrix::scalarType;
      Teuchos::SerialDenseSolver<ordinalType, scalarType> solve_for_inverseKaa;
      solve_for_inverseKaa.setMatrix(Teuchos::rcpFromRef(Kaa));
      solve_for_inverseKaa.invert();

      Core::LinAlg::SerialDenseMatrix KdaTKaa(
          NUMDOF_WEG6, soshw6_easpoisthick);  // temporary Kda^T.Kaa^{-1}
      Core::LinAlg::DenseFunctions::multiplyTN<double, NUMDOF_WEG6, soshw6_easpoisthick,
          soshw6_easpoisthick>(KdaTKaa, Kda, Kaa);

      // EAS-stiffness matrix is: Kdd - Kda^T . Kaa^-1 . Kda
      Core::LinAlg::DenseFunctions::multiply<double, NUMDOF_WEG6, soshw6_easpoisthick, NUMDOF_WEG6>(
          1.0, stiffmatrix->A(), -1.0, KdaTKaa.values(), Kda.values());

      // EAS-internal force is: fint - Kda^T . Kaa^-1 . feas
      Core::LinAlg::DenseFunctions::multiply<double, NUMDOF_WEG6, soshw6_easpoisthick, 1>(
          1.0, force->A(), -1.0, KdaTKaa.values(), feas.values());

      // store current EAS data in history
      for (int i = 0; i < soshw6_easpoisthick; ++i)
      {
        for (int j = 0; j < soshw6_easpoisthick; ++j) (*oldKaainv)(i, j) = Kaa(i, j);
        for (int j = 0; j < NUMDOF_WEG6; ++j) (*oldKda)(i, j) = Kda(i, j);
        (*oldfeas)(i, 0) = feas(i);
      }
    }  // -------------------------------------------------------------------- EAS
  }

  return;
}


/*----------------------------------------------------------------------*
 |  setup of constant ANS data (private)                       maf 05/07|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoShw6::soshw6_anssetup(
    const Core::LinAlg::Matrix<NUMNOD_WEG6, NUMDIM_WEG6>& xrefe,  // material element coords
    const Core::LinAlg::Matrix<NUMNOD_WEG6, NUMDIM_WEG6>& xcurr,  // current element coords
    std::vector<Core::LinAlg::Matrix<NUMDIM_WEG6, NUMNOD_WEG6>>**
        deriv_sp,  // derivs eval. at all sampling points
    std::vector<Core::LinAlg::Matrix<NUMDIM_WEG6, NUMDIM_WEG6>>&
        jac_sps,  // jac at all sampling points
    std::vector<Core::LinAlg::Matrix<NUMDIM_WEG6, NUMDIM_WEG6>>&
        jac_cur_sps,  // current jac at all sampling points
    Core::LinAlg::Matrix<num_ans * num_sp, NUMDOF_WEG6>& B_ans_loc)  // modified B
{
  // static matrix object of derivs at sampling points, kept in memory
  static std::vector<Core::LinAlg::Matrix<NUMDIM_WEG6, NUMNOD_WEG6>> df_sp(num_sp);
  static bool dfsp_eval;  // flag for re-evaluate everything

  if (dfsp_eval != 0)
  {                      // if true f,df already evaluated
    *deriv_sp = &df_sp;  // return adress of static object to target of pointer
  }
  else
  {
    /*====================================================================*/
    /* 6-node wedge Solid-Shell node topology
     * and location of sampling points A to E                             */
    /*--------------------------------------------------------------------*/
    /*
     *                             s
     *                   6        /
     *                 //||\\   /
     *      t        //  ||   \\
     *      ^      //    || /    \\
     *      |    //      E          \\
     *      |  //        ||            \\
     *      |//       /  ||               \\
     *      5===============================6
     *     ||      B      3                 ||
     *     ||    /      // \\               ||
     *     || /       //      \\            ||
     *   - C -  -  -// -  A  -  -\\ -  -   -D  ----> r
     *     ||     //                \\      ||
     *  /  ||   //                     \\   ||
     *     || //                          \\||
     *      1================================2
     *
     */
    /*====================================================================*/
    // (r,s,t) gp-locations of sampling points A,B,C,D,E
    // numsp = 5 here set explicitly to allow direct initializing
    std::array<double, 5> r = {0.5, 0.0, 0.0, 1.0, 0.0};
    std::array<double, 5> s = {0.0, 0.5, 0.0, 0.0, 1.0};
    std::array<double, 5> t = {0.0, 0.0, 0.0, 0.0, 0.0};

    // fill up df_sp w.r.t. rst directions (NUMDIM) at each sp
    for (int i = 0; i < num_sp; ++i)
    {
      Core::FE::shape_function_3D_deriv1(df_sp[i], r[i], s[i], t[i], Core::FE::CellType::wedge6);
    }

    // return adresses of just evaluated matrices
    *deriv_sp = &df_sp;  // return adress of static object to target of pointer
    dfsp_eval = true;    // now all arrays are filled statically
  }

  for (int sp = 0; sp < num_sp; ++sp)
  {
    // compute Jacobian matrix at all sampling points
    jac_sps[sp].Multiply(df_sp[sp], xrefe);
    // compute CURRENT Jacobian matrix at all sampling points
    jac_cur_sps[sp].Multiply(df_sp[sp], xcurr);
  }

  /*
  ** Compute modified B-operator in local(parametric) space,
  ** evaluated at all sampling points
  */
  // loop over each sampling point
  Core::LinAlg::Matrix<NUMDIM_WEG6, NUMDIM_WEG6> jac_cur;
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
    for (int inode = 0; inode < NUMNOD_WEG6; ++inode)
    {
      for (int dim = 0; dim < NUMDIM_WEG6; ++dim)
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

/*----------------------------------------------------------------------*
 |  evaluate 'T'-transformation matrix )                       maf 05/07|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoShw6::soshw6_evaluate_t(
    const Core::LinAlg::Matrix<NUMDIM_WEG6, NUMDIM_WEG6>& jac,
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
  int err = solve_for_inverseT.Invert();
  if ((err != 0) && (err2 != 0)) FOUR_C_THROW("Inversion of Tinv (Jacobian) failed");
  return;
}

/*----------------------------------------------------------------------*
 |  initialize EAS data (private)                              maf 05/08|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoShw6::soshw6_easinit()
{
  // EAS enhanced strain parameters at currently investigated load/time step
  Core::LinAlg::SerialDenseMatrix alpha(neas_, 1);
  // EAS enhanced strain parameters of last converged load/time step
  Core::LinAlg::SerialDenseMatrix alphao(neas_, 1);
  // EAS portion of internal forces, also called enhacement vector s or Rtilde
  Core::LinAlg::SerialDenseMatrix feas(neas_, 1);
  // EAS matrix K_{alpha alpha}, also called Dtilde
  Core::LinAlg::SerialDenseMatrix invKaa(neas_, neas_);
  // EAS matrix K_{d alpha}
  Core::LinAlg::SerialDenseMatrix Kda(neas_, NUMDOF_WEG6);
  // EAS increment
  Core::LinAlg::SerialDenseMatrix eas_inc(neas_, 1);

  // save EAS data into element container
  easdata_.alpha = alpha;
  easdata_.alphao = alphao;
  easdata_.feas = feas;
  easdata_.invKaa = invKaa;
  easdata_.Kda = Kda;
  easdata_.eas_inc = eas_inc;
}

/*----------------------------------------------------------------------*
 |  setup of constant EAS data (private)                       maf 05/07|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoShw6::soshw6_eassetup(
    std::vector<Core::LinAlg::SerialDenseMatrix>** M_GP,  // M-matrix evaluated at GPs
    double& detJ0,                                        // det of Jacobian at origin
    Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, Mat::NUM_STRESS_3D>&
        T0invT,                                                   // maps M(origin) local to global
    const Core::LinAlg::Matrix<NUMNOD_WEG6, NUMDIM_WEG6>& xrefe)  // material element coords
{
  // shape function derivatives, evaluated at origin (r=s=t=0.0)
  const Core::FE::IntegrationPoints3D intpoints(Core::FE::GaussRule3D::wedge_1point);
  Core::LinAlg::Matrix<NUMDIM_WEG6, NUMNOD_WEG6> df0;
  Core::FE::shape_function_3D_deriv1(df0, intpoints.qxg[0][0], intpoints.qxg[0][1],
      intpoints.qxg[0][2], Core::FE::CellType::wedge6);

  // compute Jacobian, evaluated at element origin (r=s=t=0.0)
  Core::LinAlg::Matrix<NUMDIM_WEG6, NUMDIM_WEG6> jac0;
  jac0.Multiply(df0, xrefe);

  // compute determinant of Jacobian at origin by Sarrus' rule
  detJ0 = jac0.Determinant();

  // get T-matrix at element origin
  soshw6_evaluate_t(jac0, T0invT);

  // build EAS interpolation matrix M, evaluated at the GPs of soshw6
  static std::vector<Core::LinAlg::SerialDenseMatrix> M(NUMGPT_WEG6);
  static bool M_eval;

  if (M_eval == true)
  {              // if true M already evaluated
    *M_GP = &M;  // return adress of static object to target of pointer
    return;
  }
  else
  {
    // (r,s,t) gp-locations of fully integrated linear 6-node Wedge
    const Core::FE::IntegrationPoints3D intpoints(Core::FE::GaussRule3D::wedge_6point);

    // fill up M at each gp
    if (eastype_ == soshw6_easpoisthick)
    {
      /* EAS interpolation of 1 mode corr. to linear thickness strain
      **            0
      **            0
      **    M =     t
      **            0
      **            0
      **            0
      */
      for (int i = 0; i < intpoints.nquad; ++i)
      {
        M[i].shape(Mat::NUM_STRESS_3D, soshw6_easpoisthick);
        M[i](2, 0) = intpoints.qxg[i][2];  // t at gp
        // M[i](2,1) = intpoints.qxg[i][0]*intpoints.qxg[i][2];  // r*t at gp ->not activated at all
        // due to tri M[i](2,2) = intpoints.qxg[i][1]*intpoints.qxg[i][2];  // s*t at gp ->not
        // activated at all due to tri
      }

      // return adress of just evaluated matrix
      *M_GP = &M;     // return adress of static object to target of pointer
      M_eval = true;  // now the array is filled statically
    }
    else
    {
      FOUR_C_THROW("eastype not implemented");
    }
  }
}

/*----------------------------------------------------------------------*
 |  return Cauchy stress at gp                                 maf 06/08|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoShw6::soshw6_cauchy(
    Core::LinAlg::Matrix<NUMGPT_WEG6, Mat::NUM_STRESS_3D>* elestress, const int gp,
    const Core::LinAlg::Matrix<NUMDIM_WEG6, NUMDIM_WEG6>& defgrd,
    const Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1>& glstrain,
    const Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1>& stress)
{
#if consistent_F
  // double disp1 = defgrd.NormOne();
  // double dispinf = defgrd.NormInf();

  /* to get the consistent (locking-free) F^mod, we need two spectral
   * compositions. First, find R (rotation tensor) from F=RU,
   * then from E^mod = 1/2((U^mod)^2 - 1) find U^mod,
   * and finally F^mod = RU^mod */

  // polar decomposition of displacement based F
  Core::LinAlg::SerialDenseMatrix u(NUMDIM_WEG6, NUMDIM_WEG6);
  Core::LinAlg::SerialDenseMatrix s(NUMDIM_WEG6, NUMDIM_WEG6);
  Core::LinAlg::SerialDenseMatrix v(NUMDIM_WEG6, NUMDIM_WEG6);
  SVD(defgrd, u, s, v);  // Singular Value Decomposition
  Core::LinAlg::SerialDenseMatrix rot(NUMDIM_WEG6, NUMDIM_WEG6);
  Core::LinAlg::multiply(rot, u, v);

  // get modified squared stretch (U^mod)^2 from glstrain
  Core::LinAlg::SerialDenseMatrix Usq_mod(NUMDIM_WEG6, NUMDIM_WEG6);
  for (int i = 0; i < NUMDIM_WEG6; ++i) Usq_mod(i, i) = 2.0 * glstrain(i) + 1.0;
  // off-diagonal terms are already twice in the Voigt-GLstrain-vector
  Usq_mod(0, 1) = glstrain(3);
  Usq_mod(1, 0) = glstrain(3);
  Usq_mod(1, 2) = glstrain(4);
  Usq_mod(2, 1) = glstrain(4);
  Usq_mod(0, 2) = glstrain(5);
  Usq_mod(2, 0) = glstrain(5);
  // polar decomposition of (U^mod)^2
  SVD(Usq_mod, u, s, v);  // Singular Value Decomposition
  Core::LinAlg::SerialDenseMatrix U_mod(NUMDIM_WEG6, NUMDIM_WEG6);
  for (int i = 0; i < NUMDIM_SOH8; ++i) s(i, i) = sqrt(s(i, i));
  Core::LinAlg::SerialDenseMatrix temp2(NUMDIM_WEG6, NUMDIM_WEG6);
  Core::LinAlg::multiply(temp2, u, s);
  Core::LinAlg::multiply(U_mod, temp2, v);

  // F^mod = RU^mod
  Core::LinAlg::SerialDenseMatrix defgrd_consistent(NUMDIM_WEG6, NUMDIM_WEG6);
  Core::LinAlg::multiply(defgrd_consistent, rot, U_mod);
  defgrd.SetView(defgrd_consistent.A());
#endif
  double detF = defgrd.Determinant();

  Core::LinAlg::Matrix<NUMDIM_WEG6, NUMDIM_WEG6> pkstress;
  pkstress(0, 0) = stress(0);
  pkstress(0, 1) = stress(3);
  pkstress(0, 2) = stress(5);
  pkstress(1, 0) = pkstress(0, 1);
  pkstress(1, 1) = stress(1);
  pkstress(1, 2) = stress(4);
  pkstress(2, 0) = pkstress(0, 2);
  pkstress(2, 1) = pkstress(1, 2);
  pkstress(2, 2) = stress(2);

  Core::LinAlg::Matrix<NUMDIM_WEG6, NUMDIM_WEG6> cauchystress;
  Core::LinAlg::Matrix<NUMDIM_WEG6, NUMDIM_WEG6> temp;
  temp.Multiply(1.0 / detF, defgrd, pkstress);
  cauchystress.MultiplyNT(temp, defgrd);

  (*elestress)(gp, 0) = cauchystress(0, 0);
  (*elestress)(gp, 1) = cauchystress(1, 1);
  (*elestress)(gp, 2) = cauchystress(2, 2);
  (*elestress)(gp, 3) = cauchystress(0, 1);
  (*elestress)(gp, 4) = cauchystress(1, 2);
  (*elestress)(gp, 5) = cauchystress(0, 2);
  return;
}


/*----------------------------------------------------------------------*
 | find optimal map between material space and parameter space maf 11/08|
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoShw6::soshw6_findoptparmap()
{
  // create edge vectors of lower triangle
  std::vector<std::vector<double>> edgevecs(3);
  edgevecs.at(0).resize(3);
  edgevecs.at(1).resize(3);
  edgevecs.at(2).resize(3);
  for (int i = 0; i < 3; ++i)
  {
    edgevecs.at(0)[i] = this->Nodes()[1]->X()[i] - this->Nodes()[0]->X()[i];  // a
    edgevecs.at(1)[i] = this->Nodes()[2]->X()[i] - this->Nodes()[0]->X()[i];  // b
    edgevecs.at(2)[i] = this->Nodes()[2]->X()[i] - this->Nodes()[1]->X()[i];  // c
  }

  // normalize
  for (int i = 0; i < 3; ++i)
  {
    double d = 0.;
    for (int j = 0; j < 3; ++j) d += edgevecs.at(i)[j] * edgevecs.at(i)[j];
    for (int j = 0; j < 3; ++j) edgevecs.at(i)[j] = edgevecs.at(i)[j] / sqrt(d);
  }

  // scalar product of edge vectors
  double ab = 0.;
  double bc = 0.;
  double ac = 0.;
  for (int i = 0; i < 3; ++i)
  {
    ab += edgevecs.at(0)[i] * edgevecs.at(1)[i];
    bc += edgevecs.at(1)[i] * edgevecs.at(2)[i];
    ac += edgevecs.at(0)[i] * edgevecs.at(2)[i];
  }

  // find min(ab,bc,ac)
  if ((abs(ab) <= abs(bc)) && (abs(ab) <= abs(ac)))
    return 1;  // = ab
  else if (abs(bc) <= abs(ac))
    return 3;  // = bc
  else
    return 2;  // =ac

  // impossible case
  FOUR_C_THROW("Could not find optimal map!");
  return 0;
}


/*----------------------------------------------------------------------*
 |  init the element (public)                                  maf 11/08|
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::SoShw6Type::Initialize(Discret::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<Discret::ELEMENTS::SoShw6*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_shw6* failed");

    // check whether we should align the material space optimally with the parameter space.
    // The idea is that elimination of shear-locking works best if the origin of the
    // triangle-parameter space coincides with the node where the angle between the edges
    // is closest to 90 degree.
    if ((actele->optimal_parameterspace_map_) && (!actele->nodes_rearranged_))
    {
      int originnode = actele->soshw6_findoptparmap();

      int new_nodeids[NUMNOD_WEG6];

      switch (originnode)
      {
        case 1:
        {
          // no resorting necessary, already aligned with FIRST node
          actele->nodes_rearranged_ = true;
          break;
        }
        case 2:
        {
          // resorting of nodes to align with SECOND node
          new_nodeids[0] = actele->NodeIds()[1];
          new_nodeids[1] = actele->NodeIds()[2];
          new_nodeids[2] = actele->NodeIds()[0];
          new_nodeids[3] = actele->NodeIds()[4];
          new_nodeids[4] = actele->NodeIds()[5];
          new_nodeids[5] = actele->NodeIds()[3];
          actele->SetNodeIds(NUMNOD_WEG6, new_nodeids);
          actele->nodes_rearranged_ = true;
          break;
        }
        case 3:
        {
          // resorting of nodes to align with THIRD node
          new_nodeids[0] = actele->NodeIds()[2];
          new_nodeids[1] = actele->NodeIds()[0];
          new_nodeids[2] = actele->NodeIds()[1];
          new_nodeids[3] = actele->NodeIds()[5];
          new_nodeids[4] = actele->NodeIds()[3];
          new_nodeids[5] = actele->NodeIds()[4];
          actele->SetNodeIds(NUMNOD_WEG6, new_nodeids);
          actele->nodes_rearranged_ = true;
          break;
        }
        default:
          FOUR_C_THROW("reordering of So_shw6 failed");
      }
    }
  }
  // fill complete again to reconstruct element-node pointers,
  // but without element init, etc.
  dis.fill_complete(false, false, false);

  // loop again to init Jacobian for Sosh8's
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<Discret::ELEMENTS::SoShw6*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_shw6* failed");
    actele->init_jacobian_mapping();
  }
  return 0;
}

void Discret::ELEMENTS::SoShw6::soshw6_recover(const std::vector<double>& residual)
{
  if (eastype_ == soshw6_easnone) return;

  Core::LinAlg::Matrix<NUMDOF_WEG6, 1> disi(false);
  for (int i = 0; i < NUMDOF_WEG6; ++i) disi(i) = residual[i];

  const double step_length = str_params_interface().get_step_length();

  auto* oldfeas = &easdata_.feas;
  auto* oldKda = &easdata_.Kda;
  auto* alpha = &easdata_.alpha;
  auto* eas_inc = &easdata_.eas_inc;
  auto* oldKaainv = &easdata_.invKaa;
  /* if it is a default step, we have to recover the condensed
   * solution vectors */
  if (str_params_interface().is_default_step())
  {
    // first, store the eas state of the previous accepted Newton step
    str_params_interface().sum_into_my_previous_sol_norm(
        NOX::Nln::StatusTest::quantity_eas, soshw6_easpoisthick, (*alpha)[0], Owner());

    // add Kda . res_d to feas
    Core::LinAlg::DenseFunctions::multiply<double, soshw6_easpoisthick, NUMDOF_WEG6, 1>(
        1.0, oldfeas->values(), 1.0, oldKda->values(), residual.data());
    // "new" alpha is: - Kaa^-1 . (feas + Kda . old_d), here: - Kaa^-1 . feas
    Core::LinAlg::DenseFunctions::multiply<double, soshw6_easpoisthick, soshw6_easpoisthick, 1>(
        0.0, *eas_inc, -1.0, *oldKaainv, *oldfeas);
    Core::LinAlg::DenseFunctions::update<double, soshw6_easpoisthick, 1>(1., *alpha, 1., *eas_inc);
  }
  else
  {
    Core::LinAlg::DenseFunctions::update<double, soshw6_easpoisthick, 1>(
        1., *alpha, step_length - old_step_length_, *eas_inc);
  }
  old_step_length_ = step_length;

  str_params_interface().sum_into_my_update_norm(NOX::Nln::StatusTest::quantity_eas,
      soshw6_easpoisthick, (*eas_inc)[0], (*alpha)[0], step_length, Owner());
}

FOUR_C_NAMESPACE_CLOSE

/*----------------------------------------------------------------------*/
/*! \file
\brief 'Q1P0' element in 8-node hexahedron shape

\level 2

\maintainer Christoph Meier
*/
/*----------------------------------------------------------------------*/

#include "so_hex8p1j1.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_mat/so3_material.H"

#include <Epetra_SerialComm.h>


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_Hex8P1J1::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  // Check whether the solid material PostSetup() routine has already been called and call it if not
  CheckMaterialPostSetup(params);

  // get parameter interface
  SetParamsInterfacePtr(params);

  LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8> elemat1(elemat1_epetra.A(), true);
  LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8> elemat2(elemat2_epetra.A(), true);
  LINALG::Matrix<NUMDOF_SOH8, 1> elevec1(elevec1_epetra.A(), true);
  LINALG::Matrix<NUMDOF_SOH8, 1> elevec2(elevec2_epetra.A(), true);
  LINALG::Matrix<NUMDOF_SOH8, 1> elevec3(elevec3_epetra.A(), true);

  // start with "none"
  DRT::ELEMENTS::So_hex8::ActionType act = So_hex8::none;

  // get the required action
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    dserror("No action supplied");
  else if (action == "calc_struct_linstiff")
    act = So_hex8::calc_struct_linstiff;
  else if (action == "calc_struct_nlnstiff")
    act = So_hex8::calc_struct_nlnstiff;
  else if (action == "calc_struct_internalforce")
    act = So_hex8::calc_struct_internalforce;
  else if (action == "calc_struct_linstiffmass")
    act = So_hex8::calc_struct_linstiffmass;
  else if (action == "calc_struct_nlnstiffmass")
    act = So_hex8::calc_struct_nlnstiffmass;
  else if (action == "calc_struct_nlnstifflmass")
    act = So_hex8::calc_struct_nlnstifflmass;
  else if (action == "calc_struct_stress")
    act = So_hex8::calc_struct_stress;
  else if (action == "calc_struct_eleload")
    act = So_hex8::calc_struct_eleload;
  else if (action == "calc_struct_fsiload")
    act = So_hex8::calc_struct_fsiload;
  else if (action == "calc_struct_update_istep")
    act = So_hex8::calc_struct_update_istep;
  else if (action == "calc_struct_reset_istep")
    act = So_hex8::calc_struct_reset_istep;
  else if (action == "postprocess_stress")
    act = So_hex8::postprocess_stress;
  else if (action == "multi_eas_init")
    act = So_hex8::multi_eas_init;
  else if (action == "multi_eas_set")
    act = So_hex8::multi_eas_set;
  else if (action == "multi_calc_dens")
    act = So_hex8::multi_calc_dens;
  else if (action == "multi_readrestart")
    act = So_hex8::multi_readrestart;
  else if (action == "calc_struct_recover")
    act = So_hex8::calc_recover;
  else if (action == "calc_struct_predict")
    return 0;
  else
    dserror("Unknown type of action for So_hex8");

  // what should the element do
  switch (act)
  {
    // linear business
    case calc_struct_linstiffmass:
    case calc_struct_linstiff:
      dserror("Solid Q1P0 hex8 element does not implement any linear stiffness calculation");
      break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (disp == Teuchos::null || res == Teuchos::null)
        dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);

      ForceStiffMass(lm, mydisp, myres, &elemat1, NULL, &elevec1, &elevec3, NULL, NULL, params,
          INPAR::STR::stress_none, INPAR::STR::strain_none);
    }
    break;

    // internal force vector only
    case calc_struct_internalforce:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (disp == Teuchos::null || res == Teuchos::null)
        dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8> myemat(true);

      ForceStiffMass(lm, mydisp, myres, &myemat, NULL, &elevec1, NULL, NULL, NULL, params,
          INPAR::STR::stress_none, INPAR::STR::strain_none);
    }
    break;

    // nonlinear stiffness, internal force vector, and consistent/lumped mass matrix
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (disp == Teuchos::null || res == Teuchos::null)
        dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);

      ForceStiffMass(lm, mydisp, myres, &elemat1, &elemat2, &elevec1, &elevec3, NULL, NULL, params,
          INPAR::STR::stress_none, INPAR::STR::strain_none);

      // lump mass
      if (act == calc_struct_nlnstifflmass) soh8_lumpmass(&elemat2);
    }
    break;

    // evaluate stresses and strains at gauss points
    case calc_struct_stress:
    {
      dserror(
          "Solid Q1P0 hex8 element does not yet implement writing/post-processing "
          "stresses/strains");
      // nothing to do for ghost elements
      if (discretization.Comm().MyPID() == Owner())
      {
        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
        Teuchos::RCP<std::vector<char>> stressdata =
            params.get<Teuchos::RCP<std::vector<char>>>("stress", Teuchos::null);
        Teuchos::RCP<std::vector<char>> straindata =
            params.get<Teuchos::RCP<std::vector<char>>>("strain", Teuchos::null);
        if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
        if (stressdata == Teuchos::null) dserror("Cannot get 'stress' data");
        if (straindata == Teuchos::null) dserror("Cannot get 'strain' data");
        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
        std::vector<double> myres(lm.size());
        DRT::UTILS::ExtractMyValues(*res, myres, lm);
        LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D> stress;
        LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D> strain;
        INPAR::STR::StressType iostress =
            DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
        INPAR::STR::StrainType iostrain =
            DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);
        ForceStiffMass(lm, mydisp, myres, NULL, NULL, NULL, NULL, &stress, &strain, params,
            iostress, iostrain);
        {
          DRT::PackBuffer data;
          AddtoPack(data, stress);
          data.StartPacking();
          AddtoPack(data, stress);
          std::copy(data().begin(), data().end(), std::back_inserter(*stressdata));
        }
        {
          DRT::PackBuffer data;
          AddtoPack(data, strain);
          data.StartPacking();
          AddtoPack(data, strain);
          std::copy(data().begin(), data().end(), std::back_inserter(*straindata));
        }
      }
    }
    break;

    // postprocess stresses/strains at gauss points
    //
    // note that in the following, quantities are always referred to as
    // "stresses" etc. although they might also apply to strains
    // (depending on what this routine is called for from the post filter)
    case postprocess_stress:
    {
      const Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> gpstressmap =
          params.get<Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>>>(
              "gpstressmap", Teuchos::null);
      if (gpstressmap == Teuchos::null)
        dserror("no gp stress/strain map available for postprocessing");
      std::string stresstype = params.get<std::string>("stresstype", "ndxyz");
      const int gid = Id();
      LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D> gpstress(((*gpstressmap)[gid])->A(), true);

      Teuchos::RCP<Epetra_MultiVector> poststress =
          params.get<Teuchos::RCP<Epetra_MultiVector>>("poststress", Teuchos::null);
      if (poststress == Teuchos::null) dserror("No element stress/strain vector available");

      if (stresstype == "ndxyz")
      {
        // extrapolate stresses/strains at Gauss points to nodes
        soh8_expol(gpstress, *poststress);
      }
      else if (stresstype == "cxyz")
      {
        const Epetra_BlockMap& elemap = poststress->Map();
        const int lid = elemap.LID(Id());
        if (lid != -1)
        {
          for (int i = 0; i < MAT::NUM_STRESS_3D; ++i)
          {
            double& s = (*((*poststress)(i)))[lid];  // resolve pointer for faster access
            s = 0.;
            for (unsigned j = 0; j < NUMGPT_SOH8; ++j)
            {
              s += gpstress(j, i);
            }
            s *= 1.0 / NUMGPT_SOH8;
          }
        }
        else
        {
          dserror("unknown type of stress/strain output on element level");
        }
      }
    }
    break;

    case calc_struct_eleload:
      dserror("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
      break;

    case calc_struct_fsiload:
      dserror("Case not yet implemented");
      break;

    case calc_struct_update_istep:
    {
      p_o_(0, 0) = p_(0, 0);
      t_o_(0, 0) = t_(0, 0);

      // Update of history for materials
      SolidMaterial()->Update();
    }
    break;

    case calc_struct_reset_istep:
    {
      // Reset of history (if needed)
      SolidMaterial()->ResetStep();
    }
    break;

    case multi_eas_init:
    case multi_eas_set:
    case multi_calc_dens:
    case multi_readrestart:
      dserror("multi-scale stuff not implemented for solid Q1P0 hex8 element");
      break;

    case calc_recover:
    {
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);
      soh8P1J1_recover(myres);
      break;
    }

    default:
      dserror("Unknown type of action for So_Hex8P1J1");
  }
  return 0;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_Hex8P1J1::ForceStiffMass(const std::vector<int>& lm,  // location matrix
    const std::vector<double>& disp,                             // current displacements
    const std::vector<double>& residual,                         // current residual displ
    LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* stiffmatrix,       // element stiffness matrix
    LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* massmatrix,        // element mass matrix
    LINALG::Matrix<NUMDOF_SOH8, 1>* force,                       // element internal force vector
    LINALG::Matrix<NUMDOF_SOH8, 1>* force_str,                   // structure force
    LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestress,  // stresses at GP
    LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestrain,  // strains at GP
    Teuchos::ParameterList& params,                              // algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,                       // stress output option
    const INPAR::STR::StrainType iostrain                        // strain output option
)
{
  /* ============================================================================*
  ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
  ** ============================================================================*/
  const static std::vector<LINALG::Matrix<NUMNOD_SOH8, 1>> shapefcts = soh8_shapefcts();
  const static std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = soh8_derivs();
  const static std::vector<double> gpweights = soh8_weights();
  /* ============================================================================*/

  // update element geometry
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe;  // material coord. of element
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xcurr;  // current  coord. of element
  DRT::Node** nodes = Nodes();
  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];

    xcurr(i, 0) = xrefe(i, 0) + disp[i * NODDOF_SOH8 + 0];
    xcurr(i, 1) = xrefe(i, 1) + disp[i * NODDOF_SOH8 + 1];
    xcurr(i, 2) = xrefe(i, 2) + disp[i * NODDOF_SOH8 + 2];
  }

  size_t length = residual.size();
  if (length != 24) dserror("number of element displacement dofs not equal 24");

  LINALG::Matrix<24, 1> disi(false);
  for (int i = 0; i < NUMDOF_SOH8; ++i)
  {
    disi(i, 0) = residual[i];
  }

  // check if we need to split the residuals (for Newton line search)
  // if true an additional global vector is assembled containing
  // the internal forces without the condensed EAS entries and the norm
  // of the EAS residual is calculated
  bool split_res = params.isParameter("cond_rhs_norm");

  /* ============================================================================*
  ** UPDATE OF ELEMENT VOLUME AND PRESSURE                                       *
  ** ============================================================================*/

  // this is a line search step, i.e. the direction of the eas increments
  // has been calculated by a Newton step and now it is only scaled
  if (not IsParamsInterface())
  {
    if (params.isParameter("alpha_ls"))
    {
      double alpha_ls = params.get<double>("alpha_ls");
      // undo step
      t_ -= dt_;
      p_ -= dp_;
      // scale increment
      dt_.Scale(alpha_ls);
      dp_.Scale(alpha_ls);
      // add reduced increment
      t_ += dt_;
      p_ += dp_;
    }
    else
    {
      // dt= K_tp_^-1*(-R_p_-K_pu_*du)
      dt_.MultiplyNN(-1.0 / K_pt_, K_pu_, disi);
      dt_.Update(-1.0 / K_pt_, R_p_, 1.0);

      // dp= K_tp_^-1 * (-R_t_ - K_tu_*du - K_tt*dt)
      dp_.MultiplyNN(1.0, K_tu_, disi);
      dp_.Update(K_tt_, dt_, 1.0);
      dp_.Update(1.0, R_t_, 1.0);
      dp_.Scale(-1.0 / K_pt_);

      t_ += dt_;
      p_ += dp_;
    }
  }
  K_tt_ = 0.0;
  K_pu_.PutScalar(0.0);
  K_tu_.PutScalar(0.0);
  R_t_.PutScalar(0.0);
  R_p_.PutScalar(0.0);

  //******************************************************************************************

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/

  LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrd(false);
  LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> mod_defgrd(true);
  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
  {
    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJ_[gp], derivs[gp]);
    double detJ = detJ_[gp];

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    defgrd.MultiplyTT(xcurr, N_XYZ);

    // modified deformation gradient modF = (t_/J)^1/3 * F
    const double J = defgrd.Determinant();
    // check for negative jacobian
    if ((t_(0, 0) / J) <= 0.)
    {
      // check, if errors are tolerated or should throw a dserror
      bool error_tol = false;
      if (params.isParameter("tolerate_errors")) error_tol = params.get<bool>("tolerate_errors");
      if (error_tol)
      {
        params.set<bool>("eval_error", true);
        stiffmatrix->Clear();
        force->Clear();
        return;
      }
      else
        dserror("negative jacobian determinant");
    }
    const double scalar = std::pow((t_(0, 0) / J), 1.0 / 3.0);
    mod_defgrd.Update(scalar, defgrd);

    // Modified Right Cauchy-Green tensor = modF^T * modF
    LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> mod_cauchygreen;
    mod_cauchygreen.MultiplyTN(mod_defgrd, mod_defgrd);

    // Modified Green-Lagrange strains matrix mod_E = 0.5 * (modCauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector mod_glstrain_epetra(MAT::NUM_STRESS_3D);
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> mod_glstrain(mod_glstrain_epetra.A(), true);
    mod_glstrain(0) = 0.5 * (mod_cauchygreen(0, 0) - 1.0);
    mod_glstrain(1) = 0.5 * (mod_cauchygreen(1, 1) - 1.0);
    mod_glstrain(2) = 0.5 * (mod_cauchygreen(2, 2) - 1.0);
    mod_glstrain(3) = mod_cauchygreen(0, 1);
    mod_glstrain(4) = mod_cauchygreen(1, 2);
    mod_glstrain(5) = mod_cauchygreen(2, 0);

    // return gp strains if necessary
    if (iostrain != INPAR::STR::strain_none)
      Strain(elestrain, iostrain, gp, t_(0, 0), mod_defgrd, mod_glstrain);

    // Bestimmen der N_xyz
    LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_xyz(true);
    LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> Inv_defgrd(true);
    Inv_defgrd.Invert(defgrd);

    for (int i = 0; i < NUMNOD_SOH8; ++i)
    {
      for (int j = 0; j < NUMDIM_SOH8; ++j)
      {
        for (int k = 0; k < NUMDIM_SOH8; ++k)
        {
          N_xyz(j, i) += N_XYZ(k, i) * Inv_defgrd(k, j);
        }
      }
    }

    //==========================================================================
    // ***********************    modified B-matrix   **************************
    //==========================================================================
    LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8> bopn;
    for (int i = 0; i < NUMNOD_SOH8; ++i)
    {
      bopn(0, NODDOF_SOH8 * i + 0) = N_xyz(0, i);
      bopn(0, NODDOF_SOH8 * i + 1) = 0;
      bopn(0, NODDOF_SOH8 * i + 2) = 0;
      //------------------------------
      bopn(1, NODDOF_SOH8 * i + 0) = 0;
      bopn(1, NODDOF_SOH8 * i + 1) = N_xyz(1, i);
      bopn(1, NODDOF_SOH8 * i + 2) = 0;
      //-----------------------------
      bopn(2, NODDOF_SOH8 * i + 0) = 0;
      bopn(2, NODDOF_SOH8 * i + 1) = 0;
      bopn(2, NODDOF_SOH8 * i + 2) = N_xyz(2, i);
      //-----------------------------
      bopn(3, NODDOF_SOH8 * i + 0) = N_xyz(1, i);
      bopn(3, NODDOF_SOH8 * i + 1) = N_xyz(0, i);
      bopn(3, NODDOF_SOH8 * i + 2) = 0;
      //----------------------------
      bopn(4, NODDOF_SOH8 * i + 0) = 0;
      bopn(4, NODDOF_SOH8 * i + 1) = N_xyz(2, i);
      bopn(4, NODDOF_SOH8 * i + 2) = N_xyz(1, i);
      //----------------------------
      bopn(5, NODDOF_SOH8 * i + 0) = N_xyz(2, i);
      bopn(5, NODDOF_SOH8 * i + 1) = 0;
      bopn(5, NODDOF_SOH8 * i + 2) = N_xyz(0, i);
    }

    // call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> cmat(true);
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> stress(true);
    params.set<int>("gp", gp);
    SolidMaterial()->Evaluate(&mod_defgrd, &mod_glstrain, params, &stress, &cmat, Id());

    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    LINALG::Matrix<6, 6> D_T_bar(false);
    ConvertMat(cmat, mod_defgrd, D_T_bar, t_(0, 0));

    // (secondary) Cauchy stress dependent on displacements and primary Jacobian
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> sigma_bar(false);
    {
      // Voigt indices
      std::vector<int> Index1(MAT::NUM_STRESS_3D);
      std::vector<int> Index2(MAT::NUM_STRESS_3D);
      Index1[0] = 0;
      Index2[0] = 0;  // 11
      Index1[1] = 1;
      Index2[1] = 1;  // 22
      Index1[2] = 2;
      Index2[2] = 2;  // 33
      Index1[3] = 0;
      Index2[3] = 1;  // 12
      Index1[4] = 1;
      Index2[4] = 2;  // 23
      Index1[5] = 2;
      Index2[5] = 0;  // 31

      // build
      for (int i = 0; i < MAT::NUM_STRESS_3D; ++i)
      {
        sigma_bar(i, 0) = mod_defgrd(Index1[i], 0) * stress(0, 0) * mod_defgrd(Index2[i], 0) +
                          mod_defgrd(Index1[i], 0) * stress(3, 0) * mod_defgrd(Index2[i], 1) +
                          mod_defgrd(Index1[i], 0) * stress(5, 0) * mod_defgrd(Index2[i], 2) +
                          mod_defgrd(Index1[i], 1) * stress(3, 0) * mod_defgrd(Index2[i], 0) +
                          mod_defgrd(Index1[i], 1) * stress(1, 0) * mod_defgrd(Index2[i], 1) +
                          mod_defgrd(Index1[i], 1) * stress(4, 0) * mod_defgrd(Index2[i], 2) +
                          mod_defgrd(Index1[i], 2) * stress(5, 0) * mod_defgrd(Index2[i], 0) +
                          mod_defgrd(Index1[i], 2) * stress(4, 0) * mod_defgrd(Index2[i], 1) +
                          mod_defgrd(Index1[i], 2) * stress(2, 0) * mod_defgrd(Index2[i], 2);
      }
      sigma_bar.Scale(1.0 / t_(0, 0));
    }

    // secondary pressure
    const double p_bar = 1.0 / 3.0 * (sigma_bar(0, 0) + sigma_bar(1, 0) + sigma_bar(2, 0));

    // primary pressure with respect to initial configuration (?)
    const double p_hook = p_(0, 0) * J / t_(0, 0);



    // Cauchy stress dependent on displacements, primary Jacobian and primary pressure
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> sigma_hook(sigma_bar);
    for (int i = 0; i < 3; ++i)
    {
      sigma_hook(i, 0) += p_hook - p_bar;
    }

    // deviatoric content of Cauchy stress
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> sigma_bar_dev(sigma_bar);
    for (int i = 0; i < 3; ++i)
    {
      sigma_bar_dev(i, 0) -= p_bar;
    }

    // return GP stresses if necessary
    if (iostress != INPAR::STR::stress_none)
      Stress(elestress, iostress, gp, t_(0, 0), mod_defgrd, sigma_hook);

    // integration weights
    const double detJ_w = detJ * gpweights[gp];
    const double detJ_w_t = detJ_w * t_(0, 0);
    const double detJ_w_J = detJ_w * J;

    // update of internal force vector
    if (force != NULL)
    {
      // integrate internal force vector f = f + (B^T . sigma) * Theta * detJ * w(gp)
      force->MultiplyTN(detJ_w_t, bopn, sigma_hook, 1.0);
    }
    if (split_res) force_str->MultiplyTN(detJ_w_t, bopn, sigma_hook, 1.0);

    // update of stiffness matrix
    if (stiffmatrix != NULL)
    {
      R_p_(0, 0) += (J - t_(0, 0)) * detJ_w;
      R_t_(0, 0) += (p_bar - p_(0, 0)) * detJ_w;

      LINALG::Matrix<1, 1> D_22(false);
      {
        LINALG::Matrix<1, 6> temp(false);
        temp.MultiplyTN(1.0 / 9.0, m_, D_T_bar);
        D_22.MultiplyNN(temp, m_);
        D_22(0, 0) -= 1.0 / 3.0 * p_bar;
      }

      LINALG::Matrix<6, 6> D_11(false);
      {
        LINALG::Matrix<6, 6> temp2(false);
        temp2.Multiply(I_d_, D_T_bar);
        D_11.Multiply(temp2, I_d_);
        D_11.MultiplyNT(-2.0 / 3.0, m_, sigma_bar_dev, 1.0);
        D_11.MultiplyNT(-2.0 / 3.0, sigma_bar_dev, m_, 1.0);
        double scalar = 2.0 / 3.0 * p_bar - p_hook;
        D_11.MultiplyNT(-scalar, m_, m_, 1.0);
        scalar = 2.0 * (p_bar - p_hook);
        D_11.Update(scalar, I_0_, 1.0);
      }

      LINALG::Matrix<6, 1> D_12(sigma_bar_dev);
      {
        LINALG::Matrix<6, 6> temp2(false);
        D_12.Scale(2.0 / 3.0);
        temp2.MultiplyNN(1.0 / 3.0, I_d_, D_T_bar);
        D_12.MultiplyNN(1.0, temp2, m_, 1.0);
      }

      // K_uu = (B^T . D_11 . B) *Theta *detJ * w(gp) + K_geo
      {
        LINALG::Matrix<NUMDOF_SOH8, MAT::NUM_STRESS_3D> auxmat;
        auxmat.MultiplyTN(bopn, D_11);
        stiffmatrix->MultiplyNN(detJ_w_t, auxmat, bopn, 1.0);
      }

      // integrate 'geometric' stiffness matrix and add to first part of K_uu
      // K_G = N_xyz * sigma_hook * N_xyz^T * t_ * detJ * w(gp) * I

      LINALG::Matrix<NUMNOD_SOH8, NUMNOD_SOH8> G_bar(false);
      for (int a = 0; a < NUMNOD_SOH8; ++a)
      {
        for (int b = 0; b < NUMNOD_SOH8; ++b)
        {
          G_bar(a, b) = N_xyz(0, a) * sigma_hook(0, 0) * N_xyz(0, b) +
                        N_xyz(1, a) * sigma_hook(1, 0) * N_xyz(1, b) +
                        N_xyz(2, a) * sigma_hook(2, 0) * N_xyz(2, b) +
                        N_xyz(0, a) * sigma_hook(3, 0) * N_xyz(1, b) +
                        N_xyz(1, a) * sigma_hook(3, 0) * N_xyz(0, b) +
                        N_xyz(1, a) * sigma_hook(4, 0) * N_xyz(2, b) +
                        N_xyz(2, a) * sigma_hook(4, 0) * N_xyz(1, b) +
                        N_xyz(0, a) * sigma_hook(5, 0) * N_xyz(2, b) +
                        N_xyz(2, a) * sigma_hook(5, 0) * N_xyz(0, b);
        }
      }
      G_bar.Scale(detJ_w_t);

      // add onto stiffness matrix
      for (int inod = 0; inod < NUMNOD_SOH8; ++inod)
      {
        for (int jnod = 0; jnod < NUMNOD_SOH8; ++jnod)
        {
          (*stiffmatrix)(3 * inod + 0, 3 * jnod + 0) += G_bar(inod, jnod);
          (*stiffmatrix)(3 * inod + 1, 3 * jnod + 1) += G_bar(inod, jnod);
          (*stiffmatrix)(3 * inod + 2, 3 * jnod + 2) += G_bar(inod, jnod);
        }
      }

      // K_tu = (D_12^T . B) * N_t  * detJ * w(gp) , wobei N_t=1
      K_tu_.MultiplyTN(detJ_w, D_12, bopn, 1.0);

      // K_pu = ( m^T . B) * N_p * J * detJ * w(gp) , wobei N_p=1
      K_pu_.MultiplyTN(detJ_w_J, m_, bopn, 1.0);

      // K_tt = N_t * D_22 * N_t * Theta * detJ * w(gp), wobei N_t=1
      K_tt_ += D_22(0, 0) * detJ_w / t_(0, 0);

    }  // end of stiffness matrix ++++++++++++++++++++++++++++++++++++++++++++++

    if (massmatrix != NULL)  // evaluate mass matrix +++++++++++++++++++++++++
    {
      double density = Material()->Density(gp);
      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      for (int inod = 0; inod < NUMNOD_SOH8; ++inod)
      {
        const double ifactor = shapefcts[gp](inod) * factor;
        for (int jnod = 0; jnod < NUMNOD_SOH8; ++jnod)
        {
          const double massfactor = shapefcts[gp](jnod) * ifactor;  // intermediate factor
          (*massmatrix)(NUMDIM_SOH8 * inod + 0, NUMDIM_SOH8 * jnod + 0) += massfactor;
          (*massmatrix)(NUMDIM_SOH8 * inod + 1, NUMDIM_SOH8 * jnod + 1) += massfactor;
          (*massmatrix)(NUMDIM_SOH8 * inod + 2, NUMDIM_SOH8 * jnod + 2) += massfactor;
        }
      }
    }  // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++

    /* =========================================================================*/
  } /* ==================================================== end of Loop over GP */
  /* =========================================================================*/

  //=============================================================================
  // ************   do static condensation -> reduced stiffmatrix  **************
  //=============================================================================

  // K_uu_.PutScalar(0.0);
  K_uu_.Update(1.0, *stiffmatrix);
  // F_u_.PutScalar(0.0);
  F_u_.Update(1.0, *force);

  if (split_res)
    if (params.get<int>("MyPID") == Owner())
      params.get<double>("cond_rhs_norm") += R_p_(0, 0) * R_p_(0, 0) + R_t_(0, 0) * R_t_(0, 0);

  // K_t= K_uu + K_ut * K_pt^-1 * K_pu + K_up * K_tp^-1 * K_tu + K_up * K_tp^-1 * K_tt K_pt^-1 *
  // K_pu
  const double scalar = 1.0 / K_pt_ * K_tt_ * 1.0 / K_pt_;

  stiffmatrix->MultiplyTN(-1.0 / K_pt_, K_tu_, K_pu_, 1.0);
  stiffmatrix->MultiplyTN(-1.0 / K_pt_, K_pu_, K_tu_, 1.0);
  stiffmatrix->MultiplyTN(scalar, K_pu_, K_pu_, 1.0);

  force->MultiplyTN(scalar, K_pu_, R_p_, 1.0);
  force->MultiplyTN(-1.0 / K_pt_, K_pu_, R_t_, 1.0);
  force->MultiplyTN(-1.0 / K_pt_, K_tu_, R_p_, 1.0);

  return;
}  // DRT::ELEMENTS::So_sh8::ForceStiffMass



/*----------------------------------------------------------------------------*
 |  convert constitutive tensor (material -> current configuration)   lw 02/09|
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::So_Hex8P1J1::ConvertMat(
    const LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& cmat,
    const LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& F,
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& D_T_bar, const double t)
{
  // relationship between cmat (material configuration) and D_T_bar (current configuration):

  // D_T_bar[i][j][k][l] = 1/t * F[i][I] * F[j][J] * F[k][K] * F[l][L] * cmat[I][J][K][L]
  // c^ijkl              = 1/t   F^i_I     F^j_J     F^k_K     F^l_L     C^IJKL

#if 0
  std::vector<int> Index1(MAT::NUM_STRESS_3D);
  std::vector<int> Index2(MAT::NUM_STRESS_3D);

  Index1[0] = 0; Index2[0] = 0;
  Index1[1] = 1; Index2[1] = 1;
  Index1[2] = 2; Index2[2] = 2;
  Index1[3] = 0; Index2[3] = 1;
  Index1[4] = 1; Index2[4] = 2;
  Index1[5] = 2; Index2[5] = 0;


  // matrix containing products of deformation gradient elements
  LINALG::Matrix<6,6> FxF;

  for (int i=0; i<MAT::NUM_STRESS_3D; ++i)
  {
    FxF(0,i) = F(Index1[i],0)*F(Index2[i],0);
    FxF(1,i) = F(Index1[i],1)*F(Index2[i],1);
    FxF(2,i) = F(Index1[i],2)*F(Index2[i],2);
    FxF(3,i) = F(Index1[i],0)*F(Index2[i],1) + F(Index1[i],1)*F(Index2[i],0);
    FxF(4,i) = F(Index1[i],1)*F(Index2[i],2) + F(Index1[i],2)*F(Index2[i],1);
    FxF(5,i) = F(Index1[i],0)*F(Index2[i],2) + F(Index1[i],2)*F(Index2[i],0);
  }

  // vectors containing one column of FxF each
  LINALG::Matrix<6,1> FxF1;
  LINALG::Matrix<6,1> FxF2;

  // temporary vectors needed for consecutive multiplications
  LINALG::Matrix<6,1> temp1(true);
  LINALG::Matrix<1,1> temp2(true);

  for (int i=0; i<6; ++i)
  {
    for (int k=0; k<6; ++k)
      FxF1(k,0) = FxF(k,i);

    temp1.MultiplyNN(1., cmat, FxF1, 0.);

    for (int j=0;j<6; ++j)
    {
      if (j<i)
        D_T_bar(j,i) = D_T_bar(i,j);

      else
      {
        for (int k=0; k<6; ++k)
          FxF2(k,0) = FxF(k,j);

        temp2.MultiplyTN(1./t, FxF2, temp1, 0.);
        D_T_bar(j,i) = temp2(0,0);
      }
    }
  }
#else
  LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> FxF(false);
  PushPullOperator(FxF, F, false, 1.0);
  LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> CxFxF(false);
  CxFxF.MultiplyNT(cmat, FxF);
  D_T_bar.MultiplyNN(FxF, CxFxF);
  D_T_bar.Scale(1.0 / t);
#endif

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_Hex8P1J1::Stress(LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestress,
    const INPAR::STR::StressType iostress, const int gp, const double& detdefgrd,
    const LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& defgrd,
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& stress)
{
  switch (iostress)
  {
    case INPAR::STR::stress_2pk:  // 2nd Piola-Kirchhoff stress
    {
      if (elestress == NULL) dserror("stress data not available");

      // inverse deformation gradient
      LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> invdefgrd(defgrd);
      invdefgrd.Invert();

      // pull back operator
      LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> invdefgradinvdefgradT;
      PushPullOperator(invdefgradinvdefgradT, invdefgrd, false, detdefgrd);

      // (deviatoric) Cauchy stress vector
      LINALG::Matrix<MAT::NUM_STRESS_3D, 1> pk2;
      pk2.MultiplyNN(detdefgrd, invdefgradinvdefgradT, stress);

      // store stress
      for (int i = 0; i < MAT::NUM_STRESS_3D; ++i) (*elestress)(gp, i) = pk2(i);
    }
    break;
    case INPAR::STR::stress_cauchy:  // true/Cauchy stress
    {
      if (elestress == NULL) dserror("stress data not available");

      // store stress
      for (int i = 0; i < MAT::NUM_STRESS_3D; ++i) (*elestress)(gp, i) = stress(i);
    }
    break;
    case INPAR::STR::stress_none:
      break;
    default:
      dserror("requested stress option not available");
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_Hex8P1J1::Strain(
    LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestrain,  ///< store the strain herein
    const INPAR::STR::StrainType iostrain,  ///< strain type to store for post-proc
    const int gp,                           ///< Gauss point index
    const double& detdefgrd,                ///< determinant of (assumed) deformation gradient
    const LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& defgrd,  ///< deformation gradient
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& glstrain    ///< Green-Lagrange strain vector
)
{
  switch (iostrain)
  {
    case INPAR::STR::strain_gl:  // Green-Lagrange strain
    {
      if (elestrain == NULL) dserror("strain data not available");
      // store
      for (int i = 0; i < NUMDIM_SOH8; ++i) (*elestrain)(gp, i) = glstrain(i);
      for (int i = NUMDIM_SOH8; i < MAT::NUM_STRESS_3D; ++i)
        (*elestrain)(gp, i) = 0.5 * glstrain(i);
    }
    break;
    case INPAR::STR::strain_ea:  // Euler-Almansi strain
    {
      if (elestrain == NULL) dserror("strain data not available");

      // inverse deformation gradient
      LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> invdefgrd(defgrd);
      invdefgrd.Invert();

      // create push forward 6x6 matrix
      LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> invdefgradTinvdefgrad;
      PushPullOperator(invdefgradTinvdefgrad, invdefgrd, true, 1.0);

      // push forward
      LINALG::Matrix<MAT::NUM_STRESS_3D, 1> eastrain;
      eastrain.MultiplyNN(invdefgradTinvdefgrad, glstrain);
      // store
      for (int i = 0; i < NUMDIM_SOH8; ++i) (*elestrain)(gp, i) = eastrain(i);
      for (int i = NUMDIM_SOH8; i < MAT::NUM_STRESS_3D; ++i)
        (*elestrain)(gp, i) = 0.5 * eastrain(i);
    }
    break;
    case INPAR::STR::strain_none:
      break;
    default:
      dserror("requested strain option not available");
  }

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_Hex8P1J1::PushPullOperator(
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& g,  // G_IJ^KL or G^IJ_KL
    const LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>& f,          // [F^-1]=[F^B_b] or [F]=[F^b_B]
    const bool& transpose,                                      // co-variant if true
    const double& fac                                           //  a scaling factor
)
{
  // co-variant: strain-like
  // push forward:   e^flat  = fac  F^-T . E^flat . F^-1  // fac = 1
  //                 e_ab = fac  F_a^A  E_AB  F^B_b
  //                      = fac  F_a^A  F_b^B  E_AB
  //                      = fac  G_ab^AB  E_AB
  // pull back:      E^flat  = fac  F^T . e^flat . F   // fac = 1
  //                 E_AB = fac  F_A^a  e_ab  F^b_B
  //                      = fac  G_AB^ab  e_ab
  //
  // co-variant Voigt vector notation:
  //                 [E_11 E_22 E_33 2*E_12 2*E_23 2*E_31]
  //
  // contra-variant: stress-like
  // pull back:      S^sharp  = fac  F^-1 . s^sharp . F^-T    // fac = det(J)
  //                 S^AB = fac  F^A_a  s^ab  F_b^B
  //                      = fac  F^A_a  F^B_b  s^ab
  //                      = fac  G^AB_ab  s^ab
  //
  // push forward:   s^sharp = fac  F . S^sharp . F^T  // fac = 1/det(J)
  //                 s^ab = fac  F^a_A  S^AB  F_B^b
  //                      = fac  F^a_A  F^b_B  S^AB
  //                      = fac  G^ab_AB  S^AB
  // contra-variant Voigt vector notation:
  //                 [S^11 S^22 S^33 S^12 S^23 S^31]
  if (transpose)
  {
    // G_ab^AB or G_AB^ab
    //   rows _ab (or _AB) are co-variant/stress-like
    //   cols ^AB (or ^ab) are contra-variant/strain-like
    g(0, 0) = f(0, 0) * f(0, 0);
    g(0, 1) = f(1, 0) * f(1, 0);
    g(0, 2) = f(2, 0) * f(2, 0);
    g(0, 3) = f(0, 0) * f(1, 0);
    g(0, 4) = f(1, 0) * f(2, 0);
    g(0, 5) = f(0, 0) * f(2, 0);

    g(1, 0) = f(0, 1) * f(0, 1);
    g(1, 1) = f(1, 1) * f(1, 1);
    g(1, 2) = f(2, 1) * f(2, 1);
    g(1, 3) = f(0, 1) * f(1, 1);
    g(1, 4) = f(1, 1) * f(2, 1);
    g(1, 5) = f(0, 1) * f(2, 1);

    g(2, 0) = f(0, 2) * f(0, 2);
    g(2, 1) = f(1, 2) * f(1, 2);
    g(2, 2) = f(2, 2) * f(2, 2);
    g(2, 3) = f(0, 2) * f(1, 2);
    g(2, 4) = f(1, 2) * f(2, 2);
    g(2, 5) = f(0, 2) * f(2, 2);

    g(3, 0) = 2.0 * f(0, 0) * f(0, 1);
    g(3, 1) = 2.0 * f(1, 0) * f(1, 1);
    g(3, 2) = 2.0 * f(2, 0) * f(2, 1);
    g(3, 3) = f(0, 0) * f(1, 1) + f(0, 1) * f(1, 0);
    g(3, 4) = f(1, 0) * f(2, 1) + f(1, 1) * f(2, 0);
    g(3, 5) = f(0, 0) * f(2, 1) + f(0, 1) * f(2, 0);

    g(4, 0) = 2.0 * f(0, 1) * f(0, 2);
    g(4, 1) = 2.0 * f(1, 1) * f(1, 2);
    g(4, 2) = 2.0 * f(2, 1) * f(2, 2);
    g(4, 3) = f(0, 1) * f(1, 2) + f(0, 2) * f(1, 1);
    g(4, 4) = f(1, 1) * f(2, 2) + f(1, 2) * f(2, 1);
    g(4, 5) = f(0, 1) * f(2, 2) + f(0, 2) * f(2, 1);

    g(5, 0) = 2.0 * f(0, 0) * f(0, 2);
    g(5, 1) = 2.0 * f(1, 0) * f(1, 2);
    g(5, 2) = 2.0 * f(2, 0) * f(2, 2);
    g(5, 3) = f(0, 0) * f(1, 2) + f(0, 2) * f(1, 0);
    g(5, 4) = f(1, 0) * f(2, 2) + f(1, 2) * f(2, 0);
    g(5, 5) = f(0, 0) * f(2, 2) + f(0, 2) * f(2, 0);
  }
  else
  {
    // G^ab_AB or G^AB_ab
    //   rows ^AB (or ^ab) are contra-variant/strain-like
    //   cols _ab (or _AB) are co-variant/stress-like
    g(0, 0) = f(0, 0) * f(0, 0);
    g(0, 1) = f(0, 1) * f(0, 1);
    g(0, 2) = f(0, 2) * f(0, 2);
    g(0, 3) = 2.0 * f(0, 0) * f(0, 1);
    g(0, 4) = 2.0 * f(0, 1) * f(0, 2);
    g(0, 5) = 2.0 * f(0, 0) * f(0, 2);

    g(1, 0) = f(1, 0) * f(1, 0);
    g(1, 1) = f(1, 1) * f(1, 1);
    g(1, 2) = f(1, 2) * f(1, 2);
    g(1, 3) = 2.0 * f(1, 0) * f(1, 1);
    g(1, 4) = 2.0 * f(1, 1) * f(1, 2);
    g(1, 5) = 2.0 * f(1, 0) * f(1, 2);

    g(2, 0) = f(2, 0) * f(2, 0);
    g(2, 1) = f(2, 1) * f(2, 1);
    g(2, 2) = f(2, 2) * f(2, 2);
    g(2, 3) = 2.0 * f(2, 0) * f(2, 1);
    g(2, 4) = 2.0 * f(2, 1) * f(2, 2);
    g(2, 5) = 2.0 * f(2, 0) * f(2, 2);

    g(3, 0) = f(0, 0) * f(1, 0);
    g(3, 1) = f(0, 1) * f(1, 1);
    g(3, 2) = f(0, 2) * f(1, 2);
    g(3, 3) = f(0, 0) * f(1, 1) + f(0, 1) * f(1, 0);
    g(3, 4) = f(0, 1) * f(1, 2) + f(0, 2) * f(1, 1);
    g(3, 5) = f(0, 0) * f(1, 2) + f(0, 2) * f(1, 0);

    g(4, 0) = f(1, 0) * f(2, 0);
    g(4, 1) = f(1, 1) * f(2, 1);
    g(4, 2) = f(1, 2) * f(2, 2);
    g(4, 3) = f(1, 0) * f(2, 1) + f(1, 1) * f(2, 0);
    g(4, 4) = f(1, 1) * f(2, 2) + f(1, 2) * f(2, 1);
    g(4, 5) = f(1, 0) * f(2, 2) + f(1, 2) * f(2, 0);

    g(5, 0) = f(0, 0) * f(2, 0);
    g(5, 1) = f(0, 1) * f(2, 1);
    g(5, 2) = f(0, 2) * f(2, 2);
    g(5, 3) = f(0, 0) * f(2, 1) + f(0, 1) * f(2, 0);
    g(5, 4) = f(0, 1) * f(2, 2) + f(0, 2) * f(2, 1);
    g(5, 5) = f(0, 0) * f(2, 2) + f(0, 2) * f(2, 0);
  }

  // apply scaling
  if (fac != 1.0) g.Scale(fac);

  // done
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
#if 0
void DRT::ELEMENTS::So_Hex8P1J1::test_stiffmat(
  const std::vector<int>& lm,  ///< location matrix
  const std::vector<double>& disp,  ///< current displacements
  const std::vector<double>& residual,   ///< current residual displ
  Teuchos::ParameterList& params  ///< algorithmic parameters e.g. time
  )
{
  //double epsilon = 0.0;
  double epsilon = 0.0000001;
  LINALG::Matrix<NUMDOF_SOH8,1> F(true);
  LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8> K(true);
  std::vector<double> disturbed_disp(disp);
  std::vector<double> disturbed_disi(residual);

  LINALG::Matrix<24,24> App_K(true);

  for (int disturb=0; disturb<24; ++disturb)
  {
    disturbed_disp[disturb] += epsilon;
    //disturbed_disi[disturb] += epsilon;

    /* ============================================================================*
    ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
    ** ============================================================================*/
    const static std::vector<LINALG::Matrix<NUMNOD_SOH8,1> > shapefcts = soh8_shapefcts();
    const static std::vector<LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> > derivs = soh8_derivs();
    const static std::vector<double> gpweights = soh8_weights();
    /* ============================================================================*/

    // update element geometry
    LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xrefe;  // material coord. of element
    LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xcurr;  // current  coord. of element
    DRT::Node** nodes = Nodes();
    for (int i=0; i<NUMNOD_SOH8; ++i)
    {
      const double* x = nodes[i]->X();
      xrefe(i,0) = x[0];
      xrefe(i,1) = x[1];
      xrefe(i,2) = x[2];

      xcurr(i,0) = xrefe(i,0) + disturbed_disp[i*NODDOF_SOH8+0];
      xcurr(i,1) = xrefe(i,1) + disturbed_disp[i*NODDOF_SOH8+1];
      xcurr(i,2) = xrefe(i,2) + disturbed_disp[i*NODDOF_SOH8+2];
    }

    //************************************************************************************
    size_t length = disturbed_disi.size();
    if (length != 24) dserror("number of element displacement dofs not equal 24");

    LINALG::Matrix<24,1> disi(false);
    for (int i = 0; i < NUMDOF_SOH8; ++i)
    {
      disi(i,0) = disturbed_disi[i];
    }

    // dt= K_tp_^-1*(-R_p_-K_pu_*du)
    LINALG::Matrix<1,1> dt(true);
    dt.MultiplyNN(-1.0/K_pt_, K_pu_, disi);
    dt.Update(-1.0/K_pt_, R_p_, 1.0);

    // dp= K_tp_^-1 * (-R_t_ - K_tu_*du - K_tt*dt)
    LINALG::Matrix<1,1> dp(true);
    dp.MultiplyNN(1.0, K_tu_, disi);
    dp.Update(K_tt_, dt, 1.0);
    dp.Update(1.0, R_t_, 1.0);
    dp.Scale(-1.0/K_pt_);

    // t_ und p_ neu berechnen!
    LINALG::Matrix<1,1> t(t_temp_);
    t += dt;
    LINALG::Matrix<1,1> p(p_temp_);
    p += dp;

    //******************************************************************************************

    /* =========================================================================*/
    /* ================================================= Loop over Gauss Points */
    /* =========================================================================*/

    LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> N_XYZ;
    // build deformation gradient wrt to material configuration
    // in case of prestressing, build defgrd wrt to last stored configuration
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> defgrd(false);
    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> mod_defgrd(true);
    for (unsigned gp=0; gp<NUMGPT_SOH8; ++gp)
    {
      /* get the inverse of the Jacobian matrix which looks like:
      **            [ x_,r  y_,r  z_,r ]^-1
      **     J^-1 = [ x_,s  y_,s  z_,s ]
      **            [ x_,t  y_,t  z_,t ]
      */
      // compute derivatives N_XYZ at gp w.r.t. material coordinates
      // by N_XYZ = J^-1 * N_rst
      N_XYZ.Multiply(invJ_[gp],derivs[gp]);
      double detJ = detJ_[gp];

      // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
      defgrd.MultiplyTT(xcurr,N_XYZ);

      // modified deformation gradient modF = (t_/J)^1/3 * F
      const double J = defgrd.Determinant();
      double scalar = std::pow((t(0, 0)/J), 1.0/3.0);
      mod_defgrd.Update(scalar, defgrd);

      // Modified Right Cauchy-Green tensor = modF^T * modF
      LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> mod_cauchygreen;
      mod_cauchygreen.MultiplyTN(mod_defgrd, mod_defgrd);

      // Modified Green-Lagrange strains matrix mod_E = 0.5 * (modCauchygreen - Identity)
      // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
      Epetra_SerialDenseVector mod_glstrain_epetra(MAT::NUM_STRESS_3D);
      LINALG::Matrix<MAT::NUM_STRESS_3D,1> mod_glstrain(mod_glstrain_epetra.A(),true);
      mod_glstrain(0) = 0.5 * (mod_cauchygreen(0,0) - 1.0);
      mod_glstrain(1) = 0.5 * (mod_cauchygreen(1,1) - 1.0);
      mod_glstrain(2) = 0.5 * (mod_cauchygreen(2,2) - 1.0);
      mod_glstrain(3) = mod_cauchygreen(0,1);
      mod_glstrain(4) = mod_cauchygreen(1,2);
      mod_glstrain(5) = mod_cauchygreen(2,0);

      // Bestimmen der N_xyz
      LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_xyz(true);
      LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> Inv_defgrd(true);
      Inv_defgrd.Invert(defgrd);

      for (int i=0; i<NUMNOD_SOH8; ++i)
      {
        for (int j=0; j<NUMDIM_SOH8; ++j)
        {
          for (int k=0; k<NUMDIM_SOH8; ++k)
          {
            N_xyz(j,i) += N_XYZ(k,i)*Inv_defgrd(k,j);
          }
        }
      }

      //=================================================================================================
      // ***************************    B-matrix from Zienkiewicz   *************************************
      //=================================================================================================
      LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8> bopn;
      for (int i=0; i<NUMNOD_SOH8; ++i)
      {
        bopn(0, NODDOF_SOH8*i+0) = N_xyz(0,i);
        bopn(0, NODDOF_SOH8*i+1) = 0;
        bopn(0, NODDOF_SOH8*i+2) = 0;
        //------------------------------
        bopn(1, NODDOF_SOH8*i+0) = 0;
        bopn(1, NODDOF_SOH8*i+1) = N_xyz(1,i);
        bopn(1, NODDOF_SOH8*i+2) = 0;
        //-----------------------------
        bopn(2, NODDOF_SOH8*i+0) = 0;
        bopn(2, NODDOF_SOH8*i+1) = 0;
        bopn(2, NODDOF_SOH8*i+2) = N_xyz(2,i);
        //-----------------------------
        bopn(3, NODDOF_SOH8*i+0) = N_xyz(1,i);
        bopn(3, NODDOF_SOH8*i+1) = N_xyz(0,i);
        bopn(3, NODDOF_SOH8*i+2) = 0;
        //----------------------------
        bopn(4, NODDOF_SOH8*i+0) = 0;
        bopn(4, NODDOF_SOH8*i+1) = N_xyz(2,i);
        bopn(4, NODDOF_SOH8*i+2) = N_xyz(1,i);
        //----------------------------
        bopn(5, NODDOF_SOH8*i+0) = N_xyz(2,i);
        bopn(5, NODDOF_SOH8*i+1) = 0;
        bopn(5, NODDOF_SOH8*i+2) = N_xyz(0,i);
      }

      // call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D> cmat(true);
      LINALG::Matrix<MAT::NUM_STRESS_3D,1> stress(true);
      params.set<int>("gp",gp);
      SolidMaterial()->Evaluate(&defgrd,&mod_glstrain,params,&stress,&cmat,Id());

      // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

      LINALG::Matrix<MAT::NUM_STRESS_3D,1> sigma_bar(true);
      LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> stress_matrix(false);

      stress_matrix(0,0) = stress(0,0);
      stress_matrix(1,1) = stress(1,0);
      stress_matrix(2,2) = stress(2,0);
      stress_matrix(0,1) = stress_matrix(1,0) = stress(3,0);
      stress_matrix(1,2) = stress_matrix(2,1) = stress(4,0);
      stress_matrix(0,2) = stress_matrix(2,0) = stress(5,0);

      for (int i=0; i<3; ++i)
      {
        for (int j=0; j<3; ++j)
        {
          sigma_bar(0,0) += mod_defgrd(0,i) * stress_matrix(i,j) * mod_defgrd(0,j);
          sigma_bar(1,0) += mod_defgrd(1,i) * stress_matrix(i,j) * mod_defgrd(1,j);
          sigma_bar(2,0) += mod_defgrd(2,i) * stress_matrix(i,j) * mod_defgrd(2,j);
          sigma_bar(3,0) += mod_defgrd(0,i) * stress_matrix(i,j) * mod_defgrd(1,j);
          sigma_bar(4,0) += mod_defgrd(1,i) * stress_matrix(i,j) * mod_defgrd(2,j);
          sigma_bar(5,0) += mod_defgrd(0,i) * stress_matrix(i,j) * mod_defgrd(2,j);
        }
      }
      sigma_bar.Scale(1.0/t(0,0));

      double p_hook = p(0,0) * J / t(0,0);
      double p_bar = 1.0/3.0 * (sigma_bar(0,0) + sigma_bar(1,0) + sigma_bar(2,0));

      LINALG::Matrix<MAT::NUM_STRESS_3D,1> sigma_hook(sigma_bar);

      for (int i=0; i<3; ++i)
      {
        sigma_hook(i,0) += p_hook - p_bar;
      }

      LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> sigma_hook_matrix(false);
      sigma_hook_matrix(0,0) = sigma_hook(0,0);
      sigma_hook_matrix(1,1) = sigma_hook(1,0);
      sigma_hook_matrix(2,2) = sigma_hook(2,0);
      sigma_hook_matrix(0,1) = sigma_hook_matrix(1,0) = sigma_hook(3,0);
      sigma_hook_matrix(1,2) = sigma_hook_matrix(2,1) = sigma_hook(4,0);
      sigma_hook_matrix(0,2) = sigma_hook_matrix(2,0) = sigma_hook(5,0);

      LINALG::Matrix<MAT::NUM_STRESS_3D,1> sigma_bar_dev(sigma_bar);
      for (int i=0; i<3; ++i)
      {
        sigma_bar_dev(i,0) -= p_bar;
      }

      double detJ_w = detJ * gpweights[gp];
      double detJ_w_t = detJ_w * t(0,0);

      F.MultiplyTN(detJ_w_t, bopn, sigma_hook, 1.0);


     /* =========================================================================*/
    }/* ==================================================== end of Loop over GP */
     /* =========================================================================*/

    for (int j=0; j<24; ++j)
      App_K(j, disturb) = (F(j, 0)- F_u_(j, 0));

    disturbed_disp[disturb] = disp[disturb];
    disturbed_disi[disturb] = residual[disturb];
    F.PutScalar(0.0);
    K.PutScalar(0.0);
  }

  App_K.Scale(1.0/epsilon);
  App_K.Update(-1.0, K_uu_, 1.0);
  std::cout << "A_pp:" << App_K << std::endl;


//     for (int i=0; i<24; ++i)
//     {
//       for (int j=0; j<24; ++j)
//       {
//         if (App_K(i, j) > 0.000001)
//           std::cout << "Abweichung in K_" << i << j << ": " << App_K(i, j) << std::endl;
//       }
//     }
}
#endif

/*----------------------------------------------------------------------*
 |  init the element (public)                                   lw 12/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_Hex8P1J1Type::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So_Hex8P1J1* actele =
        dynamic_cast<DRT::ELEMENTS::So_Hex8P1J1*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_Hex8P1J1* failed");
    actele->InitJacobianMapping();
    actele->InitKpt();
  }

  return 0;
}


void DRT::ELEMENTS::So_Hex8P1J1::soh8P1J1_recover(const std::vector<double>& residual)
{
  LINALG::Matrix<24, 1> disi(false);
  for (int i = 0; i < NUMDOF_SOH8; ++i) disi(i) = residual[i];

  // dt= K_tp_^-1*(-R_p_-K_pu_*du)
  dt_.MultiplyNN(-1.0 / K_pt_, K_pu_, disi);
  dt_.Update(-1.0 / K_pt_, R_p_, 1.0);

  // dp= K_tp_^-1 * (-R_t_ - K_tu_*du - K_tt*dt)
  dp_.MultiplyNN(1.0, K_tu_, disi);
  dp_.Update(K_tt_, dt_, 1.0);
  dp_.Update(1.0, R_t_, 1.0);
  dp_.Scale(-1.0 / K_pt_);

  t_ += dt_;
  p_ += dp_;
}

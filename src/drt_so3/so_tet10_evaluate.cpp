/*!----------------------------------------------------------------------*
\file so_tet10_evaluate.cpp
\brief quadratic nonlinear tetrahedron
\level 1

<pre>
\maintainer Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>

*----------------------------------------------------------------------*/
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "so_tet10.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../drt_patspec/patspec.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_contact/contact_analytical.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/so3_material.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_fem_general/drt_utils_integration.H"
#include "so_utils.H"
#include "../drt_fiber/drt_fiber_node.H"

// inverse design object
#include "inversedesign.H"
#include "prestress.H"

#include "../drt_structure_new/str_elements_paramsinterface.H"

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet10::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  LINALG::Matrix<NUMDOF_SOTET10, NUMDOF_SOTET10> elemat1(elemat1_epetra.A(), true);
  LINALG::Matrix<NUMDOF_SOTET10, NUMDOF_SOTET10> elemat2(elemat2_epetra.A(), true);
  LINALG::Matrix<NUMDOF_SOTET10, 1> elevec1(elevec1_epetra.A(), true);
  LINALG::Matrix<NUMDOF_SOTET10, 1> elevec2(elevec2_epetra.A(), true);
  LINALG::Matrix<NUMDOF_SOTET10, 1> elevec3(elevec3_epetra.A(), true);

  // start with "none"
  DRT::ELEMENTS::So_tet10::ActionType act = So_tet10::none;

  // get the required action
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    dserror("No action supplied");
  else if (action == "calc_struct_linstiff")
    act = So_tet10::calc_struct_linstiff;
  else if (action == "calc_struct_nlnstiff")
    act = So_tet10::calc_struct_nlnstiff;
  else if (action == "calc_struct_internalforce")
    act = So_tet10::calc_struct_internalforce;
  else if (action == "calc_struct_linstiffmass")
    act = So_tet10::calc_struct_linstiffmass;
  else if (action == "calc_struct_nlnstiffmass")
    act = So_tet10::calc_struct_nlnstiffmass;
  else if (action == "calc_struct_nlnstifflmass")
    act = So_tet10::calc_struct_nlnstifflmass;
  else if (action == "calc_struct_stress")
    act = So_tet10::calc_struct_stress;
  else if (action == "calc_struct_eleload")
    act = So_tet10::calc_struct_eleload;
  else if (action == "calc_struct_fsiload")
    act = So_tet10::calc_struct_fsiload;
  else if (action == "calc_struct_update_istep")
    act = So_tet10::calc_struct_update_istep;
  else if (action == "calc_struct_reset_istep")
    act = So_tet10::calc_struct_reset_istep;
  else if (action == "calc_struct_reset_all")
    act = So_tet10::calc_struct_reset_all;
  else if (action == "calc_struct_energy")
    act = So_tet10::calc_struct_energy;
  else if (action == "calc_struct_errornorms")
    act = So_tet10::calc_struct_errornorms;
  else if (action == "calc_struct_prestress_update")
    act = So_tet10::prestress_update;
  else if (action == "calc_global_gpstresses_map")
    act = So_tet10::calc_global_gpstresses_map;
  else if (action == "postprocess_stress")
    act = So_tet10::postprocess_stress;
  else if (action == "calc_struct_recover")
    return 0;
  else if (action == "calc_struct_predict")
    return 0;
  else
    dserror("Unknown type of action for So_tet10");

  // check for patient specific data
  PATSPEC::GetILTDistance(Id(), params, discretization);
  PATSPEC::GetLocalRadius(Id(), params, discretization);
  PATSPEC::GetInnerRadius(Id(), params, discretization);

  // what should the element do
  switch (act)
  {
    // linear stiffness
    case calc_struct_linstiff:
    {
      // need current displacement and residual forces
      std::vector<double> mydisp(lm.size());
      for (unsigned i = 0; i < mydisp.size(); ++i) mydisp[i] = 0.0;
      std::vector<double> myres(lm.size());
      for (unsigned i = 0; i < myres.size(); ++i) myres[i] = 0.0;
      std::vector<double> mydispmat(lm.size(), 0.0);

      so_tet10_nlnstiffmass(lm, mydisp, NULL, NULL, myres, mydispmat, &elemat1, NULL, &elevec1,
          NULL, &elevec3, NULL, NULL, params, INPAR::STR::stress_none, INPAR::STR::strain_none);
    }
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
      LINALG::Matrix<NUMDOF_SOTET10, NUMDOF_SOTET10>* matptr = NULL;
      if (elemat1.IsInitialized()) matptr = &elemat1;

      std::vector<double> mydispmat(lm.size(), 0.0);

      so_tet10_nlnstiffmass(lm, mydisp, NULL, NULL, myres, mydispmat, matptr, NULL, &elevec1, NULL,
          &elevec3, NULL, NULL, params, INPAR::STR::stress_none, INPAR::STR::strain_none);
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
      std::vector<double> mydispmat(lm.size(), 0.0);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      LINALG::Matrix<NUMDOF_SOTET10, NUMDOF_SOTET10> myemat(true);

      so_tet10_nlnstiffmass(lm, mydisp, NULL, NULL, myres, mydispmat, &myemat, NULL, &elevec1, NULL,
          NULL, NULL, NULL, params, INPAR::STR::stress_none, INPAR::STR::strain_none);
    }
    break;

    // linear stiffness and consistent mass matrix
    case calc_struct_linstiffmass:
      dserror("Case 'calc_struct_linstiffmass' not yet implemented");
      break;

    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      // need current velocities and accelerations (for non constant mass matrix)
      Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState("velocity");
      Teuchos::RCP<const Epetra_Vector> acc = discretization.GetState("acceleration");
      if (disp == Teuchos::null || res == Teuchos::null)
        dserror("Cannot get state vectors 'displacement' and/or residual");
      if (vel == Teuchos::null) dserror("Cannot get state vectors 'velocity'");
      if (acc == Teuchos::null) dserror("Cannot get state vectors 'acceleration'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myvel(lm.size());
      DRT::UTILS::ExtractMyValues(*vel, myvel, lm);
      std::vector<double> myacc(lm.size());
      DRT::UTILS::ExtractMyValues(*acc, myacc, lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);

      std::vector<double> mydispmat(lm.size(), 0.0);

      so_tet10_nlnstiffmass(lm, mydisp, &myvel, &myacc, myres, mydispmat, &elemat1, &elemat2,
          &elevec1, &elevec2, &elevec3, NULL, NULL, params, INPAR::STR::stress_none,
          INPAR::STR::strain_none);

      if (act == calc_struct_nlnstifflmass) so_tet10_lumpmass(&elemat2);
    }
    break;

    // evaluate stresses and strains at gauss points
    case calc_struct_stress:
    {
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
        LINALG::Matrix<NUMGPT_SOTET10, MAT::NUM_STRESS_3D> stress;
        LINALG::Matrix<NUMGPT_SOTET10, MAT::NUM_STRESS_3D> strain;
        INPAR::STR::StressType iostress =
            DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
        INPAR::STR::StrainType iostrain =
            DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);

        std::vector<double> mydispmat(lm.size(), 0.0);

        so_tet10_nlnstiffmass(lm, mydisp, NULL, NULL, myres, mydispmat, NULL, NULL, NULL, NULL,
            NULL, &stress, &strain, params, iostress, iostrain);
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

        // output of rotation matrix R with F = U*R
        if (DRT::INPUT::IntegralValue<bool>(DRT::Problem::Instance()->IOParams(), "OUTPUT_ROT") ==
            true)
        {
          LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> R;
          DRT::ELEMENTS::UTILS::CalcR<DRT::Element::tet10>(this, mydisp, R);

          Teuchos::RCP<std::vector<char>> rotdata =
              params.get<Teuchos::RCP<std::vector<char>>>("rotation", Teuchos::null);

          DRT::PackBuffer data;
          AddtoPack(data, R);
          data.StartPacking();
          AddtoPack(data, R);
          std::copy(data().begin(), data().end(), std::back_inserter(*rotdata));
        }
      }
    }
    break;

    // postprocess stresses/strains at gauss points

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
      int gid = Id();

      Teuchos::RCP<Epetra_MultiVector> poststress =
          params.get<Teuchos::RCP<Epetra_MultiVector>>("poststress", Teuchos::null);
      if (poststress == Teuchos::null) dserror("No element stress/strain vector available");

      if (stresstype == "ndxyz")
      {
        LINALG::Matrix<NUMGPT_SOTET10, MAT::NUM_STRESS_3D> gpstress(
            ((*gpstressmap)[gid])->A(), true);

        // extrapolate stresses/strains at Gauss points to nodes
        so_tet10_expol(gpstress, *poststress);
      }
      else if (stresstype == "cxyz")
      {
        LINALG::Matrix<NUMGPT_SOTET10, MAT::NUM_STRESS_3D> gpstress(
            ((*gpstressmap)[gid])->A(), true);

        const Epetra_BlockMap elemap = poststress->Map();
        int lid = elemap.LID(Id());
        if (lid != -1)
        {
          for (int i = 0; i < MAT::NUM_STRESS_3D; ++i)
          {
            double& s = (*((*poststress)(i)))[lid];  // resolve pointer for faster access
            s = 0.;
            for (int j = 0; j < NUMGPT_SOTET10; ++j)
            {
              s += gpstress(j, i);
            }
            s *= 1.0 / NUMGPT_SOTET10;
          }
        }
      }
      else if (stresstype == "rotation")
      {
        LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> elecenterrot(
            ((*gpstressmap)[gid])->A(), true);

        const Epetra_BlockMap elemap = poststress->Map();
        int lid = elemap.LID(Id());
        if (lid != -1)
          for (int i = 0; i < NUMDIM_SOTET10; ++i)
            for (int j = 0; j < NUMDIM_SOTET10; ++j)
              (*((*poststress)(i * NUMDIM_SOTET10 + j)))[lid] = elecenterrot(i, j);
      }
      else
      {
        dserror("unknown type of stress/strain output on element level");
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
      // Update of history for materials
      SolidMaterial()->Update();
    }
    break;

    case prestress_update:
    {
      time_ = params.get<double>("total time");
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get displacement state");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      // build incremental def gradient for every gauss point
      LINALG::SerialDenseMatrix gpdefgrd(NUMGPT_SOTET10, 9);
      DefGradient(mydisp, gpdefgrd, *prestress_);

      // update deformation gradient and put back to storage
      LINALG::Matrix<3, 3> deltaF;
      LINALG::Matrix<3, 3> Fhist;
      LINALG::Matrix<3, 3> Fnew;
      for (int gp = 0; gp < NUMGPT_SOTET10; ++gp)
      {
        prestress_->StoragetoMatrix(gp, deltaF, gpdefgrd);
        prestress_->StoragetoMatrix(gp, Fhist, prestress_->FHistory());
        Fnew.Multiply(deltaF, Fhist);
        prestress_->MatrixtoStorage(gp, Fnew, prestress_->FHistory());
      }

      // push-forward invJ for every gaussian point
      UpdateJacobianMapping(mydisp, *prestress_);

      // Update constraintmixture material
      if (Material()->MaterialType() == INPAR::MAT::m_constraintmixture)
      {
        SolidMaterial()->Update();
      }
    }
    break;

    case calc_struct_reset_istep:
    {
      // Reset of history (if needed)
      SolidMaterial()->ResetStep();
    }
    break;
    //==================================================================================
    case calc_struct_reset_all:
    {
      // Reset of history for materials
      SolidMaterial()->ResetAll(NUMGPT_SOTET10);

      // Reset prestress
      if (pstype_ == INPAR::STR::prestress_mulf)
      {
        time_ = 0.0;
        LINALG::Matrix<3, 3> Id(true);
        Id(0, 0) = Id(1, 1) = Id(2, 2) = 1.0;
        for (int gp = 0; gp < NUMGPT_SOTET10; ++gp)
        {
          prestress_->MatrixtoStorage(gp, Id, prestress_->FHistory());
          prestress_->MatrixtoStorage(gp, invJ_[gp], prestress_->JHistory());
        }
      }
      if (pstype_ == INPAR::STR::prestress_id)
        dserror("Reset of Inverse Design not yet implemented");
    }
    break;
    //==================================================================================
    case calc_struct_energy:
    {
      // check length of elevec1
      if (elevec1_epetra.Length() < 1) dserror("The given result vector is too short.");

      // initialization of internal energy
      double intenergy = 0.0;

      const static std::vector<LINALG::Matrix<NUMNOD_SOTET10, 1>> shapefcts_4gp =
          so_tet10_4gp_shapefcts();
      const static std::vector<LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs_4gp =
          so_tet10_4gp_derivs();
      const static std::vector<double> gpweights_4gp = so_tet10_4gp_weights();

      // get displacements of this processor
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state displacement vector");

      // get displacements of this element
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      // update element geometry
      LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xrefe;  // material coord. of element
      LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xcurr;  // current  coord. of element
      LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xdisp;
      DRT::Node** nodes = Nodes();
      for (int i = 0; i < NUMNOD_SOTET10; ++i)
      {
        const double* x = nodes[i]->X();
        xrefe(i, 0) = x[0];
        xrefe(i, 1) = x[1];
        xrefe(i, 2) = x[2];

        xcurr(i, 0) = xrefe(i, 0) + mydisp[i * NODDOF_SOTET10 + 0];
        xcurr(i, 1) = xrefe(i, 1) + mydisp[i * NODDOF_SOTET10 + 1];
        xcurr(i, 2) = xrefe(i, 2) + mydisp[i * NODDOF_SOTET10 + 2];

        if (pstype_ == INPAR::STR::prestress_mulf)
        {
          xdisp(i, 0) = mydisp[i * NODDOF_SOTET10 + 0];
          xdisp(i, 1) = mydisp[i * NODDOF_SOTET10 + 1];
          xdisp(i, 2) = mydisp[i * NODDOF_SOTET10 + 2];
        }
      }
      /* =========================================================================*/
      /* ================================================= Loop over Gauss Points */
      /* =========================================================================*/
      LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10> N_XYZ;
      // build deformation gradient wrt to material configuration
      // in case of prestressing, build defgrd wrt to last stored configuration
      LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> defgrd(false);
      for (int gp = 0; gp < NUMGPT_SOTET10; ++gp)
      {
        /* get the inverse of the Jacobian matrix which looks like:
        **            [ x_,r  y_,r  z_,r ]^-1
        **     J^-1 = [ x_,s  y_,s  z_,s ]
        **            [ x_,t  y_,t  z_,t ]
        */
        // compute derivatives N_XYZ at gp w.r.t. material coordinates
        // by N_XYZ = J^-1 * N_rst
        N_XYZ.Multiply(invJ_[gp], derivs_4gp[gp]);

        // Gauss weights and Jacobian determinant
        double fac = detJ_[gp] * gpweights_4gp[gp];

        if (pstype_ == INPAR::STR::prestress_mulf)
        {
          // get Jacobian mapping wrt to the stored configuration
          LINALG::Matrix<3, 3> invJdef;
          prestress_->StoragetoMatrix(gp, invJdef, prestress_->JHistory());
          // get derivatives wrt to last spatial configuration
          LINALG::Matrix<3, 10> N_xyz;
          N_xyz.Multiply(invJdef, derivs_4gp[gp]);

          // build multiplicative incremental defgrd
          defgrd.MultiplyTT(xdisp, N_xyz);
          defgrd(0, 0) += 1.0;
          defgrd(1, 1) += 1.0;
          defgrd(2, 2) += 1.0;

          // get stored old incremental F
          LINALG::Matrix<3, 3> Fhist;
          prestress_->StoragetoMatrix(gp, Fhist, prestress_->FHistory());

          // build total defgrd = delta F * F_old
          LINALG::Matrix<3, 3> Fnew;
          Fnew.Multiply(defgrd, Fhist);
          defgrd = Fnew;
        }
        else
        {
          defgrd.MultiplyTT(xcurr, N_XYZ);
        }

        // Right Cauchy-Green tensor = F^T * F
        LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> cauchygreen;
        cauchygreen.MultiplyTN(defgrd, defgrd);

        // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
        // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
        Epetra_SerialDenseVector glstrain_epetra(MAT::NUM_STRESS_3D);
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(glstrain_epetra.A(), true);
        glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
        glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
        glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
        glstrain(3) = cauchygreen(0, 1);
        glstrain(4) = cauchygreen(1, 2);
        glstrain(5) = cauchygreen(2, 0);

        // call material for evaluation of strain energy function
        double psi = 0.0;
        SolidMaterial()->StrainEnergy(glstrain, psi, Id());

        // sum up GP contribution to internal energy
        intenergy += fac * psi;
      }

      // return result
      elevec1_epetra(0) = intenergy;
    }
    break;


    case calc_struct_errornorms:
    {
      // IMPORTANT NOTES (popp 10/2010):
      // - error norms are based on a small deformation assumption (linear elasticity)
      // - extension to finite deformations would be possible without difficulties,
      //   however analytical solutions are extremely rare in the nonlinear realm
      // - only implemented for SVK material (relevant for energy norm only, L2 and
      //   H1 norms are of course valid for arbitrary materials)
      // - analytical solutions are currently stored in a repository in the CONTACT
      //   namespace, however they could (should?) be moved to a more general location

      // check length of elevec1
      if (elevec1_epetra.Length() < 3) dserror("The given result vector is too short.");

      // check material law
      Teuchos::RCP<MAT::Material> mat = Material();

      //******************************************************************
      // only for St.Venant Kirchhoff material
      //******************************************************************
      if (mat->MaterialType() == INPAR::MAT::m_stvenant)
      {
        // declaration of variables
        double l2norm = 0.0;
        double h1norm = 0.0;
        double energynorm = 0.0;

        // shape functions, derivatives and integration weights
        const static std::vector<LINALG::Matrix<NUMNOD_SOTET10, 1>> vals =
            so_tet10_11gp_shapefcts();
        const static std::vector<LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs =
            so_tet10_11gp_derivs();
        const static std::vector<double> weights = so_tet10_11gp_weights();

        // get displacements and extract values of this element
        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        if (disp == Teuchos::null) dserror("Cannot get state displacement vector");
        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

        // nodal displacement vector
        LINALG::Matrix<NUMDOF_SOTET10, 1> nodaldisp;
        for (int i = 0; i < NUMDOF_SOTET10; ++i) nodaldisp(i, 0) = mydisp[i];

        // reference geometry (nodal positions)
        LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xrefe;
        DRT::Node** nodes = Nodes();
        for (int i = 0; i < NUMNOD_SOTET10; ++i)
        {
          xrefe(i, 0) = nodes[i]->X()[0];
          xrefe(i, 1) = nodes[i]->X()[1];
          xrefe(i, 2) = nodes[i]->X()[2];
        }

        // deformation gradient = identity tensor (geometrically linear case!)
        LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> defgrd(false);
        for (int i = 0; i < NUMDIM_SOTET10; ++i) defgrd(i, i) = 1;

        // use of 11 GP for the errornorm computation

        //----------------------------------------------------------------
        // loop over all Gauss points
        //----------------------------------------------------------------
        for (int gp = 0; gp < NUMGPT_MASS_SOTET10; gp++)
        {
          // Gauss weights and Jacobian determinant
          double fac = detJ_mass_[gp] * weights[gp];

          // Gauss point in reference configuration
          LINALG::Matrix<NUMDIM_SOTET10, 1> xgp(true);
          for (int k = 0; k < NUMDIM_SOTET10; ++k)
            for (int n = 0; n < NUMNOD_SOTET10; ++n) xgp(k, 0) += (vals[gp])(n)*xrefe(n, k);

          //**************************************************************
          // get analytical solution
          LINALG::Matrix<NUMDIM_SOTET10, 1> uanalyt(true);
          LINALG::Matrix<MAT::NUM_STRESS_3D, 1> strainanalyt(true);
          LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> derivanalyt(true);

          CONTACT::AnalyticalSolutions3D(xgp, uanalyt, strainanalyt, derivanalyt);
          //**************************************************************

          //--------------------------------------------------------------
          // (1) L2 norm
          //--------------------------------------------------------------

          // compute displacements at GP
          LINALG::Matrix<NUMDIM_SOTET10, 1> ugp(true);
          for (int k = 0; k < NUMDIM_SOTET10; ++k)
            for (int n = 0; n < NUMNOD_SOTET10; ++n)
              ugp(k, 0) += (vals[gp])(n)*nodaldisp(NODDOF_SOTET10 * n + k, 0);

          // displacement error
          LINALG::Matrix<NUMDIM_SOTET10, 1> uerror(true);
          for (int k = 0; k < NUMDIM_SOTET10; ++k) uerror(k, 0) = uanalyt(k, 0) - ugp(k, 0);

          // compute GP contribution to L2 error norm
          l2norm += fac * uerror.Dot(uerror);

          //--------------------------------------------------------------
          // (2) H1 norm
          //--------------------------------------------------------------

          // compute derivatives N_XYZ at GP w.r.t. material coordinates
          // by N_XYZ = J^-1 * N_rst
          LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10> N_XYZ(true);
          N_XYZ.Multiply(invJ_mass_[gp], derivs[gp]);

          // compute partial derivatives at GP
          LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> derivgp(true);
          for (int l = 0; l < NUMDIM_SOTET10; ++l)
            for (int m = 0; m < NUMDIM_SOTET10; ++m)
              for (int k = 0; k < NUMNOD_SOTET10; ++k)
                derivgp(l, m) += N_XYZ(m, k) * nodaldisp(NODDOF_SOTET10 * k + l, 0);

          // derivative error
          LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> deriverror(true);
          for (int k = 0; k < NUMDIM_SOTET10; ++k)
            for (int m = 0; m < NUMDIM_SOTET10; ++m)
              deriverror(k, m) = derivanalyt(k, m) - derivgp(k, m);

          // compute GP contribution to H1 error norm
          h1norm += fac * deriverror.Dot(deriverror);
          h1norm += fac * uerror.Dot(uerror);

          //--------------------------------------------------------------
          // (3) Energy norm
          //--------------------------------------------------------------

          // compute linear B-operator
          LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOTET10> bop;
          for (int i = 0; i < NUMNOD_SOTET10; ++i)
          {
            bop(0, NODDOF_SOTET10 * i + 0) = N_XYZ(0, i);
            bop(0, NODDOF_SOTET10 * i + 1) = 0.0;
            bop(0, NODDOF_SOTET10 * i + 2) = 0.0;
            bop(1, NODDOF_SOTET10 * i + 0) = 0.0;
            bop(1, NODDOF_SOTET10 * i + 1) = N_XYZ(1, i);
            bop(1, NODDOF_SOTET10 * i + 2) = 0.0;
            bop(2, NODDOF_SOTET10 * i + 0) = 0.0;
            bop(2, NODDOF_SOTET10 * i + 1) = 0.0;
            bop(2, NODDOF_SOTET10 * i + 2) = N_XYZ(2, i);

            bop(3, NODDOF_SOTET10 * i + 0) = N_XYZ(1, i);
            bop(3, NODDOF_SOTET10 * i + 1) = N_XYZ(0, i);
            bop(3, NODDOF_SOTET10 * i + 2) = 0.0;
            bop(4, NODDOF_SOTET10 * i + 0) = 0.0;
            bop(4, NODDOF_SOTET10 * i + 1) = N_XYZ(2, i);
            bop(4, NODDOF_SOTET10 * i + 2) = N_XYZ(1, i);
            bop(5, NODDOF_SOTET10 * i + 0) = N_XYZ(2, i);
            bop(5, NODDOF_SOTET10 * i + 1) = 0.0;
            bop(5, NODDOF_SOTET10 * i + 2) = N_XYZ(0, i);
          }

          // compute linear strain at GP
          LINALG::Matrix<MAT::NUM_STRESS_3D, 1> straingp(true);
          straingp.Multiply(bop, nodaldisp);

          // strain error
          LINALG::Matrix<MAT::NUM_STRESS_3D, 1> strainerror(true);
          for (int k = 0; k < MAT::NUM_STRESS_3D; ++k)
            strainerror(k, 0) = strainanalyt(k, 0) - straingp(k, 0);

          // compute stress vector and constitutive matrix
          LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> cmat(true);
          LINALG::Matrix<MAT::NUM_STRESS_3D, 1> stress(true);
          params.set<int>("gp", gp);
          SolidMaterial()->Evaluate(&defgrd, &strainerror, params, &stress, &cmat, Id());

          // compute GP contribution to energy error norm
          energynorm += fac * stress.Dot(strainerror);

          // std::cout << "UAnalytical:      " << ugp << std::endl;
          // std::cout << "UDiscrete:        " << uanalyt << std::endl;
          // std::cout << "StrainAnalytical: " << strainanalyt << std::endl;
          // std::cout << "StrainDiscrete:   " << straingp << std::endl;
          // std::cout << "DerivAnalytical:  " << derivanalyt << std::endl;
          // std::cout << "DerivDiscrete:    " << derivgp << std::endl;
        }
        //----------------------------------------------------------------

        // return results
        elevec1_epetra(0) = l2norm;
        elevec1_epetra(1) = h1norm;
        elevec1_epetra(2) = energynorm;
      }
      else
        dserror("ERROR: Error norms only implemented for SVK material");

      ;  // there is nothing to do here at the moment
    }
    break;

      // evaluate stresses and strains at gauss points and store gpstresses in map <EleId,
      // gpstresses >
    case calc_global_gpstresses_map:
    {
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
        const Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> gpstressmap =
            params.get<Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>>>(
                "gpstressmap", Teuchos::null);
        if (gpstressmap == Teuchos::null)
          dserror("no gp stress map available for writing gpstresses");
        const Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> gpstrainmap =
            params.get<Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>>>(
                "gpstrainmap", Teuchos::null);
        if (gpstrainmap == Teuchos::null)
          dserror("no gp strain map available for writing gpstrains");
        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
        std::vector<double> myres(lm.size());
        DRT::UTILS::ExtractMyValues(*res, myres, lm);
        LINALG::Matrix<NUMGPT_SOTET10, MAT::NUM_STRESS_3D> stress;
        LINALG::Matrix<NUMGPT_SOTET10, MAT::NUM_STRESS_3D> strain;
        INPAR::STR::StressType iostress =
            DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
        INPAR::STR::StrainType iostrain =
            DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);

        std::vector<double> mydispmat(lm.size(), 0.0);

        // if a linear analysis is desired
        if (kintype_ == INPAR::STR::kinem_linear)
        {
          dserror("Linear case not implemented");
        }


        else
        {
          if (pstype_ == INPAR::STR::prestress_id && time_ <= pstime_)  // inverse design analysis
            dserror("Inverse Design not implemented for SOTET10");
          // invdesign_->so_tet10_nlnstiffmass(this,lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,iostress,iostrain);

          else  // standard analysis
            so_tet10_nlnstiffmass(lm, mydisp, NULL, NULL, myres, mydispmat, NULL, NULL, NULL, NULL,
                NULL, &stress, &strain, params, iostress, iostrain);
        }
        // add stresses to global map
        // get EleID Id()
        int gid = Id();
        Teuchos::RCP<Epetra_SerialDenseMatrix> gpstress =
            Teuchos::rcp(new Epetra_SerialDenseMatrix);
        gpstress->Shape(NUMGPT_SOTET10, MAT::NUM_STRESS_3D);

        // move stresses to serial dense matrix
        for (int i = 0; i < NUMGPT_SOTET10; i++)
        {
          for (int j = 0; j < MAT::NUM_STRESS_3D; j++)
          {
            (*gpstress)(i, j) = stress(i, j);
          }
        }

        // strains
        Teuchos::RCP<Epetra_SerialDenseMatrix> gpstrain =
            Teuchos::rcp(new Epetra_SerialDenseMatrix);
        gpstrain->Shape(NUMGPT_SOTET10, MAT::NUM_STRESS_3D);

        // move stresses to serial dense matrix
        for (int i = 0; i < NUMGPT_SOTET10; i++)
        {
          for (int j = 0; j < MAT::NUM_STRESS_3D; j++)
          {
            (*gpstrain)(i, j) = strain(i, j);
          }
        }

        // add to map
        (*gpstressmap)[gid] = gpstress;
        (*gpstrainmap)[gid] = gpstrain;

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


    default:
      dserror("Unknown type of action for SolidTet10");
      break;
  }
  return 0;
}



/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)              |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet10::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  // get values and switches from the condition
  const std::vector<int>* onoff = condition.Get<std::vector<int>>("onoff");
  const std::vector<double>* val = condition.Get<std::vector<double>>("val");

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  const double time = params.get("total time", -1.0);

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < NUMDIM_SOTET10)
    dserror("Fewer functions or curves defined than the element has dofs.");

  for (int checkdof = NUMDIM_SOTET10; checkdof < int(onoff->size()); ++checkdof)
  {
    if ((*onoff)[checkdof] != 0)
      dserror("Number of Dimensions in Neumann_Evalutaion is 3. Further DoFs are not considered.");
  }

  // (SPATIAL) FUNCTION BUSINESS
  const std::vector<int>* funct = condition.Get<std::vector<int>>("funct");
  LINALG::Matrix<NUMDIM_SOTET10, 1> xrefegp(false);
  bool havefunct = false;
  if (funct)
    for (int dim = 0; dim < NUMDIM_SOTET10; dim++)
      if ((*funct)[dim] > 0) havefunct = havefunct or true;

  /* ============================================================================*
  ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for TET_10 with 4 GAUSS POINTS*
  ** ============================================================================*/
  const static std::vector<LINALG::Matrix<NUMNOD_SOTET10, 1>> shapefcts = so_tet10_4gp_shapefcts();
  const static std::vector<LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs =
      so_tet10_4gp_derivs();
  const static std::vector<double> gpweights = so_tet10_4gp_weights();
  /* ============================================================================*/

  // update element geometry
  LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xrefe;  // material coord. of element
  DRT::Node** nodes = Nodes();
  for (int i = 0; i < NUMNOD_SOTET10; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }

  /* ================================================= Loop over Gauss Points */
  for (int gp = 0; gp < NUMGPT_SOTET10; ++gp)
  {
    // compute the Jacobian matrix
    LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> jac;
    jac.Multiply(derivs[gp], xrefe);

    // compute determinant of Jacobian
    const double detJ = jac.Determinant();
    if (detJ == 0.0)
      dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0)
      dserror("NEGATIVE JACOBIAN DETERMINANT");

    // material/reference co-ordinates of Gauss point
    if (havefunct)
    {
      for (int dim = 0; dim < NUMDIM_SOTET10; dim++)
      {
        xrefegp(dim) = 0.0;
        for (int nodid = 0; nodid < NUMNOD_SOTET10; ++nodid)
          xrefegp(dim) += shapefcts[gp](nodid) * xrefe(nodid, dim);
      }
    }

    // integration factor
    double fac = gpweights[gp] * detJ;
    // distribute/add over element load vector
    for (int dim = 0; dim < NUMDIM_SOTET10; dim++)
    {
      if ((*onoff)[dim])
      {
        // function evaluation
        const int functnum = (funct) ? (*funct)[dim] : -1;
        const double functfac =
            (functnum > 0)
                ? DRT::Problem::Instance()->Funct(functnum - 1).Evaluate(dim, xrefegp.A(), time)
                : 1.0;
        const double dim_fac = (*val)[dim] * fac * functfac;
        for (int nodid = 0; nodid < NUMNOD_SOTET10; ++nodid)
        {
          elevec1[nodid * NUMDIM_SOTET10 + dim] += shapefcts[gp](nodid) * dim_fac;
        }
      }
    }

  } /* ==================================================== end of Loop over GP */

  return 0;
}  // DRT::ELEMENTS::So_tet10::EvaluateNeumann


/*----------------------------------------------------------------------*
 |  init the element jacobian mapping (protected)                       |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::InitJacobianMapping()
{
  const static std::vector<LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs4gp =
      so_tet10_4gp_derivs();
  const static std::vector<LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs11gp =
      so_tet10_11gp_derivs();

  LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xrefe;
  for (int i = 0; i < NUMNOD_SOTET10; ++i)
  {
    xrefe(i, 0) = Nodes()[i]->X()[0];
    xrefe(i, 1) = Nodes()[i]->X()[1];
    xrefe(i, 2) = Nodes()[i]->X()[2];
  }
  // Initialize for stiffness integration with 4 GPs
  invJ_.resize(NUMGPT_SOTET10);
  detJ_.resize(NUMGPT_SOTET10);
  for (int gp = 0; gp < NUMGPT_SOTET10; ++gp)
  {
    invJ_[gp].Multiply(derivs4gp[gp], xrefe);
    detJ_[gp] = invJ_[gp].Invert();
    if (detJ_[gp] == 0.0)
      dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ_[gp] < 0.0)
      dserror("NEGATIVE JACOBIAN DETERMINANT");

    if (pstype_ == INPAR::STR::prestress_mulf && pstime_ >= time_)
      if (!(prestress_->IsInit()))
        prestress_->MatrixtoStorage(gp, invJ_[gp], prestress_->JHistory());
  }

  if (pstype_ == INPAR::STR::prestress_mulf && pstime_ >= time_) prestress_->IsInit() = true;

  // Initialize for mass integration with 10 GPs

  invJ_mass_.resize(NUMGPT_MASS_SOTET10);
  detJ_mass_.resize(NUMGPT_MASS_SOTET10);
  for (int gp = 0; gp < NUMGPT_MASS_SOTET10; ++gp)
  {
    invJ_mass_[gp].Multiply(derivs11gp[gp], xrefe);
    detJ_mass_[gp] = invJ_mass_[gp].Invert();
    if (detJ_mass_[gp] == 0.0)
      dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ_mass_[gp] < 0.0)
      dserror("NEGATIVE JACOBIAN DETERMINANT");
  }
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                                      |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::so_tet10_nlnstiffmass(std::vector<int>& lm,  // location matrix
    std::vector<double>& disp,                                             // current displacements
    std::vector<double>* vel,                                              // current velocities
    std::vector<double>* acc,                                              // current accelerations
    std::vector<double>& residual,                                         // current residual displ
    std::vector<double>& dispmat,                                 // current material displacements
    LINALG::Matrix<NUMDOF_SOTET10, NUMDOF_SOTET10>* stiffmatrix,  // element stiffness matrix
    LINALG::Matrix<NUMDOF_SOTET10, NUMDOF_SOTET10>* massmatrix,   // element mass matrix
    LINALG::Matrix<NUMDOF_SOTET10, 1>* force,                     // element internal force vector
    LINALG::Matrix<NUMDOF_SOTET10, 1>* forceinert,                // element inertial force vector
    LINALG::Matrix<NUMDOF_SOTET10, 1>* force_str,                 // element structural force vector
    LINALG::Matrix<NUMGPT_SOTET10, MAT::NUM_STRESS_3D>* elestress,  // stresses at GP
    LINALG::Matrix<NUMGPT_SOTET10, MAT::NUM_STRESS_3D>* elestrain,  // strains at GP
    Teuchos::ParameterList& params,         // algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,  // stress output option
    const INPAR::STR::StrainType iostrain)  // strain output option
{
  /* ============================================================================*
  ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for TET_10 with 4 GAUSS POINTS*
  ** ============================================================================*/
  const static std::vector<LINALG::Matrix<NUMNOD_SOTET10, 1>> shapefcts_4gp =
      so_tet10_4gp_shapefcts();
  const static std::vector<LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs_4gp =
      so_tet10_4gp_derivs();
  const static std::vector<double> gpweights_4gp = so_tet10_4gp_weights();
  /* ============================================================================*/

  // update element geometry
  LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xrefe;  // material coord. of element
  LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xcurr;  // current  coord. of element
  LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xdisp;
  DRT::Node** nodes = Nodes();
  for (int i = 0; i < NUMNOD_SOTET10; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];

    xcurr(i, 0) = xrefe(i, 0) + disp[i * NODDOF_SOTET10 + 0];
    xcurr(i, 1) = xrefe(i, 1) + disp[i * NODDOF_SOTET10 + 1];
    xcurr(i, 2) = xrefe(i, 2) + disp[i * NODDOF_SOTET10 + 2];

    if (pstype_ == INPAR::STR::prestress_mulf)
    {
      xdisp(i, 0) = disp[i * NODDOF_SOTET10 + 0];
      xdisp(i, 1) = disp[i * NODDOF_SOTET10 + 1];
      xdisp(i, 2) = disp[i * NODDOF_SOTET10 + 2];
    }
  }

  // interpolate fibers to gp if fiber nodes are defined
  std::vector<LINALG::Matrix<NUMDIM_SOTET10, 1>> gpfiber1(
      NUMNOD_SOTET10, LINALG::Matrix<NUMDIM_SOTET10, 1>(true));
  std::vector<LINALG::Matrix<NUMDIM_SOTET10, 1>> gpfiber2(
      NUMNOD_SOTET10, LINALG::Matrix<NUMDIM_SOTET10, 1>(true));
  DRT::FIBER::FiberNode* fnode = dynamic_cast<DRT::FIBER::FiberNode*>(nodes[0]);
  if (fnode)
    DRT::ELEMENTS::UTILS::NodalFiber<DRT::Element::tet10>(nodes, shapefcts_4gp, gpfiber1, gpfiber2);

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> defgrd(false);
  for (int gp = 0; gp < NUMGPT_SOTET10; ++gp)
  {
    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJ_[gp], derivs_4gp[gp]);
    double detJ = detJ_[gp];

    if (pstype_ == INPAR::STR::prestress_mulf)
    {
      // get Jacobian mapping wrt to the stored configuration
      LINALG::Matrix<3, 3> invJdef;
      prestress_->StoragetoMatrix(gp, invJdef, prestress_->JHistory());
      // get derivatives wrt to last spatial configuration
      LINALG::Matrix<3, 10> N_xyz;
      N_xyz.Multiply(invJdef, derivs_4gp[gp]);

      // build multiplicative incremental defgrd
      defgrd.MultiplyTT(xdisp, N_xyz);
      defgrd(0, 0) += 1.0;
      defgrd(1, 1) += 1.0;
      defgrd(2, 2) += 1.0;

      // get stored old incremental F
      LINALG::Matrix<3, 3> Fhist;
      prestress_->StoragetoMatrix(gp, Fhist, prestress_->FHistory());

      // build total defgrd = delta F * F_old
      LINALG::Matrix<3, 3> Fnew;
      Fnew.Multiply(defgrd, Fhist);
      defgrd = Fnew;
    }
    else
    {
      defgrd.MultiplyTT(xcurr, N_XYZ);
    }

    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain_epetra(MAT::NUM_STRESS_3D);
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(glstrain_epetra.A(), true);
    glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
    glstrain(3) = cauchygreen(0, 1);
    glstrain(4) = cauchygreen(1, 2);
    glstrain(5) = cauchygreen(2, 0);

    // return gp strains (only in case of stress/strain output)
    switch (iostrain)
    {
      case INPAR::STR::strain_gl:
      {
        if (elestrain == NULL) dserror("strain data not available");
        for (int i = 0; i < 3; ++i) (*elestrain)(gp, i) = glstrain(i);
        for (int i = 3; i < 6; ++i) (*elestrain)(gp, i) = 0.5 * glstrain(i);
      }
      break;
      case INPAR::STR::strain_ea:
      {
        if (elestrain == NULL) dserror("strain data not available");
        // rewriting Green-Lagrange strains in matrix format
        LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> gl;
        gl(0, 0) = glstrain(0);
        gl(0, 1) = 0.5 * glstrain(3);
        gl(0, 2) = 0.5 * glstrain(5);
        gl(1, 0) = gl(0, 1);
        gl(1, 1) = glstrain(1);
        gl(1, 2) = 0.5 * glstrain(4);
        gl(2, 0) = gl(0, 2);
        gl(2, 1) = gl(1, 2);
        gl(2, 2) = glstrain(2);

        // inverse of deformation gradient
        LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> invdefgrd;
        invdefgrd.Invert(defgrd);

        LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> temp;
        LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> euler_almansi;
        temp.Multiply(gl, invdefgrd);
        euler_almansi.MultiplyTN(invdefgrd, temp);

        (*elestrain)(gp, 0) = euler_almansi(0, 0);
        (*elestrain)(gp, 1) = euler_almansi(1, 1);
        (*elestrain)(gp, 2) = euler_almansi(2, 2);
        (*elestrain)(gp, 3) = euler_almansi(0, 1);
        (*elestrain)(gp, 4) = euler_almansi(1, 2);
        (*elestrain)(gp, 5) = euler_almansi(0, 2);
      }
      break;
      case INPAR::STR::strain_none:
        break;
      default:
        dserror("requested strain type not available");
        break;
    }

    /* non-linear B-operator (may so be called, meaning
    ** of B-operator is not so sharp in the non-linear realm) *
    ** B = F . Bl *
    **
    **      [ ... | F_11*N_{,1}^k  F_21*N_{,1}^k  F_31*N_{,1}^k | ... ]
    **      [ ... | F_12*N_{,2}^k  F_22*N_{,2}^k  F_32*N_{,2}^k | ... ]
    **      [ ... | F_13*N_{,3}^k  F_23*N_{,3}^k  F_33*N_{,3}^k | ... ]
    ** B =  [ ~~~   ~~~~~~~~~~~~~  ~~~~~~~~~~~~~  ~~~~~~~~~~~~~   ~~~ ]
    **      [       F_11*N_{,2}^k+F_12*N_{,1}^k                       ]
    **      [ ... |          F_21*N_{,2}^k+F_22*N_{,1}^k        | ... ]
    **      [                       F_31*N_{,2}^k+F_32*N_{,1}^k       ]
    **      [                                                         ]
    **      [       F_12*N_{,3}^k+F_13*N_{,2}^k                       ]
    **      [ ... |          F_22*N_{,3}^k+F_23*N_{,2}^k        | ... ]
    **      [                       F_32*N_{,3}^k+F_33*N_{,2}^k       ]
    **      [                                                         ]
    **      [       F_13*N_{,1}^k+F_11*N_{,3}^k                       ]
    **      [ ... |          F_23*N_{,1}^k+F_21*N_{,3}^k        | ... ]
    **      [                       F_33*N_{,1}^k+F_31*N_{,3}^k       ]
    */
    LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOTET10> bop;
    for (int i = 0; i < NUMNOD_SOTET10; ++i)
    {
      bop(0, NODDOF_SOTET10 * i + 0) = defgrd(0, 0) * N_XYZ(0, i);
      bop(0, NODDOF_SOTET10 * i + 1) = defgrd(1, 0) * N_XYZ(0, i);
      bop(0, NODDOF_SOTET10 * i + 2) = defgrd(2, 0) * N_XYZ(0, i);
      bop(1, NODDOF_SOTET10 * i + 0) = defgrd(0, 1) * N_XYZ(1, i);
      bop(1, NODDOF_SOTET10 * i + 1) = defgrd(1, 1) * N_XYZ(1, i);
      bop(1, NODDOF_SOTET10 * i + 2) = defgrd(2, 1) * N_XYZ(1, i);
      bop(2, NODDOF_SOTET10 * i + 0) = defgrd(0, 2) * N_XYZ(2, i);
      bop(2, NODDOF_SOTET10 * i + 1) = defgrd(1, 2) * N_XYZ(2, i);
      bop(2, NODDOF_SOTET10 * i + 2) = defgrd(2, 2) * N_XYZ(2, i);
      /* ~~~ */
      bop(3, NODDOF_SOTET10 * i + 0) = defgrd(0, 0) * N_XYZ(1, i) + defgrd(0, 1) * N_XYZ(0, i);
      bop(3, NODDOF_SOTET10 * i + 1) = defgrd(1, 0) * N_XYZ(1, i) + defgrd(1, 1) * N_XYZ(0, i);
      bop(3, NODDOF_SOTET10 * i + 2) = defgrd(2, 0) * N_XYZ(1, i) + defgrd(2, 1) * N_XYZ(0, i);
      bop(4, NODDOF_SOTET10 * i + 0) = defgrd(0, 1) * N_XYZ(2, i) + defgrd(0, 2) * N_XYZ(1, i);
      bop(4, NODDOF_SOTET10 * i + 1) = defgrd(1, 1) * N_XYZ(2, i) + defgrd(1, 2) * N_XYZ(1, i);
      bop(4, NODDOF_SOTET10 * i + 2) = defgrd(2, 1) * N_XYZ(2, i) + defgrd(2, 2) * N_XYZ(1, i);
      bop(5, NODDOF_SOTET10 * i + 0) = defgrd(0, 2) * N_XYZ(0, i) + defgrd(0, 0) * N_XYZ(2, i);
      bop(5, NODDOF_SOTET10 * i + 1) = defgrd(1, 2) * N_XYZ(0, i) + defgrd(1, 0) * N_XYZ(2, i);
      bop(5, NODDOF_SOTET10 * i + 2) = defgrd(2, 2) * N_XYZ(0, i) + defgrd(2, 0) * N_XYZ(2, i);
    }

    // call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> cmat(true);
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> stress(true);

    if (Material()->MaterialType() == INPAR::MAT::m_constraintmixture)
    {
      // gp reference coordinates
      LINALG::Matrix<NUMNOD_SOTET10, 1> funct(true);
      funct = shapefcts_4gp[gp];
      LINALG::Matrix<1, NUMDIM_SOTET10> point(true);
      point.MultiplyTN(funct, xrefe);
      params.set("gprefecoord", point);
    }

    if (fnode)
    {
      params.set("gpfiber1", gpfiber1[gp]);
      params.set("gpfiber2", gpfiber2[gp]);
    }

    params.set<int>("gp", gp);
    SolidMaterial()->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, Id());
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp stresses
    switch (iostress)
    {
      case INPAR::STR::stress_2pk:
      {
        if (elestress == NULL) dserror("stress data not available");
        for (int i = 0; i < MAT::NUM_STRESS_3D; ++i) (*elestress)(gp, i) = stress(i);
      }
      break;
      case INPAR::STR::stress_cauchy:
      {
        if (elestress == NULL) dserror("stress data not available");
        const double detF = defgrd.Determinant();

        LINALG::Matrix<3, 3> pkstress;
        pkstress(0, 0) = stress(0);
        pkstress(0, 1) = stress(3);
        pkstress(0, 2) = stress(5);
        pkstress(1, 0) = pkstress(0, 1);
        pkstress(1, 1) = stress(1);
        pkstress(1, 2) = stress(4);
        pkstress(2, 0) = pkstress(0, 2);
        pkstress(2, 1) = pkstress(1, 2);
        pkstress(2, 2) = stress(2);

        LINALG::Matrix<3, 3> temp;
        LINALG::Matrix<3, 3> cauchystress;
        temp.Multiply(1.0 / detF, defgrd, pkstress);
        cauchystress.MultiplyNT(temp, defgrd);

        (*elestress)(gp, 0) = cauchystress(0, 0);
        (*elestress)(gp, 1) = cauchystress(1, 1);
        (*elestress)(gp, 2) = cauchystress(2, 2);
        (*elestress)(gp, 3) = cauchystress(0, 1);
        (*elestress)(gp, 4) = cauchystress(1, 2);
        (*elestress)(gp, 5) = cauchystress(0, 2);
      }
      break;
      case INPAR::STR::stress_none:
        break;
      default:
        dserror("requested stress type not available");
        break;
    }

    double detJ_w = detJ * gpweights_4gp[gp];
    // update internal force vector
    if (force != NULL)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w, bop, stress, 1.0);
    }

    // update of stiffness matrix
    if (stiffmatrix != NULL)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<6, NUMDOF_SOTET10> cb;
      cb.Multiply(cmat, bop);
      stiffmatrix->MultiplyTN(detJ_w, bop, cb, 1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      LINALG::Matrix<6, 1> sfac(stress);  // auxiliary integrated stress
      sfac.Scale(detJ_w);                 // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
      std::vector<double> SmB_L(3);       // intermediate Sm.B_L
      // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
      for (int inod = 0; inod < NUMNOD_SOTET10; ++inod)
      {
        SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod) + sfac(5) * N_XYZ(2, inod);
        SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod) + sfac(4) * N_XYZ(2, inod);
        SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod) + sfac(2) * N_XYZ(2, inod);
        for (int jnod = 0; jnod < NUMNOD_SOTET10; ++jnod)
        {
          double bopstrbop = 0.0;  // intermediate value
          for (int idim = 0; idim < NUMDIM_SOTET10; ++idim)
            bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
          (*stiffmatrix)(3 * inod + 0, 3 * jnod + 0) += bopstrbop;
          (*stiffmatrix)(3 * inod + 1, 3 * jnod + 1) += bopstrbop;
          (*stiffmatrix)(3 * inod + 2, 3 * jnod + 2) += bopstrbop;
        }
      }  // end of integrate `geometric' stiffness******************************
    }
  } /* ==================================================== end of Loop over GP */

  /* ============================================================================*
  ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for TET_10 with 11 GAUSS POINTS*
  ** ============================================================================*/
  const static std::vector<LINALG::Matrix<NUMNOD_SOTET10, 1>> shapefcts_11gp =
      so_tet10_11gp_shapefcts();
  const static std::vector<LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs_11gp =
      so_tet10_11gp_derivs();
  const static std::vector<double> gpweights_11gp = so_tet10_11gp_weights();
  /* ============================================================================*/

  if (massmatrix != NULL)  // evaluate mass matrix +++++++++++++++++++++++++
  {
    // since material has only 4 GPs with density information!
    double density = 0.;
    double vol = 0.;

    for (int gp4 = 0; gp4 < NUMGPT_SOTET10; ++gp4)
    {
      vol += gpweights_4gp[gp4];
      density += gpweights_4gp[gp4] * Material()->Density(gp4);
    }
    density /= vol;

    // consistent mass matrix evaluated using a 11-point rule
    for (int gp = 0; gp < NUMGPT_MASS_SOTET10; gp++)
    {
      // integrate consistent mass matrix
      double detJ_mass = detJ_mass_[gp];
      double detJ_mass_w = detJ_mass * gpweights_11gp[gp];
      const double factor = detJ_mass_w * density;
      double ifactor, massfactor;
      for (int inod = 0; inod < NUMNOD_SOTET10; ++inod)
      {
        ifactor = shapefcts_11gp[gp](inod) * factor;
        for (int jnod = 0; jnod < NUMNOD_SOTET10; ++jnod)
        {
          massfactor = shapefcts_11gp[gp](jnod) * ifactor;  // intermediate factor
          (*massmatrix)(NUMDIM_SOTET10 * inod + 0, NUMDIM_SOTET10 * jnod + 0) += massfactor;
          (*massmatrix)(NUMDIM_SOTET10 * inod + 1, NUMDIM_SOTET10 * jnod + 1) += massfactor;
          (*massmatrix)(NUMDIM_SOTET10 * inod + 2, NUMDIM_SOTET10 * jnod + 2) += massfactor;
        }
      }

      // check for non constant mass matrix
      if (SolidMaterial()->VaryingDensity())
      {
        /*
         If the density, i.e. the mass matrix, is not constant, a linearization is neccessary.
         In general, the mass matrix can be dependent on the displacements, the velocities and the
         accelerations. We write all the additional terms into the mass matrix, hence, conversion
         from accelerations to velocities and displacements are needed. As those conversions depend
         on the time integration scheme, the factors are set within the respective time integrators
         and read from the parameter list inside the element (this is a little ugly...). */
        LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> nxyz(nxyz_);  // copy!

        double timintfac_dis = 0.0;
        double timintfac_vel = 0.0;
        if (IsParamsInterface())
        {
          timintfac_dis = StrParamsInterface().GetTimIntFactorDisp();
          timintfac_vel = StrParamsInterface().GetTimIntFactorVel();
        }
        else
        {
          timintfac_dis = params.get<double>("timintfac_dis");
          timintfac_vel = params.get<double>("timintfac_vel");
        }
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1> linmass_disp(true);
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1> linmass_vel(true);
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1> linmass(true);

        /*----------------------------------------------------------------------*
           the B-operator used is equivalent to the one used in hex8, this needs
           to be checked if it is ok, but from the mathematics point of view, the only
           thing that needed to be changed is the NUMDOF
           ----------------------------------------------------------------------*/
        /*
        ** B = F . Bl *
        **
        **      [ ... | F_11*N_{,1}^k  F_21*N_{,1}^k  F_31*N_{,1}^k | ... ]
        **      [ ... | F_12*N_{,2}^k  F_22*N_{,2}^k  F_32*N_{,2}^k | ... ]
        **      [ ... | F_13*N_{,3}^k  F_23*N_{,3}^k  F_33*N_{,3}^k | ... ]
        ** B =  [ ~~~   ~~~~~~~~~~~~~  ~~~~~~~~~~~~~  ~~~~~~~~~~~~~   ~~~ ]
        **      [       F_11*N_{,2}^k+F_12*N_{,1}^k                       ]
        **      [ ... |          F_21*N_{,2}^k+F_22*N_{,1}^k        | ... ]
        **      [                       F_31*N_{,2}^k+F_32*N_{,1}^k       ]
        **      [                                                         ]
        **      [       F_12*N_{,3}^k+F_13*N_{,2}^k                       ]
        **      [ ... |          F_22*N_{,3}^k+F_23*N_{,2}^k        | ... ]
        **      [                       F_32*N_{,3}^k+F_33*N_{,2}^k       ]
        **      [                                                         ]
        **      [       F_13*N_{,1}^k+F_11*N_{,3}^k                       ]
        **      [ ... |          F_23*N_{,1}^k+F_21*N_{,3}^k        | ... ]
        **      [                       F_33*N_{,1}^k+F_31*N_{,3}^k       ]
        */
        // size is 6x12
        LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOTET10> bop;
        for (int i = 0; i < NUMNOD_SOTET10; i++)
        {
          bop(0, NODDOF_SOTET10 * i + 0) = defgrd(0, 0) * nxyz(i, 0);
          bop(0, NODDOF_SOTET10 * i + 1) = defgrd(1, 0) * nxyz(i, 0);
          bop(0, NODDOF_SOTET10 * i + 2) = defgrd(2, 0) * nxyz(i, 0);
          bop(1, NODDOF_SOTET10 * i + 0) = defgrd(0, 1) * nxyz(i, 1);
          bop(1, NODDOF_SOTET10 * i + 1) = defgrd(1, 1) * nxyz(i, 1);
          bop(1, NODDOF_SOTET10 * i + 2) = defgrd(2, 1) * nxyz(i, 1);
          bop(2, NODDOF_SOTET10 * i + 0) = defgrd(0, 2) * nxyz(i, 2);
          bop(2, NODDOF_SOTET10 * i + 1) = defgrd(1, 2) * nxyz(i, 2);
          bop(2, NODDOF_SOTET10 * i + 2) = defgrd(2, 2) * nxyz(i, 2);
          /* ~~~ */
          bop(3, NODDOF_SOTET10 * i + 0) = defgrd(0, 0) * nxyz(i, 1) + defgrd(0, 1) * nxyz(i, 0);
          bop(3, NODDOF_SOTET10 * i + 1) = defgrd(1, 0) * nxyz(i, 1) + defgrd(1, 1) * nxyz(i, 0);
          bop(3, NODDOF_SOTET10 * i + 2) = defgrd(2, 0) * nxyz(i, 1) + defgrd(2, 1) * nxyz(i, 0);
          bop(4, NODDOF_SOTET10 * i + 0) = defgrd(0, 1) * nxyz(i, 2) + defgrd(0, 2) * nxyz(i, 1);
          bop(4, NODDOF_SOTET10 * i + 1) = defgrd(1, 1) * nxyz(i, 2) + defgrd(1, 2) * nxyz(i, 1);
          bop(4, NODDOF_SOTET10 * i + 2) = defgrd(2, 1) * nxyz(i, 2) + defgrd(2, 2) * nxyz(i, 1);
          bop(5, NODDOF_SOTET10 * i + 0) = defgrd(0, 2) * nxyz(i, 0) + defgrd(0, 0) * nxyz(i, 2);
          bop(5, NODDOF_SOTET10 * i + 1) = defgrd(1, 2) * nxyz(i, 0) + defgrd(1, 0) * nxyz(i, 2);
          bop(5, NODDOF_SOTET10 * i + 2) = defgrd(2, 2) * nxyz(i, 0) + defgrd(2, 0) * nxyz(i, 2);
        }


        // Right Cauchy-Green tensor = F^T * F
        LINALG::Matrix<NUMDIM_SOTET10, NUMDIM_SOTET10> cauchygreen;
        cauchygreen.MultiplyTN(defgrd, defgrd);

        // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
        // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
        Epetra_SerialDenseVector glstrain_epetra(MAT::NUM_STRESS_3D);
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(glstrain_epetra.A(), true);
        glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
        glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
        glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
        glstrain(3) = cauchygreen(0, 1);
        glstrain(4) = cauchygreen(1, 2);
        glstrain(5) = cauchygreen(2, 0);

        // evaluate derivative of mass w.r.t. to right cauchy green tensor
        SolidMaterial()->EvaluateNonLinMass(
            &defgrd, &glstrain, params, &linmass_disp, &linmass_vel, Id());

        // multiply by 2.0 to get derivative w.r.t green lagrange strains and multiply by time
        // integration factor
        linmass_disp.Scale(2.0 * timintfac_dis);
        linmass_vel.Scale(2.0 * timintfac_vel);
        linmass.Update(1.0, linmass_disp, 1.0, linmass_vel, 0.0);

        // evaluate accelerations at time n+1 at gauss point
        LINALG::Matrix<NUMDIM_SOTET10, 1> myacc(true);
        for (int idim = 0; idim < NUMDIM_SOTET10; ++idim)
          for (int inod = 0; inod < NUMNOD_SOTET10; ++inod)
            myacc(idim) += shapefcts_11gp[gp](inod) * (*acc)[idim + (inod * NUMDIM_SOTET10)];

        if (stiffmatrix != NULL)
        {
          // integrate linearisation of mass matrix
          //(B^T . d\rho/d disp . a) * detJ * w(gp)
          LINALG::Matrix<1, NUMDOF_SOTET10> cb;
          cb.MultiplyTN(linmass_disp, bop);
          for (int inod = 0; inod < NUMNOD_SOTET10; ++inod)
          {
            double factor = detJ_mass_w * shapefcts_11gp[gp](inod);
            for (int idim = 0; idim < NUMDIM_SOTET10; ++idim)
            {
              double massfactor = factor * myacc(idim);
              for (int jnod = 0; jnod < NUMNOD_SOTET10; ++jnod)
                for (int jdim = 0; jdim < NUMDIM_SOTET10; ++jdim)
                  (*massmatrix)(inod * NUMDIM_SOTET10 + idim, jnod * NUMDIM_SOTET10 + jdim) +=
                      massfactor * cb(jnod * NUMDIM_SOTET10 + jdim);
            }
          }
        }

        // internal force vector without EAS terms
        if (forceinert != NULL)
        {
          // integrate nonlinear inertia force term
          for (int inod = 0; inod < NUMNOD_SOTET10; ++inod)
          {
            double forcefactor = shapefcts_11gp[gp](inod) * detJ_mass_w;
            for (int idim = 0; idim < NUMDIM_SOTET10; ++idim)
              (*forceinert)(inod * NUMDIM_SOTET10 + idim) += forcefactor * density * myacc(idim);
          }
        }
      }
    }

  }  // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++

  return;
}  // DRT::ELEMENTS::So_tet10::SOTET10_nlnstiffmass



/*----------------------------------------------------------------------*
 |  lump mass matrix                                                    |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::so_tet10_lumpmass(
    LINALG::Matrix<NUMDOF_SOTET10, NUMDOF_SOTET10>* emass)
{
  // lump mass matrix
  if (emass != NULL)
  {
    // we assume #elemat2 is a square matrix
    for (unsigned c = 0; c < (*emass).N(); ++c)  // parse columns
    {
      double d = 0.0;
      for (unsigned r = 0; r < (*emass).M(); ++r)  // parse rows
      {
        d += (*emass)(r, c);  // accumulate row entries
        (*emass)(r, c) = 0.0;
      }
      (*emass)(c, c) = d;  // apply sum of row entries on diagonal
    }
  }
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet10 Shape fcts at 4 Gauss Points                         |
 *----------------------------------------------------------------------*/
const std::vector<LINALG::Matrix<NUMNOD_SOTET10, 1>>
DRT::ELEMENTS::So_tet10::so_tet10_4gp_shapefcts()
{
  static std::vector<LINALG::Matrix<NUMNOD_SOTET10, 1>> shapefcts(NUMGPT_SOTET10);
  static bool shapefcts_done = false;
  if (shapefcts_done) return shapefcts;

  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_tet_4point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int gp = 0; gp < NUMGPT_SOTET10; gp++)
  {
    const double r = intpoints.qxg[gp][0];
    const double s = intpoints.qxg[gp][1];
    const double t = intpoints.qxg[gp][2];

    DRT::UTILS::shape_function_3D(shapefcts[gp], r, s, t, DRT::Element::tet10);
  }
  shapefcts_done = true;

  return shapefcts;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet10 Shape fct derivs at 4 Gauss Points                   |
 *----------------------------------------------------------------------*/
const std::vector<LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>>&
DRT::ELEMENTS::So_tet10::so_tet10_4gp_derivs()
{
  static std::vector<LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs(NUMGPT_SOTET10);
  static bool derivs_done = false;
  if (derivs_done) return derivs;

  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_tet_4point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int gp = 0; gp < NUMGPT_SOTET10; gp++)
  {
    const double r = intpoints.qxg[gp][0];
    const double s = intpoints.qxg[gp][1];
    const double t = intpoints.qxg[gp][2];

    DRT::UTILS::shape_function_3D_deriv1(derivs[gp], r, s, t, DRT::Element::tet10);
  }
  derivs_done = true;

  return derivs;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet10 Weights at 4 Gauss Points                            |
 *----------------------------------------------------------------------*/
const std::vector<double>& DRT::ELEMENTS::So_tet10::so_tet10_4gp_weights()
{
  static std::vector<double> weights(NUMGPT_SOTET10);
  static bool weights_done = false;
  if (weights_done) return weights;

  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_tet_4point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int gp = 0; gp < NUMGPT_SOTET10; gp++) weights[gp] = intpoints.qwgt[gp];
  weights_done = true;

  return weights;
}


/*----------------------------------------------------------------------*
 |  Evaluate Tet10 Shape fcts at 10 Gauss Points                        |
 *----------------------------------------------------------------------*/
const std::vector<LINALG::Matrix<NUMNOD_SOTET10, 1>>&
DRT::ELEMENTS::So_tet10::so_tet10_11gp_shapefcts()
{
  static std::vector<LINALG::Matrix<NUMNOD_SOTET10, 1>> shapefcts(NUMGPT_MASS_SOTET10);
  static bool shapefcts_done = false;
  if (shapefcts_done) return shapefcts;

  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_tet_11point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int gp = 0; gp < NUMGPT_MASS_SOTET10; gp++)
  {
    const double r = intpoints.qxg[gp][0];
    const double s = intpoints.qxg[gp][1];
    const double t = intpoints.qxg[gp][2];

    DRT::UTILS::shape_function_3D(shapefcts[gp], r, s, t, DRT::Element::tet10);
  }
  shapefcts_done = true;

  return shapefcts;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet10 Shape fct derivs at 10 Gauss Points                  |
 *----------------------------------------------------------------------*/
const std::vector<LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>>&
DRT::ELEMENTS::So_tet10::so_tet10_11gp_derivs()
{
  static std::vector<LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs(NUMGPT_MASS_SOTET10);
  static bool derivs_done = false;
  if (derivs_done) return derivs;

  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_tet_11point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int gp = 0; gp < NUMGPT_MASS_SOTET10; gp++)
  {
    const double r = intpoints.qxg[gp][0];
    const double s = intpoints.qxg[gp][1];
    const double t = intpoints.qxg[gp][2];

    DRT::UTILS::shape_function_3D_deriv1(derivs[gp], r, s, t, DRT::Element::tet10);
  }
  derivs_done = true;

  return derivs;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet10 Weights at 10 Gauss Points                           |
 *----------------------------------------------------------------------*/
const std::vector<double>& DRT::ELEMENTS::So_tet10::so_tet10_11gp_weights()
{
  static std::vector<double> weights(NUMGPT_MASS_SOTET10);
  static bool weights_done = false;
  if (weights_done) return weights;

  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_tet_11point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int gp = 0; gp < NUMGPT_MASS_SOTET10; gp++) weights[gp] = intpoints.qwgt[gp];
  weights_done = true;

  return weights;
}


/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet10Type::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So_tet10* actele = dynamic_cast<DRT::ELEMENTS::So_tet10*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_tet10* failed");
    actele->InitJacobianMapping();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  compute def gradient at every gaussian point (protected)            |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::DefGradient(const std::vector<double>& disp,
    Epetra_SerialDenseMatrix& gpdefgrd, DRT::ELEMENTS::PreStress& prestress)
{
  const static std::vector<LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs =
      so_tet10_4gp_derivs();

  // update element geometry
  LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xdisp;  // current  coord. of element
  for (int i = 0; i < NUMNOD_SOTET10; ++i)
  {
    xdisp(i, 0) = disp[i * NODDOF_SOTET10 + 0];
    xdisp(i, 1) = disp[i * NODDOF_SOTET10 + 1];
    xdisp(i, 2) = disp[i * NODDOF_SOTET10 + 2];
  }

  for (int gp = 0; gp < NUMGPT_SOTET10; ++gp)
  {
    // get Jacobian mapping wrt to the stored deformed configuration
    LINALG::Matrix<3, 3> invJdef;
    prestress.StoragetoMatrix(gp, invJdef, prestress.JHistory());

    // by N_XYZ = J^-1 * N_rst
    LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10> N_xyz;
    N_xyz.Multiply(invJdef, derivs[gp]);

    // build defgrd (independent of xrefe!)
    LINALG::Matrix<3, 3> defgrd;
    defgrd.MultiplyTT(xdisp, N_xyz);
    defgrd(0, 0) += 1.0;
    defgrd(1, 1) += 1.0;
    defgrd(2, 2) += 1.0;

    prestress.MatrixtoStorage(gp, defgrd, gpdefgrd);
  }
  return;
}


/*----------------------------------------------------------------------*
 |  compute Jac.mapping wrt deformed configuration (protected)          |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::UpdateJacobianMapping(
    const std::vector<double>& disp, DRT::ELEMENTS::PreStress& prestress)
{
  const static std::vector<LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10>> derivs =
      so_tet10_4gp_derivs();

  // get incremental disp
  LINALG::Matrix<NUMNOD_SOTET10, NUMDIM_SOTET10> xdisp;
  for (int i = 0; i < NUMNOD_SOTET10; ++i)
  {
    xdisp(i, 0) = disp[i * NODDOF_SOTET10 + 0];
    xdisp(i, 1) = disp[i * NODDOF_SOTET10 + 1];
    xdisp(i, 2) = disp[i * NODDOF_SOTET10 + 2];
  }

  LINALG::Matrix<3, 3> invJhist;
  LINALG::Matrix<3, 3> invJ;
  LINALG::Matrix<3, 3> defgrd;
  LINALG::Matrix<NUMDIM_SOTET10, NUMNOD_SOTET10> N_xyz;
  LINALG::Matrix<3, 3> invJnew;
  for (int gp = 0; gp < NUMGPT_SOTET10; ++gp)
  {
    // get the invJ old state
    prestress.StoragetoMatrix(gp, invJhist, prestress.JHistory());
    // get derivatives wrt to invJhist
    N_xyz.Multiply(invJhist, derivs[gp]);
    // build defgrd \partial x_new / \parial x_old , where x_old != X
    defgrd.MultiplyTT(xdisp, N_xyz);
    defgrd(0, 0) += 1.0;
    defgrd(1, 1) += 1.0;
    defgrd(2, 2) += 1.0;
    // make inverse of this defgrd
    defgrd.Invert();
    // push-forward of Jinv
    invJnew.MultiplyTN(defgrd, invJhist);
    // store new reference configuration
    prestress.MatrixtoStorage(gp, invJnew, prestress.JHistory());
  }  // for (int gp=0; gp<NUMGPT_SOTET10; ++gp)

  return;
}

/*----------------------------------------------------------------------*/
/*! \file

\brief Evaluate routines for Solid Hex8 element

\level 1

\maintainer Christoph Meier

*----------------------------------------------------------------------*/

#include "so_hex8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_elements.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/voigt_notation.H"
#include "../linalg/linalg_utils_densematrix_inverse.H"
#include "../linalg/linalg_utils_densematrix_eigen.H"
#include "../linalg/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/so3_material.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_mat/growthremodel_elasthyper.H"
#include "../drt_mat/constraintmixture.H"
#include "../drt_mat/thermostvenantkirchhoff.H"
#include "../drt_mat/thermoplastichyperelast.H"
#include "../drt_mat/robinson.H"

#include "../drt_contact/contact_analytical.H"
#include "../drt_patspec/patspec.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// inverse design object
#include "inversedesign.H"
#include "prestress.H"

#include "../drt_fluid_ele/fluid_ele_parameter_timint.H"
#include "../drt_structure_new/str_elements_paramsinterface.H"
#include "../drt_structure_new/str_enum_lists.H"

#include "so3_defines.H"
#include "so_hex8_determinant_analysis.H"

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_mat/thermomech_threephase.H"
#include "so_utils.H"

using VoigtMapping = UTILS::VOIGT::IndexMappings;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  // Check whether the solid material PostSetup() routine has already been called and call it if not
  EnsureMaterialPostSetup(params);

  SetParamsInterfacePtr(params);

  LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8> elemat1(elemat1_epetra.A(), true);
  LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8> elemat2(elemat2_epetra.A(), true);
  LINALG::Matrix<NUMDOF_SOH8, 1> elevec1(elevec1_epetra.A(), true);
  LINALG::Matrix<NUMDOF_SOH8, 1> elevec2(elevec2_epetra.A(), true);
  LINALG::Matrix<NUMDOF_SOH8, 1> elevec3(elevec3_epetra.A(), true);

  // start with "none"
  ELEMENTS::ActionType act = ELEMENTS::none;

  if (IsParamsInterface())
  {
    act = ParamsInterface().GetActionType();
  }
  else
  {
    // get the required action
    std::string action = params.get<std::string>("action", "none");
    act = ELEMENTS::String2ActionType(action);
  }

  // check for patient specific data
  PATSPEC::GetILTDistance(Id(), params, discretization);
  PATSPEC::GetLocalRadius(Id(), params, discretization);
  PATSPEC::GetInnerRadius(Id(), params, discretization);

  // what should the element do
  switch (act)
  {
    //==================================================================================
    // nonlinear stiffness and internal force vector
    case ELEMENTS::struct_calc_nlnstiff:
    case ELEMENTS::struct_calc_linstiff:
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
      LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* matptr = NULL;
      if (elemat1.IsInitialized()) matptr = &elemat1;

      std::vector<double> mydispmat(lm.size(), 0.0);
      if (structale_)
      {
        Teuchos::RCP<const Epetra_Vector> dispmat =
            discretization.GetState("material_displacement");
        DRT::UTILS::ExtractMyValues(*dispmat, mydispmat, lm);
      }

      if (pstype_ == INPAR::STR::prestress_id && time_ <= pstime_)  // inverse design analysis
        invdesign_->soh8_nlnstiffmass(this, lm, mydisp, myres, matptr, NULL, &elevec1, NULL, NULL,
            params, INPAR::STR::stress_none, INPAR::STR::strain_none);

      else  // standard analysis
        nlnstiffmass(lm, mydisp, NULL, NULL, myres, mydispmat, matptr, NULL, &elevec1, NULL,
            &elevec3, NULL, NULL, NULL, params, INPAR::STR::stress_none, INPAR::STR::strain_none,
            INPAR::STR::strain_none);

      break;
    }
    //==================================================================================
    // internal force vector only
    case ELEMENTS::struct_calc_internalforce:
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

      std::vector<double> mydispmat(lm.size(), 0.0);
      if (structale_)
      {
        Teuchos::RCP<const Epetra_Vector> dispmat =
            discretization.GetState("material_displacement");
        DRT::UTILS::ExtractMyValues(*dispmat, mydispmat, lm);
      }

      nlnstiffmass(lm, mydisp, NULL, NULL, myres, mydispmat, &myemat, NULL, &elevec1, NULL, NULL,
          NULL, NULL, NULL, params, INPAR::STR::stress_none, INPAR::STR::strain_none,
          INPAR::STR::strain_none);

      break;
    }
    //==================================================================================
    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case ELEMENTS::struct_calc_nlnstiffmass:
    case ELEMENTS::struct_calc_nlnstifflmass:
    case ELEMENTS::struct_calc_linstiffmass:
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
      if (structale_)
      {
        Teuchos::RCP<const Epetra_Vector> dispmat =
            discretization.GetState("material_displacement");
        ;
        DRT::UTILS::ExtractMyValues(*dispmat, mydispmat, lm);
      }

      if (pstype_ == INPAR::STR::prestress_id && time_ <= pstime_)  // inverse design analysis
        invdesign_->soh8_nlnstiffmass(this, lm, mydisp, myres, &elemat1, &elemat2, &elevec1, NULL,
            NULL, params, INPAR::STR::stress_none, INPAR::STR::strain_none);
      else  // standard analysis
        nlnstiffmass(lm, mydisp, &myvel, &myacc, myres, mydispmat, &elemat1, &elemat2, &elevec1,
            &elevec2, &elevec3, NULL, NULL, NULL, params, INPAR::STR::stress_none,
            INPAR::STR::strain_none, INPAR::STR::strain_none);

      if (act == ELEMENTS::struct_calc_nlnstifflmass) soh8_lumpmass(&elemat2);

      break;
    }
    //==================================================================================
    case struct_calc_mass_volume:
    {
      // declaration of variables
      double volume_ref = 0.0;
      double volume_mat = 0.0;
      double volume_cur = 0.0;
      double mass_ref = 0.0;
      double mass_mat = 0.0;
      double mass_cur = 0.0;

      // get displacements and extract values of this element
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state displacement vector");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      std::vector<double> mydispmat(lm.size(), 0.0);
      if (structale_)
      {
        Teuchos::RCP<const Epetra_Vector> dispmat =
            discretization.GetState("material_displacement");
        DRT::UTILS::ExtractMyValues(*dispmat, mydispmat, lm);
      }
      // reference and current geometry (nodal positions)
      LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xref;  // reference coord. of element
      LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xcur;  // current  coord. of element
      LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xmat;  // mat  coord. of element

      for (int k = 0; k < NUMNOD_SOH8; ++k)
      {
        xref(k, 0) = Nodes()[k]->X()[0];
        xref(k, 1) = Nodes()[k]->X()[1];
        xref(k, 2) = Nodes()[k]->X()[2];
        xcur(k, 0) = xref(k, 0) + mydisp[k * NODDOF_SOH8 + 0];
        xcur(k, 1) = xref(k, 1) + mydisp[k * NODDOF_SOH8 + 1];
        xcur(k, 2) = xref(k, 2) + mydisp[k * NODDOF_SOH8 + 2];

        // material displacements for structure with ale
        if (structale_ == true)
        {
          xmat(k, 0) = xref(k, 0) + mydispmat[k * NODDOF_SOH8 + 0];
          xmat(k, 1) = xref(k, 1) + mydispmat[k * NODDOF_SOH8 + 1];
          xmat(k, 2) = xref(k, 2) + mydispmat[k * NODDOF_SOH8 + 2];
        }
      }

      // safety check before the actual evaluation starts
      const double min_detJ_curr = soh8_get_min_det_jac_at_corners(xcur);
      if (min_detJ_curr <= 0.0)
      {
        soh8_error_handling(
            min_detJ_curr, params, __LINE__, STR::ELEMENTS::ele_error_determinant_at_corner);
        elevec1_epetra(0) = 0.0;
        return 1;
      }

      const static std::vector<LINALG::Matrix<NUMNOD_SOH8, 1>> shapefcts = soh8_shapefcts();
      const static std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = soh8_derivs();
      const static std::vector<double> gpweights = soh8_weights();

      // MAT ------------------------
      // build new jacobian mapping with respect to the material configuration
      int err = 0;
      if (structale_ == true)
      {
        err = InitJacobianMapping(mydispmat);
        if (err)
        {
          // reset class variable before leaving
          InitJacobianMapping();
          return err;
        }
      }

      std::vector<double> detJmat = detJ_;
      std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>> invJmat = invJ_;

      err = InitJacobianMapping(mydisp);
      if (err)
      {
        // reset class variable before leaving
        InitJacobianMapping();
        return err;
      }

      std::vector<double> detJcur = detJ_;
      std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>> invJcur = invJ_;

      InitJacobianMapping();

      std::vector<double> detJref = detJ_;
      std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>> invJref = invJ_;


      LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_XYZ;
      LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrd(true);

      //----------------------------------------------------------------
      // loop over all Gauss points
      //----------------------------------------------------------------
      for (unsigned ip = 0; ip < NUMGPT_SOH8; ++ip)
      {
        const double density = Material()->Density(ip);
        const double wgt = gpweights[ip];

        /*------------------------------------ integration factor  -------*/
        double fac = wgt * detJref[ip];
        volume_ref += fac;
        fac = wgt * detJref[ip] * density;
        mass_ref += fac;

        // MAT ------------------------
        if (structale_)
        {
          fac = wgt * detJmat[ip];
          volume_mat += fac;
          fac = wgt * detJmat[ip] * density;
          mass_mat += fac;

          N_XYZ.Multiply(invJmat[ip], derivs[ip]);
          defgrd.MultiplyTT(xcur, N_XYZ);
          double detFmat = defgrd.Determinant();

          /*------------------------------------ integration factor  -------*/
          fac = wgt * detJcur[ip];
          volume_cur += fac;
          fac = wgt * detJcur[ip] * density * 1 / detFmat;
          mass_cur += fac;
        }
        else
        {
          N_XYZ.Multiply(invJref[ip], derivs[ip]);
          defgrd.MultiplyTT(xcur, N_XYZ);
          double detFref = defgrd.Determinant();

          /*------------------------------------ integration factor  -------*/
          fac = wgt * detJcur[ip];
          volume_cur += fac;
          fac = wgt * detJcur[ip] * density * 1 / detFref;
          mass_cur += fac;
        }
      }

      //----------------------------------------------------------------

      // return results
      if (!structale_)
      {
        volume_mat = volume_ref;
        mass_mat = mass_ref;
      }

      elevec1(0) = volume_ref;
      elevec1(1) = volume_mat;
      elevec1(2) = volume_cur;
      elevec1(3) = mass_ref;
      elevec1(4) = mass_mat;
      elevec1(5) = mass_cur;

      break;
    }

    case ELEMENTS::analyse_jacobian_determinant:
    {
      // get displacements and extract values of this element
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state displacement vector");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      std::vector<double> mydispmat(lm.size(), 0.0);
      // reference and current geometry (nodal positions)
      LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> xcurr;  // current  coord. of element

      for (int k = 0; k < NUMNOD_SOH8; ++k)
      {
        xcurr(0, k) = Nodes()[k]->X()[0] + mydisp[k * NODDOF_SOH8 + 0];
        xcurr(1, k) = Nodes()[k]->X()[1] + mydisp[k * NODDOF_SOH8 + 1];
        xcurr(2, k) = Nodes()[k]->X()[2] + mydisp[k * NODDOF_SOH8 + 2];
      }

      Teuchos::RCP<So_Hex8_Determinant_Analysis> det_analyser =
          So_Hex8_Determinant_Analysis::create();
      if (not det_analyser->isValid(xcurr))
        soh8_error_handling(-1.0, params, __LINE__, STR::ELEMENTS::ele_error_determinant_analysis);

      break;
    }

    //==================================================================================
    // nonlinear stiffness, internal force vector (GEMM)
    case ELEMENTS::struct_calc_nlnstiff_gemm:
    {
      // need old displacement, current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> dispo = discretization.GetState("old displacement");
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (dispo == Teuchos::null || disp == Teuchos::null || res == Teuchos::null)
        dserror("Cannot get state vectors '(old) displacement' and/or residual");
      std::vector<double> mydispo(lm.size());
      DRT::UTILS::ExtractMyValues(*dispo, mydispo, lm);
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);

      // default: geometrically non-linear analysis with Total Lagrangean approach
      if (kintype_ == INPAR::STR::kinem_nonlinearTotLag)
      {
        soh8_nlnstiffmass_gemm(lm, mydispo, mydisp, myres, &elemat1, NULL, &elevec1, NULL, NULL,
            NULL, params, INPAR::STR::stress_none, INPAR::STR::strain_none,
            INPAR::STR::strain_none);
      }
      // special case: geometric linear
      else  // kintype_ != INPAR::STR::kinem_nonlinearTotLag)
      {
        dserror("ERROR: Generalized EMM only makes sense in nonlinear realm");
      }

      break;
    }
    //==================================================================================
    // recover elementwise stored quantities (e.g. EAS)
    case ELEMENTS::struct_calc_recover:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");

      if (disp == Teuchos::null || res == Teuchos::null)
        dserror(
            "Cannot get state vectors \"displacement\" "
            "and/or \"residual displacement\"");

      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);

      soh8_recover(lm, myres);
      /* ToDo Probably we have to recover the history information of some special
       * materials as well.                                 hiermeier 04/2016  */

      break;
    }
    //==================================================================================
    // evaluate stresses and strains at gauss points
    case ELEMENTS::struct_calc_stress:
    {
      // nothing to do for ghost elements
      if (discretization.Comm().MyPID() == Owner())
      {
        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
        Teuchos::RCP<std::vector<char>> stressdata = Teuchos::null;
        Teuchos::RCP<std::vector<char>> straindata = Teuchos::null;
        Teuchos::RCP<std::vector<char>> plstraindata = Teuchos::null;
        INPAR::STR::StressType iostress = INPAR::STR::stress_none;
        INPAR::STR::StrainType iostrain = INPAR::STR::strain_none;
        INPAR::STR::StrainType ioplstrain = INPAR::STR::strain_none;
        if (IsParamsInterface())
        {
          stressdata = StrParamsInterface().MutableStressDataPtr();
          straindata = StrParamsInterface().MutableStrainDataPtr();
          plstraindata = StrParamsInterface().MutablePlasticStrainDataPtr();

          iostress = StrParamsInterface().GetStressOutputType();
          iostrain = StrParamsInterface().GetStrainOutputType();
          ioplstrain = StrParamsInterface().GetPlasticStrainOutputType();
        }
        else
        {
          stressdata = params.get<Teuchos::RCP<std::vector<char>>>("stress", Teuchos::null);
          straindata = params.get<Teuchos::RCP<std::vector<char>>>("strain", Teuchos::null);
          iostress =
              DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
          iostrain =
              DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);
          // in case of small strain materials calculate plastic strains for post processing
          plstraindata = params.get<Teuchos::RCP<std::vector<char>>>("plstrain", Teuchos::null);
          ioplstrain = DRT::INPUT::get<INPAR::STR::StrainType>(
              params, "ioplstrain", INPAR::STR::strain_none);
        }
        if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
        if (stressdata == Teuchos::null) dserror("Cannot get 'stress' data");
        if (straindata == Teuchos::null) dserror("Cannot get 'strain' data");
        if (plstraindata == Teuchos::null) dserror("Cannot get 'plastic strain' data");
        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
        std::vector<double> myres(lm.size());
        DRT::UTILS::ExtractMyValues(*res, myres, lm);
        LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D> stress;
        LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D> strain;
        LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D> plstrain;

        std::vector<double> mydispmat(lm.size(), 0.0);
        if (structale_)
        {
          Teuchos::RCP<const Epetra_Vector> dispmat =
              discretization.GetState("material_displacement");
          ;
          DRT::UTILS::ExtractMyValues(*dispmat, mydispmat, lm);
        }

        if (pstype_ == INPAR::STR::prestress_id && time_ <= pstime_)  // inverse design analysis
          invdesign_->soh8_nlnstiffmass(this, lm, mydisp, myres, NULL, NULL, NULL, &stress, &strain,
              params, iostress, iostrain);

        else  // standard analysis
          nlnstiffmass(lm, mydisp, NULL, NULL, myres, mydispmat, NULL, NULL, NULL, NULL, NULL,
              &stress, &strain, &plstrain, params, iostress, iostrain, ioplstrain);

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

        {
          DRT::PackBuffer data;
          AddtoPack(data, plstrain);
          data.StartPacking();
          AddtoPack(data, plstrain);
          std::copy(data().begin(), data().end(), std::back_inserter(*plstraindata));
        }
      }
    }
    break;
    //==================================================================================
    // postprocess stresses/strains at gauss points
    // note that in the following, quantities are always referred to as
    // "stresses" etc. although they might also apply to strains
    // (depending on what this routine is called for from the post filter)
    case ELEMENTS::struct_postprocess_stress:
    {
      const Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> gpstressmap =
          params.get<Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>>>(
              "gpstressmap", Teuchos::null);
      if (gpstressmap == Teuchos::null)
        dserror("no gp stress/strain map available for postprocessing");
      std::string stresstype = params.get<std::string>("stresstype", "ndxyz");
      int gid = Id();
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
        int lid = elemap.LID(Id());
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
      }
      else
      {
        dserror("unknown type of stress/strain output on element level");
      }
    }
    break;
    //==================================================================================
    case ELEMENTS::struct_calc_eleload:
      dserror("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
      break;
    //==================================================================================
    case ELEMENTS::struct_calc_fsiload:
      dserror("Case not yet implemented");
      break;
    //==================================================================================
    case ELEMENTS::struct_calc_update_istep:
    {
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      Update_element(mydisp, params, Material());
    }
    break;
    //==================================================================================
    case ELEMENTS::struct_calc_reset_istep:
    {
      // do something with internal EAS, etc parameters
      if (eastype_ != soh8_easnone)
      {
        Epetra_SerialDenseMatrix* alpha =
            data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // Alpha_{n+1}
        Epetra_SerialDenseMatrix* alphao =
            data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");  // Alpha_n
        switch (eastype_)
        {
          case DRT::ELEMENTS::So_hex8::soh8_easfull:
            LINALG::DENSEFUNCTIONS::update<double, soh8_easfull, 1>(*alpha, *alphao);
            break;
          case DRT::ELEMENTS::So_hex8::soh8_easmild:
            LINALG::DENSEFUNCTIONS::update<double, soh8_easmild, 1>(*alpha, *alphao);
            break;
          case DRT::ELEMENTS::So_hex8::soh8_eassosh8:
            LINALG::DENSEFUNCTIONS::update<double, soh8_eassosh8, 1>(*alpha, *alphao);
            break;
          case DRT::ELEMENTS::So_hex8::soh8_easnone:
            break;
          default:
            dserror("Don't know what to do with EAS type %d", eastype_);
            break;
        }
      }
      // Reset of history (if needed)
      SolidMaterial()->ResetStep();
    }
    break;
    //==================================================================================
    case ELEMENTS::struct_calc_store_istep:
    {
      int timestep = params.get<int>("timestep", -1);

      if (timestep == -1) dserror("Provide timestep number to be stored");

      // EAS
      if (eastype_ != soh8_easnone) dserror("Storage of EAS stuff must be implemented first");

      // due to the multiplicativity and futility to redo prestress steps
      // other than the last one, no need to store/recover anything
      // ... but keep in mind
      if (pstype_ != INPAR::STR::prestress_none)
      {
      }

      // Material
      SolidMaterial()->StoreHistory(timestep);
    }
    break;
    //==================================================================================
    case ELEMENTS::struct_calc_recover_istep:
    {
      int timestep = params.get<int>("timestep", -1);

      if (timestep == -1) dserror("Provide timestep number of the timestep to be recovered");

      // EAS
      if (eastype_ != soh8_easnone) dserror("Recpvery of EAS stuff must be implemented first");

      // due to the multiplicativity and futility to redo prestress steps
      // other than the last one, no need to store/recover anything
      // ... but keep in mind
      if (pstype_ != INPAR::STR::prestress_none)
      {
      }

      // Material
      SolidMaterial()->SetHistory(timestep);
    }
    break;
    //==================================================================================
    case ELEMENTS::struct_calc_reset_all:
    {
      // Reset of history for materials
      SolidMaterial()->ResetAll(NUMGPT_SOH8);

      // Reset prestress
      if (pstype_ == INPAR::STR::prestress_mulf)
      {
        time_ = 0.0;
        LINALG::Matrix<3, 3> Id(true);
        Id(0, 0) = Id(1, 1) = Id(2, 2) = 1.0;
        for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
        {
          prestress_->MatrixtoStorage(gp, Id, prestress_->FHistory());
          prestress_->MatrixtoStorage(gp, invJ_[gp], prestress_->JHistory());
        }
      }
      if (pstype_ == INPAR::STR::prestress_id)
        dserror("Reset of Inverse Design not yet implemented");

      // reset EAS parameters:
      if (eastype_ != DRT::ELEMENTS::So_hex8::soh8_easnone)
      {
        soh8_easinit();
      }
    }
    break;
    //==================================================================================
    case ELEMENTS::struct_calc_energy:
    {
      // initialization of internal energy
      double intenergy = 0.0;

      // shape functions and Gauss weights
      const static std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = soh8_derivs();
      const static std::vector<double> weights = soh8_weights();

      // get displacements of this processor
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state displacement vector");

      // get displacements of this element
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      // update element geometry
      LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe;  // material coord. of element
      LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xcurr;  // current  coord. of element
      LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xdisp;

      DRT::Node** nodes = Nodes();
      for (int i = 0; i < NUMNOD_SOH8; ++i)
      {
        xrefe(i, 0) = nodes[i]->X()[0];
        xrefe(i, 1) = nodes[i]->X()[1];
        xrefe(i, 2) = nodes[i]->X()[2];

        xcurr(i, 0) = xrefe(i, 0) + mydisp[i * NODDOF_SOH8 + 0];
        xcurr(i, 1) = xrefe(i, 1) + mydisp[i * NODDOF_SOH8 + 1];
        xcurr(i, 2) = xrefe(i, 2) + mydisp[i * NODDOF_SOH8 + 2];

        if (pstype_ == INPAR::STR::prestress_mulf)
        {
          xdisp(i, 0) = mydisp[i * NODDOF_SOH8 + 0];
          xdisp(i, 1) = mydisp[i * NODDOF_SOH8 + 1];
          xdisp(i, 2) = mydisp[i * NODDOF_SOH8 + 2];
        }
      }

      // safety check before the actual evaluation starts
      const double min_detJ_curr = soh8_get_min_det_jac_at_corners(xcurr);
      if (min_detJ_curr <= 0.0)
      {
        soh8_error_handling(
            min_detJ_curr, params, __LINE__, STR::ELEMENTS::ele_error_determinant_at_corner);
        elevec1_epetra(0) = 0.0;
        return 0;
      }

      // prepare EAS data
      Epetra_SerialDenseMatrix* alpha = NULL;                         // EAS alphas
      std::vector<Epetra_SerialDenseMatrix>* M_GP = NULL;             // EAS matrix M at all GPs
      LINALG::SerialDenseMatrix M;                                    // EAS matrix M at current GP
      double detJ0;                                                   // detJ(origin)
      LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> T0invT;  // trafo matrix
      if (eastype_ != soh8_easnone)
      {
        alpha =
            data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // get alpha of previous iteration
        soh8_eassetup(&M_GP, detJ0, T0invT, xrefe);
      }

      // loop over all Gauss points
      for (unsigned gp = 0; gp < NUMGPT_SOH8; gp++)
      {
        // Gauss weights and Jacobian determinant
        double fac = detJ_[gp] * weights[gp];

        /* get the inverse of the Jacobian matrix which looks like:
        **            [ x_,r  y_,r  z_,r ]^-1
        **     J^-1 = [ x_,s  y_,s  z_,s ]
        **            [ x_,t  y_,t  z_,t ]
        */
        // compute derivatives N_XYZ at gp w.r.t. material coordinates
        // by N_XYZ = J^-1 * N_rst
        LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_XYZ(true);
        N_XYZ.Multiply(invJ_[gp], derivs[gp]);

        // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrd(true);

        if (pstype_ == INPAR::STR::prestress_id && pstime_ < time_)
        {
          dserror("Calc Energy not implemented for prestress id");
        }

        // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
        // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(true);

        if (pstype_ == INPAR::STR::prestress_mulf)
        {
          // get Jacobian mapping wrt to the stored configuration
          LINALG::Matrix<3, 3> invJdef;
          prestress_->StoragetoMatrix(gp, invJdef, prestress_->JHistory());
          // get derivatives wrt to last spatial configuration
          LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_xyz;
          N_xyz.Multiply(invJdef, derivs[gp]);

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

          // right Cauchy-Green tensor = F^T * F
          LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> cauchygreen;
          cauchygreen.MultiplyTN(defgrd, defgrd);

          glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
          glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
          glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
          glstrain(3) = cauchygreen(0, 1);
          glstrain(4) = cauchygreen(1, 2);
          glstrain(5) = cauchygreen(2, 0);
        }
        else if (kintype_ == INPAR::STR::kinem_nonlinearTotLag)
        {
          // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
          defgrd.MultiplyTT(xcurr, N_XYZ);

          // right Cauchy-Green tensor = F^T * F
          LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> cauchygreen;
          cauchygreen.MultiplyTN(defgrd, defgrd);

          glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
          glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
          glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
          glstrain(3) = cauchygreen(0, 1);
          glstrain(4) = cauchygreen(1, 2);
          glstrain(5) = cauchygreen(2, 0);
        }
        else if (kintype_ == INPAR::STR::kinem_linear)
        {
          // in kinematically linear analysis the deformation gradient is equal to identity
          // no difference between reference and current state
          for (int i = 0; i < 3; ++i) defgrd(i, i) = 1.0;

          // nodal displacement vector
          LINALG::Matrix<NUMDOF_SOH8, 1> nodaldisp;
          for (int i = 0; i < NUMDOF_SOH8; ++i) nodaldisp(i, 0) = mydisp[i];
          // compute linear B-operator
          LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8> bop;
          for (int i = 0; i < NUMNOD_SOH8; ++i)
          {
            bop(0, NODDOF_SOH8 * i + 0) = N_XYZ(0, i);
            bop(0, NODDOF_SOH8 * i + 1) = 0.0;
            bop(0, NODDOF_SOH8 * i + 2) = 0.0;
            bop(1, NODDOF_SOH8 * i + 0) = 0.0;
            bop(1, NODDOF_SOH8 * i + 1) = N_XYZ(1, i);
            bop(1, NODDOF_SOH8 * i + 2) = 0.0;
            bop(2, NODDOF_SOH8 * i + 0) = 0.0;
            bop(2, NODDOF_SOH8 * i + 1) = 0.0;
            bop(2, NODDOF_SOH8 * i + 2) = N_XYZ(2, i);

            bop(3, NODDOF_SOH8 * i + 0) = N_XYZ(1, i);
            bop(3, NODDOF_SOH8 * i + 1) = N_XYZ(0, i);
            bop(3, NODDOF_SOH8 * i + 2) = 0.0;
            bop(4, NODDOF_SOH8 * i + 0) = 0.0;
            bop(4, NODDOF_SOH8 * i + 1) = N_XYZ(2, i);
            bop(4, NODDOF_SOH8 * i + 2) = N_XYZ(1, i);
            bop(5, NODDOF_SOH8 * i + 0) = N_XYZ(2, i);
            bop(5, NODDOF_SOH8 * i + 1) = 0.0;
            bop(5, NODDOF_SOH8 * i + 2) = N_XYZ(0, i);
          }

          // compute linear strain at GP
          glstrain.Multiply(bop, nodaldisp);
        }
        else
          dserror("unknown kinematic type for energy calculation");

        // EAS technology: "enhance the strains"  ----------------------------- EAS
        if (eastype_ != soh8_easnone)
        {
          M.LightShape(MAT::NUM_STRESS_3D, neas_);
          // map local M to global, also enhancement is referred to element origin
          // M = detJ0/detJ T0^{-T} . M
          // Epetra_SerialDenseMatrix Mtemp(M); // temp M for Matrix-Matrix-Product
          // add enhanced strains = M . alpha to GL strains to "unlock" element
          switch (eastype_)
          {
            case DRT::ELEMENTS::So_hex8::soh8_easfull:
              LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
                  soh8_easfull>(M.A(), detJ0 / detJ_[gp], T0invT.A(), (M_GP->at(gp)).A());
              LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, soh8_easfull, 1>(
                  1.0, glstrain.A(), 1.0, M.A(), alpha->A());
              break;
            case DRT::ELEMENTS::So_hex8::soh8_easmild:
              LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
                  soh8_easmild>(M.A(), detJ0 / detJ_[gp], T0invT.A(), (M_GP->at(gp)).A());
              LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, soh8_easmild, 1>(
                  1.0, glstrain.A(), 1.0, M.A(), alpha->A());
              break;
            case DRT::ELEMENTS::So_hex8::soh8_eassosh8:
              LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
                  soh8_eassosh8>(M.A(), detJ0 / detJ_[gp], T0invT.A(), (M_GP->at(gp)).A());
              LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, soh8_eassosh8, 1>(
                  1.0, glstrain.A(), 1.0, M.A(), alpha->A());
              break;
            case DRT::ELEMENTS::So_hex8::soh8_easnone:
              break;
            default:
              dserror("Don't know what to do with EAS type %d", eastype_);
              break;
          }
        }  // ------------------------------------------------------------------ EAS

        if (defgrd.Determinant() <= 0.0)
        {
          if (IsParamsInterface() and StrParamsInterface().IsTolerateErrors())
          {
            StrParamsInterface().SetEleEvalErrorFlag(
                STR::ELEMENTS::ele_error_negative_det_of_def_gradient);
            return 0;
          }
          else
          {
            dserror("Negative deformation gradient!");
          }
        }

        // call material for evaluation of strain energy function
        double psi = 0.0;
        SolidMaterial()->StrainEnergy(glstrain, psi, Id());

        // sum up GP contribution to internal energy
        intenergy += fac * psi;
      }

      if (IsParamsInterface())  // new structural time integration
      {
        StrParamsInterface().AddContributionToEnergyType(intenergy, STR::internal_energy);
      }
      else  // old structural time integration
      {
        // check length of elevec1
        if (elevec1_epetra.Length() < 1) dserror("The given result vector is too short.");

        elevec1_epetra(0) = intenergy;
      }
    }
    break;
    //==================================================================================
    case ELEMENTS::struct_calc_errornorms:
    {
      // IMPORTANT NOTES (popp 10/2010):
      // - error norms are based on a small deformation assumption (linear elasticity)
      // - extension to finite deformations would be possible without difficulties,
      //   however analytical solutions are extremely rare in the nonlinear realm
      // - only implemented for purely displacement-based version, not yet for EAS
      // - only works for materials, which implement a hyperelastic strain energy
      //   function (currently this is the case for St.-Venant-Kirchhoff and NeoHooke)
      // - analytical solutions are currently stored in a repository in the CONTACT
      //   namespace, however they could (should?) be moved to a more general location

      // check length of elevec1
      if (elevec1_epetra.Length() < 3) dserror("The given result vector is too short.");

      // not yet implemented for EAS case
      if (eastype_ != soh8_easnone) dserror("Error norms not yet implemented for EAS.");

      // declaration of variables
      double l2norm = 0.0;
      double h1norm = 0.0;
      double energynorm = 0.0;

      // shape functions, derivatives and integration weights
      const static std::vector<LINALG::Matrix<NUMNOD_SOH8, 1>> vals = soh8_shapefcts();
      const static std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = soh8_derivs();
      const static std::vector<double> weights = soh8_weights();

      // get displacements and extract values of this element
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state displacement vector");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      // nodal displacement vector
      LINALG::Matrix<NUMDOF_SOH8, 1> nodaldisp;
      for (int i = 0; i < NUMDOF_SOH8; ++i) nodaldisp(i, 0) = mydisp[i];

      // reference geometry (nodal positions)
      LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe;
      DRT::Node** nodes = Nodes();
      for (int i = 0; i < NUMNOD_SOH8; ++i)
      {
        xrefe(i, 0) = nodes[i]->X()[0];
        xrefe(i, 1) = nodes[i]->X()[1];
        xrefe(i, 2) = nodes[i]->X()[2];
      }

      // deformation gradient = identity tensor (geometrically linear case!)
      LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrd(true);
      for (int i = 0; i < NUMDIM_SOH8; ++i) defgrd(i, i) = 1.0;

      //----------------------------------------------------------------
      // loop over all Gauss points
      //----------------------------------------------------------------
      for (unsigned gp = 0; gp < NUMGPT_SOH8; gp++)
      {
        // Gauss weights and Jacobian determinant
        double fac = detJ_[gp] * weights[gp];

        // Gauss point in reference configuration
        LINALG::Matrix<NUMDIM_SOH8, 1> xgp(true);
        for (int k = 0; k < NUMDIM_SOH8; ++k)
          for (int n = 0; n < NUMNOD_SOH8; ++n) xgp(k, 0) += (vals[gp])(n)*xrefe(n, k);

        //**************************************************************
        // get analytical solution
        LINALG::Matrix<NUMDIM_SOH8, 1> uanalyt(true);
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1> strainanalyt(true);
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> derivanalyt(true);

        CONTACT::AnalyticalSolutions3D(xgp, uanalyt, strainanalyt, derivanalyt);
        //**************************************************************

        //--------------------------------------------------------------
        // (1) L2 norm
        //--------------------------------------------------------------

        // compute displacements at GP
        LINALG::Matrix<NUMDIM_SOH8, 1> ugp(true);
        for (int k = 0; k < NUMDIM_SOH8; ++k)
          for (int n = 0; n < NUMNOD_SOH8; ++n)
            ugp(k, 0) += (vals[gp])(n)*nodaldisp(NODDOF_SOH8 * n + k, 0);

        // displacement error
        LINALG::Matrix<NUMDIM_SOH8, 1> uerror(true);
        for (int k = 0; k < NUMDIM_SOH8; ++k) uerror(k, 0) = uanalyt(k, 0) - ugp(k, 0);

        // compute GP contribution to L2 error norm
        l2norm += fac * uerror.Dot(uerror);

        //--------------------------------------------------------------
        // (2) H1 norm
        //--------------------------------------------------------------

        // compute derivatives N_XYZ at GP w.r.t. material coordinates
        // by N_XYZ = J^-1 * N_rst
        LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_XYZ(true);
        N_XYZ.Multiply(invJ_[gp], derivs[gp]);

        // compute partial derivatives at GP
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> derivgp(true);
        for (int l = 0; l < NUMDIM_SOH8; ++l)
          for (int m = 0; m < NUMDIM_SOH8; ++m)
            for (int k = 0; k < NUMNOD_SOH8; ++k)
              derivgp(l, m) += N_XYZ(m, k) * nodaldisp(NODDOF_SOH8 * k + l, 0);

        // derivative error
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> deriverror(true);
        for (int k = 0; k < NUMDIM_SOH8; ++k)
          for (int m = 0; m < NUMDIM_SOH8; ++m)
            deriverror(k, m) = derivanalyt(k, m) - derivgp(k, m);

        // compute GP contribution to H1 error norm
        h1norm += fac * deriverror.Dot(deriverror);
        h1norm += fac * uerror.Dot(uerror);

        //--------------------------------------------------------------
        // (3) Energy norm
        //--------------------------------------------------------------

        // compute linear B-operator
        LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8> bop;
        for (int i = 0; i < NUMNOD_SOH8; ++i)
        {
          bop(0, NODDOF_SOH8 * i + 0) = N_XYZ(0, i);
          bop(0, NODDOF_SOH8 * i + 1) = 0.0;
          bop(0, NODDOF_SOH8 * i + 2) = 0.0;
          bop(1, NODDOF_SOH8 * i + 0) = 0.0;
          bop(1, NODDOF_SOH8 * i + 1) = N_XYZ(1, i);
          bop(1, NODDOF_SOH8 * i + 2) = 0.0;
          bop(2, NODDOF_SOH8 * i + 0) = 0.0;
          bop(2, NODDOF_SOH8 * i + 1) = 0.0;
          bop(2, NODDOF_SOH8 * i + 2) = N_XYZ(2, i);

          bop(3, NODDOF_SOH8 * i + 0) = N_XYZ(1, i);
          bop(3, NODDOF_SOH8 * i + 1) = N_XYZ(0, i);
          bop(3, NODDOF_SOH8 * i + 2) = 0.0;
          bop(4, NODDOF_SOH8 * i + 0) = 0.0;
          bop(4, NODDOF_SOH8 * i + 1) = N_XYZ(2, i);
          bop(4, NODDOF_SOH8 * i + 2) = N_XYZ(1, i);
          bop(5, NODDOF_SOH8 * i + 0) = N_XYZ(2, i);
          bop(5, NODDOF_SOH8 * i + 1) = 0.0;
          bop(5, NODDOF_SOH8 * i + 2) = N_XYZ(0, i);
        }

        // compute linear strain at GP
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1> straingp(true);
        straingp.Multiply(bop, nodaldisp);

        // strain error
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1> strainerror(true);
        for (int k = 0; k < MAT::NUM_STRESS_3D; ++k)
          strainerror(k, 0) = strainanalyt(k, 0) - straingp(k, 0);

        // compute energy error
        double psierror = 0.0;
        SolidMaterial()->StrainEnergy(strainerror, psierror, Id());

        // compute GP contribution to energy error norm
        energynorm += fac * psierror;

        // std::cout << "UAnalytical:      " << uanalyt << std::endl;
        // std::cout << "UDiscrete:        " << ugp << std::endl;
        // std::cout << "StrainAnalytical: " << strainanalyt << std::endl;
        // std::cout << "StrainDiscrete:   " << straingp << std::endl;
        // std::cout << "DerivAnalytical:  " << derivanalyt << std::endl;
        // std::cout << "DerivDiscrete:    " << derivgp << std::endl;
        // std::cout << std::endl;
      }
      //----------------------------------------------------------------

      // return results
      elevec1_epetra(0) = l2norm;
      elevec1_epetra(1) = h1norm;
      elevec1_epetra(2) = energynorm;
    }
    break;
      //==================================================================================
    case ELEMENTS::multi_calc_dens:
    {
      soh8_homog(params);
    }
    break;
      //==================================================================================
    case ELEMENTS::struct_interpolate_velocity_to_point:
    {
      // get displacements and extract values of this element (set in PrepareFluidOp())
      Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velocity");

#ifdef DEBUG
      if (dispnp == Teuchos::null) dserror("Cannot get state displacement vector");
      if (velnp == Teuchos::null) dserror("Cannot get state velocity vector");
#endif

      std::vector<double> mydispnp(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp, mydispnp, lm);

      std::vector<double> myvelnp(lm.size());
      DRT::UTILS::ExtractMyValues(*velnp, myvelnp, lm);

      // update element geometry
      LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe;  // reference coord. of element
      LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xcurr;  // current coord. of element
      DRT::Node** nodes = Nodes();
      for (int i = 0; i < NUMNOD_SOH8; ++i)
      {
        const double* X = nodes[i]->X();
        xrefe(i, 0) = X[0];
        xrefe(i, 1) = X[1];
        xrefe(i, 2) = X[2];

        xcurr(i, 0) = xrefe(i, 0) + mydispnp[i * NODDOF_SOH8 + 0];
        xcurr(i, 1) = xrefe(i, 1) + mydispnp[i * NODDOF_SOH8 + 1];
        xcurr(i, 2) = xrefe(i, 2) + mydispnp[i * NODDOF_SOH8 + 2];
      }

      // shape functions and derivatives w.r.t. r,s,t
      LINALG::Matrix<NUMNOD_SOH8, 1> shapefcts;
      LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> deriv1;
      // coordinates of given point in reference coordinates
      LINALG::Matrix<NUMDIM_SOH8, 1> xsi;
      xsi(0) = elevec2_epetra(0);
      xsi(1) = elevec2_epetra(1);
      xsi(2) = elevec2_epetra(2);
      // evaluate shape functions and derivatives at given point w.r.t r,s,t
      DRT::UTILS::shape_function<DRT::Element::hex8>(xsi, shapefcts);
      DRT::UTILS::shape_function_deriv1<DRT::Element::hex8>(xsi, deriv1);

      LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> myvelocitynp;
      for (int node = 0; node < NUMNOD_SOH8; ++node)
      {
        for (int dim = 0; dim < NUMDIM_SOH8; ++dim)
        {
          myvelocitynp(node, dim) = myvelnp[node * 3 + dim];
        }
      }

      if (params.get("calculate_velocity", 1))
      {
        //************************************************************************
        // 1.) interpolation of velocity at n+1 to given point
        //************************************************************************

        // give back velocity at given point
        LINALG::Matrix<3, 1> result;
        result.MultiplyTN(myvelocitynp, shapefcts);
        for (int i = 0; i < 3; ++i) elevec1_epetra(i) = result(i, 0);
      }
      else
      {
        //************************************************************************
        // 2.) calculation of divergence of structural velocity
        //************************************************************************

        // get Jacobian matrix
        // actually compute its transpose....
        /*
                 +-            -+ -1
                 | dx   dy   dz |
                 | --   --   -- |
                 | dr   dr   dr |
                 |              |
                 | dx   dy   dz |
        J^{-T} = | --   --   -- |
                 | ds   ds   ds |
                 |              |
                 | dx   dy   dz |
                 | --   --   -- |
                 | dt   dt   dt |
                 +-            -+
        */

        // global first derivatives of shape functions w.r.t x,y,z (by derxy1 = J^-T * deriv1)
        LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> derxy1;
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> transJ;
        transJ.Multiply(deriv1, xrefe);
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> invJ(transJ);
        invJ.Invert();
        derxy1.Multiply(invJ, deriv1);

        // build (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * derxy1^T
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrdnp(true);
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrdnp_inv(true);
        defgrdnp.MultiplyTT(xcurr, derxy1);
        defgrdnp_inv.Invert(defgrdnp);

        // evaluate divergence of structural velocity
        double velocitydivergence = 0.0;

        // div(v) = derxy1_A,I (dX_I/dx_i) v_Ai
        for (int I = 0; I < NUMDIM_SOH8; ++I)
        {
          for (int A = 0; A < NUMNOD_SOH8; ++A)
          {
            for (int i = 0; i < NUMDIM_SOH8; ++i)
            {
              velocitydivergence += derxy1(I, A) * (defgrdnp_inv(I, i) * myvelocitynp(A, i));
            }
          }
        }

        // vector to fill
        elevec1_epetra(3) = velocitydivergence;
      }
    }
    break;
    //==================================================================================
    // in case of multi-scale problems, possible EAS internal data on microscale
    // have to be stored in every macroscopic Gauss point
    // allocation and initializiation of these data arrays can only be
    // done in the elements that know the number of EAS parameters
    case ELEMENTS::multi_init_eas:
    {
      soh8_eas_init_multi(params);
    }
    break;
    //==================================================================================
    // in case of multi-scale problems, possible EAS internal data on microscale
    // have to be stored in every macroscopic Gauss point
    // before any microscale simulation, EAS internal data has to be
    // set accordingly
    case ELEMENTS::multi_set_eas:
    {
      soh8_set_eas_multi(params);
    }
    break;
    //==================================================================================
    // read restart of microscale
    case ELEMENTS::multi_readrestart:
    {
      soh8_read_restart_multi();
    }
    break;
    //==================================================================================
    case ELEMENTS::struct_update_prestress:
    {
      time_ = params.get<double>("total time");
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get displacement state");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      // build def gradient for every gauss point
      LINALG::SerialDenseMatrix gpdefgrd(NUMGPT_SOH8, 9);
      DefGradient(mydisp, gpdefgrd, *prestress_);

      // update deformation gradient and put back to storage
      LINALG::Matrix<3, 3> deltaF;
      LINALG::Matrix<3, 3> Fhist;
      LINALG::Matrix<3, 3> Fnew;
      for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
      {
        prestress_->StoragetoMatrix(gp, deltaF, gpdefgrd);
        prestress_->StoragetoMatrix(gp, Fhist, prestress_->FHistory());
        Fnew.Multiply(deltaF, Fhist);
        prestress_->MatrixtoStorage(gp, Fnew, prestress_->FHistory());
        // if(gp ==1)
        // {
        // std::cout << "Fhist  " << Fhist << std::endl;
        // std::cout << "Fhnew  " << new << std::endl;
        // }
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
    //==================================================================================
    case ELEMENTS::inversedesign_update:
    {
      time_ = params.get<double>("total time");
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get displacement state");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      invdesign_->soh8_StoreMaterialConfiguration(this, mydisp);
      invdesign_->IsInit() = true;  // this is to make the restart work
    }
    break;
    //==================================================================================
    case ELEMENTS::inversedesign_switch:
    {
      time_ = params.get<double>("total time");
    }
    break;
    //==================================================================================
    // evaluate stresses and strains at gauss points and store gpstresses in map <EleId, gpstresses
    // >
    case ELEMENTS::struct_calc_global_gpstresses_map:
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
        Teuchos::RCP<std::vector<char>> plstraindata =
            params.get<Teuchos::RCP<std::vector<char>>>("plstrain", Teuchos::null);
        if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
        if (stressdata == Teuchos::null) dserror("Cannot get 'stress' data");
        if (straindata == Teuchos::null) dserror("Cannot get 'strain' data");
        if (plstraindata == Teuchos::null) dserror("Cannot get 'plastic strain' data");
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

        std::vector<double> mydispmat(lm.size(), 0.0);
        if (structale_)
        {
          Teuchos::RCP<const Epetra_Vector> dispmat =
              discretization.GetState("material_displacement");
          ;
          DRT::UTILS::ExtractMyValues(*dispmat, mydispmat, lm);
        }

        LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D> stress;
        LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D> strain;
        LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D> plstrain;
        INPAR::STR::StressType iostress =
            DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
        INPAR::STR::StrainType iostrain =
            DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);
        INPAR::STR::StrainType ioplstrain =
            DRT::INPUT::get<INPAR::STR::StrainType>(params, "ioplstrain", INPAR::STR::strain_none);

        if (pstype_ == INPAR::STR::prestress_id && time_ <= pstime_)  // inverse design analysis
          invdesign_->soh8_nlnstiffmass(this, lm, mydisp, myres, NULL, NULL, NULL, &stress, &strain,
              params, iostress, iostrain);

        else  // standard analysis
          nlnstiffmass(lm, mydisp, NULL, NULL, myres, mydispmat, NULL, NULL, NULL, NULL, NULL,
              &stress, &strain, &plstrain, params, iostress, iostrain, ioplstrain);

        // add stresses to global map
        // get EleID Id()
        int gid = Id();
        Teuchos::RCP<Epetra_SerialDenseMatrix> gpstress =
            Teuchos::rcp(new Epetra_SerialDenseMatrix);
        gpstress->Shape(NUMGPT_SOH8, MAT::NUM_STRESS_3D);

        // move stresses to serial dense matrix
        for (unsigned i = 0; i < NUMGPT_SOH8; i++)
        {
          for (int j = 0; j < MAT::NUM_STRESS_3D; j++)
          {
            (*gpstress)(i, j) = stress(i, j);
          }
        }

        // strains
        Teuchos::RCP<Epetra_SerialDenseMatrix> gpstrain =
            Teuchos::rcp(new Epetra_SerialDenseMatrix);
        gpstrain->Shape(NUMGPT_SOH8, MAT::NUM_STRESS_3D);

        // move stresses to serial dense matrix
        for (unsigned i = 0; i < NUMGPT_SOH8; i++)
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

        {
          DRT::PackBuffer data;
          AddtoPack(data, plstrain);
          data.StartPacking();
          AddtoPack(data, plstrain);
          std::copy(data().begin(), data().end(), std::back_inserter(*plstraindata));
        }
      }
    }
    break;
    //==================================================================================
    case ELEMENTS::struct_calc_predict:
    {
      // do nothing here
      break;
    }
    //==================================================================================
    // create a backup state for all internally stored variables (e.g. EAS)
    case ELEMENTS::struct_create_backup:
    {
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (res.is_null()) dserror("Cannot get state vector \"residual displacement\"");

      // extract the part for this element
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);

      soh8_create_eas_backup_state(myres);

      break;
    }
    //==================================================================================
    /* recover internally stored state variables from a previously created backup
     * state (e.g. EAS) */
    case ELEMENTS::struct_recover_from_backup:
    {
      soh8_recover_from_eas_backup_state();

      break;
    }
    default:
      dserror("Unknown type of action for So_hex8: %s", ELEMENTS::ActionType2String(act).c_str());
      break;
  }
  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)     maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  SetParamsInterfacePtr(params);
  // get values and switches from the condition
  const std::vector<int>* onoff = condition.Get<std::vector<int>>("onoff");
  const std::vector<double>* val = condition.Get<std::vector<double>>("val");

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  double time = -1.0;
  if (IsParamsInterface())
    time = ParamsInterface().GetTotalTime();
  else
    time = params.get("total time", -1.0);

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < NUMDIM_SOH8)
    dserror("Fewer functions or curves defined than the element has dofs.");

  for (int checkdof = NUMDIM_SOH8; checkdof < int(onoff->size()); ++checkdof)
  {
    if ((*onoff)[checkdof] != 0)
      dserror("Number of Dimensions in Neumann_Evalutaion is 3. Further DoFs are not considered.");
  }

  // (SPATIAL) FUNCTION BUSINESS
  const std::vector<int>* funct = condition.Get<std::vector<int>>("funct");
  LINALG::Matrix<NUMDIM_SOH8, 1> xrefegp(false);
  bool havefunct = false;
  if (funct)
    for (int dim = 0; dim < NUMDIM_SOH8; dim++)
      if ((*funct)[dim] > 0) havefunct = havefunct or true;


  /* ============================================================================*
  ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
  ** ============================================================================*/
  const static std::vector<LINALG::Matrix<NUMNOD_SOH8, 1>> shapefcts = soh8_shapefcts();
  const static std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = soh8_derivs();
  const static std::vector<double> gpweights = soh8_weights();
  /* ============================================================================*/

  // update element geometry
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe;  // material coord. of element
  DRT::Node** nodes = Nodes();
  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }
  /* ================================================= Loop over Gauss Points */
  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
  {
    // compute the Jacobian matrix
    LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> jac;
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
      for (int dim = 0; dim < NUMDIM_SOH8; dim++)
      {
        xrefegp(dim) = 0.0;
        for (int nodid = 0; nodid < NUMNOD_SOH8; ++nodid)
          xrefegp(dim) += shapefcts[gp](nodid) * xrefe(nodid, dim);
      }
    }

    // integration factor
    const double fac = gpweights[gp] * detJ;
    // distribute/add over element load vector
    for (int dim = 0; dim < NUMDIM_SOH8; dim++)
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
        for (int nodid = 0; nodid < NUMNOD_SOH8; ++nodid)
        {
          elevec1[nodid * NUMDIM_SOH8 + dim] += shapefcts[gp](nodid) * dim_fac;
        }
      }
    }

  } /* ==================================================== end of Loop over GP */

  return 0;
}  // DRT::ELEMENTS::So_hex8::EvaluateNeumann

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const double* DRT::ELEMENTS::So_hex8::soh8_get_coordinate_of_gausspoints(const unsigned dim) const
{
  static LINALG::Matrix<NUMGPT_SOH8, NUMDIM_SOH8> coordinates_of_gps(false);
  static bool init = false;

  if (not init)
  {
    if (gp_rule_.NumPoints() != NUMGPT_SOH8)
      dserror(
          "Inconsistent number of GPs: "
          "%d != %d",
          gp_rule_.NumPoints(), NUMGPT_SOH8);

    for (unsigned gp = 0; gp < gp_rule_.NumPoints(); ++gp)
      for (unsigned d = 0; d < NUMDIM_SOH8; ++d) coordinates_of_gps(gp, d) = gp_rule_.Point(gp)[d];
    // do it only once
    init = true;
  }

  return &coordinates_of_gps(0, dim);
}

/*----------------------------------------------------------------------*
 |  init the element jacobian mapping (protected)              gee 04/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::InitJacobianMapping()
{
  const static std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = soh8_derivs();
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe;
  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    Node** nodes = Nodes();
    if (!nodes) dserror("Nodes() returned null pointer");
    xrefe(i, 0) = Nodes()[i]->X()[0];
    xrefe(i, 1) = Nodes()[i]->X()[1];
    xrefe(i, 2) = Nodes()[i]->X()[2];
  }
  invJ_.resize(NUMGPT_SOH8);
  detJ_.resize(NUMGPT_SOH8);
  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
  {
    // invJ_[gp].Shape(NUMDIM_SOH8,NUMDIM_SOH8);
    invJ_[gp].Multiply(derivs[gp], xrefe);
    detJ_[gp] = invJ_[gp].Invert();
    if (detJ_[gp] <= 0.0) dserror("Element Jacobian mapping %10.5e <= 0.0", detJ_[gp]);

    if (pstype_ == INPAR::STR::prestress_mulf && pstime_ >= time_)
      if (!(prestress_->IsInit()))
        prestress_->MatrixtoStorage(gp, invJ_[gp], prestress_->JHistory());

    if (pstype_ == INPAR::STR::prestress_id && pstime_ < time_)
      if (!(invdesign_->IsInit()))
      {
        // printf("Ele %d id use InitJacobianMapping pstime < time %10.5e <
        // %10.5e\n",Id(),pstime_,time_);
        invdesign_->MatrixtoStorage(gp, invJ_[gp], invdesign_->JHistory());
        invdesign_->DetJHistory()[gp] = detJ_[gp];
      }
  }

  if (pstype_ == INPAR::STR::prestress_mulf && pstime_ >= time_) prestress_->IsInit() = true;

  if (pstype_ == INPAR::STR::prestress_id && pstime_ < time_) invdesign_->IsInit() = true;

  return;
}
/*----------------------------------------------------------------------*
 |  init the element jacobian mapping with respect to the    farah 06/13|
 |  material configuration.                                             |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8::InitJacobianMapping(std::vector<double>& dispmat)
{
  const static std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = soh8_derivs();
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xmat;

  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    Node** nodes = Nodes();
    if (!nodes) dserror("Nodes() returned null pointer");

    xmat(i, 0) = Nodes()[i]->X()[0] + dispmat[i * NODDOF_SOH8 + 0];
    xmat(i, 1) = Nodes()[i]->X()[1] + dispmat[i * NODDOF_SOH8 + 1];
    xmat(i, 2) = Nodes()[i]->X()[2] + dispmat[i * NODDOF_SOH8 + 2];
  }
  invJ_.clear();
  detJ_.clear();
  invJ_.resize(NUMGPT_SOH8);
  detJ_.resize(NUMGPT_SOH8);
  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
  {
    // invJ_[gp].Shape(NUMDIM_SOH8,NUMDIM_SOH8);
    invJ_[gp].Multiply(derivs[gp], xmat);
    detJ_[gp] = invJ_[gp].Invert();

    if (detJ_[gp] <= 0.0)
    {
      if (IsParamsInterface() and StrParamsInterface().IsTolerateErrors())
      {
        StrParamsInterface().SetEleEvalErrorFlag(
            STR::ELEMENTS::ele_error_negative_det_of_def_gradient);
        return 1;
      }
      else
        dserror("Element Jacobian mapping %10.5e <= 0.0", detJ_[gp]);
    }
  }

  return 0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double DRT::ELEMENTS::So_hex8::soh8_get_min_det_jac_at_corners(
    const LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8>& xcurr) const
{
  LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> xcurr_t(false);
  xcurr_t.UpdateT(xcurr);
  return DRT::UTILS::GetMinimalJacDeterminantAtNodes<DRT::Element::hex8>(xcurr_t);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_error_handling(const double& det_curr,
    Teuchos::ParameterList& params, const int line_id, const STR::ELEMENTS::EvalErrorFlag flag)
{
  ErrorHandling(det_curr, params, line_id, flag);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_computeEASInc(
    const std::vector<double>& residual, Epetra_SerialDenseMatrix* const eas_inc)
{
  Epetra_SerialDenseMatrix* oldKaainv = data_.GetMutable<Epetra_SerialDenseMatrix>("invKaa");
  Epetra_SerialDenseMatrix* oldKda = data_.GetMutable<Epetra_SerialDenseMatrix>("Kda");
  Epetra_SerialDenseMatrix* oldfeas = data_.GetMutable<Epetra_SerialDenseMatrix>("feas");
  if (!oldKaainv || !oldKda || !oldfeas) dserror("Missing EAS history data");

  // we need the (residual) displacement at the previous step
  LINALG::SerialDenseVector res_d_eas(NUMDOF_SOH8);
  for (int i = 0; i < NUMDOF_SOH8; ++i) res_d_eas(i) = residual[i];
  // --- EAS default update ---------------------------
  Epetra_SerialDenseMatrix eashelp(neas_, 1);
  /*----------- make multiplication eashelp = oldLt * disp_incr[kstep] */
  oldKda->Multiply(false, res_d_eas, eashelp);
  /*---------------------------------------- add old Rtilde to eashelp */
  eashelp += *oldfeas;
  /*--------- make multiplication alpha_inc = - old Dtildinv * eashelp */
  oldKaainv->Multiply(false, eashelp, *eas_inc);
  eas_inc->Scale(-1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_recover(
    const std::vector<int>& lm, const std::vector<double>& residual)
{
  // for eas
  Epetra_SerialDenseMatrix* alpha = NULL;
  Epetra_SerialDenseMatrix* eas_inc = NULL;
  // get access to the interface parameters
  const double step_length = StrParamsInterface().GetStepLength();
  const bool iseas = (eastype_ != soh8_easnone);

  // have eas?
  if (iseas)
  {
    // access general eas history stuff stored in element
    // get alpha of previous iteration
    alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");
    // get the old eas increment
    eas_inc = data_.GetMutable<Epetra_SerialDenseMatrix>("eas_inc");
    if (!alpha || !eas_inc) dserror("Missing EAS history data (eas_inc and/or alpha)");
  }

  /* if it is a default step, we have to recover the condensed
   * solution vectors */
  if (StrParamsInterface().IsDefaultStep())
  {
    /* recovery of the enhanced assumed strain increment and
     * update of the eas dofs. */
    if (iseas)
    {
      // first, store the eas state of the previous accepted Newton step
      StrParamsInterface().SumIntoMyPreviousSolNorm(
          NOX::NLN::StatusTest::quantity_eas, neas_, (*alpha)[0], Owner());

      // compute the eas increment
      soh8_computeEASInc(residual, eas_inc);

      /*--------------------------- update alpha += step_length * alfa_inc */
      for (int i = 0; i < neas_; ++i) (*alpha)(i, 0) += step_length * (*eas_inc)(i, 0);
    }  // if (iseas)
  }    // if (*isdefault_step_ptr_)
  /* if it is no default step, we can correct the update and the current eas
   * state without the need for any matrix-vector products. */
  else
  {
    // The first step has to be a default step!
    if (old_step_length_ < 0.0) dserror("The old step length was not defined!");
    /* if this is no full step, we have to adjust the length of the
     * enhanced assumed strain incremental step. */
    if (iseas)
    {
      /* undo the previous step:
       *            alpha_new = alpha_old - old_step * alpha_inc
       * and update the solution variable with the new step length:
       *            alpha_new = alpha_new + new_step * alpha_inc */
      for (int i = 0; i < neas_; ++i)
        (*alpha)(i, 0) += (step_length - old_step_length_) * (*eas_inc)(i, 0);

      //      {
      //        std::cout << "EAS #" << Id() << ":\n";
      //        alpha->Print( std::cout );
      //        eas_inc->Print( std::cout );
      //        std::cout << "\n";
      //      }
    }  // if (nhyb_)
  }    // else

  // save the old step length
  old_step_length_ = step_length;

  // Check if the eas incr is tested and if yes, calculate the element
  // contribution to the norm
  if (iseas)
    StrParamsInterface().SumIntoMyUpdateNorm(NOX::NLN::StatusTest::quantity_eas, neas_,
        (*eas_inc)[0], (*alpha)[0], step_length, Owner());

  // the element internal stuff should be up-to-date for now...
  return;
}



/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::nlnstiffmass(std::vector<int>& lm,  // location matrix
    std::vector<double>& disp,                                   // current displacements
    std::vector<double>* vel,                                    // current velocities
    std::vector<double>* acc,                                    // current accelerations
    std::vector<double>& residual,                               // current residual displ
    std::vector<double>& dispmat,                                // current material displacements
    LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* stiffmatrix,       // element stiffness matrix
    LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* massmatrix,        // element mass matrix
    LINALG::Matrix<NUMDOF_SOH8, 1>* force,                       // element internal force vector
    LINALG::Matrix<NUMDOF_SOH8, 1>* forceinert,                  // element inertial force vector
    LINALG::Matrix<NUMDOF_SOH8, 1>* force_str,                   // element structural force vector
    LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestress,  // stresses at GP
    LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestrain,  // strains at GP
    LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* eleplstrain,  // plastic strains at GP
    Teuchos::ParameterList& params,           // algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,    // stress output option
    const INPAR::STR::StrainType iostrain,    // strain output option
    const INPAR::STR::StrainType ioplstrain)  // plastic strain output option
{
  /* ============================================================================*
  ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
  ** ============================================================================*/
  const static std::vector<LINALG::Matrix<NUMNOD_SOH8, 1>> shapefcts = soh8_shapefcts();
  const static std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = soh8_derivs();
  const static std::vector<double> gpweights = soh8_weights();
  /* ============================================================================*/

  // check for prestressing
  if (pstype_ != INPAR::STR::prestress_none && eastype_ != soh8_easnone)
    dserror("No way you can do mulf or id prestressing with EAS turned on!");

  // update element geometry
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe;  // reference coord. of element
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xcurr;  // current  coord. of element
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xdisp;

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

    if (pstype_ == INPAR::STR::prestress_mulf)
    {
      xdisp(i, 0) = disp[i * NODDOF_SOH8 + 0];
      xdisp(i, 1) = disp[i * NODDOF_SOH8 + 1];
      xdisp(i, 2) = disp[i * NODDOF_SOH8 + 2];
    }
  }

  // safety check before the actual evaluation starts
  const double min_detJ_curr = soh8_get_min_det_jac_at_corners(xcurr);
  if (min_detJ_curr <= 0.0)
  {
    soh8_error_handling(
        min_detJ_curr, params, __LINE__, STR::ELEMENTS::ele_error_determinant_at_corner);
    return;
  }

  double elediagonallength = 0.0;
  if (!analyticalmaterialtangent_)
    elediagonallength = sqrt(pow(xrefe(0, 0) - xrefe(6, 0), 2) + pow(xrefe(0, 1) - xrefe(6, 1), 2) +
                             pow(xrefe(0, 2) - xrefe(6, 2), 2));

  // we need the (residual) displacement at the previous step
  LINALG::Matrix<NUMDOF_SOH8, 1> nodaldisp;
  for (int i = 0; i < NUMDOF_SOH8; ++i) nodaldisp(i, 0) = disp[i];

  /*
  ** EAS Technology: declare, intialize, set up, and alpha history -------- EAS
  */
  // in any case declare variables, sizes etc. only in eascase
  Epetra_SerialDenseMatrix* alpha = NULL;              // EAS alphas
  std::vector<Epetra_SerialDenseMatrix>* M_GP = NULL;  // EAS matrix M at all GPs
  LINALG::SerialDenseMatrix M;                         // EAS matrix M at current GP
  Epetra_SerialDenseVector feas;                       // EAS portion of internal forces
  Epetra_SerialDenseMatrix Kaa;                        // EAS matrix Kaa
  Epetra_SerialDenseMatrix Kda;                        // EAS matrix Kda
  double detJ0;                                        // detJ(origin)
  Epetra_SerialDenseMatrix* oldfeas = NULL;            // EAS history
  Epetra_SerialDenseMatrix* oldKaainv = NULL;          // EAS history
  Epetra_SerialDenseMatrix* oldKda = NULL;             // EAS history
  Epetra_SerialDenseMatrix* eas_inc = NULL;            // EAS increment

  // transformation matrix T0, maps M-matrix evaluated at origin
  // between local element coords and global coords
  // here we already get the inverse transposed T0
  LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> T0invT;  // trafo matrix

  if (eastype_ != soh8_easnone)
  {
    /*
    ** EAS Update of alphas:
    ** the current alphas are (re-)evaluated out of
    ** Kaa and Kda of previous step to avoid additional element call.
    ** This corresponds to the (innermost) element update loop
    ** in the nonlinear FE-Skript page 120 (load-control alg. with EAS)
    */
    alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // get alpha of previous iteration
    oldfeas = data_.GetMutable<Epetra_SerialDenseMatrix>("feas");
    oldKaainv = data_.GetMutable<Epetra_SerialDenseMatrix>("invKaa");
    oldKda = data_.GetMutable<Epetra_SerialDenseMatrix>("Kda");
    eas_inc = data_.GetMutable<Epetra_SerialDenseMatrix>("eas_inc");

    if (!alpha || !oldKaainv || !oldKda || !oldfeas || !eas_inc)
      dserror("Missing EAS history-data");

    // ============================== DEPRECATED ==============================
    // FixMe deprecated implementation
    if (not IsParamsInterface())
    {
      // we need the (residual) displacement at the previous step
      LINALG::SerialDenseVector res_d_eas(NUMDOF_SOH8);
      for (int i = 0; i < NUMDOF_SOH8; ++i) res_d_eas(i) = residual[i];

      // this is a line search step, i.e. the direction of the eas increments
      // has been calculated by a Newton step and now it is only scaled
      if (params.isParameter("alpha_ls"))
      {
        double alpha_ls = params.get<double>("alpha_ls");
        // undo step
        eas_inc->Scale(-1.);
        alpha->operator+=(*eas_inc);
        // scale increment
        eas_inc->Scale(-1. * alpha_ls);
        // add reduced increment
        alpha->operator+=(*eas_inc);
      }
      // add Kda . res_d to feas
      // new alpha is: - Kaa^-1 . (feas + Kda . old_d), here: - Kaa^-1 . feas
      else
      {
        switch (eastype_)
        {
          case DRT::ELEMENTS::So_hex8::soh8_easfull:
            LINALG::DENSEFUNCTIONS::multiply<double, soh8_easfull, NUMDOF_SOH8, 1>(
                1.0, *oldfeas, 1.0, *oldKda, res_d_eas);
            LINALG::DENSEFUNCTIONS::multiply<double, soh8_easfull, soh8_easfull, 1>(
                0.0, *eas_inc, -1.0, *oldKaainv, *oldfeas);
            LINALG::DENSEFUNCTIONS::update<double, soh8_easfull, 1>(1., *alpha, 1., *eas_inc);
            break;
          case DRT::ELEMENTS::So_hex8::soh8_easmild:
            LINALG::DENSEFUNCTIONS::multiply<double, soh8_easmild, NUMDOF_SOH8, 1>(
                1.0, *oldfeas, 1.0, *oldKda, res_d_eas);
            LINALG::DENSEFUNCTIONS::multiply<double, soh8_easmild, soh8_easmild, 1>(
                0.0, *eas_inc, -1.0, *oldKaainv, *oldfeas);
            LINALG::DENSEFUNCTIONS::update<double, soh8_easmild, 1>(1., *alpha, 1., *eas_inc);
            break;
          case DRT::ELEMENTS::So_hex8::soh8_eassosh8:
            LINALG::DENSEFUNCTIONS::multiply<double, soh8_eassosh8, NUMDOF_SOH8, 1>(
                1.0, *oldfeas, 1.0, *oldKda, res_d_eas);
            LINALG::DENSEFUNCTIONS::multiply<double, soh8_eassosh8, soh8_eassosh8, 1>(
                0.0, *eas_inc, -1.0, *oldKaainv, *oldfeas);
            LINALG::DENSEFUNCTIONS::update<double, soh8_eassosh8, 1>(1., *alpha, 1., *eas_inc);
            break;
          case DRT::ELEMENTS::So_hex8::soh8_easnone:
            break;
          default:
            dserror("Don't know what to do with EAS type %d", eastype_);
            break;
        }
      }
    }  // if (not IsInterface())
    // ============================== DEPRECATED ==============================

    /* end of EAS Update ******************/

    // EAS portion of internal forces, also called enhacement vector s or Rtilde
    feas.Size(neas_);

    // EAS matrix K_{alpha alpha}, also called Dtilde
    Kaa.Shape(neas_, neas_);

    // EAS matrix K_{d alpha}
    Kda.Shape(neas_, NUMDOF_SOH8);


    /* evaluation of EAS variables (which are constant for the following):
    ** -> M defining interpolation of enhanced strains alpha, evaluated at GPs
    ** -> determinant of Jacobi matrix at element origin (r=s=t=0.0)
    ** -> T0^{-T}
    */
    soh8_eassetup(&M_GP, detJ0, T0invT, xrefe);
  }  // -------------------------------------------------------------------- EAS

  // build new jacobian mapping with respect to the material configuration
  if (structale_ == true) InitJacobianMapping(dispmat);

  // check if we need to split the residuals (for Newton line search)
  // if true an additional global vector is assembled containing
  // the internal forces without the condensed EAS entries and the norm
  // of the EAS residual is calculated
  bool split_res = params.isParameter("cond_rhs_norm");

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_XYZ;

  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrd(true);
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

    if (pstype_ == INPAR::STR::prestress_mulf)
    {
      // get Jacobian mapping wrt to the stored configuration
      LINALG::Matrix<3, 3> invJdef;
      prestress_->StoragetoMatrix(gp, invJdef, prestress_->JHistory());
      // get derivatives wrt to last spatial configuration
      LINALG::Matrix<3, 8> N_xyz;
      N_xyz.Multiply(invJdef, derivs[gp]);

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
    else if (kintype_ == INPAR::STR::kinem_nonlinearTotLag)
    {
      // standard kinematically nonlinear analysis
      defgrd.MultiplyTT(xcurr, N_XYZ);
    }
    else
    {
      // in kinematically linear analysis the deformation gradient is equal to identity
      // no difference between reference and current state
      for (int i = 0; i < 3; ++i) defgrd(i, i) = 1.0;
    }

    if (IsParamsInterface())
    {
      double det_defgrd = defgrd.Determinant();
      if (det_defgrd < 0.0)
      {
        if (StrParamsInterface().IsTolerateErrors())
        {
          StrParamsInterface().SetEleEvalErrorFlag(
              STR::ELEMENTS::ele_error_negative_det_of_def_gradient);
          stiffmatrix->Clear();
          force->Clear();
          return;
        }
        else
          dserror("negative deformation gradient determinant");
      }  // if (det_defgrd<0.0)
    }

    if (pstype_ == INPAR::STR::prestress_id && pstime_ < time_)
    {
      // printf("Ele %d entering id poststress\n",Id());
      // make the multiplicative update so that defgrd refers to
      // the reference configuration that resulted from the inverse
      // design analysis
      LINALG::Matrix<3, 3> Fhist;
      invdesign_->StoragetoMatrix(gp, Fhist, invdesign_->FHistory());
      LINALG::Matrix<3, 3> tmp3x3;
      tmp3x3.Multiply(defgrd, Fhist);
      defgrd = tmp3x3;

      // make detJ and invJ refer to the ref. configuration that resulted from
      // the inverse design analysis
      detJ = invdesign_->DetJHistory()[gp];
      invdesign_->StoragetoMatrix(gp, tmp3x3, invdesign_->JHistory());
      N_XYZ.Multiply(tmp3x3, derivs[gp]);
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
    LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8> bop;
    for (int i = 0; i < NUMNOD_SOH8; ++i)
    {
      bop(0, NODDOF_SOH8 * i + 0) = defgrd(0, 0) * N_XYZ(0, i);
      bop(0, NODDOF_SOH8 * i + 1) = defgrd(1, 0) * N_XYZ(0, i);
      bop(0, NODDOF_SOH8 * i + 2) = defgrd(2, 0) * N_XYZ(0, i);
      bop(1, NODDOF_SOH8 * i + 0) = defgrd(0, 1) * N_XYZ(1, i);
      bop(1, NODDOF_SOH8 * i + 1) = defgrd(1, 1) * N_XYZ(1, i);
      bop(1, NODDOF_SOH8 * i + 2) = defgrd(2, 1) * N_XYZ(1, i);
      bop(2, NODDOF_SOH8 * i + 0) = defgrd(0, 2) * N_XYZ(2, i);
      bop(2, NODDOF_SOH8 * i + 1) = defgrd(1, 2) * N_XYZ(2, i);
      bop(2, NODDOF_SOH8 * i + 2) = defgrd(2, 2) * N_XYZ(2, i);
      /* ~~~ */
      bop(3, NODDOF_SOH8 * i + 0) = defgrd(0, 0) * N_XYZ(1, i) + defgrd(0, 1) * N_XYZ(0, i);
      bop(3, NODDOF_SOH8 * i + 1) = defgrd(1, 0) * N_XYZ(1, i) + defgrd(1, 1) * N_XYZ(0, i);
      bop(3, NODDOF_SOH8 * i + 2) = defgrd(2, 0) * N_XYZ(1, i) + defgrd(2, 1) * N_XYZ(0, i);
      bop(4, NODDOF_SOH8 * i + 0) = defgrd(0, 1) * N_XYZ(2, i) + defgrd(0, 2) * N_XYZ(1, i);
      bop(4, NODDOF_SOH8 * i + 1) = defgrd(1, 1) * N_XYZ(2, i) + defgrd(1, 2) * N_XYZ(1, i);
      bop(4, NODDOF_SOH8 * i + 2) = defgrd(2, 1) * N_XYZ(2, i) + defgrd(2, 2) * N_XYZ(1, i);
      bop(5, NODDOF_SOH8 * i + 0) = defgrd(0, 2) * N_XYZ(0, i) + defgrd(0, 0) * N_XYZ(2, i);
      bop(5, NODDOF_SOH8 * i + 1) = defgrd(1, 2) * N_XYZ(0, i) + defgrd(1, 0) * N_XYZ(2, i);
      bop(5, NODDOF_SOH8 * i + 2) = defgrd(2, 2) * N_XYZ(0, i) + defgrd(2, 0) * N_XYZ(2, i);
    }

    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain_epetra(MAT::NUM_STRESS_3D);
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(glstrain_epetra.A(), true);
    if (kintype_ == INPAR::STR::kinem_nonlinearTotLag)
    {
      // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
      glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
      glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
      glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
      glstrain(3) = cauchygreen(0, 1);
      glstrain(4) = cauchygreen(1, 2);
      glstrain(5) = cauchygreen(2, 0);
    }
    else
    {
      // build the linearised strain epsilon = B_L . d
      glstrain.Multiply(bop, nodaldisp);
    }

    // deformation gradient consistent with (potentially EAS-modified) GL strains
    // without eas this is equal to the regular defgrd.
    LINALG::Matrix<3, 3> defgrd_mod(defgrd);

    // EAS technology: "enhance the strains"  ----------------------------- EAS
    if (eastype_ != soh8_easnone)
    {
      M.LightShape(MAT::NUM_STRESS_3D, neas_);
      // map local M to global, also enhancement is refered to element origin
      // M = detJ0/detJ T0^{-T} . M
      // Epetra_SerialDenseMatrix Mtemp(M); // temp M for Matrix-Matrix-Product
      // add enhanced strains = M . alpha to GL strains to "unlock" element
      switch (eastype_)
      {
        case DRT::ELEMENTS::So_hex8::soh8_easfull:
          LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
              soh8_easfull>(M.A(), detJ0 / detJ, T0invT.A(), (M_GP->at(gp)).A());
          LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, soh8_easfull, 1>(
              1.0, glstrain.A(), 1.0, M.A(), alpha->A());
          break;
        case DRT::ELEMENTS::So_hex8::soh8_easmild:
          LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
              soh8_easmild>(M.A(), detJ0 / detJ, T0invT.A(), (M_GP->at(gp)).A());
          LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, soh8_easmild, 1>(
              1.0, glstrain.A(), 1.0, M.A(), alpha->A());
          break;
        case DRT::ELEMENTS::So_hex8::soh8_eassosh8:
          LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
              soh8_eassosh8>(M.A(), detJ0 / detJ, T0invT.A(), (M_GP->at(gp)).A());
          LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, soh8_eassosh8, 1>(
              1.0, glstrain.A(), 1.0, M.A(), alpha->A());
          break;
        case DRT::ELEMENTS::So_hex8::soh8_easnone:
          break;
        default:
          dserror("Don't know what to do with EAS type %d", eastype_);
          break;
      }

      // calculate deformation gradient consistent with modified GL strain tensor
      if (Teuchos::rcp_static_cast<MAT::So3Material>(Material())->NeedsDefgrd())
        CalcConsistentDefgrd(defgrd, glstrain, defgrd_mod);

      const double det_defgrd_mod = defgrd_mod.Determinant();
      if (det_defgrd_mod <= 0.0)
      {
        soh8_error_handling(det_defgrd_mod, params, __LINE__,
            STR::ELEMENTS::ele_error_negative_det_of_def_gradient);
        return;
      }
    }  // ------------------------------------------------------------------ EAS

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
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> gl;
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
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> invdefgrd;
        invdefgrd.Invert(defgrd);

        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> temp;
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> euler_almansi;
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
      case INPAR::STR::strain_log:
      {
        if (elestrain == NULL) dserror("strain data not available");

        /// the Eularian logarithmic strain is defined as the natural logarithm of the left stretch
        /// tensor [1,2]: \f[
        ///    e_{log} = e_{hencky} = ln (\mathbf{V}) = \sum_{i=1}^3 (ln \lambda_i) \mathbf{n}_i
        ///    \otimes \mathbf{n}_i
        /// \f]
        ///< h3>References</h3>
        /// <ul>
        /// <li> [1] H. Xiao, Beijing, China, O. T. Bruhns and A. Meyers (1997) Logarithmic strain,
        /// logarithmic spin and logarithmic rate, Eq. 5 <li> [2] Caminero et al. (2011) Modeling
        /// large strain anisotropic elasto-plasticity with logarithmic strain and stress measures,
        /// Eq. 70
        /// </ul>
        ///
        /// \author HdV
        /// \date 08/13

        // eigenvalue decomposition (from elasthyper.cpp)
        LINALG::Matrix<3, 3> prstr2(true);  // squared principal stretches
        LINALG::Matrix<3, 1> prstr(true);   // principal stretch
        LINALG::Matrix<3, 3> prdir(true);   // principal directions
        LINALG::SYEV(cauchygreen, prstr2, prdir);

        // THE principal stretches
        for (int al = 0; al < 3; ++al) prstr(al) = std::sqrt(prstr2(al, al));

        // populating the logarithmic strain matrix
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> lnv(true);

        // checking if cauchy green is correctly determined to ensure eigen vectors in correct
        // direction i.e. a flipped eigenvector is also a valid solution C = \sum_{i=1}^3
        // (\lambda_i^2) \mathbf{n}_i \otimes \mathbf{n}_i
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> tempCG(true);

        for (int k = 0; k < 3; ++k)
        {
          double n_00, n_01, n_02, n_11, n_12, n_22 = 0.0;

          n_00 = prdir(0, k) * prdir(0, k);
          n_01 = prdir(0, k) * prdir(1, k);
          n_02 = prdir(0, k) * prdir(2, k);
          n_11 = prdir(1, k) * prdir(1, k);
          n_12 = prdir(1, k) * prdir(2, k);
          n_22 = prdir(2, k) * prdir(2, k);

          // only compute the symmetric components from a single eigenvector,
          // because eigenvalue directions are not consistent (it can be flipped)
          tempCG(0, 0) += (prstr(k)) * (prstr(k)) * n_00;
          tempCG(0, 1) += (prstr(k)) * (prstr(k)) * n_01;
          tempCG(0, 2) += (prstr(k)) * (prstr(k)) * n_02;
          tempCG(1, 0) += (prstr(k)) * (prstr(k)) * n_01;  // symmetry
          tempCG(1, 1) += (prstr(k)) * (prstr(k)) * n_11;
          tempCG(1, 2) += (prstr(k)) * (prstr(k)) * n_12;
          tempCG(2, 0) += (prstr(k)) * (prstr(k)) * n_02;  // symmetry
          tempCG(2, 1) += (prstr(k)) * (prstr(k)) * n_12;  // symmetry
          tempCG(2, 2) += (prstr(k)) * (prstr(k)) * n_22;

          // Computation of the Logarithmic strain tensor

          lnv(0, 0) += (std::log(prstr(k))) * n_00;
          lnv(0, 1) += (std::log(prstr(k))) * n_01;
          lnv(0, 2) += (std::log(prstr(k))) * n_02;
          lnv(1, 0) += (std::log(prstr(k))) * n_01;  // symmetry
          lnv(1, 1) += (std::log(prstr(k))) * n_11;
          lnv(1, 2) += (std::log(prstr(k))) * n_12;
          lnv(2, 0) += (std::log(prstr(k))) * n_02;  // symmetry
          lnv(2, 1) += (std::log(prstr(k))) * n_12;  // symmetry
          lnv(2, 2) += (std::log(prstr(k))) * n_22;
        }

        // compare CG computed with deformation gradient with CG computed
        // with eigenvalues and -vectors to determine/ensure the correct
        // orientation of the eigen vectors
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> diffCG(true);

        for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 3; ++j)
          {
            diffCG(i, j) = cauchygreen(i, j) - tempCG(i, j);
            // the solution to this problem is to evaluate the cauchygreen tensor with tempCG
            // computed with every combination of eigenvector orientations -- up to nine comparisons
            if (diffCG(i, j) > 1e-10)
              dserror(
                  "eigenvector orientation error with the diffCG giving problems: %10.5e \n BUILD "
                  "SOLUTION TO FIX IT",
                  diffCG(i, j));
          }
        }

        (*elestrain)(gp, 0) = lnv(0, 0);
        (*elestrain)(gp, 1) = lnv(1, 1);
        (*elestrain)(gp, 2) = lnv(2, 2);
        (*elestrain)(gp, 3) = lnv(0, 1);
        (*elestrain)(gp, 4) = lnv(1, 2);
        (*elestrain)(gp, 5) = lnv(0, 2);
      }
      break;
      case INPAR::STR::strain_none:
        break;
      default:
        dserror("requested strain type not available");
        break;
    }


    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix must be retrieved,
    ** all necessary data must be passed.
    */
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> cmat(true);
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> stress(true);

    UTILS::GetTemperatureForStructuralMaterial<hex8>(shapefcts[gp], params);

    if (Material()->MaterialType() == INPAR::MAT::m_constraintmixture ||
        Material()->MaterialType() == INPAR::MAT::m_growthremodel_elasthyper ||
        Material()->MaterialType() == INPAR::MAT::m_mixture_elasthyper)
    {
      LINALG::Matrix<1, NUMDIM_SOH8> point(true);
      soh8_GaussPointRefeCoords(point, xrefe, gp);
      params.set("gprefecoord", point);

      // center of element in reference configuration
      point.Clear();
      soh8_ElementCenterRefeCoords(point, xrefe);
      params.set("elecenter", point);
    }

    params.set<int>("gp", gp);

    // if output is requested only active stresses are written.
    params.set<int>("iostress", iostress);

    Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_static_cast<MAT::So3Material>(Material());
    so3mat->Evaluate(&defgrd_mod, &glstrain, params, &stress, &cmat, Id());

    // stop if the material evaluation fails
    if (IsParamsInterface() and StrParamsInterface().IsTolerateErrors())
      if (StrParamsInterface().GetEleEvalErrorFlag() != STR::ELEMENTS::ele_error_none) return;

    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp plastic strains (only in case of plastic strain output)
    switch (ioplstrain)
    {
      case INPAR::STR::strain_gl:
      {
        if (eleplstrain == NULL) dserror("plastic strain data not available");
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1> plglstrain =
            params.get<LINALG::Matrix<MAT::NUM_STRESS_3D, 1>>("plglstrain");
        for (int i = 0; i < 3; ++i) (*eleplstrain)(gp, i) = plglstrain(i);
        for (int i = 3; i < 6; ++i) (*eleplstrain)(gp, i) = 0.5 * plglstrain(i);
        break;
      }
      case INPAR::STR::strain_ea:
      {
        if (eleplstrain == NULL) dserror("plastic strain data not available");
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1> plglstrain =
            params.get<LINALG::Matrix<MAT::NUM_STRESS_3D, 1>>("plglstrain");
        // rewriting Green-Lagrange strains in matrix format
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> gl;
        gl(0, 0) = plglstrain(0);
        gl(0, 1) = 0.5 * plglstrain(3);
        gl(0, 2) = 0.5 * plglstrain(5);
        gl(1, 0) = gl(0, 1);
        gl(1, 1) = plglstrain(1);
        gl(1, 2) = 0.5 * plglstrain(4);
        gl(2, 0) = gl(0, 2);
        gl(2, 1) = gl(1, 2);
        gl(2, 2) = plglstrain(2);

        // inverse of deformation gradient
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> invdefgrd;
        invdefgrd.Invert(defgrd);

        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> temp;
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> euler_almansi;
        temp.Multiply(gl, invdefgrd);
        euler_almansi.MultiplyTN(invdefgrd, temp);

        (*eleplstrain)(gp, 0) = euler_almansi(0, 0);
        (*eleplstrain)(gp, 1) = euler_almansi(1, 1);
        (*eleplstrain)(gp, 2) = euler_almansi(2, 2);
        (*eleplstrain)(gp, 3) = euler_almansi(0, 1);
        (*eleplstrain)(gp, 4) = euler_almansi(1, 2);
        (*eleplstrain)(gp, 5) = euler_almansi(0, 2);
        break;
      }
      case INPAR::STR::strain_none:
        break;
      default:
        dserror("requested plastic strain type not available");
        break;
    }

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
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> cauchystress(false);
        PK2toCauchy(&stress, &defgrd, &cauchystress);

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

    const double detJ_w = detJ * gpweights[gp];
    // update internal force vector
    if (force != NULL)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w, bop, stress, 1.0);
    }

    // structural force vector
    if (split_res && force_str != NULL) force_str->MultiplyTN(detJ_w, bop, stress, 1.);

    // update stiffness matrix
    if (stiffmatrix != NULL)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<6, NUMDOF_SOH8> cb;
      cb.Multiply(cmat, bop);

      if (analyticalmaterialtangent_)
        stiffmatrix->MultiplyTN(detJ_w, bop, cb, 1.0);  // standard hex8 evaluation
      else
      {
        EvaluateFiniteDifferenceMaterialTangent(stiffmatrix, stress, disp, detJ_w, detJ, detJ0,
            elediagonallength, bop, cb, N_XYZ, T0invT, M_GP, alpha, M, gp, params);
      }


      if (kintype_ == INPAR::STR::kinem_nonlinearTotLag)
      {
        // integrate `geometric' stiffness matrix and add to keu *****************
        LINALG::Matrix<6, 1> sfac(stress);  // auxiliary integrated stress
        sfac.Scale(detJ_w);                 // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
        std::vector<double> SmB_L(3);       // intermediate Sm.B_L
        // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
        for (int inod = 0; inod < NUMNOD_SOH8; ++inod)
        {
          SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod) + sfac(5) * N_XYZ(2, inod);
          SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod) + sfac(4) * N_XYZ(2, inod);
          SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod) + sfac(2) * N_XYZ(2, inod);
          for (int jnod = 0; jnod < NUMNOD_SOH8; ++jnod)
          {
            double bopstrbop = 0.0;  // intermediate value
            for (int idim = 0; idim < NUMDIM_SOH8; ++idim)
              bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
            (*stiffmatrix)(3 * inod + 0, 3 * jnod + 0) += bopstrbop;
            (*stiffmatrix)(3 * inod + 1, 3 * jnod + 1) += bopstrbop;
            (*stiffmatrix)(3 * inod + 2, 3 * jnod + 2) += bopstrbop;
          }

        }  // end of integrate `geometric' stiffness******************************
      }

      // EAS technology: integrate matrices --------------------------------- EAS
      if (eastype_ != soh8_easnone)
      {
        // integrate Kaa: Kaa += (M^T . cmat . M) * detJ * w(gp)
        // integrate Kda: Kda += (M^T . cmat . B) * detJ * w(gp)
        // integrate feas: feas += (M^T . sigma) * detJ *wp(gp)
        LINALG::SerialDenseMatrix cM(MAT::NUM_STRESS_3D, neas_);  // temporary c . M
        switch (eastype_)
        {
          case DRT::ELEMENTS::So_hex8::soh8_easfull:
            LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
                soh8_easfull>(cM.A(), cmat.A(), M.A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_easfull, MAT::NUM_STRESS_3D,
                soh8_easfull>(1.0, Kaa, detJ_w, M, cM);
            LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_easfull, MAT::NUM_STRESS_3D,
                NUMDOF_SOH8>(1.0, Kda.A(), detJ_w, M.A(), cb.A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_easfull, MAT::NUM_STRESS_3D, 1>(
                1.0, feas.A(), detJ_w, M.A(), stress.A());
            break;
          case DRT::ELEMENTS::So_hex8::soh8_easmild:
            LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
                soh8_easmild>(cM.A(), cmat.A(), M.A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_easmild, MAT::NUM_STRESS_3D,
                soh8_easmild>(1.0, Kaa, detJ_w, M, cM);
            LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_easmild, MAT::NUM_STRESS_3D,
                NUMDOF_SOH8>(1.0, Kda.A(), detJ_w, M.A(), cb.A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_easmild, MAT::NUM_STRESS_3D, 1>(
                1.0, feas.A(), detJ_w, M.A(), stress.A());
            break;
          case DRT::ELEMENTS::So_hex8::soh8_eassosh8:
            LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
                soh8_eassosh8>(cM.A(), cmat.A(), M.A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_eassosh8, MAT::NUM_STRESS_3D,
                soh8_eassosh8>(1.0, Kaa, detJ_w, M, cM);
            LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_eassosh8, MAT::NUM_STRESS_3D,
                NUMDOF_SOH8>(1.0, Kda.A(), detJ_w, M.A(), cb.A());
            LINALG::DENSEFUNCTIONS::multiplyTN<double, soh8_eassosh8, MAT::NUM_STRESS_3D, 1>(
                1.0, feas.A(), detJ_w, M.A(), stress.A());
            break;
          case DRT::ELEMENTS::So_hex8::soh8_easnone:
            break;
          default:
            dserror("Don't know what to do with EAS type %d", eastype_);
            break;
        }
      }  // ---------------------------------------------------------------- EAS
    }

    if (massmatrix != NULL)  // evaluate mass matrix +++++++++++++++++++++++++
    {
      const double density = Material()->Density(gp);

      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      double ifactor, massfactor;
      for (int inod = 0; inod < NUMNOD_SOH8; ++inod)
      {
        ifactor = shapefcts[gp](inod) * factor;
        for (int jnod = 0; jnod < NUMNOD_SOH8; ++jnod)
        {
          massfactor = shapefcts[gp](jnod) * ifactor;  // intermediate factor
          (*massmatrix)(NUMDIM_SOH8 * inod + 0, NUMDIM_SOH8 * jnod + 0) += massfactor;
          (*massmatrix)(NUMDIM_SOH8 * inod + 1, NUMDIM_SOH8 * jnod + 1) += massfactor;
          (*massmatrix)(NUMDIM_SOH8 * inod + 2, NUMDIM_SOH8 * jnod + 2) += massfactor;
        }
      }

      // check for non constant mass matrix
      if (so3mat->VaryingDensity())
      {
        /*
         If the density, i.e. the mass matrix, is not constant, a linearization is neccessary.
         In general, the mass matrix can be dependent on the displacements, the velocities and the
         accelerations. We write all the additional terms into the mass matrix, hence, conversion
         from accelerations to velocities and displacements are needed. As those conversions depend
         on the time integration scheme, the factors are set within the respective time integrators
         and read from the parameter list inside the element (this is a little ugly...). */
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

        // evaluate derivative of mass w.r.t. to right cauchy green tensor
        so3mat->EvaluateNonLinMass(&defgrd, &glstrain, params, &linmass_disp, &linmass_vel, Id());

        // multiply by 2.0 to get derivative w.r.t green lagrange strains and multiply by time
        // integration factor
        linmass_disp.Scale(2.0 * timintfac_dis);
        linmass_vel.Scale(2.0 * timintfac_vel);
        linmass.Update(1.0, linmass_disp, 1.0, linmass_vel, 0.0);

        // evaluate accelerations at time n+1 at gauss point
        LINALG::Matrix<NUMDIM_SOH8, 1> myacc(true);
        for (int idim = 0; idim < NUMDIM_SOH8; ++idim)
          for (int inod = 0; inod < NUMNOD_SOH8; ++inod)
            myacc(idim) += shapefcts[gp](inod) * (*acc)[idim + (inod * NUMDIM_SOH8)];

        if (stiffmatrix != NULL)
        {
          // integrate linearisation of mass matrix
          //(B^T . d\rho/d disp . a) * detJ * w(gp)
          LINALG::Matrix<1, NUMDOF_SOH8> cb;
          cb.MultiplyTN(linmass_disp, bop);
          for (int inod = 0; inod < NUMNOD_SOH8; ++inod)
          {
            double factor = detJ_w * shapefcts[gp](inod);
            for (int idim = 0; idim < NUMDIM_SOH8; ++idim)
            {
              double massfactor = factor * myacc(idim);
              for (int jnod = 0; jnod < NUMNOD_SOH8; ++jnod)
                for (int jdim = 0; jdim < NUMDIM_SOH8; ++jdim)
                  (*massmatrix)(inod * NUMDIM_SOH8 + idim, jnod * NUMDIM_SOH8 + jdim) +=
                      massfactor * cb(jnod * NUMDIM_SOH8 + jdim);
            }
          }
        }

        // internal force vector without EAS terms
        if (forceinert != NULL)
        {
          // integrate nonlinear inertia force term
          for (int inod = 0; inod < NUMNOD_SOH8; ++inod)
          {
            double forcefactor = shapefcts[gp](inod) * detJ_w;
            for (int idim = 0; idim < NUMDIM_SOH8; ++idim)
              (*forceinert)(inod * NUMDIM_SOH8 + idim) += forcefactor * density * myacc(idim);
          }
        }
      }

    }  // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
       /* =========================================================================*/
  }    /* ==================================================== end of Loop over GP */
       /* =========================================================================*/

  // rhs norm of eas equations
  if (eastype_ != soh8_easnone && split_res && force != NULL)
    // only add for row-map elements
    if (params.get<int>("MyPID") == Owner())
      params.get<double>("cond_rhs_norm") += pow(feas.Norm2(), 2.);

  if (force != NULL && stiffmatrix != NULL)
  {
    // EAS technology: ------------------------------------------------------ EAS
    // subtract EAS matrices from disp-based Kdd to "soften" element
    if (eastype_ != soh8_easnone)
    {
      // we need the inverse of Kaa. Catch Inf/NaN case
      const double norm1 = Kaa.OneNorm();
      if (std::isnan(norm1) || std::isinf(norm1) || norm1 == 0.)
      {
        for (int i = 0; i < Kaa.N(); ++i)
          for (int j = 0; j < Kaa.M(); ++j) Kaa(j, i) = std::numeric_limits<double>::quiet_NaN();
      }
      else
      {
        Epetra_SerialDenseSolver solve_for_inverseKaa;
        solve_for_inverseKaa.SetMatrix(Kaa);
        solve_for_inverseKaa.Invert();
      }

      // EAS-stiffness matrix is: Kdd - Kda^T . Kaa^-1 . Kda
      // EAS-internal force is: fint - Kda^T . Kaa^-1 . feas

      LINALG::SerialDenseMatrix KdaKaa(NUMDOF_SOH8, neas_);  // temporary Kda.Kaa^{-1}
      switch (eastype_)
      {
        case DRT::ELEMENTS::So_hex8::soh8_easfull:
          LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMDOF_SOH8, soh8_easfull, soh8_easfull>(
              KdaKaa, Kda, Kaa);
          LINALG::DENSEFUNCTIONS::multiply<double, NUMDOF_SOH8, soh8_easfull, NUMDOF_SOH8>(
              1.0, stiffmatrix->A(), -1.0, KdaKaa.A(), Kda.A());
          LINALG::DENSEFUNCTIONS::multiply<double, NUMDOF_SOH8, soh8_easfull, 1>(
              1.0, force->A(), -1.0, KdaKaa.A(), feas.A());
          break;
        case DRT::ELEMENTS::So_hex8::soh8_easmild:
          LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMDOF_SOH8, soh8_easmild, soh8_easmild>(
              KdaKaa, Kda, Kaa);
          LINALG::DENSEFUNCTIONS::multiply<double, NUMDOF_SOH8, soh8_easmild, NUMDOF_SOH8>(
              1.0, stiffmatrix->A(), -1.0, KdaKaa.A(), Kda.A());
          LINALG::DENSEFUNCTIONS::multiply<double, NUMDOF_SOH8, soh8_easmild, 1>(
              1.0, force->A(), -1.0, KdaKaa.A(), feas.A());
          break;
        case DRT::ELEMENTS::So_hex8::soh8_eassosh8:
          LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMDOF_SOH8, soh8_eassosh8, soh8_eassosh8>(
              KdaKaa, Kda, Kaa);
          LINALG::DENSEFUNCTIONS::multiply<double, NUMDOF_SOH8, soh8_eassosh8, NUMDOF_SOH8>(
              1.0, stiffmatrix->A(), -1.0, KdaKaa.A(), Kda.A());
          LINALG::DENSEFUNCTIONS::multiply<double, NUMDOF_SOH8, soh8_eassosh8, 1>(
              1.0, force->A(), -1.0, KdaKaa.A(), feas.A());
          break;
        case DRT::ELEMENTS::So_hex8::soh8_easnone:
          break;
        default:
          dserror("Don't know what to do with EAS type %d", eastype_);
          break;
      }

      // store current EAS data in history
      for (int i = 0; i < neas_; ++i)
      {
        for (int j = 0; j < neas_; ++j) (*oldKaainv)(i, j) = Kaa(i, j);
        for (int j = 0; j < NUMDOF_SOH8; ++j) (*oldKda)(i, j) = Kda(i, j);
        (*oldfeas)(i, 0) = feas(i);
      }
    }  // -------------------------------------------------------------------- EAS
  }
  return;
}  // DRT::ELEMENTS::So_hex8::nlnstiffmass


/*----------------------------------------------------------------------*
 |  evaluate the element for GEMM (private)                   popp 09/11|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_nlnstiffmass_gemm(std::vector<int>& lm,  // location matrix
    std::vector<double>& dispo,                                            // old displacements
    std::vector<double>& disp,                                             // current displacements
    std::vector<double>& residual,                                         // current residual displ
    LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* stiffmatrix,         // element stiffness matrix
    LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* massmatrix,          // element mass matrix
    LINALG::Matrix<NUMDOF_SOH8, 1>* force,                         // element internal force vector
    LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestress,    // stresses at GP
    LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestrain,    // strains at GP
    LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* eleplstrain,  // plastic strains at GP
    Teuchos::ParameterList& params,           // algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,    // stress output option
    const INPAR::STR::StrainType iostrain,    // strain output option
    const INPAR::STR::StrainType ioplstrain)  // plastic strain output option
{
  /* ============================================================================*
  ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
  ** ============================================================================*/
  const static std::vector<LINALG::Matrix<NUMNOD_SOH8, 1>> shapefcts = soh8_shapefcts();
  const static std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = soh8_derivs();
  const static std::vector<double> gpweights = soh8_weights();
  /* ============================================================================*/

  // check for prestressing or EAS
  if (pstype_ != INPAR::STR::prestress_none || eastype_ != soh8_easnone)
    dserror("GEMM for Sohex8 not (yet) compatible with EAS / prestressing!");

  // GEMM coefficients
  const double gemmalphaf = params.get<double>("alpha f");
  const double gemmxi = params.get<double>("xi");

  // update element geometry
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe;   // material coord. of element
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xcurr;   // current  coord. of element
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xcurro;  // old  coord. of element

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

    xcurro(i, 0) = xrefe(i, 0) + dispo[i * NODDOF_SOH8 + 0];
    xcurro(i, 1) = xrefe(i, 1) + dispo[i * NODDOF_SOH8 + 1];
    xcurro(i, 2) = xrefe(i, 2) + dispo[i * NODDOF_SOH8 + 2];
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_XYZ;
  // build deformation gradient wrt to material configuration
  LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrd(false);
  LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrdo(false);
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
    defgrdo.MultiplyTT(xcurro, N_XYZ);

    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> rcg;
    LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> rcgo;
    rcg.MultiplyTN(defgrd, defgrd);
    rcgo.MultiplyTN(defgrdo, defgrdo);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain_epetra(MAT::NUM_STRESS_3D);
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(glstrain_epetra.A(), true);
    Epetra_SerialDenseVector glstrain_epetrao(MAT::NUM_STRESS_3D);
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstraino(glstrain_epetrao.A(), true);
    glstrain(0) = 0.5 * (rcg(0, 0) - 1.0);
    glstrain(1) = 0.5 * (rcg(1, 1) - 1.0);
    glstrain(2) = 0.5 * (rcg(2, 2) - 1.0);
    glstrain(3) = rcg(0, 1);
    glstrain(4) = rcg(1, 2);
    glstrain(5) = rcg(2, 0);
    glstraino(0) = 0.5 * (rcgo(0, 0) - 1.0);
    glstraino(1) = 0.5 * (rcgo(1, 1) - 1.0);
    glstraino(2) = 0.5 * (rcgo(2, 2) - 1.0);
    glstraino(3) = rcgo(0, 1);
    glstraino(4) = rcgo(1, 2);
    glstraino(5) = rcgo(2, 0);

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
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> gl;
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
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> invdefgrd;
        invdefgrd.Invert(defgrd);

        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> temp;
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> euler_almansi;
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
    LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8> bop;
    for (int i = 0; i < NUMNOD_SOH8; ++i)
    {
      bop(0, NODDOF_SOH8 * i + 0) = defgrd(0, 0) * N_XYZ(0, i);
      bop(0, NODDOF_SOH8 * i + 1) = defgrd(1, 0) * N_XYZ(0, i);
      bop(0, NODDOF_SOH8 * i + 2) = defgrd(2, 0) * N_XYZ(0, i);
      bop(1, NODDOF_SOH8 * i + 0) = defgrd(0, 1) * N_XYZ(1, i);
      bop(1, NODDOF_SOH8 * i + 1) = defgrd(1, 1) * N_XYZ(1, i);
      bop(1, NODDOF_SOH8 * i + 2) = defgrd(2, 1) * N_XYZ(1, i);
      bop(2, NODDOF_SOH8 * i + 0) = defgrd(0, 2) * N_XYZ(2, i);
      bop(2, NODDOF_SOH8 * i + 1) = defgrd(1, 2) * N_XYZ(2, i);
      bop(2, NODDOF_SOH8 * i + 2) = defgrd(2, 2) * N_XYZ(2, i);
      /* ~~~ */
      bop(3, NODDOF_SOH8 * i + 0) = defgrd(0, 0) * N_XYZ(1, i) + defgrd(0, 1) * N_XYZ(0, i);
      bop(3, NODDOF_SOH8 * i + 1) = defgrd(1, 0) * N_XYZ(1, i) + defgrd(1, 1) * N_XYZ(0, i);
      bop(3, NODDOF_SOH8 * i + 2) = defgrd(2, 0) * N_XYZ(1, i) + defgrd(2, 1) * N_XYZ(0, i);
      bop(4, NODDOF_SOH8 * i + 0) = defgrd(0, 1) * N_XYZ(2, i) + defgrd(0, 2) * N_XYZ(1, i);
      bop(4, NODDOF_SOH8 * i + 1) = defgrd(1, 1) * N_XYZ(2, i) + defgrd(1, 2) * N_XYZ(1, i);
      bop(4, NODDOF_SOH8 * i + 2) = defgrd(2, 1) * N_XYZ(2, i) + defgrd(2, 2) * N_XYZ(1, i);
      bop(5, NODDOF_SOH8 * i + 0) = defgrd(0, 2) * N_XYZ(0, i) + defgrd(0, 0) * N_XYZ(2, i);
      bop(5, NODDOF_SOH8 * i + 1) = defgrd(1, 2) * N_XYZ(0, i) + defgrd(1, 0) * N_XYZ(2, i);
      bop(5, NODDOF_SOH8 * i + 2) = defgrd(2, 2) * N_XYZ(0, i) + defgrd(2, 0) * N_XYZ(2, i);
    }
    LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8> bopo;
    for (int i = 0; i < NUMNOD_SOH8; ++i)
    {
      bopo(0, NODDOF_SOH8 * i + 0) = defgrdo(0, 0) * N_XYZ(0, i);
      bopo(0, NODDOF_SOH8 * i + 1) = defgrdo(1, 0) * N_XYZ(0, i);
      bopo(0, NODDOF_SOH8 * i + 2) = defgrdo(2, 0) * N_XYZ(0, i);
      bopo(1, NODDOF_SOH8 * i + 0) = defgrdo(0, 1) * N_XYZ(1, i);
      bopo(1, NODDOF_SOH8 * i + 1) = defgrdo(1, 1) * N_XYZ(1, i);
      bopo(1, NODDOF_SOH8 * i + 2) = defgrdo(2, 1) * N_XYZ(1, i);
      bopo(2, NODDOF_SOH8 * i + 0) = defgrdo(0, 2) * N_XYZ(2, i);
      bopo(2, NODDOF_SOH8 * i + 1) = defgrdo(1, 2) * N_XYZ(2, i);
      bopo(2, NODDOF_SOH8 * i + 2) = defgrdo(2, 2) * N_XYZ(2, i);
      /* ~~~ */
      bopo(3, NODDOF_SOH8 * i + 0) = defgrdo(0, 0) * N_XYZ(1, i) + defgrdo(0, 1) * N_XYZ(0, i);
      bopo(3, NODDOF_SOH8 * i + 1) = defgrdo(1, 0) * N_XYZ(1, i) + defgrdo(1, 1) * N_XYZ(0, i);
      bopo(3, NODDOF_SOH8 * i + 2) = defgrdo(2, 0) * N_XYZ(1, i) + defgrdo(2, 1) * N_XYZ(0, i);
      bopo(4, NODDOF_SOH8 * i + 0) = defgrdo(0, 1) * N_XYZ(2, i) + defgrdo(0, 2) * N_XYZ(1, i);
      bopo(4, NODDOF_SOH8 * i + 1) = defgrdo(1, 1) * N_XYZ(2, i) + defgrdo(1, 2) * N_XYZ(1, i);
      bopo(4, NODDOF_SOH8 * i + 2) = defgrdo(2, 1) * N_XYZ(2, i) + defgrdo(2, 2) * N_XYZ(1, i);
      bopo(5, NODDOF_SOH8 * i + 0) = defgrdo(0, 2) * N_XYZ(0, i) + defgrdo(0, 0) * N_XYZ(2, i);
      bopo(5, NODDOF_SOH8 * i + 1) = defgrdo(1, 2) * N_XYZ(0, i) + defgrdo(1, 0) * N_XYZ(2, i);
      bopo(5, NODDOF_SOH8 * i + 2) = defgrdo(2, 2) * N_XYZ(0, i) + defgrdo(2, 0) * N_XYZ(2, i);
    }

    // GEMM: computed averaged mid-point quantities

    // non-linear mid-B-operator
    // B_m = (1.0-gemmalphaf)*B_{n+1} + gemmalphaf*B_{n}
    LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8> bopm;
    bopm.Update(1.0 - gemmalphaf, bop, gemmalphaf, bopo);

    // mid-strain GL vector
    // E_m = (1.0-gemmalphaf+gemmxi)*E_{n+1} + (gemmalphaf-gemmxi)*E_n
    Epetra_SerialDenseVector glstrain_epetram(MAT::NUM_STRESS_3D);
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrainm(glstrain_epetram.A(), true);
    glstrainm.Update(1.0 - gemmalphaf + gemmxi, glstrain, gemmalphaf - gemmxi, glstraino);

    // call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> cmat(true);
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> stress(true);
    double density;
    SolidMaterial()->EvaluateGEMM(
        &stress, &cmat, &density, &glstrainm, &glstrain, &glstraino, &rcg, &rcgo, Id());
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

    double detJ_w = detJ * gpweights[gp];
    // update internal force vector
    if (force != NULL)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w, bopm, stress, 1.0);
    }

    // update stiffness matrix
    if (stiffmatrix != NULL)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      const double faceu = (1.0 - gemmalphaf + gemmxi) * detJ_w;
      LINALG::Matrix<6, NUMDOF_SOH8> cb;
      cb.Multiply(cmat, bop);                         // B_{n+1} here!!!
      stiffmatrix->MultiplyTN(faceu, bopm, cb, 1.0);  // B_m here!!!

      // integrate `geometric' stiffness matrix
      const double facg = (1.0 - gemmalphaf) * detJ_w;
      LINALG::Matrix<6, 1> sfac(stress);  // auxiliary integrated stress
      sfac.Scale(facg);                   // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
      std::vector<double> SmB_L(3);       // intermediate Sm.B_L
      // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
      for (int inod = 0; inod < NUMNOD_SOH8; ++inod)
      {
        SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod) + sfac(5) * N_XYZ(2, inod);
        SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod) + sfac(4) * N_XYZ(2, inod);
        SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod) + sfac(2) * N_XYZ(2, inod);
        for (int jnod = 0; jnod < NUMNOD_SOH8; ++jnod)
        {
          double bopstrbop = 0.0;  // intermediate value
          for (int idim = 0; idim < NUMDIM_SOH8; ++idim)
            bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
          (*stiffmatrix)(3 * inod + 0, 3 * jnod + 0) += bopstrbop;
          (*stiffmatrix)(3 * inod + 1, 3 * jnod + 1) += bopstrbop;
          (*stiffmatrix)(3 * inod + 2, 3 * jnod + 2) += bopstrbop;
        }
      }  // end of integrate `geometric' stiffness******************************
    }

    if (massmatrix != NULL)  // evaluate mass matrix +++++++++++++++++++++++++
    {
      double density = Material()->Density(gp);

      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      double ifactor, massfactor;
      for (int inod = 0; inod < NUMNOD_SOH8; ++inod)
      {
        ifactor = shapefcts[gp](inod) * factor;
        for (int jnod = 0; jnod < NUMNOD_SOH8; ++jnod)
        {
          massfactor = shapefcts[gp](jnod) * ifactor;  // intermediate factor
          (*massmatrix)(NUMDIM_SOH8 * inod + 0, NUMDIM_SOH8 * jnod + 0) += massfactor;
          (*massmatrix)(NUMDIM_SOH8 * inod + 1, NUMDIM_SOH8 * jnod + 1) += massfactor;
          (*massmatrix)(NUMDIM_SOH8 * inod + 2, NUMDIM_SOH8 * jnod + 2) += massfactor;
        }
      }

    }  // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
       /* =========================================================================*/
  }    /* ==================================================== end of Loop over GP */
       /* =========================================================================*/

  return;
}  // DRT::ELEMENTS::So_hex8::soh8_nlnstiffmass_gemm


/*----------------------------------------------------------------------*
 |  lump mass matrix (private)                               bborn 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_lumpmass(LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* emass)
{
  // lump mass matrix
  if (emass != NULL)
  {
    // we assume #elemat2 is a square matrix
    for (unsigned int c = 0; c < (*emass).N(); ++c)  // parse columns
    {
      double d = 0.0;
      for (unsigned int r = 0; r < (*emass).M(); ++r)  // parse rows
      {
        d += (*emass)(r, c);  // accumulate row entries
        (*emass)(r, c) = 0.0;
      }
      (*emass)(c, c) = d;  // apply sum of row entries on diagonal
    }
  }
}

/*----------------------------------------------------------------------*
 |  Evaluate Hex8 Shape fcts at all 8 Gauss Points             maf 05/08|
 *----------------------------------------------------------------------*/
const std::vector<LINALG::Matrix<NUMNOD_SOH8, 1>> DRT::ELEMENTS::So_hex8::soh8_shapefcts() const
{
  std::vector<LINALG::Matrix<NUMNOD_SOH8, 1>> shapefcts(NUMGPT_SOH8);

  // fill up nodal f at each gp
  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
  {
    const LINALG::Matrix<NUMDIM_SOH8, 1> rst_gp(gp_rule_.Point(gp), true);
    DRT::UTILS::shape_function<DRT::Element::hex8>(rst_gp, shapefcts[gp]);
  }

  return shapefcts;
}


/*----------------------------------------------------------------------*
 |  Evaluate Hex8 Shape fct derivs at all 8 Gauss Points       maf 05/08|
 *----------------------------------------------------------------------*/
const std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> DRT::ELEMENTS::So_hex8::soh8_derivs()
    const
{
  std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs(NUMGPT_SOH8);

  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
  {
    const LINALG::Matrix<NUMDIM_SOH8, 1> rst_gp(gp_rule_.Point(gp), true);
    DRT::UTILS::shape_function_deriv1<DRT::Element::hex8>(rst_gp, derivs[gp]);
  }

  return derivs;
}

/*----------------------------------------------------------------------*
 |  Evaluate Hex8 Weights at all 8 Gauss Points                maf 05/08|
 *----------------------------------------------------------------------*/
const std::vector<double> DRT::ELEMENTS::So_hex8::soh8_weights() const
{
  std::vector<double> weights(NUMGPT_SOH8);
  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp) weights[gp] = gp_rule_.Weight(gp);

  return weights;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_create_eas_backup_state(const std::vector<double>& displ_incr)
{
  if (eastype_ == soh8_easnone) return;

  // --- create EAS state backup ----------------------------------------------
  {
    const Epetra_SerialDenseMatrix* alpha = data_.Get<Epetra_SerialDenseMatrix>("alpha");
    if (not alpha) dserror("Can't access the current enhanced strain state.");

    Epetra_SerialDenseMatrix* alpha_backup_ptr =
        data_.GetMutable<Epetra_SerialDenseMatrix>("alpha_backup");
    if (alpha_backup_ptr)
      *alpha_backup_ptr = *alpha;
    else
      data_.Add("alpha_backup", *alpha);
  }

  // --- create EAS increment backup ------------------------------------------
  {
    // compute the current eas increment
    Epetra_SerialDenseMatrix eas_inc(neas_, 1);
    soh8_computeEASInc(displ_incr, &eas_inc);

    Epetra_SerialDenseMatrix* eas_inc_backup_ptr =
        data_.GetMutable<Epetra_SerialDenseMatrix>("eas_inc_backup");
    if (eas_inc_backup_ptr)
      *eas_inc_backup_ptr = eas_inc;
    else
      data_.Add("eas_inc_backup", eas_inc);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_recover_from_eas_backup_state()
{
  if (eastype_ == soh8_easnone) return;

  Epetra_SerialDenseMatrix* alpha = NULL;
  Epetra_SerialDenseMatrix* eas_inc = NULL;

  // --- recover state from EAS backup ----------------------------------------
  {
    const Epetra_SerialDenseMatrix* alpha_backup =
        data_.Get<Epetra_SerialDenseMatrix>("alpha_backup");
    if (not alpha_backup)
      dserror(
          "Can't access the enhanced strain backup state. Did you "
          "create a backup? See soh8_create_eas_backup_state().");

    alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");
    if (not alpha) dserror("Can't access the enhanced strain state.");

    *alpha = *alpha_backup;
  }

  // --- recover increment from EAS backup ------------------------------------
  {
    const Epetra_SerialDenseMatrix* eas_inc_backup =
        data_.Get<Epetra_SerialDenseMatrix>("eas_inc_backup");
    if (not eas_inc_backup)
      dserror(
          "Can't access the enhanced strain increment backup. Did you "
          "create a backup? See soh8_create_eas_backup_state().");

    eas_inc = data_.GetMutable<Epetra_SerialDenseMatrix>("eas_inc");
    if (not eas_inc) dserror("Can't access the enhanced strain increment.");

    *eas_inc = *eas_inc_backup;
  }

  // Finally, we have to update the backup state, otherwise a follow-up
  // step length adaption will lead to a wrong eas state.
  {
    old_step_length_ = 0.0;
  }
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                  gee 04/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8Type::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So_hex8* actele = dynamic_cast<DRT::ELEMENTS::So_hex8*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_hex8* failed");
    actele->InitJacobianMapping();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  compute def gradient at every gaussian point (protected)   gee 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::DefGradient(const std::vector<double>& disp,
    Epetra_SerialDenseMatrix& gpdefgrd, DRT::ELEMENTS::PreStress& prestress)
{
  const static std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = soh8_derivs();

  // update element geometry
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xdisp;  // current  coord. of element
  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    xdisp(i, 0) = disp[i * NODDOF_SOH8 + 0];
    xdisp(i, 1) = disp[i * NODDOF_SOH8 + 1];
    xdisp(i, 2) = disp[i * NODDOF_SOH8 + 2];
  }

  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
  {
    // get Jacobian mapping wrt to the stored deformed configuration
    LINALG::Matrix<3, 3> invJdef;
    prestress.StoragetoMatrix(gp, invJdef, prestress.JHistory());

    // by N_XYZ = J^-1 * N_rst
    LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_xyz;
    N_xyz.Multiply(invJdef, derivs[gp]);

    // build defgrd (independent of xrefe!)
    LINALG::Matrix<3, 3> defgrd(true);
    if (kintype_ == INPAR::STR::kinem_nonlinearTotLag)
    {
      defgrd.MultiplyTT(xdisp, N_xyz);
    }
    defgrd(0, 0) += 1.0;
    defgrd(1, 1) += 1.0;
    defgrd(2, 2) += 1.0;

    prestress.MatrixtoStorage(gp, defgrd, gpdefgrd);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  compute Jac.mapping wrt deformed configuration (protected) gee 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::UpdateJacobianMapping(
    const std::vector<double>& disp, DRT::ELEMENTS::PreStress& prestress)
{
  const static std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = soh8_derivs();

  // get incremental disp
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xdisp;
  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    xdisp(i, 0) = disp[i * NODDOF_SOH8 + 0];
    xdisp(i, 1) = disp[i * NODDOF_SOH8 + 1];
    xdisp(i, 2) = disp[i * NODDOF_SOH8 + 2];
  }

  LINALG::Matrix<3, 3> invJhist;
  LINALG::Matrix<3, 3> invJ;
  LINALG::Matrix<3, 3> defgrd(true);
  LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_xyz;
  LINALG::Matrix<3, 3> invJnew;
  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
  {
    // get the invJ old state
    prestress.StoragetoMatrix(gp, invJhist, prestress.JHistory());
    // get derivatives wrt to invJhist
    N_xyz.Multiply(invJhist, derivs[gp]);
    // build defgrd \partial x_new / \parial x_old , where x_old != X
    if (kintype_ == INPAR::STR::kinem_nonlinearTotLag)
    {
      defgrd.MultiplyTT(xdisp, N_xyz);
    }
    defgrd(0, 0) += 1.0;
    defgrd(1, 1) += 1.0;
    defgrd(2, 2) += 1.0;
    // make inverse of this defgrd
    defgrd.Invert();
    // push-forward of Jinv
    invJnew.MultiplyTN(defgrd, invJhist);
    // store new reference configuration
    prestress.MatrixtoStorage(gp, invJnew, prestress.JHistory());
  }  // for (int gp=0; gp<NUMGPT_SOH8; ++gp)

  return;
}

/*---------------------------------------------------------------------------------------------*
 |  Update history variables (e.g. remodeling of fiber directions) (protected)      braeu 07/16|
 *---------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::Update_element(
    std::vector<double>& disp, Teuchos::ParameterList& params, Teuchos::RCP<MAT::Material> mat)
{
  bool remodel;
  const Teuchos::ParameterList& patspec = DRT::Problem::Instance()->PatSpecParams();
  remodel = DRT::INPUT::IntegralValue<int>(patspec, "REMODEL");

  // Calculate current deformation gradient
  if ((mat->MaterialType() == INPAR::MAT::m_constraintmixture) ||
      (mat->MaterialType() == INPAR::MAT::m_elasthyper) ||
      (mat->MaterialType() == INPAR::MAT::m_growthremodel_elasthyper) ||
      (SolidMaterial()->UsesExtendedUpdate()))
  {
    // in a first step ommit everything with prestress and EAS!!

    const static std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = soh8_derivs();

    // update element geometry
    LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe;  // reference coord. of element
    LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xcurr;  // current  coord. of element
    LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xdisp;
    DRT::Node** nodes = Nodes();
    for (int i = 0; i < NUMNOD_SOH8; ++i)
    {
      const double* x = nodes[i]->X();
      xrefe(i, 0) = x[0];
      xrefe(i, 1) = x[1];
      xrefe(i, 2) = x[2];
      xcurr(i, 0) = x[0] + disp[i * NODDOF_SOH8 + 0];
      xcurr(i, 1) = x[1] + disp[i * NODDOF_SOH8 + 1];
      xcurr(i, 2) = x[2] + disp[i * NODDOF_SOH8 + 2];

      if (pstype_ == INPAR::STR::prestress_mulf)
      {
        xdisp(i, 0) = disp[i * NODDOF_SOH8 + 0];
        xdisp(i, 1) = disp[i * NODDOF_SOH8 + 1];
        xdisp(i, 2) = disp[i * NODDOF_SOH8 + 2];
      }
    }
    /* =========================================================================*/
    /* ================================================= Loop over Gauss Points */
    /* =========================================================================*/
    LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_XYZ;
    // interpolated values of stress and defgrd for remodeling
    LINALG::Matrix<3, 3> avg_stress(true);
    LINALG::Matrix<3, 3> avg_defgrd(true);

    // build deformation gradient wrt to material configuration
    LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrd(false);
    params.set<int>("numgp", static_cast<int>(NUMGPT_SOH8));

    // center of element in reference configuration
    LINALG::Matrix<1, NUMDIM_SOH8> point(false);
    point.Clear();
    soh8_ElementCenterRefeCoords(point, xrefe);
    params.set("elecenter", point);

    for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
    {
      soh8_GaussPointRefeCoords(point, xrefe, gp);
      params.set("gprefecoord", point);


      /* get the inverse of the Jacobian matrix which looks like:
       **            [ x_,r  y_,r  z_,r ]^-1
       **     J^-1 = [ x_,s  y_,s  z_,s ]
       **            [ x_,t  y_,t  z_,t ]
       */
      // compute derivatives N_XYZ at gp w.r.t. material coordinates
      // by N_XYZ = J^-1 * N_rst
      N_XYZ.Multiply(invJ_[gp], derivs[gp]);

      if (pstype_ == INPAR::STR::prestress_mulf)
      {
        // get Jacobian mapping wrt to the stored configuration
        LINALG::Matrix<3, 3> invJdef;
        prestress_->StoragetoMatrix(gp, invJdef, prestress_->JHistory());
        // get derivatives wrt to last spatial configuration
        LINALG::Matrix<3, 8> N_xyz;
        N_xyz.Multiply(invJdef, derivs[gp]);

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
        // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
        defgrd.MultiplyTT(xcurr, N_XYZ);


      // call material update if material = m_growthremodel_elasthyper (calculate and update
      // inelastic deformation gradient)
      if (SolidMaterial()->UsesExtendedUpdate())
      {
        SolidMaterial()->Update(defgrd, gp, params, Id());
      }

      // determine new fiber directions
      if (remodel)
      {
        // determine new fiber directions
        // Right Cauchy-Green tensor = F^T * F
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> cauchygreen;
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

        // call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
        LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> cmat(true);
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1> stress(true);
        params.set<int>("gp", gp);
        SolidMaterial()->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, Id());
        // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

        // Cauchy stress
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

        LINALG::Matrix<3, 3> temp(true);
        LINALG::Matrix<3, 3> cauchystress(true);
        temp.Multiply(1.0 / detF, defgrd, pkstress);
        cauchystress.MultiplyNT(temp, defgrd);

        // evaluate eigenproblem based on stress of previous step
        LINALG::Matrix<3, 3> lambda(true);
        LINALG::Matrix<3, 3> locsys(true);
        LINALG::SYEV(cauchystress, lambda, locsys);

        if (mat->MaterialType() == INPAR::MAT::m_constraintmixture)
        {
          MAT::ConstraintMixture* comi = static_cast<MAT::ConstraintMixture*>(mat.get());
          comi->EvaluateFiberVecs(gp, locsys, defgrd);
        }
        else if (mat->MaterialType() == INPAR::MAT::m_elasthyper)
        {
          // we only have fibers at element center, thus we interpolate stress and defgrd
          avg_stress.Update(1.0 / NUMGPT_SOH8, cauchystress, 1.0);
          avg_defgrd.Update(1.0 / NUMGPT_SOH8, defgrd, 1.0);
        }
        else
          dserror("material not implemented for remodeling");
      }
    }  // end loop over gauss points

    // determine new fiber directions
    if (remodel)
    {
      if (mat->MaterialType() == INPAR::MAT::m_elasthyper)
      {
        // evaluate eigenproblem based on stress of previous step
        LINALG::Matrix<3, 3> lambda(true);
        LINALG::Matrix<3, 3> locsys(true);
        LINALG::SYEV(avg_stress, lambda, locsys);

        // modulation function acc. Hariton: tan g = 2nd max lambda / max lambda
        double newgamma = atan2(lambda(1, 1), lambda(2, 2));
        // compression in 2nd max direction, thus fibers are alligned to max principal direction
        if (lambda(1, 1) < 0) newgamma = 0.0;

        // new fiber vectors
        MAT::ElastHyper* elast = static_cast<MAT::ElastHyper*>(mat.get());
        elast->EvaluateFiberVecs(newgamma, locsys, avg_defgrd);
      }
    }
  }

  // do something with internal EAS, etc parameters
  if (eastype_ != soh8_easnone)
  {
    Epetra_SerialDenseMatrix* alpha =
        data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // Alpha_{n+1}
    Epetra_SerialDenseMatrix* alphao =
        data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");  // Alpha_n
    // alphao := alpha
    switch (eastype_)
    {
      case DRT::ELEMENTS::So_hex8::soh8_easfull:
        LINALG::DENSEFUNCTIONS::update<double, soh8_easfull, 1>(*alphao, *alpha);
        break;
      case DRT::ELEMENTS::So_hex8::soh8_easmild:
        LINALG::DENSEFUNCTIONS::update<double, soh8_easmild, 1>(*alphao, *alpha);
        break;
      case DRT::ELEMENTS::So_hex8::soh8_eassosh8:
        LINALG::DENSEFUNCTIONS::update<double, soh8_eassosh8, 1>(*alphao, *alpha);
        break;
      case DRT::ELEMENTS::So_hex8::soh8_easnone:
        break;
      default:
        dserror("Don't know what to do with EAS type %d", eastype_);
        break;
    }
  }
  SolidMaterial()->Update();

  return;
}


/*----------------------------------------------------------------------*
 | push forward of material to spatial stresses              dano 11/12 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::GLtoEA(LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain,
    LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>* defgrd,
    LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>* euler_almansi)
{
  // e = F^{T-1} . E . F^{-1}

  // rewrite Green-Lagrange strain in tensor notation
  LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> gl;
  gl(0, 0) = (*glstrain)(0);
  gl(0, 1) = 0.5 * (*glstrain)(3);
  gl(0, 2) = 0.5 * (*glstrain)(5);
  gl(1, 0) = gl(0, 1);
  gl(1, 1) = (*glstrain)(1);
  gl(1, 2) = 0.5 * (*glstrain)(4);
  gl(2, 0) = gl(0, 2);
  gl(2, 1) = gl(1, 2);
  gl(2, 2) = (*glstrain)(2);

  // inverse of deformation gradient
  LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> invdefgrd;
  invdefgrd.Invert((*defgrd));

  // (3x3) = (3x3) (3x3) (3x3)
  LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> temp;
  temp.Multiply(gl, invdefgrd);
  (*euler_almansi).MultiplyTN(invdefgrd, temp);

}  // GLtoEA()


/*----------------------------------------------------------------------*
 | push forward of material to spatial stresses              dano 11/12 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::PK2toCauchy(LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* stress,
    LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>* defgrd,
    LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>* cauchystress)
{
  // calculate the Jacobi-deterinant
  const double detF = (*defgrd).Determinant();

  // sigma = 1/J . F . S . F^T
  LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> pkstress;
  pkstress(0, 0) = (*stress)(0);
  pkstress(0, 1) = (*stress)(3);
  pkstress(0, 2) = (*stress)(5);
  pkstress(1, 0) = pkstress(0, 1);
  pkstress(1, 1) = (*stress)(1);
  pkstress(1, 2) = (*stress)(4);
  pkstress(2, 0) = pkstress(0, 2);
  pkstress(2, 1) = pkstress(1, 2);
  pkstress(2, 2) = (*stress)(2);

  LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> temp;
  temp.Multiply((1.0 / detF), (*defgrd), pkstress);
  (*cauchystress).MultiplyNT(temp, (*defgrd));

}  // PK2toCauchy()

/*----------------------------------------------------------------------*
 |  Calculate consistent deformation gradient               seitz 04/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::CalcConsistentDefgrd(LINALG::Matrix<3, 3> defgrd_disp,
    LINALG::Matrix<6, 1> glstrain_mod, LINALG::Matrix<3, 3>& defgrd_mod) const
{
  LINALG::Matrix<3, 3> R;       // rotation tensor
  LINALG::Matrix<3, 3> U_mod;   // modified right stretch tensor
  LINALG::Matrix<3, 3> U_disp;  // displacement-based right stretch tensor
  LINALG::Matrix<3, 3> EW;      // temporarily store eigenvalues
  LINALG::Matrix<3, 3> tmp;     // temporary matrix for matrix matrix matrix products
  LINALG::Matrix<3, 3> tmp2;    // temporary matrix for matrix matrix matrix products

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

  LINALG::SYEV(U_mod, EW, U_mod);
  for (int i = 0; i < 3; ++i) EW(i, i) = sqrt(EW(i, i));
  tmp.Multiply(U_mod, EW);
  tmp2.MultiplyNT(tmp, U_mod);
  U_mod.Update(tmp2);

  // ******************************************************************
  // calculate displacement-based right stretch tensor
  // ******************************************************************
  U_disp.MultiplyTN(defgrd_disp, defgrd_disp);

  LINALG::SYEV(U_disp, EW, U_disp);
  for (int i = 0; i < 3; ++i) EW(i, i) = sqrt(EW(i, i));
  tmp.Multiply(U_disp, EW);
  tmp2.MultiplyNT(tmp, U_disp);
  U_disp.Update(tmp2);

  // ******************************************************************
  // compose consistent deformation gradient
  // ******************************************************************
  U_disp.Invert();
  R.Multiply(defgrd_disp, U_disp);
  defgrd_mod.Multiply(R, U_mod);

  // you're done here
  return;
}

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | check the constitutive tensor and/or use the approximation as        |
 | elastic stiffness matrix                                  rauch 07/13|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::EvaluateFiniteDifferenceMaterialTangent(
    LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* stiffmatrix,
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& stress, std::vector<double>& disp,
    const double detJ_w, const double detJ, const double detJ0, const double charelelength,
    const LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8>& bop,
    const LINALG::Matrix<6, NUMDOF_SOH8>& cb, const LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>& N_XYZ,
    const LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& T0invT,
    const std::vector<Epetra_SerialDenseMatrix>* M_GP, const Epetra_SerialDenseMatrix* alpha,
    LINALG::SerialDenseMatrix& M, const int gp, Teuchos::ParameterList& params)
{
  // build elastic stiffness matrix directly by finite differences

  Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_static_cast<MAT::So3Material>(Material());

#ifdef MATERIALFDCHECK
  static LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8> stiffmatrix_analytical;
  static LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8> stiffmatrix_fd;
  bool success = true;
  if (gp == 0)
  {
    stiffmatrix_analytical.PutScalar(0.0);
    stiffmatrix_fd.PutScalar(0.0);
  }
  stiffmatrix_analytical.MultiplyTN(detJ_w, bop, cb, 1.0);
#endif

  const double delta = charelelength * 1.0e-08;

  // matrices and vectors
  LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8> cb_fd(true);
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> stress_fd(true);
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> finitedifference(true);

  // update element geometry
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe;  // reference coord. of element
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xcurr;  // current   coord. of element
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xdisp;

  // get nodes
  DRT::Node** nodes = Nodes();

  //////////////////////////////////////////////////////////////////////////////
  ////// evaluate partial derivatives of stress (S(d_n+delta) - S(d_n))/delta
  /////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////
  // loop over columns
  /////////////////////////////////////////////
  for (int i = 0; i < NUMDOF_SOH8; ++i)
  {
    // undo disturbance for disp[i-1]
    if (i > 0) disp[i - 1] -= delta;

    // disturb displacements
    disp[i] += delta;

    for (int k = 0; k < NUMNOD_SOH8; ++k)
    {
      const double* x = nodes[k]->X();
      xrefe(k, 0) = x[0];
      xrefe(k, 1) = x[1];
      xrefe(k, 2) = x[2];

      xcurr(k, 0) = xrefe(k, 0) + disp[k * NODDOF_SOH8 + 0];
      xcurr(k, 1) = xrefe(k, 1) + disp[k * NODDOF_SOH8 + 1];
      xcurr(k, 2) = xrefe(k, 2) + disp[k * NODDOF_SOH8 + 2];
    }

    // build deformation gradient wrt to material configuration

    LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrd_fd(true);
    defgrd_fd.MultiplyTT(xcurr, N_XYZ);


    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> cauchygreen_fd(true);
    cauchygreen_fd.MultiplyTN(defgrd_fd, defgrd_fd);

    // Green-Lagrange strains matrix E = 0.5 * (cauchygreen_fd - Identity)
    // GL strain vector glstrain_fd={E11,E22,E33,2*E12,2*E23,2*E31}
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain_fd(true);
    glstrain_fd(0) = 0.5 * (cauchygreen_fd(0, 0) - 1.0);
    glstrain_fd(1) = 0.5 * (cauchygreen_fd(1, 1) - 1.0);
    glstrain_fd(2) = 0.5 * (cauchygreen_fd(2, 2) - 1.0);
    glstrain_fd(3) = cauchygreen_fd(0, 1);
    glstrain_fd(4) = cauchygreen_fd(1, 2);
    glstrain_fd(5) = cauchygreen_fd(2, 0);

    // deformation gradient consistent with (potentially EAS-modified) GL strains
    // without eas this is equal to the regular defgrd_fd.
    LINALG::Matrix<3, 3> defgrd_fd_mod(defgrd_fd);

    // EAS technology: "enhance the strains"  ----------------------------- EAS
    if (eastype_ != soh8_easnone)
    {
      dserror("be careful ! fdcheck has not been tested with EAS, yet! ");
      M.LightShape(MAT::NUM_STRESS_3D, neas_);
      // map local M to global, also enhancement is refered to element origin
      // M = detJ0/detJ T0^{-T} . M
      // Epetra_SerialDenseMatrix Mtemp(M); // temp M for Matrix-Matrix-Product
      // add enhanced strains = M . alpha to GL strains to "unlock" element
      switch (eastype_)
      {
        case DRT::ELEMENTS::So_hex8::soh8_easfull:
          LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
              soh8_easfull>(M.A(), detJ0 / detJ, T0invT.A(), (M_GP->at(gp)).A());
          LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, soh8_easfull, 1>(
              1.0, glstrain_fd.A(), 1.0, M.A(), alpha->A());
          break;
        case DRT::ELEMENTS::So_hex8::soh8_easmild:
          LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
              soh8_easmild>(M.A(), detJ0 / detJ, T0invT.A(), (M_GP->at(gp)).A());
          LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, soh8_easmild, 1>(
              1.0, glstrain_fd.A(), 1.0, M.A(), alpha->A());
          break;
        case DRT::ELEMENTS::So_hex8::soh8_eassosh8:
          LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
              soh8_eassosh8>(M.A(), detJ0 / detJ, T0invT.A(), (M_GP->at(gp)).A());
          LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, soh8_eassosh8, 1>(
              1.0, glstrain_fd.A(), 1.0, M.A(), alpha->A());
          break;
        case DRT::ELEMENTS::So_hex8::soh8_easnone:
          break;
        default:
          dserror("Don't know what to do with EAS type %d", eastype_);
          break;
      }

      // calculate deformation gradient consistent with modified GL strain tensor
      if (Teuchos::rcp_static_cast<MAT::So3Material>(Material())->NeedsDefgrd())
        CalcConsistentDefgrd(defgrd_fd, glstrain_fd, defgrd_fd_mod);
    }  // ------------------------------------------------------------------ EAS

    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> cmat_fd;
    so3mat->Evaluate(&defgrd_fd_mod, &glstrain_fd, params, &stress_fd, &cmat_fd, Id());

    // finite difference approximation of partial derivative
    //
    //      d S_ij,columnindex
    //
    //
    finitedifference.Update(1.0, stress_fd, -1.0, stress);
    finitedifference.Scale(1.0 / delta);

    /////////////////////////
    // loop over rows
    ////////////////////////
    for (int j = 0; j < MAT::NUM_STRESS_3D; ++j)
    {
      cb_fd(j, i) = finitedifference(j, 0);
    }  // j-loop (rows)


    // reset disp at last loop execution
    if (i == (NUMDOF_SOH8 - 1))
    {
      disp[i] -= delta;

      // reset xcurr
      for (int k = 0; k < NUMNOD_SOH8; ++k)
      {
        const double* x = nodes[k]->X();
        xrefe(k, 0) = x[0];
        xrefe(k, 1) = x[1];
        xrefe(k, 2) = x[2];

        xcurr(k, 0) = xrefe(k, 0) + disp[k * NODDOF_SOH8 + 0];
        xcurr(k, 1) = xrefe(k, 1) + disp[k * NODDOF_SOH8 + 1];
        xcurr(k, 2) = xrefe(k, 2) + disp[k * NODDOF_SOH8 + 2];
      }
    }

  }  // i-loop (columns)
  /////////////////////////////// FD LOOP


  ///////////////////////////////////////
  // build approximated stiffness matrix
  ///////////////////////////////////////
  stiffmatrix->MultiplyTN(detJ_w, bop, cb_fd, 1.0);
#ifdef MATERIALFDCHECK
  stiffmatrix_fd.MultiplyTN(detJ_w, bop, cb_fd, 1.0);
#endif

  ///////////////////////////////////////
#ifdef MATERIALFDCHECK
  // after last gp was evaluated
  if (gp == (NUMGPT_SOH8 - 1))
  {
    LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8> errormatrix(true);

    // calc error (subtraction stiffmatrix - stiffmatrix_analytical)
    errormatrix.Update(1.0, stiffmatrix_fd, -1.0, stiffmatrix_analytical);

    for (int i = 0; i < NUMDOF_SOH8; ++i)
    {
      for (int j = 0; j < NUMDOF_SOH8; ++j)
      {
        if (1)
        {
          double relerror = abs(errormatrix(i, j)) / abs((stiffmatrix_analytical)(i, j));
          if (std::min(abs(errormatrix(i, j)), relerror) > delta * 1000.0)
          {
            std::cout << "ELEGID:" << this->Id() << "  gp: " << gp << "  ROW: " << i
                      << "  COL: " << j << "    REL. ERROR: " << relerror
                      << "    ABS. ERROR: " << abs(errormatrix(i, j))
                      << "    stiff. val: " << stiffmatrix_analytical(i, j)
                      << "    approx. val: " << stiffmatrix_fd(i, j) << std::endl;
            success = false;
          }
        }
      }
    }  // check errors

    /*  if (!success)
      {
        std::cout<<"FDCHECK FAILED!"<<std::endl;
        dserror("encountered errors checking the elastic stiffness matrix");
      }*/

  }  // if last gp of element is reached
#endif


  return;
}


/*----------------------------------------------------------------------*
 |  evaluate cauchy stress tensor                           seitz 11/16|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::GetCauchyAtXi(const LINALG::Matrix<3, 1>& xi,
    const std::vector<double>& disp, const LINALG::Matrix<3, 1>& n, const LINALG::Matrix<3, 1>& t,
    double& sigma_nt, Epetra_SerialDenseMatrix* DsntDd, Epetra_SerialDenseMatrix* D2sntDd2,
    Epetra_SerialDenseMatrix* D2sntDdDn, Epetra_SerialDenseMatrix* D2sntDdDt,
    Epetra_SerialDenseMatrix* D2sntDdDxi, LINALG::Matrix<3, 1>* DsntDn,
    LINALG::Matrix<3, 1>* DsntDt, LINALG::Matrix<3, 1>* DsntDxi, const std::vector<double>* temp,
    Epetra_SerialDenseMatrix* DsntDT, Epetra_SerialDenseMatrix* D2sntDdDT)
{
  if (temp || DsntDT || D2sntDdDT)
    dserror("Thermo-elastic Nitsche contact not yet implemented in so hex8");

  sigma_nt = 0.0;

  static LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe(true);  // reference coord. of element
  static LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xcurr(true);  // current  coord. of element
  xrefe.Clear();
  xcurr.Clear();
  DRT::Node** nodes = Nodes();

  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    const double* x = nodes[i]->X();
    for (int d = 0; d < NUMDIM_SOH8; ++d)
    {
      xrefe(i, d) = x[d];
      xcurr(i, d) = xrefe(i, d) + disp[i * NODDOF_SOH8 + d];
    }
  }

  static LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> deriv(true);
  deriv.Clear();
  DRT::UTILS::shape_function_deriv1<DRT::Element::hex8>(xi, deriv);

  static LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> N_XYZ(true);
  static LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> invJ(true);
  invJ.Multiply(1.0, deriv, xrefe, 0.0);
  invJ.Invert();
  N_XYZ.Multiply(1.0, invJ, deriv, 0.0);
  static LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> defgrd(true);
  defgrd.MultiplyTT(1.0, xcurr, N_XYZ, 0.0);

  // linearization of deformation gradient F w.r.t. displacements
  static LINALG::Matrix<9, NUMDOF_SOH8> DFDd(true);
  DFDd.Clear();
  if (DsntDd || D2sntDdDn || D2sntDdDt || D2sntDd2 || D2sntDdDxi)
  {
    for (int i = 0; i < NUMNOD_SOH8; ++i)
    {
      DFDd(0, NODDOF_SOH8 * i + 0) = N_XYZ(0, i);
      DFDd(1, NODDOF_SOH8 * i + 1) = N_XYZ(1, i);
      DFDd(2, NODDOF_SOH8 * i + 2) = N_XYZ(2, i);
      DFDd(3, NODDOF_SOH8 * i + 0) = N_XYZ(1, i);
      DFDd(4, NODDOF_SOH8 * i + 1) = N_XYZ(2, i);
      DFDd(5, NODDOF_SOH8 * i + 0) = N_XYZ(2, i);
      DFDd(6, NODDOF_SOH8 * i + 1) = N_XYZ(0, i);
      DFDd(7, NODDOF_SOH8 * i + 2) = N_XYZ(1, i);
      DFDd(8, NODDOF_SOH8 * i + 2) = N_XYZ(0, i);
    }
  }

  static LINALG::Matrix<9, 1> DsntDF(true);
  static LINALG::Matrix<9, 9> D2sntDF2(true);
  static LINALG::Matrix<9, NUMDIM_SOH8> D2sntDFDn(true);
  static LINALG::Matrix<9, NUMDIM_SOH8> D2sntDFDt(true);

  SolidMaterial()->EvaluateCauchy(
      defgrd, n, t, sigma_nt, DsntDn, DsntDt, &DsntDF, &D2sntDF2, &D2sntDFDn, &D2sntDFDt, 0);

  if (DsntDd)
  {
    DsntDd->Reshape(NUMDOF_SOH8, 1);
    LINALG::Matrix<NUMDOF_SOH8, 1> DsntDd_m(DsntDd->A(), true);
    DsntDd_m.MultiplyTN(1.0, DFDd, DsntDF, 0.0);
  }

  if (D2sntDdDn)
  {
    D2sntDdDn->Reshape(NUMDOF_SOH8, NUMDIM_SOH8);
    LINALG::Matrix<NUMDOF_SOH8, NUMDIM_SOH8> D2sntDdDn_m(D2sntDdDn->A(), true);
    D2sntDdDn_m.MultiplyTN(1.0, DFDd, D2sntDFDn, 0.0);
  }

  if (D2sntDdDt)
  {
    D2sntDdDt->Reshape(NUMDOF_SOH8, NUMDIM_SOH8);
    LINALG::Matrix<NUMDOF_SOH8, NUMDIM_SOH8> D2sntDdDt_m(D2sntDdDt->A(), true);
    D2sntDdDt_m.MultiplyTN(1.0, DFDd, D2sntDFDt, 0.0);
  }

  if (D2sntDd2)
  {
    D2sntDd2->Reshape(NUMDOF_SOH8, NUMDOF_SOH8);
    LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8> D2sntDd2_m(D2sntDd2->A(), true);
    static LINALG::Matrix<9, NUMDOF_SOH8> D2sntDF2DFDd(true);
    D2sntDF2DFDd.Multiply(1.0, D2sntDF2, DFDd, 0.0);
    D2sntDd2_m.MultiplyTN(1.0, DFDd, D2sntDF2DFDd, 0.0);
  }


  if (D2sntDdDxi)
  {
    D2sntDdDxi->Reshape(NUMDOF_SOH8, NUMDIM_SOH8);
    LINALG::Matrix<NUMDOF_SOH8, NUMDIM_SOH8> D2sntDdDxi_m(D2sntDdDxi->A(), true);

    static LINALG::Matrix<DRT::UTILS::DisTypeToNumDeriv2<DRT::Element::hex8>::numderiv2,
        NUMNOD_SOH8>
        deriv2(true);
    deriv2.Clear();
    DRT::UTILS::shape_function_deriv2<DRT::Element::hex8>(xi, deriv2);

    static LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xXF(true);
    static LINALG::Matrix<NUMDIM_SOH8,
        DRT::UTILS::DisTypeToNumDeriv2<DRT::Element::hex8>::numderiv2>
        xXFsec(true);
    xXF.Update(1.0, xcurr, 0.0);
    xXF.MultiplyNT(-1.0, xrefe, defgrd, 1.0);
    xXFsec.MultiplyTT(1.0, xXF, deriv2, 0.0);

    static LINALG::Matrix<9, NUMDIM_SOH8> DFDxi(true);
    DFDxi.Clear();
    for (int a = 0; a < NUMDIM_SOH8; ++a)
      for (int b = 0; b < NUMDIM_SOH8; ++b)
      {
        DFDxi(VoigtMapping::NonSymToVoigt9(a, b), 0) +=
            xXFsec(a, 0) * invJ(b, 0) + xXFsec(a, 3) * invJ(b, 1) + xXFsec(a, 4) * invJ(b, 2);
        DFDxi(VoigtMapping::NonSymToVoigt9(a, b), 1) +=
            xXFsec(a, 3) * invJ(b, 0) + xXFsec(a, 1) * invJ(b, 1) + xXFsec(a, 5) * invJ(b, 2);
        DFDxi(VoigtMapping::NonSymToVoigt9(a, b), 2) +=
            xXFsec(a, 4) * invJ(b, 0) + xXFsec(a, 5) * invJ(b, 1) + xXFsec(a, 2) * invJ(b, 2);
      }

    DsntDxi->MultiplyTN(1.0, DFDxi, DsntDF, 0.0);

    static LINALG::Matrix<DRT::UTILS::DisTypeToNumDeriv2<DRT::Element::hex8>::numderiv2,
        NUMDIM_SOH8>
        Xsec(true);
    static LINALG::Matrix<NUMNOD_SOH8,
        DRT::UTILS::DisTypeToNumDeriv2<DRT::Element::hex8>::numderiv2>
        N_XYZ_Xsec(true);
    Xsec.Multiply(1.0, deriv2, xrefe, 0.0);
    N_XYZ_Xsec.MultiplyTT(1.0, N_XYZ, Xsec, 0.0);

    static LINALG::Matrix<9, NUMDOF_SOH8> D2sntDF2DFDd(true);
    D2sntDF2DFDd.Multiply(1.0, D2sntDF2, DFDd, 0.0);
    D2sntDdDxi_m.MultiplyTN(1.0, D2sntDF2DFDd, DFDxi, 0.0);

    static LINALG::Matrix<9, NUMDIM_SOH8 * NUMDOF_SOH8> D2FDxiDd(true);
    D2FDxiDd.Clear();
    for (int m = 0; m < NUMDIM_SOH8; ++m)
      for (int n = 0; n < NUMDIM_SOH8; ++n)
        for (int k = 0; k < NUMNOD_SOH8; ++k)
        {
          D2FDxiDd(VoigtMapping::NonSymToVoigt9(m, n), NODDOF_SOH8 * (NODDOF_SOH8 * k + m) + 0) +=
              deriv2(0, k) * invJ(n, 0) + deriv2(3, k) * invJ(n, 1) + deriv2(4, k) * invJ(n, 2) -
              N_XYZ_Xsec(k, 0) * invJ(n, 0) - N_XYZ_Xsec(k, 3) * invJ(n, 1) -
              N_XYZ_Xsec(k, 4) * invJ(n, 2);

          D2FDxiDd(VoigtMapping::NonSymToVoigt9(m, n), NODDOF_SOH8 * (NODDOF_SOH8 * k + m) + 1) +=
              deriv2(3, k) * invJ(n, 0) + deriv2(1, k) * invJ(n, 1) + deriv2(5, k) * invJ(n, 2) -
              N_XYZ_Xsec(k, 3) * invJ(n, 0) - N_XYZ_Xsec(k, 1) * invJ(n, 1) -
              N_XYZ_Xsec(k, 5) * invJ(n, 2);

          D2FDxiDd(VoigtMapping::NonSymToVoigt9(m, n), NODDOF_SOH8 * (NODDOF_SOH8 * k + m) + 2) +=
              deriv2(4, k) * invJ(n, 0) + deriv2(5, k) * invJ(n, 1) + deriv2(2, k) * invJ(n, 2) -
              N_XYZ_Xsec(k, 4) * invJ(n, 0) - N_XYZ_Xsec(k, 5) * invJ(n, 1) -
              N_XYZ_Xsec(k, 2) * invJ(n, 2);

          for (int l = 0; l < NUMDIM_SOH8; ++l)
            D2sntDdDxi_m(k * 3 + m, l) += DsntDF(VoigtMapping::NonSymToVoigt9(m, n), 0) *
                                          D2FDxiDd(VoigtMapping::NonSymToVoigt9(m, n),
                                              NODDOF_SOH8 * (NODDOF_SOH8 * k + m) + l);
        }
  }

  return;
}

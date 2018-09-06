/*----------------------------------------------------------------------*/
/*!
\file so_sh8p8_evaluate.cpp
\brief some element evaluate
\level 3
\maintainer Martin Kronbichler

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
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_io/io_gmsh.H"
#include "Epetra_Time.h"
#include "Teuchos_TimeMonitor.hpp"
#include "../drt_mat/material_service.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_mat/neohooke.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_mat/aaaraghavanvorp_damage.H"
#include "../drt_mat/visconeohooke.H"
#include "../drt_mat/viscoanisotropic.H"
#include "../drt_mat/viscoelasthyper.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_mat/micromaterial.H"

#include "../headers/definitions.h"

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            bborn 03/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_sh8p8::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  LINALG::Matrix<NUMDOF_, NUMDOF_> elemat1(elemat1_epetra.A(), true);
  LINALG::Matrix<NUMDOF_, NUMDOF_> elemat2(elemat2_epetra.A(), true);
  LINALG::Matrix<NUMDOF_, 1> elevec1(elevec1_epetra.A(), true);
  LINALG::Matrix<NUMDOF_, 1> elevec2(elevec2_epetra.A(), true);
  // elevec3 is not used anyway

  // start with "none"
  DRT::ELEMENTS::So_hex8::ActionType act = So_hex8::none;

  // get the required action
  const std::string action = params.get<std::string>("action", "none");
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
  else if (action == "calc_stc_matrix")
    act = So_hex8::calc_stc_matrix;
  else if (action == "calc_stc_matrix_inverse")
    act = So_hex8::calc_stc_matrix_inverse;
  else if (action == "calc_struct_predict")
    return 0;
  else if (action == "calc_struct_recover")
    return 0;
  else
  {
    std::cout << action << std::endl;
    dserror("Unknown type of action for So_hex8");
  }
  // what should the element do
  switch (act)
  {
    // linear stiffness
    case calc_struct_linstiff:
    {
      // need zero current displacement and residual forces
      LINALG::Matrix<NUMDISP_, 1> mydisp(true);
      LINALG::Matrix<NUMPRES_, 1> mypres(true);
      LINALG::Matrix<NUMDISP_, 1> mydispi(true);
      LINALG::Matrix<NUMPRES_, 1> mypresi(true);
      LINALG::Matrix<NUMDISP_, NUMDISP_> stiffmatrix(true);
      LINALG::Matrix<NUMDISP_, NUMPRES_> gradmatrix(true);
      LINALG::Matrix<NUMPRES_, NUMPRES_> stabmatrix(true);
      LINALG::Matrix<NUMDISP_, 1> force(true);
      LINALG::Matrix<NUMPRES_, 1> incomp(true);
      if ((stab_ == stab_spatial) or (stab_ == stab_spatialaffine))
      {
        LINALG::Matrix<NUMPRES_, NUMDISP_> dargmatrix;
        ForceStiffMass(lm, mydisp, mypres, mydispi, mypresi, NULL, &stiffmatrix, &gradmatrix, NULL,
            &stabmatrix, &force, &incomp, NULL, NULL, NULL, params, INPAR::STR::stress_none,
            INPAR::STR::strain_none);
        BuildElementMatrix(&elemat1, &stiffmatrix, &gradmatrix, &dargmatrix, &stabmatrix);
      }
      else
      {
        ForceStiffMass(lm, mydisp, mypres, mydispi, mypresi, NULL, &stiffmatrix, NULL, NULL,
            &stabmatrix, &force, &incomp, NULL, NULL, NULL, params, INPAR::STR::stress_none,
            INPAR::STR::strain_none);
        BuildElementMatrix(&elemat1, &stiffmatrix, NULL, NULL, &stabmatrix);
      }
      BuildElementVector(&elevec1, &force, &incomp);
    }
    break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    {
      // need current displacement
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mystat(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mystat, lm);
      LINALG::Matrix<NUMDISP_, 1> mydisp;
      LINALG::Matrix<NUMPRES_, 1> mypres;
      ExtractDispAndPres(mystat, mydisp, mypres);
      // residual displacements
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (res == Teuchos::null) dserror("Didn't get \"residual displacement\"");
      std::vector<double> mystati(lm.size());
      DRT::UTILS::ExtractMyValues(*res, mystati, lm);
      LINALG::Matrix<NUMDISP_, 1> mydispi;
      LINALG::Matrix<NUMPRES_, 1> mypresi;
      ExtractDispAndPres(mystati, mydispi, mypresi);
      // allocate element quantities
      LINALG::Matrix<NUMDISP_, NUMDISP_> stiffmatrix(true);
      LINALG::Matrix<NUMDISP_, NUMPRES_> gradmatrix(true);
      LINALG::Matrix<NUMPRES_, NUMPRES_> stabmatrix(true);
      LINALG::Matrix<NUMDISP_, 1> force(true);
      LINALG::Matrix<NUMPRES_, 1> incomp(true);
      double volume = 0.0;
      if ((stab_ == stab_spatial) or (stab_ == stab_spatialaffine))
      {
        LINALG::Matrix<NUMPRES_, NUMDISP_> dargmatrix;
        ForceStiffMass(lm, mydisp, mypres, mydispi, mypresi, NULL, &stiffmatrix, &gradmatrix,
            &dargmatrix, &stabmatrix, &force, &incomp, NULL, NULL, &volume, params,
            INPAR::STR::stress_none, INPAR::STR::strain_none);
        BuildElementMatrix(&elemat1, &stiffmatrix, &gradmatrix, &dargmatrix, &stabmatrix);
      }
      else
      {
        ForceStiffMass(lm, mydisp, mypres, mydispi, mypresi, NULL, &stiffmatrix, &gradmatrix, NULL,
            &stabmatrix, &force, &incomp, NULL, NULL, &volume, params, INPAR::STR::stress_none,
            INPAR::STR::strain_none);
        BuildElementMatrix(&elemat1, &stiffmatrix, &gradmatrix, NULL, &stabmatrix);
      }
      BuildElementVector(&elevec1, &force, &incomp);
      //      AssembleVolume(params,volume);
      // a desperate call for help
      // GnuplotOut(params,mystat,elevec1,elemat1);
    }
    break;

    // internal force vector only
    case calc_struct_internalforce:
    {
      // need current displacement
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mystat(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mystat, lm);
      LINALG::Matrix<NUMDISP_, 1> mydisp;
      LINALG::Matrix<NUMPRES_, 1> mypres;
      ExtractDispAndPres(mystat, mydisp, mypres);
      // residual displacements
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (res == Teuchos::null) dserror("Didn't get \"residual displacement\"");
      std::vector<double> mystati(lm.size());
      DRT::UTILS::ExtractMyValues(*res, mystati, lm);
      LINALG::Matrix<NUMDISP_, 1> mydispi;
      LINALG::Matrix<NUMPRES_, 1> mypresi;
      ExtractDispAndPres(mystati, mydispi, mypresi);
      // allocate element quantities
      LINALG::Matrix<NUMDISP_, NUMDISP_> stiffmatrix(true);
      LINALG::Matrix<NUMDISP_, NUMPRES_> gradmatrix(true);
      LINALG::Matrix<NUMPRES_, NUMPRES_> stabmatrix(true);
      LINALG::Matrix<NUMDISP_, 1> force(true);
      LINALG::Matrix<NUMPRES_, 1> incomp(true);
      double volume = 0.0;
      ForceStiffMass(lm, mydisp, mypres, mydispi, mypresi, NULL, &stiffmatrix, &gradmatrix, NULL,
          &stabmatrix, &force, &incomp, NULL, NULL, &volume, params, INPAR::STR::stress_none,
          INPAR::STR::strain_none);
      BuildElementVector(&elevec1, &force, &incomp);
      //      AssembleVolume(params,volume);
    }
    break;

    // linear stiffness and consistent mass matrix
    case calc_struct_linstiffmass:
      dserror("Case 'calc_struct_linstiffmass' not yet implemented");
      break;

    // nonlinear stiffness, internal force vector, and consistent/lumped mass matrix
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    {
      // need current displacement
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mystat(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mystat, lm);
      LINALG::Matrix<NUMDISP_, 1> mydisp;
      LINALG::Matrix<NUMPRES_, 1> mypres;
      ExtractDispAndPres(mystat, mydisp, mypres);
      // residual displacements
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (res == Teuchos::null) dserror("Didn't get \"residual displacement\"");
      std::vector<double> mystati(lm.size());
      DRT::UTILS::ExtractMyValues(*res, mystati, lm);
      LINALG::Matrix<NUMDISP_, 1> mydispi;
      LINALG::Matrix<NUMPRES_, 1> mypresi;
      ExtractDispAndPres(mystati, mydispi, mypresi);
      // allocate element quantities
      LINALG::Matrix<NUMDISP_, NUMDISP_> massmatrix(true);
      LINALG::Matrix<NUMDISP_, NUMDISP_> stiffmatrix(true);
      LINALG::Matrix<NUMDISP_, NUMPRES_> gradmatrix(true);
      LINALG::Matrix<NUMPRES_, NUMPRES_> stabmatrix(true);
      LINALG::Matrix<NUMDISP_, 1> force(true);
      LINALG::Matrix<NUMPRES_, 1> incomp(true);
      double volume = 0.0;
      if ((stab_ == stab_spatial) or (stab_ == stab_spatialaffine))
      {
        LINALG::Matrix<NUMPRES_, NUMDISP_> dargmatrix;
        ForceStiffMass(lm, mydisp, mypres, mydispi, mypresi, &massmatrix, &stiffmatrix, &gradmatrix,
            &dargmatrix, &stabmatrix, &force, &incomp, NULL, NULL, &volume, params,
            INPAR::STR::stress_none, INPAR::STR::strain_none);
        BuildElementMatrix(&elemat1, &stiffmatrix, &gradmatrix, &dargmatrix, &stabmatrix);
      }
      else
      {
        ForceStiffMass(lm, mydisp, mypres, mydispi, mypresi, &massmatrix, &stiffmatrix, &gradmatrix,
            NULL, &stabmatrix, &force, &incomp, NULL, NULL, &volume, params,
            INPAR::STR::stress_none, INPAR::STR::strain_none);
        BuildElementMatrix(&elemat1, &stiffmatrix, &gradmatrix, NULL, &stabmatrix);
      }
      // lump mass
      if (act == calc_struct_nlnstifflmass) soh8_lumpmass(&massmatrix);
      // assemble displacement pressure parts
      BuildElementMatrix(&elemat2, &massmatrix, NULL, NULL, NULL);
      BuildElementVector(&elevec1, &force, &incomp);
      //      AssembleVolume(params,volume);
    }
    break;

    // evaluate stresses and strains at gauss points
    case calc_struct_stress:
    {
      // nothing to do for ghost elements
      if (discretization.Comm().MyPID() == Owner())
      {
        // data
        Teuchos::RCP<std::vector<char>> stressdata =
            params.get<Teuchos::RCP<std::vector<char>>>("stress", Teuchos::null);
        Teuchos::RCP<std::vector<char>> straindata =
            params.get<Teuchos::RCP<std::vector<char>>>("strain", Teuchos::null);
        if (stressdata == Teuchos::null) dserror("Cannot get stress 'data'");
        if (straindata == Teuchos::null) dserror("Cannot get strain 'data'");
        // current displacements
        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
        std::vector<double> mystat(lm.size());
        DRT::UTILS::ExtractMyValues(*disp, mystat, lm);
        LINALG::Matrix<NUMDISP_, 1> mydisp;
        LINALG::Matrix<NUMPRES_, 1> mypres;
        ExtractDispAndPres(mystat, mydisp, mypres);
        // residual displacements
        Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
        if (res == Teuchos::null) dserror("Didn't get \"residual displacement\"");
        std::vector<double> mystati(lm.size());
        DRT::UTILS::ExtractMyValues(*res, mystati, lm);
        LINALG::Matrix<NUMDISP_, 1> mydispi;
        LINALG::Matrix<NUMPRES_, 1> mypresi;
        ExtractDispAndPres(mystati, mydispi, mypresi);
        // types
        LINALG::Matrix<NUMGPT_, MAT::NUM_STRESS_3D> stress;
        LINALG::Matrix<NUMGPT_, MAT::NUM_STRESS_3D> strain;
        INPAR::STR::StressType iostress =
            DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
        INPAR::STR::StrainType iostrain =
            DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);
        ForceStiffMass(lm, mydisp, mypres, mydispi, mypresi, NULL, NULL, NULL, NULL, NULL, NULL,
            NULL, &stress, &strain, NULL, params, iostress, iostrain);
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
      LINALG::Matrix<NUMGPT_, MAT::NUM_STRESS_3D> gpstress(((*gpstressmap)[gid])->A(), true);

      Teuchos::RCP<Epetra_MultiVector> poststress =
          params.get<Teuchos::RCP<Epetra_MultiVector>>("poststress", Teuchos::null);
      if (poststress == Teuchos::null) dserror("No element stress/strain vector available");

      if (stresstype == "ndxyz")
      {
        // extrapolate stresses/strains at Gauss points to nodes
        sosh8p8_expol(gpstress, *poststress);
      }
      else if (stresstype == "cxyz")
      {
        const Epetra_BlockMap elemap = poststress->Map();
        int lid = elemap.LID(Id());
        if (lid != -1)
        {
          for (int i = 0; i < MAT::NUM_STRESS_3D; ++i)
          {
            double& s = (*((*poststress)(i)))[lid];  // resolve pointer for faster access
            s = 0.;
            for (int j = 0; j < NUMGPT_; ++j)
            {
              s += gpstress(j, i);
            }
            s *= 1.0 / NUMGPT_;
          }
        }
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
      // do something with internal EAS, etc parameters
      if (eastype_ != soh8_easnone)
      {
        Epetra_SerialDenseMatrix* alpha =
            data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // Alpha_{n+1}
        Epetra_SerialDenseMatrix* alphao =
            data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");  // Alpha_n
        // alphao := alpha
        if (eastype_ == soh8_eassosh8)
        {
          LINALG::DENSEFUNCTIONS::update<double, NUMEAS_SOSH8_, 1>(*alphao, *alpha);
        }
        else if (eastype_ == soh8_easa)
        {
          LINALG::DENSEFUNCTIONS::update<double, NUMEAS_A_, 1>(*alphao, *alpha);
        }
        else
        {
          dserror("Not impl.");
        }
      }
      // Update of history for materials
      SolidMaterial()->Update();
    }
    break;

    case calc_struct_reset_istep:
    {
      // do something with internal EAS, etc parameters
      if (eastype_ != soh8_easnone)
      {
        Epetra_SerialDenseMatrix* alpha =
            data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // Alpha_{n+1}
        Epetra_SerialDenseMatrix* alphao =
            data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");  // Alpha_n
        // alpha := alphao
        if (eastype_ == soh8_eassosh8)
        {
          LINALG::DENSEFUNCTIONS::update<double, NUMEAS_SOSH8_, 1>(*alpha, *alphao);
        }
        else if (eastype_ == soh8_easa)
        {
          LINALG::DENSEFUNCTIONS::update<double, NUMEAS_A_, 1>(*alpha, *alphao);
        }
        else
        {
          dserror("Not impl.");
        }
      }
      // Reset of history (if needed)
      SolidMaterial()->ResetStep();
    }
    break;

    case multi_calc_dens:
    {
      soh8_homog(params);
    }
    break;

    // in case of multi-scale problems, possible EAS internal data on microscale
    // have to be stored in every macroscopic Gauss point
    // allocation and initializiation of these data arrays can only be
    // done in the elements that know the number of EAS parameters
    case multi_eas_init:
    {
      dserror("Meaningful?");
      if (eastype_ != soh8_easnone)
      {
        soh8_eas_init_multi(params);
      }
    }
    break;

    // in case of multi-scale problems, possible EAS internal data on microscale
    // have to be stored in every macroscopic Gauss point
    // before any microscale simulation, EAS internal data has to be
    // set accordingly
    case multi_eas_set:
    {
      dserror("Meaningful?");
      if (eastype_ != soh8_easnone)
      {
        soh8_set_eas_multi(params);
      }
    }
    break;

    // read restart of microscale
    case multi_readrestart:
    {
      Teuchos::RCP<MAT::Material> mat = Material();

      if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale) soh8_read_restart_multi();
    }
    break;

    case calc_stc_matrix:
    {
      const INPAR::STR::STC_Scale stc_scaling =
          DRT::INPUT::get<INPAR::STR::STC_Scale>(params, "stc_scaling");
      if (stc_scaling == INPAR::STR::stc_none)
        dserror("To scale or not to scale, that's the querry!");
      else
      {
        CalcSTCMatrix(
            elemat1, stc_scaling, params.get<int>("stc_layer"), lm, discretization, false);
      }
    }
    break;
    case calc_stc_matrix_inverse:
    {
      const INPAR::STR::STC_Scale stc_scaling =
          DRT::INPUT::get<INPAR::STR::STC_Scale>(params, "stc_scaling");
      if (stc_scaling == INPAR::STR::stc_none)
        dserror("To scale or not to scale, that's the query!");
      else
      {
        CalcSTCMatrix(elemat1, stc_scaling, params.get<int>("stc_layer"), lm, discretization, true);
      }
    }
    break;
    default:
    {
      dserror("Unknown type of action for So_sh8p8");
      break;
    }
  }
  return 0;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_sh8p8::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  // build 24x1 location vector
  std::vector<int> lm_d(NUMDISP_);
  for (int inode = 0; inode < NUMNOD_; ++inode)
    for (int idof = 0; idof < NODDISP_; ++idof)
      lm_d[inode * NODDISP_ + idof] = lm[inode * NODDOF_ + idof];

  // 24x1 load vector
  Epetra_SerialDenseVector elevec1_d(NUMDISP_);
  // 24x24 tangent
  Epetra_SerialDenseMatrix* elemat1_dd = NULL;
  if (elemat1 != NULL) elemat1_dd = new Epetra_SerialDenseMatrix(NUMDISP_, NUMDISP_);
  // determine (displacement) load vector (and tangent)
  const int rv =
      So_hex8::EvaluateNeumann(params, discretization, condition, lm_d, elevec1_d, elemat1_dd);
  // copy vector
  LINALG::Matrix<NUMDOF_, 1> ev1(elevec1, true);       // only view
  LINALG::Matrix<NUMDISP_, 1> ev1_d(elevec1_d, true);  // only view
  BuildElementVector(&ev1, &ev1_d, NULL);
  // copy matrix
  if ((elemat1_dd != NULL) and (elemat1 != NULL))
  {
    LINALG::Matrix<NUMDOF_, NUMDOF_> em1(*elemat1, true);          // only view
    LINALG::Matrix<NUMDISP_, NUMDISP_> em1_dd(*elemat1_dd, true);  // only view
    BuildElementMatrix(&em1, &em1_dd, NULL, NULL, NULL);
    delete elemat1_dd;
  }

  return rv;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::ForceStiffMass(const std::vector<int>& lm,  // location matrix
    const LINALG::Matrix<NUMDISP_, 1>& disp,                              // current displacements
    const LINALG::Matrix<NUMPRES_, 1>& pres,                              // current pressures
    const LINALG::Matrix<NUMDISP_, 1>& dispi,                // current residual displacements
    const LINALG::Matrix<NUMPRES_, 1>& presi,                // current residual pressures
    LINALG::Matrix<NUMDISP_, NUMDISP_>* massmatrix,          // element mass matrix
    LINALG::Matrix<NUMDISP_, NUMDISP_>* stiffmatrix,         // element stiffness matrix
    LINALG::Matrix<NUMDISP_, NUMPRES_>* gradmatrix,          // element gradient matrix
    LINALG::Matrix<NUMPRES_, NUMDISP_>* dargmatrix,          // 'transposed' element gradient matrix
    LINALG::Matrix<NUMPRES_, NUMPRES_>* stabmatrix,          // element stabilisation matrix
    LINALG::Matrix<NUMDISP_, 1>* force,                      // element internal force vector
    LINALG::Matrix<NUMPRES_, 1>* incomp,                     // incompressibility residual
    LINALG::Matrix<NUMGPT_, MAT::NUM_STRESS_3D>* elestress,  // stresses at GP
    LINALG::Matrix<NUMGPT_, MAT::NUM_STRESS_3D>* elestrain,  // strains at GP
    double* volume,                                          // element volume
    Teuchos::ParameterList& params,                          // algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,                   // stress output option
    const INPAR::STR::StrainType iostrain                    // strain output option
)
{
  /* ============================================================================*
  ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
  ** ============================================================================*/
  const static std::vector<LINALG::Matrix<NUMNOD_, 1>> shapefcts = soh8_shapefcts();
  const static std::vector<LINALG::Matrix<NUMDIM_, NUMNOD_>> derivs = soh8_derivs();
  const static std::vector<double> gpweights = soh8_weights();
  /* ============================================================================*/

  // update element geometry
  LINALG::Matrix<NUMNOD_, NUMDIM_> xrefe;  // material coord. of element
  LINALG::Matrix<NUMNOD_, NUMDIM_> xcurr;  // current  coord. of element
  DRT::Node** nodes = Nodes();
  for (int i = 0; i < NUMNOD_; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];

    xcurr(i, 0) = xrefe(i, 0) + disp(i * NODDISP_ + 0, 0);
    xcurr(i, 1) = xrefe(i, 1) + disp(i * NODDISP_ + 1, 0);
    xcurr(i, 2) = xrefe(i, 2) + disp(i * NODDISP_ + 2, 0);
  }

  // EAS Technology: declare, intialize, set up, and alpha history
  //
  // the current alphas are (re-)evaluated out of
  // Kaa and Kda, Kpa of previous step to avoid additional element call.
  // This corresponds to the (innermost) element update loop
  // in the nonlinear FE-Skript page 120 (load-control alg. with EAS)
  //
  // in any case declare variables, sizes etc. only in EAS case
  Epetra_SerialDenseMatrix* oldfeas = NULL;                     // EAS history
  Epetra_SerialDenseMatrix* oldKaainv = NULL;                   // EAS history
  Epetra_SerialDenseMatrix* oldKad = NULL;                      // EAS history
  Epetra_SerialDenseMatrix* oldKap = NULL;                      // EAS history
  Teuchos::RCP<Epetra_SerialDenseVector> feas = Teuchos::null;  // EAS portion of internal forces
  Teuchos::RCP<Epetra_SerialDenseMatrix> Kaa = Teuchos::null;   // EAS matrix Kaa
  Teuchos::RCP<Epetra_SerialDenseMatrix> Kap = Teuchos::null;   // EAS matrix Kpa
  Teuchos::RCP<Epetra_SerialDenseMatrix> Kad = Teuchos::null;   // EAS matrix Kda
  Epetra_SerialDenseMatrix* alpha = NULL;                       // EAS alphas
  Teuchos::RCP<Epetra_SerialDenseMatrix> M =
      Teuchos::null;  // EAS matrix M at current GP, fixed for sosh8
  // make update
  if (eastype_ == soh8_easnone)
  {
    ;  // continue
  }
  else if (eastype_ == soh8_eassosh8)
  {
    EasUpdateIncrementally<NUMEAS_SOSH8_>(
        oldfeas, oldKaainv, oldKad, oldKap, feas, Kaa, Kad, Kap, alpha, M, data_, dispi, presi);
  }
  else if (eastype_ == soh8_easa)
  {
    EasUpdateIncrementally<NUMEAS_A_>(
        oldfeas, oldKaainv, oldKad, oldKap, feas, Kaa, Kad, Kap, alpha, M, data_, dispi, presi);
  }
  else
  {
    dserror("EAS ain't easy");
  }

  // evaluation of EAS variables (which are constant for the following):
  // -> M defining interpolation of enhanced strains alpha, evaluated at GPs
  // -> determinant of Jacobi matrix at element origin (r=s=t=0.0)
  // -> T0^{-T}
  std::vector<Epetra_SerialDenseMatrix>* M_GP = NULL;             // EAS matrix M at all GPs
  double detJ0;                                                   // detJ(origin)
  LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> T0invT;  // trafo matrix
  if (eastype_ != soh8_easnone)
  {
    soh8_eassetup(&M_GP, detJ0, T0invT, xrefe);
  }

  //
  // ANS Element technology to remedy
  //  - transverse-shear locking E_rt and E_st
  //  - trapezoidal (curvature-thickness) locking E_tt
  //

  // modified B-operator in local(parameter) element space

  // ANS modified rows of bop in local(parameter) coords
  LINALG::Matrix<NUMANS_ * NUMSP_, NUMDISP_> B_ans_loc;
  // Jacobian evaluated at all ANS sampling points
  std::vector<LINALG::Matrix<NUMDIM_, NUMDIM_>> jac_sps(NUMSP_);
  // CURRENT Jacobian evaluated at all ANS sampling points
  std::vector<LINALG::Matrix<NUMDIM_, NUMDIM_>> jac_cur_sps(NUMSP_);
  // pointer to derivs evaluated at all sampling points
  std::vector<LINALG::Matrix<NUMDIM_, NUMNOD_>>* deriv_sp = NULL;
  // evaluate all necessary variables for ANS
  sosh8_anssetup(xrefe, xcurr, &deriv_sp, jac_sps, jac_cur_sps, B_ans_loc);
  // (r,s) gp-locations of fully integrated linear 8-node Hex
  // necessary for ANS interpolation
  const double gploc = 1.0 / sqrt(3.0);  // gp sampling point value for linear fct
  const double r[NUMGPT_] = {-gploc, gploc, gploc, -gploc, -gploc, gploc, gploc, -gploc};
  const double s[NUMGPT_] = {-gploc, -gploc, gploc, gploc, -gploc, -gploc, gploc, gploc};
  // const double t[NUMGPT_] = {-gploc,-gploc,-gploc,-gploc, gploc, gploc, gploc, gploc};

  // proportions of element
  Teuchos::RCP<LINALG::Matrix<NUMDIM_, NUMDIM_>> jac0 = Teuchos::null;  // Jacobian at origin
  Teuchos::RCP<LINALG::Matrix<NUMDIM_, 1>> axmetr0 = Teuchos::null;     // axial metrics at origin
  Teuchos::RCP<const double> hths = Teuchos::null;
  Teuchos::RCP<const double> hthr = Teuchos::null;
#if 0
  std::vector<LINALG::Matrix<NUMDIM_,NUMDIM_> > metr_sps;
#endif
  if (ans_ == ans_onspot)
  {
    jac0 = Teuchos::rcp(new LINALG::Matrix<NUMDIM_, NUMDIM_>(false));
    axmetr0 = Teuchos::rcp(new LINALG::Matrix<NUMDIM_, 1>(false));
    AxialMetricsAtOrigin(xrefe, *jac0, *axmetr0);
    hths = Teuchos::rcp(new double((*axmetr0)(2) / (*axmetr0)(1)));
    hthr = Teuchos::rcp(new double((*axmetr0)(2) / (*axmetr0)(0)));
//    std::cout << "hths=" << *hths << ", hthr=" << *hthr << std::endl;
#if 0
    metr_sps.resize(NUMSP_);
    for (int sp=0; sp<NUMSP_; ++sp)
      LocalMetrics(jac_sps[sp],metr_sps[sp]);
#endif
  }


  // ---------------------------------------------------------------------
  // first loop over Gauss point
  // stabilisation matrices
  LINALG::Matrix<NUMPRESBRO_, NUMPRESBRO_> stabAA(true);  // element volume
  LINALG::Matrix<NUMPRES_, NUMPRESBRO_> stabHA(
      true);  // integral of pressure shape functions across element
  LINALG::Matrix<NUMPRES_, NUMPRES_> stabHH(true);  // mass-like matrix
  std::vector<LINALG::Matrix<NUMPRESBRO_, NUMPRESBRO_>> stabA(
      NUMGPT_);  // shape functions for projected Q0/constant pressure
  for (int gp = 0; gp < NUMGPT_; ++gp)
  {
    // integration weight
    double wdetJ = 0.0;
    if ((stab_ == stab_spatial) or (stab_ == stab_spatialaffine))
    {
      // (transposed) spatial-to-parametric Jacobian j = (x_{,xi})^T
      LINALG::Matrix<NUMDIM_, NUMDIM_> jac;
      jac.Multiply(derivs[gp], xcurr);
      const double detj = jac.Determinant();
      wdetJ = detj * gpweights[gp];
    }
    else
    {
      // (transposed) material-to-parametric Jacobian J = (X_{,xi})^T
      LINALG::Matrix<NUMDIM_, NUMDIM_> Jac;
      Jac.Multiply(derivs[gp], xrefe);
      const double detJ = Jac.Determinant();
      wdetJ = detJ * gpweights[gp];
    }

    stabA[gp].PutScalar(1.0);

    stabAA.MultiplyTN(wdetJ, stabA[gp], stabA[gp], 1.0);
    stabHA.MultiplyNT(wdetJ, shapefcts[gp], stabA[gp], 1.0);
    stabHH.MultiplyNT(wdetJ, shapefcts[gp], shapefcts[gp], 1.0);
  }

  // stabilisation matrix
  if (stabmatrix != NULL)
  {
    // shear modulus
    const double shearmod = ShearMod();
    // (-Cem) = -1./shearmod*( Mem - Eem'*inv(Dem)*Eem );
    stabmatrix->Update(stabHH);
    stabmatrix->MultiplyNT(-1.0 / stabAA(0, 0), stabHA, stabHA, 1.0);
    stabmatrix->Scale(-1.0 / shearmod);
  }

  // extra matrices to deal with linearisation of stabilisation
  // matrix in case of spatial version
  Teuchos::RCP<LINALG::Matrix<NUMPRES_, NUMDISP_>> stabHHbydisp = Teuchos::null;
  Teuchos::RCP<LINALG::Matrix<NUMPRES_ * NUMPRESBRO_, NUMDISP_>> stabHAbydisp = Teuchos::null;
  Teuchos::RCP<LINALG::Matrix<NUMPRESBRO_ * NUMPRESBRO_, NUMDISP_>> stabAAbydisp = Teuchos::null;
  double voldev = 0.0;  // volumetric deviation of element
  Teuchos::RCP<LINALG::Matrix<1, NUMDISP_>> detdefgradbydisp = Teuchos::null;
  if ((stab_ == stab_spatial) or (stab_ == stab_spatialaffine))
  {
    stabHHbydisp = Teuchos::rcp(new LINALG::Matrix<NUMNOD_, NUMDISP_>(true));
    stabHAbydisp = Teuchos::rcp(new LINALG::Matrix<NUMNOD_ * NUMPRESBRO_, NUMDISP_>(true));
    stabAAbydisp = Teuchos::rcp(new LINALG::Matrix<NUMPRESBRO_ * NUMPRESBRO_, NUMDISP_>(true));
    detdefgradbydisp = Teuchos::rcp(new LINALG::Matrix<1, NUMDISP_>(true));
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp = 0; gp < NUMGPT_; ++gp)
  {
    /* compute the Jacobian matrix which looks like:
    **         [ x_,r  y_,r  z_,r ]
    **     J = [ x_,s  y_,s  z_,s ]
    **         [ x_,t  y_,t  z_,t ]
    */
    LINALG::Matrix<NUMDIM_, NUMDIM_> jac;
    jac.Multiply(derivs[gp], xrefe);

    // compute determinant of Jacobian by Sarrus' rule
    double detJ = jac.Determinant();
    if (fabs(detJ) <= EPS10)
      dserror("JACOBIAN DETERMINANT CLOSE TO ZERO");
    else if (detJ < 0.0)
      dserror("NEGATIVE JACOBIAN DETERMINANT");

    // integration factor
    const double detJ_w = detJ * gpweights[gp];

    /* compute the CURRENT Jacobian matrix which looks like:
    **         [ xcurr_,r  ycurr_,r  zcurr_,r ]
    **  Jcur = [ xcurr_,s  ycurr_,s  zcurr_,s ]
    **         [ xcurr_,t  ycurr_,t  zcurr_,t ]
    ** Used to transform the global displacements into parametric space
    */
    LINALG::Matrix<NUMDIM_, NUMDIM_> jac_cur;
    jac_cur.Multiply(derivs[gp], xcurr);

    // compute determinant of Jacobian by Sarrus' rule
    double detJ_cur = jac_cur.Determinant();
    if (detJ_cur == 0.0)
      dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ_cur < 0.0)
      dserror("NEGATIVE JACOBIAN DETERMINANT");

    // set up B-Operator in local(parameter) element space including ANS
    LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDISP_> bop_loc;
    for (int inode = 0; inode < NUMNOD_; ++inode)
    {
      for (int dim = 0; dim < NUMDIM_; ++dim)
      {
        // B_loc_rr = N_r.X_r
        bop_loc(0, inode * 3 + dim) = derivs[gp](0, inode) * jac_cur(0, dim);
        // B_loc_ss = N_s.X_s
        bop_loc(1, inode * 3 + dim) = derivs[gp](1, inode) * jac_cur(1, dim);
        // B_loc_rs = N_r.X_s + N_s.X_r
        bop_loc(3, inode * 3 + dim) =
            derivs[gp](0, inode) * jac_cur(1, dim) + derivs[gp](1, inode) * jac_cur(0, dim);
        if (ans_ == ans_none)
        {
          // B_loc_tt = N_t.X_t
          bop_loc(2, inode * 3 + dim) = derivs[gp](2, inode) * jac_cur(2, dim);
          // B_loc_st = N_s.X_t + N_t.X_s
          bop_loc(4, inode * 3 + dim) =
              derivs[gp](1, inode) * jac_cur(2, dim) + derivs[gp](2, inode) * jac_cur(1, dim);
          // B_loc_rt = N_r.X_t + N_t.X_r
          bop_loc(5, inode * 3 + dim) =
              derivs[gp](0, inode) * jac_cur(2, dim) + derivs[gp](2, inode) * jac_cur(0, dim);
        }
        else
        {
          // B_loc_tt = interpolation along (r x s) of ANS B_loc_tt
          //          = (1-r)(1-s)/4 * B_ans(SP E) + (1+r)(1-s)/4 * B_ans(SP F)
          //           +(1+r)(1+s)/4 * B_ans(SP G) + (1-r)(1+s)/4 * B_ans(SP H)
          bop_loc(2, inode * 3 + dim) =
              0.25 * (1 - r[gp]) * (1 - s[gp]) * B_ans_loc(0 + 4 * NUMANS_, inode * 3 + dim) +
              0.25 * (1 + r[gp]) * (1 - s[gp]) * B_ans_loc(0 + 5 * NUMANS_, inode * 3 + dim) +
              0.25 * (1 + r[gp]) * (1 + s[gp]) * B_ans_loc(0 + 6 * NUMANS_, inode * 3 + dim) +
              0.25 * (1 - r[gp]) * (1 + s[gp]) * B_ans_loc(0 + 7 * NUMANS_, inode * 3 + dim);
          // B_loc_st and B_loc_rt
          if (ans_ == ans_lateral)
          {
            // B_loc_st = interpolation along r of ANS B_loc_st
            //          = (1+r)/2 * B_ans(SP B) + (1-r)/2 * B_ans(SP D)
            bop_loc(4, inode * 3 + dim) =
                0.5 * (1.0 + r[gp]) * B_ans_loc(1 + 1 * NUMANS_, inode * 3 + dim) +
                0.5 * (1.0 - r[gp]) * B_ans_loc(1 + 3 * NUMANS_, inode * 3 + dim);
            // B_loc_rt = interpolation along s of ANS B_loc_rt
            //          = (1-s)/2 * B_ans(SP A) + (1+s)/2 * B_ans(SP C)
            bop_loc(5, inode * 3 + dim) =
                0.5 * (1.0 - s[gp]) * B_ans_loc(2 + 0 * NUMANS_, inode * 3 + dim) +
                0.5 * (1.0 + s[gp]) * B_ans_loc(2 + 2 * NUMANS_, inode * 3 + dim);
          }
          else if (ans_ == ans_onspot)
          {
            // B_loc_st = interpolation along r of ANS B_loc_st
            //          = (1+r)/2 * B_ans(SP B) + (1-r)/2 * B_ans(SP D)
            bop_loc(4, inode * 3 + dim) = 0.5 * (1.0 + r[gp]) * 0.5 * (1 - s[gp]) *
                                              B_ans_loc(1 + 1 * NUMANS_, inode * 3 + dim)  // B
                                          + 0.5 * (1.0 + r[gp]) * 0.5 * (1 - s[gp]) * (*hths) *
                                                B_ans_loc(1 + 5 * NUMANS_, inode * 3 + dim)  // F
                                          + 0.5 * (1.0 + r[gp]) * 0.5 * (1 + s[gp]) *
                                                B_ans_loc(1 + 1 * NUMANS_, inode * 3 + dim)  // B
                                          + 0.5 * (1.0 + r[gp]) * 0.5 * (1 + s[gp]) * (*hths) *
                                                B_ans_loc(1 + 6 * NUMANS_, inode * 3 + dim)  // G
                                          + 0.5 * (1.0 - r[gp]) * 0.5 * (1 - s[gp]) *
                                                B_ans_loc(1 + 3 * NUMANS_, inode * 3 + dim)  // D
                                          + 0.5 * (1.0 - r[gp]) * 0.5 * (1 - s[gp]) * (*hths) *
                                                B_ans_loc(1 + 4 * NUMANS_, inode * 3 + dim)  // E
                                          + 0.5 * (1.0 - r[gp]) * 0.5 * (1 + s[gp]) *
                                                B_ans_loc(1 + 3 * NUMANS_, inode * 3 + dim)  // D
                                          + 0.5 * (1.0 - r[gp]) * 0.5 * (1 + s[gp]) * (*hths) *
                                                B_ans_loc(1 + 7 * NUMANS_, inode * 3 + dim);  // H
            // B_loc_rt = interpolation along s of ANS B_loc_rt
            //          = (1-s)/2 * B_ans(SP A) + (1+s)/2 * B_ans(SP C)
            bop_loc(5, inode * 3 + dim) = 0.5 * (1.0 - s[gp]) * 0.5 * (1 - r[gp]) *
                                              B_ans_loc(2 + 0 * NUMANS_, inode * 3 + dim)  // A
                                          + 0.5 * (1.0 - s[gp]) * 0.5 * (1 - r[gp]) * (*hthr) *
                                                B_ans_loc(2 + 4 * NUMANS_, inode * 3 + dim)  // E
                                          + 0.5 * (1.0 - s[gp]) * 0.5 * (1 + r[gp]) *
                                                B_ans_loc(2 + 0 * NUMANS_, inode * 3 + dim)  // A
                                          + 0.5 * (1.0 - s[gp]) * 0.5 * (1 + r[gp]) * (*hthr) *
                                                B_ans_loc(2 + 5 * NUMANS_, inode * 3 + dim)  // F
                                          + 0.5 * (1.0 + s[gp]) * 0.5 * (1 - r[gp]) *
                                                B_ans_loc(2 + 2 * NUMANS_, inode * 3 + dim)  // C
                                          + 0.5 * (1.0 + s[gp]) * 0.5 * (1 - r[gp]) * (*hthr) *
                                                B_ans_loc(2 + 7 * NUMANS_, inode * 3 + dim)  // H
                                          + 0.5 * (1.0 + s[gp]) * 0.5 * (1 + r[gp]) *
                                                B_ans_loc(2 + 2 * NUMANS_, inode * 3 + dim)  // C
                                          + 0.5 * (1.0 + s[gp]) * 0.5 * (1 + r[gp]) * (*hthr) *
                                                B_ans_loc(2 + 6 * NUMANS_, inode * 3 + dim);  // G
          }
          else
          {
            dserror("You should not turn arb here.");
          }
        }
      }
    }

    // transformation from local (parameter) element space to global(material) space
    // with famous 'T'-matrix already used for EAS but now evaluated at each gp
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> TinvT;
    sosh8_evaluateT(jac, TinvT);
    LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDISP_> bop;
    bop.Multiply(TinvT, bop_loc);

    // local GL strain vector lstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    // but with modified ANS strains E33, E23 and E13
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> lstrain;
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
    // remaining natural strains
    if (ans_ == ans_none)
    {
      // Ett = 0.5 * (dz/dt * dz/dt^T - dZ/dt * dZ/dt^T)
      lstrain(2) =
          0.5 * (+(jac_cur(2, 0) * jac_cur(2, 0) + jac_cur(2, 1) * jac_cur(2, 1) +
                     jac_cur(2, 2) * jac_cur(2, 2)) -
                    (jac(2, 0) * jac(2, 0) + jac(2, 1) * jac(2, 1) + jac(2, 2) * jac(2, 2)));
      // Est = (dx/dt * dy/ds^T - dX/dt * dY/ds^T)
      lstrain(4) = (+(jac_cur(1, 0) * jac_cur(2, 0) + jac_cur(1, 1) * jac_cur(2, 1) +
                        jac_cur(1, 2) * jac_cur(2, 2)) -
                    (jac(1, 0) * jac(2, 0) + jac(1, 1) * jac(2, 1) + jac(1, 2) * jac(2, 2)));
      // Etr = (dx/dr * dy/dt^T - dX/dr * dY/dt^T)
      lstrain(5) = (+(jac_cur(2, 0) * jac_cur(0, 0) + jac_cur(2, 1) * jac_cur(0, 1) +
                        jac_cur(2, 2) * jac_cur(0, 2)) -
                    (jac(2, 0) * jac(0, 0) + jac(2, 1) * jac(0, 1) + jac(2, 2) * jac(0, 2)));
    }
    // ANS modification of strains ************************************** ANS
    else
    {
      double dxdt_A = 0.0;
      double dXdt_A = 0.0;
      double dydt_B = 0.0;
      double dYdt_B = 0.0;
      double dxdt_C = 0.0;
      double dXdt_C = 0.0;
      double dydt_D = 0.0;
      double dYdt_D = 0.0;

      double dxdt_E = 0.0;
      double dXdt_E = 0.0;
      double dydt_E = 0.0;
      double dYdt_E = 0.0;
      double dxdt_F = 0.0;
      double dXdt_F = 0.0;
      double dydt_F = 0.0;
      double dYdt_F = 0.0;
      double dxdt_G = 0.0;
      double dXdt_G = 0.0;
      double dydt_G = 0.0;
      double dYdt_G = 0.0;
      double dxdt_H = 0.0;
      double dXdt_H = 0.0;
      double dydt_H = 0.0;
      double dYdt_H = 0.0;

      double dzdt_E = 0.0;
      double dZdt_E = 0.0;
      double dzdt_F = 0.0;
      double dZdt_F = 0.0;
      double dzdt_G = 0.0;
      double dZdt_G = 0.0;
      double dzdt_H = 0.0;
      double dZdt_H = 0.0;

      // vector product of rows of jacobians at corresponding sampling point    std::cout <<
      // jac_cur_sps;
      for (int dim = 0; dim < NUMDIM_; ++dim)
      {
        dxdt_A += jac_cur_sps[0](0, dim) * jac_cur_sps[0](2, dim);  // g_13^A
        dXdt_A += jac_sps[0](0, dim) * jac_sps[0](2, dim);          // G_13^A
        dydt_B += jac_cur_sps[1](1, dim) * jac_cur_sps[1](2, dim);  // g_23^B
        dYdt_B += jac_sps[1](1, dim) * jac_sps[1](2, dim);          // G_23^B
        dxdt_C += jac_cur_sps[2](0, dim) * jac_cur_sps[2](2, dim);  // g_13^C
        dXdt_C += jac_sps[2](0, dim) * jac_sps[2](2, dim);          // G_13^C
        dydt_D += jac_cur_sps[3](1, dim) * jac_cur_sps[3](2, dim);  // g_23^D
        dYdt_D += jac_sps[3](1, dim) * jac_sps[3](2, dim);          // G_23^D

        if (ans_ == ans_onspot)
        {
          dxdt_E += jac_cur_sps[4](0, dim) * jac_cur_sps[4](2, dim);  // g_13^E
          dXdt_E += jac_sps[4](0, dim) * jac_sps[4](2, dim);          // G_13^E
          dydt_E += jac_cur_sps[4](1, dim) * jac_cur_sps[4](2, dim);  // g_23^E
          dYdt_E += jac_sps[4](1, dim) * jac_sps[4](2, dim);          // G_23^E
          dxdt_F += jac_cur_sps[5](0, dim) * jac_cur_sps[5](2, dim);  // g_13^F
          dXdt_F += jac_sps[5](0, dim) * jac_sps[5](2, dim);          // G_13^F
          dydt_F += jac_cur_sps[5](1, dim) * jac_cur_sps[5](2, dim);  // g_23^F
          dYdt_F += jac_sps[5](1, dim) * jac_sps[5](2, dim);          // G_23^F
          dxdt_G += jac_cur_sps[6](0, dim) * jac_cur_sps[6](2, dim);  // g_13^G
          dXdt_G += jac_sps[6](0, dim) * jac_sps[6](2, dim);          // G_13^G
          dydt_G += jac_cur_sps[6](1, dim) * jac_cur_sps[6](2, dim);  // g_23^G
          dYdt_G += jac_sps[6](1, dim) * jac_sps[6](2, dim);          // G_23^G
          dxdt_H += jac_cur_sps[7](0, dim) * jac_cur_sps[7](2, dim);  // g_13^H
          dXdt_H += jac_sps[7](0, dim) * jac_sps[7](2, dim);          // G_13^H
          dydt_H += jac_cur_sps[7](1, dim) * jac_cur_sps[7](2, dim);  // g_23^H
          dYdt_H += jac_sps[7](1, dim) * jac_sps[7](2, dim);          // G_23^H
        }

        dzdt_E += jac_cur_sps[4](2, dim) * jac_cur_sps[4](2, dim);
        dZdt_E += jac_sps[4](2, dim) * jac_sps[4](2, dim);
        dzdt_F += jac_cur_sps[5](2, dim) * jac_cur_sps[5](2, dim);
        dZdt_F += jac_sps[5](2, dim) * jac_sps[5](2, dim);
        dzdt_G += jac_cur_sps[6](2, dim) * jac_cur_sps[6](2, dim);
        dZdt_G += jac_sps[6](2, dim) * jac_sps[6](2, dim);
        dzdt_H += jac_cur_sps[7](2, dim) * jac_cur_sps[7](2, dim);
        dZdt_H += jac_sps[7](2, dim) * jac_sps[7](2, dim);
      }
      // E33: remedy of curvature thickness locking
      // Ett = 0.5* ( (1-r)(1-s)/4 * Ett(SP E) + ... + (1-r)(1+s)/4 * Ett(SP H) )
      lstrain(2) = 0.5 * (0.25 * (1.0 - r[gp]) * (1.0 - s[gp]) * (dzdt_E - dZdt_E) +
                             0.25 * (1.0 + r[gp]) * (1.0 - s[gp]) * (dzdt_F - dZdt_F) +
                             0.25 * (1.0 + r[gp]) * (1.0 + s[gp]) * (dzdt_G - dZdt_G) +
                             0.25 * (1.0 - r[gp]) * (1.0 + s[gp]) * (dzdt_H - dZdt_H));
      // E23 and E31
      if (ans_ == ans_lateral)
      {
        // E23: remedy of transverse shear locking
        // Est = (1+r)/2 * Est(SP B) + (1-r)/2 * Est(SP D)
        lstrain(4) = 0.5 * (1 + r[gp]) * (dydt_B - dYdt_B) + 0.5 * (1 - r[gp]) * (dydt_D - dYdt_D);
        // E13: remedy of transverse shear locking
        // Ert = (1-s)/2 * Ert(SP A) + (1+s)/2 * Ert(SP C)
        lstrain(5) = 0.5 * (1 - s[gp]) * (dxdt_A - dXdt_A) + 0.5 * (1 + s[gp]) * (dxdt_C - dXdt_C);
      }
      else if (ans_ == ans_onspot)
      {
#if 0
        const double hths_B = std::sqrt(metr_sps[1](2,2)/metr_sps[1](1,1));
        const double hths_D = std::sqrt(metr_sps[3](2,2)/metr_sps[3](1,1));
        const double hthr_A = std::sqrt(metr_sps[0](2,2)/metr_sps[0](0,0));
        const double hthr_C = std::sqrt(metr_sps[2](2,2)/metr_sps[2](0,0));
        const double atas_B = 1;///(1-metr_sps[1](0,2)/std::sqrt(metr_sps[1](0,0)*metr_sps[1](2,2)));
        const double atas_D = 1;///(1-metr_sps[3](0,2)/std::sqrt(metr_sps[3](0,0)*metr_sps[3](2,2)));
        const double atar_A = 1;///(1-metr_sps[0](1,2)/std::sqrt(metr_sps[0](1,1)*metr_sps[0](2,2)));
        const double atar_C = 1;///(1-metr_sps[2](1,2)/std::sqrt(metr_sps[2](1,1)*metr_sps[2](2,2)));
//        std::cout << "atas_B="  << atas_B << ", atas_D" << atas_D << ", atar_A=" << atar_A << "atar_C=" << atar_C << std::endl;
        // E23: remedy of transverse shear locking
        // Est = (1+r)/2 * Est(SP B) + (1-r)/2 * Est(SP D)
        lstrain(4)
          = 0.5*(1+r[gp]) * 0.5*(1-s[gp]) * ( (dydt_B - dYdt_B) + atas_B*hths_B*(dydt_F - dYdt_F) )
          + 0.5*(1+r[gp]) * 0.5*(1+s[gp]) * ( (dydt_B - dYdt_B) + atas_B*hths_B*(dydt_G - dYdt_G) )
          + 0.5*(1-r[gp]) * 0.5*(1-s[gp]) * ( (dydt_D - dYdt_D) + atas_D*hths_D*(dydt_E - dYdt_E) )
          + 0.5*(1-r[gp]) * 0.5*(1+s[gp]) * ( (dydt_D - dYdt_D) + atas_D*hths_D*(dydt_H - dYdt_H) );
        // E13: remedy of transverse shear locking
        // Ert = (1-s)/2 * Ert(SP A) + (1+s)/2 * Ert(SP C)
        lstrain(5)
          = 0.5*(1-s[gp]) * 0.5*(1-r[gp]) * ( (dxdt_A - dXdt_A) + atar_A*hthr_A*(dxdt_E - dXdt_E) )
          + 0.5*(1-s[gp]) * 0.5*(1+r[gp]) * ( (dxdt_A - dXdt_A) + atar_A*hthr_A*(dxdt_F - dXdt_F) )
          + 0.5*(1+s[gp]) * 0.5*(1-r[gp]) * ( (dxdt_C - dXdt_C) + atar_C*hthr_C*(dxdt_H - dXdt_H) )
          + 0.5*(1+s[gp]) * 0.5*(1+r[gp]) * ( (dxdt_C - dXdt_C) + atar_C*hthr_C*(dxdt_G - dXdt_G) );
#else
        // E23: remedy of transverse shear locking
        // Est = (1+r)/2 * Est(SP B) + (1-r)/2 * Est(SP D)
        lstrain(4) = 0.5 * (1 + r[gp]) * 0.5 * (1 - s[gp]) *
                         ((dydt_B - dYdt_B) + (*hths) * (dydt_F - dYdt_F)) +
                     0.5 * (1 + r[gp]) * 0.5 * (1 + s[gp]) *
                         ((dydt_B - dYdt_B) + (*hths) * (dydt_G - dYdt_G)) +
                     0.5 * (1 - r[gp]) * 0.5 * (1 - s[gp]) *
                         ((dydt_D - dYdt_D) + (*hths) * (dydt_E - dYdt_E)) +
                     0.5 * (1 - r[gp]) * 0.5 * (1 + s[gp]) *
                         ((dydt_D - dYdt_D) + (*hths) * (dydt_H - dYdt_H));
        // E13: remedy of transverse shear locking
        // Ert = (1-s)/2 * Ert(SP A) + (1+s)/2 * Ert(SP C)
        lstrain(5) = 0.5 * (1 - s[gp]) * 0.5 * (1 - r[gp]) *
                         ((dxdt_A - dXdt_A) + (*hthr) * (dxdt_E - dXdt_E)) +
                     0.5 * (1 - s[gp]) * 0.5 * (1 + r[gp]) *
                         ((dxdt_A - dXdt_A) + (*hthr) * (dxdt_F - dXdt_F)) +
                     0.5 * (1 + s[gp]) * 0.5 * (1 - r[gp]) *
                         ((dxdt_C - dXdt_C) + (*hthr) * (dxdt_H - dXdt_H)) +
                     0.5 * (1 + s[gp]) * 0.5 * (1 + r[gp]) *
                         ((dxdt_C - dXdt_C) + (*hthr) * (dxdt_G - dXdt_G));
#endif
      }
      else
      {
        dserror("Your should not turn up here.");
      }
    }
    // end of ANS modification of strains

    // push local/natural/parametric glstrains forward to global/material space
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain;
    glstrain.Multiply(TinvT, lstrain);

    // EAS technology: "enhance the strains"
    if (eastype_ == soh8_easnone)
    {
      ;  // continue
    }
    else if (eastype_ == soh8_eassosh8)
    {
      // map local M to global, also enhancement is refered to element origin
      EasMaterialiseShapeFcts<NUMEAS_SOSH8_>(M, detJ0, detJ, T0invT, M_GP->at(gp));
      // add enhanced strains = M . alpha to GL strains to "unlock" element
      EasAddStrain<NUMEAS_SOSH8_>(glstrain, M, alpha);
    }
    else if (eastype_ == soh8_easa)
    {
      // map local M to global, also enhancement is refered to element origin
      EasMaterialiseShapeFcts<NUMEAS_A_>(M, detJ0, detJ, T0invT, M_GP->at(gp));
      // add enhanced strains = M . alpha to GL strains to "unlock" element
      EasAddStrain<NUMEAS_A_>(glstrain, M, alpha);
    }
    else
    {
      dserror("EAS ain't easy.");
    }  // end of EAS modification of strains

    // recover deformation gradient incoperating assumed GL strain
    double detdefgrad;                         // determinant of assumed def.grad.
    LINALG::Matrix<NUMDIM_, NUMDIM_> defgrad;  // assumed def.grad.
    LINALG::Matrix<NUMDIM_, NUMDIM_>
        invdefgrad;                           // inverse of deformation gradient and its derivative
    LINALG::Matrix<NUMDIM_, NUMDIM_> rgtstr;  // assumed material stretch
    LINALG::Matrix<NUMDIM_, NUMDIM_> defgradD;  // pure disp-based def.grad.
    LINALG::Matrix<NUMDIM_, NUMDIM_> rgtstrD;   // pure disp.-based material stretch
    LINALG::Matrix<NUMDIM_, NUMDIM_>
        invrgtstrD;  // inverse of pure disp-based material stretch tensor U^{d;-1}
    AssDefGrad(detdefgrad, defgrad, invdefgrad, rgtstr, defgradD, rgtstrD, invrgtstrD, invJ_[gp],
        jac, jac_cur, glstrain);

    // Gauss weighted Jacobian mapping volume of spatial to parameter configuration
    const double vol_w = detdefgrad * detJ_w;

    // integrate volumetric deviation of element
    voldev += (1.0 - detdefgrad) * stabA[gp](0, 0) * detJ_w;

    // return gp strains if necessary
    if (iostrain != INPAR::STR::strain_none)
      Strain(elestrain, iostrain, gp, detdefgrad, defgrad, invdefgrad, glstrain);

    // call material law
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> stress(true);  // 2nd PK stress
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> cmat(true);
    params.set<int>("gp", gp);
    SolidMaterial()->Evaluate(&defgrad, &glstrain, params, &stress, &cmat, Id());
    if (iso_ == iso_enforced)
    {
      LINALG::Matrix<MAT::NUM_STRESS_3D, 1> pk2gen(
          stress);  // may contain non-isochoric material response
      LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> cgen(
          cmat);  // may contain non-isochoric material response
      MAT::VolumetrifyAndIsochorify(NULL, NULL, &stress, &cmat, glstrain, pk2gen, cgen);
    }
    // end of call material law

    // linearly interpolated pressure at Gauss point
    const double pressure = (shapefcts[gp]).Dot(pres);

    // return Gauss point stresses if necessary
    if (iostress != INPAR::STR::stress_none)
      Stress(elestress, iostress, gp, detdefgrad, defgrad, glstrain, stress, pressure);

    // effective shape function of scalar pressure field at current Gauss point
    LINALG::Matrix<NUMPRES_, 1> prshfct(true);
    if ((stab_ == stab_nonaffine) or (stab_ == stab_spatial))
      prshfct.MultiplyNT(1.0 / stabAA(0, 0), stabHA, stabA[gp]);
    else if ((stab_ == stab_affine) or (stab_ == stab_spatialaffine) or (stab_ == stab_puredisp))
      prshfct.Update(shapefcts[gp]);
    else
      dserror("Cannot handle requested stabilisation type %d", stab_);

    // update of internal force vector
    if (force != NULL)
    {
      // integrate internal force vector
      // fint := fint
      //      + (B^T . sigma) * detJ * w(gp)
      //      + (-G) . ep   // will be done _after_ Gauss point loop
      force->MultiplyTN(detJ_w, bop, stress, 1.0);
    }

    // incompressiblity equation
    if (incomp != NULL)
    {
      // pint := pint
      //       - He . (Fdet - 1.0) * detJ * wp(gp)
      //       + (-Ce) . ep   // will be done _after_ Gauss point loop
      incomp->Update(-(detdefgrad - 1.0) * detJ_w, prshfct, 1.0);
    }

    // update of stiffness matrix
    if (stiffmatrix != NULL)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDISP_> cb;
      cb.Multiply(cmat, bop);  // temporary C . B
      stiffmatrix->MultiplyTN(detJ_w, bop, cb, 1.0);

      // integrate `geometric' stiffness matrix and add to keu
      // here also the ANS interpolation comes into play
      Teuchos::RCP<LINALG::Matrix<MAT::NUM_STRESS_3D, NUMNOD_* NUMNOD_>> bopbydisp = Teuchos::null;
      if (lin_ > lin_sixth)
        bopbydisp = Teuchos::rcp(new LINALG::Matrix<MAT::NUM_STRESS_3D, NUMNOD_ * NUMNOD_>());
      for (int jnod = 0; jnod < NUMNOD_; ++jnod)
      {
        for (int inod = 0; inod < NUMNOD_; ++inod)
        {
          LINALG::Matrix<MAT::NUM_STRESS_3D, 1> G_ij;
          G_ij(0) = derivs[gp](0, inod) * derivs[gp](0, jnod);  // rr-dir
          G_ij(1) = derivs[gp](1, inod) * derivs[gp](1, jnod);  // ss-dir
          G_ij(3) = derivs[gp](0, inod) * derivs[gp](1, jnod) +
                    derivs[gp](1, inod) * derivs[gp](0, jnod);  // rs-dir
          if (ans_ == ans_none)
          {
            G_ij(2) = derivs[gp](2, inod) * derivs[gp](2, jnod);  // tt-dir
            G_ij(4) = derivs[gp](1, inod) * derivs[gp](2, jnod) +
                      derivs[gp](2, inod) * derivs[gp](1, jnod);  // st-dir
            G_ij(5) = derivs[gp](2, inod) * derivs[gp](0, jnod) +
                      derivs[gp](0, inod) * derivs[gp](2, jnod);  // tr-dir
          }
          else
          {
            // ANS modification in tt-dir
            G_ij(2) = 0.25 * (1.0 - r[gp]) * (1.0 - s[gp]) * (*deriv_sp)[4](2, inod) *
                          (*deriv_sp)[4](2, jnod) +
                      0.25 * (1.0 + r[gp]) * (1.0 - s[gp]) * (*deriv_sp)[5](2, inod) *
                          (*deriv_sp)[5](2, jnod) +
                      0.25 * (1.0 + r[gp]) * (1.0 + s[gp]) * (*deriv_sp)[6](2, inod) *
                          (*deriv_sp)[6](2, jnod) +
                      0.25 * (1.0 - r[gp]) * (1.0 + s[gp]) * (*deriv_sp)[7](2, inod) *
                          (*deriv_sp)[7](2, jnod);
            // st-dir and rt-dir
            if (ans_ == ans_lateral)
            {
              // ANS modification in st-dir
              G_ij(4) =
                  0.5 * ((1.0 + r[gp]) * ((*deriv_sp)[1](1, inod) * (*deriv_sp)[1](2, jnod) +
                                             (*deriv_sp)[1](2, inod) * (*deriv_sp)[1](1, jnod)) +
                            (1.0 - r[gp]) * ((*deriv_sp)[3](1, inod) * (*deriv_sp)[3](2, jnod) +
                                                (*deriv_sp)[3](2, inod) * (*deriv_sp)[3](1, jnod)));
              // ANS modification in rt-dir
              G_ij(5) =
                  0.5 * ((1.0 - s[gp]) * ((*deriv_sp)[0](0, inod) * (*deriv_sp)[0](2, jnod) +
                                             (*deriv_sp)[0](2, inod) * (*deriv_sp)[0](0, jnod)) +
                            (1.0 + s[gp]) * ((*deriv_sp)[2](0, inod) * (*deriv_sp)[2](2, jnod) +
                                                (*deriv_sp)[2](2, inod) * (*deriv_sp)[2](0, jnod)));
            }
            else if (ans_ == ans_onspot)
            {
              // ANS modification in st-dir
              G_ij(4) = 0.5 * (1.0 + r[gp]) * 0.5 * (1 - s[gp]) *
                            ((*deriv_sp)[1](1, inod) * (*deriv_sp)[1](2, jnod) +
                                (*deriv_sp)[1](2, inod) * (*deriv_sp)[1](1, jnod))  // B
                        + 0.5 * (1.0 + r[gp]) * 0.5 * (1 - s[gp]) * (*hths) *
                              ((*deriv_sp)[5](1, inod) * (*deriv_sp)[5](2, jnod) +
                                  (*deriv_sp)[5](2, inod) * (*deriv_sp)[5](1, jnod))  // F
                        + 0.5 * (1.0 + r[gp]) * 0.5 * (1 + s[gp]) *
                              ((*deriv_sp)[1](1, inod) * (*deriv_sp)[1](2, jnod) +
                                  (*deriv_sp)[1](2, inod) * (*deriv_sp)[1](1, jnod))  // B
                        + 0.5 * (1.0 + r[gp]) * 0.5 * (1 + s[gp]) * (*hths) *
                              ((*deriv_sp)[6](1, inod) * (*deriv_sp)[6](2, jnod) +
                                  (*deriv_sp)[6](2, inod) * (*deriv_sp)[6](1, jnod))  // G
                        + 0.5 * (1.0 - r[gp]) * 0.5 * (1 - s[gp]) *
                              ((*deriv_sp)[3](1, inod) * (*deriv_sp)[3](2, jnod) +
                                  (*deriv_sp)[3](2, inod) * (*deriv_sp)[3](1, jnod))  // D
                        + 0.5 * (1.0 - r[gp]) * 0.5 * (1 - s[gp]) * (*hths) *
                              ((*deriv_sp)[4](1, inod) * (*deriv_sp)[4](2, jnod) +
                                  (*deriv_sp)[4](2, inod) * (*deriv_sp)[4](1, jnod))  // E
                        + 0.5 * (1.0 - r[gp]) * 0.5 * (1 + s[gp]) *
                              ((*deriv_sp)[3](1, inod) * (*deriv_sp)[3](2, jnod) +
                                  (*deriv_sp)[3](2, inod) * (*deriv_sp)[3](1, jnod))  // D
                        + 0.5 * (1.0 - r[gp]) * 0.5 * (1 + s[gp]) * (*hths) *
                              ((*deriv_sp)[7](1, inod) * (*deriv_sp)[7](2, jnod) +
                                  (*deriv_sp)[7](2, inod) * (*deriv_sp)[7](1, jnod));  // H
              // ANS modification in rt-dir
              G_ij(5) = 0.5 * (1.0 - s[gp]) * 0.5 * (1 - r[gp]) *
                            ((*deriv_sp)[0](0, inod) * (*deriv_sp)[0](2, jnod) +
                                (*deriv_sp)[0](2, inod) * (*deriv_sp)[0](0, jnod))  // A
                        + 0.5 * (1.0 - s[gp]) * 0.5 * (1 - r[gp]) * (*hthr) *
                              ((*deriv_sp)[4](0, inod) * (*deriv_sp)[4](2, jnod) +
                                  (*deriv_sp)[4](2, inod) * (*deriv_sp)[4](0, jnod))  // E
                        + 0.5 * (1.0 - s[gp]) * 0.5 * (1 + r[gp]) *
                              ((*deriv_sp)[0](0, inod) * (*deriv_sp)[0](2, jnod) +
                                  (*deriv_sp)[0](2, inod) * (*deriv_sp)[0](0, jnod))  // A
                        + 0.5 * (1.0 - s[gp]) * 0.5 * (1 + r[gp]) * (*hthr) *
                              ((*deriv_sp)[5](0, inod) * (*deriv_sp)[5](2, jnod) +
                                  (*deriv_sp)[5](2, inod) * (*deriv_sp)[5](0, jnod))  // F
                        + 0.5 * (1.0 + s[gp]) * 0.5 * (1 - r[gp]) *
                              ((*deriv_sp)[2](0, inod) * (*deriv_sp)[2](2, jnod) +
                                  (*deriv_sp)[2](2, inod) * (*deriv_sp)[2](0, jnod))  // C
                        + 0.5 * (1.0 + s[gp]) * 0.5 * (1 - r[gp]) * (*hthr) *
                              ((*deriv_sp)[7](0, inod) * (*deriv_sp)[7](2, jnod) +
                                  (*deriv_sp)[7](2, inod) * (*deriv_sp)[7](0, jnod))  // H
                        + 0.5 * (1.0 + s[gp]) * 0.5 * (1 + r[gp]) *
                              ((*deriv_sp)[2](0, inod) * (*deriv_sp)[2](2, jnod) +
                                  (*deriv_sp)[2](2, inod) * (*deriv_sp)[2](0, jnod))  // C
                        + 0.5 * (1.0 + s[gp]) * 0.5 * (1 + r[gp]) * (*hthr) *
                              ((*deriv_sp)[6](0, inod) * (*deriv_sp)[6](2, jnod) +
                                  (*deriv_sp)[6](2, inod) * (*deriv_sp)[6](0, jnod));  // G
            }
            else
            {
              dserror("You should not turn up here.");
            }
          }

          // transformation of local(parameter) space 'back' to global(material) space
          LINALG::Matrix<MAT::NUM_STRESS_3D, 1> G_ij_glob;
          G_ij_glob.Multiply(TinvT, G_ij);

          // store B_{aBd,k}
          if (lin_ > lin_sixth)
            for (int istr = 0; istr < MAT::NUM_STRESS_3D; ++istr)
              (*bopbydisp)(istr, NUMNOD_ * inod + jnod) = G_ij_glob(istr);

          // Scalar Gij results from product of G_ij with stress, scaled with detJ*weights
          double Gij = detJ_w * stress.Dot(G_ij_glob);

          // add "geometric part" Gij times detJ*weights to stiffness matrix
          (*stiffmatrix)(NUMDIM_ * inod + 0, NUMDIM_ * jnod + 0) += Gij;
          (*stiffmatrix)(NUMDIM_ * inod + 1, NUMDIM_ * jnod + 1) += Gij;
          (*stiffmatrix)(NUMDIM_ * inod + 2, NUMDIM_ * jnod + 2) += Gij;
        }
      }  // end of integrate `geometric' stiffness

      // add (incomplete) derivative of pressure-proportional force w.r.t. displacements
      // Kp = (-Gm*ep')_,d
      //    = (-dFv'*fvT*pN * Fdet*detJ*wp(i))_,d * ep'
      // Kp = Kp + dFv'*fvT*(pN*ep')*fvT'*dFv * Fdet*detJ*wp(i)  // due to Fdet_,d = fvT'*dFv * Fdet
      //         + (pN*ep')*dFv'*WmT*dFv * Fdet*detJ*wp(i)       // due to fvT_,d = WmT*dFv
      //         + ddFv * Fdet*detJ*wp(i)                        // due to ddFv = dFv'_,d*fvT =
      //         (fv'*dFv_,d)'
      // Ke = Keu + Kg - Kp;
      {
        // effective pressure at Gauss point
        const double effpressure = prshfct.Dot(pres);  // pN*ep'
        // Voigt 9-vector of transposed & inverted deformation gradient fvT := F^{-T}
        LINALG::Matrix<NUMDFGR_, 1> tinvdefgrad;
        Matrix2TensorToVector9Voigt(tinvdefgrad, invdefgrad, true);

        // derivative of WmT := F^{-T}_{,F} in Voigt vector notation
        LINALG::Matrix<NUMDFGR_, NUMDFGR_> WmT;
        InvVector9VoigtDiffByItself(WmT, invdefgrad, true);
        // WmT := WmT + fvT*fvT'
        WmT.MultiplyNT(1.0, tinvdefgrad, tinvdefgrad, 1.0);

        // material derivatives of shape functions
        LINALG::Matrix<NUMDIM_, NUMNOD_> derivsmat;
        derivsmat.MultiplyNN(invJ_[gp], derivs[gp]);

        // linear B-op
        // derivative of displ-based def.grad with respect to nodal displacements
        // F^d_{aC,d}
        int iboplin[NUMNOD_][NUMDFGR_];  // index entries which are non-equal zero
        LINALG::Matrix<NUMDFGR_, NUMDISP_> boplin(true);
        for (int m = 0; m < NUMNOD_; ++m)
        {
          for (int ij = 0; ij < NUMDFGR_; ++ij)
          {
            const int i = VOIGT9ROW_[ij];
            const int j = VOIGT9COL_[ij];
            const int k = m * NODDISP_ + i;
            iboplin[m][ij] = k;
            boplin(ij, k) = derivsmat(j, m);
          }
        }

        // derivative of def.grad. w.r.t. displacements Fv_{,d}
        LINALG::Matrix<NUMDFGR_, NUMDISP_> defgradbydisp;
        if (ans_ == ans_none)
        {
          defgradbydisp.Update(boplin);
        }
        else
        {
          // derivative of pure-disp. inverse material stretch tensor with respect to nodal
          // displacements U^{d;-1}_{,d} = U^{d;-1}_{,U} . (C^d_{,U^d})^{-1} . C^d_{,d} on exit of
          // this block the following variables are going to hold ...
          LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDISP_> invrgtstrDbydisp;  // ...U^{d;-1}_{,d}
          LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>
              invrgtstrDbyrgtstrD;  // ...U^{d;-1}_{,U}
          LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>
              rcgDbyrgtstrD;                                           // ...(C^d_{,U^d})^{-1}
          LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDISP_> rgtstrDbydisp;  // ...U^d_{,d}
          {
            // U^{d;-1}_{,U}
            InvVector6VoigtDiffByItself(invrgtstrDbyrgtstrD, invrgtstrD);

            // C^d_{,U^d} = (U^d . U^d)_{,U^d}
            SqVector6VoigtDiffByItself(rcgDbyrgtstrD, rgtstrD);

            // displ-based deformation gradient as Voigt matrix
            LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDFGR_> defgradDm;
            Matrix2TensorToMatrix6x9Voigt(defgradDm, defgradD, true);

            // C^d_{,d} = 2 * Fm^d * Boplin, 6x24
            LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDISP_> rcgDbydisp;
            rcgDbydisp.MultiplyNN(2.0, defgradDm, boplin);

            // U^d_{,d} = (C^d_{,U^d})^{-1} . C^d_{,d}
            {
              LINALG::FixedSizeSerialDenseSolver<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D, NUMDISP_>
                  rcgDbyrgtstrDsolver;
              rcgDbyrgtstrDsolver.SetMatrix(rcgDbyrgtstrD);               // LHS
              rcgDbyrgtstrDsolver.SetVectors(rgtstrDbydisp, rcgDbydisp);  // SOL, RHS
              const int err = rcgDbyrgtstrDsolver.Solve();
              if (err != 0) dserror("Failed to solve, error=%d", err);
              if (lin_ >= lin_one)
              {
                const int err = rcgDbyrgtstrDsolver.Invert();
                if (err != 0) dserror("Failed to invert, error=%d", err);
              }
            }

            // U^{d;-1}_{,d} = U^{d;-1}_{,U} . U^d_{,d}
            invrgtstrDbydisp.MultiplyNN(invrgtstrDbyrgtstrD, rgtstrDbydisp);
          }

          // derivative of ass. mat. stretch tensor with respect to nodal displacements
          // U^{ass}_{,d} = (C^{ass}_{,U^{ass}})^{-1} . C^{ass}_{,d}
          {
            // derivative of ass. right Cauchy-Green with respect to ass. material stretch tensor
            // C^{ass}_{,U^{ass}}
            LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> rcgbyrgtstr;
            SqVector6VoigtDiffByItself(rcgbyrgtstr, rgtstr);

            // C^{ass}_{,d} = 2 * bop
            // C^{ass}_{AB,k} = 2 B_{ABk}

            // derivative of ass. mat. stretch tensor with respect to nodal displacements
            // U^{ass}_{,d} = (C^{ass}_{,U^{ass}})^{-1} . C^{ass}_{,d}
            LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDISP_> rgtstrbydisp;  // ... U^{ass}_{,d}
            LINALG::FixedSizeSerialDenseSolver<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D, NUMDISP_>
                rcgbyrgtstrsolver;
            {
              rcgbyrgtstrsolver.SetMatrix(rcgbyrgtstr);         // LHS
              rcgbyrgtstrsolver.SetVectors(rgtstrbydisp, bop);  // SOL, RHS
              const int err = rcgbyrgtstrsolver.Solve();
              if (err != 0) dserror("Failed to solve, error=%d", err);
              if (lin_ > lin_half)
              {
                const int err = rcgbyrgtstrsolver.Invert();
                if (err != 0) dserror("Failed to invert, error=%d", err);
              }
            }
            rgtstrbydisp.Scale(2.0);

            // derivative of def.grad. with respect to k nodal displacements d^k
            {
              // pseudo identity 2-tensor
              // I^{assd}_{CB} = U^{d;-1}_{CD} . U^{ass}_{DB}
              LINALG::Matrix<NUMDIM_, NUMDIM_> invrgtstrDtimesrgtstr;
              invrgtstrDtimesrgtstr.MultiplyNN(invrgtstrD, rgtstr);

              // derivative of pseudo identity with respect to displacements
              // I^{assd}_{CB,k} = U^{d;-1}_{CD,k} . U^{ass}_{DB}
              //                 + U^{d;-1}_{CD} . U^{ass}_{DB,k}
              // WARNING: I^{assd}_{CB} and I^{assd}_{CB,k} might be non-symmetric in CB
              LINALG::Matrix<NUMDFGR_, NUMDISP_> pseudoidenity;
              for (int k = 0; k < NUMDISP_; ++k)
              {
                for (int CB = 0; CB < NUMDFGR_; ++CB)
                {
                  const int C = VOIGT9ROW_[CB];
                  const int B = VOIGT9COL_[CB];
                  double pseudoidenity_CBk = 0.0;
                  for (int D = 0; D < NUMDIM_; ++D)
                  {
                    const int CD = VOIGT3X3SYM_[NUMDIM_ * C + D];
                    const int DB = VOIGT3X3SYM_[NUMDIM_ * D + B];
                    const double CDfact = (C == D) ? 1.0 : 0.5;
                    const double DBfact = (D == B) ? 1.0 : 0.5;
                    pseudoidenity_CBk += CDfact * invrgtstrDbydisp(CD, k) * rgtstr(D, B) +
                                         invrgtstrD(C, D) * DBfact * rgtstrbydisp(DB, k);
                  }
                  pseudoidenity(CB, k) = pseudoidenity_CBk;
                }
              }

              // derivative of def.grad. with respect to k nodal displacements d^k
              // F_{aB,k} = F^d_{aC,k} . U^{d;-1}_{CD} . U^{ass}_{DB}
              //          + F^d_{aC} . U^{d;-1}_{CD,k} . U^{ass}_{DB}
              //          + F^d_{aC} . U^{d;-1}_{CD} . U^{ass}_{DB,k}
              //          = F^d_{aC,k} . I^{assd}_{CB}
              //          + F^d_{aC} . I^{assd}_{CB,k}
              for (int k = 0; k < NUMDISP_; ++k)
              {
                const int m = k / NODDISP_;
                for (int aB = 0; aB < NUMDFGR_; ++aB)
                {
                  const int a = VOIGT9ROW_[aB];
                  const int B = VOIGT9COL_[aB];
                  double defgradbydisp_aBk = 0.0;
                  for (int C = 0; C < NUMDIM_; ++C)
                  {
                    const int CB = VOIGT3X3_[NUMDIM_ * C + B];
                    defgradbydisp_aBk += defgradD(a, C) * pseudoidenity(CB, k);
                  }
                  for (int C = 0; C < NUMDIM_; ++C)
                  {
                    const int aC = VOIGT3X3_[NUMDIM_ * a + C];
                    if (k == iboplin[m][aC])
                      defgradbydisp_aBk += boplin(aC, k) * invrgtstrDtimesrgtstr(C, B);
                  }
                  defgradbydisp(aB, k) = defgradbydisp_aBk;
                }
              }
            }  // end block

            // EAS technology: integrate matrices
            if (eastype_ == soh8_easnone)
            {
              ;  // continue
            }
            else if (eastype_ == soh8_eassosh8)
            {
              EasConstraintAndTangent<NUMEAS_SOSH8_>(feas, Kaa, Kad, Kap, defgradD, invrgtstrD,
                  rcgbyrgtstr, detdefgrad, tinvdefgrad, WmT, cmat, stress, effpressure, detJ_w, cb,
                  defgradbydisp, prshfct, M);
            }
            else if (eastype_ == soh8_easa)
            {
              EasConstraintAndTangent<NUMEAS_A_>(feas, Kaa, Kad, Kap, defgradD, invrgtstrD,
                  rcgbyrgtstr, detdefgrad, tinvdefgrad, WmT, cmat, stress, effpressure, detJ_w, cb,
                  defgradbydisp, prshfct, M);
            }
            else
            {
              dserror("EAS ain't easy.");
            }  // if (eastype_ == ... )

            // ext(p)ensive computation to achieve full tangent
            if (lin_ > lin_half)
            {
              // REMARK:
              // on #rcgDbyrgtstrD is stored the inverse of C^{d}_{,U^{d}}
              // on #rcgbyrgtstr is stored the inverse of C^{ass}_{,U^{ass}}

              // second derivative of pure-disp right Cauchy-Green tensor w.r.t. displacements
              // TEMP^ass = C^{ass}_{EF,dk} - C^{ass}_{,UU} . U^{ass}_{,d} . U^{ass}_{,k}
              // TEMP^d = C^{d}_{EF,dk} - C^{d}_{,UU} . U^{d}_{,d} . U^{d}_{,k}
              LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDISPSQSYM_> temp;
              LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDISPSQSYM_> tempD;
              {
                // second derivative of assumed right Cauchy-Green tensor
                // w.r.t. to right stretch tensor
                // C^{ass}_{,U^{ass} U^{ass}} = const
                int ircgbybyrgtstr[MAT::NUM_STRESS_3D * 6];  // for sparse access
                LINALG::Matrix<MAT::NUM_STRESS_3D, 6> rcgbybyrgtstr;
                SqVector6VoigtTwiceDiffByItself(ircgbybyrgtstr, rcgbybyrgtstr);

                // second derivative of disp-based right Cauchy-Green tensor
                // w.r.t. to right stretch tensor
                // C^{d}_{,U^{d} U^{d}} = const
                // MARK: an extra variable is not needed as same as for assumed right CG tensor
                // (above)

                // second derivative of pure-disp right Cauchy-Green tensor w.r.t. displacements
                // TEMP^ass = C^{ass}_{EF,dk} - C^{ass}_{,UU} . U^{ass}_{,d} . U^{ass}_{,d}
                // TEMP^d = C^{d}_{EF,dk} - C^{d}_{,UU} . U^{d}_{,d} . U^{d}_{,d}
                for (int d = 0; d < NUMDISP_; ++d)
                {
                  for (int k = d; k < NUMDISP_; ++k)
                  {  // symmetric in d,k
                    const int dk = NUMDISP_ * d - d * (d + 1) / 2 + k;
                    int ndnk = -1;
                    if (d % NODDISP_ == k % NODDISP_)
                    {
                      const int nd = d / NODDISP_;
                      const int nk = k / NODDISP_;
                      ndnk = nd * NUMNOD_ + nk;
                    }
                    for (int EF = 0; EF < MAT::NUM_STRESS_3D; ++EF)
                    {
                      const int E = VOIGT6ROW_[EF];
                      const int F = VOIGT6COL_[EF];
                      // C^{ass}_{,UU} . U^{ass}_{,d} . U^{ass}_{,d}
                      // and
                      // C^{d}_{,UU} . U^{d}_{,d} . U^{d}_{,d}
                      double temp_EFdk = 0.0;
                      double tempD_EFdk = 0.0;
                      // this looks terrible, but speed is everything
                      for (int GHIJ = 0; GHIJ < 6; ++GHIJ)
                      {
                        if (ircgbybyrgtstr[MAT::NUM_STRESS_3D * EF + GHIJ] != -1)
                        {
                          const int GH =
                              ircgbybyrgtstr[MAT::NUM_STRESS_3D * EF + GHIJ] / MAT::NUM_STRESS_3D;
                          const int IJ =
                              ircgbybyrgtstr[MAT::NUM_STRESS_3D * EF + GHIJ] % MAT::NUM_STRESS_3D;
                          // C^{ass}_{,UU} . U^{ass}_{,d} . U^{ass}_{,d}
                          temp_EFdk +=
                              rcgbybyrgtstr(EF, GHIJ) * rgtstrbydisp(IJ, d) * rgtstrbydisp(GH, k);
                          // C^{d}_{,UU} . U^{d}_{,d} . U^{d}_{,d}
                          // C^{d}_{,U^d U^d} = C^{ass}_{,U^ass U^ass} = const
                          tempD_EFdk +=
                              rcgbybyrgtstr(EF, GHIJ) * rgtstrDbydisp(IJ, d) * rgtstrDbydisp(GH, k);
                        }
                      }
                      // U^{ass}_{DB,dk}
                      double rcgbybydisp_EFdk = 0.0;
                      if (ndnk != -1) rcgbybydisp_EFdk = 2.0 * (*bopbydisp)(EF, ndnk);
                      double rcgDbybydisp_EFdk = 0.0;
                      for (int a = 0; a < NUMDIM_; ++a)
                      {
                        const int aE = VOIGT3X3_[NUMDIM_ * a + E];
                        const int aF = VOIGT3X3_[NUMDIM_ * a + F];
                        if (E == F)  // make strain-like 6-Voigt vector
                          rcgDbybydisp_EFdk += 2.0 * boplin(aE, d) * boplin(aF, k);
                        else  // thus setting  V_EF + V_FE if E!=F
                          rcgDbybydisp_EFdk += 2.0 * boplin(aE, d) * boplin(aF, k) +
                                               2.0 * boplin(aF, d) * boplin(aE, k);
                      }
                      temp(EF, dk) = rcgbybydisp_EFdk - temp_EFdk;
                      tempD(EF, dk) = rcgDbybydisp_EFdk - tempD_EFdk;
                    }
                  }
                }
              }

              // second derivative of assumed right stretch tensor w.r.t. displacements
              // U^{ass}_{DB,dk} = (C^{ass}_{,U^{ass}})_{DBEF}^{-1}
              //                 . ( C^{ass}_{EF,dk} - C^{ass}_{EF,GHIJ}  U^{ass}_{IJ,d}
              //                 U^{ass}_{GH,k} )
              // and
              // second derivative of pure-disp right stretch tensor w.r.t. displacements
              // U^{d}_{DB,dk} = (C^{d}_{,U^{d}})_{DBEF}^{-1}
              //               . ( C^{d}_{EF,dk} - C^{d}_{EF,GHIJ}  U^{d}_{IJ,d}  U^{d}_{GH,k} )
              LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDISPSQSYM_> rgtstrbybydisp;
              LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDISPSQSYM_> rgtstrDbybydisp;
              for (int d = 0; d < NUMDISP_; ++d)
              {
                for (int k = d; k < NUMDISP_; ++k)
                {  // symmetric in d,k
                  const int dk = NUMDISP_ * d - d * (d + 1) / 2 + k;
                  for (int DB = 0; DB < MAT::NUM_STRESS_3D; ++DB)
                  {
                    double rgtstrbybydisp_DBdk = 0.0;
                    double rgtstrDbybydisp_DBdk = 0.0;
                    for (int EF = 0; EF < MAT::NUM_STRESS_3D; ++EF)
                    {
                      // TEMP^ass = C^{ass}_{EF,dk} - C^{ass}_{,UU} . U^{ass}_{,d} . U^{ass}_{,k}
                      const double temp_EFdk = temp(EF, dk);
                      // (C^{ass}_{,U^{ass}})_{DBEF}^{-1}
                      const double rcgbyrgtstr_DBEF = rcgbyrgtstr(DB, EF);
                      // U^{ass}_{DB,dk}
                      rgtstrbybydisp_DBdk += rcgbyrgtstr_DBEF * temp_EFdk;
                      // TEMP^d = C^{d}_{EF,dk} - C^{d}_{,UU} . U^{d}_{,d} . U^{d}_{,k}
                      const double tempD_EFdk = tempD(EF, dk);
                      // (C^{d}_{,U^{d}})_{DBEF}^{-1}
                      const double rcgDbyrgtstrD_DBEF = rcgDbyrgtstrD(DB, EF);
                      // U^{d}_{DB,dk}
                      rgtstrDbybydisp_DBdk += rcgDbyrgtstrD_DBEF * tempD_EFdk;
                    }
                    rgtstrbybydisp(DB, dk) = rgtstrbybydisp_DBdk;
                    rgtstrDbybydisp(DB, dk) = rgtstrDbybydisp_DBdk;
                  }
                }
              }

              // add (very pricy) contribution of second derivative with respect
              // to element displacements of (assumed) deformation gradient
              if (lin_ >= lin_one)
              {
                // second derivative of disp-based inverse right stretch tensor
                // w.r.t. disp-based right stretch tensor
                // U^{d-1}_{,U^d U^d}
                LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D * MAT::NUM_STRESS_3D>
                    invrgtstrDbybyrgtstrD;
                InvVector6VoigtTwiceDiffByItself(invrgtstrDbybyrgtstrD, invrgtstrD);

                // second derivative of pure-disp inverse right stretch tensor w.r.t. displacements
                // U^{d-1}_{CD,dk} = U^{d-1}_{CD,EFGH} U^{d}_{GH,k} U^{d}_{EF,d}
                //                 + U^{d-1}_{CD,EF} U^{d}_{EF,dk}
                LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDISPSQSYM_>
                    invrgtstrDbybydisp;  // ... U^{d-1}_{DB,dk}
                for (int d = 0; d < NUMDISP_; ++d)
                {
                  for (int k = d; k < NUMDISP_; ++k)
                  {  // symmetric in d,k
                    const int dk = NUMDISP_ * d - d * (d + 1) / 2 + k;
                    for (int CD = 0; CD < MAT::NUM_STRESS_3D; ++CD)
                    {
                      double invrgtstrDbybydisp_CDdk = 0.0;
                      for (int EF = 0; EF < MAT::NUM_STRESS_3D; ++EF)
                      {
                        const double rgtstrDbybydisp_EFdk = rgtstrDbybydisp(EF, dk);
                        invrgtstrDbybydisp_CDdk +=
                            invrgtstrDbyrgtstrD(CD, EF) * rgtstrDbybydisp_EFdk;
                        for (int GH = 0; GH < MAT::NUM_STRESS_3D; ++GH)
                        {
                          const int EFGH = MAT::NUM_STRESS_3D * EF + GH;
                          invrgtstrDbybydisp_CDdk +=
                              invrgtstrDbybyrgtstrD(CD, EFGH)  // col are strain-like 6-Voigt
                              * rgtstrDbydisp(GH, k)           // row are strain-like 6-Voigt too
                              * rgtstrDbydisp(EF, d);          // row are strain-like 6-Voigt too
                        }
                      }
                      invrgtstrDbybydisp(CD, dk) = invrgtstrDbybydisp_CDdk;
                    }
                  }
                }


                // inverse assumed deformation gradient times #boplin
                // Z^{assd}_{BCd} = F^{-1}_{Ba}  F^d_{aC,d}
                // WARNING: Z^{assd}_{BCd} might be non-symmetric in BC
                LINALG::Matrix<NUMDFGR_, NUMDISP_> invdefgradtimesboplin;
                for (int d = 0; d < NUMDISP_; ++d)
                {
                  for (int BC = 0; BC < NUMDFGR_; ++BC)
                  {
                    const int B = VOIGT9ROW_[BC];
                    const int C = VOIGT9COL_[BC];
                    double invdefgradtimesboplin_BCd = 0.0;
                    for (int a = 0; a < NUMDIM_; ++a)
                    {
                      const int aC = VOIGT3X3_[NUMDIM_ * a + C];
                      invdefgradtimesboplin_BCd += invdefgrad(B, a) * boplin(aC, d);
                    }
                    invdefgradtimesboplin(BC, d) = invdefgradtimesboplin_BCd;
                  }
                }

                // pseudo identity matrix
                // I^{assd}_{BC} = F^{-1}_{Ba} . F^{d}_{aC}
                //               = (R_{aE} . U_{EB})^{-1}_{Ba} . (R_{aD} . U_{DC})^{d}_{aC}
                //               = U^{-1}_{BE} . (R^{-1})_{Ea} . R^{d}_{aD} . U^{d}_{DC}
                // because we have only one rotation matrix R^{d}=R
                // I^{assd}_{BC} = U^{-1}_{BE} . (R^{-1})_{Ea} . R_{aD} . U^{d}_{DC}
                //               = U^{-1}_{BE} . I_{ED} . U^{d}_{DC}
                //               = U^{-1}_{BD} . U^{d}_{DC}
                LINALG::Matrix<NUMDIM_, NUMDIM_> invdefgradtimesdefgradD;
                invdefgradtimesdefgradD.MultiplyNN(invdefgrad, defgradD);

                // pseudo inverse disp-based right stretch tensor
                // U^{assd-1}_{BD} = ( F^{-1}_{Ba} . F^{d}_{aC} ) . U^{d;-1}_{CD}
                //                 = I^{assd}_{BC} . U^{d;-1}_{CD}
                // By the way, we have the identity
                // U^{assd-1}_{BD} = U^{-1}_{BE} . U^{d}_{EC} . U^{d;-1}_{CD}
                //                 = U^{-1}_{BE} . I_{ED}
                //                 = U^{-1}_{BD}
                LINALG::Matrix<NUMDIM_, NUMDIM_> invdefgradtimesdefgradDtimesinvrgtstrD;
                invdefgradtimesdefgradDtimesinvrgtstrD.MultiplyNN(
                    invdefgradtimesdefgradD, invrgtstrD);

                // contribute stuff containing second derivatives in displacements
                // F^{-T}_{aB} F_{aB,dk}
                // = F^{-1}_{Ba} F_{aB,dk}
                // = F^{-1}_{Ba}  ( F^d_{aC} U^{d;-1}_{CD} U^{ass}_{DB} )_{,dk}
                //
                // = F^{-1}_{Ba}  F^d_{aC,dk} U^{d;-1}_{CD} U^{ass}_{DB}         |  = 0
                // + F^{-1}_{Ba}  F^d_{aC} U^{d;-1}_{CD,dk} U^{ass}_{DB}         |  # 0, very pricy
                // + F^{-1}_{Ba}  F^d_{aC} U^{d;-1}_{CD} U^{ass}_{DB,dk}         |  # 0, pricy
                // + F^{-1}_{Ba}  F^d_{aC,d} U^{d;-1}_{CD,k} U^{ass}_{DB}        |  # 0, okay
                // + F^{-1}_{Ba}  F^d_{aC,d} U^{d;-1}_{CD} U^{ass}_{DB,k}        |  # 0, okay
                // + F^{-1}_{Ba}  F^d_{aC,k} U^{d;-1}_{CD,d} U^{ass}_{DB}        |  # 0, okay
                // + F^{-1}_{Ba}  F^d_{aC} U^{d;-1}_{CD,d} U^{ass}_{DB,k}        |  # 0, okay
                // + F^{-1}_{Ba}  F^d_{aC,k} U^{d;-1}_{CD} U^{ass}_{DB,d}        |  # 0, okay
                // + F^{-1}_{Ba}  F^d_{aC} U^{d;-1}_{CD,k} U^{ass}_{DB,d}        |  # 0, okay
                //
                // = F^{-1}_{Ba}  F^d_{aC,dk} U^{d;-1}_{CD} U^{ass}_{DB}         // zero
                // + I^{assd}_{BC} U^{d;-1}_{CD,dk} U^{ass}_{DB}                 //    (7)
                // + U^{assd-1}_{BD} U^{ass}_{DB,dk}                             //       (8)
                // + Z^{assd}_{BCd} U^{d;-1}_{CD,k} U^{ass}_{DB}                 // (1)
                // + Z^{assd}_{BCd} U^{d;-1}_{CD} U^{ass}_{DB,k}                 // (2)
                // + Z^{assd}_{BCk} U^{d;-1}_{CD,d} U^{ass}_{DB}                 // (3)
                // + I^{assd}_{BC} U^{d;-1}_{CD,d} U^{ass}_{DB,k}                //    (5)
                // + Z^{assd}_{BCk} U^{d;-1}_{CD} U^{ass}_{DB,d}                 // (4)
                // + I^{assd}_{BC} U^{d;-1}_{CD,k} U^{ass}_{DB,d}                //    (6)
                for (int d = 0; d < NUMDISP_; ++d)
                {
                  for (int k = d; k < NUMDISP_; ++k)
                  {  // symmetric matrix : only upper right triangle is computed
                    const int dk = NUMDISP_ * d - d * (d + 1) / 2 + k;
                    double defgradbybydisp_dk = 0.0;
                    for (int B = 0; B < NUMDIM_; ++B)
                    {
                      for (int D = 0; D < NUMDIM_; ++D)
                      {
                        const int DB = VOIGT3X3SYM_[NUMDIM_ * D + B];
                        const double DBfact = (D == B) ? 1.0 : 0.5;
                        for (int C = 0; C < NUMDIM_; ++C)
                        {
                          const int BC = VOIGT3X3_[NUMDIM_ * B + C];
                          const int CD = VOIGT3X3SYM_[NUMDIM_ * C + D];
                          const double CDfact = (C == D) ? 1.0 : 0.5;
                          defgradbybydisp_dk += invdefgradtimesboplin(BC, d) * CDfact *
                                                    invrgtstrDbydisp(CD, k) * rgtstr(D, B)  // (1)
                                                + invdefgradtimesboplin(BC, d) * invrgtstrD(C, D) *
                                                      DBfact * rgtstrbydisp(DB, k);  // (2)
                          defgradbybydisp_dk += invdefgradtimesboplin(BC, k) * CDfact *
                                                    invrgtstrDbydisp(CD, d) * rgtstr(D, B)  // (3)
                                                + invdefgradtimesboplin(BC, k) * invrgtstrD(C, D) *
                                                      DBfact * rgtstrbydisp(DB, d);  // (4)
                          defgradbybydisp_dk +=
                              invdefgradtimesdefgradD(B, C) * CDfact * invrgtstrDbydisp(CD, d) *
                                  DBfact * rgtstrbydisp(DB, k)  // (5)
                              + invdefgradtimesdefgradD(B, C) * CDfact * invrgtstrDbydisp(CD, k) *
                                    DBfact * rgtstrbydisp(DB, d);  // (6)
                          defgradbybydisp_dk += invdefgradtimesdefgradD(B, C) * CDfact *
                                                invrgtstrDbybydisp(CD, dk) * rgtstr(D, B);  // (7)
                        }
                        defgradbybydisp_dk += invdefgradtimesdefgradDtimesinvrgtstrD(B, D) *
                                              DBfact * rgtstrbybydisp(DB, dk);  // (8)
                      }
                    }
                    (*stiffmatrix)(d, k) -= defgradbybydisp_dk * effpressure * vol_w;
                    if (k != d) (*stiffmatrix)(k, d) -= defgradbybydisp_dk * effpressure * vol_w;
                  }
                }

              }  // if (lin_ >= lin_one)
            }    //  if (lin_ > lin_half)
          }      // end block

        }  // if (ans_ == ans_none) else

        // contribute stuff containing first derivatives in displacements
        if (stab_ != stab_puredisp)
        {
          // AUX = (WmT + fvT*fvT') * dFv
          LINALG::Matrix<NUMDFGR_, NUMDISP_> aux;
          aux.MultiplyNN(WmT, defgradbydisp);
          // K -= dFv' * AUX * detJ * w(gp)
          stiffmatrix->MultiplyTN(-effpressure * vol_w, defgradbydisp, aux, 1.0);
        }

        // deformation gradient differentiated with respect to displacements
        // double contracted with transposed inverted deformation gradient
        // (dFv'*fvT) = dFv'*fvT
        LINALG::Matrix<NUMDISP_, 1> defgradbydisptimestinvdefgrad;
        if ((gradmatrix != NULL) or (dargmatrix != NULL))
          defgradbydisptimestinvdefgrad.MultiplyTN(defgradbydisp, tinvdefgrad);

        // derivative of pressure contribution in internal force with respect to pressure
        // (-G) = -(dFv'*fvT)*pN * Fdet*detJ*wp(i);
        if (gradmatrix != NULL)
        {
          // contribute to G-op
          gradmatrix->MultiplyNT(-vol_w, defgradbydisptimestinvdefgrad, prshfct, 1.0);
        }

        // derivatives of spatially evaluated stabilisation matrix
        if ((stab_ == stab_spatial) or (stab_ == stab_spatialaffine))
        {
          // Mem_,d . ep
          // (h.p)_,d
          stabHHbydisp->MultiplyNT(
              pressure * vol_w, shapefcts[gp], defgradbydisptimestinvdefgrad, 1.0);
          // Eem_,d
          // wT_,d
          stabHAbydisp->MultiplyNT(
              stabA[gp](0, 0) * vol_w, shapefcts[gp], defgradbydisptimestinvdefgrad, 1.0);
          // Dem_,d
          // a_,d
          stabAAbydisp->MultiplyTT(
              stabA[gp](0, 0) * vol_w, stabA[gp], defgradbydisptimestinvdefgrad, 1.0);
          // J_,d
          detdefgradbydisp->UpdateT(stabA[gp](0, 0) * vol_w, defgradbydisptimestinvdefgrad, 1.0);
        }

      }  // end block
    }    // if (stiffmatrix != NULL)

    // mass matrix
    if (massmatrix != NULL)
    {
      double density = Material()->Density(gp);
      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      for (int inod = 0; inod < NUMNOD_; ++inod)
      {
        const double ifactor = shapefcts[gp](inod) * factor;
        for (int jnod = 0; jnod < NUMNOD_; ++jnod)
        {
          const double massfactor = shapefcts[gp](jnod) * ifactor;  // intermediate factor
          (*massmatrix)(NUMDIM_ * inod + 0, NUMDIM_ * jnod + 0) += massfactor;
          (*massmatrix)(NUMDIM_ * inod + 1, NUMDIM_ * jnod + 1) += massfactor;
          (*massmatrix)(NUMDIM_ * inod + 2, NUMDIM_ * jnod + 2) += massfactor;
        }
      }
    }  // end of mass matrix

    // store volume
    if (volume != NULL) *volume += vol_w;

    /* =========================================================================*/
  } /* ==================================================== end of Loop over GP */
  /* =========================================================================*/

  // finish with internal force etc
  // internal force
  if ((force != NULL) and (stabmatrix != NULL))
  {
    // integrate internal force vector
    // fint := fint
    //      + (B^T . sigma) * detJ * w(gp)   // already done
    //      + (-G) . ep
    if (stab_ != stab_puredisp) force->MultiplyNN(1.0, *gradmatrix, pres, 1.0);
  }
  // incompressiblity equation
  if ((incomp != NULL) and (stabmatrix != NULL))
  {
    // pint := pint
    //       - H . (Fdet - 1.0) * detJ * wp(gp)  // already done
    //       + (-Ce) . ep
    if (stab_ != stab_puredisp)
      incomp->MultiplyNN(1.0, *stabmatrix, pres, 1.0);
    else
      incomp->Clear();
  }

  // copy grad-matrix onto drag-matrix
  // (-G^T) = -pN * fvT'*dFv * Fdet*detJ*wp(i)
  if ((dargmatrix != NULL) and (gradmatrix != NULL))
  {
    // copy ordinary part from gradmatrix by transposing
    dargmatrix->UpdateT(*gradmatrix);
  }

  // derivative of incompressibility residual with respect to displacements
  // (-G^T) = -pN * fvT'*dFv * Fdet*detJ*wp(i)  // = gradmatrix^T // already done
  //        -pN_,d * fvT'*dFv * Fdet*detJ*wp(i)  // done here (k_ss)
  //        + ((-Ce) . ep)_{,d}  // done here (k_s)
  if ((stab_ == stab_spatial) or (stab_ == stab_spatialaffine))
  {
    // displacement derivative of spatially described stabilisation matrix (k_ss)
    if (dargmatrix != NULL)
    {
      // shear modulus
      const double shearmod = ShearMod();
      // integral pressure
      // (w.p)
      const double intpressure = stabHA.Dot(pres);
      // (-Cem . ep)_{,d} = -1./shearmod*( (Mem . ep)_{,d} - (Eem.inv(Dem).Eem')_{,d} . ep )
      // which is
      // -1./shearmod * (Mem . ep)_{,d}
      // -1/shearmod * (h . p)_,d
      dargmatrix->Update(-1.0 / shearmod, *stabHHbydisp, 1.0);
      // -1./shearmod * Eem_{,d} . (-inv(Dem).Eem'.ep)
      // +1/shearmod * a^{-1} * (w.p) * w^T_,d
      dargmatrix->Update(+1.0 / shearmod / stabAA(0, 0) * intpressure, *stabHAbydisp, 1.0);
      // Eem . inv(Dem)_,d = Eem . (-inv(Dem)^{-2}) . Dem_,d
      // -1./shearmod * (Eem_ . inv(Dem)_,d) . (-Eem' . ep)
      // +1/shearmod * (w.p) * (-a^{-2}) * w^T . a_,d
      dargmatrix->MultiplyNN(
          -intpressure / (shearmod * stabAA(0, 0) * stabAA(0, 0)), stabHA, *stabAAbydisp, 1.0);
      // (ep^T . Eem_{,d})^T = Eem_{,d}^T . ep
      // (p^T w^T_,d)
      LINALG::Matrix<1, NUMDISP_> prestimesstabHAbydisp;
      prestimesstabHAbydisp.MultiplyTN(pres, *stabHAbydisp);
      // -1./shearmod * (-Eem . inv(Dem)) . (ep^T . Eem_{,d})^T
      // +1/shearmod * a^{-1} * w^T . (p^T . w^T_,d)
      dargmatrix->MultiplyNN(+1.0 / (shearmod * stabAA(0, 0)), stabHA, prestimesstabHAbydisp, 1.0);
      //
      if (stab_ == stab_spatial)
      {
        if (stiffmatrix != NULL)
        {
          // -(w.p) * (-a^{-2}) * J_,d . a_,d
          stiffmatrix->MultiplyTN(
              intpressure / (stabAA(0, 0) * stabAA(0, 0)), *detdefgradbydisp, *stabAAbydisp, 1.0);
          // -a^{-1} * J_,d . w^t_,d
          stiffmatrix->MultiplyTN(
              -1.0 / stabAA(0, 0), *detdefgradbydisp, prestimesstabHAbydisp, 1.0);
        }
      }
    }
  }
  // displacement derivative of projected pressure shape functions (k_s)
  if (stab_ == stab_spatial)
  {
    if (dargmatrix != NULL)
    {
      // int((1-J)*A) * ( -a^{-2} * w^T . a_,d )
      dargmatrix->MultiplyNN(-voldev / (stabAA(0, 0) * stabAA(0, 0)), stabHA, *stabAAbydisp, 1.0);
      // int((1-J)*A) * ( a^{-1} * w^T_,d )
      dargmatrix->Update(voldev / stabAA(0, 0), *stabHAbydisp, 1.0);
    }
  }

  // EAS static condensation
  // subtract EAS matrices from disp-based Kdd to "soften" element
  if (eastype_ == soh8_easnone)
  {
    ;  // continue
  }
  else if (eastype_ == soh8_eassosh8)
  {
    EasCondensation<NUMEAS_SOSH8_>(force, stiffmatrix, gradmatrix, incomp, dargmatrix, stabmatrix,
        oldfeas, oldKaainv, oldKad, oldKap, feas, Kaa, Kad, Kap);
  }
  else if (eastype_ == soh8_easa)
  {
    EasCondensation<NUMEAS_A_>(force, stiffmatrix, gradmatrix, incomp, dargmatrix, stabmatrix,
        oldfeas, oldKaainv, oldKad, oldKap, feas, Kaa, Kad, Kap);
  }
  else
  {
    dserror("EAS ain't easy.");
  }

  // fake pure-disp based approach (ANS might be active)
  if (stab_ == stab_puredisp)
  {
    if (gradmatrix != NULL) gradmatrix->Clear();
    if (dargmatrix != NULL) dargmatrix->Clear();
    if (stabmatrix != NULL)
    {
      stabmatrix->Clear();
      for (int i = 0; i < NUMPRES_; ++i) (*stabmatrix)(i, i) = 1.0;
    }
  }

  // get away from here
  return;
}  // DRT::ELEMENTS::So_sh8p8::ForceStiffMass


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Stress(LINALG::Matrix<NUMGPT_, MAT::NUM_STRESS_3D>* elestress,
    const INPAR::STR::StressType iostress, const int gp, const double& detdefgrd,
    const LINALG::Matrix<NUMDIM_, NUMDIM_>& defgrd,
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& glstrain,
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& stress, const double& pressure)
{
  switch (iostress)
  {
    case INPAR::STR::stress_2pk:  // 2nd Piola-Kirchhoff stress
    {
      if (elestress == NULL) dserror("stress data not available");
      // determine stress
      if (stab_ == stab_puredisp)
      {
        // store stress
        for (int i = 0; i < MAT::NUM_STRESS_3D; ++i) (*elestress)(gp, i) = stress(i);
      }
      else
      {
        // inverted right Cauchy-Green strain tensor
        LINALG::Matrix<NUMDIM_, NUMDIM_> invcg;
        invcg.MultiplyTN(defgrd, defgrd);
        invcg.Invert();
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1> invcgv;
        Matrix2TensorToVector6Voigt(invcgv, invcg, voigt6_stress);
        // store stress
        for (int i = 0; i < MAT::NUM_STRESS_3D; ++i)
        {
          (*elestress)(gp, i) = stress(i) - pressure * detdefgrd * invcgv(i);
        }
      }
    }
    break;
    case INPAR::STR::stress_cauchy:  // true/Cauchy stress
    {
      if (elestress == NULL) dserror("stress data not available");
      // push forward
      LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> defgraddefgradT;
      Matrix2TensorToLeftRightProductMatrix6x6Voigt(
          defgraddefgradT, defgrd, true, voigt6_stress, voigt6_strain);
      // (deviatoric/isochoric) Cauchy stress vector
      LINALG::Matrix<MAT::NUM_STRESS_3D, 1> cauchyv;
      cauchyv.MultiplyNN(1.0 / detdefgrd, defgraddefgradT, stress);

      // determine total stress
      if (stab_ != stab_puredisp)
      {
        // above computed #cauchyv is deviatoric/isochoric true stress
        // isochoric Cauchy stress vector
        // the volumetric part is added here
        for (int i = 0; i < NUMDIM_; ++i) cauchyv(i) -= pressure;
      }

      // store stress
      for (int i = 0; i < MAT::NUM_STRESS_3D; ++i) (*elestress)(gp, i) = cauchyv(i);
    }
    break;
    case INPAR::STR::stress_none:
      break;
    default:
    {
      dserror("requested stress option not available");
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Strain(
    LINALG::Matrix<NUMGPT_, MAT::NUM_STRESS_3D>* elestrain,  ///< store the strain herein
    const INPAR::STR::StrainType iostrain,
    const int gp,             ///< Gauss point index
    const double& detdefgrd,  ///< determinant of (assumed) deformation gradient
    const LINALG::Matrix<NUMDIM_, NUMDIM_>& defgrd,     ///< (assumed) deformation gradient
    const LINALG::Matrix<NUMDIM_, NUMDIM_>& invdefgrd,  ///< (assumed) inverted deformation gradient
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& glstrain  ///< Green-Lagrange strain vector
)
{
  switch (iostrain)
  {
    case INPAR::STR::strain_gl:  // Green-Lagrange strain
    {
      if (elestrain == NULL) dserror("strain data not available");
      // store
      for (int i = 0; i < NUMDIM_; ++i) (*elestrain)(gp, i) = glstrain(i);
      for (int i = NUMDIM_; i < MAT::NUM_STRESS_3D; ++i) (*elestrain)(gp, i) = 0.5 * glstrain(i);
    }
    break;
    case INPAR::STR::strain_ea:  // Euler-Almansi strain
    {
      if (elestrain == NULL) dserror("strain data not available");
      // create push forward 6x6 matrix
      LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> invdefgradTdefgrad;
      Matrix2TensorToLeftRightProductMatrix6x6Voigt(
          invdefgradTdefgrad, invdefgrd, false, voigt6_strain, voigt6_stress);
      // push forward
      LINALG::Matrix<MAT::NUM_STRESS_3D, 1> eastrain;
      eastrain.MultiplyNN(invdefgradTdefgrad, glstrain);
      // store
      for (int i = 0; i < NUMDIM_; ++i) (*elestrain)(gp, i) = eastrain(i);
      for (int i = NUMDIM_; i < MAT::NUM_STRESS_3D; ++i) (*elestrain)(gp, i) = 0.5 * eastrain(i);
    }
    break;
    case INPAR::STR::strain_none:
      break;
    default:
    {
      dserror("requested strain option not available");
      break;
    }
  }

  // bye
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::AssDefGrad(double& detdefgrad,
    LINALG::Matrix<NUMDIM_, NUMDIM_>& defgrad, LINALG::Matrix<NUMDIM_, NUMDIM_>& invdefgrad,
    LINALG::Matrix<NUMDIM_, NUMDIM_>& rgtstr, LINALG::Matrix<NUMDIM_, NUMDIM_>& defgradD,
    LINALG::Matrix<NUMDIM_, NUMDIM_>& rgtstrD, LINALG::Matrix<NUMDIM_, NUMDIM_>& invrgtstrD,
    const LINALG::Matrix<NUMDIM_, NUMDIM_>& Jinv, const LINALG::Matrix<NUMDIM_, NUMDIM_>& Jac,
    const LINALG::Matrix<NUMDIM_, NUMDIM_>& jac,
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& glstrain)
{
  // pure displacement-based deformation gradient
  // F = x_{,X} = x_{,xi} . xi_{,X} = x_{,xi} . (X_{,xi})^{-1} = jac^T . Jinv^T
  defgradD.MultiplyTT(jac, Jinv);

  // pure displacement-based right Cauchy-Green strain
  LINALG::Matrix<NUMDIM_, NUMDIM_> cgD;
  cgD.MultiplyTN(defgradD, defgradD);

  // rotation matrix in pure displacement based deformation gradient
  // and pure disp-based material stretch tensor
  LINALG::Matrix<NUMDIM_, NUMDIM_> rot;
  {
#if 0
    StretchTensor(NULL,&rgtstrD,&invrgtstrD,cgD);
    // rotation matrix
    rot.MultiplyNN(defgradD,invrgtstrD);
#else
#if 0
    LINALG::Matrix<NUMDIM_,NUMDIM_> nd;
    LINALG::Matrix<NUMDIM_,NUMDIM_> lamd;
    LINALG::Matrix<NUMDIM_,NUMDIM_> NdT;
    LINALG::SVD(defgradD,nd,lamd,NdT);
    rot.MultiplyNN(nd,NdT);
    // pure disp-based material stretch tensor
    LINALG::Matrix<NUMDIM_,NUMDIM_> aux;
    aux.MultiplyTN(NdT,lamd);
    rgtstrD.MultiplyNN(aux,NdT);
#else
    // spectral decomposition of disp-based right Cauchy-Green tensor
    LINALG::Matrix<NUMDIM_, NUMDIM_> NdT;
    LINALG::Matrix<NUMDIM_, NUMDIM_> lamd;
#if 0
#if 0
    const int err = SymSpectralDecompJacIter(lamd,NdT,cgD,EPS12,12);
    if (err!=0) dserror("spectral decomposition failed");
    // spectral composition of disp-based right stretch tensor
    for (int i=0; i<NUMDIM_; ++i) lamd(i,i) = sqrt(lamd(i,i));
    LINALG::Matrix<NUMDIM_,NUMDIM_> aux;
    aux.MultiplyNN(NdT,lamd);
    rgtstrD.MultiplyNT(aux,NdT);
#else
    LINALG::Matrix<NUMDIM_,NUMDIM_> Nd;
    LINALG::SVD(cgD,NdT,lamd,Nd);
    // spectral composition of disp-based right stretch tensor
    for (int i=0; i<NUMDIM_; ++i) lamd(i,i) = sqrt(lamd(i,i));
    LINALG::Matrix<NUMDIM_,NUMDIM_> aux;
    aux.MultiplyNN(NdT,lamd);
    rgtstrD.MultiplyNN(aux,Nd);
#endif
#else
    LINALG::SYEV(cgD, lamd, NdT);
    // spectral composition of disp-based right stretch tensor
    for (int al = 0; al < NUMDIM_; ++al) lamd(al, al) = sqrt(lamd(al, al));
    LINALG::Matrix<NUMDIM_, NUMDIM_> aux;
    aux.MultiplyNN(NdT, lamd);
    rgtstrD.MultiplyNT(aux, NdT);
#endif
    // inverse disp-based right stretch tensor
#if 0
    invrgtstrD.Update(rgtstrD);
    const double detrgtstrD = invrgtstrD.Invert();
    if (detrgtstrD < 0.0) dserror("Trouble during inversion of right stretch tensor");
#else
    invrgtstrD.Clear();
    // double detdefgradD = 1.0;
    for (int al = 0; al < NUMDIM_; ++al)
    {
      // detdefgradD *= lamd(al,al);
      for (int j = 0; j < NUMDIM_; ++j)
      {
        const double NdT_jal_by_lamd_alal = NdT(j, al) / lamd(al, al);
        for (int i = 0; i < NUMDIM_; ++i) invrgtstrD(i, j) += NdT(i, al) * NdT_jal_by_lamd_alal;
      }
    }
#endif
    // rotation matrix
    rot.MultiplyNN(defgradD, invrgtstrD);
#endif
#endif

    // // U_{,C}
    //     // correct, but same speed like solution with Lapack and inaccurate
    //     {
    //       for (int kl=0; kl<MAT::NUM_STRESS_3D; ++kl) {
    //         const int k = VOIGT6ROW_[kl];
    //         const int l = VOIGT6COL_[kl];
    //         for (int ij=0; ij<MAT::NUM_STRESS_3D; ++ij) {
    //           const int i = VOIGT6ROW_[ij];
    //           const int j = VOIGT6COL_[ij];
    //           double rgtstrDbyrcgD_ijkl = 0.0;
    //           for (int al=0; al<NUMDIM_; ++al) {
    //             for (int be=0; be<NUMDIM_; ++be) {
    //               rgtstrDbyrcgD_ijkl
    //                 += 0.5*( NdT(i,al)*NdT(k,al)*NdT(l,be)*NdT(j,be)
    //                          + NdT(j,al)*NdT(l,al)*NdT(i,be)*NdT(k,be)
    //                   ) / (2.0*lamd(be,be));
    //               if (ij >= NUMDIM_)
    //                 rgtstrDbyrcgD_ijkl
    //                   +=  0.5*( NdT(j,al)*NdT(k,al)*NdT(l,be)*NdT(i,be)
    //                             + NdT(i,al)*NdT(l,al)*NdT(j,be)*NdT(k,be)
    //                     ) / (2.0*lamd(be,be));
    //             }
    //           }
    //           rgtstrDbyrcgD(ij,kl) = rgtstrDbyrcgD_ijkl;
    //         }
    //       }
    //     }
  }

  // assumed material stretch tensor
  LINALG::Matrix<NUMDIM_, NUMDIM_> invrgtstr;
  {
    LINALG::Matrix<NUMDIM_, NUMDIM_> cga;
    {
      for (int i = 0; i < NUMDIM_; ++i) cga(i, i) = 2.0 * glstrain(i, 0) + 1.0;
      // off-diagonal terms are already twice in the Voigt-GLstrain-vector
      cga(0, 1) = glstrain(3);
      cga(1, 0) = glstrain(3);
      cga(1, 2) = glstrain(4);
      cga(2, 1) = glstrain(4);
      cga(0, 2) = glstrain(5);
      cga(2, 0) = glstrain(5);
    }
    LINALG::Matrix<NUMDIM_, NUMDIM_> lama;
    LINALG::Matrix<NUMDIM_, NUMDIM_> NaT;
#if 0
    LINALG::Matrix<NUMDIM_,NUMDIM_> Na;
    LINALG::SVD(cga,NaT,lama,Na);
    for (int al=0; al<NUMDIM_; ++al) lama(al,al) = sqrt(lama(al,al));
    LINALG::Matrix<NUMDIM_,NUMDIM_> aux;
    aux.MultiplyNN(NaT,lama);
    rgtstr.MultiplyNN(aux,Na);
#else
    LINALG::SYEV(cga, lama, NaT);
    for (int al = 0; al < NUMDIM_; ++al) lama(al, al) = sqrt(lama(al, al));
#if 0
    rgtstr.Clear();
    for (int al=0; al<NUMDIM_; ++al) {
      for (int j=0; j<NUMDIM_; ++j) {
        const double NaT_jal_times_lama_alal = NaT(j,al)*lama(al,al);
        for (int i=0; i<NUMDIM_; ++i)
          rgtstr(i,j) += NaT(i,al)*NaT_jal_times_lama_alal;
      }
    }
#else
    LINALG::Matrix<NUMDIM_, NUMDIM_> aux;
    aux.MultiplyNN(NaT, lama);
    rgtstr.MultiplyNT(aux, NaT);
#endif
#endif
#if 0
    invrgtstr.Update(rgtstr);
    detdefgrad = invrgtstr.Invert();
#else
    invrgtstr.Clear();
    detdefgrad = 1.0;
    for (int al = 0; al < NUMDIM_; ++al)
    {
      detdefgrad *= lama(al, al);
      for (int j = 0; j < NUMDIM_; ++j)
      {
        const double NaT_jal_by_lama_alal = NaT(j, al) / lama(al, al);
        for (int i = 0; i < NUMDIM_; ++i) invrgtstr(i, j) += NaT(i, al) * NaT_jal_by_lama_alal;
      }
    }
#endif
  }

  // assumed deformation gradient
  defgrad.MultiplyNN(rot, rgtstr);

  // inverse of assumed deformation gradient
  invdefgrad.MultiplyNT(invrgtstr, rot);

  // done
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double DRT::ELEMENTS::So_sh8p8::ShearMod() const
{
  // All materials that have a pure LINALG::Matrix
  // interface go to the material law here.
  // the old interface does not exist anymore...
  Teuchos::RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case INPAR::MAT::m_stvenant: /*-------- st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast<MAT::StVenantKirchhoff*>(mat.get());
      return stvk->ShearMod();
      break;
    }
    case INPAR::MAT::m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast<MAT::NeoHooke*>(mat.get());
      return neo->ShearMod();
      break;
    }
    case INPAR::MAT::m_aaaneohooke: /*-----------AAA NeoHookean Material */
    {
      MAT::AAAneohooke* aaaneo = static_cast<MAT::AAAneohooke*>(mat.get());
      return aaaneo->ShearMod();
      break;
    }
    case INPAR::MAT::m_aaaraghavanvorp_damage: /*-AAA RaghavanVorp Material with damage*/
    {
      MAT::AAAraghavanvorp_damage* aaadam = static_cast<MAT::AAAraghavanvorp_damage*>(mat.get());
      return aaadam->ShearMod();
      break;
    }
    case INPAR::MAT::m_viscoanisotropic: /*-----------Anisotropic Viscous Material */
    {
      MAT::ViscoAnisotropic* visco = static_cast<MAT::ViscoAnisotropic*>(Material().get());
      return visco->ShearMod();
      break;
    }
    case INPAR::MAT::m_visconeohooke: /*-----------Viscous NeoHookean Material */
    {
      MAT::ViscoNeoHooke* visconeo = static_cast<MAT::ViscoNeoHooke*>(Material().get());
      return visconeo->ShearMod();
      break;
    }
    case INPAR::MAT::m_viscoelasthyper: /*-----------general visco-hyperelastic material */
    {
      MAT::ViscoElastHyper* viscoelasthyper = static_cast<MAT::ViscoElastHyper*>(Material().get());
      return viscoelasthyper->ShearMod();
      break;
    }
    case INPAR::MAT::m_elasthyper: /*----------- general hyperelastic matrial */
    {
      MAT::ElastHyper* hyper = static_cast<MAT::ElastHyper*>(Material().get());
      return hyper->ShearMod();
      break;
    }
    default:
      dserror("Cannot ask material for shear modulus");
      break;
  }  // switch (mat->MaterialType())

  return 0;
}

void DRT::ELEMENTS::So_sh8p8::CalcSTCMatrix(LINALG::Matrix<NUMDOF_, NUMDOF_>& elemat1,
    const INPAR::STR::STC_Scale stc_scaling, const int stc_layer, std::vector<int>& lm,
    DRT::Discretization& discretization, bool calcinverse)
{
  double stc_fact = 0.0;
  if (stc_scaling == INPAR::STR::stc_currsym)
  {
    // stc_fact = sqrt(sosh8_calcaspectratio());
    stc_fact = sosh8_calcaspectratio();
  }
  else
  {
    // stc_fact = sosh8_calcaspectratio();
    stc_fact = sosh8_calcaspectratio() * sosh8_calcaspectratio();
  }

  // Compute different scaling factors for STC or Inv(STC)
  const double factor1 = (stc_fact + 1.0) / (2.0 * stc_fact);
  const double factor2 = (stc_fact - 1.0) / (2.0 * stc_fact);
  const double factor3 = (1.0 / stc_fact);
  const double factor4 = (1.0 - 1.0 / stc_fact);

  if (stc_scaling == INPAR::STR::stc_curr or stc_scaling == INPAR::STR::stc_currsym)
  {
    LINALG::Matrix<NUMDOF_, 1> adjele(true);
    DRT::Node** nodes = Nodes();

    std::vector<DRT::Condition*> cond0;
    int condnum0 = 1000;    // minimun STCid of layer with nodes 0..3
    bool current0 = false;  // layer with nodes 0..4 to be scaled
    (nodes[0])->GetCondition("STC Layer", cond0);
    std::vector<DRT::Condition*> cond1;
    int condnum1 = 1000;    // minimun STCid of layer with nodes 4..7
    bool current1 = false;  // minimun STCid of layer with nodes 4..7
    (nodes[NUMNOD_ / 2])->GetCondition("STC Layer", cond1);

    for (unsigned int conu = 0; conu < cond0.size(); ++conu)
    {
      int tmp = cond0[conu]->GetInt("ConditionID");
      if (tmp < condnum0) condnum0 = tmp;
    }
    if (condnum0 == stc_layer) current0 = true;


    for (unsigned int conu = 0; conu < cond1.size(); ++conu)
    {
      int tmp = cond1[conu]->GetInt("ConditionID");
      if (tmp < condnum1) condnum1 = tmp;
    }
    if (condnum1 == stc_layer) current1 = true;

    // both surfaces are to be scaled
    if (current0 and current1)
    {
      // only valid for first round
      if (condnum0 != 1)
        dserror("STC error: non-initial layer is not connected to a smaller id");
      else
      {
        for (int i = 0; i < NUMNOD_; i++)
        {
          adjele(NODDOF_ * i + 0, 0) = nodes[i]->NumElement();
          adjele(NODDOF_ * i + 1, 0) = nodes[i]->NumElement();
          adjele(NODDOF_ * i + 2, 0) = nodes[i]->NumElement();
          adjele(NODDOF_ * i + 3, 0) = nodes[i]->NumElement();
        }

        for (int ind1 = 0; ind1 < NUMNOD_ / 2; ind1++)
        {
          for (int ind2 = 0; ind2 < NUMDIM_; ind2++)
          {
            elemat1(NODDOF_ * ind1 + ind2, NODDOF_ * ind1 + ind2) +=
                factor1 / adjele(NODDOF_ * ind1 + ind2, 0) * cond0.size();
            elemat1(NODDOF_ * ind1 + ind2 + NUMDOF_ / 2, NODDOF_ * ind1 + ind2 + NUMDOF_ / 2) +=
                factor1 / adjele(NODDOF_ * ind1 + ind2 + NUMDOF_ / 2, 0) * cond1.size();
            elemat1(NODDOF_ * ind1 + ind2, NODDOF_ * ind1 + ind2 + NUMDOF_ / 2) +=
                factor2 / adjele(NODDOF_ * ind1 + ind2, 0) * cond0.size();
            elemat1(NODDOF_ * ind1 + ind2 + NUMDOF_ / 2, NODDOF_ * ind1 + ind2) +=
                factor2 / adjele(NODDOF_ * ind1 + ind2 + NUMDOF_ / 2, 0) * cond1.size();
          }
          elemat1(NODDOF_ * ind1 + NUMDIM_, NODDOF_ * ind1 + NUMDIM_) +=
              1.0 / adjele(NODDOF_ * ind1 + NUMDIM_, 0) * cond0.size();
          elemat1(NODDOF_ * ind1 + NUMDIM_ + NUMDOF_ / 2, NODDOF_ * ind1 + NUMDIM_ + NUMDOF_ / 2) +=
              1.0 / adjele(NODDOF_ * ind1 + NUMDIM_ + NUMDOF_ / 2, 0) * cond1.size();
        }
      }
    }
    // surface with nodes 0..3 is to be scaled
    else if (current0)
    {
      // but not by this element
      if (condnum1 > condnum0)
      {
        for (int i = 0; i < NUMNOD_; i++)
        {
          adjele(NODDOF_ * i + 0, 0) = nodes[i]->NumElement();
          adjele(NODDOF_ * i + 1, 0) = nodes[i]->NumElement();
          adjele(NODDOF_ * i + 2, 0) = nodes[i]->NumElement();
          adjele(NODDOF_ * i + 3, 0) = nodes[i]->NumElement();
        }

        for (int ind1 = NUMNOD_ / 2; ind1 < NUMNOD_; ind1++)
        {
          for (int ind2 = 0; ind2 < NODDOF_; ind2++)
          {
            elemat1(NODDOF_ * ind1 + ind2, NODDOF_ * ind1 + ind2) +=
                1.0 / adjele(NODDOF_ * ind1 + ind2, 0);
          }
        }
      }
      // this element has to do the whole scaling
      else if (condnum1 < condnum0)
      {
        for (int i = 0; i < NUMNOD_; i++)
        {
          adjele(NODDOF_ * i + 0, 0) = nodes[i]->NumElement();
          adjele(NODDOF_ * i + 1, 0) = nodes[i]->NumElement();
          adjele(NODDOF_ * i + 2, 0) = nodes[i]->NumElement();
          adjele(NODDOF_ * i + 3, 0) = nodes[i]->NumElement();
        }

        for (int ind1 = 0; ind1 < NUMNOD_; ind1++)
        {
          for (int ind2 = 0; ind2 < NUMDIM_; ind2++)
          {
            if (ind1 < NUMNOD_ / 2)
            {
              elemat1(NODDOF_ * ind1 + ind2, NODDOF_ * ind1 + ind2) +=
                  factor3 / adjele(NODDOF_ * ind1 + ind2, 0) * cond0.size();
              elemat1(NODDOF_ * ind1 + ind2, NODDOF_ * ind1 + ind2 + NUMDOF_ / 2) +=
                  factor4 / adjele(NODDOF_ * ind1 + ind2, 0) * cond0.size();
            }
            else
            {
              elemat1(NODDOF_ * ind1 + ind2, NODDOF_ * ind1 + ind2) +=
                  1.0 / adjele(NODDOF_ * ind1 + NUMDIM_, 0);
            }
          }
          elemat1(NODDOF_ * ind1 + NUMDIM_, NODDOF_ * ind1 + NUMDIM_) +=
              1.0 / adjele(NODDOF_ * ind1 + NUMDIM_, 0);
        }
      }
    }
    // surface with nodes 4..7 is to be scaled
    else if (current1)
    {
      // but not by this element
      if (condnum0 > condnum1)
      {
        for (int i = 0; i < NUMNOD_; i++)
        {
          adjele(NODDOF_ * i + 0, 0) = nodes[i]->NumElement();
          adjele(NODDOF_ * i + 1, 0) = nodes[i]->NumElement();
          adjele(NODDOF_ * i + 2, 0) = nodes[i]->NumElement();
          adjele(NODDOF_ * i + 3, 0) = nodes[i]->NumElement();
        }

        for (int ind1 = 0; ind1 < NUMNOD_ / 2; ind1++)
        {
          for (int ind2 = 0; ind2 < NODDOF_; ind2++)
          {
            elemat1(NODDOF_ * ind1 + ind2, NODDOF_ * ind1 + ind2) +=
                1.0 / adjele(NODDOF_ * ind1 + ind2, 0);
          }
        }
      }
      // this element has to do the whole scaling
      else if (condnum0 < condnum1)
      {
        for (int i = 0; i < NUMNOD_; i++)
        {
          adjele(NODDOF_ * i + 0, 0) = nodes[i]->NumElement();
          adjele(NODDOF_ * i + 1, 0) = nodes[i]->NumElement();
          adjele(NODDOF_ * i + 2, 0) = nodes[i]->NumElement();
          adjele(NODDOF_ * i + 3, 0) = nodes[i]->NumElement();
        }

        for (int ind1 = 0; ind1 < NUMNOD_; ind1++)
        {
          for (int ind2 = 0; ind2 < NUMDIM_; ind2++)
          {
            if (ind1 >= NUMNOD_ / 2)
            {
              elemat1(NODDOF_ * ind1 + ind2, NODDOF_ * ind1 + ind2) +=
                  factor3 / adjele(NODDOF_ * ind1 + ind2, 0) * cond1.size();
              elemat1(NODDOF_ * ind1 + ind2, NODDOF_ * ind1 + ind2 - NUMDOF_ / 2) +=
                  factor4 / adjele(NODDOF_ * ind1 + ind2, 0) * cond1.size();
            }
            else
            {
              elemat1(NODDOF_ * ind1 + ind2, NODDOF_ * ind1 + ind2) +=
                  1.0 / adjele(NODDOF_ * ind1 + NUMDIM_, 0);
            }
          }
          elemat1(NODDOF_ * ind1 + NUMDIM_, NODDOF_ * ind1 + NUMDIM_) +=
              1.0 / adjele(NODDOF_ * ind1 + NUMDIM_, 0);
        }
      }
    }
    // nothing to be done by this element
    else
    {
      for (int i = 0; i < NUMNOD_; i++)
      {
        adjele(NODDOF_ * i + 0, 0) = nodes[i]->NumElement();
        adjele(NODDOF_ * i + 1, 0) = nodes[i]->NumElement();
        adjele(NODDOF_ * i + 2, 0) = nodes[i]->NumElement();
        adjele(NODDOF_ * i + 3, 0) = nodes[i]->NumElement();
      }
      for (int ind1 = 0; ind1 < NUMDOF_; ind1++)
      {
        elemat1(ind1, ind1) += 1.0 / adjele(ind1, 0);
      }
    }
  }
  else
    dserror("Chosen STC_SCALING not supported!");
}


/*======================================================================*/
/*======================================================================*/
/*======================================================================*/
/*======================================================================*/


/*----------------------------------------------------------------------*
 |  init the element (public)                                  maf 07/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_sh8p8Type::Initialize(DRT::Discretization& dis)
{
  // sosh8_gmshplotdis(dis);

  int num_morphed_so_hex8_easnone = 0;

  // Loop through all elements
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    // get the actual element
    if (dis.lColElement(i)->ElementType() != *this) continue;
    // go on for So_sh8p elements
    DRT::ELEMENTS::So_sh8p8* actele = dynamic_cast<DRT::ELEMENTS::So_sh8p8*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_sh8p8* failed");

    if (!actele->nodes_rearranged_)
    {
      bool altered = false;
      switch (actele->thickdir_)
      {
        // check for automatic definition of thickness direction
        case DRT::ELEMENTS::So_sh8::autoj:
        {
          actele->thickdir_ = actele->sosh8_findthickdir();
          altered = true;
          break;
        }
        // check for enforced definition of thickness direction
        case DRT::ELEMENTS::So_sh8::globx:
        {
          LINALG::Matrix<DRT::ELEMENTS::So_sh8p8::NUMDIM_, 1> thickdirglo(true);
          thickdirglo(0) = 1.0;
          actele->thickdir_ = actele->sosh8_enfthickdir(thickdirglo);
          altered = true;
          break;
        }
        case DRT::ELEMENTS::So_sh8::globy:
        {
          LINALG::Matrix<DRT::ELEMENTS::So_sh8p8::NUMDIM_, 1> thickdirglo(true);
          thickdirglo(1) = 1.0;
          actele->thickdir_ = actele->sosh8_enfthickdir(thickdirglo);
          altered = true;
          break;
        }
        case DRT::ELEMENTS::So_sh8::globz:
        {
          LINALG::Matrix<DRT::ELEMENTS::So_sh8p8::NUMDIM_, 1> thickdirglo(true);
          thickdirglo(2) = 1.0;
          actele->thickdir_ = actele->sosh8_enfthickdir(thickdirglo);
          altered = true;
          break;
        }
        default:
          break;
      }

      if (altered and (actele->thickdir_ != DRT::ELEMENTS::So_sh8p8::undefined))
      {
        // special element-dependent input of material parameters
        if (actele->Material()->MaterialType() == INPAR::MAT::m_viscoanisotropic)
        {
          MAT::ViscoAnisotropic* visco =
              static_cast<MAT::ViscoAnisotropic*>(actele->Material().get());
          visco->Setup(DRT::ELEMENTS::So_sh8p8::NUMGPT_, actele->thickvec_);
          if (actele->thickvec_.size() == 0)
            dserror("zero size thickness vector for element %d", actele->Id());
        }
      }

      int new_nodeids[DRT::ELEMENTS::So_sh8p8::NUMNOD_];

      switch (actele->thickdir_)
      {
        case DRT::ELEMENTS::So_sh8::globx:
        case DRT::ELEMENTS::So_sh8::globy:
        case DRT::ELEMENTS::So_sh8::globz:
        {
          dserror("This should have been replaced by auto(r|s|t)");
          break;
        }
        case DRT::ELEMENTS::So_sh8::autor:
        case DRT::ELEMENTS::So_sh8::enfor:
        {
          // resorting of nodes to arrive at local t-dir for global x-dir
          new_nodeids[0] = actele->NodeIds()[7];
          new_nodeids[1] = actele->NodeIds()[4];
          new_nodeids[2] = actele->NodeIds()[0];
          new_nodeids[3] = actele->NodeIds()[3];
          new_nodeids[4] = actele->NodeIds()[6];
          new_nodeids[5] = actele->NodeIds()[5];
          new_nodeids[6] = actele->NodeIds()[1];
          new_nodeids[7] = actele->NodeIds()[2];
          //        actele->sosh8_gmshplotlabeledelement(actele->NodeIds());
          //        actele->sosh8_gmshplotlabeledelement(new_nodeids);
          actele->SetNodeIds(DRT::ELEMENTS::So_sh8p8::NUMNOD_, new_nodeids);
          actele->nodes_rearranged_ = true;
          break;
        }
        case DRT::ELEMENTS::So_sh8::autos:
        case DRT::ELEMENTS::So_sh8::enfos:
        {
          // resorting of nodes to arrive at local t-dir for global y-dir
          new_nodeids[0] = actele->NodeIds()[4];
          new_nodeids[1] = actele->NodeIds()[5];
          new_nodeids[2] = actele->NodeIds()[1];
          new_nodeids[3] = actele->NodeIds()[0];
          new_nodeids[4] = actele->NodeIds()[7];
          new_nodeids[5] = actele->NodeIds()[6];
          new_nodeids[6] = actele->NodeIds()[2];
          new_nodeids[7] = actele->NodeIds()[3];
          actele->SetNodeIds(DRT::ELEMENTS::So_sh8p8::NUMNOD_, new_nodeids);
          actele->nodes_rearranged_ = true;
          break;
        }
        case DRT::ELEMENTS::So_sh8::autot:
        case DRT::ELEMENTS::So_sh8::enfot:
        {
          // no resorting necessary
          for (int node = 0; node < 8; ++node)
          {
            new_nodeids[node] = actele->NodeIds()[node];
          }
          actele->SetNodeIds(DRT::ELEMENTS::So_sh8p8::NUMNOD_, new_nodeids);
          actele->nodes_rearranged_ = true;
          break;
        }
        case DRT::ELEMENTS::So_sh8::undefined:
        {
          // here comes plan B: switch off ANS
          actele->SetANS(DRT::ELEMENTS::So_sh8p8::ans_none);
          num_morphed_so_hex8_easnone++;
          break;
        }
        case DRT::ELEMENTS::So_sh8::none:
          break;
        default:
        {
          dserror("no thickness direction for So_sh8p8");
          break;
        }
      }
      // actele->sosh8p8_gmshplotlabeledelement(actele->NodeIds());
    }
  }

  if (num_morphed_so_hex8_easnone > 0)
  {
    std::cout << std::endl
              << num_morphed_so_hex8_easnone
              << " Sosh8p8-Elements have no clear 'thin' direction and ANS is disabled!"
              << std::endl;
  }

  // fill complete again to reconstruct element-node pointers,
  // but without element init, etc.
  dis.FillComplete(false, false, false);

  // loop again to init Jacobian for Sosh8p8's
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So_sh8p8* actele = dynamic_cast<DRT::ELEMENTS::So_sh8p8*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_sh8p8* failed");
    actele->InitJacobianMapping();  // this sets #invJ_ in So_hex8
  }

  // **************** debug printout ot gmesh **********************************
  // sosh8_gmshplotdis(dis);

  return 0;
}

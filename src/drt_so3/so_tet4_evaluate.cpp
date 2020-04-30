/*----------------------------------------------------------------------*/
/*! \file
\brief quadratic nonlinear tetrahedron
\level 1
\maintainer Christoph Meier
*----------------------------------------------------------------------*/
#include "so_tet4.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/voigt_notation.H"
#include "../linalg/linalg_utils_densematrix_inverse.H"
#include "../linalg/linalg_utils_densematrix_eigen.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_patspec/patspec.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/elasthyper.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_mat/constraintmixture.H"
#include "../drt_contact/contact_analytical.H"
#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <boost/static_assert.hpp>

// inverse design object
#include "inversedesign.H"
#include "prestress.H"

#include "../drt_structure_new/str_elements_paramsinterface.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/thermostvenantkirchhoff.H"
#include "../drt_mat/thermoplastichyperelast.H"
#include "../drt_mat/robinson.H"
#include "so_utils.H"

//#define PRINT_DEBUG
#ifdef PRINT_DEBUG
#include <string>
#include <sstream>
#include <cstd::string>
template <class T>
void writeArray(const T& mat, std::string name = "unnamed")
{
  std::stringstream header;
  header << 'M' << name << ':' << mat.M() << 'x' << mat.N() << ':';
  unsigned int s = header.str().size() + mat.M() * mat.N() * sizeof(double);
  std::cerr.write(reinterpret_cast<const char*>(&s), sizeof(unsigned int));
  std::cerr << header.str();
  for (int i = 0; i < mat.M() * mat.N(); ++i)
  {
    std::cerr.write(reinterpret_cast<const char*>(&(mat.A()[i])), sizeof(double));
  }
}

void writeComment(const std::string v)
{
  unsigned int s = v.size() + 1;
  std::cerr.write(reinterpret_cast<const char*>(&s), sizeof(unsigned int));
  std::cerr << 'C' << v;
}
#endif  // PRINT_DEBUG

using VoigtMapping = UTILS::VOIGT::IndexMappings;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              vlf 06/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet4::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  // Check whether the solid material PostSetup() routine has already been called and call it if not
  EnsureMaterialPostSetup(params);

  LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4> elemat1(elemat1_epetra.A(), true);
  LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4> elemat2(elemat2_epetra.A(), true);
  LINALG::Matrix<NUMDOF_SOTET4, 1> elevec1(elevec1_epetra.A(), true);
  LINALG::Matrix<NUMDOF_SOTET4, 1> elevec2(elevec2_epetra.A(), true);
  LINALG::Matrix<NUMDOF_SOTET4, 1> elevec3(elevec3_epetra.A(), true);

  // start with "none"
  DRT::ELEMENTS::So_tet4::ActionType act = So_tet4::none;

  // get the required action
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    dserror("No action supplied");
  else if (action == "calc_struct_linstiff")
    act = So_tet4::calc_struct_linstiff;
  else if (action == "calc_struct_nlnstiff")
    act = So_tet4::calc_struct_nlnstiff;
  else if (action == "calc_struct_internalforce")
    act = So_tet4::calc_struct_internalforce;
  else if (action == "calc_struct_linstiffmass")
    act = So_tet4::calc_struct_linstiffmass;
  else if (action == "calc_struct_nlnstiffmass")
    act = So_tet4::calc_struct_nlnstiffmass;
  else if (action == "calc_struct_nlnstifflmass")
    act = So_tet4::calc_struct_nlnstifflmass;
  else if (action == "calc_struct_stress")
    act = So_tet4::calc_struct_stress;
  else if (action == "postprocess_stress")
    act = So_tet4::postprocess_stress;
  else if (action == "calc_struct_eleload")
    act = So_tet4::calc_struct_eleload;
  else if (action == "calc_struct_fsiload")
    act = So_tet4::calc_struct_fsiload;
  else if (action == "calc_struct_store_istep")
    act = So_tet4::struct_calc_store_istep;
  else if (action == "calc_struct_recover_istep")
    act = So_tet4::struct_calc_recover_istep;
  else if (action == "calc_struct_update_istep")
    act = So_tet4::calc_struct_update_istep;
  else if (action == "calc_struct_reset_istep")
    act = So_tet4::calc_struct_reset_istep;
  else if (action == "calc_struct_reset_all")
    act = So_tet4::calc_struct_reset_all;
  else if (action == "calc_struct_errornorms")
    act = So_tet4::calc_struct_errornorms;
  else if (action == "calc_struct_prestress_update")
    act = So_tet4::prestress_update;
  else if (action == "calc_global_gpstresses_map")
    act = So_tet4::calc_global_gpstresses_map;
  else if (action == "calc_struct_energy")
    act = So_tet4::calc_struct_energy;
  else if (action == "calc_struct_output_E")
    act = So_tet4::calc_struct_output_E;
  else if (action == "calc_struct_inversedesign_update")
    act = So_tet4::inversedesign_update;
  else if (action == "calc_struct_inversedesign_switch")
    act = So_tet4::inversedesign_switch;
  else if (action == "multi_calc_dens")
    act = So_tet4::multi_calc_dens;
  else if (action == "multi_readrestart")
    act = So_tet4::multi_readrestart;
  else if (action == "calc_struct_recover")
    return 0;
  else if (action == "calc_struct_predict")
    return 0;
  else
    dserror("Unknown type of action for So_tet4");

  // check for patient specific data
  PATSPEC::GetILTDistance(Id(), params, discretization);
  PATSPEC::GetLocalRadius(Id(), params, discretization);
  PATSPEC::GetInnerRadius(Id(), params, discretization);

  // what should the element do
  switch (act)
  {
    //==================================================================================
    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    case calc_struct_linstiff:
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
      //      LINALG::Matrix<NUMDOF_SOTET4,NUMDOF_SOTET4>* matptr = NULL;
      //      if (elemat1.IsInitialized()) matptr = &elemat1;

      std::vector<double> mydispmat(lm.size(), 0.0);

      if (pstype_ == INPAR::STR::PreStress::id && time_ <= pstime_)  // inverse design analysis
        invdesign_->so_tet4_nlnstiffmass(params, this, lm, mydisp, myres, &elemat1, NULL, &elevec1,
            NULL, NULL, INPAR::STR::stress_none, INPAR::STR::strain_none);
      else
        nlnstiffmass(lm, mydisp, NULL, NULL, myres, mydispmat, &elemat1, NULL, &elevec1, NULL,
            &elevec3, NULL, NULL, NULL, params, INPAR::STR::stress_none, INPAR::STR::strain_none,
            INPAR::STR::strain_none);
    }
    break;

    //==================================================================================
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
      LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4> myemat(true);  // to zero

      nlnstiffmass(lm, mydisp, NULL, NULL, myres, mydispmat, &myemat, NULL, &elevec1, NULL, NULL,
          NULL, NULL, NULL, params, INPAR::STR::stress_none, INPAR::STR::strain_none,
          INPAR::STR::strain_none);
    }
    break;

    //==================================================================================
    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    case calc_struct_linstiffmass:
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

      if (pstype_ == INPAR::STR::PreStress::id && time_ <= pstime_)  // inverse design analysis
        invdesign_->so_tet4_nlnstiffmass(params, this, lm, mydisp, myres, &elemat1, &elemat2,
            &elevec1, NULL, NULL, INPAR::STR::stress_none, INPAR::STR::strain_none);
      else
        nlnstiffmass(lm, mydisp, &myvel, &myacc, myres, mydispmat, &elemat1, &elemat2, &elevec1,
            &elevec2, &elevec3, NULL, NULL, NULL, params, INPAR::STR::stress_none,
            INPAR::STR::strain_none, INPAR::STR::strain_none);

      if (act == calc_struct_nlnstifflmass) so_tet4_lumpmass(&elemat2);
    }
    break;

    //==================================================================================
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
        LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D> stress(true);  // set to zero
        LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D> strain(true);
        INPAR::STR::StressType iostress =
            DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
        INPAR::STR::StrainType iostrain =
            DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);

        std::vector<double> mydispmat(lm.size(), 0.0);

        if (pstype_ == INPAR::STR::PreStress::id && time_ <= pstime_)  // inverse design analysis
          invdesign_->so_tet4_nlnstiffmass(params, this, lm, mydisp, myres, NULL, NULL, NULL,
              &stress, &strain, iostress, iostrain);
        else
          nlnstiffmass(lm, mydisp, NULL, NULL, myres, mydispmat, NULL, NULL, NULL, NULL, NULL,
              &stress, &strain, NULL, params, iostress, iostrain, INPAR::STR::strain_none);

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

    //==================================================================================
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
      LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D> gpstress(((*gpstressmap)[gid])->A(), true);

      Teuchos::RCP<Epetra_MultiVector> poststress =
          params.get<Teuchos::RCP<Epetra_MultiVector>>("poststress", Teuchos::null);
      if (poststress == Teuchos::null) dserror("No element stress/strain vector available");

      if (stresstype == "ndxyz")
      {
        // extrapolate stresses/strains at Gauss points to nodes
        so_tet4_expol(gpstress, *poststress);
      }
      else if (stresstype == "cxyz")
      {
        const Epetra_BlockMap elemap = poststress->Map();
        int lid = elemap.LID(Id());
        if (lid != -1)
        {
          for (int i = 0; i < MAT::NUM_STRESS_3D; ++i)
          {
            double& s = (*((*poststress)(i)))[lid];
            s = 0.;
            for (int j = 0; j < NUMGPT_SOTET4; ++j)
            {
              s += gpstress(j, i);
            }
            s /= NUMGPT_SOTET4;
          }
        }
      }
      else
        dserror("unknown type of stress/strain output on element level");
    }
    break;

    //==================================================================================
    case prestress_update:
    {
      time_ = params.get<double>("total time");
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get displacement state");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      // build incremental def gradient for every gauss point
      LINALG::SerialDenseMatrix gpdefgrd(NUMGPT_SOTET4, 9);
      DefGradient(mydisp, gpdefgrd, *prestress_);

      // update deformation gradient and put back to storage
      LINALG::Matrix<3, 3> deltaF;
      LINALG::Matrix<3, 3> Fhist;
      LINALG::Matrix<3, 3> Fnew;
      for (int gp = 0; gp < NUMGPT_SOTET4; ++gp)
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
    //==================================================================================
    // this is a dummy output for strain energy
    case calc_struct_energy:
    {
      // check length of elevec1
      if (elevec1_epetra.Length() < 1) dserror("The given result vector is too short.");

      // initialization of internal energy
      double intenergy = 0.0;

      const static std::vector<LINALG::Matrix<NUMNOD_SOTET4, 1>> shapefcts =
          so_tet4_1gp_shapefcts();
      const static std::vector<LINALG::Matrix<NUMDIM_SOTET4 + 1, NUMNOD_SOTET4>> derivs =
          so_tet4_1gp_derivs();
      const static std::vector<double> gpweights = so_tet4_1gp_weights();

      // get displacements of this processor
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state displacement vector");

      // get displacements of this element
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      /* ============================================================================*/
      // element geometry
      /* structure of xrefe:
       **             [  X_1   Y_1   Z_1  ]
       **     xrefe = [  X_2   Y_2   Z_2  ]
       **             [   |     |     |   ]
       **             [  X_4   Y_4   Z_4  ]
       */
      /* structure of xcurr:
       **             [  x_1   y_1   z_1  ]
       **     xcurr = [  x_2   y_2   z_2  ]
       **             [   |     |     |   ]
       **             [  x_4   y_4   z_4  ]
       */
      // current  displacements of element
      LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xrefe;
      LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xdisp;
      DRT::Node** nodes = Nodes();
      for (int i = 0; i < NUMNOD_SOTET4; ++i)
      {
        const double* x = nodes[i]->X();
        xrefe(i, 0) = x[0];
        xrefe(i, 1) = x[1];
        xrefe(i, 2) = x[2];

        xdisp(i, 0) = mydisp[i * NODDOF_SOTET4 + 0];
        xdisp(i, 1) = mydisp[i * NODDOF_SOTET4 + 1];
        xdisp(i, 2) = mydisp[i * NODDOF_SOTET4 + 2];
      }


      // volume of a tetrahedra
      double detJ = V_;

      /* =========================================================================*/
      /* ============================================== Loop over Gauss Points ===*/
      /* =========================================================================*/
      for (int gp = 0; gp < NUMGPT_SOTET4; gp++)
      {
        LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> nxyz(nxyz_);  // copy!

        //                                      d xcurr
        // (material) deformation gradient F = --------- = xcurr^T * nxyz^T
        //                                      d xrefe

        /*structure of F
        **             [    dx       dy       dz    ]
        **             [  ------   ------   ------  ]
        **             [    dX       dX       dX    ]
        **             [                            ]
        **      F   =  [    dx       dy       dz    ]
        **             [  ------   ------   ------  ]
        **             [    dY       dY       dY    ]
        **             [                            ]
        **             [    dx       dy       dz    ]
        **             [  ------   ------   ------  ]
        **             [    dZ       dZ       dZ    ]
        */

        // size is 3x3
        LINALG::Matrix<3, 3> defgrd(true);
        // Gauss weights and Jacobian determinant
        double fac = detJ * gpweights[gp];

        if (pstype_ == INPAR::STR::PreStress::mulf)
        {
          // get derivatives wrt to last spatial configuration
          LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> N_xyz;
          prestress_->StoragetoMatrix(gp, N_xyz, prestress_->JHistory());

          // build multiplicative incremental defgrd
          // defgrd.Multiply('T','N',1.0,xdisp,N_xyz,0.0);
          if (kintype_ == INPAR::STR::kinem_nonlinearTotLag)
          {
            defgrd.MultiplyTN(xdisp, N_xyz);
          }
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
          // in kinematically linear analysis the deformation gradient is equal to identity
          if (kintype_ == INPAR::STR::kinem_nonlinearTotLag)
          {
            defgrd.MultiplyTN(xdisp, nxyz);
          }
          defgrd(0, 0) += 1.0;
          defgrd(1, 1) += 1.0;
          defgrd(2, 2) += 1.0;
        }


        if (pstype_ == INPAR::STR::PreStress::id && pstime_ < time_)
        {
          dserror("Calc Energy not implemented for prestress id");
        }

        // Right Cauchy-Green tensor = F^T * F
        // size is 3x3
        LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> cauchygreen;
        cauchygreen.MultiplyTN(defgrd, defgrd);

        // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
        // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
        LINALG::Matrix<6, 1> glstrain(false);
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
    // this is needed by bone topology optimization
    case calc_struct_output_E:
    {
      Teuchos::RCP<MAT::Material> mat = Material();
      // check length of elevec1
      if (elevec1_epetra.Length() < 1) dserror("The given result vector is too short.");
      MAT::Material* rawmat = mat.get();
      MAT::StVenantKirchhoff* stvk = dynamic_cast<MAT::StVenantKirchhoff*>(rawmat);
      if (!stvk) dserror("dynamic cast to stvenant failed");
      double E = stvk->Youngs();
      elevec1_epetra(0) = E;
    }
    break;

    //==================================================================================
    case inversedesign_update:
    {
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get displacement state");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      invdesign_->sot4_StoreMaterialConfiguration(this, mydisp);
      invdesign_->IsInit() = true;  // this is to make the restart work
    }
    break;

    //==================================================================================
    case inversedesign_switch:
    {
      time_ = params.get<double>("total time");
    }
    break;

    //==================================================================================
    case calc_struct_eleload:
      dserror("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
      break;

    //==================================================================================
    case calc_struct_fsiload:
      dserror("Case not yet implemented");
      break;

    //==================================================================================
    case struct_calc_store_istep:
    {
      int timestep = params.get<int>("timestep", -1);

      if (timestep == -1) dserror("Provide timestep number to be stored");

      // due to the multiplicativity and futility to redo prestress steps
      // other than the last one, no need to store/recover anything
      // ... but keep in mind
      if (pstype_ != INPAR::STR::PreStress::none)
      {
      }

      // Material
      SolidMaterial()->StoreHistory(timestep);
    }
    break;

    //==================================================================================
    case struct_calc_recover_istep:
    {
      int timestep = params.get<int>("timestep", -1);

      if (timestep == -1) dserror("Provide timestep number of the timestep to be recovered");

      // due to the multiplicativity and futility to redo prestress steps
      // other than the last one, no need to store/recover anything
      // ... but keep in mind
      if (pstype_ != INPAR::STR::PreStress::none)
      {
      }

      // Material
      SolidMaterial()->SetHistory(timestep);
    }
    break;

    //==================================================================================
    case calc_struct_update_istep:
    {
      // determine new fiber directions
      bool remodel;
      const Teuchos::ParameterList& patspec = DRT::Problem::Instance()->PatSpecParams();
      remodel = DRT::INPUT::IntegralValue<int>(patspec, "REMODEL");
      if (remodel)
      {
        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
        so_tet4_remodel(lm, mydisp, params, Material());
      }
      // Update of history for materials
      SolidMaterial()->Update();
    }
    break;

    //==================================================================================
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
      SolidMaterial()->ResetAll(NUMGPT_SOTET4);

      // Reset prestress
      if (pstype_ == INPAR::STR::PreStress::mulf)
      {
        time_ = 0.0;
        LINALG::Matrix<3, 3> Id(true);
        Id(0, 0) = Id(1, 1) = Id(2, 2) = 1.0;
        for (int gp = 0; gp < NUMGPT_SOTET4; ++gp)
        {
          prestress_->MatrixtoStorage(gp, Id, prestress_->FHistory());
          prestress_->MatrixtoStorage(gp, nxyz_, prestress_->JHistory());
        }
      }
      if (pstype_ == INPAR::STR::PreStress::id)
        dserror("Reset of Inverse Design not yet implemented");
    }
    break;

    //==================================================================================
    case calc_struct_errornorms:
    {
      // IMPORTANT NOTES (popp 10/2010):
      // - error norms are based on a small deformation assumption (linear elasticity)
      // - extension to finite deformations would be possible without difficulties,
      //   however analytical solutions are extremely rare in the nonlinear realm
      // - 4 Gauss point rule is used for integration of error norms
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

        // use 4 Gauss points for integration, not only 1
        int ngp = 4;

        // shape functions, derivatives and integration weights
        const static std::vector<LINALG::Matrix<NUMNOD_SOTET4, 1>> vals = so_tet4_4gp_shapefcts();
        const static std::vector<double> weights = so_tet4_4gp_weights();

        // get displacements and extract values of this element
        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        if (disp == Teuchos::null) dserror("Cannot get state displacement vector");
        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

        // nodal displacement vector
        LINALG::Matrix<NUMDOF_SOTET4, 1> nodaldisp;
        for (int i = 0; i < NUMDOF_SOTET4; ++i) nodaldisp(i, 0) = mydisp[i];

        // reference geometry (nodal positions)
        LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xrefe;
        DRT::Node** nodes = Nodes();
        for (int i = 0; i < NUMNOD_SOTET4; ++i)
        {
          xrefe(i, 0) = nodes[i]->X()[0];
          xrefe(i, 1) = nodes[i]->X()[1];
          xrefe(i, 2) = nodes[i]->X()[2];
        }

        // deformation gradient = identity tensor (geometrically linear case!)
        LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> defgrd(true);
        for (int i = 0; i < NUMDIM_SOTET4; ++i) defgrd(i, i) = 1.0;

        //----------------------------------------------------------------
        // loop over all Gauss points
        //----------------------------------------------------------------
        for (int gp = 0; gp < ngp; gp++)
        {
          // Gauss weights and Jacobian determinant
          double fac = V_ * weights[gp];

          // Gauss point in reference configuration
          LINALG::Matrix<NUMDIM_SOTET4, 1> xgp(true);
          for (int k = 0; k < NUMDIM_SOTET4; ++k)
            for (int n = 0; n < NUMNOD_SOTET4; ++n) xgp(k, 0) += (vals[gp])(n)*xrefe(n, k);

          //**************************************************************
          // get analytical solution
          LINALG::Matrix<NUMDIM_SOTET4, 1> uanalyt(true);
          LINALG::Matrix<MAT::NUM_STRESS_3D, 1> strainanalyt(true);
          LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> derivanalyt(true);

          CONTACT::AnalyticalSolutions3D(xgp, uanalyt, strainanalyt, derivanalyt);
          //**************************************************************

          //--------------------------------------------------------------
          // (1) L2 norm
          //--------------------------------------------------------------

          // compute displacements at GP
          LINALG::Matrix<NUMDIM_SOTET4, 1> ugp(true);
          for (int k = 0; k < NUMDIM_SOTET4; ++k)
            for (int n = 0; n < NUMNOD_SOTET4; ++n)
              ugp(k, 0) += (vals[gp])(n)*nodaldisp(NODDOF_SOTET4 * n + k, 0);

          // displacement error
          LINALG::Matrix<NUMDIM_SOTET4, 1> uerror(true);
          for (int k = 0; k < NUMDIM_SOTET4; ++k) uerror(k, 0) = uanalyt(k, 0) - ugp(k, 0);

          // compute GP contribution to L2 error norm
          l2norm += fac * uerror.Dot(uerror);

          //--------------------------------------------------------------
          // (2) H1 norm
          //--------------------------------------------------------------

          // compute derivatives N_XYZ at GP w.r.t. material coordinates
          // (we have the transposed stored in nxyz_ and then switch indices)
          LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> tr_N_XYZ(nxyz_);  // copy!
          LINALG::Matrix<NUMDIM_SOTET4, NUMNOD_SOTET4> N_XYZ(true);
          for (int k = 0; k < NUMNOD_SOTET4; ++k)
            for (int m = 0; m < NUMDIM_SOTET4; ++m) N_XYZ(m, k) = tr_N_XYZ(k, m);

          // compute partial derivatives at GP
          LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> derivgp(true);
          for (int l = 0; l < NUMDIM_SOTET4; ++l)
            for (int m = 0; m < NUMDIM_SOTET4; ++m)
              for (int k = 0; k < NUMNOD_SOTET4; ++k)
                derivgp(l, m) += N_XYZ(m, k) * nodaldisp(NODDOF_SOTET4 * k + l, 0);

          // derivative error
          LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> deriverror(true);
          for (int k = 0; k < NUMDIM_SOTET4; ++k)
            for (int m = 0; m < NUMDIM_SOTET4; ++m)
              deriverror(k, m) = derivanalyt(k, m) - derivgp(k, m);

          // compute GP contribution to H1 error norm
          h1norm += fac * deriverror.Dot(deriverror);
          h1norm += fac * uerror.Dot(uerror);

          //--------------------------------------------------------------
          // (3) Energy norm
          //--------------------------------------------------------------

          // compute linear B-operator
          LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOTET4> bop;
          for (int i = 0; i < NUMNOD_SOTET4; ++i)
          {
            bop(0, NODDOF_SOTET4 * i + 0) = N_XYZ(0, i);
            bop(0, NODDOF_SOTET4 * i + 1) = 0.0;
            bop(0, NODDOF_SOTET4 * i + 2) = 0.0;
            bop(1, NODDOF_SOTET4 * i + 0) = 0.0;
            bop(1, NODDOF_SOTET4 * i + 1) = N_XYZ(1, i);
            bop(1, NODDOF_SOTET4 * i + 2) = 0.0;
            bop(2, NODDOF_SOTET4 * i + 0) = 0.0;
            bop(2, NODDOF_SOTET4 * i + 1) = 0.0;
            bop(2, NODDOF_SOTET4 * i + 2) = N_XYZ(2, i);

            bop(3, NODDOF_SOTET4 * i + 0) = N_XYZ(1, i);
            bop(3, NODDOF_SOTET4 * i + 1) = N_XYZ(0, i);
            bop(3, NODDOF_SOTET4 * i + 2) = 0.0;
            bop(4, NODDOF_SOTET4 * i + 0) = 0.0;
            bop(4, NODDOF_SOTET4 * i + 1) = N_XYZ(2, i);
            bop(4, NODDOF_SOTET4 * i + 2) = N_XYZ(1, i);
            bop(5, NODDOF_SOTET4 * i + 0) = N_XYZ(2, i);
            bop(5, NODDOF_SOTET4 * i + 1) = 0.0;
            bop(5, NODDOF_SOTET4 * i + 2) = N_XYZ(0, i);
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
    }
    break;

    //==================================================================================
    case multi_calc_dens:
    {
      sotet4_homog(params);
    }
    break;

    //==================================================================================
    // read restart of microscale
    case multi_readrestart:
    {
      sotet4_read_restart_multi();
    }
    break;
    // evaluate stresses and strains at gauss points and store gpstresses in map <EleId, gpstresses
    // >
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
        LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D> stress;
        LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D> strain;
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
          if (pstype_ == INPAR::STR::PreStress::id && time_ <= pstime_)  // inverse design analysis
            invdesign_->so_tet4_nlnstiffmass(params, this, lm, mydisp, myres, NULL, NULL, NULL,
                &stress, &strain, iostress, iostrain);
          else
            nlnstiffmass(lm, mydisp, NULL, NULL, myres, mydispmat, NULL, NULL, NULL, NULL, NULL,
                &stress, &strain, NULL, params, iostress, iostrain, INPAR::STR::strain_none);
        }
        // add stresses to global map
        // get EleID Id()
        int gid = Id();
        Teuchos::RCP<Epetra_SerialDenseMatrix> gpstress =
            Teuchos::rcp(new Epetra_SerialDenseMatrix);
        gpstress->Shape(NUMGPT_SOTET4, MAT::NUM_STRESS_3D);

        // move stresses to serial dense matrix
        for (int i = 0; i < NUMGPT_SOTET4; i++)
        {
          for (int j = 0; j < MAT::NUM_STRESS_3D; j++)
          {
            (*gpstress)(i, j) = stress(i, j);
          }
        }

        // strains
        Teuchos::RCP<Epetra_SerialDenseMatrix> gpstrain =
            Teuchos::rcp(new Epetra_SerialDenseMatrix);
        gpstrain->Shape(NUMGPT_SOTET4, MAT::NUM_STRESS_3D);

        // move stresses to serial dense matrix
        for (int i = 0; i < NUMGPT_SOTET4; i++)
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
      dserror("Unknown type of action for so_tet4");
      break;
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)     maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet4::EvaluateNeumann(Teuchos::ParameterList& params,
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
  if (int(onoff->size()) < NUMDIM_SOTET4)
    dserror("Fewer functions or curves defined than the element has dofs.");

  for (int checkdof = NUMDIM_SOTET4; checkdof < int(onoff->size()); ++checkdof)
  {
    if ((*onoff)[checkdof] != 0)
      dserror("Number of Dimensions in Neumann_Evalutaion is 3. Further DoFs are not considered.");
  }

  // (SPATIAL) FUNCTION BUSINESS
  BOOST_STATIC_ASSERT((NUMGPT_SOTET4 == 1)) BACI_ATTRIBUTE_UNUSED;
  const std::vector<int>* funct = condition.Get<std::vector<int>>("funct");
  LINALG::Matrix<NUMDIM_SOTET4, 1> xrefegp(false);
  bool havefunct = false;
  if (funct)
    for (int dim = 0; dim < NUMDIM_SOTET4; dim++)
      if ((*funct)[dim] > 0)
      {
        havefunct = true;
        break;
      }


  /* =============================================================================*
   * CONST SHAPE FUNCTIONS and WEIGHTS for TET_4 with 1 GAUSS POINTS              *
   * =============================================================================*/
  const static std::vector<LINALG::Matrix<NUMNOD_SOTET4, 1>> shapefcts = so_tet4_1gp_shapefcts();
  const static std::vector<double> gpweights = so_tet4_1gp_weights();
  /* ============================================================================*/

  // update element geometry
  LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xrefe;
  DRT::Node** nodes = Nodes();
  for (int i = 0; i < NUMNOD_SOTET4; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }

  /* get the matrix of the coordinates of nodes needed to compute the volume,
  ** which is used here as detJ in the quadrature rule.
  ** ("Jacobian matrix") for the quadrature rule:
  **             [  1    1    1    1  ]
  **         J = [ X_1  X_2  X_3  X_4 ]
  **             [ Y_1  Y_2  Y_3  Y_4 ]
  **             [ Z_1  Z_2  Z_3  Z_4 ]
  */
  LINALG::Matrix<NUMCOORD_SOTET4, NUMCOORD_SOTET4> jac;
  for (int i = 0; i < 4; i++) jac(0, i) = 1.0;
  for (int row = 0; row < 3; row++)
    for (int col = 0; col < 4; col++) jac(row + 1, col) = xrefe(col, row);

  // compute determinant of Jacobian once outside Gauss point loop since it is constant
  // jac.Determinant() delivers six times the reference volume of the tet
  const double detJ = jac.Determinant() * (1.0 / 6.0);

  if (detJ == 0.0)
    dserror("ZERO JACOBIAN DETERMINANT");
  else if (detJ < 0.0)
    dserror("NEGATIVE JACOBIAN DETERMINANT");

  /* ================================================= Loop over Gauss Points */
  for (int gp = 0; gp < NUMGPT_SOTET4; gp++)
  {
    // material/reference co-ordinates of Gauss point
    if (havefunct)
    {
      for (int dim = 0; dim < NUMDIM_SOTET4; dim++)
      {
        xrefegp(dim) = 0.0;
        for (int nodid = 0; nodid < NUMNOD_SOTET4; ++nodid)
          xrefegp(dim) += shapefcts[gp](nodid) * xrefe(nodid, dim);
      }
    }

    // integration factor
    double fac = gpweights[gp] * detJ;
    // distribute/add over element load vector
    for (int dim = 0; dim < NUMDIM_SOTET4; dim++)
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
        for (int nodid = 0; nodid < NUMNOD_SOTET4; ++nodid)
        {
          elevec1[nodid * NUMDIM_SOTET4 + dim] += shapefcts[gp](nodid) * dim_fac;
        }
      }
    }


  } /* ==================================================== end of Loop over GP */

  return 0;
}  // DRT::ELEMENTS::So_tet4::EvaluateNeumann


/*----------------------------------------------------------------------*
 |  init the element jacobian mapping (protected)              gee 05/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::InitJacobianMapping()
{
  LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xrefe;
  DRT::Node** nodes = Nodes();
  for (int i = 0; i < NUMNOD_SOTET4; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }
  /* get the matrix of the coordinates of nodes needed to compute the volume,
  ** which is used here as detJ in the quadrature rule.
  ** ("Jacobian matrix") for the quadrarture rule:
  **             [  1    1    1    1  ]
  **         J = [ X_1  X_2  X_3  X_4 ]
  **             [ Y_1  Y_2  Y_3  Y_4 ]
  **             [ Z_1  Z_2  Z_3  Z_4 ]
  */
  LINALG::Matrix<NUMCOORD_SOTET4, NUMCOORD_SOTET4> jac;
  for (int i = 0; i < 4; i++) jac(0, i) = 1;
  for (int row = 0; row < 3; row++)
    for (int col = 0; col < 4; col++) jac(row + 1, col) = xrefe(col, row);
  // volume of the element
  V_ = jac.Determinant() / 6.0;
  if (V_ <= 0.0) dserror("Element volume %10.5e <= 0.0", V_);

  // nxyz_.resize(NUMGPT_SOTET4);
  const static std::vector<LINALG::Matrix<NUMDIM_SOTET4 + 1, NUMNOD_SOTET4>> derivs =
      so_tet4_1gp_derivs();
  LINALG::Matrix<NUMCOORD_SOTET4 - 1, NUMCOORD_SOTET4> tmp;
  for (int gp = 0; gp < NUMGPT_SOTET4; ++gp)
  {
    tmp.MultiplyTN(xrefe, derivs[gp]);
    for (int i = 0; i < 4; i++) jac(0, i) = 1;
    for (int row = 0; row < 3; row++)
      for (int col = 0; col < 4; col++) jac(row + 1, col) = tmp(row, col);
    // size is 4x3
    LINALG::Matrix<NUMCOORD_SOTET4, NUMDIM_SOTET4> I_aug(true);
    // size is 4x3
    LINALG::Matrix<NUMCOORD_SOTET4, NUMDIM_SOTET4> partials(true);
    I_aug(1, 0) = 1;
    I_aug(2, 1) = 1;
    I_aug(3, 2) = 1;

    // solve A.X=B
    LINALG::FixedSizeSerialDenseSolver<NUMCOORD_SOTET4, NUMCOORD_SOTET4, NUMDIM_SOTET4>
        solve_for_inverseJac;
    solve_for_inverseJac.SetMatrix(jac);               // set A=jac
    solve_for_inverseJac.SetVectors(partials, I_aug);  // set X=partials, B=I_aug
    solve_for_inverseJac.FactorWithEquilibration(true);
    int err2 = solve_for_inverseJac.Factor();
    int err = solve_for_inverseJac.Solve();  // partials = jac^-1.I_aug
    if ((err != 0) || (err2 != 0)) dserror("Inversion of Jacobian failed");

    // nxyz_[gp] = N_xsi_k*partials
    nxyz_.Multiply(derivs[gp], partials);
    /* structure of N_XYZ:
    **             [   dN_1     dN_1     dN_1   ]
    **             [  ------   ------   ------  ]
    **             [    dX       dY       dZ    ]
    **    N_XYZ =  [     |        |        |    ]
    **             [                            ]
    **             [   dN_4     dN_4     dN_4   ]
    **             [  -------  -------  ------- ]
    **             [    dX       dY       dZ    ]
    */

    if (pstype_ == INPAR::STR::PreStress::mulf && pstime_ >= time_)
      if (!(prestress_->IsInit())) prestress_->MatrixtoStorage(gp, nxyz_, prestress_->JHistory());

    if (pstype_ == INPAR::STR::PreStress::id && pstime_ < time_)
      if (!(invdesign_->IsInit()))
      {
        invdesign_->MatrixtoStorage(gp, nxyz_, invdesign_->JHistory());
        invdesign_->DetJHistory()[gp] = V_;
      }

  }  // for (int gp=0; gp<NUMGPT_SOTET4; ++gp)

  if (pstype_ == INPAR::STR::PreStress::mulf && pstime_ >= time_) prestress_->IsInit() = true;

  if (pstype_ == INPAR::STR::PreStress::id && pstime_ < time_) invdesign_->IsInit() = true;

  return;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (private)                            vlf 08/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::nlnstiffmass(std::vector<int>& lm,  // location matrix
    std::vector<double>& disp,                                   // current displacements
    std::vector<double>* vel,                                    // current velocities
    std::vector<double>* acc,                                    // current accelerations
    std::vector<double>& residual,                               // current residual displ
    std::vector<double>& dispmat,                                // current material displacements
    LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4>* stiffmatrix,   // element stiffness matrix
    LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4>* massmatrix,    // element mass matrix
    LINALG::Matrix<NUMDOF_SOTET4, 1>* force,                     // element internal force vector
    LINALG::Matrix<NUMDOF_SOTET4, 1>* forceinert,                // element inertial force vector
    LINALG::Matrix<NUMDOF_SOTET4, 1>* force_str,                 // element structural force vector
    LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D>* elestress,    // stresses at GP
    LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D>* elestrain,    // strains at GP
    LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D>* eleplstrain,  // plastic strains at GP
    Teuchos::ParameterList& params,          // algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,   // stress output option
    const INPAR::STR::StrainType iostrain,   // strain output option
    const INPAR::STR::StrainType ioplstrain  // plastic strain output option
)
{
  /* =============================================================================*
  ** CONST DERIVATIVES and WEIGHTS for TET_4  with 1 GAUSS POINTS*
  ** =============================================================================*/
  const static std::vector<LINALG::Matrix<NUMNOD_SOTET4, 1>> shapefcts = so_tet4_1gp_shapefcts();
  const static std::vector<LINALG::Matrix<NUMDIM_SOTET4 + 1, NUMNOD_SOTET4>> derivs =
      so_tet4_1gp_derivs();
  const static std::vector<double> gpweights = so_tet4_1gp_weights();
  /* ============================================================================*/
  // element geometry
  /* structure of xrefe:
   **             [  X_1   Y_1   Z_1  ]
   **     xrefe = [  X_2   Y_2   Z_2  ]
   **             [   |     |     |   ]
   **             [  X_4   Y_4   Z_4  ]
   */
  /* structure of xcurr:
   **             [  x_1   y_1   z_1  ]
   **     xcurr = [  x_2   y_2   z_2  ]
   **             [   |     |     |   ]
   **             [  x_4   y_4   z_4  ]
   */
  // current  displacements of element
  LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xrefe;
  LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xdisp;
  DRT::Node** nodes = Nodes();
  for (int i = 0; i < NUMNOD_SOTET4; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];

    xdisp(i, 0) = disp[i * NODDOF_SOTET4 + 0];
    xdisp(i, 1) = disp[i * NODDOF_SOTET4 + 1];
    xdisp(i, 2) = disp[i * NODDOF_SOTET4 + 2];
  }


  // volume of a tetrahedra
  double detJ = V_;


  // size is 3x3
  LINALG::Matrix<3, 3> defgrd(true);
  /* =========================================================================*/
  /* ============================================== Loop over Gauss Points ===*/
  /* =========================================================================*/
  for (int gp = 0; gp < NUMGPT_SOTET4; gp++)
  {
    LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> nxyz(nxyz_);  // copy!

    //                                      d xcurr
    // (material) deformation gradient F = --------- = xcurr^T * nxyz^T
    //                                      d xrefe

    /*structure of F
    **             [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dX       dX       dX    ]
    **             [                            ]
    **      F   =  [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dY       dY       dY    ]
    **             [                            ]
    **             [    dx       dy       dz    ]
    **             [  ------   ------   ------  ]
    **             [    dZ       dZ       dZ    ]
    */

    if (pstype_ == INPAR::STR::PreStress::mulf)
    {
      // get derivatives wrt to last spatial configuration
      LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> N_xyz;
      prestress_->StoragetoMatrix(gp, N_xyz, prestress_->JHistory());

      // build multiplicative incremental defgrd
      if (kintype_ == INPAR::STR::kinem_nonlinearTotLag)
      {
        // defgrd.Multiply('T','N',1.0,xdisp,N_xyz,0.0);
        defgrd.MultiplyTN(xdisp, N_xyz);
      }
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
      // in kinematically linear analysis the deformation gradient is equal to identity
      if (kintype_ == INPAR::STR::kinem_nonlinearTotLag)
      {
        defgrd.MultiplyTN(xdisp, nxyz);
      }
      defgrd(0, 0) += 1.0;
      defgrd(1, 1) += 1.0;
      defgrd(2, 2) += 1.0;
    }

    if (pstype_ == INPAR::STR::PreStress::id && pstime_ < time_)
    {
      // make the multiplicative update so that defgrd refers to
      // the reference configuration that resulted from the inverse
      // design analysis
      LINALG::Matrix<3, 3> Fhist;
      invdesign_->StoragetoMatrix(gp, Fhist, invdesign_->FHistory());
      LINALG::Matrix<3, 3> tmp3x3;
      tmp3x3.Multiply(defgrd, Fhist);
      defgrd = tmp3x3;

      // make detJ and nxyzmat refer to the ref. configuration that resulted from
      // the inverse design analysis
      detJ = invdesign_->DetJHistory()[gp];
      invdesign_->StoragetoMatrix(gp, nxyz, invdesign_->JHistory());
    }

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
    LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOTET4> bop;
    for (int i = 0; i < NUMNOD_SOTET4; i++)
    {
      bop(0, NODDOF_SOTET4 * i + 0) = defgrd(0, 0) * nxyz(i, 0);
      bop(0, NODDOF_SOTET4 * i + 1) = defgrd(1, 0) * nxyz(i, 0);
      bop(0, NODDOF_SOTET4 * i + 2) = defgrd(2, 0) * nxyz(i, 0);
      bop(1, NODDOF_SOTET4 * i + 0) = defgrd(0, 1) * nxyz(i, 1);
      bop(1, NODDOF_SOTET4 * i + 1) = defgrd(1, 1) * nxyz(i, 1);
      bop(1, NODDOF_SOTET4 * i + 2) = defgrd(2, 1) * nxyz(i, 1);
      bop(2, NODDOF_SOTET4 * i + 0) = defgrd(0, 2) * nxyz(i, 2);
      bop(2, NODDOF_SOTET4 * i + 1) = defgrd(1, 2) * nxyz(i, 2);
      bop(2, NODDOF_SOTET4 * i + 2) = defgrd(2, 2) * nxyz(i, 2);
      /* ~~~ */
      bop(3, NODDOF_SOTET4 * i + 0) = defgrd(0, 0) * nxyz(i, 1) + defgrd(0, 1) * nxyz(i, 0);
      bop(3, NODDOF_SOTET4 * i + 1) = defgrd(1, 0) * nxyz(i, 1) + defgrd(1, 1) * nxyz(i, 0);
      bop(3, NODDOF_SOTET4 * i + 2) = defgrd(2, 0) * nxyz(i, 1) + defgrd(2, 1) * nxyz(i, 0);
      bop(4, NODDOF_SOTET4 * i + 0) = defgrd(0, 1) * nxyz(i, 2) + defgrd(0, 2) * nxyz(i, 1);
      bop(4, NODDOF_SOTET4 * i + 1) = defgrd(1, 1) * nxyz(i, 2) + defgrd(1, 2) * nxyz(i, 1);
      bop(4, NODDOF_SOTET4 * i + 2) = defgrd(2, 1) * nxyz(i, 2) + defgrd(2, 2) * nxyz(i, 1);
      bop(5, NODDOF_SOTET4 * i + 0) = defgrd(0, 2) * nxyz(i, 0) + defgrd(0, 0) * nxyz(i, 2);
      bop(5, NODDOF_SOTET4 * i + 1) = defgrd(1, 2) * nxyz(i, 0) + defgrd(1, 0) * nxyz(i, 2);
      bop(5, NODDOF_SOTET4 * i + 2) = defgrd(2, 2) * nxyz(i, 0) + defgrd(2, 0) * nxyz(i, 2);
    }

    // Right Cauchy-Green tensor = F^T * F
    // size is 3x3
    LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(false);
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
      LINALG::Matrix<NUMDOF_SOTET4, 1> nodaldisp;
      for (int i = 0; i < NUMNOD_SOTET4; ++i)
      {
        nodaldisp(3 * i, 0) = xdisp(i, 0);
        nodaldisp(3 * i + 1, 0) = xdisp(i, 1);
        nodaldisp(3 * i + 2, 0) = xdisp(i, 2);
      }

      // build the linearised strain epsilon = B_L . d
      glstrain.Multiply(bop, nodaldisp);
    }

    // return gp strains (only in case of stress/strain output)
    switch (iostrain)
    {
      case INPAR::STR::strain_gl:
      {
        if (elestrain == NULL) dserror("no strain data available");
        for (int i = 0; i < 3; ++i) (*elestrain)(gp, i) = glstrain(i);
        for (int i = 3; i < 6; ++i) (*elestrain)(gp, i) = 0.5 * glstrain(i);
      }
      break;
      case INPAR::STR::strain_ea:
      {
        if (elestrain == NULL) dserror("no strain data available");

        // rewriting Green-Lagrange strains in matrix format
        LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> gl;
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
        // Epetra_SerialDenseMatrix invdefgrd(defgrd); // make a copy here otherwise defgrd is
        // destroyed! LINALG::NonsymInverse3x3(invdefgrd);
        LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> invdefgrd;
        invdefgrd.Invert(defgrd);

        LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> temp;
        LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> euler_almansi;
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
        LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> lnv(true);

        // checking if cauchy green is correctly determined to ensure eigen vectors in correct
        // direction i.e. a flipped eigenvector is also a valid solution C = \sum_{i=1}^3
        // (\lambda_i^2) \mathbf{n}_i \otimes \mathbf{n}_i
        LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> tempCG(true);

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
        LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> diffCG(true);

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
        dserror("requested strain option not available");
        break;
    }

    // call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> cmat(true);
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> stress(true);

    if (Material()->MaterialType() == INPAR::MAT::m_constraintmixture ||
        Material()->MaterialType() == INPAR::MAT::m_mixture_elasthyper)
    {
      // gp reference coordinates
      LINALG::Matrix<NUMNOD_SOTET4, 1> funct(true);
      funct = shapefcts[gp];
      LINALG::Matrix<1, NUMDIM_SOTET4> point(true);
      point.MultiplyTN(funct, xrefe);
      params.set("gprefecoord", point);
    }

    params.set<int>("gp", gp);


    UTILS::GetTemperatureForStructuralMaterial<tet4>(shapefcts[gp], params);

    SolidMaterial()->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, Id());

    // return gp stresses
    switch (iostress)
    {
      case INPAR::STR::stress_2pk:
      {
        if (elestress == NULL) dserror("no stress data available");
        for (int i = 0; i < MAT::NUM_STRESS_3D; ++i) (*elestress)(gp, i) = stress(i);
      }
      break;
      case INPAR::STR::stress_cauchy:
      {
        if (elestress == NULL) dserror("no stress data available");
        double detF = defgrd.Determinant();

        LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> pkstress;
        pkstress(0, 0) = stress(0);
        pkstress(0, 1) = stress(3);
        pkstress(0, 2) = stress(5);
        pkstress(1, 0) = pkstress(0, 1);
        pkstress(1, 1) = stress(1);
        pkstress(1, 2) = stress(4);
        pkstress(2, 0) = pkstress(0, 2);
        pkstress(2, 1) = pkstress(1, 2);
        pkstress(2, 2) = stress(2);

        LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> temp;
        LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> cauchystress;
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

    double detJ_w = detJ * (gpweights)[gp];

    // update of internal force vector
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
      // size is 6x12
      LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOTET4> cb;
      cb.Multiply(cmat, bop);  // temporary C . B
      // size is 12x12
      stiffmatrix->MultiplyTN(detJ_w, bop, cb, 1.0);

      if (kintype_ == INPAR::STR::kinem_nonlinearTotLag)
      {
        // integrate `geometric' stiffness matrix and add to keu
        // auxiliary integrated stress
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1> sfac(stress);
        // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
        sfac.Scale(detJ_w);
        // intermediate Sm.B_L
        double SmB_L[NUMDIM_SOTET4];
        // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)
        // with B_L = Ni,Xj see NiliFEM-Skript
        for (int inod = 0; inod < NUMNOD_SOTET4; ++inod)
        {
          SmB_L[0] = sfac(0) * nxyz(inod, 0) + sfac(3) * nxyz(inod, 1) + sfac(5) * nxyz(inod, 2);
          SmB_L[1] = sfac(3) * nxyz(inod, 0) + sfac(1) * nxyz(inod, 1) + sfac(4) * nxyz(inod, 2);
          SmB_L[2] = sfac(5) * nxyz(inod, 0) + sfac(4) * nxyz(inod, 1) + sfac(2) * nxyz(inod, 2);
          for (int jnod = 0; jnod < NUMNOD_SOTET4; ++jnod)
          {
            double bopstrbop = 0.0;  // intermediate value
            for (int idim = 0; idim < NUMDIM_SOTET4; ++idim)
              bopstrbop += nxyz(jnod, idim) * SmB_L[idim];
            (*stiffmatrix)(NUMDIM_SOTET4 * inod + 0, NUMDIM_SOTET4 * jnod + 0) += bopstrbop;
            (*stiffmatrix)(NUMDIM_SOTET4 * inod + 1, NUMDIM_SOTET4 * jnod + 1) += bopstrbop;
            (*stiffmatrix)(NUMDIM_SOTET4 * inod + 2, NUMDIM_SOTET4 * jnod + 2) += bopstrbop;
          }
        }
      }
    }
    /* =========================================================================*/
  } /* ==================================================== end of Loop over GP */
  /* =========================================================================*/


  // static integrator created in any case to safe "if-case"
  const static std::vector<LINALG::Matrix<NUMNOD_SOTET4, 1>> shapefcts4gp = so_tet4_4gp_shapefcts();
  const static std::vector<double> gpweights4gp = so_tet4_4gp_weights();
  // evaluate mass matrix
  if (massmatrix != NULL)
  {
    double density = Material()->Density(0);  // density at the only Gauss point the material has!
    // consistent mass matrix evaluated using a 4-point rule
    for (int gp = 0; gp < 4; gp++)
    {
      double factor = density * detJ * gpweights4gp[gp];
      double ifactor, massfactor;
      for (int inod = 0; inod < NUMNOD_SOTET4; ++inod)
      {
        ifactor = (shapefcts4gp[gp])(inod)*factor;
        for (int jnod = 0; jnod < NUMNOD_SOTET4; ++jnod)
        {
          massfactor = (shapefcts4gp[gp])(jnod)*ifactor;
          (*massmatrix)(NUMDIM_SOTET4 * inod + 0, NUMDIM_SOTET4 * jnod + 0) += massfactor;
          (*massmatrix)(NUMDIM_SOTET4 * inod + 1, NUMDIM_SOTET4 * jnod + 1) += massfactor;
          (*massmatrix)(NUMDIM_SOTET4 * inod + 2, NUMDIM_SOTET4 * jnod + 2) += massfactor;
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

        LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> nxyz(nxyz_);  // copy!

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
        LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOTET4> bop;
        for (int i = 0; i < NUMNOD_SOTET4; i++)
        {
          bop(0, NODDOF_SOTET4 * i + 0) = defgrd(0, 0) * nxyz(i, 0);
          bop(0, NODDOF_SOTET4 * i + 1) = defgrd(1, 0) * nxyz(i, 0);
          bop(0, NODDOF_SOTET4 * i + 2) = defgrd(2, 0) * nxyz(i, 0);
          bop(1, NODDOF_SOTET4 * i + 0) = defgrd(0, 1) * nxyz(i, 1);
          bop(1, NODDOF_SOTET4 * i + 1) = defgrd(1, 1) * nxyz(i, 1);
          bop(1, NODDOF_SOTET4 * i + 2) = defgrd(2, 1) * nxyz(i, 1);
          bop(2, NODDOF_SOTET4 * i + 0) = defgrd(0, 2) * nxyz(i, 2);
          bop(2, NODDOF_SOTET4 * i + 1) = defgrd(1, 2) * nxyz(i, 2);
          bop(2, NODDOF_SOTET4 * i + 2) = defgrd(2, 2) * nxyz(i, 2);
          /* ~~~ */
          bop(3, NODDOF_SOTET4 * i + 0) = defgrd(0, 0) * nxyz(i, 1) + defgrd(0, 1) * nxyz(i, 0);
          bop(3, NODDOF_SOTET4 * i + 1) = defgrd(1, 0) * nxyz(i, 1) + defgrd(1, 1) * nxyz(i, 0);
          bop(3, NODDOF_SOTET4 * i + 2) = defgrd(2, 0) * nxyz(i, 1) + defgrd(2, 1) * nxyz(i, 0);
          bop(4, NODDOF_SOTET4 * i + 0) = defgrd(0, 1) * nxyz(i, 2) + defgrd(0, 2) * nxyz(i, 1);
          bop(4, NODDOF_SOTET4 * i + 1) = defgrd(1, 1) * nxyz(i, 2) + defgrd(1, 2) * nxyz(i, 1);
          bop(4, NODDOF_SOTET4 * i + 2) = defgrd(2, 1) * nxyz(i, 2) + defgrd(2, 2) * nxyz(i, 1);
          bop(5, NODDOF_SOTET4 * i + 0) = defgrd(0, 2) * nxyz(i, 0) + defgrd(0, 0) * nxyz(i, 2);
          bop(5, NODDOF_SOTET4 * i + 1) = defgrd(1, 2) * nxyz(i, 0) + defgrd(1, 0) * nxyz(i, 2);
          bop(5, NODDOF_SOTET4 * i + 2) = defgrd(2, 2) * nxyz(i, 0) + defgrd(2, 0) * nxyz(i, 2);
        }

        // Right Cauchy-Green tensor = F^T * F
        // size is 3x3
        LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> cauchygreen;
        cauchygreen.MultiplyTN(defgrd, defgrd);

        // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
        // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(false);
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
          LINALG::Matrix<NUMDOF_SOTET4, 1> nodaldisp;
          for (int i = 0; i < NUMNOD_SOTET4; ++i)
          {
            nodaldisp(3 * i, 0) = xdisp(i, 0);
            nodaldisp(3 * i + 1, 0) = xdisp(i, 1);
            nodaldisp(3 * i + 2, 0) = xdisp(i, 2);
          }

          // build the linearised strain epsilon = B_L . d
          glstrain.Multiply(bop, nodaldisp);
        }

        // evaluate derivative of mass w.r.t. to right cauchy green tensor
        SolidMaterial()->EvaluateNonLinMass(
            &defgrd, &glstrain, params, &linmass_disp, &linmass_vel, Id());

        // multiply by 2.0 to get derivative w.r.t green lagrange strains and multiply by time
        // integration factor
        linmass_disp.Scale(2.0 * timintfac_dis);
        linmass_vel.Scale(2.0 * timintfac_vel);
        linmass.Update(1.0, linmass_disp, 1.0, linmass_vel, 0.0);

        // evaluate accelerations at time n+1 at gauss point
        LINALG::Matrix<NUMDIM_SOTET4, 1> myacc(true);
        for (int idim = 0; idim < NUMDIM_SOTET4; ++idim)
          for (int inod = 0; inod < NUMNOD_SOTET4; ++inod)
            myacc(idim) += shapefcts4gp[gp](inod) * (*acc)[idim + (inod * NUMDIM_SOTET4)];

        if (stiffmatrix != NULL)
        {
          // integrate linearisation of mass matrix
          //(B^T . d\rho/d disp . a) * detJ * w(gp)
          LINALG::Matrix<1, NUMDOF_SOTET4> cb;
          cb.MultiplyTN(linmass_disp, bop);
          for (int inod = 0; inod < NUMNOD_SOTET4; ++inod)
          {
            double factor = detJ * gpweights4gp[gp] * shapefcts4gp[gp](inod);
            for (int idim = 0; idim < NUMDIM_SOTET4; ++idim)
            {
              double massfactor = factor * myacc(idim);
              for (int jnod = 0; jnod < NUMNOD_SOTET4; ++jnod)
                for (int jdim = 0; jdim < NUMDIM_SOTET4; ++jdim)
                  (*massmatrix)(inod * NUMDIM_SOTET4 + idim, jnod * NUMDIM_SOTET4 + jdim) +=
                      massfactor * cb(jnod * NUMDIM_SOTET4 + jdim);
            }
          }
        }

        // internal force vector
        if (forceinert != NULL)
        {
          // integrate nonlinear inertia force term
          for (int inod = 0; inod < NUMNOD_SOTET4; ++inod)
          {
            double forcefactor = shapefcts4gp[gp](inod) * detJ * gpweights4gp[gp];
            for (int idim = 0; idim < NUMDIM_SOTET4; ++idim)
              (*forceinert)(inod * NUMDIM_SOTET4 + idim) += forcefactor * density * myacc(idim);
          }
        }
      }

    }  // end loop over mass matrix Gauss points

  }  // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++


  return;
}  // DRT::ELEMENTS::So_tet4::nlnstiffmass


/*----------------------------------------------------------------------*
 |  lump mass matrix (private)                               bborn 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::so_tet4_lumpmass(LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4>* emass)
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
 |  init the element (public)                                  gee 05/08|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_tet4Type::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So_tet4* actele = dynamic_cast<DRT::ELEMENTS::So_tet4*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_tet4* failed");
    actele->InitJacobianMapping();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet4 Shape fcts at 1 Gauss Point                           |
 *----------------------------------------------------------------------*/
const std::vector<LINALG::Matrix<NUMNOD_SOTET4, 1>> DRT::ELEMENTS::So_tet4::so_tet4_1gp_shapefcts()
{
  std::vector<LINALG::Matrix<NUMNOD_SOTET4, 1>> shapefcts(NUMGPT_SOTET4);

  // There is only one gausspoint, so the loop (and the vector) is not really needed.
  for (int gp = 0; gp < NUMGPT_SOTET4; gp++)
  {
    (shapefcts[gp])(0) = 0.25;
    (shapefcts[gp])(1) = 0.25;
    (shapefcts[gp])(2) = 0.25;
    (shapefcts[gp])(3) = 0.25;
  }

  return shapefcts;
}


/*----------------------------------------------------------------------*
 |  Evaluate Tet4 Shape fct derivs at 1 Gauss Point                     |
 *----------------------------------------------------------------------*/
const std::vector<LINALG::Matrix<NUMDIM_SOTET4 + 1, NUMNOD_SOTET4>>
DRT::ELEMENTS::So_tet4::so_tet4_1gp_derivs()
{
  std::vector<LINALG::Matrix<NUMDIM_SOTET4 + 1, NUMNOD_SOTET4>> derivs(NUMGPT_SOTET4);
  // There is only one gausspoint, so the loop (and the vector) is not really needed.
  for (int gp = 0; gp < NUMGPT_SOTET4; gp++)
  {
    (derivs[gp])(0, 0) = 1.0;
    (derivs[gp])(1, 0) = 0.0;
    (derivs[gp])(2, 0) = 0.0;
    (derivs[gp])(3, 0) = 0.0;

    (derivs[gp])(0, 1) = 0.0;
    (derivs[gp])(1, 1) = 1.0;
    (derivs[gp])(2, 1) = 0.0;
    (derivs[gp])(3, 1) = 0.0;

    (derivs[gp])(0, 2) = 0.0;
    (derivs[gp])(1, 2) = 0.0;
    (derivs[gp])(2, 2) = 1.0;
    (derivs[gp])(3, 2) = 0.0;

    (derivs[gp])(0, 3) = 0.0;
    (derivs[gp])(1, 3) = 0.0;
    (derivs[gp])(2, 3) = 0.0;
    (derivs[gp])(3, 3) = 1.0;
  }
  return derivs;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet4 Weights at 1 Gauss Point                              |
 *----------------------------------------------------------------------*/
const std::vector<double> DRT::ELEMENTS::So_tet4::so_tet4_1gp_weights()
{
  std::vector<double> weights(NUMGPT_SOTET4);
  // There is only one gausspoint, so the loop (and the vector) is not really needed.
  for (int i = 0; i < NUMGPT_SOTET4; ++i) weights[i] = 1.0;
  return weights;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet4 Shape fcts at 4 Gauss Points                          |
 *----------------------------------------------------------------------*/
const std::vector<LINALG::Matrix<NUMNOD_SOTET4, 1>> DRT::ELEMENTS::So_tet4::so_tet4_4gp_shapefcts()
{
  std::vector<LINALG::Matrix<NUMNOD_SOTET4, 1>> shapefcts(4);

  const double gploc_alpha =
      (5.0 + 3.0 * sqrt(5.0)) / 20.0;  // gp sampling point value for quadr. fct
  const double gploc_beta = (5.0 - sqrt(5.0)) / 20.0;

  const double xsi1[4] = {gploc_alpha, gploc_beta, gploc_beta, gploc_beta};
  const double xsi2[4] = {gploc_beta, gploc_alpha, gploc_beta, gploc_beta};
  const double xsi3[4] = {gploc_beta, gploc_beta, gploc_alpha, gploc_beta};
  const double xsi4[4] = {gploc_beta, gploc_beta, gploc_beta, gploc_alpha};

  for (int gp = 0; gp < 4; gp++)
  {
    (shapefcts[gp])(0) = xsi1[gp];
    (shapefcts[gp])(1) = xsi2[gp];
    (shapefcts[gp])(2) = xsi3[gp];
    (shapefcts[gp])(3) = xsi4[gp];
  }

  return shapefcts;
}


/*----------------------------------------------------------------------*
 |  Evaluate Tet4 Shape fct derivs at 4 Gauss Points                    |
 *----------------------------------------------------------------------*/
const std::vector<LINALG::Matrix<NUMDIM_SOTET4 + 1, NUMNOD_SOTET4>>
DRT::ELEMENTS::So_tet4::so_tet4_4gp_derivs()
{
  std::vector<LINALG::Matrix<NUMDIM_SOTET4 + 1, NUMNOD_SOTET4>> derivs(4);

  for (int gp = 0; gp < 4; gp++)
  {
    (derivs[gp])(0, 0) = 1.0;
    (derivs[gp])(1, 0) = 0.0;
    (derivs[gp])(2, 0) = 0.0;
    (derivs[gp])(3, 0) = 0.0;

    (derivs[gp])(0, 1) = 0.0;
    (derivs[gp])(1, 1) = 1.0;
    (derivs[gp])(2, 1) = 0.0;
    (derivs[gp])(3, 1) = 0.0;

    (derivs[gp])(0, 2) = 0.0;
    (derivs[gp])(1, 2) = 0.0;
    (derivs[gp])(2, 2) = 1.0;
    (derivs[gp])(3, 2) = 0.0;

    (derivs[gp])(0, 3) = 0.0;
    (derivs[gp])(1, 3) = 0.0;
    (derivs[gp])(2, 3) = 0.0;
    (derivs[gp])(3, 3) = 1.0;
  }
  return derivs;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet4 Weights at 4 Gauss Points                             |
 *----------------------------------------------------------------------*/
const std::vector<double> DRT::ELEMENTS::So_tet4::so_tet4_4gp_weights()
{
  std::vector<double> weights(4);
  for (int i = 0; i < 4; ++i)
  {
    weights[i] = 0.25;
  }
  return weights;
}


/*----------------------------------------------------------------------*
 |  compute def gradient at every gaussian point (protected)   gee 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::DefGradient(const std::vector<double>& disp,
    Epetra_SerialDenseMatrix& gpdefgrd, DRT::ELEMENTS::PreStress& prestress)
{
  // update element geometry
  LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xdisp;
  for (int i = 0; i < NUMNOD_SOTET4; ++i)
  {
    xdisp(i, 0) = disp[i * NODDOF_SOTET4 + 0];
    xdisp(i, 1) = disp[i * NODDOF_SOTET4 + 1];
    xdisp(i, 2) = disp[i * NODDOF_SOTET4 + 2];
  }

  for (int gp = 0; gp < NUMGPT_SOTET4; ++gp)
  {
    // get derivatives wrt to last spatial configuration
    LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> N_xyz;
    prestress_->StoragetoMatrix(gp, N_xyz, prestress_->JHistory());

    // build multiplicative incremental defgrd
    LINALG::Matrix<3, 3> defgrd(true);
    if (kintype_ == INPAR::STR::kinem_nonlinearTotLag)
    {
      defgrd.MultiplyTN(xdisp, N_xyz);
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
void DRT::ELEMENTS::So_tet4::UpdateJacobianMapping(
    const std::vector<double>& disp, DRT::ELEMENTS::PreStress& prestress)
{
  // get incremental disp
  LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xdisp;
  for (int i = 0; i < NUMNOD_SOTET4; ++i)
  {
    xdisp(i, 0) = disp[i * NODDOF_SOTET4 + 0];
    xdisp(i, 1) = disp[i * NODDOF_SOTET4 + 1];
    xdisp(i, 2) = disp[i * NODDOF_SOTET4 + 2];
  }

  LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> nxyzhist;
  LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> nxyznew;
  LINALG::Matrix<3, 3> defgrd(true);

  for (int gp = 0; gp < NUMGPT_SOTET4; ++gp)
  {
    // get the nxyz old state
    prestress.StoragetoMatrix(gp, nxyzhist, prestress.JHistory());
    // build multiplicative incremental defgrd
    if (kintype_ == INPAR::STR::kinem_nonlinearTotLag)
    {
      defgrd.MultiplyTN(xdisp, nxyzhist);
    }
    defgrd(0, 0) += 1.0;
    defgrd(1, 1) += 1.0;
    defgrd(2, 2) += 1.0;
    // make inverse of this defgrd
    defgrd.Invert();

    // push-forward of nxyz
    nxyznew.Multiply(nxyzhist, defgrd);
    // store new reference configuration
    prestress.MatrixtoStorage(gp, nxyznew, prestress.JHistory());

  }  // for (int gp=0; gp<NUMGPT_SOTET4; ++gp)

  return;
}

/*----------------------------------------------------------------------*
  |  remodeling of fiber directions (protected)               tinkl 01/10|
  *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::so_tet4_remodel(std::vector<int>& lm,  // location matrix
    std::vector<double>& disp,                                      // current displacements
    Teuchos::ParameterList& params,   // algorithmic parameters e.g. time
    Teuchos::RCP<MAT::Material> mat)  // material
{
  if ((Material()->MaterialType() == INPAR::MAT::m_constraintmixture) ||
      (Material()->MaterialType() == INPAR::MAT::m_elasthyper))
  {
    // in a first step ommit everything with prestress

    // current  displacements of element
    LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xdisp;
    for (int i = 0; i < NUMNOD_SOTET4; ++i)
    {
      xdisp(i, 0) = disp[i * NODDOF_SOTET4 + 0];
      xdisp(i, 1) = disp[i * NODDOF_SOTET4 + 1];
      xdisp(i, 2) = disp[i * NODDOF_SOTET4 + 2];
    }

    /* =========================================================================*/
    /* ============================================== Loop over Gauss Points ===*/
    /* =========================================================================*/
    // interpolated values of stress and defgrd for remodeling
    LINALG::Matrix<3, 3> avg_stress(true);
    LINALG::Matrix<3, 3> avg_defgrd(true);

    for (int gp = 0; gp < NUMGPT_SOTET4; gp++)
    {
      const LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4>& nxyz = nxyz_;

      //                                      d xcurr
      // (material) deformation gradient F = --------- = xcurr^T * nxyz^T
      //                                      d xrefe

      /*structure of F
      **             [    dx       dy       dz    ]
      **             [  ------   ------   ------  ]
      **             [    dX       dX       dX    ]
      **             [                            ]
      **      F   =  [    dx       dy       dz    ]
      **             [  ------   ------   ------  ]
      **             [    dY       dY       dY    ]
      **             [                            ]
      **             [    dx       dy       dz    ]
      **             [  ------   ------   ------  ]
      **             [    dZ       dZ       dZ    ]
      */

      // size is 3x3
      LINALG::Matrix<3, 3> defgrd(false);

      if (pstype_ == INPAR::STR::PreStress::mulf)
      {
        // get derivatives wrt to last spatial configuration
        LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> N_xyz;
        prestress_->StoragetoMatrix(gp, N_xyz, prestress_->JHistory());

        // build multiplicative incremental defgrd
        // defgrd.Multiply('T','N',1.0,xdisp,N_xyz,0.0);
        defgrd.MultiplyTN(xdisp, N_xyz);
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
        defgrd.MultiplyTN(xdisp, nxyz);
        defgrd(0, 0) += 1;
        defgrd(1, 1) += 1;
        defgrd(2, 2) += 1;
      }

      // Right Cauchy-Green tensor = F^T * F
      // size is 3x3
      LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> cauchygreen;
      cauchygreen.MultiplyTN(defgrd, defgrd);

      // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
      // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
      LINALG::Matrix<6, 1> glstrain(false);
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
        avg_stress.Update(1.0 / NUMGPT_SOTET4, cauchystress, 1.0);
        avg_defgrd.Update(1.0 / NUMGPT_SOTET4, defgrd, 1.0);
      }
      else
        dserror("material not implemented for remodeling");

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
    }  // end loop over gauss points
  }
}

/*----------------------------------------------------------------------*
 |  evaluate cauchy stress tensor                            seitz 11/16|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::GetCauchyAtXi(const LINALG::Matrix<3, 1>& xi,
    const std::vector<double>& disp, const LINALG::Matrix<3, 1>& n, const LINALG::Matrix<3, 1>& t,
    double& sigma_nt, Epetra_SerialDenseMatrix* DsntDd, Epetra_SerialDenseMatrix* D2sntDd2,
    Epetra_SerialDenseMatrix* D2sntDdDn, Epetra_SerialDenseMatrix* D2sntDdDt,
    Epetra_SerialDenseMatrix* D2sntDdDxi, LINALG::Matrix<3, 1>* DsntDn,
    LINALG::Matrix<3, 1>* DsntDt, LINALG::Matrix<3, 1>* DsntDxi, const std::vector<double>* temp,
    Epetra_SerialDenseMatrix* DsntDT, Epetra_SerialDenseMatrix* D2sntDdDT)
{
  if (temp || DsntDT || D2sntDdDT)
    dserror("Thermo-elastic Nitsche contact not yet implemented in so tet4");

  sigma_nt = 0.0;

  static LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xrefe(true);  // reference coord. of element
  static LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xcurr(true);  // current  coord. of element
  xrefe.Clear();
  xcurr.Clear();
  DRT::Node** nodes = Nodes();

  for (int i = 0; i < NUMNOD_SOTET4; ++i)
  {
    const double* x = nodes[i]->X();
    for (int d = 0; d < NUMDIM_SOTET4; ++d)
    {
      xrefe(i, d) = x[d];
      xcurr(i, d) = xrefe(i, d) + disp[i * NODDOF_SOTET4 + d];
    }
  }

  static LINALG::Matrix<NUMDIM_SOTET4, NUMNOD_SOTET4> deriv(true);
  deriv.Clear();
  DRT::UTILS::shape_function_deriv1<DRT::Element::tet4>(xi, deriv);

  static LINALG::Matrix<NUMDIM_SOTET4, NUMNOD_SOTET4> N_XYZ(true);
  static LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> invJ(true);
  invJ.Multiply(1.0, deriv, xrefe, 0.0);
  invJ.Invert();
  N_XYZ.Multiply(1.0, invJ, deriv, 0.0);
  static LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> defgrd(true);
  defgrd.MultiplyTT(1.0, xcurr, N_XYZ, 0.0);

  // linearization of deformation gradient F w.r.t. displacements
  static LINALG::Matrix<9, NUMDOF_SOTET4> DFDd(true);
  DFDd.Clear();
  if (DsntDd || D2sntDdDn || D2sntDdDt || D2sntDd2 || D2sntDdDxi)
  {
    for (int i = 0; i < NUMNOD_SOTET4; ++i)
    {
      DFDd(0, NODDOF_SOTET4 * i + 0) = N_XYZ(0, i);
      DFDd(1, NODDOF_SOTET4 * i + 1) = N_XYZ(1, i);
      DFDd(2, NODDOF_SOTET4 * i + 2) = N_XYZ(2, i);
      DFDd(3, NODDOF_SOTET4 * i + 0) = N_XYZ(1, i);
      DFDd(4, NODDOF_SOTET4 * i + 1) = N_XYZ(2, i);
      DFDd(5, NODDOF_SOTET4 * i + 0) = N_XYZ(2, i);
      DFDd(6, NODDOF_SOTET4 * i + 1) = N_XYZ(0, i);
      DFDd(7, NODDOF_SOTET4 * i + 2) = N_XYZ(1, i);
      DFDd(8, NODDOF_SOTET4 * i + 2) = N_XYZ(0, i);
    }
  }

  static LINALG::Matrix<9, 1> DsntDF(true);
  static LINALG::Matrix<9, 9> D2sntDF2(true);
  static LINALG::Matrix<9, NUMDIM_SOTET4> D2sntDFDn(true);
  static LINALG::Matrix<9, NUMDIM_SOTET4> D2sntDFDt(true);

  SolidMaterial()->EvaluateCauchy(
      defgrd, n, t, sigma_nt, DsntDn, DsntDt, &DsntDF, &D2sntDF2, &D2sntDFDn, &D2sntDFDt, 0);

  if (DsntDd)
  {
    DsntDd->Reshape(NUMDOF_SOTET4, 1);
    LINALG::Matrix<NUMDOF_SOTET4, 1> DsntDd_m(DsntDd->A(), true);
    DsntDd_m.MultiplyTN(1.0, DFDd, DsntDF, 0.0);
  }

  if (D2sntDdDn)
  {
    D2sntDdDn->Reshape(NUMDOF_SOTET4, NUMDIM_SOTET4);
    LINALG::Matrix<NUMDOF_SOTET4, NUMDIM_SOTET4> D2sntDdDn_m(D2sntDdDn->A(), true);
    D2sntDdDn_m.MultiplyTN(1.0, DFDd, D2sntDFDn, 0.0);
  }

  if (D2sntDdDt)
  {
    D2sntDdDt->Reshape(NUMDOF_SOTET4, NUMDIM_SOTET4);
    LINALG::Matrix<NUMDOF_SOTET4, NUMDIM_SOTET4> D2sntDdDt_m(D2sntDdDt->A(), true);
    D2sntDdDt_m.MultiplyTN(1.0, DFDd, D2sntDFDt, 0.0);
  }

  if (D2sntDd2)
  {
    D2sntDd2->Reshape(NUMDOF_SOTET4, NUMDOF_SOTET4);
    LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4> D2sntDd2_m(D2sntDd2->A(), true);
    static LINALG::Matrix<9, NUMDOF_SOTET4> D2sntDF2DFDd(true);
    D2sntDF2DFDd.Multiply(1.0, D2sntDF2, DFDd, 0.0);
    D2sntDd2_m.MultiplyTN(1.0, DFDd, D2sntDF2DFDd, 0.0);
  }

  if (D2sntDdDxi)
  {
    D2sntDdDxi->Reshape(NUMDOF_SOTET4, NUMDIM_SOTET4);
    LINALG::Matrix<NUMDOF_SOTET4, NUMDIM_SOTET4> D2sntDdDxi_m(D2sntDdDxi->A(), true);

    static LINALG::Matrix<DRT::UTILS::DisTypeToNumDeriv2<DRT::Element::tet4>::numderiv2,
        NUMNOD_SOTET4>
        deriv2(true);
    deriv2.Clear();
    DRT::UTILS::shape_function_deriv2<DRT::Element::tet4>(xi, deriv2);

    static LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xXF(true);
    static LINALG::Matrix<NUMDIM_SOTET4,
        DRT::UTILS::DisTypeToNumDeriv2<DRT::Element::tet4>::numderiv2>
        xXFsec(true);
    xXF.Update(1.0, xcurr, 0.0);
    xXF.MultiplyNT(-1.0, xrefe, defgrd, 1.0);
    xXFsec.MultiplyTT(1.0, xXF, deriv2, 0.0);

    static LINALG::Matrix<9, NUMDIM_SOTET4> DFDxi(true);
    DFDxi.Clear();
    for (int a = 0; a < NUMDIM_SOTET4; ++a)
      for (int b = 0; b < NUMDIM_SOTET4; ++b)
      {
        DFDxi(VoigtMapping::NonSymToVoigt9(a, b), 0) +=
            xXFsec(a, 0) * invJ(b, 0) + xXFsec(a, 3) * invJ(b, 1) + xXFsec(a, 4) * invJ(b, 2);
        DFDxi(VoigtMapping::NonSymToVoigt9(a, b), 1) +=
            xXFsec(a, 3) * invJ(b, 0) + xXFsec(a, 1) * invJ(b, 1) + xXFsec(a, 5) * invJ(b, 2);
        DFDxi(VoigtMapping::NonSymToVoigt9(a, b), 2) +=
            xXFsec(a, 4) * invJ(b, 0) + xXFsec(a, 5) * invJ(b, 1) + xXFsec(a, 2) * invJ(b, 2);
      }

    DsntDxi->MultiplyTN(1.0, DFDxi, DsntDF, 0.0);

    static LINALG::Matrix<DRT::UTILS::DisTypeToNumDeriv2<DRT::Element::tet4>::numderiv2,
        NUMDIM_SOTET4>
        Xsec(true);
    static LINALG::Matrix<NUMNOD_SOTET4,
        DRT::UTILS::DisTypeToNumDeriv2<DRT::Element::tet4>::numderiv2>
        N_XYZ_Xsec(true);
    Xsec.Multiply(1.0, deriv2, xrefe, 0.0);
    N_XYZ_Xsec.MultiplyTT(1.0, N_XYZ, Xsec, 0.0);

    static LINALG::Matrix<9, NUMDOF_SOTET4> D2sntDF2DFDd(true);
    D2sntDF2DFDd.Multiply(1.0, D2sntDF2, DFDd, 0.0);
    D2sntDdDxi_m.MultiplyTN(1.0, D2sntDF2DFDd, DFDxi, 0.0);

    static LINALG::Matrix<9, NUMDIM_SOTET4 * NUMDOF_SOTET4> D2FDxiDd(true);
    D2FDxiDd.Clear();
    for (int m = 0; m < NUMDIM_SOTET4; ++m)
      for (int n = 0; n < NUMDIM_SOTET4; ++n)
        for (int k = 0; k < NUMNOD_SOTET4; ++k)
        {
          D2FDxiDd(
              VoigtMapping::NonSymToVoigt9(m, n), NODDOF_SOTET4 * (NODDOF_SOTET4 * k + m) + 0) +=
              deriv2(0, k) * invJ(n, 0) + deriv2(3, k) * invJ(n, 1) + deriv2(4, k) * invJ(n, 2) -
              N_XYZ_Xsec(k, 0) * invJ(n, 0) - N_XYZ_Xsec(k, 3) * invJ(n, 1) -
              N_XYZ_Xsec(k, 4) * invJ(n, 2);

          D2FDxiDd(
              VoigtMapping::NonSymToVoigt9(m, n), NODDOF_SOTET4 * (NODDOF_SOTET4 * k + m) + 1) +=
              deriv2(3, k) * invJ(n, 0) + deriv2(1, k) * invJ(n, 1) + deriv2(5, k) * invJ(n, 2) -
              N_XYZ_Xsec(k, 3) * invJ(n, 0) - N_XYZ_Xsec(k, 1) * invJ(n, 1) -
              N_XYZ_Xsec(k, 5) * invJ(n, 2);

          D2FDxiDd(
              VoigtMapping::NonSymToVoigt9(m, n), NODDOF_SOTET4 * (NODDOF_SOTET4 * k + m) + 2) +=
              deriv2(4, k) * invJ(n, 0) + deriv2(5, k) * invJ(n, 1) + deriv2(2, k) * invJ(n, 2) -
              N_XYZ_Xsec(k, 4) * invJ(n, 0) - N_XYZ_Xsec(k, 5) * invJ(n, 1) -
              N_XYZ_Xsec(k, 2) * invJ(n, 2);

          for (int l = 0; l < NUMDIM_SOTET4; ++l)
            D2sntDdDxi_m(k * 3 + m, l) += DsntDF(VoigtMapping::NonSymToVoigt9(m, n), 0) *
                                          D2FDxiDd(VoigtMapping::NonSymToVoigt9(m, n),
                                              NODDOF_SOTET4 * (NODDOF_SOTET4 * k + m) + l);
        }
  }

  return;
}

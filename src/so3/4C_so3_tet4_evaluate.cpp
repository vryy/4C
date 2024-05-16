/*----------------------------------------------------------------------*/
/*! \file
\brief quadratic nonlinear tetrahedron
\level 1
*----------------------------------------------------------------------*/

#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_utils_densematrix_eigen.hpp"
#include "4C_linalg_utils_densematrix_inverse.hpp"
#include "4C_mat_constraintmixture.hpp"
#include "4C_mat_elasthyper.hpp"
#include "4C_mat_robinson.hpp"
#include "4C_mat_stvenantkirchhoff.hpp"
#include "4C_mat_thermoplastichyperelast.hpp"
#include "4C_mat_thermostvenantkirchhoff.hpp"
#include "4C_so3_prestress.hpp"
#include "4C_so3_prestress_service.hpp"
#include "4C_so3_tet4.hpp"
#include "4C_so3_utils.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_function.hpp"

#include <Teuchos_SerialDenseSolver.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

// #define PRINT_DEBUG
#ifdef PRINT_DEBUG
#include <sstream>
#include <string>

#include <cstd::string>

template <class T>
void writeArray(const T& mat, std::string name = "unnamed")
{
  std::stringstream header;
  header << 'M' << name << ':' << mat.numRows() << 'x' << mat.numCols() << ':';
  unsigned int s = header.str().size() + mat.numRows() * mat.numCols() * sizeof(double);
  std::cerr.write(reinterpret_cast<const char*>(&s), sizeof(unsigned int));
  std::cerr << header.str();
  for (int i = 0; i < mat.numRows() * mat.numCols(); ++i)
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
#endif

using VoigtMapping = CORE::LINALG::VOIGT::IndexMappings;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              vlf 06/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoTet4::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  // Check whether the solid material PostSetup() routine has already been called and call it if not
  EnsureMaterialPostSetup(params);

  SetParamsInterfacePtr(params);

  CORE::LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4> elemat1(elemat1_epetra.values(), true);
  CORE::LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4> elemat2(elemat2_epetra.values(), true);
  CORE::LINALG::Matrix<NUMDOF_SOTET4, 1> elevec1(elevec1_epetra.values(), true);
  CORE::LINALG::Matrix<NUMDOF_SOTET4, 1> elevec2(elevec2_epetra.values(), true);
  CORE::LINALG::Matrix<NUMDOF_SOTET4, 1> elevec3(elevec3_epetra.values(), true);

  // start with "none"
  DRT::ELEMENTS::SoTet4::ActionType act = SoTet4::none;

  // get the required action
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    FOUR_C_THROW("No action supplied");
  else if (action == "calc_struct_linstiff")
    act = SoTet4::calc_struct_linstiff;
  else if (action == "calc_struct_nlnstiff")
    act = SoTet4::calc_struct_nlnstiff;
  else if (action == "calc_struct_internalforce")
    act = SoTet4::calc_struct_internalforce;
  else if (action == "calc_struct_linstiffmass")
    act = SoTet4::calc_struct_linstiffmass;
  else if (action == "calc_struct_nlnstiffmass")
    act = SoTet4::calc_struct_nlnstiffmass;
  else if (action == "calc_struct_nlnstifflmass")
    act = SoTet4::calc_struct_nlnstifflmass;
  else if (action == "calc_struct_stress")
    act = SoTet4::calc_struct_stress;
  else if (action == "calc_struct_eleload")
    act = SoTet4::calc_struct_eleload;
  else if (action == "calc_struct_fsiload")
    act = SoTet4::calc_struct_fsiload;
  else if (action == "calc_struct_store_istep")
    act = SoTet4::struct_calc_store_istep;
  else if (action == "calc_struct_recover_istep")
    act = SoTet4::struct_calc_recover_istep;
  else if (action == "calc_struct_update_istep")
    act = SoTet4::calc_struct_update_istep;
  else if (action == "calc_struct_reset_istep")
    act = SoTet4::calc_struct_reset_istep;
  else if (action == "calc_struct_reset_all")
    act = SoTet4::calc_struct_reset_all;
  else if (action == "calc_struct_prestress_update")
    act = SoTet4::prestress_update;
  else if (action == "calc_global_gpstresses_map")
    act = SoTet4::calc_global_gpstresses_map;
  else if (action == "calc_struct_energy")
    act = SoTet4::calc_struct_energy;
  else if (action == "calc_struct_output_E")
    act = SoTet4::calc_struct_output_E;
  else if (action == "multi_calc_dens")
    act = SoTet4::multi_calc_dens;
  else if (action == "multi_readrestart")
    act = SoTet4::multi_readrestart;
  else if (action == "calc_struct_recover")
    return 0;
  else if (action == "calc_struct_predict")
    return 0;
  else
    FOUR_C_THROW("Unknown type of action for So_tet4");

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
        FOUR_C_THROW("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      CORE::FE::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      CORE::FE::ExtractMyValues(*res, myres, lm);
      //      CORE::LINALG::Matrix<NUMDOF_SOTET4,NUMDOF_SOTET4>* matptr = nullptr;
      //      if (elemat1.IsInitialized()) matptr = &elemat1;

      std::vector<double> mydispmat(lm.size(), 0.0);

      nlnstiffmass(lm, mydisp, nullptr, nullptr, myres, mydispmat, &elemat1, nullptr, &elevec1,
          nullptr, &elevec3, nullptr, nullptr, nullptr, params, INPAR::STR::stress_none,
          INPAR::STR::strain_none, INPAR::STR::strain_none);
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
        FOUR_C_THROW("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      CORE::FE::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      CORE::FE::ExtractMyValues(*res, myres, lm);
      std::vector<double> mydispmat(lm.size(), 0.0);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      CORE::LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4> myemat(true);  // to zero

      nlnstiffmass(lm, mydisp, nullptr, nullptr, myres, mydispmat, &myemat, nullptr, &elevec1,
          nullptr, nullptr, nullptr, nullptr, nullptr, params, INPAR::STR::stress_none,
          INPAR::STR::strain_none, INPAR::STR::strain_none);
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
        FOUR_C_THROW("Cannot get state vectors 'displacement' and/or residual");
      if (vel == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'velocity'");
      if (acc == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'acceleration'");
      std::vector<double> mydisp(lm.size());
      CORE::FE::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myvel(lm.size());
      CORE::FE::ExtractMyValues(*vel, myvel, lm);
      std::vector<double> myacc(lm.size());
      CORE::FE::ExtractMyValues(*acc, myacc, lm);
      std::vector<double> myres(lm.size());
      CORE::FE::ExtractMyValues(*res, myres, lm);

      std::vector<double> mydispmat(lm.size(), 0.0);

      nlnstiffmass(lm, mydisp, &myvel, &myacc, myres, mydispmat, &elemat1, &elemat2, &elevec1,
          &elevec2, &elevec3, nullptr, nullptr, nullptr, params, INPAR::STR::stress_none,
          INPAR::STR::strain_none, INPAR::STR::strain_none);

      if (act == calc_struct_nlnstifflmass) so_tet4_lumpmass(&elemat2);
    }
    break;

    //==================================================================================
    // evaluate stresses and strains at gauss points
    case calc_struct_stress:
    {
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      Teuchos::RCP<std::vector<char>> stressdata =
          params.get<Teuchos::RCP<std::vector<char>>>("stress", Teuchos::null);
      Teuchos::RCP<std::vector<char>> straindata =
          params.get<Teuchos::RCP<std::vector<char>>>("strain", Teuchos::null);
      if (disp == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'displacement'");
      if (stressdata == Teuchos::null) FOUR_C_THROW("Cannot get 'stress' data");
      if (straindata == Teuchos::null) FOUR_C_THROW("Cannot get 'strain' data");
      std::vector<double> mydisp(lm.size());
      CORE::FE::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      CORE::FE::ExtractMyValues(*res, myres, lm);
      CORE::LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D> stress(true);  // set to zero
      CORE::LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D> strain(true);
      auto iostress = CORE::UTILS::GetAsEnum<INPAR::STR::StressType>(
          params, "iostress", INPAR::STR::stress_none);
      auto iostrain = CORE::UTILS::GetAsEnum<INPAR::STR::StrainType>(
          params, "iostrain", INPAR::STR::strain_none);

      std::vector<double> mydispmat(lm.size(), 0.0);

      nlnstiffmass(lm, mydisp, nullptr, nullptr, myres, mydispmat, nullptr, nullptr, nullptr,
          nullptr, nullptr, &stress, &strain, nullptr, params, iostress, iostrain,
          INPAR::STR::strain_none);

      {
        CORE::COMM::PackBuffer data;
        AddtoPack(data, stress);
        data.StartPacking();
        AddtoPack(data, stress);
        std::copy(data().begin(), data().end(), std::back_inserter(*stressdata));
      }
      {
        CORE::COMM::PackBuffer data;
        AddtoPack(data, strain);
        data.StartPacking();
        AddtoPack(data, strain);
        std::copy(data().begin(), data().end(), std::back_inserter(*straindata));
      }
    }
    break;

    //==================================================================================
    case prestress_update:
    {
      time_ = params.get<double>("total time");
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) FOUR_C_THROW("Cannot get displacement state");
      std::vector<double> mydisp(lm.size());
      CORE::FE::ExtractMyValues(*disp, mydisp, lm);

      switch (pstype_)
      {
        case INPAR::STR::PreStress::mulf:
        {
          // build incremental def gradient for every gauss point
          CORE::LINALG::SerialDenseMatrix gpdefgrd(NUMGPT_SOTET4, 9);
          DefGradient(mydisp, gpdefgrd, *prestress_);

          // update deformation gradient and put back to storage
          CORE::LINALG::Matrix<3, 3> deltaF;
          CORE::LINALG::Matrix<3, 3> Fhist;
          CORE::LINALG::Matrix<3, 3> Fnew;
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
          if (Material()->MaterialType() == CORE::Materials::m_constraintmixture)
          {
            SolidMaterial()->Update();
          }
          break;
        }
        default:
          FOUR_C_THROW(
              "You should either not be here, or the prestressing method you are using is not "
              "implemented for TET4 elements!");
      }
    }
    break;
    //==================================================================================
    // this is a dummy output for strain energy
    case calc_struct_energy:
    {
      // initialization of internal energy
      double intenergy = 0.0;

      const static std::vector<CORE::LINALG::Matrix<NUMNOD_SOTET4, 1>> shapefcts =
          so_tet4_1gp_shapefcts();
      const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET4 + 1, NUMNOD_SOTET4>> derivs =
          so_tet4_1gp_derivs();
      const static std::vector<double> gpweights = so_tet4_1gp_weights();

      // get displacements of this processor
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) FOUR_C_THROW("Cannot get state displacement vector");

      // get displacements of this element
      std::vector<double> mydisp(lm.size());
      CORE::FE::ExtractMyValues(*disp, mydisp, lm);

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
      CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xrefe;
      CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xdisp;
      DRT::Node** nodes = Nodes();
      for (int i = 0; i < NUMNOD_SOTET4; ++i)
      {
        const auto& x = nodes[i]->X();
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
        CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> nxyz(nxyz_);  // copy!

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
        CORE::LINALG::Matrix<3, 3> defgrd(true);
        // Gauss weights and Jacobian determinant
        double fac = detJ * gpweights[gp];

        if (PRESTRESS::IsMulf(pstype_))
        {
          // get derivatives wrt to last spatial configuration
          CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> N_xyz;
          prestress_->StoragetoMatrix(gp, N_xyz, prestress_->JHistory());

          // build multiplicative incremental defgrd
          // defgrd.Multiply('T','N',1.0,xdisp,N_xyz,0.0);
          if (kintype_ == INPAR::STR::KinemType::nonlinearTotLag)
          {
            defgrd.MultiplyTN(xdisp, N_xyz);
          }
          defgrd(0, 0) += 1.0;
          defgrd(1, 1) += 1.0;
          defgrd(2, 2) += 1.0;

          // get stored old incremental F
          CORE::LINALG::Matrix<3, 3> Fhist;
          prestress_->StoragetoMatrix(gp, Fhist, prestress_->FHistory());

          // build total defgrd = delta F * F_old
          CORE::LINALG::Matrix<3, 3> Fnew;
          Fnew.Multiply(defgrd, Fhist);
          defgrd = Fnew;
        }
        else
        {
          // in kinematically linear analysis the deformation gradient is equal to identity
          if (kintype_ == INPAR::STR::KinemType::nonlinearTotLag)
          {
            defgrd.MultiplyTN(xdisp, nxyz);
          }
          defgrd(0, 0) += 1.0;
          defgrd(1, 1) += 1.0;
          defgrd(2, 2) += 1.0;
        }

        // Right Cauchy-Green tensor = F^T * F
        // size is 3x3
        CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> cauchygreen;
        cauchygreen.MultiplyTN(defgrd, defgrd);

        // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
        // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
        CORE::LINALG::Matrix<6, 1> glstrain(false);
        glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
        glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
        glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
        glstrain(3) = cauchygreen(0, 1);
        glstrain(4) = cauchygreen(1, 2);
        glstrain(5) = cauchygreen(2, 0);

        // call material for evaluation of strain energy function
        double psi = 0.0;
        SolidMaterial()->StrainEnergy(glstrain, psi, gp, Id());

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
        if (elevec1_epetra.length() < 1) FOUR_C_THROW("The given result vector is too short.");

        elevec1_epetra(0) = intenergy;
      }
    }
    break;
    // this is needed by bone topology optimization
    case calc_struct_output_E:
    {
      Teuchos::RCP<CORE::MAT::Material> mat = Material();
      // check length of elevec1
      if (elevec1_epetra.length() < 1) FOUR_C_THROW("The given result vector is too short.");
      CORE::MAT::Material* rawmat = mat.get();
      auto* stvk = dynamic_cast<MAT::StVenantKirchhoff*>(rawmat);
      if (!stvk) FOUR_C_THROW("dynamic cast to stvenant failed");
      double E = stvk->Youngs();
      elevec1_epetra(0) = E;
    }
    break;

    //==================================================================================
    case calc_struct_eleload:
      FOUR_C_THROW("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
      break;

    //==================================================================================
    case calc_struct_fsiload:
      FOUR_C_THROW("Case not yet implemented");
      break;

    //==================================================================================
    case struct_calc_store_istep:
    {
      int timestep = params.get<int>("timestep", -1);

      if (timestep == -1) FOUR_C_THROW("Provide timestep number to be stored");

      // due to the multiplicativity and futility to redo prestress steps
      // other than the last one, no need to store/recover anything
      // ... but keep in mind
      if (PRESTRESS::IsAny(pstype_))
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

      if (timestep == -1) FOUR_C_THROW("Provide timestep number of the timestep to be recovered");

      // due to the multiplicativity and futility to redo prestress steps
      // other than the last one, no need to store/recover anything
      // ... but keep in mind
      if (PRESTRESS::IsAny(pstype_))
      {
      }

      // Material
      SolidMaterial()->SetHistory(timestep);
    }
    break;

    //==================================================================================
    case calc_struct_update_istep:
    {
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      CORE::FE::ExtractMyValues(*disp, mydisp, lm);

      if (SolidMaterial()->UsesExtendedUpdate())
      {
        CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMNOD_SOTET4> xdispT(mydisp.data());

        // build deformation gradient wrt to material configuration
        CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> defgrd(false);
        for (unsigned gp = 0; gp < NUMGPT_SOTET4; ++gp)
        {
          // Compute deformation gradient
          ComputeDeformationGradient(defgrd, xdispT, gp);

          // call material update if material = m_growthremodel_elasthyper (calculate and update
          // inelastic deformation gradient)
          if (SolidMaterial()->UsesExtendedUpdate())
          {
            SolidMaterial()->Update(defgrd, gp, params, Id());
          }
        }
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
        if (disp == Teuchos::null) FOUR_C_THROW("Cannot get state vectors 'displacement'");
        if (stressdata == Teuchos::null) FOUR_C_THROW("Cannot get 'stress' data");
        if (straindata == Teuchos::null) FOUR_C_THROW("Cannot get 'strain' data");
        const Teuchos::RCP<std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>>
            gpstressmap = params.get<
                Teuchos::RCP<std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>>>(
                "gpstressmap", Teuchos::null);
        if (gpstressmap == Teuchos::null)
          FOUR_C_THROW("no gp stress map available for writing gpstresses");
        const Teuchos::RCP<std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>>
            gpstrainmap = params.get<
                Teuchos::RCP<std::map<int, Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>>>>(
                "gpstrainmap", Teuchos::null);
        if (gpstrainmap == Teuchos::null)
          FOUR_C_THROW("no gp strain map available for writing gpstrains");
        std::vector<double> mydisp(lm.size());
        CORE::FE::ExtractMyValues(*disp, mydisp, lm);
        std::vector<double> myres(lm.size());
        CORE::FE::ExtractMyValues(*res, myres, lm);
        CORE::LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D> stress;
        CORE::LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D> strain;
        auto iostress = CORE::UTILS::GetAsEnum<INPAR::STR::StressType>(
            params, "iostress", INPAR::STR::stress_none);
        auto iostrain = CORE::UTILS::GetAsEnum<INPAR::STR::StrainType>(
            params, "iostrain", INPAR::STR::strain_none);

        std::vector<double> mydispmat(lm.size(), 0.0);

        // if a linear analysis is desired
        if (kintype_ == INPAR::STR::KinemType::linear)
        {
          FOUR_C_THROW("Linear case not implemented");
        }

        else
        {
          nlnstiffmass(lm, mydisp, nullptr, nullptr, myres, mydispmat, nullptr, nullptr, nullptr,
              nullptr, nullptr, &stress, &strain, nullptr, params, iostress, iostrain,
              INPAR::STR::strain_none);
        }
        // add stresses to global map
        // get EleID Id()
        int gid = Id();
        Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> gpstress =
            Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix);
        gpstress->shape(NUMGPT_SOTET4, MAT::NUM_STRESS_3D);

        // move stresses to serial dense matrix
        for (int i = 0; i < NUMGPT_SOTET4; i++)
        {
          for (int j = 0; j < MAT::NUM_STRESS_3D; j++)
          {
            (*gpstress)(i, j) = stress(i, j);
          }
        }

        // strains
        Teuchos::RCP<CORE::LINALG::SerialDenseMatrix> gpstrain =
            Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix);
        gpstrain->shape(NUMGPT_SOTET4, MAT::NUM_STRESS_3D);

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
          CORE::COMM::PackBuffer data;
          AddtoPack(data, stress);
          data.StartPacking();
          AddtoPack(data, stress);
          std::copy(data().begin(), data().end(), std::back_inserter(*stressdata));
        }

        {
          CORE::COMM::PackBuffer data;
          AddtoPack(data, strain);
          data.StartPacking();
          AddtoPack(data, strain);
          std::copy(data().begin(), data().end(), std::back_inserter(*straindata));
        }
      }
    }
    break;

    default:
      FOUR_C_THROW("Unknown type of action for so_tet4");
      break;
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)     maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::SoTet4::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, CORE::Conditions::Condition& condition,
    std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
    CORE::LINALG::SerialDenseMatrix* elemat1)
{
  // get values and switches from the condition
  const auto* onoff = &condition.parameters().Get<std::vector<int>>("onoff");
  const auto* val = &condition.parameters().Get<std::vector<double>>("val");

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  const double time = std::invoke(
      [&]()
      {
        if (IsParamsInterface())
          return StrParamsInterface().GetTotalTime();
        else
          return params.get("total time", -1.0);
      });

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < NUMDIM_SOTET4)
    FOUR_C_THROW("Fewer functions or curves defined than the element has dofs.");

  for (int checkdof = NUMDIM_SOTET4; checkdof < int(onoff->size()); ++checkdof)
  {
    if ((*onoff)[checkdof] != 0)
      FOUR_C_THROW(
          "Number of Dimensions in Neumann_Evalutaion is 3. Further DoFs are not considered.");
  }

  // (SPATIAL) FUNCTION BUSINESS
  static_assert((NUMGPT_SOTET4 == 1));
  const auto* funct = &condition.parameters().Get<std::vector<int>>("funct");
  CORE::LINALG::Matrix<NUMDIM_SOTET4, 1> xrefegp(false);
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
  const static std::vector<CORE::LINALG::Matrix<NUMNOD_SOTET4, 1>> shapefcts =
      so_tet4_1gp_shapefcts();
  const static std::vector<double> gpweights = so_tet4_1gp_weights();
  /* ============================================================================*/

  // update element geometry
  CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xrefe;
  DRT::Node** nodes = Nodes();
  for (int i = 0; i < NUMNOD_SOTET4; ++i)
  {
    const auto& x = nodes[i]->X();
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
  CORE::LINALG::Matrix<NUMCOORD_SOTET4, NUMCOORD_SOTET4> jac;
  for (int i = 0; i < 4; i++) jac(0, i) = 1.0;
  for (int row = 0; row < 3; row++)
    for (int col = 0; col < 4; col++) jac(row + 1, col) = xrefe(col, row);

  // compute determinant of Jacobian once outside Gauss point loop since it is constant
  // jac.Determinant() delivers six times the reference volume of the tet
  const double detJ = jac.Determinant() * (1.0 / 6.0);

  if (detJ == 0.0)
    FOUR_C_THROW("ZERO JACOBIAN DETERMINANT");
  else if (detJ < 0.0)
    FOUR_C_THROW("NEGATIVE JACOBIAN DETERMINANT");

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
            (functnum > 0) ? GLOBAL::Problem::Instance()
                                 ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(functnum - 1)
                                 .Evaluate(xrefegp.A(), time, dim)
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
void DRT::ELEMENTS::SoTet4::InitJacobianMapping()
{
  CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xrefe;
  DRT::Node** nodes = Nodes();
  for (int i = 0; i < NUMNOD_SOTET4; ++i)
  {
    const auto& x = nodes[i]->X();
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
  CORE::LINALG::Matrix<NUMCOORD_SOTET4, NUMCOORD_SOTET4> jac;
  for (int i = 0; i < 4; i++) jac(0, i) = 1;
  for (int row = 0; row < 3; row++)
    for (int col = 0; col < 4; col++) jac(row + 1, col) = xrefe(col, row);
  // volume of the element
  V_ = jac.Determinant() / 6.0;
  if (V_ <= 0.0) FOUR_C_THROW("Element volume %10.5e <= 0.0 (Id: %i)", V_, Id());

  // nxyz_.resize(NUMGPT_SOTET4);
  const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET4 + 1, NUMNOD_SOTET4>> derivs =
      so_tet4_1gp_derivs();
  CORE::LINALG::Matrix<NUMCOORD_SOTET4 - 1, NUMCOORD_SOTET4> tmp;
  for (int gp = 0; gp < NUMGPT_SOTET4; ++gp)
  {
    tmp.MultiplyTN(xrefe, derivs[gp]);
    for (int i = 0; i < 4; i++) jac(0, i) = 1;
    for (int row = 0; row < 3; row++)
      for (int col = 0; col < 4; col++) jac(row + 1, col) = tmp(row, col);
    // size is 4x3
    CORE::LINALG::Matrix<NUMCOORD_SOTET4, NUMDIM_SOTET4> I_aug(true);
    // size is 4x3
    CORE::LINALG::Matrix<NUMCOORD_SOTET4, NUMDIM_SOTET4> partials(true);
    I_aug(1, 0) = 1;
    I_aug(2, 1) = 1;
    I_aug(3, 2) = 1;

    // solve A.X=B
    CORE::LINALG::FixedSizeSerialDenseSolver<NUMCOORD_SOTET4, NUMCOORD_SOTET4, NUMDIM_SOTET4>
        solve_for_inverseJac;
    solve_for_inverseJac.SetMatrix(jac);               // set A=jac
    solve_for_inverseJac.SetVectors(partials, I_aug);  // set X=partials, B=I_aug
    solve_for_inverseJac.FactorWithEquilibration(true);
    int err2 = solve_for_inverseJac.Factor();
    int err = solve_for_inverseJac.Solve();  // partials = jac^-1.I_aug
    if ((err != 0) || (err2 != 0)) FOUR_C_THROW("Inversion of Jacobian failed");

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

    if (PRESTRESS::IsMulfActive(time_, pstype_, pstime_))
    {
      if (!(prestress_->IsInit())) prestress_->MatrixtoStorage(gp, nxyz_, prestress_->JHistory());
    }

  }  // for (int gp=0; gp<NUMGPT_SOTET4; ++gp)

  if (PRESTRESS::IsMulfActive(time_, pstype_, pstime_)) prestress_->IsInit() = true;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (private)                            vlf 08/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoTet4::nlnstiffmass(std::vector<int>& lm,  // location matrix
    std::vector<double>& disp,                                  // current displacements
    std::vector<double>* vel,                                   // current velocities
    std::vector<double>* acc,                                   // current accelerations
    std::vector<double>& residual,                              // current residual displ
    std::vector<double>& dispmat,                               // current material displacements
    CORE::LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4>* stiffmatrix,  // element stiffness matrix
    CORE::LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4>* massmatrix,   // element mass matrix
    CORE::LINALG::Matrix<NUMDOF_SOTET4, 1>* force,       // element internal force vector
    CORE::LINALG::Matrix<NUMDOF_SOTET4, 1>* forceinert,  // element inertial force vector
    CORE::LINALG::Matrix<NUMDOF_SOTET4, 1>* force_str,   // element structural force vector
    CORE::LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D>* elestress,    // stresses at GP
    CORE::LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D>* elestrain,    // strains at GP
    CORE::LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D>* eleplstrain,  // plastic strains at GP
    Teuchos::ParameterList& params,          // algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,   // stress output option
    const INPAR::STR::StrainType iostrain,   // strain output option
    const INPAR::STR::StrainType ioplstrain  // plastic strain output option
)
{
  /* =============================================================================*
  ** CONST DERIVATIVES and WEIGHTS for TET_4  with 1 GAUSS POINTS*
  ** =============================================================================*/
  const static std::vector<CORE::LINALG::Matrix<NUMNOD_SOTET4, 1>> shapefcts =
      so_tet4_1gp_shapefcts();
  const static std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET4 + 1, NUMNOD_SOTET4>> derivs =
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
  CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xrefe;
  DRT::Node** nodes = Nodes();
  for (int i = 0; i < NUMNOD_SOTET4; ++i)
  {
    const auto& x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }
  CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMNOD_SOTET4> xdispT(disp.data());


  // volume of a tetrahedra
  double detJ = V_;


  // size is 3x3
  CORE::LINALG::Matrix<3, 3> defgrd(true);
  /* =========================================================================*/
  /* ============================================== Loop over Gauss Points ===*/
  /* =========================================================================*/
  for (int gp = 0; gp < NUMGPT_SOTET4; gp++)
  {
    CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> nxyz(nxyz_);  // copy!

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

    // Evaluate deformation gradient
    ComputeDeformationGradient(defgrd, xdispT, gp);

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
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOTET4> bop;
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
    CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(false);
    if (kintype_ == INPAR::STR::KinemType::nonlinearTotLag)
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
      CORE::LINALG::Matrix<NUMDOF_SOTET4, 1> nodaldisp;
      for (int i = 0; i < NUMNOD_SOTET4; ++i)
      {
        nodaldisp(3 * i, 0) = xdispT(0, i);
        nodaldisp(3 * i + 1, 0) = xdispT(1, i);
        nodaldisp(3 * i + 2, 0) = xdispT(2, i);
      }

      // build the linearised strain epsilon = B_L . d
      glstrain.Multiply(bop, nodaldisp);
    }

    // return gp strains (only in case of stress/strain output)
    switch (iostrain)
    {
      case INPAR::STR::strain_gl:
      {
        if (elestrain == nullptr) FOUR_C_THROW("no strain data available");
        for (int i = 0; i < 3; ++i) (*elestrain)(gp, i) = glstrain(i);
        for (int i = 3; i < 6; ++i) (*elestrain)(gp, i) = 0.5 * glstrain(i);
      }
      break;
      case INPAR::STR::strain_ea:
      {
        if (elestrain == nullptr) FOUR_C_THROW("no strain data available");

        // rewriting Green-Lagrange strains in matrix format
        CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> gl;
        gl(0, 0) = glstrain(0);
        gl(0, 1) = 0.5 * glstrain(3);
        gl(0, 2) = 0.5 * glstrain(5);
        gl(1, 0) = gl(0, 1);
        gl(1, 1) = glstrain(1);
        gl(1, 2) = 0.5 * glstrain(4);
        gl(2, 0) = gl(0, 2);
        gl(2, 1) = gl(1, 2);
        gl(2, 2) = glstrain(2);

        // Inverse of deformation gradient (make a copy here otherwise defgrd is destroyed).
        CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> invdefgrd;
        invdefgrd.Invert(defgrd);

        CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> temp;
        CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> euler_almansi;
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
        if (elestrain == nullptr) FOUR_C_THROW("strain data not available");

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
        CORE::LINALG::Matrix<3, 3> prstr2(true);  // squared principal stretches
        CORE::LINALG::Matrix<3, 1> prstr(true);   // principal stretch
        CORE::LINALG::Matrix<3, 3> prdir(true);   // principal directions
        CORE::LINALG::SYEV(cauchygreen, prstr2, prdir);

        // THE principal stretches
        for (int al = 0; al < 3; ++al) prstr(al) = std::sqrt(prstr2(al, al));

        // populating the logarithmic strain matrix
        CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> lnv(true);

        // checking if cauchy green is correctly determined to ensure eigen vectors in correct
        // direction i.e. a flipped eigenvector is also a valid solution C = \sum_{i=1}^3
        // (\lambda_i^2) \mathbf{n}_i \otimes \mathbf{n}_i
        CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> tempCG(true);

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
        CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> diffCG(true);

        for (int i = 0; i < 3; ++i)
        {
          for (int j = 0; j < 3; ++j)
          {
            diffCG(i, j) = cauchygreen(i, j) - tempCG(i, j);
            // the solution to this problem is to evaluate the cauchygreen tensor with tempCG
            // computed with every combination of eigenvector orientations -- up to nine comparisons
            if (diffCG(i, j) > 1e-10)
              FOUR_C_THROW(
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
        FOUR_C_THROW("requested strain option not available");
        break;
    }

    // call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> cmat(true);
    CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> stress(true);

    if (Material()->MaterialType() == CORE::Materials::m_constraintmixture ||
        Material()->MaterialType() == CORE::Materials::m_mixture)
    {
      // gp reference coordinates
      CORE::LINALG::Matrix<NUMNOD_SOTET4, 1> funct(true);
      funct = shapefcts[gp];
      CORE::LINALG::Matrix<NUMDIM_SOTET4, 1> point(true);
      point.MultiplyTN(xrefe, funct);
      params.set("gp_coords_ref", point);
    }

    UTILS::GetTemperatureForStructuralMaterial<CORE::FE::CellType::tet4>(shapefcts[gp], params);

    SolidMaterial()->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, gp, Id());

    // return gp stresses
    switch (iostress)
    {
      case INPAR::STR::stress_2pk:
      {
        if (elestress == nullptr) FOUR_C_THROW("no stress data available");
        for (int i = 0; i < MAT::NUM_STRESS_3D; ++i) (*elestress)(gp, i) = stress(i);
      }
      break;
      case INPAR::STR::stress_cauchy:
      {
        if (elestress == nullptr) FOUR_C_THROW("no stress data available");
        double detF = defgrd.Determinant();

        CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> pkstress;
        pkstress(0, 0) = stress(0);
        pkstress(0, 1) = stress(3);
        pkstress(0, 2) = stress(5);
        pkstress(1, 0) = pkstress(0, 1);
        pkstress(1, 1) = stress(1);
        pkstress(1, 2) = stress(4);
        pkstress(2, 0) = pkstress(0, 2);
        pkstress(2, 1) = pkstress(1, 2);
        pkstress(2, 2) = stress(2);

        CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> temp;
        CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> cauchystress;
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
        FOUR_C_THROW("requested stress type not available");
        break;
    }

    double detJ_w = detJ * (gpweights)[gp];

    // update of internal force vector
    if (force != nullptr)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w, bop, stress, 1.0);
    }

    // update of stiffness matrix
    if (stiffmatrix != nullptr)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      // size is 6x12
      CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOTET4> cb;
      cb.Multiply(cmat, bop);  // temporary C . B
      // size is 12x12
      stiffmatrix->MultiplyTN(detJ_w, bop, cb, 1.0);

      if (kintype_ == INPAR::STR::KinemType::nonlinearTotLag)
      {
        // integrate `geometric' stiffness matrix and add to keu
        // auxiliary integrated stress
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> sfac(stress);
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
  const static std::vector<CORE::LINALG::Matrix<NUMNOD_SOTET4, 1>> shapefcts4gp =
      so_tet4_4gp_shapefcts();
  const static std::vector<double> gpweights4gp = so_tet4_4gp_weights();
  // evaluate mass matrix
  if (massmatrix != nullptr)
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

        CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> nxyz(nxyz_);  // copy!

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
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> linmass_disp(true);
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> linmass_vel(true);
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> linmass(true);

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
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOTET4> bop;
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
        CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> cauchygreen;
        cauchygreen.MultiplyTN(defgrd, defgrd);

        // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
        // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
        CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(false);
        if (kintype_ == INPAR::STR::KinemType::nonlinearTotLag)
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
          CORE::LINALG::Matrix<NUMDOF_SOTET4, 1> nodaldisp;
          for (int i = 0; i < NUMNOD_SOTET4; ++i)
          {
            nodaldisp(3 * i, 0) = xdispT(0, i);
            nodaldisp(3 * i + 1, 0) = xdispT(1, i);
            nodaldisp(3 * i + 2, 0) = xdispT(2, i);
          }

          // build the linearised strain epsilon = B_L . d
          glstrain.Multiply(bop, nodaldisp);
        }

        // evaluate derivative of mass w.r.t. to right cauchy green tensor
        SolidMaterial()->EvaluateNonLinMass(
            &defgrd, &glstrain, params, &linmass_disp, &linmass_vel, gp, Id());

        // multiply by 2.0 to get derivative w.r.t green lagrange strains and multiply by time
        // integration factor
        linmass_disp.Scale(2.0 * timintfac_dis);
        linmass_vel.Scale(2.0 * timintfac_vel);
        linmass.Update(1.0, linmass_disp, 1.0, linmass_vel, 0.0);

        // evaluate accelerations at time n+1 at gauss point
        CORE::LINALG::Matrix<NUMDIM_SOTET4, 1> myacc(true);
        for (int idim = 0; idim < NUMDIM_SOTET4; ++idim)
          for (int inod = 0; inod < NUMNOD_SOTET4; ++inod)
            myacc(idim) += shapefcts4gp[gp](inod) * (*acc)[idim + (inod * NUMDIM_SOTET4)];

        if (stiffmatrix != nullptr)
        {
          // integrate linearisation of mass matrix
          //(B^T . d\rho/d disp . a) * detJ * w(gp)
          CORE::LINALG::Matrix<1, NUMDOF_SOTET4> cb;
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
        if (forceinert != nullptr)
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
void DRT::ELEMENTS::SoTet4::so_tet4_lumpmass(
    CORE::LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4>* emass)
{
  // lump mass matrix
  if (emass != nullptr)
  {
    // we assume #elemat2 is a square matrix
    for (unsigned c = 0; c < (*emass).numCols(); ++c)  // parse columns
    {
      double d = 0.0;
      for (unsigned r = 0; r < (*emass).numRows(); ++r)  // parse rows
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
int DRT::ELEMENTS::SoTet4Type::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    auto* actele = dynamic_cast<DRT::ELEMENTS::SoTet4*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to So_tet4* failed");
    actele->InitJacobianMapping();
  }
  return 0;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet4 Shape fcts at 1 Gauss Point                           |
 *----------------------------------------------------------------------*/
std::vector<CORE::LINALG::Matrix<NUMNOD_SOTET4, 1>> DRT::ELEMENTS::SoTet4::so_tet4_1gp_shapefcts()
{
  std::vector<CORE::LINALG::Matrix<NUMNOD_SOTET4, 1>> shapefcts(NUMGPT_SOTET4);

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
std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET4 + 1, NUMNOD_SOTET4>>
DRT::ELEMENTS::SoTet4::so_tet4_1gp_derivs()
{
  std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET4 + 1, NUMNOD_SOTET4>> derivs(NUMGPT_SOTET4);
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
std::vector<double> DRT::ELEMENTS::SoTet4::so_tet4_1gp_weights()
{
  std::vector<double> weights(NUMGPT_SOTET4);
  // There is only one gausspoint, so the loop (and the vector) is not really needed.
  for (int i = 0; i < NUMGPT_SOTET4; ++i) weights[i] = 1.0;
  return weights;
}

/*----------------------------------------------------------------------*
 |  Evaluate Tet4 Shape fcts at 4 Gauss Points                          |
 *----------------------------------------------------------------------*/
std::vector<CORE::LINALG::Matrix<NUMNOD_SOTET4, 1>> DRT::ELEMENTS::SoTet4::so_tet4_4gp_shapefcts()
{
  std::vector<CORE::LINALG::Matrix<NUMNOD_SOTET4, 1>> shapefcts(4);

  const double gploc_alpha =
      (5.0 + 3.0 * sqrt(5.0)) / 20.0;  // gp sampling point value for quadr. fct
  const double gploc_beta = (5.0 - sqrt(5.0)) / 20.0;

  const std::array<double, 4> xsi1 = {gploc_alpha, gploc_beta, gploc_beta, gploc_beta};
  const std::array<double, 4> xsi2 = {gploc_beta, gploc_alpha, gploc_beta, gploc_beta};
  const std::array<double, 4> xsi3 = {gploc_beta, gploc_beta, gploc_alpha, gploc_beta};
  const std::array<double, 4> xsi4 = {gploc_beta, gploc_beta, gploc_beta, gploc_alpha};

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
std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET4 + 1, NUMNOD_SOTET4>>
DRT::ELEMENTS::SoTet4::so_tet4_4gp_derivs()
{
  std::vector<CORE::LINALG::Matrix<NUMDIM_SOTET4 + 1, NUMNOD_SOTET4>> derivs(4);

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
std::vector<double> DRT::ELEMENTS::SoTet4::so_tet4_4gp_weights()
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
void DRT::ELEMENTS::SoTet4::DefGradient(const std::vector<double>& disp,
    CORE::LINALG::SerialDenseMatrix& gpdefgrd, DRT::ELEMENTS::PreStress& prestress)
{
  // update element geometry
  CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xdisp;
  for (int i = 0; i < NUMNOD_SOTET4; ++i)
  {
    xdisp(i, 0) = disp[i * NODDOF_SOTET4 + 0];
    xdisp(i, 1) = disp[i * NODDOF_SOTET4 + 1];
    xdisp(i, 2) = disp[i * NODDOF_SOTET4 + 2];
  }

  for (int gp = 0; gp < NUMGPT_SOTET4; ++gp)
  {
    // get derivatives wrt to last spatial configuration
    CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> N_xyz;
    prestress_->StoragetoMatrix(gp, N_xyz, prestress_->JHistory());

    // build multiplicative incremental defgrd
    CORE::LINALG::Matrix<3, 3> defgrd(true);
    if (kintype_ == INPAR::STR::KinemType::nonlinearTotLag)
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

void DRT::ELEMENTS::SoTet4::ComputeDeformationGradient(
    CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4>& defgrd,
    const CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMNOD_SOTET4>& xdisp, const int gp) const
{
  if (kintype_ == INPAR::STR::KinemType::linear)
  {
    // in the linear case, the deformation gradient is the identity matrix
    defgrd.Clear();
    defgrd(0, 0) += 1.0;
    defgrd(1, 1) += 1.0;
    defgrd(2, 2) += 1.0;

    return;
  }

  if (pstype_ == INPAR::STR::PreStress::mulf)
  {
    // get derivatives wrt to last spatial configuration
    CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> N_xyz;
    prestress_->StoragetoMatrix(gp, N_xyz, prestress_->JHistory());

    // build multiplicative incremental defgrd
    CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> Finc;

    Finc.MultiplyNN(xdisp, N_xyz);

    // build multiplicative incremental defgrd
    Finc(0, 0) += 1.0;
    Finc(1, 1) += 1.0;
    Finc(2, 2) += 1.0;

    // get stored old incremental F
    CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> Fhist;
    prestress_->StoragetoMatrix(gp, Fhist, prestress_->FHistory());

    // build total defgrd = delta F * F_old
    defgrd.Multiply(Finc, Fhist);
  }
  else
  {
    // (material) deformation gradient F = d xcurr / d xrefe = I + xdisp * N_XYZ^T
    defgrd.MultiplyNN(xdisp, nxyz_);
    defgrd(0, 0) += 1.0;
    defgrd(1, 1) += 1.0;
    defgrd(2, 2) += 1.0;
  }
}

/*----------------------------------------------------------------------*
 |  compute Jac.mapping wrt deformed configuration (protected) gee 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoTet4::UpdateJacobianMapping(
    const std::vector<double>& disp, DRT::ELEMENTS::PreStress& prestress)
{
  // get incremental disp
  CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xdisp;
  for (int i = 0; i < NUMNOD_SOTET4; ++i)
  {
    xdisp(i, 0) = disp[i * NODDOF_SOTET4 + 0];
    xdisp(i, 1) = disp[i * NODDOF_SOTET4 + 1];
    xdisp(i, 2) = disp[i * NODDOF_SOTET4 + 2];
  }

  CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> nxyzhist;
  CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> nxyznew;
  CORE::LINALG::Matrix<3, 3> defgrd(true);

  for (int gp = 0; gp < NUMGPT_SOTET4; ++gp)
  {
    // get the nxyz old state
    prestress.StoragetoMatrix(gp, nxyzhist, prestress.JHistory());
    // build multiplicative incremental defgrd
    if (kintype_ == INPAR::STR::KinemType::nonlinearTotLag)
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
void DRT::ELEMENTS::SoTet4::so_tet4_remodel(std::vector<int>& lm,  // location matrix
    std::vector<double>& disp,                                     // current displacements
    Teuchos::ParameterList& params,                // algorithmic parameters e.g. time
    const Teuchos::RCP<CORE::MAT::Material>& mat)  // material
{
  if ((Material()->MaterialType() == CORE::Materials::m_constraintmixture) ||
      (Material()->MaterialType() == CORE::Materials::m_elasthyper))
  {
    // in a first step ommit everything with prestress

    // current  displacements of element
    CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xdisp;
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
    CORE::LINALG::Matrix<3, 3> avg_stress(true);
    CORE::LINALG::Matrix<3, 3> avg_defgrd(true);

    for (int gp = 0; gp < NUMGPT_SOTET4; gp++)
    {
      const CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4>& nxyz = nxyz_;

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
      CORE::LINALG::Matrix<3, 3> defgrd(false);

      if (PRESTRESS::IsMulf(pstype_))
      {
        // get derivatives wrt to last spatial configuration
        CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> N_xyz;
        prestress_->StoragetoMatrix(gp, N_xyz, prestress_->JHistory());

        // build multiplicative incremental defgrd
        // defgrd.Multiply('T','N',1.0,xdisp,N_xyz,0.0);
        defgrd.MultiplyTN(xdisp, N_xyz);
        defgrd(0, 0) += 1.0;
        defgrd(1, 1) += 1.0;
        defgrd(2, 2) += 1.0;

        // get stored old incremental F
        CORE::LINALG::Matrix<3, 3> Fhist;
        prestress_->StoragetoMatrix(gp, Fhist, prestress_->FHistory());

        // build total defgrd = delta F * F_old
        CORE::LINALG::Matrix<3, 3> Fnew;
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
      CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> cauchygreen;
      cauchygreen.MultiplyTN(defgrd, defgrd);

      // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
      // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
      CORE::LINALG::Matrix<6, 1> glstrain(false);
      glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
      glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
      glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
      glstrain(3) = cauchygreen(0, 1);
      glstrain(4) = cauchygreen(1, 2);
      glstrain(5) = cauchygreen(2, 0);

      // call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> cmat(true);
      CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1> stress(true);

      SolidMaterial()->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, gp, Id());
      // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

      // Cauchy stress
      const double detF = defgrd.Determinant();

      CORE::LINALG::Matrix<3, 3> pkstress;
      pkstress(0, 0) = stress(0);
      pkstress(0, 1) = stress(3);
      pkstress(0, 2) = stress(5);
      pkstress(1, 0) = pkstress(0, 1);
      pkstress(1, 1) = stress(1);
      pkstress(1, 2) = stress(4);
      pkstress(2, 0) = pkstress(0, 2);
      pkstress(2, 1) = pkstress(1, 2);
      pkstress(2, 2) = stress(2);

      CORE::LINALG::Matrix<3, 3> temp(true);
      CORE::LINALG::Matrix<3, 3> cauchystress(true);
      temp.Multiply(1.0 / detF, defgrd, pkstress);
      cauchystress.MultiplyNT(temp, defgrd);

      // evaluate eigenproblem based on stress of previous step
      CORE::LINALG::Matrix<3, 3> lambda(true);
      CORE::LINALG::Matrix<3, 3> locsys(true);
      CORE::LINALG::SYEV(cauchystress, lambda, locsys);

      if (mat->MaterialType() == CORE::Materials::m_constraintmixture)
      {
        auto* comi = dynamic_cast<MAT::ConstraintMixture*>(mat.get());
        comi->EvaluateFiberVecs(gp, locsys, defgrd);
      }
      else if (mat->MaterialType() == CORE::Materials::m_elasthyper)
      {
        // we only have fibers at element center, thus we interpolate stress and defgrd
        avg_stress.Update(1.0 / NUMGPT_SOTET4, cauchystress, 1.0);
        avg_defgrd.Update(1.0 / NUMGPT_SOTET4, defgrd, 1.0);
      }
      else
        FOUR_C_THROW("material not implemented for remodeling");

      if (mat->MaterialType() == CORE::Materials::m_elasthyper)
      {
        // evaluate eigenproblem based on stress of previous step
        CORE::LINALG::Matrix<3, 3> lambda(true);
        CORE::LINALG::Matrix<3, 3> locsys(true);
        CORE::LINALG::SYEV(avg_stress, lambda, locsys);

        // modulation function acc. Hariton: tan g = 2nd max lambda / max lambda
        double newgamma = atan2(lambda(1, 1), lambda(2, 2));
        // compression in 2nd max direction, thus fibers are alligned to max principal direction
        if (lambda(1, 1) < 0) newgamma = 0.0;

        // new fiber vectors
        auto* elast = dynamic_cast<MAT::ElastHyper*>(mat.get());
        elast->EvaluateFiberVecs(newgamma, locsys, avg_defgrd);
      }
    }  // end loop over gauss points
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoTet4::GetCauchyNDirAndDerivativesAtXi(const CORE::LINALG::Matrix<3, 1>& xi,
    const std::vector<double>& disp, const CORE::LINALG::Matrix<3, 1>& n,
    const CORE::LINALG::Matrix<3, 1>& dir, double& cauchy_n_dir,
    CORE::LINALG::SerialDenseMatrix* d_cauchyndir_dd,
    CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd2,
    CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd_dn,
    CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd_ddir,
    CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd_dxi,
    CORE::LINALG::Matrix<3, 1>* d_cauchyndir_dn, CORE::LINALG::Matrix<3, 1>* d_cauchyndir_ddir,
    CORE::LINALG::Matrix<3, 1>* d_cauchyndir_dxi, const std::vector<double>* temp,
    CORE::LINALG::SerialDenseMatrix* d_cauchyndir_dT,
    CORE::LINALG::SerialDenseMatrix* d2_cauchyndir_dd_dT, const double* concentration,
    double* d_cauchyndir_dc)
{
  if (temp || d_cauchyndir_dT || d2_cauchyndir_dd_dT)
    FOUR_C_THROW("Thermo-elastic Nitsche contact not yet implemented in so tet4");

  cauchy_n_dir = 0.0;

  static CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xrefe(
      true);  // reference coord. of element
  static CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xcurr(
      true);  // current  coord. of element
  xrefe.Clear();
  xcurr.Clear();
  DRT::Node** nodes = Nodes();

  for (int i = 0; i < NUMNOD_SOTET4; ++i)
  {
    const auto& x = nodes[i]->X();
    for (int d = 0; d < NUMDIM_SOTET4; ++d)
    {
      xrefe(i, d) = x[d];
      xcurr(i, d) = xrefe(i, d) + disp[i * NODDOF_SOTET4 + d];
    }
  }

  static CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMNOD_SOTET4> deriv(true);
  deriv.Clear();
  CORE::FE::shape_function_deriv1<CORE::FE::CellType::tet4>(xi, deriv);

  static CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMNOD_SOTET4> N_XYZ(true);
  static CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> invJ(true);
  invJ.Multiply(1.0, deriv, xrefe, 0.0);
  invJ.Invert();
  N_XYZ.Multiply(1.0, invJ, deriv, 0.0);
  static CORE::LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> defgrd(true);
  defgrd.MultiplyTT(1.0, xcurr, N_XYZ, 0.0);

  // linearization of deformation gradient F w.r.t. displacements
  static CORE::LINALG::Matrix<9, NUMDOF_SOTET4> d_F_dd(true);
  d_F_dd.Clear();
  if (d_cauchyndir_dd || d2_cauchyndir_dd_dn || d2_cauchyndir_dd_ddir || d2_cauchyndir_dd2 ||
      d2_cauchyndir_dd_dxi)
  {
    for (int i = 0; i < NUMNOD_SOTET4; ++i)
    {
      d_F_dd(0, NODDOF_SOTET4 * i + 0) = N_XYZ(0, i);
      d_F_dd(1, NODDOF_SOTET4 * i + 1) = N_XYZ(1, i);
      d_F_dd(2, NODDOF_SOTET4 * i + 2) = N_XYZ(2, i);
      d_F_dd(3, NODDOF_SOTET4 * i + 0) = N_XYZ(1, i);
      d_F_dd(4, NODDOF_SOTET4 * i + 1) = N_XYZ(2, i);
      d_F_dd(5, NODDOF_SOTET4 * i + 0) = N_XYZ(2, i);
      d_F_dd(6, NODDOF_SOTET4 * i + 1) = N_XYZ(0, i);
      d_F_dd(7, NODDOF_SOTET4 * i + 2) = N_XYZ(1, i);
      d_F_dd(8, NODDOF_SOTET4 * i + 2) = N_XYZ(0, i);
    }
  }

  static CORE::LINALG::Matrix<9, 1> d_cauchyndir_dF(true);
  static CORE::LINALG::Matrix<9, 9> d2_cauchyndir_dF2(true);
  static CORE::LINALG::Matrix<9, NUMDIM_SOTET4> d2_cauchyndir_dF_dn(true);
  static CORE::LINALG::Matrix<9, NUMDIM_SOTET4> d2_cauchyndir_dF_ddir(true);

  SolidMaterial()->EvaluateCauchyNDirAndDerivatives(defgrd, n, dir, cauchy_n_dir, d_cauchyndir_dn,
      d_cauchyndir_ddir, &d_cauchyndir_dF, &d2_cauchyndir_dF2, &d2_cauchyndir_dF_dn,
      &d2_cauchyndir_dF_ddir, -1, Id(), concentration, nullptr, nullptr, nullptr);

  if (d_cauchyndir_dd)
  {
    d_cauchyndir_dd->reshape(NUMDOF_SOTET4, 1);
    CORE::LINALG::Matrix<NUMDOF_SOTET4, 1> d_cauchyndir_dd_mat(d_cauchyndir_dd->values(), true);
    d_cauchyndir_dd_mat.MultiplyTN(1.0, d_F_dd, d_cauchyndir_dF, 0.0);
  }

  if (d2_cauchyndir_dd_dn)
  {
    d2_cauchyndir_dd_dn->reshape(NUMDOF_SOTET4, NUMDIM_SOTET4);
    CORE::LINALG::Matrix<NUMDOF_SOTET4, NUMDIM_SOTET4> d2_cauchyndir_dd_dn_mat(
        d2_cauchyndir_dd_dn->values(), true);
    d2_cauchyndir_dd_dn_mat.MultiplyTN(1.0, d_F_dd, d2_cauchyndir_dF_dn, 0.0);
  }

  if (d2_cauchyndir_dd_ddir)
  {
    d2_cauchyndir_dd_ddir->reshape(NUMDOF_SOTET4, NUMDIM_SOTET4);
    CORE::LINALG::Matrix<NUMDOF_SOTET4, NUMDIM_SOTET4> d2_cauchyndir_dd_ddir_mat(
        d2_cauchyndir_dd_ddir->values(), true);
    d2_cauchyndir_dd_ddir_mat.MultiplyTN(1.0, d_F_dd, d2_cauchyndir_dF_ddir, 0.0);
  }

  if (d2_cauchyndir_dd2)
  {
    d2_cauchyndir_dd2->reshape(NUMDOF_SOTET4, NUMDOF_SOTET4);
    CORE::LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4> d2_cauchyndir_dd2_mat(
        d2_cauchyndir_dd2->values(), true);
    static CORE::LINALG::Matrix<9, NUMDOF_SOTET4> d2_cauchyndir_dd2_d_F_dd(true);
    d2_cauchyndir_dd2_d_F_dd.Multiply(1.0, d2_cauchyndir_dF2, d_F_dd, 0.0);
    d2_cauchyndir_dd2_mat.MultiplyTN(1.0, d_F_dd, d2_cauchyndir_dd2_d_F_dd, 0.0);
  }

  // prepare evaluation of d_cauchyndir_dxi or d2_cauchyndir_dd_dxi
  static CORE::LINALG::Matrix<CORE::FE::DisTypeToNumDeriv2<CORE::FE::CellType::tet4>::numderiv2,
      NUMNOD_SOTET4>
      deriv2(true);
  static CORE::LINALG::Matrix<9, NUMDIM_SOTET4> d_F_dxi(true);
  deriv2.Clear();
  d_F_dxi.Clear();

  if (d_cauchyndir_dxi or d2_cauchyndir_dd_dxi)
  {
    CORE::FE::shape_function_deriv2<CORE::FE::CellType::tet4>(xi, deriv2);

    static CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xXF(true);
    static CORE::LINALG::Matrix<NUMDIM_SOTET4,
        CORE::FE::DisTypeToNumDeriv2<CORE::FE::CellType::tet4>::numderiv2>
        xXFsec(true);
    xXF.Update(1.0, xcurr, 0.0);
    xXF.MultiplyNT(-1.0, xrefe, defgrd, 1.0);
    xXFsec.MultiplyTT(1.0, xXF, deriv2, 0.0);

    for (int a = 0; a < NUMDIM_SOTET4; ++a)
    {
      for (int b = 0; b < NUMDIM_SOTET4; ++b)
      {
        d_F_dxi(VoigtMapping::NonSymToVoigt9(a, b), 0) +=
            xXFsec(a, 0) * invJ(b, 0) + xXFsec(a, 3) * invJ(b, 1) + xXFsec(a, 4) * invJ(b, 2);
        d_F_dxi(VoigtMapping::NonSymToVoigt9(a, b), 1) +=
            xXFsec(a, 3) * invJ(b, 0) + xXFsec(a, 1) * invJ(b, 1) + xXFsec(a, 5) * invJ(b, 2);
        d_F_dxi(VoigtMapping::NonSymToVoigt9(a, b), 2) +=
            xXFsec(a, 4) * invJ(b, 0) + xXFsec(a, 5) * invJ(b, 1) + xXFsec(a, 2) * invJ(b, 2);
      }
    }
  }

  if (d_cauchyndir_dxi)
  {
    d_cauchyndir_dxi->MultiplyTN(1.0, d_F_dxi, d_cauchyndir_dF, 0.0);
  }

  if (d2_cauchyndir_dd_dxi)
  {
    d2_cauchyndir_dd_dxi->reshape(NUMDOF_SOTET4, NUMDIM_SOTET4);
    CORE::LINALG::Matrix<NUMDOF_SOTET4, NUMDIM_SOTET4> d2_cauchyndir_dd_dxi_mat(
        d2_cauchyndir_dd_dxi->values(), true);

    static CORE::LINALG::Matrix<CORE::FE::DisTypeToNumDeriv2<CORE::FE::CellType::tet4>::numderiv2,
        NUMDIM_SOTET4>
        Xsec(true);
    static CORE::LINALG::Matrix<NUMNOD_SOTET4,
        CORE::FE::DisTypeToNumDeriv2<CORE::FE::CellType::tet4>::numderiv2>
        N_XYZ_Xsec(true);
    Xsec.Multiply(1.0, deriv2, xrefe, 0.0);
    N_XYZ_Xsec.MultiplyTT(1.0, N_XYZ, Xsec, 0.0);

    static CORE::LINALG::Matrix<9, NUMDOF_SOTET4> d2_cauchyndir_dF2_d_F_dd(true);
    d2_cauchyndir_dF2_d_F_dd.Multiply(1.0, d2_cauchyndir_dF2, d_F_dd, 0.0);
    d2_cauchyndir_dd_dxi_mat.MultiplyTN(1.0, d2_cauchyndir_dF2_d_F_dd, d_F_dxi, 0.0);

    static CORE::LINALG::Matrix<9, NUMDIM_SOTET4 * NUMDOF_SOTET4> d2_F_dxi_dd(true);
    d2_F_dxi_dd.Clear();
    for (int i = 0; i < NUMDIM_SOTET4; ++i)
    {
      for (int j = 0; j < NUMDIM_SOTET4; ++j)
      {
        for (int k = 0; k < NUMNOD_SOTET4; ++k)
        {
          d2_F_dxi_dd(
              VoigtMapping::NonSymToVoigt9(i, j), NODDOF_SOTET4 * (NODDOF_SOTET4 * k + i) + 0) +=
              deriv2(0, k) * invJ(j, 0) + deriv2(3, k) * invJ(j, 1) + deriv2(4, k) * invJ(j, 2) -
              N_XYZ_Xsec(k, 0) * invJ(j, 0) - N_XYZ_Xsec(k, 3) * invJ(j, 1) -
              N_XYZ_Xsec(k, 4) * invJ(j, 2);

          d2_F_dxi_dd(
              VoigtMapping::NonSymToVoigt9(i, j), NODDOF_SOTET4 * (NODDOF_SOTET4 * k + i) + 1) +=
              deriv2(3, k) * invJ(j, 0) + deriv2(1, k) * invJ(j, 1) + deriv2(5, k) * invJ(j, 2) -
              N_XYZ_Xsec(k, 3) * invJ(j, 0) - N_XYZ_Xsec(k, 1) * invJ(j, 1) -
              N_XYZ_Xsec(k, 5) * invJ(j, 2);

          d2_F_dxi_dd(
              VoigtMapping::NonSymToVoigt9(i, j), NODDOF_SOTET4 * (NODDOF_SOTET4 * k + i) + 2) +=
              deriv2(4, k) * invJ(j, 0) + deriv2(5, k) * invJ(j, 1) + deriv2(2, k) * invJ(j, 2) -
              N_XYZ_Xsec(k, 4) * invJ(j, 0) - N_XYZ_Xsec(k, 5) * invJ(j, 1) -
              N_XYZ_Xsec(k, 2) * invJ(j, 2);

          for (int l = 0; l < NUMDIM_SOTET4; ++l)
          {
            d2_cauchyndir_dd_dxi_mat(k * 3 + i, l) +=
                d_cauchyndir_dF(VoigtMapping::NonSymToVoigt9(i, j), 0) *
                d2_F_dxi_dd(VoigtMapping::NonSymToVoigt9(i, j),
                    NODDOF_SOTET4 * (NODDOF_SOTET4 * k + i) + l);
          }
        }
      }
    }
  }

  if (d_cauchyndir_dc != nullptr)
  {
    static CORE::LINALG::Matrix<9, 1> d_F_dc(true);
    SolidMaterial()->EvaluateLinearizationOD(defgrd, *concentration, &d_F_dc);
    *d_cauchyndir_dc = d_cauchyndir_dF.Dot(d_F_dc);
  }
}

FOUR_C_NAMESPACE_CLOSE

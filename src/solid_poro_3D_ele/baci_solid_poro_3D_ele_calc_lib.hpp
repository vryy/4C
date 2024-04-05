/*! \file

\brief A library of free functions for a solid-poro element

\level 1
*/

#ifndef FOUR_C_SOLID_PORO_3D_ELE_CALC_LIB_HPP
#define FOUR_C_SOLID_PORO_3D_ELE_CALC_LIB_HPP


#include "baci_config.hpp"

#include "baci_discretization_fem_general_extract_values.hpp"
#include "baci_lib_element.hpp"
#include "baci_lib_element_integration_select.hpp"
#include "baci_mat_fluidporo_multiphase.hpp"
#include "baci_mat_structporo.hpp"
#include "baci_solid_3D_ele_calc_lib.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

#include <numeric>

BACI_NAMESPACE_OPEN



namespace DRT::ELEMENTS
{

  template <CORE::FE::CellType celltype>
  constexpr auto GetGaussRuleStiffnessMatrixPoro()
  {
    return DRT::ELEMENTS::DisTypeToOptGaussRule<celltype>::rule;
  }

  //! extract element data from global vector
  template <CORE::FE::CellType celltype>
  void ExtractValuesFromGlobalVector(const DRT::Discretization& discretization, const int& dofset,
      const std::vector<int>& lm,
      CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_nodes<celltype>>* matrixtofill,
      CORE::LINALG::Matrix<DETAIL::num_nodes<celltype>, 1>* vectortofill, const std::string& state,
      const DRT::Element& ele)
  {
    // get state of the global vector
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(dofset, state);
    if (matrix_state == Teuchos::null) dserror("Cannot get state vector %s", state.c_str());

    // ask for the number of dofs of dofset
    const int numdofpernode = discretization.NumDof(dofset, ele.Nodes()[0]);

    // extract local values of the global vectors
    std::vector<double> mymatrix(lm.size());
    CORE::FE::ExtractMyValues(*matrix_state, mymatrix, lm);

    if (numdofpernode == DETAIL::num_dim<celltype> + 1)
    {
      for (int inode = 0; inode < DETAIL::num_nodes<celltype>; ++inode)  // number of nodes
      {
        // fill a vector field via a pointer
        if (matrixtofill != nullptr)
        {
          for (int idim = 0; idim < DETAIL::num_dim<celltype>; ++idim)  // number of dimensions
          {
            (*matrixtofill)(idim, inode) = mymatrix[idim + (inode * numdofpernode)];
          }
        }
        // fill a scalar field via a pointer
        if (vectortofill != nullptr)
          (*vectortofill)(inode, 0) = mymatrix[DETAIL::num_dim<celltype> + (inode * numdofpernode)];
      }
    }
    else if (numdofpernode == DETAIL::num_dim<celltype>)
    {
      for (int inode = 0; inode < DETAIL::num_nodes<celltype>; ++inode)  // number of nodes
      {
        // fill a vector field via a pointer
        if (matrixtofill != nullptr)
        {
          for (int idim = 0; idim < DETAIL::num_dim<celltype>; ++idim)  // number of dimensions
          {
            (*matrixtofill)(idim, inode) = mymatrix[idim + (inode * numdofpernode)];
          }
        }
      }
    }
    else if (numdofpernode == 1)
    {
      for (std::size_t inode = 0; inode < DETAIL::num_nodes<celltype>; ++inode)  // number of nodes
      {
        if (vectortofill != nullptr) (*vectortofill)(inode, 0) = mymatrix[inode * numdofpernode];
      }
    }
    else
    {
      dserror(
          "Unknown degrees of freedom per node. Currently, only dim+1, dim and 1 are supported. "
          "You have %d dofs per node.",
          numdofpernode);
    }
  }

  /*!
   * @brief Calculate volume change
   *
   * @tparam celltype: Cell type
   * @param spatial_material_mapping (in) : An object holding quantities of the spatial material
   * mapping (deformation_gradient, inverse_deformation_gradient,
   * determinant_deformation_gradient)
   * @param jacobian_mapping (in) : An object holding quantities of the jacobian mapping
   * (inverse Jacobian, determinant, derivatives of the shape functions w.r.t. XYZ)
   * @param ele (in) : Element
   * @param discretization (in) : Discretization
   * @param lm (in) : Location vector of the element, i.e., global dof numbers of elemental dofs
   * @param kinematictype (in): kinematic type of element
   * @return volchange: volume change
   */
  template <CORE::FE::CellType celltype>
  double ComputeVolumeChange(const SpatialMaterialMapping<celltype>& spatial_material_mapping,
      const JacobianMapping<celltype>& jacobian_mapping, const DRT::Element& ele,
      const DRT::Discretization& discretization, const std::vector<int>& lm,
      const INPAR::STR::KinemType& kinematictype)
  {
    if (kinematictype == INPAR::STR::KinemType::linear)
    {
      // for linear kinematics the volume change is the trace of the linearized strains

      // gradient of displacements
      CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_nodes<celltype>> mydisp(true);
      ExtractValuesFromGlobalVector<celltype>(
          discretization, 0, lm, &mydisp, nullptr, "displacement", ele);

      CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>> dispgrad;
      dispgrad.Clear();
      // gradient of displacements
      dispgrad.MultiplyNT(mydisp, jacobian_mapping.N_XYZ_);

      double volchange = 1.0;
      // volchange = 1 + trace of the linearized strains (= trace of displacement gradient)
      for (int i = 0; i < DETAIL::num_dim<celltype>; ++i) volchange += dispgrad(i, i);

      return volchange;
    }
    else
    {
      return spatial_material_mapping.determinant_deformation_gradient_;
    }
  }

  /*!
   * @brief Calculate derivative of determinant of deformation gradient w.r.t.
   * the displacements
   *
   * @tparam celltype: Cell type
   * @param spatial_material_mapping (in) : An object holding quantities of the spatial material
   * mapping (deformation_gradient, inverse_deformation_gradient,
   * determinant_deformation_gradient)
   * @param jacobian_mapping (in) : An object holding quantities of the jacobian mapping
   * (inverse Jacobian, determinant, derivatives of the shape functions w.r.t. XYZ)
   * @param kinematictype (in): kinematic type of element
   * @return dDetDefGrad_dDisp: derivative of determinant of deformation gradient w.r.t.
   * the displacements
   */
  template <CORE::FE::CellType celltype, std::enable_if_t<DETAIL::num_dim<celltype> == 3, int> = 0>
  CORE::LINALG::Matrix<1, DETAIL::num_dof_per_ele<celltype>>
  ComputeLinearizationOfDetDefGradWrtDisp(
      const SpatialMaterialMapping<celltype> spatial_material_mapping,
      const JacobianMapping<celltype> jacobian_mapping, const INPAR::STR::KinemType& kinematictype)
  {
    CORE::LINALG::Matrix<1, DETAIL::num_dof_per_ele<celltype>> dDetDefGrad_dDisp;

    if (kinematictype == INPAR::STR::KinemType::linear)
    {
      dDetDefGrad_dDisp.Clear();
      return dDetDefGrad_dDisp;
    }
    else
    {
      CORE::LINALG::Matrix<DETAIL::num_dim<celltype> * DETAIL::num_dim<celltype>, 1> defgrd_inv_vec;
      defgrd_inv_vec(0) = spatial_material_mapping.inverse_deformation_gradient_(0, 0);
      defgrd_inv_vec(1) = spatial_material_mapping.inverse_deformation_gradient_(0, 1);
      defgrd_inv_vec(2) = spatial_material_mapping.inverse_deformation_gradient_(0, 2);
      defgrd_inv_vec(3) = spatial_material_mapping.inverse_deformation_gradient_(1, 0);
      defgrd_inv_vec(4) = spatial_material_mapping.inverse_deformation_gradient_(1, 1);
      defgrd_inv_vec(5) = spatial_material_mapping.inverse_deformation_gradient_(1, 2);
      defgrd_inv_vec(6) = spatial_material_mapping.inverse_deformation_gradient_(2, 0);
      defgrd_inv_vec(7) = spatial_material_mapping.inverse_deformation_gradient_(2, 1);
      defgrd_inv_vec(8) = spatial_material_mapping.inverse_deformation_gradient_(2, 2);

      // build N_X operator (w.r.t. material configuration)
      CORE::LINALG::Matrix<9, DETAIL::num_dof_per_ele<celltype>> N_X(true);  // set to zero
      for (int i = 0; i < DETAIL::num_nodes<celltype>; ++i)
      {
        N_X(0, 3 * i + 0) = jacobian_mapping.N_XYZ_(0, i);
        N_X(1, 3 * i + 1) = jacobian_mapping.N_XYZ_(0, i);
        N_X(2, 3 * i + 2) = jacobian_mapping.N_XYZ_(0, i);

        N_X(3, 3 * i + 0) = jacobian_mapping.N_XYZ_(1, i);
        N_X(4, 3 * i + 1) = jacobian_mapping.N_XYZ_(1, i);
        N_X(5, 3 * i + 2) = jacobian_mapping.N_XYZ_(1, i);

        N_X(6, 3 * i + 0) = jacobian_mapping.N_XYZ_(2, i);
        N_X(7, 3 * i + 1) = jacobian_mapping.N_XYZ_(2, i);
        N_X(8, 3 * i + 2) = jacobian_mapping.N_XYZ_(2, i);
      }

      // linearization of determinant of deformation gradient detF w.r.t. structural displacement
      // u_s dDetF/du_s = dDetF/dF : dF/du_s = DetF * F^-T * N,X
      dDetDefGrad_dDisp.MultiplyTN(
          spatial_material_mapping.determinant_deformation_gradient_, defgrd_inv_vec, N_X);
      return dDetDefGrad_dDisp;
    }
  }

  /*!
   * @brief Calculate derivative of volume change w.r.t. the displacements
   *
   * @tparam celltype: Cell type
   * @param dDetDefGrad_dDisp (in) : derivative of determinant of deformation gradient w.r.t.
   * the displacements
   * @param jacobian_mapping (in) : An object holding quantities of the jacobian mapping
   * (inverse Jacobian, determinant, derivatives of the shape functions w.r.t. XYZ)
   * @param kinematictype (in): kinematic type of element
   * @return dVolchange_dDisp: derivative of volume change w.r.t. the displacements
   */
  template <CORE::FE::CellType celltype, std::enable_if_t<DETAIL::num_dim<celltype> == 3, int> = 0>
  CORE::LINALG::Matrix<1, DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>
  ComputeLinearizationOfVolchangeWrtDisp(
      const CORE::LINALG::Matrix<1, DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>
          dDetDefGrad_dDisp,
      const JacobianMapping<celltype>& jacobian_mapping, const INPAR::STR::KinemType& kinematictype)
  {
    if (kinematictype == INPAR::STR::KinemType::linear)
    {
      CORE::LINALG::Matrix<1, DETAIL::num_dof_per_ele<celltype>> dVolchange_dDisp;

      for (int i = 0; i < DETAIL::num_dim<celltype>; ++i)
        for (int j = 0; j < DETAIL::num_nodes<celltype>; ++j)
          dVolchange_dDisp(DETAIL::num_dim<celltype> * j + i) = jacobian_mapping.N_XYZ_(i, j);

      return dVolchange_dDisp;
    }
    else
    {
      return dDetDefGrad_dDisp;
    }
  }

  /*!
   * @brief Calculate porosity depending on poro law given in input file and derivatives of
   * multiphase fluid primary variables w.r.t. the displacements
   *
   * @tparam celltype: Cell type
   * @param porostructmat (in) : material of skeleton (solid phase of porous domain)
   * @param params (in) : List of additional parameter to pass quantities from the time integrator
   * to the material
   * @param solidpressure (in): solid pressure
   * @param gp (in): Gauss point
   * @param volchange (in): volume change
   * @param porosity (in/out): porosity (volfrac of multiphase porspace + volfracs of additional
   *  @param dDetDefGrad_dDisp (in): derivative of determinan of deformation gradient w.r.t. the
   * displacements
   * @param dPorosity_dDisp (in/out): derivative of porosity w.r.t. the displacements
   */
  template <CORE::FE::CellType celltype>
  void ComputePorosityAndLinearization(MAT::StructPoro& porostructmat,
      Teuchos::ParameterList& params, const double solidpressure, const int gp,
      const double volchange, double& porosity,
      const CORE::LINALG::Matrix<1, DETAIL::num_dof_per_ele<celltype>>& dDetDefGrad_dDisp,
      CORE::LINALG::Matrix<1, DETAIL::num_dof_per_ele<celltype>>& dPorosity_dDisp)
  {
    double dphi_dJ = 0.0;

    porostructmat.ComputePorosity(params, solidpressure, volchange, gp, porosity,
        nullptr,  // dphi_dp not needed
        &dphi_dJ,
        nullptr,  // dphi_dJdp not needed
        nullptr,  // dphi_dJJ not needed
        nullptr   // dphi_dpp not needed
    );

    dPorosity_dDisp.Update(dphi_dJ, dDetDefGrad_dDisp);
  }

  /*!
   * @brief Calculate porosity depending on poro law given in input file
   *
   * @tparam celltype: Cell type
   * @param porostructmat (in) : material of skeleton (solid phase of porous domain)
   * @param params (in) : List of additional parameter to pass quantities from the time integrator
   * to the material
   * @param solidpressure (in): solid pressure
   * @param volchange (in): volume change
   * @param gp (in): Gauss point
   * @return porosity (volfrac of multiphase porspace + volfracs of additional porous networks)
   */
  template <CORE::FE::CellType celltype>
  double ComputePorosity(MAT::StructPoro& porostructmat, Teuchos::ParameterList& params,
      const double solidpressure, const double volchange, const int gp)
  {
    double porosity = 0.0;
    porostructmat.ComputePorosity(params, solidpressure, volchange, gp, porosity, nullptr,
        nullptr,  // dphi_dJ not needed
        nullptr,  // dphi_dJdp not needed
        nullptr,  // dphi_dJJ not needed
        nullptr   // dphi_dpp not needed
    );
    return porosity;
  }

  /*!
   * @brief Recalculate derivative of solidpressure w.r.t. the displacements in case of volfracs
   *
   * @tparam celltype: Cell type
   * @param fluidpress (in) : solid pressure contribution coming from the multiphase fluid S_i*p_i
   * @param porosity (in) : porosity = volumefraction in multiphase porespace + volfracs from
   * additional porous networks
   * @param nummultifluiddofpernode (in): number of fluid multiphase dofs per node
   * @param numfluidphases (in): number of fluidphases in multiphase porespace
   * @param numvolfrac (in): number of volfracs
   * @param fluidmultiphase_phiAtGP (in): fluid multiphase primary variables at GP
   * @param dPorosity_dDisp (in): derivative of porosity w.r.t. the displacements
   * @param dSolidpressure_dDisp (in/out): derivative of solidpressure w.r.t. the displacements
   */
  template <CORE::FE::CellType celltype>
  void RecalculateLinearizationOfSolPressWrtDisp(const double fluidpress, const double porosity,
      const int nummultifluiddofpernode, const int numfluidphases, const int numvolfrac,
      const std::vector<double>& fluidmultiphase_phiAtGP,
      const CORE::LINALG::Matrix<1, DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>&
          dPorosity_dDisp,
      CORE::LINALG::Matrix<1, DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>&
          dSolidpressure_dDisp)
  {
    // get volume fraction primary variables
    std::vector<double> volfracphi(&fluidmultiphase_phiAtGP[numfluidphases],
        &fluidmultiphase_phiAtGP[numfluidphases + numvolfrac]);
    double sumaddvolfrac = 0.0;
    for (int ivolfrac = 0; ivolfrac < numvolfrac; ivolfrac++) sumaddvolfrac += volfracphi[ivolfrac];

    // get volume fraction pressure at [numfluidphases+numvolfrac...nummultifluiddofpernode-1]
    std::vector<double> volfracpressure(&fluidmultiphase_phiAtGP[numfluidphases + numvolfrac],
        &fluidmultiphase_phiAtGP[nummultifluiddofpernode]);

    // p_s = (porosity - sumaddvolfrac)/porosity * fluidpress
    //       + 1.0 / porosity sum_i=1^numvolfrac (volfrac_i*pressure_i)
    // d (p_s) / d porosity = + sumaddvolfrac/porosity/porosity * fluidpress
    double dps_dphi = sumaddvolfrac / (porosity * porosity) * fluidpress;

    // ... + 1.0 / porosity / porosity sum_i=1^numvolfrac (volfrac_i*pressure_i)
    for (int ivolfrac = 0; ivolfrac < numvolfrac; ivolfrac++)
      dps_dphi -= volfracphi[ivolfrac] * volfracpressure[ivolfrac] / (porosity * porosity);

    // d (p_s) / d u_s = d (p_s) / d porosity * d porosity / d u_s
    dSolidpressure_dDisp.Update(dps_dphi, dPorosity_dDisp);
  }

  /*!
   * @brief Recalculate solidpressure in case of volfracs
   *
   * @param press (in) : current solid pressure
   * @param porosity (in) : porosity = volumefraction in multiphase porespace + volfracs from
   * additional porous networks
   * @param nummultifluiddofpernode (in): number of fluid multiphase dofs per node
   * @param numfluidphases (in): number of fluidphases in multiphase porespace
   * @param numvolfrac (in): number of volfracs
   * @param fluidmultiphase_phiAtGP (in): fluid multiphase primary variables at GP
   * @return solid pressure
   */
  double RecalculateSolPressureAtGP(double press, const double porosity,
      const int nummultifluiddofpernode, const int numfluidphases, const int numvolfrac,
      const std::vector<double>& fluidmultiphase_phiAtGP)
  {
    // get volume fraction primary variables at [numfluidphases-1...numfluidphase-1+numvolfrac]
    std::vector<double> volfracphi(&fluidmultiphase_phiAtGP[numfluidphases],
        &fluidmultiphase_phiAtGP[numfluidphases + numvolfrac]);
    double sumaddvolfrac = 0.0;
    for (int ivolfrac = 0; ivolfrac < numvolfrac; ivolfrac++) sumaddvolfrac += volfracphi[ivolfrac];

    // p_s = (porosity - sumaddvolfrac)/porosity * fluidpress
    //      + 1.0 / porosity * sum_i=1^numvolfrac (volfrac_i*pressure_i)
    // first part
    press *= (porosity - sumaddvolfrac) / porosity;

    // get volfrac pressures at [numfluidphases+numvolfrac...nummultifluiddofpernode-1]
    std::vector<double> volfracpressure(&fluidmultiphase_phiAtGP[numfluidphases + numvolfrac],
        &fluidmultiphase_phiAtGP[nummultifluiddofpernode]);

    // second part
    for (int ivolfrac = 0; ivolfrac < numvolfrac; ivolfrac++)
      press += volfracphi[ivolfrac] / porosity * volfracpressure[ivolfrac];

    // note: in RecalculateSolidPressure in porofluid_phasemanager calculation is performed a bit
    //       differently since we already pass porosity = porosity - sumaddvolfrac, but result is
    //       equivalent

    return press;
  }

  /*!
   * @brief Update the internal force vector with poroelasticity contribution of one Gauss point
   *
   * @tparam celltype: Cell type
   * @param detJ_w (in) : integration factor (Gauss point weight times the determinant of
   * the jacobian)
   * @param solidpressure (in) : solid pressure
   * @param DetDefGrad (in) : Determinant of deformation gradient
   * * @param BopCinv (in) : B^T . C^-1
   * @param force_vector (in/out) : Force vector where the local contribution is added to
   */
  template <CORE::FE::CellType celltype>
  void UpdateInternalForceVectorMultiPhasePressureBased(const double detJ_w,
      const double solidpressure, const double DetDefGrad,
      const CORE::LINALG::Matrix<DETAIL::num_dof_per_ele<celltype>, 1>& BopCinv,
      CORE::LINALG::Matrix<DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>, 1>&
          force_vector)
  {
    // additional fluid stress- stiffness term RHS -(B^T .  C^-1  * J * p^f * detJ * w(gp))
    double factor = -detJ_w * solidpressure * DetDefGrad;
    force_vector.Update(factor, BopCinv, 1.0);
  }

  /*!
   * @brief Update elastic stiffness matrix with poroelasticity contribution of one Gauss point
   *
   * @tparam celltype: Cell type
   * @param detJ_w (in) : integration factor (Gauss point weight times the determinant of
   * the jacobian)
   * @param solidpressure (in) : solidpressure
   * @param DetDefGrad (in) : Determinant of deformation gradient
   * @param BopCinv (in) : B^T . C^-1
   * @param Bop (in) : Strain gradient (B-Operator)
   * @param dDetDefGrad_dDisp (in) : derivative of determinant of derformation gradient w.r.t. the
   * displacements
   * @param dSolidpressure_dDisp (in) : derivative of solidpressure w.r.t. the
   * displacements
   * @param dInverseRightCauchyGreen_dDisp (in) : derivatives of inverse right cauchy green tensor
   * w.r.t. the displacements
   * @param stiffness_matrix (in/out) : stiffness matrix where the local contribution is added to
   */
  template <CORE::FE::CellType celltype>
  void UpdateElasticStiffnessMatrixMultiPhasePressureBased(const double detJ_w,
      const double solidpressure, const double DetDefGrad,
      const CORE::LINALG::Matrix<DETAIL::num_dof_per_ele<celltype>, 1>& BopCinv,
      const CORE::LINALG::Matrix<DETAIL::num_str<celltype>, DETAIL::num_dof_per_ele<celltype>>& Bop,
      const CORE::LINALG::Matrix<1, DETAIL::num_dof_per_ele<celltype>>& dDetDefGrad_dDisp,
      const CORE::LINALG::Matrix<1, DETAIL::num_dof_per_ele<celltype>>& dSolidpressure_dDisp,
      const CORE::LINALG::Matrix<DETAIL::num_str<celltype>, DETAIL::num_dof_per_ele<celltype>>&
          dInverseRightCauchyGreen_dDisp,
      CORE::LINALG::Matrix<DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>,
          DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>& stiffness_matrix)
  {
    CORE::LINALG::Matrix<DETAIL::num_dof_per_ele<celltype>, DETAIL::num_dof_per_ele<celltype>> tmp;

    // additional fluid stress- stiffness term -(B^T . C^-1 . dJ/d(us) * p^f * detJ * w(gp))
    tmp.Multiply((-detJ_w * solidpressure), BopCinv, dDetDefGrad_dDisp);
    stiffness_matrix.Update(1.0, tmp, 1.0);

    // additional fluid stress- stiffness term -(B^T .  dC^-1/d(us) * J * p^f * detJ * w(gp))
    tmp.MultiplyTN((-detJ_w * solidpressure * DetDefGrad), Bop, dInverseRightCauchyGreen_dDisp);
    stiffness_matrix.Update(1.0, tmp, 1.0);

    // additional fluid stress- stiffness term -(B^T .  dC^-1 * J * dp^s/d(us) * detJ * w(gp))
    tmp.Multiply((-detJ_w * DetDefGrad), BopCinv, dSolidpressure_dDisp);
    stiffness_matrix.Update(1.0, tmp, 1.0);
  }

  /*!
   * @brief Update geometric stiffness matrix with poroelasticity contribution of one Gauss point
   *
   * @tparam celltype: Cell type
   * @param detJ_w (in) : integration factor (Gauss point weight times the determinant of
   * the jacobian)
   * @param solidpressure (in) : solidpressure
   * @param DetDefGrad (in) : Determinant of deformation gradient
   * @param C_inv_vec (in) : inverse Right Cauchy-Green tensor as vector in voigt notation
   * @param n_xyz (in) : derivatives of the shape functions w.r.t. XYZ
   * @param stiffness_matrix (in/out) : stiffness matrix where the local contribution is added to
   */
  template <CORE::FE::CellType celltype>
  void UpdateGeometricStiffnessMatrixMultiPhasePressureBased(const double detJ_w,
      const double solidpressure, const double DetDefGrad,
      const CORE::LINALG::Matrix<DETAIL::num_str<celltype>, 1>& C_inv_vec,
      const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_nodes<celltype>>& N_XYZ_,
      CORE::LINALG::Matrix<DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>,
          DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>& stiffness_matrix)
  {
    // integrate `geometric' stiffness matrix and add to keu *****************
    CORE::LINALG::Matrix<DETAIL::num_str<celltype>, 1> sfac(
        C_inv_vec);  // auxiliary integrated stress

    // scale
    sfac.Scale((-detJ_w * solidpressure * DetDefGrad));

    std::vector<double> SmB_L(3);  // intermediate Sm.B_L
    // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
    for (int inod = 0; inod < DETAIL::num_nodes<celltype>; ++inod)
    {
      SmB_L[0] = sfac(0) * N_XYZ_(0, inod) + sfac(3) * N_XYZ_(1, inod) + sfac(5) * N_XYZ_(2, inod);
      SmB_L[1] = sfac(3) * N_XYZ_(0, inod) + sfac(1) * N_XYZ_(1, inod) + sfac(4) * N_XYZ_(2, inod);
      SmB_L[2] = sfac(5) * N_XYZ_(0, inod) + sfac(4) * N_XYZ_(1, inod) + sfac(2) * N_XYZ_(2, inod);
      for (int jnod = 0; jnod < DETAIL::num_nodes<celltype>; ++jnod)
      {
        double bopstrbop = 0.0;  // intermediate value
        for (int idim = 0; idim < DETAIL::num_dim<celltype>; ++idim)
          bopstrbop += N_XYZ_(idim, jnod) * SmB_L[idim];
        stiffness_matrix(DETAIL::num_dim<celltype> * inod + 0,
            DETAIL::num_dim<celltype> * jnod + 0) += bopstrbop;
        stiffness_matrix(DETAIL::num_dim<celltype> * inod + 1,
            DETAIL::num_dim<celltype> * jnod + 1) += bopstrbop;
        stiffness_matrix(DETAIL::num_dim<celltype> * inod + 2,
            DETAIL::num_dim<celltype> * jnod + 2) += bopstrbop;
      }
    }
  }

  /*!
   * @brief Add coupling contribution (poroelasticity OD entries) of one Gauss point to stiffness
   * matrix
   *
   * @tparam celltype: Cell type
   * @param detJ_w (in) : integration factor (Gauss point weight times the determinant of
   * the jacobian)
   * @param solidpressurederiv (in) : derivative of solidpressure w.r.t. fluid multiphase
   * @param BopCinv (in) : B^T * C^-1
   * @param shape_functions (in) : Shape function
   * * @param DetDefGrad (in) : Determinant of deformation gradient
   * @param nummultifluiddofpernode (in) : number of fluid multiphase dofs per node
   * @param stiffness_matrix (in/out) : stiffness matrix where the local contribution is added to
   */
  template <CORE::FE::CellType celltype>
  void UpdateStiffnessMatrixCouplingMultiPhasePressureBased(const double detJ_w,
      const std::vector<double>& solidpressurederiv,
      const CORE::LINALG::Matrix<DETAIL::num_dof_per_ele<celltype>, 1>& BopCinv,
      ShapeFunctionsAndDerivatives<celltype> shape_functions, const double DetDefGrad,
      const int nummultifluiddofpernode, CORE::LINALG::SerialDenseMatrix& stiffness_matrix)
  {
    for (int i = 0; i < DETAIL::num_nodes<celltype>; i++)
    {
      const int fi = DETAIL::num_dim<celltype> * i;

      for (int j = 0; j < DETAIL::num_dim<celltype>; j++)
      {
        for (int k = 0; k < DETAIL::num_nodes<celltype>; k++)
        {
          for (int iphase = 0; iphase < nummultifluiddofpernode; iphase++)
          {
            int fk_press = k * nummultifluiddofpernode + iphase;

            /*-------structure- fluid pressure coupling: "stress term"
             -B^T . ( -1*J*C^-1 ) * Dp
             */
            stiffness_matrix(fi + j, fk_press) += detJ_w * BopCinv(fi + j) * (-1.0) * DetDefGrad *
                                                  shape_functions.shapefunctions_(k) *
                                                  solidpressurederiv[iphase];
          }
        }
      }
    }
  }

  /*!
   * @brief Calculates fluid mulltiphase primary variables at GP
   *
   * @tparam celltype: Cell type
   * @param fluidmultiphase_ephi (in) : primary variables of multiphase porous medium flow
   * @param nummultifluiddofpernode (in) : number of fluid multiphase dofs per node
   * @param: shape_functions (in): Shape functions
   * @returns: fluidmultiphase_phiAtGP: fluid multiphase primary variables at GP
   */
  template <CORE::FE::CellType celltype>
  std::vector<double> ComputeFluidMultiPhasePrimaryVariablesAtGP(
      const std::vector<double>& fluidmultiphase_ephi, const int nummultifluiddofpernode,
      const ShapeFunctionsAndDerivatives<celltype>& shape_functions)
  {
    std::vector<double> fluidmultiphase_phiAtGP(nummultifluiddofpernode);
    // compute phi at GP = phi * shapefunction
    for (int i = 0; i < DETAIL::num_nodes<celltype>; i++)
    {
      for (int j = 0; j < nummultifluiddofpernode; j++)
      {
        fluidmultiphase_phiAtGP[j] += shape_functions.shapefunctions_(i) *
                                      fluidmultiphase_ephi[i * nummultifluiddofpernode + j];
      }
    }

    return fluidmultiphase_phiAtGP;
  }

  /*!
   * @brief Calculates solid pressure derivatives w.r.t. primary variables of fluid phases in
   * multiphase porespace
   *
   * @tparam celltype: Cell type
   * @param porofluidmat (in) : material of multiphase fluid
   * @param fluidmultiphase_phiAtGP (in) : fluid multiphase primary variables at GP
   * @param numfluidphases (in): number of fluidphases in multiphase porespace
   * @returns solidpressurederiv: derivative of solidpressure w.r.t. fluid multiphase
   * primary variables
   */
  template <CORE::FE::CellType celltype>
  std::vector<double> ComputeSolidPressureDeriv(MAT::FluidPoroMultiPhase& porofluidmat,
      const std::vector<double>& fluidmultiphase_phiAtGP, const int numfluidphases)
  {
    // zero out everything
    std::vector<double> solidpressurederiv(fluidmultiphase_phiAtGP.size());

    // initialize auxiliary variables
    std::vector<double> genpress(numfluidphases);
    std::vector<double> press(numfluidphases);
    std::vector<double> sat(numfluidphases);
    CORE::LINALG::SerialDenseMatrix helpderiv(numfluidphases, numfluidphases, true);
    CORE::LINALG::SerialDenseMatrix satderiv(numfluidphases, numfluidphases, true);
    CORE::LINALG::SerialDenseMatrix pressderiv(numfluidphases, numfluidphases, true);
    std::vector<double> fluidphi(
        &fluidmultiphase_phiAtGP[0], &fluidmultiphase_phiAtGP[numfluidphases]);

    // evaluate the pressures
    porofluidmat.EvaluateGenPressure(genpress, fluidphi);

    // transform generalized pressures to true pressure values
    porofluidmat.TransformGenPresToTruePres(genpress, press);

    // explicit evaluation of saturation
    porofluidmat.EvaluateSaturation(sat, fluidphi, press);

    // calculate the derivative of the pressure (actually first its inverse)
    porofluidmat.EvaluateDerivOfDofWrtPressure(pressderiv, fluidphi);

    // now invert the derivatives of the dofs w.r.t. pressure to get the derivatives
    // of the pressure w.r.t. the dofs
    {
      Teuchos::SerialDenseSolver<int, double> inverse;

      inverse.setMatrix(Teuchos::rcpFromRef(pressderiv));
      int err = inverse.invert();
      if (err != 0)
        dserror("Inversion of matrix for pressure derivative failed with error code %d.", err);
    }

    // calculate derivatives of saturation w.r.t. pressure
    porofluidmat.EvaluateDerivOfSaturationWrtPressure(helpderiv, press);

    // chain rule: the derivative of saturation w.r.t. dof =
    // (derivative of saturation w.r.t. pressure) * (derivative of pressure w.r.t. dof)
    satderiv.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, helpderiv, pressderiv, 0.0);

    // compute derivative of solid pressure w.r.t. dofs with product rule
    // standard derivative: no volume fractions present
    for (int iphase = 0; iphase < numfluidphases; iphase++)
    {
      for (int jphase = 0; jphase < numfluidphases; jphase++)
        solidpressurederiv[iphase] +=
            pressderiv(jphase, iphase) * sat[jphase] + satderiv(jphase, iphase) * press[jphase];
    }

    return solidpressurederiv;
  }

  /*!
   * @brief Calculates solidpressure
   *
   * @tparam celltype: Cell type
   * @param nummultifluiddofpernode (in) : number of fluid multiphase dofs per node
   * @param numfluidphases (in) : number of fluidphases in multiphase porespace
   * @param: fluidmultiphase_phiAtGP (in): fluid multiphase primary variables at GP
   * @param: porofluidmat (in): material of multiphase fluid
   * @return solidpressure
   */
  template <CORE::FE::CellType celltype>
  double ComputeSolPressureAtGP(const int nummultifluiddofpernode, const int numfluidphases,
      const std::vector<double>& fluidmultiphase_phiAtGP, MAT::FluidPoroMultiPhase& porofluidmat)
  {
    // initialize auxiliary variables
    std::vector<double> genpress(numfluidphases, 0.0);
    std::vector<double> sat(numfluidphases, 0.0);
    std::vector<double> press(numfluidphases, 0.0);
    std::vector<double> fluidphi(
        &fluidmultiphase_phiAtGP[0], &fluidmultiphase_phiAtGP[numfluidphases]);

    // evaluate the pressures
    porofluidmat.EvaluateGenPressure(genpress, fluidphi);

    //! transform generalized pressures to true pressure values
    porofluidmat.TransformGenPresToTruePres(genpress, press);

    // explicit evaluation of saturation
    porofluidmat.EvaluateSaturation(sat, fluidphi, press);

    // solid pressure = sum (S_i*p_i)
    const double solidpressure = std::inner_product(sat.begin(), sat.end(), press.begin(), 0.0);

    return solidpressure;
  }

  /*!
   * @brief Recalculates solid pressure derivative in case of volfracs
   *
   * @param fluidmultiphase_phiAtGP (in) : fluid multiphase primary variables at GP
   * @param nummultifluiddofpernode (in) : number of fluid multiphase dofs per node
   * @param: numfluidphases (in): number of fluidphases in multiphase porespace
   * @param: numvolfrac (in): number of volfracs from additionalporous networks
   * @param: solidpressure (in): solidpressue
   * @param: porosity (in): porosity = volumefraction in multiphase porespace + volfracs from
   * additional porous networks
   * @param: solidpressurederiv (in/out): derivative of solidpressure w.r.t. fluid multiphase
   * primary variables and volfracs
   */
  void RecalculateSolPressureDeriv(const std::vector<double>& fluidmultiphase_phiAtGP,
      const int nummultifluiddofpernode, const int numfluidphases, const int numvolfrac,
      const double solidpressure, const double porosity, std::vector<double>& solidpressurederiv)
  {
    // get volume fraction primary variables
    std::vector<double> volfracphi(&fluidmultiphase_phiAtGP[numfluidphases],
        &fluidmultiphase_phiAtGP[numfluidphases + numvolfrac]);
    double sumaddvolfrac = 0.0;
    for (int ivolfrac = 0; ivolfrac < numvolfrac; ivolfrac++) sumaddvolfrac += volfracphi[ivolfrac];

    // p_s = (porosity - sumaddvolfrac)/porosity * fluidpress
    //      + 1.0 / porosity sum_i=1^numvolfrac (volfrac_i*pressure_i)
    const double scale = (porosity - sumaddvolfrac) / porosity;

    // scale original fluid press deriv with (porosity - sumaddvolfrac)/porosity
    for (int iphase = 0; iphase < numfluidphases; iphase++) solidpressurederiv[iphase] *= scale;

    // get volfrac pressures at [numfluidphases+numvolfrac...nummultifluiddofpernode-1]
    std::vector<double> volfracpressure(&fluidmultiphase_phiAtGP[numfluidphases + numvolfrac],
        &fluidmultiphase_phiAtGP[nummultifluiddofpernode]);

    for (int ivolfrac = 0; ivolfrac < numvolfrac; ivolfrac++)
    {
      // d p_s / d volfrac = - fluidpress/porosity + volfracpressure/porosity
      solidpressurederiv[ivolfrac + numfluidphases] =
          -1.0 / porosity * solidpressure + 1.0 / porosity * volfracpressure[ivolfrac];
      // d p_s / d volfracpress = + volfracphi/porosity
      solidpressurederiv[ivolfrac + numfluidphases + numvolfrac] = volfracphi[ivolfrac] / porosity;
    }
  }

  template <CORE::FE::CellType celltype>
  struct CauchyGreenAndInverse
  {
    CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>> right_cauchy_green_;
    CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>
        inverse_right_cauchy_green_;
  };

  /*!
   * @brief Evaluates right Cauchy-Green deformation tensor and its inverse
   *
   * @tparam celltype: Cell type
   * @param spatial_material_mapping (in) : An object holding quantities of the spatial material
   * mapping (deformation_gradient, inverse_deformation_gradient,
   * determinant_deformation_gradient)
   * @return CauchyGreenAndInverse<celltype> : An object holding the right Cauchy-Green deformation
   * tensor and its inverse
   */
  template <CORE::FE::CellType celltype, std::enable_if_t<DETAIL::num_dim<celltype> == 3, int> = 0>
  CauchyGreenAndInverse<celltype> EvaluateCauchyGreenAndInverse(
      const SpatialMaterialMapping<celltype>& spatial_material_mapping)
  {
    CauchyGreenAndInverse<celltype> cauchygreen;

    cauchygreen.right_cauchy_green_ = DRT::ELEMENTS::EvaluateCauchyGreen(spatial_material_mapping);
    cauchygreen.inverse_right_cauchy_green_.Invert(cauchygreen.right_cauchy_green_);

    return cauchygreen;
  }

  /*!
   * @brief Evaluates the derivative of the inverse right Cauchy-Green deformation tensor w.r.t. the
   * displacements
   *
   * @tparam celltype: Cell type
   * @param cauchygreen (in) : An object holding the right Cauchy-Green deformation tensor and its
   * inverse
   * @param jacobian_mapping (in) : n object holding quantities of the jacobian mapping
   * (inverse Jacobian, determinant, derivatives of the shape functions w.r.t. XYZ)
   * @param: spatial_material_mapping (in): An object holding quantities of the spatial material
   * mapping (deformation_gradient, inverse_deformation_gradient,
   * determinant_deformation_gradient)
   * @return dInverseCauchyGreen_dDisp : derivative of the inverse right Cauchy-Green deformation
   * tensor w.r.t. the displacements
   */
  template <CORE::FE::CellType celltype, std::enable_if_t<DETAIL::num_dim<celltype> == 3, int> = 0>
  CORE::LINALG::Matrix<DETAIL::num_str<celltype>,
      DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>
  EvaluateInverseCauchyGreenLinearization(CauchyGreenAndInverse<celltype> cauchygreen,
      JacobianMapping<celltype> jacobian_mapping,
      SpatialMaterialMapping<celltype> spatial_material_mapping)
  {
    // dC^-1/dDisp
    CORE::LINALG::Matrix<DETAIL::num_str<celltype>,
        DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>
        dInverseCauchyGreen_dDisp(true);

    for (int n = 0; n < DETAIL::num_nodes<celltype>; ++n)
    {
      for (int k = 0; k < DETAIL::num_dim<celltype>; ++k)
      {
        const int gid = n * DETAIL::num_dim<celltype> + k;
        for (int i = 0; i < DETAIL::num_dim<celltype>; ++i)
        {
          dInverseCauchyGreen_dDisp(0, gid) +=
              -2 * cauchygreen.inverse_right_cauchy_green_(0, i) * jacobian_mapping.N_XYZ_(i, n) *
              spatial_material_mapping.inverse_deformation_gradient_(0, k);
          dInverseCauchyGreen_dDisp(1, gid) +=
              -2 * cauchygreen.inverse_right_cauchy_green_(1, i) * jacobian_mapping.N_XYZ_(i, n) *
              spatial_material_mapping.inverse_deformation_gradient_(1, k);
          dInverseCauchyGreen_dDisp(2, gid) +=
              -2 * cauchygreen.inverse_right_cauchy_green_(2, i) * jacobian_mapping.N_XYZ_(i, n) *
              spatial_material_mapping.inverse_deformation_gradient_(2, k);
          /* ~~~ */
          dInverseCauchyGreen_dDisp(3, gid) +=
              -cauchygreen.inverse_right_cauchy_green_(0, i) * jacobian_mapping.N_XYZ_(i, n) *
                  spatial_material_mapping.inverse_deformation_gradient_(1, k) -
              spatial_material_mapping.inverse_deformation_gradient_(0, k) *
                  jacobian_mapping.N_XYZ_(i, n) * cauchygreen.inverse_right_cauchy_green_(1, i);
          dInverseCauchyGreen_dDisp(4, gid) +=
              -cauchygreen.inverse_right_cauchy_green_(1, i) * jacobian_mapping.N_XYZ_(i, n) *
                  spatial_material_mapping.inverse_deformation_gradient_(2, k) -
              spatial_material_mapping.inverse_deformation_gradient_(1, k) *
                  jacobian_mapping.N_XYZ_(i, n) * cauchygreen.inverse_right_cauchy_green_(2, i);
          dInverseCauchyGreen_dDisp(5, gid) +=
              -cauchygreen.inverse_right_cauchy_green_(2, i) * jacobian_mapping.N_XYZ_(i, n) *
                  spatial_material_mapping.inverse_deformation_gradient_(0, k) -
              spatial_material_mapping.inverse_deformation_gradient_(2, k) *
                  jacobian_mapping.N_XYZ_(i, n) * cauchygreen.inverse_right_cauchy_green_(0, i);
        }
      }
    }
    return dInverseCauchyGreen_dDisp;
  }
}  // namespace DRT::ELEMENTS

BACI_NAMESPACE_CLOSE

#endif
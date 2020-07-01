/*----------------------------------------------------------------------------*/
/*! \file
\brief (augmented) contact integration policies

\level 3

*/
/*----------------------------------------------------------------------------*/

//#define CONTACT_AUG_JACOBIAN_DEBUG_OUTPUT

#include "contact_aug_integrator_policy.H"

#include "../drt_mortar/mortar_element.H"
#include "../drt_contact/contact_node.H"
#include "../drt_io/io_pstream.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::IncompleteIntPolicy<probdim, slavetype, mastertype>::Get_Deriv2nd_Jacobian(
    const MORTAR::MortarElement& sele,
    const LINALG::Matrix<probdim, my::SLAVENUMNODE, int>& nodal_dofs,
    const LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& deriv,
    const LINALG::Matrix<probdim, 1>& unit_normal, const double length_n_inv,
    const Deriv1stVecMap& d_non_unit_normal, Deriv2ndMap& dd_jac) const
{
  // do nothing
  GEN_DATA::reset(0, dd_jac);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::CompleteIntPolicy<probdim, slavetype, mastertype>::Get_Deriv2nd_Jacobian(
    const MORTAR::MortarElement& sele,
    const LINALG::Matrix<probdim, my::SLAVENUMNODE, int>& nodal_dofs,
    const LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& deriv,
    const LINALG::Matrix<probdim, 1>& unit_normal, const double length_n_inv,
    const Deriv1stVecMap& d_non_unit_normal, Deriv2ndMap& dd_jac) const
{
  /*----------------------------------------------------------------------*/
  // non-unit normal vector: 2-nd order derivative
  // 1-st int: vector index corresponds to the normal component index
  // 2-nd int: 1-st paired vector key corresponds to varied dof GID
  // 3-rd int: 2-nd paired vector key corresponds to linearized dof GID
  Deriv2ndVecMap dd_non_unit_normal(probdim);
  this->Deriv2nd_NonUnitSlaveElementNormal(sele, nodal_dofs, deriv, dd_non_unit_normal);

  /*----------------------------------------------------------------------*/
  // unit normal vector: 1-st order derivative
  // 1-st int: vector index corresponds to the normal component index
  // 2-nd int: paired vector key corresponds to varied dof GID
  Deriv1stVecMap d_unit_normal(probdim);
  this->Deriv1st_UnitSlaveElementNormal(
      unit_normal, length_n_inv, d_non_unit_normal, d_unit_normal);

  /*----------------------------------------------------------------------*/
  // jacobian determinant: 2-nd order derivative
  // 1-st int: 1-st paired vector key corresponds to varied dof GID
  // 2-nd int: 2-nd paired vector key corresponds to linearized dof GID
  this->Deriv2nd_Jacobian(
      d_unit_normal, d_non_unit_normal, unit_normal, dd_non_unit_normal, dd_jac);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype>
void CONTACT::AUG::BaseSlaveIntPolicy<probdim, slavetype>::Deriv1st_Jacobian(
    const LINALG::Matrix<probdim, 1>& unit_normal, const Deriv1stVecMap& d_non_unit_normal,
    Deriv1stMap& d_jac) const
{
  GEN_DATA::reset(probdim * SLAVENUMNODE, d_jac);

  for (unsigned n_dof = 0; n_dof < probdim; ++n_dof)
  {
    for (auto& d_nun_pair : d_non_unit_normal[n_dof])
    {
      d_jac[d_nun_pair.first] += d_nun_pair.second * unit_normal(n_dof, 0);
    }
  }

#ifdef CONTACT_AUG_JACOBIAN_DEBUG_OUTPUT
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  for (auto d_jac_pair : d_jac)
    std::cout << "varied GID (" << d_jac_pair.first << ") = " << d_jac_pair.second << std::endl;
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype>
void CONTACT::AUG::BaseSlaveIntPolicy<probdim, slavetype>::Deriv2nd_Jacobian(
    const Deriv1stVecMap& d_unit_normal, const Deriv1stVecMap& d_non_unit_normal,
    const LINALG::Matrix<probdim, 1>& unit_normal, const Deriv2ndVecMap& dd_non_unit_normal,
    Deriv2ndMap& dd_jac) const
{
  this->timer_.start(TimeID::Deriv2nd_Jacobian);

#ifdef CONTACT_AUG_JACOBIAN_DEBUG_OUTPUT
  std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif
  GEN_DATA::reset(probdim * SLAVENUMNODE, dd_jac);

  // (1) Inner product of varied non-unit normal and linearized unit-normal
  for (unsigned n_dof = 0; n_dof < probdim; ++n_dof)
  {
    for (auto& d_nun_pair : d_non_unit_normal[n_dof])
    {
      Deriv1stMap& dd_jac_var_dof = dd_jac[d_nun_pair.first];

      for (auto& d_un_pair : d_unit_normal[n_dof])
      {
        dd_jac_var_dof[d_un_pair.first] += d_nun_pair.second * d_un_pair.second;
      }
    }
  }

  // (2) Inner product of linearized varied non-unit normal and the unit-normal
  for (unsigned n_dof = 0; n_dof < probdim; ++n_dof)
  {
    for (auto& dd_nun_pair_var : dd_non_unit_normal[n_dof])
    {
      Deriv1stMap& dd_jac_var_dof = dd_jac[dd_nun_pair_var.first];
      for (auto& dd_nun_pair_lin : dd_nun_pair_var.second)
      {
        dd_jac_var_dof[dd_nun_pair_lin.first] += unit_normal(n_dof, 0) * dd_nun_pair_lin.second;
      }
    }
  }

#ifdef CONTACT_AUG_JACOBIAN_DEBUG_OUTPUT
  for (auto& dd_jac_pair_var : dd_jac)
  {
    for (auto& dd_jac_pair_lin : dd_jac_pair_var.second)
    {
      std::cout << "varied GID (" << dd_jac_pair_var.first << "), "
                << "lin GID (" << dd_jac_pair_lin.first << ") = " << dd_jac_pair_lin.second << "\n";
    }
  }
#endif

  this->timer_.stop(TimeID::Deriv2nd_Jacobian);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype>
void CONTACT::AUG::BaseSlaveIntPolicy<probdim, slavetype>::Deriv1st_UnitSlaveElementNormal(
    const LINALG::Matrix<probdim, 1>& unit_normal, const double length_n_inv,
    const Deriv1stVecMap& d_non_unit_normal, Deriv1stVecMap& d_unit_normal, const bool reset) const
{
#ifdef CONTACT_AUG_JACOBIAN_DEBUG_OUTPUT
  std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

  if (reset) GEN_DATA::reset(probdim, SLAVENUMNODE * probdim, d_unit_normal);
  // safety check
  else
  {
    if (d_unit_normal.size() < probdim) dserror("Given vector has the wrong size.");
    for (auto& d_un_ndof : d_unit_normal)
      if (d_un_ndof.capacity() < SLAVENUMNODE * probdim)
        dserror("Given pairedvector provides an insufficient capacity.");
  }

  // (0) calculate the projection matrix into the tangential plain
  LINALG::Matrix<probdim, probdim> tproj_mat(true);
  ProjectionIntoTangentialPlain(unit_normal, tproj_mat);

  // (1) calculate the product of projection matrix and the first derivative
  //     of the non-unit normal. Additionally, the result is scaled by the
  //     reciprocal length of the non-unit normal.
  for (unsigned n_dof = 0; n_dof < probdim; ++n_dof)
  {
    Deriv1stMap& d_un_dof = d_unit_normal[n_dof];

    for (unsigned nn_dof = 0; nn_dof < probdim; ++nn_dof)
    {
      const Deriv1stMap& d_nun_dof = d_non_unit_normal[nn_dof];
      for (auto& d_nun_dof_var : d_nun_dof)
      {
        d_un_dof(d_nun_dof_var.first) +=
            length_n_inv * tproj_mat(n_dof, nn_dof) * d_nun_dof_var.second;
      }
    }
  }

  GEN_DATA::complete(d_unit_normal);

#ifdef CONTACT_AUG_JACOBIAN_DEBUG_OUTPUT
  for (unsigned n_dof = 0; n_dof < probdim; ++n_dof)
    for (auto& pair : d_unit_normal[n_dof])
    {
      std::cout << " n-dof (" << n_dof << "), varied GID (" << pair.first << ") = " << pair.second
                << "\n";
    }
#endif
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype>
double CONTACT::AUG::BaseSlaveIntPolicy<probdim, slavetype>::UnitSlaveElementNormal(
    const MORTAR::MortarElement& sele, const LINALG::Matrix<3, 2>& tau,
    LINALG::Matrix<probdim, 1>& unit_normal) const
{
  // 1-st column view
  const LINALG::Matrix<3, 1> tau_1(&tau(0, 0), true);
  // 2-nd column view
  const LINALG::Matrix<3, 1> tau_2(&tau(0, 1), true);
  // necessary for the 2-d case
  LINALG::Matrix<3, 1> non_unit_normal_extended;
  non_unit_normal_extended.CrossProduct(tau_1, tau_2);

  std::copy(non_unit_normal_extended.A(), non_unit_normal_extended.A() + probdim, unit_normal.A());

  double length_n_inv = unit_normal.Norm2();
  if (length_n_inv == 0.0) dserror("The length of the slave element normal is equal to 0.0!");

  length_n_inv = 1.0 / length_n_inv;
  const double normal_fac = sele.NormalFac();
  unit_normal.Scale(normal_fac * length_n_inv);

  return (length_n_inv);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype>
void CONTACT::AUG::BaseSlaveIntPolicy<probdim, slavetype>::Deriv1st_NonUnitSlaveElementNormal(
    const MORTAR::MortarElement& sele, const LINALG::Matrix<probdim, SLAVENUMNODE, int>& nodal_dofs,
    const LINALG::Matrix<SLAVEDIM, SLAVENUMNODE>& deriv, const LINALG::Matrix<3, 2>& tau,
    Deriv1stVecMap& d_non_unit_normal) const
{
#ifdef CONTACT_AUG_JACOBIAN_DEBUG_OUTPUT
  std::cout << __PRETTY_FUNCTION__ << std::endl;
#endif

  this->timer_.start(TimeID::Deriv1st_NonUnitSlaveElementNormal);

  const double normal_fac = sele.NormalFac();
  GEN_DATA::reset(probdim, SLAVENUMNODE * probdim, d_non_unit_normal);

  // loop over all components of the normal vector
  for (unsigned n_dof = 0; n_dof < probdim; ++n_dof)
  {
    Deriv1stMap& d_n_dof = d_non_unit_normal[n_dof];

    // loop over all nodes of the current element for the variation
    for (unsigned n_j = 0; n_j < SLAVENUMNODE; ++n_j)
    {
      // loop over all varied degrees of freedom
      for (unsigned j = 0; j < probdim; ++j)
      {
        // skip equal indices since Levi Civita symbol would be equal to 0
        if (j == n_dof) continue;

        const int var_jdof = nodal_dofs(j, n_j);

        double& d_nn_var_jdof = d_n_dof(var_jdof);

        // loop over the full 3-component convective base vectors, since in 2-D
        // the used cross product makes the extension necessary:
        //                tau(k,1) = { 0.0, 0.0, 1.0 }
        for (unsigned k = 0; k < 3; ++k)
        {
          // skip equal indices since Levi Civita symbol would be equal to 0
          if (j == k or n_dof == k) continue;

          switch (SLAVEDIM)
          {
            // --- 3-D case only ----------------------------------------------
            // Note: deriv(n_j,1) is equal to zero in the 2-D case.
            case 2:
            {
              const double e_ikj = INTEGRATOR::LeviCivitaSymbol(n_dof, k, j);
              d_nn_var_jdof += e_ikj * normal_fac * deriv(1, n_j) * tau(k, 0);
            }
            // no break
            // --- 3-D and 2-D case -------------------------------------------
            case 1:
            {
              const double e_ijk = INTEGRATOR::LeviCivitaSymbol(n_dof, j, k);
              d_nn_var_jdof += e_ijk * normal_fac * deriv(0, n_j) * tau(k, 1);
              break;
            }
            default:
            {
              dserror("Unsupported slave element dimension! (dim=%d)", SLAVEDIM);
              exit(EXIT_FAILURE);
            }
          }
#ifdef CONTACT_AUG_JACOBIAN_DEBUG_OUTPUT
          std::cout << "n-dof (" << n_dof << "), "
                    << "varied GID (" << var_jdof << ") | "
                    << "e_ijk = " << e_ijk << ", e_ikj = " << e_ikj << " = " << d_nn_var_jdof
                    << std::endl;
#endif
        }
      }
    }
  }

  GEN_DATA::complete(d_non_unit_normal);

  this->timer_.stop(TimeID::Deriv1st_NonUnitSlaveElementNormal);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype>
void CONTACT::AUG::BaseSlaveIntPolicy<probdim, slavetype>::Deriv2nd_NonUnitSlaveElementNormal(
    const MORTAR::MortarElement& sele, const LINALG::Matrix<probdim, SLAVENUMNODE, int>& nodal_dofs,
    const LINALG::Matrix<SLAVEDIM, SLAVENUMNODE>& deriv, Deriv2ndVecMap& dd_non_unit_normal) const
{
#ifdef CONTACT_AUG_JACOBIAN_DEBUG_OUTPUT
  std::cout << CONTACT_FUNC_NAME << std::endl;
#endif

  this->timer_.start(TimeID::Deriv2nd_NonUnitSlaveElementNormal);

  const double normal_fac = sele.NormalFac();
  // loop over all components of the normal vector
  for (unsigned n_dof = 0; n_dof < probdim; ++n_dof)
  {
    Deriv2ndMap& dd_n_dof = dd_non_unit_normal[n_dof];

    GEN_DATA::reset(SLAVENUMNODE * probdim, dd_n_dof);

    // loop over all nodes of the current element for the variation
    for (unsigned n_j = 0; n_j < SLAVENUMNODE; ++n_j)
    {
      // loop over varied degrees of freedom
      for (unsigned j = 0; j < probdim; ++j)
      {
        // skip equal indices since Levi Civita symbol e_ijk would be equal to 0
        if (j == n_dof) continue;

        // GID of the varied dof
        const int var_jdof = nodal_dofs(j, n_j);

        Deriv1stMap& dd_nn_var_jdof = dd_n_dof[var_jdof];

        for (unsigned n_k = 0; n_k < SLAVENUMNODE; ++n_k)
        {
          // loop over linearization degrees of freedom
          for (unsigned k = 0; k < probdim; ++k)
          {
            // skip equal indices since Levi Civita symbol would be equal to 0
            if (k == j or k == n_dof) continue;

            // GID of the linearized DOF
            const int lin_kdof = nodal_dofs(k, n_k);

            const double e_ijk = INTEGRATOR::LeviCivitaSymbol(n_dof, j, k);
            dd_nn_var_jdof(lin_kdof) += e_ijk * normal_fac * deriv(0, n_j) * deriv(1, n_k);

            const double e_ikj = INTEGRATOR::LeviCivitaSymbol(n_dof, k, j);
            dd_nn_var_jdof(lin_kdof) += e_ikj * normal_fac * deriv(1, n_j) * deriv(0, n_k);

#ifdef CONTACT_AUG_JACOBIAN_DEBUG_OUTPUT
            std::cout << n_dof << " | " << j << " (" << var_jdof << ") | " << k << " (" << lin_kdof
                      << ") : " << e_ijk << ", " << e_ikj << "; " << dd_nn_var_jdof[lin_kdof]
                      << std::endl;
#endif
          }
        }
      }
    }
  }

  GEN_DATA::complete(dd_non_unit_normal);

  this->timer_.stop(TimeID::Deriv2nd_NonUnitSlaveElementNormal);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype>
void CONTACT::AUG::BaseSlaveIntPolicy<probdim, slavetype>::Deriv2nd_UnitSlaveElementNormal(
    const LINALG::Matrix<probdim, 1>& unit_normal, const double length_n_inv,
    const Deriv1stVecMap& d_non_unit_normal, const Deriv1stVecMap& d_unit_normal,
    const Deriv2ndVecMap& dd_non_unit_normal, Deriv2ndVecMap& dd_unit_normal) const
{
  this->timer_.start(TimeID::Deriv2nd_UnitSlaveElementNormal);

  // temporal data structures
  Deriv1stMap dn_n(0);
  GEN_DATA::reset(d_non_unit_normal[0].capacity(), dn_n);
  Deriv2ndMap dn_dn(0);
  GEN_DATA::reset(d_non_unit_normal[0].capacity(), d_unit_normal[0].capacity(), dn_dn);

  /*--------------------------------------------------------------------------*/
  // (0-0) evaluate the scalar product of unit-normal vector and the first
  //       derivative of the non-unit normal vector
  // loop over the normal vector components
  InnerProductOfVectorAndDeriv1stVector(unit_normal, d_non_unit_normal, dn_n);

  // (0-1) Multiply the scalar product calculated in (0-0) with the varied
  //       unit normal and scale the result by the negative inverse normal
  //       length.
  for (unsigned n_dof = 0; n_dof < probdim; ++n_dof)
  {
    const Deriv1stMap& d_un_ndof = d_unit_normal[n_dof];
    Deriv2ndMap& dd_un_ndof = dd_unit_normal[n_dof];

    for (auto& d_un_ndof_var : d_un_ndof)
    {
      Deriv1stMap& dd_un_ndof_var = dd_un_ndof.repetitive_access(d_un_ndof_var.first, gp_id_);

      for (auto& dn_n_lin : dn_n)
      {
        dd_un_ndof_var.repetitive_access(dn_n_lin.first, gp_id_) -=
            length_n_inv * dn_n_lin.second * d_un_ndof_var.second;
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  // (1-0) evaluate the scalar product of unit normal vector and the first
  //       derivative of the non-unit normal vector
  // reuse the previously calculated quantity ( see (0-0) )

  // (1-1) Multiply the scalar product calculated in (1-0) with the linearized
  //       unit normal and scale the result by the negative inverse normal
  //       length.
  for (auto& dn_n_var : dn_n)
  {
    for (unsigned n_dof = 0; n_dof < probdim; ++n_dof)
    {
      const Deriv1stMap& d_un_ndof = d_unit_normal[n_dof];
      Deriv2ndMap& dd_un_ndof = dd_unit_normal[n_dof];
      Deriv1stMap& dd_un_ndof_var = dd_un_ndof.repetitive_access(dn_n_var.first, gp_id_);

      for (auto& d_un_ndof_lin : d_un_ndof)
      {
        dd_un_ndof_var.repetitive_access(d_un_ndof_lin.first, gp_id_) -=
            length_n_inv * dn_n_var.second * d_un_ndof_lin.second;
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  // (2-0) evaluate the scalar product of the linearized unit normal vector and
  //       the varied non-unit normal vector
  for (unsigned n_dof = 0; n_dof < probdim; ++n_dof)
  {
    const Deriv1stMap& d_nun_ndof = d_non_unit_normal[n_dof];
    const Deriv1stMap& d_un_ndof = d_unit_normal[n_dof];

    for (auto& d_nun_ndof_var : d_nun_ndof)
    {
      Deriv1stMap& dn_dn_var = dn_dn[d_nun_ndof_var.first];

      for (auto& d_un_ndof_lin : d_un_ndof)
      {
        dn_dn_var(d_un_ndof_lin.first) += d_nun_ndof_var.second * d_un_ndof_lin.second;
      }
    }
  }

  // complete the auxiliary variable
  dn_dn.complete();

  // (2-1) Multiply the scalar product calculated in (2-0) with the unit normal
  //       vector and scale the result by the negative inverse normal length
  for (unsigned n_dof = 0; n_dof < probdim; ++n_dof)
  {
    Deriv2ndMap& dd_un_ndof = dd_unit_normal[n_dof];
    const double un_ndof = unit_normal(n_dof, 0);

    for (auto& dn_dn_var : dn_dn)
    {
      Deriv1stMap& dd_un_ndof_var = dd_un_ndof.repetitive_access(dn_dn_var.first, gp_id_);
      for (auto& dn_dn_var_lin : dn_dn_var.second)
      {
        dd_un_ndof_var.repetitive_access(dn_dn_var_lin.first, gp_id_) -=
            length_n_inv * dn_dn_var_lin.second * un_ndof;
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  // (3-0) Calculate the projection matrix into the tangential plain
  LINALG::Matrix<probdim, probdim> tproj_mat(false);
  ProjectionIntoTangentialPlain(unit_normal, tproj_mat);

  // (3-1) Multiply the tangential projection matrix with the second derviatives
  //       of the non-unit normal and scale the result with the inverse
  //       normal length
  for (unsigned n_dof = 0; n_dof < probdim; ++n_dof)
  {
    Deriv2ndMap& dd_un_ndof = dd_unit_normal[n_dof];

    for (unsigned nn_dof = 0; nn_dof < probdim; ++nn_dof)
    {
      const Deriv2ndMap& dd_nun_ndof = dd_non_unit_normal[nn_dof];
      for (auto& dd_nun_ndof_var : dd_nun_ndof)
      {
        Deriv1stMap& dd_un_ndof_var = dd_un_ndof.repetitive_access(dd_nun_ndof_var.first, gp_id_);

        for (auto& dd_nun_ndof_var_lin : dd_nun_ndof_var.second)
        {
          dd_un_ndof_var.repetitive_access(dd_nun_ndof_var_lin.first, gp_id_) +=
              length_n_inv * tproj_mat(n_dof, nn_dof) * dd_nun_ndof_var_lin.second;
        }
      }
    }
  }

  GEN_DATA::complete(dd_unit_normal);
  this->timer_.stop(TimeID::Deriv2nd_UnitSlaveElementNormal);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype>
void CONTACT::AUG::BaseSlaveIntPolicy<probdim, slavetype>::InnerProductOfVectorAndDeriv1stVector(
    const LINALG::Matrix<probdim, 1>& vec, const Deriv1stVecMap& d_vec, Deriv1stMap& dvec_vec) const
{
  // loop over the normal vector components
  for (unsigned dof = 0; dof < probdim; ++dof)
  {
    const Deriv1stMap& d_v_dof = d_vec[dof];
    const double v_dof = vec(dof, 0);

    for (auto& d_v_dof_pair : d_v_dof)
    {
      dvec_vec(d_v_dof_pair.first) += v_dof * d_v_dof_pair.second;
    }
  }

  dvec_vec.complete();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype>
void CONTACT::AUG::BaseSlaveIntPolicy<probdim, slavetype>::ProjectionIntoTangentialPlain(
    const LINALG::Matrix<probdim, 1>& unit_normal,
    LINALG::Matrix<probdim, probdim>& tproj_mat) const
{
  tproj_mat.MultiplyNT(-1.0, unit_normal, unit_normal);

  for (unsigned i = 0; i < probdim; ++i) tproj_mat(i, i) += 1.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::BaseIntPolicy<probdim, slavetype, mastertype>::LMatrixInverse(
    const LINALG::Matrix<3, 2>& mtau, const LINALG::Matrix<probdim, 1>& snormal,
    LINALG::Matrix<probdim, probdim>& lmat_inv) const
{
  for (unsigned c = 0; c < probdim; ++c)
  {
    if (c < MASTERDIM)
      std::copy(&mtau(0, c), &mtau(0, c) + probdim, &lmat_inv(0, c));
    else
    {
      LINALG::Matrix<probdim, 1> last_col(&lmat_inv(0, c), true);
      last_col.Update(-1.0, snormal);
    }
  }

  if (std::abs(lmat_inv.Determinant()) < 1.0e-15)
  {
    lmat_inv.Print(std::cout);
    dserror("L-matrix is almost singular!");
  }

  lmat_inv.Invert();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype>
void CONTACT::AUG::BaseSlaveIntPolicy<probdim, slavetype>::AveragedNormalAtXi(
    MORTAR::MortarElement& sele, const LINALG::Matrix<SLAVENUMNODE, 1>& sval,
    LINALG::Matrix<probdim, 1>& snormal) const
{
  std::fill(snormal.A(), snormal.A() + probdim, 0.0);

  const DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < SLAVENUMNODE; ++i)
  {
    const CoNode& cnode = static_cast<const CoNode&>(*snodes[i]);
    const LINALG::Matrix<probdim, 1> nn(cnode.MoData().n(), true);
    snormal.Update(sval(i, 0), nn, 1.0);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype>
void CONTACT::AUG::BaseSlaveIntPolicy<probdim, slavetype>::CompleteNodeData(
    MORTAR::MortarElement& sele) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < SLAVENUMNODE; ++i)
  {
    CoNode& cnode = dynamic_cast<CoNode&>(*snodes[i]);
    cnode.AugData().Complete();
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::BaseIntPolicy<probdim, slavetype, mastertype>::Deriv1st_MXiGP(
    const LINALG::Matrix<probdim, probdim>& lmat_inv, MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
    const LINALG::Matrix<MASTERNUMNODE, 1>& mval, const double alpha, Deriv1stVecMap& d_mxi,
    Deriv1stMap& d_alpha) const
{
  this->timer_.start(TimeID::Deriv1st_MXiGP);

  const DRT::Node* const* snodes = sele.Nodes();

  LINALG::Matrix<probdim, my::SLAVENUMNODE, int> snodal_dofs;
  CONTACT::INTEGRATOR::GetElementNodalDofs(sele, snodal_dofs);

  LINALG::Matrix<probdim, MASTERNUMNODE, int> mnodal_dofs;
  CONTACT::INTEGRATOR::GetElementNodalDofs(mele, mnodal_dofs);

  for (unsigned i = 0; i < probdim; ++i)
  {
    // switch between the 1-st order derivative of the master parametric
    // coordinates and the 1-st order derivative of the projection length
    // parameter alpha.
    Deriv1stMap* deriv_i_ptr = NULL;
    if (i < MASTERDIM)
      deriv_i_ptr = &d_mxi[i];
    else
      deriv_i_ptr = &d_alpha;

    Deriv1stMap& deriv_i = *deriv_i_ptr;

    // (0) varied slave dofs
    for (unsigned j = 0; j < my::SLAVENUMNODE; ++j)
    {
      for (unsigned k = 0; k < probdim; ++k)
      {
        deriv_i.repetitive_access(snodal_dofs(k, j), my::gp_id_) += lmat_inv(i, k) * sval(j);
      }
    }

    // (1) varied master dofs
    for (unsigned j = 0; j < MASTERNUMNODE; ++j)
    {
      for (unsigned k = 0; k < probdim; ++k)
      {
        deriv_i.repetitive_access(mnodal_dofs(k, j), my::gp_id_) -= lmat_inv(i, k) * mval(j);
      }
    }

    // (2) varied smooth normal
    for (unsigned j = 0; j < my::SLAVENUMNODE; ++j)
    {
      const CoNode& cnode = static_cast<const CoNode&>(*snodes[j]);
      const Deriv1stVecMap& d_n = cnode.AugData().GetDeriv1st_N();

      for (unsigned k = 0; k < probdim; ++k)
      {
        const double tmp = lmat_inv(i, k) * sval(j) * alpha;

        for (auto& d_n_k : d_n[k])
          deriv_i.repetitive_access(d_n_k.first, my::gp_id_) += tmp * d_n_k.second;
      }
    }
  }

  GEN_DATA::complete(d_mxi);
  d_alpha.complete();

  this->timer_.stop(TimeID::Deriv1st_MXiGP);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::CompleteIntPolicy<probdim, slavetype, mastertype>::Get_Deriv2nd_MXiGP(
    const LINALG::Matrix<probdim, probdim>& lmat_inv, MORTAR::MortarElement& sele,
    MORTAR::MortarElement& mele, const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
    const LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
    const LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
    const LINALG::Matrix<3, my::MASTERNUMNODE>& mderiv2, const LINALG::Matrix<3, 2>& mtau,
    const double* mxi, const double alpha, const Deriv1stVecMap& d_mxigp,
    const Deriv1stMap& d_alpha, Deriv2ndVecMap& dd_mxigp) const
{
  LINALG::Matrix<probdim, my::MASTERNUMNODE, int> mnodal_dofs;
  CONTACT::INTEGRATOR::GetElementNodalDofs(mele, mnodal_dofs);

  LINALG::Matrix<3, my::MASTERNUMNODE> mcoord;
  mele.GetNodalCoords(mcoord);

  // --- contributions from the master side -----------------------------------
  this->Add_Deriv2nd_MaDispl(lmat_inv, mnodal_dofs, mderiv, d_mxigp, dd_mxigp);

  this->Add_Deriv1st_MaMetric(
      lmat_inv, mnodal_dofs, mtau, mderiv, mderiv2, mcoord, d_mxigp, dd_mxigp);

  // --- contributions from the slave side -------------------------------------
  this->Add_Deriv1st_Alpha_Deriv1st_Normal(lmat_inv, d_alpha, sele, sval, dd_mxigp);

  this->Add_Alpha_Deriv2nd_Normal(lmat_inv, alpha, sele, sval, dd_mxigp);

  GEN_DATA::complete(dd_mxigp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::BaseIntPolicy<probdim, slavetype, mastertype>::Add_Deriv2nd_MaDispl(
    const LINALG::Matrix<probdim, probdim>& lmat_inv,
    const LINALG::Matrix<probdim, MASTERNUMNODE, int>& mnodal_dofs,
    const LINALG::Matrix<MASTERDIM, MASTERNUMNODE>& mderiv, const Deriv1stVecMap& d_mxigp,
    Deriv2ndVecMap& dd_mxigp) const
{
  // we consider only the second derivatives of the master parametric coordinate
  // here, and neglect the projection length alpha for now.

  // loop over the parametric coordinates of the master side
  for (unsigned i = 0; i < MASTERDIM; ++i)
  {
    Deriv2ndMap& dd_mxigp_i = dd_mxigp[i];

    for (unsigned k = 0; k < MASTERNUMNODE; ++k)
    {
      for (unsigned j = 0; j < probdim; ++j)
      {
        Deriv1stMap& dd_mxigp_i_var = dd_mxigp_i.repetitive_access(mnodal_dofs(j, k), my::gp_id_);

        for (unsigned l = 0; l < MASTERDIM; ++l)
        {
          for (auto& d_mxigp_l_pair : d_mxigp[l])
          {
            dd_mxigp_i_var.repetitive_access(d_mxigp_l_pair.first, my::gp_id_) -=
                lmat_inv(i, j) * mderiv(l, k) * d_mxigp_l_pair.second;
          }
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::BaseIntPolicy<probdim, slavetype, mastertype>::Add_Deriv1st_MaMetric(
    const LINALG::Matrix<probdim, probdim>& lmat_inv,
    const LINALG::Matrix<probdim, MASTERNUMNODE, int>& mnodal_dofs,
    const LINALG::Matrix<3, 2>& mtau, const LINALG::Matrix<MASTERDIM, MASTERNUMNODE>& mderiv,
    const LINALG::Matrix<3, MASTERNUMNODE>& mderiv2nd,
    const LINALG::Matrix<3, MASTERNUMNODE>& mcoord, const Deriv1stVecMap& d_mxigp,
    Deriv2ndVecMap& dd_mxigp) const
{
  static const unsigned dd_2nd_indices[2][2] = {{0, 2}, {2, 1}};

  for (unsigned i = 0; i < MASTERDIM; ++i)
  {
    Deriv2ndMap& dd_mxigp_i = dd_mxigp[i];

    for (unsigned j = 0; j < MASTERDIM; ++j)
    {
      const Deriv1stMap& d_mxigp_j = d_mxigp[j];

      for (auto& d_mxi_gp_j_pair : d_mxigp_j)
      {
        Deriv1stMap& dd_mxigp_i_var =
            dd_mxigp_i.repetitive_access(d_mxi_gp_j_pair.first, my::gp_id_);

        // linearization of the master position multiplied by the varied
        // master parametric coordinates
        for (unsigned k = 0; k < MASTERNUMNODE; ++k)
        {
          const double tmp = mderiv(j, k) * d_mxi_gp_j_pair.second;

          for (unsigned d = 0; d < probdim; ++d)
          {
            dd_mxigp_i_var.repetitive_access(mnodal_dofs(d, k), my::gp_id_) -= lmat_inv(i, d) * tmp;
          }
        }

        // 1st derivative of the convective base vectors multiplied by the
        // linearized and varied parametric coordinates
        for (unsigned l = 0; l < MASTERDIM; ++l)
        {
          // index for the 2-nd derivative:
          //        N_{,11}    | N_{,12}    | N_{,21}    | N_{,21} (== N_{,22})
          // (j,l): (0,0) -> 0 | (0,1) -> 2 | (1,0) -> 2 | (1,1) -> 1
          const unsigned j_2nd = dd_2nd_indices[j][l];
          //          std::cout << "(" << j << "," << l << ") -> " << j_2nd << "\n";

          const Deriv1stMap& d_mxigp_l = d_mxigp[l];

          double tmp = 0.0;

          for (unsigned d = 0; d < probdim; ++d)
          {
            for (unsigned k = 0; k < MASTERNUMNODE; ++k)
            {
              tmp += lmat_inv(i, d) * mderiv2nd(j_2nd, k) * mcoord(d, k);
            }
          }

          tmp *= d_mxi_gp_j_pair.second;

          for (auto& d_mxigp_l_pair : d_mxigp_l)
          {
            dd_mxigp_i_var.repetitive_access(d_mxigp_l_pair.first, my::gp_id_) -=
                tmp * d_mxigp_l_pair.second;
          }
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::BaseIntPolicy<probdim, slavetype,
    mastertype>::Add_Deriv1st_Alpha_Deriv1st_Normal(const LINALG::Matrix<probdim, probdim>&
                                                        lmat_inv,
    const Deriv1stMap& d_alpha, const MORTAR::MortarElement& sele,
    const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval, Deriv2ndVecMap& dd_mxigp) const
{
  const DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < MASTERDIM; ++i)
  {
    Deriv2ndMap& dd_mxigp_i = dd_mxigp[i];

    // (0) varied alpha multiplied with the linearized smooth normal vector
    for (auto& d_alpha_var : d_alpha)
    {
      Deriv1stMap& dd_mxigp_i_var = dd_mxigp_i.repetitive_access(d_alpha_var.first, my::gp_id_);

      for (unsigned j = 0; j < my::SLAVENUMNODE; ++j)
      {
        const double tmp1 = d_alpha_var.second * sval(j, 0);

        const CoNode& cnode = static_cast<const CoNode&>(*snodes[j]);
        const Deriv1stVecMap& d_n = cnode.AugData().GetDeriv1st_N();

        // ToDo Probably, it is possible to reduce the index insert look-up.
        for (unsigned k = 0; k < probdim; ++k)
        {
          const double tmp2 = tmp1 * lmat_inv(i, k);

          const Deriv1stMap& d_n_k = d_n[k];

          for (auto& d_n_k_lin : d_n_k)
          {
            //            std::cout << k << ", " << d_n_k_lin.first << std::endl;
            dd_mxigp_i_var.repetitive_access(d_n_k_lin.first, my::gp_id_) +=
                tmp2 * d_n_k_lin.second;
          }
        }
      }
    }

    // (1) linearized alpha multiplied with the varied smooth normal vector
    for (unsigned j = 0; j < my::SLAVENUMNODE; ++j)
    {
      const CoNode& cnode = static_cast<const CoNode&>(*snodes[j]);
      const Deriv1stVecMap& d_n = cnode.AugData().GetDeriv1st_N();

      for (unsigned k = 0; k < probdim; ++k)
      {
        const Deriv1stMap& d_n_k = d_n[k];

        const double tmp1 = sval(j, 0) * lmat_inv(i, k);

        for (auto& d_n_k_var : d_n_k)
        {
          const double tmp2 = tmp1 * d_n_k_var.second;

          Deriv1stMap& dd_mxigp_i_var = dd_mxigp_i.repetitive_access(d_n_k_var.first, my::gp_id_);

          for (auto& d_alpha_lin : d_alpha)
          {
            dd_mxigp_i_var.repetitive_access(d_alpha_lin.first, my::gp_id_) +=
                tmp2 * d_alpha_lin.second;
          }
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::BaseIntPolicy<probdim, slavetype, mastertype>::Add_Alpha_Deriv2nd_Normal(
    const LINALG::Matrix<probdim, probdim>& lmat_inv, const double alpha,
    const MORTAR::MortarElement& sele, const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
    Deriv2ndVecMap& dd_mxigp) const
{
  const DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < MASTERDIM; ++i)
  {
    Deriv2ndMap& dd_mxigp_i = dd_mxigp[i];

    for (unsigned j = 0; j < my::SLAVENUMNODE; ++j)
    {
      const CoNode& cnode = static_cast<const CoNode&>(*snodes[j]);
      const Deriv2ndVecMap& dd_n = cnode.AugData().GetDeriv2nd_N();

      double tmp1 = alpha * sval(j, 0);

      for (unsigned k = 0; k < probdim; ++k)
      {
        const double tmp2 = tmp1 * lmat_inv(i, k);

        for (auto& dd_n_var : dd_n[k])
        {
          Deriv1stMap& dd_mxigp_i_var = dd_mxigp_i.repetitive_access(dd_n_var.first, my::gp_id_);

          for (auto& dd_n_var_lin : dd_n_var.second)
          {
            dd_mxigp_i_var.repetitive_access(dd_n_var_lin.first, my::gp_id_) +=
                tmp2 * dd_n_var_lin.second;
          }
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::IncompleteIntPolicy<probdim, slavetype, mastertype>::Get_Deriv1st_GapN(
    MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
    const LINALG::Matrix<my::MASTERNUMNODE, 1>& mval, const double* gpn,
    const LINALG::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp, Deriv1stMap& deriv_gapn_sl,
    Deriv1stMap& deriv_gapn_ma) const
{
  // get slave element nodes
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** mnodes = mele.Nodes();

  this->Deriv1st_GapN_Sl(snodes, sval, gpn, deriv_gapn_sl);

  // NOTE: skip the part related to the projected parametric coordinates
  const Deriv1stVecMap empty_d_mxigp(my::MASTERDIM, Deriv1stMap(0));
  this->Deriv1st_GapN_Ma(mnodes, mval, gpn, mtau, empty_d_mxigp, deriv_gapn_ma);

  deriv_gapn_sl.complete();
  deriv_gapn_ma.complete();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::CompleteIntPolicy<probdim, slavetype, mastertype>::Get_Deriv1st_GapN(
    MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
    const LINALG::Matrix<my::MASTERNUMNODE, 1>& mval, const double* gpn,
    const LINALG::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp, Deriv1stMap& deriv_gapn_sl,
    Deriv1stMap& deriv_gapn_ma) const
{
  // get slave element nodes
  DRT::Node** snodes = sele.Nodes();
  DRT::Node** mnodes = mele.Nodes();

  this->Deriv1st_GapN_Sl(snodes, sval, gpn, deriv_gapn_sl);

  this->Deriv1st_GapN_Ma(mnodes, mval, gpn, mtau, d_mxigp, deriv_gapn_ma);

  deriv_gapn_sl.complete();
  deriv_gapn_ma.complete();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::BaseIntPolicy<probdim, slavetype, mastertype>::Deriv1st_GapN_Sl(
    const DRT::Node* const* snodes, const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
    const double* gpn, Deriv1stMap& deriv_gapn_sl) const
{
  for (unsigned k = 0; k < my::SLAVENUMNODE; ++k)
  {
    const CoNode& snode = static_cast<const CoNode&>(*snodes[k]);
    const int* sdof = snode.Dofs();

    // variation of the slave position
    for (unsigned d = 0; d < probdim; ++d) deriv_gapn_sl(sdof[d]) += sval(k, 0) * gpn[d];
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::BaseIntPolicy<probdim, slavetype, mastertype>::Deriv1st_GapN_Ma(
    const DRT::Node* const* mnodes, const LINALG::Matrix<MASTERNUMNODE, 1>& mval, const double* gpn,
    const LINALG::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp,
    Deriv1stMap& deriv_gapn_ma) const
{
  for (unsigned k = 0; k < MASTERNUMNODE; ++k)
  {
    const CoNode& mnode = static_cast<const CoNode&>(*mnodes[k]);
    const int* mdof = mnode.Dofs();

    // (0) variation of the master position
    for (unsigned d = 0; d < probdim; ++d)
    {
      deriv_gapn_ma(mdof[d]) += mval(k, 0) * gpn[d];
    }
  }

  // (1) derivative of the master convective parameter coordinate
  //     ( these contributions are skipped for the incomplete variant )
  for (unsigned j = 0; j < MASTERDIM; ++j)
  {
    for (auto& d_mxigp_j : d_mxigp[j])
    {
      double& deriv_gapn_ma_j = deriv_gapn_ma(d_mxigp_j.first);

      for (unsigned d = 0; d < probdim; ++d)
      {
        deriv_gapn_ma_j += gpn[d] * mtau(d, j) * d_mxigp_j.second;
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::IncompleteIntPolicy<probdim, slavetype, mastertype>::Get_Deriv1st_WGap(
    MORTAR::MortarElement& sele, const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const double gapn_sl, const double gapn_ma, const double wgt, const double jac,
    const Deriv1stMap& d_jac, const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma) const
{
  // get slave element nodes
  DRT::Node** snodes = sele.Nodes();

  this->Add_Deriv1stGapNContributions(snodes, wgt * jac, lmval, d_gapn_sl, d_gapn_ma);

  // NOTE: skip contributions related to the 1-st order derivative of the
  // element jacobian
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::CompleteIntPolicy<probdim, slavetype, mastertype>::Get_Deriv1st_WGap(
    MORTAR::MortarElement& sele, const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const double gapn_sl, const double gapn_ma, const double wgt, const double jac,
    const Deriv1stMap& d_jac, const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma) const
{
  // get slave element nodes
  DRT::Node** snodes = sele.Nodes();

  this->Add_Deriv1stGapNContributions(snodes, wgt * jac, lmval, d_gapn_sl, d_gapn_ma);

  this->Add_Deriv1stJacobianContributions(snodes, wgt, lmval, gapn_sl, gapn_ma, d_jac);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::IncompleteIntPolicy<probdim, slavetype, mastertype>::Get_Deriv1st_WGap_Complete(
    const int linsize, MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
    const LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
    const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double* gpn,
    const LINALG::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp, const double gapn_sl,
    const double gapn_ma, const double wgt, const double jac, const Deriv1stMap& d_jac) const
{
  // create an instance of the complete integration strategy, such that we
  // can reuse the already existing implementations
  const CompleteIntPolicy<probdim, slavetype, mastertype> complete_policy;

  Deriv1stMap deriv_gapn_sl(probdim * my::SLAVENUMNODE);
  Deriv1stMap deriv_gapn_ma(linsize + probdim * my::MASTERNUMNODE);

  /*--------------------------------------------------------------------------*/
  // (0) evaluate the discrete gap ( complete )
  complete_policy.Get_Deriv1st_GapN(
      sele, mele, sval, mval, gpn, mtau, d_mxigp, deriv_gapn_sl, deriv_gapn_ma);

  /*--------------------------------------------------------------------------*/
  // (1) swap the paired vector content
  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    DRT::Node* const* snodes = sele.Nodes();

    CoNode& cnode = static_cast<CoNode&>(*snodes[i]);

    // swap slave contributions
    {
      Deriv1stMap& wgap_sl = cnode.AugData().GetDeriv1st_WGapSl();
      Deriv1stMap& wgap_sl_complete = cnode.AugData().GetDeriv1st_WGapSl_Complete();
      wgap_sl_complete.swap(wgap_sl);
    }

    // swap master contributions
    {
      Deriv1stMap& wgap_ma = cnode.AugData().GetDeriv1st_WGapMa();
      Deriv1stMap& wgap_ma_complete = cnode.AugData().GetDeriv1st_WGapMa_Complete();
      wgap_ma_complete.swap(wgap_ma);
    }
  }

  /*--------------------------------------------------------------------------*/
  // (2) call the evaluate and assemble routines of the CompleteIntPolicy
  complete_policy.Get_Deriv1st_WGap(
      sele, lmval, gapn_sl, gapn_ma, wgt, jac, d_jac, deriv_gapn_sl, deriv_gapn_ma);

  /*--------------------------------------------------------------------------*/
  // (3) swap the paired vector content back
  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    DRT::Node* const* snodes = sele.Nodes();

    CoNode& cnode = static_cast<CoNode&>(*snodes[i]);

    // swap slave contributions
    {
      Deriv1stMap& wgap_sl = cnode.AugData().GetDeriv1st_WGapSl();
      Deriv1stMap& wgap_sl_complete = cnode.AugData().GetDeriv1st_WGapSl_Complete();
      wgap_sl_complete.swap(wgap_sl);
    }

    // swap master contributions
    {
      Deriv1stMap& wgap_ma = cnode.AugData().GetDeriv1st_WGapMa();
      Deriv1stMap& wgap_ma_complete = cnode.AugData().GetDeriv1st_WGapMa_Complete();
      wgap_ma_complete.swap(wgap_ma);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::CompleteIntPolicy<probdim, slavetype, mastertype>::Get_Deriv1st_WGap_Complete(
    const int linsize, MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
    const LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
    const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double* gpn,
    const LINALG::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp, const double gapn_sl,
    const double gapn_ma, const double wgt, const double jac, const Deriv1stMap& d_jac) const
{
  // get slave element nodes
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    CoNode& cnode = static_cast<CoNode&>(*snodes[i]);

    Deriv1stMap& d_wgap_sl = cnode.AugData().GetDeriv1st_WGapSl();
    Deriv1stMap& d_wgap_ma = cnode.AugData().GetDeriv1st_WGapMa();

    // reset the pointers to the already calculated complete first order
    // derivative values

    // slave
    Teuchos::RCP<Deriv1stMap>& d_wgap_sl_complete_ptr =
        cnode.AugData().GetDeriv1st_WGapSl_Complete_Ptr();
    d_wgap_sl_complete_ptr = Teuchos::rcpFromRef(d_wgap_sl);

    // master
    Teuchos::RCP<Deriv1stMap>& d_wgap_ma_complete_ptr =
        cnode.AugData().GetDeriv1st_WGapMa_Complete_Ptr();
    d_wgap_ma_complete_ptr = Teuchos::rcpFromRef(d_wgap_ma);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::BaseIntPolicy<probdim, slavetype, mastertype>::Add_Deriv1stGapNContributions(
    DRT::Node* const* snodes, const double scale, const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma) const
{
  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    const double tmp = scale * lmval(i, 0);
    CoNode& cnode = static_cast<CoNode&>(*snodes[i]);

    /*------------------------------------------------------------------------*/
    // slave contributions
    {
      Deriv1stMap& d_wgap_sl = cnode.AugData().GetDeriv1st_WGapSl();
      for (auto& d_gapn_sl_var : d_gapn_sl)
      {
        d_wgap_sl(d_gapn_sl_var.first) += tmp * d_gapn_sl_var.second;
      }
    }

    /*------------------------------------------------------------------------*/
    // master contributions
    {
      Deriv1stMap& d_wgap_ma = cnode.AugData().GetDeriv1st_WGapMa();
      for (auto& d_gapn_ma_var : d_gapn_ma)
      {
        GEN_DATA::increaseCapacity(d_wgap_ma);
        d_wgap_ma(d_gapn_ma_var.first) += tmp * d_gapn_ma_var.second;
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::BaseIntPolicy<probdim, slavetype, mastertype>::Add_Deriv1stJacobianContributions(
    DRT::Node* const* snodes, const double wgt, const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const double gapn_sl, const double gapn_ma, const Deriv1stMap& d_jac) const
{
  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    CoNode& cnode = static_cast<CoNode&>(*snodes[i]);

    const double tmp_sl_ma = lmval(i, 0) * wgt * (gapn_sl - gapn_ma);

    // slave contributions
    Deriv1stMap& d_wgap_sl = cnode.AugData().GetDeriv1st_WGapSl();

    for (auto& d_jac_var : d_jac) d_wgap_sl(d_jac_var.first) += tmp_sl_ma * d_jac_var.second;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::IncompleteIntPolicy<probdim, slavetype, mastertype>::Get_Deriv2nd_WGap(
    MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
    const LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
    const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
    const LINALG::Matrix<3, my::MASTERNUMNODE>& mderiv2nd, const LINALG::Matrix<3, 2>& mtau,
    const double* gpn, const double wgt, const double gapn_sl, const double gapn_ma,
    const double jac, const Deriv1stMap& d_jac, const Deriv2ndMap& dd_jac,
    const Deriv1stVecMap& d_mxigp, const Deriv2ndVecMap& dd_mxigp, const Deriv1stVecMap& d_n_unit,
    const Deriv2ndVecMap& dd_n_unit, const Deriv1stMap& d_gapn_sl,
    const Deriv1stMap& d_gapn_ma) const
{
  Add_Jac_Deriv2nd_GapN(
      sele, mele, sval, mval, lmval, mderiv, mtau, gpn, wgt, jac, d_mxigp, d_n_unit, dd_n_unit);

  Add_Deriv1st_GapN_Deriv1st_Jac(sele, lmval, wgt, d_gapn_sl, d_gapn_ma, d_jac);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::IncompleteIntPolicy<probdim, slavetype, mastertype>::Add_Jac_Deriv2nd_GapN(
    MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
    const LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
    const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
    const LINALG::Matrix<3, 2>& mtau, const double* gpn, const double wgt, const double jac,
    const Deriv1stVecMap& d_mxigp, const Deriv1stVecMap& d_n_unit,
    const Deriv2ndVecMap& dd_n_unit) const
{
  this->timer_.start(TimeID::INCOMPLETE_Add_Jac_Deriv2nd_GapN);

  DRT::Node* const* snodes = sele.Nodes();
  const DRT::Node* const* mnodes = mele.Nodes();

  LINALG::Matrix<3, my::SLAVENUMNODE> scoord;
  sele.GetNodalCoords(scoord);

  LINALG::Matrix<3, 1> xs(true);
  xs.MultiplyNN(scoord, sval);

  LINALG::Matrix<3, my::MASTERNUMNODE> mcoord;
  mele.GetNodalCoords(mcoord);

  LINALG::Matrix<3, 1> xm(true);
  xm.MultiplyNN(mcoord, mval);

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    CoNode& cnode = static_cast<CoNode&>(*snodes[i]);

    // constant scaling term
    const double tmp = wgt * jac * lmval(i, 0);

    /*------------------------------------------------------------------------*/
    // slave contributions
    Deriv2ndMap& dd_wgap_sl = cnode.AugData().GetDeriv2nd_WGapSl();

    for (unsigned k = 0; k < my::SLAVENUMNODE; ++k)
    {
      const CoNode& snode = static_cast<const CoNode&>(*snodes[k]);
      const int* sdof = snode.Dofs();

      // variation of the slave position multiplied with the linearized
      // smooth normal
      for (unsigned d = 0; d < probdim; ++d)
      {
        GEN_DATA::increaseCapacity(dd_wgap_sl);
        Deriv1stMap& dd_wgap_sl_var = dd_wgap_sl.repetitive_access(sdof[d], my::gp_id_);

        for (auto& d_n_unit_lin : d_n_unit[d])
        {
          GEN_DATA::increaseCapacity(dd_wgap_sl_var);
          dd_wgap_sl_var.repetitive_access(d_n_unit_lin.first, my::gp_id_) +=
              tmp * sval(k, 0) * d_n_unit_lin.second;
        }
      }
    }

    // varied smooth normal multiplied with the linearized slave position
    for (unsigned d = 0; d < probdim; ++d)
    {
      for (auto& d_n_unit_var : d_n_unit[d])
      {
        GEN_DATA::increaseCapacity(dd_wgap_sl);
        Deriv1stMap& dd_wgap_sl_var = dd_wgap_sl.repetitive_access(d_n_unit_var.first, my::gp_id_);

        const double val = tmp * d_n_unit_var.second;

        for (unsigned k = 0; k < my::SLAVENUMNODE; ++k)
        {
          const CoNode& snode = static_cast<const CoNode&>(*snodes[k]);
          const int* sdof = snode.Dofs();

          GEN_DATA::increaseCapacity(dd_wgap_sl_var);
          dd_wgap_sl_var.repetitive_access(sdof[d], my::gp_id_) += sval(k, 0) * val;
        }
      }
    }

    /*------------------------------------------------------------------------*/
    // master contributions
    Deriv2ndMap& dd_wgap_ma = cnode.AugData().GetDeriv2nd_WGapMa();

    for (unsigned k = 0; k < my::MASTERNUMNODE; ++k)
    {
      const CoNode& mnode = static_cast<const CoNode&>(*mnodes[k]);
      const int* mdof = mnode.Dofs();

      for (unsigned d = 0; d < probdim; ++d)
      {
        GEN_DATA::increaseCapacity(dd_wgap_ma);
        Deriv1stMap& dd_wgap_ma_var = dd_wgap_ma.repetitive_access(mdof[d], my::gp_id_);

        // variation of the master position multiplied with the linearized
        // smooth normal
        for (auto& d_n_unit_lin : d_n_unit[d])
        {
          GEN_DATA::increaseCapacity(dd_wgap_ma_var);
          dd_wgap_ma_var.repetitive_access(d_n_unit_lin.first, my::gp_id_) +=
              tmp * mval(k, 0) * d_n_unit_lin.second;
        }

        // 2-nd order derivative of the master position multiplied by the smooth
        // normal
        for (unsigned j = 0; j < my::MASTERDIM; ++j)
        {
          const double val = tmp * mderiv(j, k) * gpn[d];
          //          if ( val == 0 )
          //            continue;

          for (auto& d_mxigp_j_lin : d_mxigp[j])
          {
            GEN_DATA::increaseCapacity(dd_wgap_ma_var);
            dd_wgap_ma_var.repetitive_access(d_mxigp_j_lin.first, my::gp_id_) +=
                val * d_mxigp_j_lin.second;
          }
        }
      }
    }

    // varied smooth normal multiplied with the linearized master position
    for (unsigned d = 0; d < probdim; ++d)
    {
      for (auto& d_n_unit_var : d_n_unit[d])
      {
        GEN_DATA::increaseCapacity(dd_wgap_ma);
        Deriv1stMap& dd_wgap_ma_var = dd_wgap_ma.repetitive_access(d_n_unit_var.first, my::gp_id_);

        const double val = tmp * d_n_unit_var.second;

        // variation of the master position multiplied with the linearized
        // smooth normal ( 1-st part )
        for (unsigned k = 0; k < my::MASTERNUMNODE; ++k)
        {
          const CoNode& mnode = static_cast<const CoNode&>(*mnodes[k]);
          const int* mdof = mnode.Dofs();

          GEN_DATA::increaseCapacity(dd_wgap_ma_var);
          dd_wgap_ma_var.repetitive_access(mdof[d], my::gp_id_) += val * mval(k, 0);
        }

        /* linearization of the master convective parametric coordinate
         * multiplied by the master convective base vectors multiplied with
         * the varied smooth unit normal */
        for (unsigned j = 0; j < my::MASTERDIM; ++j)
        {
          const double val_1 = val * mtau(d, j);

          for (auto& d_mxigp_j_lin : d_mxigp[j])
          {
            GEN_DATA::increaseCapacity(dd_wgap_ma_var);
            dd_wgap_ma_var.repetitive_access(d_mxigp_j_lin.first, my::gp_id_) +=
                d_mxigp_j_lin.second * val_1;
          }
        }
      }
    }

    /*------------------------------------------------------------------------*/
    // 2-nd order derivative of the smooth unit normal multiplied by the
    // slave and master position, respectively.
#ifndef SWITCH_2ND_DERIV_NORMAL_OFF
    for (unsigned d = 0; d < probdim; ++d)
    {
      const double val_sl = tmp * xs(d, 0);
      const double val_ma = tmp * xm(d, 0);

      for (auto& dd_n_unit_var : dd_n_unit[d])
      {
        GEN_DATA::increaseCapacity(dd_wgap_sl);
        Deriv1stMap& dd_wgap_sl_var = dd_wgap_sl.repetitive_access(dd_n_unit_var.first, my::gp_id_);

        GEN_DATA::increaseCapacity(dd_wgap_ma);
        Deriv1stMap& dd_wgap_ma_var = dd_wgap_ma.repetitive_access(dd_n_unit_var.first, my::gp_id_);

        for (auto& dd_n_unit_var_lin : dd_n_unit_var.second)
        {
          const int gid_lin = dd_n_unit_var_lin.first;

          GEN_DATA::increaseCapacity(dd_wgap_sl_var);
          dd_wgap_sl_var.repetitive_access(gid_lin, my::gp_id_) +=
              val_sl * dd_n_unit_var_lin.second;

          GEN_DATA::increaseCapacity(dd_wgap_ma_var);
          dd_wgap_ma_var.repetitive_access(gid_lin, my::gp_id_) +=
              val_ma * dd_n_unit_var_lin.second;
        }
      }
    }
#endif
  }

  this->timer_.stop(TimeID::INCOMPLETE_Add_Jac_Deriv2nd_GapN);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::CompleteIntPolicy<probdim, slavetype, mastertype>::Get_Deriv2nd_WGap(
    MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
    const LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
    const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
    const LINALG::Matrix<3, my::MASTERNUMNODE>& mderiv2nd, const LINALG::Matrix<3, 2>& mtau,
    const double* gpn, const double wgt, const double gapn_sl, const double gapn_ma,
    const double jac, const Deriv1stMap& d_jac, const Deriv2ndMap& dd_jac,
    const Deriv1stVecMap& d_mxigp, const Deriv2ndVecMap& dd_mxigp, const Deriv1stVecMap& d_n_unit,
    const Deriv2ndVecMap& dd_n_unit, const Deriv1stMap& d_gapn_sl,
    const Deriv1stMap& d_gapn_ma) const
{
  // evaluate the complete second order derivative of the discrete normal gap
  Add_Jac_Deriv2nd_GapN(sele, mele, sval, mval, lmval, mderiv, mderiv2nd, mtau, gpn, wgt, jac,
      d_mxigp, dd_mxigp, d_n_unit, dd_n_unit);

  this->Add_Deriv1st_GapN_Deriv1st_Jac(sele, lmval, wgt, d_gapn_sl, d_gapn_ma, d_jac);

  this->Add_GapN_Deriv2nd_Jac(sele, lmval, wgt, gapn_sl, gapn_ma, dd_jac);

  DRT::Node* const* snodes = sele.Nodes();
  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    CoNode& cnode = static_cast<CoNode&>(*snodes[i]);

    // slave contributions
    cnode.AugData().GetDeriv2nd_WGapSl().complete();
    // master contributions
    cnode.AugData().GetDeriv2nd_WGapMa().complete();
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::BaseIntPolicy<probdim, slavetype, mastertype>::Add_GapN_Deriv2nd_Jac(
    MORTAR::MortarElement& sele, const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
    const double gapn_sl, const double gapn_ma, const Deriv2ndMap& dd_jac) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    CoNode& cnode = static_cast<CoNode&>(*snodes[i]);

    const double tmp_sl_ma = wgt * lmval(i, 0) * (gapn_sl - gapn_ma);

    /*------------------------------------------------------------------------*/
    // slave + master contributions
    Deriv2ndMap& dd_wgap_sl = cnode.AugData().GetDeriv2nd_WGapSl();

    for (auto& dd_jac_var : dd_jac)
    {
      const int gid_var = dd_jac_var.first;
      Deriv1stMap& dd_wgap_sl_var = dd_wgap_sl.repetitive_access(gid_var, my::gp_id_);

      for (auto& dd_jac_var_lin : dd_jac_var.second)
      {
        const int gid_lin = dd_jac_var_lin.first;
        dd_wgap_sl_var.repetitive_access(gid_lin, my::gp_id_) += tmp_sl_ma * dd_jac_var_lin.second;
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::BaseIntPolicy<probdim, slavetype, mastertype>::Add_Var_GapN_Lin_Jac(
    MORTAR::MortarElement& sele, const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
    const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma, const Deriv1stMap& d_jac) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    CoNode& cnode = static_cast<CoNode&>(*snodes[i]);

    const double tmp = wgt * lmval(i, 0);

    /*------------------------------------------------------------------------*/
    // slave contributions
    Deriv2ndMap& dd_wgap_sl = cnode.AugData().GetDeriv2nd_WGapSl();

    /* varied slave part of the normal gap multiplied with the linearized
     * jacobian */
    for (auto& d_gapn_sl_var : d_gapn_sl)
    {
      GEN_DATA::increaseCapacity(dd_wgap_sl);
      Deriv1stMap& dd_wgap_sl_var = dd_wgap_sl.repetitive_access(d_gapn_sl_var.first, my::gp_id_);

      const double val = tmp * d_gapn_sl_var.second;

      for (auto& d_jac_lin : d_jac)
      {
        GEN_DATA::increaseCapacity(dd_wgap_sl_var);
        dd_wgap_sl_var.repetitive_access(d_jac_lin.first, my::gp_id_) += d_jac_lin.second * val;
      }
    }

    /*------------------------------------------------------------------------*/
    // master contributions
    Deriv2ndMap& dd_wgap_ma = cnode.AugData().GetDeriv2nd_WGapMa();

    /* varied master part of the normal gap multiplied with the linearized
     * jacobian */
    for (auto& d_gapn_ma_var : d_gapn_ma)
    {
      GEN_DATA::increaseCapacity(dd_wgap_ma);
      Deriv1stMap& dd_wgap_ma_var = dd_wgap_ma.repetitive_access(d_gapn_ma_var.first, my::gp_id_);

      const double val = tmp * d_gapn_ma_var.second;

      for (auto& d_jac_lin : d_jac)
      {
        GEN_DATA::increaseCapacity(dd_wgap_ma_var);
        dd_wgap_ma_var.repetitive_access(d_jac_lin.first, my::gp_id_) += d_jac_lin.second * val;
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::BaseIntPolicy<probdim, slavetype, mastertype>::Add_Var_Jac_Lin_GapN(
    MORTAR::MortarElement& sele, const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
    const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma, const Deriv1stMap& d_jac) const
{
  DRT::Node* const* snodes = sele.Nodes();

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    CoNode& cnode = static_cast<CoNode&>(*snodes[i]);

    const double tmp = wgt * lmval(i, 0);

    // slave contributions
    Deriv2ndMap& dd_wgap_sl = cnode.AugData().GetDeriv2nd_WGapSl();
    // master contributions
    Deriv2ndMap& dd_wgap_ma = cnode.AugData().GetDeriv2nd_WGapMa();

    /*------------------------------------------------------------------------*/
    // slave + master contributions
    /* varied jacobian */
    for (auto& d_jac_var : d_jac)
    {
      const double val = tmp * d_jac_var.second;

      // linearized slave part of the normal gap
      GEN_DATA::increaseCapacity(dd_wgap_sl);
      Deriv1stMap& dd_wgap_sl_var = dd_wgap_sl.repetitive_access(d_jac_var.first, my::gp_id_);
      for (auto& d_gapn_sl_lin : d_gapn_sl)
      {
        GEN_DATA::increaseCapacity(dd_wgap_sl_var);
        dd_wgap_sl_var.repetitive_access(d_gapn_sl_lin.first, my::gp_id_) +=
            d_gapn_sl_lin.second * val;
      }

      // linearized master part of the normal gap
      Deriv1stMap& dd_wgap_ma_var = dd_wgap_ma.repetitive_access(d_jac_var.first, my::gp_id_);
      for (auto& d_gapn_ma_lin : d_gapn_ma)
      {
        GEN_DATA::increaseCapacity(dd_wgap_ma_var);
        dd_wgap_ma_var.repetitive_access(d_gapn_ma_lin.first, my::gp_id_) +=
            d_gapn_ma_lin.second * val;
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::IncompleteIntPolicy<probdim, slavetype,
    mastertype>::Add_Deriv1st_GapN_Deriv1st_Jac(MORTAR::MortarElement& sele,
    const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
    const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma, const Deriv1stMap& d_jac) const
{
  this->timer_.start(TimeID::INCOMPLETE_Add_Deriv1st_GapN_Deriv1st_Jac);
  this->Add_Var_GapN_Lin_Jac(sele, lmval, wgt, d_gapn_sl, d_gapn_ma, d_jac);
  this->timer_.stop(TimeID::INCOMPLETE_Add_Deriv1st_GapN_Deriv1st_Jac);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::CompleteIntPolicy<probdim, slavetype,
    mastertype>::Add_Deriv1st_GapN_Deriv1st_Jac(MORTAR::MortarElement& sele,
    const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
    const Deriv1stMap& d_gapn_sl, const Deriv1stMap& d_gapn_ma, const Deriv1stMap& d_jac) const
{
  this->Add_Var_GapN_Lin_Jac(sele, lmval, wgt, d_gapn_sl, d_gapn_ma, d_jac);
  this->Add_Var_Jac_Lin_GapN(sele, lmval, wgt, d_gapn_sl, d_gapn_ma, d_jac);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::CompleteIntPolicy<probdim, slavetype, mastertype>::Add_Jac_Deriv2nd_GapN(
    MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
    const LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
    const LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
    const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
    const LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
    const LINALG::Matrix<3, my::MASTERNUMNODE>& mderiv2nd, const LINALG::Matrix<3, 2>& mtau,
    const double* gpn, const double wgt, const double jac, const Deriv1stVecMap& d_mxigp,
    const Deriv2ndVecMap& dd_mxigp, const Deriv1stVecMap& d_n_unit,
    const Deriv2ndVecMap& dd_n_unit) const
{
  DRT::Node* const* snodes = sele.Nodes();
  const DRT::Node* const* mnodes = mele.Nodes();

  static const unsigned dd_2nd_indices[2][2] = {{0, 2}, {2, 1}};

  LINALG::Matrix<3, my::SLAVENUMNODE> scoord;
  sele.GetNodalCoords(scoord);

  LINALG::Matrix<3, 1> xs(true);
  xs.MultiplyNN(scoord, sval);

  LINALG::Matrix<3, my::MASTERNUMNODE> mcoord;
  mele.GetNodalCoords(mcoord);

  LINALG::Matrix<3, 1> xm(true);
  xm.MultiplyNN(mcoord, mval);

  for (unsigned i = 0; i < my::SLAVENUMNODE; ++i)
  {
    CoNode& cnode = static_cast<CoNode&>(*snodes[i]);

    // constant scaling term
    const double tmp = wgt * jac * lmval(i, 0);

    /*------------------------------------------------------------------------*/
    // slave contributions
    Deriv2ndMap& dd_wgap_sl = cnode.AugData().GetDeriv2nd_WGapSl();

    for (unsigned k = 0; k < my::SLAVENUMNODE; ++k)
    {
      const CoNode& snode = static_cast<const CoNode&>(*snodes[k]);
      const int* sdof = snode.Dofs();

      // variation of the slave position multiplied with the linearized
      // smooth normal
      for (unsigned d = 0; d < probdim; ++d)
      {
        GEN_DATA::increaseCapacity(dd_wgap_sl);
        Deriv1stMap& dd_wgap_sl_var = dd_wgap_sl.repetitive_access(sdof[d], my::gp_id_);

        for (auto& d_n_unit_lin : d_n_unit[d])
        {
          GEN_DATA::increaseCapacity(dd_wgap_sl_var);
          dd_wgap_sl_var.repetitive_access(d_n_unit_lin.first, my::gp_id_) +=
              tmp * sval(k, 0) * d_n_unit_lin.second;
        }
      }
    }

    // varied smooth normal multiplied with the linearized slave position
    for (unsigned d = 0; d < probdim; ++d)
    {
      for (auto& d_n_unit_var : d_n_unit[d])
      {
        GEN_DATA::increaseCapacity(dd_wgap_sl);
        Deriv1stMap& dd_wgap_sl_var = dd_wgap_sl.repetitive_access(d_n_unit_var.first, my::gp_id_);

        const double val = tmp * d_n_unit_var.second;

        for (unsigned k = 0; k < my::SLAVENUMNODE; ++k)
        {
          const CoNode& snode = static_cast<const CoNode&>(*snodes[k]);
          const int* sdof = snode.Dofs();

          GEN_DATA::increaseCapacity(dd_wgap_sl_var);
          dd_wgap_sl_var.repetitive_access(sdof[d], my::gp_id_) += sval(k, 0) * val;
        }
      }
    }

    /*------------------------------------------------------------------------*/
    // master contributions
    Deriv2ndMap& dd_wgap_ma = cnode.AugData().GetDeriv2nd_WGapMa();
    for (unsigned k = 0; k < my::MASTERNUMNODE; ++k)
    {
      const CoNode& mnode = static_cast<const CoNode&>(*mnodes[k]);
      const int* mdof = mnode.Dofs();

      for (unsigned d = 0; d < probdim; ++d)
      {
        GEN_DATA::increaseCapacity(dd_wgap_ma);
        Deriv1stMap& dd_wgap_ma_var = dd_wgap_ma.repetitive_access(mdof[d], my::gp_id_);

        // variation of the master position multiplied with the linearized
        // smooth normal
        for (auto& d_n_unit_lin : d_n_unit[d])
        {
          GEN_DATA::increaseCapacity(dd_wgap_ma_var);
          dd_wgap_ma_var.repetitive_access(d_n_unit_lin.first, my::gp_id_) +=
              tmp * mval(k, 0) * d_n_unit_lin.second;
        }

        // 2-nd order derivative of the master position multiplied by the smooth
        // normal
        for (unsigned j = 0; j < my::MASTERDIM; ++j)
        {
          const double val = tmp * mderiv(j, k) * gpn[d];
          //          if ( val == 0 )
          //            continue;

          for (auto& d_mxigp_j_lin : d_mxigp[j])
          {
            GEN_DATA::increaseCapacity(dd_wgap_ma_var);
            dd_wgap_ma_var.repetitive_access(d_mxigp_j_lin.first, my::gp_id_) +=
                val * d_mxigp_j_lin.second;
          }
        }
      }

      for (unsigned j = 0; j < my::MASTERDIM; ++j)
      {
        for (auto& d_mxigp_j_var : d_mxigp[j])
        {
          GEN_DATA::increaseCapacity(dd_wgap_ma);
          Deriv1stMap& dd_wgap_ma_var =
              dd_wgap_ma.repetitive_access(d_mxigp_j_var.first, my::gp_id_);

          /* varied master convective parametric coordinate multiplied by the
           * linearized master convective base vectors and the smooth unit
           * normal */
          // 1-st part
          const double val = mderiv(j, k) * tmp * d_mxigp_j_var.second;
          for (unsigned d = 0; d < probdim; ++d)
          {
            //            if ( val == 0)
            //              continue;

            GEN_DATA::increaseCapacity(dd_wgap_ma_var);
            dd_wgap_ma_var.repetitive_access(mdof[d], my::gp_id_) += gpn[d] * val;
          }

          // 2-nd part
          const double val1 = tmp * d_mxigp_j_var.second;
          //          if ( val1 == 0 )
          //            continue;

          for (unsigned l = 0; l < my::MASTERDIM; ++l)
          {
            // index for the 2-nd derivative:
            //        N_{,11}    | N_{,12}    | N_{,21}    | N_{,22}
            // (j,l): (0,0) -> 0 | (0,1) -> 2 | (1,0) -> 2 | (1,1) -> 1
            const unsigned j_2nd = dd_2nd_indices[j][l];

            for (auto& d_mxigp_l_lin : d_mxigp[l])
            {
              GEN_DATA::increaseCapacity(dd_wgap_ma_var);
              double& dd_wgap_ma_var_lin =
                  dd_wgap_ma_var.repetitive_access(d_mxigp_l_lin.first, my::gp_id_);

              const double val2 = mderiv2nd(j_2nd, k) * d_mxigp_l_lin.second * val1;

              //              if ( val2 == 0)
              //                continue;

              for (unsigned d = 0; d < probdim; ++d)
              {
                dd_wgap_ma_var_lin += mcoord(d, k) * gpn[d] * val2;
              }
            }
          }
        }
      }
    }

    /* 2-nd derivative of the master parametric coordinates multiplied by
     * the convective base vectors and the smooth unit normal */
    for (unsigned j = 0; j < my::MASTERDIM; ++j)
    {
      for (auto& dd_mxigp_j_var : dd_mxigp[j])
      {
        GEN_DATA::increaseCapacity(dd_wgap_ma);
        Deriv1stMap& dd_wgap_ma_var =
            dd_wgap_ma.repetitive_access(dd_mxigp_j_var.first, my::gp_id_);

        for (auto& dd_mxigp_j_var_lin : dd_mxigp_j_var.second)
        {
          GEN_DATA::increaseCapacity(dd_wgap_ma_var);
          double& dd_wgap_ma_var_lin =
              dd_wgap_ma_var.repetitive_access(dd_mxigp_j_var_lin.first, my::gp_id_);

          const double val = dd_mxigp_j_var_lin.second * tmp;
          //          if ( val == 0 )
          //            continue;

          for (unsigned d = 0; d < probdim; ++d)
          {
            dd_wgap_ma_var_lin += gpn[d] * mtau(d, j) * val;
          }
        }
      }
    }

    /* variation of the master convective parametric coordinate
     * multiplied by the master convective base vectors multiplied with
     * the linearized smooth unit normal */
    for (unsigned j = 0; j < my::MASTERDIM; ++j)
    {
      for (auto& d_mxigp_j_var : d_mxigp[j])
      {
        GEN_DATA::increaseCapacity(dd_wgap_ma);
        Deriv1stMap& dd_wgap_ma_var = dd_wgap_ma.repetitive_access(d_mxigp_j_var.first, my::gp_id_);

        for (unsigned d = 0; d < probdim; ++d)
        {
          const double val = mtau(d, j) * d_mxigp_j_var.second * tmp;

          for (auto& d_n_unit_lin : d_n_unit[d])
          {
            GEN_DATA::increaseCapacity(dd_wgap_ma_var);
            dd_wgap_ma_var.repetitive_access(d_n_unit_lin.first, my::gp_id_) +=
                val * d_n_unit_lin.second;
          }
        }
      }
    }

    // varied smooth normal multiplied with the linearized master position
    for (unsigned d = 0; d < probdim; ++d)
    {
      for (auto& d_n_unit_var : d_n_unit[d])
      {
        GEN_DATA::increaseCapacity(dd_wgap_ma);
        Deriv1stMap& dd_wgap_ma_var = dd_wgap_ma.repetitive_access(d_n_unit_var.first, my::gp_id_);

        const double val = tmp * d_n_unit_var.second;

        // variation of the master position multiplied with the linearized
        // smooth normal ( 1-st part )
        for (unsigned k = 0; k < my::MASTERNUMNODE; ++k)
        {
          const CoNode& mnode = static_cast<const CoNode&>(*mnodes[k]);
          const int* mdof = mnode.Dofs();

          GEN_DATA::increaseCapacity(dd_wgap_ma_var);
          dd_wgap_ma_var.repetitive_access(mdof[d], my::gp_id_) += val * mval(k, 0);
        }

        /* linearization of the master convective parametric coordinate
         * multiplied by the master convective base vectors multiplied with
         * the varied smooth unit normal */
        for (unsigned j = 0; j < my::MASTERDIM; ++j)
        {
          const double val_1 = val * mtau(d, j);

          for (auto& d_mxigp_j_lin : d_mxigp[j])
          {
            GEN_DATA::increaseCapacity(dd_wgap_ma_var);
            dd_wgap_ma_var.repetitive_access(d_mxigp_j_lin.first, my::gp_id_) +=
                d_mxigp_j_lin.second * val_1;
          }
        }
      }
    }

    /*------------------------------------------------------------------------*/
    // 2-nd order derivative of the smooth unit normal multiplied by the
    // slave and master position, respectively.
    for (unsigned d = 0; d < probdim; ++d)
    {
      const double val_sl_ma = tmp * (xs(d, 0) - xm(d, 0));

      for (auto& dd_n_unit_var : dd_n_unit[d])
      {
        GEN_DATA::increaseCapacity(dd_wgap_sl);
        Deriv1stMap& dd_wgap_sl_var = dd_wgap_sl.repetitive_access(dd_n_unit_var.first, my::gp_id_);

        for (auto& dd_n_unit_var_lin : dd_n_unit_var.second)
        {
          const int gid_lin = dd_n_unit_var_lin.first;

          GEN_DATA::increaseCapacity(dd_wgap_sl_var);
          dd_wgap_sl_var.repetitive_access(gid_lin, my::gp_id_) +=
              val_sl_ma * dd_n_unit_var_lin.second;
        }
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::IncompleteIntPolicy<probdim, slavetype, mastertype>::Get_Deriv1st_WGapN_Error(
    const MORTAR::MortarElement& sele, const std::vector<unsigned>& active_nlids,
    const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double* gpn, const double gapn_sl,
    const double gapn_ma, const double wgt, const double jac, const Deriv1stMap& d_jac,
    const LINALG::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp, Deriv1stMap& d_gapn_ma,
    std::unordered_map<int, Deriv1stMap>& error_ma,
    std::unordered_map<int, Deriv1stMap>& error_jac) const
{
  const int* sl_node_gids = sele.NodeIds();

  // (0) derivative of the master convective parameter coordinate
  for (unsigned j = 0; j < my::MASTERDIM; ++j)
  {
    for (auto& d_mxigp_j : d_mxigp[j])
    {
      double& deriv_gapn_ma_j = d_gapn_ma(d_mxigp_j.first);

      for (unsigned d = 0; d < probdim; ++d)
      {
        deriv_gapn_ma_j += gpn[d] * mtau(d, j) * d_mxigp_j.second;
      }
    }
  }

  d_gapn_ma.complete();

  for (const unsigned snlid : active_nlids)
  {
    const int sngid = sl_node_gids[snlid];

    // --- error due to the missing variation of the master parametric coordinate
    Deriv1stMap& error_ma_i = error_ma[sngid];
    const double tmp_ma = wgt * lmval(snlid, 0) * jac;

    for (const auto& d_gap_ma_j : d_gapn_ma)
    {
      GEN_DATA::increaseCapacity(error_ma_i);
      error_ma_i[d_gap_ma_j.first] += tmp_ma * d_gap_ma_j.second;
    }

    // --- error due to the missing variation of the jacobian -----------------
    const double tmp_jac = lmval(snlid, 0) * wgt * (gapn_ma - gapn_sl);

    // slave contributions
    Deriv1stMap& error_jac_i = error_jac[sngid];

    for (auto& d_jac_var : d_jac)
    {
      GEN_DATA::increaseCapacity(error_jac_i);
      error_jac_i(d_jac_var.first) += tmp_jac * d_jac_var.second;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned probdim, DRT::Element::DiscretizationType slavetype,
    DRT::Element::DiscretizationType mastertype>
void CONTACT::AUG::CompleteIntPolicy<probdim, slavetype, mastertype>::Get_Deriv1st_WGapN_Error(
    const MORTAR::MortarElement& sele, const std::vector<unsigned>& active_nlids,
    const LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double* gpn, const double gapn_sl,
    const double gapn_ma, const double wgt, const double jacslave, const Deriv1stMap& d_jac,
    const LINALG::Matrix<3, 2>& mtau, const Deriv1stVecMap& d_mxigp, Deriv1stMap& d_gapn_ma,
    std::unordered_map<int, Deriv1stMap>& error_ma,
    std::unordered_map<int, Deriv1stMap>& error_jac) const
{
  // for the complete policy the error is zero
  IO::cout << "LINE " << __LINE__ << " -- " << __FUNCTION__
           << ": "
              "There is no error for the complete variational approach."
           << IO::endl;
  error_ma.clear();
  error_jac.clear();
}

/*----------------------------------------------------------------------------*/
template class CONTACT::AUG::BaseSlaveIntPolicy<2, DRT::Element::line2>;

template class CONTACT::AUG::BaseSlaveIntPolicy<2, DRT::Element::nurbs2>;

template class CONTACT::AUG::BaseSlaveIntPolicy<2, DRT::Element::nurbs3>;

template class CONTACT::AUG::BaseSlaveIntPolicy<3, DRT::Element::quad4>;

template class CONTACT::AUG::BaseSlaveIntPolicy<3, DRT::Element::tri3>;

template class CONTACT::AUG::BaseSlaveIntPolicy<3, DRT::Element::nurbs4>;

template class CONTACT::AUG::BaseSlaveIntPolicy<3, DRT::Element::nurbs9>;

/*----------------------------------------------------------------------------*/
template class CONTACT::AUG::BaseIntPolicy<2, DRT::Element::line2, DRT::Element::line2>;
template class CONTACT::AUG::BaseIntPolicy<2, DRT::Element::line2, DRT::Element::nurbs2>;
template class CONTACT::AUG::BaseIntPolicy<2, DRT::Element::line2, DRT::Element::nurbs3>;

template class CONTACT::AUG::BaseIntPolicy<2, DRT::Element::nurbs2, DRT::Element::nurbs2>;
template class CONTACT::AUG::BaseIntPolicy<2, DRT::Element::nurbs2, DRT::Element::line2>;
template class CONTACT::AUG::BaseIntPolicy<2, DRT::Element::nurbs2, DRT::Element::nurbs3>;

template class CONTACT::AUG::BaseIntPolicy<2, DRT::Element::nurbs3, DRT::Element::nurbs3>;
template class CONTACT::AUG::BaseIntPolicy<2, DRT::Element::nurbs3, DRT::Element::line2>;
template class CONTACT::AUG::BaseIntPolicy<2, DRT::Element::nurbs3, DRT::Element::nurbs2>;

template class CONTACT::AUG::BaseIntPolicy<3, DRT::Element::quad4, DRT::Element::quad4>;
template class CONTACT::AUG::BaseIntPolicy<3, DRT::Element::quad4, DRT::Element::tri3>;
template class CONTACT::AUG::BaseIntPolicy<3, DRT::Element::quad4, DRT::Element::nurbs4>;
template class CONTACT::AUG::BaseIntPolicy<3, DRT::Element::quad4, DRT::Element::nurbs9>;

template class CONTACT::AUG::BaseIntPolicy<3, DRT::Element::tri3, DRT::Element::quad4>;
template class CONTACT::AUG::BaseIntPolicy<3, DRT::Element::tri3, DRT::Element::tri3>;
template class CONTACT::AUG::BaseIntPolicy<3, DRT::Element::tri3, DRT::Element::nurbs4>;
template class CONTACT::AUG::BaseIntPolicy<3, DRT::Element::tri3, DRT::Element::nurbs9>;

template class CONTACT::AUG::BaseIntPolicy<3, DRT::Element::nurbs4, DRT::Element::nurbs4>;
template class CONTACT::AUG::BaseIntPolicy<3, DRT::Element::nurbs4, DRT::Element::quad4>;
template class CONTACT::AUG::BaseIntPolicy<3, DRT::Element::nurbs4, DRT::Element::tri3>;
template class CONTACT::AUG::BaseIntPolicy<3, DRT::Element::nurbs4, DRT::Element::nurbs9>;

template class CONTACT::AUG::BaseIntPolicy<3, DRT::Element::nurbs9, DRT::Element::nurbs9>;
template class CONTACT::AUG::BaseIntPolicy<3, DRT::Element::nurbs9, DRT::Element::quad4>;
template class CONTACT::AUG::BaseIntPolicy<3, DRT::Element::nurbs9, DRT::Element::tri3>;
template class CONTACT::AUG::BaseIntPolicy<3, DRT::Element::nurbs9, DRT::Element::nurbs4>;

/*----------------------------------------------------------------------------*/
template class CONTACT::AUG::IncompleteIntPolicy<2, DRT::Element::line2, DRT::Element::line2>;
template class CONTACT::AUG::IncompleteIntPolicy<2, DRT::Element::line2, DRT::Element::nurbs2>;
template class CONTACT::AUG::IncompleteIntPolicy<2, DRT::Element::line2, DRT::Element::nurbs3>;

template class CONTACT::AUG::IncompleteIntPolicy<2, DRT::Element::nurbs2, DRT::Element::nurbs2>;
template class CONTACT::AUG::IncompleteIntPolicy<2, DRT::Element::nurbs2, DRT::Element::line2>;
template class CONTACT::AUG::IncompleteIntPolicy<2, DRT::Element::nurbs2, DRT::Element::nurbs3>;

template class CONTACT::AUG::IncompleteIntPolicy<2, DRT::Element::nurbs3, DRT::Element::nurbs3>;
template class CONTACT::AUG::IncompleteIntPolicy<2, DRT::Element::nurbs3, DRT::Element::line2>;
template class CONTACT::AUG::IncompleteIntPolicy<2, DRT::Element::nurbs3, DRT::Element::nurbs2>;

template class CONTACT::AUG::IncompleteIntPolicy<3, DRT::Element::quad4, DRT::Element::quad4>;
template class CONTACT::AUG::IncompleteIntPolicy<3, DRT::Element::quad4, DRT::Element::tri3>;
template class CONTACT::AUG::IncompleteIntPolicy<3, DRT::Element::quad4, DRT::Element::nurbs4>;
template class CONTACT::AUG::IncompleteIntPolicy<3, DRT::Element::quad4, DRT::Element::nurbs9>;

template class CONTACT::AUG::IncompleteIntPolicy<3, DRT::Element::tri3, DRT::Element::quad4>;
template class CONTACT::AUG::IncompleteIntPolicy<3, DRT::Element::tri3, DRT::Element::tri3>;
template class CONTACT::AUG::IncompleteIntPolicy<3, DRT::Element::tri3, DRT::Element::nurbs4>;
template class CONTACT::AUG::IncompleteIntPolicy<3, DRT::Element::tri3, DRT::Element::nurbs9>;

template class CONTACT::AUG::IncompleteIntPolicy<3, DRT::Element::nurbs4, DRT::Element::nurbs4>;
template class CONTACT::AUG::IncompleteIntPolicy<3, DRT::Element::nurbs4, DRT::Element::quad4>;
template class CONTACT::AUG::IncompleteIntPolicy<3, DRT::Element::nurbs4, DRT::Element::tri3>;
template class CONTACT::AUG::IncompleteIntPolicy<3, DRT::Element::nurbs4, DRT::Element::nurbs9>;

template class CONTACT::AUG::IncompleteIntPolicy<3, DRT::Element::nurbs9, DRT::Element::nurbs9>;
template class CONTACT::AUG::IncompleteIntPolicy<3, DRT::Element::nurbs9, DRT::Element::quad4>;
template class CONTACT::AUG::IncompleteIntPolicy<3, DRT::Element::nurbs9, DRT::Element::tri3>;
template class CONTACT::AUG::IncompleteIntPolicy<3, DRT::Element::nurbs9, DRT::Element::nurbs4>;

/*----------------------------------------------------------------------------*/
template class CONTACT::AUG::CompleteIntPolicy<2, DRT::Element::line2, DRT::Element::line2>;
template class CONTACT::AUG::CompleteIntPolicy<2, DRT::Element::line2, DRT::Element::nurbs2>;
template class CONTACT::AUG::CompleteIntPolicy<2, DRT::Element::line2, DRT::Element::nurbs3>;

template class CONTACT::AUG::CompleteIntPolicy<2, DRT::Element::nurbs2, DRT::Element::nurbs2>;
template class CONTACT::AUG::CompleteIntPolicy<2, DRT::Element::nurbs2, DRT::Element::line2>;
template class CONTACT::AUG::CompleteIntPolicy<2, DRT::Element::nurbs2, DRT::Element::nurbs3>;

template class CONTACT::AUG::CompleteIntPolicy<2, DRT::Element::nurbs3, DRT::Element::nurbs3>;
template class CONTACT::AUG::CompleteIntPolicy<2, DRT::Element::nurbs3, DRT::Element::line2>;
template class CONTACT::AUG::CompleteIntPolicy<2, DRT::Element::nurbs3, DRT::Element::nurbs2>;

template class CONTACT::AUG::CompleteIntPolicy<3, DRT::Element::quad4, DRT::Element::quad4>;
template class CONTACT::AUG::CompleteIntPolicy<3, DRT::Element::quad4, DRT::Element::tri3>;
template class CONTACT::AUG::CompleteIntPolicy<3, DRT::Element::quad4, DRT::Element::nurbs4>;
template class CONTACT::AUG::CompleteIntPolicy<3, DRT::Element::quad4, DRT::Element::nurbs9>;

template class CONTACT::AUG::CompleteIntPolicy<3, DRT::Element::tri3, DRT::Element::quad4>;
template class CONTACT::AUG::CompleteIntPolicy<3, DRT::Element::tri3, DRT::Element::tri3>;
template class CONTACT::AUG::CompleteIntPolicy<3, DRT::Element::tri3, DRT::Element::nurbs4>;
template class CONTACT::AUG::CompleteIntPolicy<3, DRT::Element::tri3, DRT::Element::nurbs9>;

template class CONTACT::AUG::CompleteIntPolicy<3, DRT::Element::nurbs4, DRT::Element::nurbs4>;
template class CONTACT::AUG::CompleteIntPolicy<3, DRT::Element::nurbs4, DRT::Element::quad4>;
template class CONTACT::AUG::CompleteIntPolicy<3, DRT::Element::nurbs4, DRT::Element::tri3>;
template class CONTACT::AUG::CompleteIntPolicy<3, DRT::Element::nurbs4, DRT::Element::nurbs9>;

template class CONTACT::AUG::CompleteIntPolicy<3, DRT::Element::nurbs9, DRT::Element::nurbs9>;
template class CONTACT::AUG::CompleteIntPolicy<3, DRT::Element::nurbs9, DRT::Element::quad4>;
template class CONTACT::AUG::CompleteIntPolicy<3, DRT::Element::nurbs9, DRT::Element::tri3>;
template class CONTACT::AUG::CompleteIntPolicy<3, DRT::Element::nurbs9, DRT::Element::nurbs4>;

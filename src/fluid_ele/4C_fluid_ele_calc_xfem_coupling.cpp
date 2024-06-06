/*----------------------------------------------------------------------*/
/*! \file

\brief Factory class for providing an implementation for coupling with
       Mixed/Stress/Hybrid methods or Nitsche's method to enforce
       interface conditions in the XFEM weakly

\level 2


*/
/*----------------------------------------------------------------------*/


#include "4C_fluid_ele_calc_xfem_coupling.hpp"

#include "4C_bele_bele3.hpp"
#include "4C_fluid_ele_calc_xfem_coupling_impl.hpp"

FOUR_C_NAMESPACE_OPEN

using namespace Discret::ELEMENTS::XFLUID;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void SlaveElementInterface<distype>::DefineStateNames(
    Core::FE::CellType slave_distype,  ///< coupling slave discretization type
    std::string& disp_statename,       ///< name of displacement state at current step
    std::string& vel_statename,        ///< name of velocity state at current step
    std::string& veln_statename        ///< name of velocity state at previous step
)
{
  switch (slave_distype)
  {
      //  case Core::FE::CellType::tri3:
      //  case Core::FE::CellType::tri6:
    case Core::FE::CellType::quad4:
    case Core::FE::CellType::quad8:
    case Core::FE::CellType::quad9:
    {
      disp_statename = std::string("idispnp");
      vel_statename = std::string("ivelnp");
      veln_statename = std::string("iveln");
      break;
    }
    case Core::FE::CellType::hex8:
    case Core::FE::CellType::hex20:
    case Core::FE::CellType::hex27:
    {
      disp_statename = std::string("dispnp");
      vel_statename = std::string("velaf");
      veln_statename = std::string("veln");
      break;
    }
    default:
      FOUR_C_THROW("Unsupported element shape %d", slave_distype);
      break;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Teuchos::RCP<SlaveElementInterface<distype>>
SlaveElementInterface<distype>::create_slave_element_representation(
    Core::Elements::Element* slave_ele,  ///< coupling slave element
    Core::LinAlg::SerialDenseMatrix&
        slave_xyz  ///< global node coordinates of coupling slave element
)
{
  SlaveElementInterface* sla = nullptr;

  // get number of dofs for this slave element
  const unsigned numdofpernode = slave_ele->NumDofPerNode(*slave_ele->Nodes()[0]);

  if (numdofpernode == 3)
  {
    switch (slave_ele->Shape())
    {
        //      case Core::FE::CellType::tri3:
        //      {
        //        typedef
        //        SlaveElementRepresentation<distype,Core::FE::CellType::tri3,3>
        //        SlaveEleType; sla = new SlaveEleType(slave_xyz); break;
        //      }
        //      case Core::FE::CellType::tri6:
        //      {
        //        typedef
        //        SlaveElementRepresentation<distype,Core::FE::CellType::tri6,3>
        //        SlaveEleType; sla = new SlaveEleType(slave_xyz); break;
        //      }
      case Core::FE::CellType::quad4:
      {
        typedef SlaveElementRepresentation<distype, Core::FE::CellType::quad4, 3> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        typedef SlaveElementRepresentation<distype, Core::FE::CellType::quad8, 3> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        typedef SlaveElementRepresentation<distype, Core::FE::CellType::quad9, 3> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported boundary element shape %d", slave_ele->Shape());
        break;
    }
  }
  else if (numdofpernode ==
           4)  // volumetric coupling partners only required for fluid-fluid coupling
  {
    switch (slave_ele->Shape())
    {
        //      case Core::FE::CellType::tri3:
        //      {
        //        typedef
        //        SlaveElementRepresentation<distype,Core::FE::CellType::tri3,4>
        //        SlaveEleType; sla = new SlaveEleType(slave_xyz); break;
        //      }
        //      case Core::FE::CellType::tri6:
        //      {
        //        typedef
        //        SlaveElementRepresentation<distype,Core::FE::CellType::tri6,4>
        //        SlaveEleType; sla = new SlaveEleType(slave_xyz); break;
        //      }
      case Core::FE::CellType::quad4:
      {
        typedef SlaveElementRepresentation<distype, Core::FE::CellType::quad4, 4> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        typedef SlaveElementRepresentation<distype, Core::FE::CellType::quad8, 4> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        typedef SlaveElementRepresentation<distype, Core::FE::CellType::quad9, 4> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case Core::FE::CellType::hex8:
      {
        typedef SlaveElementRepresentation<distype, Core::FE::CellType::hex8, 4> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case Core::FE::CellType::hex20:
      {
        typedef SlaveElementRepresentation<distype, Core::FE::CellType::hex20, 4> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case Core::FE::CellType::hex27:
      {
        typedef SlaveElementRepresentation<distype, Core::FE::CellType::hex27, 4> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      default:
        FOUR_C_THROW("unsupported boundary element shape %d", slave_ele->Shape());
        break;
    }
  }

  return Teuchos::rcp(sla);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Teuchos::RCP<NitscheInterface<distype>>
NitscheInterface<distype>::create_nitsche_coupling_x_fluid_wdbc(
    Core::LinAlg::SerialDenseMatrix::Base& C_umum, Core::LinAlg::SerialDenseMatrix::Base& rhC_um,
    const Discret::ELEMENTS::FluidEleParameterXFEM& fldparaxfem)
{
  NitscheInterface* nit = nullptr;
  typedef NitscheCoupling<distype, Core::FE::CellType::dis_none, 3> NitscheCouplType;
  nit = new NitscheCouplType(C_umum, rhC_um, fldparaxfem);

  return Teuchos::rcp(nit);
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Teuchos::RCP<NitscheInterface<distype>>
NitscheInterface<distype>::create_nitsche_coupling_x_fluid_wdbc(Core::Elements::Element* bele,
    Core::LinAlg::SerialDenseMatrix::Base& bele_xyz, Core::LinAlg::SerialDenseMatrix::Base& C_umum,
    Core::LinAlg::SerialDenseMatrix::Base& rhC_um,
    const Discret::ELEMENTS::FluidEleParameterXFEM& fldparaxfem)
{
  NitscheInterface* nit = nullptr;

  // get number of dofs for this boundary element
  const unsigned numdofpernode = bele->NumDofPerNode(*bele->Nodes()[0]);

  // three dofs per node, for standard Dirichlet coupling, four dofs per node for background
  // geometry coupling
  if (numdofpernode == 3)
  {
    switch (bele->Shape())
    {
      //      case Core::FE::CellType::tri3:
      //      {
      //        typedef NitscheCoupling<distype,Core::FE::CellType::tri3,3>
      //        NitscheCouplType; nit = new NitscheCouplType(bele_xyz,C_umum,rhC_um,); break;
      //      }
      //      case Core::FE::CellType::tri6:
      //      {
      //        typedef NitscheCoupling<distype,Core::FE::CellType::tri6,3>
      //        NitscheCouplType; nit = new NitscheCouplType(bele_xyz,C_umum,rhC_um,); break;
      //      }
      case Core::FE::CellType::quad4:
      {
        typedef NitscheCoupling<distype, Core::FE::CellType::quad4, 3> NitscheCouplType;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        typedef NitscheCoupling<distype, Core::FE::CellType::quad8, 3> NitscheCouplType;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        typedef NitscheCoupling<distype, Core::FE::CellType::quad9, 3> NitscheCouplType;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported boundary element shape %d", bele->Shape());
        break;
    }
  }
  else if (numdofpernode == 4)
  {
    switch (bele->Shape())
    {
      //      case Core::FE::CellType::tri3:
      //      {
      //        typedef NitscheCoupling<distype,Core::FE::CellType::tri3,3>
      //        NitscheCouplType; nit = new NitscheCouplType(bele_xyz,C_umum,rhC_um,); break;
      //      }
      //      case Core::FE::CellType::tri6:
      //      {
      //        typedef NitscheCoupling<distype,Core::FE::CellType::tri6,3>
      //        NitscheCouplType; nit = new NitscheCouplType(bele_xyz,C_umum,rhC_um,); break;
      //      }
      case Core::FE::CellType::quad4:
      {
        typedef NitscheCoupling<distype, Core::FE::CellType::quad4, 4> NitscheCouplType;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        typedef NitscheCoupling<distype, Core::FE::CellType::quad8, 4> NitscheCouplType;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        typedef NitscheCoupling<distype, Core::FE::CellType::quad9, 4> NitscheCouplType;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported boundary element shape %d", bele->Shape());
        break;
    }
  }
  else
    FOUR_C_THROW("Unsupported number of %d nodes for coupling slave element.", numdofpernode);

  return Teuchos::rcp(nit);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Teuchos::RCP<NitscheInterface<distype>>
NitscheInterface<distype>::create_nitsche_coupling_x_fluid_sided(Core::Elements::Element* bele,
    Core::LinAlg::SerialDenseMatrix::Base& bele_xyz, Core::LinAlg::SerialDenseMatrix::Base& C_umum,
    Core::LinAlg::SerialDenseMatrix::Base& C_usum, Core::LinAlg::SerialDenseMatrix::Base& C_umus,
    Core::LinAlg::SerialDenseMatrix::Base& C_usus, Core::LinAlg::SerialDenseMatrix::Base& rhC_um,
    Core::LinAlg::SerialDenseMatrix::Base& rhC_us,
    const Discret::ELEMENTS::FluidEleParameterXFEM& fldparaxfem)
{
  NitscheInterface* nit = nullptr;

  // get number of dofs for this boundary element
  const unsigned numdofpernode = bele->NumDofPerNode(*bele->Nodes()[0]);

  // three dofs per node, for monolithic XFSI, four dofs per node for background-sided fluid-fluid
  // coupling
  if (numdofpernode == 3)
  {
    switch (bele->Shape())
    {
        //    case Core::FE::CellType::tri3:
        //    {
        //      typedef NitscheCoupling<distype,Core::FE::CellType::tri3,3>
        //      NitscheCouplType; nit = new
        //      NitscheCouplType(bele_xyz,C_umum,C_usum,C_umus,C_usus,rhC_um,rhC_us,is_viscAdjointSymmetric);
        //      break;
        //    }
        //    case Core::FE::CellType::tri6:
        //    {
        //      typedef NitscheCoupling<distype,Core::FE::CellType::tri6,3>
        //      NitscheCouplType; nit = new
        //      NitscheCouplType(bele_xyz,C_umum,C_usum,C_umus,C_usus,rhC_um,rhC_us,is_viscAdjointSymmetric);
        //      break;
        //    }
      case Core::FE::CellType::quad4:
      {
        typedef NitscheCoupling<distype, Core::FE::CellType::quad4, 3> NitscheCouplType;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        typedef NitscheCoupling<distype, Core::FE::CellType::quad8, 3> NitscheCouplType;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        typedef NitscheCoupling<distype, Core::FE::CellType::quad9, 3> NitscheCouplType;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported boundary element shape %d", bele->Shape());
        break;
    }
  }
  else if (numdofpernode == 4)
  {
    switch (bele->Shape())
    {
      //    case Core::FE::CellType::tri3:
      //    {
      //      typedef NitscheCoupling<distype,Core::FE::CellType::tri3,4>
      //      NitscheCouplType; nit = new
      //      NitscheCouplType(bele_xyz,C_umum,C_usum,C_umus,C_usus,rhC_um,rhC_us,is_viscAdjointSymmetric);
      //      break;
      //    }
      //    case Core::FE::CellType::tri6:
      //    {
      //      typedef NitscheCoupling<distype,Core::FE::CellType::tri6,4>
      //      NitscheCouplType; nit = new
      //      NitscheCouplType(bele_xyz,C_umum,C_usum,C_umus,C_usus,rhC_um,rhC_us,is_viscAdjointSymmetric);
      //      break;
      //    }
      case Core::FE::CellType::quad4:
      {
        typedef NitscheCoupling<distype, Core::FE::CellType::quad4, 4> NitscheCouplType;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        typedef NitscheCoupling<distype, Core::FE::CellType::quad8, 4> NitscheCouplType;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        typedef NitscheCoupling<distype, Core::FE::CellType::quad9, 4> NitscheCouplType;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported boundary element shape %d", bele->Shape());
        break;
    }
  }
  else
    FOUR_C_THROW("Unsupported number of %d nodes for coupling slave element.", numdofpernode);

  return Teuchos::rcp(nit);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Teuchos::RCP<NitscheInterface<distype>>
NitscheInterface<distype>::create_nitsche_coupling_two_sided(Core::Elements::Element* vele,
    Core::LinAlg::SerialDenseMatrix::Base& vele_xyz, Core::LinAlg::SerialDenseMatrix::Base& C_umum,
    Core::LinAlg::SerialDenseMatrix::Base& C_usum, Core::LinAlg::SerialDenseMatrix::Base& C_umus,
    Core::LinAlg::SerialDenseMatrix::Base& C_usus, Core::LinAlg::SerialDenseMatrix::Base& rhC_um,
    Core::LinAlg::SerialDenseMatrix::Base& rhC_us,
    const Discret::ELEMENTS::FluidEleParameterXFEM& fldparaxfem)
{
  NitscheInterface* nit = nullptr;

  // get number of dofs for the embedded element
  const unsigned numdofpernode = vele->NumDofPerNode(*vele->Nodes()[0]);

  if (numdofpernode == 4)
  {
    switch (vele->Shape())
    {
        //    case Core::FE::CellType::tet4:
        //    {
        //      typedef NitscheCoupling<distype,Core::FE::CellType::tet4,4>
        //      NitscheCouplType; nit = new
        //      NitscheCouplType(vele_xyz,C_umum,C_usum,C_umus,C_usus,rhC_um,rhC_us,is_viscAdjointSymmetric);
        //      break;
        //    }
        //    case Core::FE::CellType::tet10:
        //    {
        //      typedef NitscheCoupling<distype,Core::FE::CellType::tet10,4>
        //      NitscheCouplType; nit = new
        //      NitscheCouplType(vele_xyz,C_umum,C_usum,C_umus,C_usus,rhC_um,rhC_us,is_viscAdjointSymmetric);
        //      break;
        //    }
      case Core::FE::CellType::hex8:
      {
        typedef NitscheCoupling<distype, Core::FE::CellType::hex8, 4> NitscheCouplType;
        nit = new NitscheCouplType(
            vele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case Core::FE::CellType::hex20:
      {
        typedef NitscheCoupling<distype, Core::FE::CellType::hex20, 4> NitscheCouplType;
        nit = new NitscheCouplType(
            vele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case Core::FE::CellType::hex27:
      {
        typedef NitscheCoupling<distype, Core::FE::CellType::hex27, 4> NitscheCouplType;
        nit = new NitscheCouplType(
            vele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported volume element shape %d.", vele->Shape());
        break;
    }
  }
  else if (numdofpernode == 3)
  {
    switch (vele->Shape())
    {
      case Core::FE::CellType::hex8:
      {
        typedef NitscheCoupling<distype, Core::FE::CellType::hex8, 3> NitscheCouplType;
        nit = new NitscheCouplType(
            vele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      default:
        // expecting 3 dofs per slave element node as this is fluid-solid coupling
        FOUR_C_THROW("Unsupported volume element shape %d.", vele->Shape());
        break;
    }
  }

  return Teuchos::rcp(nit);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Teuchos::RCP<HybridLMInterface<distype>>
HybridLMInterface<distype>::create_hybrid_lm_coupling_x_fluid_wdbc(
    bool is_viscAdjointSymmetric  ///< flag that indicates equal signs of Nitsche's standard &
                                  ///< adjoint viscous term
)
{
  HybridLMInterface* hybridlm = nullptr;
  typedef HybridLMCoupling<distype, Core::FE::CellType::dis_none, 3> HybridLMCouplType;
  hybridlm = new HybridLMCouplType(is_viscAdjointSymmetric);

  return Teuchos::rcp(hybridlm);
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Teuchos::RCP<HybridLMInterface<distype>>
HybridLMInterface<distype>::create_hybrid_lm_coupling_x_fluid_wdbc(
    Core::Elements::Element* bele,              ///< boundary element
    Core::LinAlg::SerialDenseMatrix& bele_xyz,  ///< global node coordinates of boundary element
    bool is_viscAdjointSymmetric  ///< flag that indicates equal signs of Nitsche's standard &
                                  ///< adjoint viscous term
)
{
  // get number of dofs for this boundary element
  const unsigned numdofpernode = bele->NumDofPerNode(*bele->Nodes()[0]);

  // three dofs per node, for standard Dirichlet coupling, four dofs per node for fluid-fluid
  // coupling
  if (numdofpernode != 3)
  {
    FOUR_C_THROW("Unsupported number of %d nodes for standard Dirichlet coupling.", numdofpernode);
  }

  switch (bele->Shape())
  {
      //    case Core::FE::CellType::tri3:
      //    {
      //      typedef HybridLMCoupling<distype,Core::FE::CellType::tri3,3>
      //      HybridLMCouplType; return Teuchos::rcp(new HybridLMCouplType(bele_xyz)); break;
      //    }
      //    case Core::FE::CellType::tri6:
      //    {
      //      typedef HybridLMCoupling<distype,Core::FE::CellType::tri6,3>
      //      HybridLMCouplType; return Teuchos::rcp(new HybridLMCouplType(bele_xyz)); break;
      //    }
    case Core::FE::CellType::quad4:
    {
      typedef HybridLMCoupling<distype, Core::FE::CellType::quad4, 3> HybridLMCouplType;
      return Teuchos::rcp(new HybridLMCouplType(bele_xyz, is_viscAdjointSymmetric));
      break;
    }
    case Core::FE::CellType::quad8:
    {
      typedef HybridLMCoupling<distype, Core::FE::CellType::quad8, 3> HybridLMCouplType;
      return Teuchos::rcp(new HybridLMCouplType(bele_xyz, is_viscAdjointSymmetric));
      break;
    }
    case Core::FE::CellType::quad9:
    {
      typedef HybridLMCoupling<distype, Core::FE::CellType::quad9, 3> HybridLMCouplType;
      return Teuchos::rcp(new HybridLMCouplType(bele_xyz, is_viscAdjointSymmetric));
      break;
    }
    default:
      FOUR_C_THROW("Unsupported boundary element shape %d", bele->Shape());
      break;
  }

  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Teuchos::RCP<HybridLMInterface<distype>>
HybridLMInterface<distype>::create_hybrid_lm_coupling_x_fluid_sided(
    Core::Elements::Element* bele,              ///< boundary element
    Core::LinAlg::SerialDenseMatrix& bele_xyz,  ///< global node coordinates of boundary element
    Core::LinAlg::SerialDenseMatrix& C_usum,    ///< C_usum coupling matrix
    Core::LinAlg::SerialDenseMatrix& C_umus,    ///< C_umus coupling matrix
    Core::LinAlg::SerialDenseMatrix& rhC_us,    ///< C_us coupling rhs
    Core::LinAlg::SerialDenseMatrix& G_s_us,    ///< \f$G_{u^s \sigma}\f$ coupling matrix
    Core::LinAlg::SerialDenseMatrix& G_us_s,    ///< \f$G_{\sigma u^s}\f$ coupling matrix
    bool is_viscAdjointSymmetric  ///< flag that indicates equal signs of Nitsche's standard &
                                  ///< adjoint viscous term
)
{
  HybridLMInterface* hlm = nullptr;

  // get number of dofs for this boundary element
  const unsigned numdofpernode = bele->NumDofPerNode(*bele->Nodes()[0]);

  // three dofs per node, for monolithic XFSI, four dofs per node for backgrounds-sided fluid-fluid
  // coupling
  if (numdofpernode == 3)
  {
    switch (bele->Shape())
    {
        //    case Core::FE::CellType::tri3:
        //    {
        //      typedef HybridLMCoupling<distype,Core::FE::CellType::tri3,3>
        //      HybridLMCouplType; hlm = new
        //      HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,is_viscAdjointSymmetric);
        //      break;
        //    }
        //    case Core::FE::CellType::tri6:
        //    {
        //      typedef HybridLMCoupling<distype,Core::FE::CellType::tri6,3>
        //      HybridLMCouplType; hlm = new
        //      HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,is_viscAdjointSymmetric);
        //      break;
        //    }
      case Core::FE::CellType::quad4:
      {
        typedef HybridLMCoupling<distype, Core::FE::CellType::quad4, 3> HybridLMCouplType;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        typedef HybridLMCoupling<distype, Core::FE::CellType::quad8, 3> HybridLMCouplType;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        typedef HybridLMCoupling<distype, Core::FE::CellType::quad9, 3> HybridLMCouplType;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported boundary element shape %d", bele->Shape());
        break;
    }
  }
  else if (numdofpernode == 4)
  {
    switch (bele->Shape())
    {
        //      case Core::FE::CellType::tri3:
        //      {
        //        typedef HybridLMCoupling<distype,Core::FE::CellType::tri3,4>
        //        HybridLMCouplType; hlm = new
        //        HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,is_viscAdjointSymmetric);
        //        break;
        //      }
        //      case Core::FE::CellType::tri6:
        //      {
        //        typedef HybridLMCoupling<distype,Core::FE::CellType::tri6,4>
        //        HybridLMCouplType; hlm = new
        //        HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,is_viscAdjointSymmetric);
        //        break;
        //      }
      case Core::FE::CellType::quad4:
      {
        typedef HybridLMCoupling<distype, Core::FE::CellType::quad4, 4> HybridLMCouplType;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      case Core::FE::CellType::quad8:
      {
        typedef HybridLMCoupling<distype, Core::FE::CellType::quad8, 4> HybridLMCouplType;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      case Core::FE::CellType::quad9:
      {
        typedef HybridLMCoupling<distype, Core::FE::CellType::quad9, 4> HybridLMCouplType;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      default:
        FOUR_C_THROW("Unsupported boundary element shape %d", bele->Shape());
        break;
    }
  }
  else
    FOUR_C_THROW("Unsupported number of %d nodes for coupling slave element.", numdofpernode);

  return Teuchos::rcp(hlm);
}

template class Discret::ELEMENTS::XFLUID::SlaveElementInterface<Core::FE::CellType::hex8>;
template class Discret::ELEMENTS::XFLUID::SlaveElementInterface<Core::FE::CellType::hex20>;
template class Discret::ELEMENTS::XFLUID::SlaveElementInterface<Core::FE::CellType::hex27>;
template class Discret::ELEMENTS::XFLUID::SlaveElementInterface<Core::FE::CellType::tet4>;
template class Discret::ELEMENTS::XFLUID::SlaveElementInterface<Core::FE::CellType::tet10>;
template class Discret::ELEMENTS::XFLUID::SlaveElementInterface<Core::FE::CellType::wedge6>;
template class Discret::ELEMENTS::XFLUID::SlaveElementInterface<Core::FE::CellType::wedge15>;
// template class
// Discret::ELEMENTS::XFLUID::SlaveElementInterface<Core::FE::CellType::pyramid5>;

template class Discret::ELEMENTS::XFLUID::NitscheInterface<Core::FE::CellType::hex8>;
template class Discret::ELEMENTS::XFLUID::NitscheInterface<Core::FE::CellType::hex20>;
template class Discret::ELEMENTS::XFLUID::NitscheInterface<Core::FE::CellType::hex27>;
template class Discret::ELEMENTS::XFLUID::NitscheInterface<Core::FE::CellType::tet4>;
template class Discret::ELEMENTS::XFLUID::NitscheInterface<Core::FE::CellType::tet10>;
template class Discret::ELEMENTS::XFLUID::NitscheInterface<Core::FE::CellType::wedge6>;
template class Discret::ELEMENTS::XFLUID::NitscheInterface<Core::FE::CellType::wedge15>;
// template class
// Discret::ELEMENTS::XFLUID::NitscheInterface<Core::FE::CellType::pyramid5>;

template class Discret::ELEMENTS::XFLUID::HybridLMInterface<Core::FE::CellType::hex8>;
template class Discret::ELEMENTS::XFLUID::HybridLMInterface<Core::FE::CellType::hex20>;
template class Discret::ELEMENTS::XFLUID::HybridLMInterface<Core::FE::CellType::hex27>;
template class Discret::ELEMENTS::XFLUID::HybridLMInterface<Core::FE::CellType::tet4>;
template class Discret::ELEMENTS::XFLUID::HybridLMInterface<Core::FE::CellType::tet10>;
template class Discret::ELEMENTS::XFLUID::HybridLMInterface<Core::FE::CellType::wedge6>;
template class Discret::ELEMENTS::XFLUID::HybridLMInterface<Core::FE::CellType::wedge15>;
// template class
// Discret::ELEMENTS::XFLUID::HybridLMInterface<Core::FE::CellType::pyramid5>;

FOUR_C_NAMESPACE_CLOSE

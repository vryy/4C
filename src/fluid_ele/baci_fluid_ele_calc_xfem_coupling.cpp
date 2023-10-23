/*----------------------------------------------------------------------*/
/*! \file

\brief Factory class for providing an implementation for coupling with
       Mixed/Stress/Hybrid methods or Nitsche's method to enforce
       interface conditions in the XFEM weakly

\level 2


*/
/*----------------------------------------------------------------------*/


#include "baci_fluid_ele_calc_xfem_coupling.H"

#include "baci_bele_bele3.H"
#include "baci_fluid_ele_calc_xfem_coupling_impl.H"

using namespace DRT::ELEMENTS::XFLUID;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void SlaveElementInterface<distype>::DefineStateNames(
    DRT::Element::DiscretizationType slave_distype,  ///< coupling slave discretization type
    std::string& disp_statename,                     ///< name of displacement state at current step
    std::string& vel_statename,                      ///< name of velocity state at current step
    std::string& veln_statename                      ///< name of velocity state at previous step
)
{
  switch (slave_distype)
  {
      //  case DRT::Element::DiscretizationType::tri3:
      //  case DRT::Element::DiscretizationType::tri6:
    case DRT::Element::DiscretizationType::quad4:
    case DRT::Element::DiscretizationType::quad8:
    case DRT::Element::DiscretizationType::quad9:
    {
      disp_statename = std::string("idispnp");
      vel_statename = std::string("ivelnp");
      veln_statename = std::string("iveln");
      break;
    }
    case DRT::Element::DiscretizationType::hex8:
    case DRT::Element::DiscretizationType::hex20:
    case DRT::Element::DiscretizationType::hex27:
    {
      disp_statename = std::string("dispnp");
      vel_statename = std::string("velaf");
      veln_statename = std::string("veln");
      break;
    }
    default:
      dserror("Unsupported element shape %d", slave_distype);
      break;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
Teuchos::RCP<SlaveElementInterface<distype>>
SlaveElementInterface<distype>::CreateSlaveElementRepresentation(
    DRT::Element* slave_ele,  ///< coupling slave element
    CORE::LINALG::SerialDenseMatrix&
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
        //      case DRT::Element::DiscretizationType::tri3:
        //      {
        //        typedef
        //        SlaveElementRepresentation<distype,DRT::Element::DiscretizationType::tri3,3>
        //        SlaveEleType; sla = new SlaveEleType(slave_xyz); break;
        //      }
        //      case DRT::Element::DiscretizationType::tri6:
        //      {
        //        typedef
        //        SlaveElementRepresentation<distype,DRT::Element::DiscretizationType::tri6,3>
        //        SlaveEleType; sla = new SlaveEleType(slave_xyz); break;
        //      }
      case DRT::Element::DiscretizationType::quad4:
      {
        typedef SlaveElementRepresentation<distype, DRT::Element::DiscretizationType::quad4, 3>
            SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::DiscretizationType::quad8:
      {
        typedef SlaveElementRepresentation<distype, DRT::Element::DiscretizationType::quad8, 3>
            SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::DiscretizationType::quad9:
      {
        typedef SlaveElementRepresentation<distype, DRT::Element::DiscretizationType::quad9, 3>
            SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      default:
        dserror("Unsupported boundary element shape %d", slave_ele->Shape());
        break;
    }
  }
  else if (numdofpernode ==
           4)  // volumetric coupling partners only required for fluid-fluid coupling
  {
    switch (slave_ele->Shape())
    {
        //      case DRT::Element::DiscretizationType::tri3:
        //      {
        //        typedef
        //        SlaveElementRepresentation<distype,DRT::Element::DiscretizationType::tri3,4>
        //        SlaveEleType; sla = new SlaveEleType(slave_xyz); break;
        //      }
        //      case DRT::Element::DiscretizationType::tri6:
        //      {
        //        typedef
        //        SlaveElementRepresentation<distype,DRT::Element::DiscretizationType::tri6,4>
        //        SlaveEleType; sla = new SlaveEleType(slave_xyz); break;
        //      }
      case DRT::Element::DiscretizationType::quad4:
      {
        typedef SlaveElementRepresentation<distype, DRT::Element::DiscretizationType::quad4, 4>
            SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::DiscretizationType::quad8:
      {
        typedef SlaveElementRepresentation<distype, DRT::Element::DiscretizationType::quad8, 4>
            SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::DiscretizationType::quad9:
      {
        typedef SlaveElementRepresentation<distype, DRT::Element::DiscretizationType::quad9, 4>
            SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::DiscretizationType::hex8:
      {
        typedef SlaveElementRepresentation<distype, DRT::Element::DiscretizationType::hex8, 4>
            SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::DiscretizationType::hex20:
      {
        typedef SlaveElementRepresentation<distype, DRT::Element::DiscretizationType::hex20, 4>
            SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::DiscretizationType::hex27:
      {
        typedef SlaveElementRepresentation<distype, DRT::Element::DiscretizationType::hex27, 4>
            SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      default:
        dserror("unsupported boundary element shape %d", slave_ele->Shape());
        break;
    }
  }

  return Teuchos::rcp(sla);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
Teuchos::RCP<NitscheInterface<distype>> NitscheInterface<distype>::CreateNitscheCoupling_XFluidWDBC(
    CORE::LINALG::SerialDenseMatrix::Base& C_umum, CORE::LINALG::SerialDenseMatrix::Base& rhC_um,
    const DRT::ELEMENTS::FluidEleParameterXFEM& fldparaxfem)
{
  NitscheInterface* nit = nullptr;
  typedef NitscheCoupling<distype, DRT::Element::DiscretizationType::dis_none, 3> NitscheCouplType;
  nit = new NitscheCouplType(C_umum, rhC_um, fldparaxfem);

  return Teuchos::rcp(nit);
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
Teuchos::RCP<NitscheInterface<distype>> NitscheInterface<distype>::CreateNitscheCoupling_XFluidWDBC(
    DRT::Element* bele, CORE::LINALG::SerialDenseMatrix::Base& bele_xyz,
    CORE::LINALG::SerialDenseMatrix::Base& C_umum, CORE::LINALG::SerialDenseMatrix::Base& rhC_um,
    const DRT::ELEMENTS::FluidEleParameterXFEM& fldparaxfem)
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
      //      case DRT::Element::DiscretizationType::tri3:
      //      {
      //        typedef NitscheCoupling<distype,DRT::Element::DiscretizationType::tri3,3>
      //        NitscheCouplType; nit = new NitscheCouplType(bele_xyz,C_umum,rhC_um,); break;
      //      }
      //      case DRT::Element::DiscretizationType::tri6:
      //      {
      //        typedef NitscheCoupling<distype,DRT::Element::DiscretizationType::tri6,3>
      //        NitscheCouplType; nit = new NitscheCouplType(bele_xyz,C_umum,rhC_um,); break;
      //      }
      case DRT::Element::DiscretizationType::quad4:
      {
        typedef NitscheCoupling<distype, DRT::Element::DiscretizationType::quad4, 3>
            NitscheCouplType;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      case DRT::Element::DiscretizationType::quad8:
      {
        typedef NitscheCoupling<distype, DRT::Element::DiscretizationType::quad8, 3>
            NitscheCouplType;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      case DRT::Element::DiscretizationType::quad9:
      {
        typedef NitscheCoupling<distype, DRT::Element::DiscretizationType::quad9, 3>
            NitscheCouplType;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      default:
        dserror("Unsupported boundary element shape %d", bele->Shape());
        break;
    }
  }
  else if (numdofpernode == 4)
  {
    switch (bele->Shape())
    {
      //      case DRT::Element::DiscretizationType::tri3:
      //      {
      //        typedef NitscheCoupling<distype,DRT::Element::DiscretizationType::tri3,3>
      //        NitscheCouplType; nit = new NitscheCouplType(bele_xyz,C_umum,rhC_um,); break;
      //      }
      //      case DRT::Element::DiscretizationType::tri6:
      //      {
      //        typedef NitscheCoupling<distype,DRT::Element::DiscretizationType::tri6,3>
      //        NitscheCouplType; nit = new NitscheCouplType(bele_xyz,C_umum,rhC_um,); break;
      //      }
      case DRT::Element::DiscretizationType::quad4:
      {
        typedef NitscheCoupling<distype, DRT::Element::DiscretizationType::quad4, 4>
            NitscheCouplType;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      case DRT::Element::DiscretizationType::quad8:
      {
        typedef NitscheCoupling<distype, DRT::Element::DiscretizationType::quad8, 4>
            NitscheCouplType;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      case DRT::Element::DiscretizationType::quad9:
      {
        typedef NitscheCoupling<distype, DRT::Element::DiscretizationType::quad9, 4>
            NitscheCouplType;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      default:
        dserror("Unsupported boundary element shape %d", bele->Shape());
        break;
    }
  }
  else
    dserror("Unsupported number of %d nodes for coupling slave element.", numdofpernode);

  return Teuchos::rcp(nit);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
Teuchos::RCP<NitscheInterface<distype>>
NitscheInterface<distype>::CreateNitscheCoupling_XFluidSided(DRT::Element* bele,
    CORE::LINALG::SerialDenseMatrix::Base& bele_xyz, CORE::LINALG::SerialDenseMatrix::Base& C_umum,
    CORE::LINALG::SerialDenseMatrix::Base& C_usum, CORE::LINALG::SerialDenseMatrix::Base& C_umus,
    CORE::LINALG::SerialDenseMatrix::Base& C_usus, CORE::LINALG::SerialDenseMatrix::Base& rhC_um,
    CORE::LINALG::SerialDenseMatrix::Base& rhC_us,
    const DRT::ELEMENTS::FluidEleParameterXFEM& fldparaxfem)
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
        //    case DRT::Element::DiscretizationType::tri3:
        //    {
        //      typedef NitscheCoupling<distype,DRT::Element::DiscretizationType::tri3,3>
        //      NitscheCouplType; nit = new
        //      NitscheCouplType(bele_xyz,C_umum,C_usum,C_umus,C_usus,rhC_um,rhC_us,is_viscAdjointSymmetric);
        //      break;
        //    }
        //    case DRT::Element::DiscretizationType::tri6:
        //    {
        //      typedef NitscheCoupling<distype,DRT::Element::DiscretizationType::tri6,3>
        //      NitscheCouplType; nit = new
        //      NitscheCouplType(bele_xyz,C_umum,C_usum,C_umus,C_usus,rhC_um,rhC_us,is_viscAdjointSymmetric);
        //      break;
        //    }
      case DRT::Element::DiscretizationType::quad4:
      {
        typedef NitscheCoupling<distype, DRT::Element::DiscretizationType::quad4, 3>
            NitscheCouplType;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case DRT::Element::DiscretizationType::quad8:
      {
        typedef NitscheCoupling<distype, DRT::Element::DiscretizationType::quad8, 3>
            NitscheCouplType;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case DRT::Element::DiscretizationType::quad9:
      {
        typedef NitscheCoupling<distype, DRT::Element::DiscretizationType::quad9, 3>
            NitscheCouplType;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      default:
        dserror("Unsupported boundary element shape %d", bele->Shape());
        break;
    }
  }
  else if (numdofpernode == 4)
  {
    switch (bele->Shape())
    {
      //    case DRT::Element::DiscretizationType::tri3:
      //    {
      //      typedef NitscheCoupling<distype,DRT::Element::DiscretizationType::tri3,4>
      //      NitscheCouplType; nit = new
      //      NitscheCouplType(bele_xyz,C_umum,C_usum,C_umus,C_usus,rhC_um,rhC_us,is_viscAdjointSymmetric);
      //      break;
      //    }
      //    case DRT::Element::DiscretizationType::tri6:
      //    {
      //      typedef NitscheCoupling<distype,DRT::Element::DiscretizationType::tri6,4>
      //      NitscheCouplType; nit = new
      //      NitscheCouplType(bele_xyz,C_umum,C_usum,C_umus,C_usus,rhC_um,rhC_us,is_viscAdjointSymmetric);
      //      break;
      //    }
      case DRT::Element::DiscretizationType::quad4:
      {
        typedef NitscheCoupling<distype, DRT::Element::DiscretizationType::quad4, 4>
            NitscheCouplType;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case DRT::Element::DiscretizationType::quad8:
      {
        typedef NitscheCoupling<distype, DRT::Element::DiscretizationType::quad8, 4>
            NitscheCouplType;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case DRT::Element::DiscretizationType::quad9:
      {
        typedef NitscheCoupling<distype, DRT::Element::DiscretizationType::quad9, 4>
            NitscheCouplType;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      default:
        dserror("Unsupported boundary element shape %d", bele->Shape());
        break;
    }
  }
  else
    dserror("Unsupported number of %d nodes for coupling slave element.", numdofpernode);

  return Teuchos::rcp(nit);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
Teuchos::RCP<NitscheInterface<distype>> NitscheInterface<distype>::CreateNitscheCoupling_TwoSided(
    DRT::Element* vele, CORE::LINALG::SerialDenseMatrix::Base& vele_xyz,
    CORE::LINALG::SerialDenseMatrix::Base& C_umum, CORE::LINALG::SerialDenseMatrix::Base& C_usum,
    CORE::LINALG::SerialDenseMatrix::Base& C_umus, CORE::LINALG::SerialDenseMatrix::Base& C_usus,
    CORE::LINALG::SerialDenseMatrix::Base& rhC_um, CORE::LINALG::SerialDenseMatrix::Base& rhC_us,
    const DRT::ELEMENTS::FluidEleParameterXFEM& fldparaxfem)
{
  NitscheInterface* nit = nullptr;

  // get number of dofs for the embedded element
  const unsigned numdofpernode = vele->NumDofPerNode(*vele->Nodes()[0]);

  if (numdofpernode == 4)
  {
    switch (vele->Shape())
    {
        //    case DRT::Element::DiscretizationType::tet4:
        //    {
        //      typedef NitscheCoupling<distype,DRT::Element::DiscretizationType::tet4,4>
        //      NitscheCouplType; nit = new
        //      NitscheCouplType(vele_xyz,C_umum,C_usum,C_umus,C_usus,rhC_um,rhC_us,is_viscAdjointSymmetric);
        //      break;
        //    }
        //    case DRT::Element::DiscretizationType::tet10:
        //    {
        //      typedef NitscheCoupling<distype,DRT::Element::DiscretizationType::tet10,4>
        //      NitscheCouplType; nit = new
        //      NitscheCouplType(vele_xyz,C_umum,C_usum,C_umus,C_usus,rhC_um,rhC_us,is_viscAdjointSymmetric);
        //      break;
        //    }
      case DRT::Element::DiscretizationType::hex8:
      {
        typedef NitscheCoupling<distype, DRT::Element::DiscretizationType::hex8, 4>
            NitscheCouplType;
        nit = new NitscheCouplType(
            vele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case DRT::Element::DiscretizationType::hex20:
      {
        typedef NitscheCoupling<distype, DRT::Element::DiscretizationType::hex20, 4>
            NitscheCouplType;
        nit = new NitscheCouplType(
            vele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case DRT::Element::DiscretizationType::hex27:
      {
        typedef NitscheCoupling<distype, DRT::Element::DiscretizationType::hex27, 4>
            NitscheCouplType;
        nit = new NitscheCouplType(
            vele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      default:
        dserror("Unsupported volume element shape %d.", vele->Shape());
        break;
    }
  }
  else if (numdofpernode == 3)
  {
    switch (vele->Shape())
    {
      case DRT::Element::DiscretizationType::hex8:
      {
        typedef NitscheCoupling<distype, DRT::Element::DiscretizationType::hex8, 3>
            NitscheCouplType;
        nit = new NitscheCouplType(
            vele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      default:
        // expecting 3 dofs per slave element node as this is fluid-solid coupling
        dserror("Unsupported volume element shape %d.", vele->Shape());
        break;
    }
  }

  return Teuchos::rcp(nit);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
Teuchos::RCP<HybridLMInterface<distype>>
HybridLMInterface<distype>::CreateHybridLMCoupling_XFluidWDBC(
    bool is_viscAdjointSymmetric  ///< flag that indicates equal signs of Nitsche's standard &
                                  ///< adjoint viscous term
)
{
  HybridLMInterface* hybridlm = nullptr;
  typedef HybridLMCoupling<distype, DRT::Element::DiscretizationType::dis_none, 3>
      HybridLMCouplType;
  hybridlm = new HybridLMCouplType(is_viscAdjointSymmetric);

  return Teuchos::rcp(hybridlm);
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
Teuchos::RCP<HybridLMInterface<distype>>
HybridLMInterface<distype>::CreateHybridLMCoupling_XFluidWDBC(
    DRT::Element* bele,                         ///< boundary element
    CORE::LINALG::SerialDenseMatrix& bele_xyz,  ///< global node coordinates of boundary element
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
    dserror("Unsupported number of %d nodes for standard Dirichlet coupling.", numdofpernode);
  }

  switch (bele->Shape())
  {
      //    case DRT::Element::DiscretizationType::tri3:
      //    {
      //      typedef HybridLMCoupling<distype,DRT::Element::DiscretizationType::tri3,3>
      //      HybridLMCouplType; return Teuchos::rcp(new HybridLMCouplType(bele_xyz)); break;
      //    }
      //    case DRT::Element::DiscretizationType::tri6:
      //    {
      //      typedef HybridLMCoupling<distype,DRT::Element::DiscretizationType::tri6,3>
      //      HybridLMCouplType; return Teuchos::rcp(new HybridLMCouplType(bele_xyz)); break;
      //    }
    case DRT::Element::DiscretizationType::quad4:
    {
      typedef HybridLMCoupling<distype, DRT::Element::DiscretizationType::quad4, 3>
          HybridLMCouplType;
      return Teuchos::rcp(new HybridLMCouplType(bele_xyz, is_viscAdjointSymmetric));
      break;
    }
    case DRT::Element::DiscretizationType::quad8:
    {
      typedef HybridLMCoupling<distype, DRT::Element::DiscretizationType::quad8, 3>
          HybridLMCouplType;
      return Teuchos::rcp(new HybridLMCouplType(bele_xyz, is_viscAdjointSymmetric));
      break;
    }
    case DRT::Element::DiscretizationType::quad9:
    {
      typedef HybridLMCoupling<distype, DRT::Element::DiscretizationType::quad9, 3>
          HybridLMCouplType;
      return Teuchos::rcp(new HybridLMCouplType(bele_xyz, is_viscAdjointSymmetric));
      break;
    }
    default:
      dserror("Unsupported boundary element shape %d", bele->Shape());
      break;
  }

  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
Teuchos::RCP<HybridLMInterface<distype>>
HybridLMInterface<distype>::CreateHybridLMCoupling_XFluidSided(
    DRT::Element* bele,                         ///< boundary element
    CORE::LINALG::SerialDenseMatrix& bele_xyz,  ///< global node coordinates of boundary element
    CORE::LINALG::SerialDenseMatrix& C_usum,    ///< C_usum coupling matrix
    CORE::LINALG::SerialDenseMatrix& C_umus,    ///< C_umus coupling matrix
    CORE::LINALG::SerialDenseMatrix& rhC_us,    ///< C_us coupling rhs
    CORE::LINALG::SerialDenseMatrix& G_s_us,    ///< \f$G_{u^s \sigma}\f$ coupling matrix
    CORE::LINALG::SerialDenseMatrix& G_us_s,    ///< \f$G_{\sigma u^s}\f$ coupling matrix
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
        //    case DRT::Element::DiscretizationType::tri3:
        //    {
        //      typedef HybridLMCoupling<distype,DRT::Element::DiscretizationType::tri3,3>
        //      HybridLMCouplType; hlm = new
        //      HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,is_viscAdjointSymmetric);
        //      break;
        //    }
        //    case DRT::Element::DiscretizationType::tri6:
        //    {
        //      typedef HybridLMCoupling<distype,DRT::Element::DiscretizationType::tri6,3>
        //      HybridLMCouplType; hlm = new
        //      HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,is_viscAdjointSymmetric);
        //      break;
        //    }
      case DRT::Element::DiscretizationType::quad4:
      {
        typedef HybridLMCoupling<distype, DRT::Element::DiscretizationType::quad4, 3>
            HybridLMCouplType;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      case DRT::Element::DiscretizationType::quad8:
      {
        typedef HybridLMCoupling<distype, DRT::Element::DiscretizationType::quad8, 3>
            HybridLMCouplType;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      case DRT::Element::DiscretizationType::quad9:
      {
        typedef HybridLMCoupling<distype, DRT::Element::DiscretizationType::quad9, 3>
            HybridLMCouplType;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      default:
        dserror("Unsupported boundary element shape %d", bele->Shape());
        break;
    }
  }
  else if (numdofpernode == 4)
  {
    switch (bele->Shape())
    {
        //      case DRT::Element::DiscretizationType::tri3:
        //      {
        //        typedef HybridLMCoupling<distype,DRT::Element::DiscretizationType::tri3,4>
        //        HybridLMCouplType; hlm = new
        //        HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,is_viscAdjointSymmetric);
        //        break;
        //      }
        //      case DRT::Element::DiscretizationType::tri6:
        //      {
        //        typedef HybridLMCoupling<distype,DRT::Element::DiscretizationType::tri6,4>
        //        HybridLMCouplType; hlm = new
        //        HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,is_viscAdjointSymmetric);
        //        break;
        //      }
      case DRT::Element::DiscretizationType::quad4:
      {
        typedef HybridLMCoupling<distype, DRT::Element::DiscretizationType::quad4, 4>
            HybridLMCouplType;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      case DRT::Element::DiscretizationType::quad8:
      {
        typedef HybridLMCoupling<distype, DRT::Element::DiscretizationType::quad8, 4>
            HybridLMCouplType;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      case DRT::Element::DiscretizationType::quad9:
      {
        typedef HybridLMCoupling<distype, DRT::Element::DiscretizationType::quad9, 4>
            HybridLMCouplType;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      default:
        dserror("Unsupported boundary element shape %d", bele->Shape());
        break;
    }
  }
  else
    dserror("Unsupported number of %d nodes for coupling slave element.", numdofpernode);

  return Teuchos::rcp(hlm);
}

template class DRT::ELEMENTS::XFLUID::SlaveElementInterface<DRT::Element::DiscretizationType::hex8>;
template class DRT::ELEMENTS::XFLUID::SlaveElementInterface<
    DRT::Element::DiscretizationType::hex20>;
template class DRT::ELEMENTS::XFLUID::SlaveElementInterface<
    DRT::Element::DiscretizationType::hex27>;
template class DRT::ELEMENTS::XFLUID::SlaveElementInterface<DRT::Element::DiscretizationType::tet4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementInterface<
    DRT::Element::DiscretizationType::tet10>;
template class DRT::ELEMENTS::XFLUID::SlaveElementInterface<
    DRT::Element::DiscretizationType::wedge6>;
template class DRT::ELEMENTS::XFLUID::SlaveElementInterface<
    DRT::Element::DiscretizationType::wedge15>;
// template class
// DRT::ELEMENTS::XFLUID::SlaveElementInterface<DRT::Element::DiscretizationType::pyramid5>;

template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::DiscretizationType::hex8>;
template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::DiscretizationType::hex20>;
template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::DiscretizationType::hex27>;
template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::DiscretizationType::tet4>;
template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::DiscretizationType::tet10>;
template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::DiscretizationType::wedge6>;
template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::DiscretizationType::wedge15>;
// template class
// DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::DiscretizationType::pyramid5>;

template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::DiscretizationType::hex8>;
template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::DiscretizationType::hex20>;
template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::DiscretizationType::hex27>;
template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::DiscretizationType::tet4>;
template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::DiscretizationType::tet10>;
template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::DiscretizationType::wedge6>;
template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::DiscretizationType::wedge15>;
// template class
// DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::DiscretizationType::pyramid5>;

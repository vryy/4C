/*----------------------------------------------------------------------*/
/*!

\brief Factory class for providing an implementation for coupling with
       Mixed/Stress/Hybrid methods or Nitsche's method to enforce
       interface conditions in the XFEM weakly

\level 2

\maintainer  Christoph Ager

*/
/*----------------------------------------------------------------------*/


#include "fluid_ele_calc_xfem_coupling.H"
#include "fluid_ele_calc_xfem_coupling_impl.H"

#include "../drt_bele3/bele3.H"

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
      //  case DRT::Element::tri3:
      //  case DRT::Element::tri6:
    case DRT::Element::quad4:
    case DRT::Element::quad8:
    case DRT::Element::quad9:
    {
      disp_statename = std::string("idispnp");
      vel_statename = std::string("ivelnp");
      veln_statename = std::string("iveln");
      break;
    }
    case DRT::Element::hex8:
    case DRT::Element::hex20:
    case DRT::Element::hex27:
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
    DRT::Element* slave_ele,             ///< coupling slave element
    Epetra_SerialDenseMatrix& slave_xyz  ///< global node coordinates of coupling slave element
)
{
  SlaveElementInterface* sla = NULL;

  // get number of dofs for this slave element
  const unsigned numdofpernode = slave_ele->NumDofPerNode(*slave_ele->Nodes()[0]);

  if (numdofpernode == 3)
  {
    switch (slave_ele->Shape())
    {
        //      case DRT::Element::tri3:
        //      {
        //        typedef SlaveElementRepresentation<distype,DRT::Element::tri3,3> SlaveEleType;
        //        sla = new SlaveEleType(slave_xyz);
        //        break;
        //      }
        //      case DRT::Element::tri6:
        //      {
        //        typedef SlaveElementRepresentation<distype,DRT::Element::tri6,3> SlaveEleType;
        //        sla = new SlaveEleType(slave_xyz);
        //        break;
        //      }
      case DRT::Element::quad4:
      {
        typedef SlaveElementRepresentation<distype, DRT::Element::quad4, 3> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::quad8:
      {
        typedef SlaveElementRepresentation<distype, DRT::Element::quad8, 3> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::quad9:
      {
        typedef SlaveElementRepresentation<distype, DRT::Element::quad9, 3> SlaveEleType;
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
        //      case DRT::Element::tri3:
        //      {
        //        typedef SlaveElementRepresentation<distype,DRT::Element::tri3,4> SlaveEleType;
        //        sla = new SlaveEleType(slave_xyz);
        //        break;
        //      }
        //      case DRT::Element::tri6:
        //      {
        //        typedef SlaveElementRepresentation<distype,DRT::Element::tri6,4> SlaveEleType;
        //        sla = new SlaveEleType(slave_xyz);
        //        break;
        //      }
      case DRT::Element::quad4:
      {
        typedef SlaveElementRepresentation<distype, DRT::Element::quad4, 4> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::quad8:
      {
        typedef SlaveElementRepresentation<distype, DRT::Element::quad8, 4> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::quad9:
      {
        typedef SlaveElementRepresentation<distype, DRT::Element::quad9, 4> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::hex8:
      {
        typedef SlaveElementRepresentation<distype, DRT::Element::hex8, 4> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::hex20:
      {
        typedef SlaveElementRepresentation<distype, DRT::Element::hex20, 4> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::hex27:
      {
        typedef SlaveElementRepresentation<distype, DRT::Element::hex27, 4> SlaveEleType;
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
    Epetra_SerialDenseMatrix& C_umum, Epetra_SerialDenseMatrix& rhC_um,
    const DRT::ELEMENTS::FluidEleParameterXFEM& fldparaxfem)
{
  NitscheInterface* nit = NULL;
  typedef NitscheCoupling<distype, DRT::Element::dis_none, 3> NitscheCouplType;
  nit = new NitscheCouplType(C_umum, rhC_um, fldparaxfem);

  return Teuchos::rcp(nit);
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
Teuchos::RCP<NitscheInterface<distype>> NitscheInterface<distype>::CreateNitscheCoupling_XFluidWDBC(
    DRT::Element* bele, Epetra_SerialDenseMatrix& bele_xyz, Epetra_SerialDenseMatrix& C_umum,
    Epetra_SerialDenseMatrix& rhC_um, const DRT::ELEMENTS::FluidEleParameterXFEM& fldparaxfem)
{
  NitscheInterface* nit = NULL;

  // get number of dofs for this boundary element
  const unsigned numdofpernode = bele->NumDofPerNode(*bele->Nodes()[0]);

  // three dofs per node, for standard Dirichlet coupling, four dofs per node for background
  // geometry coupling
  if (numdofpernode == 3)
  {
    switch (bele->Shape())
    {
      //      case DRT::Element::tri3:
      //      {
      //        typedef NitscheCoupling<distype,DRT::Element::tri3,3> NitscheCouplType;
      //        nit = new NitscheCouplType(bele_xyz,C_umum,rhC_um,);
      //        break;
      //      }
      //      case DRT::Element::tri6:
      //      {
      //        typedef NitscheCoupling<distype,DRT::Element::tri6,3> NitscheCouplType;
      //        nit = new NitscheCouplType(bele_xyz,C_umum,rhC_um,);
      //        break;
      //      }
      case DRT::Element::quad4:
      {
        typedef NitscheCoupling<distype, DRT::Element::quad4, 3> NitscheCouplType;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      case DRT::Element::quad8:
      {
        typedef NitscheCoupling<distype, DRT::Element::quad8, 3> NitscheCouplType;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      case DRT::Element::quad9:
      {
        typedef NitscheCoupling<distype, DRT::Element::quad9, 3> NitscheCouplType;
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
      //      case DRT::Element::tri3:
      //      {
      //        typedef NitscheCoupling<distype,DRT::Element::tri3,3> NitscheCouplType;
      //        nit = new NitscheCouplType(bele_xyz,C_umum,rhC_um,);
      //        break;
      //      }
      //      case DRT::Element::tri6:
      //      {
      //        typedef NitscheCoupling<distype,DRT::Element::tri6,3> NitscheCouplType;
      //        nit = new NitscheCouplType(bele_xyz,C_umum,rhC_um,);
      //        break;
      //      }
      case DRT::Element::quad4:
      {
        typedef NitscheCoupling<distype, DRT::Element::quad4, 4> NitscheCouplType;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      case DRT::Element::quad8:
      {
        typedef NitscheCoupling<distype, DRT::Element::quad8, 4> NitscheCouplType;
        nit = new NitscheCouplType(bele_xyz, C_umum, rhC_um, fldparaxfem);
        break;
      }
      case DRT::Element::quad9:
      {
        typedef NitscheCoupling<distype, DRT::Element::quad9, 4> NitscheCouplType;
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
    Epetra_SerialDenseMatrix& bele_xyz, Epetra_SerialDenseMatrix& C_umum,
    Epetra_SerialDenseMatrix& C_usum, Epetra_SerialDenseMatrix& C_umus,
    Epetra_SerialDenseMatrix& C_usus, Epetra_SerialDenseMatrix& rhC_um,
    Epetra_SerialDenseMatrix& rhC_us, const DRT::ELEMENTS::FluidEleParameterXFEM& fldparaxfem)
{
  NitscheInterface* nit = NULL;

  // get number of dofs for this boundary element
  const unsigned numdofpernode = bele->NumDofPerNode(*bele->Nodes()[0]);

  // three dofs per node, for monolithic XFSI, four dofs per node for background-sided fluid-fluid
  // coupling
  if (numdofpernode == 3)
  {
    switch (bele->Shape())
    {
        //    case DRT::Element::tri3:
        //    {
        //      typedef NitscheCoupling<distype,DRT::Element::tri3,3> NitscheCouplType;
        //      nit = new
        //      NitscheCouplType(bele_xyz,C_umum,C_usum,C_umus,C_usus,rhC_um,rhC_us,is_viscAdjointSymmetric);
        //      break;
        //    }
        //    case DRT::Element::tri6:
        //    {
        //      typedef NitscheCoupling<distype,DRT::Element::tri6,3> NitscheCouplType;
        //      nit = new
        //      NitscheCouplType(bele_xyz,C_umum,C_usum,C_umus,C_usus,rhC_um,rhC_us,is_viscAdjointSymmetric);
        //      break;
        //    }
      case DRT::Element::quad4:
      {
        typedef NitscheCoupling<distype, DRT::Element::quad4, 3> NitscheCouplType;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case DRT::Element::quad8:
      {
        typedef NitscheCoupling<distype, DRT::Element::quad8, 3> NitscheCouplType;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case DRT::Element::quad9:
      {
        typedef NitscheCoupling<distype, DRT::Element::quad9, 3> NitscheCouplType;
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
      //    case DRT::Element::tri3:
      //    {
      //      typedef NitscheCoupling<distype,DRT::Element::tri3,4> NitscheCouplType;
      //      nit = new
      //      NitscheCouplType(bele_xyz,C_umum,C_usum,C_umus,C_usus,rhC_um,rhC_us,is_viscAdjointSymmetric);
      //      break;
      //    }
      //    case DRT::Element::tri6:
      //    {
      //      typedef NitscheCoupling<distype,DRT::Element::tri6,4> NitscheCouplType;
      //      nit = new
      //      NitscheCouplType(bele_xyz,C_umum,C_usum,C_umus,C_usus,rhC_um,rhC_us,is_viscAdjointSymmetric);
      //      break;
      //    }
      case DRT::Element::quad4:
      {
        typedef NitscheCoupling<distype, DRT::Element::quad4, 4> NitscheCouplType;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case DRT::Element::quad8:
      {
        typedef NitscheCoupling<distype, DRT::Element::quad8, 4> NitscheCouplType;
        nit = new NitscheCouplType(
            bele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case DRT::Element::quad9:
      {
        typedef NitscheCoupling<distype, DRT::Element::quad9, 4> NitscheCouplType;
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
    DRT::Element* vele, Epetra_SerialDenseMatrix& vele_xyz, Epetra_SerialDenseMatrix& C_umum,
    Epetra_SerialDenseMatrix& C_usum, Epetra_SerialDenseMatrix& C_umus,
    Epetra_SerialDenseMatrix& C_usus, Epetra_SerialDenseMatrix& rhC_um,
    Epetra_SerialDenseMatrix& rhC_us, const DRT::ELEMENTS::FluidEleParameterXFEM& fldparaxfem)
{
  NitscheInterface* nit = NULL;

  // get number of dofs for the embedded element
  const unsigned numdofpernode = vele->NumDofPerNode(*vele->Nodes()[0]);

  if (numdofpernode == 4)
  {
    switch (vele->Shape())
    {
        //    case DRT::Element::tet4:
        //    {
        //      typedef NitscheCoupling<distype,DRT::Element::tet4,4> NitscheCouplType;
        //      nit = new
        //      NitscheCouplType(vele_xyz,C_umum,C_usum,C_umus,C_usus,rhC_um,rhC_us,is_viscAdjointSymmetric);
        //      break;
        //    }
        //    case DRT::Element::tet10:
        //    {
        //      typedef NitscheCoupling<distype,DRT::Element::tet10,4> NitscheCouplType;
        //      nit = new
        //      NitscheCouplType(vele_xyz,C_umum,C_usum,C_umus,C_usus,rhC_um,rhC_us,is_viscAdjointSymmetric);
        //      break;
        //    }
      case DRT::Element::hex8:
      {
        typedef NitscheCoupling<distype, DRT::Element::hex8, 4> NitscheCouplType;
        nit = new NitscheCouplType(
            vele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case DRT::Element::hex20:
      {
        typedef NitscheCoupling<distype, DRT::Element::hex20, 4> NitscheCouplType;
        nit = new NitscheCouplType(
            vele_xyz, C_umum, C_usum, C_umus, C_usus, rhC_um, rhC_us, fldparaxfem);
        break;
      }
      case DRT::Element::hex27:
      {
        typedef NitscheCoupling<distype, DRT::Element::hex27, 4> NitscheCouplType;
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
      case DRT::Element::hex8:
      {
        typedef NitscheCoupling<distype, DRT::Element::hex8, 3> NitscheCouplType;
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
  HybridLMInterface* hybridlm = NULL;
  typedef HybridLMCoupling<distype, DRT::Element::dis_none, 3> HybridLMCouplType;
  hybridlm = new HybridLMCouplType(is_viscAdjointSymmetric);

  return Teuchos::rcp(hybridlm);
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
Teuchos::RCP<HybridLMInterface<distype>>
HybridLMInterface<distype>::CreateHybridLMCoupling_XFluidWDBC(
    DRT::Element* bele,                  ///< boundary element
    Epetra_SerialDenseMatrix& bele_xyz,  ///< global node coordinates of boundary element
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
      //    case DRT::Element::tri3:
      //    {
      //      typedef HybridLMCoupling<distype,DRT::Element::tri3,3> HybridLMCouplType;
      //      return Teuchos::rcp(new HybridLMCouplType(bele_xyz));
      //      break;
      //    }
      //    case DRT::Element::tri6:
      //    {
      //      typedef HybridLMCoupling<distype,DRT::Element::tri6,3> HybridLMCouplType;
      //      return Teuchos::rcp(new HybridLMCouplType(bele_xyz));
      //      break;
      //    }
    case DRT::Element::quad4:
    {
      typedef HybridLMCoupling<distype, DRT::Element::quad4, 3> HybridLMCouplType;
      return Teuchos::rcp(new HybridLMCouplType(bele_xyz, is_viscAdjointSymmetric));
      break;
    }
    case DRT::Element::quad8:
    {
      typedef HybridLMCoupling<distype, DRT::Element::quad8, 3> HybridLMCouplType;
      return Teuchos::rcp(new HybridLMCouplType(bele_xyz, is_viscAdjointSymmetric));
      break;
    }
    case DRT::Element::quad9:
    {
      typedef HybridLMCoupling<distype, DRT::Element::quad9, 3> HybridLMCouplType;
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
    DRT::Element* bele,                  ///< boundary element
    Epetra_SerialDenseMatrix& bele_xyz,  ///< global node coordinates of boundary element
    Epetra_SerialDenseMatrix& C_usum,    ///< C_usum coupling matrix
    Epetra_SerialDenseMatrix& C_umus,    ///< C_umus coupling matrix
    Epetra_SerialDenseMatrix& rhC_us,    ///< C_us coupling rhs
    Epetra_SerialDenseMatrix& G_s_us,    ///< \f$G_{u^s \sigma}\f$ coupling matrix
    Epetra_SerialDenseMatrix& G_us_s,    ///< \f$G_{\sigma u^s}\f$ coupling matrix
    bool is_viscAdjointSymmetric  ///< flag that indicates equal signs of Nitsche's standard &
                                  ///< adjoint viscous term
)
{
  HybridLMInterface* hlm = NULL;

  // get number of dofs for this boundary element
  const unsigned numdofpernode = bele->NumDofPerNode(*bele->Nodes()[0]);

  // three dofs per node, for monolithic XFSI, four dofs per node for backgrounds-sided fluid-fluid
  // coupling
  if (numdofpernode == 3)
  {
    switch (bele->Shape())
    {
        //    case DRT::Element::tri3:
        //    {
        //      typedef HybridLMCoupling<distype,DRT::Element::tri3,3> HybridLMCouplType;
        //      hlm = new
        //      HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,is_viscAdjointSymmetric);
        //      break;
        //    }
        //    case DRT::Element::tri6:
        //    {
        //      typedef HybridLMCoupling<distype,DRT::Element::tri6,3> HybridLMCouplType;
        //      hlm = new
        //      HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,is_viscAdjointSymmetric);
        //      break;
        //    }
      case DRT::Element::quad4:
      {
        typedef HybridLMCoupling<distype, DRT::Element::quad4, 3> HybridLMCouplType;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      case DRT::Element::quad8:
      {
        typedef HybridLMCoupling<distype, DRT::Element::quad8, 3> HybridLMCouplType;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      case DRT::Element::quad9:
      {
        typedef HybridLMCoupling<distype, DRT::Element::quad9, 3> HybridLMCouplType;
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
        //      case DRT::Element::tri3:
        //      {
        //        typedef HybridLMCoupling<distype,DRT::Element::tri3,4> HybridLMCouplType;
        //        hlm = new
        //        HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,is_viscAdjointSymmetric);
        //        break;
        //      }
        //      case DRT::Element::tri6:
        //      {
        //        typedef HybridLMCoupling<distype,DRT::Element::tri6,4> HybridLMCouplType;
        //        hlm = new
        //        HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,is_viscAdjointSymmetric);
        //        break;
        //      }
      case DRT::Element::quad4:
      {
        typedef HybridLMCoupling<distype, DRT::Element::quad4, 4> HybridLMCouplType;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      case DRT::Element::quad8:
      {
        typedef HybridLMCoupling<distype, DRT::Element::quad8, 4> HybridLMCouplType;
        hlm = new HybridLMCouplType(
            bele_xyz, C_usum, C_umus, rhC_us, G_s_us, G_us_s, is_viscAdjointSymmetric);
        break;
      }
      case DRT::Element::quad9:
      {
        typedef HybridLMCoupling<distype, DRT::Element::quad9, 4> HybridLMCouplType;
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

template class DRT::ELEMENTS::XFLUID::SlaveElementInterface<DRT::Element::hex8>;
template class DRT::ELEMENTS::XFLUID::SlaveElementInterface<DRT::Element::hex20>;
template class DRT::ELEMENTS::XFLUID::SlaveElementInterface<DRT::Element::hex27>;
template class DRT::ELEMENTS::XFLUID::SlaveElementInterface<DRT::Element::tet4>;
template class DRT::ELEMENTS::XFLUID::SlaveElementInterface<DRT::Element::tet10>;
template class DRT::ELEMENTS::XFLUID::SlaveElementInterface<DRT::Element::wedge6>;
template class DRT::ELEMENTS::XFLUID::SlaveElementInterface<DRT::Element::wedge15>;
// template class DRT::ELEMENTS::XFLUID::SlaveElementInterface<DRT::Element::pyramid5>;

template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::hex8>;
template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::hex20>;
template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::hex27>;
template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::tet4>;
template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::tet10>;
template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::wedge6>;
template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::wedge15>;
// template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::pyramid5>;

template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::hex8>;
template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::hex20>;
template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::hex27>;
template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::tet4>;
template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::tet10>;
template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::wedge6>;
template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::wedge15>;
// template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::pyramid5>;

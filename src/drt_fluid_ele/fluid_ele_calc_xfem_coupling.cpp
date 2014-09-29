/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_xfem_coupling.cpp

\brief Factory class for providing an implementation for coupling with
       Mixed/Stress/Hybrid methods or Nitsche's method to enforce
       interface conditions in the XFEM weakly
<pre>
Maintainer: Raffaela Kruse /Benedikt Schott
            kruse@lnm.mw.tum.de
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>

#include <fstream>

#include "../drt_bele3/bele3.H"

#include "fluid_ele_calc_xfem_coupling.H"
#include "fluid_ele_calc_xfem_coupling_impl.H"

using namespace DRT::ELEMENTS::XFLUID;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
Teuchos::RCP<SlaveElementInterface<distype> > SlaveElementInterface<distype>::CreateSlaveElementRepresentation(
  DRT::Element *              slave_ele,       ///< coupling slave element
  Epetra_SerialDenseMatrix &  slave_xyz        ///< global node coordinates of coupling slave element
)
{
  SlaveElementInterface * sla = NULL;

  // get number of dofs for this slave element
  const unsigned numdofpernode = slave_ele->NumDofPerNode(*slave_ele->Nodes()[0]);

  if (numdofpernode == 3)
  {
    switch ( slave_ele->Shape() )
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
        typedef SlaveElementRepresentation<distype,DRT::Element::quad4,3> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::quad8:
      {
        typedef SlaveElementRepresentation<distype,DRT::Element::quad8,3> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::quad9:
      {
        typedef SlaveElementRepresentation<distype,DRT::Element::quad9,3> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      default:
        dserror( "Unsupported boundary element shape %d", slave_ele->Shape() ); break;
    }
  }
  else if (numdofpernode == 4) // volumetric coupling partners only required for fluid-fluid coupling
  {

    switch ( slave_ele->Shape() )
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
        typedef SlaveElementRepresentation<distype,DRT::Element::quad4,4> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::quad8:
      {
        typedef SlaveElementRepresentation<distype,DRT::Element::quad8,4> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::quad9:
      {
        typedef SlaveElementRepresentation<distype,DRT::Element::quad9,4> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::hex8:
      {
        typedef SlaveElementRepresentation<distype,DRT::Element::hex8,4> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::hex20:
      {
        typedef SlaveElementRepresentation<distype,DRT::Element::hex20,4> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      case DRT::Element::hex27:
      {
        typedef SlaveElementRepresentation<distype,DRT::Element::hex27,4> SlaveEleType;
        sla = new SlaveEleType(slave_xyz);
        break;
      }
      default:
        dserror( "unsupported boundary element shape %d", slave_ele->Shape() ); break;
    }
  }

  return Teuchos::rcp(sla);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
Teuchos::RCP<NitscheInterface<distype> > NitscheInterface<distype>::CreateNitscheCoupling_XFluidWDBC(
  bool  is_viscAdjointSymmetric ///< flag that indicates equal signs of Nitsche's standard & adjoint viscous term
)
{
  NitscheInterface * nit = NULL;
  typedef NitscheCoupling<distype,DRT::Element::dis_none,3> NitscheCouplType;
  nit = new NitscheCouplType(is_viscAdjointSymmetric);

  return Teuchos::rcp(nit);
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
Teuchos::RCP<NitscheInterface<distype> > NitscheInterface<distype>::CreateNitscheCoupling_XFluidWDBC(
  DRT::Element *              bele,           ///< boundary element
  Epetra_SerialDenseMatrix &  bele_xyz,       ///< global node coordinates of boundary element
  bool                        is_viscAdjointSymmetric ///< flag that indicates equal signs of Nitsche's standard & adjoint viscous term
)
{
  NitscheInterface * nit = NULL;

  // get number of dofs for this boundary element
  const unsigned numdofpernode = bele->NumDofPerNode(*bele->Nodes()[0]);

  // three dofs per node, for standard Dirichlet coupling, four dofs per node for fluid-fluid coupling
  if (numdofpernode != 3)
  {
    dserror("Unsupported number of %d nodes for standard Dirichlet coupling.", numdofpernode);
  }

  switch ( bele->Shape() )
  {
//      case DRT::Element::tri3:
//      {
//        typedef NitscheCoupling<distype,DRT::Element::tri3,3> NitscheCouplType;
//        nit = new NitscheCouplType(slave_xyz);
//        break;
//      }
//      case DRT::Element::tri6:
//      {
//        typedef NitscheCoupling<distype,DRT::Element::tri6,3> NitscheCouplType;
//        nit = new NitscheCouplType(slave_xyz);
//        break;
//      }
    case DRT::Element::quad4:
    {
      typedef NitscheCoupling<distype,DRT::Element::quad4,3> NitscheCouplType;
      nit = new NitscheCouplType(bele_xyz,is_viscAdjointSymmetric);
      break;
    }
    case DRT::Element::quad8:
    {
      typedef NitscheCoupling<distype,DRT::Element::quad8,3> NitscheCouplType;
      nit = new NitscheCouplType(bele_xyz,is_viscAdjointSymmetric);
      break;
    }
    case DRT::Element::quad9:
    {
      typedef NitscheCoupling<distype,DRT::Element::quad9,3> NitscheCouplType;
      nit = new NitscheCouplType(bele_xyz,is_viscAdjointSymmetric);
      break;
    }
    default:
      dserror( "Unsupported boundary element shape %d", bele->Shape() ); break;
  }

  return Teuchos::rcp(nit);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
Teuchos::RCP<NitscheInterface<distype> > NitscheInterface<distype>::CreateNitscheCoupling_XFluidSided(
  DRT::Element *              bele,                   ///< boundary element
  Epetra_SerialDenseMatrix &  bele_xyz,               ///< global node coordinates of boundary element
  Epetra_SerialDenseMatrix &  C_usum,                 ///< C_usum coupling matrix
  Epetra_SerialDenseMatrix &  C_umus,                 ///< C_umus coupling matrix
  Epetra_SerialDenseMatrix &  rhC_us,                 ///< C_us coupling rhs
  Epetra_SerialDenseMatrix &  C_usus,                 ///< C_usus coupling matrix
  bool                        is_viscAdjointSymmetric ///< flag that indicates equal signs of Nitsche's standard & adjoint viscous term
)
{
  NitscheInterface * nit = NULL;

  // get number of dofs for this boundary element
  const unsigned numdofpernode = bele->NumDofPerNode(*bele->Nodes()[0]);

  // three dofs per node, for monolithic XFSI, four dofs per node for background-sided fluid-fluid coupling
  if ( numdofpernode == 3 )
  {
    switch ( bele->Shape() )
    {
//    case DRT::Element::tri3:
//    {
//      typedef NitscheCoupling<distype,DRT::Element::tri3,3> NitscheCouplType;
//      nit = new NitscheCouplType(bele_xyz,C_usum,C_umus,rhC_us,C_usus,
//          SlaveElementInterface<distype>::MonolithicXFSI,is_viscAdjointSymmetric);
//      break;
//    }
//    case DRT::Element::tri6:
//    {
//      typedef NitscheCoupling<distype,DRT::Element::tri6,3> NitscheCouplType;
//      nit = new NitscheCouplType(bele_xyz,C_usum,C_umus,rhC_us,C_usus,
//          SlaveElementInterface<distype>::MonolithicXFSI,is_viscAdjointSymmetric);
//      break;
//    }
    case DRT::Element::quad4:
    {
      typedef NitscheCoupling<distype,DRT::Element::quad4,3> NitscheCouplType;
      nit = new NitscheCouplType(bele_xyz,C_usum,C_umus,rhC_us,C_usus,
          SlaveElementInterface<distype>::MonolithicXFSI,is_viscAdjointSymmetric);
      break;
    }
    case DRT::Element::quad8:
    {
      typedef NitscheCoupling<distype,DRT::Element::quad8,3> NitscheCouplType;
      nit = new NitscheCouplType(bele_xyz,C_usum,C_umus,rhC_us,C_usus,
          SlaveElementInterface<distype>::MonolithicXFSI,is_viscAdjointSymmetric);
      break;
    }
    case DRT::Element::quad9:
    {
      typedef NitscheCoupling<distype,DRT::Element::quad9,3> NitscheCouplType;
      nit = new NitscheCouplType(bele_xyz,C_usum,C_umus,rhC_us,C_usus,
          SlaveElementInterface<distype>::MonolithicXFSI,is_viscAdjointSymmetric);
      break;
    }
    default:
      dserror( "Unsupported boundary element shape %d", bele->Shape() ); break;
    }
  }
  else if ( numdofpernode == 4 )
  {
    switch ( bele->Shape() )
    {
    //    case DRT::Element::tri3:
    //    {
    //      typedef NitscheCoupling<distype,DRT::Element::tri3,4> NitscheCouplType;
    //      nit = new NitscheCouplType(bele_xyz,C_usum,C_umus,rhC_us,C_usus,
    //          SlaveElementInterface<distype>::MonolithicXFluidFluid,is_viscAdjointSymmetric);
    //      break;
    //    }
    //    case DRT::Element::tri6:
    //    {
    //      typedef NitscheCoupling<distype,DRT::Element::tri6,4> NitscheCouplType;
    //      nit = new NitscheCouplType(bele_xyz,C_usum,C_umus,rhC_us,C_usus,
    //          SlaveElementInterface<distype>::MonolithicXFluidFluid,is_viscAdjointSymmetric);
    //      break;
    //    }
    case DRT::Element::quad4:
    {
      typedef NitscheCoupling<distype,DRT::Element::quad4,4> NitscheCouplType;
      nit = new NitscheCouplType(bele_xyz,C_usum,C_umus,rhC_us,C_usus,
          SlaveElementInterface<distype>::XFluidFluid,is_viscAdjointSymmetric);
      break;
    }
    case DRT::Element::quad8:
    {
      typedef NitscheCoupling<distype,DRT::Element::quad8,4> NitscheCouplType;
      nit = new NitscheCouplType(bele_xyz,C_usum,C_umus,rhC_us,C_usus,
          SlaveElementInterface<distype>::XFluidFluid,is_viscAdjointSymmetric);
      break;
    }
    case DRT::Element::quad9:
    {
      typedef NitscheCoupling<distype,DRT::Element::quad9,4> NitscheCouplType;
      nit = new NitscheCouplType(bele_xyz,C_usum,C_umus,rhC_us,C_usus,
          SlaveElementInterface<distype>::XFluidFluid,is_viscAdjointSymmetric);
      break;
    }
    default:
      dserror( "Unsupported boundary element shape %d", bele->Shape() ); break;
    }
  }
  else
    dserror("Unsupported number of %d nodes for coupling slave element.", numdofpernode);

  return Teuchos::rcp(nit);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
Teuchos::RCP<NitscheInterface<distype> > NitscheInterface<distype>::CreateNitscheCoupling_TwoSided(
  DRT::Element *              vele,           ///< volumetric element to couple with
  Epetra_SerialDenseMatrix &  vele_xyz,       ///< global node coordinates of volumetric element
  Epetra_SerialDenseMatrix &  C_usum,         ///< C_usum coupling matrix
  Epetra_SerialDenseMatrix &  C_umus,         ///< C_umus coupling matrix
  Epetra_SerialDenseMatrix &  rhC_us,         ///< C_us coupling rhs
  Epetra_SerialDenseMatrix &  C_usus,         ///< C_usus coupling matrix
  bool                        is_viscAdjointSymmetric ///< flag that indicates equal signs of Nitsche's standard & adjoint viscous term
)
{
  NitscheInterface * nit = NULL;

  // get number of dofs for the embedded element
  const unsigned numdofpernode = vele->NumDofPerNode(*vele->Nodes()[0]);

  // expecting 4 dofs per slave element node as this is fluid-fluid coupling
  if (numdofpernode != 4)
  {
    dserror("Unsupported number of %d nodes for fluid-fluid coupling of slave element.", numdofpernode);
    return Teuchos::null;
  }

  switch ( vele->Shape() )
  {
//    case DRT::Element::tet4:
//    {
//      typedef NitscheCoupling<distype,DRT::Element::tet4,4> NitscheCouplType;
//      nit = new NitscheCouplType(vele_xyz,C_usum,C_umus,rhC_us,C_usus,
//          SlaveElementInterface<distype>::XFluidFluid,is_viscAdjointSymmetric);
//      break;
//    }
//    case DRT::Element::tet10:
//    {
//      typedef NitscheCoupling<distype,DRT::Element::tet10,4> NitscheCouplType;
//      nit = new NitscheCouplType(vele_xyz,C_usum,C_umus,rhC_us,C_usus,
//          SlaveElementInterface<distype>::XFluidFluid,is_viscAdjointSymmetric);
//      break;
//    }
    case DRT::Element::hex8:
    {
      typedef NitscheCoupling<distype,DRT::Element::hex8,4> NitscheCouplType;
      nit = new NitscheCouplType(vele_xyz,C_usum,C_umus,rhC_us,C_usus,
          SlaveElementInterface<distype>::XFluidFluid,is_viscAdjointSymmetric);
      break;
    }
    case DRT::Element::hex20:
    {
      typedef NitscheCoupling<distype,DRT::Element::hex20,4> NitscheCouplType;
      nit = new NitscheCouplType(vele_xyz,C_usum,C_umus,rhC_us,C_usus,
          SlaveElementInterface<distype>::XFluidFluid,is_viscAdjointSymmetric);
      break;
    }
    case DRT::Element::hex27:
    {
      typedef NitscheCoupling<distype,DRT::Element::hex27,4> NitscheCouplType;
      nit = new NitscheCouplType(vele_xyz,C_usum,C_umus,rhC_us,C_usus,
          SlaveElementInterface<distype>::XFluidFluid,is_viscAdjointSymmetric);
      break;
    }
    default:
      dserror( "Unsupported volume element shape %d.", vele->Shape() ); break;
  }

  return Teuchos::rcp(nit);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
Teuchos::RCP<HybridLMInterface<distype> > HybridLMInterface<distype>::CreateHybridLMCoupling_XFluidWDBC(
  DRT::Element *              bele,           ///< boundary element
  Epetra_SerialDenseMatrix &  bele_xyz,       ///< global node coordinates of boundary element
  bool                        is_viscAdjointSymmetric ///< flag that indicates equal signs of Nitsche's standard & adjoint viscous term
)
{
  // get number of dofs for this boundary element
  const unsigned numdofpernode = bele->NumDofPerNode(*bele->Nodes()[0]);

  // three dofs per node, for standard Dirichlet coupling, four dofs per node for fluid-fluid coupling
  if (numdofpernode != 3)
  {
    dserror("Unsupported number of %d nodes for standard Dirichlet coupling.", numdofpernode);
  }

  switch ( bele->Shape() )
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
      typedef HybridLMCoupling<distype,DRT::Element::quad4,3> HybridLMCouplType;
      return Teuchos::rcp(new HybridLMCouplType(bele_xyz,is_viscAdjointSymmetric));
      break;
    }
    case DRT::Element::quad8:
    {
      typedef HybridLMCoupling<distype,DRT::Element::quad8,3> HybridLMCouplType;
      return Teuchos::rcp(new HybridLMCouplType(bele_xyz,is_viscAdjointSymmetric));
      break;
    }
    case DRT::Element::quad9:
    {
      typedef HybridLMCoupling<distype,DRT::Element::quad9,3> HybridLMCouplType;
      return Teuchos::rcp(new HybridLMCouplType(bele_xyz,is_viscAdjointSymmetric));
      break;
    }
    default:
      dserror( "Unsupported boundary element shape %d", bele->Shape() ); break;
  }

  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
Teuchos::RCP<HybridLMInterface<distype> > HybridLMInterface<distype>::CreateHybridLMCoupling_XFluidSided(
  DRT::Element *              bele,           ///< boundary element
  Epetra_SerialDenseMatrix &  bele_xyz,       ///< global node coordinates of boundary element
  Epetra_SerialDenseMatrix &  C_usum,         ///< C_usum coupling matrix
  Epetra_SerialDenseMatrix &  C_umus,         ///< C_umus coupling matrix
  Epetra_SerialDenseMatrix &  rhC_us,         ///< C_us coupling rhs
  Epetra_SerialDenseMatrix &  G_s_us,         ///< \f$G_{u^s \sigma}\f$ coupling matrix
  Epetra_SerialDenseMatrix &  G_us_s,         ///< \f$G_{\sigma u^s}\f$ coupling matrix
  bool                        is_viscAdjointSymmetric ///< flag that indicates equal signs of Nitsche's standard & adjoint viscous term
)
{
  HybridLMInterface * hlm = NULL;

  // get number of dofs for this boundary element
  const unsigned numdofpernode = bele->NumDofPerNode(*bele->Nodes()[0]);

  // three dofs per node, for monolithic XFSI, four dofs per node for backgrounds-sided fluid-fluid coupling
  if ( numdofpernode == 3 )
  {
    switch ( bele->Shape() )
    {
//    case DRT::Element::tri3:
//    {
//      typedef HybridLMCoupling<distype,DRT::Element::tri3,3> HybridLMCouplType;
//      hlm = new HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,
//          SlaveElementInterface<distype>::MonolithicXFSI,is_viscAdjointSymmetric);
//      break;
//    }
//    case DRT::Element::tri6:
//    {
//      typedef HybridLMCoupling<distype,DRT::Element::tri6,3> HybridLMCouplType;
//      hlm = new HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,
//          SlaveElementInterface<distype>::MonolithicXFSI,is_viscAdjointSymmetric);
//      break;
//    }
    case DRT::Element::quad4:
    {
      typedef HybridLMCoupling<distype,DRT::Element::quad4,3> HybridLMCouplType;
      hlm = new HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,
          SlaveElementInterface<distype>::MonolithicXFSI,is_viscAdjointSymmetric);
      break;
    }
    case DRT::Element::quad8:
    {
      typedef HybridLMCoupling<distype,DRT::Element::quad8,3> HybridLMCouplType;
      hlm = new HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,
          SlaveElementInterface<distype>::MonolithicXFSI,is_viscAdjointSymmetric);
      break;
    }
    case DRT::Element::quad9:
    {
      typedef HybridLMCoupling<distype,DRT::Element::quad9,3> HybridLMCouplType;
      hlm = new HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,
          SlaveElementInterface<distype>::MonolithicXFSI,is_viscAdjointSymmetric);
      break;
    }
    default:
      dserror( "Unsupported boundary element shape %d", bele->Shape() ); break;
    }
  }
  else if ( numdofpernode == 4 )
  {
    switch ( bele->Shape() )
    {
//      case DRT::Element::tri3:
//      {
//        typedef HybridLMCoupling<distype,DRT::Element::tri3,4> HybridLMCouplType;
//        hlm = new HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,
//            SlaveElementInterface<distype>::XFluidFluid,is_viscAdjointSymmetric);
//        break;
//      }
//      case DRT::Element::tri6:
//      {
//        typedef HybridLMCoupling<distype,DRT::Element::tri6,4> HybridLMCouplType;
//        hlm = new HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,
//            SlaveElementInterface<distype>::XFluidFluid,is_viscAdjointSymmetric);
//        break;
//      }
      case DRT::Element::quad4:
      {
        typedef HybridLMCoupling<distype,DRT::Element::quad4,4> HybridLMCouplType;
        hlm = new HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,
            SlaveElementInterface<distype>::XFluidFluid,is_viscAdjointSymmetric);
        break;
      }
      case DRT::Element::quad8:
      {
        typedef HybridLMCoupling<distype,DRT::Element::quad8,4> HybridLMCouplType;
        hlm = new HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,
            SlaveElementInterface<distype>::XFluidFluid,is_viscAdjointSymmetric);
        break;
      }
      case DRT::Element::quad9:
      {
        typedef HybridLMCoupling<distype,DRT::Element::quad9,4> HybridLMCouplType;
        hlm = new HybridLMCouplType(bele_xyz,C_usum,C_umus,rhC_us,G_s_us,G_us_s,
            SlaveElementInterface<distype>::XFluidFluid,is_viscAdjointSymmetric);
        break;
      }
      default:
        dserror( "Unsupported boundary element shape %d", bele->Shape() ); break;
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
//template class DRT::ELEMENTS::XFLUID::SlaveElementInterface<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::XFLUID::SlaveElementInterface<DRT::Element::pyramid5>;

template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::hex8>;
template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::hex20>;
template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::hex27>;
template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::tet4>;
template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::tet10>;
//template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::XFLUID::NitscheInterface<DRT::Element::pyramid5>;

template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::hex8>;
template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::hex20>;
template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::hex27>;
template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::tet4>;
template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::tet10>;
//template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::XFLUID::HybridLMInterface<DRT::Element::pyramid5>;

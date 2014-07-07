/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_xfem_coupling.cpp

\brief Interface class for coupling of sides and elements from two different meshes

<pre>
Maintainer: Shadan Shahmiri /Benedikt Schott
            shahmiri@lnm.mw.tum.de
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



/*--------------------------------------------------------------------------------
 * Mixed/Stress/Hybrid (MSH) coupling
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
Teuchos::RCP<SideInterface<distype> > SideInterface<distype>::Impl(
    DRT::Element *              side,            ///< side element
    Epetra_SerialDenseMatrix &  C_uiu,           ///< C_uiu coupling matrix
    Epetra_SerialDenseMatrix &  C_uui,           ///< C_uui coupling matrix
    Epetra_SerialDenseMatrix &  rhC_ui,          ///< C_ui coupling rhs
    Epetra_SerialDenseMatrix &  Gsui,            ///< interface sigma-u_interface coupling
    Epetra_SerialDenseMatrix &  Guis,            ///< interface u_interface-sigma coupling
    Epetra_SerialDenseMatrix &  side_xyze,       ///< global node coordinates
    bool isViscAdjointSymmetric                  ///< symmetric or skew-symmetric formulation of MHVS adjoint viscous term
)
{
  SideInterface * si = NULL;

  // get number of dofs for this side element
  const int numdofpernode = side->NumDofPerNode(*side->Nodes()[0]);

  if (numdofpernode == 3) // three dofs per node, for standard Dirichlet coupling
  {
    switch ( side->Shape() )
    {
//    case DRT::Element::tri3:
//    {
//      typedef SideImpl<distype,DRT::Element::tri3, 3> SideImplType;
//      si = new SideImplType(side,C_uiu,C_uui,rhC_ui,Gsui,Guis,side_xyze,isViscAdjointSymmetric);
//      break;
//    }
//    case DRT::Element::tri6:
//    {
//      typedef SideImpl<distype,DRT::Element::tri6, 3> SideImplType;
//      si = new SideImplType(side,C_uiu,C_uui,rhC_ui,Gsui,Guis,side_xyze,isViscAdjointSymmetric);
//      break;
//    }
    case DRT::Element::quad4:
    {
      typedef SideImpl<distype,DRT::Element::quad4, 3> SideImplType;
      si = new SideImplType(side,C_uiu,C_uui,rhC_ui,Gsui,Guis,side_xyze,isViscAdjointSymmetric);
      break;
    }
    case DRT::Element::quad8:
    {
      typedef SideImpl<distype,DRT::Element::quad8, 3> SideImplType;
      si = new SideImplType(side,C_uiu,C_uui,rhC_ui,Gsui,Guis,side_xyze,isViscAdjointSymmetric);
      break;
    }
    case DRT::Element::quad9:
    {
      typedef SideImpl<distype,DRT::Element::quad9, 3> SideImplType;
      si = new SideImplType(side,C_uiu,C_uui,rhC_ui,Gsui,Guis,side_xyze,isViscAdjointSymmetric);
      break;
    }
    default:
      dserror( "unsupported side shape %d", side->Shape() ); break;
    }
  }
  else if (numdofpernode == 4) // four dofs per node, for standard Dirichlet coupling
  {
    switch ( side->Shape() )
    {
//    case DRT::Element::tri3:
//    {
//      typedef SideImpl<distype,DRT::Element::tri3, 4> SideImplType;
//      si = new SideImplType(side,C_uiu,C_uui,rhC_ui,Gsui,Guis,side_xyze,isViscAdjointSymmetric);
//      break;
//    }
//    case DRT::Element::tri6:
//    {
//      typedef SideImpl<distype,DRT::Element::tri6, 4> SideImplType;
//      si = new SideImplType(side,C_uiu,C_uui,rhC_ui,Gsui,Guis,side_xyze,isViscAdjointSymmetric);
//      break;
//    }
    case DRT::Element::quad4:
    {
      typedef SideImpl<distype,DRT::Element::quad4, 4> SideImplType;
      si = new SideImplType(side,C_uiu,C_uui,rhC_ui,Gsui,Guis,side_xyze,isViscAdjointSymmetric);
      break;
    }
    case DRT::Element::quad8:
    {
      typedef SideImpl<distype,DRT::Element::quad8, 4> SideImplType;
      si = new SideImplType(side,C_uiu,C_uui,rhC_ui,Gsui,Guis,side_xyze,isViscAdjointSymmetric);
      break;
    }
    case DRT::Element::quad9:
    {
      typedef SideImpl<distype,DRT::Element::quad9, 4> SideImplType;
      si = new SideImplType(side,C_uiu,C_uui,rhC_ui,Gsui,Guis,side_xyze,isViscAdjointSymmetric);
      break;
    }
    default:
      dserror( "unsupported side shape %d", side->Shape() ); break;
    }
  }

  return Teuchos::rcp(si);
}


/*--------------------------------------------------------------------------------
 * Nitsche (NIT) coupling between background element and side element
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
Teuchos::RCP<SideInterface<distype> > SideInterface<distype>::Impl(
    DRT::Element *              side,            ///< side element
    Epetra_SerialDenseMatrix &  C_uiu,           ///< C_uiu coupling matrix
    Epetra_SerialDenseMatrix &  C_uui,           ///< C_uui coupling matrix
    Epetra_SerialDenseMatrix &  rhC_ui,          ///< C_ui coupling rhs
    Epetra_SerialDenseMatrix &  C_uiui,          ///< Cuiui coupling matrix
    Epetra_SerialDenseMatrix &  side_xyze,       ///< global node coordinates
    bool isViscAdjointSymmetric                  ///< symmetric or skew-symmetric formulation of Nitsche's adjoint viscous term
)
{
  SideInterface * si = NULL;

  // get number of dofs for this side element
  const int numdofpernode = side->NumDofPerNode(*side->Nodes()[0]);

  if (numdofpernode == 3) // three dofs per node, for coupling
  {
    switch ( side->Shape() )
    {
//    case DRT::Element::tri3:
//    {
//      typedef SideImpl<distype,DRT::Element::tri3, 3> SideImplType;
//      si = new SideImplType(side,C_uiu,C_uui,rhC_ui,C_uiui,side_xyze,isViscAdjointSymmetric);
//      break;
//    }
//    case DRT::Element::tri6:
//    {
//      typedef SideImpl<distype,DRT::Element::tri6, 3> SideImplType;
//      si = new SideImplType(side,C_uiu,C_uui,rhC_ui,C_uiui,side_xyze,isViscAdjointSymmetric);
//      break;
//    }
    case DRT::Element::quad4:
    {
      typedef SideImpl<distype,DRT::Element::quad4, 3> SideImplType;
      si = new SideImplType(side,C_uiu,C_uui,rhC_ui,C_uiui,side_xyze,isViscAdjointSymmetric);
      break;
    }
    case DRT::Element::quad8:
    {
      typedef SideImpl<distype,DRT::Element::quad8, 3> SideImplType;
      si = new SideImplType(side,C_uiu,C_uui,rhC_ui,C_uiui,side_xyze,isViscAdjointSymmetric);
      break;
    }
    case DRT::Element::quad9:
    {
      typedef SideImpl<distype,DRT::Element::quad9, 3> SideImplType;
      si = new SideImplType(side,C_uiu,C_uui,rhC_ui,C_uiui,side_xyze,isViscAdjointSymmetric);
      break;
    }
    default:
      dserror( "unsupported side shape %d", side->Shape() ); break;
    }
  }
  else if (numdofpernode == 4) // four dofs per node, for coupling
  {
    switch ( side->Shape() )
    {
//    case DRT::Element::tri3:
//    {
//      typedef SideImpl<distype,DRT::Element::tri3, 4> SideImplType;
//      si = new SideImplType(side,C_uiu,C_uui,rhC_ui,C_uiui,side_xyze,isViscAdjointSymmetric);
//      break;
//    }
//    case DRT::Element::tri6:
//    {
//      typedef SideImpl<distype,DRT::Element::tri6, 4> SideImplType;
//      si = new SideImplType(side,C_uiu,C_uui,rhC_ui,C_uiui,side_xyze,isViscAdjointSymmetric);
//      break;
//    }
    case DRT::Element::quad4:
    {
      typedef SideImpl<distype,DRT::Element::quad4, 4> SideImplType;
      si = new SideImplType(side,C_uiu,C_uui,rhC_ui,C_uiui,side_xyze,isViscAdjointSymmetric);
      break;
    }
    case DRT::Element::quad8:
    {
      typedef SideImpl<distype,DRT::Element::quad8, 4> SideImplType;
      si = new SideImplType(side,C_uiu,C_uui,rhC_ui,C_uiui,side_xyze,isViscAdjointSymmetric);
      break;
    }
    case DRT::Element::quad9:
    {
      typedef SideImpl<distype,DRT::Element::quad9, 4> SideImplType;
      si = new SideImplType(side,C_uiu,C_uui,rhC_ui,C_uiui,side_xyze,isViscAdjointSymmetric);
      break;
    }
    default:
      dserror( "unsupported side shape %d", side->Shape() ); break;
    }
  }

  return Teuchos::rcp(si);
}


/*--------------------------------------------------------------------------------
 * simple side coupling together with embedded Nitsche coupling
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
Teuchos::RCP<SideInterface<distype> > SideInterface<distype>::Impl(
    DRT::Element *              side,            ///< side element
    Epetra_SerialDenseMatrix &  side_xyze        ///< global node coordinates
)
{
  SideInterface * si = NULL;

  // get number of dofs for this side element
  const int numdofpernode = side->NumDofPerNode(*side->Nodes()[0]);

  if (numdofpernode == 3) // three dofs per node, for standard Dirichlet coupling
  {
    switch ( side->Shape() )
    {
//    case DRT::Element::tri3:
//    {
//      typedef SideImpl<distype,DRT::Element::tri3, 3> SideImplType;
//      si = new SideImplType(side,side_xyze);
//      break;
//    }
//    case DRT::Element::tri6:
//    {
//      typedef SideImpl<distype,DRT::Element::tri6, 3> SideImplType;
//      si = new SideImplType(side,side_xyze);
//      break;
//    }
    case DRT::Element::quad4:
    {
      typedef SideImpl<distype,DRT::Element::quad4, 3> SideImplType;
      si = new SideImplType(side,side_xyze);
      break;
    }
    case DRT::Element::quad8:
    {
      typedef SideImpl<distype,DRT::Element::quad8, 3> SideImplType;
      si = new SideImplType(side,side_xyze);
      break;
    }
    case DRT::Element::quad9:
    {
      typedef SideImpl<distype,DRT::Element::quad9, 3> SideImplType;
      si = new SideImplType(side,side_xyze);
      break;
    }
    default:
      dserror( "unsupported side shape %d", side->Shape() ); break;
    }
  }
  else if (numdofpernode == 4) // three dofs per node, for standard Dirichlet coupling
  {
    switch ( side->Shape() )
    {
//    case DRT::Element::tri3:
//    {
//      typedef SideImpl<distype,DRT::Element::tri3, 4> SideImplType;
//      si = new SideImplType(side,side_xyze);
//      break;
//    }
//    case DRT::Element::tri6:
//    {
//      typedef SideImpl<distype,DRT::Element::tri6, 4> SideImplType;
//      si = new SideImplType(side,side_xyze);
//      break;
//    }
    case DRT::Element::quad4:
    {
      typedef SideImpl<distype,DRT::Element::quad4, 4> SideImplType;
      si = new SideImplType(side,side_xyze);
      break;
    }
    case DRT::Element::quad8:
    {
      typedef SideImpl<distype,DRT::Element::quad8, 4> SideImplType;
      si = new SideImplType(side,side_xyze);
      break;
    }
    case DRT::Element::quad9:
    {
      typedef SideImpl<distype,DRT::Element::quad9, 4> SideImplType;
      si = new SideImplType(side,side_xyze);
      break;
    }
    default:
      dserror( "unsupported side shape %d", side->Shape() ); break;
    }
  }

  return Teuchos::rcp(si);
}


/*--------------------------------------------------------------------------------
 * Nitsche (NIT) coupling between background element and embedded element
 *--------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
Teuchos::RCP<EmbCoupling<distype> > EmbCoupling<distype>::TwoSidedImpl(
    DRT::Element *              emb_ele,         ///< side element
    Epetra_SerialDenseMatrix &  C_uiu,           ///< C_uiu coupling matrix
    Epetra_SerialDenseMatrix &  C_uui,           ///< C_uui coupling matrix
    Epetra_SerialDenseMatrix &  rhC_ui,          ///< C_ui coupling rhs
    Epetra_SerialDenseMatrix &  C_uiui,          ///< Cuiui coupling matrix
    Epetra_SerialDenseMatrix &  emb_xyze,        ///< global node coordinates
    bool isViscAdjointSymmetric                  ///< symmetric or skew-symmetric formulation of Nitsche's adjoint viscous term
)
{

  EmbCoupling * emb = NULL;

  switch ( emb_ele->Shape() )
  {
//  case DRT::Element::tet4:
//  {
//    typedef EmbImpl<distype,DRT::Element::tet4> EmbImplType;
//    emb = new EmbImplType(emb_ele,C_uiu,C_uui,rhC_ui,C_uiui,emb_xyze,isViscAdjointSymmetric);
//    break;
//  }
//  case DRT::Element::tet10:
//  {
//    typedef EmbImpl<distype,DRT::Element::tet10> EmbImplType;
//    emb = new EmbImplType(emb_ele,C_uiu,C_uui,rhC_ui,C_uiui,emb_xyze,isViscAdjointSymmetric);
//    break;
//  }
  case DRT::Element::hex8:
  {
    typedef EmbImpl<distype,DRT::Element::hex8> EmbImplType;
    emb = new EmbImplType(emb_ele,C_uiu,C_uui,rhC_ui,C_uiui,emb_xyze,isViscAdjointSymmetric);
    break;
  }
  case DRT::Element::hex20:
  {
    typedef EmbImpl<distype,DRT::Element::hex20> EmbImplType;
    emb = new EmbImplType(emb_ele,C_uiu,C_uui,rhC_ui,C_uiui,emb_xyze,isViscAdjointSymmetric);
    break;
  }
  case DRT::Element::hex27:
  {
    typedef EmbImpl<distype,DRT::Element::hex27> EmbImplType;
    emb = new EmbImplType(emb_ele,C_uiu,C_uui,rhC_ui,C_uiui,emb_xyze,isViscAdjointSymmetric);
    break;
  }
  default:
    dserror( "unsupported side shape %d", emb_ele->Shape() ); break;
  }
  return Teuchos::rcp(emb);

}





// create coupling objects just for 3D elements
template class DRT::ELEMENTS::XFLUID::SideInterface<DRT::Element::hex8>;
template class DRT::ELEMENTS::XFLUID::SideInterface<DRT::Element::hex20>;
template class DRT::ELEMENTS::XFLUID::SideInterface<DRT::Element::hex27>;
template class DRT::ELEMENTS::XFLUID::SideInterface<DRT::Element::tet4>;
template class DRT::ELEMENTS::XFLUID::SideInterface<DRT::Element::tet10>;
//template class DRT::ELEMENTS::XFLUID::SideInterface<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::XFLUID::SideInterface<DRT::Element::pyramid5>;


// create coupling objects just for 3D elements
template class DRT::ELEMENTS::XFLUID::EmbCoupling<DRT::Element::hex8>;
template class DRT::ELEMENTS::XFLUID::EmbCoupling<DRT::Element::hex20>;
template class DRT::ELEMENTS::XFLUID::EmbCoupling<DRT::Element::hex27>;
template class DRT::ELEMENTS::XFLUID::EmbCoupling<DRT::Element::tet4>;
template class DRT::ELEMENTS::XFLUID::EmbCoupling<DRT::Element::tet10>;
//template class DRT::ELEMENTS::XFLUID::EmbCoupling<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::XFLUID::EmbCoupling<DRT::Element::pyramid5>;



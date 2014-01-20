/*!
\file acou_ele_boundary_calc.cpp
\brief

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
*/
/*--------------------------------------------------------------------------*/


#include "acou_ele_boundary_calc.H"
#include "acou_ele.H"
#include "../drt_lib/drt_node.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::AcouBoundaryImplInterface* DRT::ELEMENTS::AcouBoundaryImplInterface::Impl(const DRT::Element* ele)
{
  switch (ele->Shape())
  {
  case DRT::Element::quad4:
  {
    return AcouBoundaryImpl<DRT::Element::quad4>::Instance();
  }
  case DRT::Element::quad8:
  {
    return AcouBoundaryImpl<DRT::Element::quad8>::Instance();
  }
  case DRT::Element::quad9:
  {
    return AcouBoundaryImpl<DRT::Element::quad9>::Instance();
  }
  case DRT::Element::tri3:
  {
    return AcouBoundaryImpl<DRT::Element::tri3>::Instance();
  }
  case DRT::Element::tri6:
  {
    return AcouBoundaryImpl<DRT::Element::tri6>::Instance();
  }
  case DRT::Element::line2:
  {
    return AcouBoundaryImpl<DRT::Element::line2>::Instance();
  }
  case DRT::Element::line3:
  {
    return AcouBoundaryImpl<DRT::Element::line3>::Instance();
  }
  case DRT::Element::nurbs2:    // 1D nurbs boundary element
  {
    return AcouBoundaryImpl<DRT::Element::nurbs2>::Instance();
  }
  case DRT::Element::nurbs3:    // 1D nurbs boundary element
  {
    return AcouBoundaryImpl<DRT::Element::nurbs3>::Instance();
  }
  case DRT::Element::nurbs4:    // 2D nurbs boundary element
  {
    return AcouBoundaryImpl<DRT::Element::nurbs4>::Instance();
  }
  case DRT::Element::nurbs9:    // 2D nurbs boundary element
  {
    return AcouBoundaryImpl<DRT::Element::nurbs9>::Instance();
  }
  default:
    dserror("Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
    break;
  }
  return NULL;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/


template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AcouBoundaryImpl<distype> * DRT::ELEMENTS::AcouBoundaryImpl<distype>::Instance(bool create)
{
  static AcouBoundaryImpl<distype> * instance;
  if (create)
  {
    if (instance==NULL)
      instance = new AcouBoundaryImpl<distype>();
  }
  else
  {
    if (instance!=NULL)
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouBoundaryImpl<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AcouBoundaryImpl<distype>::AcouBoundaryImpl()
  : xyze_(true),
    funct_(true),
    deriv_(true),
    unitnormal_(true),
    velint_(true),
    drs_(0.0),
    fac_(0.0)
{
  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouBoundaryImpl<distype>::EvaluateNeumann(
                              DRT::ELEMENTS::AcouBoundary* ele,
                              Teuchos::ParameterList&        params,
                              DRT::Discretization&           discretization,
                              DRT::Condition&                condition,
                              std::vector<int>&              lm,
                              Epetra_SerialDenseVector&      elevec1_epetra,
                              Epetra_SerialDenseMatrix*      elemat1_epetra)
{

  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouBoundaryImpl<distype>::Absorbing(
                              DRT::ELEMENTS::AcouBoundary*   ele,
                              Teuchos::ParameterList&        params,
                              DRT::Discretization&           discretization,
                              std::vector<int>&              lm,
                              Epetra_SerialDenseMatrix&      elemat1_epetra,
                              Epetra_SerialDenseMatrix&      elemat2_epetra,
                              Epetra_SerialDenseVector&  elevec1_epetra,
                              Epetra_SerialDenseVector&  elevec2_epetra,
                              Epetra_SerialDenseVector&  elevec3_epetra)
{
  /* the term representing absorbing first order boundary conditions for the
   * here given problem looks like < lambda, mu > over Gamma_ext, hence it belongs
   * to the matrix Gmat evaluated at Gamma_ext. When condensing the local
   * unknowns we build K with G as summand. Hence, we can just add the terms
   * resulting from this boundary condition to K (and hence G)
   */

  const int* nodeids = ele->NodeIds();
  DRT::Element* parent = ele->ParentElement();
  DRT::Element** faces = parent->Faces();
  bool same = false;
  for(int i=0; i<parent->NumFace(); ++i)
  {
    const int* nodeidsfaces = faces[i]->NodeIds();
    if( faces[i]->NumNode() != ele->NumNode() ) dserror("error");

    for(int j=0; j<ele->NumNode(); ++j)
    {
      if(nodeidsfaces[j]==nodeids[j])
        same = true;
      else
      {
        same = false;
        break;
      }
    }
    if(same == true)
    {
      // i is the number we were searching for!!!!
      params.set<int>("face",i);
      break;
    }
  }
  if (same == false)
    dserror("either nodeids are sorted differently or a boundary element does not know to whom it belongs");

  ele->ParentElement()->Evaluate(params,discretization,lm,elemat1_epetra,elemat2_epetra,elevec1_epetra,elevec2_epetra,elevec3_epetra);

  return 0;
}


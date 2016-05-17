/*!
\file scatra_ele_hdg_boundary_calc.cpp

\brief Routines for ScaTraHDG boundary elements

<pre>
\level 3

\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
*/
/*--------------------------------------------------------------------------*/


#include "scatra_ele_hdg_boundary_calc.H"
#include "scatra_ele_hdg.H"
#include "../drt_lib/drt_node.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraHDGBoundaryImplInterface* DRT::ELEMENTS::ScaTraHDGBoundaryImplInterface::Impl(const DRT::Element* ele)
{
  switch (ele->Shape())
  {
  case DRT::Element::quad4:
  {
    return ScaTraHDGBoundaryImpl<DRT::Element::quad4>::Instance();
  }
  case DRT::Element::quad8:
  {
    return ScaTraHDGBoundaryImpl<DRT::Element::quad8>::Instance();
  }
  case DRT::Element::quad9:
  {
    return ScaTraHDGBoundaryImpl<DRT::Element::quad9>::Instance();
  }
  case DRT::Element::tri3:
  {
    return ScaTraHDGBoundaryImpl<DRT::Element::tri3>::Instance();
  }
  case DRT::Element::tri6:
  {
    return ScaTraHDGBoundaryImpl<DRT::Element::tri6>::Instance();
  }
  case DRT::Element::line2:
  {
    return ScaTraHDGBoundaryImpl<DRT::Element::line2>::Instance();
  }
  case DRT::Element::line3:
  {
    return ScaTraHDGBoundaryImpl<DRT::Element::line3>::Instance();
  }
  case DRT::Element::nurbs2:    // 1D nurbs boundary element
  {
    return ScaTraHDGBoundaryImpl<DRT::Element::nurbs2>::Instance();
  }
  case DRT::Element::nurbs3:    // 1D nurbs boundary element
  {
    return ScaTraHDGBoundaryImpl<DRT::Element::nurbs3>::Instance();
  }
  case DRT::Element::nurbs4:    // 2D nurbs boundary element
  {
    return ScaTraHDGBoundaryImpl<DRT::Element::nurbs4>::Instance();
  }
  case DRT::Element::nurbs9:    // 2D nurbs boundary element
  {
    return ScaTraHDGBoundaryImpl<DRT::Element::nurbs9>::Instance();
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
DRT::ELEMENTS::ScaTraHDGBoundaryImpl<distype> * DRT::ELEMENTS::ScaTraHDGBoundaryImpl<distype>::Instance(bool create)
{
  static ScaTraHDGBoundaryImpl<distype> * instance;
  if (create)
  {
    if (instance==NULL)
      instance = new ScaTraHDGBoundaryImpl<distype>();
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
void DRT::ELEMENTS::ScaTraHDGBoundaryImpl<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraHDGBoundaryImpl<distype>::ScaTraHDGBoundaryImpl()
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
int DRT::ELEMENTS::ScaTraHDGBoundaryImpl<distype>::EvaluateNeumann(
    DRT::ELEMENTS::ScaTraHDGBoundary*   ele,
    Teuchos::ParameterList&             params,
    DRT::Discretization&                discretization,
    DRT::Element::LocationArray&        la,
    Epetra_SerialDenseMatrix&           elemat1_epetra,
    Epetra_SerialDenseVector&           elevec1_epetra
    )
{

  Epetra_SerialDenseVector dummy_vec2, dummy_vec3;
  Epetra_SerialDenseMatrix dummy_mat2;

  const int* nodeids = ele->NodeIds();

  DRT::Element* parent = ele->ParentElement();
  DRT::FaceElement** faces = parent->Faces();
  bool same = false;
  for(int i=0; i<parent->NumFace(); ++i)
  {
    const int* nodeidsfaces = faces[i]->NodeIds();

    if( faces[i]->NumNode() != ele->NumNode() ) break;

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
      ele->ParentElement()->Evaluate(params,discretization,la,elemat1_epetra,dummy_mat2,elevec1_epetra,dummy_vec2,dummy_vec3);
      //break;
    }
  }
  if(same == false && ( faces[0]->NumNode() != ele->NumNode() ) )
    dserror("Neumann boundary condition implemented only for surface elements");

  return 0;
}


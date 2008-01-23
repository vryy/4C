/*!----------------------------------------------------------------------
\file drt_xelement.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

//#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils.H"
#include "drt_xelement.H"


//using namespace DRT::ELEMENTS;

/*----------------------------------------------------------------------*
 |  ctor                                                     mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::XElement::XElement(const int id, const DRT::Element::ElementType etype, const int owner) :
Element(id, etype, owner)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::XElement::XElement(const DRT::ELEMENTS::XElement& old) :
Element(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::XElement::~XElement()
{
  return;
}


/*----------------------------------------------------------------------*
 |  evaluate element dummy (public)                          mwgee 12/06|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::XElement::EvaluateShape(const ParameterList& params,
                                 const DRT::Discretization&      discretization,
                               const vector<int>&              lm,
                               Epetra_SerialDenseMatrix& elemat1,
                               Epetra_SerialDenseMatrix& elemat2,
                               Epetra_SerialDenseVector& elevec1,
                               Epetra_SerialDenseVector& elevec2,
                               Epetra_SerialDenseVector& elevec3)
{
  cout << "DRT::Element::Evaluate:\n"
       << "Base class dummy routine DRT::Element::Evaluate(...) called\n"
       << __FILE__ << ":" << __LINE__ << endl;
  return -1;
}


int DRT::ELEMENTS::XElement::NumDofPerNode(const DRT::Node& node) const 
{
    return 4;
}

#endif  // #ifdef CCADISCRET

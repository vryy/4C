/*!----------------------------------------------------------------------
\file beam3tosolidmeshtying.cpp

\brief meshtying element for meshtying between a 3D beam end a 2D surface (belonging to a 3D solid) element

\level 3

\maintainer Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238

*----------------------------------------------------------------------*/

#include "beam3tosolidmeshtying.H"
#include "beam3contact_defines.H"
#include "beam3contact_utils.H"


/*----------------------------------------------------------------------*
 | Constructor (public)                                      popp 05/16 |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
CONTACT::Beam3tosolidmeshtying<numnodessol, numnodes, numnodalvalues>::Beam3tosolidmeshtying(
    const DRT::Discretization& pdiscret,
    const DRT::Discretization& cdiscret,
    const std::map<int,int>& dofoffsetmap,
    DRT::Element* element1,
    DRT::Element* element2,
    Teuchos::ParameterList beamcontactparams):
pdiscret_(pdiscret),
cdiscret_(cdiscret),
dofoffsetmap_(dofoffsetmap),
element1_(element1),
element2_(element2),
contactflag_(false)
{
  // MASTER THESIS JK 2016
  // (...)
  // still empty constructor body

  return;
}
/*----------------------------------------------------------------------*
 | End: constructor                                                     |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | copy-constructor (public)                               popp 05/2016 |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
CONTACT::Beam3tosolidmeshtying<numnodessol, numnodes, numnodalvalues>::Beam3tosolidmeshtying(const Beam3tosolidmeshtying& old):
pdiscret_(old.pdiscret_),
cdiscret_(old.cdiscret_),
dofoffsetmap_(old.dofoffsetmap_),
element1_(old.element1_),
element2_(old.element2_),
contactflag_(old.contactflag_)
{
  dserror("ERROR: Copy constructor incomplete");
  return;
}
/*----------------------------------------------------------------------*
 | End: copy-constructor                                                |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Evaluate the element (public)                           popp 05/2016 |
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes, const int numnodalvalues>
bool CONTACT::Beam3tosolidmeshtying<numnodessol, numnodes, numnodalvalues>::Evaluate(
     LINALG::SparseMatrix& stiffmatrix,
     Epetra_Vector& fint,
     const double& pp)
{
  // MASTER THESIS JK 2016
  // (...)
  // still to be implemented

  return true;
}
/*----------------------------------------------------------------------*
 | End: Evaluate the element                                            |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Template specification (public)                         popp 05/2016 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::Beam3tosolidmeshtyinginterface> CONTACT::Beam3tosolidmeshtyinginterface::Impl(
    const int numnodessol,
    const int numnodes,
    const int numnodalvalues,
    const DRT::Discretization& pdiscret,
    const DRT::Discretization& cdiscret,
    const std::map<int,int>& dofoffsetmap,
    DRT::Element* element1,
    DRT::Element* element2,
    Teuchos::ParameterList beamcontactparams)
{

  if (numnodalvalues!=1 and numnodalvalues!=2)
    dserror("Only the values 1 and 2 are valid for numnodalvalues!");

  if (numnodes!=2 and numnodes!=3 and numnodes!=4 and numnodes!=5)
    dserror("Only the values 2, 3, 4 and 5 are valid for numnodes!");

  if (numnodessol!=3 and numnodessol!=6 and numnodessol!=4 and numnodessol!=8 and numnodessol!=9)
    dserror("Only the values 3, 4, 6, 8 and 9 are valid for numnodessol!");


  switch (numnodessol)
  {
    case 4:
    {
      switch (numnodalvalues)
      {
        case 2:
        {
          return Teuchos::rcp (new CONTACT::Beam3tosolidmeshtying<4,2,2>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
          break;
        }
        default:
        {
          dserror("ERROR: Case numnodalvalues != 2 not yet implemented for Beam-to-solid meshtying)");
          break;
        }
      }
      break;
    }
    default:
    {
      dserror("ERROR: Case numnodessol != 4 not yet implemented for Beam-to-solid meshtying)");
      break;
    }
  }

  return Teuchos::null;
}
/*----------------------------------------------------------------------*
 | End: Template specification                                          |
 *----------------------------------------------------------------------*/

// Possible template cases: this is necessary for the compiler
template class CONTACT::Beam3tosolidmeshtying<3,2,2>;     // Hermite beam element, quad4 surface element
template class CONTACT::Beam3tosolidmeshtying<4,2,2>;     // Hermite beam element, quad4 surface element
//template class CONTACT::Beam3tosolidmeshtying<6,2,2>;     // Hermite beam element, tri6 surface element
//template class CONTACT::Beam3tosolidmeshtying<8,2,2>;     // Hermite beam element, quad8 surface element
//template class CONTACT::Beam3tosolidmeshtying<9,2,2>;     // Hermite beam element, quad9 surface element

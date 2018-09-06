/*!----------------------------------------------------------------------
\file beam3tosolidmeshtying.cpp

\brief meshtying element for meshtying between a 3D beam end a 2D surface (belonging to a 3D solid)
element

\level 3

\maintainer Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238

*----------------------------------------------------------------------*/

#include "beam3tosolidmeshtying.H"
#include "../drt_beaminteraction/beam3contact_utils.H"

#include "../drt_beam3/beam3.H"
#include "../drt_beam3/beam3r.H"
#include "../drt_beam3/beam3k.H"
#include "../drt_beam3/beam3eb.H"
#include "../drt_beaminteraction/beam3contact_defines.H"


/*----------------------------------------------------------------------*
 | Constructor (public)                                      popp 05/16 |
 *----------------------------------------------------------------------*/
template <const int numnodessol, const int numnodes, const int numnodalvalues>
CONTACT::Beam3tosolidmeshtying<numnodessol, numnodes, numnodalvalues>::Beam3tosolidmeshtying(
    const DRT::Discretization& pdiscret, const DRT::Discretization& cdiscret,
    const std::map<int, int>& dofoffsetmap, DRT::Element* element1,
    std::vector<DRT::Element*> element2, Teuchos::ParameterList beamcontactparams)
    : pdiscret_(pdiscret),
      cdiscret_(cdiscret),
      dofoffsetmap_(dofoffsetmap),
      element1_(element1),
      element2_(element2),
      contactflag_(false)
{
  // (1) initialize reference nodal positions (and tangents) for beam element
  for (int i = 0; i < 3 * numnodes * numnodalvalues; i++) ele1posref_(i) = 0.0;

  // (2) set reference nodal positions (and tangents) for beam element
  for (int n = 0; n < numnodes; ++n)
  {
    const DRT::Node* node = element1_->Nodes()[n];
    for (int d = 0; d < 3; ++d) ele1posref_(3 * numnodalvalues * n + d) = node->X()[d];

    // tangents
    if (numnodalvalues == 2)
    {
      LINALG::Matrix<3, 1> tan;
      const DRT::ElementType& eot = element1_->ElementType();
      if (eot == DRT::ELEMENTS::Beam3Type::Instance())
      {
        dserror("ERROR: Beam3tosolidmeshtying: numnodalvalues=2 detected for beam3 element");
      }
      else if (eot == DRT::ELEMENTS::Beam3rType::Instance())
      {
        const DRT::ELEMENTS::Beam3r* ele = dynamic_cast<const DRT::ELEMENTS::Beam3r*>(element1_);
        if (ele->HermiteCenterlineInterpolation())
          tan = ele->Tref()[n];
        else
          dserror(
              "ERROR: Beam3tosolidmeshtying: numnodalvalues=2 detected for beam3r element w/o "
              "Hermite CL");
      }
      else if (eot == DRT::ELEMENTS::Beam3kType::Instance())
      {
        const DRT::ELEMENTS::Beam3k* ele = dynamic_cast<const DRT::ELEMENTS::Beam3k*>(element1_);
        tan = ele->Tref()[n];
      }
      else if (eot == DRT::ELEMENTS::Beam3ebType::Instance())
      {
        const DRT::ELEMENTS::Beam3eb* ele = dynamic_cast<const DRT::ELEMENTS::Beam3eb*>(element1_);
        tan = ele->Tref()[n];
      }
      else
      {
        dserror("ERROR: Beam3tosolidmeshtying: Invalid beam element type");
      }

      for (int d = 0; d < 3; ++d) ele1posref_(3 * numnodalvalues * n + d + 3) = tan(d, 0);
    }
  }

  // (3) initialize reference nodal positions for all solid surface elements
  ele2posref_.resize((int)element2_.size());
  for (int e = 0; e < (int)element2_.size(); e++)
    for (int i = 0; i < 3 * numnodessol; i++) (ele2posref_[e])(i) = 0.0;

  // (4) set reference nodal positions for all solid surface elements
  for (int e = 0; e < (int)element2_.size(); e++)
  {
    for (int n = 0; n < numnodessol; ++n)
    {
      const DRT::Node* node = element2_[e]->Nodes()[n];
      for (int d = 0; d < 3; ++d) (ele2posref_[e])(3 * n + d) = node->X()[d];
    }
  }

  // (5) initialize current nodal positions (and tangents) for beam element
  for (int i = 0; i < 3 * numnodes * numnodalvalues; i++) ele1pos_(i) = 0.0;

  // (6) initialize current nodal positions for all solid surface elements
  ele2pos_.resize((int)element2_.size());
  for (int e = 0; e < (int)element2_.size(); e++)
    for (int i = 0; i < 3 * numnodessol; i++) (ele2pos_[e])(i) = 0.0;

  // MASTER THESIS JK 2016
  // (...)

  return;
}
/*----------------------------------------------------------------------*
 | End: constructor                                                     |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | copy-constructor (public)                               popp 05/2016 |
 *----------------------------------------------------------------------*/
template <const int numnodessol, const int numnodes, const int numnodalvalues>
CONTACT::Beam3tosolidmeshtying<numnodessol, numnodes, numnodalvalues>::Beam3tosolidmeshtying(
    const Beam3tosolidmeshtying& old)
    : pdiscret_(old.pdiscret_),
      cdiscret_(old.cdiscret_),
      dofoffsetmap_(old.dofoffsetmap_),
      element1_(old.element1_),
      element2_(old.element2_),
      ele1pos_(old.ele1pos_),
      ele2pos_(old.ele2pos_),
      ele1posref_(old.ele1posref_),
      ele2posref_(old.ele2posref_),
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
template <const int numnodessol, const int numnodes, const int numnodalvalues>
bool CONTACT::Beam3tosolidmeshtying<numnodessol, numnodes, numnodalvalues>::Evaluate(
    LINALG::SparseMatrix& stiffmatrix, Epetra_Vector& fint, const double& pp)
{
  // MASTER THESIS JK 2016
  // (...)
  // still to be implemented
  std::cout << "-->FOUND NEW MESHTYING GROUP" << std::endl;
  std::cout << "   Slave ID:          " << element1_->Id() << std::endl;
  std::cout << "   # Master Elements: " << (int)element2_.size() << std::endl;
  for (int m = 0; m < (int)element2_.size(); ++m)
    std::cout << "   Master ID:         " << element2_[m]->Id() << std::endl;
  std::cout << std::endl;

  return true;
}
/*----------------------------------------------------------------------*
 | End: Evaluate the element                                            |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Update nodal coordinates (public)                         popp 05/16 |
 *----------------------------------------------------------------------*/
template <const int numnodessol, const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tosolidmeshtying<numnodessol, numnodes, numnodalvalues>::UpdateElePos(
    Epetra_SerialDenseMatrix& newele1pos, std::vector<Epetra_SerialDenseMatrix>& newele2pos)
{
  // Beam element positions
  for (int i = 0; i < 3 * numnodalvalues; i++)
    for (int j = 0; j < numnodes; j++) ele1pos_(3 * numnodalvalues * j + i) = newele1pos(i, j);

  // Solid element positions
  for (int e = 0; e < (int)element2_.size(); e++)  // Loop over solid elements
    for (int i = 0; i < 3; i++)                    // Loop over nodal dofs
      for (int j = 0; j < numnodessol; j++)        // Loop over nodes
        (ele2pos_[e])(3 * j + i) = (newele2pos[e])(i, j);

  return;
}
/*----------------------------------------------------------------------*
 | End: Update nodal coordinates (public)                               |
 *----------------------------------------------------------------------*/



/*----------------------------------------------------------------------*
 | Template specification (public)                         popp 05/2016 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::Beam3tosolidmeshtyinginterface> CONTACT::Beam3tosolidmeshtyinginterface::Impl(
    const int numnodessol, const int numnodes, const int numnodalvalues,
    const DRT::Discretization& pdiscret, const DRT::Discretization& cdiscret,
    const std::map<int, int>& dofoffsetmap, DRT::Element* element1,
    std::vector<DRT::Element*> element2, Teuchos::ParameterList beamcontactparams)
{
  if (numnodalvalues != 1 and numnodalvalues != 2)
    dserror("Only the values 1 and 2 are valid for numnodalvalues!");

  if (numnodes != 2 and numnodes != 3 and numnodes != 4 and numnodes != 5)
    dserror("Only the values 2, 3, 4 and 5 are valid for numnodes!");

  if (numnodessol != 3 and numnodessol != 6 and numnodessol != 4 and numnodessol != 8 and
      numnodessol != 9)
    dserror("Only the values 3, 4, 6, 8 and 9 are valid for numnodessol!");


  switch (numnodessol)
  {
    case 4:
    {
      switch (numnodalvalues)
      {
        case 2:
        {
          return Teuchos::rcp(new CONTACT::Beam3tosolidmeshtying<4, 2, 2>(
              pdiscret, cdiscret, dofoffsetmap, element1, element2, beamcontactparams));
          break;
        }
        default:
        {
          dserror(
              "ERROR: Case numnodalvalues != 2 not yet implemented for Beam-to-solid meshtying)");
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
template class CONTACT::Beam3tosolidmeshtying<3, 2,
    2>;  // Hermite beam element, quad4 surface element
template class CONTACT::Beam3tosolidmeshtying<4, 2,
    2>;  // Hermite beam element, quad4 surface element
// template class CONTACT::Beam3tosolidmeshtying<6,2,2>;     // Hermite beam element, tri6 surface
// element template class CONTACT::Beam3tosolidmeshtying<8,2,2>;     // Hermite beam element, quad8
// surface element template class CONTACT::Beam3tosolidmeshtying<9,2,2>;     // Hermite beam
// element, quad9 surface element

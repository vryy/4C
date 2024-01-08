/*-----------------------------------------------------------------------*/
/*! \file
\brief Assembly routines for contact integration

\level 1

*/
/*-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Header                                                    farah 09/13|
 *----------------------------------------------------------------------*/
#include "baci_contact_defines.H"
#include "baci_contact_element.H"
#include "baci_contact_friction_node.H"
#include "baci_contact_integrator.H"
#include "baci_contact_node.H"
#include "baci_inpar_contact.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_linalg_serialdensevector.H"
#include "baci_mortar_coupling3d_classes.H"
#include "baci_mortar_defines.H"
#include "baci_mortar_projector.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  Assemble g~ contribution (2D / 3D)                        popp 01/08|
 |  This method assembles the contribution of a 1D/2D slave and master  |
 |  overlap pair to the weighted gap of the adjacent slave nodes.       |
 *----------------------------------------------------------------------*/
bool CONTACT::Integrator::AssembleG(
    const Epetra_Comm& comm, MORTAR::MortarElement& sele, CORE::LINALG::SerialDenseVector& gseg)
{
  // get adjacent slave nodes to assemble to
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("AssembleG: Null pointer for snodes!");

  // loop over all slave nodes
  for (int slave = 0; slave < sele.NumNode(); ++slave)
  {
    CONTACT::Node* snode = dynamic_cast<CONTACT::Node*>(snodes[slave]);

    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID()) continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound()) continue;

    double val = gseg(slave);
    snode->AddgValue(val);

    /*
#ifdef DEBUG
    std::cout << "Node: " << snode->Id() << "  Owner: " << snode->Owner() << std::endl;
    std::cout << "Weighted gap: " << snode->Getg() << std::endl;
#endif // #ifdef DEBUG
    */
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  Assemble g~ contribution (2D / 3D)                        popp 02/10|
 |  PIECEWISE LINEAR LM INTERPOLATION VERSION                           |
 *----------------------------------------------------------------------*/
bool CONTACT::Integrator::AssembleG(
    const Epetra_Comm& comm, MORTAR::IntElement& sintele, CORE::LINALG::SerialDenseVector& gseg)
{
  // get adjacent slave int nodes to assemble to
  DRT::Node** snodes = sintele.Nodes();
  if (!snodes) dserror("AssembleG: Null pointer for sintnodes!");

  // loop over all slave nodes
  for (int slave = 0; slave < sintele.NumNode(); ++slave)
  {
    CONTACT::Node* snode = dynamic_cast<CONTACT::Node*>(snodes[slave]);

    // only process slave node rows that belong to this proc
    if (snode->Owner() != comm.MyPID()) continue;

    // do not process slave side boundary nodes
    // (their row entries would be zero anyway!)
    if (snode->IsOnBound()) continue;

    double val = gseg(slave);
    snode->AddgValue(val);
  }

  return true;
}

BACI_NAMESPACE_CLOSE

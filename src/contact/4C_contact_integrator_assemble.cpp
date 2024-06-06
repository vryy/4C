/*-----------------------------------------------------------------------*/
/*! \file
\brief Assembly routines for contact integration

\level 1

*/
/*-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Header                                                    farah 09/13|
 *----------------------------------------------------------------------*/
#include "4C_contact_defines.hpp"
#include "4C_contact_element.hpp"
#include "4C_contact_friction_node.hpp"
#include "4C_contact_integrator.hpp"
#include "4C_contact_node.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_mortar_coupling3d_classes.hpp"
#include "4C_mortar_defines.hpp"
#include "4C_mortar_projector.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  Assemble g~ contribution (2D / 3D)                        popp 01/08|
 |  This method assembles the contribution of a 1D/2D slave and master  |
 |  overlap pair to the weighted gap of the adjacent slave nodes.       |
 *----------------------------------------------------------------------*/
bool CONTACT::Integrator::AssembleG(
    const Epetra_Comm& comm, Mortar::Element& sele, Core::LinAlg::SerialDenseVector& gseg)
{
  // get adjacent slave nodes to assemble to
  Core::Nodes::Node** snodes = sele.Nodes();
  if (!snodes) FOUR_C_THROW("AssembleG: Null pointer for snodes!");

  // loop over all slave nodes
  for (int slave = 0; slave < sele.num_node(); ++slave)
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
#ifdef FOUR_C_ENABLE_ASSERTIONS
    std::cout << "Node: " << snode->Id() << "  Owner: " << snode->Owner() << std::endl;
    std::cout << "Weighted gap: " << snode->Getg() << std::endl;
#endif
    */
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  Assemble g~ contribution (2D / 3D)                        popp 02/10|
 |  PIECEWISE LINEAR LM INTERPOLATION VERSION                           |
 *----------------------------------------------------------------------*/
bool CONTACT::Integrator::AssembleG(
    const Epetra_Comm& comm, Mortar::IntElement& sintele, Core::LinAlg::SerialDenseVector& gseg)
{
  // get adjacent slave int nodes to assemble to
  Core::Nodes::Node** snodes = sintele.Nodes();
  if (!snodes) FOUR_C_THROW("AssembleG: Null pointer for sintnodes!");

  // loop over all slave nodes
  for (int slave = 0; slave < sintele.num_node(); ++slave)
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

FOUR_C_NAMESPACE_CLOSE

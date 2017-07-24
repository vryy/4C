/*!----------------------------------------------------------------------
\file contact_integrator_assemble.cpp

\level 1

<pre>
\maintainer Alexander Popp

</pre>

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Header                                                    farah 09/13|
 *----------------------------------------------------------------------*/
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "contact_integrator.H"
#include "contact_node.H"
#include "contact_element.H"
#include "contact_defines.H"
#include "friction_node.H"
#include "../drt_mortar/mortar_defines.H"

#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"

#include "../drt_inpar/inpar_contact.H"
#include "../drt_mortar/mortar_projector.H"
#include "../drt_mortar/mortar_coupling3d_classes.H"


/*----------------------------------------------------------------------*
 |  Assemble g~ contribution (2D / 3D)                        popp 01/08|
 |  This method assembles the contribution of a 1D/2D slave and master  |
 |  overlap pair to the weighted gap of the adjacent slave nodes.       |
 *----------------------------------------------------------------------*/
bool CONTACT::CoIntegrator::AssembleG(const Epetra_Comm& comm,
                                      MORTAR::MortarElement& sele,
                                      Epetra_SerialDenseVector& gseg)
{
  // get adjacent slave nodes to assemble to
  DRT::Node** snodes = sele.Nodes();
  if (!snodes) dserror("ERROR: AssembleG: Null pointer for snodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sele.NumNode();++slave)
  {
    CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(snodes[slave]);

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
bool CONTACT::CoIntegrator::AssembleG(const Epetra_Comm& comm,
                                      MORTAR::IntElement& sintele,
                                      Epetra_SerialDenseVector& gseg)
{
  // get adjacent slave int nodes to assemble to
  DRT::Node** snodes = sintele.Nodes();
  if (!snodes) dserror("ERROR: AssembleG: Null pointer for sintnodes!");

  // loop over all slave nodes
  for (int slave=0;slave<sintele.NumNode();++slave)
  {
    CONTACT::CoNode* snode = dynamic_cast<CONTACT::CoNode*>(snodes[slave]);

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

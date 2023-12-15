/*----------------------------------------------------------------------*/
/*! \file
\brief a class to manage an enhanced discretization including varying number
of dofs per node on a fluid discretization for xwall

\level 2


*/
/*---------------------------------------------------------------------*/


#include "baci_lib_discret_xwall.H"

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  Constructor (public)                                        bk 10/14|
 *----------------------------------------------------------------------*/
DRT::DiscretizationXWall::DiscretizationXWall(
    const std::string name, Teuchos::RCP<Epetra_Comm> comm)
    : DiscretizationFaces(name, comm)  // use base class constructor
      {};

BACI_NAMESPACE_CLOSE

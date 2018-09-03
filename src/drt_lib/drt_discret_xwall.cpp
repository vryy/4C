/*----------------------------------------------------------------------*/
/*!
\file drt_discret_xwall.cpp
\brief a class to manage an enhanced discretization including varying number of dofs per node on a
fluid discretization for xwall

<pre>
\level 2

\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>


*----------------------------------------------------------------------*/


#include "drt_discret_xwall.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                        bk 10/14|
 *----------------------------------------------------------------------*/
DRT::DiscretizationXWall::DiscretizationXWall(
    const std::string name, Teuchos::RCP<Epetra_Comm> comm)
    : DiscretizationFaces(name, comm)  // use base class constructor
      {};

/*!----------------------------------------------------------------------
\file drt_dirichletextractor.H
\brief Declaration

<pre>
\brief Declaration
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef DRT_DIRICHLETEXTRACTOR_H
#define DRT_DIRICHLETEXTRACTOR_H

#include "../linalg/linalg_mapextractor.H"

namespace DRT
{
  class Discretization;
}

namespace DRT
{

  /// MapExtractor that covers all Dirichlet dofs
  class DirichletExtractor : public LINALG::MultiMapExtractor
  {
  public:

    enum {
      cond_nondirichlet = 0,
      cond_dirichlet = 1
    };

    using LINALG::MultiMapExtractor::Setup;

    /// setup the whole thing
    void Setup(const DRT::Discretization& dis);

    /// put zero to all Dirichlet values
    void ZeroDirichlets( Teuchos::RCP<Epetra_Vector> residual ) const;

    MAP_EXTRACTOR_VECTOR_METHODS( NonDirichlet, cond_nondirichlet )
    MAP_EXTRACTOR_VECTOR_METHODS( Dirichlet   , cond_dirichlet )
  };
}

#endif

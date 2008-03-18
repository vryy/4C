/*----------------------------------------------------------------------*/
/*!
\file elch_algorithm.H

\brief Basis of all ELCH algorithms

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef ELCH_ALGORITHM_H
#define ELCH_ALGORITHM_H

#include <Epetra_Comm.h>

//#include <Teuchos_RCP.hpp>
//#include <Teuchos_ParameterList.hpp>

#include "../drt_fsi/fsi_fluid.H"
//#include "fsi_ale.H"



namespace ELCH
{

  /// ELCH algorithm base
  /*!

    Base class of ELCH algorithms. Derives from fluid BaseAlgorithm.
    There can (and will) be different subclasses that implement 
    different coupling schemes.

    \author gjb
    \date 03/08
   */
  class Algorithm : public FSI::FluidBaseAlgorithm
  {
  public:

    /// create using a Epetra_Comm
    explicit Algorithm(Epetra_Comm& comm);

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// communicator
    const Epetra_Comm& Comm() const { return comm_; }

  protected:

  private:

    /// communication (mainly for screen output)
    Epetra_Comm& comm_;

  };

} // namespace ELCH

#endif // ELCH_ALGORITHM_H
#endif // CCADISCRET

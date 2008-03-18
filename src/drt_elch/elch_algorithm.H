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

    /// outer level ELCH time loop
    void TimeLoop();

    /// communicator
    const Epetra_Comm& Comm() const { return comm_; }

  protected:

    /// time step size
    double Dt() const { return dt_; }

    /// step number
    int Step() const { return step_; }

    /// start a new time step
    void PrepareTimeStep();

    /// take current results for converged and save for next time step
    void Update();

    /// write output
    void Output();

  private:

    /// communication (mainly for screen output)
    Epetra_Comm& comm_;

    //! @name Time stepping variables
    int step_;
    int nstep_;
    double time_;
    double maxtime_;
    double dt_;
    //@}

  };

} // namespace ELCH

#endif // ELCH_ALGORITHM_H
#endif // CCADISCRET

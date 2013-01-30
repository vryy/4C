/*----------------------------------------------------------------------*/
/*!
\file cavitation_algorithm.H

\brief Algorithm to control cavitation simulations

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-152537
</pre>
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              ghamm 11/12 |
 *----------------------------------------------------------------------*/
#ifndef CAVITATION_ALGORITHM_H
#define CAVITATION_ALGORITHM_H

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 11/12 |
 *----------------------------------------------------------------------*/
#include "particle_algorithm.H"

#include "Epetra_MpiComm.h"
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>

/*----------------------------------------------------------------------*
 | forward declarations                                     ghamm 11/12 |
 *----------------------------------------------------------------------*/
namespace ADAPTER
{
  class Fluid;
}

/*----------------------------------------------------------------------*
 | cavitation algorithm                                     ghamm 11/12 |
 *----------------------------------------------------------------------*/
/// CAVITATION: Cavitation problem
namespace CAVITATION
{

  class Algorithm : public PARTICLE::Algorithm
  {
  public:

    Algorithm(
      const Epetra_Comm& lcomm,  ///< local epetra communicator
      const Teuchos::ParameterList& params   ///< problem parameters
      );

    virtual ~Algorithm(){};

    /// do initialization of problem
    virtual void Init();

    /// read restart
    virtual void ReadRestart(
      int restart  ///< time step from which to restart
      );

    /// setup of the system
    virtual void SetupSystem();

    /// outer level TFSI time loop
    virtual void Timeloop();

    /// test results (if necessary)
    virtual void TestResults(
      const Epetra_Comm& comm  ///< local epetra communicator
      );


  protected:

    /// prepare time step
    virtual void PrepareTimeStep();

    /// solve the current particle time step
    virtual void Integrate();

    /// force on particles is calculated and applied
    virtual void CalculateAndApplyFluidForcesToParticles();

    /// update the current time step
    virtual void Update();

    /// output particle time step
    virtual void Output();

    /// find XAABB and divide into bins using the cutoff radius
    virtual void CreateBins();

    /// bins are distributed to the processors
    virtual Teuchos::RCP<Epetra_Map> DistributeBinsToProcs();

    /// setup ghosting of bins and particles
    virtual void SetupGhosting(
      Teuchos::RCP<Epetra_Map> binrowmap  ///< rowmap of bins
      );


  private:
    //! @name fluid related stuff
    //@{

    Teuchos::RCP<DRT::Discretization> fluiddis_;     ///< underlying fluid discretization

    Teuchos::RCP<ADAPTER::Fluid> fluid_;     ///< underlying fluid time integrator

    //! Info: Currently it is assumed that bins and fluid elements are fixed in space (initial content is stored here)
    //! Future: Dynamical registering of fluid elements in bins when they are moving
    std::map<int, std::set<int> > extendedfluidghosting_;   ///< map that contains all fluid element ids of each bin

    //@}

  }; // Algorithm

} // namespace CAVITATION


/*----------------------------------------------------------------------*/
#endif  // CAVITATION_ALGORITHM_H

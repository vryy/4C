/*----------------------------------------------------------------------*/
/*!
\file fsi_algorithm.H

\brief Basis of all FSI algorithms

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef FSI_ALGORITHM_H
#define FSI_ALGORITHM_H

#include <Epetra_Comm.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "fsi_structure.H"
#include "fsi_fluid.H"
#include "fsi_ale.H"

//#include "../drt_mfsi/mfsi_structure.H"


namespace FSI
{

  /// FSI algorithm base
  /*!

    Base class of FSI algorithms with ALE field. There can (and will) be
    different subclasses that implement different coupling schemes.

    \note Right now we have partitioned FSI implemented using NOX::Epetra and
    monolithic FSI implemented using NOX::Thyra. These are closely linked and
    will be even more so in the future.

    \author u.kue
    \date 02/08
   */
  class Algorithm
  {
  public:

    /// create using a Epetra_Comm
    explicit Algorithm(Epetra_Comm& comm);

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// read restart data
    void ReadRestart(int step);

    /// communicator
    const Epetra_Comm& Comm() const { return comm_; }

    /// structural field solver
    Structure& StructureField() { return *structure_; }

    /// fluid field solver
    Fluid& FluidField() { return *fluid_; }

    /// ale field solver
    Ale& AleField() { return *ale_; }

  protected:

    /// time step size
    double Dt() const { return dt_; }

    /// step number
    int Step() const { return step_; }

    //! @name Time loop building blocks

    /// tests if there are more time steps to do
    bool NotFinished() { return step_ < nstep_ and time_ <= maxtime_; }

    /// start a new time step
    void PrepareTimeStep();

    /// take current results for converged and save for next time step
    void Update();

    /// write output
    void Output();

    /// set method name for screen output
    void SetMethod(std::string method) { method_ = method; }

    //@}

  private:

    /// setup structure algorithm
    void SetupStructure();

    /// setup fluid algorithm
    void SetupFluid();

    /// setup ale algorithm
    void SetupAle();

    /// structural field solver
    Teuchos::RCP<Structure> structure_;

    /// fluid field solver
    Teuchos::RCP<Fluid>     fluid_;

    /// ale field solver
    Teuchos::RCP<Ale>       ale_;

    /// comunication (mainly for screen output)
    Epetra_Comm& comm_;

    /// method name prepared for output
    std::string method_;

    //! @name Time stepping variables
    int step_;
    int nstep_;
    double time_;
    double maxtime_;
    double dt_;
    //@}
  };
}

#endif
#endif

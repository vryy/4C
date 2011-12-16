/*----------------------------------------------------------------------*/
/*!
\file tsi_algorithm.H

\brief Basis of all TSI algorithms that perform a coupling between the
       structural field equation and temperature field equations

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 12/09 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef TSI_ALGORITHM_H
#define TSI_ALGORITHM_H

/*----------------------------------------------------------------------*
 | headers                                                   dano 12/09 |
 *----------------------------------------------------------------------*/
//#include <Epetra_Comm.h>
//
//#include <Teuchos_RCP.hpp>
//#include <Teuchos_ParameterList.hpp>

#include "../drt_adapter/adapter_structure.H"
#include "../drt_adapter/adapter_thermo.H"
#include "../drt_adapter/adapter_algorithmbase.H"


/*----------------------------------------------------------------------*
 |                                                           dano 12/09 |
 *----------------------------------------------------------------------*/
/// TSI: Thermo-Structure Interaction
namespace TSI
{
  /// TSI algorithm base
  //!
  //!
  //!  Base class of TSI algorithms. Derives from StructureBaseAlgorithm and
  //!  ThermoBaseAlgorithm with temperature field.
  //!  There can (and will) be different subclasses that implement different
  //!  coupling schemes.
  //!
  //!  \warning The order of calling the two BaseAlgorithm-constructors (that
  //!  is the order in which we list the base classes) is important here! In the
  //!  constructors control file entries are written. And these entries define
  //!  the order in which the filters handle the Discretizations, which in turn
  //!  defines the dof number ordering of the Discretizations... Don't get
  //!  confused. Just always list structure, thermo. In that order.
  //!
  //!  \author u.kue
  //!  \date 02/08
  class Algorithm : public ADAPTER::AlgorithmBase,
                    public ADAPTER::StructureBaseAlgorithm,
                    public ADAPTER::ThermoBaseAlgorithm
  {
  public:

    /// create using a Epetra_Comm
    explicit Algorithm(const Epetra_Comm& comm);

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// outer level time loop (to be implemented by deriving classes)
    void TimeLoop();

    /// time loop for TSI algorithm with one-way coupling
    void TimeLoopSequStagg();

    /// time loop for TSI algorithm with one-way coupling
    void TimeLoopOneWay();

    /// time loop for TSI algorithm with iteration between fields (full coupling)
    void TimeLoopFull();

    /// read restart data
    virtual void ReadRestart(
      int step  ///< step number where the calculation is continued
      );

  protected:

    //! @name Time loop building blocks

    /// calculate stresses, strains, energies
    void PrepareOutput()
    {
      StructureField().PrepareOutput();
    }

    /// take current results for converged and save for next time step
    virtual void Update();

    /// write output
    virtual void Output();

    /// communicate displacement vector to thermal field to enable their
    /// visualisation on the deformed body
    void OutputDeformationInThr(
      Teuchos::RCP<const Epetra_Vector> dispnp,
      Teuchos::RCP<DRT::Discretization> structdis
      );

    /// displacements at time n+1 for thermal output
    Teuchos::RCP<Epetra_MultiVector> dispnp_;

    //@}

    //! @name Solve

    /// solve temperature equations for current time step
    void DoThermoStep();

    /// solve displacement equations for current time step
    void DoStructureStep();

    // two-way coupling (iterative staggered)

    /// outer iteration loop
    void OuterIterationLoop();

    //@}

    //! @name Access methods

    /// velocity calculation given the displacements (like in FSI)
    Teuchos::RCP<Epetra_Vector> CalcVelocity(
      Teuchos::RCP<const Epetra_Vector> dispnp
      ) const;

    //@}

    /// convergence check for iterative staggered TSI solver
    bool ConvergenceCheck(int itnum, int itmax, double ittol);

    /// maximum iteration steps
    int itmax_;
    /// convergence tolerance
    double ittol_;

    /// temperature increment of the outer loop
    Teuchos::RCP<Epetra_Vector> tempincnp_;
    /// displacement increment of the outer loop
    Teuchos::RCP<Epetra_Vector> dispincnp_;

  private:

    /// displacements from time step begin (t_n)
    Teuchos::RCP<Epetra_Vector> dispn_;
    /// velocities from time step begin (t_n)
    Teuchos::RCP<Epetra_Vector> veln_;
    /// velocities from time step begin (t_n+1)
    Teuchos::RCP<Epetra_Vector> velnp_;

    /// temperature from time step begin (t_n)
    Teuchos::RCP<Epetra_Vector> tempn_;

    /// we couple based on displacements
    bool displacementcoupling_;
    /// quasi-static solution of the mechancial equation
    bool quasistatic_;

  };  // Algorithm
}  // namespace TSI

/*----------------------------------------------------------------------*/
#endif  // TSI_ALGORITHM_H
#endif  // CCADISCRET

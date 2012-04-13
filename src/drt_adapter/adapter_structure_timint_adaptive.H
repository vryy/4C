/*----------------------------------------------------------------------*/
/*!
\file adapter_structure_timint_adaptive.H

\brief Structural time integration field adapter

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef ADAPTER_STRUCTURE_TIMINT_ADAPT_H
#define ADAPTER_STRUCTURE_TIMINT_ADAPT_H

/*----------------------------------------------------------------------*/
/* headers */
#include "FSIStructureWrapper.H"


// forward declarations:
namespace STR
{
  class TimInt;
  class TimAda;
}


/*----------------------------------------------------------------------*/
/* adapting adapter */
namespace ADAPTER {

  /*====================================================================*/
  /*!
   * \brief Adapter to general implicit structural time integration
   *
   * \date 08/08
   */
  class StructureTimIntAda : public FSIStructureWrapper
  {
  public:

    /// Constructor
    StructureTimIntAda(
      Teuchos::RCP<STR::TimAda> sta,
      Teuchos::RCP<STR::TimInt> sti
    );

    /// recalculate structural matrices for tsi simulations
    virtual void TSIMatrix() { dserror("not implemented"); return; }

    // access to contact/meshtying manager
    virtual Teuchos::RCP<MORTAR::ManagerBase> ContactManager() {dserror("not implemented"); return Teuchos::null; }

    // access to locsys manager
    virtual Teuchos::RCP<DRT::UTILS::LocsysManager> LocsysManager() {dserror("not implemented"); return Teuchos::null; }

    /// @name Time step helpers
    //@{

    /// Integrate from t1 to t2
    virtual void Integrate();

    /// set time step size
    virtual void SetTimeStepSize(double timestepsize){ dserror("not implemented"); return; }

    /// poroelasticity specific method
    virtual void PoroInitForceStiffResidual() { dserror("not implemented"); return; }

    /// start new time step
    virtual void PrepareTimeStep() { dserror("not implemented"); return; }

    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc) { dserror("not implemented"); return; }

    /// update at time step end
    virtual void Update() { dserror("not implemented"); return; }

    /// prepare output (i.e. calculate stresses, strains, energies)
    virtual void PrepareOutput();

    /// output results
    virtual void Output();

    //@}


    /// evaluate reference state (for frictional contact)
    virtual void EvaluateReferenceState() { dserror("not implemented"); return; }

    /// get constraint manager defined in the structure
    virtual const Teuchos::RCP<UTILS::ConstrManager> GetConstraintManager() { dserror("not implemented"); return Teuchos::null; }


    /// @name Solver calls
    //@{

    /// non-linear solve
    ///
    /// Do the nonlinear solve for the time step. All boundary conditions have
    /// been set.
    virtual void Solve() { dserror("not implemented"); return; }


    //@}

    //! @name TSI specific methods
    //@{

    //! Apply temperatures  for TSI

    /// apply temperatures to structural solver
    ///
    /// This prepares a new solve of the structural field within one time
    /// step. The middle values are newly created.
    ///
    /// \note This is not yet the most efficient implementation.
    virtual void ApplyTemperatures(Teuchos::RCP<const Epetra_Vector> temp) { dserror("not implemented"); return; }

    // poroelasticity specific method
    virtual void ApplyVelAndPress(Teuchos::RCP<const Epetra_Vector> velnp) { dserror("not implemented"); return; }

    //@}

    //! @name Structure with ale specific methods
    //@{

    /// material displacements (structure with ale)
    virtual Teuchos::RCP<Epetra_Vector> DispMat() {dserror("not implemented"); return Teuchos::null; }

    //@}


  private:

    /// the actual structure algorithm
    Teuchos::RCP<STR::TimAda> structure_;

  };  // class StructureTimIntAda

}  // namespace ADAPTER

/*----------------------------------------------------------------------*/
#endif

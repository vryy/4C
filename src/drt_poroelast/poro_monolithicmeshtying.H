/*----------------------------------------------------------------------*/
/*!
 \file poro_monolithicmeshtying.H

 \brief base for porous media monolithic meshtying method

// Masterthesis of h.Willmann under supervision of Anh-Tu Vuong and Christoph Ager
// Originates from poro_monolithic

\level 2

\maintainer Ager Christoph
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
 *----------------------------------------------------------------------*/


#ifndef PORO_MONOLITHICMESHTYING_H_
#define PORO_MONOLITHICMESHTYING_H_

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/

#include "poroelast_monolithic.H"

namespace ADAPTER
{
  class CouplingPoroMortar;
}

namespace POROELAST
{
  //
  class MonolithicMeshtying: public Monolithic
  {
  public:

    //! create using a Epetra_Comm
    explicit MonolithicMeshtying(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    //! virtual destructor to support polymorph destruction
    virtual ~MonolithicMeshtying(){};

    //! Setup the monolithic system
    virtual void SetupSystem();

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> iterinc, //!< increment between iteration i and i+1
                  bool firstiter = false);

    //use monolithic update and set old meshtying quantities at the end of a timestep
    virtual void Update();

    //Recover Lagrange Multiplier after Newton step
    virtual void RecoverLagrangeMultiplierAfterNewtonStep(Teuchos::RCP<const Epetra_Vector> x);

    //build meshtying specific norms where meshtying constraint residuals are evaluated separately
    virtual void BuildConvergenceNorms();

    //extractor to split fluid RHS vector for convergence check
    //should be named fluidvelocityactiverowdofmap_
    Teuchos::RCP<const LINALG::MultiMapExtractor> FluidVelActiveDofExtractor() const
    { return fvelactiverowdofmap_;}

    //setup meshtying activedof extractors
    void SetupExtractor();

    //decide convergence with additional evaluation of meshtying constraint residuals
    virtual bool Converged();

    //setup solver with additional residual tolerances for meshtying
    virtual bool SetupSolver();

    //contains header to PrintNewtonIter with meshtying solver tolerance
    virtual void PrintNewtonIterHeaderStream(std::ostringstream& oss);

    //contains text to PrintNewtonIter with meshtying residuals
    virtual void PrintNewtonIterTextStream(std::ostringstream& oss);

  private:
    //nonlinear mortar adapter used to evaluate meshtying
    Teuchos::RCP<ADAPTER::CouplingPoroMortar> mortar_adapter_;

    // fluid velocity dof row map splitted in active dofs and the rest (no pressures)
    Teuchos::RCP<LINALG::MultiMapExtractor> fvelactiverowdofmap_; //should be named fluidvelocityactiverowdofmap_, but kept shorter

    double normrhsfactiven_;     //!< norm of coupling part of residual forces (fluid )

    double tolfres_ncoup_;      //residuum tolerance for porofluid normal coupling condition

  };
}

#endif /* PORO_MONOLITHICMESHTYING_H_ */

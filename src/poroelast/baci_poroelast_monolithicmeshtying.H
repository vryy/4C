/*----------------------------------------------------------------------*/
/*! \file

 \brief base for porous media monolithic meshtying method

// Masterthesis of h.Willmann under supervision of Anh-Tu Vuong and Johannes Kremheller
// Originates from poro_monolithic

\level 2

 *----------------------------------------------------------------------*/


#ifndef BACI_POROELAST_MONOLITHICMESHTYING_H
#define BACI_POROELAST_MONOLITHICMESHTYING_H


#include "baci_config.H"

#include "baci_poroelast_monolithic.H"

BACI_NAMESPACE_OPEN

namespace ADAPTER
{
  class CouplingPoroMortar;
}

namespace POROELAST
{
  class MonolithicMeshtying : public Monolithic
  {
   public:
    //! create using a Epetra_Comm
    explicit MonolithicMeshtying(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams,
        Teuchos::RCP<CORE::LINALG::MapExtractor> porosity_splitter);

    //! Setup the monolithic system
    void SetupSystem() override;

    //! evaluate all fields at x^n+1_i+1 with x^n+1_i+1 = x_n+1_i + iterinc
    void Evaluate(
        Teuchos::RCP<const Epetra_Vector> iterinc,  //!< increment between iteration i and i+1
        bool firstiter = false) override;

    //! use monolithic update and set old meshtying quantities at the end of a timestep
    void Update() override;

    //! Recover Lagrange Multiplier after Newton step
    void RecoverLagrangeMultiplierAfterNewtonStep(
        Teuchos::RCP<const Epetra_Vector> iterinc) override;

    //! build meshtying specific norms where meshtying constraint residuals are evaluated separately
    void BuildConvergenceNorms() override;

    //! extractor to split fluid RHS vector for convergence check
    //! should be named fluidvelocityactiverowdofmap_
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> FluidVelActiveDofExtractor() const
    {
      return fvelactiverowdofmap_;
    }

    //! setup meshtying activedof extractors
    void SetupExtractor();

    //! decide convergence with additional evaluation of meshtying constraint residuals
    bool Converged() override;

    //! setup solver with additional residual tolerances for meshtying
    bool SetupSolver() override;

    //! contains header to PrintNewtonIter with meshtying solver tolerance
    void PrintNewtonIterHeaderStream(std::ostringstream& oss) override;

    //! contains text to PrintNewtonIter with meshtying residuals
    void PrintNewtonIterTextStream(std::ostringstream& oss) override;

   private:
    //! nonlinear mortar adapter used to evaluate meshtying
    Teuchos::RCP<ADAPTER::CouplingPoroMortar> mortar_adapter_;

    //! fluid velocity dof row map splitted in active dofs and the rest (no pressures)
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor>
        fvelactiverowdofmap_;  //!< should be named fluidvelocityactiverowdofmap_, but kept shorter

    double normrhsfactiven_;  //!< norm of coupling part of residual forces (fluid )

    double tolfres_ncoup_;  //!< residuum tolerance for porofluid normal coupling condition
  };
}  // namespace POROELAST

BACI_NAMESPACE_CLOSE

#endif  // POROELAST_MONOLITHICMESHTYING_H

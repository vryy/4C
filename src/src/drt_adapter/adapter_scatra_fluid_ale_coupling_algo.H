/*----------------------------------------------------------------------*/
/*!
\file adapter_scatra_fluid_ale_coupling_algo.H

\brief Basis of all algorithms that perform a coupling between Navier-Stokes
       and (active or passive) scalar transport equations

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef ADAPTER_SCATRA_FLUID_ALE_COUPLING_ALGO_H
#define ADAPTER_SCATRA_FLUID_ALE_COUPLING_ALGO_H

#include "../drt_ale/ale.H"
#include "../drt_adapter/adapter_scatra_fluid_coupling_algorithm.H"
#include "adapter_coupling.H"

namespace ADAPTER
{
  /// basis coupling algorithm for scalar transport with Navier-Stokes on moving meshes
  /*!

    Base class for scalar transport problems coupled to Navier-Stokes velocity field on
    deforming meshes.
    Derives from ScaTraFluidCouplingAlgorithm and AleBaseAlgorithm and establishes
    the fluid-ale coupling.
    Different application coupling algorithms will inherit from this base class
    (at the moment only electrochemistry applications).

    \author gjb
    \date 05/09
   */
  class ScaTraFluidAleCouplingAlgorithm :
        public ADAPTER::ScaTraFluidCouplingAlgorithm,
        public ALE::AleBaseAlgorithm
  {
  public:

    /// constructor using a Epetra_Comm
    explicit ScaTraFluidAleCouplingAlgorithm(
        Epetra_Comm& comm,              ///< communicator
        const Teuchos::ParameterList& prbdyn, ///< problem-specific parameters
        const std::string condname ///< name of condition that defines fluid-ale coupling
        );

    /// virtual destructor to support polymorph destruction
    virtual ~ScaTraFluidAleCouplingAlgorithm();

    /// read restart data
    virtual void ReadRestart(
        int step ///< step number where the calculation is continued
        );

    /// solve fluid-ale
    virtual void FluidAleNonlinearSolve(
        Teuchos::RCP<Epetra_Vector> idisp,
        Teuchos::RCP<Epetra_Vector> ivel
        );

  protected:

    //! @name Transfer helpers

    /// field transform
    virtual Teuchos::RCP<Epetra_Vector> AleToFluidField(Teuchos::RCP<Epetra_Vector> iv) const;

    /// field transform
    virtual Teuchos::RCP<Epetra_Vector> AleToFluidField(Teuchos::RCP<const Epetra_Vector> iv) const;

    /// interface transform
    virtual Teuchos::RCP<Epetra_Vector> FluidToAle(Teuchos::RCP<Epetra_Vector> iv) const;

    /// interface transform
    virtual Teuchos::RCP<Epetra_Vector> FluidToAle(Teuchos::RCP<const Epetra_Vector> iv) const;

  private:

    /// coupling of fluid and ale (whole field)
    Coupling               coupfa_;

    /// coupling of fluid and ale (interface only)
    Coupling               icoupfa_;

    /// coupling of fluid and ale at the free surface
    Coupling               fscoupfa_;

  };

} // namespace ADAPTER

#endif // ADAPTER_SCATRA_FLUID_ALE_COUPLING_ALGO_H
#endif // CCADISCRET

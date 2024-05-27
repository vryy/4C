/*----------------------------------------------------------------------*/
/*! \file

\brief Basis of all algorithms that perform a coupling between Navier-Stokes
       and (active or passive) scalar transport equations
\level 2


*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_ADAPTER_SCATRA_FLUID_ALE_COUPLING_ALGO_HPP
#define FOUR_C_ADAPTER_SCATRA_FLUID_ALE_COUPLING_ALGO_HPP

#include "4C_config.hpp"

#include "4C_adapter_ale_fluid.hpp"
#include "4C_adapter_scatra_fluid_coupling_algorithm.hpp"
#include "4C_coupling_adapter.hpp"

FOUR_C_NAMESPACE_OPEN


namespace ADAPTER
{
  // forward declarations

  class Coupling;


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
  class ScaTraFluidAleCouplingAlgorithm : public ADAPTER::ScaTraFluidCouplingAlgorithm,
                                          public ADAPTER::AleBaseAlgorithm
  {
   public:
    /// constructor using a Epetra_Comm
    ScaTraFluidAleCouplingAlgorithm(const Epetra_Comm& comm,  ///< communicator
        const Teuchos::ParameterList& prbdyn,                 ///< problem-specific parameters
        const std::string condname,  ///< name of condition that defines fluid-ale coupling
        const Teuchos::ParameterList& solverparams);


    /// setup
    void Setup() override;

    /// init
    void Init() override;

    /// read restart data (pure virtual)
    void read_restart(int step  ///< step number where the calculation is continued
        ) override = 0;

    /// solve fluid-ale
    virtual void fluid_ale_nonlinear_solve(Teuchos::RCP<Epetra_Vector> idisp,
        Teuchos::RCP<Epetra_Vector> ivel, const bool pseudotransient);

    /// access to ale field
    const Teuchos::RCP<ADAPTER::AleFluidWrapper>& ale_field() { return ale_; }

   protected:
    //! @name Transfer helpers

    /// field transform
    virtual Teuchos::RCP<Epetra_Vector> ale_to_fluid_field(Teuchos::RCP<Epetra_Vector> iv) const;

    /// field transform
    virtual Teuchos::RCP<Epetra_Vector> ale_to_fluid_field(
        Teuchos::RCP<const Epetra_Vector> iv) const;

    /// interface transform
    virtual Teuchos::RCP<Epetra_Vector> fluid_to_ale(Teuchos::RCP<Epetra_Vector> iv) const;

    /// interface transform
    virtual Teuchos::RCP<Epetra_Vector> fluid_to_ale(Teuchos::RCP<const Epetra_Vector> iv) const;

   private:
    /// ALE-fluid wrapper
    Teuchos::RCP<AleFluidWrapper> ale_;

    /// coupling of fluid and ale (whole field)
    Teuchos::RCP<CORE::ADAPTER::Coupling> coupfa_;

    /// coupling of fluid and ale (interface only)
    Teuchos::RCP<CORE::ADAPTER::Coupling> icoupfa_;

    /// coupling of fluid and ale at the free surface
    Teuchos::RCP<CORE::ADAPTER::Coupling> fscoupfa_;

    /// condition name
    const std::string condname_;
  };

}  // namespace ADAPTER

FOUR_C_NAMESPACE_CLOSE

#endif

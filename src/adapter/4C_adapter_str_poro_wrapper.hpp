/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for structure or poro time integration

\level 2


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_ADAPTER_STR_PORO_WRAPPER_HPP
#define FOUR_C_ADAPTER_STR_PORO_WRAPPER_HPP

#include "4C_config.hpp"

#include "4C_adapter_field_wrapper.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_structure_aux.hpp"

FOUR_C_NAMESPACE_OPEN

namespace GLOBAL
{
  class Problem;
}

namespace POROELAST
{
  class Monolithic;
}

namespace ADAPTER
{
  class FluidPoro;

  /// Just wrap, do nothing new, provides methods which are not available for Base Class
  /// ADAPTER::Field!
  class StructurePoroWrapper : public FieldWrapper
  {
   public:
    /// constructor
    explicit StructurePoroWrapper(
        Teuchos::RCP<Field> field, FieldWrapper::Fieldtype type, bool NOXCorrection = false);

    /// setup
    void Setup();

    /// communication object at the interface
    virtual Teuchos::RCP<const STR::MapExtractor> Interface() const
    {
      return structure_->Interface();
    }

    /// direct access to discretization
    virtual Teuchos::RCP<DRT::Discretization> Discretization()
    {
      return structure_->Discretization();
    }

    /// return time integration factor
    virtual double TimIntParam() const { return structure_->TimIntParam(); }

    /// Access to output object
    virtual Teuchos::RCP<IO::DiscretizationWriter> DiscWriter() { return structure_->DiscWriter(); }

    /// unknown displacements at \f$t_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp() const { return structure_->Dispnp(); }

    /// unknown displacements at \f$t_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessDispnp() const
    {
      return structure_->WriteAccessDispnp();
    }

    /// get constraint manager defined in the structure
    virtual Teuchos::RCP<CONSTRAINTS::ConstrManager> get_constraint_manager()
    {
      return structure_->get_constraint_manager();
    }

    // access to contact/meshtying bridge
    virtual Teuchos::RCP<CONTACT::MeshtyingContactBridge> meshtying_contact_bridge()
    {
      return structure_->meshtying_contact_bridge();
    }

    /// extract interface displacements at \f$t_{n}\f$
    virtual Teuchos::RCP<Epetra_Vector> extract_interface_dispn()
    {
      return structure_->extract_interface_dispn();
    }

    /// extract interface displacements at \f$t_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> extract_interface_dispnp()
    {
      return structure_->extract_interface_dispnp();
    }

    /// extract interface displacements at \f$t_{n}\f$
    virtual Teuchos::RCP<Epetra_Vector> extract_fpsi_interface_dispn()
    {
      return structure_->Interface()->ExtractFPSICondVector(structure_->Dispn());
    }

    /// extract interface displacements at \f$t_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> extract_fpsi_interface_dispnp()
    {
      return structure_->Interface()->ExtractFPSICondVector(structure_->Dispnp());
    }

    //! unique map of all dofs that should be constrained with DBC
    virtual Teuchos::RCP<const Epetra_Map> CombinedDBCMap();

    //! perform result test
    void TestResults(GLOBAL::Problem* problem);

    //! return poro PoroField
    const Teuchos::RCP<POROELAST::Monolithic>& PoroField();

    //! return poro StructureField
    const Teuchos::RCP<FSIStructureWrapper>& StructureField();

    //! return poro FluidField
    const Teuchos::RCP<ADAPTER::FluidPoro>& FluidField();

    //! Insert FSI Condition Vector
    Teuchos::RCP<Epetra_Vector> InsertFSICondVector(Teuchos::RCP<const Epetra_Vector> cond);

    //! Recover Lagrange Multiplier during iteration (does nothing for structure)
    void recover_lagrange_multiplier_after_newton_step(Teuchos::RCP<Epetra_Vector> iterinc);

    bool isPoro() { return (type_ == FieldWrapper::type_PoroField); }

   protected:
    Teuchos::RCP<POROELAST::Monolithic> poro_;     ///< underlying poro time integration
    Teuchos::RCP<FSIStructureWrapper> structure_;  ///< underlying structural time integration
  };
}  // namespace ADAPTER

FOUR_C_NAMESPACE_CLOSE

#endif

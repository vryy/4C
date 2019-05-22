/*----------------------------------------------------------------------*/
/*!

\brief Wrapper for structure or poro time integration

\level 2

\maintainer Christoph Ager

*/
/*----------------------------------------------------------------------*/

#ifndef ADAPTER_STRUCTPORO_WRAPPER_H
#define ADAPTER_STRUCTPORO_WRAPPER_H

#include "ad_field_wrapper.H"
#include "ad_str_fsiwrapper.H"
#include "../drt_structure/stru_aux.H"

namespace DRT
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
    virtual Teuchos::RCP<const STR::AUX::MapExtractor> Interface() const
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
    virtual const Teuchos::RCP<UTILS::ConstrManager> GetConstraintManager()
    {
      return structure_->GetConstraintManager();
    }

    // access to contact/meshtying bridge
    virtual Teuchos::RCP<CONTACT::MeshtyingContactBridge> MeshtyingContactBridge()
    {
      return structure_->MeshtyingContactBridge();
    }

    /// extract interface displacements at \f$t_{n}\f$
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispn()
    {
      return structure_->ExtractInterfaceDispn();
    }

    /// extract interface displacements at \f$t_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispnp()
    {
      return structure_->ExtractInterfaceDispnp();
    }

    /// extract interface displacements at \f$t_{n}\f$
    virtual Teuchos::RCP<Epetra_Vector> ExtractFPSIInterfaceDispn()
    {
      return structure_->Interface()->ExtractFPSICondVector(structure_->Dispn());
    }

    /// extract interface displacements at \f$t_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> ExtractFPSIInterfaceDispnp()
    {
      return structure_->Interface()->ExtractFPSICondVector(structure_->Dispnp());
    }

    //! unique map of all dofs that should be constrained with DBC
    virtual Teuchos::RCP<const Epetra_Map> CombinedDBCMap();

    //! perform result test
    void TestResults(DRT::Problem* problem);

    //! return poro PoroField
    const Teuchos::RCP<POROELAST::Monolithic>& PoroField();

    //! return poro StructureField
    const Teuchos::RCP<FSIStructureWrapper>& StructureField();

    //! return poro FluidField
    const Teuchos::RCP<ADAPTER::FluidPoro>& FluidField();

    //! Insert FSI Condition Vector
    Teuchos::RCP<Epetra_Vector> InsertFSICondVector(Teuchos::RCP<const Epetra_Vector> cond);

    //! Recover Lagrange Multiplier during iteration (does nothing for structure)
    void RecoverLagrangeMultiplierAfterNewtonStep(Teuchos::RCP<Epetra_Vector> iterinc);

    bool isPoro() { return (type_ == FieldWrapper::type_PoroField); }

   protected:
    Teuchos::RCP<POROELAST::Monolithic> poro_;     ///< underlying poro time integration
    Teuchos::RCP<FSIStructureWrapper> structure_;  ///< underlying structural time integration
  };
}  // namespace ADAPTER

#endif

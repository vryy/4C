/*----------------------------------------------------------------------*/
/*!
\file ad_str_poro_wrapper.H

\brief Wrapper for structure or poro time integration

<pre>
Maintainer: Ager Christoph
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef ADAPTER_STRUCTPORO_WRAPPER_H
#define ADAPTER_STRUCTPORO_WRAPPER_H

#include "ad_field_wrapper.H"

#include "ad_str_fpsiwrapper.H"
#include "../drt_poroelast/poroelast_monolithic.H"
#include "ad_fld_poro.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_structure/stru_aux.H"


namespace ADAPTER
{
  /// Just wrap, do nothing new, provides methods which are not available for Base Class ADAPTER::Field!
  class StructurePoroWrapper : public FieldWrapper
  {
  public:

    /// constructor
    explicit StructurePoroWrapper(Teuchos::RCP<Field> field, FieldWrapper::Fieldtype type, bool NOXCorrection = false)
    :FieldWrapper(field,type,NOXCorrection)
    {
      switch (type_)
      {
        case FieldWrapper::type_StructureField:
          structure_ = Teuchos::rcp_dynamic_cast<FSIStructureWrapper>(field_);
          if (structure_ == Teuchos::null)
            dserror("StructurePoroWrapper: Cast from Field to FSIStructureWrapper failed!");
          poro_ = Teuchos::null;
        break;
        case FieldWrapper::type_PoroField:
          poro_ = Teuchos::rcp_dynamic_cast<POROELAST::Monolithic>(field_);
          if (poro_ == Teuchos::null)
                      dserror("StructurePoroWrapper: Cast from Field to PoroBase failed!");
          structure_ = poro_->StructureField();
        break;
        default:
          dserror("StructurePoroWrapper - FieldWrapper::Fieldtype not available for this wrapper!");
        break;
      }
    }

    /// communication object at the interface
    virtual Teuchos::RCP<const STR::AUX::MapExtractor> Interface() const
    { return structure_->Interface(); }

    /// direct access to discretization
    virtual Teuchos::RCP<DRT::Discretization> Discretization()
    { return structure_->Discretization(); }

    /// return time integration factor
    virtual double TimIntParam() const
    { return structure_->TimIntParam(); }

    /// Access to output object
    virtual Teuchos::RCP<IO::DiscretizationWriter> DiscWriter()
    { return structure_->DiscWriter();   }

    /// unknown displacements at \f$t_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp() const
    { return structure_->Dispnp(); }

    /// get constraint manager defined in the structure
    virtual const Teuchos::RCP< UTILS::ConstrManager> GetConstraintManager()
    { return structure_->GetConstraintManager(); }

    // access to contact/meshtying bridge
    virtual Teuchos::RCP<CONTACT::MeshtyingContactBridge> MeshtyingContactBridge()
    { return structure_->MeshtyingContactBridge(); }

    /// extract interface displacements at \f$t_{n}\f$
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispn()
    { return structure_->ExtractInterfaceDispn();}

    /// extract interface displacements at \f$t_{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceDispnp()
    { return structure_->ExtractInterfaceDispnp();}

    //! unique map of all dofs that should be constrained with DBC
    virtual Teuchos::RCP<const Epetra_Map> CombinedDBCMap()
    {
     switch(type_)
     {
     case FieldWrapper::type_StructureField:
       return structure_->GetDBCMapExtractor()->CondMap();
       break;
     case FieldWrapper::type_PoroField:
       return poro_->CombinedDBCMap();
       break;
     default:
       dserror("StructurePoroWrapper: type for this wrapper not considered!");
       return Teuchos::null;
       break;
     }
   }

//   //! perform result test
   void TestResults(DRT::Problem* problem)
   {
     problem->AddFieldTest(structure_->CreateFieldTest());

     if (type_ == FieldWrapper::type_PoroField)
       problem->AddFieldTest(poro_->FluidField()->CreateFieldTest());
   }

   const Teuchos::RCP<POROELAST::Monolithic>& PoroField()
    {
     if (type_ == FieldWrapper::type_PoroField)
       return poro_;
     else
       dserror("StructurePoroWrapper - Field not PoroField!");
     return poro_; //do not remove dserror!!! - return just to make complier happy :-)
    }

   const Teuchos::RCP<FSIStructureWrapper>& StructureField()
    {
     if (type_ == FieldWrapper::type_PoroField || type_ == FieldWrapper::type_StructureField)
       return structure_;
     else
       dserror("StructurePoroWrapper - Field not Structural- or Poro-Field!");
     return structure_; //do not remove dserror!!! - return just to make complier happy :-)
    }

   //! return poro FluidField
   const Teuchos::RCP<ADAPTER::FluidPoro>& FluidField()
    {
     if (type_ == FieldWrapper::type_PoroField)
       return poro_->FluidField();
     else
       dserror("StructurePoroWrapper - Field not PoroField (no poro fluid field!");
     return poro_->FluidField(); //do not remove dserror!!! - return just to make complier happy :-)
    }

   //!Insert FSI Condition Vector
   Teuchos::RCP<Epetra_Vector> InsertFSICondVector(Teuchos::RCP<const Epetra_Vector> cond)
    {
       Teuchos::RCP<Epetra_Vector> tmpcond;
       switch(type_)
       {
       case FieldWrapper::type_StructureField:
         return Interface()->InsertFSICondVector(cond);
         break;
       case FieldWrapper::type_PoroField:
         tmpcond = Interface()->InsertFSICondVector(cond);
         return poro_->Extractor()->InsertVector(tmpcond,0); //into structural part = 0
         break;
       default:
         dserror("StructurePoroWrapper: type for this wrapper not considered!");
         return Teuchos::null;
         break;
       }
    }

   //!Recover Lagrange Multiplier during iteration (does nothing for structure)
   void RecoverLagrangeMultiplierAfterNewtonStep(Teuchos::RCP<Epetra_Vector> iterinc)
   {
     if (type_ == FieldWrapper::type_PoroField)
       poro_->RecoverLagrangeMultiplierAfterNewtonStep(iterinc);
   }

  protected:

    Teuchos::RCP<POROELAST::Monolithic> poro_;   ///< underlying poro time integration
    Teuchos::RCP<FSIStructureWrapper> structure_;   ///< underlying structural time integration
  };
}

#endif

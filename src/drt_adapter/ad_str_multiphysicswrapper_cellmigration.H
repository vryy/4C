/*----------------------------------------------------------------------*/
/*! \file

\brief Adapter wrapper for multiphysics in cell migration.

\maintainer Jonas Eichinger

\level 1
*/
/*----------------------------------------------------------------------*/


#ifndef MULTIPHYSICS_STRUCTURE_WRAPPER_CELLMIGRATION_H
#define MULTIPHYSICS_STRUCTURE_WRAPPER_CELLMIGRATION_H


#include "ad_str_wrapper.H"
#include "../drt_immersed_problem/str_model_evaluator_multiphysics_cellmigration.H"

namespace ADAPTER
{
  class StructAleWrapper;
  class SSIStructureWrapper;
  class FSIStructureWrapperImmersed;


  /*! \brief This Multiphysics adapter wraps the subproblem specific structure wrappers.
   *
   *
   * \date 12/2016
   * \author rauch (rauch@lnm.mw.tum.de)          */
  class MultiphysicsStructureWrapperCellMigration : public StructureWrapper
  {
   public:
    /// constructor
    explicit MultiphysicsStructureWrapperCellMigration(Teuchos::RCP<Structure> structure);

    /// Destructor
    virtual ~MultiphysicsStructureWrapperCellMigration() { ; }

    //! @name Access to field wrappers
    //@{
    /// return the structure wrapper for SSI
    const Teuchos::RCP<::ADAPTER::SSIStructureWrapper>& GetSSIStructureWrapperPtr()
    {
      return ssi_structure_wrapper_;
    };
    /// return the structure wrapper for immersed FSI
    const Teuchos::RCP<::ADAPTER::FSIStructureWrapperImmersed>& GetFSIStructureWrapperPtr()
    {
      return fsi_structure_wrapper_;
    };
    /// return the structure wrapper for immersed FSI
    const Teuchos::RCP<::ADAPTER::StructAleWrapper>& GetStructAleStructureWrapperPtr()
    {
      return struct_ale_structure_wrapper_;
    };
    //@}

    /*! \brief Initialize this object

    Hand in all objects/parameters/etc. from outside.
    Construct and manipulate internal objects.

    \return int
    \date 08/16
    \author rauch    */
    virtual int Init(Teuchos::RCP<Structure> ti_strategy);

    /*! \brief Setup all class internal objects and members.

    \return void
    \date 11/16
    \author rauch  */
    virtual void Setup();

    /// set pointer to cellmigration model evaluator
    void SetModelEvaluatorPtr(Teuchos::RCP<::STR::MODELEVALUATOR::CellMigration> me)
    {
      cellmigration_model_evaluator_ = me;
      return;
    }

    /// access the cellmigration model evaluator
    Teuchos::RCP<::STR::MODELEVALUATOR::CellMigration> CellmigrationModelEvaluator()
    {
      return cellmigration_model_evaluator_;
    };

   private:
    //! @name Field wrappers
    //@{
    /// underlying structure of the SSI problem
    Teuchos::RCP<::ADAPTER::SSIStructureWrapper> ssi_structure_wrapper_;
    /// underlying structure of the SSI problem
    Teuchos::RCP<::ADAPTER::FSIStructureWrapperImmersed> fsi_structure_wrapper_;
    /// underlying structure-ale of the SSI problem
    Teuchos::RCP<::ADAPTER::StructAleWrapper> struct_ale_structure_wrapper_;
    //@}


    /// the cellmigration model evaluator object
    Teuchos::RCP<::STR::MODELEVALUATOR::CellMigration> cellmigration_model_evaluator_;

    //! flag indicating if class is setup
    bool issetup_;

    //! flag indicating if class is initialized
    bool isinit_;

   protected:
    //! returns true if Setup() was called and is still valid
    bool IsSetup() { return issetup_; };

    //! returns true if Init(..) was called and is still valid
    bool IsInit() { return isinit_; };

    //! check if \ref Setup() was called
    void CheckIsSetup()
    {
      if (not IsSetup()) dserror("Setup() was not called.");
    };

    //! check if \ref Init() was called
    void CheckIsInit()
    {
      if (not IsInit()) dserror("Init(...) was not called.");
    };

    //! set flag true after setup or false if setup became invalid
    void SetIsSetup(bool trueorfalse) { issetup_ = trueorfalse; };

    //! set flag true after init or false if init became invalid
    void SetIsInit(bool trueorfalse) { isinit_ = trueorfalse; };

  };  // class MultiphysicsStructureWrapperCellMigration
}  // namespace ADAPTER


#endif /* MULTIPHYSICS_STRUCTURE_WRAPPER_CELLMIGRATION_H */

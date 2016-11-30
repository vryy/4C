/*----------------------------------------------------------------------*/
/*!
\file ad_str_multiphysicswrapper_cellmigration.H

\brief Adapter wrapper for multiphysics in cell migration.

\maintainer Andreas Rauch

\level 1
*/
/*----------------------------------------------------------------------*/


#ifndef MULTIPHYSICS_STRUCTURE_WRAPPER_CELLMIGRATION_H
#define MULTIPHYSICS_STRUCTURE_WRAPPER_CELLMIGRATION_H


#include "ad_str_wrapper.H"

namespace ADAPTER
{

class SSIStructureWrapper;
class FSIStructureWrapperImmersed;

/*!
\brief This Multiphysics adapter wraps the subproblem specific structure wrappers.

*/
  class MultiphysicsStructureWrapperCellMigration : public StructureWrapper
  {
  public:

    /// constructor
    explicit MultiphysicsStructureWrapperCellMigration(Teuchos::RCP<Structure> structure);

    /// Destructor
    virtual ~MultiphysicsStructureWrapperCellMigration() { ; }

    /// return the structure wrapper for SSI
    const Teuchos::RCP< ::ADAPTER::SSIStructureWrapper>& GetSSIStructureWrapperPtr(){return ssi_structure_wrapper_;};

    /// return the structure wrapper for immersed FSI
    const Teuchos::RCP< ::ADAPTER::FSIStructureWrapperImmersed>& GetFSIStructureWrapperPtr(){return fsi_structure_wrapper_;};

    /*! \brief Initialize this object

    Hand in all objects/parameters/etc. from outside.
    Construct and manipulate internal objects.

    \note Try to only perform actions in Init(), which are still valid
          after parallel redistribution of discretizations.
          If you have to perform an action depending on the parallel
          distribution, make sure you adapt the affected objects after
          parallel redistribution.
          Example: cloning a discretization from another discretization is
          OK in Init(...). However, after redistribution of the source
          discretization do not forget to also redistribute the cloned
          discretization.
          All objects relying on the parallel distribution are supposed to
          the constructed in \ref Setup().

    \warning none
    \return int
    \date 08/16
    \author rauch  */
    virtual int Init(Teuchos::RCP<Structure> ti_strategy);

    /*! \brief Setup all class internal objects and members

     Setup() is not supposed to have any input arguments !

     Must only be called after Init().

     Construct all objects depending on the parallel distribution and
     relying on valid maps like, e.g. the state vectors, system matrices, etc.

     Call all Setup() routines on previously initialized internal objects and members.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, e.g. vectors may have wrong maps.

    \warning none
    \return void
    \date 11/16
    \author rauch  */
    virtual void Setup();

  private:

    /// underlying structure of the SSI problem
    Teuchos::RCP< ::ADAPTER::SSIStructureWrapper> ssi_structure_wrapper_;

    /// underlying structure of the SSI problem
    Teuchos::RCP< ::ADAPTER::FSIStructureWrapperImmersed> fsi_structure_wrapper_;

  private:

    //! flag indicating if class is setup
    bool issetup_;

    //! flag indicating if class is initialized
    bool isinit_;

  protected:

    //! returns true if Setup() was called and is still valid
    bool IsSetup()
    {return issetup_;};

    //! returns true if Init(..) was called and is still valid
    bool IsInit()
    {return isinit_;};

    //! check if \ref Setup() was called
    void CheckIsSetup()
    {if(not IsSetup()) dserror("Setup() was not called.");};

    //! check if \ref Init() was called
    void CheckIsInit()
    {if(not IsInit()) dserror("Init(...) was not called.");};

    //! set flag true after setup or false if setup became invalid
    void SetIsSetup(bool trueorfalse){issetup_ = trueorfalse;};

    //! set flag true after init or false if init became invalid
    void SetIsInit(bool trueorfalse){isinit_ = trueorfalse;};

  }; // class MultiphysicsStructureWrapperCellMigration
} // namespace ADAPTER


#endif /* MULTIPHYSICS_STRUCTURE_WRAPPER_CELLMIGRATION_H */

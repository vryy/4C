/*----------------------------------------------------------------------------*/
/*! \file

\brief Abstract interface for adapters to couple two discretizations

\level 1

*/
/*----------------------------------------------------------------------------*/

#ifndef COUPLING_ADAPTER_BASE_H
#define COUPLING_ADAPTER_BASE_H

/*----------------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace ADAPTER
{
  /*! \class CouplingBase
   *  \brief Abstract interface to coupling managers
   *
   *  This class is an interface to the coupling adapters.
   */
  class CouplingBase
  {
   public:
    /// empty constructor
    CouplingBase(){};

    /// virtual destructor
    virtual ~CouplingBase(){};

    /// @name Conversion between master and slave
    //!@{

    /// transfer a dof vector from master to slave
    virtual Teuchos::RCP<Epetra_Vector> MasterToSlave(
        Teuchos::RCP<Epetra_Vector> mv  ///< master vector (to be transferred)
    ) const = 0;

    /// transfer a dof vector from slave to master
    virtual Teuchos::RCP<Epetra_Vector> SlaveToMaster(
        Teuchos::RCP<Epetra_Vector> sv  ///< slave vector (to be transferred)
    ) const = 0;

    /// transfer a dof vector from master to slave
    virtual Teuchos::RCP<Epetra_MultiVector> MasterToSlave(
        Teuchos::RCP<Epetra_MultiVector> mv  ///< master vector (to be transferred)
    ) const = 0;

    /// transfer a dof vector from slave to master
    virtual Teuchos::RCP<Epetra_MultiVector> SlaveToMaster(
        Teuchos::RCP<Epetra_MultiVector> sv  ///< slave vector (to be transferred)
    ) const = 0;

    /// transfer a dof vector from master to slave
    virtual Teuchos::RCP<Epetra_Vector> MasterToSlave(
        Teuchos::RCP<const Epetra_Vector> mv  ///< master vector (to be transferred)
    ) const = 0;

    /// transfer a dof vector from slave to master
    virtual Teuchos::RCP<Epetra_Vector> SlaveToMaster(
        Teuchos::RCP<const Epetra_Vector> sv  ///< slave vector (to be transferred)
    ) const = 0;

    /// transfer a dof vector from master to slave
    virtual Teuchos::RCP<Epetra_MultiVector> MasterToSlave(
        Teuchos::RCP<const Epetra_MultiVector> mv  ///< master vector (to be transferred)
    ) const = 0;

    /// transfer a dof vector from slave to master
    virtual Teuchos::RCP<Epetra_MultiVector> SlaveToMaster(
        Teuchos::RCP<const Epetra_MultiVector> sv  ///< slave vector (to be transferred)
    ) const = 0;

    /// transfer a dof vector from master to slave
    virtual void MasterToSlave(
        Teuchos::RCP<const Epetra_MultiVector> mv,  ///< master vector (to be transferred)
        Teuchos::RCP<Epetra_MultiVector> sv         ///< slave vector (containing result)
    ) const = 0;

    /// transfer a dof vector from slave to master
    virtual void SlaveToMaster(
        Teuchos::RCP<const Epetra_MultiVector> sv,  ///< slave vector (to be transferred)
        Teuchos::RCP<Epetra_MultiVector> mv         ///< master vector (containing result)
    ) const = 0;

    //!@}

    //! @name Coupled maps
    //!@{

    /// the interface dof map of the master side
    virtual Teuchos::RCP<const Epetra_Map> MasterDofMap() const = 0;

    /// the interface dof map of the slave side
    virtual Teuchos::RCP<const Epetra_Map> SlaveDofMap() const = 0;

    //!@}
  };
}  // namespace ADAPTER

#endif  // COUPLING_ADAPTER_BASE_H

/*----------------------------------------------------------------------------*/
/*! \file

\brief Abstract interface for adapters to couple two discretizations

\level 1

*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_COUPLING_ADAPTER_BASE_HPP
#define FOUR_C_COUPLING_ADAPTER_BASE_HPP

/*----------------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace Core::Adapter
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
    virtual ~CouplingBase() = default;

    /// @name Conversion between master and slave
    //!@{

    /// transfer a dof vector from master to slave
    virtual Teuchos::RCP<Epetra_Vector> master_to_slave(
        Teuchos::RCP<Epetra_Vector> mv  ///< master vector (to be transferred)
    ) const = 0;

    /// transfer a dof vector from slave to master
    virtual Teuchos::RCP<Epetra_Vector> slave_to_master(
        Teuchos::RCP<Epetra_Vector> sv  ///< slave vector (to be transferred)
    ) const = 0;

    /// transfer a dof vector from master to slave
    virtual Teuchos::RCP<Epetra_MultiVector> master_to_slave(
        Teuchos::RCP<Epetra_MultiVector> mv  ///< master vector (to be transferred)
    ) const = 0;

    /// transfer a dof vector from slave to master
    virtual Teuchos::RCP<Epetra_MultiVector> slave_to_master(
        Teuchos::RCP<Epetra_MultiVector> sv  ///< slave vector (to be transferred)
    ) const = 0;

    /// transfer a dof vector from master to slave
    virtual Teuchos::RCP<Epetra_Vector> master_to_slave(
        Teuchos::RCP<const Epetra_Vector> mv  ///< master vector (to be transferred)
    ) const = 0;

    /// transfer a dof vector from slave to master
    virtual Teuchos::RCP<Epetra_Vector> slave_to_master(
        Teuchos::RCP<const Epetra_Vector> sv  ///< slave vector (to be transferred)
    ) const = 0;

    /// transfer a dof vector from master to slave
    virtual Teuchos::RCP<Epetra_MultiVector> master_to_slave(
        Teuchos::RCP<const Epetra_MultiVector> mv  ///< master vector (to be transferred)
    ) const = 0;

    /// transfer a dof vector from slave to master
    virtual Teuchos::RCP<Epetra_MultiVector> slave_to_master(
        Teuchos::RCP<const Epetra_MultiVector> sv  ///< slave vector (to be transferred)
    ) const = 0;

    /// transfer a dof vector from master to slave
    virtual void master_to_slave(
        Teuchos::RCP<const Epetra_MultiVector> mv,  ///< master vector (to be transferred)
        Teuchos::RCP<Epetra_MultiVector> sv         ///< slave vector (containing result)
    ) const = 0;

    /// transfer a dof vector from slave to master
    virtual void slave_to_master(
        Teuchos::RCP<const Epetra_MultiVector> sv,  ///< slave vector (to be transferred)
        Teuchos::RCP<Epetra_MultiVector> mv         ///< master vector (containing result)
    ) const = 0;

    //!@}

    //! @name Coupled maps
    //!@{

    /// the interface dof map of the master side
    virtual Teuchos::RCP<const Epetra_Map> master_dof_map() const = 0;

    /// the interface dof map of the slave side
    virtual Teuchos::RCP<const Epetra_Map> slave_dof_map() const = 0;

    //!@}
  };
}  // namespace Core::Adapter

FOUR_C_NAMESPACE_CLOSE

#endif

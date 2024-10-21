#ifndef FOUR_C_COUPLING_ADAPTER_BASE_HPP
#define FOUR_C_COUPLING_ADAPTER_BASE_HPP

/*----------------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"

#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace Coupling::Adapter
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
    virtual Teuchos::RCP<Core::LinAlg::Vector<double>> master_to_slave(
        const Core::LinAlg::Vector<double>& mv  ///< master vector (to be transferred)
    ) const = 0;

    /// transfer a dof vector from slave to master
    virtual Teuchos::RCP<Core::LinAlg::Vector<double>> slave_to_master(
        const Core::LinAlg::Vector<double>& sv  ///< slave vector (to be transferred)
    ) const = 0;

    /// transfer a dof vector from master to slave
    virtual Teuchos::RCP<Core::LinAlg::MultiVector<double>> master_to_slave(
        const Core::LinAlg::MultiVector<double>& mv  ///< master vector (to be transferred)
    ) const = 0;

    /// transfer a dof vector from slave to master
    virtual Teuchos::RCP<Core::LinAlg::MultiVector<double>> slave_to_master(
        const Core::LinAlg::MultiVector<double>& sv  ///< slave vector (to be transferred)
    ) const = 0;

    /// transfer a dof vector from master to slave
    virtual void master_to_slave(
        const Core::LinAlg::MultiVector<double>& mv,  ///< master vector (to be transferred)
        Core::LinAlg::MultiVector<double>& sv         ///< slave vector (containing result)
    ) const = 0;

    /// transfer a dof vector from slave to master
    virtual void slave_to_master(
        const Core::LinAlg::MultiVector<double>& sv,  ///< slave vector (to be transferred)
        Core::LinAlg::MultiVector<double>& mv         ///< master vector (containing result)
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
}  // namespace Coupling::Adapter

FOUR_C_NAMESPACE_CLOSE

#endif

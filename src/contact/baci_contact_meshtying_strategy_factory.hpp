/*---------------------------------------------------------------------*/
/*! \file
\brief Factory to create the desired meshtying strategy.


\level 2
*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_MESHTYING_STRATEGY_FACTORY_HPP
#define FOUR_C_CONTACT_MESHTYING_STRATEGY_FACTORY_HPP

#include "baci_config.hpp"

#include "baci_inpar_contact.hpp"
#include "baci_mortar_element.hpp"
#include "baci_mortar_strategy_factory.hpp"

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

BACI_NAMESPACE_OPEN

// forward declarations
namespace CONTACT
{
  class MtAbstractStrategy;
}  // namespace CONTACT

namespace MORTAR
{
  class InterfaceDataContainer;
  class MortarAbstractStrategy;
  class Element;
  class Interface;
  class StratDataContainer;

  namespace STRATEGY
  {
    /*! \brief Factory to create meshtying strategy
     *
     */
    class FactoryMT : public MORTAR::STRATEGY::Factory
    {
     public:
      //! constructor
      FactoryMT() = default;

      //! derived
      void Setup() override;

      /*! \brief Read and check meshtying/contact input parameters
       *
       * All specified contact-related input parameters are read from the
       * GLOBAL::Problem::Instance() and stored into a local variable of
       * type Teuchos::ParameterList. Invalid parameter combinations are
       * sorted out and throw a dserror.
       *
       * \param[in/out] params ParameterList with meshtying/contact parameters from input file
       *
       * \author Popp
       */
      void ReadAndCheckInput(Teuchos::ParameterList& params) const;

      /** \brief Create the meshtying interfaces
       *
       * \param[in] params ParameterList with mortar meshtying/contact parameters from input file
       * \param[in/out] interfaces Collection of all mortar interfaces
       * \param poroslave
       * \param poromaster
       *
       * \todo ToDo Get rid of poroslave and poromaster parameters.
       *
       * \author Popp
       */
      void BuildInterfaces(const Teuchos::ParameterList& params,
          std::vector<Teuchos::RCP<MORTAR::Interface>>& interfaces, bool& poroslave,
          bool& poromaster) const;

      /*! \brief Create the solver strategy object and pass all necessary data to it
       *
       * \param[in] params ParameterList with mortar meshtying/contact parameters from input file
       * \param[in] poroslave
       * \param[in] poromaster
       * \param[in] dof_offset
       * \param interfaces Collection of all mortar interfaces
       *
       * \todo ToDo Get rid of poroslave and poromaster parameters.
       *
       * \author Popp */
      Teuchos::RCP<CONTACT::MtAbstractStrategy> BuildStrategy(const Teuchos::ParameterList& params,
          const bool& poroslave, const bool& poromaster, const int& dof_offset,
          std::vector<Teuchos::RCP<MORTAR::Interface>>& interfaces) const;

      /*! \brief Create the solver strategy object and pass all necessary data to it
       *
       * \param[in] stype Type of solution strategy
       * \param[in] params ParameterList with mortar meshtying/contact parameters from input file
       * \param[in] poroslave
       * \param[in] poromaster
       * \param[in] dof_offset
       * \param[in] dof_row_map Dof row map
       * \param[in] node_row_map Node row map
       * \param[in] dim Spatial dimension
       * \param[in] comm_ptr Communicator
       * \param data_ptr Strategy data container
       *
       * \todo ToDo Get rid of poroslave and poromaster parameters.
       *
       * \note This routine can be used like a non-member function. If you need
       * access to the class members, use the alternative call.
       *
       * \author hiermeier \date 03/17 */
      static Teuchos::RCP<CONTACT::MtAbstractStrategy> BuildStrategy(
          const INPAR::CONTACT::SolvingStrategy stype, const Teuchos::ParameterList& params,
          const bool& poroslave, const bool& poromaster, const int& dof_offset,
          std::vector<Teuchos::RCP<MORTAR::Interface>>& interfaces, const Epetra_Map* dof_row_map,
          const Epetra_Map* node_row_map, const int dim,
          const Teuchos::RCP<const Epetra_Comm>& comm_ptr,
          Teuchos::RCP<MORTAR::StratDataContainer> data_ptr);

     protected:
     private:
    };  // class FactoryMT
  }     // namespace STRATEGY
}  // namespace MORTAR

BACI_NAMESPACE_CLOSE

#endif

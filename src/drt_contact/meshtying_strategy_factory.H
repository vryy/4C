/*---------------------------------------------------------------------*/
/*!
\brief Factory to create the desired meshtying strategy.

\maintainer Matthias Mayr

\level 2
*/
/*---------------------------------------------------------------------*/
#ifndef SRC_DRT_CONTACT_MESHTYING_STRATEGY_FACTORY_H_
#define SRC_DRT_CONTACT_MESHTYING_STRATEGY_FACTORY_H_

#include "../drt_inpar/inpar_contact.H"

#include "../drt_mortar/mortar_element.H"
#include "../drt_mortar/mortar_strategy_factory.H"

#include "../drt_xfem/xfem_enum_lists.H"

// forward declarations
namespace CONTACT
{
  class MtAbstractStrategy;
}  // namespace CONTACT

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

namespace MORTAR
{
  class InterfaceDataContainer;
  class MortarAbstractStrategy;
  class MortarElement;
  class MortarInterface;
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
      FactoryMT(){};

      //! destructor
      virtual ~FactoryMT(){};

      //! derived
      virtual void Setup();

      /*! \brief Read and check meshtying/contact input parameters
       *
       * All specified contact-related input parameters are read from the
       * DRT::Problem::Instance() and stored into a local variable of
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
          std::vector<Teuchos::RCP<MORTAR::MortarInterface>>& interfaces, bool& poroslave,
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
          std::vector<Teuchos::RCP<MORTAR::MortarInterface>>& interfaces) const;

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
          std::vector<Teuchos::RCP<MORTAR::MortarInterface>>& interfaces,
          const Epetra_Map* dof_row_map, const Epetra_Map* node_row_map, const int dim,
          const Teuchos::RCP<const Epetra_Comm>& comm_ptr,
          Teuchos::RCP<MORTAR::StratDataContainer> data_ptr);

     protected:
     private:
    };  // class FactoryMT
  }     // namespace STRATEGY
}  // namespace MORTAR

#endif /* SRC_DRT_CONTACT_MESHTYING_STRATEGY_FACTORY_H_ */

/*---------------------------------------------------------------------*/
/*!
\file meshtying_strategy_factory.H

\brief Factory to create the desired meshtying strategy.

\maintainer Alexander Seitz

\level 3

*/
/*---------------------------------------------------------------------*/


#ifndef SRC_DRT_CONTACT_MESHTYING_STRATEGY_FACTORY_H_
#define SRC_DRT_CONTACT_MESHTYING_STRATEGY_FACTORY_H_

#include "../drt_mortar/mortar_element.H"
#include "../drt_mortar/mortar_strategy_factory.H"
#include "../drt_xfem/xfem_enum_lists.H"
#include "../drt_inpar/inpar_contact.H"

// forward declarations...
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos
namespace CONTACT
{
  class MtAbstractStrategy;
}
namespace MORTAR
{
  class MortarInterface;
  class MortarElement;
  class MortarAbstractStrategy;
  class StratDataContainer;
  class IDataContainer;
  namespace STRATEGY
  {
    class FactoryMT : public MORTAR::STRATEGY::Factory
    {
     public:
      //! constructor
      FactoryMT(){};

      //! destructor
      virtual ~FactoryMT(){};

      //! derived
      virtual void Setup();

      /*! \brief Read and check contact input parameters
       *
       *  All specified contact-related input parameters are read from the
       *  DRT::Problem::Instance() and stored into a local variable of
       *  type Teuchos::ParameterList. Invalid parameter combinations are
       *  sorted out and throw a dserror.
       *
       *  \author Popp */
      void ReadAndCheckInput(Teuchos::ParameterList& cparams) const;

      /** \brief Create the contact interfaces
       *
       *  \author Popp */
      void BuildInterfaces(const Teuchos::ParameterList& mtparams,
          std::vector<Teuchos::RCP<MORTAR::MortarInterface>>& interfaces, bool& poroslave,
          bool& poromaster) const;

      /*! \brief Create the solver strategy object and pass all necessary data to it
       *
       * \author Popp */
      Teuchos::RCP<CONTACT::MtAbstractStrategy> BuildStrategy(const Teuchos::ParameterList& cparams,
          const bool& poroslave, const bool& poromaster, const int& dof_offset,
          std::vector<Teuchos::RCP<MORTAR::MortarInterface>>& interfaces) const;

      /*! \brief Create the solver strategy object and pass all necessary data to it
       *
       *  \note This routine can be used like a non-member function. If you need
       *  access to the class members, use the alternative call.
       *
       * \author hiermeier \date 03/17 */
      static Teuchos::RCP<CONTACT::MtAbstractStrategy> BuildStrategy(
          const INPAR::CONTACT::SolvingStrategy stype, const Teuchos::ParameterList& cparams,
          const bool& poroslave, const bool& poromaster, const int& dof_offset,
          std::vector<Teuchos::RCP<MORTAR::MortarInterface>>& interfaces,
          const Epetra_Map* dof_row_map, const Epetra_Map* node_row_map, const int dim,
          const Teuchos::RCP<const Epetra_Comm>& comm_ptr,
          Teuchos::RCP<MORTAR::StratDataContainer> data_ptr);

     protected:
    };  // class FactoryMT
  }     // namespace STRATEGY
}  // namespace MORTAR



#endif /* SRC_DRT_CONTACT_MESHTYING_STRATEGY_FACTORY_H_ */

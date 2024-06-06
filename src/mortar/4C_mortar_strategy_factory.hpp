/*---------------------------------------------------------------------*/
/*! \file
\brief Base class for the CONTACT/MESHTYING factories.

\level 2

*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_MORTAR_STRATEGY_FACTORY_HPP
#define FOUR_C_MORTAR_STRATEGY_FACTORY_HPP

#include "4C_config.hpp"

#include "4C_inpar_contact.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  class Discretization;
}  // namespace Discret

namespace Core::Nodes
{
  class Node;
}

namespace Core::Elements
{
  class Element;
}

namespace Core::IO
{
  class DiscretizationReader;
  class DiscretizationWriter;
}  // namespace Core::IO

namespace STR
{
  namespace TimeInt
  {
    class BaseDataGlobalState;
  }  // namespace TimeInt
}  // namespace STR

namespace Mortar
{
  class Element;
  class Interface;
  class Node;
  class StrategyBase;

  namespace STRATEGY
  {
    /*! \brief Base class for mortar meshtying/contact factories
     *
     */
    class Factory
    {
     public:
      /**
       * Virtual destructor.
       */
      virtual ~Factory() = default;

      //! constructor
      Factory();

      //! initialization of class variables
      void Init(const Teuchos::RCP<STR::TimeInt::BaseDataGlobalState>& gstate_ptr);
      void Init(Teuchos::RCP<Discret::Discretization> dis);

      /*! \brief Setup of class variables
       *
       * \note Since this is an abstract class, the setup flag stays false. It has to be set by
       * the Setup routing of the derived class.
       */
      virtual void Setup();

      /*! \brief print strategy banner
       *
       *  \param[in] soltype contact solving strategy type
       */
      static void PrintStrategyBanner(const enum Inpar::CONTACT::SolvingStrategy soltype);

      /*! \brief Create the desired search tree object
       *
       * We loop over all \c interfaces and build the search tree for every interface.
       *
       * @param[in] interfaces All meshtying interfaces
       */
      void BuildSearchTree(const std::vector<Teuchos::RCP<Mortar::Interface>>& interfaces) const;

      /*! \brief Check the problem dimension
       *
       * Mortar meshtying/contact problems are implemented for 2D and 3D only.
       * Throw an error in case of any other spatial dimension.
       */
      void CheckDimension() const;

     protected:
      /*! \brief Set #issetup_ flag with user-given value
       *
       * @param[in] issetup Setup flag (default = true)
       */
      inline void set_is_setup(const bool issetup = true) { issetup_ = issetup; };

      //! Returns true, if Init() has been called
      inline const bool& is_init() const { return isinit_; };

      //! Returns true, if Setup() has been called
      inline const bool& is_setup() const { return issetup_; };

      //! Checks, if Init() and Setup() have been called
      void check_init_setup() const;

      //! Checks if Init() has been called
      void check_init() const;

      //! @name NURBS related stuff
      //!@{

      /*! \brief Prepare mortar element for NURBS case
       *
       *  Stores knot vector, zerosized information and normal factor
       *
       *  \author Farah */
      void prepare_nurbs_element(const Discret::Discretization& discret,
          Teuchos::RCP<Core::Elements::Element> ele, Teuchos::RCP<Mortar::Element> cele) const;

      /*! \brief Prepare mortar node for NURBS case
       *
       *  Stores control point weight
       *
       *  \author Farah */
      void prepare_nurbs_node(
          const Core::Nodes::Node* node, Teuchos::RCP<Mortar::Node> mnode) const;

      //!@}

      //! @name Accessors
      //!@{

      //! Returns the global state data container (read-only, do not change this!!!)
      const STR::TimeInt::BaseDataGlobalState& g_state() const;

      //! Returns the (structural) discretization
      Discret::Discretization& discret();
      const Discret::Discretization& discret() const;

      //! returns a reference to a copy of the structural communicator
      Epetra_Comm& comm();
      const Epetra_Comm& comm() const;
      Teuchos::RCP<Epetra_Comm> comm_ptr();
      Teuchos::RCP<const Epetra_Comm> comm_ptr() const;

      //! returns the problem dimension
      const int& dim() const;

      //!@}

      //! pointer to the structural problem discretization
      Teuchos::RCP<Discret::Discretization> discret_ptr_;

     private:
      //! @name Status flags
      //!@{

      //! init flag
      bool isinit_;

      //! setup flag
      bool issetup_;

      //!@}

      //! pointer to the structural global state data container
      Teuchos::RCP<STR::TimeInt::BaseDataGlobalState> gstate_ptr_;

      //! pointer to a COPY of the structural communicator
      Teuchos::RCP<Epetra_Comm> comm_ptr_;

      int dim_;
    };  // namespace STRATEGY
  }     // namespace STRATEGY
}  // namespace Mortar

FOUR_C_NAMESPACE_CLOSE

#endif

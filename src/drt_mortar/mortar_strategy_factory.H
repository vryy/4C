/*---------------------------------------------------------------------*/
/*!
\file mortar_strategy_factory.H

\brief Base class for the CONTACT/MESHTYING factories.

\level 3

\maintainer Michael Hiermeier

\date Feb 4, 2016

*/
/*---------------------------------------------------------------------*/
#ifndef SRC_DRT_MORTAR_MORTAR_STRATEGY_FACTORY_H_
#define SRC_DRT_MORTAR_MORTAR_STRATEGY_FACTORY_H_

#include <Teuchos_RCP.hpp>

// forward declarations
class Epetra_Vector;
class Epetra_Comm;
namespace STR
{
  namespace TIMINT
  {
    class BaseDataGlobalState;
  }  // namespace TIMINT
}  // namespace STR
namespace IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
}  // namespace IO
namespace DRT
{
  class Node;
  class Element;
  class DiscretizationInterface;
}  // namespace DRT
namespace MORTAR
{
  class MortarNode;
  class MortarElement;
  class StrategyBase;
  namespace STRATEGY
  {
    class Factory
    {
     public:
      //! constructor
      Factory();

      //! destructor
      virtual ~Factory(){};

      //! initialization of class variables
      void Init(const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr);
      void Init(Teuchos::RCP<DRT::DiscretizationInterface> dis);

      //! setup of class variables
      void Setup();

     protected:
      //! Returns true, if Init() has been called
      inline const bool& IsInit() const { return isinit_; };

      //! Returns true, if Setup() has been called
      inline const bool& IsSetup() const { return issetup_; };

      //! Checks, if Init() and Setup() have been called
      void CheckInitSetup() const;

      //! Checks if Init() has been called
      void CheckInit() const;

      //! @name NURBS related stuff
      //! @{
      /*! \brief Prepare mortar element for NURBS case
       *
       *  Stores knot vector, zerosized information and normal factor
       *
       *  \author Farah */
      void PrepareNURBSElement(const DRT::DiscretizationInterface& discret,
          Teuchos::RCP<DRT::Element> ele, Teuchos::RCP<MORTAR::MortarElement> cele) const;

      /*! \brief Prepare mortar node for NURBS case
       *
       *  Stores control point weight
       *
       *  \author Farah */
      void PrepareNURBSNode(const DRT::Node* node, Teuchos::RCP<MORTAR::MortarNode> mnode) const;
      //! @}
      //! @}

      //! @name accessors
      //! @{
      //! Returns the global state data container (read-only, do not change this!!!)
      const STR::TIMINT::BaseDataGlobalState& GState() const;

      //! Returns the (structural) discretization
      DRT::DiscretizationInterface& Discret();
      const DRT::DiscretizationInterface& Discret() const;

      //! returns a reference to a copy of the structural communicator
      Epetra_Comm& Comm();
      const Epetra_Comm& Comm() const;
      Teuchos::RCP<Epetra_Comm> CommPtr();
      Teuchos::RCP<const Epetra_Comm> CommPtr() const;

      //! returns the problem dimension
      const int& Dim() const;
      //! @}Manager
     protected:
      //! init flag
      bool isinit_;

      //! setup flag
      bool issetup_;

     private:
      //! pointer to the structural global state data container
      Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> gstate_ptr_;

      //! pointer to the structural problem discretization
      Teuchos::RCP<DRT::DiscretizationInterface> discret_ptr_;

      //! pointer to a COPY of the structural communicator
      Teuchos::RCP<Epetra_Comm> comm_ptr_;

      int dim_;
    };  // class Factory
  }     // namespace STRATEGY
}  // namespace MORTAR


#endif /* SRC_DRT_MORTAR_MORTAR_STRATEGY_FACTORY_H_ */

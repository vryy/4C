/*---------------------------------------------------------------------*/
/*!
\brief Base class for the CONTACT/MESHTYING factories.

\level 2

\maintainer Matthias Mayr
*/
/*---------------------------------------------------------------------*/
#ifndef SRC_DRT_MORTAR_MORTAR_STRATEGY_FACTORY_H_
#define SRC_DRT_MORTAR_MORTAR_STRATEGY_FACTORY_H_

#include "../drt_inpar/inpar_contact.H"

#include <Teuchos_RCP.hpp>

// forward declarations
class Epetra_Comm;
class Epetra_Vector;

namespace DRT
{
  class DiscretizationInterface;
  class Element;
  class Node;
}  // namespace DRT

namespace IO
{
  class DiscretizationReader;
  class DiscretizationWriter;
}  // namespace IO

namespace STR
{
  namespace TIMINT
  {
    class BaseDataGlobalState;
  }  // namespace TIMINT
}  // namespace STR

namespace MORTAR
{
  class MortarElement;
  class MortarInterface;
  class MortarNode;
  class StrategyBase;

  namespace STRATEGY
  {
    /*! \brief Base class for mortar meshtying/contact factories
     *
     */
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

      /*! \brief Setup of class variables
       *
       * \note Since this is an abstract class, the setup flag stays false. It has to be set by
       * the Setup routing of the derived class.
       */
      void Setup();

      /*! \brief print strategy banner
       *
       *  \param[in] soltype contact solving strategy type
       */
      static void PrintStrategyBanner(const enum INPAR::CONTACT::SolvingStrategy soltype);

      /*! \brief Create the desired search tree object
       *
       * We loop over all \c interfaces and build the search tree for every interface.
       *
       * @param[in] interfaces All meshtying interfaces
       */
      void BuildSearchTree(
          const std::vector<Teuchos::RCP<MORTAR::MortarInterface>>& interfaces) const;

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
      inline void SetIsSetup(const bool issetup = true) { issetup_ = issetup; };

      //! Returns true, if Init() has been called
      inline const bool& IsInit() const { return isinit_; };

      //! Returns true, if Setup() has been called
      inline const bool& IsSetup() const { return issetup_; };

      //! Checks, if Init() and Setup() have been called
      void CheckInitSetup() const;

      //! Checks if Init() has been called
      void CheckInit() const;

      //! @name NURBS related stuff
      //!@{

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

      //!@}

      //! @name Accessors
      //!@{

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

      //!@}

      //! pointer to the structural problem discretization
      Teuchos::RCP<DRT::DiscretizationInterface> discret_ptr_;

     private:
      //! @name Status flags
      //!@{

      //! init flag
      bool isinit_;

      //! setup flag
      bool issetup_;

      //!@}

      //! pointer to the structural global state data container
      Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> gstate_ptr_;

      //! pointer to a COPY of the structural communicator
      Teuchos::RCP<Epetra_Comm> comm_ptr_;

      int dim_;
    };  // namespace STRATEGY
  }     // namespace STRATEGY
}  // namespace MORTAR

#endif /* SRC_DRT_MORTAR_MORTAR_STRATEGY_FACTORY_H_ */

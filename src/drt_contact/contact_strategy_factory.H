/*---------------------------------------------------------------------*/
/*!
\file contact_strategy_factory.H

\brief Factory to create the desired contact strategy.

\maintainer Michael Hiermeier

\date Feb 4, 2016

\level 3

*/
/*---------------------------------------------------------------------*/


#ifndef SRC_DRT_CONTACT_CONTACT_STRATEGY_FACTORY_H_
#define SRC_DRT_CONTACT_CONTACT_STRATEGY_FACTORY_H_

#include "../drt_mortar/mortar_element.H"
#include "../drt_mortar/mortar_strategy_factory.H"

// forward declarations...
namespace Teuchos {
  class ParameterList;
} // namespace Teuchos
namespace CONTACT {
  class CoInterface;
  class CoElement;
  class CoAbstractStrategy;
namespace STRATEGY {
  class Factory : public MORTAR::STRATEGY::Factory
  {
  public:
    //! constructor
    Factory() {};

    //! destructor
    virtual ~Factory() {};

    virtual void Setup();

    //! Check the problem dimension
    void CheckDimension() const;

    /*! \brief Read and check contact input parameters
     *
     *  All specified contact-related input parameters are read from the
     *  DRT::Problem::Instance() and stored into a local variable of
     *  type Teuchos::ParameterList. Invalid parameter combinations are
     *  sorted out and throw a dserror.
     *
     *  \author Popp */
    void ReadAndCheckInput(Teuchos::ParameterList& cparams) const;

    /* \brief Create the contact interfaces
     *
     * \author Popp */
    void BuildInterfaces(const Teuchos::ParameterList& cparams,
        std::vector<Teuchos::RCP<CONTACT::CoInterface> >& interfaces,
        bool& poroslave,
        bool& poromaster) const;

    /*! \brief Create the solver strategy object and pass all necessary data to it
     *
     * \author Popp */
    Teuchos::RCP<CONTACT::CoAbstractStrategy> BuildStrategy(
        const Teuchos::ParameterList& cparams,
        const bool& poroslave,
        const bool& poromaster,
        const int& dof_offset,
        std::vector<Teuchos::RCP<CONTACT::CoInterface> >& interfaces) const;

    //! Create the desired search tree object
    void BuildSearchTree(
        std::vector<Teuchos::RCP<CONTACT::CoInterface> >& interfaces) const;

    //! print some final screen output
    void Print(
        std::vector<Teuchos::RCP<CONTACT::CoInterface> >& interfaces,
        const Teuchos::RCP<CONTACT::CoAbstractStrategy>& strategy_ptr,
        const Teuchos::ParameterList& cparams) const;

  protected:
    /*! \brief Set Parent Elements for Poro Face Elements
     *
     *  \author Ager */
    void SetPoroParentElement(
        enum MORTAR::MortarElement::PhysicalType& slavetype,
        enum MORTAR::MortarElement::PhysicalType& mastertype,
        Teuchos::RCP<CONTACT::CoElement>& cele,
        Teuchos::RCP<DRT::Element>& ele) const;

    /*! \brief Find Physical Type (Poro or Structure) of Poro Interface
     *
     *  \author Ager */
    void FindPoroInterfaceTypes(
        bool& poromaster,
        bool& poroslave,
        bool& structmaster,
        bool& structslave,
        enum MORTAR::MortarElement::PhysicalType& slavetype,
        enum MORTAR::MortarElement::PhysicalType& mastertype) const;

  };  // class Factory
} // namespace STRATEGY
} // namespace CONTACT



#endif /* SRC_DRT_CONTACT_CONTACT_STRATEGY_FACTORY_H_ */

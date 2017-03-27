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
#include "../drt_xfem/xfem_enum_lists.H"
#include "../drt_inpar/inpar_contact.H"

// forward declarations...
namespace Teuchos {
  class ParameterList;
} // namespace Teuchos
namespace CONTACT {
  class CoInterface;
  class CoElement;
  class CoAbstractStrategy;
  class AbstractStratDataContainer;
  class IDataContainer;
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

    /** \brief Create the contact interfaces
     *
     *  \author Popp */
    void BuildInterfaces(const Teuchos::ParameterList& cparams,
        std::vector<Teuchos::RCP<CONTACT::CoInterface> >& interfaces,
        bool& poroslave,
        bool& poromaster) const;

    /** \brief Create a contact interface object based on the given information
     *
     *  \author hiermeier \date 03/17 */
    static Teuchos::RCP< ::CONTACT::CoInterface> CreateInterface(
        const int id,
        const Epetra_Comm& comm,
        const int dim,
        Teuchos::ParameterList& icparams,
        const bool selfcontact,
        const enum INPAR::MORTAR::RedundantStorage redundant,
        const Teuchos::RCP<std::pair<enum XFEM::FieldName,
            Teuchos::RCP<const DRT::DiscretizationInterface> > >& parent_dis_pair,
        Teuchos::RCP<CONTACT::IDataContainer> idata_ptr = Teuchos::null );

    /** \brief Create a contact interface object based on the given information
     *
     *  \author hiermeier \date 03/17 */
    static Teuchos::RCP< ::CONTACT::CoInterface> CreateInterface(
        const enum INPAR::CONTACT::SolvingStrategy stype,
        const int id,
        const Epetra_Comm& comm,
        const int dim,
        Teuchos::ParameterList& icparams,
        const bool selfcontact,
        const enum INPAR::MORTAR::RedundantStorage redundant,
        const Teuchos::RCP<std::pair<enum XFEM::FieldName,
            Teuchos::RCP<const DRT::DiscretizationInterface> > >& parent_dis_pair,
        Teuchos::RCP<CONTACT::IDataContainer> idata_ptr );

    /*! \brief Create the solver strategy object and pass all necessary data to it
     *
     * \author Popp */
    Teuchos::RCP<CONTACT::CoAbstractStrategy> BuildStrategy(
        const Teuchos::ParameterList& cparams,
        const bool& poroslave,
        const bool& poromaster,
        const int& dof_offset,
        std::vector<Teuchos::RCP<CONTACT::CoInterface> >& interfaces) const;

    /*! \brief Create the solver strategy object and pass all necessary data to it
     *
     *  \note This routine can be used like a non-member function. If you need
     *  access to the class members, use the alternative call.
     *
     * \author hiermeier \date 03/17 */
    static Teuchos::RCP<CONTACT::CoAbstractStrategy> BuildStrategy(
        const INPAR::CONTACT::SolvingStrategy stype,
        const Teuchos::ParameterList& cparams,
        const bool& poroslave,
        const bool& poromaster,
        const int& dof_offset,
        std::vector<Teuchos::RCP<CONTACT::CoInterface> >& interfaces,
        const Epetra_Map* dof_row_map,
        const Epetra_Map* node_row_map,
        const int dim,
        const Teuchos::RCP<const Epetra_Comm>& comm_ptr,
        Teuchos::RCP<CONTACT::AbstractStratDataContainer> data_ptr );

    //! Create the desired search tree object
    void BuildSearchTree(
        const std::vector<Teuchos::RCP<CONTACT::CoInterface> >& interfaces) const;

    //! print some final screen output
    void Print(
        std::vector<Teuchos::RCP<CONTACT::CoInterface> >& interfaces,
        const Teuchos::RCP<CONTACT::CoAbstractStrategy>& strategy_ptr,
        const Teuchos::ParameterList& cparams) const;

    /*! \brief print strategy banner
     *
     *  \param soltype (in) : contact solving strategy type */
    static void PrintStrategyBanner(
        const enum INPAR::CONTACT::SolvingStrategy soltype );

  protected:

    /*! \brief Set Parent Elements for Poro Face Elements
     *
     *  \author Ager */
    void SetPoroParentElement(
        enum MORTAR::MortarElement::PhysicalType& slavetype,
        enum MORTAR::MortarElement::PhysicalType& mastertype,
        Teuchos::RCP<CONTACT::CoElement>& cele,
        Teuchos::RCP<DRT::Element>& ele,
        const DRT::DiscretizationInterface & discret) const;

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

    /*! \brief Extract the discretization belonging to the given contact
     *   condition pair/group
     *
     *  \param full_discret (in) : full structural discretization
     *  \param given_cgroup (in) : given contact condition group/pair
     *
     *  \author hiermeier
     *  \date 10/16 */
    void ExtractParentDiscretization(
        const DRT::DiscretizationInterface & full_discret,
        const std::vector<DRT::Condition*> & given_ccgroup,
        std::pair<enum XFEM::FieldName, Teuchos::RCP<const DRT::DiscretizationInterface> > & xdispair) const;
  };  // class Factory
} // namespace STRATEGY
} // namespace CONTACT



#endif /* SRC_DRT_CONTACT_CONTACT_STRATEGY_FACTORY_H_ */

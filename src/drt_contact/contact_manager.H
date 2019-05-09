/*----------------------------------------------------------------------*/
/*!
\brief BACI implementation of main class to control all contact

\level 1

\maintainer Matthias Mayr

*/
/*----------------------------------------------------------------------*/
#ifndef CONTACT_MANAGER_H
#define CONTACT_MANAGER_H

#include "../drt_mortar/mortar_manager_base.H"


// forward declarations
namespace DRT
{
  class Discretization;
}

namespace CONTACT
{
  // forward declarations
  class CoElement;

  /*!
  \brief BACI implementation of main class to control all contact

  \author popp (popp@lnm.mw.tum.de)
  */

  class CoManager : public MORTAR::ManagerBase
  {
   public:
    /*!
    \brief Standard Constructor

    The constructor takes a discretization  that is expected to have at least
    two contact boundary conditions. It extracts all contact boundary conditions
    and constructs one or multiple contact interfaces from them and stores them.

    It calls CONTACT::CoInterface::FillComplete() on all contact interfaces which
    makes the nodes and elements of a contact interfaces redundant on all
    processors that either own a node or an element on the interfaces in the
    input discretization.

    In addition, it creates the necessary solver strategy object which handles
    the whole contact evaluation.

    \param discret (in): A discretization containing contact boundary conditions
    \param alphaf (in): Generalized-alpha parameter (set to 0.0 by default)

    */
    CoManager(DRT::Discretization& discret, double alphaf = 0.0);

    /*!
    \brief Destructor

    */
    virtual ~CoManager() = default;

    //! @name Access methods
    //! @{

    /*!
    \brief Get discretization

    */
    const DRT::Discretization& Discret() const { return discret_; };

    //! @}

    //! @name Evaluation methods
    //! @{

    /*!
    \brief Write restart information for contact

    The additionally necessary restart information in the contact
    case are the current Lagrange multiplier values and the current
    active set status of each slave node.

    \param output (in): IO::Discretization writer for restart

    */
    void WriteRestart(IO::DiscretizationWriter& output, bool forcedrestart = false) final;

    /*!
    \brief Read restart information for contact

    This method has the inverse functionality of WriteRestart, as
    it reads the activetoggle / lmold vectors and stores the restart
    status into each slave node. Moreover, all global maps concerning
    the active set and the old mortar matrices D,M are rebuilt based
    on the restart information.

    \param reader (in): IO::Discretization reader for restart
    \param dis (in)   : global dof displacement vector
    \param zero (in)  : global dof zero vector

    */
    void ReadRestart(IO::DiscretizationReader& reader, Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<Epetra_Vector> zero) final;

    /*!
    \brief Write interface quantities for postprocessing

    \param output (in): IO::Discretization writer for restart

    */
    void PostprocessQuantities(IO::DiscretizationWriter& output) final;

    /*!
    \brief Reconnect Contact Element -- Parent Element Pointers

    (As during the Restart the initial created structural elements are destructed and created again,
    the pointer of these elements changes and therefore needs to be reconnected)
    */
    void ReconnectParentElements();

    /*!
    \brief Set Parent Elements for Poro Face Elements

    \param output (in): int& type of slave/master elements --> = (-1; //1 poro, 0 struct, -1
    default) Teuchos::RCP<CONTACT::CoElement>& Reference to pointer of contact face element!
                        Teuchos::RCP<DRT::Element>& Reference to pointer of contact parent element!

    */
    void SetPoroParentElement(int& slavetype, int& mastertype,
        Teuchos::RCP<CONTACT::CoElement>& cele, Teuchos::RCP<DRT::Element>& ele);

    /*!
    \brief Find Physical Type (Poro or Structure) of Poro Interface

    \param output (in): bool& poromaster, bool& poroslave,
        bool& structmaster, bool& structslave, int& slavetype, int& mastertype

    */
    void FindPoroInterfaceTypes(bool& poromaster, bool& poroslave, bool& structmaster,
        bool& structslave, int& slavetype, int& mastertype);

    //! @}

   private:
    /*!
    \brief Read and check contact input parameters

    All specified contact-related input parameters are read from the
    DRT::Problem::Instance() and stored into a local variable of
    type Teuchos::ParameterList. Invalid parameter combinations are
    sorted out and throw a dserror.

    */
    bool ReadAndCheckInput(Teuchos::ParameterList& cparams);

   protected:
    // don't want = operator and cctor
    CoManager operator=(const CoManager& old) = delete;
    CoManager(const CoManager& old) = delete;

    //! the underlying problem discretization
    DRT::Discretization& discret_;

  };  // class CoManager
}  // namespace CONTACT

#endif  // #ifndef CONTACT_MANAGER_H

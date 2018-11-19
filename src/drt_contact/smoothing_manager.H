/*---------------------------------------------------------------------*/
/*!
\brief Manager for the smoothing contact approach.

\level 3

\maintainer Alexander Popp

*/
/*---------------------------------------------------------------------*/
#ifndef SMOOTHING_MANAGER_H_
#define SMOOTHING_MANAGER_H_

#include "../drt_mortar/mortar_manager_base.H"

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace CONTACT
{
  class SmoothingManager : public MORTAR::ManagerBase
  {
   public:
    /*!
    \brief Constructor

    \param discret (in): A discretization containing contact boundary conditions
    \param alphaf (in): Generalized-alpha parameter (set to 0.0 by default)
    */
    SmoothingManager(DRT::Discretization& discret, double alphaf = 0.0);

    /*!
    \brief Destructor

    */
    virtual ~SmoothingManager() = default;

    //! @name Access methods

    /*!
    \brief Get discretization

    */
    const DRT::Discretization& Discret() const { return discret_; };

    //@}

    //! @name Evaluation methods

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
    \brief Write interface tractions for postprocessing

    \param output (in): IO::Discretization writer for restart

    */
    void PostprocessQuantities(IO::DiscretizationWriter& output) final;

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

    // don't want = operator and cctor
    SmoothingManager operator=(const SmoothingManager& old) = delete;
    SmoothingManager(const SmoothingManager& old) = delete;

    //! Underlying problem discretization
    DRT::Discretization& discret_;

  };  // class SmoothingManager
}  // namespace CONTACT

#endif /* SMOOTHING_MANAGER_H_ */

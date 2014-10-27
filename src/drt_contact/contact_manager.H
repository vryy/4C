/*!----------------------------------------------------------------------
\file contact_manager.H

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/
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
  virtual ~CoManager() {};

  //! @name Access methods

  /*!
  \brief Get discretization

  */
  DRT::Discretization& Discret() const { return discret_;};

  //@}

  //! @name Evaluation methods

  /*!
  \brief Read and check contact input parameters

  All specified contact-related input parameters are read from the
  DRT::Problem::Instance() and stored into a local variable of
  type Teuchos::ParameterList. Invalid parameter combinations are
  sorted out and throw a dserror.

  */
  bool ReadAndCheckInput(Teuchos::ParameterList& cparams);

  /*!
  \brief Write restart information for contact

  The additionally necessary restart information in the contact
  case are the current Lagrange multiplier values and the current
  active set status of each slave node.

  \param output (in): IO::Discretization writer for restart

  */
  void WriteRestart(IO::DiscretizationWriter& output, bool forcedrestart = false);

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
  void ReadRestart(IO::DiscretizationReader& reader,
                   Teuchos::RCP<Epetra_Vector> dis,
                   Teuchos::RCP<Epetra_Vector> zero);

  /*!
  \brief Write interface tractions for postprocessing

  \param output (in): IO::Discretization writer for restart

  */
  void PostprocessTractions(IO::DiscretizationWriter& output);

  //@}


protected:

  // don't want = operator and cctor
  CoManager operator = (const CoManager& old);
  CoManager(const CoManager& old);

  DRT::Discretization&  discret_; // the underlying problem discretization

}; // class CoManager
} // namespace CONTACT

#endif  // #ifndef CONTACT_MANAGER_H

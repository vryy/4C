/*!----------------------------------------------------------------------
\file meshtying_manager.H

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*-----------------------------------------------------------------------*/
#ifndef MESHTYING_MANAGER_H
#define MESHTYING_MANAGER_H

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
  \brief BACI implementation of main class to control all meshtying

  \author popp (popp@lnm.mw.tum.de)
  */
  class MtManager : public MORTAR::ManagerBase
  {
   public:
    /*!
    \brief Standard Constructor

    The constructor takes a discretization that is expected to have at least
    two meshtying boundary conditions. It extracts all meshtying boundary conditions
    and constructs one or multiple meshtying interfaces from them and stores them.

    It calls MORTAR::MortarInterface::FillComplete() on all meshtying interfaces which
    makes the nodes and elements of a meshtying interfaces redundant on all
    processors that either own a node or an element on the interfaces in the
    input discretization.

    In addition, it creates the necessary solver strategy object which handles
    the whole meshtying evaluation.

    \param discret (in): A discretization containing meshtying boundary conditions
    \param alphaf (in): Generalized-alpha parameter (set to 0.0 by default)

    */
    MtManager(DRT::Discretization& discret, double alphaf = 0.0);

    /*!
    \brief Destructor

    */
    virtual ~MtManager(){};

    //! @name Access methods

    /*!
    \brief Get discretization

    */
    DRT::Discretization& Discret() const { return discret_; };

    //@}

    //! @name Evaluation methods

    /*!
    \brief Read and check input parameters

    All specified meshtying-related input parameters are read from the
    DRT::Problem::Instance() and stored into a local variable of
    type Teuchos::ParameterList. Invalid parameter combinations are
    sorted out and throw a dserror.

    */
    bool ReadAndCheckInput(Teuchos::ParameterList& mtparams);

    /*!
    \brief Write restart information for meshtying

    The additionally necessary restart information in the meshtying
    case are the current Lagrange multiplier values.

    \param output (in): IO::Discretization writer for restart

    */
    void WriteRestart(IO::DiscretizationWriter& output, bool forcedrestart = false);

    /*!
    \brief Read restart information for contact

    This method has the inverse functionality of WriteRestart, as
    it reads the restart Lagrange mulitplier vectors. Moreover,
    all mortar coupling quantities (e.g. D and M) have to be
    re-computed upon restart..

    \param reader (in): IO::Discretization reader for restart
    \param dis (in)   : global dof displacement vector
    \param zero (in)  : global dof zero vector

    */
    void ReadRestart(IO::DiscretizationReader& reader, Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<Epetra_Vector> zero);

    /*!
    \brief Write interface tractions for postprocessing

    \param output (in): IO::Discretization writer for restart

    */
    void PostprocessQuantities(IO::DiscretizationWriter& output);

    //@}


   protected:
    // don't want = operator and cctor
    MtManager operator=(const MtManager& old);
    MtManager(const MtManager& old);

    DRT::Discretization& discret_;  // the underlying problem discretization

  };  // class MtManager
}  // namespace CONTACT

#endif  // #ifndef MESHTYING_MANAGER_H

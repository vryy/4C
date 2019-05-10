/*!----------------------------------------------------------------------
\file meshtying_contact_bridge.H

\level 2

<pre>
\maintainer Alexander Popp

</pre>

*-----------------------------------------------------------------------*/
#ifndef MESHTYING_CONTACT_BRIDGE_H_
#define MESHTYING_CONTACT_BRIDGE_H_

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "../drt_inpar/inpar_contact.H"
#include "Epetra_Comm.h"
#include "../drt_mortar/mortar_manager_base.H"  // TODO forward declaration!!!
#include "../drt_lib/drt_dserror.H"

// forward declarations
namespace DRT
{
  class Discretization;
  class Condition;
}  // namespace DRT

namespace LINALG
{
  class MapExtractor;
}

namespace MORTAR
{
  class StrategyBase;
}

namespace CONTACT
{
  // class CoManager;
  // class MtManager;

  class MeshtyingContactBridge
  {
   public:
    //@}

    MeshtyingContactBridge(DRT::Discretization& dis, std::vector<DRT::Condition*>& mtcond,
        std::vector<DRT::Condition*>& ccond, double alphaf, bool smoothing = false);

    /*!
    \brief Destructor

    */
    virtual ~MeshtyingContactBridge(){};

    //! @name Access methods

    /*!
    \brief Get Epetra communicator

    */
    const Epetra_Comm& Comm() const
    {
      if (cman_ != Teuchos::null) return cman_->Comm();
      if (mtman_ != Teuchos::null) return mtman_->Comm();
      dserror("can't get comm()");
      return cman_->Comm();
    }

    /*!
    \brief Get contact manager

    */
    Teuchos::RCP<MORTAR::ManagerBase> ContactManager()
    {
      if (sman_ != Teuchos::null)
        return sman_;
      else
        return cman_;
    }

    /*!
    \brief Get meshtying manager

    */
    Teuchos::RCP<MORTAR::ManagerBase> MtManager() { return mtman_; }

    /*!
    \brief Get smoothing manager

    */
    Teuchos::RCP<MORTAR::ManagerBase> SManager() { return sman_; }

    /*!
    \brief Get strategy of mt/contact problem

    */
    MORTAR::StrategyBase& GetStrategy();

    /*!
    \brief return bool indicating if contact is defined

    */
    bool HaveContact()
    {
      if (sman_ != Teuchos::null)
        return true;
      else
        return (cman_ != Teuchos::null);
    }

    /*!
    \brief return bool indicating if meshtying is defined

    */
    bool HaveMeshtying() { return (mtman_ != Teuchos::null); }

    /*!
    \brief return bool indicating if smoothing is defined

    */
    bool HaveSmoothing() { return (sman_ != Teuchos::null); }

    /*!
    \brief postprocess tractions for mesht/contact

    */
    void PostprocessQuantities(Teuchos::RCP<IO::DiscretizationWriter>& output);

    /*!
    \brief read restart

    */
    void ReadRestart(IO::DiscretizationReader& reader, Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<Epetra_Vector> zero);
    /*!
    \brief recover lagr. mult. for contact/meshtying and slave displ for mesht.

    */
    void Recover(Teuchos::RCP<Epetra_Vector> disi);

    /*!
    \brief set state vector

    */
    void SetState(Teuchos::RCP<Epetra_Vector> zeros);

    /*!
    \brief store dirichlet status

    */
    void StoreDirichletStatus(Teuchos::RCP<LINALG::MapExtractor> dbcmaps);

    /*!
    \brief update

    */
    void Update(Teuchos::RCP<Epetra_Vector> dis);

    /*!
    \brief visualize stuff with gmsh

    */
    void VisualizeGmsh(const int istep, const int iter = -1);

    /*!
    \brief write restart

    */
    void WriteRestart(Teuchos::RCP<IO::DiscretizationWriter>& output, bool forcedrestart = false);
    //@}

   protected:
    // don't want cctor (= operator impossible anyway for abstract class)
    MeshtyingContactBridge(const MeshtyingContactBridge& old) = delete;

    // internal variables:

    Teuchos::RCP<MORTAR::ManagerBase> cman_;   //! contact manager
    Teuchos::RCP<MORTAR::ManagerBase> mtman_;  //! meshtying manager
    Teuchos::RCP<MORTAR::ManagerBase> sman_;   //! smoothing manager

  };  // class MeshtyingContactBridge
}  // namespace CONTACT

#endif /* MESHTYING_CONTACT_BRIDGE_H_ */

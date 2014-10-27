/*!----------------------------------------------------------------------
\file meshtying_contact_bridge.H

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*-----------------------------------------------------------------------*/
#ifndef MESHTYING_CONTACT_BRIDGE_H_
#define MESHTYING_CONTACT_BRIDGE_H_

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "../drt_inpar/inpar_contact.H"
#include "Epetra_Comm.h"
#include "../drt_mortar/mortar_manager_base.H" // TODO forward declaration!!!
//#include "Epetra_Map.h"


// forward declarations
//namespace IO
//{
// class DiscretizationWriter;
// //class DiscretizationReader;
//}

namespace DRT
{
 class Discretization;
 class Condition;
}

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
//class CoManager;
//class MtManager;

class MeshtyingContactBridge
{
public:

  //@}

  MeshtyingContactBridge(DRT::Discretization& dis,
                         std::vector<DRT::Condition*>& mtcond,
                         std::vector<DRT::Condition*>& ccond,
                         double alphaf);

  /*!
  \brief Destructor

  */
  virtual ~MeshtyingContactBridge() {};

  //! @name Access methods

  /*!
  \brief Get Epetra communicator

  */
  const Epetra_Comm& Comm() const { return cman_->Comm();}

  /*!
  \brief Get contact manager

  */
  Teuchos::RCP<MORTAR::ManagerBase> ContactManager() {return cman_;}

  /*!
  \brief Get meshtying manager

  */
  Teuchos::RCP<MORTAR::ManagerBase> MtManager() {return mtman_;}

  /*!
  \brief Get strategy of mt/contact problem

  */
  MORTAR::StrategyBase& GetStrategy();

  /*!
  \brief return bool indicating if contact is defined

  */
  bool HaveContact() { return (cman_ != Teuchos::null); }

  /*!
  \brief return bool indicating if meshtying is defined

  */
  bool HaveMeshtying() { return (mtman_ != Teuchos::null); }

  /*!
  \brief postprocess tractions for mesht/contact

  */
  void PostprocessTractions(Teuchos::RCP<IO::DiscretizationWriter>& output);

  /*!
  \brief read restart

  */
  void ReadRestart(IO::DiscretizationReader& reader,
                   Teuchos::RCP<Epetra_Vector> dis,
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
  void Update(int istep, Teuchos::RCP<Epetra_Vector> dis);

  /*!
  \brief write restart

  */
  void WriteRestart(Teuchos::RCP<IO::DiscretizationWriter>& output,
                    bool forcedrestart = false);
  //@}

protected:

  // don't want cctor (= operator impossible anyway for abstract class)
  MeshtyingContactBridge(const MeshtyingContactBridge& old);

  // internal variables:

  Teuchos::RCP<MORTAR::ManagerBase>   cman_;   //! contact manager
  Teuchos::RCP<MORTAR::ManagerBase>   mtman_;  //! meshtying manager

}; // class MeshtyingContactBridge
} // namespace CONTACT

#endif /* MESHTYING_CONTACT_BRIDGE_H_ */

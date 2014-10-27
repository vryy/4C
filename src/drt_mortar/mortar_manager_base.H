/*!----------------------------------------------------------------------
\file mortar_manager_base.H

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/
#ifndef MORTAR_MANAGER_BASE_H
#define MORTAR_MANAGER_BASE_H

#include "Teuchos_RCP.hpp"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"


// forward declarations
namespace IO
{
 class DiscretizationWriter;
 class DiscretizationReader;
}

namespace MORTAR
{

// forward declarations
class StrategyBase;

/*!
\brief Abstract base class to control all mortar coupling

\author popp (popp@lnm.mw.tum.de)
*/
class ManagerBase
{
public:

  //! @name Enums and Friends
  //@}

  /*!
  \brief Standard Constructor

  The base class constructor is empty.

  One needs a derived class for a concrete implementation of the Manager
  class into a given FE code environment (see e.g. contact_manager.H and
  contact_manager.cpp for the BACI mortar contact implementation or
  meshtying_manager.H and meshtying_manager.coo for the BACI mortar
  meshtying implementation).

  This constructor then has to be fed with a discretization that is expected
  to carry at least two mortar boundary conditions (one is only sufficient
  in the case of self contact simulations). It extracts all mortar boundary
  conditions, constructs one or multiple mortar interfaces and stores them.

  It also builds the corresponding strategy solver object and stores a
  reference in the strategy_ member variable.

  */
  ManagerBase();

  /*!
  \brief Destructor

  */
  virtual ~ManagerBase() {};

  //! @name Access methods

  /*!
  \brief Get Epetra communicator

  */
  const Epetra_Comm& Comm() const { return *comm_;}

  /*!
  \brief Return the object for the solving strategy.

  All necessary steps for the computation algorithm
  have to be specialized in subclasses of StrategyBase

  */
  MORTAR::StrategyBase& GetStrategy() { return *strategy_;}

  //@}

  //! @name Purely virtual functions

  // All these functions are defined in one or more specific derived classes,
  // such as CONTACT::CManager or CONTACT::MtManager. As the base class
  // MORTAR::ManagerBase is always called from the control routine (time
  // integrator), these functions need to be defined purely virtual here.

  virtual void PostprocessTractions(IO::DiscretizationWriter& output) = 0;
  virtual void ReadRestart(IO::DiscretizationReader& reader, Teuchos::RCP<Epetra_Vector> dis, Teuchos::RCP<Epetra_Vector> zero) = 0;
  virtual void WriteRestart(IO::DiscretizationWriter& output, bool forcedrestart=false) = 0;

  //@}

protected:

  // don't want cctor (= operator impossible anyway for abstract class)
  ManagerBase(const ManagerBase& old);

  Teuchos::RCP<Epetra_Comm>               comm_;     // communicator
  Teuchos::RCP<MORTAR::StrategyBase>      strategy_; // strategy object

}; // class ManagerBase
} // namespace MORTAR


#endif  // #ifndef MORTAR_MANAGER_BASE_H

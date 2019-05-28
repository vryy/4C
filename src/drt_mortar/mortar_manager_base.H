/*-----------------------------------------------------------------------*/
/*!
\maintainer Matthias Mayr

\level 2

\brief Base class for contact and meshtying managers (structural problems only)
       all other problem types use mortar adapters
*/
/*---------------------------------------------------------------------*/
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
}  // namespace IO

namespace DRT
{
  class Node;
  class Discretization;
  class Element;
}  // namespace DRT

namespace Teuchos
{
  class ParameterList;
}

namespace MORTAR
{
  // forward declarations
  class StrategyBase;
  class MortarNode;
  class MortarElement;

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
    virtual ~ManagerBase() = default;

    //! @name Access methods

    /*!
    \brief Get Epetra communicator

    */
    const Epetra_Comm& Comm() const { return *comm_; }

    /*!
    \brief Return the object for the solving strategy.

    All necessary steps for the computation algorithm
    have to be specialized in subclasses of StrategyBase

    */
    MORTAR::StrategyBase& GetStrategy() { return *strategy_; }

    //@}

    //! @name Purely virtual functions
    //! @{

    //! Write interface quantities for postprocessing
    virtual void PostprocessQuantities(IO::DiscretizationWriter& output) = 0;

    /*!
    \brief Write results for visualization separately for each meshtying/contact interface

    Call each interface, such that each interface can handle its own output of results.

    \param[in] outputParams Parameter list with stuff required by interfaces to write output
    */
    virtual void PostprocessQuantitiesPerInterface(
        Teuchos::RCP<Teuchos::ParameterList> outputParams) = 0;

    //! Read restart data from disk
    virtual void ReadRestart(IO::DiscretizationReader& reader, Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<Epetra_Vector> zero) = 0;

    //! Write restart data to disk
    virtual void WriteRestart(IO::DiscretizationWriter& output, bool forcedrestart = false) = 0;

    //! @}

   protected:
    // don't want cctor (= operator impossible anyway for abstract class)
    ManagerBase(const ManagerBase& old) = delete;

    //! Communicator
    Teuchos::RCP<Epetra_Comm> comm_;

    //! Strategy object
    Teuchos::RCP<MORTAR::StrategyBase> strategy_;

  };  // class ManagerBase
}  // namespace MORTAR


#endif  // #ifndef MORTAR_MANAGER_BASE_H

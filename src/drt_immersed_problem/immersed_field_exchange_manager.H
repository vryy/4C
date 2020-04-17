/*----------------------------------------------------------------------*/
/*! \file

\brief manage access to and provide data globally in immersed problems

\level 3

\maintainer Jonas Eichinger

*----------------------------------------------------------------------*/

#ifndef IMMERSED_FIELD_EXCHANGE_MANAGER_H
#define IMMERSED_FIELD_EXCHANGE_MANAGER_H

#include <map>
#include <set>
#include <cstdlib>
#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"

#include "../drt_lib/drt_singletondestruction.H"
#include "../linalg/linalg_fixedsizematrix.H"

namespace GEO
{
  class SearchTree;
}

namespace ADAPTER
{
  class FSIStructureWrapperImmersed;
}

namespace DRT
{
  class Discretization;

  class ImmersedFieldExchangeManager : public DRT::SingletonDestruction
  {
   public:
    /// virtual empty destructor (dtor)
    virtual ~ImmersedFieldExchangeManager() {}

    /// singleton access method
    static ImmersedFieldExchangeManager* Instance(bool create = true);

    /// called upon destruction
    virtual void Done();

    /// pointer to adapter FSIStructureWrapperImmersed
    Teuchos::RCP<ADAPTER::FSIStructureWrapperImmersed> FSIStructureWrapperImmersed_;

    void SetAdapter(Teuchos::RCP<ADAPTER::FSIStructureWrapperImmersed> adapter)
    {
      FSIStructureWrapperImmersed_ = adapter;
      return;
    };

    Teuchos::RCP<ADAPTER::FSIStructureWrapperImmersed> GetAdapter()
    {
      return FSIStructureWrapperImmersed_;
    };

    //@}

   protected:
    /// protected constructor since we are a Singleton.
    ImmersedFieldExchangeManager() = default;

    /// declaration of the static instance
    static DRT::ImmersedFieldExchangeManager* instance_;


  };  // class ImmersedFieldExchangeManager

}  // namespace DRT


#endif /*! IMMERSED_FIELD_EXCHANGE_MANAGER_H */

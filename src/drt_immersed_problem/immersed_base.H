/*!----------------------------------------------------------------------
\file immersed_base.cpp

\brief base class for all immersed algorithms

<pre>
Maintainers: Andreas Rauch & Anh-Tu Vuong
             {rauch,vuong}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289--15240 / 15264
</pre>
*----------------------------------------------------------------------*/
#ifndef IMMERSED_BASE_H
#define IMMERSED_BASE_H

//#include "../drt_lib/drt_discret.H" // check if really necessary or far too much
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "../drt_lib/drt_condition.H"

namespace DRT{
  class Discretization;
}

namespace IMMERSED{

class ImmersedBase
{
public:

  /// setup base class for all immersed algorithms
  explicit ImmersedBase();

  std::vector<int> DetermineImmersionDomain(Teuchos::RCP<DRT::Discretization> backgrounddis, Teuchos::RCP<DRT::Discretization> immerseddis, bool firstcall=false);

  std::vector<int> DetermineImmersionBoundaryDomain(Teuchos::RCP<DRT::Discretization> backgrounddis, Teuchos::RCP<DRT::Discretization> immerseddis,  const std::string& condname, bool firstcall=false);

  void CreateGhosting(Teuchos::RCP<DRT::Discretization> distobeghosted);

  void CreateVolumeCondition(Teuchos::RCP<DRT::Discretization> dis, std::vector<int> dvol_fenode, DRT::Condition::ConditionType condtype, std::string condname);

  void UpdateVolumeCondition(Teuchos::RCP<DRT::Discretization> dis, std::vector<int> dvol_fenode, DRT::Condition::ConditionType condtype, std::string condname);

  Teuchos::RCP<Epetra_Vector> InterpolateToImmersedNodes(Teuchos::RCP<DRT::Discretization> dis, std::vector<int> dvol_fenode, std::string condname);


};//class ImmersedBase
}// namsepace IMMERSED

#endif

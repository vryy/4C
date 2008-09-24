/*!----------------------------------------------------------------------
\file statmech_time.H

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15234
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STATMECH_H
#define STATMECH_H

#include "../drt_structure/strugenalpha.H"
#include "../drt_statmech/statmech_manager.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <random/normal.h>


// forward declarations

/*!
\brief time integration for structural problems with statistical mechanics

\author cyron (cyron@lnm.mw.tum.de)
*/
class StatMechTime : public virtual StruGenAlpha
{
public:

  //! \brief Standard Constructor
  explicit StatMechTime(ParameterList& params,
                   DRT::Discretization& dis,
                   LINALG::Solver& solver,
                   IO::DiscretizationWriter& output);

  //! \brief Destructor
  virtual ~StatMechTime() {};

  //! \brief Do time integration
  void Integrate();
    

protected:

  // don't want = operator and cctor
  StatMechTime operator = (const StatMechTime& old);
  StatMechTime(const StatMechTime& old);
  
  //! \brief manager for administration of statistical mechanics stuff
  RCP<StatMechManager>         statmechmanager_;


}; // class StatMechTime

#endif  // #ifndef STATMECH_H
#endif  // #ifdef CCADISCRET

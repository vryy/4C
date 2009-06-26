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
#ifndef STATMECH_TIME_H
#define STATMECH_TIME_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_structure/strugenalpha.H"
#include "../drt_lib/drt_globalproblem.H"
/*!
\brief time integration for structural problems with statistical mechanics

\author cyron (cyron@lnm.mw.tum.de)
*/

// forward declarations
class StatMechManager;


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
  
  //! \brief do consistent predictor step
  virtual void ConsistentPredictor();    

  //!\brief do full Newton iteration
  virtual void FullNewton();
  
  //!\brief solve non-linear problem by means of PTC method
  virtual void PTC();
  
  //!\brief writing output including restart information
  virtual void Output();
  
  //!\brief read restart information
  virtual void ReadRestart(int step);
    

protected:

  // don't want = operator and cctor
  StatMechTime operator = (const StatMechTime& old);
  StatMechTime(const StatMechTime& old);
  
  //! \brief manager for administration of statistical mechanics stuff
  RCP<StatMechManager>         statmechmanager_;
  
  //! \brief class variable to make statistical forces globally accessible in time integration
  RCP<Epetra_Vector>    browniancol_;


}; // class StatMechTime

#endif  // #ifndef STATMECH_H
#endif  // #ifdef CCADISCRET

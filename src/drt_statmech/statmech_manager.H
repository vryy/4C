/*!----------------------------------------------------------------------
\file statmech_manager.H

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15234
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STATMECH_MANAGER_H
#define STATMECH_MANAGER_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/drt_globalproblem.H"
#include <iostream>

// forward declarations

/*!
\brief management and auxiliary functions for statistical mechanics

\author cyron (cyron@lnm.mw.tum.de)
*/
class StatMechManager
{
public:

  //! \brief Standard Constructor
  explicit StatMechManager(ParameterList& params);

  //! \brief Destructor
  virtual ~StatMechManager() {};

  //! \brief write special output for statistical mechanics into files
  void StatMechOutput(const double& time,const int& num_dof,const int& istep, const double& dt, const Epetra_Vector& dis);

protected:

  // don't want = operator and cctor
  StatMechManager operator = (const StatMechManager& old);
  StatMechManager(const StatMechManager& old);
  
  Teuchos::ParameterList statmechparams_; // containing statistical mechanics input parameters
  double maxtime_; //maximal simulation time
  double starttimeoutput_; //point of time at which writing special statistical output is started
  double endtoendref_; //end to end distance of a single filament at the beginning of statistical output
  int istart_; //time step number at the beginning of statistical output
  std::ostringstream outputfilename_; //name of file into which output is written
  

}; // class StatMechManager

#endif /*STATMECH_MANAGER_H*/
#endif  // #ifdef CCADISCRET

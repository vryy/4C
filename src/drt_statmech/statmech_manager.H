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
#include "../drt_lib/linalg_utils.H"

// forward declarations

/*!
\brief management and auxiliary functions for statistical mechanics

\author cyron (cyron@lnm.mw.tum.de)
*/
class StatMechManager
{
public:

  //! \brief Standard Constructor
  explicit StatMechManager(ParameterList& params, DRT::Discretization& discret);

  //! \brief Destructor
  virtual ~StatMechManager() {};

  //! \brief write special output for statistical mechanics into files
  void StatMechOutput(const double& time,const int& num_dof,const int& istep, const double& dt, const Epetra_Vector& dis);
  
  //! \brief special update for statistical mechanics including setting and deleting crosslinkers in the discretization
  void StatMechUpdate();

protected:

  // don't want = operator and cctor
  StatMechManager operator = (const StatMechManager& old);
  StatMechManager(const StatMechManager& old);
  
  //! \brief containing statistical mechanics input parameters
  Teuchos::ParameterList statmechparams_; 
  
  //! \brief maximal simulation time
  double maxtime_; 
  
  //! \brief point of time at which writing special statistical output is started
  double starttimeoutput_; 
  
  //! \brief end to end distance of a single filament at the beginning of statistical output
  double endtoendref_; 
  
  //! \brief time step number at the beginning of statistical output
  int istart_; 
  
  //! \brief maximal radius within which crosslinkers are established
  double rlink_;
  
  //! \brief vector storing the currently unused global element Ids (due to former deletion of crosslinker elements)
  vector<int> unusedids_;
  
  //! \brief reference to current discretization
  DRT::Discretization& discret_; 
  
  //! \brief name of file into which output is written
  std::ostringstream outputfilename_; 
  
  //! \brief storing for each node global Id of node with which it is crosslinked; -1 means: no crosslink
  RCP<Epetra_Vector>  crosslinkerpartner_;
  
  //! \brief storing for each node global Id of crosslinking Element; -1 means: no crosslink
  RCP<Epetra_Vector> crosslinkerelements_;
  

}; // class StatMechManager

#endif /*STATMECH_MANAGER_H*/
#endif  // #ifdef CCADISCRET

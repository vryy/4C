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


/*!
\brief management and auxiliary functions for statistical mechanics
\author cyron (cyron@lnm.mw.tum.de)
*/

class StatMechManager
{
public:

  //! \brief Standard Constructor
  explicit StatMechManager(ParameterList& params, DRT::Discretization& discret, RCP<LINALG::SparseMatrix>& stiff, RCP<LINALG::SparseMatrix>& damp);

  //! \brief Destructor
  virtual ~StatMechManager() {};

  //! \brief write special output for statistical mechanics into files
  void StatMechOutput(const double& time,const int& num_dof,const int& istep, const double& dt, const Epetra_Vector& dis);
  
  //! \brief special update for statistical mechanics including setting and deleting crosslinkers in the discretization
  void StatMechUpdate(const double dt, const Epetra_Vector& dis);
  
  //! \brief updates system damping matrix and external force vector according to influence of thermal bath
  void StatMechForceDamp(ParameterList& params, RCP<Epetra_Vector> dis, RCP<Epetra_Vector> fext, RCP<LINALG::SparseMatrix> damp);
  
  //! \brief containing statistical mechanics input parameters
  Teuchos::ParameterList statmechparams_;

protected:

  // don't want = operator and cctor
  StatMechManager operator = (const StatMechManager& old);
  StatMechManager(const StatMechManager& old);
   
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
  
  //! \brief number of elements which are part of the discretization already before starting adding crosslinkers
  int basiselements_;
  
  //! \brief reference to current discretization
  DRT::Discretization& discret_; 
  
  //! \brief reference to current stiffness matrix
  RCP<LINALG::SparseMatrix>& stiff_;
  
  //! \brief reference to current damping matrix
  RCP<LINALG::SparseMatrix>& damp_;
  
  //! \brief name of file into which output is written
  std::ostringstream outputfilename_; 
  
  //! \brief storing for each node global Id of node with which it is crosslinked; -1 means: no crosslink
  RCP<Epetra_Vector>  crosslinkerpartner_;
  
  //! \brief storing for each node number of filament to which it belongs
  RCP<Epetra_Vector> filamentnumber_;

  

}; // class StatMechManager

#endif /*STATMECH_MANAGER_H*/
#endif  // #ifdef CCADISCRET

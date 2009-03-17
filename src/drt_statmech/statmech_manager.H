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
#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/intersection_service.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"


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
  void StatMechOutput(ParameterList& params, const int ndim, const double& time,const int& istep, const double& dt, const Epetra_Vector& dis, const Epetra_Vector& fint);
   
  //! \brief saving all necessary class variables for restart
  void StatMechWriteRestart(IO::DiscretizationWriter& output);
  
  //! \brief reading all necessary information after a restart
  void StatMechReadRestart(IO::DiscretizationReader& reader);
    
  //! \brief special update for statistical mechanics including setting and deleting crosslinkers in the discretization
  void StatMechUpdate(const double dt, const Epetra_Vector& disrow);

  //! \brief updates system damping matrix and external force vector according to influence of thermal bath
  void StatMechBrownian(ParameterList& params, RCP<Epetra_Vector> dis, RCP<Epetra_Vector> fext, RCP<LINALG::SparseMatrix> damp);
  
  //! \brief set crosslinkers according to setcrosslinkercol (if entry for a node equals 1 a crosslinker may be set for this node)
  void SetCrosslinkers(const Epetra_Vector& setcrosslinkercol,const Epetra_Map& nodecolmap,const Epetra_Vector& nearestneighbour,const std::map<int,LINALG::Matrix<3,1> >& currentpositions,const std::map<int,LINALG::Matrix<3,1> >& currentrotations);
  
  //! \brief initializing special output for statistical mechanics by looking for a suitable name of the outputfile and setting up an empty file with this name
  void StatMechInitOutput(const int ndim);
  
  //! \brief write output for Gmsh postprocessing
  void GmshOutput(const Epetra_Vector& disrow, const std::ostringstream& filename, const int& step);
  
  //! \brief containing statistical mechanics input parameters
  Teuchos::ParameterList statmechparams_;

protected:

  // don't want = operator and cctor
  StatMechManager operator = (const StatMechManager& old);
  StatMechManager(const StatMechManager& old); 
  
  //! \brief delete crosslinkers according to delcrosslinkercol (if entry for a node equals 1 its crosslinker is deleted)
  void DelCrosslinkers(const Epetra_Vector& delcrosslinkercol,const Epetra_Map& nodecolmap);
   
  //! \brief maximal simulation time
  double maxtime_; 
  
  //! \brief point of time at which writing special statistical output is started
  double starttimeoutput_; 
  
  //! \brief end to end distance of a single filament at the beginning of statistical output
  double endtoendref_; 
  
  //! \brief time step number at the beginning of statistical output
  int istart_; 
  
  //! \brief position of first node in last time step
  LINALG::Matrix<3,1> beginold_;
  
  //! \brief position of last nodein last time step
  LINALG::Matrix<3,1> endold_;
  
  //! \brief maximal radius within which crosslinkers are established
  double rlink_;
  
  //! \brief number of elements which are part of the discretization already before starting adding crosslinkers
  int basiselements_;
  
  //! \brief number of file into which output is written (only important for single filament dynamics)
  int outputfilenumber_;
  
  //! \brief reference to current discretization
  DRT::Discretization& discret_; 
  
  //! \brief reference to current stiffness matrix
  RCP<LINALG::SparseMatrix>& stiff_;
  
  //! \brief reference to current damping matrix
  RCP<LINALG::SparseMatrix>& damp_;
  
  //! \brief storing for each node global Id of node with which it is crosslinked; -1 means: no crosslink
  RCP<Epetra_Vector>  crosslinkerpartner_;
  
  //! \brief storing for each node number of filament to which it belongs
  RCP<Epetra_Vector> filamentnumber_;
  
  //! \brief storing for degree of freedom whether force has to be measured
  RCP<Epetra_Vector> forcesensor_;
  
  //! SearchTree which allows accelerated search for neighbouring nodes in case of dynamic crosslinkers
  Teuchos::RCP<GEO::SearchTree>  octTree_;

  

}; // class StatMechManager

#endif /*STATMECH_MANAGER_H*/
#endif  // #ifdef CCADISCRET

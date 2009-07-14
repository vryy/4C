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
  explicit StatMechManager(ParameterList& params, DRT::Discretization& discret, RCP<LINALG::SparseMatrix>& stiff);

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
  void StatMechBrownian(ParameterList& params, RCP<Epetra_Vector> dis, RCP<Epetra_Vector> browniancol);
  
  //! \brief initializing special output for statistical mechanics by looking for a suitable name of the outputfile and setting up an empty file with this name
  void StatMechInitOutput(const int ndim, const double& dt);
  
  //! \brief write output for Gmsh postprocessing
  void GmshOutput(const Epetra_Vector& disrow, const std::ostringstream& filename, const int& step);
  
  //! \brief set crosslinker if two nodes are close enough and probability check is passed
  void SetCrosslinkers(const double& dt, const Epetra_Map& noderowmap, const Epetra_Map& nodecolmap,const std::map<int,LINALG::Matrix<3,1> >& currentpositions,const std::map<int,LINALG::Matrix<3,1> >& currentrotations);
  
  //! \brief containing statistical mechanics input parameters
  Teuchos::ParameterList statmechparams_;

protected:

  // don't want = operator and cctor
  StatMechManager operator = (const StatMechManager& old);
  StatMechManager(const StatMechManager& old); 
   
  //! \brief delete crosslinkers according to delcrosslinkercol (if entry for a node equals 1 its crosslinker is deleted)
  void DelCrosslinkers(const double& dt, const Epetra_Map& noderowmap);
  
  //! \brief Searches and saves in variable  crosslinkerneighbours_ neighbours for each row node
  void SearchNeighbours(const std::map<int,LINALG::Matrix<3,1> > currentpositions);
   
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
  int basisnodes_;
  
  //! \brief number of file into which output is written (only important for single filament dynamics)
  int outputfilenumber_;
  
  //! \brief reference to current discretization
  DRT::Discretization& discret_; 
  
  //! \brief reference to current stiffness matrix
  RCP<LINALG::SparseMatrix>& stiff_;
  
  //! \brief storing for each row map node the column map LID of all proper neighbour nodes
  std::vector< vector<int> >  crosslinkerneighbours_;
  
  //! \brief storing for each node GIDs of nodes to which a crosslinker has been established
  std::vector< vector<int> >  crosslinkerpartner_;
  
  //! \brief storing for each node LID number of filament to which it belongs
  RCP<Epetra_Vector> filamentnumber_;
  
  //! \brief storing for degree of freedom whether force has to be measured
  RCP<Epetra_Vector> forcesensor_;
  
  //! SearchTree which allows accelerated search for neighbouring nodes in case of dynamic crosslinkers
  Teuchos::RCP<GEO::SearchTree>  octTree_;

  
}; // class StatMechManager

#endif /*STATMECH_MANAGER_H*/
#endif  // #ifdef CCADISCRET

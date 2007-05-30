/*!----------------------------------------------------------------------
\file drt_periodicbc.cpp

\brief Methods to apply periodic boundary conditions


<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef PERIODICBC_H
#define PERIODICBC_H

#include "Teuchos_RefCountPtr.hpp"
#include "drt_discret.H"
#include "drt_dofset.H"
#include "drt_pbcdofset.H"
#include "drt_nodematchingoctree.H"
#include <Teuchos_TimeMonitor.hpp>

using namespace Teuchos;

/*!
\class DRT::PeriodicBoundaryConditions

\brief

\note 
\author gammi
*/
class PeriodicBoundaryConditions
{
public:

  /*!
  \brief Standard Constructor

  \param  dis (i) Discretisation which contains the nodes of the
                  periodc boundary condition

                  
  \return void
  
  */
  PeriodicBoundaryConditions(RefCountPtr<DRT::Discretization> dis);
  
  /*!
  \brief Destructor

  */
  virtual ~PeriodicBoundaryConditions();


  //! @name Framework

  /*!
  \brief Proceeds all pairs of periodic boundary conditions
  
  \return void

  */
  void UpdateDofsForPeriodicBoundaryConditions();
  //@}


  //! @name Methods to process one pair of pbcs

  
  /*!
  \brief Couple nodes for specific pair of periodic boundary conditions

  \param  midtosid        (o) map from master to slavenodes
  \param  masternodeids   (i) all master node ids
  \param  slavenodeids    (i) all slave node ids
  \param  dofsforpbcplane (i) periodic boundary is in this plane

  \return void
  
  */
  void CreateNodeCouplingForSinglePBC(
    map<int,int>       &midtosid,
    const vector <int>  masternodeids,
    const vector <int>  slavenodeids,
    const vector <int>  dofsforpbcplane);

  /*!
  \brief Add the connectivity from this condition to the connectivity
  of all previously processed periodic boundary conditions.
  Redistribute the nodes and assign the dofs to the
  current distribution of nodes

  \param  midtosid        (i) map from master to slavenodes
  \param  pbcid           (i) the id of the periodic boundary condition
  
  \return void

  */
  void AddConnectivityRedistributeAndCreateDofCoupling(
    map<int,int> &midtosid,const int pbcid);
  
  //@}

  
protected:
  //!\brief the discretisation
  RefCountPtr<DRT::Discretization> discret_;
  
  //!\brief number of pairs of periodic boundary conditions
  int                     numpbcpairs_;

  //!\brief vector of periodic surface boundary conditions
  vector<DRT::Condition*> mysurfpbcs_;

  //!\brief map connecting master to slave nodes owned by this proc
  //       master node -> list of his slave node(s)
  RefCountPtr<map<int,vector<int> > > allcoupledrownodes_;

  //!\brief map connecting master to slave nodes owned or ghosted by this proc
  //       master node -> list of his slave node(s)
  RefCountPtr<map<int,vector<int> > > allcoupledcolnodes_;



  //!\brief time measurement (total)
  RefCountPtr<Time> timepbctot_     ;
  //!\brief time measurement (create master slave matching for pairs)
  RefCountPtr<Time> timepbcmidtosid_;
  //!\brief time measurement (create octree)
  RefCountPtr<Time> timepbcmidoct_  ;
  //!\brief time measurement (search in octree)
  RefCountPtr<Time> timepbcmidmatch_;
  //!\brief time measurement (redistribute nodes)
  RefCountPtr<Time> timepbcreddis_  ;
  //!\brief time measurement (add connectivity to previous pbcs)
  RefCountPtr<Time> timepbcaddcon_  ;
  //!\brief time measurement (fetch slave nodes to masterproc)
  RefCountPtr<Time> timepbcfetchs_  ;
  //!\brief time measurement (repair ghosting)
  RefCountPtr<Time> timepbcghost_   ;
  //!\brief time measurement (make colmap for ghosting) 
  RefCountPtr<Time> timepbcmakeghostmap_   ;
  //!\brief time measurement (discret->Redistribute)
  RefCountPtr<Time> timepbcrenumdofs_   ;

  
  RefCountPtr<TimeMonitor> tm0_ref_ ;
  RefCountPtr<TimeMonitor> tm1_ref_ ;
  RefCountPtr<TimeMonitor> tm2_ref_ ;
  RefCountPtr<TimeMonitor> tm3_ref_ ;
  RefCountPtr<TimeMonitor> tm4_ref_ ;
  RefCountPtr<TimeMonitor> tm5_ref_ ;
  RefCountPtr<TimeMonitor> tm6_ref_ ;
  RefCountPtr<TimeMonitor> tm7_ref_ ;
  RefCountPtr<TimeMonitor> tm8_ref_ ;
  RefCountPtr<TimeMonitor> tm9_ref_ ;


  
};

#endif  // #ifndef PERIODICBC_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET

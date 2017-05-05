/*----------------------------------------------------------------------*/
/*!
\file xfem_neumann.cpp

\brief base xfem Neumann boundary conditions

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning think about removing these routines!!!

*/
/*----------------------------------------------------------------------*/


#include "../drt_xfem/xfem_neumann.H"

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fluid_ele/fluid_ele.H"
#include "../linalg/linalg_utils.H"



/*----------------------------------------------------------------------*
 |  evaluate Neumann conditions (public)                    schott 08/11|
 *----------------------------------------------------------------------*/
void XFEM::EvaluateNeumann(                    Teuchos::ParameterList&              params,
                                               Teuchos::RCP<DRT::Discretization>    discret,
                                               Teuchos::RCP<Epetra_Vector>          systemvector,
                                               Teuchos::RCP<LINALG::SparseOperator> systemmatrix)
{
  if (systemmatrix==Teuchos::null)
    EvaluateNeumann(params,discret,*systemvector);
  else
    EvaluateNeumann(params,discret,*systemvector,systemmatrix.get());
  return;
}


/*----------------------------------------------------------------------*
 |  evaluate Neumann conditions (public)                    schott 08/11|
 *----------------------------------------------------------------------*/
void XFEM::EvaluateNeumann(  Teuchos::ParameterList&              params,
                             Teuchos::RCP<DRT::Discretization>    discret,
                             Epetra_Vector&                       systemvector,
                             LINALG::SparseOperator*              systemmatrix)
{

  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 5) EvaluateNeumann" );


  if (!discret->Filled()) dserror("FillComplete() was not called");
  if (!discret->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  bool assemblemat = (systemmatrix != NULL);

  // get the current time
  const double time = params.get("total time",-1.0);

  std::multimap<std::string,DRT::Condition* >::iterator fool;
  std::multimap<std::string,DRT::Condition* > condition;

  // vector for conditions of one special type
  std::vector<DRT::Condition *> condition_vec;

  //================================================
  // Load Neumann conditions from discretization
  // REMARK: ° standard volume Neumann conditions are not loaded -> evaluated in Evaluate
  //         ° for XFEM Neumann boundaries: we assumme only XFEM Surface(!) Neumann conditions
  //================================================

  // get standard Point Neumann conditions
  condition_vec.clear();
  discret->GetCondition("PointNeumann", condition_vec);
  // copy conditions to a condition multimap
  for(size_t i=0; i< condition_vec.size(); i++)
  {
    condition.insert( std::pair<std::string,DRT::Condition* >(std::string("PointNeumann"),condition_vec[i]));
  }

  // get standard Surface Neumann conditions
  condition_vec.clear();
  discret->GetCondition("LineNeumann", condition_vec);
  for(size_t i=0; i< condition_vec.size(); i++)
  {
    condition.insert( std::pair<std::string,DRT::Condition* >(std::string("LineNeumann"),condition_vec[i]));
  }

  // get standard Surface Neumann conditions
  condition_vec.clear();
  discret->GetCondition("SurfaceNeumann", condition_vec);
  for(size_t i=0; i< condition_vec.size(); i++)
  {
    condition.insert( std::pair<std::string,DRT::Condition* >(std::string("SurfaceNeumann"),condition_vec[i]));
  }

  // NOTE: WE SKIP VolumeNeumann conditions -> there are evaluated at the fluid element level (Bodyforce!)

  //TODO: we should introduce safety checks, when the Neumann boundary is cut!
  //TODO: shift the functionality to XFEM Neumann conditions and automatically detect, whether the element is cut by an interface
  // similar to Weak Dirichlet conditions!

  // evaluate standard Neumann conditions
  EvaluateNeumannStandard(condition,
      time,
      assemblemat,
      params,
      discret,
      systemvector,
      systemmatrix);

  return;

}



/*----------------------------------------------------------------------*
 |  evaluate Neumann for standard conditions (public)       schott 08/11|
 *----------------------------------------------------------------------*/
void XFEM::EvaluateNeumannStandard( std::multimap<std::string,DRT::Condition* > &   condition,
                                    const double                          time,
                                    bool                                  assemblemat,
                                    Teuchos::ParameterList&               params,
                                    Teuchos::RCP<DRT::Discretization>     discret,
                                    Epetra_Vector&                        systemvector,
                                    LINALG::SparseOperator*               systemmatrix)
{
  //TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::EvaluateNeumannStandard" );

  std::multimap<std::string,DRT::Condition* >::iterator fool;

  //--------------------------------------------------------
  // loop through Point Neumann conditions and evaluate them
  //--------------------------------------------------------
  for (fool=condition.begin(); fool!=condition.end(); ++fool)
  {
    if (fool->first != (std::string)"PointNeumann") continue;
    if (assemblemat && !systemvector.Comm().MyPID())
      std::cout << "WARNING: No linearization of PointNeumann conditions" << std::endl;
    DRT::Condition& cond = *(fool->second);
    const std::vector<int>* nodeids = cond.Nodes();
    if (!nodeids) dserror("PointNeumann condition does not have nodal cloud");
    const int nnode = (*nodeids).size();
    const std::vector<int>*    funct  = cond.Get<std::vector<int> >("funct");
    const std::vector<int>*    onoff  = cond.Get<std::vector<int> >("onoff");
    const std::vector<double>* val    = cond.Get<std::vector<double> >("val");
    // Neumann BCs for some historic reason only have one curve
    int functnum = -1;
    if (funct) functnum = (*funct)[0];
    double functfac = 1.0;
    if (functnum>=0)
      functfac = DRT::Problem::Instance()->Funct(functnum).EvaluateTime(time);
    for (int i=0; i<nnode; ++i)
    {
      // do only nodes in my row map
      if (!discret->NodeRowMap()->MyGID((*nodeids)[i])) continue;
      DRT::Node* actnode = discret->gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node %d",(*nodeids)[i]);
      // call explicitly the main dofset, i.e. the first column
      std::vector<int> dofs = discret->Dof(0,actnode);
      const unsigned numdf = dofs.size();
      for (unsigned j=0; j<numdf; ++j)
      {
        if ((*onoff)[j]==0) continue;
        const int gid = dofs[j];
        double value  = (*val)[j];
        value *= functfac;
        const int lid = systemvector.Map().LID(gid);
        if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
        systemvector[lid] += value;
      }
    }
  }

  //--------------------------------------------------------
  // loop through line/surface Neumann BCs and evaluate them
  // ATTENTION: VolumeNeumann conditions (bodyforces) are evaluated in Evaluate
  //--------------------------------------------------------
  for (fool=condition.begin(); fool!=condition.end(); ++fool)
    if (fool->first == (std::string)"LineNeumann" ||
        fool->first == (std::string)"SurfaceNeumann"
    )
    {
      DRT::Condition& cond = *(fool->second);
      std::map<int,Teuchos::RCP<DRT::Element> >& geom = cond.Geometry();
      std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
      Epetra_SerialDenseVector elevector;
      Epetra_SerialDenseMatrix elematrix;
      for (curr=geom.begin(); curr!=geom.end(); ++curr)
      {
        // get element location vector, dirichlet flags and ownerships
        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        curr->second->LocationVector(*discret,lm,lmowner,lmstride);
        elevector.Size((int)lm.size());
        if (!assemblemat)
        {
          curr->second->EvaluateNeumann(params,*discret,cond,lm,elevector);
          LINALG::Assemble(systemvector,elevector,lm,lmowner);
        }
        else
        {
          const int size = (int)lm.size();
          if (elematrix.M() != size) elematrix.Shape(size,size);
          else memset(elematrix.A(),0,size*size*sizeof(double));
          curr->second->EvaluateNeumann(params,*discret,cond,lm,elevector,&elematrix);
          LINALG::Assemble(systemvector,elevector,lm,lmowner);
          systemmatrix->Assemble(curr->second->Id(),lmstride,elematrix,lm,lmowner);
        }
      }
    }

  return;
}

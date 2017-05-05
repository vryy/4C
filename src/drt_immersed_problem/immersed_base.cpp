/*!----------------------------------------------------------------------
\file immersed_base.cpp

\brief base class for all immersed algorithms

\level 2

\maintainer Andreas Rauch
            rauch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 -15240

*----------------------------------------------------------------------*/
#include "immersed_base.H"
#include "../linalg/linalg_utils.H"
#include "../drt_adapter/ad_fld_wrapper.H"
#include "../drt_adapter/ad_str_fsiwrapper_immersed.H"


IMMERSED::ImmersedBase::ImmersedBase():
issetup_(false),
isinit_(false)
{
  // empty
  return;
} // ImmersedBase


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::CreateVolumeCondition(const Teuchos::RCP<DRT::Discretization>& dis,
                                                   const std::vector<int> dvol_fenode,
                                                   const DRT::Condition::ConditionType condtype,
                                                   const std::string condname)
{
  // determine id of condition
  std::multimap<std::string,Teuchos::RCP<DRT::Condition> > allconditions;
  allconditions = dis->GetAllConditions();
  int id = (int)allconditions.size();
  id += 1;

  // build condition
  bool buildgeometry = false;
  Teuchos::RCP<DRT::Condition> condition =
          Teuchos::rcp(new DRT::Condition(id,condtype,buildgeometry,DRT::Condition::Volume));

  // add nodes to conditions
   condition->Add("Node Ids",dvol_fenode);

   // add condition to discretization
   dis->SetCondition(condname,condition);

   // fill complete if necessary
   if (!dis->Filled())
     dis -> FillComplete(false,false,false);

   return;
} // CreateVolumeCondition


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::BuildConditionDofMap(
  const Teuchos::RCP<const DRT::Discretization>&  dis,
  const std::string                         condname,
  const Teuchos::RCP<const Epetra_Map>&     cond_dofmap_orig,
  const int                                 numdof,
  Teuchos::RCP<Epetra_Map> &                cond_dofmap)
{
  // declare dof vector
  std::vector<int> mydirichdofs(0);

  // get condition and conditioned nodes
  DRT::Condition* condition = dis->GetCondition(condname);
  const std::vector<int>* cond_nodes = condition->Nodes();
  int cond_nodes_size = cond_nodes->size();

  if(cond_nodes_size==0)
    dserror("No nodes in nodal cloud of condition %s",condname.c_str());

  // loop over all conditioned nodes
  for(int node=0;node<cond_nodes_size;node++)
  {
    // get node id
    int nodeid = cond_nodes->at(node);
    // get node pointer
    DRT::Node* node_ptr = dis->gNode(nodeid);
    if(node_ptr==NULL)
      dserror("Could not get node with id %d",nodeid);

    // get dofs
    std::vector<int> dofs = dis->Dof(0,node_ptr);

    for (int dim=0;dim<numdof;++dim)
    {
      // if not already in original dirich map
      if(cond_dofmap_orig->LID(dofs[dim]) == -1)
        mydirichdofs.push_back(dofs[dim]);
    }

  } // loop over all conditioned nodes

  int nummydirichvals = mydirichdofs.size();
  cond_dofmap = Teuchos::rcp( new Epetra_Map(-1,nummydirichvals,&(mydirichdofs[0]),0,dis->Comm()) );

  return;
} // BuildConditionDofRowMap


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::DoDirichletCond(
    const Teuchos::RCP<Epetra_Vector>&        statevector,
    const Teuchos::RCP<const Epetra_Vector>&  dirichvals,
    const Teuchos::RCP<const Epetra_Map>&     dbcmap)
{
  int mynumvals = dbcmap->NumMyElements();
  double* myvals = dirichvals->Values();

  for(int i=0;i<mynumvals;++i)
  {
    int gid = dbcmap->GID(i);

#ifdef DEBUG
    if(mynumvals==0)
      dserror("dbcmap empty!");
    int err = -2;
    int lid = dirichvals->Map().LID(gid);
    err = statevector -> ReplaceGlobalValue(gid,0,myvals[lid]);
    if(err==-1)
      dserror("VectorIndex >= NumVectors()");
    else if (err==1)
        dserror("GlobalRow not associated with calling processor");
    else if (err != -1 and err != 1 and err != 0)
      dserror("Trouble using ReplaceGlobalValue on fluid state vector. ErrorCode = %d",err);
#else
    int lid = dirichvals->Map().LID(gid);
    statevector -> ReplaceGlobalValue(gid,0,myvals[lid]);
#endif

  }
  return;
} // DoDirichletCond


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::ApplyDirichlet(
    const Teuchos::RCP< ::ADAPTER::StructureWrapper>& field_wrapper,
    const Teuchos::RCP<DRT::Discretization>&  dis,
    const std::string                         condname,
    Teuchos::RCP<Epetra_Map>&                 cond_dofrowmap,
    const int                                 numdof,
    const Teuchos::RCP<const Epetra_Vector>&  dirichvals)
{
  // build map of dofs subjected to Dirichlet condition
  BuildConditionDofMap(
      dis,
      condname,
      field_wrapper->GetDBCMapExtractor()->CondMap(),
      numdof,
      cond_dofrowmap);

  // add adhesion dofs to dbc map
  field_wrapper->AddDirichDofs(cond_dofrowmap);

  // write Dirichlet values to systemvector
  DoDirichletCond(
      field_wrapper->WriteAccessDispnp(),
      dirichvals,
      field_wrapper->GetDBCMapExtractor()->CondMap());

  return;
} // ApplyDirichlet


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::ApplyDirichletToFluid(
    const Teuchos::RCP< ::ADAPTER::FluidWrapper>& field_wrapper,
    const Teuchos::RCP<DRT::Discretization>&  dis,
    const std::string                         condname,
    Teuchos::RCP<Epetra_Map>&                 cond_dofrowmap,
    const int                                 numdof,
    const Teuchos::RCP<const Epetra_Vector>&  dirichvals)
{
  // build map of dofs subjected to Dirichlet condition
  BuildConditionDofMap(
      dis,
      condname,
      field_wrapper->GetDBCMapExtractor()->CondMap(),
      numdof,
      cond_dofrowmap);

  // add adhesion dofs to dbc map
  field_wrapper->AddDirichCond(cond_dofrowmap);

  // write Dirichlet values to systemvector
  DoDirichletCond(
      field_wrapper->WriteAccessVelnp(),
      dirichvals,
      field_wrapper->GetDBCMapExtractor()->CondMap());

  return;
} // ApplyDirichlet


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::RemoveDirichlet(
    const Teuchos::RCP<const Epetra_Map>& cond_dofmap,
    const Teuchos::RCP< ::ADAPTER::StructureWrapper>& field_wrapper)
{
  field_wrapper->RemoveDirichDofs(cond_dofmap);
  return;
} // RemoveDirichlet


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::RemoveDirichletFromFluid(
    const Teuchos::RCP<const Epetra_Map>& cond_dofmap,
    const Teuchos::RCP< ::ADAPTER::FluidWrapper>& field_wrapper)
{
  field_wrapper->RemoveDirichCond(cond_dofmap);
  return;
} // RemoveDirichlet


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::EvaluateImmersed(Teuchos::ParameterList& params,
                                              Teuchos::RCP<DRT::Discretization> dis,
                                              DRT::AssembleStrategy* strategy,
                                              std::map<int,std::set<int> >* elementstoeval,
                                              Teuchos::RCP<GEO::SearchTree> structsearchtree,
                                              std::map<int,LINALG::Matrix<3,1> >* currpositions_struct,
                                              int action,
                                              bool evaluateonlyboundary)
{
  // pointer to element
  DRT::Element* ele;

  for(std::map<int, std::set<int> >::const_iterator closele = elementstoeval->begin(); closele != elementstoeval->end(); closele++)
  {
    for(std::set<int>::const_iterator eleIter = (closele->second).begin(); eleIter != (closele->second).end(); eleIter++)
    {
      ele=dis->gElement(*eleIter);

      DRT::ELEMENTS::FluidImmersedBase* immersedelebase = dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(ele);
      if(immersedelebase==NULL)
        dserror("dynamic cast from DRT::Element* to DRT::ELEMENTS::FluidImmersedBase* failed");

      // evaluate this element and fill vector with immersed dirichlets
      int row = strategy->FirstDofSet();
      int col = strategy->SecondDofSet();

      params.set<int>("action",action);
      params.set<Teuchos::RCP<GEO::SearchTree> >("structsearchtree_rcp",structsearchtree);
      params.set<std::map<int,LINALG::Matrix<3,1> >* >("currpositions_struct",currpositions_struct);
      params.set<int>("Physical Type",INPAR::FLUID::poro_p1);

      DRT::Element::LocationArray la(1);
      immersedelebase->LocationVector(*dis,la,false);
      strategy->ClearElementStorage( la[row].Size(), la[col].Size() );

      if(!evaluateonlyboundary)
        immersedelebase->Evaluate(params,*dis,la[0].lm_,
            strategy->Elematrix1(),
            strategy->Elematrix2(),
            strategy->Elevector1(),
            strategy->Elevector2(),
            strategy->Elevector3());
      else
      {
        if(immersedelebase->IsBoundaryImmersed())
          immersedelebase->Evaluate(params,*dis,la[0].lm_,
              strategy->Elematrix1(),
              strategy->Elematrix2(),
              strategy->Elevector1(),
              strategy->Elevector2(),
              strategy->Elevector3());
      }

      strategy->AssembleVector1( la[row].lm_, la[row].lmowner_ );
    }
  }
  return;
} // EvaluateImmersed

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::EvaluateImmersedNoAssembly(Teuchos::ParameterList& params,
                                              Teuchos::RCP<DRT::Discretization> dis,
                                              std::map<int,std::set<int> >* elementstoeval,
                                              Teuchos::RCP<GEO::SearchTree> structsearchtree,
                                              std::map<int,LINALG::Matrix<3,1> >* currpositions_struct,
                                              int action
                                                        )
{
  // pointer to element
  DRT::Element* ele;

  for(std::map<int, std::set<int> >::const_iterator closele = elementstoeval->begin(); closele != elementstoeval->end(); closele++)
  {
    for(std::set<int>::const_iterator eleIter = (closele->second).begin(); eleIter != (closele->second).end(); eleIter++)
    {
      ele=dis->gElement(*eleIter);

      DRT::ELEMENTS::FluidImmersedBase* immersedelebase = dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(ele);
      if(immersedelebase==NULL)
        dserror("dynamic cast from DRT::Element* to DRT::ELEMENTS::FluidImmersedBase* failed");

      // provide important objects to ParameterList
      params.set<int>("action",action);
      params.set<Teuchos::RCP<GEO::SearchTree> >("structsearchtree_rcp",structsearchtree);
      params.set<std::map<int,LINALG::Matrix<3,1> >* >("currpositions_struct",currpositions_struct);
      params.set<int>("Physical Type",INPAR::FLUID::poro_p1);
      if(dis->Name()=="fluid")
        params.set<std::string>("immerseddisname","structure");
      else if (dis->Name()=="porofluid")
        params.set<std::string>("immerseddisname","cell");
      else
        dserror("no corresponding immerseddisname set for this type of backgrounddis!");

      // evaluate the element
      Epetra_SerialDenseMatrix dummymat;
      Epetra_SerialDenseVector dummyvec;

      DRT::Element::LocationArray la(1);
      immersedelebase->LocationVector(*dis,la,false);

      immersedelebase->Evaluate(params,*dis,la[0].lm_,dummymat,dummymat,dummyvec,dummyvec,dummyvec);
    }
  }
  return;
} // EvaluateImmersedNoAssembly

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::EvaluateScaTraWithInternalCommunication(Teuchos::RCP<DRT::Discretization> dis,
                                                                     const Teuchos::RCP<const DRT::Discretization> idis,
                                                                     DRT::AssembleStrategy* strategy,
                                                                     std::map<int,std::set<int> >* elementstoeval,
                                                                     Teuchos::RCP<GEO::SearchTree> structsearchtree,
                                                                     std::map<int,LINALG::Matrix<3,1> >* currpositions_struct,
                                                                     Teuchos::ParameterList& params,
                                                                     bool evaluateonlyboundary)
{
  // pointer to element
  DRT::Element* ele;
  DRT::Element* iele;

  for(std::map<int, std::set<int> >::const_iterator closele = elementstoeval->begin(); closele != elementstoeval->end(); closele++)
  {
    for(std::set<int>::const_iterator eleIter = (closele->second).begin(); eleIter != (closele->second).end(); eleIter++)
    {
      ele=dis->gElement(*eleIter);
      iele=idis->gElement(*eleIter);

      DRT::ELEMENTS::FluidImmersedBase* immersedelebase = dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(iele);
      if(immersedelebase==NULL)
        dserror("dynamic cast from DRT::Element* to DRT::ELEMENTS::FluidImmersedBase* failed");

      // evaluate this element and fill vector with immersed dirichlets
      int row = strategy->FirstDofSet();
      int col = strategy->SecondDofSet();

      params.set<Teuchos::RCP<GEO::SearchTree> >("structsearchtree_rcp",structsearchtree);
      params.set<std::map<int,LINALG::Matrix<3,1> >* >("currpositions_struct",currpositions_struct);
      params.set<int>("Physical Type",INPAR::FLUID::poro_p1);

      DRT::Element::LocationArray la(dis->NumDofSets());
      ele->LocationVector(*dis,la,false);
      strategy->ClearElementStorage( la[row].Size(), la[col].Size() );

      if(!evaluateonlyboundary)
        ele->Evaluate(params,*dis,la,
            strategy->Elematrix1(),
            strategy->Elematrix2(),
            strategy->Elevector1(),
            strategy->Elevector2(),
            strategy->Elevector3());
      else
      {
        if(immersedelebase->IsBoundaryImmersed())
          ele->Evaluate(params,*dis,la,
              strategy->Elematrix1(),
              strategy->Elematrix2(),
              strategy->Elevector1(),
              strategy->Elevector2(),
              strategy->Elevector3());
      }

      strategy->AssembleVector1( la[row].lm_, la[row].lmowner_ );
    }
  }
} // EvaluateWithInternalCommunication

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/// Reduces to standard EvaluateCondition on one proc.
/// Evaluate a specific condition using assemble strategy allowing communication at element level
/// until every conditioned element is evaluated. Needed especially during interpolation from an
/// other discretization to the conditioned elements (e.g. in immersed method).
/// The integration point of a conditioned element requesting a quantity may be owned by another
/// proc as the interpolating element providing this quantity.  rauch 05/14
void IMMERSED::ImmersedBase::EvaluateInterpolationCondition
(
    Teuchos::RCP<DRT::Discretization> evaldis,
    Teuchos::ParameterList& params,
    DRT::AssembleStrategy & strategy,
    const std::string& condstring,
    const int condid
)
{
# ifdef DEBUG
  if (!(evaldis->Filled()) ) dserror("FillComplete() was not called");
  if (!(evaldis->HaveDofs())) dserror("AssignDegreesOfFreedom() was not called");
# endif

  int row = strategy.FirstDofSet();
  int col = strategy.SecondDofSet();

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  params.set<int>("dummy_call",0);

  DRT::Element::LocationArray la(evaldis->NumDofSets());

  std::multimap<std::string,Teuchos::RCP<DRT::Condition> >::iterator fool;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (fool=evaldis->GetAllConditions().begin(); fool!=evaldis->GetAllConditions().end(); ++fool)
  {
    if (fool->first == condstring)
    {
      DRT::Condition& cond = *(fool->second);
      if (condid == -1 || condid ==cond.GetInt("ConditionID"))
      {
        std::map<int,Teuchos::RCP<DRT::Element> >& geom = cond.Geometry();
        if (geom.empty()) dserror("evaluation of condition with empty geometry on proc %d",evaldis->Comm().MyPID());

        std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;

        // Evaluate Loadcurve if defined. Put current load factor in parameterlist
        const std::vector<int>* curve  = cond.Get<std::vector<int> >("curve");
        int curvenum = -1;
        if (curve) curvenum = (*curve)[0];
        double curvefac = 1.0;
        if (curvenum>=0 && usetime)
          curvefac = DRT::Problem::Instance()->Funct(curvenum).EvaluateTime(time);

        // Get ConditionID of current condition if defined and write value in parameterlist
        const std::vector<int>*    CondIDVec  = cond.Get<std::vector<int> >("ConditionID");
        if (CondIDVec)
        {
          params.set("ConditionID",(*CondIDVec)[0]);
          char factorname[30];
          sprintf(factorname,"LoadCurveFactor %d",(*CondIDVec)[0]);
          params.set(factorname,curvefac);
        }
        else
        {
          params.set("LoadCurveFactor",curvefac);
        }
        params.set<Teuchos::RCP<DRT::Condition> >("condition", fool->second);

        int mygeometrysize=-1234;
        if(geom.empty()==true)
          mygeometrysize=0;
        else
          mygeometrysize=geom.size();
        int maxgeometrysize=-1234;
        evaldis->Comm().MaxAll(&mygeometrysize,&maxgeometrysize,1);
        curr=geom.begin();

#ifdef DEBUG
        std::cout<<"PROC "<<evaldis->Comm().MyPID()<<": mygeometrysize = "<<mygeometrysize<<" maxgeometrysize = "<<maxgeometrysize<<std::endl;
#endif


        // enter loop on every proc until the last proc evaluated his last geometry element
        // because there is communication happening inside
        for (int i=0;i<maxgeometrysize;++i)
        {
          if(i>=mygeometrysize)
            params.set<int>("dummy_call",1);

          // get element location vector and ownerships
          // the LocationVector method will return the the location vector
          // of the dofs this condition is meant to assemble into.
          // These dofs do not need to be the same as the dofs of the element
          // (this is the standard case, though). Special boundary conditions,
          // like weak dirichlet conditions, assemble into the dofs of the parent element.
          curr->second->LocationVector(*evaldis,la,false,condstring,params);

          // get dimension of element matrices and vectors
          // Reshape element matrices and vectors and init to zero

          strategy.ClearElementStorage( la[row].Size(), la[col].Size() );

          // call the element specific evaluate method
          int err = curr->second->Evaluate(params,*evaldis,la[0].lm_,
              strategy.Elematrix1(),
              strategy.Elematrix2(),
              strategy.Elevector1(),
              strategy.Elevector2(),
              strategy.Elevector3());
          if (err) dserror("error while evaluating elements");

          // assemble every element contribution only once
          // do not assemble after dummy call for internal communication
          if(i<mygeometrysize)
          {
            // assembly
            int eid = curr->second->Id();
            strategy.AssembleMatrix1( eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_ );
            strategy.AssembleMatrix2( eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_ );
            strategy.AssembleVector1( la[row].lm_, la[row].lmowner_ );
            strategy.AssembleVector2( la[row].lm_, la[row].lmowner_ );
            strategy.AssembleVector3( la[row].lm_, la[row].lmowner_ );
          }

          // go to next element
          if (i<(mygeometrysize-1))
            ++curr;

        } // for 0 to max. geometrysize over all procs
      } // if check of condid successful
    } // if condstring found
  } //for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  return;
} // EvaluateInterpolationCondition

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::SearchPotentiallyCoveredBackgrdElements(
        std::map<int,std::set<int> >* current_subset_tofill,
        Teuchos::RCP<GEO::SearchTree> backgrd_SearchTree,
        const DRT::Discretization& dis,
        const std::map<int, LINALG::Matrix<3, 1> >& currentpositions,
        const LINALG::Matrix<3, 1>& point,
        const double radius,
        const int label)
{
  *current_subset_tofill = backgrd_SearchTree->searchElementsInRadius(dis,currentpositions,point,radius,label);
  return;
} // SearchPotentiallyCoveredBackgrdElements


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::EvaluateSubsetElements(Teuchos::ParameterList& params,
                                                    Teuchos::RCP<DRT::Discretization> dis,
                                                    std::map<int,std::set<int> >& elementstoeval,
                                                    int action)
{
  // pointer to element
  DRT::Element* ele;

  // initialize location array
  DRT::Element::LocationArray la(1);

  for(std::map<int, std::set<int> >::const_iterator closele = elementstoeval.begin(); closele != elementstoeval.end(); closele++)
  {
    for(std::set<int>::const_iterator eleIter = (closele->second).begin(); eleIter != (closele->second).end(); eleIter++)
    {
        ele = dis->gElement(*eleIter);

        Epetra_SerialDenseMatrix dummymatrix;
        Epetra_SerialDenseVector dummyvector;
        ele->Evaluate(params,*dis,la,dummymatrix,dummymatrix,dummyvector,dummyvector,dummyvector);
    }
  }

  return;
} // EvaluateSubsetElements


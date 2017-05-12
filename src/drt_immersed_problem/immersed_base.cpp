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
#include "../drt_io/io_control.H"
#include "../linalg/linalg_utils.H"
#include "../drt_adapter/ad_fld_wrapper.H"
#include "../drt_adapter/ad_str_fsiwrapper_immersed.H"
#include "../drt_fluid_ele/fluid_ele_poro_immersed.H"


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
                                                   const std::string condname,
                                                   bool buildgeometry)
{
  // determine id of condition
  std::multimap<std::string,Teuchos::RCP<DRT::Condition> > allconditions;
  allconditions = dis->GetAllConditions();
  int id = (int)allconditions.size();
  id += 1;

  // build condition
  Teuchos::RCP<DRT::Condition> condition =
          Teuchos::rcp(new DRT::Condition(id,condtype,buildgeometry,DRT::Condition::Volume));

  // add nodes to conditions
   condition->Add("Node Ids",dvol_fenode);

   // add condition to discretization
   dis->SetCondition(condname,condition);

   // fill complete if necessary
   if (!dis->Filled())
     dis -> FillComplete(false,false,buildgeometry);

   std::map<int,Teuchos::RCP<DRT::Element> >& geom = dis->GetCondition(condname)->Geometry();
   std::map<int,Teuchos::RCP<DRT::Element> >::iterator it;
   for(it=geom.begin();it!=geom.end();it++)
   {
     int id = it->second->Id();
     dis->gElement(id)->SetCondition(condname,condition);
   }

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
    const Teuchos::RCP<const Epetra_Map>&     dbcmap_new)
{
  int mynumvals = dbcmap_new->NumMyElements();
  double* myvals = dirichvals->Values();

  for(int i=0;i<mynumvals;++i)
  {
    int gid = dbcmap_new->GID(i);

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
void IMMERSED::ImmersedBase::DoDirichletCond(
    const Teuchos::RCP<Epetra_Vector>&        statevector,
    const Teuchos::RCP<const Epetra_Vector>&  dirichvals,
    const Teuchos::RCP<const Epetra_Map>&     dbcmap_new,
    const Teuchos::RCP<const Epetra_Map>&     dbcmap_orig)
{
  int mynumvals = dbcmap_new->NumMyElements();
  double* myvals = dirichvals->Values();

  for(int i=0;i<mynumvals;++i)
  {
    int gid = dbcmap_new->GID(i);

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

    // we do not want to overwrite original values
    if(dbcmap_orig->LID(gid) == -1)
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
  const Teuchos::RCP<const Epetra_Map> condmap_orig =
      field_wrapper->GetDBCMapExtractor()->CondMap();

  // build map of dofs subjected to Dirichlet condition
  BuildConditionDofMap(
      dis,
      condname,
      condmap_orig,
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
  // save the original condition map
  const Teuchos::RCP<const Epetra_Map> condmap_orig =
      field_wrapper->GetDBCMapExtractor()->CondMap();

  // build map of dofs subjected to Dirichlet condition
  BuildConditionDofMap(
      dis,
      condname,
      condmap_orig,
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
void IMMERSED::ImmersedBase::ApplyDirichletToArtificialECM(
    const Teuchos::RCP< ::ADAPTER::StructureWrapper>& field_wrapper,
    const Teuchos::RCP<DRT::Discretization>&  dis,
    const std::string                         condname,
    Teuchos::RCP<Epetra_Map>&                 cond_dofrowmap,
    const int                                 numdof,
    const Teuchos::RCP<const Epetra_Vector>&  dirichvals)
{
  const Teuchos::RCP<const Epetra_Map> cond_dofmap_orig =
      field_wrapper->GetDBCMapExtractor()->CondMap();

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
    bool attachedtoimmersedelement = false;
    // get node id
    int nodeid = cond_nodes->at(node);
    // is node on this proc ?
    if( DRT::Problem::Instance()->GetDis("porofluid")->NodeColMap()->MyGID(nodeid) )
    {
      // get node pointer
      DRT::Node* node_ptr = DRT::Problem::Instance()->GetDis("porofluid")->gNode(nodeid);

      if(node_ptr->Owner()==dis->Comm().MyPID())
      {
        // get adjacent elements
        DRT::Element** element = node_ptr->Elements();
        int numadjacentelements = node_ptr->NumElement();
        for(int aele=0; aele<numadjacentelements;aele++)
        {
          DRT::Element* iele = DRT::Problem::Instance()->GetDis("porofluid")->gElement(element[aele]->Id());
          DRT::ELEMENTS::FluidPoroImmersed* immersedele=dynamic_cast<DRT::ELEMENTS::FluidPoroImmersed* >(iele);
          if(immersedele==NULL)
            dserror("Dynamic cast from DRT::Element* to DRT::ELEMENTS::FluidPoroImmersed* failed");

          if(immersedele->IsImmersed())
          {
            if(immersedele->IsBoundaryImmersed())
              dserror("Element must not be labeled IsImmersed and IsImmersedBoundary at once!");

            attachedtoimmersedelement=true;
            break;
          }
        }// loop over all adjacent elements

        IMMERSED::ImmersedNode* immersed_node_ptr =
            dynamic_cast<IMMERSED::ImmersedNode* >(node_ptr);
        if(node_ptr==NULL)
          dserror("Dynamic cast from DRT::Node* to IMMERSED::ImmersedNode* failed");

        // we only set Dirichlet values to fully artificial nodes,
        // which do not belong to pseudo-boundary
        if((immersed_node_ptr->IsBoundaryImmersed()==false) and attachedtoimmersedelement)
        {
          // get dofs
          std::vector<int> dofs = dis->Dof(0,node_ptr);

          for (int dim=0;dim<numdof;++dim)
          {
            // if not already in original dirich map
            if(cond_dofmap_orig->LID(dofs[dim]) == -1)
              mydirichdofs.push_back(dofs[dim]);
          }
        }
      }// if proc owns node
    }// if proc has node
  } // loop over all conditioned nodes

  int nummydirichvals = mydirichdofs.size();

  // construct new dof row map
  cond_dofrowmap = Teuchos::rcp( new Epetra_Map(-1,nummydirichvals,&(mydirichdofs[0]),0,dis->Comm()) );

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
void IMMERSED::ImmersedBase::ApplyDirichletToRealECM(
    const Teuchos::RCP< ::ADAPTER::StructureWrapper>& field_wrapper,
    const Teuchos::RCP<DRT::Discretization>&  dis,
    const std::string                         condname,
    Teuchos::RCP<Epetra_Map>&                 cond_dofrowmap,
    const int                                 numdof,
    const Teuchos::RCP<const Epetra_Vector>&  dirichvals)
{
  const Teuchos::RCP<const Epetra_Map> cond_dofmap_orig =
      field_wrapper->GetDBCMapExtractor()->CondMap();

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
    bool attachedtoimmersedelement = false;
    bool attachedtoimmersedbdryelement = false;
    // get node id
    int nodeid = cond_nodes->at(node);
    // is node on this proc ?
    if( DRT::Problem::Instance()->GetDis("porofluid")->NodeColMap()->MyGID(nodeid) )
    {
      // get node pointer
      DRT::Node* node_ptr = DRT::Problem::Instance()->GetDis("porofluid")->gNode(nodeid);
      if(node_ptr==NULL)
        dserror("Could not get node with id %d",nodeid);
      if(node_ptr->Owner()==dis->Comm().MyPID())
      {
        // get adjacent elements
        DRT::Element** element = node_ptr->Elements();
        int numadjacentelements = node_ptr->NumElement();
        for(int aele=0; aele<numadjacentelements;aele++)
        {
          DRT::Element* iele = DRT::Problem::Instance()->GetDis("porofluid")->gElement(element[aele]->Id());
          DRT::ELEMENTS::FluidPoroImmersed* immersedele=dynamic_cast<DRT::ELEMENTS::FluidPoroImmersed* >(iele);
          if(immersedele==NULL)
            dserror("Dynamic cast from DRT::Element* to DRT::ELEMENTS::FluidPoroImmersed* failed");

          if(immersedele->IsImmersed())
          {
            attachedtoimmersedelement=true;
          }

          if(immersedele->IsBoundaryImmersed())
          {
            attachedtoimmersedbdryelement=true;
          }

        } // loop over all adjacent elements

        IMMERSED::ImmersedNode* immersed_node_ptr =
            dynamic_cast<IMMERSED::ImmersedNode* >(node_ptr);
        if(immersed_node_ptr==NULL)
          dserror("Dynamic cast from DRT::Node* to IMMERSED::ImmersedNode* failed");

        // we only set Dirichlet values to fully artificial nodes,
        // which do not belong to pseudo-boundary
        if(attachedtoimmersedbdryelement or (attachedtoimmersedbdryelement == false and attachedtoimmersedelement == false))
        {
          // get dofs
          std::vector<int> dofs = dis->Dof(0,node_ptr);

          for (int dim=0;dim<numdof;++dim)
          {
            // if not already in original dirich map
            if(cond_dofmap_orig->LID(dofs[dim]) == -1)
              mydirichdofs.push_back(dofs[dim]);
          }
        }
      }// if node is owned by proc
    }// if proc has node
  }// loop over all conditioned nodes

  int nummydirichvals = mydirichdofs.size();
  if(nummydirichvals == 0)
    dserror("Proc %i does not add any dof to Dirichlet map.",dis->Comm().MyPID());

  // construct new dof row map
  cond_dofrowmap = Teuchos::rcp( new Epetra_Map(-1,nummydirichvals,&(mydirichdofs[0]),0,dis->Comm()) );

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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedBase::WriteExtraOutput(
    const Epetra_Comm& comm,
    const double time,
    const std::string filenameending,
    const std::vector<double> valuetowrite,
    const std::vector<double> valuetowrite2,
    const std::vector<double> valuetowrite3)
{
  // append values to output file
  if (comm.MyPID()==0)
  {
    const std::string fname1 = DRT::Problem::Instance()->OutputControlFile()->FileName()+"."+filenameending;

    std::ofstream f1;
    f1.open(fname1.c_str(),std::fstream::ate | std::fstream::app);

    f1 << time << " " << valuetowrite[0] << " " << valuetowrite[1] << " " << valuetowrite[2] << " " << valuetowrite[3]<< " " <<
                         valuetowrite2[0] << " " << valuetowrite2[1] << " " << valuetowrite2[2] << " " << valuetowrite2[3] << " " <<
                         valuetowrite3[0] << " " << valuetowrite3[1] << " " << valuetowrite3[2] << " " << valuetowrite3[3] << "   ";

    f1 << "\n";
    f1.flush();
    f1.close();

  } // only proc 0 writes

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<double> IMMERSED::ImmersedBase::
    CalcGlobalResultantfromEpetraVector(
        const Epetra_Comm& comm,
        const Teuchos::RCP<const DRT::Discretization>& dis,
        const Teuchos::RCP<const Epetra_Vector>& vec_epetra)
{
  double summyrowentriesx= 0.0;
  double summyrowentriesy= 0.0;
  double summyrowentriesz= 0.0;
  double result_globalx  = 0.0;
  double result_globaly  = 0.0;
  double result_globalz  = 0.0;
  double result_L2norm   = 0.0;

  std::vector<double> result;

  const int nummyrownodes = dis->NumMyRowNodes();
  const int myveclength = vec_epetra->MyLength();

  if(myveclength!=nummyrownodes*3)
    dserror("local vector length is invalid!");

  for(int i=0; i<nummyrownodes; i++)
  {
    summyrowentriesx += vec_epetra->Values()[i*3+0];
    summyrowentriesy += vec_epetra->Values()[i*3+1];
    summyrowentriesz += vec_epetra->Values()[i*3+2];
  }

  comm.Barrier();
  comm.SumAll(&summyrowentriesx,&result_globalx,1);
  comm.SumAll(&summyrowentriesy,&result_globaly,1);
  comm.SumAll(&summyrowentriesz,&result_globalz,1);

  result_L2norm = sqrt(pow(result_globalx,2)+pow(result_globaly,2)+pow(result_globalz,2));

  result.push_back(result_globalx);
  result.push_back(result_globaly);
  result.push_back(result_globalz);
  result.push_back(result_L2norm);

  return result;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::vector<int> IMMERSED::ImmersedBase::DetermineElesInFirstImmersedRow(
    const Teuchos::RCP<const DRT::Discretization>& dis,
    std::vector<int>* eleids)
{
  std::vector<int> returnvar;
  bool isconnectedtointersectedelement=false;

  const Epetra_Map* elecolmapptr = dis->ElementColMap();
  const int nummyelements = elecolmapptr->NumMyElements();

  // loop over all row elements
  for(int lid=0;lid<nummyelements;lid++)
  {
    isconnectedtointersectedelement=false;

    DRT::Element* eleptr = dis->gElement(elecolmapptr->GID(lid));
    DRT::ELEMENTS::FluidPoroImmersed* immersedeleptr =
        dynamic_cast<DRT::ELEMENTS::FluidPoroImmersed* >(eleptr);
    if(immersedeleptr==NULL)
      dserror("dynamic cast to FluidPoroImmersed failed!");

    // if element is labeled as IsImmersed
    if(immersedeleptr->IsImmersed())
    {
      // safety check
      if(immersedeleptr->IsBoundaryImmersed())
        dserror("Element cannot be IsImmersed and IsBoundaryImmersed at once!");

      // loop over all nodes of immersed element
      DRT::Node** nodesptr = immersedeleptr->Nodes();
      for (int node=0;node<immersedeleptr->NumNode();node++)
      {
        DRT::Node* nodeptr = nodesptr[node];
        IMMERSED::ImmersedNode* immersednodeptr = dynamic_cast<IMMERSED::ImmersedNode* >(nodeptr);
        if(immersednodeptr==NULL)
          dserror("dynamic cast to ImmersedNode failed!");

        if(immersednodeptr->IsBoundaryImmersed())
        {
          eleids->push_back(immersedeleptr->Id());
          immersedeleptr->SetIsImmersedFirstRow(true);
          isconnectedtointersectedelement = true;
          break;
        }
      } // loop over nodes of immersed element

      // add nodeids to vector
      if(isconnectedtointersectedelement)
      {
        for (int node=0;node<immersedeleptr->NumNode();node++)
        {
          dynamic_cast<IMMERSED::ImmersedNode* >(nodesptr[node])->SetIsPseudoBoundary(true);
          returnvar.push_back(nodesptr[node]->Id());
        }
      }

    } // if immersed element
  } // loop over all row elements

  return returnvar;
}

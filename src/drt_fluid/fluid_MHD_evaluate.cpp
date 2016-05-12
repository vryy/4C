/*----------------------------------------------------------------------*/
/*!
\file fluid_MHD_evaluate.cpp

\brief Class FLD::FluidMHDEvaluate

<pre>

\maintainer Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236

\level 1
</pre>

*/

#include "../drt_fluid/fluid_MHD_evaluate.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_fluid/drt_periodicbc.H"
#include "../drt_lib/drt_dofset_transparent.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_lib/drt_utils_parmetis.H"
#include "../drt_lib/drt_element.H"

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Destructor dtor (public)                                  gammi 02/11|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::FluidMHDEvaluate::~FluidMHDEvaluate()
{
  return;
}// FluidMHDEvaluate::~FluidMHDEvaluate()


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                     gammi 02/11|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
FLD::FluidMHDEvaluate::FluidMHDEvaluate(
  Teuchos::RCP<DRT::Discretization>    actdis
  )
  :
  pdiscret_(actdis)
{
  std::vector<DRT::Condition*> MHDcnd;
  pdiscret_->GetCondition("SurfaceMixHybDirichlet",MHDcnd);

  std::vector <int> allcnd_MHDnodeids;

  if(MHDcnd.size()!=0)
  {
    if(pdiscret_->Comm().MyPID()==0)
    {
      printf("+----------------\n");
      printf("|\n");
      printf("| Generating a boundary discretisation for all elements next ");
      printf("to a mixed/hybrid\n");
      printf("| Dirichlet boundary\n");
      printf("|\n");
    }

    // generate an empty boundary discretisation
    bnd_discret_
      = Teuchos::rcp(new DRT::Discretization((std::string)"boundary discretisation",
                                    Teuchos::rcp(pdiscret_->Comm().Clone())));

    // make the condition known to the boundary discretisation
    for (unsigned numcond=0;numcond<MHDcnd.size();++numcond)
    {
      // We use the same nodal ids and therefore we can just copy the
      // conditions.
      bnd_discret_->SetCondition("SurfaceMixHybDirichlet",
                                 Teuchos::rcp(new DRT::Condition(*MHDcnd[numcond])));
    }

    // get set of ids of all MHD nodes
    std::set<int> MHDnodeset;
    {

      for (unsigned numcond=0;numcond<MHDcnd.size();++numcond)
      {
        const std::vector <int>* MHDnodeids = (*MHDcnd[numcond]).Nodes();

        allcnd_MHDnodeids.reserve( allcnd_MHDnodeids.size() + MHDnodeids->size());
        allcnd_MHDnodeids.insert ( allcnd_MHDnodeids.end(), MHDnodeids->begin(), MHDnodeids->end());
      }

      for(std::vector<int>::iterator id =allcnd_MHDnodeids.begin();
          id!=allcnd_MHDnodeids.end();
          ++id)
      {
        MHDnodeset.insert(*id);
      }
    }

    // determine sets of nodes next to MHD nodes
    std::set<int> adjacent_row;
    std::set<int> adjacent_col;

    // loop all column elements and label all row nodes next to a MHD node
    for (int i=0; i<pdiscret_->NumMyColElements(); ++i)
    {
      DRT::Element* actele = pdiscret_->lColElement(i);

      // get the node ids of this element
      const int  numnode = actele->NumNode();
      const int* nodeids = actele->NodeIds();

      bool found=false;

      // loop nodeids, check if a MHD condition is active
      for(int rr=0;rr<numnode;++rr)
      {
        int gid=nodeids[rr];

        std::set<int>::iterator curr=MHDnodeset.find(gid);
        if(curr!=MHDnodeset.end())
        {
          found=true;
        }
      }

      // yes, we have a MHD condition
      if(found==true)
      {
        // loop nodeids
        for(int rr=0;rr<numnode;++rr)
        {
          int gid=nodeids[rr];

          if ((pdiscret_->NodeRowMap())->LID(gid)>-1)
          {
            adjacent_row.insert(gid);
          }
          adjacent_col.insert(gid);
        }
      }
    }

    // all row nodes next to a MHD node are now contained in the bndydis
    for(std::set<int>::iterator id = adjacent_row.begin();
        id!=adjacent_row.end();
        ++id)
    {
      DRT::Node* actnode=pdiscret_->gNode(*id);

      Teuchos::RCP<DRT::Node> bndnode =Teuchos::rcp(actnode->Clone());

      bnd_discret_->AddNode(bndnode);
    }

    // loop all row elements and add all elements with a MHD node
    for (int i=0; i<pdiscret_->NumMyRowElements(); ++i)
    {
      DRT::Element* actele = pdiscret_->lRowElement(i);

      // get the node ids of this element
      const int  numnode = actele->NumNode();
      const int* nodeids = actele->NodeIds();

      bool found=false;

      // loop nodeids, check if a MHD condition is active
      for(int rr=0;rr<numnode;++rr)
      {
        int gid=nodeids[rr];

        std::set<int>::iterator curr=MHDnodeset.find(gid);
        if(curr!=MHDnodeset.end())
        {
          found=true;
        }
      }

      // yes, we have a MHD condition
      if(found==true)
      {
        Teuchos::RCP<DRT::Element> bndele =Teuchos::rcp(actele->Clone());

        bnd_discret_->AddElement(bndele);
      }
    }

    // bndydis needs a full NodeRowMap and a NodeColMap
    Teuchos::RCP<Epetra_Map> newrownodemap;
    Teuchos::RCP<Epetra_Map> newcolnodemap;

    {

      std::vector<int> rownodes;

      // all row nodes next to a MHD node are now contained in the bndydis
      for(std::set<int>::iterator id = adjacent_row.begin();
          id!=adjacent_row.end();
          ++id)
      {
        rownodes.push_back(*id);
      }

      // build noderowmap for new distribution of nodes
      newrownodemap = Teuchos::rcp(new Epetra_Map(-1,
                                         rownodes.size(),
                                         &rownodes[0],
                                         0,
                                         bnd_discret_->Comm()));

      std::vector<int> colnodes;

      for(std::set<int>::iterator id = adjacent_col.begin();
          id!=adjacent_col.end();
          ++id)
      {
        colnodes.push_back(*id);
      }
      // build nodecolmap for new distribution of nodes
      newcolnodemap = Teuchos::rcp(new Epetra_Map(-1,
                                         colnodes.size(),
                                         &colnodes[0],
                                         0,
                                         bnd_discret_->Comm()));
    }

    if(bnd_discret_->Comm().MyPID()==0)
    {
      printf("| Redistribute according to the initial nodemaps\n");
    }

    bnd_discret_->Redistribute(*newrownodemap,*newcolnodemap,false,false,false);

    if(bnd_discret_->Comm().MyPID()==0)
    {
      std::cout << "| ... done.\n";
    }

    {
      if(bnd_discret_->Comm().MyPID()==0)
      {
        std::cout << "| Inherit periodic boundary conditions, redistribute again ";
        std::cout << "to fetch slave nodes\n";
        std::cout << "| to the master's proc\n";
      }

      // make the pbc condition known to the boundary discretisation
      std::vector<DRT::Condition*> mysurfpbcs;

      // get periodic surface boundary conditions
      pdiscret_->GetCondition("SurfacePeriodic",mysurfpbcs);

      for (unsigned numcond=0;numcond<mysurfpbcs.size();++numcond)
      {
        // We use the same nodal ids --- nevertheless, we just use a subset
        // of the node ids and thus cannot copy the conditions completely.
        std::vector<int> reduced_ids;

        const std::vector<int>* candidates = (*mysurfpbcs[numcond]).Nodes();

        std::vector<int> mytoggle(candidates->size(),0);
        std::vector<int> toggle(candidates->size(),0);

        for(unsigned rr=0;rr<candidates->size();++rr)
        {
          if(newrownodemap->LID((*candidates)[rr])>-1)
          {
            mytoggle[rr]=1;
          }
        }

        bnd_discret_->Comm().SumAll(&mytoggle[0],&toggle[0],toggle.size());

        for(unsigned rr=0;rr<candidates->size();++rr)
        {
          if(toggle[rr]>0)
          {
            reduced_ids.push_back((*candidates)[rr]);
          }
        }

        (*mysurfpbcs[numcond]).Delete("Node Ids");
        (*mysurfpbcs[numcond]).Add("Node Ids",reduced_ids);

        bnd_discret_->SetCondition("SurfacePeriodic", Teuchos::rcp(new DRT::Condition(*mysurfpbcs[numcond])));
      }

      PeriodicBoundaryConditions pbc(bnd_discret_,false);
      pbc.UpdateDofsForPeriodicBoundaryConditions();

      if(bnd_discret_->Comm().MyPID()==0)
      {
        std::cout << "| ... done.\n";
      }
    }

    if(bnd_discret_->Comm().MyPID()==0)
    {
      std::cout << "| Replace dofset by a transparent dofset that copies ";
      std::cout << "the dofs of the original\n";
      std::cout << "| (parent) discretisation. At this place";
      std::cout << " a sub-dofrowmap (identical layout) of\n";
    }

    // idea: use a transparent dofset and hand through the dof numbering

    bnd_discret_->ReplaceDofSet(Teuchos::rcp(new DRT::TransparentDofSet(pdiscret_,true)));

    bnd_discret_->Redistribute(*newrownodemap,*newcolnodemap,true,true,true);

    if(bnd_discret_->Comm().MyPID()==0)
    {
      std::cout << "| the parent discretisation is generated. It is used to define a system\n";
      std::cout << "| matrix for the boundary ";
      std::cout << "dofs, which is filled and assembled into the global\n";
      std::cout << "| matrix later on.\n";
    }

    subdofrowmap_=Teuchos::rcp(new Epetra_Map(*bnd_discret_->DofRowMap()));
    Epetra_Map subdofcolmap(*bnd_discret_->DofColMap());


    bndmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*subdofrowmap_,
                                                   500,
                                                   false,
                                                   true,
                                                   LINALG::SparseMatrix::FE_MATRIX));

    if(bnd_discret_->Comm().MyPID()==0)
    {
      std::cout << "| ... done.\n";
    }

    if(bnd_discret_->Comm().MyPID()==0)
    {
      std::cout << "| Call PARMETIS on the boundary discretisation and ";
      std::cout << "redistribute according to\n";
      std::cout << "| the new maps\n";
    }

    std::vector<int> bndnids;
    std::vector<int> bndnidslocal(bnd_discret_->NodeRowMap()->NumMyElements());

    for (int i=0; i<bnd_discret_->NodeRowMap()->NumMyElements(); ++i)
    bndnidslocal[i] = bnd_discret_->NodeRowMap()->GID(i);

    const int numproc = pdiscret_->Comm().NumProc();

    // vector containing all proc ids
    std::vector<int> allproc(numproc);
    for (int i=0; i<numproc; ++i) allproc[i] = i;

    LINALG::Gather<int>(bndnidslocal,bndnids,numproc,&allproc[0],pdiscret_->Comm());

    //**********************************************************************
    // call PARMETIS (again with #ifdef to be on the safe side)
#if defined(PARALLEL) && defined(HAVE_PARMETIS)

    Teuchos::RCP<Epetra_Map> bndrownodes;
    Teuchos::RCP<Epetra_Map> bndcolnodes;

    Teuchos::RCP<Epetra_Map> belemap = Teuchos::rcp( new Epetra_Map(*bnd_discret_->ElementRowMap()));
    Epetra_Time time(pdiscret_->Comm());
    Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(pdiscret_->Comm().Clone());

    DRT::UTILS::PartUsingParMetis(bnd_discret_,
                                  belemap,
                                  bndrownodes,
                                  bndcolnodes,
                                  comm,
                                  false);

#else
#if defined(PARALLEL)
    dserror("require PARMETIS not METIS");
#else
    bndrownodes=Teuchos::rcp(new Epetra_Map(*newrownodemap));
    bndcolnodes=Teuchos::rcp(new Epetra_Map(*newcolnodemap));
#endif
#endif
    if(bnd_discret_->Comm().MyPID()==0)
    {
      std::cout << "| Redistributing .";
    }
    bnd_discret_->Redistribute(*bndrownodes,*bndcolnodes,false,false);

    if(bnd_discret_->Comm().MyPID()==0)
    {
      std::cout << ".. done.\n";
    }

    if(bnd_discret_->Comm().MyPID()==0)
    {
      std::cout << "| Apply periodic boundary conditions to the redistributed";
      std::cout << " discretisation to\n";
      std::cout << "| fetch slave nodes to the master's proc\n";
    }

    {
      PeriodicBoundaryConditions pbc(bnd_discret_,false);
      pbc.UpdateDofsForPeriodicBoundaryConditions();
    }

    if(bnd_discret_->Comm().MyPID()==0)
    {
      std::cout << "| Assign the dofs for the redistributed layout, again using ";
      std::cout << "a parallel version\n";
      std::cout << "| of the transparent dofset\n";
    }

    // idea: use a transparent dofset and hand through the dof numbering
    bnd_discret_->ReplaceDofSet(Teuchos::rcp(new DRT::TransparentDofSet(pdiscret_,true)));

    bnd_discret_->FillComplete();

    if(bnd_discret_->Comm().MyPID()==0)
    {
      std::cout << "| ... done.\n";
      printf("|\n");
      printf("+----------------\n\n");

    }

    {
      std::vector<int> my_n_nodes   (numproc,0);
      std::vector<int>    n_nodes   (numproc,0);
      std::vector<int> my_n_elements(numproc,0);
      std::vector<int>    n_elements(numproc,0);
      std::vector<int> my_n_ghostele(numproc,0);
      std::vector<int>    n_ghostele(numproc,0);
      std::vector<int> my_n_dof     (numproc,0);
      std::vector<int>    n_dof     (numproc,0);

      int myrank=bnd_discret_->Comm().MyPID();

      my_n_nodes   [myrank]=bnd_discret_->NodeRowMap()->NumMyElements ();
      my_n_elements[myrank]=bnd_discret_->NumMyColElements();
      my_n_ghostele[myrank]=bnd_discret_->NumMyColElements()-bnd_discret_->NumMyRowElements();
      my_n_dof     [myrank]=bnd_discret_->DofRowMap()->NumMyElements();

      bnd_discret_->Comm().SumAll(&my_n_nodes[0]   ,&n_nodes[0]   ,numproc);
      bnd_discret_->Comm().SumAll(&my_n_elements[0],&n_elements[0],numproc);
      bnd_discret_->Comm().SumAll(&my_n_ghostele[0],&n_ghostele[0],numproc);
      bnd_discret_->Comm().SumAll(&my_n_dof[0]     ,&n_dof[0]     ,numproc);

      if(bnd_discret_->Comm().MyPID()==0)
      {
        printf("   +-----+---------------+-----------------+----------------+-----------------+\n");
        printf("   +                       boundary discretisation                            +\n");
        printf("   +-----+---------------+-----------------+----------------+-----------------+\n");
        printf("   | PID |    n_nodes    |    n_elements   |   n_ghostele   |      n_dof      |\n");
        printf("   +-----+---------------+-----------------+----------------+-----------------+\n");
        for(int npid=0;npid<numproc ;++npid)
        {
          printf("   | %3d | %13d | %15d | %14d | %15d |\n",npid,n_nodes[npid],n_elements[npid],n_ghostele[npid],n_dof[npid]);
          printf("   +-----+---------------+-----------------+----------------+-----------------+\n");
        }
        std::cout << std::endl <<std::endl;
      }
    }

    // ---------------------------------------------------------------
    // The remaining part are just sanity checks for the redistributed discretisation

    bool insane=false;

    // loop all column eles, check dofs for each node
    for (int i=0; i<bnd_discret_->NumMyColElements(); ++i)
    {
      DRT::Element* actele = bnd_discret_->lColElement(i);

      // get the node ids of this element
      const int  numnode = actele->NumNode();
      const int* nodeids = actele->NodeIds();

      // loop nodeids, check if a MHD condition is active
      for(int rr=0;rr<numnode;++rr)
      {
        DRT::Node*  node = bnd_discret_->gNode(nodeids[rr]);
        std::vector<int> nodedofset = bnd_discret_->Dof(node);

        for(unsigned index=0;index<nodedofset.size();++index)
        {
          int gid = nodedofset[index];

          if(bnd_discret_->DofColMap()->LID(gid)<0)
          {
            insane=true;
            printf("myrank %d dof %d not in colmap\n",bnd_discret_->Comm().MyPID(),gid);
          }
        }
      }
    }
    if(insane) dserror("invalid dof col map");

    {
      std::set<int> testset;
      for(int rr=0;rr<bnd_discret_->DofRowMap()->NumMyElements();++rr)
      {
        int id=bnd_discret_->DofRowMap()->MyGlobalElements()[rr];

        std::set<int>::iterator curr=testset.find(id);
        if(curr!=testset.end())
        {
          dserror("DofRowMap of bnd dis is not unique on this proc");
        }
        testset.insert(id);
      }

      if (!bnd_discret_->DofRowMap()->UniqueGIDs())
      {
        std::cout  << *bnd_discret_->DofRowMap();

        dserror("DofRowMap  of bnd dis is not unique (global)");
      }
    }
  }

  return;
}







void FLD::FluidMHDEvaluate::BoundaryElementLoop(
  Teuchos::ParameterList&              mhdbcparams,
  Teuchos::RCP<Epetra_Vector>          velaf_     ,
  Teuchos::RCP<Epetra_Vector>          velnp_     ,
  Teuchos::RCP<Epetra_Vector>          residual_  ,
  Teuchos::RCP<LINALG::SparseMatrix>   sysmat_
  )
{
    // set the required state vectors
    Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*bnd_discret_->DofColMap(),true);


    LINALG::Export(*velaf_,*tmp);
    bnd_discret_->SetState("u and p (trial)",tmp);
    bnd_discret_->SetState("velaf",tmp); // was missing. gjb 12/12

    {
      Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*bnd_discret_->DofColMap(),true);
      LINALG::Export(*velnp_,*tmp);
      bnd_discret_->SetState("u and p (trial,n+1)",tmp);
      bnd_discret_->SetState("velnp",tmp); // was missing. gjb 12/12
    }

    // small sysmat and residual
    //bndmat_->Zero(); // please tell me why zero doesn't work in parallel
    bndmat_->Reset();

    Teuchos::RCP<Epetra_Vector> bndres = LINALG::CreateVector(*bnd_discret_->DofRowMap(),true);

    std::vector<DRT::Condition*> bndMHDcnd;
    const std::string condstring = "SurfaceMixHybDirichlet";
    bnd_discret_->GetCondition(condstring,bndMHDcnd);

    // evaluate all mixed hybrid Dirichlet boundary conditions
    {
      // loop all MHD conditions
      for (unsigned numcond=0;numcond<bndMHDcnd.size();++numcond)
      {
        std::map<int,Teuchos::RCP<DRT::Element> >& geom = (*bndMHDcnd[numcond]).Geometry();

        Teuchos::RCP<DRT::Condition> thiscond = Teuchos::rcp(bndMHDcnd[numcond],false);

        mhdbcparams.set<Teuchos::RCP<DRT::Condition> >("condition",thiscond);

        // define element matrices and vectors
        Epetra_SerialDenseMatrix elematrix1;
        Epetra_SerialDenseMatrix dummymat;
        Epetra_SerialDenseVector elevector1;
        Epetra_SerialDenseVector dummyvec;

        // element matrices and vectors will be reshaped
        // during the element call!

        for (std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr=geom.begin(); curr!=geom.end(); ++curr)
        {
          // get element location vector and ownerships
          // the LocationVector method will return the the location vector
          // of the dofs this condition is meant to assemble into.
          // These dofs do not need to be the same as the dofs of the element
          // (this is the standard case, though). Special boundary conditions,
          // like weak dirichlet conditions, assemble into the dofs of the parent element.
          DRT::Element::LocationArray la(1);

          curr->second->LocationVector(*bnd_discret_,la,false,condstring,mhdbcparams);

          // call the element specific evaluate method
          int err = curr->second->Evaluate(mhdbcparams,*bnd_discret_,la[0].lm_,elematrix1,dummymat,
                                           elevector1,dummyvec,dummyvec);
          if (err) dserror("error while evaluating elements");

          // assembly to all parent dofs even if we just integrated
          // over a boundary element
          int eid = curr->second->Id();

          bndmat_->FEAssemble(eid,elematrix1,la[0].lm_,la[0].lmowner_,la[0].lm_);
          LINALG::Assemble(*bndres,elevector1,la[0].lm_,la[0].lmowner_);
        } // end loop geometry elements of this conditions
      }
    }

    // complete system matrix --- do all communication internally
    bndmat_->Complete();

    // loop all local entries of my boundary matrix and add them to sysmat in
    // the same position
    // should be OK since bndmat_ is constructed on a subset of dofs of sysmat
    // in this parallel layout

    Epetra_CrsMatrix* Epetra_Crs_bndmat
      =
      dynamic_cast<Epetra_CrsMatrix*>(bndmat_->EpetraOperator().get());

    if(Epetra_Crs_bndmat==NULL)
    {
      dserror("NULL Epetra_Crs_bndmat\n");
    }

    LINALG::Add(*Epetra_Crs_bndmat,
                false,
                1.0,
                *(sysmat_->EpetraMatrix()),
                1.0);

    {
      Teuchos::RCP<Epetra_Vector> tmp=LINALG::CreateVector(*pdiscret_->DofRowMap(),true);

      Epetra_Export exporter(bndres->Map(),tmp->Map());
      int err = tmp->Export(*bndres,exporter,Add);
      if (err) dserror("Export using exporter returned err=%d",err);

      residual_->Update(1.0,*tmp,1.0);
    }
    return;
}

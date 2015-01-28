/*!----------------------------------------------------------------------
\file xfluidfluidnew.cpp
\brief Control routine for fluid-fluid (in)stationary solvers with XFEM

ATTENTION! Class is still a prototype. Does not provide full (and correct)
functionality yet!

<pre>
Maintainer:  Raffaela Kruse
             kruse@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
*----------------------------------------------------------------------*/

#include "xfluidfluidnew.H"
#include "xfluid_state_creator.H"

#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid/fluid_utils.H"

#include "../drt_lib/drt_discret_faces.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_dofset_transparent_independent.H"
#include "../drt_lib/drt_utils_parallel.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

#include "../drt_fluid_ele/fluid_ele_parameter_xfem.H"

#include "../drt_xfem/xfem_edgestab.H"
#include "../drt_xfem/xfem_dofset.H"
#include "../drt_xfem/xfluidfluid_timeInt.H"


FLD::XFluidFluidNew::XFluidFluidNew(
  const Teuchos::RCP<FLD::FluidImplicitTimeInt> & embedded_fluid,     ///< embedded fluid
  const Teuchos::RCP<DRT::Discretization>&        xfluiddis,          ///< background fluid discretization
  const Teuchos::RCP<LINALG::Solver>&             solver,             ///< fluid solver
  const Teuchos::RCP<Teuchos::ParameterList>&     params,             ///< xfluid params
  const Teuchos::RCP<IO::DiscretizationWriter>&   output,             ///< discretization writer for paraview output
  bool                                            ale_xfluid,         ///< background (XFEM) fluid in ALE-formulation
  bool                                            ale_fluid           ///< embedded fluid in ALE-formulation
) : XFluid(
  xfluiddis,
  embedded_fluid->Discretization(),
  solver,
  params,
  output,
  ale_xfluid),
  embedded_fluid_(embedded_fluid),
  alefluid_(ale_fluid)
{
  // Write background fluid mesh
 // output_->WriteMesh(0,0.0);
}

FLD::XFluidFluidNew::~XFluidFluidNew()
{
}

void FLD::XFluidFluidNew::Init()
{
  // base class init
  XFluid::Init();

  // set parameters specific for fluid-fluid coupling
  SetXFluidFluidParams();

  if(meshcoupl_dis_.size() != 1) dserror("we expect exact one mesh coupling discretization for Xfluidfluid at the moment!");

  soliddis_ = meshcoupl_dis_[0];

  // make the dofset of boundarydis be a subset of the embedded dis
  Teuchos::RCP<Epetra_Map> newcolnodemap = DRT::UTILS::ComputeNodeColMap(
      soliddis_,boundarydis_);
  soliddis_->Redistribute(*(embedded_fluid_->Discretization()->NodeRowMap()), *newcolnodemap);
  Teuchos::RCP<DRT::DofSet> newdofset=Teuchos::rcp(new DRT::TransparentIndependentDofSet(
      soliddis_,false));
  boundarydis_->ReplaceDofSet(newdofset); // do not call this with true!!
  boundarydis_->FillComplete();

  //DRT::UTILS::PrintParallelDistribution(*boundarydis_);

  INPAR::FLUID::CalcError calcerr = DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_,"calculate error");
  if ((coupling_strategy_ != INPAR::XFEM::Xfluid_Sided_Coupling) or (calcerr != INPAR::FLUID::no_error_calculation))
  {
    PrepareEmbeddedDistribution();
    CreateBoundaryEmbeddedMap();
  }

  embedded_fluid_->Init();

  //-------------------------------------------------------------------

  if (alefluid_)
    dispnpoldstate_ = Teuchos::rcp(new Epetra_Vector(*embedded_fluid_->Dispnp()));

  //--------------------------------------------------
  // Create XFluidFluid State
  //-----------------------------------------------
  const int restart = DRT::Problem::Instance()->Restart();

  if (restart)
  {
    embedded_fluid_->ReadRestart(restart);
  }

  state_ = this->GetNewState();

  //----------------------------------------------------------------
  // create auxiliary discretization of outer embedded fluid elements
  // Nitsche-Eigenvalue-Problem and for application of EOS
  // pressure term to outer embedded element
  //----------------------------------------------------------------
  // Todo: Shift to XFEM::MeshCouplingFluidFluid
  if (nitsche_evp_ || xff_eos_pres_emb_layer_)
    CreateEmbeddedBoundaryDiscretization();

  // non-stationary fluid-fluid interface requires
  // proper time integration approach
  if (alefluid_)
    xfluidfluid_timeint_ =  Teuchos::rcp(new XFEM::XFluidFluidTimeIntegration(
        XFluid::discret_,
        embedded_fluid_->Discretization(),
        state_->Wizard(), step_,
        xfem_timeintapproach_,*params_));

}

void FLD::XFluidFluidNew::SetXFluidFluidParams()
{
  Teuchos::ParameterList&   params_xf_stab = params_->sublist("XFLUID DYNAMIC/STABILIZATION");

  hybrid_lm_l2_proj_     = DRT::INPUT::IntegralValue<INPAR::XFEM::Hybrid_LM_L2_Proj>(params_xf_stab, "HYBRID_LM_L2_PROJ");

  xff_conv_stab_scaling_ = DRT::INPUT::IntegralValue<INPAR::XFEM::XFF_ConvStabScaling>(params_xf_stab,"XFF_CONV_STAB_SCALING");

  // Flag, whether any face-based terms are active
  eval_eos_ = edge_based_ || ghost_penalty_;

  // additional eos pressure stabilization on the elements of the embedded discretization,
  // that contribute to the interface
  xff_eos_pres_emb_layer_ = DRT::INPUT::IntegralValue<bool>(params_xf_stab,"XFF_EOS_PRES_EMB_LAYER");

  // whether an eigenvalue problem has to be solved to estimate Nitsche's parameter
  nitsche_evp_ = (DRT::INPUT::IntegralValue<INPAR::XFEM::ViscStab_TraceEstimate>(params_xf_stab,"VISC_STAB_TRACE_ESTIMATE")
                  == INPAR::XFEM::ViscStab_TraceEstimate_eigenvalue);

  // get general XFEM/XFFSI specific parameters
  monolithic_approach_  = DRT::INPUT::IntegralValue<INPAR::XFEM::Monolithic_xffsi_Approach>(params_->sublist("XFLUID DYNAMIC/GENERAL"),"MONOLITHIC_XFFSI_APPROACH");
  xfem_timeintapproach_ = DRT::INPUT::IntegralValue<INPAR::XFEM::XFluidFluidTimeInt>(params_->sublist("XFLUID DYNAMIC/GENERAL"),"XFLUIDFLUID_TIMEINT");

  // get information about active shape derivatives
  active_shapederivatives_ = alefluid_ && params_->get<bool>("shape derivatives");
}

void FLD::XFluidFluidNew::TimeUpdate()
{
  XFluid::Update();
  embedded_fluid_->Update();
}

void FLD::XFluidFluidNew::CutAndSetStateVectors()
{
  if (!alefluid_)
   return;

  // save the old state
  staten_ = state_;

  const int restart = DRT::Problem::Instance()->Restart();
  // if restart
  if (restart and ((restart+1) == step_))
  {
    xfluidfluid_timeint_->CreateBgNodeMapsForRestart(staten_->Wizard());
  }
}

Teuchos::RCP<LINALG::BlockSparseMatrixBase> FLD::XFluidFluidNew::BlockSystemMatrix(
    Teuchos::RCP<Epetra_Map> innermap,
    Teuchos::RCP<Epetra_Map> condmap)
{
  //Map of fluid FSI DOFs: condmap
  //Map of inner fluid DOFs: innermap

  //Get the fluid-fluid system matrix as sparse matrix
  Teuchos::RCP<LINALG::SparseMatrix> sparsesysmat = SystemMatrix();

  //F_{II}, F_{I\Gamma}, F_{\GammaI}, F_{\Gamma\Gamma}
  Teuchos::RCP<LINALG::SparseMatrix> fii, fig, fgi, fgg;
  // Split sparse system matrix into blocks according to the given maps
  LINALG::SplitMatrix2x2(sparsesysmat,innermap,condmap,innermap,condmap,fii,fig,fgi,fgg);
  // create a new block matrix out of the 4 blocks
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockmat = LINALG::BlockMatrix2x2(*fii,*fig,*fgi,*fgg);

  if ( blockmat == Teuchos::null )
    dserror("Creation of fluid-fluid block matrix failed.");

  return blockmat;
}

Teuchos::RCP<const Epetra_Map> FLD::XFluidFluidNew::PressureRowMap()
{
  return state_->xffluidvelpressplitter_->CondMap();
}

Teuchos::RCP<const Epetra_Map> FLD::XFluidFluidNew::VelocityRowMap()
{
  return state_->xffluidvelpressplitter_->OtherMap();
}

void FLD::XFluidFluidNew::CreateState()
{
  // new cut for this time step
  state_ = this->GetNewState();
  // map of background-fluid's standard and enriched node-ids and
  // their dof-gids for new cut
  xfluidfluid_timeint_->SaveBgNodeMapsAndCreateNew(state_->Wizard());
}

Teuchos::RCP<FLD::XFluidFluidState> FLD::XFluidFluidNew::GetNewState()
{
  if (alefluid_)
  {
    LINALG::Export(*dispnp_,*idispnp_);
  }

  Teuchos::RCP<FLD::XFluidFluidState> state = state_creator_->Create(
    xdiscret_,
    embedded_fluid_->Discretization(),
    Teuchos::null, //!< col vector holding background ALE displacements for backdis
    solver_->Params(),
    step_,
    time_);

  // increment vector for merged background & embedded fluid
  // (not the classical Newton increment but the difference to
  // the value at the last time step)
  stepinc_ = LINALG::CreateVector(*state->xffluiddofrowmap_,true);

  // build a merged map from fluid-fluid dbc-maps
  state->CreateMergedDBCMapExtractor(embedded_fluid_->GetDBCMapExtractor());

  // create object for edgebased stabilization
  if (eval_eos_)
    edgestab_ =  Teuchos::rcp(new XFEM::XFEM_EdgeStab(state->Wizard(), discret_, include_inner_));

  return state;
}

void FLD::XFluidFluidNew::PrepareEmbeddedDistribution()
{
  DRT::UTILS::ConditionSelector conds(*embedded_fluid_->Discretization(), "XFEMSurfFluidFluid");
  std::vector<int> embnode_outer;
  std::vector<int> embele_outer;

  // select all outer embedded nodes and elements
  for (int  inode=0; inode<embedded_fluid_->Discretization()->NumMyRowNodes(); inode++)
  {
   DRT::Node* embnode = embedded_fluid_->Discretization()->lRowNode(inode);
   if (conds.ContainsNode(embnode->Id()))
   {
     embnode_outer.push_back(embnode->Id());
     const size_t numele = embnode->NumElement();
     // get list of adjacent elements of this node
     DRT::Element** adjeles = embnode->Elements();
     int mypid =embedded_fluid_->Discretization()->Comm().MyPID();

     for (size_t j=0; j<numele; ++j)
     {
       DRT::Element* adjele = adjeles[j];

       // if the element belongs to this processor insert it
       if (adjele->Owner() == mypid)
         embele_outer.push_back(adjele->Id());

     }
   }
  }

  // embnode_outer and embele_outer on all processors
  std::vector<int> embnode_outer_all;
  std::vector<int> embele_outer_all;

  // information how many processors work at all
  std::vector<int> allproc(embedded_fluid_->Discretization()->Comm().NumProc());

  // in case of n processors allproc becomes a vector with entries (0,1,...,n-1)
  for (int i=0; i<embedded_fluid_->Discretization()->Comm().NumProc(); ++i) allproc[i] = i;

  // gathers information of all processors
  LINALG::Gather<int>(embnode_outer,embnode_outer_all,(int)embedded_fluid_->Discretization()->Comm().NumProc(),&allproc[0],embedded_fluid_->Discretization()->Comm());
  LINALG::Gather<int>(embele_outer,embele_outer_all,(int)embedded_fluid_->Discretization()->Comm().NumProc(),&allproc[0],embedded_fluid_->Discretization()->Comm());

  // combine the colmap of every processor with the embedded outer map
  for (int  cnode=0; cnode<embedded_fluid_->Discretization()->NumMyColNodes(); cnode++)
  {
   DRT::Node* embnode = embedded_fluid_->Discretization()->lColNode(cnode);
   embnode_outer_all.push_back(embnode->Id());
  }
  for (int  cele=0; cele<embedded_fluid_->Discretization()->NumMyColElements(); cele++)
  {
   DRT::Element* ele = embedded_fluid_->Discretization()->lColElement(cele);
   embele_outer_all.push_back(ele->Id());
  }

  // create node and element distribution of outer layer elements and nodes
  // of the embedded discretization ghosted on all processors
  Teuchos::RCP<const Epetra_Map> embnodeoutermap = Teuchos::rcp(new Epetra_Map(-1, embnode_outer_all.size(), &embnode_outer_all[0], 0, embedded_fluid_->Discretization()->Comm()));
  Teuchos::RCP<const Epetra_Map> embeleoutermap = Teuchos::rcp(new Epetra_Map(-1, embele_outer_all.size(), &embele_outer_all[0], 0, embedded_fluid_->Discretization()->Comm()));

  const Epetra_Map embnodeoutercolmap = *LINALG::AllreduceOverlappingEMap(*embnodeoutermap);
  const Epetra_Map embelemoutercolmap = *LINALG::AllreduceOverlappingEMap(*embeleoutermap);

  // redistribute nodes and elements to column (ghost) map
  embedded_fluid_->Discretization()->ExportColumnNodes(embnodeoutercolmap);
  embedded_fluid_->Discretization()->ExportColumnElements(embelemoutercolmap);

  embedded_fluid_->Discretization()->FillComplete();
}

void FLD::XFluidFluidNew::CreateBoundaryEmbeddedMap()
{
  // Todo: Shift to XFEM::MeshCouplingFluidFluid and re-implement!

  // map to store the local id (value) of the boundary element (key) w.r.t to the embedded element
  std::map<int,int> boundary_emb_face_lid_map;
  std::map<int,int> boundary_emb_gid_map;
  // fill boundary_embedded_mapdmap between boundary element id and its corresponding embedded element id
  for (int iele=0; iele< boundarydis_->NumMyColElements(); ++iele)
  {
    // boundary element and its nodes
    DRT::Element* bele = boundarydis_->lColElement(iele);
    const int * belenodeIds = bele->NodeIds();

    bool bele_found = false;

    // ask all conditioned embedded elements for this boundary element
    for(int it=0; it< embedded_fluid_->Discretization()->NumMyColElements(); ++it)
    {
      DRT::Element* ele = embedded_fluid_->Discretization()->lColElement(it);
      const int * elenodeIds = ele->NodeIds();

      // get the surface-element map for the embedded element
      std::vector<std::vector<int> > face_node_map = DRT::UTILS::getEleNodeNumberingFaces(ele->Shape());

      // loop the faces of the element
      for(int f=0; f< ele->NumFace(); f++)
      {
        // assume the element has been found
        bele_found = true;

        const int face_numnode = face_node_map[f].size();

        if(bele->NumNode() != face_numnode) continue; // this face cannot be the right one

        // check all nodes of the boundary element
        for(int inode=0; inode<bele->NumNode();  ++inode)
        {
          // boundary node
          const int belenodeId = belenodeIds[inode];

          bool node_found = false;
          for (int fnode=0; fnode<face_numnode; ++fnode)
          {
            const int facenodeId = elenodeIds[face_node_map[f][fnode]];

            if(facenodeId == belenodeId)
            {
              // nodes are the same
              node_found = true;
              break;
            }
          } // loop nodes of element's face
          if(node_found==false) // this node is not contained in this face
          {
            bele_found = false; // element not the right one, if at least one boundary node is not found
            break; // node not found
          }
        } // loop nodes of boundary element

        if(bele_found==true)
        {
          boundary_emb_gid_map.insert(std::pair<int,int>(bele->Id(),ele->Id()));
          boundary_emb_face_lid_map.insert(std::pair<int,int>(bele->Id(),f));
          break;
        }
      } // loop element faces
      if(bele_found) break; // do not continue the search

    }

    if(bele_found == false) dserror("corresponding embele for boundary element with boundary id %i not found on proc %i ! Please ghost corresponding embedded elements on all procs!", bele->Id(), myrank_);
  }

  // update the estimate of the maximal eigenvalues in the parameter list to access on element level
  DRT::ELEMENTS::FluidEleParameterXFEM::Instance()->Set_boundary_emb_face_lid_map(boundary_emb_face_lid_map);
  DRT::ELEMENTS::FluidEleParameterXFEM::Instance()->Set_boundary_emb_gid_map(boundary_emb_gid_map);
}

void FLD::XFluidFluidNew::CreateEmbeddedBoundaryDiscretization()
{
  // Todo: Shift to XFEM::MeshCouplingFluidFluid and re-implement!

  // generate an empty boundary discretization
  embboundarydis_ = Teuchos::rcp(new DRT::Discretization(std::string("boundary discretization"),
      Teuchos::rcp(embedded_fluid_->Discretization()->Comm().Clone())));

  std::vector<DRT::Condition*> xfemcnd;
  embedded_fluid_->Discretization()->GetCondition("XFEMSurfFluidFluid",xfemcnd);

  // make the condition known to the boundary discretization
  for (unsigned cond=0; cond<xfemcnd.size(); ++cond)
  {
   // We use the same nodal ids and therefore we can just copy the conditions.
    embboundarydis_->SetCondition("XFEMSurfFluidFluid",Teuchos::rcp(new DRT::Condition(*xfemcnd[cond])));
  }

  // get the set of ids of all xfem nodes
  std::set<int> xfemnodeset;//(xfemcnd.size()*(xfemcnd[0]->Nodes()->size()));
  {
    for (unsigned cond=0; cond<xfemcnd.size(); ++cond)
    {
      // conditioned node ids
      const std::vector<int>* nodeids_cnd = xfemcnd[cond]->Nodes();
      for (std::vector<int>::const_iterator c = nodeids_cnd->begin();
           c != nodeids_cnd->end(); ++c)
        xfemnodeset.insert(*c);
    }
  }

  // determine sets of nodes next to xfem nodes
  std::set<int> adjacent_row;
  std::set<int> adjacent_col;

  // loop all column elements and label all row nodes next to a xfem node
  for (int i=0; i<embedded_fluid_->Discretization()->NumMyColElements(); ++i)
  {
    DRT::Element* actele = embedded_fluid_->Discretization()->lColElement(i);

    // get the node ids of this element
    const int  numnode = actele->NumNode();
    const int* nodeids = actele->NodeIds();

    bool found=false;

    // loop the element's nodes, check if a xfem condition is active
    for (int n=0; n<numnode; ++n)
    {
      const int node_gid(nodeids[n]);
      std::set<int>::iterator curr = xfemnodeset.find(node_gid);
      found = (curr!=xfemnodeset.end());
      if (found) break;
    }

    if (!found) continue;

    // if at least one of the element's nodes holds a xfem condition,
    // add all node gids to the adjecent node sets
    for (int n=0; n<numnode; ++n)
    {
      const int node_gid(nodeids[n]);
      // yes, we have a xfem condition:
      // node stored on this proc? add to the set of row nodes!
      if (embedded_fluid_->Discretization()->NodeRowMap()->LID(node_gid) > -1)
        adjacent_row.insert(node_gid);

      // always add to set of col nodes
      adjacent_col.insert(node_gid);
    }

    // add the element to the discretization
    if (embedded_fluid_->Discretization()->ElementRowMap()->LID(actele->Id()) > -1)
    {
      Teuchos::RCP<DRT::Element> bndele =Teuchos::rcp(actele->Clone());
      embboundarydis_->AddElement(bndele);
    }
  } // end loop over column elements

  // all row nodes next to a xfem node are now added to the embedded boundary discretization
  for (std::set<int>::iterator id=adjacent_row.begin();
       id!=adjacent_row.end(); ++id)
  {
    DRT::Node* actnode=embedded_fluid_->Discretization()->gNode(*id);
    Teuchos::RCP<DRT::Node> bndnode =Teuchos::rcp(actnode->Clone());
    embboundarydis_->AddNode(bndnode);
  }

  // build nodal row & col maps to redistribute the discretization
  Teuchos::RCP<Epetra_Map> newrownodemap;
  Teuchos::RCP<Epetra_Map> newcolnodemap;

  {
    // copy row/col node gids to std::vector
    // (expected by Epetra_Map ctor)
    std::vector<int> rownodes(adjacent_row.begin(),adjacent_row.end());
    // build noderowmap for new distribution of nodes
    newrownodemap = Teuchos::rcp(new Epetra_Map(-1,
                                                rownodes.size(),
                                                &rownodes[0],
                                                0,
                                                embboundarydis_->Comm()));

    std::vector<int> colnodes(adjacent_col.begin(),adjacent_col.end());

    // build nodecolmap for new distribution of nodes
    newcolnodemap = Teuchos::rcp(new Epetra_Map(-1,
                                                colnodes.size(),
                                                &colnodes[0],
                                                0,
                                                embboundarydis_->Comm()));

    embboundarydis_->Redistribute(*newrownodemap,*newcolnodemap,false,false,false);

    // make boundary discretization have the same dofs as the embedded fluid
    Teuchos::RCP<DRT::DofSet> newdofset = Teuchos::rcp(new DRT::TransparentIndependentDofSet(soliddis_,true));
    embboundarydis_->ReplaceDofSet(newdofset); // do not call this with true (no replacement in static dofsets intended)
    embboundarydis_->FillComplete();
  }
}

void FLD::XFluidFluidNew::AssembleMatAndRHS(
    int itnum                           ///< iteration number
)
{
  // Todo: doesn't work yet...
  XFluid::AssembleMatAndRHS(itnum);
  //embedded_fluid_->AssembleMatAndRHS();
}

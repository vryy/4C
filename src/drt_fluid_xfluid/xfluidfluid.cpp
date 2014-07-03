/*!----------------------------------------------------------------------
\file xfluidfluid.cpp
\brief Control routine for fluid-fluid (in)stationary solvers with XFEM

<pre>
Maintainer:  Shadan Shahmiri
             shahmiri@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15265
</pre>

*----------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>

#include "../drt_fluid_ele/fluid_ele_action.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_discret_faces.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_dofset_transparent_independent.H"
#include "../drt_lib/drt_dofset.H"
#include "../drt_lib/drt_dofset_independent.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_parallel.H"

#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include "../drt_geometry/geo_intersection.H"

#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_sidehandle.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_cut/cut_integrationcell.H"
#include "../drt_cut/cut_point.H"
#include "../drt_cut/cut_meshintersection.H"
#include "../drt_cut/cut_position.H"

#include "../drt_fluid_ele/fluid_ele.H"

#include "../drt_fluid_ele/fluid_ele_interface.H"
#include "../drt_fluid_ele/fluid_ele_factory.H"

#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_fluid_ele/fluid_ele_parameter_xfem.H"

#include "../drt_io/io.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_xfem/xfem_edgestab.H"
#include "../drt_xfem/xfem_neumann.H"
#include "../drt_xfem/xfem_fluiddofset.H"
#include "../drt_xfem/xfem_fluidwizard.H"
#include "../drt_xfem/xfluidfluid_timeInt.H"

#include "../drt_combust/combust_utils_time_integration.H"
#include "xfluidfluidresulttest.H"

#include "../drt_fluid/fluid_utils.H"

#include "xfluid_defines.H"

#include "xfluidfluid.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matpar_bundle.H"

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
FLD::XFluidFluid::XFluidFluidState::XFluidFluidState( XFluidFluid & xfluid, Epetra_Vector & idispcol )
  : xfluid_( xfluid ),
    wizard_( Teuchos::rcp( new XFEM::FluidWizard(*xfluid.bgdis_, *xfluid.boundarydis_ )))
{
  // do the (parallel!) cut for the 0 timestep and find the fluid dofset
  wizard_->Cut( false,                                 // include_inner
                idispcol,                              // interface displacements
                xfluid_.VolumeCellGaussPointBy_,       // how to create volume cell Gauss points?
                xfluid_.BoundCellGaussPointBy_,        // how to create boundary cell Gauss points?
                true,                                  // use parallel cut framework
                xfluid_.gmsh_cut_out_,                 // gmsh output for cut library
                true                                   // find point positions
                );

  // estimate (and reserve) the maximum number of dof
  // (nodes)*4*(max. number of xfem-dofsets p. node)
  int maxNumMyReservedDofs = xfluid.bgdis_->NumGlobalNodes()*(xfluid.maxnumdofsets_)*4;

  // ask fluid wizard for new background fluid dofset
  dofset_ = wizard_->DofSet(maxNumMyReservedDofs);

  // get minimal GID
  const int restart = DRT::Problem::Instance()->Restart();
  if ((xfluid.step_ < 1) or restart)
  {
    xfluid.minnumdofsets_ = xfluid.bgdis_->DofRowMap()->MinAllGID();
  }
  // set the minimal GID of xfem dis
  dofset_->MinGID(xfluid.minnumdofsets_);
  // set new dofset
  xfluid_.bgdis_->ReplaceDofSet( dofset_, true);
  xfluid_.bgdis_->FillComplete();

  // print all dofsets
  //---ALE---|---BG---|---Emb---|---Str---|
  xfluid_.bgdis_->GetDofSetProxy()->PrintAllDofsets(xfluid_.bgdis_->Comm());

  // ---------------------------------------------------------------------------
  // TODO:
  // recompute nullspace based on new number of dofs per node
  // REMARK: this has to be done after replacing the discret_' dofset (via xfluid.discret_->ReplaceDofSet)
  // we a build a nullspace for both discrretizations and set it into ML-parameter list.
  if (xfluid_.solver_->Params().isSublist("ML Parameters"))
  {
	  // build a nullspace for bg-dis
    xfluid_.bgdis_->ComputeNullSpaceIfNecessary(xfluid_.solver_->Params(),true);
	  Teuchos::ParameterList& ml_params = (xfluid_.solver_->Params()).sublist("ML Parameters");
	  Teuchos::RCP<std::vector<double> > nullspace = ml_params.get<Teuchos::RCP<std::vector<double> > >("nullspace");

	  // build a nullspace for emb-dis
	  xfluid_.embdis_->ComputeNullSpaceIfNecessary(xfluid_.solver_->Params(),true);
	  Teuchos::RCP<std::vector<double> > nullspace_embdis = ml_params.get<Teuchos::RCP<std::vector<double> > >("nullspace");

	  // build a full nullspace
	  nullspace->insert(nullspace->end(), nullspace_embdis->begin(), nullspace_embdis->end());
	  ml_params.set<Teuchos::RCP<std::vector<double> > >("nullspace",nullspace);
	  Teuchos::RCP<std::vector<double> > nullspaceend = ml_params.get<Teuchos::RCP<std::vector<double> > >("nullspace");
  }
  // ---------------------------------------------------------------------------

  // split dof row map into velocity and pressure dof
  // (btw this has nothing to do with monolithic fluidsplit FSI)
  FLD::UTILS::SetupFluidSplit(*xfluid.bgdis_,xfluid.numdim_, 1,velpressplitter_);

  // get the background fluid dof-rowmap
  fluiddofrowmap_ = Teuchos::rcp(xfluid.bgdis_->DofRowMap(), false);

  // create an EpetraFECrs matrix that does communication for non-local rows and columns
  // * this enables to do the evaluate loop over just row elements instead of col elements
  // * time consuming assemble for cut elements is done only once on a unique row processor
  // REMARK: call the SparseMatrix: * explicitdirichlet = true (is used in ApplyDirichlet, false uses Epetra memory based operations
  //                                                            that are not ensured to be always compatible with baci)
  //                                * savegraph = true, store the matrix graph (pattern for non-zero entries)
  //                                * with FE_MATRIX flag

  // REMARK: The max. number of entries to be reserved for the sparse matrix is set together as follows:
  // for a standard hex8 element one would obtain 27(nodes)*4(dof)=108.
  // Additionally, for edge-based stabilization, one obtains additional 6(faces)*9(nodes)*1(dof)=54 entries per row.
  // (Why 1(dof)? --> edge-based stabilization couples v_x->u_x, v_y->u_y, v_z->u_z, q->p)

  // setup of background fluid system matrix
  // savegraph-option set to true, as the matrix graph will not change throughout the lifetime of
  // this state-object (cut happens in ctor)
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*fluiddofrowmap_,108+54,true,true,LINALG::SparseMatrix::FE_MATRIX));

  // ---------------------------------------------------------------------------
  // background fluid vectors
  // ---------------------------------------------------------------------------

  // background fluid velocity/pressure at time n+1, n and n-1
  velnp_ = LINALG::CreateVector(*fluiddofrowmap_,true);
  veln_  = LINALG::CreateVector(*fluiddofrowmap_,true);
  velnm_ = LINALG::CreateVector(*fluiddofrowmap_,true);

  // acceleration/(scalar time derivative) at time n+1 and n
  accnp_ = LINALG::CreateVector(*fluiddofrowmap_,true);
  accn_  = LINALG::CreateVector(*fluiddofrowmap_,true);

  // velocity/pressure at time n+alpha_F
  velaf_ = LINALG::CreateVector(*fluiddofrowmap_,true);

  // acceleration/(scalar time derivative) at time n+alpha_M/(n+alpha_M/n)
  accam_ = LINALG::CreateVector(*fluiddofrowmap_,true);

  // scalar at time n+alpha_F/n+1 and n+alpha_M/n
  // (only required for low-Mach-number case)
  scaaf_ = LINALG::CreateVector(*fluiddofrowmap_,true);
  scaam_ = LINALG::CreateVector(*fluiddofrowmap_,true);

  // history vector
  hist_ = LINALG::CreateVector(*fluiddofrowmap_,true);

  // the vector containing body and surface forces
  neumann_loads_= LINALG::CreateVector(*fluiddofrowmap_,true);

  // rhs: standard (stabilized) residual vector (rhs for the incremental form)
  residual_     = LINALG::CreateVector(*fluiddofrowmap_,true);
  trueresidual_ = LINALG::CreateVector(*fluiddofrowmap_,true);

  // right hand side vector for linearised solution;
  rhs_ = LINALG::CreateVector(*fluiddofrowmap_,true);

  // Nonlinear iteration increment vector
  incvel_ = LINALG::CreateVector(*fluiddofrowmap_,true);

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_   = LINALG::CreateVector(*fluiddofrowmap_,true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());

  {
    Teuchos::ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time",xfluid.time_);
    xfluid.bgdis_->EvaluateDirichlet(eleparams, zeros_, Teuchos::null, Teuchos::null,
                                     Teuchos::null, dbcmaps_);

    zeros_->PutScalar(0.0); // just in case of change
  }

  // ---------------------------------------------------------------------------
  // combined dof-maps & map extractors for merged background & embedded fluid
  // ---------------------------------------------------------------------------

  // merge the background and embedded fluid dof-maps
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  Teuchos::RCP<const Epetra_Map> alefluiddofrowmap = Teuchos::rcp(xfluid.embdis_->DofRowMap(),false);
  maps.push_back(fluiddofrowmap_);

  if ((fluiddofrowmap_->MaxAllGID()) > (alefluiddofrowmap->MinAllGID()))
		  dserror("Maximum number of reserved dofs is not enough! Change MAX_NUM_DOFSETS in you dat-file.");

  maps.push_back(alefluiddofrowmap);
  fluidfluiddofrowmap_ = LINALG::MultiMapExtractor::MergeMaps(maps);

  // setup the embedded fluid/background fluid map extractor
  fluidfluidsplitter_ = Teuchos::rcp(new FLD::UTILS::FluidXFluidMapExtractor());
  fluidfluidsplitter_->Setup(*fluidfluiddofrowmap_,alefluiddofrowmap,fluiddofrowmap_);

  FLD::UTILS::SetupFluidFluidVelPresSplit(*xfluid.bgdis_,xfluid.numdim_,*xfluid.embdis_,fluidfluidvelpressplitter_,
                                          fluidfluiddofrowmap_);

  // ---------------------------------------------------------------------------
  // matrices & vectors for merged background & embedded fluid
  // ---------------------------------------------------------------------------

  // the combined fluid system matrix is not of FECrs-type - it is solely composed out of
  // fully assembled submatrices
  fluidfluidsysmat_   = Teuchos::rcp(new LINALG::SparseMatrix(*fluidfluiddofrowmap_,108,false,true));

  fluidfluidresidual_ = LINALG::CreateVector(*fluidfluiddofrowmap_,true);
  fluidfluidincvel_   = LINALG::CreateVector(*fluidfluiddofrowmap_,true);
  fluidfluidvelnp_    = LINALG::CreateVector(*fluidfluiddofrowmap_,true);
  fluidfluidveln_     = LINALG::CreateVector(*fluidfluiddofrowmap_,true);
  fluidfluidzeros_    = LINALG::CreateVector(*fluidfluiddofrowmap_,true);

  // increment vector for merged background & embedded fluid
  // (not the classical Newton increment but the difference to
  // the value at the last time step)
  stepinc_ = LINALG::CreateVector(*fluidfluiddofrowmap_,true);

  // create object for edgebased stabilization
  if(xfluid_.edge_based_ or xfluid_.ghost_penalty_ )
    edgestab_ =  Teuchos::rcp(new XFEM::XFEM_EdgeStab(wizard_, xfluid.bgdis_));

  // allocate memory for fluid-fluid coupling matrices
  PrepareCouplingMatrices(true);

  // build merged fluid dirichlet maps
  CreateFluidFluidDBCMaps();
}

/*------------------------------------------------------------------------------------------------*
 * estimate the scaling from the trace inequality via solving local eigenvalue problems for Nitsche's method
 *------------------------------------------------------------------------------------------------*/
void FLD::XFluidFluid::EstimateNitscheTraceMaxEigenvalue()
{
	Teuchos::ParameterList params;

  // set action for elements
  params.set<int>("action",FLD::estimate_Nitsche_trace_maxeigenvalue_);

  Teuchos::RCP<Epetra_Vector> embboundarydispnp = LINALG::CreateVector(*embboundarydis_->DofRowMap(),true);
  LINALG::Export(*aledispnp_,*embboundarydispnp);

  embboundarydis_->SetState("dispnp", embboundarydispnp);

  std::map<int,double>  boundaryeleid_to_maxeigenvalue;

  /// map of embedded element ID to the value of it's Nitsche parameter
  Teuchos::RCP<std::map<int,double> > trace_estimate_max_eigenvalue_map;
  trace_estimate_max_eigenvalue_map =  Teuchos::rcp(new std::map<int,double> (boundaryeleid_to_maxeigenvalue));
  params.set<Teuchos::RCP<std::map<int,double > > >("trace_estimate_max_eigenvalue_map", trace_estimate_max_eigenvalue_map);

  Teuchos::RCP<LINALG::SparseOperator> systemmatrixA;
  Teuchos::RCP<LINALG::SparseOperator> systemmatrixB;

  // Evaluate the general eigen value problem Ax = lambda Bx for local for the elements of embboundarydis
  embboundarydis_->EvaluateCondition(params,
                                     systemmatrixA,
                                     systemmatrixB,
                                     Teuchos::null,
                                     Teuchos::null,
                                     Teuchos::null,
                                     "XFEMCoupling");

  // gather the information form all processors
  Teuchos::RCP<std::map<int,double> > tmp_map = params.get<Teuchos::RCP<std::map<int,double > > >("trace_estimate_max_eigenvalue_map");

  // information how many processors work at all
  std::vector<int> allproc(embboundarydis_->Comm().NumProc());

  // in case of n processors allproc becomes a vector with entries (0,1,...,n-1)
  for (int i=0; i<embboundarydis_->Comm().NumProc(); ++i) allproc[i] = i;

  // gather the information from all procs
  LINALG::Gather<int,double>(*tmp_map,*trace_estimate_max_eigenvalue_map,(int)bgdis_->Comm().NumProc(),&allproc[0], bgdis_->Comm());

  // update the estimate of the maximal eigenvalues in the parameter list to access on element level
  DRT::ELEMENTS::FluidEleParameterXFEM::Instance()->Update_TraceEstimate_MaxEigenvalue(trace_estimate_max_eigenvalue_map);
}

/*------------------------------------------------------------------------------------------------*
 | create embedded-boundary dis, which includes the fluid elements of the embedded discretization
 | which are contribute to the fluid-fluid interface (this is not a boundary discretization)
 *------------------------------------------------------------------------------------------------*/
void FLD::XFluidFluid::CreateEmbeddedBoundarydis()
{
  // generate an empty boundary discretization
	embboundarydis_ = Teuchos::rcp(new DRT::Discretization(std::string("boundary discretization"),
		  Teuchos::rcp(embdis_->Comm().Clone())));

  std::vector<DRT::Condition*> xfemcnd;
  embdis_->GetCondition("XFEMCoupling",xfemcnd);

  std::vector<int> allcnd_xfemnodeids;

  // make the condition known to the boundary discretization
  for (unsigned numcond=0;numcond<xfemcnd.size();++numcond)
  {
   // We use the same nodal ids and therefore we can just copy the conditions.
	  embboundarydis_->SetCondition("XFEMCoupling",Teuchos::rcp(new DRT::Condition(*xfemcnd[numcond])));
  }

  // get the set of ids of all xfem nodes
  std::set<int> xfemnodeset;
  {

    for (unsigned numcond=0;numcond<xfemcnd.size();++numcond)
    {
      const std::vector <int>* xfemnodeids = (*xfemcnd[numcond]).Nodes();

      allcnd_xfemnodeids.reserve( allcnd_xfemnodeids.size() + xfemnodeids->size());
      allcnd_xfemnodeids.insert ( allcnd_xfemnodeids.end(), xfemnodeids->begin(), xfemnodeids->end());
    }

    for(std::vector<int>::iterator id =allcnd_xfemnodeids.begin();
        id!=allcnd_xfemnodeids.end();
        ++id)
    {
      xfemnodeset.insert(*id);
    }
  }

  // determine sets of nodes next to xfem nodes
  std::set<int> adjacent_row;
  std::set<int> adjacent_col;

  // loop all column elements and label all row nodes next to a xfem node
  for (int i=0; i<embdis_->NumMyColElements(); ++i)
  {
    DRT::Element* actele = embdis_->lColElement(i);

    // get the node ids of this element
    const int  numnode = actele->NumNode();
    const int* nodeids = actele->NodeIds();

    bool found=false;

    // loop nodeids, check if a xfem condition is active
    for(int rr=0;rr<numnode;++rr)
    {
      int gid=nodeids[rr];

      std::set<int>::iterator curr=xfemnodeset.find(gid);
      if(curr!=xfemnodeset.end())
      {
        found=true;
      }
    }

    // yes, we have a xfem condition
    if(found==true)
    {
      // loop nodeids
      for(int rr=0;rr<numnode;++rr)
      {
        int gid=nodeids[rr];

        if ((embdis_->NodeRowMap())->LID(gid)>-1)
        {
          adjacent_row.insert(gid);
        }
        adjacent_col.insert(gid);
      }
    }
  }

  // all row nodes next to a xfem nodes are now contained in the discretization
  for(std::set<int>::iterator id = adjacent_row.begin();
      id!=adjacent_row.end();
      ++id)
  {
    DRT::Node* actnode=embdis_->gNode(*id);

    Teuchos::RCP<DRT::Node> bndnode =Teuchos::rcp(actnode->Clone());

    embboundarydis_->AddNode(bndnode);
  }

  // loop all row elements and add all elements with a xfem node
  for (int i=0; i<embdis_->NumMyRowElements(); ++i)
  {
    DRT::Element* actele = embdis_->lRowElement(i);

    // get the node ids of this element
    const int  numnode = actele->NumNode();
    const int* nodeids = actele->NodeIds();

    bool found=false;

    // loop nodeids, check if a MHD condition is active
    for(int rr=0;rr<numnode;++rr)
    {
      int gid=nodeids[rr];

      std::set<int>::iterator curr=xfemnodeset.find(gid);
      if(curr!=xfemnodeset.end())
      {
        found=true;
      }
    }

    // yes, we have a xfem condition
    if(found==true)
    {
      Teuchos::RCP<DRT::Element> bndele =Teuchos::rcp(actele->Clone());

      embboundarydis_->AddElement(bndele);
    }
  }

  // the discretization needs a full NodeRowMap and a NodeColMap
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
                                                embboundarydis_->Comm()));

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
                                                embboundarydis_->Comm()));

    embboundarydis_->Redistribute(*newrownodemap,*newcolnodemap,false,false,false);

    // make embboundarydis have the same dofs as embedded-dis
    Teuchos::RCP<DRT::DofSet> newdofset=Teuchos::rcp(new DRT::TransparentIndependentDofSet(embdis_,false,Teuchos::null));
    embboundarydis_->ReplaceDofSet(newdofset); // do not call this with true!!
    embboundarydis_->FillComplete();
  }
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void FLD::XFluidFluid::XFluidFluidState::EvaluateFluidFluid(const Teuchos::RCP<DRT::Discretization>    & bgdis,
                                                            const Teuchos::RCP<DRT::Discretization>    & cutdiscret,
                                                            const Teuchos::RCP<DRT::Discretization>    & embdis)
{
  TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluidFluid::XFluidFluidState::EvaluateFluidFluid" );

  // set background fluid matrix to zero
  sysmat_->Zero();

  // set embedded fluid matrix to zero
  xfluid_.alesysmat_->Zero();

  // add Neumann loads
  residual_->Update(1.0,*neumann_loads_,0.0);
  xfluid_.aleresidual_->Update(1.0,*xfluid_.aleneumann_loads_,0.0);

  // create an column residual vector (background fluid residual)
  // for assembly over row elements, that has to be communicated at the end
  Teuchos::RCP<Epetra_Vector> residual_col = LINALG::CreateVector(*(bgdis->DofColMap()),true);
  // create an column fluid-fluid coupling rhs-vector for assembly over row elements,
  // that has to be communicated at the end
  Teuchos::RCP<Epetra_Vector> rhC_ui_col;

  // ---------------------------------------------------------------------------

  // set general vector values needed by background elements
  bgdis->ClearState();
  bgdis->SetState("hist" ,hist_ );
  bgdis->SetState("veln" ,veln_ );
  bgdis->SetState("accam",accam_);
  bgdis->SetState("scaaf",scaaf_);
  bgdis->SetState("scaam",scaam_);

  // set general vector values needed by embedded elements
  embdis->ClearState();
  embdis->SetState("hist" ,xfluid_.alehist_ );
  embdis->SetState("veln" ,xfluid_.aleveln_ );
  embdis->SetState("accam",xfluid_.aleaccam_);
  embdis->SetState("scaaf",xfluid_.alescaaf_);
  embdis->SetState("scaam",xfluid_.alescaam_);

  // if we have a moving embedded fluid or embedded-sided coupling,
  // set the embedded fluid displacement
  if (xfluid_.alefluid_ or xfluid_.coupling_approach_ == CouplingNitsche_EmbFluid or
                           xfluid_.coupling_approach_ == CouplingNitsche_TwoSided)
    embdis->SetState("dispnp", xfluid_.aledispnp_);

  if (xfluid_.alefluid_)
    embdis->SetState("gridv", xfluid_.gridv_);

  // set general vector values of boundarydis needed by elements
  LINALG::Export(*(xfluid_.alevelnp_),*(xfluid_.ivelnp_));

  cutdiscret->SetState("ivelnp",xfluid_.ivelnp_);

  // set interface dispnp needed for the elements
  if (xfluid_.alefluid_)
    LINALG::Export(*(xfluid_.aledispnp_),*(xfluid_.idispnp_));

  cutdiscret->SetState("idispnp",xfluid_.idispnp_);

  // set scheme-specific element parameters and vector values
  if (xfluid_.timealgo_==INPAR::FLUID::timeint_afgenalpha)
  {
    dserror("no genalpha for fluid-fluid!!");
    bgdis->SetState("velaf",velaf_);
    embdis->SetState("velaf",xfluid_.alevelaf_);
  }
  else
  {
    bgdis->SetState("velaf",velnp_);
    embdis->SetState("velaf",xfluid_.alevelnp_);
  }

  DRT::AssembleStrategy strategy(0, 0, sysmat_,Teuchos::null,residual_col,Teuchos::null,Teuchos::null);
  DRT::AssembleStrategy alestrategy(0, 0, xfluid_.alesysmat_,xfluid_.shapederivatives_, xfluid_.aleresidual_,Teuchos::null,Teuchos::null);

  PrepareCouplingMatrices();

  if (xfluid_.coupling_approach_ == CouplingMHCS_XFluid or
      xfluid_.coupling_approach_ == CouplingNitsche_XFluid or
      xfluid_.coupling_approach_ == CouplingMHVS_XFluid)
  {
    rhC_ui_col= LINALG::CreateVector(*cutdiscret->DofColMap(),true);
  }
  else if (xfluid_.coupling_approach_ == CouplingNitsche_EmbFluid or xfluid_.coupling_approach_ == CouplingNitsche_TwoSided)
  {
    rhC_ui_col= LINALG::CreateVector(*embdis->DofColMap(),true);
  }

  DRT::Element::LocationArray la( 1 );
  DRT::Element::LocationArray alela( 1 );
  DRT::Element::LocationArray ila ( 1 );

  // dummy
  Teuchos::ParameterList eleparams;

  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // loop over row elements
  // ---------------------------------------------------------------------------

  const int numrowele = bgdis->NumMyRowElements();

  // REMARK: in this XFEM framework the whole evaluate routine uses only row elements and
  // assembles into EpetraFECrs matrix
  // this is baci-unusual but more efficient in all XFEM applications
  for (int i=0; i<numrowele; ++i)
  {
    DRT::Element * actele = bgdis->lRowElement(i);
    Teuchos::RCP<MAT::Material> mat = actele->Material();

    DRT::ELEMENTS::Fluid * ele = dynamic_cast<DRT::ELEMENTS::Fluid *>( actele );
    if ( ele==NULL )
    {
      dserror("Failed to cast element %d from background discretization to fluid element.", actele->Id());
    }

    DRT::ELEMENTS::FluidEleInterface * impl = DRT::ELEMENTS::FluidFactory::ProvideImplXFEM(actele->Shape(), "xfem");

    // wizard returns NULL, if element is not intersected
    GEO::CUT::ElementHandle * e = wizard_->GetElement( actele );

    // evaluate an intersected background fluid element
    if ( e!=NULL )
    {

#ifdef DOFSETS_NEW
      // sets of volume-cells associated with current element
      std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
      // sets of nodal dofsets associated with current element
      std::vector< std::vector<int> > nds_sets;
      // sets of integration points associated with current element
      std::vector<std::vector< DRT::UTILS::GaussIntegration > > intpoints_sets;

      // the volume-cell set at position i in cell_sets is associated with
      // the nodal dofset vector at position i in nds_sets
      // and with the set of gauss points at position i in intpoints_sets
      e->GetCellSets_DofSets_GaussPoints( cell_sets, nds_sets, intpoints_sets, xfluid_.VolumeCellGaussPointBy_ );

      if ( cell_sets.size() != intpoints_sets.size() )
        dserror("Non-matching number of volume-cell sets (%d) and integration point sets (%d).", cell_sets.size(), intpoints_sets.size());
      if ( cell_sets.size() != nds_sets.size() )
        dserror("Non-matching number of volume-cell sets (%d) and sets of nodal dofsets (%d).", cell_sets.size(), nds_sets.size());

      int set_counter = 0;

      // run through the volume-cell sets
      for( std::vector< GEO::CUT::plain_volumecell_set>::iterator s=cell_sets.begin();
           s!=cell_sets.end();
           s++ )
      {
        // for each side that is involved in the cut for this element, the coupling matrices Cuiu, Cuui and the rhs has to be built
        std::map<int, std::vector<Epetra_SerialDenseMatrix> > side_coupling;
        GEO::CUT::plain_volumecell_set & cells = *s;
        const std::vector<int> & nds = nds_sets[set_counter];

        // we have to assembly all volume cells of this set
        // for linear elements, there should be only one volumecell for each set
        // for quadratic elements, there are some volumecells with respect to subelements, that have to be assembled at once

        // get element location vector, dirichlet flags and ownerships
        actele->LocationVector(*bgdis,nds,la,false);

        // get dimension of element matrices and vectors
        // Reshapelement matrices and vectors and init to zero
        strategy.ClearElementStorage( la[0].Size(), la[0].Size() );

        {
          // ---------------------------------------------------------------------------
          // Evaluate domain integrals
          TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate cut domain" );

          // call the element evaluate method
          int err = impl->EvaluateXFEM( ele, *bgdis, la[0].lm_, eleparams, mat,
                                        strategy.Elematrix1(),
                                        strategy.Elematrix2(),
                                        strategy.Elevector1(),
                                        strategy.Elevector2(),
                                        strategy.Elevector3(),
                                        intpoints_sets[set_counter],
                                        cells);

          if (err)
            dserror("Proc %d: Element %d returned err=%d",bgdis->Comm().MyPID(),actele->Id(),err);
          // ---------------------------------------------------------------------------
        }

        // do cut interface condition
        // maps of sid and corresponding boundary cells (for quadratic elements: collected via volumecells of subelements)
        std::map<int, std::vector<GEO::CUT::BoundaryCell*> > bcells;
        std::map<int, std::vector<DRT::UTILS::GaussIntegration> > bintpoints;

        for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
        {
          GEO::CUT::VolumeCell * vc = *i;
          if ( vc->Position()==GEO::CUT::Point::outside )
          {
              vc->GetBoundaryCells( bcells );
          }
        }

        if ( bcells.size() > 0 )
        {
          TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate boundary" );

          // Attention: switch also the flag in fluid3_impl.cpp
#ifdef BOUNDARYCELL_TRANSFORMATION_OLD
          // original Axel's transformation
          e->BoundaryCellGaussPoints( wizard_->CutWizard().Mesh(), 0, bcells, bintpoints );
#else
          // new Benedikt's transformation
          e->BoundaryCellGaussPointsLin( wizard_->CutWizard().Mesh(), 0, bcells, bintpoints );
#endif


          std::vector<int> patchelementslm;
          std::vector<int> patchelementslmowner;

          // initialize the coupling matrices for each side and the current element
          for ( std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
                bc!=bcells.end(); ++bc )
          {
            int sid = bc->first; // all boundary cells within the current iterator belong to the same side
            DRT::Element * side = cutdiscret->gElement( sid );

            std::vector<int> patchlm;
            std::vector<int> patchlmowner;
            std::vector<int> patchlmstride;
            // for nitsche embedded and two-sided we couple with the whole embedded element not only with its side
            if (xfluid_.coupling_approach_ == CouplingMHCS_XFluid or
                xfluid_.coupling_approach_ == CouplingNitsche_XFluid or
                xfluid_.coupling_approach_ == CouplingMHVS_XFluid)
              side->LocationVector(*cutdiscret, patchlm, patchlmowner, patchlmstride);
            else if(xfluid_.coupling_approach_ == CouplingNitsche_EmbFluid or xfluid_.coupling_approach_ == CouplingNitsche_TwoSided)
            {
              // get the corresponding embedded element for nitsche
              // embedded and two-sided
              int emb_eid = xfluid_.boundary_emb_gid_map_.find(sid)->second;
              DRT::Element * emb_ele = embdis->gElement( emb_eid );
              emb_ele->LocationVector(*embdis, patchlm, patchlmowner, patchlmstride);
            }

            patchelementslm.reserve( patchelementslm.size() + patchlm.size());
            patchelementslm.insert(patchelementslm.end(), patchlm.begin(), patchlm.end());

            patchelementslmowner.reserve( patchelementslmowner.size() + patchlmowner.size());
            patchelementslmowner.insert( patchelementslmowner.end(), patchlmowner.begin(), patchlmowner.end());

            const size_t ndof_i = patchlm.size();     // sum over number of dofs of all sides
            const size_t ndof   = la[0].lm_.size();   // number of dofs for background element

            std::vector<Epetra_SerialDenseMatrix> & couplingmatrices = side_coupling[sid]; // the function inserts a new element with that key and returns a reference to its mapped value
            if ( couplingmatrices.size()!=0 )
              dserror("zero sized vector expected");

            couplingmatrices.resize(3);

            // no coupling for pressure in stress based method, but the coupling matrices include entries for pressure coupling
            couplingmatrices[0].Reshape(ndof_i,ndof);  //C_uiu
            couplingmatrices[1].Reshape(ndof,ndof_i);  //C_uui
            couplingmatrices[2].Reshape(ndof_i,1);     //rhC_ui

          }

          const size_t nui = patchelementslm.size();
          Epetra_SerialDenseMatrix  Cuiui(nui,nui);

          if (xfluid_.coupling_approach_ == CouplingMHCS_XFluid or xfluid_.coupling_approach_ == CouplingMHVS_XFluid)
            impl->ElementXfemInterfaceHybridLM(
                                              ele,
                                              *bgdis,
                                              la[0].lm_,
                                              intpoints_sets[set_counter],
                                              *cutdiscret,
                                              bcells,
                                              bintpoints,
                                              side_coupling,
                                              eleparams,
                                              strategy.Elematrix1(),
                                              strategy.Elevector1(),
                                              Cuiui,
                                              cells,
                                              true);
          else if (xfluid_.coupling_approach_ == CouplingNitsche_XFluid)
            impl->ElementXfemInterfaceNIT(    ele,
                                              *bgdis,
                                              la[0].lm_,
                                              *cutdiscret,
                                              bcells,
                                              bintpoints,
                                              side_coupling,
                                              eleparams,
                                              strategy.Elematrix1(),
                                              strategy.Elevector1(),
                                              Cuiui,
                                              cells,
                                              true);
          else if (xfluid_.coupling_approach_ == CouplingNitsche_EmbFluid or xfluid_.coupling_approach_ == CouplingNitsche_TwoSided)
            impl->ElementXfemInterfaceNIT2(   ele,
                                              *bgdis,
                                              la[0].lm_,
                                              *cutdiscret,
                                              bcells,
                                              bintpoints,
                                              side_coupling,
                                              eleparams,
                                              *embdis,
                                              xfluid_.boundary_emb_gid_map_,
                                              strategy.Elematrix1(),
                                              strategy.Elevector1(),
                                              Cuiui,
                                              cells);
          else
            dserror("The coupling method you have chosen is not (yet) implemented.");



          for ( std::map<int, std::vector<Epetra_SerialDenseMatrix> >::const_iterator sc=side_coupling.begin();
                sc!=side_coupling.end(); ++sc )
          {
            std::vector<Epetra_SerialDenseMatrix>  couplingmatrices = sc->second;

            int sid = sc->first;

            if ( cutdiscret->HaveGlobalElement(sid) )
            {
              std::vector<int> patchlm;
              std::vector<int> patchlmowner;
              std::vector<int> patchlmstride;
              if (xfluid_.coupling_approach_ == CouplingMHCS_XFluid or
                  xfluid_.coupling_approach_ == CouplingNitsche_XFluid or
                  xfluid_.coupling_approach_ == CouplingMHVS_XFluid)
              {
                DRT::Element * side = cutdiscret->gElement( sid );
                side->LocationVector(*cutdiscret, patchlm, patchlmowner, patchlmstride);
              }
              else if (xfluid_.coupling_approach_ == CouplingNitsche_EmbFluid or xfluid_.coupling_approach_ == CouplingNitsche_TwoSided)
              {
                int emb_eid = xfluid_.boundary_emb_gid_map_.find(sid)->second;
                DRT::Element * emb_ele = embdis->gElement( emb_eid );
                emb_ele->LocationVector(*embdis, patchlm, patchlmowner, patchlmstride);
              }

              // assemble Cuiu
              //create a dummy mypatchlmowner that assembles also non-local rows and communicates the required data
              std::vector<int> mypatchlmowner;
              for(size_t index=0; index<patchlmowner.size(); index++) mypatchlmowner.push_back(xfluid_.myrank_);

              Cuiu_->FEAssemble(-1, couplingmatrices[0],patchlm,mypatchlmowner,la[0].lm_);

              // assemble Cuui
              std::vector<int> mylmowner;
              for(size_t index=0; index<la[0].lmowner_.size(); index++) mylmowner.push_back(xfluid_.myrank_);

              Cuui_->FEAssemble(-1, couplingmatrices[1],la[0].lm_,mylmowner, patchlm);


              // assemble rhC_ui_col
              Epetra_SerialDenseVector rhC_ui_eptvec(::View,couplingmatrices[2].A(),patchlm.size());
              LINALG::Assemble(*rhC_ui_col, rhC_ui_eptvec, patchlm, mypatchlmowner);
            }
          }

          // assemble Cuiui
          std::vector<int> mypatchelementslmowner;
          for(size_t index=0; index<patchelementslm.size(); index++) mypatchelementslmowner.push_back(xfluid_.myrank_);

          Cuiui_->FEAssemble(-1,Cuiui, patchelementslm, mypatchelementslmowner, patchelementslm );

        }

        int eid = actele->Id();

        // introduce an vector containing the rows for that values have to be communicated
        // REMARK: when assembling row elements also non-row rows have to be communicated
        std::vector<int> myowner;
        for(size_t index=0; index<la[0].lmowner_.size(); index++) myowner.push_back(xfluid_.myrank_);

        // calls the Assemble function for EpetraFECrs matrices including communication of non-row entries
        sysmat_->FEAssemble(eid, strategy.Elematrix1(), la[0].lm_,myowner,la[0].lm_);

        // REMARK:: call Assemble without lmowner
        // to assemble the residual_col vector on only row elements also column nodes have to be assembled
        // do not exclude non-row nodes (modify the real owner to myowner)
        // after assembly the col vector it has to be exported to the row residual_ vector
        // using the 'Add' flag to get the right value for shared nodes
        LINALG::Assemble(*strategy.Systemvector1(),strategy.Elevector1(),la[0].lm_,myowner);

        set_counter += 1;

      } // end of loop over cellsets // end of assembly for each set of cells
#else
      GEO::CUT::plain_volumecell_set cells;
      std::vector<DRT::UTILS::GaussIntegration> intpoints;
      std::vector<std::vector<double> > refEqns;
      e->VolumeCellGaussPoints( cells, intpoints, refEqns, xfluid_.VolumeCellGaussPointBy_);//modify gauss type

      int count = 0;
      for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        GEO::CUT::VolumeCell * vc = *i;
        std::map<int, std::vector<Epetra_SerialDenseMatrix> > side_coupling;
        if ( vc->Position()==GEO::CUT::Point::outside )
        {
          //const std::vector<int> & nds = vc->NodalDofSet();

          // one set of dofsets
          std::vector<int>  ndstest;
          for (int t=0;t<8; ++t)
            ndstest.push_back(0);

          // actele->LocationVector(discret,nds,la,false);
          actele->LocationVector(discret,ndstest,la,false);

          // get dimension of element matrices and vectors
          // Reshape element matrices and vectors and init to zero
          strategy.ClearElementStorage( la[0].Size(), la[0].Size() );

          {
            TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluidFluid::XFluidFluidState::Evaluate cut domain" );

            // call the element evaluate methods
            int err = impl->Evaluate( ele, discret, la[0].lm_, eleparams, mat,
                                      strategy.Elematrix1(),
                                      strategy.Elematrix2(),
                                      strategy.Elevector1(),
                                      strategy.Elevector2(),
                                      strategy.Elevector3(),
                                      intpoints[count] );
            if (err)
              dserror("Proc %d: Element %d returned err=%d",discret.Comm().MyPID(),actele->Id(),err);
          }

          // do cut interface condition

          std::map<int, std::vector<GEO::CUT::BoundaryCell*> > bcells;
          vc->GetBoundaryCells( bcells );

          if ( bcells.size() > 0 )
          {
            TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluidFluid::XFluidFluidState::Evaluate boundary" );

            std::map<int, std::vector<DRT::UTILS::GaussIntegration> > bintpoints;

            // Attention: switch also the flag in fluid3_impl.cpp
#ifdef BOUNDARYCELL_TRANSFORMATION_OLD
            // original Axel's transformation
            e->BoundaryCellGaussPoints( wizard_->CutWizard().Mesh(), 0, bcells, bintpoints );
#else
            // new Benedikt's transformation
            e->BoundaryCellGaussPointsLin( wizard_->CutWizard().Mesh(), 0, bcells, bintpoints );
#endif

            //std::map<int, std::vector<Epetra_SerialDenseMatrix> > side_coupling;

            std::set<int> begids;
            for (std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
                     bc!=bcells.end(); ++bc )
            {
              int sid = bc->first;
              begids.insert(sid);
            }


            std::vector<int> patchelementslm;
            std::vector<int> patchelementslmowner;
            for ( std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
                  bc!=bcells.end(); ++bc )
            {
              int sid = bc->first;
              DRT::Element * side = cutdiscret.gElement( sid );

              std::vector<int> patchlm;
              std::vector<int> patchlmowner;
              std::vector<int> patchlmstride;
              side->LocationVector(cutdiscret, patchlm, patchlmowner, patchlmstride);

              patchelementslm.reserve( patchelementslm.size() + patchlm.size());
              patchelementslm.insert(patchelementslm.end(), patchlm.begin(), patchlm.end());

              patchelementslmowner.reserve( patchelementslmowner.size() + patchlmowner.size());
              patchelementslmowner.insert( patchelementslmowner.end(), patchlmowner.begin(), patchlmowner.end());

              const size_t ndof_i = patchlm.size();
              const size_t ndof   = la[0].lm_.size();

              std::vector<Epetra_SerialDenseMatrix> & couplingmatrices = side_coupling[sid];
              if ( couplingmatrices.size()!=0 )
                dserror("zero sized vector expected");

              couplingmatrices.resize(3);
              couplingmatrices[0].Reshape(ndof_i,ndof);  //C_uiu
              couplingmatrices[1].Reshape(ndof,ndof_i);  //C_uui
              couplingmatrices[2].Reshape(ndof_i,1);     //rhC_ui
            }

            const size_t nui = patchelementslm.size();
            Epetra_SerialDenseMatrix  Cuiui(nui,nui);

            // all boundary cells that belong to one cut element
            impl->ElementXfemInterfaceMHCS(    ele,
                                              discret,
                                              la[0].lm_,
                                              intpoints[count],
                                              cutdiscret,
                                              bcells,
                                              bintpoints,
                                              side_coupling,
                                              eleparams,
                                              strategy.Elematrix1(),
                                              strategy.Elevector1(),
                                              Cuiui,
                                              cells,
                                              true);

            for ( std::map<int, std::vector<Epetra_SerialDenseMatrix> >::const_iterator sc=side_coupling.begin();
                  sc!=side_coupling.end(); ++sc )
            {
              std::vector<Epetra_SerialDenseMatrix>  couplingmatrices = sc->second;

              int sid = sc->first;

              if ( cutdiscret.HaveGlobalElement(sid) )
              {
                DRT::Element * side = cutdiscret.gElement( sid );
                std::vector<int> patchlm;
                std::vector<int> patchlmowner;
                std::vector<int> patchlmstride;
                side->LocationVector(cutdiscret, patchlm, patchlmowner, patchlmstride);

                // assemble Cuiu
                //create a dummy mypatchlmowner that assembles also non-local rows and communicates the required data
                std::vector<int> mypatchlmowner;
                for(size_t index=0; index<patchlmowner.size(); index++) mypatchlmowner.push_back(xfluid_.myrank_);

                Cuiu_->FEAssemble(-1, couplingmatrices[0],patchlm,mypatchlmowner,la[0].lm_);

                // assemble Cuui
                std::vector<int> mylmowner;
                for(size_t index=0; index<la[0].lmowner_.size(); index++) mylmowner.push_back(xfluid_.myrank_);

                Cuui_->FEAssemble(-1, couplingmatrices[1],la[0].lm_,mylmowner, patchlm);


                // assemble rhC_ui_col
                Epetra_SerialDenseVector rhC_ui_eptvec(::View,couplingmatrices[2].A(),patchlm.size());
                LINALG::Assemble(*rhC_ui_col, rhC_ui_eptvec, patchlm, mypatchlmowner);
              }
            }

            // assemble Cuiui
            std::vector<int> mypatchelementslmowner;
            for(size_t index=0; index<patchelementslm.size(); index++) mypatchelementslmowner.push_back(xfluid_.myrank_);

            Cuiui_->FEAssemble(-1,Cuiui, patchelementslm, mypatchelementslmowner, patchelementslm );
          }

          int eid = actele->Id();

          // introduce an vector containing the rows for that values have to be communicated
          // REMARK: when assembling row elements also non-row rows have to be communicated
          std::vector<int> myowner;
          for(size_t index=0; index<la[0].lmowner_.size(); index++) myowner.push_back(xfluid_.myrank_);

          // calls the Assemble function for EpetraFECrs matrices including communication of non-row entries
          sysmat_->FEAssemble(eid, strategy.Elematrix1(), la[0].lm_,myowner,la[0].lm_);

          // REMARK:: call Assemble without lmowner
          // to assemble the residual_col vector on only row elements also column nodes have to be assembled
          // do not exclude non-row nodes (modify the real owner to myowner)
          // after assembly the col vector it has to be exported to the row residual_ vector
          // using the 'Add' flag to get the right value for shared nodes
          LINALG::Assemble(*strategy.Systemvector1(),strategy.Elevector1(),la[0].lm_,myowner);

        }
        count += 1;
      }
#endif
    }
    else // evaluation of a non-intersected background element
    {
      TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluidFluid::XFluidFluidState::Evaluate normal" );
      // get element location vector, dirichlet flags and ownerships
      actele->LocationVector(*bgdis,la,false);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      strategy.ClearElementStorage( la[0].Size(), la[0].Size() );

      // call the element evaluate method
      int err = impl->Evaluate( ele, *bgdis, la[0].lm_, eleparams, mat,
                                strategy.Elematrix1(),
                                strategy.Elematrix2(),
                                strategy.Elevector1(),
                                strategy.Elevector2(),
                                strategy.Elevector3());

      if (err) dserror("Proc %d: Element %d returned err=%d",bgdis->Comm().MyPID(),actele->Id(),err);

      int eid = actele->Id();

      // introduce an vector containing the rows for that values have to be communicated
      // REMARK: when assembling row elements also non-row rows have to be communicated
      std::vector<int> myowner;
      for(size_t index=0; index<la[0].lmowner_.size(); index++) myowner.push_back(xfluid_.myrank_);

      // calls the Assemble function for EpetraFECrs matrices including communication of non-row entries
      sysmat_->FEAssemble(eid, strategy.Elematrix1(), la[0].lm_,myowner,la[0].lm_);

      // REMARK:: call Assemble without lmowner
      // to assemble the residual_col vector on only row elements also column nodes have to be assembled
      // do not exclude non-row nodes (modify the real owner to myowner)
      // after assembly the col vector it has to be exported to the row residual_ vector
      // using the 'Add' flag to get the right value for shared nodes
      LINALG::Assemble(*strategy.Systemvector1(),strategy.Elevector1(),la[0].lm_,myowner);

    }
  } // end of loop over bgdis

  // call edge stabilization
  if ( xfluid_.edge_based_ or xfluid_.ghost_penalty_ )
  {
    TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 4) EOS" );

    Teuchos::ParameterList faceparams;

    // set additional faceparams according to ghost-penalty terms due to Nitsche's method
    faceparams.set("ghost_penalty_reconstruct", false); // no XFEM timeintegration reconstruction call

    //------------------------------------------------------------
    // loop over row faces

    const Teuchos::RCP<DRT::DiscretizationFaces> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(bgdis, true);

    const int numrowintfaces = xdiscret->NumMyRowFaces();

    // REMARK: in this XFEM framework the whole evaluate routine uses only row internal faces
    // and assembles into EpetraFECrs matrix
    // this is baci-unusual but more efficient in all XFEM applications

    for (int i=0; i<numrowintfaces; ++i)
    {
      DRT::Element* actface = xdiscret->lRowFace(i);

      DRT::ELEMENTS::FluidIntFace * ele = dynamic_cast<DRT::ELEMENTS::FluidIntFace *>( actface );
      if ( ele==NULL ) dserror( "expect FluidIntFace element" );

      edgestab_->EvaluateEdgeStabGhostPenalty(faceparams, bgdis, ele, sysmat_, strategy.Systemvector1(), xfluid_.gmsh_EOS_out_);
    }
  }

  bgdis->ClearState();

  // finalize the complete matrices
  if (xfluid_.coupling_approach_ == CouplingMHCS_XFluid or
      xfluid_.coupling_approach_ == CouplingNitsche_XFluid or
      xfluid_.coupling_approach_ == CouplingMHVS_XFluid)
  {
    // REMARK: for EpetraFECrs matrices Complete() calls the GlobalAssemble() routine to gather entries from all processors
    Cuui_->Complete(*xfluid_.boundarydofrowmap_,*fluiddofrowmap_);
    Cuiu_->Complete(*fluiddofrowmap_,*xfluid_.boundarydofrowmap_);
    Cuiui_->Complete(*xfluid_.boundarydofrowmap_,*xfluid_.boundarydofrowmap_);
  }
  else if (xfluid_.coupling_approach_ == CouplingNitsche_EmbFluid or xfluid_.coupling_approach_ == CouplingNitsche_TwoSided)
  {
    // REMARK: for EpetraFECrs matrices Complete() calls the GlobalAssemble() routine to gather entries from all processors
    Cuui_->Complete(*xfluid_.aledofrowmap_,*fluiddofrowmap_);
    Cuiu_->Complete(*fluiddofrowmap_,*xfluid_.aledofrowmap_);
    Cuiui_->Complete(*xfluid_.aledofrowmap_,*xfluid_.aledofrowmap_);
  }

  //-------------------------------------------------------------------------------
  // export the rhs coupling vector to a row vector
  Epetra_Vector rhC_ui_tmp(rhC_ui_->Map(),true);
  Epetra_Export exporter_rhC_ui_col(rhC_ui_col->Map(),rhC_ui_tmp.Map());
  int err3 = rhC_ui_tmp.Export(*rhC_ui_col,exporter_rhC_ui_col,Add);
  if (err3) dserror("Export using exporter returned err=%d",err3);
  rhC_ui_->Update(1.0,rhC_ui_tmp,0.0);

  //-------------------------------------------------------------------------------
  // need to export residual_col to systemvector1 (residual_)
  Epetra_Vector res_tmp(residual_->Map(),true);
  Epetra_Export exporter(strategy.Systemvector1()->Map(),res_tmp.Map());
  int err2 = res_tmp.Export(*strategy.Systemvector1(),exporter,Add);
  if (err2) dserror("Export using exporter returned err=%d",err2);
  residual_->Update(1.0,res_tmp,1.0);


  //-------------------------------------------------------------------------------
  // finalize the complete matrix
  // REMARK: for EpetraFECrs matrices Complete() calls the GlobalAssemble() routine to gather entries from all processors
  sysmat_->Complete();

  //////////////////////////////////////////////////////////////////////////////////////////
  //
  // loop over column elements of fluid-ale discretization
  //
  ////////////////////////////////////////////////////////////////////////////////////////
  const int numcolaleele = embdis->NumMyColElements();
  for (int i=0; i<numcolaleele; ++i)
  {
    DRT::Element* actaleele = embdis->lColElement(i);
    Teuchos::RCP<MAT::Material> mat = actaleele->Material();

    DRT::ELEMENTS::Fluid * aleele = dynamic_cast<DRT::ELEMENTS::Fluid *>( actaleele );
    if ( aleele==NULL )
    {
      dserror( "expect fluid element" );
    }

    DRT::ELEMENTS::FluidEleInterface * impl = DRT::ELEMENTS::FluidFactory::ProvideImplXFEM(actaleele->Shape(), "xfem");

    GEO::CUT::ElementHandle * e = wizard_->GetElement( actaleele );
    if ( e!=NULL )
    {
      dserror("ALE element geschnitten?!!!!");
    }
    else
    {
      TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluidFluid::XFluidFluidState::Evaluate normal" );

      // get element location vector, dirichlet flags and ownerships
      actaleele->LocationVector(*embdis,alela,false);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      alestrategy.ClearElementStorage( alela[0].Size(), alela[0].Size() );

      // call the element evaluate method
      int err = impl->Evaluate( aleele, *embdis, alela[0].lm_, eleparams, mat,
                                alestrategy.Elematrix1(),
                                alestrategy.Elematrix2(),
                                alestrategy.Elevector1(),
                                alestrategy.Elevector2(),
                                alestrategy.Elevector3() );

      if (err) dserror("Proc %d: Element %d returned err=%d",embdis->Comm().MyPID(),actaleele->Id(),err);

      int eid = actaleele->Id();
      alestrategy.AssembleMatrix1(eid,alela[0].lm_,alela[0].lm_,alela[0].lmowner_,alela[0].stride_);
      alestrategy.AssembleMatrix2(eid,alela[0].lm_,alela[0].lm_,alela[0].lmowner_,alela[0].stride_);
      alestrategy.AssembleVector1(alela[0].lm_,alela[0].lmowner_);

    }
  } // end of loop over embedded discretization

  // call edge stabilization
  if( xfluid_.edge_based_ or xfluid_.ghost_penalty_ )
  {
    TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 4) EOS" );

    Teuchos::ParameterList faceparams;

    // set additional faceparams according to ghost-penalty terms due to Nitsche's method
    faceparams.set("ghost_penalty_reconstruct", false); // no XFEM timeintegration reconstruction call

    //------------------------------------------------------------
    Teuchos::RCP<Epetra_Vector> ale_residual_col = LINALG::CreateVector(*embdis->DofColMap(),true);

    //------------------------------------------------------------
    const Epetra_Map* rmap = NULL;

    Teuchos::RCP<Epetra_FECrsMatrix> sysmat_FE;
    if (xfluid_.alesysmat_ != Teuchos::null)
    {
      rmap = &(xfluid_.alesysmat_->OperatorRangeMap());
      sysmat_FE = Teuchos::rcp(new Epetra_FECrsMatrix(::Copy,*rmap,256,false));
    }
    else dserror("alesysmat is NULL!");

    Teuchos::RCP<LINALG::SparseMatrix> sysmat_linalg = Teuchos::rcp(new LINALG::SparseMatrix(sysmat_FE,true,true,LINALG::SparseMatrix::FE_MATRIX));

    //------------------------------------------------------------
    // loop over row faces

    const Teuchos::RCP<DRT::DiscretizationFaces> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(xfluid_.embdis_, true);

    const int numrowintfaces = xdiscret->NumMyRowFaces();

    // REMARK: in this XFEM framework the whole evaluate routine uses only row internal faces
    // and assembles into EpetraFECrs matrix
    // this is baci-unusual but more efficient in all XFEM applications
    for (int i=0; i<numrowintfaces; ++i)
    {
      DRT::Element* actface = xdiscret->lRowFace(i);
      DRT::ELEMENTS::FluidIntFace * ele = dynamic_cast<DRT::ELEMENTS::FluidIntFace *>( actface );
      if ( ele==NULL ) dserror( "expect FluidIntFace element" );
      edgestab_->EvaluateEdgeStabStd(faceparams, xfluid_.embdis_, ele, sysmat_linalg, ale_residual_col);
    }

    //------------------------------------------------------------
    sysmat_linalg->Complete();

    (xfluid_.alesysmat_)->Add(*sysmat_linalg, false, 1.0, 1.0);

    //------------------------------------------------------------
    // need to export ale_residual_col to systemvector1 (aleresidual_)
    Epetra_Vector res_tmp(xfluid_.aleresidual_->Map(),true);
    Epetra_Export exporter(ale_residual_col->Map(),res_tmp.Map());
    int err2 = res_tmp.Export(*ale_residual_col,exporter,Add);
    if (err2) dserror("Export using exporter returned err=%d",err2);
    xfluid_.aleresidual_->Update(1.0,res_tmp,1.0);

    //------------------------------------------------------------
  }

  cutdiscret->ClearState();
  embdis->ClearState();

  // finalize the complete matrices
  xfluid_.alesysmat_->Complete();

  // adding rhC_ui_ to fluidale residual
  for (int iter=0; iter<rhC_ui_->MyLength();++iter)
  {
    int rhsdgid = rhC_ui_->Map().GID(iter);
    if (rhC_ui_->Map().MyGID(rhsdgid) == false) dserror("rhsd_ should be on all prossesors");
    if (xfluid_.aleresidual_->Map().MyGID(rhsdgid))
      (*xfluid_.aleresidual_)[xfluid_.aleresidual_->Map().LID(rhsdgid)]=(*xfluid_.aleresidual_)[xfluid_.aleresidual_->Map().LID(rhsdgid)] +
                                                                        (*rhC_ui_)[rhC_ui_->Map().LID(rhsdgid)];
    else dserror("cut dof not on ale discret available");
  }


}

// -------------------------------------------------------------------
// -------------------------------------------------------------------
void FLD::XFluidFluid::XFluidFluidState::GmshOutput( DRT::Discretization & discret,
                                                     DRT::Discretization & alefluiddis,
                                                     DRT::Discretization & cutdiscret,
                                                     const std::string & filename_base,
                                                     int countiter,
                                                     int step,
                                                     Teuchos::RCP<Epetra_Vector> vel,
                                                     Teuchos::RCP<Epetra_Vector> alevel,
                                                     Teuchos::RCP<Epetra_Vector> dispntotal)
{
  Teuchos::RCP<const Epetra_Vector> col_vel =
    DRT::UTILS::GetColVersionOfRowVector(xfluid_.bgdis_, vel);

  Teuchos::RCP<const Epetra_Vector> col_alevel =
    DRT::UTILS::GetColVersionOfRowVector(xfluid_.embdis_, alevel);

  Teuchos::RCP<const Epetra_Vector> col_dis;

  if (xfluid_.alefluid_)
    col_dis = DRT::UTILS::GetColVersionOfRowVector(xfluid_.embdis_, dispntotal);

  const int step_diff = 1;
  const bool screen_out = 0;

  // output for Element and Node IDs
  std::ostringstream filename_base_vel;
  if(countiter > -1) filename_base_vel << filename_base << "_" << countiter << "_"<< step << "_vel";
  else           filename_base_vel << filename_base << "_"<< step << "_vel";
  const std::string filename_vel = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_vel.str(), step, step_diff, screen_out, discret.Comm().MyPID());
  IO::cout << IO::endl;
  std::ofstream gmshfilecontent_vel(filename_vel.c_str());
  gmshfilecontent_vel.setf(std::ios::scientific,std::ios::floatfield);
  gmshfilecontent_vel.precision(16);

  std::ostringstream filename_base_press;
  if(countiter > -1) filename_base_press << filename_base << "_" << countiter << "_" << step << "_press";
  else           filename_base_press << filename_base << "_"<< step << "_press";
  const std::string filename_press = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filename_base_press.str(), step, step_diff, screen_out, discret.Comm().MyPID());
  IO::cout << IO::endl;
  std::ofstream gmshfilecontent_press(filename_press.c_str());
  gmshfilecontent_press.setf(std::ios::scientific,std::ios::floatfield);
  gmshfilecontent_press.precision(16);

  if(countiter > -1) // for residual output
  {
    gmshfilecontent_vel   << "View \"" << "SOL " << "vel "   << countiter << "_" << step <<   "\" {\n";
    gmshfilecontent_press << "View \"" << "SOL " << "press " << countiter << "_" << step << "\" {\n";
  }
  else
  {
    gmshfilecontent_vel   << "View \"" << "SOL " << "vel "  << "_" << step << "\" {\n";
    gmshfilecontent_press << "View \"" << "SOL " << "press " << "_" << step << "\" {\n";
  }


  const int numrowele = discret.NumMyRowElements();
  for (int i=0; i<numrowele; ++i)
  {
    DRT::Element* actele = discret.lRowElement(i);

    GEO::CUT::ElementHandle * e = wizard_->GetElement( actele );
    if ( e!=NULL )
    {
#ifdef DOFSETS_NEW

      std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
      std::vector< std::vector<int> > nds_sets;

      e->GetVolumeCellsDofSets( cell_sets, nds_sets );

      int set_counter = 0;

      for( std::vector< GEO::CUT::plain_volumecell_set>::iterator s=cell_sets.begin();
           s!=cell_sets.end();
           s++)
      {
        GEO::CUT::plain_volumecell_set & cells = *s;

        const std::vector<int> & nds = nds_sets[set_counter];

        for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
        {
          GEO::CUT::VolumeCell * vc = *i;
          if ( vc->Position()==GEO::CUT::Point::outside )
          {
            if ( e->IsCut() )
            {
              GmshOutputVolumeCell( discret, gmshfilecontent_vel, gmshfilecontent_press, actele, e, vc, col_vel, nds );
            }
            else
            {
              GmshOutputElement( discret, gmshfilecontent_vel, gmshfilecontent_press, actele, col_vel );
            }
          }
        }
        set_counter += 1;
      }

#else
      GEO::CUT::plain_volumecell_set cells;
      std::vector<DRT::UTILS::GaussIntegration> intpoints;
      e->VolumeCellGaussPoints( cells, intpoints,xfluid_.VolumeCellGaussPointBy_);//modify gauss type
      int count = 0;
      for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        GEO::CUT::VolumeCell * vc = *i;
        if ( vc->Position()==GEO::CUT::Point::outside )
        {
#ifdef DOFSETS_NEW
          const std::vector<int> & nds = vc->NodalDofSet();
          if ( e->IsCut() )
            GmshOutputVolumeCell( discret, gmshfilecontent_vel, gmshfilecontent_press, actele, e, vc, col_vel, nds );
#else
          std::vector<int>  ndstest;
          for (int t=0;t<8; ++t)
            ndstest.push_back(0);
          if ( e->IsCut() )
             GmshOutputVolumeCell( discret, gmshfilecontent_vel, gmshfilecontent_press, actele, e, vc, col_vel, ndstest );
#endif
          else
          {
            GmshOutputElement( discret, gmshfilecontent_vel, gmshfilecontent_press, actele, col_vel );
          }
        }
      }
      count += 1;
#endif
    }
    else
    {
      GmshOutputElement( discret, gmshfilecontent_vel, gmshfilecontent_press, actele, col_vel);
    }
  }


  // comment out to have sepperate views of embedded and background fluid
  gmshfilecontent_vel << "};\n";
  gmshfilecontent_press << "};\n";

  if(countiter > -1) // for residual output
  {
    gmshfilecontent_vel   << "View \"" << "SOL " << "embedded " << countiter << "_" << step << "\" {\n";
    gmshfilecontent_press << "View \"" << "SOL " << "embedded " << countiter << "_" << step << "\" {\n";
  }
  else
  {
    gmshfilecontent_vel   << "View \"" << "SOL " << "embedded " << "_" << step << "\" {\n";
    gmshfilecontent_press << "View \"" << "SOL " << "embedded " << "_" << step << "\" {\n";
  }


  const int numalerowele = alefluiddis.NumMyRowElements();
  for (int i=0; i<numalerowele; ++i)
  {
    DRT::Element* actele = alefluiddis.lRowElement(i);
    GmshOutputElementEmb( alefluiddis, gmshfilecontent_vel, gmshfilecontent_press, actele,col_alevel,col_dis );
  }

  gmshfilecontent_vel << "};\n";
  gmshfilecontent_press << "};\n";

  if(countiter > -1) // for residual output
  {
    gmshfilecontent_vel   << "View \"" << "SOL " << "void " << countiter << "_" << step  << "\" {\n";
    gmshfilecontent_press << "View \"" << "SOL " << "void " << countiter << "_" << step  << "\" {\n";
  }
  else
  {
    gmshfilecontent_vel   << "View \"" << "SOL " << "void " << "_" << step  << "\" {\n";
    gmshfilecontent_press << "View \"" << "SOL " << "void " << "_" << step  << "\" {\n";
  }

  for (int i=0; i<numrowele; ++i)
  {
    DRT::Element* actele = discret.lRowElement(i);

    GEO::CUT::ElementHandle * e = wizard_->GetElement( actele );
    if ( e!=NULL )
    {
      GEO::CUT::plain_volumecell_set cells;
      std::vector<DRT::UTILS::GaussIntegration> intpoints;

      e->GetVolumeCells( cells );
      e->VolumeCellGaussPoints( cells, intpoints,xfluid_.VolumeCellGaussPointBy_);//modify gauss type

      int count = 0;
      for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
      {
        GEO::CUT::VolumeCell * vc = *i;
        if ( vc->Position()==GEO::CUT::Point::outside )
        {
          if ( e->IsCut() )
          {
            GmshOutputElement( discret, gmshfilecontent_vel,  gmshfilecontent_press, actele, col_vel);
          }
        }
      }
      count += 1;
    }
  }
  gmshfilecontent_vel << "};\n";
  gmshfilecontent_press << "};\n";
}
// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::XFluidFluidState::GmshOutputElement( DRT::Discretization & discret,
                                                  std::ofstream & vel_f,
                                                  std::ofstream & press_f,
                                                  DRT::Element * actele,
                                                  Teuchos::RCP<const Epetra_Vector> vel )
{
  DRT::Element::LocationArray la( 1 );

  // get element location vector, dirichlet flags and ownerships
  actele->LocationVector(discret,la,false);

  std::vector<double> m(la[0].lm_.size());
  DRT::UTILS::ExtractMyValues(*vel,m,la[0].lm_);

  switch ( actele->Shape() )
  {
  case DRT::Element::hex8:
  case DRT::Element::hex20:
    vel_f << "VH(";
    press_f << "SH(";
    break;
  default:
    dserror( "unsupported shape" ); break;
  }

//  for ( int i=0; i<actele->NumNode(); ++i )
  for ( int i=0; i<8; ++i )
  {
    if ( i > 0 )
    {
      vel_f << ",";
      press_f << ",";
    }
    const double * x = actele->Nodes()[i]->X();
    vel_f   << x[0] << "," << x[1] << "," << x[2];
    press_f << x[0] << "," << x[1] << "," << x[2];
  }
  vel_f << "){";
  press_f << "){";

//  for ( int i=0; i<actele->NumNode(); ++i )
  for ( int i=0; i<8; ++i )
  {
    if ( i > 0 )
    {
      vel_f << ",";
      press_f << ",";
    }
    int j = 4*i;
    vel_f   << m[j] << "," << m[j+1] << "," << m[j+2];
    press_f << m[j+3];
  }

  vel_f << "};\n";
  press_f << "};\n";
}

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::XFluidFluidState::GmshOutputElementEmb( DRT::Discretization & discret,
                                                               std::ofstream & vel_f,
                                                               std::ofstream & press_f,
                                                               DRT::Element * actele,
                                                               Teuchos::RCP<const Epetra_Vector> vel,
                                                               Teuchos::RCP<const Epetra_Vector> disp)
{
  DRT::Element::LocationArray la( 1 );

  // get element location vector, dirichlet flags and ownerships
  actele->LocationVector(discret,la,false);

  std::vector<double> m(la[0].lm_.size());
  DRT::UTILS::ExtractMyValues(*vel,m,la[0].lm_);

  std::vector<double> dis(la[0].lm_.size());
  if (xfluid_.alefluid_)
    DRT::UTILS::ExtractMyValues(*disp,dis,la[0].lm_);

  switch ( actele->Shape() )
  {
  case DRT::Element::hex8:
  case DRT::Element::hex20:
    vel_f << "VH(";
    press_f << "SH(";
    break;
  default:
    dserror( "unsupported shape" ); break;
  }

  for ( int i=0; i<8; ++i )
  {
    if ( i > 0 )
    {
      vel_f << ",";
      press_f << ",";
    }
    const double * x = actele->Nodes()[i]->X();
    int k = 4*i;
    if (xfluid_.alefluid_)
    {
      vel_f   << x[0]+dis[k] << "," << x[1]+dis[k+1] << "," << x[2]+dis[k+2];
      press_f << x[0]+dis[k] << "," << x[1]+dis[k+1]<< "," << x[2]+dis[k+2];
    }
    else
    {
      vel_f   << x[0] << "," << x[1] << "," << x[2];
      press_f << x[0] << "," << x[1] << "," << x[2];
    }
  }
  vel_f << "){";
  press_f << "){";

  for ( int i=0; i<8; ++i )
  {
    if ( i > 0 )
    {
      vel_f << ",";
      press_f << ",";
    }
    int j = 4*i;
    vel_f   << m[j] << "," << m[j+1] << "," << m[j+2];
    press_f << m[j+3];
  }

  vel_f << "};\n";
  press_f << "};\n";
}
// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::XFluidFluidState::GmshOutputVolumeCell( DRT::Discretization & discret,
                                                     std::ofstream & vel_f,
                                                     std::ofstream & press_f,
                                                     DRT::Element * actele,
                                                     GEO::CUT::ElementHandle * e,
                                                     GEO::CUT::VolumeCell * vc,
                                                     Teuchos::RCP<const Epetra_Vector> velvec,
                                                     const std::vector<int> & nds)
{

  DRT::Element::LocationArray la( 1 );

  // get element location vector, dirichlet flags and ownerships
  actele->LocationVector(discret,nds,la,false);
  //actele->LocationVector(discret,nds,la,false);

   std::vector<double> m(la[0].lm_.size());

  DRT::UTILS::ExtractMyValues(*velvec,m,la[0].lm_);

  Epetra_SerialDenseMatrix vel( 3, actele->NumNode() );
  Epetra_SerialDenseMatrix press( 1, actele->NumNode() );

  for ( int i=0; i<actele->NumNode(); ++i )
  {
    vel( 0, i ) = m[4*i+0];
    vel( 1, i ) = m[4*i+1];
    vel( 2, i ) = m[4*i+2];
    press( 0, i ) = m[4*i+3];
  }

  // facet based output for cut volumes
  // integrationcells are not available because tessellation is not used
  if( xfluid_.VolumeCellGaussPointBy_!=INPAR::CUT::VCellGaussPts_Tessellation )
  {
    const GEO::CUT::plain_facet_set & facete = vc->Facets();
    for(GEO::CUT::plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
    {
      // split facet into tri and quad cell
      GEO::CUT::Facet *fe = *i;
      std::vector<std::vector<GEO::CUT::Point*> > split;
      std::vector<GEO::CUT::Point*> corners = fe->CornerPoints();

      if( corners.size()==3 ) // only Tri can be used directly. Quad may be concave
        split.push_back( corners );
      else
      {
        if( !fe->IsFacetSplit() )
          fe->SplitFacet( fe->CornerPoints() );
         split = fe->GetSplitCells();
      }

      for( std::vector<std::vector<GEO::CUT::Point*> >::const_iterator j=split.begin();
                                                                       j!=split.end();j++ )
      {
        std::vector<GEO::CUT::Point*> cell = *j;

        switch ( cell.size() )
        {
        case 3:
          vel_f << "VT(";
          press_f << "ST(";
          break;
        case 4:
          vel_f << "VQ(";
          press_f << "SQ(";
          break;
        default:
          dserror( "splitting facets failed" ); break;
        }

        for ( unsigned k=0; k<cell.size(); ++k )
        {
          if ( k > 0 )
          {
            vel_f << ",";
            press_f << ",";
          }
          const double * x = cell[k]->X();
          vel_f   << x[0] << "," << x[1] << "," << x[2];
          press_f << x[0] << "," << x[1] << "," << x[2];
        }
        vel_f << "){";
        press_f << "){";

        for ( unsigned k=0; k<cell.size(); ++k )
        {
          LINALG::Matrix<3,1> v( true );
          LINALG::Matrix<1,1> p( true );

          GEO::CUT::Point * point = cell[k];
          const LINALG::Matrix<3,1> & rst = e->LocalCoordinates( point );

          switch ( actele->Shape() )
          {
          case DRT::Element::hex8:
          {
            const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
            LINALG::Matrix<numnodes,1> funct;
            DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::hex8 );
            LINALG::Matrix<3,numnodes> velocity( vel, true );
            LINALG::Matrix<1,numnodes> pressure( press, true );

            v.Multiply( 1, velocity, funct, 1 );
            p.Multiply( 1, pressure, funct, 1 );
            break;
          }
          case DRT::Element::hex20:
          {
            // TODO: check the output for hex20
            const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
            LINALG::Matrix<numnodes,1> funct;
            DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::hex20 );
            LINALG::Matrix<3,numnodes> velocity( vel, true );
            LINALG::Matrix<1,numnodes> pressure( press, true );

            v.Multiply( 1, velocity, funct, 1 );
            p.Multiply( 1, pressure, funct, 1 );
            break;
          }
          default:
            dserror( "unsupported shape" ); break;
          }

          if ( k > 0 )
          {
            vel_f << ",";
            press_f << ",";
          }
          vel_f   << v( 0 ) << "," << v( 1 ) << "," << v( 2 );
          press_f << p( 0 );
        }

        vel_f << "};\n";
        press_f << "};\n";
      }
    }
  }
  // integrationcells based output for tessellation
  else
  {
    const GEO::CUT::plain_integrationcell_set & intcells = vc->IntegrationCells();
    for ( GEO::CUT::plain_integrationcell_set::const_iterator i=intcells.begin();
          i!=intcells.end();
          ++i )
    {
      GEO::CUT::IntegrationCell * ic = *i;

      const std::vector<GEO::CUT::Point*> & points = ic->Points();
      Epetra_SerialDenseMatrix values( 4, points.size() );

      switch ( ic->Shape() )
      {
      case DRT::Element::hex8:
        vel_f << "VH(";
        press_f << "SH(";
        break;
      case DRT::Element::tet4:
        vel_f << "VS(";
        press_f << "SS(";
        break;
      default:
        dserror( "unsupported shape" ); break;
      }

      for ( unsigned i=0; i<points.size(); ++i )
      {
        if ( i > 0 )
        {
          vel_f << ",";
          press_f << ",";
        }
        const double * x = points[i]->X();
        vel_f   << x[0] << "," << x[1] << "," << x[2];
        press_f << x[0] << "," << x[1] << "," << x[2];
      }
      vel_f << "){";
      press_f << "){";

      for ( unsigned i=0; i<points.size(); ++i )
      {
        LINALG::Matrix<3,1> v( true );
        LINALG::Matrix<1,1> p( true );

         GEO::CUT::Point * point = points[i];
        const LINALG::Matrix<3,1> & rst = e->LocalCoordinates( point );

        switch ( actele->Shape() )
        {
        case DRT::Element::hex8:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
          LINALG::Matrix<numnodes,1> funct;
          DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::hex8 );
          LINALG::Matrix<3,numnodes> velocity( vel, true );
          LINALG::Matrix<1,numnodes> pressure( press, true );

          v.Multiply( 1, velocity, funct, 1 );
          p.Multiply( 1, pressure, funct, 1 );
          break;
        }
        case DRT::Element::hex20:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
          LINALG::Matrix<numnodes,1> funct;
          DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::hex20 );
          LINALG::Matrix<3,numnodes> velocity( vel, true );
          LINALG::Matrix<1,numnodes> pressure( press, true );

          v.Multiply( 1, velocity, funct, 1 );
          p.Multiply( 1, pressure, funct, 1 );
          break;
        }
        default:
          dserror( "unsupported shape" ); break;
        }

        if ( i > 0 )
        {
          vel_f << ",";
          press_f << ",";
        }
        vel_f   << v( 0 ) << "," << v( 1 ) << "," << v( 2 );
        press_f << p( 0 );
      }

      vel_f << "};\n";
      press_f << "};\n";
    }
  }
}

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::XFluidFluidState::GmshOutputBoundaryCell( DRT::Discretization & discret,
                                                                 DRT::Discretization & cutdiscret,
                                                                 std::ofstream & bound_f,
                                                                 DRT::Element * actele,
                                                                 GEO::CUT::ElementHandle * e,
                                                                 GEO::CUT::VolumeCell * vc )
{
  LINALG::Matrix<3,1> normal;
  LINALG::Matrix<2,2> metrictensor;
  double drs;

  GEO::CUT::MeshIntersection & mesh = wizard_->CutWizard().Mesh();

  std::map<int, std::vector<GEO::CUT::BoundaryCell*> > bcells;
  vc->GetBoundaryCells( bcells );
  for ( std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::iterator i=bcells.begin();
        i!=bcells.end();
        ++i )
  {
    int sid = i->first;
    std::vector<GEO::CUT::BoundaryCell*> & bcs = i->second;

    DRT::Element * side = cutdiscret.gElement( sid );
    GEO::CUT::SideHandle * s = mesh.GetCutSide( sid, 0 );

    const int numnodes = side->NumNode();
    DRT::Node ** nodes = side->Nodes();
    Epetra_SerialDenseMatrix side_xyze( 3, numnodes );
    for ( int i=0; i<numnodes; ++i )
    {
      const double * x = nodes[i]->X();
      std::copy( x, x+3, &side_xyze( 0, i ) );
    }

    for ( std::vector<GEO::CUT::BoundaryCell*>::iterator i=bcs.begin();
          i!=bcs.end();
          ++i )
    {
      GEO::CUT::BoundaryCell * bc = *i;

      switch ( bc->Shape() )
      {
      case DRT::Element::quad4:
        bound_f << "VQ(";
        break;
      case DRT::Element::tri3:
        bound_f << "VT(";
        break;
      default:
        dserror( "unsupported shape" ); break;
      }

      const std::vector<GEO::CUT::Point*> & points = bc->Points();
      for ( std::vector<GEO::CUT::Point*>::const_iterator i=points.begin();
            i!=points.end();
            ++i )
      {
        GEO::CUT::Point * p = *i;

        if ( i!=points.begin() )
          bound_f << ",";

        const double * x = p->X();
        bound_f << x[0] << "," << x[1] << "," << x[2];
      }

      bound_f << "){";

      for ( std::vector<GEO::CUT::Point*>::const_iterator i=points.begin();
            i!=points.end();
            ++i )
      {
        GEO::CUT::Point * p = *i;

        const LINALG::Matrix<2,1> & eta = s->LocalCoordinates( p );

        switch ( side->Shape() )
        {
        case DRT::Element::tri3:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement;
          LINALG::Matrix<3,numnodes> xyze( side_xyze, true );
          //LINALG::Matrix<numnodes,1> funct;
          LINALG::Matrix<2,numnodes> deriv;
          //DRT::UTILS::shape_function_2D( funct, eta( 0 ), eta( 1 ), DRT::Element::quad4 );
          DRT::UTILS::shape_function_2D_deriv1( deriv, eta( 0 ), eta( 1 ), DRT::Element::tri3 );
          DRT::UTILS::ComputeMetricTensorForBoundaryEle<DRT::Element::tri3>( xyze, deriv, metrictensor, drs, &normal );
          //x.Multiply( xyze, funct );
          break;
        }
        case DRT::Element::quad4:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement;
          LINALG::Matrix<3,numnodes> xyze( side_xyze, true );
          //LINALG::Matrix<numnodes,1> funct;
          LINALG::Matrix<2,numnodes> deriv;
          //DRT::UTILS::shape_function_2D( funct, eta( 0 ), eta( 1 ), DRT::Element::quad4 );
          DRT::UTILS::shape_function_2D_deriv1( deriv, eta( 0 ), eta( 1 ), DRT::Element::quad4 );
          DRT::UTILS::ComputeMetricTensorForBoundaryEle<DRT::Element::quad4>( xyze, deriv, metrictensor, drs, &normal );
          //x.Multiply( xyze, funct );
          break;
        }
        case DRT::Element::quad8:
        {
          const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad8>::numNodePerElement;
          LINALG::Matrix<3,numnodes> xyze( side_xyze, true );
          //LINALG::Matrix<numnodes,1> funct;
          LINALG::Matrix<2,numnodes> deriv;
          //DRT::UTILS::shape_function_2D( funct, eta( 0 ), eta( 1 ), DRT::Element::quad4 );
          DRT::UTILS::shape_function_2D_deriv1( deriv, eta( 0 ), eta( 1 ), DRT::Element::quad8 );
          DRT::UTILS::ComputeMetricTensorForBoundaryEle<DRT::Element::quad8>( xyze, deriv, metrictensor, drs, &normal );
          //x.Multiply( xyze, funct );
          break;
        }
        default:
          dserror( "unsupported side shape %d", side->Shape() ); break;
        }

        if ( i!=points.begin() )
          bound_f << ",";

        bound_f << normal( 0 ) << "," << normal( 1 ) << "," << normal( 2 );
      }
      bound_f << "};\n";
    }
  }
}

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
FLD::XFluidFluid::XFluidFluid(
    const Teuchos::RCP<DRT::Discretization>&      actdis,
    const Teuchos::RCP<DRT::Discretization>&      embdis,
    const Teuchos::RCP<LINALG::Solver>&           solver,
    const Teuchos::RCP<Teuchos::ParameterList>&   params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output,
    bool                                          alefluid ,
    bool                                          monolithicfluidfluidfsi
):TimInt(actdis, solver, params, output),
  bgdis_(discret_),
  embdis_(embdis),
  alefluid_(alefluid),
  monolithicfluidfluidfsi_(monolithicfluidfluidfsi)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluid::Init()
{
  numdim_            = DRT::Problem::Instance()->NDim(); //params_->get<int>("DIM");
  dtp_               = params_->get<double>("time step size");
  theta_             = params_->get<double>("theta");
  newton_            = DRT::INPUT::get<INPAR::FLUID::LinearisationAction>(*params_, "Linearisation");
  convform_          = params_->get<string>("form of convective term","convective");

  Teuchos::ParameterList&   params_xfem    = params_->sublist("XFEM");
  Teuchos::ParameterList&   params_xf_gen  = params_->sublist("XFLUID DYNAMIC/GENERAL");
  Teuchos::ParameterList&   params_xf_stab = params_->sublist("XFLUID DYNAMIC/STABILIZATION");

  // get general XFEM specific parameters
  maxnumdofsets_           = params_->sublist("XFEM").get<int>("MAX_NUM_DOFSETS");
  VolumeCellGaussPointBy_  = DRT::INPUT::IntegralValue<INPAR::CUT::VCellGaussPts>(params_xfem, "VOLUME_GAUSS_POINTS_BY");
  BoundCellGaussPointBy_   = DRT::INPUT::IntegralValue<INPAR::CUT::BCellGaussPts>(params_xfem, "BOUNDARY_GAUSS_POINTS_BY");

  // get interface stabilization specific parameters
  coupling_method_       = DRT::INPUT::IntegralValue<INPAR::XFEM::CouplingMethod>(params_xf_stab,"COUPLING_METHOD");
  coupling_strategy_  = DRT::INPUT::IntegralValue<INPAR::XFEM::CouplingStrategy>(params_xf_stab,"COUPLING_STRATEGY");

  hybrid_lm_l2_proj_ = DRT::INPUT::IntegralValue<INPAR::XFEM::Hybrid_LM_L2_Proj>(params_xf_stab, "HYBRID_LM_L2_PROJ");

  xff_conv_stab_scaling_ = DRT::INPUT::IntegralValue<INPAR::XFEM::XFF_ConvStabScaling>(params_xf_stab,"XFF_CONV_STAB_SCALING");
  conv_stab_scaling_     = DRT::INPUT::IntegralValue<INPAR::XFEM::ConvStabScaling>(params_xf_stab,"CONV_STAB_SCALING");

  // set flag if any edge-based fluid stabilization has to integrated as std or gp stabilization
  edge_based_        = (params_->sublist("RESIDUAL-BASED STABILIZATION").get<string>("STABTYPE")=="edge_based"
                        or params_->sublist("EDGE-BASED STABILIZATION").get<string>("EOS_PRES")        != "none"
                        or params_->sublist("EDGE-BASED STABILIZATION").get<string>("EOS_CONV_STREAM") != "none"
                        or params_->sublist("EDGE-BASED STABILIZATION").get<string>("EOS_CONV_CROSS")  != "none"
                        or params_->sublist("EDGE-BASED STABILIZATION").get<string>("EOS_DIV")         != "none");

  // set flag if a viscous or transient (1st or 2nd order) ghost-penalty stabiliation due to Nitsche's method has to be integrated
  ghost_penalty_                    = (    (bool)DRT::INPUT::IntegralValue<int>(params_xf_stab,"GHOST_PENALTY_STAB")
                                        or (bool)DRT::INPUT::IntegralValue<int>(params_xf_stab,"GHOST_PENALTY_TRANSIENT_STAB")
                                        or (bool)DRT::INPUT::IntegralValue<int>(params_xf_stab,"GHOST_PENALTY_2nd_STAB") );

  ghost_penalty_fac_ = params_xf_stab.get<double>("GHOST_PENALTY_FAC", 0.0);

  velgrad_interface_stab_ =  (bool)DRT::INPUT::IntegralValue<int>(params_xf_stab,"VELGRAD_INTERFACE_STAB");

  presscoupling_interface_stab_ = (bool)DRT::INPUT::IntegralValue<int>(params_xf_stab,"PRESSCOUPLING_INTERFACE_STAB");
  presscoupling_interface_fac_  = params_xf_stab.get<double>("PRESSCOUPLING_INTERFACE_FAC", 0.0);

  nitsche_evp_ = (DRT::INPUT::IntegralValue<INPAR::XFEM::ViscStab_TraceEstimate>(params_xf_stab,"VISC_STAB_TRACE_ESTIMATE")
       == INPAR::XFEM::ViscStab_TraceEstimate_eigenvalue);

  // get general XFEM specific parameters

  monolithic_approach_= DRT::INPUT::IntegralValue<INPAR::XFEM::Monolithic_xffsi_Approach>(params_->sublist("XFLUID DYNAMIC/GENERAL"),"MONOLITHIC_XFFSI_APPROACH");
  xfem_timeintapproach_= DRT::INPUT::IntegralValue<INPAR::XFEM::XFluidFluidTimeInt>(params_->sublist("XFLUID DYNAMIC/GENERAL"),"XFLUIDFLUID_TIMEINT");
  relaxing_ale_ = (bool)DRT::INPUT::IntegralValue<int>(params_xf_gen,"RELAXING_ALE");
  relaxing_ale_every_ = params_xf_gen.get<int>("RELAXING_ALE_EVERY", 0);

  gmsh_count_ = 0;

  // load GMSH output flags
  gmsh_sol_out_          = (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_SOL_OUT");
  gmsh_debug_out_        = (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_DEBUG_OUT");
  gmsh_debug_out_screen_ = (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_DEBUG_OUT_SCREEN");
  gmsh_EOS_out_          = ((bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_EOS_OUT") && (edge_based_ or ghost_penalty_));
  gmsh_discret_out_      = (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_DISCRET_OUT");
  gmsh_step_diff_        = 500;
  gmsh_cut_out_          = (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_CUT_OUT");

  // set the element name for boundary elements BELE3_4
  std::string element_name;

  // check input parameters set for boundary/interface coupling and
  // set the action

  if(coupling_strategy_ == INPAR::XFEM::Xfluid_Sided_weak_DBC)
    dserror("do not choose Xfluid_Sided_weak_DBC for fluid-fluid-couplings!");

  switch (coupling_method_)
  {
  case INPAR::XFEM::Hybrid_LM_Cauchy_stress:
    // method needs just 3 dofs, but we use 4 to be consistent with Nitsche & MHVS
    element_name = "BELE3_4";
    IO::cout << "Coupling Mixed/Hybrid Cauchy Stress-Based LM Xfluid-Sided" << IO::endl;
    if (coupling_strategy_ != INPAR::XFEM::Xfluid_Sided_Coupling)
      dserror("Choose Xfluid_Sided_Coupling for MHCS");
    coupling_approach_ = CouplingMHCS_XFluid;
    break;
  case INPAR::XFEM::Hybrid_LM_viscous_stress:
    element_name = "BELE3_4"; // use 4 dofs
    IO::cout << "XFEM interface coupling method: TypeMHVS" << IO::endl;
    if (coupling_strategy_ != INPAR::XFEM::Xfluid_Sided_Coupling)
      dserror("Embedded-sided or two-sided MHVS coupling not supported yet.");
    coupling_approach_ = CouplingMHVS_XFluid;
    IO::cout << "Coupling Mixed/Hybrid Viscous Stress-Based LM Xfluid-Sided" << IO::endl;
    break;
  case INPAR::XFEM::Nitsche:
    element_name = "BELE3_4"; // use 4 dofs
    IO::cout << "XFEM interface coupling method: ";
    if (coupling_strategy_ == INPAR::XFEM::Two_Sided_Coupling)
    {
      coupling_approach_ = CouplingNitsche_TwoSided;
      IO::cout << "Coupling Nitsche Two-Sided" << IO::endl;
      IO::cout << "ATTENTION: choose reasonable weights (k1,k2) for mortaring" << IO::endl;
    }
    if (coupling_strategy_ == INPAR::XFEM::Xfluid_Sided_Coupling)
    {
      coupling_approach_ = CouplingNitsche_XFluid;
      IO::cout << "Coupling Nitsche Xfluid-Sided" << IO::endl;
    }
    if (coupling_strategy_ == INPAR::XFEM::Embedded_Sided_Coupling)
    {
      coupling_approach_ = CouplingNitsche_EmbFluid;
      IO::cout << "Coupling Nitsche Embedded-Sided" << IO::endl;
    }
    break;
  default:
    dserror("Unknown fluid-fluid coupling type."); break;
  }

  if (monolithicfluidfluidfsi_)
  {
    switch (monolithic_approach_)
    {
    case INPAR::XFEM::XFFSI_FixedALE_Partitioned:
      monotype_ = FixedALEPartitioned;
      break;
    case INPAR::XFEM::XFFSI_Full_Newton:
      monotype_ = FullyNewton;
      break;
    case INPAR::XFEM::XFFSI_FixedALE_Interpolation:
      monotype_ = FixedALEInterpolation;
      break;
    default:
      dserror("Unknown monolithic XFFSI approach."); break;
    }
  }
  else
    monotype_ = NoMonolithicXFFSI;

  // check if fluidfluidcoupling is set for the right type. fluidfluidcoupling
  // is just related to fixedale monolithic approach
  std::vector<int> condnodes;
  std::vector<int> condnodesglobal;
  DRT::UTILS::FindConditionedNodes(*embdis_,"FluidFluidCoupling",condnodes);

  //information how many processors work at all
  std::vector<int> allproc(embdis_->Comm().NumProc());

  //in case of n processors allproc becomes a vector with entries (0,1,...,n-1)
  for (int i=0; i<embdis_->Comm().NumProc(); ++i) allproc[i] = i;

  //gathers information of condnodes of all processors
  LINALG::Gather<int>(condnodes,condnodesglobal,(int)embdis_->Comm().NumProc(),&allproc[0],embdis_->Comm());

  if (((monotype_ == FixedALEPartitioned) or (monotype_ == FixedALEInterpolation)) and (condnodesglobal.size()==0))
    dserror("FluidFluidCoupling condition is missing!");
  else if (((monotype_ == FullyNewton) or (monotype_ == NoMonolithicXFFSI)) and (condnodesglobal.size()!=0))
    dserror("FluidFluidCoupling condition is not related here!");

  // check xfluid input params
  CheckXFluidFluidParams(params_xfem,params_xf_gen,params_xf_stab);

  PrintStabilizationParams();

  // compute or set 1.0 - theta for time-integration schemes
  if (timealgo_ == INPAR::FLUID::timeint_one_step_theta)  omtheta_ = 1.0 - theta_;
  else                                      omtheta_ = 0.0;

  // parameter for linearization scheme (fixed-point-like or Newton)
  newton_ = DRT::INPUT::get<INPAR::FLUID::LinearisationAction>(*params_, "Linearisation");

  predictor_ = params_->get<std::string>("predictor","steady_state_predictor");

  // form of convective term
  convform_ = params_->get<std::string>("form of convective term","convective");

  emboutput_ = embdis_->Writer();
  emboutput_->WriteMesh(0,0.0);

  // ensure that degrees of freedom in the discretization have been set
  if ( not bgdis_->Filled() or not discret_->HaveDofs() )
    bgdis_->FillComplete();

  std::vector<std::string> conditions_to_copy;
  conditions_to_copy.push_back("XFEMCoupling");
  boundarydis_ = DRT::UTILS::CreateDiscretizationFromCondition(embdis_, "XFEMCoupling", "boundary", element_name, conditions_to_copy);

  if (boundarydis_->NumGlobalNodes() == 0)
  {
    dserror("Empty XFEM-boundary discretization detected!");
  }

  // create node and element distribution with boundarydis elements and nodes
  // ghosted on all processors
  const Epetra_Map noderowmap = *boundarydis_->NodeRowMap();
  const Epetra_Map elemrowmap = *boundarydis_->ElementRowMap();

  // put all boundary nodes and elements onto all processors
  const Epetra_Map nodecolmap = *LINALG::AllreduceEMap(noderowmap);
  const Epetra_Map elemcolmap = *LINALG::AllreduceEMap(elemrowmap);

  // redistribute nodes and elements to column (ghost) map
  boundarydis_->ExportColumnNodes(nodecolmap);
  boundarydis_->ExportColumnElements(elemcolmap);

  boundarydis_->FillComplete();

  // make the dofset of boundarydis be a subset of the embedded dis
  Teuchos::RCP<Epetra_Map> newcolnodemap = DRT::UTILS::ComputeNodeColMap(embdis_,boundarydis_);
  embdis_->Redistribute(*(embdis_->NodeRowMap()), *newcolnodemap);
  Teuchos::RCP<DRT::DofSet> newdofset=Teuchos::rcp(new DRT::TransparentIndependentDofSet(embdis_,false,Teuchos::null));
  boundarydis_->ReplaceDofSet(newdofset); // do not call this with true!!
  boundarydis_->FillComplete();

  //DRT::UTILS::PrintParallelDistribution(*boundarydis_);

  // prepare embedded dis for Nitsche-Coupling-Type Ale-Sided
  // this also needed for error calculation method
  INPAR::FLUID::CalcError calcerr = DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_,"calculate error");
  if ((coupling_strategy_ != INPAR::XFEM::Xfluid_Sided_Coupling) or (calcerr != INPAR::FLUID::no_error_calculation))
  {
    PrepareEmbeddedDistribution();
    CreateBoundaryEmbeddedMap();
  }

  // store a dofset with the complete fluid unknowns
  dofset_out_ = Teuchos::rcp(new DRT::IndependentDofSet());
  dofset_out_->Reset();
  dofset_out_->AssignDegreesOfFreedom(*bgdis_,0,0);
  // split based on complete fluid field (standard splitter that handles one dofset)
  FLD::UTILS::SetupFluidSplit(*bgdis_,*dofset_out_,numdim_,velpressplitterForOutput_);

  // create vector according to the dofset_out row map holding all standard fluid unknowns
  outvec_fluid_ = LINALG::CreateVector(*dofset_out_->DofRowMap(),true);

  // create fluid output object
  output_ = bgdis_->Writer();
  output_->WriteMesh(0,0.0);


  // used to write out owner of elements just once
  firstoutputofrun_ = true;

  // counter for number of written restarts, used to decide when we have to clear the MapStack (explanation see Output() )
  restart_count_ = 0;


  //-------------------------------------------------------------------
  // create internal faces extension for edge based stabilization
  if(edge_based_ or ghost_penalty_)
  {
    Teuchos::RCP<DRT::DiscretizationFaces> actembdis = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(embdis_, true);
    actembdis->CreateInternalFacesExtension();

    Teuchos::RCP<DRT::DiscretizationFaces> actbgdis = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(bgdis_, true);
    actbgdis->CreateInternalFacesExtension();
  }
  //-------------------------------------------------------------------


//   output_ = bgdis_->Writer();
//   output_->WriteMesh(0,0.0);

  // embedded fluid state vectors
  FLD::UTILS::SetupFluidSplit(*embdis_,numdim_,alevelpressplitter_);

  aledofrowmap_ = embdis_->DofRowMap();

  // explicitdirichlet=false, savegraph=true
  alesysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*aledofrowmap_,108,false,true));

  // Vectors passed to the element
  // -----------------------------
  // velocity/pressure at time n+1, n and n-1
  alevelnp_ = LINALG::CreateVector(*aledofrowmap_,true);
  aleveln_  = LINALG::CreateVector(*aledofrowmap_,true);
  alevelnm_ = LINALG::CreateVector(*aledofrowmap_,true);

  // we need the displacement vector of an ALE element if alefluid_ or when we do not use Xfluid-sided-mortaring
  //if(alefluid_ || coupling_strategy_ != INPAR::XFEM::Xfluid_Sided_Mortaring)
  aledispnp_ = LINALG::CreateVector(*aledofrowmap_,true);

  if (alefluid_)
  {
    aledispn_  = LINALG::CreateVector(*aledofrowmap_,true);
    aledispnm_ = LINALG::CreateVector(*aledofrowmap_,true);
    gridv_  = LINALG::CreateVector(*aledofrowmap_,true);
  }

  if (monolithicfluidfluidfsi_)
    aledispnpoldstate_ = LINALG::CreateVector(*aledofrowmap_,true);

  aletotaldispnp_ = LINALG::CreateVector(*aledofrowmap_,true);
  aletotaldispn_ = LINALG::CreateVector(*aledofrowmap_,true);

  // acceleration/(scalar time derivative) at time n+1 and n
  aleaccnp_ = LINALG::CreateVector(*aledofrowmap_,true);
  aleaccn_  = LINALG::CreateVector(*aledofrowmap_,true);

  // velocity/pressure at time n+alpha_F
  alevelaf_ = LINALG::CreateVector(*aledofrowmap_,true);

  // rhs: standard (stabilized) residual vector (rhs for the incremental form)
  aleresidual_     = LINALG::CreateVector(*aledofrowmap_,true);
  aletrueresidual_ = LINALG::CreateVector(*aledofrowmap_,true);

  aleneumann_loads_= LINALG::CreateVector(*aledofrowmap_,true);

  // acceleration/(scalar time derivative) at time n+alpha_M/(n+alpha_M/n)
  aleaccam_ = LINALG::CreateVector(*aledofrowmap_,true);

  // scalar at time n+alpha_F/n+1 and n+alpha_M/n
  // (only required for low-Mach-number case)
  alescaaf_ = LINALG::CreateVector(*aledofrowmap_,true);
  alescaam_ = LINALG::CreateVector(*aledofrowmap_,true);

  // history vector
  alehist_ = LINALG::CreateVector(*aledofrowmap_,true);

  // right hand side vector for linearised solution;
  alerhs_ = LINALG::CreateVector(*aledofrowmap_,true);

  // Nonlinear iteration increment vector
  aleincvel_ = LINALG::CreateVector(*aledofrowmap_,true);

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  alezeros_   = LINALG::CreateVector(*aledofrowmap_,true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  aledbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());

  {
    Teuchos::ParameterList eleparams;
    // other parameters needed by the elements
    embdis_->EvaluateDirichlet(eleparams, alezeros_, Teuchos::null, Teuchos::null,
                               Teuchos::null, aledbcmaps_);

    alezeros_->PutScalar(0.0); // just in case of change
  }

  //--------------------------------------------------------
  // FluidFluid-Boundary Vectros passes to element
  // -------------------------------------------------------
  boundarydofrowmap_ = boundarydis_->DofRowMap();
  ivelnp_ = LINALG::CreateVector(*boundarydofrowmap_,true);
  idispnp_ = LINALG::CreateVector(*boundarydofrowmap_,true);


  //----------------------------------------------------------------------
  turbmodel_ = INPAR::FLUID::dynamic_smagorinsky;
  //----------------------------------------------------------------------

  // -----------------------------------------------------------------
  // set general fluid parameter
  // -----------------------------------------------------------------
  SetElementGeneralFluidXFEMParameter();
  SetFaceGeneralFluidXFEMParameter();

  //--------------------------------------------------
  // XFluidFluid State
  //-----------------------------------------------
  const int restart = DRT::Problem::Instance()->Restart();
  Epetra_Vector idispcol( *boundarydis_->DofColMap() );
  idispcol.PutScalar( 0.0 );
  if (restart)
  {
    ReadRestartEmb(restart);
    LINALG::Export(*aledispn_,idispcol);
  }
  state_ = Teuchos::rcp( new XFluidFluidState( *this, idispcol ) );


  if ( not bgdis_->Filled() or not bgdis_->HaveDofs() )
    bgdis_->FillComplete();

  //----------------------------------------------------------------
  // create embedded-boundary discretization needed for solving
  // Nistche-Eigenvalue-Problem
  //----------------------------------------------------------------
  if (nitsche_evp_)
    CreateEmbeddedBoundarydis();

  //gmsh discretization output
  if(gmsh_discret_out_)
  {
    OutputDiscret();
  }

  if (alefluid_)
    xfluidfluid_timeint_ =  Teuchos::rcp(new XFEM::XFluidFluidTimeIntegration(bgdis_, embdis_, state_->wizard_, step_,
                                                                              xfem_timeintapproach_,*params_));

  return;
}

// ------------------------------------------------------------------------
// prepare embedded discretization for Ale-sided-coupling
// Ghost the outer elements of the embedded discretization on all processors
// -----------------------------------------------------------------------
void FLD::XFluidFluid::PrepareEmbeddedDistribution()
{
#ifdef PARALLEL

  DRT::UTILS::ConditionSelector conds(*embdis_, "XFEMCoupling");
  std::vector<int> embnode_outer;
  std::vector<int> embele_outer;

  // select all outer embedded nodes and elements
  for (int  inode=0; inode<embdis_->NumMyRowNodes(); inode++)
  {
    DRT::Node* embnode = embdis_->lRowNode(inode);
    if (conds.ContainsNode(embnode->Id()))
    {
      embnode_outer.push_back(embnode->Id());
      const size_t numele = embnode->NumElement();
      // get list of adjacent elements of this node
      DRT::Element** adjeles = embnode->Elements();
      int mypid =embdis_->Comm().MyPID();

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
  std::vector<int> allproc(embdis_->Comm().NumProc());

  // in case of n processors allproc becomes a vector with entries (0,1,...,n-1)
  for (int i=0; i<embdis_->Comm().NumProc(); ++i) allproc[i] = i;

  // gathers information of all processors
  LINALG::Gather<int>(embnode_outer,embnode_outer_all,(int)embdis_->Comm().NumProc(),&allproc[0],embdis_->Comm());
  LINALG::Gather<int>(embele_outer,embele_outer_all,(int)embdis_->Comm().NumProc(),&allproc[0],embdis_->Comm());

  // combine the colmap of every processor with the embedded outer map
  for (int  cnode=0; cnode<embdis_->NumMyColNodes(); cnode++)
  {
    DRT::Node* embnode = embdis_->lColNode(cnode);
    embnode_outer_all.push_back(embnode->Id());
  }
  for (int  cele=0; cele<embdis_->NumMyColElements(); cele++)
  {
    DRT::Element* ele = embdis_->lColElement(cele);
    embele_outer_all.push_back(ele->Id());
  }

  // create node and element distribution of outer layer elements and nodes
  // of the embedded discretization ghosted on all processors
  Teuchos::RCP<const Epetra_Map> embnodeoutermap = Teuchos::rcp(new Epetra_Map(-1, embnode_outer_all.size(), &embnode_outer_all[0], 0, embdis_->Comm()));
  Teuchos::RCP<const Epetra_Map> embeleoutermap = Teuchos::rcp(new Epetra_Map(-1, embele_outer_all.size(), &embele_outer_all[0], 0, embdis_->Comm()));

  const Epetra_Map embnodeoutercolmap = *LINALG::AllreduceOverlappingEMap(*embnodeoutermap);
  const Epetra_Map embelemoutercolmap = *LINALG::AllreduceOverlappingEMap(*embeleoutermap);

  // redistribute nodes and elements to column (ghost) map
  embdis_->ExportColumnNodes(embnodeoutercolmap);
  embdis_->ExportColumnElements(embelemoutercolmap);

  embdis_->FillComplete();

#endif
}//FLD::XFluidFluid::PrepareEmbeddedDistribution()

// -------------------------------------------------------------------
// create boundary-embedded Map for Ale-sided-coupling
// -------------------------------------------------------------------
void FLD::XFluidFluid::CreateBoundaryEmbeddedMap()
{
  // fill boundary_embedded_mapdmap between boundary element id and its corresponding embedded element id
  for (int iele=0; iele< boundarydis_->NumMyColElements(); ++iele)
  {
    // boundary element and its nodes
    DRT::Element* bele = boundarydis_->lColElement(iele);
    const int* inodes = bele->NodeIds();

    bool bele_found = false;

    // ask all conditioned embedded elements for this boundary element
    for(int it=0; it< embdis_->NumMyColElements(); ++it)
    {
      DRT::Element* ele = embdis_->lColElement(it);
      const int* elenodes = (ele)->NodeIds();

      // assume the element has been founduntied
      bele_found = true;

      // check all nodes of the boundary element
      for(int inode=0; inode<bele->NumNode();  ++inode)
      {
        // boundary node
        const int inode_ID = inodes[inode];

        bool node_found = false;
        for (int enode=0; enode<ele->NumNode(); ++enode)
        {
          const int enode_ID = elenodes[enode];

          if(enode_ID == inode_ID)
          {
            node_found = true;
            break; // breaks the element nodes loop
          }
        }
        if(node_found==false) // this node is not contained in this element
        {
          bele_found = false; // element not the right one, if at least one boundary node is not found
          break; // node not found
        }
      }

      if(bele_found==true)
      {
        boundary_emb_gid_map_.insert(std::pair<int,int>(bele->Id(),ele->Id()));
        break;
      }

    }

    if(bele_found == false) dserror("corresponding embele for boundary element with boundary id %i not found on proc %i ! Please ghost corresponding embedded elements on all procs!", bele->Id(), myrank_);
  }

}//FLD::XFluidFluid::CreateBoundaryEmbeddedMap()

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::IntegrateFluidFluid()
{
  // distinguish stationary and instationary case
  if (timealgo_==INPAR::FLUID::timeint_stationary)
    SolveStationaryProblemFluidFluid();
  else
    TimeLoop();

  // print the results of time measurements
  Teuchos::TimeMonitor::summarize();

}

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::TimeLoop()
{
  TEUCHOS_FUNC_TIME_MONITOR("Time Loop");
  while (step_<stepmax_ and time_<maxtime_)
  {
    PrepareTimeStep();

    // -------------------------------------------------------------------
    //                       output to screen
    // -------------------------------------------------------------------
    if (myrank_==0)
    {
      switch (timealgo_)
      {
      case INPAR::FLUID::timeint_one_step_theta:
        IO::cout << "TIME: " << time_ << "/"<< maxtime_ << " " << "DT =  " << dta_ << "   One-Step-Theta    STEP =  "
        		 << step_<< "/"<< stepmax_ << IO::endl;
        break;
      case INPAR::FLUID::timeint_afgenalpha:
          IO::cout << "TIME: " << time_ << "/"<< maxtime_ << " " << "DT =  " << dta_ << "   Generalized-Alpha    STEP =  "
          		 << step_<< "/"<< stepmax_ << IO::endl;
        break;
      case INPAR::FLUID::timeint_bdf2:
          IO::cout << "TIME: " << time_ << "/"<< maxtime_ << " " << "DT =  " << dta_ << "   BDF2      STEP =  "
          		 << step_<< "/"<< stepmax_ << IO::endl;
        break;
      default:
        dserror("parameter out of range: IOP\n"); break;
      } /* end of switch(timealgo) */
    }

    // -----------------------------------------------------------------
    //        prepare nonlinear solve (used for Solve())
    // -----------------------------------------------------------------
    PrepareNonlinearSolve();

    // -----------------------------------------------------------------
    //                     solve nonlinear equation
    // -----------------------------------------------------------------
    Solve();

    // -------------------------------------------------------------------
    //                         update solution
    //        current solution becomes old solution of next timestep
    // -------------------------------------------------------------------
    TimeUpdate();


    // -------------------------------------------------------------------
    //  lift'n'drag forces, statistics time sample and output of solution
    //  and statistics
    // -------------------------------------------------------------------
    StatisticsAndOutput();

    // -------------------------------------------------------------------
    // evaluate error for test flows with analytical solutions
    // -------------------------------------------------------------------
    EvaluateErrorComparedToAnalyticalSol();

    // -------------------------------------------------------------------
    //                       update time step sizes
    // -------------------------------------------------------------------
    dtp_ = dta_;

    // -------------------------------------------------------------------
    //                    stop criterium for timeloop
    // -------------------------------------------------------------------
  }
}//FLD::XFluidFluid::IntegrateFluidFluid()

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::SolveStationaryProblemFluidFluid()
{
  // -------------------------------------------------------------------
  // pseudo time loop (continuation loop)
  // -------------------------------------------------------------------
  // slightly increasing b.c. values by given (pseudo-)timecurves to reach
  // convergence also for higher Reynolds number flows
  // as a side effect, you can do parameter studies for different Reynolds
  // numbers within only ONE simulation when you apply a proper
  // (pseudo-)timecurve

  while (step_< stepmax_)
  {
    // -------------------------------------------------------------------
    //              set (pseudo-)time dependent parameters
    // -------------------------------------------------------------------
    step_ += 1;
    time_ += dta_;
    // -------------------------------------------------------------------
    //                         out to screen
    // -------------------------------------------------------------------
    if (myrank_==0)
    {
    	IO::cout <<  "Stationary Fluid Solver - STEP = " << step_ << "/" << stepmax_ << IO::endl;
    }

    SetElementTimeParameter();

    SetDirichletNeumannBC();

    if (coupling_approach_ == CouplingNitsche_EmbFluid and nitsche_evp_)
      EstimateNitscheTraceMaxEigenvalue();

    // -------------------------------------------------------------------
    //                     solve nonlinear equation system
    // -------------------------------------------------------------------
    Solve();

    // -------------------------------------------------------------------
    //         calculate lift'n'drag forces from the residual
    // -------------------------------------------------------------------
    LiftDrag();

    // -------------------------------------------------------------------
    //                        compute flow rates
    // -------------------------------------------------------------------
    //ComputeFlowRates();

    // -------------------------------------------------------------------
    // evaluate error for test flows with analytical solutions
    // -------------------------------------------------------------------
    EvaluateErrorComparedToAnalyticalSol();


    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    Output();


  }
}// FLD::XFluidFluid::SolveStationaryProblemFluidFluid()

/*----------------------------------------------------------------------*
 |  check xfluid input parameters/ safety checks           schott 05/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluid::CheckXFluidFluidParams( Teuchos::ParameterList& params_xfem,
                                               Teuchos::ParameterList& params_xf_gen,
                                               Teuchos::ParameterList& params_xf_stab)
{
  if (myrank_==0)
  {
    // ----------------------------------------------------------------------
    // check XFEM GENERAL parameter list
    // ----------------------------------------------------------------------

    // ----------------------------------------------------------------------
    // check XFLUID DYNAMIC/GENERAL parameter list
    // ----------------------------------------------------------------------


    // ----------------------------------------------------------------------
    // check XFLUID DYNAMIC/STABILIZATION parameter list
    // ----------------------------------------------------------------------

    // convective stabilization parameter (scaling factor and stabilization factor)
    if (conv_stab_scaling_ != INPAR::XFEM::ConvStabScaling_none)
      dserror("INPAR::XFEM::ConvStabScaling should be set to ConvStabScaling_none for XFluidFluid!");

    if (velgrad_interface_stab_ and coupling_approach_ != CouplingNitsche_EmbFluid )
      dserror("VELGRAD_INTERFACE_STAB just for embedded sided Nitsche-Coupling!");

    if (presscoupling_interface_stab_ and coupling_approach_ != CouplingNitsche_EmbFluid)
      dserror("PRESSCOUPLING_INTERFACE_STAB just for embedded sided Nitsche-Coupling!");

    if (nitsche_evp_ and coupling_approach_ != CouplingNitsche_EmbFluid)
      dserror("NITSCHE_EVP just for embedded sided Nitsche-Coupling!");

  return;
  }
}


/*----------------------------------------------------------------------*
 |  Print fluid stabilization parameters                   schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluid::PrintStabilizationParams()
{
  // output of stabilization details
  if (myrank_==0)
  {
    Teuchos::ParameterList *  stabparams                =&(params_->sublist("RESIDUAL-BASED STABILIZATION"));
    Teuchos::ParameterList *  stabparams_edgebased      =&(params_->sublist("EDGE-BASED STABILIZATION"));
    Teuchos::ParameterList *  interfstabparams          =&(params_->sublist("XFLUID DYNAMIC/STABILIZATION"));


    IO::cout << "+------------------------------------------------------------------------------------+" << IO::endl;
    IO::cout << "                              FLUID-STABILIZATION                      \n " << IO::endl;

    IO::cout << "Stabilization type: " << stabparams->get<std::string>("STABTYPE") << "\n\n";

    //---------------------------------------------------------------------------------------------
    // output for residual-based fluid stabilization
    if(DRT::INPUT::IntegralValue<INPAR::FLUID::StabType>(*stabparams, "STABTYPE") == INPAR::FLUID::stabtype_residualbased)
    {
      IO::cout << "RESIDUAL-BASED fluid stabilization " << "\n";
      IO::cout << "                    " << stabparams->get<std::string>("TDS")<< "\n";
      IO::cout << "\n";

      std::string def_tau = stabparams->get<std::string>("DEFINITION_TAU");

//      if(    def_tau == "Franca_Barrenechea_Valentin_Frey_Wall"
//          or def_tau == "Franca_Barrenechea_Valentin_Frey_Wall_wo_dt") dserror("do not use Franca_Barrenechea_Valentin_Frey_Wall stabilization for XFEM -no stable results!");


      // instationary case
      if (timealgo_!=INPAR::FLUID::timeint_stationary)
      {
        IO::cout <<  "                    " << "Tau Type        = " << def_tau <<"\n";


        // check for instationary version of tau definitions
        if(def_tau != "Taylor_Hughes_Zarins" and
            def_tau != "Taylor_Hughes_Zarins_Whiting_Jansen" and
            def_tau != "Taylor_Hughes_Zarins_scaled" and
            def_tau != "Franca_Barrenechea_Valentin_Frey_Wall" and
            def_tau != "Shakib_Hughes_Codina" and
            def_tau != "Codina" and
            def_tau != "Franca_Madureira_Valentin_Badia_Codina" and
            def_tau != "Smoothed_FBVW")
        {
          IO::cout << RED_LIGHT
              << "Are you sure that you want to use stationary version of stabilization parameters "
              << "for instationary computations (just reasonable for small time steps dt)"
              << END_COLOR << IO::endl;
        }
      }
      else // stationary case
      {
        if(def_tau != "Taylor_Hughes_Zarins_wo_dt" and
            def_tau != "Taylor_Hughes_Zarins_Whiting_Jansen_wo_dt" and
            def_tau != "Taylor_Hughes_Zarins_scaled_wo_dt" and
            def_tau != "Franca_Barrenechea_Valentin_Frey_Wall_wo_dt" and
            def_tau != "Shakib_Hughes_Codina_wo_dt" and
            def_tau != "Codina_wo_dt" and
            def_tau != "Franca_Madureira_Valentin_Badia_Codina_wo_dt")
        {
//          dserror("not a valid tau definition (DEFINITION_TAU) for stationary problems");
        }
      }
      IO::cout << "\n";

      if(stabparams->get<std::string>("TDS") == "quasistatic")
      {
        if(stabparams->get<std::string>("TRANSIENT")=="yes_transient")
        {
          dserror("The quasistatic version of the residual-based stabilization currently does not support the incorporation of the transient term.");
        }
      }
      IO::cout <<  "                    " << "TRANSIENT       = " << stabparams->get<std::string>("TRANSIENT")      <<"\n";
      IO::cout <<  "                    " << "SUPG            = " << stabparams->get<std::string>("SUPG")           <<"\n";
      IO::cout <<  "                    " << "PSPG            = " << stabparams->get<std::string>("PSPG")           <<"\n";
      IO::cout <<  "                    " << "VSTAB           = " << stabparams->get<std::string>("VSTAB")          <<"\n";
      IO::cout <<  "                    " << "GRAD_DIV        = " << stabparams->get<std::string>("GRAD_DIV")       <<"\n";
      IO::cout <<  "                    " << "CROSS-STRESS    = " << stabparams->get<std::string>("CROSS-STRESS")   <<"\n";
      IO::cout <<  "                    " << "REYNOLDS-STRESS = " << stabparams->get<std::string>("REYNOLDS-STRESS")<<"\n";

      if(stabparams->get<std::string>("VSTAB")           != "no_vstab")    dserror("check VSTAB for XFEM");
      if(stabparams->get<std::string>("CROSS-STRESS")    != "no_cross")    dserror("check CROSS-STRESS for XFEM");
      if(stabparams->get<std::string>("REYNOLDS-STRESS") != "no_reynolds") dserror("check REYNOLDS-STRESS for XFEM");

    }
    else if(DRT::INPUT::IntegralValue<INPAR::FLUID::StabType>(*stabparams, "STABTYPE") == INPAR::FLUID::stabtype_edgebased)
    {
      // safety check for combinations of edge-based and residual-based stabilizations
      if((DRT::INPUT::IntegralValue<int>(*stabparams,"PSPG") != false)     or
         (DRT::INPUT::IntegralValue<int>(*stabparams,"SUPG") != false)     or
         (DRT::INPUT::IntegralValue<int>(*stabparams,"GRAD_DIV") != false)    or
         (stabparams->get<std::string>("VSTAB")          != "no_vstab")    or
         (stabparams->get<std::string>("CROSS-STRESS")   != "no_cross")    or
         (stabparams->get<std::string>("REYNOLDS-STRESS")!= "no_reynolds"))
         {
           dserror("if you want to combine residual-based stabilizations with edgebased-ghost-penalty stabilizations, please choose STABTYPE = residualbased");
         }
    }

    // check for non-valid combinations of residual-based and edge-based fluid stabilizations in the XFEM
    if( (DRT::INPUT::IntegralValue<int>(*stabparams,"PSPG") != false) and (stabparams_edgebased->get<std::string>("EOS_PRES") == "std_eos") )
      dserror("combine PSPG only with ghost-penalty variant of EOS_PRES ! ");
    if( (DRT::INPUT::IntegralValue<int>(*stabparams,"SUPG") != false) and (stabparams_edgebased->get<std::string>("EOS_CONV_STREAM") == "std_eos") )
      dserror("combine SUPG only with ghost-penalty variant of EOS_CONV_STREAM ! ");
    if( (DRT::INPUT::IntegralValue<int>(*stabparams,"SUPG") != false) and (stabparams_edgebased->get<std::string>("EOS_CONV_CROSS") == "std_eos") )
      dserror("combine SUPG only with ghost-penalty variant of EOS_CONV_CROSS ! ");

    //---------------------------------------------------------------------------------------------
    IO::cout << "\n\nEDGE-BASED (EOS) fluid stabilizations " << "\n";

    IO::cout <<  "                    " << "EOS_PRES             = " << stabparams_edgebased->get<std::string>("EOS_PRES")      <<"\n";
    IO::cout <<  "                    " << "EOS_CONV_STREAM      = " << stabparams_edgebased->get<std::string>("EOS_CONV_STREAM")      <<"\n";
    IO::cout <<  "                    " << "EOS_CONV_CROSS       = " << stabparams_edgebased->get<std::string>("EOS_CONV_CROSS")      <<"\n";
    IO::cout <<  "                    " << "EOS_DIV              = " << stabparams_edgebased->get<std::string>("EOS_DIV")      <<"\n";
    IO::cout <<  "                    " << "EOS_DEFINITION_TAU   = " << stabparams_edgebased->get<std::string>("EOS_DEFINITION_TAU")      <<"\n";
    IO::cout <<  "                    " << "EOS_H_DEFINITION     = " << stabparams_edgebased->get<std::string>("EOS_H_DEFINITION")      <<"\n";
    IO::cout << "+------------------------------------------------------------------------------------+" << IO::endl;
    IO::cout << "\n";


    //---------------------------------------------------------------------------------------------

    IO::cout << "+------------------------------------------------------------------------------------+" << IO::endl;
    IO::cout << "                              INTERFACE-STABILIZATION                       \n" << IO::endl;
    IO::cout << "Stabilization type          : " << interfstabparams->get<std::string>("COUPLING_METHOD") << "\n";
    IO::cout << "Coupling strategy           : " << interfstabparams->get<std::string>("COUPLING_STRATEGY") << "\n";

    if(coupling_method_ == INPAR::XFEM::Hybrid_LM_Cauchy_stress or coupling_method_ == INPAR::XFEM::Hybrid_LM_viscous_stress)
      IO::cout << "HYBRID_LM_L2_PROJ                 : " << interfstabparams->get<std::string>("HYBRID_LM_L2_PROJ") << "\n";

    IO::cout << "GHOST_PENALTY_STAB:           " << interfstabparams->get<std::string>("GHOST_PENALTY_STAB") << "\n";
    IO::cout << "GHOST_PENALTY_TRANSIENT_STAB: " << interfstabparams->get<std::string>("GHOST_PENALTY_TRANSIENT_STAB") << "\n";
    IO::cout << "GHOST_PENALTY_FAC:            " << interfstabparams->get<double>("GHOST_PENALTY_FAC") << "\n";
    IO::cout << "GHOST_PENALTY_TRANSIENT_FAC:  " << interfstabparams->get<double>("GHOST_PENALTY_TRANSIENT_FAC") << "\n";
    IO::cout << "GHOST_PENALTY_2nd_STAB:       " << interfstabparams->get<std::string>("GHOST_PENALTY_2nd_STAB") << "\n";
    IO::cout << "INFLOW_CONV_STAB_STRATEGY   : " << interfstabparams->get<std::string>("XFF_CONV_STAB_SCALING") << IO::endl;
    IO::cout << "VELGRAD_INTERFACE_STAB      : " << interfstabparams->get<std::string>("VELGRAD_INTERFACE_STAB")<< IO::endl;
    IO::cout << "PRESSCOUPLING_INTERFACE_STAB: " << interfstabparams->get<std::string>("PRESSCOUPLING_INTERFACE_STAB")<< IO::endl;
    IO::cout << "PRESSCOUPLING_INTERFACE_FAC : " << presscoupling_interface_fac_ << IO::endl;
    IO::cout << "VISC_STAB_FAC:              : " << interfstabparams->get<double>("VISC_STAB_FAC") << "\n";
    IO::cout << "VISC_STAB_TRACE_ESTIMATE:   : " << interfstabparams->get<std::string>("VISC_STAB_TRACE_ESTIMATE") << "\n";
    IO::cout << "VISC_STAB_HK                : " << interfstabparams->get<std::string>("VISC_STAB_HK")  << IO::endl;
    if (coupling_method_ != INPAR::XFEM::Hybrid_LM_viscous_stress)
      IO::cout << "VISC_ADJOINT_SYMMETRY       : " << interfstabparams->get<std::string>("VISC_ADJOINT_SYMMETRY") << IO::endl;
    IO::cout << IO::endl;

    IO::cout << "+------------------------------------------------------------------------------------+" << IO::endl;
    IO::cout << "\n";

  }

}


// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::PrepareTimeStep()
{
  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  step_ += 1;
  time_ += dta_;

  gmsh_count_ = 0;

  // -------------------------------------------------------------------
  // set time parameters dependent on time integration scheme and step
  // -------------------------------------------------------------------
  if (timealgo_ == INPAR::FLUID::timeint_stationary)
  {
    theta_ = 1.0;
  }
  else
  {
    // do a backward Euler step for the first timestep
    if (step_==1)
    {
      theta_ = params_->get<double>("start theta");
    }
    else if (step_ > 1)
    {
      // for OST
      if(timealgo_ == INPAR::FLUID::timeint_one_step_theta) theta_ = params_->get<double>("theta");

      // for BDF2, theta is set by the time-step sizes, 2/3 for const. dt
      if (timealgo_==INPAR::FLUID::timeint_bdf2) theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);
    }
    else dserror("number of time step is wrong");
  }

  // -------------------------------------------------------------------
  //  Set time parameter for element call
  // -------------------------------------------------------------------
  SetElementTimeParameter();

  if (monolithicfluidfluidfsi_)
  {
    SetHistoryValues();

    SetDirichletNeumannBC();
  }

  if (coupling_approach_ == CouplingNitsche_EmbFluid and nitsche_evp_)
    EstimateNitscheTraceMaxEigenvalue();

}//FLD::XFluidFluid::PrepareTimeStep()

// ----------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::PrepareNonlinearSolve()
{
  // do the cut for this timestep, if the embedded fluid is an ALE-fluid
  if (alefluid_)
    CutAndSaveBgFluidStatus();

  if (not monolithicfluidfluidfsi_)
    SetBgStateVectors(aledispn_);
  else if(monotype_ == FixedALEPartitioned or monotype_ == FixedALEInterpolation)
    SetBgStateVectors(aledispnpoldstate_);
  else if(monotype_ == FullyNewton)
    SetBgStateVectors(aledispnp_);
  else
    dserror("Unknown monolithic approach.");

  SetHistoryValues();
  SetDirichletNeumannBC();

}//FLD::XFluidFluid::PrepareNonlinearSolve()

// ----------------------------------------------------------------
// Prepare monolithic step (called in TimeUpdate)
// - set time parameters
// - do the cut and xfluidfluid time integration
// - do the nonlinearsolve with fsi-dofs as dirichlet values
// -------------------------------------------------------------------
void FLD::XFluidFluid::PrepareMonolithicFixedAle()
{
  if (bgdis_->Comm().MyPID() == 0)
    IO::cout << "Update monolithic fluid solution " << IO::endl;

  // for BDF2, theta is set by the time-step sizes, 2/3 for const. dtp_
  if (timealgo_==INPAR::FLUID::timeint_bdf2) theta_ = (dta_+dtp_)/(2.0*dta_ + dtp_);

  // -------------------------------------------------------------------
  //  Set time parameter for element call
  // -------------------------------------------------------------------
  SetElementTimeParameter();

  // cut and do xfluidfluid time integration.
  PrepareNonlinearSolve();

  if(monotype_ == FixedALEPartitioned)
    UpdateMonolithicFluidSolution();
}

// ----------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::Solve()
{
  // ---------------------------------------------- nonlinear iteration
  // ------------------------------- stop nonlinear iteration when both
  //                                 increment-norms are below this bound
  const double  ittol     =params_->get<double>("tolerance for nonlin iter");

  //------------------------------ turn adaptive solver tolerance on/off
  const bool   isadapttol    = params_->get<bool>("ADAPTCONV",true);
  const double adaptolbetter = params_->get<double>("ADAPTCONV_BETTER",0.01);

  int  itnum = 0;
  bool stopnonliniter = false;

  int itemax  = params_->get<int>   ("max nonlin iter steps");

  dtsolve_  = 0.0;
  dtele_    = 0.0;
  dtfilter_ = 0.0;

  if (myrank_ == 0)
  {
    IO::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+" << IO::endl
             << "|- step/max -|- tol      [norm] -|-- vel-res ---|-- pre-res ---|-- vel-inc ---|-- pre-inc ---|" << IO::endl;
  }

//    const int numcolele = bgdis_->NumMyColElements();
//    for (int j=0; j<numcolele; ++j)
//    {
//      DRT::Element* actele = bgdis_->lColElement(j);
//      cout << "actele " << actele->Id() << endl;
//      const DRT::Node*const* elenodes = actele->Nodes();
//      for (int inode=0; inode<actele->NumNode(); ++inode)
//      {
//        cout << " node ids: " << elenodes[inode]->Id()  ;
//        std::vector<int> gdofs = bgdis_->Dof(elenodes[inode]);
//        cout << " dofs: " ;
//        for (int d=0; d<gdofs.size();++d)
//          cout << " " << gdofs.at(d) << " ";
//        cout << endl;
//      }
//    }

  while (stopnonliniter==false)
  {
    // Insert fluid and xfluid vectors to fluidxfluid
    state_->fluidfluidsplitter_->InsertXFluidVector(state_->velnp_,state_->fluidfluidvelnp_);
    state_->fluidfluidsplitter_->InsertFluidVector(alevelnp_,state_->fluidfluidvelnp_);

    itnum++;

    // -------------------------------------------------------------------
    // Call elements to calculate system matrix and RHS
    // -------------------------------------------------------------------
    {
      // get cpu time
      const double tcpu=Teuchos::Time::wallTime();

      // create the parameters for the discretization
      Teuchos::ParameterList eleparams;

      // set vector values needed by elements
      bgdis_->ClearState();
      bgdis_->SetState("velaf",state_->velnp_);

      embdis_->ClearState();
      embdis_->SetState("velaf",alevelnp_);


      int itemax  = params_->get<int>("max nonlin iter steps");

      if (itnum != itemax)
      {
        state_->EvaluateFluidFluid( bgdis_,
                                    boundarydis_,
                                    embdis_);
      }

      // end time measurement for element
      dtele_=Teuchos::Time::wallTime()-tcpu;

    }


    // scaling to get true residual vector
    state_->trueresidual_->Update(ResidualScaling(),*state_->residual_,0.0);
    aletrueresidual_->Update(ResidualScaling(),*aleresidual_,0.0);

    // blank residual DOFs which are on Dirichlet BC
    // We can do this because the values at the dirichlet positions
    // are not used anyway.
    // We could avoid this though, if velrowmap_ and prerowmap_ would
    // not include the dirichlet values as well. But it is expensive
    // to avoid that.

    state_->dbcmaps_->InsertCondVector(state_->dbcmaps_->ExtractCondVector(state_->zeros_), state_->residual_);
    aledbcmaps_->InsertCondVector(aledbcmaps_->ExtractCondVector(alezeros_), aleresidual_);

    if(monotype_ == FixedALEPartitioned)
    {
      // set the aleresidual values to zeros at the fsi-interface
      Teuchos::RCP<Epetra_Vector> fixedfsizeros = LINALG::CreateVector(*fixedfsidofmap_,true);
      LINALG::Export(*fixedfsizeros,*aleresidual_);
    }

    // insert fluid and alefluid residuals to fluidfluidresidual
    state_->fluidfluidsplitter_->InsertXFluidVector(state_->residual_,state_->fluidfluidresidual_);
    state_->fluidfluidsplitter_->InsertFluidVector(aleresidual_,state_->fluidfluidresidual_);

    double incvelnorm_L2 = 0.0;
    double incprenorm_L2 = 0.0;

    double velnorm_L2 = 0.0;
    double prenorm_L2 = 0.0;

    double vresnorm = 0.0;
    double presnorm = 0.0;


    Teuchos::RCP<Epetra_Vector> onlyvel = state_->fluidfluidvelpressplitter_.ExtractOtherVector(state_->fluidfluidresidual_);
    onlyvel->Norm2(&vresnorm);

    state_->fluidfluidvelpressplitter_.ExtractOtherVector(state_->fluidfluidincvel_,onlyvel);
    onlyvel->Norm2(&incvelnorm_L2);

    state_->fluidfluidvelpressplitter_.ExtractOtherVector(state_->fluidfluidvelnp_,onlyvel);;
    onlyvel->Norm2(&velnorm_L2);


    Teuchos::RCP<Epetra_Vector> onlypre = state_->fluidfluidvelpressplitter_.ExtractCondVector(state_->fluidfluidresidual_);
    onlypre->Norm2(&presnorm);

    state_->fluidfluidvelpressplitter_.ExtractCondVector(state_->fluidfluidincvel_,onlypre);
    onlypre->Norm2(&incprenorm_L2);

    state_->fluidfluidvelpressplitter_.ExtractCondVector(state_->fluidfluidvelnp_,onlypre);
    onlypre->Norm2(&prenorm_L2);


    // care for the case that nothing really happens in the velocity
    // or pressure field
    if (velnorm_L2 < 1e-5) velnorm_L2 = 1.0;
    if (prenorm_L2 < 1e-5) prenorm_L2 = 1.0;

    //-------------------------------------------------- output to screen
    /* special case of very first iteration step:
        - solution increment is not yet available
        - convergence check is not required (we solve at least once!)    */
    if (itnum == 1)
    {
      double min = 1.0e19;
      double max = 0.0;
      bgdis_->Comm().MinAll(&dtele_,&min,1);
      bgdis_->Comm().MaxAll(&dtele_,&max,1);
      if (myrank_ == 0)
      {
        IO::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itemax_ << "   | "
        		 << std::setw(10) << std::setprecision(3) << std::scientific << ittol << "[L_2 ]  | "
                 << std::setw(10) << std::setprecision(3) << std::scientific << vresnorm << "   | "
                 << std::setw(10) << std::setprecision(3) << std::scientific << presnorm << "   |      --      |      --      | (      --     ,te_min="
                 << std::setw(10) << std::setprecision(3) << std::scientific << min << ",te_max="
                 << std::setw(10) << std::setprecision(3) << std::scientific << max << ")";
        IO::cout << IO::endl;
      }
    }
    /* ordinary case later iteration steps:
        - solution increment can be printed
        - convergence check should be done*/
    else
    {
    // this is the convergence check
    // We always require at least one solve. Otherwise the
    // perturbation at the FSI interface might get by unnoticed.
      if (vresnorm <= ittol and presnorm <= ittol and
          incvelnorm_L2/velnorm_L2 <= ittol and incprenorm_L2/prenorm_L2 <= ittol)
      {
        stopnonliniter=true;
        double min = 1.0e19;
        double max = 0.0;
        if (myrank_ == 0)
        {
          IO::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itemax_ << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << ittol << "[L_2 ]  | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << vresnorm << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << presnorm << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << incvelnorm_L2/velnorm_L2 << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << incprenorm_L2/prenorm_L2 << "   | (ts="
                   << std::setw(10) << std::setprecision(3) << std::scientific << dtsolve_ << ",te_min="
                   << std::setw(10) << std::setprecision(3) << std::scientific << min << ",te_max="
                   << std::setw(10) << std::setprecision(3) << std::scientific << max << ")" << IO::endl;
          IO::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+" << IO::endl;

          FILE* errfile = params_->get<FILE*>("err file",NULL);
          if (errfile!=NULL)
          {
            fprintf(errfile,"fluid solve:   %3d/%3d  tol=%10.3E[L_2 ]  vres=%10.3E  pres=%10.3E  vinc=%10.3E  pinc=%10.3E\n",
                    itnum,itemax,ittol,vresnorm,presnorm,
                    incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
          }
        }
        break;
      }
      else // if not yet converged
      {
    	double min = 1.0e19;
        double max = 0.0;
        bgdis_->Comm().MinAll(&dtele_,&min,1);
        bgdis_->Comm().MaxAll(&dtele_,&max,1);

        if (myrank_ == 0)
        {
          IO::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itemax_ << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << ittol << "[L_2 ]  | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << vresnorm << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << presnorm << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << incvelnorm_L2/velnorm_L2 << "   | "
                   << std::setw(10) << std::setprecision(3) << std::scientific << incprenorm_L2/prenorm_L2 << "   | (ts="
                   << std::setw(10) << std::setprecision(3) << std::scientific << dtsolve_ << ",te_min="
                   << std::setw(10) << std::setprecision(3) << std::scientific << min << ",te_max="
                   << std::setw(10) << std::setprecision(3) << std::scientific << max << ")" << IO::endl;
        }
      }
    }

    // warn if itemax is reached without convergence, but proceed to
    // next timestep...
    if ((itnum == itemax) and (vresnorm > ittol or presnorm > ittol or
                             incvelnorm_L2/velnorm_L2 > ittol or
                             incprenorm_L2/prenorm_L2 > ittol))
    {
      stopnonliniter=true;
      if (myrank_ == 0)
      {
        printf("+---------------------------------------------------------------+\n");
        printf("|            >>>>>> not converged in itemax steps!              |\n");
        printf("+---------------------------------------------------------------+\n");

        FILE* errfile = params_->get<FILE*>("err file",NULL);
        if (errfile!=NULL)
        {
          fprintf(errfile,"fluid unconverged solve:   %3d/%3d  tol=%10.3E[L_2 ]  vres=%10.3E  pres=%10.3E  vinc=%10.3E  pinc=%10.3E\n",
                  itnum,itemax,ittol,vresnorm,presnorm,
                  incvelnorm_L2/velnorm_L2,incprenorm_L2/prenorm_L2);
        }
      }
      break;
    }

    //--------- Apply Dirichlet boundary conditions to system of equations
    //          residual displacements are supposed to be zero at
    //          boundary conditions
    state_->fluidfluidincvel_->PutScalar(0.0);

    // Add the fluid & xfluid & couple-matrices to fluidxfluidsysmat;
    state_->fluidfluidsysmat_->Zero();
    state_->fluidfluidsysmat_->Add(*state_->sysmat_,false,1.0,0.0);
    state_->fluidfluidsysmat_->Add(*alesysmat_,false,1.0,1.0);
    state_->fluidfluidsysmat_->Add(*state_->Cuui_,false,1.0,1.0);
    state_->fluidfluidsysmat_->Add(*state_->Cuiu_,false,1.0,1.0);
    state_->fluidfluidsysmat_->Add(*state_->Cuiui_,false,1.0,1.0);
    state_->fluidfluidsysmat_->Complete();

    // build a merged map from fluid-fluid dirichlet-maps
    state_->CreateFluidFluidDBCMaps();

    LINALG::ApplyDirichlettoSystem(state_->fluidfluidsysmat_,state_->fluidfluidincvel_,state_->fluidfluidresidual_,state_->fluidfluidzeros_,*state_->fluidfluiddbcmaps_);

    // set the fsi dirichlet values for monolithic_fixedale_partitioned
    if (monotype_ == FixedALEPartitioned)
    {
      LINALG::ApplyDirichlettoSystem(state_->fluidfluidsysmat_,state_->fluidfluidincvel_,state_->fluidfluidresidual_,state_->fluidfluidzeros_,toggle_);
    }

    //-------solve for residual displacements to correct incremental displacements
    {
      // get cpu time
      const double tcpusolve=Teuchos::Time::wallTime();

      // do adaptive linear solver tolerance (not in first solve)
      if (isadapttol && itnum>1)
      {
        double currresidual = std::max(vresnorm,presnorm);
        currresidual = std::max(currresidual,incvelnorm_L2/velnorm_L2);
        currresidual = std::max(currresidual,incprenorm_L2/prenorm_L2);
        solver_->AdaptTolerance(ittol,currresidual,adaptolbetter);
      }


      Teuchos::RCP<LINALG::SparseMatrix> sysmatmatrixmatlab = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(state_->fluidfluidsysmat_);
      solver_->Solve(state_->fluidfluidsysmat_->EpetraOperator(),state_->fluidfluidincvel_,state_->fluidfluidresidual_,true,itnum==1);
      solver_->ResetTolerance();

      // end time measurement for solver
      dtsolve_ = Teuchos::Time::wallTime()-tcpusolve;
    }

    // -------------------------------------------------------------------
    // update velocity and pressure values by increments
    //
    // -------------------------------------------------------------------
    state_->fluidfluidvelnp_->Update(1.0,*state_->fluidfluidincvel_,1.0);
    // extract velnp_
    state_->velnp_ = state_->fluidfluidsplitter_->ExtractXFluidVector(state_->fluidfluidvelnp_);
    alevelnp_ = state_->fluidfluidsplitter_->ExtractFluidVector(state_->fluidfluidvelnp_);

    // extract residual
    state_->residual_ = state_->fluidfluidsplitter_->ExtractXFluidVector(state_->fluidfluidresidual_);
    aleresidual_ = state_->fluidfluidsplitter_->ExtractFluidVector(state_->fluidfluidresidual_);

    // Update the fluid material velocity along the interface (ivelnp_), source (in): state_.alevelnp_
    LINALG::Export(*(alevelnp_),*(ivelnp_));
    boundarydis_->SetState("ivelnp",ivelnp_);

    // -------------------------------------------------------------------
    // For af-generalized-alpha: update accelerations
    // Furthermore, calculate velocities, pressures, scalars and
    // accelerations at intermediate time steps n+alpha_F and n+alpha_M,
    // respectively, for next iteration.
    // This has to be done at the end of the iteration, since we might
    // need the velocities at n+alpha_F in a potential coupling
    // algorithm, for instance.
    // -------------------------------------------------------------------
    if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
    {
      GenAlphaUpdateAcceleration();

      GenAlphaIntermediateValues();
    }
  }

  if (alefluid_)
    aletotaldispn_->Update(1.0,*aledispn_,1.0);

}//FLD::XFluidFluid::Solve()

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::LinearSolve()
{

}

// -------------------------------------------------------------------
// evaluate method for monolithic fluid-fluid-fsi
// -------------------------------------------------------------------
void FLD::XFluidFluid::Evaluate(
  Teuchos::RCP<const Epetra_Vector> stepinc ///< solution increment between time step n and n+1
  )
{
  // set embedded fluid system matrix to 0
  alesysmat_->Zero();

  // set shapederivatives to 0, if activated
  if (shapederivatives_ != Teuchos::null)
    shapederivatives_->Zero();

  gmsh_count_++;

  // set the new solution we just got. Note: the solution we got here
  // is the step increment which means the sum of all iteration
  // increments of the time step.
  if (stepinc!=Teuchos::null)
  {
    // Take Dirichlet values from velnp and add vel to veln for non-Dirichlet
    // values.

    Teuchos::RCP<Epetra_Vector> stepinc_bg = LINALG::CreateVector(*state_->fluiddofrowmap_,true);
    Teuchos::RCP<Epetra_Vector> stepinc_emb = LINALG::CreateVector(*aledofrowmap_,true);

    stepinc_bg = state_->fluidfluidsplitter_->ExtractXFluidVector(stepinc);
    stepinc_emb = state_->fluidfluidsplitter_->ExtractFluidVector(stepinc);

    Teuchos::RCP<Epetra_Vector> stepinc_bg_tmp = LINALG::CreateVector(*state_->fluiddofrowmap_,true);
    Teuchos::RCP<Epetra_Vector> stepinc_emb_tmp = LINALG::CreateVector(*aledofrowmap_,true);

    stepinc_bg_tmp->Update(1.0, *state_->veln_, 1.0, *stepinc_bg, 0.0);
    stepinc_emb_tmp->Update(1.0, *aleveln_, 1.0, *stepinc_emb, 0.0);

    // Apply DBCs to stepinc
    state_->dbcmaps_->InsertCondVector(state_->dbcmaps_->ExtractCondVector(state_->velnp_), stepinc_bg_tmp );
    aledbcmaps_->InsertCondVector(aledbcmaps_->ExtractCondVector(alevelnp_), stepinc_emb_tmp );

    *state_->velnp_ = *stepinc_bg_tmp;
    *alevelnp_ = *stepinc_emb_tmp;

    // prepare new iteration for fully_newton approach
    if(monotype_ == FullyNewton)
    {
      // cut and set state vectors
      PrepareNonlinearSolve();

      stepinc_bg = LINALG::CreateVector(*state_->fluiddofrowmap_,true);

      // calculate the stepinc for both fluids. This is needed in
      // monolithic fsi to sum up the iteration steps
      stepinc_bg->Update(1.0,*state_->velnp_,-1.0,*state_->veln_,0.0);
      stepinc_emb->Update(1.0,*alevelnp_,-1.0,*aleveln_,0.0);

      state_->fluidfluidsplitter_->InsertXFluidVector(stepinc_bg,state_->stepinc_);
      state_->fluidfluidsplitter_->InsertFluidVector(stepinc_emb,state_->stepinc_);
    }

    state_->fluidfluidsplitter_->InsertXFluidVector(state_->velnp_,state_->fluidfluidvelnp_);
    state_->fluidfluidsplitter_->InsertFluidVector(alevelnp_,state_->fluidfluidvelnp_);

    // mit iterinc--------------------------------------
//     state_->fluidfluidvelnp_->Update(1.0,*stepinc,1.0);
//     // extract velnp_
//     state_->velnp_ = state_->fluidfluidsplitter_->ExtractXFluidVector(state_->fluidfluidvelnp_);
//     alevelnp_ = state_->fluidfluidsplitter_->ExtractFluidVector(state_->fluidfluidvelnp_);

    // extract residual
    state_->residual_ = state_->fluidfluidsplitter_->ExtractXFluidVector(state_->fluidfluidresidual_);
    aleresidual_ = state_->fluidfluidsplitter_->ExtractFluidVector(state_->fluidfluidresidual_);

    // Update the fluid material velocity along the interface (ivelnp_), source (in): state_.alevelnp_
    LINALG::Export(*(alevelnp_),*(ivelnp_));
  }

  state_->EvaluateFluidFluid( bgdis_,
                              boundarydis_,
                              embdis_);

  // scaling to get true residual vector
  state_->trueresidual_->Update(ResidualScaling(),*state_->residual_,0.0);
  aletrueresidual_->Update(ResidualScaling(),*aleresidual_,0.0);

  // Add the fluid & xfluid & couple-matrices to fluidxfluidsysmat
  state_->fluidfluidsysmat_->Add(*state_->sysmat_,false,1.0,0.0);
  state_->fluidfluidsysmat_->Add(*alesysmat_,false,1.0,1.0);
  state_->fluidfluidsysmat_->Add(*state_->Cuui_,false,1.0,1.0);
  state_->fluidfluidsysmat_->Add(*state_->Cuiu_,false,1.0,1.0);
  state_->fluidfluidsysmat_->Add(*state_->Cuiui_,false,1.0,1.0);
  state_->fluidfluidsysmat_->Complete();

  if (shapederivatives_ != Teuchos::null)
  {
    shapederivatives_->Complete();
    // apply Dirichlet conditions to a non-diagonal matrix
    // (The Dirichlet rows will become all zero, no diagonal one.)
    shapederivatives_->ApplyDirichlet(*(aledbcmaps_->CondMap()),false);
  }

  state_->fluidfluidincvel_->PutScalar(0.0);

  state_->dbcmaps_->InsertCondVector(state_->dbcmaps_->ExtractCondVector(state_->zeros_), state_->residual_);
  aledbcmaps_->InsertCondVector(aledbcmaps_->ExtractCondVector(alezeros_), aleresidual_);

  // insert fluid and alefluid residuals to fluidfluidresidual
  state_->fluidfluidsplitter_->InsertXFluidVector(state_->residual_,state_->fluidfluidresidual_);
  state_->fluidfluidsplitter_->InsertFluidVector(aleresidual_,state_->fluidfluidresidual_);

  // build a merged map from fluid-fluid dbc-maps
  state_->CreateFluidFluidDBCMaps();

  LINALG::ApplyDirichlettoSystem(state_->fluidfluidsysmat_,state_->fluidfluidincvel_,state_->fluidfluidresidual_,
                                 state_->fluidfluidzeros_,*state_->fluidfluiddbcmaps_);

  if(gmsh_debug_out_)
  {
    if(monotype_ == FullyNewton)
      state_->GmshOutput(*bgdis_,*embdis_,*boundarydis_, "result_evaluate", gmsh_count_, step_, state_->velnp_, alevelnp_,
                         aledispnp_);
    else
      state_->GmshOutput(*bgdis_,*embdis_,*boundarydis_, "res_evaluate", -1, step_, state_->velnp_, alevelnp_,
                         aledispn_);
  }

}//FLD::XFluidFluid::Evaluate

// -------------------------------------------------------------------
// Read Restart data for background discretization
// -------------------------------------------------------------------
void FLD::XFluidFluid::ReadRestart(int step)
{
  //-------- background discretization
  //  ART_exp_timeInt_->ReadRestart(step);
  IO::DiscretizationReader reader(bgdis_,step);

  reader.ReadVector(state_->velnp_,"velnp_bg");
  reader.ReadVector(state_->velnm_,"velnm_bg");
  reader.ReadVector(state_->veln_,"veln_bg");
  reader.ReadVector(state_->accnp_,"accnp_bg");
  reader.ReadVector(state_->accn_ ,"accn_bg");

  // set element time parameter after restart:
  // Here it is already needed by AVM3 and impedance boundary condition!!
  SetElementTimeParameter();

  // ensure that the overall dof numbering is identical to the one
  // that was used when the restart data was written. Especially
  // in case of multiphysics problems & periodic boundary conditions
  // it is better to check the consistency of the maps here:
  if (not (bgdis_->DofRowMap())->SameAs(state_->velnp_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (bgdis_->DofRowMap())->SameAs(state_->veln_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (bgdis_->DofRowMap())->SameAs(state_->accn_->Map()))
    dserror("Global dof numbering in maps does not match");
}

// -------------------------------------------------------------------
// Read Restart data for embedded discretization
// -------------------------------------------------------------------
void FLD::XFluidFluid::ReadRestartEmb(int step)
{
  //-------- embedded discretization
  IO::DiscretizationReader embreader(embdis_,step);
  time_ = embreader.ReadDouble("time");
  step_ = embreader.ReadInt("step");

  embreader.ReadVector(alevelnp_,"velnp_emb");
  embreader.ReadVector(aleveln_, "veln_emb");
  embreader.ReadVector(alevelnm_,"velnm_emb");
  embreader.ReadVector(aleaccnp_,"accnp_emb");
  embreader.ReadVector(aleaccn_ ,"accn_emb");
  // set element time parameter after restart:
  // Here it is already needed by AVM3 and impedance boundary condition!!
  SetElementTimeParameter();

  if (alefluid_)
  {
    embreader.ReadVector(aledispnp_,"dispnp_emb");
    embreader.ReadVector(aledispn_ , "dispn_emb");
    embreader.ReadVector(aledispnm_,"dispnm_emb");
    embreader.ReadVector(gridv_,"gridv_emb");
  }

  // ensure that the overall dof numbering is identical to the one
  // that was used when the restart data was written. Especially
  // in case of multiphysics problems & periodic boundary conditions
  // it is better to check the consistency of the maps here:
  if (not (embdis_->DofRowMap())->SameAs(alevelnp_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (embdis_->DofRowMap())->SameAs(aleveln_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (embdis_->DofRowMap())->SameAs(aleaccn_->Map()))
    dserror("Global dof numbering in maps does not match");
}

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::UpdateGridv()
{
  // get order of accuracy of grid velocity determination
  // from input file data (updated to BE, BDF2 inpars rauch 09/13)
  const Teuchos::ParameterList& fluiddynparams =  DRT::Problem::Instance()->FluidDynamicParams();
  const int order = DRT::INPUT::IntegralValue<INPAR::FLUID::Gridvel>(fluiddynparams, "GRIDVEL");

  switch (order)
  {
    case INPAR::FLUID::BE:
      /* get gridvelocity from BE time discretisation of mesh motion:
           -> cheap
           -> easy
           -> limits FSI algorithm to first order accuracy in time

                  x^n+1 - x^n
             uG = -----------
                    Delta t                        */
      gridv_->Update(1/dta_, *aledispnp_, -1/dta_, *aledispn_, 0.0);
    break;
    case INPAR::FLUID::BDF2:
      /* get gridvelocity from BDF2 time discretisation of mesh motion:
           -> requires one more previous mesh position or displacement
           -> somewhat more complicated
           -> allows second order accuracy for the overall flow solution  */
      gridv_->Update(1.5/dta_, *aledispnp_, -2.0/dta_, *aledispn_, 0.0);
      gridv_->Update(0.5/dta_, *aledispnm_, 1.0);
    break;
    case INPAR::FLUID::OST:
      /* get gridvelocity from OST time discretisation of mesh motion:
         -> needed to allow consistent linearization of FPSI problem  */
      dserror("One-Step-Theta gridvelocity determination not implemented for xfluids. fix your input file! GRIDVEL = BE or BDF2 ");
    break;
    default:
      dserror("Unknown or invalid type of grid velocity determination. Fix GRIDVEL section of your input file.");
      break;
  }
}//FLD::XFluidFluid::UpdateGridv()

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd)
{
  std::vector<Teuchos::RCP<const Epetra_Map> > condmaps;
  condmaps.push_back(maptoadd);
  condmaps.push_back(aledbcmaps_->CondMap());
  Teuchos::RCP<Epetra_Map> condmerged = LINALG::MultiMapExtractor::MergeMaps(condmaps);
  *(aledbcmaps_) = LINALG::MapExtractor(*(embdis_->DofRowMap()), condmerged);

  // update the combined fluid dirichlet maps immediately
  state_->CreateFluidFluidDBCMaps();
  return;
}//FLD::XFluidFluid::AddDirichCond

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::TimeUpdate()
{
  // parameter for monolithic_approach_fixedale, do not move
  // the embedded mesh in this monolithic step. (No Cut and ale-relaxing;
  // just update the solution)
  bool RelaxingAleInthisTimestep = relaxing_ale_;

  if (step_%relaxing_ale_every_!=0)
    RelaxingAleInthisTimestep = false;

  if ((monotype_ == FixedALEPartitioned or monotype_ == FixedALEInterpolation)
      and RelaxingAleInthisTimestep)
  {
    PrepareMonolithicFixedAle();
  }

  Teuchos::ParameterList *  stabparams=&(params_->sublist("RESIDUAL-BASED STABILIZATION"));

  if(stabparams->get<string>("TDS") == "time_dependent")
  {
    const double tcpu=Teuchos::Time::wallTime();

    if(myrank_==0)
    {
      IO::cout << "time update for subscales";
    }

    // call elements to calculate system matrix and rhs and assemble
    // this is required for the time update of the subgrid scales and
    // makes sure that the current subgrid scales correspond to the
    // current residual
    AssembleMatAndRHS();

    // create the parameters for the discretization
    Teuchos::ParameterList eleparams;

    // update time paramters
    if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
    {
      eleparams.set("gamma"  ,gamma_);
    }
    else if (timealgo_==INPAR::FLUID::timeint_one_step_theta)
    {
      eleparams.set("gamma"  ,theta_);
    }
    else if (timealgo_==INPAR::FLUID::timeint_bdf2)
    {
      eleparams.set("gamma"  ,1.0);
    }
    else
    {
      dserror("Unknown time integration approach.");
    }

    eleparams.set("dt"     ,dta_    );

    // call loop over elements to update subgrid scales
    bgdis_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
    embdis_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

    if(myrank_==0)
    {
      IO::cout << "("<<Teuchos::Time::wallTime()-tcpu<<")\n";
    }
  }

  // Compute accelerations
  {
    Teuchos::RCP<Epetra_Vector> onlyaccn  = state_->velpressplitter_.ExtractOtherVector(state_->accn_ );
    Teuchos::RCP<Epetra_Vector> onlyaccnp = state_->velpressplitter_.ExtractOtherVector(state_->accnp_);
    Teuchos::RCP<Epetra_Vector> onlyvelnm = state_->velpressplitter_.ExtractOtherVector(state_->velnm_);
    Teuchos::RCP<Epetra_Vector> onlyveln  = state_->velpressplitter_.ExtractOtherVector(state_->veln_ );
    Teuchos::RCP<Epetra_Vector> onlyvelnp = state_->velpressplitter_.ExtractOtherVector(state_->velnp_);

    COMBUST::UTILS::CalculateAcceleration(onlyvelnp,
                                              onlyveln ,
                                              onlyvelnm,
                                              onlyaccn ,
                                              timealgo_,
                                              step_    ,
                                              theta_   ,
                                              dta_     ,
                                              dtp_     ,
                                              onlyaccnp);

    // copy back into global vector
    LINALG::Export(*onlyaccnp,*state_->accnp_);

    Teuchos::RCP<Epetra_Vector> aleonlyaccn  = alevelpressplitter_.ExtractOtherVector(aleaccn_ );
    Teuchos::RCP<Epetra_Vector> aleonlyaccnp = alevelpressplitter_.ExtractOtherVector(aleaccnp_);
    Teuchos::RCP<Epetra_Vector> aleonlyvelnm = alevelpressplitter_.ExtractOtherVector(alevelnm_);
    Teuchos::RCP<Epetra_Vector> aleonlyveln  = alevelpressplitter_.ExtractOtherVector(aleveln_ );
    Teuchos::RCP<Epetra_Vector> aleonlyvelnp = alevelpressplitter_.ExtractOtherVector(alevelnp_);

    COMBUST::UTILS::CalculateAcceleration(aleonlyvelnp,
                                              aleonlyveln ,
                                              aleonlyvelnm,
                                              aleonlyaccn ,
                                              timealgo_,
                                              step_    ,
                                              theta_   ,
                                              dta_     ,
                                              dtp_     ,
                                              aleonlyaccnp);

    // copy back into global vector
    LINALG::Export(*aleonlyaccnp,*aleaccnp_);
  }

  // update old acceleration
  state_->accn_->Update(1.0,*state_->accnp_,0.0);
  aleaccn_->Update(1.0,*aleaccnp_,0.0);

  // velocities/pressures of this step become most recent
  // velocities/pressures of the last step
  state_->velnm_->Update(1.0,*state_->veln_ ,0.0);
  state_->veln_ ->Update(1.0,*state_->velnp_,0.0);

  alevelnm_->Update(1.0,*aleveln_ ,0.0);
  aleveln_ ->Update(1.0,*alevelnp_,0.0);

  if (alefluid_)
  {
    aledispnm_->Update(1.0,*aledispn_,0.0);
    aledispn_->Update(1.0,*aledispnp_,0.0);
  }

} //XFluidFluid::TimeUpdate()

// -------------------------------------------------------------------
// In this function:
//
// - Save the old state_ and old status of bg nodes (std/enriched/void)
// - New cut with the new ale displacement
// - Save the status of new bg nodes
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::CutAndSaveBgFluidStatus()
{
  // save the old state vector
  staten_ = state_;

  const int restart = DRT::Problem::Instance()->Restart();
  // if restart
  if (restart and ((restart+1) == step_))
  {
    xfluidfluid_timeint_->CreateBgNodeMapsForRestart(bgdis_,staten_->wizard_);
  }

  // new cut for this time step
  Epetra_Vector idispcol( *boundarydis_->DofColMap() );
  idispcol.PutScalar( 0.0 );
  LINALG::Export(*aledispnp_,idispcol);
  state_ = Teuchos::rcp( new XFluidFluidState( *this, idispcol ) );

  // map of background-fluid's standard and enriched node-ids and
  // their dof-gids for new cut
  samemaps_ = xfluidfluid_timeint_->SaveAndCreateNewBgNodeMaps(bgdis_,state_->wizard_);

}//CutAndSaveBgFluidStatus()

// -------------------------------------------------------------------
// In this function:
// - Set new bg state vectors: veln_, velnm_, accn_
// - Set velnp_ and alevelnp_ to the veln_ and aleveln_ as start values
// - disp is the displacment of the embedded fluid at the time we want to
//   interpolate values from it
// -------------------------------------------------------------------
void FLD::XFluidFluid::SetBgStateVectors(Teuchos::RCP<Epetra_Vector>    disp)
{
// create patch boxes of embedded elements
//  std::map<int,GEO::CUT::BoundingBox> patchboxes;
//   if (alefluid_)
//     CreatePatchBoxes(patchboxes);

  const Teuchos::ParameterList& fdyn  = DRT::Problem::Instance()->FluidDynamicParams();
  INPAR::FLUID::InitialField initfield = DRT::INPUT::IntegralValue<INPAR::FLUID::InitialField>(fdyn,"INITIALFIELD");
  int startfuncno = fdyn.get<int>("STARTFUNCNO");
  if (monolithicfluidfluidfsi_ or
      (not monolithicfluidfluidfsi_ and step_>1 and alefluid_))
  {
    if (xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_FullProj or
        xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_KeepGhostValues or
        (xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_ProjIfMoved and (not samemaps_)) or
         xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_IncompProj)
    {
      // export the vectors to the column distribution map
      Teuchos::RCP<Epetra_Vector> alevelnpcol = LINALG::CreateVector(*embdis_->DofColMap(),true);
      Teuchos::RCP<Epetra_Vector> alevelncol = LINALG::CreateVector(*embdis_->DofColMap(),true);
      Teuchos::RCP<Epetra_Vector> alevelnmcol = LINALG::CreateVector(*embdis_->DofColMap(),true);
      Teuchos::RCP<Epetra_Vector> aleaccncol = LINALG::CreateVector(*embdis_->DofColMap(),true);
      Teuchos::RCP<Epetra_Vector> aleaccnpcol = LINALG::CreateVector(*embdis_->DofColMap(),true);
      Teuchos::RCP<Epetra_Vector> aledispcol = LINALG::CreateVector(*embdis_->DofColMap(),true);

      LINALG::Export(*alevelnp_,*alevelnpcol);
      LINALG::Export(*aleveln_,*alevelncol);
      LINALG::Export(*alevelnm_,*alevelnmcol);
      LINALG::Export(*aleaccn_,*aleaccncol);
      LINALG::Export(*aleaccnp_,*aleaccnpcol);
      LINALG::Export(*disp,*aledispcol);

      // we have five state vectors, which need values from the last time step
      xfluidfluid_timeint_->SetNewBgStatevectorAndProjectEmbToBg(bgdis_,
    		                                                     staten_->velnp_,state_->velnp_,alevelnpcol,
    		                                                     staten_->veln_,state_->veln_,alevelncol,
    		                                                     staten_->velnm_,state_->velnm_,alevelnmcol,
    		                                                     staten_->accn_,state_->accn_,aleaccncol,
    		                                                     staten_->accnp_,state_->accnp_,aleaccnpcol,
    		                                                     aledispcol);

      //-------------------------
      // Enforce incompressibility
      if (xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_IncompProj)
      {
        // find the incompressibility patch
        xfluidfluid_timeint_->PatchelementForIncompressibility(bgdis_,staten_->wizard_,state_->wizard_,state_->dbcmaps_);

        // prepare the incompressibility discretization and evaluate
        // incompressibility condition
        xfluidfluid_timeint_->PrepareIncompDiscret(bgdis_,state_->wizard_);
        xfluidfluid_timeint_->EvaluateIncompressibility(bgdis_,state_->wizard_);

        //solve the optimization problem for the state vectors
        xfluidfluid_timeint_->SolveIncompOptProb(state_->velnp_);
        xfluidfluid_timeint_->SolveIncompOptProb(state_->veln_);
        xfluidfluid_timeint_->SolveIncompOptProb(state_->velnm_);

        if (gmsh_debug_out_)
          state_->GmshOutput( *bgdis_, *embdis_, *boundarydis_, "incomvel", 0,  step_,  state_->velnp_, alevelnp_, aledispnp_);
      }

    }
    // Note: if Xff_TimeInt_ProjIfMoved is chosen and the maps remain the same
    // (TODO: they remain the same just for one dofset)
    // the enriched values are not projected from the embedded fluid anymore.
    else if(xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_ProjIfMoved and samemaps_)
    {
      // we use the old velocity as start value
      IO::cout << "samemaps ..." << IO::endl;
      state_->velnp_->Update(1.0,*staten_->velnp_,0.0);
      state_->veln_->Update(1.0,*staten_->veln_,0.0);
      state_->velnm_->Update(1.0,*staten_->velnm_,0.0);
      state_->accn_->Update(1.0,*staten_->accn_,0.0);
      state_->accnp_->Update(1.0,*staten_->accnp_,0.0);
    }
  }
  else if(step_==1 and initfield != INPAR::FLUID::initfield_zero_field){
    SetInitialFlowField(initfield,startfuncno);
  }

  if (monotype_ == FixedALEInterpolation)
  {
    // set the embedded state vectors to the new ale displacement
    // before updating the vectors of the old time step

    IO::cout << "Interpolate the embedded state vectors ... " << IO::endl;

    if (embdis_->Comm().NumProc() == 1)
    {
      xfluidfluid_timeint_->SetNewEmbStatevector(bgdis_,staten_->velnp_, alevelnp_, alevelnp_, aledispnp_, disp);
      xfluidfluid_timeint_->SetNewEmbStatevector(bgdis_,staten_->veln_ , aleveln_ , aleveln_ , aledispnp_, disp);
      xfluidfluid_timeint_->SetNewEmbStatevector(bgdis_, staten_->accnp_, aleaccnp_, aleaccnp_, aledispnp_, disp);
    }
    else
    {
      // export the vectors to the column distribution map
      Teuchos::RCP<Epetra_Vector> alevelnpcol = LINALG::CreateVector(*embdis_->DofColMap(),true);
      Teuchos::RCP<Epetra_Vector> alevelncol = LINALG::CreateVector(*embdis_->DofColMap(),true);
      Teuchos::RCP<Epetra_Vector> alevelnmcol = LINALG::CreateVector(*embdis_->DofColMap(),true);
      Teuchos::RCP<Epetra_Vector> aleaccncol = LINALG::CreateVector(*embdis_->DofColMap(),true);
      Teuchos::RCP<Epetra_Vector> aleaccnpcol = LINALG::CreateVector(*embdis_->DofColMap(),true);
      Teuchos::RCP<Epetra_Vector> aledispcol = LINALG::CreateVector(*embdis_->DofColMap(),true);
      Teuchos::RCP<Epetra_Vector> aledispnpcol = LINALG::CreateVector(*embdis_->DofColMap(),true);

      LINALG::Export(*alevelnp_,*alevelnpcol);
      LINALG::Export(*aleveln_,*alevelncol);
      LINALG::Export(*aleaccnp_,*aleaccnpcol);
      LINALG::Export(*disp,*aledispcol);
      LINALG::Export(*aledispnp_,*aledispnpcol);

      xfluidfluid_timeint_->SetNewEmbStatevector(bgdis_,staten_->velnp_, alevelnpcol, alevelnpcol, aledispnpcol, aledispcol);
      xfluidfluid_timeint_->SetNewEmbStatevector(bgdis_,staten_->veln_ , alevelncol , alevelncol , aledispnpcol, aledispcol);
      xfluidfluid_timeint_->SetNewEmbStatevector(bgdis_, staten_->accnp_, aleaccnpcol, aleaccnpcol, aledispnpcol, aledispcol);
    }
  }
}//SetBgStateVectors()

// -------------------------------------------------------------------
// Update the fluid solution for the monolithic timestep
// -------------------------------------------------------------------
void FLD::XFluidFluid::UpdateMonolithicFluidSolution()
{

  DRT::UTILS::ConditionSelector conds(*embdis_, "FSICoupling");
  std::vector<int> conddofs;
  for (int lnid=0; lnid<embdis_->NumMyRowNodes(); lnid++)
  {
    DRT::Node* embnode = embdis_->lRowNode(lnid);
    if(conds.ContainsNode(embnode->Id()))
    {
      std::vector<int> pgdofs(4);
      embdis_->Dof(embnode,0,pgdofs);

      conddofs.push_back(pgdofs[0]);
      conddofs.push_back(pgdofs[1]);
      conddofs.push_back(pgdofs[2]);
      //conddofs.push_back(pgdofs[3]); //dirichlet for pressure?
    }
  }

  fixedfsidofmap_ = Teuchos::rcp(new Epetra_Map(-1, conddofs.size(), &conddofs[0], 0, embdis_->Comm()));

  Teuchos::RCP<Epetra_Vector> tmpvec = Teuchos::rcp(new Epetra_Vector(*fixedfsidofmap_));
  tmpvec->PutScalar(1.0);

  // the toggle vector with values 1 and 0
  toggle_ = LINALG::CreateVector(*state_->fluidfluiddofrowmap_,true);
  LINALG::Export(*tmpvec,*toggle_);

  if (gmsh_debug_out_)
  {
    int count = 1;
    Teuchos::RCP<Epetra_Vector> testbg = state_->fluidfluidsplitter_->ExtractXFluidVector(toggle_);
    Teuchos::RCP<Epetra_Vector> testemb = state_->fluidfluidsplitter_->ExtractFluidVector(toggle_);
    state_->GmshOutput( *bgdis_, *embdis_, *boundarydis_, "toggle", count,  step_, testbg , testemb, aledispnp_);
    state_->GmshOutput(*bgdis_,*embdis_,*boundarydis_, "res_nachcut", gmsh_count_, step_, state_->velnp_, alevelnp_,
                       aledispnp_);
  }

  Solve();

}//FLD::XFluidFluid::UpdateMonolithicFluidSolution()

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::SetDirichletNeumannBC()
{
  // -------------------------------------------------------------------
  //         evaluate Dirichlet and Neumann boundary conditions
  // -------------------------------------------------------------------
  {
    Teuchos::ParameterList eleparams;

    // other parameters needed by the elements
    eleparams.set("total time",time_);

    // set vector values needed by elements
    bgdis_->ClearState();
    bgdis_->SetState("velaf",state_->velnp_);
    // predicted dirichlet values
    // velnp then also holds prescribed new dirichlet values
    bgdis_->EvaluateDirichlet(eleparams,state_->velnp_,Teuchos::null,Teuchos::null,Teuchos::null,state_->dbcmaps_);
    bgdis_->ClearState();

    embdis_->ClearState();
    embdis_->SetState("velaf",alevelnp_);
    //don't call this with the mapextractor. Otherwise the Mapextractor will
    //be built again.
    embdis_->EvaluateDirichlet(eleparams,alevelnp_,Teuchos::null,Teuchos::null,Teuchos::null);
    embdis_->ClearState();

    // set thermodynamic pressure
    eleparams.set("thermodynamic pressure",thermpressaf_);

    // Neumann
    state_->neumann_loads_->PutScalar(0.0);
    bgdis_->SetState("scaaf",state_->scaaf_);
//    bgdis_->EvaluateNeumann(eleparams,*state_->neumann_loads_);
    XFEM::EvaluateNeumann(state_->wizard_, eleparams, *bgdis_, *boundarydis_, state_->neumann_loads_);
    bgdis_->ClearState();

    aleneumann_loads_->PutScalar(0.0);
    embdis_->SetState("scaaf",alescaaf_);
    embdis_->EvaluateNeumann(eleparams,*aleneumann_loads_);
    embdis_->ClearState();
  }
}//FLD::XFluidFluid::SetDirichletNeumannBC()

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::SetHistoryValues()
{

  // -------------------------------------------------------------------
  // set part(s) of the rhs vector(s) belonging to the old timestep
  // (only meaningful for momentum part)
  //
  // stationary/af-generalized-alpha: hist_ = 0.0
  //
  // one-step-Theta:                  hist_ = veln_  + dt*(1-Theta)*accn_
  //
  // BDF2: for constant time step:    hist_ = 4/3 veln_  - 1/3 velnm_
  //
  //
  // ------------------------------------------------------------------
  COMBUST::UTILS::SetOldPartOfRighthandside(state_->veln_,state_->velnm_, state_->accn_,
                                                timealgo_, dta_, theta_, state_->hist_);
  COMBUST::UTILS::SetOldPartOfRighthandside(aleveln_,alevelnm_, aleaccn_,
                                                timealgo_, dta_, theta_, alehist_);

}//FLD::XFluidFluid::SetHistoryValues()

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::StatisticsAndOutput()
{
  // time measurement: output and statistics
  TEUCHOS_FUNC_TIME_MONITOR("      + output and statistics");

  // -------------------------------------------------------------------
  //          calculate lift'n'drag forces from the residual
  // -------------------------------------------------------------------
  LiftDrag();

  // compute equation-of-state factor
  //const double eosfac = thermpressaf_/gasconstant_;
  // store subfilter stresses for additional output
//   if (scale_similarity_)
//   {
//     Teuchos::RCP<Epetra_Vector> stress12 = CalcSFS(1,2);
//     statisticsmanager_->StoreNodalValues(step_, stress12);
//   }
  // -------------------------------------------------------------------
  //   add calculated velocity to mean value calculation (statistics)
  // -------------------------------------------------------------------
//  statisticsmanager_->DoTimeSample(step_,time_,eosfac);

  // -------------------------------------------------------------------
  //                        compute flow rates
  // -------------------------------------------------------------------
//  ComputeFlowRates();

  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  Output();

  // -------------------------------------------------------------------
  //          dumping of turbulence statistics if required
  // -------------------------------------------------------------------
  //statisticsmanager_->DoOutput(output_,step_,eosfac);

  return;
}//FLD::XFluidFluid::StatisticsAndOutput()


/*----------------------------------------------------------------------*
 |  write discretization output                            schott 03/12 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluid::OutputDiscret()
{
  // compute the current solid and boundary position
  std::map<int,LINALG::Matrix<3,1> >      curralepositions;
  std::map<int,LINALG::Matrix<3,1> >      currinterfacepositions;

  //---------------------------------- GMSH DISCRET OUTPUT (element and node ids for all discretizations) ------------------------
  if(gmsh_discret_out_)
  {
    ExtractNodeVectors(*embdis_, aledispnp_, curralepositions);
    ExtractNodeVectors(*boundarydis_, idispnp_, currinterfacepositions);

    // cast to DiscretizationXFEM
    Teuchos::RCP<DRT::DiscretizationFaces> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret_, true);
    if (xdiscret == Teuchos::null)
      dserror("Failed to cast DRT::Discretization to DRT::DiscretizationFaces.");

    // output for Element and Node IDs
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("element_node_id", step_, gmsh_step_diff_, gmsh_debug_out_screen_, discret_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    gmshfilecontent.setf(std::ios::scientific,std::ios::floatfield);
    gmshfilecontent.precision(16);

    disToStream(discret_,     "bg ",     true, false, true, false, true,  false, gmshfilecontent);
    disToStream(embdis_,      "emb",     true, false, true, false, true,  false, gmshfilecontent, &curralepositions);
    disToStream(boundarydis_, "cut",     true, true,  true, true,  false, false, gmshfilecontent, &currinterfacepositions);

    if (nitsche_evp_)
    	disToStream(embboundarydis_, "emb_bound", true, true, true, true, false,  false, gmshfilecontent);

    gmshfilecontent.close();

  } // end if gmsh_discret_out_
}


// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::Output()
{
  const bool write_visualization_data = step_%upres_ == 0;
  const bool write_restart_data = step_!=0 and uprestart_ != 0 and step_%uprestart_ == 0;

  const int step_diff = 500;
  bool screen_out = false;

  if( !write_visualization_data && !write_restart_data )
    return;

  OutputDiscret();

  std::map<int,LINALG::Matrix<3,1> >      curralepositions;
  std::map<int,LINALG::Matrix<3,1> >      currinterfacepositions;

  ExtractNodeVectors(*embdis_, aledispnp_, curralepositions);
  ExtractNodeVectors(*boundarydis_, idispnp_, currinterfacepositions);

  //  ART_exp_timeInt_->Output();
  // output of solution


  int count = -1;
  if (gmsh_sol_out_)
    state_->GmshOutput( *bgdis_, *embdis_, *boundarydis_, "result", count,  step_, state_->velnp_ , alevelnp_, aledispnp_);




  //---------------------------------- GMSH DISCRET OUTPUT (element and node ids for all discretizations) ------------------------
  if(gmsh_EOS_out_)
  {

    // cast to DiscretizationXFEM
    Teuchos::RCP<DRT::DiscretizationFaces> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret_, true);
    if (xdiscret == Teuchos::null)
      dserror("Failed to cast DRT::Discretization to DRT::DiscretizationFaces.");

    Teuchos::RCP<DRT::DiscretizationFaces> embxdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(embdis_, true);
    if (embxdiscret == Teuchos::null)
      dserror("Failed to cast DRT::Discretization to DRT::DiscretizationFaces.");


    // output for Element and Node IDs
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("EOS", step_, step_diff, screen_out, discret_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());
    gmshfilecontent.setf(std::ios::scientific,std::ios::floatfield);
    gmshfilecontent.precision(16);

    if( xdiscret->FilledExtension() == true && ghost_penalty_ ) // stabilization output
    {
      // draw internal faces elements with associated face's gid
      gmshfilecontent << "View \" " << "ghost penalty stabilized \" {\n";

      for (int i=0; i<xdiscret->NumMyRowFaces(); ++i)
      {
        const DRT::Element* actele = xdiscret->lRowFace(i);
        std::map<int,int> & ghost_penalty_map = state_->EdgeStab()->GetGhostPenaltyMap();

        std::map<int,int>::iterator it = ghost_penalty_map.find(actele->Id());
        if(it != ghost_penalty_map.end())
        {
          int ghost_penalty = it->second;

          if(ghost_penalty) IO::GMSH::elementAtInitialPositionToStream(double(ghost_penalty),actele, gmshfilecontent);
        }
        else dserror("face %d in map not found", actele->Id());
      }
      gmshfilecontent << "};\n";
    }
    if(xdiscret->FilledExtension() == true && edge_based_)
    {
      // draw internal faces elements with associated face's gid
      gmshfilecontent << "View \" " << "edgebased stabilized \" {\n";

      for (int i=0; i<xdiscret->NumMyRowFaces(); ++i)
      {
        const DRT::Element* actele = xdiscret->lRowFace(i);
        std::map<int,int> & edge_based_map = state_->EdgeStab()->GetEdgeBasedMap();
        std::map<int,int>::iterator it = edge_based_map.find(actele->Id());

        if(it != edge_based_map.end())
        {
          int edge_stab =it->second;

          if(edge_stab) IO::GMSH::elementAtInitialPositionToStream(double(edge_stab),actele, gmshfilecontent);
        }
      }
      gmshfilecontent << "};\n";
    } // end stabilization output

    gmshfilecontent.close();

  } // end if gmsh_discret_out_



  if (write_visualization_data)
  {
    //Velnp()->Print(cout);
    /*vector<int> lm;
    lm.push_back(17041);
    lm.push_back(17042);
    lm.push_back(17043);
    lm.push_back(17044);
    Epetra_SerialDenseVector result(lm.size());
    DRT::UTILS::ExtractMyValues(*Velnp(),result,lm);
    result.Print(cout<<std::scientific<<setprecision(10));*/

    // step number and time
    output_->NewStep(step_,time_);

//     for (int iter=0; iter<state_->velnp_->MyLength();++iter)
//     {
//       int gid = state_->velnp_->Map().GID(iter);
//       if (state_->velnpoutput_->Map().MyGID(gid))
//         (*state_->velnpoutput_)[state_->velnpoutput_->Map().LID(gid)]=(*state_->velnpoutput_)[state_->velnpoutput_->Map().LID(gid)] +                                                                (*state_->velnp_)[state_->velnp_->Map().LID(gid)];
//     }

    const Epetra_Map* dofrowmap = dofset_out_->DofRowMap(); // original fluid unknowns
    const Epetra_Map* xdofrowmap = bgdis_->DofRowMap();    // fluid unknown for current cut

    for (int i=0; i<bgdis_->NumMyRowNodes(); ++i)
    {
      // get row node via local id
      const DRT::Node* xfemnode = bgdis_->lRowNode(i);

      // the dofset_out_ contains the original dofs for each row node
      const std::vector<int> gdofs_original(dofset_out_->Dof(xfemnode));

      //cout << "node->Id() " << gid << "gdofs " << gdofs_original[0] << " " << gdofs_original[1] << "etc" << endl;

      // if the dofs for this node do not exist in the xdofrowmap, then a hole is given
      // else copy the right nodes
      const std::vector<int> gdofs_current(bgdis_->Dof(xfemnode));

      if(gdofs_current.size() == 0)
      {
        // IO::cout << "no dofs available->hole" << IO::endl;
      }
      else if(gdofs_current.size() == gdofs_original.size())
      {
        //IO::cout << "same number of dofs available" << IO::endl;
      }
      else if(gdofs_current.size() > gdofs_original.size())
      {
        //IO::cout << "more dofs available->decide" << IO::endl;
      }
      else IO::cout << "decide which dofs can be copied and which have to be set to zero" << IO::endl;


      if(gdofs_current.size() == 0) //void
      {
        size_t numdof = gdofs_original.size();

#ifdef INTERPOLATEFOROUTPUT
        LINALG::Matrix<3,1> bgnodecords(true);
        bgnodecords(0,0) = xfemnode->X()[0];
        bgnodecords(1,0) = xfemnode->X()[1];
        bgnodecords(2,0) = xfemnode->X()[2];

        // take the values of embedded fluid is available
        bool insideelement = false;
        int count = 0;
        // check all embedded elements to find the right one
        for (int e=0; e<embdis_->NumMyColElements(); e++)
        {
          DRT::Element* pele = embdis_->lColElement(e);
          LINALG::Matrix<4,1> interpolatedvec(true);
          insideelement = ComputeSpacialToElementCoordAndProject(pele,bgnodecords,interpolatedvec,*alevelnp_,aledispnp_,embdis_);
          if (insideelement)
          {
            // hier set state
            (*outvec_fluid_)[dofrowmap->LID(gdofs_original[0])] = interpolatedvec(0);
            (*outvec_fluid_)[dofrowmap->LID(gdofs_original[1])] = interpolatedvec(1);
            (*outvec_fluid_)[dofrowmap->LID(gdofs_original[2])] = interpolatedvec(2);
            (*outvec_fluid_)[dofrowmap->LID(gdofs_original[3])] = interpolatedvec(3);
            break;
          }
          count ++;
          if (count == embdis_->NumMyColElements())
          {
            for (std::size_t idof = 0; idof < numdof; ++idof)
            {
              //cout << dofrowmap->LID(gdofs[idof]) << endl;
              (*outvec_fluid_)[dofrowmap->LID(gdofs_original[idof])] = 0.0;
            }
          }
        }
#else
        for (std::size_t idof = 0; idof < numdof; ++idof)
        {
          //cout << dofrowmap->LID(gdofs[idof]) << endl;
          (*outvec_fluid_)[dofrowmap->LID(gdofs_original[idof])] = 0.0;
        }
#endif
      }
      else if(gdofs_current.size() == gdofs_original.size())
      {
        size_t numdof = gdofs_original.size();
        // copy all values
        for (std::size_t idof = 0; idof < numdof; ++idof)
        {
          //cout << dofrowmap->LID(gdofs[idof]) << endl;
          (*outvec_fluid_)[dofrowmap->LID(gdofs_original[idof])] = (*state_->velnp_)[xdofrowmap->LID(gdofs_current[idof])];
        }
      }
      // choose the std-dofset for the nodes which has std-dofs. If a node
      // doesn't have a std-dofset the first ghost-set is taken.
      else if(gdofs_current.size() % gdofs_original.size() == 0) //multiple dofsets
      {
        // if there are multiple dofsets we write output for the standard dofset
        GEO::CUT::Node* node = state_->wizard_->GetNode(xfemnode->Id());

        GEO::CUT::Point* p = node->point();

        const std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp> > & dofcellsets = node->DofCellSets();

        int nds = 0;
        bool is_std_set = false;

        // find the standard dofset
        for(std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp> >::const_iterator cellsets= dofcellsets.begin();
            cellsets!=dofcellsets.end();
            cellsets++)
        {
          // at least one vc has to contain the node
          for(std::set<GEO::CUT::plain_volumecell_set>::const_iterator sets=cellsets->begin(); sets!=cellsets->end(); sets++)
          {
            const GEO::CUT::plain_volumecell_set& set = *sets;

            for(GEO::CUT::plain_volumecell_set::const_iterator vcs=set.begin(); vcs!=set.end(); vcs++)
            {
              if((*vcs)->Contains(p))
              {
                // return if at least one vc contains this point
                is_std_set=true;
                break;
              }
            }
            // break the outer loop if at least one vc contains this point
            if(is_std_set == true) break;
          }
          if(is_std_set == true) break;
          nds++;
        }

        size_t numdof = gdofs_original.size();
        size_t offset = 0;  // no offset in case of no std-dofset means take the first dofset for output

        if(is_std_set)
        {
          offset = numdof*nds;   // offset to start the loop over dofs at the right nds position
        }

        // copy all values
        for (std::size_t idof = 0; idof < numdof; ++idof)
        {
          (*outvec_fluid_)[dofrowmap->LID(gdofs_original[idof])] = (*state_->velnp_)[xdofrowmap->LID(gdofs_current[offset+idof])];
        }

        // if there are just ghost values, then we write output for the set with largest physical fluid volume
      }
      else dserror("unknown number of dofs for output");

    };


    // velocity/pressure vector
    //output_->WriteVector("velnp",state_->velnpoutput_);
    output_->WriteVector("velnp", outvec_fluid_);

    // output (hydrodynamic) pressure for visualization
    Teuchos::RCP<Epetra_Vector> pressure = velpressplitterForOutput_.ExtractCondVector(outvec_fluid_);
    output_->WriteVector("pressure", pressure);

    // write domain decomposition for visualization (only once!)
    if (step_==upres_) output_->WriteElementData(true);

  }


  // write restart
  if (write_restart_data)
  {
    IO::cout << "---  write restart... " << IO::endl;

    restart_count_++;

    // velocity/pressure vector
    output_->WriteVector("velnp_bg",state_->velnp_);

    // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
    output_->WriteVector("accnp_bg",state_->accnp_);
    output_->WriteVector("accn_bg",state_->accn_);
    output_->WriteVector("veln_bg",state_->veln_);
    output_->WriteVector("velnm_bg",state_->velnm_);
  }


  //-----------------------------------------------------------
  // REMARK on "Why to clear the MapCache" for restarts
  //-----------------------------------------------------------
  // every time, when output-vectors are written based on a new(!), still unknown map
  // the map is stored in a mapstack (io.cpp, WriteVector-routine) for efficiency in standard applications.
  // However, in the XFEM for each timestep or FSI-iteration we have to cut and the map is created newly
  // (done by the FillComplete call).
  // This is the reason why the MapStack increases and the storage is overwritten for large problems.
  // Hence, we have clear the MapCache in regular intervals of written restarts.
  // In case of writing paraview-output, here, we use a standard map which does not change over time, that's okay.
  // For the moment, restart_count = 5 is set quite arbitrary, in case that we need for storage, we have to reduce this number

  if(restart_count_ == 5)
  {
    if(myrank_ == 0) IO::cout << "\t... Clear MapCache after " << restart_count_ << " written restarts." << IO::endl;

    output_->ClearMapCache(); // clear the output's map-cache
    restart_count_ = 0;
  }


  // embedded fluid output
  if (write_visualization_data)
  {
    // step number and time
    emboutput_->NewStep(step_,time_);

    // velocity/pressure vector
    emboutput_->WriteVector("velnp",alevelnp_);

    // (hydrodynamic) pressure
    Teuchos::RCP<Epetra_Vector> pressure = alevelpressplitter_.ExtractCondVector(alevelnp_);
    emboutput_->WriteVector("pressure", pressure);

    //output_.WriteVector("residual", trueresidual_);

    if (alefluid_) emboutput_->WriteVector("dispnp", aledispnp_);

    if (step_==upres_) emboutput_->WriteElementData(true);

  }

  // write restart
  if (write_restart_data)
  {
    // velocity/pressure vector
    emboutput_->WriteVector("velnp_emb",alevelnp_);

    //output_->WriteVector("residual", trueresidual_);
    if (alefluid_)
    {
      emboutput_->WriteVector("dispnp_emb", aledispnp_);
      emboutput_->WriteVector("dispn_emb", aledispn_);
      emboutput_->WriteVector("dispnm_emb",aledispnm_);
    }

    if(alefluid_)
      emboutput_->WriteVector("gridv_emb", gridv_);

    // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
    emboutput_->WriteVector("accnp_emb",aleaccnp_);
    emboutput_->WriteVector("accn_emb", aleaccn_);
    emboutput_->WriteVector("veln_emb", aleveln_);
    emboutput_->WriteVector("velnm_emb", alevelnm_);
    // emboutput_->WriteVector("neumann_loads",aleneumann_loads_);

  }

   return;
}// XFluidFluid::Output


/*----------------------------------------------------------------------*
 | print discretization to gmsh stream                     schott 01/13 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluid::disToStream(Teuchos::RCP<DRT::Discretization> dis,
                                   const std::string& disname,
                                   const bool elements,
                                   const bool elecol,
                                   const bool nodes,
                                   const bool nodecol,
                                   const bool faces,
                                   const bool facecol,
                                   std::ostream& s,
                                   std::map<int, LINALG::Matrix<3,1> >* curr_pos)
{
  if(elements)
  {
    // draw bg elements with associated gid
    s << "View \" " << disname;
    if(elecol)
    {
      s << " col e->Id() \" {\n";
      for (int i=0; i<dis->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = dis->lColElement(i);
        if(curr_pos == NULL)
          IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, s);
        else
          IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, *curr_pos, s);
      };
    }
    else
    {
      s << " row e->Id() \" {\n";
      for (int i=0; i<dis->NumMyRowElements(); ++i)
      {
        const DRT::Element* actele = dis->lRowElement(i);
        if(curr_pos == NULL)
          IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, s);
        else
          IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, *curr_pos, s);
      };
    }
    s << "};\n";
  }

  if(nodes)
  {
    s << "View \" " << disname;
    if(nodecol)
    {
      s << " col n->Id() \" {\n";
      for (int i=0; i<dis->NumMyColNodes(); ++i)
      {
        const DRT::Node* actnode = dis->lColNode(i);
        LINALG::Matrix<3,1> pos(true);

        if(curr_pos != NULL)
        {
          const LINALG::Matrix<3,1>& curr_x = curr_pos->find(actnode->Id())->second;
          pos(0) = curr_x(0);
          pos(1) = curr_x(1);
          pos(2) = curr_x(2);
        }
        else
        {
          const LINALG::Matrix<3,1> x(actnode->X());
          pos(0) = x(0);
          pos(1) = x(1);
          pos(2) = x(2);
        }
        IO::GMSH::cellWithScalarToStream(DRT::Element::point1, actnode->Id(), pos, s);
      }
    }
    else
    {
      s << " row n->Id() \" {\n";
      for (int i=0; i<dis->NumMyRowNodes(); ++i)
      {
        const DRT::Node* actnode = dis->lRowNode(i);
        LINALG::Matrix<3,1> pos(true);

        if(curr_pos != NULL)
        {
          const LINALG::Matrix<3,1>& curr_x = curr_pos->find(actnode->Id())->second;
          pos(0) = curr_x(0);
          pos(1) = curr_x(1);
          pos(2) = curr_x(2);
        }
        else
        {
          const LINALG::Matrix<3,1> x(actnode->X());
          pos(0) = x(0);
          pos(1) = x(1);
          pos(2) = x(2);
        }
        IO::GMSH::cellWithScalarToStream(DRT::Element::point1, actnode->Id(), pos, s);
      }
    }
    s << "};\n";
  }

  if(faces)
  {
    // cast to DiscretizationXFEM
    Teuchos::RCP<DRT::DiscretizationFaces> xdis = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(dis, true);
    if (xdis == Teuchos::null)
      dserror("Failed to cast DRT::Discretization to DRT::DiscretizationFaces.");

    if( xdis->FilledExtension() == true )     // faces output
    {
      s << "View \" " << disname;
      if(facecol)
      {
        s << " col f->Id() \" {\n";
        for (int i=0; i<xdis->NumMyColFaces(); ++i)
        {
          const DRT::Element* actele = xdis->lColFace(i);
          if(curr_pos == NULL)
            IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, s);
          else
            IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, *curr_pos, s);
        };
      }
      else
      {
        s << " row f->Id() \" {\n";
        for (int i=0; i<xdis->NumMyRowFaces(); ++i)
        {
          const DRT::Element* actele = xdis->lRowFace(i);
          if(curr_pos == NULL)
            IO::GMSH::elementAtInitialPositionToStream(double(actele->Id()), actele, s);
          else
            IO::GMSH::elementAtCurrentPositionToStream(double(actele->Id()), actele, *curr_pos, s);
        };
      }
      s << "};\n";
    }
  }
}


/*--------------------------------------------------------------------------*
 | extract the nodal vectors and store them in node-vector-map schott 01/13 |
 *--------------------------------------------------------------------------*/
void FLD::XFluidFluid::ExtractNodeVectors(DRT::Discretization & dis,
                                          Teuchos::RCP<Epetra_Vector> dofrowvec,
                                          std::map<int, LINALG::Matrix<3,1> >& nodevecmap)
{
  Epetra_Vector dispcol( *dis.DofColMap() );
  dispcol.PutScalar( 0. );

  LINALG::Export(*dofrowvec,dispcol);

  nodevecmap.clear();

  for (int lid = 0; lid < dis.NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = dis.lColNode(lid);
    std::vector<int> lm;
    dis.Dof(node, lm);
    std::vector<double> mydisp;
    DRT::UTILS::ExtractMyValues(dispcol,mydisp,lm);
    if (mydisp.size() < 3)
      dserror("we need at least 3 dofs here");

    LINALG::Matrix<3,1> currpos;
    currpos(0) = node->X()[0] + mydisp[0];
    currpos(1) = node->X()[1] + mydisp[1];
    currpos(2) = node->X()[2] + mydisp[2];
    nodevecmap.insert(std::make_pair(node->Id(),currpos));
  }
}


// -------------------------------------------------------------------
// set general face fluid parameter (BS 06/2014)
// -------------------------------------------------------------------
void FLD::XFluidFluid::SetElementGeneralFluidXFEMParameter()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action",FLD::set_general_fluid_xfem_parameter); // do not call another action as then another object of the std-class will be created

  //------------------------------------------------------------------------------------------------------
  // set general element parameters
  eleparams.set("form of convective term",convform_);
  eleparams.set<int>("Linearisation",newton_);
  eleparams.set<int>("Physical Type", physicaltype_);

  // parameter for stabilization
  eleparams.sublist("RESIDUAL-BASED STABILIZATION") = params_->sublist("RESIDUAL-BASED STABILIZATION");

  //set time integration scheme
  eleparams.set<int>("TimeIntegrationScheme", timealgo_);

  //------------------------------------------------------------------------------------------------------
  // set general parameters for turbulent flow
  eleparams.sublist("TURBULENCE MODEL") = params_->sublist("TURBULENCE MODEL");

  // set model-dependent parameters
  eleparams.sublist("SUBGRID VISCOSITY") = params_->sublist("SUBGRID VISCOSITY");
  eleparams.sublist("MULTIFRACTAL SUBGRID SCALES") = params_->sublist("MULTIFRACTAL SUBGRID SCALES");


  //------------------------------------------------------------------------------------------------------
  // set general XFEM element parameters

  eleparams.sublist("XFEM")                         = params_->sublist("XFEM");
  eleparams.sublist("XFLUID DYNAMIC/GENERAL")       = params_->sublist("XFLUID DYNAMIC/GENERAL");
  eleparams.sublist("XFLUID DYNAMIC/STABILIZATION") = params_->sublist("XFLUID DYNAMIC/STABILIZATION");


  //------------------------------------------------------------------------------------------------------
  // set the params in the XFEM-parameter-list class
  DRT::ELEMENTS::FluidType::Instance().PreEvaluate(*discret_,eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  return;
}

// -------------------------------------------------------------------
// set general face fluid parameter (BS 06/2014)
// -------------------------------------------------------------------
void FLD::XFluidFluid::SetFaceGeneralFluidXFEMParameter()
{

  //------------------------------------------------------------------------------------------------------
  // set general fluid stabilization parameter for faces
  {
    Teuchos::ParameterList faceparams;

    faceparams.set<int>("action",FLD::set_general_face_fluid_parameter);

    faceparams.sublist("EDGE-BASED STABILIZATION")     = params_->sublist("EDGE-BASED STABILIZATION");

    faceparams.set<int>("STABTYPE", DRT::INPUT::IntegralValue<INPAR::FLUID::StabType>( params_->sublist("RESIDUAL-BASED STABILIZATION"), "STABTYPE"));

    DRT::ELEMENTS::FluidIntFaceType::Instance().PreEvaluate(*discret_,faceparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  }

  //------------------------------------------------------------------------------------------------------
  // set XFEM specific parameter for faces
  {
    Teuchos::ParameterList faceparams;

    faceparams.set<int>("action",FLD::set_general_face_xfem_parameter);

    // set general fluid face parameters are contained in the following two sublists
    faceparams.sublist("XFLUID DYNAMIC/STABILIZATION") = params_->sublist("XFLUID DYNAMIC/STABILIZATION");

    DRT::ELEMENTS::FluidIntFaceType::Instance().PreEvaluate(*discret_,faceparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  }

  return;
}

// -------------------------------------------------------------------
// set general time parameter (AE 01/2011)
// -------------------------------------------------------------------
void FLD::XFluidFluid::SetElementTimeParameter()
{
  Teuchos::ParameterList eleparams;

  // set action
  eleparams.set<int>("action",FLD::set_time_parameter);
  // set time integration scheme
  eleparams.set<int>("TimeIntegrationScheme", timealgo_);
  // set general element parameters
  eleparams.set("dt",dta_);
  eleparams.set("theta",theta_);
  eleparams.set("omtheta",omtheta_);

  // set scheme-specific element parameters and vector values
  if (timealgo_==INPAR::FLUID::timeint_stationary)
  {
    eleparams.set("total time",time_);
  }
  else if (timealgo_==INPAR::FLUID::timeint_afgenalpha)
  {
    eleparams.set("total time",time_-(1-alphaF_)*dta_);
    eleparams.set("alphaF",alphaF_);
    eleparams.set("alphaM",alphaM_);
    eleparams.set("gamma",gamma_);
  }
  else
  {
    eleparams.set("total time",time_);
  }

  DRT::ELEMENTS::FluidType::Instance().PreEvaluate(*bgdis_,eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  // not necessary
  //DRT::ELEMENTS::FluidType::Instance().PreEvaluate(*embdis_,eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

}

// -------------------------------------------------------------------
// return time integration factor
// -------------------------------------------------------------------
const double FLD::XFluidFluid::TimIntParam() const
{
  double retval = 0.0;
  switch (TimIntScheme())
  {
  case INPAR::FLUID::timeint_afgenalpha:
  case INPAR::FLUID::timeint_npgenalpha:
    // this is the interpolation weight for quantities from last time step
    retval = 1.0 - alphaF_;
  break;
  case INPAR::FLUID::timeint_one_step_theta:
    // this is the interpolation weight for quantities from last time step
    retval = 0.0;
  break;
  case INPAR::FLUID::timeint_bdf2:
    // this is the interpolation weight for quantities from last time step
    retval = 0.0;
  break;
  case INPAR::FLUID::timeint_stationary:
    // this is the interpolation weight for quantities from last time step
    retval = 0.0;
  break;
  default:
    dserror("Unknown time integration scheme");
  break;
  }
  return retval;
}


//----------------------------------------------------------------------
// LiftDrag                                               rasthofer
//----------------------------------------------------------------------
//calculate lift&drag forces and angular moments
//
//Lift and drag forces are based upon the right hand side true-residual entities
//of the corresponding nodes. The contribution of the end node of a line is entirely
//added to a present L&D force.
//--------------------------------------------------------------------
void FLD::XFluidFluid::LiftDrag() const
{
  // initially check whether computation of lift and drag values is required
  if (params_->get<bool>("LIFTDRAG")) 
  {
    // in this map, the results of the lift drag calculation are stored
    Teuchos::RCP<std::map<int,std::vector<double> > > liftdragvals;

    FLD::UTILS::LiftDrag(embdis_,aletrueresidual_,aledispnp_,numdim_,liftdragvals,alefluid_);

    if (liftdragvals!=Teuchos::null and embdis_->Comm().MyPID() == 0)
      FLD::UTILS::WriteLiftDragToFile(time_, step_, *liftdragvals);
  }

  return;
}


//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------
void FLD::XFluidFluid::GenAlphaIntermediateValues()
{
  state_->GenAlphaIntermediateValues();
}
//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------
void FLD::XFluidFluid::AssembleMatAndRHS()
{

}
//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------
void FLD::XFluidFluid::GenAlphaUpdateAcceleration()
{
  state_->GenAlphaUpdateAcceleration();
}

//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------
void FLD::XFluidFluid::XFluidFluidState::GenAlphaIntermediateValues()
{
  //       n+alphaM                n+1                      n
  //    acc         = alpha_M * acc     + (1-alpha_M) *  acc
  //       (i)                     (i)
  {
    // extract the degrees of freedom associated with velocities
    // only these are allowed to be updated, otherwise you will
    // run into trouble in loma, where the 'pressure' component
    // is used to store the acceleration of the temperature
    Teuchos::RCP<Epetra_Vector> onlyaccn  = velpressplitter_.ExtractOtherVector(accn_ );
    Teuchos::RCP<Epetra_Vector> onlyaccnp = velpressplitter_.ExtractOtherVector(accnp_);

    Teuchos::RCP<Epetra_Vector> onlyaccam = Teuchos::rcp(new Epetra_Vector(onlyaccnp->Map()));

    onlyaccam->Update((xfluid_.alphaM_),*onlyaccnp,(1.0-xfluid_.alphaM_),*onlyaccn,0.0);

    // copy back into global vector
    LINALG::Export(*onlyaccam,*accam_);
  }

  // set intermediate values for velocity
  //
  //       n+alphaF              n+1                   n
  //      u         = alpha_F * u     + (1-alpha_F) * u
  //       (i)                   (i)
  //
  // and pressure
  //
  //       n+alphaF              n+1                   n
  //      p         = alpha_F * p     + (1-alpha_F) * p
  //       (i)                   (i)
  //
  // note that its af-genalpha with mid-point treatment of the pressure,
  // not implicit treatment as for the genalpha according to Whiting
  velaf_->Update((xfluid_.alphaF_),*velnp_,(1.0-xfluid_.alphaF_),*veln_,0.0);
}
//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------
void FLD::XFluidFluid::XFluidFluidState::GenAlphaUpdateAcceleration()
{
  //                                  n+1     n
  //                               vel   - vel
  //       n+1      n  gamma-1.0      (i)
  //    acc    = acc * --------- + ------------
  //       (i)           gamma      gamma * dt
  //

  // extract the degrees of freedom associated with velocities
  // only these are allowed to be updated, otherwise you will
  // run into trouble in loma, where the 'pressure' component
  // is used to store the acceleration of the temperature
  Teuchos::RCP<Epetra_Vector> onlyaccn  = velpressplitter_.ExtractOtherVector(accn_ );
  Teuchos::RCP<Epetra_Vector> onlyveln  = velpressplitter_.ExtractOtherVector(veln_ );
  Teuchos::RCP<Epetra_Vector> onlyvelnp = velpressplitter_.ExtractOtherVector(velnp_);

  Teuchos::RCP<Epetra_Vector> onlyaccnp = Teuchos::rcp(new Epetra_Vector(onlyaccn->Map()));

  const double fact1 = 1.0/(xfluid_.gamma_*xfluid_.dta_);
  const double fact2 = 1.0 - (1.0/xfluid_.gamma_);
  onlyaccnp->Update(fact2,*onlyaccn,0.0);
  onlyaccnp->Update(fact1,*onlyvelnp,-fact1,*onlyveln,1.0);

  // copy back into global vector
  LINALG::Export(*onlyaccnp,*accnp_);

}
//---------------------------------------------------------------
// Teuchos::RCP<Epetra_Vector> FLD::XFluidFluid::calcStresses()
// {
//   string condstring("FluidStressCalc");

//   Teuchos::RCP<Epetra_Vector> integratedshapefunc = IntegrateInterfaceShape(condstring);

//   // compute traction values at specified nodes; otherwise do not touch the zero values
//   for (int i=0;i<integratedshapefunc->MyLength();i++)
//   {
//     if ((*integratedshapefunc)[i] != 0.0)
//     {
//       // overwrite integratedshapefunc values with the calculated traction coefficients,
//       // which are reconstructed out of the nodal forces (trueresidual_) using the
//       // same shape functions on the boundary as for velocity and pressure.
//       (*integratedshapefunc)[i] = (*fluidtrueresidual_)[i]/(*integratedshapefunc)[i];
//     }
//   }

//   return integratedshapefunc;
//}

// -------------------------------------------------------------------
//
// -------------------------------------------------------------------
void FLD::XFluidFluid::EvaluateErrorComparedToAnalyticalSol()
{
  // this functions provides a general implementation for calculating error norms between computed solutions
  // and an analytical solution which is implemented or given by a function in the input file

  INPAR::FLUID::CalcError calcerr = DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_,"calculate error");

  if(calcerr != INPAR::FLUID::no_error_calculation)
  {
    //TODO: decide between absolute and relative errors

    // set the time to evaluate errors
    //

    // define the norms that have to be computed

    //-------------------------------------------------------------------------------------------------------------------
    // domain error norms w.r.t incompressible Navier-Stokes equations
    //
    //--------------------------------------
    // background domain
    //--------------------------------------
    // standard domain errors
    // 1.   || u - u_b ||_L2(Omega)            =   standard L2-norm for velocity
    // 2.   || grad( u - u_b ) ||_L2(Omega)    =   standard H1-seminorm for velocity
    // 3.   || u - u_b ||_H1(Omega)            =   standard H1-norm for velocity
    //                                         =   sqrt( || u - u_b ||^2_L2(Omega) + || grad( u - u_b ) ||^2_L2(Omega) )
    // 4.   || p - p_b ||_L2(Omega)            =   standard L2-norm for for pressure
    //
    // viscosity-scaled domain errors
    // 5.   || nu^(+1/2) grad( u - u_b ) ||_L2(Omega)   =   visc-scaled H1-seminorm for velocity
    //                                                  =   nu^(+1/2) * || grad( u - u_b ) ||_L2(Omega) (for homogeneous visc)
    // 6.   || nu^(-1/2) ( p - p_b ) ||_L2(Omega)       =   visc-scaled L2-norm for for pressure
    //                                                  =   nu^(-1/2) * || p - p_b ||_L2(Omega) (for homogeneous visc)
    //
    // Error on Functionals from solution (Sudhakar)
    // 7.   | sin(x) ( u,x - u,x exact ) | (background)
    //
    //
    //--------------------------------------
    // embedded domain
    //--------------------------------------
    // 1.   || u - u_e ||_L2(Omega)            =   standard L2-norm for velocity
    // 2.   || grad( u - u_e ) ||_L2(Omega)    =   standard H1-seminorm for velocity
    // 3.   || u - u_e ||_H1(Omega)            =   standard H1-norm for velocity
    //                                         =   sqrt( || u - u_e ||^2_L2(Omega) + || grad( u - u_e ) ||^2_L2(Omega) )
    // 4.  || p - p_e ||_L2(Omega)            =   standard L2-norm for for pressure
    //
    // viscosity-scaled domain errors
    // 5.  || nu^(+1/2) grad( u - u_e ) ||_L2(Omega)   =   visc-scaled H1-seminorm for velocity
    //                                                  =   nu^(+1/2) * || grad( u - u_e ) ||_L2(Omega) (for homogeneous visc)
    // 6.  || nu^(-1/2) ( p - p_e ) ||_L2(Omega)       =   visc-scaled L2-norm for for pressure
    //                                                  =   nu^(-1/2) * || p - p_e ||_L2(Omega) (for homogeneous visc)
    // Error on Functionals from solution (Sudhakar)
    // 7.  | sin(x) ( u,x - u,x exact ) | (embedded)
    //-------------------------------------------------------------------------------------------------------------------
    // interface/boundary error norms at the XFEM-interface, boundary
    // w.r.t Nitsche's method to enforce interface/boundary conditions
    //
    // 1.   || nu^(+1/2) (u_b - u_e) ||_H1/2(Gamma)          =  broken H1/2 Sobolev norm for boundary/coupling condition
    // 2.   || nu^(+1/2) grad( u_b - u_e )*n ||_H-1/2(Gamma) =  standard H-1/2 Sobolev norm for normal flux (velocity part)
    // 3.   || nu^(-1/2) ( p_b - p_e )*n ||_H-1/2(Gamma)     =  standard H-1/2 Sobolev norm for normal flux (pressure part)
    //
    //-------------------------------------------------------------------------------------------------------------------
    // errors introduced by stabilizations (edge-based fluid stabilizations and ghost-penalty stabilizations)
    //
    // ...
    //-------------------------------------------------------------------------------------------------------------------

    // number of norms that have to be calculated
    int num_dom_norms    = 7;
    int num_interf_norms = 5;
    int num_stab_norms   = 3;

    Epetra_SerialDenseVector cpu_dom_norms_bg(num_dom_norms);
    Epetra_SerialDenseVector cpu_dom_norms_emb(num_dom_norms);
    Epetra_SerialDenseVector cpu_interf_norms(num_interf_norms);
    Epetra_SerialDenseVector cpu_stab_norms(num_stab_norms);

    Teuchos::RCP<Epetra_SerialDenseVector> glob_dom_norms_bg  = Teuchos::rcp(new Epetra_SerialDenseVector(num_dom_norms));
    Teuchos::RCP<Epetra_SerialDenseVector> glob_dom_norms_emb = Teuchos::rcp(new Epetra_SerialDenseVector(num_dom_norms));
    Teuchos::RCP<Epetra_SerialDenseVector> glob_interf_norms  = Teuchos::rcp(new Epetra_SerialDenseVector(num_interf_norms));
    Teuchos::RCP<Epetra_SerialDenseVector> glob_stab_norms    = Teuchos::rcp(new Epetra_SerialDenseVector(num_stab_norms));

    // set vector values needed by elements
    bgdis_->ClearState();
    bgdis_->SetState("u and p at time n+1 (converged)", state_->velnp_);

    embdis_->ClearState();
    embdis_->SetState("velaf", alevelnp_);
    //embdis_->SetState("dispnp", aledispnp_);

    boundarydis_->ClearState();
    boundarydis_->SetState("ivelnp", ivelnp_);
    boundarydis_->SetState("idispnp",idispnp_);

    // evaluate domain error norms and interface/boundary error norms at XFEM-interface
    // loop row elements
    const int numrowele = bgdis_->NumMyRowElements();
    for (int i=0; i<numrowele; ++i)
    {
      // local element-wise squared error norms
      Epetra_SerialDenseVector ele_dom_norms_bg(num_dom_norms);
      Epetra_SerialDenseVector ele_interf_norms(num_interf_norms);


      // pointer to current element
      DRT::Element* actele = bgdis_->lRowElement(i);

      Teuchos::RCP<MAT::Material> mat = actele->Material();

      DRT::ELEMENTS::Fluid * ele = dynamic_cast<DRT::ELEMENTS::Fluid *>( actele );

      GEO::CUT::ElementHandle * e = state_->wizard_->GetElement( actele );
      DRT::Element::LocationArray la( 1 );

      DRT::ELEMENTS::FluidEleInterface * impl = DRT::ELEMENTS::FluidFactory::ProvideImplXFEM( actele->Shape(), "xfem");

      // xfem element
      if ( e!=NULL )
      {
#ifdef DOFSETS_NEW

        std::vector< GEO::CUT::plain_volumecell_set > cell_sets;
        std::vector< std::vector<int> > nds_sets;
        std::vector<std::vector< DRT::UTILS::GaussIntegration > >intpoints_sets;

        e->GetCellSets_DofSets_GaussPoints( cell_sets, nds_sets, intpoints_sets, VolumeCellGaussPointBy_ );

        if(cell_sets.size() != intpoints_sets.size()) dserror("number of cell_sets and intpoints_sets not equal!");
        if(cell_sets.size() != nds_sets.size()) dserror("number of cell_sets and nds_sets not equal!");

        int set_counter = 0;


        // loop over volume cells
        for( std::vector< GEO::CUT::plain_volumecell_set>::iterator s=cell_sets.begin();
            s!=cell_sets.end();
            s++)
        {
          // needed for fluid-fluid Coupling
          std::map<int, std::vector<Epetra_SerialDenseMatrix> >  side_coupling;

          GEO::CUT::plain_volumecell_set & cells = *s;
          const std::vector<int> & nds = nds_sets[set_counter];

          // get element location vector, dirichlet flags and ownerships
          actele->LocationVector(*bgdis_,nds,la,false);

          for( unsigned cellcount=0;cellcount!=cell_sets[set_counter].size();cellcount++)
          {
            //------------------------------------------------------------
            // Evaluate domain integral errors
            impl->ComputeError(ele,
                               *params_,
                               mat,
                               *bgdis_,
                               la[0].lm_,
                               ele_dom_norms_bg,
                               intpoints_sets[set_counter][cellcount]
            );

            // sum up (on each processor)
            cpu_dom_norms_bg += ele_dom_norms_bg;


            //------------------------------------------------------------
            // Evaluate interface integral errors
            // do cut interface condition

            // maps of sid and corresponding boundary cells ( for quadratic elements: collected via volumecells of subelements)
            std::map<int, std::vector<GEO::CUT::BoundaryCell*> > bcells;
            std::map<int, std::vector<DRT::UTILS::GaussIntegration> > bintpoints;

            for ( GEO::CUT::plain_volumecell_set::iterator i=cells.begin(); i!=cells.end(); ++i )
            {
              GEO::CUT::VolumeCell * vc = *i;
              if ( vc->Position()==GEO::CUT::Point::outside )
              {
                  vc->GetBoundaryCells( bcells );
              }
            }

            if ( bcells.size() > 0 )
            {
                TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluid::XFluidState::Evaluate 2) interface" );

                // Attention: switch also the flag in fluid_ele_calc_xfem.cpp
#ifdef BOUNDARYCELL_TRANSFORMATION_OLD
                // original Axel's transformation
                e->BoundaryCellGaussPoints( wizard_->CutWizard().Mesh(), 0, bcells, bintpoints );
#else
                // new Benedikt's transformation
                e->BoundaryCellGaussPointsLin( state_->wizard_->CutWizard().Mesh(), 0, bcells, bintpoints );
#endif

                std::vector<int> patchelementslm;
                std::vector<int> patchelementslmowner;

                // initialize the coupling matrices for each side and the current element
                for ( std::map<int,  std::vector<GEO::CUT::BoundaryCell*> >::const_iterator bc=bcells.begin();
                      bc!=bcells.end(); ++bc )
                {
                  int sid = bc->first; // all boundary cells within the current iterator belong to the same side
                  DRT::Element * side = boundarydis_->gElement( sid );

                  std::vector<int> patchlm;
                  std::vector<int> patchlmowner;
                  std::vector<int> patchlmstride;
                  // for nitsche embedded and two-sided we couple with the whole embedded element not only with its side
                  if (coupling_approach_ == CouplingMHCS_XFluid or
                      coupling_approach_ == CouplingNitsche_XFluid or
                      coupling_approach_ == CouplingMHVS_XFluid)
                    side->LocationVector(*boundarydis_, patchlm, patchlmowner, patchlmstride);
                  else if(coupling_approach_ == CouplingNitsche_EmbFluid or coupling_approach_ == CouplingNitsche_TwoSided)
                  {
                    // get the corresponding embedded element for nitsche
                    // embedded and two-sided
                    int emb_eid = boundary_emb_gid_map_.find(sid)->second;
                    DRT::Element * emb_ele = embdis_->gElement( emb_eid );
                    emb_ele->LocationVector(*embdis_, patchlm, patchlmowner, patchlmstride);
                  }

                  patchelementslm.reserve( patchelementslm.size() + patchlm.size());
                  patchelementslm.insert(patchelementslm.end(), patchlm.begin(), patchlm.end());

                  patchelementslmowner.reserve( patchelementslmowner.size() + patchlmowner.size());
                  patchelementslmowner.insert( patchelementslmowner.end(), patchlmowner.begin(), patchlmowner.end());

                  const size_t ndof_i = patchlm.size();     // sum over number of dofs of all sides
                  const size_t ndof   = la[0].lm_.size();   // number of dofs for background element

                  std::vector<Epetra_SerialDenseMatrix> & couplingmatrices = side_coupling[sid];
                  if ( couplingmatrices.size()!=0 )
                    dserror("zero sized vector expected");
                  couplingmatrices.resize(3);

                  // no coupling for pressure in stress based method, but the coupling matrices include entries for pressure coupling
                  couplingmatrices[0].Shape(ndof_i,ndof);  //C_uiu
                  couplingmatrices[1].Shape(ndof,ndof_i);  //C_uui
                  couplingmatrices[2].Shape(ndof_i,1);     //rhC_ui
                }

                const size_t nui = patchelementslm.size();
                Epetra_SerialDenseMatrix  Cuiui(nui,nui);


                if(coupling_method_ == INPAR::XFEM::Hybrid_LM_Cauchy_stress or
                   coupling_method_ == INPAR::XFEM::Hybrid_LM_viscous_stress or
                   coupling_method_ == INPAR::XFEM::Nitsche)
                {
                  impl->ComputeErrorInterfacefluidfluidcoupling(
                      ele,
                      *bgdis_,
                      la[0].lm_,
                      mat,
                      ele_interf_norms,
                      *boundarydis_,
                      *embdis_,
                      bcells,
                      bintpoints,
                      side_coupling,
                      *params_,
                      cells,
                      boundary_emb_gid_map_);
                }
            } // bcells
          }

          set_counter += 1;
        }

#else
        GEO::CUT::plain_volumecell_set cells;
        std::vector<DRT::UTILS::GaussIntegration> intpoints;
        std::vector<std::vector<double> > refEqns;
        e->VolumeCellGaussPoints( cells, intpoints ,refEqns, VolumeCellGaussPointBy_);//modify gauss type

        int count = 0;
        for ( GEO::CUT::plain_volumecell_set::iterator s=cells.begin(); s!=cells.end(); ++s )
        {
          GEO::CUT::VolumeCell * vc = *s;
          if ( vc->Position()==GEO::CUT::Point::outside )
          {
            //             // one set of dofs
            //             std::vector<int>  ndstest;
            //             for (int t=0;t<8; ++t)
            //               ndstest.push_back(0);

            const std::vector<int> & nds = vc->NodalDofSet();
            actele->LocationVector(*bgdis_,nds,la,false);

            impl->ComputeError(ele,
                               *params_,
                               mat,
                               *bgdis_
                               la[0].lm_,
                               ele_dom_norms_bg,
                               intpoints[count]
            );

            // sum up (on each processor)
            cpu_dom_norms_bg += ele_dom_norms_bg;
          }
          count += 1;
        }

#endif

        // sum up (on each processor)
        cpu_interf_norms += ele_interf_norms;
      }
      // standard (no xfem) element
      else
      {
        TEUCHOS_FUNC_TIME_MONITOR( "FLD::XFluidFluid::EvaluateErrorComparedToAnalyticalSol::Evaluate normal" );

        // get element location vector, dirichlet flags and ownerships
        actele->LocationVector(*bgdis_,la,false);

        Epetra_SerialDenseMatrix elemat1;
        Epetra_SerialDenseMatrix elemat2;
        Epetra_SerialDenseVector elevec2;
        Epetra_SerialDenseVector elevec3;
        params_->set<int>("action",FLD::calc_fluid_error);

        DRT::ELEMENTS::FluidFactory::ProvideImplXFEM( actele->Shape(), "xfem")->EvaluateService(ele,
                                                                                             *params_,
                                                                                             mat,
                                                                                             *bgdis_,
                                                                                             la[0].lm_,
                                                                                             elemat1,
                                                                                             elemat2,
                                                                                             ele_dom_norms_bg,
                                                                                             elevec2,
                                                                                             elevec3);

        // sum up (on each processor)
        cpu_dom_norms_bg += ele_dom_norms_bg;

        // no interface norms on non-xfem elements
      }
    }//end loop over bg-fluid elements

    //-----------------------------------------------
    // Embedded discretization
    //---------------------------------------------
    // set vector values needed by elements
    embdis_->ClearState();
    embdis_->SetState("u and p at time n+1 (converged)", alevelnp_);
    //embdis_->SetState("velaf", alevelnp_);

    // evaluate domain error norms and interface/boundary error norms at XFEM-interface
    // loop row elements
    const int numrowele_emb = embdis_->NumMyRowElements();
    for (int i=0; i<numrowele_emb; ++i)
    {

      // local element-wise squared error norms
      Epetra_SerialDenseVector ele_dom_norms_emb(num_dom_norms);

      // pointer to current element
      DRT::Element* actele = embdis_->lRowElement(i);

      Teuchos::RCP<MAT::Material> mat = actele->Material();

      DRT::ELEMENTS::Fluid * ele = dynamic_cast<DRT::ELEMENTS::Fluid *>( actele );

      DRT::Element::LocationArray la( 1 );

//      DRT::ELEMENTS::FluidEleInterface * impl = DRT::ELEMENTS::FluidFactory::ProvideImplXFEM( actele->Shape(), "xfem");

      // get element location vector, dirichlet flags and ownerships
      actele->LocationVector(*embdis_,la,false);

      Epetra_SerialDenseMatrix elemat1;
      Epetra_SerialDenseMatrix elemat2;
      Epetra_SerialDenseVector elevec2;
      Epetra_SerialDenseVector elevec3;
      params_->set<int>("action",FLD::calc_fluid_error);

      DRT::ELEMENTS::FluidFactory::ProvideImplXFEM( actele->Shape(), "xfem")->EvaluateService(ele,
                                                                                           *params_,
                                                                                           mat,
                                                                                           *embdis_,
                                                                                           la[0].lm_,
                                                                                           elemat1,
                                                                                           elemat2,
                                                                                           ele_dom_norms_emb,
                                                                                           elevec2,
                                                                                           elevec3);

      // sum up (on each processor)
      cpu_dom_norms_emb += ele_dom_norms_emb;

    } // end loop over embedded fluid elements

    //--------------------------------------------------------
    // reduce and sum over all procs
    for (int i=0; i<num_dom_norms; ++i) (*glob_dom_norms_bg)(i) = 0.0;
    bgdis_->Comm().SumAll(cpu_dom_norms_bg.Values(), glob_dom_norms_bg->Values(), num_dom_norms);

    for (int i=0; i<num_dom_norms; ++i) (*glob_dom_norms_emb)(i) = 0.0;
    embdis_->Comm().SumAll(cpu_dom_norms_emb.Values(), glob_dom_norms_emb->Values(), num_dom_norms);

    for (int i=0; i<num_interf_norms; ++i) (*glob_interf_norms)(i) = 0.0;
    bgdis_->Comm().SumAll(cpu_interf_norms.Values(), glob_interf_norms->Values(), num_interf_norms);

    // standard domain errors bg-dis
    double dom_bg_err_vel_L2      = 0.0;         //  || u - u_b ||_L2(Omega)           =   standard L2-norm for velocity
    double dom_bg_err_vel_H1_semi = 0.0;         //  || grad( u - u_b ) ||_L2(Omega)   =   standard H1-seminorm for velocity
    double dom_bg_err_vel_H1      = 0.0;         //  || u - u_b ||_H1(Omega)           =   standard H1-norm for velocity
    double dom_bg_err_pre_L2      = 0.0;         //  || p - p_b ||_L2(Omega)           =   standard L2-norm for for pressure

    // viscosity-scaled domain errors
    double dom_bg_err_vel_H1_semi_nu_scaled = 0.0;  //  || nu^(+1/2) grad( u - u_b ) ||_L2(Omega)  =   visc-scaled H1-seminorm for velocity
    double dom_bg_err_pre_L2_nu_scaled      = 0.0;  //  || nu^(-1/2) (p - p_b) ||_L2(Omega)        =   visc-scaled L2-norm for for pressure

    // standard domain errors bg-dis
    double dom_emb_err_vel_L2      = 0.0;         //  || u - u_e ||_L2(Omega)           =   standard L2-norm for velocity
    double dom_emb_err_vel_H1_semi = 0.0;         //  || grad( u - u_e ) ||_L2(Omega)   =   standard H1-seminorm for velocity
    double dom_emb_err_vel_H1      = 0.0;         //  || u - u_e ||_H1(Omega)           =   standard H1-norm for velocity
    double dom_emb_err_pre_L2      = 0.0;         //  || p - p_e ||_L2(Omega)           =   standard L2-norm for for pressure

    // viscosity-scaled domain errors
    double dom_emb_err_vel_H1_semi_nu_scaled = 0.0;  //  || nu^(+1/2) grad( u - u_e ) ||_L2(Omega)  =   visc-scaled H1-seminorm for velocity
    double dom_emb_err_pre_L2_nu_scaled      = 0.0;  //  || nu^(-1/2) (p - p_e) ||_L2(Omega)        =   visc-scaled L2-norm for for pressure

    // interface errors
    double interf_err_Honehalf    = 0.0;         //  || nu^(+1/2) (u_b - u_e) ||_H1/2(Gamma)          =  broken H1/2 Sobolev norm for boundary/coupling condition
    double interf_err_Hmonehalf_u = 0.0;         //  || nu^(+1/2) grad( u_b - u_e )*n ||_H-1/2(Gamma) =  broken H-1/2 Sobolev norm for normal flux (velocity part)
    double interf_err_Hmonehalf_p = 0.0;         //  || nu^(-1/2) (p_b - p_e)*n ||_H-1/2(Gamma)         =  broken H-1/2 Sobolev norm for normal flux (pressure part)

    dom_bg_err_vel_L2             = sqrt((*glob_dom_norms_bg)[0]);
    dom_bg_err_vel_H1_semi        = sqrt((*glob_dom_norms_bg)[1]);
    dom_bg_err_vel_H1             = sqrt((*glob_dom_norms_bg)[2]);
    dom_bg_err_pre_L2             = sqrt((*glob_dom_norms_bg)[3]);

    dom_bg_err_vel_H1_semi_nu_scaled = sqrt((*glob_dom_norms_bg)[4]);
    dom_bg_err_pre_L2_nu_scaled      = sqrt((*glob_dom_norms_bg)[5]);

    dom_emb_err_vel_L2             = sqrt((*glob_dom_norms_emb)[0]);
    dom_emb_err_vel_H1_semi        = sqrt((*glob_dom_norms_emb)[1]);
    dom_emb_err_vel_H1             = sqrt((*glob_dom_norms_emb)[2]);
    dom_emb_err_pre_L2             = sqrt((*glob_dom_norms_emb)[3]);

    dom_emb_err_vel_H1_semi_nu_scaled = sqrt((*glob_dom_norms_emb)[4]);
    dom_emb_err_pre_L2_nu_scaled      = sqrt((*glob_dom_norms_emb)[5]);

    interf_err_Honehalf           = sqrt((*glob_interf_norms)[0]);
    interf_err_Hmonehalf_u        = sqrt((*glob_interf_norms)[1]);
    interf_err_Hmonehalf_p        = sqrt((*glob_interf_norms)[2]);

    if (myrank_ == 0)
    {
      {
        cout.precision(8);
        IO::cout << IO::endl << "---- error norm for analytical solution Nr. "
             <<  DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_,"calculate error")
             <<  " ----------" << IO::endl;
        IO::cout << "-------------- domain error norms (background)------------"        << IO::endl;
        IO::cout << "|| u - u_b ||_L2(Omega)                        =  " << dom_bg_err_vel_L2                << IO::endl;
        IO::cout << "|| grad( u - u_b ) ||_L2(Omega)                =  " << dom_bg_err_vel_H1_semi           << IO::endl;
        IO::cout << "|| u - u_b ||_H1(Omega)                        =  " << dom_bg_err_vel_H1                << IO::endl;
        IO::cout << "|| p - p_b ||_L2(Omega)                        =  " << dom_bg_err_pre_L2                << IO::endl;
        IO::cout << "-------------- domain error norms (embedded)  ------------"       << IO::endl;
        IO::cout << "|| u - u_e ||_L2(Omega)                        =  " << dom_emb_err_vel_L2               << IO::endl;
        IO::cout << "|| grad( u_ - u_h ) ||_L2(Omega)               =  " << dom_emb_err_vel_H1_semi          << IO::endl;
        IO::cout << "|| u - u_e ||_H1(Omega)                        =  " << dom_emb_err_vel_H1               << IO::endl;
        IO::cout << "|| p - p_e ||_L2(Omega)                        =  " << dom_emb_err_pre_L2               << IO::endl;
        IO::cout << "----viscosity-scaled domain error norms (background)------"       << IO::endl;
        IO::cout << "|| nu^(+1/2) grad( u - u_b ) ||_L2(Omega)      =  " << dom_bg_err_vel_H1_semi_nu_scaled << IO::endl;
        IO::cout << "|| nu^(-1/2) (p - p_b) ||_L2(Omega)            =  " << dom_bg_err_pre_L2_nu_scaled      << IO::endl;
        IO::cout << "----viscosity-scaled domain error norms (embedded) ------"       << IO::endl;
        IO::cout << "|| nu^(+1/2) grad( u - u_e ) ||_L2(Omega)      =  " << dom_emb_err_vel_H1_semi_nu_scaled<< IO::endl;
        IO::cout << "|| nu^(-1/2) (p - p_e) ||_L2(Omega)            =  " << dom_emb_err_pre_L2_nu_scaled     << IO::endl;
        IO::cout << "---------------------------------------------------------"       << IO::endl;
        IO::cout << "-------------- interface/boundary error norms -----------"       << IO::endl;
        IO::cout << "|| nu^(+1/2) (u_b - u_e) ||_H1/2(Gamma)            =  " << interf_err_Honehalf          << IO::endl;
        IO::cout << "|| nu^(+1/2) grad( u_b - u_e )*n ||_H-1/2(Gamma)   =  " << interf_err_Hmonehalf_u       << IO::endl;
        IO::cout << "|| nu^(-1/2) (p_b - p_e)*n ||_H-1/2(Gamma)         =  " << interf_err_Hmonehalf_p       << IO::endl;
        IO::cout << "---------------------------------------------------------"       << IO::endl;
        IO::cout << "-------------- Error on Functionals from solution  ------------"       << IO::endl;
        IO::cout << " | sin(x) ( u,x - u,x exact ) | (background)        = " << (*glob_dom_norms_bg)[6]      << IO::endl;
        IO::cout << " | sin(x) ( u,x - u,x exact ) | (embedded)          = " << (*glob_dom_norms_emb)[6]     << IO::endl;
      }

      // append error of the last time step to the error file
      if ((step_==stepmax_) or (time_==maxtime_))// write results to file
      {
        std::ostringstream temp;
        const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
        const std::string fname = simulation+".xfem_abserror";

        std::ofstream f;
        f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
        f << "#| " << simulation << "\n";
        f << "#| Step"
          << " | Time"
          << " | || u - u_b ||_L2(Omega)"
          << " | || grad( u - u_b ) ||_L2(Omega)"
          << " | || u - u_b ||_H1(Omega)"
          << " | || p - p_b ||_L2(Omega)"
          << " | || nu^(+1/2) grad( u - u_b ) ||_L2(Omega)"
          << " | || nu^(-1/2) ( p - p_b ) ||_L2(Omega)"
          << " | || u - u_e ||_L2(Omega)"
          << " | || grad( u - u_e ) ||_L2(Omega)"
          << " | || u - u_e ||_H1(Omega)"
          << " | || p - p_e ||_L2(Omega)"
          << " | || nu^(+1/2) grad( u - u_e ) ||_L2(Omega)"
          << " | || nu^(-1/2) ( p - p_e) ||_L2(Omega)"
          << " | || nu^(+1/2) (u_b - u_e) ||_H1/2(Gamma)"
          << " | || nu^(+1/2) grad( u_b - u_e )*n ||_H-1/2(Gamma)"
          << " | || nu^(-1/2) (p_b - p_e)*n |_H-1/2(Gamma)"
          << " |  | sin(x) ( u,x - u,x exact ) | (background)"
          << " |  | sin(x) ( u,x - u,x exact ) | (embedded)"
          << " |\n";
        f << step_ << " "
          << time_ << " "
          << dom_bg_err_vel_L2 << " "
          << dom_bg_err_vel_H1_semi << " "
          << dom_bg_err_vel_H1 << " "
          << dom_bg_err_pre_L2 << " "
          << dom_bg_err_vel_H1_semi_nu_scaled << " "
          << dom_bg_err_pre_L2_nu_scaled << " "
          << dom_emb_err_vel_L2 << " "
          << dom_emb_err_vel_H1_semi << " "
          << dom_emb_err_vel_H1 << " "
          << dom_emb_err_pre_L2 << " "
          << dom_emb_err_vel_H1_semi_nu_scaled << " "
          << dom_emb_err_pre_L2_nu_scaled << " "
          << interf_err_Honehalf << " "
          << interf_err_Hmonehalf_u << " "
          << interf_err_Hmonehalf_p << " "
          << (*glob_dom_norms_bg)[6] << " "
          << (*glob_dom_norms_emb)[6] << " "
          <<"\n";
        f.flush();
        f.close();
        }
      std::ostringstream temp;
      const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
      const std::string fname = simulation+"_time.xfem_abserror";

      if(step_==1)
      {
        std::ofstream f;
        f.open(fname.c_str());
        f << "#| " << simulation << "\n";
        f << "#| Step"
          << " | Time"
          << " | || u - u_b ||_L2(Omega)"
          << " | || grad( u - u_b ) ||_L2(Omega)"
          << " | || u - u_b ||_H1(Omega)"
          << " | || p - p_b ||_L2(Omega)"
          << " | || nu^(+1/2) grad( u - u_b ) ||_L2(Omega)"
          << " | || nu^(-1/2) ( p - p_b ) ||_L2(Omega)"
          << " | || u - u_e ||_L2(Omega)"
          << " | || grad( u - u_e ) ||_L2(Omega)"
          << " | || u - u_e ||_H1(Omega)"
          << " | || p - p_e ||_L2(Omega)"
          << " | || nu^(+1/2) grad( u - u_e ) ||_L2(Omega)"
          << " | || nu^(-1/2) ( p - p_e) ||_L2(Omega)"
          << " | || nu^(+1/2) (u_b - u_e) ||_H1/2(Gamma)"
          << " | || nu^(+1/2) grad( u_b - u_e )*n ||_H-1/2(Gamma)"
          << " | || nu^(-1/2) (p_b - p_e)*n |_H-1/2(Gamma)"
          << " |  | sin(x) ( u,x - u,x exact ) | (background)"
          << " |  | sin(x) ( u,x - u,x exact ) | (embedded)"
          << " |\n";
        f << step_ << " "
          << time_ << " "
          << dom_bg_err_vel_L2 << " "
          << dom_bg_err_vel_H1_semi << " "
          << dom_bg_err_vel_H1 << " "
          << dom_bg_err_pre_L2 << " "
          << dom_bg_err_vel_H1_semi_nu_scaled << " "
          << dom_bg_err_pre_L2_nu_scaled << " "
          << dom_emb_err_vel_L2 << " "
          << dom_emb_err_vel_H1_semi << " "
          << dom_emb_err_vel_H1 << " "
          << dom_emb_err_pre_L2 << " "
          << dom_emb_err_vel_H1_semi_nu_scaled << " "
          << dom_emb_err_pre_L2_nu_scaled << " "
          << interf_err_Honehalf << " "
          << interf_err_Hmonehalf_u << " "
          << interf_err_Hmonehalf_p << " "
          << (*glob_dom_norms_bg)[6] << " "
          << (*glob_dom_norms_emb)[6] << " "
          <<"\n";

        f.flush();
        f.close();
      }
      else
      {
        std::ofstream f;
        f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
        f << step_ << " "
          << time_ << " "
          << dom_bg_err_vel_L2 << " "
          << dom_bg_err_vel_H1_semi << " "
          << dom_bg_err_vel_H1 << " "
          << dom_bg_err_pre_L2 << " "
          << dom_bg_err_vel_H1_semi_nu_scaled << " "
          << dom_bg_err_pre_L2_nu_scaled << " "
          << dom_emb_err_vel_L2 << " "
          << dom_emb_err_vel_H1_semi << " "
          << dom_emb_err_vel_H1 << " "
          << dom_emb_err_pre_L2 << " "
          << dom_emb_err_vel_H1_semi_nu_scaled << " "
          << dom_emb_err_pre_L2_nu_scaled << " "
          << interf_err_Honehalf << " "
          << interf_err_Hmonehalf_u << " "
          << interf_err_Hmonehalf_p << " "
          << (*glob_dom_norms_bg)[6] << " "
          << (*glob_dom_norms_emb)[6] << " "
          <<"\n";

        f.flush();
        f.close();
      }
    } // myrank = 0

  }
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void FLD::XFluidFluid::SetInitialFlowField(
  const INPAR::FLUID::InitialField initfield,
  const int startfuncno
  )
{
  // initial field by (undisturbed) function (init==2)
  // or disturbed function (init==3)
  if (initfield == INPAR::FLUID::initfield_field_by_function or
      initfield == INPAR::FLUID::initfield_disturbed_field_from_function)
  {
    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<bgdis_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = bgdis_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      const std::vector<int> nodedofset = bgdis_->Dof(lnode);

      if (nodedofset.size()!=0)
      {
        for(int index=0;index<numdim_+1;++index)
        {
          int gid = nodedofset[index];

          double initialval=DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(index,lnode->X(),time_,NULL);
          state_->velnp_->ReplaceGlobalValues(1,&initialval,&gid);
        }
      }
    }

    // initialize veln_ as well.
    state_->veln_->Update(1.0,*state_->velnp_ ,0.0);

    // loop all nodes of embedded fluid on the processor
    for(int lnodeid=0;lnodeid<embdis_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = embdis_->lRowNode(lnodeid);
      // the set of degrees of freedom associated with the node
      const std::vector<int> nodedofset = embdis_->Dof(lnode);

      for(int index=0;index<numdim_+1;++index)
      {
        int gid = nodedofset[index];

        double initialval=DRT::Problem::Instance()->Funct(startfuncno-1).Evaluate(index,lnode->X(),time_,NULL);

        alevelnp_->ReplaceGlobalValues(1,&initialval,&gid);
      }
    }

    // initialize veln_ as well.
    aleveln_->Update(1.0,*alevelnp_ ,0.0);
    LINALG::Export(*(alevelnp_),*(ivelnp_));
  }

  // special initial function: Beltrami flow (3-D)
  else if (initfield == INPAR::FLUID::initfield_beltrami_flow)
  {
    const Epetra_Map* dofrowmap = bgdis_->DofRowMap();

    int err = 0;

    const int npredof = numdim_;

    double         p;
    std::vector<double> u  (numdim_);
    std::vector<double> xyz(numdim_);

    // check whether present flow is indeed three-dimensional
    if (numdim_!=3) dserror("Beltrami flow is a three-dimensional flow!");

    // set constants for analytical solution
    const double a = M_PI/4.0;
    const double d = M_PI/2.0;

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<bgdis_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = bgdis_->lRowNode(lnodeid);

      // the set of degrees of freedom associated with the node
      std::vector<int> nodedofset = bgdis_->Dof(lnode);

      //there wont be any dof for nodes which are inside the structure
      //the cut algorithm erases these dofs
      if(nodedofset.size()==0)
        continue;

      // set node coordinates
      for(int dim=0;dim<numdim_;dim++)
      {
        xyz[dim]=lnode->X()[dim];
      }

      // compute initial velocity components
      u[0] = -a * ( exp(a*xyz[0]) * sin(a*xyz[1] + d*xyz[2]) +
                    exp(a*xyz[2]) * cos(a*xyz[0] + d*xyz[1]) );
      u[1] = -a * ( exp(a*xyz[1]) * sin(a*xyz[2] + d*xyz[0]) +
                    exp(a*xyz[0]) * cos(a*xyz[1] + d*xyz[2]) );
      u[2] = -a * ( exp(a*xyz[2]) * sin(a*xyz[0] + d*xyz[1]) +
                    exp(a*xyz[1]) * cos(a*xyz[2] + d*xyz[0]) );

      // compute initial pressure
      int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_fluid);
      if (id==-1) dserror("Newtonian fluid material could not be found");
      const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
      const MAT::PAR::NewtonianFluid* actmat = static_cast<const MAT::PAR::NewtonianFluid*>(mat);
      double dens = actmat->density_;
      p = -a*a/2.0 * dens *
        ( exp(2.0*a*xyz[0])
          + exp(2.0*a*xyz[1])
          + exp(2.0*a*xyz[2])
          + 2.0 * sin(a*xyz[0] + d*xyz[1]) * cos(a*xyz[2] + d*xyz[0]) * exp(a*(xyz[1]+xyz[2]))
          + 2.0 * sin(a*xyz[1] + d*xyz[2]) * cos(a*xyz[0] + d*xyz[1]) * exp(a*(xyz[2]+xyz[0]))
          + 2.0 * sin(a*xyz[2] + d*xyz[0]) * cos(a*xyz[1] + d*xyz[2]) * exp(a*(xyz[0]+xyz[1]))
          );

      // set initial velocity components
      for(int nveldof=0;nveldof<numdim_;nveldof++)
      {
        const int gid = nodedofset[nveldof];
        int lid = dofrowmap->LID(gid);
        err += state_->velnp_->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += state_->veln_ ->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err += state_->velnm_->ReplaceMyValues(1,&(u[nveldof]),&lid);
      }

      // set initial pressure
      const int gid = nodedofset[npredof];
      int lid = dofrowmap->LID(gid);
      err += state_->velnp_->ReplaceMyValues(1,&p,&lid);
      err += state_->veln_ ->ReplaceMyValues(1,&p,&lid);
      err += state_->velnm_->ReplaceMyValues(1,&p,&lid);
    } // end loop nodes lnodeid

    if (err!=0) dserror("dof not on proc");

    const Epetra_Map* dofrowmap1 = embdis_->DofRowMap();

    int err1 = 0;

    // check whether present flow is indeed three-dimensional
    if (numdim_!=3) dserror("Beltrami flow is a three-dimensional flow!");

    // loop all nodes on the processor
    for(int lnodeid=0;lnodeid<embdis_->NumMyRowNodes();lnodeid++)
    {
      // get the processor local node
      DRT::Node*  lnode      = embdis_->lRowNode(lnodeid);

      // the set of degrees of freedom associated with the node
      std::vector<int> nodedofset = embdis_->Dof(lnode);

      //there wont be any dof for nodes which are inside the structure
      //the cut algorithm erases these dofs
      if(nodedofset.size()==0)
        continue;

      // set node coordinates
      for(int dim=0;dim<numdim_;dim++)
      {
        xyz[dim]=lnode->X()[dim];
      }

          // compute initial velocity components
      u[0] = -a * ( exp(a*xyz[0]) * sin(a*xyz[1] + d*xyz[2]) +
                    exp(a*xyz[2]) * cos(a*xyz[0] + d*xyz[1]) );
      u[1] = -a * ( exp(a*xyz[1]) * sin(a*xyz[2] + d*xyz[0]) +
                        exp(a*xyz[0]) * cos(a*xyz[1] + d*xyz[2]) );
      u[2] = -a * ( exp(a*xyz[2]) * sin(a*xyz[0] + d*xyz[1]) +
                    exp(a*xyz[1]) * cos(a*xyz[2] + d*xyz[0]) );

      // compute initial pressure
      p = -a*a/2.0 *
          ( exp(2.0*a*xyz[0])
            + exp(2.0*a*xyz[1])
            + exp(2.0*a*xyz[2])
            + 2.0 * sin(a*xyz[0] + d*xyz[1]) * cos(a*xyz[2] + d*xyz[0]) * exp(a*(xyz[1]+xyz[2]))
            + 2.0 * sin(a*xyz[1] + d*xyz[2]) * cos(a*xyz[0] + d*xyz[1]) * exp(a*(xyz[2]+xyz[0]))
            + 2.0 * sin(a*xyz[2] + d*xyz[0]) * cos(a*xyz[1] + d*xyz[2]) * exp(a*(xyz[0]+xyz[1]))
            );

      // set initial velocity components
      for(int nveldof=0;nveldof<numdim_;nveldof++)
      {
        const int gid = nodedofset[nveldof];
        int lid = dofrowmap1->LID(gid);
        err1 += alevelnp_->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err1 += alevelaf_ ->ReplaceMyValues(1,&(u[nveldof]),&lid);
        err1 += aleveln_->ReplaceMyValues(1,&(u[nveldof]),&lid);
      }

      // set initial pressure
      const int gid = nodedofset[npredof];
      int lid = dofrowmap1->LID(gid);
      err1 += alevelnp_->ReplaceMyValues(1,&p,&lid);
      err1 += alevelaf_ ->ReplaceMyValues(1,&p,&lid);
      err1 += aleveln_->ReplaceMyValues(1,&p,&lid);

    } // end loop nodes lnodeid

    //    state_->veln_->Print(std::cout);
    //    state_->GmshOutput(*bgdis_, *embdis_, *boundarydis_, "initial", -1, step_, state_->veln_, aleveln_, aledispnp_);

    if (err1!=0) dserror("dof not on proc");
  }

  else
  {
    dserror("Only initial fields auch as a zero field, initial fields by (un-)disturbed functions and  Beltrami flow!");
  }

  return;
} // end SetInitialFlowField

/*------------------------------------------------------------------------------------------------*
 | use shapederivatives as block matrix
 | for the case of monolithic fluid-fluid fsi with fluidsplit
 *------------------------------------------------------------------------------------------------*/
void FLD::XFluidFluid::UseBlockMatrix(Teuchos::RCP<std::set<int> >     condelements,
                                      const LINALG::MultiMapExtractor& shapederivdomainmaps,
                                      const LINALG::MultiMapExtractor& shapederivrangemaps)
{
  // no shapederivatives active
  if (! params_->get<bool>("shape derivatives"))
    return;

  // here we initialize the shapederivates
  // REMARK: in case of monolithic fluid-fluid fsi the passed map extractors contain
  // background fluid dof. that's ok, as the entire matrix is first set to zero before we reach
  // the evaluation loop for subsequent embedded elements. there is no chance to accidentally
  // manipulate a background fluid entry...
  // kruse 05/2014
  Teuchos::RCP<LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy> > mat =
    Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>(shapederivdomainmaps,shapederivrangemaps,108,false,true));
  mat->SetCondElements(condelements);
  shapederivatives_ = mat;
}

/*------------------------------------------------------------------------------------------------*
 | create a result test
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> FLD::XFluidFluid::CreateFieldTest()
{
  return Teuchos::rcp(new FLD::XFluidFluidResultTest(*this));
}

/*------------------------------------------------------------------------------------------------*
 | get coupled fluid-fluid system matrix
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> FLD::XFluidFluid::SystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(state_->fluidfluidsysmat_);
}

// -------------------------------------------------------------------
// extrapolate from time mid-point to end-point         (mayr 12/2011)
// -------------------------------------------------------------------
Teuchos::RCP<Epetra_Vector> FLD::XFluidFluid::ExtrapolateEndPoint
(
  Teuchos::RCP<Epetra_Vector> vecn,
  Teuchos::RCP<Epetra_Vector> vecm
)
{
  Teuchos::RCP<Epetra_Vector> vecnp = Teuchos::rcp(new Epetra_Vector(*vecm));

//   // For gen-alpha extrapolate mid-point quantities to end-point.
//   // Otherwise, equilibrium time level is already end-point.
//   if (is_genalpha_)
//     vecnp->Update((alphaF_-1.0)/alphaF_,*vecn,1.0/alphaF_);

  return vecnp;
}

/*------------------------------------------------------------------------------------------------*
 | get fluid-fluid system matrix as block matrix
 | for the case of monolithic fluid-fluid fsi with fluidsplit
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> FLD::XFluidFluid::BlockSystemMatrix(
    const Teuchos::RCP<const LINALG::MultiMapExtractor> fsiextractor,
    Teuchos::RCP<Epetra_Map> innermap,
    Teuchos::RCP<Epetra_Map> condmap)
{
  //Map of fluid FSI DOFs: condmap
  //Map of inner fluid DOFs: innermap

  //Get the fluid-fluid system matrix as sparse matrix
  Teuchos::RCP<LINALG::SparseMatrix> sparsesysmat = SystemMatrix();

  //F_{II}
  Teuchos::RCP<LINALG::SparseMatrix> fii;
  //F_{I\Gamma}
  Teuchos::RCP<LINALG::SparseMatrix> fig;
  //F_{\GammaI}
  Teuchos::RCP<LINALG::SparseMatrix> fgi;
  //F_{\Gamma\Gamma}
  Teuchos::RCP<LINALG::SparseMatrix> fgg;

  // The method expects non-const objects
  // Filling the empty block matrix with entries from sparsesysmat:
  LINALG::SplitMatrix2x2(sparsesysmat,condmap,innermap,condmap,innermap,fgg,fgi,fig,fii);

  blockmat_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*fsiextractor,*fsiextractor,108,false,true));

  blockmat_->Assign(0,0,View,*fii);
  blockmat_->Assign(0,1,View,*fig);
  blockmat_->Assign(1,0,View,*fgi);
  blockmat_->Assign(1,1,View,*fgg);

  blockmat_->Complete();

  if ( blockmat_ == Teuchos::null )
    dserror("Creation of fluid-fluid block matrix failed.");

  return blockmat_;
}

/*------------------------------------------------------------------------------------------------*
 | remove dirichlet-constrained dof from dirichlet map
 *------------------------------------------------------------------------------------------------*/
void FLD::XFluidFluid::RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove)
{
  // the removal works as follows:
  // add the map of dirichlet dof that shall be removed
  // to othermap and setup with the new othermap and iscondmap=false
	std::vector<Teuchos::RCP<const Epetra_Map> > othermaps;
	othermaps.push_back(maptoremove);
	othermaps.push_back(aledbcmaps_->OtherMap());
	// the new map of non-dirichlet dof
	Teuchos::RCP<Epetra_Map> new_othermap = LINALG::MultiMapExtractor::MergeMaps(othermaps);
	*(aledbcmaps_)=LINALG::MapExtractor(*(embdis_->DofRowMap()),new_othermap,false);

  // update the combined fluid dirichlet maps immediately
  state_->CreateFluidFluidDBCMaps();
}

/*------------------------------------------------------------------------------------------------*
 | prepare fluid-fluid coupling matrices for evaluation routine
 *------------------------------------------------------------------------------------------------*/
void FLD::XFluidFluid::XFluidFluidState::PrepareCouplingMatrices(bool initial_call)
{
  // allocate, if this is the first call after construction.
  // the boundary and embedded fluid dof maps will never change,
  // and the background fluid map also won't change during the lifetime of this state object
  if (initial_call)
  {
    if (xfluid_.coupling_approach_ == CouplingMHCS_XFluid or
        xfluid_.coupling_approach_ == CouplingNitsche_XFluid or
        xfluid_.coupling_approach_ == CouplingMHVS_XFluid)
    {
      Cuui_  = Teuchos::rcp(new LINALG::SparseMatrix(*fluiddofrowmap_,0,true,true,LINALG::SparseMatrix::FE_MATRIX));
      Cuiu_  = Teuchos::rcp(new LINALG::SparseMatrix(*xfluid_.boundarydofrowmap_,0,true,true,LINALG::SparseMatrix::FE_MATRIX));
      Cuiui_ = Teuchos::rcp(new LINALG::SparseMatrix(*xfluid_.boundarydofrowmap_,0,true,true,LINALG::SparseMatrix::FE_MATRIX));
      rhC_ui_= LINALG::CreateVector(*xfluid_.boundarydofrowmap_,true);
    }
    else if (xfluid_.coupling_approach_ == CouplingNitsche_EmbFluid or xfluid_.coupling_approach_ == CouplingNitsche_TwoSided)
    {
      Cuui_  = Teuchos::rcp(new LINALG::SparseMatrix(*fluiddofrowmap_,0,true,true,LINALG::SparseMatrix::FE_MATRIX));
      Cuiu_  = Teuchos::rcp(new LINALG::SparseMatrix(*xfluid_.aledofrowmap_,0,true,true,LINALG::SparseMatrix::FE_MATRIX));
      Cuiui_ = Teuchos::rcp(new LINALG::SparseMatrix(*xfluid_.aledofrowmap_,0,true,true,LINALG::SparseMatrix::FE_MATRIX));
      rhC_ui_= LINALG::CreateVector(*xfluid_.aledofrowmap_,true);
    }
    else
      dserror("Unknown fluid-fluid coupling type. Failed to prepare coupling matrices.");

    return;
  }

  // if this is not the initial call, set everything to zero
  Cuui_->Zero();
  Cuiu_->Zero();
  Cuiui_->Zero();
  rhC_ui_->Scale(0.0);
}

/*------------------------------------------------------------------------------------------------*
 | create map of dirichlet-constrained dof from both fluids
 *------------------------------------------------------------------------------------------------*/
void FLD::XFluidFluid::XFluidFluidState::CreateFluidFluidDBCMaps()
{
  // create merged dbc map from both fluids
  std::vector<Teuchos::RCP<const Epetra_Map> > dbcmaps;
  dbcmaps.push_back(dbcmaps_->CondMap());
  dbcmaps.push_back(xfluid_.aledbcmaps_->CondMap());
  fluidfluiddbcmaps_ = LINALG::MultiMapExtractor::MergeMaps(dbcmaps);
}

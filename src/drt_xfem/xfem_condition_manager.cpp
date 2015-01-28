/*!
\file xfem_condition_manager.cpp

\brief manages the different types of mesh and level-set based coupling conditions and thereby builds the bridge between the
xfluid class and the cut-library

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/



#include <Teuchos_TimeMonitor.hpp>

#include "xfem_condition_manager.H"

#include "../drt_inpar/inpar_xfem.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_discret_faces.H"
#include "../drt_lib/drt_discret_xfem.H"
#include "../drt_lib/drt_dofset_transparent_independent.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_parallel.H"


#include "../linalg/linalg_utils.H"

#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_fluid.H"

#include "../drt_crack/crackUtils.H"

#include "../drt_io/io.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

void XFEM::CouplingBase::SetElementConditions()
{
  // number of column cutter boundary elements
  int nummycolele = cutter_dis_->NumMyColElements();

  cutterele_conds_.clear();
  cutterele_conds_.reserve(nummycolele);

  // initialize the vector invalid coupling-condition type "NONE"
  EleCoupCond init_pair = EleCoupCond(INPAR::XFEM::CouplingCond_NONE,NULL);
  for(int lid=0; lid<nummycolele; lid++) cutterele_conds_.push_back(init_pair);

  //-----------------------------------------------------------------------------------
  // loop all column cutting elements on this processor
  for(int lid=0; lid<nummycolele; lid++)
  {
    DRT::Element* cutele = cutter_dis_->lColElement(lid);

    // loop all possible XFEM-coupling conditions
    for(size_t cond=0; cond < conditions_to_copy_.size(); cond++)
    {
      // get all conditions with given condition name
      std::vector<DRT::Condition*> mycond;
      DRT::UTILS::FindElementConditions(cutele, conditions_to_copy_[cond], mycond);

      // safety checks
      if (mycond.size()>1)
      {
//        dserror("%i conditions of the same name for element %i! %s coupling-condition not unique!", mycond.size(), cutele->Id(), conditions_to_copy_[cond].c_str());
      }
      else if(mycond.size() == 0)
      {
        continue; // try the next condition type
      }
      else
      {
//        std::cout << "unique condition found!" << std::endl;
      }

      INPAR::XFEM::EleCouplingCondType cond_type = CondType_stringToEnum(conditions_to_copy_[cond]);

      // non-coupling condition found (e.g. FSI coupling)
      if(cond_type == INPAR::XFEM::CouplingCond_NONE) continue;

      // non-unique conditions for one cutter element
      if( cutterele_conds_[lid].first != INPAR::XFEM::CouplingCond_NONE )
      {
        dserror("There are two different condition types for the same cutter dis element with id %i: 1st %i, 2nd %i. Make the XFEM coupling conditions unique!",
            cutele->Id(), cutterele_conds_[lid].first, cond_type);
      }

      // store the unique condition pointer to the cutting element
      cutterele_conds_[lid] = EleCoupCond(cond_type, mycond[0]);
    }
  }

  //-----------------------------------------------------------------------------------
  // check if all column cutter elements have a valid condition type
  // loop all column cutting elements on this processor
  for(int lid=0; lid<nummycolele; lid++)
  {
    if(cutterele_conds_[lid].first == INPAR::XFEM::CouplingCond_NONE)
      dserror("cutter element with local id %i has no valid coupling-condition", lid);
  }

}

void XFEM::CouplingBase::SetCouplingStrategy()
{
  const INPAR::XFEM::EleCouplingCondType cond_type = CondType_stringToEnum(cond_name_);

  switch(cond_type)
  {
  case INPAR::XFEM::CouplingCond_SURF_FSI_MONO:
  {
    coupling_strategy_ = INPAR::XFEM::Xfluid_Sided_Coupling; //TODO: rename to xfluid (only one for coupling and wdbc/Neumann)
    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID:
  {
    // ask the first cutter element
    const int lid=0;
    const int val = cutterele_conds_[lid].second->GetInt("COUPSTRATEGY");
    coupling_strategy_ = static_cast<INPAR::XFEM::CouplingStrategy>(val);
    break;
  }
  case INPAR::XFEM::CouplingCond_LEVELSET_TWOPHASE:
  case INPAR::XFEM::CouplingCond_LEVELSET_COMBUSTION:
  {
    coupling_strategy_ = INPAR::XFEM::Harmonic;
    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_FSI_PART:
  case INPAR::XFEM::CouplingCond_SURF_CRACK_FSI_PART:
  case INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET: // set this to Teuchos::null when the values are read from the function instead of the ivelnp vector
  case INPAR::XFEM::CouplingCond_SURF_NEUMANN:
  case INPAR::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET:
  case INPAR::XFEM::CouplingCond_LEVELSET_NEUMANN:
  {
    coupling_strategy_ = INPAR::XFEM::Xfluid_Sided_weak_DBC; //TODO: rename to xfluid
    break;
  }
  default: dserror("which is the coupling discretization for this type of coupling %i?", cond_type); break;
  }
}


void XFEM::CouplingBase::SetCouplingDiscretization()
{
  const INPAR::XFEM::EleCouplingCondType cond_type = CondType_stringToEnum(cond_name_);

  switch(cond_type)
  {
  case INPAR::XFEM::CouplingCond_SURF_FSI_MONO:
  {
    coupl_dis_ = cutter_dis_; break;
    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID:
  {
    // depending on the weighting strategy
    if(coupling_strategy_==INPAR::XFEM::Xfluid_Sided_Coupling)
    {
      coupl_dis_ = cutter_dis_;
    }
    else if(coupling_strategy_==INPAR::XFEM::Embedded_Sided_Coupling or
        coupling_strategy_==INPAR::XFEM::Two_Sided_Coupling )
    {
      coupl_dis_ = cond_dis_;
    }
    else dserror("invalid coupling strategy for fluid-fluid application");
    break;
  }
  case INPAR::XFEM::CouplingCond_LEVELSET_TWOPHASE:
  case INPAR::XFEM::CouplingCond_LEVELSET_COMBUSTION:
  {
    coupl_dis_ = bg_dis_;
    break;
  }
  case INPAR::XFEM::CouplingCond_SURF_FSI_PART:
  case INPAR::XFEM::CouplingCond_SURF_CRACK_FSI_PART:
  case INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET: // set this to Teuchos::null when the values are read from the function instead of the ivelnp vector
  case INPAR::XFEM::CouplingCond_SURF_NEUMANN:
  {
    coupl_dis_ = cutter_dis_;
    break;
  }
  case INPAR::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET:
  case INPAR::XFEM::CouplingCond_LEVELSET_NEUMANN:
  {
    coupl_dis_ = Teuchos::null;
    break;
  }
  default: dserror("which is the coupling discretization for this type of coupling %i?", cond_type); break;
  }
}



XFEM::MeshCoupling::MeshCoupling(
    Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
    const std::string &                 cond_name,///< name of the condition, by which the derived cutter discretization is identified
    Teuchos::RCP<DRT::Discretization>&  cond_dis  ///< discretization from which the cutter discretization is derived
) : CouplingBase(bg_dis, cond_name, cond_dis)
{

  // set list of conditions that will be copied to the new cutter discretization
  SetConditionsToCopy();

  // create a cutter discretization from conditioned nodes of the given coupling discretization
  CreateCutterDisFromCondition();

  std::cout << *cutter_dis_ << std::endl;

  // set unique element conditions
  SetElementConditions();

  // set the coupling strategy
  SetCouplingStrategy();

  // set coupling discretization
  SetCouplingDiscretization();

  // initialize state vectors based on cutter discretization
  InitStateVectors();

}


void XFEM::MeshCoupling::SetConditionsToCopy()
{
  // fill list of conditions that will be copied to the new cutter discretization
  conditions_to_copy_.push_back(cond_name_);

  // additional conditions required for the new boundary conditions
  conditions_to_copy_.push_back("FSICoupling");  // for partitioned and monolithic XFSI
}


/*--------------------------------------------------------------------------*
 | Create the cutter discretization                                       |
 *--------------------------------------------------------------------------*/
void XFEM::MeshCoupling::CreateCutterDisFromCondition()
{
  // create name string for new cutter discretization (e.g, "boundary_of_struct" or "boundary_of_fluid")
  std::string cutterdis_name ("boundary_of_");
  cutterdis_name += cond_dis_->Name();


  //--------------------------------
  // create the new cutter discretization form the conditioned coupling discretization
  cutter_dis_ = DRT::UTILS::CreateDiscretizationFromCondition(
      cond_dis_,                ///< discretization with condition
      cond_name_,               ///< name of the condition, by which the derived discretization is identified
      cutterdis_name,           ///< name of the new discretization
      GetBELEName(cond_dis_),   ///< name/type of the elements to be created
      conditions_to_copy_       ///< list of conditions that will be copied to the new discretization
  );
  //--------------------------------


  if (cutter_dis_->NumGlobalNodes() == 0)
  {
    dserror("Empty cutter discretization detected. No coupling can be performed...");
  }

  // for parallel jobs we have to call TransparentDofSet with additional flag true
  bool parallel = cond_dis_->Comm().NumProc() > 1;
  Teuchos::RCP<DRT::DofSet> newdofset = Teuchos::rcp(new
      DRT::TransparentIndependentDofSet(cond_dis_,parallel));

  cutter_dis_->ReplaceDofSet(newdofset); //do not call this with true!!
  cutter_dis_->FillComplete();

  // create node and element distribution with elements and nodes ghosted on all processors
  const Epetra_Map noderowmap = *cutter_dis_->NodeRowMap();
  const Epetra_Map elemrowmap = *cutter_dis_->ElementRowMap();

  // put all boundary nodes and elements onto all processors
  const Epetra_Map nodecolmap = *LINALG::AllreduceEMap(noderowmap);
  const Epetra_Map elemcolmap = *LINALG::AllreduceEMap(elemrowmap);

  // redistribute nodes and elements to column (ghost) map
  cutter_dis_->ExportColumnNodes(nodecolmap);
  cutter_dis_->ExportColumnElements(elemcolmap);

  cutter_dis_->FillComplete();
}


void XFEM::MeshCoupling::InitStateVectors()
{
  const Epetra_Map* cutterdofrowmap = cutter_dis_->DofRowMap();

  ivelnp_  = LINALG::CreateVector(*cutterdofrowmap,true);
  iveln_   = LINALG::CreateVector(*cutterdofrowmap,true);
  ivelnm_  = LINALG::CreateVector(*cutterdofrowmap,true);

  idispnp_ = LINALG::CreateVector(*cutterdofrowmap,true);
  idispn_  = LINALG::CreateVector(*cutterdofrowmap,true);
}

void XFEM::MeshCoupling::SetState()
{ std::cout << "set state!" << std::endl;
  // set general vector values of cutterdis needed by background element evaluate routine
  cutter_dis_->ClearState();

  cutter_dis_->SetState("ivelnp", ivelnp_ );
  cutter_dis_->SetState("iveln",  iveln_  );
  cutter_dis_->SetState("idispnp",idispnp_);
}


Teuchos::RCP<const Epetra_Vector> XFEM::MeshCoupling::GetCutterDispCol()
{
  // export cut-discretization mesh displacements
  Teuchos::RCP<Epetra_Vector> idispcol = LINALG::CreateVector(*cutter_dis_->DofColMap(),true);
  LINALG::Export(*idispnp_,*idispcol);

  return idispcol;
}



//! constructor
XFEM::MeshCouplingFSI::MeshCouplingFSI(
    Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
    const std::string &                 cond_name,///< name of the condition, by which the derived cutter discretization is identified
    Teuchos::RCP<DRT::Discretization>&  cond_dis  ///< discretization from which cutter discretization can be derived
) : MeshCoupling(bg_dis,cond_name,cond_dis), firstoutputofrun_(true)
{
  InitStateVectors_FSI();
  PrepareCutterOutput();
}



void XFEM::MeshCouplingFSI::InitStateVectors_FSI()
{
  const Epetra_Map* cutterdofrowmap = cutter_dis_->DofRowMap();
  const Epetra_Map* cutterdofcolmap = cutter_dis_->DofColMap();

  itrueresidual_ = LINALG::CreateVector(*cutterdofrowmap,true);
  iforcecol_     = LINALG::CreateVector(*cutterdofcolmap,true);
}


void XFEM::MeshCouplingFSI::CompleteStateVectors()
{
  //-------------------------------------------------------------------------------
  // finalize itrueresidual vector

  // need to export the interface forces
  Epetra_Vector iforce_tmp(itrueresidual_->Map(),true);
  Epetra_Export exporter_iforce(iforcecol_->Map(),iforce_tmp.Map());
  int err1 = iforce_tmp.Export(*iforcecol_,exporter_iforce,Add);
  if (err1) dserror("Export using exporter returned err=%d",err1);

  // scale the interface trueresidual with -1.0 to get the forces acting on structural side (no residual-scaling!)
  itrueresidual_->Update(-1.0,iforce_tmp,0.0);
}


void XFEM::MeshCouplingFSI::ZeroStateVectors_FSI()
{
  itrueresidual_->PutScalar(0.0);
  iforcecol_->PutScalar(0.0);
}

// -------------------------------------------------------------------
// Read Restart data for cutter discretization
// -------------------------------------------------------------------
void XFEM::MeshCouplingFSI::ReadRestart(
    const int step
)
{
  const int myrank = cutter_dis_->Comm().MyPID();

  if(myrank) IO::cout << "ReadRestart for boundary discretization " << IO::endl;

  //-------- boundary discretization
  IO::DiscretizationReader boundaryreader(cutter_dis_, step);

  const double time = boundaryreader.ReadDouble("time");
//  const int    step = boundaryreader.ReadInt("step");

  if(myrank == 0)
  {
    IO::cout << "time: " << time << IO::endl;
    IO::cout << "step: " << step << IO::endl;
  }

  boundaryreader.ReadVector(iveln_,   "iveln_res");
  boundaryreader.ReadVector(idispn_,  "idispn_res");

  // REMARK: ivelnp_ and idispnp_ are set again for the new time step in PrepareSolve()
  boundaryreader.ReadVector(ivelnp_,  "ivelnp_res");
  boundaryreader.ReadVector(idispnp_, "idispnp_res");

  if (not (cutter_dis_->DofRowMap())->SameAs(ivelnp_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (cutter_dis_->DofRowMap())->SameAs(iveln_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (cutter_dis_->DofRowMap())->SameAs(idispnp_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (cutter_dis_->DofRowMap())->SameAs(idispn_->Map()))
    dserror("Global dof numbering in maps does not match");


}


void XFEM::MeshCouplingFSI::GmshOutput(
    const std::string & filename_base,
    const int step,
    const int gmsh_step_diff,
    const bool gmsh_debug_out_screen
)
{
  std::ostringstream filename_base_fsi;
  filename_base_fsi << filename_base << "_force";

  // compute the current boundary position
  std::map<int,LINALG::Matrix<3,1> > currinterfacepositions;
  ExtractNodeVectors(cutter_dis_, idispnp_, currinterfacepositions);


  const std::string filename =
      IO::GMSH::GetNewFileNameAndDeleteOldFiles(
          filename_base_fsi.str(),
          step,
          gmsh_step_diff,
          gmsh_debug_out_screen,
          cutter_dis_->Comm().MyPID()
      );

  std::ofstream gmshfilecontent(filename.c_str());

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "iforce \" {" << std::endl;
    // draw vector field 'force' for every node
    IO::GMSH::SurfaceVectorFieldDofBasedToGmsh(cutter_dis_,itrueresidual_,currinterfacepositions,gmshfilecontent,3,3);
    gmshfilecontent << "};" << std::endl;
  }

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "idispnp \" {" << std::endl;
    // draw vector field 'idispnp' for every node
    IO::GMSH::SurfaceVectorFieldDofBasedToGmsh(cutter_dis_,idispnp_,currinterfacepositions,gmshfilecontent,3,3);
    gmshfilecontent << "};" << std::endl;
  }

  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "ivelnp \" {" << std::endl;
    // draw vector field 'ivelnp' for every node
    IO::GMSH::SurfaceVectorFieldDofBasedToGmsh(cutter_dis_,ivelnp_,currinterfacepositions,gmshfilecontent,3,3);
    gmshfilecontent << "};" << std::endl;
  }

  gmshfilecontent.close();
}

void XFEM::MeshCouplingFSI::PrepareCutterOutput()
{
  // -------------------------------------------------------------------
  // prepare output
  // -------------------------------------------------------------------

  cutter_dis_->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(cutter_dis_)));
  cutter_output_ = cutter_dis_->Writer();
  cutter_output_->WriteMesh(0,0.0);
}

void XFEM::MeshCouplingFSI::Output(
    const int step,
    const double time,
    const bool write_restart_data
)
{
  // output for interface
  cutter_output_->NewStep(step,time);

  cutter_output_->WriteVector("ivelnp", ivelnp_);
  cutter_output_->WriteVector("idispnp", idispnp_);
  cutter_output_->WriteVector("itrueresnp", itrueresidual_);

  cutter_output_->WriteElementData(firstoutputofrun_);
  firstoutputofrun_ = false;

  // write restart
  if (write_restart_data)
  {
    cutter_output_->WriteVector("iveln_res",   iveln_);
    cutter_output_->WriteVector("idispn_res",  idispn_);
    cutter_output_->WriteVector("ivelnp_res",  ivelnp_);
    cutter_output_->WriteVector("idispnp_res", idispnp_);
  }
}


/*--------------------------------------------------------------------------*
 | extract the nodal vectors and store them in node-vector-map schott 01/13 |
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::ExtractNodeVectors(
    Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<Epetra_Vector> dofrowvec,
    std::map<int, LINALG::Matrix<3,1> >& nodevecmap
)
{
  Teuchos::RCP<const Epetra_Vector> dispcol = DRT::UTILS::GetColVersionOfRowVector(dis, dofrowvec);

  nodevecmap.clear();

  for (int lid = 0; lid < dis->NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = dis->lColNode(lid);
    std::vector<int> lm;
    dis->Dof(node, lm);
    std::vector<double> mydisp;
    DRT::UTILS::ExtractMyValues(*dispcol,mydisp,lm);
    if (mydisp.size() < 3)
      dserror("we need at least 3 dofs here");

    LINALG::Matrix<3,1> currpos;
    currpos(0) = node->X()[0] + mydisp[0];
    currpos(1) = node->X()[1] + mydisp[1];
    currpos(2) = node->X()[2] + mydisp[2];
    nodevecmap.insert(std::make_pair(node->Id(),currpos));
  }
}


//----------------------------------------------------------------------
// LiftDrag                                                  chfoe 11/07
//----------------------------------------------------------------------
//calculate lift&drag forces
//
//Lift and drag forces are based upon the right hand side true-residual entities
//of the corresponding nodes. The contribution of the end node of a line is entirely
//added to a present L&D force.
/*----------------------------------------------------------------------*/
void XFEM::MeshCouplingFSI::LiftDrag(
    const int step,
    const double time
) const
{
  // get forces on all procs
  // create interface DOF vectors using the fluid parallel distribution
  Teuchos::RCP<const Epetra_Vector> iforcecol = DRT::UTILS::GetColVersionOfRowVector(cutter_dis_, itrueresidual_);

  if (cutter_dis_->Comm().MyPID() == 0)
  {
    // compute force components
    const int nsd = 3;
    const Epetra_Map* dofcolmap = cutter_dis_->DofColMap();
    LINALG::Matrix<3,1> c(true);
    for (int inode = 0; inode < cutter_dis_->NumMyColNodes(); ++inode)
    {
      const DRT::Node* node = cutter_dis_->lColNode(inode);
      const std::vector<int> dof = cutter_dis_->Dof(node);
      for (int isd = 0; isd < nsd; ++isd)
      {
        // [// minus to get correct sign of lift and drag (force acting on the body) ]
        c(isd) += (*iforcecol)[dofcolmap->LID(dof[isd])];
      }
    }

    // print to file
    std::ostringstream s;
    std::ostringstream header;

    header << std::left  << std::setw(10) << "Time"
        << std::right << std::setw(16) << "F_x"
        << std::right << std::setw(16) << "F_y"
        << std::right << std::setw(16) << "F_z";
    s << std::left  << std::setw(10) << std::scientific << time
        << std::right << std::setw(16) << std::scientific << c(0)
        << std::right << std::setw(16) << std::scientific << c(1)
        << std::right << std::setw(16) << std::scientific << c(2);

    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                                        + ".liftdrag."
                                        + cond_name_
                                        + ".txt";
    if (step <= 1)
    {
      f.open(fname.c_str(),std::fstream::trunc);
      f << header.str() << std::endl;
    }
    else
    {
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
    }
    f << s.str() << "\n";
    f.close();

    std::cout << header.str() << std::endl << s.str() << std::endl;
  }
}



//! constructor
XFEM::MeshCouplingFSICrack::MeshCouplingFSICrack(
    Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
    const std::string &                 cond_name,///< name of the condition, by which the derived cutter discretization is identified
    Teuchos::RCP<DRT::Discretization>&  cond_dis  ///< discretization from which cutter discretization can be derived
) : MeshCouplingFSI(bg_dis,cond_name,cond_dis)
{
  InitCrackInitiationsPoints();

  {
    // @ Sudhakar
    // keep a pointer to the original boundary discretization
    // note: for crack problems, the discretization is replaced by new ones during the simulation.
    // Paraview output based on changing discretizations is not possible so far.
    // To enable at least restarts, the IO::DiscretizationWriter(boundarydis_) has to be kept alive,
    // however, in case that the initial boundary dis, used for creating the Writer, is replaced, it will be deleted,
    // as now other RCP points to it anymore. Then the functionality of the Writer breaks down. Therefore, we artificially
    // hold second pointer to the original boundary dis for Crack-problems.
    cutterdis_init_output_ = cutter_dis_;
  }
}


void XFEM::MeshCouplingFSICrack::SetCutterDis(Teuchos::RCP<DRT::Discretization> cutter_dis_new)
{
  cutter_dis_ = cutter_dis_new;

  // update the Coupling object

  // set unique element conditions
  SetElementConditions();

  // set the coupling strategy
  SetCouplingStrategy();

  // set coupling discretization
  SetCouplingDiscretization();
}


void XFEM::MeshCouplingFSICrack::InitCrackInitiationsPoints()
{
  tip_nodes_.clear();

  DRT::Condition* crackpts = cond_dis_->GetCondition( "CrackInitiationPoints" );

  const std::vector<int>* crackpt_nodes = const_cast<std::vector<int>* >(crackpts->Nodes());


  for(std::vector<int>::const_iterator inod=crackpt_nodes->begin(); inod!=crackpt_nodes->end();inod++ )
  {
    const int nodid = *inod;
    LINALG::Matrix<3, 1> xnod( true );

    tip_nodes_[nodid] = xnod;
  }

  if( tip_nodes_.size() == 0 )
    dserror("crack initiation points unspecified\n");

/*---------------------- POSSIBILITY 2 --- adding crack tip elements ----------------------------*/
/*
{
DRT::Condition* crackpts = soliddis_->GetCondition( "CrackInitiationPoints" );

const std::vector<int>* tipnodes = const_cast<std::vector<int>* >(crackpts->Nodes());

if( tipnodes->size() == 0 )
  dserror("crack initiation points unspecified\n");

addCrackTipElements( tipnodes );
}*/

/*  Teuchos::RCP<DRT::DofSet> newdofset1 = Teuchos::rcp(new DRT::TransparentIndependentDofSet(soliddis_,true));

boundarydis_->ReplaceDofSet(newdofset1);//do not call this with true!!
boundarydis_->FillComplete();*/
}


void XFEM::MeshCouplingFSICrack::UpdateBoundaryValuesAfterCrack(
    const std::map<int,int>& oldnewIds
)
{
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( cutter_dis_, ivelnp_, oldnewIds );
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( cutter_dis_, iveln_,  oldnewIds );
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( cutter_dis_, ivelnm_, oldnewIds );

  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( cutter_dis_, idispnp_, oldnewIds );
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( cutter_dis_, idispn_,  oldnewIds );

  //itrueresidual_ = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);
  DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( cutter_dis_, itrueresidual_, oldnewIds );

  //TODO: I guess the following lines are unnecessary (Sudhakar)
  {
    //iforcenp_ = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);
    //LINALG::Export( *itrueresidual_, *iforcenp_ );

  }

  //TODO: Check whether the output in case of crack-FSI work properly (Sudhakar)
  //boundarydis_->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(boundarydis_)));
  //boundary_output_ = boundarydis_->Writer();
}


//void XFEM::MeshCouplingFSICrack::PrepareCutterOutputCrack()
//{
//  // -------------------------------------------------------------------
//  // prepare output
//  // -------------------------------------------------------------------
//
//  cutter_dis_->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(cutter_dis_)));
//
//
//
//  cutter_output_ = cutter_dis_->Writer();
//  cutter_output_->WriteMesh(0,0.0);
//}



XFEM::LevelSetCoupling::LevelSetCoupling(
    Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
    const std::string &                 cond_name ///< name of the condition, by which the derived cutter discretization is identified
) : CouplingBase(bg_dis, cond_name, bg_dis)
{
  /// level-set field is given w.r.t background mesh
  /// NOTE: more generally it would be possible cutterdis != bg_dis for the single LevelSetCoupling,
  /// however, the unique bg_phinp vector stored in the ConditionManager has to be given w.r.t bgdis
  cutter_dis_ = bg_dis;

  SetConditionsToCopy();

  SetElementConditions();

  // set the coupling strategy
  SetCouplingStrategy();

  // set coupling discretization
  SetCouplingDiscretization();

  // create node-based vector with level-set values
  phinp_ = Teuchos::rcp(new Epetra_Vector(*cutter_dis_->NodeRowMap()));

}


void XFEM::LevelSetCoupling::SetConditionsToCopy()
{
  // set only the unique given condition name
  conditions_to_copy_.push_back(cond_name_);
}


/*----------------------------------------------------------------------*
 | ... |
 *----------------------------------------------------------------------*/
void XFEM::LevelSetCoupling::SetLevelSetField(
   Teuchos::RCP<const Epetra_Vector> scalaraf,
   Teuchos::RCP<DRT::Discretization> scatradis
   )
{
  // initializations
  int err(0);
  double value(0.0);
  std::vector<int> nodedofs;


  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<cutter_dis_->NumMyRowNodes();lnodeid++)
  {
    // get the processor's local scatra node
    DRT::Node* lscatranode = scatradis->lRowNode(lnodeid);

    // find out the global dof id of the last(!) dof at the scatra node
    const int numscatradof = scatradis->NumDof(0,lscatranode);
    const int globalscatradofid = scatradis->Dof(0,lscatranode,numscatradof-1);
    const int localscatradofid = scalaraf->Map().LID(globalscatradofid);
    if (localscatradofid < 0)
      dserror("localdofid not found in map for given globaldofid");

    // now copy the values
    value = (*scalaraf)[localscatradofid];
    err = phinp_->ReplaceMyValue(lnodeid,0,value);
    if (err != 0) dserror("error while inserting value into phinp_");

  }

  return;
}


/*----------------------------------------------------------------------*
 | ... |
 *----------------------------------------------------------------------*/
void XFEM::LevelSetCoupling::SetLevelSetField(const double time)
{
  // initializations
  int err(0);
  double value(0.0);
  std::vector<int> nodedofs;

  // get the function from the first element
  const int lid=0;
  DRT::Condition* cond = cutterele_conds_[lid].second;
  const int func_no = cond->GetInt("levelsetfieldno");

  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<cutter_dis_->NumMyRowNodes();lnodeid++)
  {
    // get the processor's local scatra node
    DRT::Node* lnode = cutter_dis_->lRowNode(lnodeid);

    // get value
    value=DRT::Problem::Instance()->Funct(func_no-1).Evaluate(0,lnode->X(),time,NULL);

    // now copy the values
    err = phinp_->ReplaceMyValue(lnodeid,0,value);
    if (err != 0) dserror("error while inserting value into phinp_");
  }

  return;
}



//constructor
XFEM::ConditionManager::ConditionManager(
      Teuchos::RCP<DRT::Discretization> &               bg_dis,          ///< background discretization
      std::vector<Teuchos::RCP<DRT::Discretization> > & meshcoupl_dis    ///< mesh coupling discretizations
  ) :
  bg_dis_(bg_dis),
  levelset_gid_(-1),
  time_(0.0),
  is_levelset_uptodate_(false),
  bg_phinp_(Teuchos::null)
{

  // create Levelset Coupling objects
  {
    std::vector<std::string> conditions_to_check;

    conditions_to_check.push_back("XFEMLevelsetWeakDirichlet");
    conditions_to_check.push_back("XFEMLevelsetNeumann");
    conditions_to_check.push_back("XFEMLevelsetTwophase");
    conditions_to_check.push_back("XFEMLevelsetCombustion");

    std::vector<std::string> names;
    bg_dis_->GetConditionNames( names );

    // check if background discretization has relevant conditioned nodes
    // create new coupling object for each type of condition
    for(size_t c=0; c<conditions_to_check.size(); c++)
    {
      if(std::find(names.begin(), names.end(), conditions_to_check[c]) != names.end())
        CreateNewLevelSetCoupling(conditions_to_check[c]);

    }
  }

  // create Mesh Coupling objects
  {
    std::vector<std::string> conditions_to_check;
    conditions_to_check.push_back("XFEMSurfFSIPart");
    conditions_to_check.push_back("XFEMSurfFSIMono");
    conditions_to_check.push_back("XFEMSurfCrackFSIPart");
    conditions_to_check.push_back("XFEMSurfFluidFluid");
    conditions_to_check.push_back("XFEMSurfWeakDirichlet");
    conditions_to_check.push_back("XFEMSurfNeumann");


    // check if a coupling discretization has relevant conditioned nodes
    // create new coupling object for each type of condition and each coupling discretization
    for(size_t mc_idx=0; mc_idx<meshcoupl_dis.size(); mc_idx++) // loop all specified mesh coupling discretizations
    {
      if(meshcoupl_dis[mc_idx] == Teuchos::null) continue;

      std::vector<std::string> names;
      meshcoupl_dis[mc_idx]->GetConditionNames( names );

      for(size_t c=0; c<conditions_to_check.size(); c++)
      {
        if(std::find(names.begin(), names.end(), conditions_to_check[c]) != names.end())
          CreateNewMeshCoupling(conditions_to_check[c], meshcoupl_dis[mc_idx]);
      }
    }
  }

}

void XFEM::ConditionManager::CreateNewLevelSetCoupling(const std::string& cond_name)
{
  AddLevelSetCoupling( cond_name );
}


void XFEM::ConditionManager::CreateNewMeshCoupling(
    const std::string& cond_name,
    Teuchos::RCP<DRT::Discretization> cond_dis         ///< discretization from which the cutter discretization can be derived
)
{
  AddMeshCoupling( cond_name, cond_dis );
}


void XFEM::ConditionManager::Create(const double time)
{
  time_ = time;
  numglobal_coupling_sides = 0;
  mesh_coupl_start_gid_.reserve(mesh_coupl_.size());
  levelset_gid_ = -1;

  // set global side Ids for all Mesh coupling discretizations and level-set sides

  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    Teuchos::RCP<DRT::Discretization> mc_cutdis = mesh_coupl_[mc]->GetCutterDis();
    if(mc_cutdis == Teuchos::null) dserror("cutter dis is Teuchos::null");

    // set current number of global coupling sides as start index for global id this coupling object
    mesh_coupl_start_gid_[mc] = numglobal_coupling_sides;

    std::cout << "mesh coupling object " << mc << " starts with global side index " << mesh_coupl_start_gid_[mc] << std::endl;

    // increase total number of sides with number of global side elements of this mesh coupling object
    numglobal_coupling_sides += mc_cutdis->NumGlobalElements();
  }

  // TODO: unify the level-set coupling objects to one unique level-set field
  // combine the level-set values

  if(levelset_coupl_.size() > 0)
  {

    // add one global levelset side used in the cut library
    levelset_gid_ = numglobal_coupling_sides;
    numglobal_coupling_sides+=1;

    std::cout << "levelset coupling object " << " has global side id index " << levelset_gid_ << std::endl;

    bg_phinp_= LINALG::CreateVector(*bg_dis_->NodeRowMap(), true);

    //TODO: note: information about the coupling condition for level-sets is obtained via the background element
    // for which we store the index of the level-set coupling object
    // we allow for multiple level-set coupling objects however only for one level-set side
  }



  //TODO: Status-routine to print to screen

}


void XFEM::ConditionManager::SetLevelSetField( const double time )
{
  if(levelset_coupl_.size() != 1)
    dserror("level-set field is not unique, which level-set field to be set?");

  is_levelset_uptodate_ = false;

  // update the unique level-set field
  levelset_coupl_[0]->SetLevelSetField(time);
}


void XFEM::ConditionManager::SetLevelSetField(
   Teuchos::RCP<const Epetra_Vector> scalaraf,
   Teuchos::RCP<DRT::Discretization> scatradis
   )
{
  if(levelset_coupl_.size() != 1)
    dserror("level-set field is not unique, which level-set field to be set by given vector?");

  is_levelset_uptodate_ = false;

  // update the unique level-set field
  levelset_coupl_[0]->SetLevelSetField(scalaraf, scatradis);
}


Teuchos::RCP<const Epetra_Vector> XFEM::ConditionManager::GetLevelSetFieldCol()
{
  if(levelset_coupl_.size()==0) return Teuchos::null;

  // export nodal level-set values to node column map
  Teuchos::RCP<Epetra_Vector> bg_phinp_col = Teuchos::rcp(new Epetra_Vector(*bg_dis_->NodeColMap()));
  LINALG::Export(*GetLevelSetField(),*bg_phinp_col);

  return bg_phinp_col;
}

void XFEM::ConditionManager::UpdateLevelSetField()
{

  if(levelset_coupl_.size() != 1)
      dserror("level-set field is not unique, update of the global bg_phinp by more than one phinp not implemented yet");

  // assume same maps between background fluid dis and the cutterdis (scatra dis)
  bg_phinp_->Update(1.0, *(levelset_coupl_[0]->GetLevelSetField()), 0.0);

  is_levelset_uptodate_ = true;
}

void XFEM::ConditionManager::SetState()
{
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->SetState();
  }
}

void XFEM::ConditionManager::CompleteStateVectors()
{
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->CompleteStateVectors();
  }
}

void XFEM::ConditionManager::ZeroStateVectors_FSI()
{
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->ZeroStateVectors_FSI();
  }
}

void XFEM::ConditionManager::GmshOutput(
    const std::string & filename_base,
    const int step,
    const int gmsh_step_diff,
    const bool gmsh_debug_out_screen
)
{
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->GmshOutput(filename_base, step, gmsh_step_diff, gmsh_debug_out_screen);
  }
}

void XFEM::ConditionManager::Output(
    const int step,
    const double time,
    const bool write_restart_data
)
{
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->Output(step, time, write_restart_data);
  }
}

void XFEM::ConditionManager::LiftDrag(
    const int step,
    const double time
)
{
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->LiftDrag(step, time);
  }
}

void XFEM::ConditionManager::ReadRestart(
    const int step
)
{
  // loop all mesh coupling objects
  for(int mc=0; mc<(int)mesh_coupl_.size(); mc++)
  {
    mesh_coupl_[mc]->ReadRestart(step);
  }
}

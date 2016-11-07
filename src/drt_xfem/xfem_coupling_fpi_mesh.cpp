/*!
\file xfem_coupling_fpi_mesh.cpp

\brief manages mesh based coupling of fluid and porous media
xfluid class and the cut-library

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>
*/

#include <Teuchos_TimeMonitor.hpp>

#include "xfem_coupling_fpi_mesh.H"

#include "xfem_utils.H"
#include "xfem_discretization_utils.H"

#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_dofset_transparent_independent.H"

#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_ele/fluid_ele_parameter_xfem.H"

#include "../linalg/linalg_utils.H"

#include "../drt_io/io.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

//! constructor
XFEM::MeshCouplingFPI::MeshCouplingFPI(
    Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
    const std::string &                 cond_name,///< name of the condition, by which the derived cutter discretization is identified
    Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< discretization from which cutter discretization can be derived
    const int                           coupling_id,///< id of composite of coupling conditions
    const double                        time,      ///< time
    const int                           step,       ///< time step
    MeshCouplingFPI::coupled_field     field      ///< which field is coupled to the fluid
    ) : MeshCoupling(bg_dis,cond_name,cond_dis,coupling_id, time, step,
        (field==MeshCouplingFPI::ps_ps ? "_ps_ps" : (field==MeshCouplingFPI::ps_pf ? "_ps_pf" : (field==MeshCouplingFPI::pf_ps ? "_pf_ps" : "_pf_pf")))),
        coupled_field_(field)
{
  //Todo Calculation of the Porosity in the gausspoint will folllow here soon ...
  std::cout << RED << "###########################################################################" << END_COLOR << std::endl;
  std::cout << RED << "###########################################################################" << END_COLOR << std::endl;
  std::cout << RED << "## WARNING: At the moment we suppose to have a constant porosity of 0.5! ##" << END_COLOR << std::endl;
  std::cout << RED << "###########################################################################" << END_COLOR << std::endl;
  std::cout << RED << "###########################################################################" << END_COLOR << std::endl;
  //Set Condname dependent on the field type ... to be able to differentiate between the two fpi coupling objects!
  std::stringstream str;
  str << cond_name << (field==MeshCouplingFPI::ps_ps ? "_ps_ps" : (field==MeshCouplingFPI::ps_pf ? "_ps_pf" : (field==MeshCouplingFPI::pf_ps ? "_pf_ps" : "_pf_pf")));

  SetConditionSpecificParameters(cond_name);

  cond_name_ = str.str();

  full_BJ_ = true; //Todo Ager: Will create an Input parameter in the *.dat-file ...

  InitStateVectors_FPI();

  // create map from side to embedded element ID
  CreateCuttingToEmbeddedElementMap();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFPI::InitStateVectors_FPI()
{
  const Epetra_Map* cutterdofrowmap = cutter_dis_->DofRowMap();
  const Epetra_Map* cutterdofcolmap = cutter_dis_->DofColMap();

  itrueresidual_ = LINALG::CreateVector(*cutterdofrowmap,true);
  iforcecol_     = LINALG::CreateVector(*cutterdofcolmap,true);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFPI::CompleteStateVectors()
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

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFPI::InitConfigurationMap()
{

  //for the first tests not done ...
  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFPI::UpdateConfigurationMap_GP(
    double& kappa_m,
    double& visc_m,
    double& visc_s,
    double& visc_stab,
    double& full_stab,
    const LINALG::Matrix<3,1>& x,
    const DRT::Condition* cond)
{
  double dynvisc   = (kappa_m*visc_m + (1.0-kappa_m)*visc_s);

  double stabnit = 0.0;
  double stabadj = 0.0;
  XFEM::UTILS::GetNavierSlipStabilizationParameters(full_stab,visc_stab,dynvisc,BJ_coeff_,stabnit,stabadj);

  // Calculation on Porosity will follow soon ... todo
  double porosity;
  //const Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<DISTYPE> > & si,              ///< associated slave element coupling object
  //int coup_sid                                                             ///< id of the coupling side
  // porosity = CalcPorosity<DISTYPE>(si,coup_sid);
  // if (fabs(porosity - 0.5) > 1e-14)
  //    std::cout << "porosity: " << porosity << std::endl;
  porosity = 0.5;

  //Overall there are 9 coupling blocks to evaluate for fpi:
  // 1 - ps_ps --> ff,fps,fps,psf
  // 2 - ps_pf --> fpf,pspf
  // 3 - pf_ps --> pff, pfps
  // 4 - pf_pf --> pfpf
switch(coupled_field_)
{
  case MeshCouplingFPI::ps_ps:
  {
    //Configuration of Consistency Terms
    configuration_map_[INPAR::XFEM::F_Con_Row] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::F_Con_Col] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::X_Con_Row] = std::pair<bool,double>(true,1.0);
    //Configuration of Adjount Consistency Terms
    configuration_map_[INPAR::XFEM::F_Adj_n_Row] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::F_Adj_n_Col] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::X_Adj_n_Col] = std::pair<bool,double>(true,1-porosity);
    configuration_map_[INPAR::XFEM::F_Adj_t_Row] = std::pair<bool,double>(true,stabadj);
    configuration_map_[INPAR::XFEM::F_Adj_t_Col] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::X_Adj_t_Col] = std::pair<bool,double>(true,1-full_BJ_*porosity);
    configuration_map_[INPAR::XFEM::FStr_Adj_t_Col] = std::pair<bool,double>(true, BJ_coeff_); //Here we need alpha BJ finally!!! Todo
    //Configuration of Penalty Terms
    configuration_map_[INPAR::XFEM::F_Pen_n_Row] = std::pair<bool,double>(true,visc_stab);
    configuration_map_[INPAR::XFEM::F_Pen_n_Col] = std::pair<bool,double>(true,1);
    configuration_map_[INPAR::XFEM::X_Pen_n_Row] = std::pair<bool,double>(true,visc_stab);
    configuration_map_[INPAR::XFEM::X_Pen_n_Col] = std::pair<bool,double>(true,1-porosity);
    configuration_map_[INPAR::XFEM::F_Pen_t_Row] = std::pair<bool,double>(true,stabnit);
    configuration_map_[INPAR::XFEM::F_Pen_t_Col] = std::pair<bool,double>(true,1);
    configuration_map_[INPAR::XFEM::X_Pen_t_Row] = std::pair<bool,double>(true,stabnit);
    configuration_map_[INPAR::XFEM::X_Pen_t_Col] = std::pair<bool,double>(true,1-full_BJ_*porosity);
    configuration_map_[INPAR::XFEM::F_Con_t_Row] = std::pair<bool,double>(true,-stabnit); //+sign for penalty!
    configuration_map_[INPAR::XFEM::F_Con_t_Col] = std::pair<bool,double>(true,BJ_coeff_/dynvisc);
    configuration_map_[INPAR::XFEM::X_Con_t_Row] = std::pair<bool,double>(true,-stabnit); //+sign for penalty!
    break;
  }
  case MeshCouplingFPI::ps_pf:
  {
    //Configuration of Adjount Consistency Terms
    configuration_map_[INPAR::XFEM::F_Adj_n_Row] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::X_Adj_n_Col] = std::pair<bool,double>(true,porosity);
    configuration_map_[INPAR::XFEM::F_Adj_t_Row] = std::pair<bool,double>(true,stabadj);
    configuration_map_[INPAR::XFEM::X_Adj_t_Col] = std::pair<bool,double>(full_BJ_,full_BJ_*porosity);
//    //Configuration of Penalty Terms
    configuration_map_[INPAR::XFEM::F_Pen_n_Row] = std::pair<bool,double>(true,visc_stab);
    configuration_map_[INPAR::XFEM::X_Pen_n_Row] = std::pair<bool,double>(true,visc_stab);
    configuration_map_[INPAR::XFEM::X_Pen_n_Col] = std::pair<bool,double>(true,porosity);
    configuration_map_[INPAR::XFEM::F_Pen_t_Row] = std::pair<bool,double>(true,stabnit);
    configuration_map_[INPAR::XFEM::X_Pen_t_Row] = std::pair<bool,double>(true,stabnit);
    configuration_map_[INPAR::XFEM::X_Pen_t_Col] = std::pair<bool,double>(full_BJ_,full_BJ_*porosity);
    break;
  }
  case MeshCouplingFPI::pf_ps:
  {
    //Configuration of Consistency Terms
   configuration_map_[INPAR::XFEM::X_Con_n_Row] = std::pair<bool,double>(true,1.0);
   configuration_map_[INPAR::XFEM::F_Con_n_Col] = std::pair<bool,double>(true,1.0);
   //Configuration of Penalty Terms
    configuration_map_[INPAR::XFEM::X_Pen_n_Row] = std::pair<bool,double>(true,visc_stab);
    configuration_map_[INPAR::XFEM::F_Pen_n_Col] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::X_Pen_n_Col] = std::pair<bool,double>(true,1-porosity);

    //does nothing but should just be done in case we don't use the adjoint
    configuration_map_[INPAR::XFEM::F_Adj_n_Col].second = configuration_map_[INPAR::XFEM::F_Pen_n_Col].second;
    configuration_map_[INPAR::XFEM::X_Adj_n_Col].second = configuration_map_[INPAR::XFEM::X_Pen_n_Col].second;
    break;
  }
  case MeshCouplingFPI::pf_pf:
  {
    //Configuration of Penalty Terms
     configuration_map_[INPAR::XFEM::X_Pen_n_Row] = std::pair<bool,double>(true,visc_stab);
     configuration_map_[INPAR::XFEM::X_Pen_n_Col] = std::pair<bool,double>(true,porosity);

     //does nothing but should just be done in case we don't use the adjoint
     configuration_map_[INPAR::XFEM::X_Adj_n_Col].second = configuration_map_[INPAR::XFEM::X_Pen_n_Col].second;
    break;
  }
}

  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFPI::ZeroStateVectors_FPI()
{
  itrueresidual_->PutScalar(0.0);
  iforcecol_->PutScalar(0.0);
}

// -------------------------------------------------------------------
// Read Restart data for cutter discretization
// -------------------------------------------------------------------
void XFEM::MeshCouplingFPI::ReadRestart(
    const int step
)
{
  if(myrank_) IO::cout << "ReadRestart for boundary discretization " << IO::endl;

  //-------- boundary discretization
  IO::DiscretizationReader boundaryreader(cutter_dis_, step);

  const double time = boundaryreader.ReadDouble("time");
//  const int    step = boundaryreader.ReadInt("step");

  if(myrank_ == 0)
  {
    IO::cout << "time: " << time << IO::endl;
    IO::cout << "step: " << step << IO::endl;
  }

  boundaryreader.ReadVector(iveln_,   "iveln_res");
  boundaryreader.ReadVector(idispn_,  "idispn_res");

  // REMARK: ivelnp_ and idispnp_ are set again for the new time step in PrepareSolve()
  boundaryreader.ReadVector(ivelnp_,  "ivelnp_res");
  boundaryreader.ReadVector(idispnp_, "idispnp_res");
  boundaryreader.ReadVector(idispnpi_, "idispnpi_res");

  if (not (cutter_dis_->DofRowMap())->SameAs(ivelnp_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (cutter_dis_->DofRowMap())->SameAs(iveln_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (cutter_dis_->DofRowMap())->SameAs(idispnp_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (cutter_dis_->DofRowMap())->SameAs(idispn_->Map()))
    dserror("Global dof numbering in maps does not match");
  if (not (cutter_dis_->DofRowMap())->SameAs(idispnpi_->Map()))
    dserror("Global dof numbering in maps does not match");


}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFPI::GmshOutput(
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
  XFEM::UTILS::ExtractNodeVectors(cutter_dis_, currinterfacepositions,idispnp_);


  const std::string filename =
      IO::GMSH::GetNewFileNameAndDeleteOldFiles(
          filename_base_fsi.str(),
          step,
          gmsh_step_diff,
          gmsh_debug_out_screen,
          myrank_
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

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFPI::GmshOutputDiscretization(
  std::ostream& gmshfilecontent
)
{
  // print surface discretization
  XFEM::MeshCoupling::GmshOutputDiscretization(gmshfilecontent);

  // compute the current solid and boundary position
  std::map<int,LINALG::Matrix<3,1> > currsolidpositions;

  // write dis with zero solid displacements here!
  Teuchos::RCP<Epetra_Vector> solid_dispnp = LINALG::CreateVector(*cond_dis_->DofRowMap(), true);

  XFEM::UTILS::ExtractNodeVectors(cond_dis_, currsolidpositions, solid_dispnp);

  XFEM::UTILS::PrintDiscretizationToStream(cond_dis_,
      cond_dis_->Name(), true, false, true, false, false, false, gmshfilecontent, &currsolidpositions);
}

void XFEM::MeshCouplingFPI::Output(
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
    cutter_output_->WriteVector("idispnpi_res", idispnpi_);
  }
}

void XFEM::MeshCouplingFPI::SetConditionSpecificParameters(  const std::string &  cond_name)
{

  std::vector< DRT::Condition* >  conditions_XFPI;
  cutter_dis_->GetCondition(cond_name, conditions_XFPI);

  // Create maps for easy extraction at gausspoint level
  for(std::vector<DRT::Condition*>::iterator i=conditions_XFPI.begin();
      i!=conditions_XFPI.end();
      ++i)
  {
    (*i)->Print(std::cout);

    if (i!=conditions_XFPI.begin())
      if (fabs(BJ_coeff_-(*i)->GetDouble("bj_coeff")) > 1e-16)
        dserror("XFEM::MeshCouplingFPI::SetConditionSpecificParameters: You defined two FPI conditions, with different BJ_coeff!");

    BJ_coeff_ = (*i)->GetDouble("bj_coeff");

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
void XFEM::MeshCouplingFPI::LiftDrag(
    const int step,
    const double time
) const
{
  // get forces on all procs
  // create interface DOF vectors using the fluid parallel distribution
  Teuchos::RCP<const Epetra_Vector> iforcecol = DRT::UTILS::GetColVersionOfRowVector(cutter_dis_, itrueresidual_);

  if (myrank_ == 0)
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

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::MeshCouplingFPI::CreateCuttingToEmbeddedElementMap()
{
  // fill map between boundary (cutting) element id and its corresponding embedded (coupling) element id
  for (int ibele=0; ibele< cutter_dis_->NumMyColElements(); ++ ibele)
  {
    // boundary element and its nodes
    DRT::Element* bele = cutter_dis_->lColElement(ibele);
    const int * bele_node_ids = bele->NodeIds();

    bool bele_found = false;

    // ask all conditioned embedded elements for this boundary element
    for(int iele = 0; iele< cond_dis_->NumMyColElements(); ++iele)
    {
      DRT::Element* ele = cond_dis_->lColElement(iele);
      const int * ele_node_ids = ele->NodeIds();

      // get nodes for every face of the embedded element
      std::vector<std::vector<int> > face_node_map = DRT::UTILS::getEleNodeNumberingFaces(ele->Shape());

      // loop the faces of the element and check node equality for every boundary element
      // Todo: Efficiency?
      for (int f = 0; f < ele->NumFace(); f++)
      {
        bele_found = true;

        const int face_numnode = face_node_map[f].size();

        if(bele->NumNode() != face_numnode) continue; // this face cannot be the right one

        // check all nodes of the boundary element
        for(int inode=0; inode<bele->NumNode();  ++inode)
        {
          // boundary node
          const int belenodeId = bele_node_ids[inode];

          bool node_found = false;
          for (int fnode=0; fnode<face_numnode; ++fnode)
          {
            const int facenodeId = ele_node_ids[face_node_map[f][fnode]];

            if(facenodeId == belenodeId)
            {
              // nodes are the same
              node_found = true;
              break;
            }
          } // loop nodes of element's face
          if (node_found==false) // this node is not contained in this face
          {
            bele_found = false; // element not the right one, if at least one boundary node is not found
            break; // node not found
          }
        } // loop nodes of boundary element


        if (bele_found)
        {
          cutting_emb_gid_map_.insert(std::pair<int,int>(bele->Id(),ele->Id()));
          cutting_emb_face_lid_map_.insert(std::pair<int,int>(bele->Id(),f));
          break;
        }
      } // loop element faces
      if (bele_found) break; // do not continue the search

    }

    if(bele_found == false)
      dserror("Corresponding embedded element for boundary element id %i not found on proc %i ! Please ghost corresponding embedded elements on all procs!",
        bele->Id(), cond_dis_->Comm().MyPID());
  }
}

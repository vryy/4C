/*----------------------------------------------------------------------*/
/*!
\file xfem_coupling_levelset.cpp

\brief manages the different types of level-set based coupling conditions and thereby builds the bridge between the
xfluid class and the cut-library

\level 2

<pre>
\maintainer Benedikt Schott & Magnus Winter
            {schott, winter}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>

#include "xfem_coupling_levelset.H"
#include "xfem_utils.H"
#include "xfem_interface_utils.H"

#include "../drt_cut/cut_cutwizard.H"
#include "../drt_cut/cut_node.H"
#include "../drt_cut/cut_point.H"


#include "../linalg/linalg_utils.H"

#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_fluid.H"

#include "../drt_io/io.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_lib/drt_discret_xfem.H"

//Needed to find element conditions
#include "../drt_lib/drt_condition_utils.H"

#include "../drt_mat/newtonianfluid.H"

//TODO: CouplingBase should become abstract class


XFEM::LevelSetCoupling::LevelSetCoupling(
    Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
    const std::string &                 cond_name, ///< name of the condition, by which the derived cutter discretization is identified
    Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< full discretization from which the cutter discretization is derived
    const int                           coupling_id,///< id of composite of coupling conditions
    const double                        time,      ///< time
    const int                           step       ///< time step
) : CouplingBase(bg_dis, cond_name, cond_dis, coupling_id, time, step),
    bg_nds_phi_(-1),
    cutter_nds_phi_(-1),
    normal_orientation_(-1.0), // level set gradient is directed from inside to outside -- normal in xfluid points from outside to inside
    have_nodematching_dis_(false)
{
}

void XFEM::LevelSetCoupling::SetCouplingDofsets()
{
  bg_nds_phi_ = GetCouplingDofsetNds("phi_scatra_proxy_in_fluid");
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCoupling::SetConditionsToCopy()
{
  // set only the unique given condition name
  conditions_to_copy_.push_back(cond_name_);

  // additional conditions required for the levelset field based on the cutter (background) mesh
  conditions_to_copy_.push_back("XFEMSurfDisplacement");
}


//TODO: needs to be generalized!
/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCoupling::SetCutterDiscretization()
{
  if(&*cond_dis_ != &*bg_dis_) // check not required for two-phase anymore
    dserror("for non-twophase couplings, we have not checked the functionality of using an extra scatra this -- lets try!");

  /// level-set field is given w.r.t background mesh
  /// NOTE: more generally it would be possible cutterdis != bg_dis for the single LevelSetCoupling,
  /// however, the unique bg_phinp vector stored in the ConditionManager has to be given w.r.t bgdis

  // TODO: shall we allow to use subsets of the cond_dis_ as a cutterdis?
  cutter_dis_ = cond_dis_;

  // do we have node-matching disretizations? otherwise we need to somehow project quantities between the meshes
  have_nodematching_dis_ = HaveMatchingNodes(cutter_dis_, bg_dis_);


  std::string dofset_name = "";

  if(cutter_dis_->Name()=="scatra")
    dofset_name = "phi_in_scatra";
  else if(cutter_dis_->Name()=="fluid")
    dofset_name = "phi_scatra_proxy_in_fluid";
  else
    dserror("unsupported cutter dis!");

  if(not(dofset_coupling_map_.count(dofset_name) > 0))
      dserror("dofset not set in dofset_coupling_map for cutter dis!");

  cutter_nds_phi_ = dofset_coupling_map_[dofset_name]; // dofset id for scalar field
}

//TODO: shift to DRT::Utils...
/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
bool XFEM::LevelSetCoupling::HaveMatchingNodes(
    const Teuchos::RCP<DRT::Discretization> & dis_A,
    const Teuchos::RCP<DRT::Discretization> & dis_B
)
{
  // check for equal node row maps
  const Epetra_Map * noderowmap_A = dis_A->NodeRowMap();
  const Epetra_Map * noderowmap_B = dis_B->NodeRowMap();

  if(!(noderowmap_A->SameAs(*noderowmap_B)))
    return false;

  // check for equal node coordinates
  for(int lid=0; lid<noderowmap_A->NumMyElements(); ++lid)
  {
    const DRT::Node * node_A = dis_A->lRowNode(lid);
    const DRT::Node * node_B = dis_B->lRowNode(lid);

    const int nsd = node_A->Dim();

    Epetra_SerialDenseVector X_A(nsd);
    Epetra_SerialDenseVector X_B(nsd);

    std::copy(node_A->X(), node_A->X()+nsd, X_A.A());
    std::copy(node_B->X(), node_B->X()+nsd, X_B.A());

    Epetra_SerialDenseVector diff(X_A);
    diff.Scale(-1.0);
    diff+=X_B;

    if(diff.Norm2()> 1e-14)
      return false;
  }

  return true;
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCoupling::InitStateVectors()
{
  // initialize state vectors w.r.t. background discretization
  InitStateVectors_Bg();
  // initialize state vectors w.r.t. cutter (potential subset of scatra) discretization
  InitStateVectors_Cutter();
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCoupling::InitStateVectors_Bg()
{
  // background-dis (fluid) related state vectors
  const Epetra_Map* bg_dofrowmap = bg_dis_->DofRowMap(bg_nds_phi_);

  phinp_ = LINALG::CreateVector(*bg_dofrowmap, true);
}


//TODO: check if all vectors are really used and implement a save export strategy... (also used in mesh coupling?)
/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCoupling::InitStateVectors_Cutter()
{
  // cutter-dis related state vectors
  const Epetra_Map* cutter_dofrowmap  = cutter_dis_->DofRowMap(cutter_nds_phi_); // used for level set field and its derivatives
  const Epetra_Map* cutter_dofcolmap  = cutter_dis_->DofColMap(cutter_nds_phi_); // used for level set field and its derivatives

  cutter_phinp_                = LINALG::CreateVector(*cutter_dofrowmap, true);
  gradphinp_smoothed_node_     = LINALG::CreateMultiVector(*cutter_dofrowmap,nsd_, true);
  gradphinp_smoothed_node_col_ = LINALG::CreateMultiVector(*cutter_dofcolmap,nsd_, true);
  curvaturenp_node_            = LINALG::CreateVector(*cutter_dofrowmap, true);
  curvaturenp_node_col_        = LINALG::CreateVector(*cutter_dofcolmap, true);
}


//TODO: check output functionality
/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCoupling::PrepareCutterOutput()
{
  // -------------------------------------------------------------------
  // prepare output
  // -------------------------------------------------------------------
  cutter_output_ = cutter_dis_->Writer();

  if(cutter_output_ == Teuchos::null)
      cutter_dis_->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(cutter_dis_)));

  bg_output_ = bg_dis_->Writer();
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCoupling::DoConditionSpecificSetup()
{
  //TODO: remove the smoothed normal stuff from this function!!!
  // read initial level-set field
  SetLevelSetField(time_);

  // set level-boolean type (may be overwritten in constructors of derived class
  SetLevelSetBooleanType();
}

//TODO: store the first element condition to access cutterele_conds_[0] at several places
/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCoupling::SetLevelSetBooleanType()
{

  if(cutterele_conds_.size() == 0)
    dserror("no element condition for LevelSetCouplingBC set. Not possible to extract BOOLEANTYPE!");

  DRT::Condition* cond = (cutterele_conds_[0]).second;
  const std::string* booleantype = cond->Get<std::string>("booleantype");

  if(*booleantype == "none")
    ls_boolean_type_ = ls_none;
  else if(*booleantype == "cut")
    ls_boolean_type_ = ls_cut;
  else if(*booleantype == "union")
    ls_boolean_type_ = ls_union;
  else if(*booleantype == "difference")
    ls_boolean_type_ = ls_difference;
  else if(*booleantype == "sym_difference")
    ls_boolean_type_ = ls_sym_difference;
  else
    dserror("not a valid boolean type %s: ", booleantype->c_str());
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
bool XFEM::LevelSetCoupling::ApplyComplementaryOperator()
{
  if(cutterele_conds_.size() == 0)
    dserror("no element condition for LevelSetCouplingBC set. Not possible to extract BOOLEANTYPE!");

  DRT::Condition* cond = (cutterele_conds_[0]).second;
  bool complementary = (bool)cond->GetInt("complementary");

  return complementary;
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCoupling::Output(
    const int step,
    const double time,
    const bool write_restart_data,
    const int lsc_idx
)
{
  // output for level-set interface
  // bg_output_->NewStep(step,time); // not required, as already called for the bgdis when output is written for fluid fields

  std::ostringstream temp;
  temp << lsc_idx;
  std::string name = "phinp_"+temp.str();

  bg_output_->WriteVector(name, phinp_);

  // write restart
  if (write_restart_data)
  {
    std::string name_restart = "phinp_res_"+temp.str();

    bg_output_->WriteVector(name_restart, phinp_);
  }

  cutter_output_->NewStep(step,time); // required, as already called for the bgdis when output is written for fluid fields

  // write restart
  if (write_restart_data)
  {
    std::ostringstream temp2;
    temp2 << lsc_idx;
    std::string name_restart = "cutter_phinp_res_"+temp.str();

    cutter_output_->WriteVector(name_restart, cutter_phinp_);
  }

}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCoupling::GmshOutput(
    const std::string & filename_base,
    const int step,
    const int gmsh_step_diff,
    const bool gmsh_debug_out_screen
)
{
  //TODO: adapt!!!
  // gmsh output of geometric quantities based on the cutterdis (=fluid dis) (mapped quantities from scatradis to fluiddis)
  // cutterdis=fluiddis within the condition manager

  std::ostringstream filename_base_fsi;
  filename_base_fsi << filename_base << "_levelset";

//  // compute the current boundary position
//  std::map<int,LINALG::Matrix<3,1> > currinterfacepositions;
//  XFEM::UTILS::ExtractNodeVectors(cutter_dis_, currinterfacepositions,idispnp_);


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
    gmshfilecontent << "View \" " << "cutter_phinp_ \" {" << std::endl;
    // draw vector field 'force' for every node
    IO::GMSH::ScalarFieldNodeBasedToGmsh(cutter_dis_,cutter_phinp_,gmshfilecontent);
    gmshfilecontent << "};" << std::endl;
  }

  //TODO: add further output for other state vectors!

//  {
//    // add 'View' to Gmsh postprocessing file
//    gmshfilecontent << "View \" " << "smoothed gradphi \" {" << std::endl;
//    // draw vector field 'idispnp' for every node
////    IO::GMSH::VectorFieldNodeBasedToGmsh(cutter_dis_,gradphinp_smoothed_node_,gmshfilecontent);
//    IO::GMSH::VectorFieldDofBasedToGmsh(cutter_dis_,gradphinp_smoothed_node_,gmshfilecontent,3);
//    gmshfilecontent << "};" << std::endl;
//  }
//
//  {
//    // add 'View' to Gmsh postprocessing file
//    gmshfilecontent << "View \" " << "smoothed curvature \" {" << std::endl;
//    // draw vector field 'idispnp' for every node
//    IO::GMSH::ScalarFieldNodeBasedToGmsh(cutter_dis_,curvaturenp_node_,gmshfilecontent);
//    gmshfilecontent << "};" << std::endl;
//  }

//  {
//    // add 'View' to Gmsh postprocessing file
//    gmshfilecontent << "View \" " << "transport velocity \" {" << std::endl;
//    // draw vector field 'force' for every node
//    IO::GMSH::ScalarFieldNodeBasedToGmsh(cutter_dis_,cutter_phinp_,gmshfilecontent);
//    gmshfilecontent << "};" << std::endl;
//  }
//
//  gmshfilecontent.close();
//  ;

  // gmsh output for geometric quantities
  // gmsh output for transport velocity...
}

// -------------------------------------------------------------------
// Read Restart data for cutter discretization
// -------------------------------------------------------------------
void XFEM::LevelSetCoupling::ReadRestart(
    const int step,
    const int lsc_idx
)
{

//  dserror("Not tested Level Set restart from file. Should most likely work though if this dserror is removed.");

  //-------- boundary discretization
  IO::DiscretizationReader boundaryreader(cutter_dis_, step);

  const double time = boundaryreader.ReadDouble("time");

  if(myrank_ == 0)
  {
    IO::cout << "            RESTART IS PERFORMED FROM FUNCTION IN INPUT FILE!                  " << IO::endl;
    IO::cout << "ReadRestart for Level Set Cut in Xfluid (time="<< time <<" ; step="<< step <<")" << IO::endl;
  }

  SetLevelSetField(time);

}

//TODO: remove the Navier-Slip stuff
/*---------------------------------------------------------------------------*
 | Set the level set field and if smoothed gradients are needed create these |
 |                                                                           |
 *---------------------------------------------------------------------------*/
bool XFEM::LevelSetCoupling::SetLevelSetField(const double time)
{
  //TODO: clean this routine!!!

  // make a copy of last time step

  Teuchos::RCP<Epetra_Vector> delta_phi = LINALG::CreateVector(cutter_phinp_->Map(),true);
  delta_phi->Update(1.0, *cutter_phinp_, 0.0);

  // initializations
  int err(0);
  double value(0.0);
  std::vector<int> nodedofs;

  // get the function from the first element
  const int lid=0;
  DRT::Condition* cond = cutterele_conds_[lid].second;
  const int func_no = cond->GetInt("levelsetfieldno");

  // check for potential time curve
  const int curvenum  = cond->GetInt("levelsetcurve");

  // initialization of time-curve factor
  double curvefac = 0.0;

  // compute potential time curve or set time-curve factor to one
  if (curvenum >= 0)
  {
    // time factor (negative time indicating error)
    if (time >= 0.0)
      curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
    else dserror("Negative time in function evaluation: time = %f", time);
  }
  else curvefac = 1.0;

  // loop all nodes on the processor
  for(int lnodeid=0;lnodeid<cutter_dis_->NumMyRowNodes();lnodeid++)
  {
    // get the processor's local scatra node
    DRT::Node* lnode = cutter_dis_->lRowNode(lnodeid);

    // get value
    if(func_no < 0)
      value = FunctImplementation(func_no, lnode->X(),time);
    else if(func_no >= 1)
      value=DRT::Problem::Instance()->Funct(func_no-1).Evaluate(0,lnode->X(),time,NULL);
    else
      dserror("invalid function no. to set level-set field!");

    double final_val = curvefac*value;

    // now copy the values
    err = cutter_phinp_->ReplaceMyValue(lnodeid,0,final_val);
    if (err != 0) dserror("error while inserting value into cutter_phinp_");
  }


  //TODO: remove this part from this function!!!

  // Might make this available for other condition types!
  const INPAR::XFEM::EleCouplingCondType cond_type = CondType_stringToEnum(cond_name_);

  if(cond_type == INPAR::XFEM::CouplingCond_LEVELSET_NAVIER_SLIP)
  {

    //Do we need smoothed gradients? I.e. what type is it?
    const int val = cutterele_conds_[lid].second->GetInt("SURFACE_PROJECTION");
    projtosurf_ = static_cast<INPAR::XFEM::ProjToSurface>(val);

    if(projtosurf_!=INPAR::XFEM::Proj_normal) //and projtosurf_!=INPAR::XFEM::Proj_normal_phi
    {

      // check for potential L2_Projection smoothing
      const int l2_proj_num  = ( cond->GetInt("l2projsolv") + 1 );
      if(l2_proj_num<1)
        dserror("Issue with L2_PROJECTION_SOLVER, smaller than 1!!!");

      // SMOOTHED GRAD PHI!!!!!! (Create from nodal map on Xfluid discretization)
      // This method might be a bit too complicated as we need to save modphinp (size of fluid discretization),
      // as well as gradphinp_smoothed_rownode as the function ComputeNodalL2Projection returns the row node map.
      //----------------------------------------------
      // create the parameters for the discretization
      Teuchos::ParameterList eleparams;
      // action for elements
      eleparams.set<int>("action",FLD::presgradient_projection);

      //To get phi nodal values into pressure dofs in the fluid discretization!!! - any idea for nice implementation?
      const Epetra_Map* modphinp_dofrowmap = Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(cutter_dis_)->InitialDofRowMap();
      Teuchos::RCP<Epetra_Vector> modphinp = Teuchos::rcp(new Epetra_Vector(*modphinp_dofrowmap,true));

      double* val = cutter_phinp_->Values();

      int numrows = cutter_dis_->NumMyRowNodes();
      // loop all column nodes on the processor
      for(int lnodeid=0;lnodeid<numrows;++lnodeid)
      {
        // get the processor's local node
        DRT::Node* lsnode = cutter_dis_->lRowNode(lnodeid);
        if (lsnode == NULL)
          dserror("Returned node is null-pointer.");

        std::vector<int> initialdof = Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(cutter_dis_)->InitialDof(lsnode);

        if(initialdof.size()!=4)
          dserror("Initial Dof Size is not 4! Size: %d",initialdof.size());

        const int gid = initialdof[3];

        int err = modphinp->ReplaceGlobalValues(1,&val[lnodeid],&gid);
        if(err!=0)
          dserror("Something went wrong when replacing the values.");

      } // Loop over all nodes

      //SAFETY check
      // dependent on the desired projection, just remove this line
      if(not modphinp->Map().SameAs(*Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(cutter_dis_)->InitialDofRowMap()))
        dserror("input map is not a dof row map of the fluid");

      // set given state for element evaluation
      cutter_dis_->ClearState();
      Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(cutter_dis_)->SetInitialState(0,"pres",modphinp);

//      Teuchos::RCP<Epetra_MultiVector> gradphinp_smoothed_rownode= DRT::UTILS::ComputeNodalL2Projection(cutter_dis_,modphinp,"pres",3,eleparams,l2_proj_num);
      Teuchos::RCP<Epetra_MultiVector> gradphinp_smoothed_rownode= DRT::UTILS::ComputeNodalL2Projection(cutter_dis_,"pres",3,eleparams,l2_proj_num);
      if(gradphinp_smoothed_rownode==Teuchos::null)
        dserror("A smoothed grad phi is required, but an empty one is provided!");

      gradphinp_smoothed_node_ = Teuchos::rcp(new Epetra_MultiVector(*cutter_dis_->NodeColMap(),3));
      LINALG::Export(*gradphinp_smoothed_rownode,*gradphinp_smoothed_node_);

      //---------------------------------------------- // SMOOTHED GRAD PHI END
    }
  }

  // map the cutterdis-based phinp to the bgdis-noderowmap based phinp
  MapCutterToBgVector(
      cutter_dis_,
      cutter_phinp_,
      cutter_nds_phi_,
      bg_dis_,
      phinp_,
      bg_nds_phi_);

  // check if boundary position changed from the last step

  delta_phi->Update(1.0, *cutter_phinp_, -1.0); // phinp - phin

  double norm = 0.0;
  delta_phi->Norm2(&norm);

  return (norm > 1e-14); // did interface change?
}


//TODO: generalization in DRT::UTILS???
/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void XFEM::LevelSetCoupling::MapCutterToBgVector(
    const Teuchos::RCP<DRT::Discretization> & source_dis,
    const Teuchos::RCP<Epetra_Vector> &       source_vec_dofbased,
    const int                                 source_nds,
    const Teuchos::RCP<DRT::Discretization> & target_dis,
    const Teuchos::RCP<Epetra_Vector> &       target_vec_dofbased,
    const int                                 target_nds
)
{
  if( HaveMatchingNodes(source_dis, target_dis) ) // check for equal node positions
  {
    // here we assume that source_dis and target_dis are equal!

    // loop the nodes
    for(int lnodeid=0; lnodeid<target_dis->NumMyRowNodes(); ++lnodeid)
    {
      DRT::Node* node_source = source_dis->lRowNode(lnodeid);
      DRT::Node* node_target = target_dis->lRowNode(lnodeid);

      // get the set of source dof IDs for this node
      std::vector<int> lm_source;
      source_dis->Dof( source_nds, node_source,lm_source );

      std::vector<int> lm_target;
      target_dis->Dof( target_nds, node_target, lm_target );

      if(static_cast<int>(lm_source.size())!=1)
        dserror("we expect a unique dof per node here!");

      if(static_cast<int>(lm_target.size())!=1)
        dserror("we expect a unique dof per node here!");

      std::vector<double> val_source;
      DRT::UTILS::ExtractMyValues(*source_vec_dofbased, val_source, lm_source);

      // set to a dofrowmap based vector!
      const int lid_target = target_vec_dofbased->Map().LID(lm_target[0]);
      const int err = target_vec_dofbased->ReplaceMyValues(1,&val_source[0],&lid_target);
      if (err) dserror("could not replace values for convective velocity");
    }
  }
  else
  {
    dserror("nonmatching discretizations not supported so far! - Implement a mesh projector?");
  }
}



/*----------------------------------------------------------------------*
 | set interface level set field at current time           schott 02/15 |
 *----------------------------------------------------------------------*/
double XFEM::LevelSetCoupling::FunctImplementation(
    const int      func_no,
    const double * coords,
    const double t
)
{
//  dserror("you try to evaluate an implemented function for level-set field! Which one?");
  // WARNING!

  double x = coords[0];
  double y = coords[1];
  double z = coords[2];

  double val = 0.0;

  double R = 0.2;
  double r = 0.1;

  const double alpha = 0.6;


  if(func_no == -1)
  {


    // level set field for a helical pipe
    // sqrt( (x-Rcos(2 pi t))^2 + (y-Rsin(2 pi t))^2 + (z-alpha t)^2 -r = 0
    // with t(x,y,z) solution of minimization problem
    // dist((x,y,z), curve(t(x,y,z))) = min!

    // with curve(t) a parametrized curve (e.g. a helix)



    // NEWTON SYSTEM FOR SOLVING FOR t
    // d''(t) Delta_t = -d'(t)
    // t_n+1 = t_n + Delta_t

    // HELICAL CURVE z=alpha*t


    const double two_alpha_squared = 2.0*alpha*alpha;
    double two_PI = 2.0*PI;

    double t_0 = z/alpha;

    double Jac = 0.0;
    double rhs = 1.0;

    int armijo_steps = 50;

    int maxiter = 50;

    for(int i=0; i< maxiter; i++)
    {
      if(fabs(rhs)<1e-13) break;

      double arc= two_PI*t_0;
      double cosine = cos(arc);
      double sine   = sin(arc);
      Jac = 4.0*PI*R*(two_PI*x*cosine + two_PI*y*sine)+two_alpha_squared;
      rhs = 4.0*PI*R*(x*sine-y*cosine) + two_alpha_squared * t_0 - 2.0*alpha*z;


      double dt = -rhs/Jac;


      double armijo = 1.0;

      if(i<armijo_steps)
      {
        // it may happen, that the Newton direction is not a descent direction, then change the search direction
        // grad(f(x))^T * searchdir < 0 !   <=>  d'(t)*dt < 0   <=>  rhs*dt < 0
        if(dt*rhs > 0.0)
          dt*=-1.0;

        for(int l=0; l< 5; ++l)
        {
          if( l>0)
            armijo *= 0.5;

          // d(t+armijo*dt) < d(t) !!! and armijo (0,1] möglichst nahe an 1
          double t_new = t_0+armijo*dt;
          double arc_new = two_PI*t_new ;
          double cosine_new  = cos(arc_new);
          double sine_new    = sin(arc_new);

          double tmpx_new = x-R*cosine_new;
          double tmpy_new = y-R*sine_new;
          double tmpz_new = z-alpha*t_new;
          double norm1_squared = tmpx_new*tmpx_new+tmpy_new*tmpy_new+tmpz_new*tmpz_new;

          double tmpx = x-R*cosine;
          double tmpy = y-R*sine;
          double tmpz = z-alpha*t_0;
          double norm2_squared = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;

          if(norm1_squared < norm2_squared)
            break;
        }
      }

      t_0 += dt*armijo;

      if(i > maxiter-1)
      {
        std::cout << "Jac: " << Jac << std::endl;
        std::cout << "i: " << i << " rhs " << rhs << std::endl;
        std::cout << "armijo: " << armijo << std::endl;

        dserror("did not converge properly, intial guess not good enough - increase helixal height alpha!");
      }

    }


    double curve = alpha*t_0;

    double angle = two_PI*t_0;

    double cosine = cos(angle);
    double sine   = sin(angle);

    double tmp1 = x-R*cosine;
    double tmp2 = y-R*sine;
    double tmp3 = z-curve;

    double val_helix = sqrt(tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3)-r;

    return val_helix;
  }
  else if(func_no == -2)
  {
    double n1 = 0.0;
    double n2 = 2.0*PI*R;
    double n3 = alpha;

    double norm = sqrt(n1*n1+n2*n2+n3*n3);
    n1/=norm;
    n2/=norm;
    n3/=norm;

    // inflow region
    // point_on_plane_x
    double pop_x = 0.0; // here arbitrary
    double pop_y = 0.0;
    double pop_z = -2.0*alpha;

    double dist =  n1*pop_x + n2*pop_y + n3*pop_z;

    double val_plane_inflow = n1*x + n2*y + n3*z - dist;

    double val_inflow = std::max(val_plane_inflow, z-(pop_z+r*1.1));
    val_inflow = std::min(val_inflow, (z-(pop_z-r*1.1)));

    return -val_inflow;
  }
  else if(func_no == -3)
  {
    // outflow region

    double n1_out = 0.0;
    double n2_out = -2.0*PI*R;
    double n3_out = -alpha;

    double norm_out = sqrt(n1_out*n1_out+n2_out*n2_out+n3_out*n3_out);
    n1_out/=norm_out;
    n2_out/=norm_out;
    n3_out/=norm_out;

    // point_on_plane_x
    double pop_x_out = 0.0; // here arbitrary
    double pop_y_out = 0.0;
    double pop_z_out = +2.0*alpha;

    double dist_out =  n1_out*pop_x_out + n2_out*pop_y_out + n3_out*pop_z_out;

    double val_plane_outflow = n1_out*x + n2_out*y + n3_out*z - dist_out;

    double val_outflow = std::max(val_plane_outflow, -(z-(pop_z_out-r*1.1)));
    val_outflow = std::min(val_outflow, -(z-(pop_z_out+r*1.1)));

    return -val_outflow;
  }
  else if(func_no == -4)
  {
    double val_inner_ring_cyl = sqrt((x+0.2)*(x+0.2)+y*y)-0.14;
    double val_z_limit_inner = z+0.9;
    return std::max(val_inner_ring_cyl,val_z_limit_inner);
  }
  else if(func_no == -5)
  {
    double val_outer_ring_cyl = sqrt((x+0.2)*(x+0.2)+y*y)-0.22;
    double val_inner_ring_cyl = sqrt((x+0.2)*(x+0.2)+y*y)-0.14;
    double val_ring = std::max(val_outer_ring_cyl,-val_inner_ring_cyl);
    double val_z_limit_inner = z+0.9;
    double val_cylinder_ring_half = std::max(val_ring,val_z_limit_inner);
    return val_cylinder_ring_half;
  }
  else if(func_no == -6) // cylinder at inflow of a helix
  {
    double n1 = 0.0;
    double n2 = 2.0*PI*R;
    double n3 = alpha;

    double norm = sqrt(n1*n1+n2*n2+n3*n3);

    n1/=norm;
    n2/=norm;
    n3/=norm;

    // inflow region
    // point_on_plane_x
    double pop_x = 0.2; // here arbitrary
    double pop_y = 0.0;
    double pop_z = -2.0*alpha;

    double dist =  n1*pop_x + n2*pop_y + n3*pop_z;

    double val_plane_inflow = n1*x + n2*y + n3*z - dist;

    double coord_x_center = x-0.2;
    double coord_y_center = y;
    double coord_z_center = z+alpha*2.0;

    double coord_dot_n = coord_x_center*n1+coord_y_center*n2+coord_z_center*n3;

    double tmp1 = (coord_x_center-n1*coord_dot_n);
    tmp1*=tmp1;
    double tmp2 = (coord_y_center-n2*coord_dot_n);
    tmp2*=tmp2;
    double tmp3 = (coord_z_center-n3*coord_dot_n);
    tmp3*=tmp3;

    double val_cylinder = sqrt(tmp1+tmp2+tmp3)-r;
    val_cylinder = std::max(val_cylinder, val_plane_inflow);

    return val_cylinder;
  }
  else if(func_no == -7) // box for Oseen
  {
    //return -(std::max( (fabs(x-0.5+0.013))/0.3, (fabs(y-0.5+0.013))/0.3)-1.0);
    return -(std::max( (fabs(x-1.0))/0.45, std::max((fabs(y-0.5))/0.45, (fabs(z-0.5))/0.45) )-1.0);
  }


  //val = std::max(val_helix, std::max(val_inflow, val_outflow) );


////  double z_centerline = 1.0;
//  double alpha = 0.21;
//
//  double arc=0.0;
//  int n=0;
//
//  double z_max_at_arczero = r;
//
////  if(fabs(x)<1e-14 or fabs(y)<1e-14)
//////    if(fabs(x)<1e-14 and fabs(y)<1e-14)
////  {
////    val = 0.5;
////    return val;
////  }
//
//  double sgn_y = 1.0;
//
//  if(y>1e-14)
//    sgn_y= 1.0;
//  else if(y<1e-14)
//    sgn_y= -1.0;
//  else
//    sgn_y = 0.0;
//
//  double sgn_x = 1.0;
//
//  if(x>1e-14)
//    sgn_x= 1.0;
//  else if(x<1e-14)
//    sgn_x= -1.0;
//  else
//    sgn_x = 0.0;
//
//  n = 0;
//
//  if(z>=0.0)
//  {
//    // look for the first
//    for(int i = 0; ; i++)
//    {
//      if(i*alpha <= z and z<(i+1)*alpha)
//      {
//        n=i;
//        break;
//      }
//    }
//  }
//  else // z<0.0
//  {
//    // look for the first
//    for(int i = 0; ; i--)
//    {
//      if(i*alpha <= z and z<(i+1)*alpha)
//      {
//        n=i;
//        break;
//      }
//    }
//  }
//
//  // three possible i's
////  if(fabs(z-(n-1)*alpha) < fabs(z-n*alpha))
////    n--;
//
//  double arc_1 = 0.0;
//  double arc_2 = 0.0;
//  double arc_3 = 0.0;
//
//
//  if(fabs(x)>1e-14 and fabs(y)>1e-14)
//  {
//    arc_1 = sgn_y* 1.0/(2.0*PI)*acos(sgn_x/(sqrt(1+yy/xx))) + 0.5*(1.0 - sgn_y) +n;
//    arc_2 = sgn_y* 1.0/(2.0*PI)*acos(sgn_x/(sqrt(1+yy/xx))) + 0.5*(1.0 - sgn_y) +(n-1);
////    arc_3 = sgn_y* 1.0/(2.0*PI)*acos(sgn_x/(sqrt(1+yy/xx))) + 0.5*(1.0 - sgn_y) +(n-2);
//
//    //arc_3 = sgn_y* 1.0/(2.0*PI)*acos(sgn_x/(sqrt(1+yy/xx))) + 0.5*(1.0 - sgn_y) +(n+1);
//  }
//  arc= std::max(arc_1,arc_2);
// //   arc= std::max(arc_1,std::max(arc_2,arc_3)); //sgn_y* 1.0/(2.0*PI)*acos(sgn_x/(sqrt(1+yy/xx))) + 0.5*(1.0 - sgn_y) +n;
////  else
////  {
////    if(fabs(x)>1e-14 and fabs(y)<1e-14)
////    {
////      if(x>0)
////        arc = n;
////      else
////        arc = 0.5+n;
////    }
////    else if(fabs(x)<1e-14 and fabs(y)>1e-14)
////    {
////      if(y>0)
////        arc = 0.25+n;
////      else
////        arc = 0.75+n;
////    }
////    else
////    {
////      arc = 0.5 + n;
////    }
////  }
//
////  if(y>1e-14)
////  {
////    arc= sgn* 1.0/(2.0*PI)*acos(1.0/(sqrt(1+yy/xx))) +n;
////  }
////  else if(y < 1e-14)
////  {
////    arc= 1.0-1.0/(2.0*PI)*acos(1.0/(sqrt(1+yy/xx))) +n-1;
////  }
////  else
////  {
////    if(x > 0.0)
////      arc=n;
////    else if(x< 0.0)
////      arc=0.5+n;
////    else
////      arc=0.0;
////
////  }
//
//  double z_centerline = alpha * arc;
//
//  double tmp_1 = sqrt(xx+yy)-R;
//  double tmp_2 = z-z_centerline;
//
//  val = tmp_1*tmp_1 + tmp_2*tmp_2-r*r;

  return val;
}

//TODO: has_interface_moved_ checks its functionality and there is another flag in meshcoupling i think
XFEM::LevelSetCouplingBC::LevelSetCouplingBC(
     Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
     const std::string &                 cond_name,///< name of the condition, by which the derived cutter discretization is identified
     Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< full discretization from which the cutter discretization is derived
     const int                           coupling_id,///< id of composite of coupling conditions
     const double                        time,      ///< time
     const int                           step       ///< time step
 ) : LevelSetCoupling(bg_dis, cond_name, cond_dis, coupling_id, time, step),
     has_interface_moved_(true)
{
}

/*----------------------------------------------------------------------*
 | set interface level set field at current time           schott 02/15 |
 *----------------------------------------------------------------------*/
void XFEM::LevelSetCouplingBC::PrepareSolve()
{

  if(myrank_ == 0) IO::cout << "\t set level-set field, time " << time_ << IO::endl;

  has_interface_moved_ = SetLevelSetField(time_);
  return;
}



bool XFEM::LevelSetCouplingBC::HasMovingInterface()
{
  return has_interface_moved_;
}



void XFEM::LevelSetCouplingWeakDirichlet::EvaluateCouplingConditions(
    LINALG::Matrix<3,1>& ivel,
    LINALG::Matrix<3,1>& itraction,
    const LINALG::Matrix<3,1>& x,
    const DRT::Condition* cond
)
{
  // evaluate interface velocity (given by weak Dirichlet condition)
  EvaluateDirichletFunction(ivel, x, cond, time_);

  // no interface traction to be evaluated
  itraction.Clear();
}

// TODO: remove old state implementation?!
void XFEM::LevelSetCouplingWeakDirichlet::EvaluateCouplingConditionsOldState(
    LINALG::Matrix<3,1>& ivel,
    LINALG::Matrix<3,1>& itraction,
    const LINALG::Matrix<3,1>& x,
    const DRT::Condition* cond
)
{
  // evaluate interface velocity (given by weak Dirichlet condition)
  EvaluateDirichletFunction(ivel, x, cond, time_-dt_);

  // no interface traction to be evaluated
  itraction.Clear();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingWeakDirichlet::SetupConfigurationMap()
{
  //Configuration of Consistency Terms
  configuration_map_[INPAR::XFEM::F_Con_Row] = std::pair<bool,double>(true,1.0);
  configuration_map_[INPAR::XFEM::F_Con_Col] = std::pair<bool,double>(true,1.0);

  //Configuration of Adjount Consistency Terms
  configuration_map_[INPAR::XFEM::F_Adj_Row] = std::pair<bool,double>(true,1.0);
  configuration_map_[INPAR::XFEM::F_Adj_Col] = std::pair<bool,double>(true,1.0);

  //Configuration of Penalty Terms
  configuration_map_[INPAR::XFEM::F_Pen_Row] = std::pair<bool,double>(true,1.0);
  configuration_map_[INPAR::XFEM::F_Pen_Col] = std::pair<bool,double>(true,1.0);
  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingWeakDirichlet::UpdateConfigurationMap_GP(
    double& kappa_m,
    double& visc_m,
    double& visc_s,
    double& visc_stab,
    double& full_stab,
    const LINALG::Matrix<3,1>& x,
    const DRT::Condition* cond)
{
  //Configuration of Penalty Terms
  configuration_map_[INPAR::XFEM::F_Pen_Row].second = full_stab;

  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingNeumann::EvaluateCouplingConditions(
    LINALG::Matrix<3,1>& ivel,
    LINALG::Matrix<3,1>& itraction,
    const LINALG::Matrix<3,1>& x,
    const DRT::Condition* cond
)
{
  // no interface velocity to be evaluated
  ivel.Clear();

  // evaluate interface traction (given by Neumann condition)
  EvaluateNeumannFunction(itraction, x, cond, time_);
}

//TODO: combine it with the function before with optional time parameter?!
/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingNeumann::EvaluateCouplingConditionsOldState(
    LINALG::Matrix<3,1>& ivel,
    LINALG::Matrix<3,1>& itraction,
    const LINALG::Matrix<3,1>& x,
    const DRT::Condition* cond
)
{
  // no interface velocity to be evaluated
  ivel.Clear();

  // evaluate interface traction (given by Neumann condition)
  EvaluateNeumannFunction(itraction, x, cond, time_-dt_);
}

// TODO: inheritance from the bc seems to be the wrong concept, more delegate functionality to a Dirichlet object and a Neumann object...
// the same is implemented in mesh coupling object again! wrong concept!
/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
XFEM::LevelSetCouplingNavierSlip::LevelSetCouplingNavierSlip(
    Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
    const std::string &                 cond_name, ///< name of the condition, by which the derived cutter discretization is identified
    Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< full discretization from which the cutter discretization is derived
    const int                           coupling_id,///< id of composite of coupling conditions
    const double                        time,      ///< time
    const int                           step       ///< time step
) : LevelSetCouplingBC(bg_dis, cond_name, cond_dis, coupling_id, time, step)
{
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingNavierSlip::SetElementConditions()
{
  XFEM::LevelSetCouplingBC::SetElementConditions();

  if(cutterele_conds_.size()==0)
    dserror("call SetElementConditions() first!");

  DRT::Condition* cond = cutterele_conds_[0].second; // get condition of first element

  // Get robin coupling IDs
  robin_dirichlet_id_ = cond->GetInt("robin_id_dirch");
  robin_neumann_id_   = cond->GetInt("robin_id_neumann");

  has_neumann_jump_ = (robin_neumann_id_ < 0) ? false : true;

  if (has_neumann_jump_)
  {
    std::cout << "#########################################################################################################\n";
    std::cout << "#########################################################################################################\n";
    std::cout << "### WARNING:: XFEM::LevelSetCouplingNavierSlip                              The traction jump is      ###\n";
    std::cout << "### divided by the dynviscosity on Gausspoint Level, this might be expensed and not really necessary! ###\n";
    std::cout << "#########################################################################################################\n";
    std::cout << "#########################################################################################################" << std::endl;
  }

  // set the navier-slip specific element conditions
  SetElementSpecificConditions(cutterele_cond_robin_dirichlet_,"XFEMRobinDirichletVol",robin_dirichlet_id_);
  if(has_neumann_jump_)
    SetElementSpecificConditions(cutterele_cond_robin_neumann_,"XFEMRobinNeumannVol",robin_neumann_id_);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingNavierSlip::SetElementSpecificConditions(
    std::vector<DRT::Condition* > & cutterele_cond,
    const std::string & cond_name,
    const int & robin_id)
{

  //TODO: can we combine this function with SetElementConditions in the coupling base routine!

  // number of column cutter boundary elements
  int nummycolele = cutter_dis_->NumMyColElements();

  cutterele_cond.clear();
  cutterele_cond.reserve(nummycolele);

  //// initialize the vector invalid coupling-condition type "NONE"
  //EleCoupCond init_pair = EleCoupCond(INPAR::XFEM::CouplingCond_NONE,NULL);
  for(int lid=0; lid<nummycolele; ++lid) cutterele_cond.push_back(NULL);

  //-----------------------------------------------------------------------------------
  // loop all column cutting elements on this processor
  for(int lid=0; lid<nummycolele; ++lid)
  {
    DRT::Element* cutele = cutter_dis_->lColElement(lid);

    // get all conditions with given condition name
    std::vector<DRT::Condition*> mycond;
    DRT::UTILS::FindElementConditions(cutele, cond_name, mycond);

    std::vector<DRT::Condition*> mynewcond;
    GetConditionByRobinId(mycond, robin_id, mynewcond);

    DRT::Condition* cond_unique = NULL;

    // safety checks
    if(mynewcond.size() != 1)
    {
      dserror("%i conditions of the same name with robin id %i, for element %i! %s coupling-condition not unique!", mynewcond.size(), (robin_id+1), cutele->Id(), cond_name.c_str());

    }
    else if(mynewcond.size() == 1) // unique condition found
    {
      cond_unique = mynewcond[0];
    }

    // store the unique condition pointer to the cutting element
    cutterele_cond[lid] = cond_unique;
  }
  //  //-----------------------------------------------------------------------------------
  //  // check if all column cutter elements have a valid condition type
  //  // loop all column cutting elements on this processor
  for(int lid=0; lid<nummycolele; ++lid)
  {
    if(cutterele_cond[lid] == NULL)
      dserror("cutter element with local id %i has no Robin-condition!!!", lid);
  }

}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingNavierSlip::EvaluateCouplingConditions(
    LINALG::Matrix<3,1>& ivel,
    LINALG::Matrix<3,1>& itraction,
    const LINALG::Matrix<3,1>& x,
    const DRT::Condition* cond
)
{
  if(cutterele_cond_robin_dirichlet_.size()==0)
    dserror("initialize cutterele_cond_robin_dirichlet_ first!");

  // evaluate interface velocity (given by weak Dirichlet condition)
  EvaluateDirichletFunction(ivel, x, cutterele_cond_robin_dirichlet_[0], time_);

  // evaluate interface traction (given by Neumann condition)
  if(has_neumann_jump_)
  {
    if(cutterele_cond_robin_neumann_.size()==0)
      dserror("initialize cutterele_cond_robin_neumann_ first!");

    EvaluateNeumannFunction(itraction, x, cutterele_cond_robin_neumann_[0], time_);
  }

}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingNavierSlip::EvaluateCouplingConditionsOldState(
    LINALG::Matrix<3,1>& ivel,
    LINALG::Matrix<3,1>& itraction,
    const LINALG::Matrix<3,1>& x,
    const DRT::Condition* cond
)
{
  if(cutterele_cond_robin_dirichlet_.size()==0)
    dserror("initialize cutterele_cond_robin_dirichlet_ first!");

  // evaluate interface velocity (given by weak Dirichlet condition)
  EvaluateDirichletFunction(ivel, x, cutterele_cond_robin_dirichlet_[0], time_-dt_);

  // evaluate interface traction (given by Neumann condition)
  if(has_neumann_jump_)
  {
    if(cutterele_cond_robin_neumann_.size()==0)
      dserror("initialize cutterele_cond_robin_neumann_ first!");

    EvaluateNeumannFunction(itraction, x, cutterele_cond_robin_neumann_[0], time_-dt_);
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingNavierSlip::GetSlipCoefficient(
    double& slipcoeff,
    const LINALG::Matrix<3,1>& x,
    const DRT::Condition* cond
)
{
  if(is_constant_sliplength_)
    slipcoeff = sliplength_;
  else
    EvaluateScalarFunction(slipcoeff,x.A(),sliplength_,cond,time_);

}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingNavierSlip::SetConditionSpecificParameters()
{
  if(cutterele_conds_.size()==0)
    dserror("call SetElementConditions() first!");

  DRT::Condition* cond = cutterele_conds_[0].second; // get condition of first element

  // Get the scaling factor for the slip length
  sliplength_ = cond->GetDouble("slipcoeff");

  //Temporary variable for readability.
  bool tmp_bool;

  // Is the slip length constant? Don't call functions at GP-level unnecessary.
  tmp_bool = (cond->GetInt("curve") < 0 and cond->GetInt("funct") < 1);
  is_constant_sliplength_ = (tmp_bool) ? true : false;

  // Project the prescribed velocity in tangential direction, to remove "spurious velocities"
  //  from the geometry approximation.
  tmp_bool=( (cond->GetInt("force_tang_vel")) ==0 );
  forcetangvel_ = (tmp_bool) ? false : true;

}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingNavierSlip::GetConditionByRobinId(
    const std::vector<DRT::Condition*> & mycond,
    const int coupling_id,
    std::vector<DRT::Condition*> & mynewcond
)
{
  mynewcond.clear();

  // select the conditions with specified "couplingID"
  for(size_t i=0; i< mycond.size(); ++i)
  {
    DRT::Condition* cond = mycond[i];
    const int id = cond->GetInt("robin_id");

    if(id == coupling_id)
      mynewcond.push_back(cond);
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingNavierSlip::SetupConfigurationMap()
{
  if (GetAveragingStrategy() == INPAR::XFEM::Xfluid_Sided)
  {
    //Configuration of Consistency Terms
    configuration_map_[INPAR::XFEM::F_Con_Row] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::F_Con_Col] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::F_Con_t_Row] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::F_Con_t_Col] = std::pair<bool,double>(true,1.0);

    //Configuration of Adjount Consistency Terms
    configuration_map_[INPAR::XFEM::F_Adj_n_Row] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::F_Adj_n_Col] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::F_Adj_t_Row] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::F_Adj_t_Col] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::FStr_Adj_t_Col] = std::pair<bool,double>(true,1.0);

    //Configuration of Penalty Terms
    configuration_map_[INPAR::XFEM::F_Pen_n_Row] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::F_Pen_n_Col] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::F_Pen_t_Row] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::F_Pen_t_Col] = std::pair<bool,double>(true,1.0);
  }
  else if (GetAveragingStrategy() == INPAR::XFEM::invalid)
    dserror("XFEM::LevelSetCouplingNavierSlip: Averaging Strategy not set!");
  else
    dserror("XFEM::LevelSetCouplingNavierSlip: You want to initialize another strategy than Xfluid_Sided?");
  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingNavierSlip::UpdateConfigurationMap_GP(
    double& kappa_m,
    double& visc_m,
    double& visc_s,
    double& visc_stab,
    double& full_stab,
    const LINALG::Matrix<3,1>& x,
    const DRT::Condition* cond)
{

  double dynvisc   = (kappa_m*visc_m + (1.0-kappa_m)*visc_s);
  double sliplength = 0.0;
  GetSlipCoefficient(sliplength,x,cond);

  if(sliplength < 0.0)
    dserror("The slip length can not be negative.");

  if ( sliplength != 0.0 )
  {
    double stabnit = 0.0;
    double stabadj = 0.0;
    XFEM::UTILS::GetNavierSlipStabilizationParameters(full_stab,visc_stab,dynvisc,sliplength,stabnit,stabadj);
    configuration_map_[INPAR::XFEM::F_Pen_t_Row].second = stabnit;
    configuration_map_[INPAR::XFEM::F_Con_t_Row] = std::pair<bool,double>(true,-stabnit); //+sign for penalty!
    configuration_map_[INPAR::XFEM::F_Con_t_Col] = std::pair<bool,double>(true,sliplength/dynvisc);
    configuration_map_[INPAR::XFEM::F_Adj_t_Row].second = stabadj;
    configuration_map_[INPAR::XFEM::FStr_Adj_t_Col] = std::pair<bool,double>(true,sliplength);
  }
  else
  {
    configuration_map_[INPAR::XFEM::F_Pen_t_Row].second = visc_stab;
    configuration_map_[INPAR::XFEM::F_Con_t_Row] = std::pair<bool,double>(false,0.0);
    configuration_map_[INPAR::XFEM::F_Con_t_Col] = std::pair<bool,double>(false,0.0);
    configuration_map_[INPAR::XFEM::F_Adj_t_Row].second = 1.0;
    configuration_map_[INPAR::XFEM::FStr_Adj_t_Col] = std::pair<bool,double>(false,0.0);
  }

  //Configuration of Penalty Terms
  configuration_map_[INPAR::XFEM::F_Pen_n_Row].second = visc_stab; //full_stab <-- to keep results!

  return;
}


/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
XFEM::LevelSetCouplingTwoPhase::LevelSetCouplingTwoPhase(
    Teuchos::RCP<DRT::Discretization>&  bg_dis,    ///< background discretization
    const std::string &                 cond_name, ///< name of the condition, by which the derived cutter discretization is identified
    Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< discretization from which the cutter discretization is derived
    const int                           coupling_id,///< id of composite of coupling conditions
    const double                        time,      ///< time
    const int                           step       ///< time step
) : LevelSetCoupling(bg_dis, cond_name,  cond_dis, coupling_id, time, step),
    twophasedyn_(DRT::Problem::Instance()->TwoPhaseFlowParams()), // access parameter for two phase flow TODO: remove Problem
    surftensapprox_(INPAR::TWOPHASE::surface_tension_approx_none),
    laplacebeltrami_(INPAR::TWOPHASE::matrix_mixed_smoothed),
    surfacetension_init_(false),
    require_smoothedgradphi_(false),
    require_nodalcurvature_(false),
    transport_direction_(INPAR::TWOPHASE::transport_dir_all),
    cutter_nds_vel_(-1),
    col_vectors_valid_(false)
{
}


//TODO: needs to be generalized!
//TODO: allow to create a cutter dis from condition set on the conddis = full scatra-dis
/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingTwoPhase::SetCutterDiscretization()
{
  //TODO: for the moment we assume matching! introduce a check here, maybe directly clone the dis here, or extract one
  //if(matching) then clone strategy, otherwise create Dis from condition...
  cutter_dis_ = cond_dis_;

  have_nodematching_dis_ = HaveMatchingNodes(cutter_dis_, bg_dis_);

  if(not(dofset_coupling_map_.count("phi_in_scatra") == 1))
      dserror("phi_in_scatra-dofset not set in dofset_coupling_map for scatra dis!");

  if(not(dofset_coupling_map_.count("vel_fluid_proxy_in_scatra") == 1))
      dserror("vel_fluid_proxy_in_scatra-dofset not set in dofset_coupling_map for scatra dis!");

  cutter_nds_phi_ = dofset_coupling_map_["phi_in_scatra"];             // dofset id for scalar field
  cutter_nds_vel_ = dofset_coupling_map_["vel_fluid_proxy_in_scatra"]; // dofset id for transport velocity

}


/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingTwoPhase::InitStateVectors()
{
  XFEM::LevelSetCoupling::InitStateVectors();

  const Epetra_Map * cutter_dofrowmap = cutter_dis_->DofRowMap(cutter_nds_vel_);
  cutter_transport_vel_ = LINALG::CreateVector( *cutter_dofrowmap, true);
}


/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingTwoPhase::SetConditionSpecificParameters()
{

  DRT::Condition* cond = cutterele_conds_[0].second; // get condition of first element

  // set parameters
  SetParameters_Physical(cond);

  SetParameters_SurfaceTension();

  // set flags
  SetFlags_InterfaceTransport(cond);

  SetFlags_GeometricQuantities();
}


/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingTwoPhase::SetParameters_SurfaceTension()
{
  const Teuchos::ParameterList & params = twophasedyn_.sublist("SURFACE TENSION");

  surftensapprox_  = DRT::INPUT::IntegralValue<INPAR::TWOPHASE::SurfaceTensionApprox>(params,"SURFTENSAPPROX");
  laplacebeltrami_ = DRT::INPUT::IntegralValue<INPAR::TWOPHASE::LaplaceBeltramiCalc>(params,"LAPLACE_BELTRAMI");

  //SAFETY-CHECKS
  if(DRT::INPUT::IntegralValue<bool>(params,"L2_PROJECTION_SECOND_DERIVATIVES"))
    dserror("Second L2-projected derivatives can not be calculated as of now for the Level Set.");

  if(DRT::INPUT::IntegralValue<INPAR::TWOPHASE::SmoothGradPhi>(params,"SMOOTHGRADPHI")!=INPAR::TWOPHASE::smooth_grad_phi_l2_projection)
    dserror("No other smoothing for the gradient of the level set other than L2 is allowed for now.");

  if(DRT::INPUT::IntegralValue<INPAR::TWOPHASE::NodalCurvatureCalc>(params,"NODAL_CURVATURE")!=INPAR::TWOPHASE::l2_projected)
    dserror("No other way to calculate the nodal curvature than L2.");

  if(twophasedyn_.sublist("SURFACE TENSION").get<double>("SMOOTHING_PARAMETER")!=0.0)
    dserror("No smoothing available for now.");


  surfacetension_init_ = true;

}


/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingTwoPhase::SetFlags_InterfaceTransport(
    DRT::Condition * cond
)
{
  //TODO: think about alternatives for twophase flow! more accurate interface transport needed
  // normal transport is implemented, however not fully tested so far! (see combustion coupling object!)
  transport_direction_ = INPAR::TWOPHASE::transport_dir_all;
}


/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingTwoPhase::SetFlags_GeometricQuantities()
{
  CheckInit_SurfaceTension();

  // options for classical two-phase flow
  if(surftensapprox_==INPAR::TWOPHASE::surface_tension_approx_divgrad_normal)
  {
    require_smoothedgradphi_ = true;
  }
  else if(surftensapprox_==INPAR::TWOPHASE::surface_tension_approx_nodal_curvature)
  {
    require_nodalcurvature_ = true;
  }
  else if(surftensapprox_==INPAR::TWOPHASE::surface_tension_approx_laplacebeltrami)
  {
    if(laplacebeltrami_ != INPAR::TWOPHASE::matrix_non_smoothed)
      require_smoothedgradphi_ = true;
    else
      require_smoothedgradphi_ = false;
  }
  else
  {
    require_smoothedgradphi_ = false;
    require_nodalcurvature_  = false;
  }

  if(transport_direction_==INPAR::TWOPHASE::transport_dir_normal)
    require_smoothedgradphi_ = true;
}


/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingTwoPhase::WriteAccess_GeometricQuantities(
    Teuchos::RCP<Epetra_Vector> &      scalaraf,
    Teuchos::RCP<Epetra_MultiVector> & smoothed_gradphiaf,
    Teuchos::RCP<Epetra_Vector> &      curvatureaf
   )
{
  if(cutter_phinp_==Teuchos::null) dserror("cutter_phinp null pointer");

  scalaraf           = cutter_phinp_;
  smoothed_gradphiaf = (require_smoothedgradphi_) ? gradphinp_smoothed_node_ : Teuchos::null;
  curvatureaf        = (require_nodalcurvature_)  ? curvaturenp_node_        : Teuchos::null;

  // column vectors might not be valid anymore
  col_vectors_valid_ =  false;
}


/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingTwoPhase::ExportGeometricQuantities()
{
  std::cout << "CALLED EXPORT ROUTINE!!!" << std::endl;

  if(require_smoothedgradphi_)
      LINALG::Export(*gradphinp_smoothed_node_, *gradphinp_smoothed_node_col_);

  if(require_nodalcurvature_)
    LINALG::Export(*curvaturenp_node_, *curvaturenp_node_col_);

  // column vectors are uptodate gain
  col_vectors_valid_=true;
}


/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Vector> XFEM::LevelSetCouplingTwoPhase::ComputeTransportVelocity(
    const Teuchos::RCP<GEO::CutWizard> & wizard,                   ///< the cut wizard
    const Teuchos::RCP<const Epetra_Vector>& convective_velocity   ///< the convective fluid velocity based on initial dofmap (just velocity)
)
{
  CheckForValidVectors();

  // Check for Matching nodes
  if(!have_nodematching_dis_)
    dserror("introduce mappings between non-node-matching discretizations here!");

  //------------------------
  // compute the flame velocity at nodes of the cutter-dis
  //------------------------

  // loop over nodes on this processor
  for(int lnodeid=0; lnodeid < cutter_dis_->NumMyRowNodes(); ++lnodeid)
  {
    // get the current node
    DRT::Node* cutter_node = cutter_dis_->lRowNode(lnodeid);

    // compute the flame velocity
    // u_flame = u_convective + u_relative_flame
    // u_convective     = flvelconv =    the discontinuous fluid velocity u           (OPTION 1)
    //                                OR more precisely the normal part (u*n^ij)*n^ij (OPTION 2)
    // u_relative_flame = flvelrel  = the discontinuous relative flame velocity ( M/rho^i in Omega^i and M/rho^j in Omega^j
    //                             where the mass flow rate M = -rho^i*sl with sl the laminar flamespeed w.r.t Omega

    //-------------------------------------------------------------
    // get (smoothed) gradient of the G-function field at this node
    //-------------------------------------------------------------
    // smoothed normal vector at this node
    Epetra_SerialDenseVector nvec(nsd_);
    double curv = 0.0;

    GetSmoothedQuantitiesAtNode(nvec, curv, cutter_node, gradphinp_smoothed_node_col_, curvaturenp_node_col_, cutter_dis_, cutter_nds_phi_);


    if(require_smoothedgradphi_)
    {
      // do normalization or a manipulation of the smoothed normal vector
      if(!RescaleNormal(nvec))
      {
        std::cout << "\n/!\\ phi gradient too small at node "
            << cutter_node->Id() << " -> interface velocity is only the convective velocity" << std::endl;
      }

      // possible change in normal direction compared to phi-gradient
      nvec.Scale(normal_orientation_);
    }

    //------------------------
    // get fluid data at current scatra's node coordinates
    //------------------------

    // get the positioning of the node
    GEO::CUT::Point::PointPosition pos = GEO::CUT::Point::undecided;

    DRT::Node * bg_node = NULL;

    if(have_nodematching_dis_)
    {
      const int bgnode_id = cutter_node->Id(); // this only holds for matching discretizations!!!

      // get the cut node and its position
      GEO::CUT::Node * n = wizard->GetNode( bgnode_id );

      // get the points position!
      pos = n->Position();

      bg_node = bg_dis_->gNode(bgnode_id);
    }
    else
    {
      // const double * cutter_node_xyz = cutter_node->X();
      // search for the fluid element, the current scatra node lies in to obtain the node's GEO::CUT::Point::PointPosition
      // interpolate fluid quantities within the element afterwards

      dserror("ask a mesh projector with a search tree for the right fluid phase!");
    }

    Epetra_SerialDenseVector flvelconv(nsd_);  // the convective fluid velocity (Navier-Stokes) at this node
    Epetra_SerialDenseVector flvelrel (nsd_);  // the relative interface(flame) velocity at this node
    Epetra_SerialDenseVector flvelabs (nsd_);  // the absolute interface(flame) velocity (transport velocity for interface)

    ComputeRelativeTransportVelocity(flvelrel, pos, nvec, curv);


    //-----------------------------------------------
    // compute (absolute) flame velocity at this node
    //-----------------------------------------------

    // get the set of dof IDs for this node (nsd x vel + 1 x pressure) from enriched cutFEM dofset
    //TODO: get this dof-index from map!
    const int nds = 0;
    std::vector<int> lm_fld;
    bg_dis_->Dof( lm_fld, bg_node, 0, nds );

    if(lm_fld.size()!=(nsd_+1))
      dserror("we expect fluid and pressure dofs here!");

    // remove the pressure dof
    lm_fld.erase(lm_fld.end()-1);

    DRT::UTILS::ExtractMyValues(*convective_velocity, flvelconv, lm_fld);


    if(transport_direction_==INPAR::TWOPHASE::transport_dir_normal) // OPTION 2: just the normal part (u*n^ij)*n^ij (OPTION 2)
    {
      const double normal_vel = flvelconv.Dot(nvec);
      flvelconv=nvec;
      flvelconv.Scale(normal_vel);
    }
    else if(transport_direction_==INPAR::TWOPHASE::transport_dir_all)
    {
      // ELSE: nothing to do -- // OPTION 1: use the full convective velocity (normal and tangential parts)
      // just continue in this loop
    }
    else
      dserror("unsupported type of transport option for interface");


    flvelabs+=flvelconv;
    flvelabs+=flvelrel;


    const std::vector<int> lm_vel = cutter_dis_->Dof(cutter_nds_vel_, cutter_node);

    if(lm_vel.size()!=(nsd_+1))
      dserror("assume nsd_ dofs in cutterdis-Dofset for transport velocity");

    // add fluid velocity (Navier Stokes solution) and relative flame velocity, in addition set the pressure dof to zero!
    for (size_t icomp=0; icomp<nsd_+1; ++icomp)
    {
      double tmp = 0.0;
      if(icomp!=nsd_)
        tmp = flvelabs(icomp);

      const int gid = lm_vel[icomp];
      const int lid = cutter_transport_vel_->Map().LID(gid);
      const int err = cutter_transport_vel_->ReplaceMyValues(1,&tmp,&lid);
      if (err) dserror("could not replace values for convective velocity");

    }

#if(0)
    std::cout << "flvelconv" << flvelconv << std::endl;
    std::cout << "flvelrel"  << flvelrel  << std::endl;
    std::cout << "computed flvelabs"  << flvelabs  << std::endl;
#endif
  }

  return cutter_transport_vel_;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingTwoPhase::ComputeRelativeTransportVelocity(
    Epetra_SerialDenseVector & flvelrel,
    const GEO::CUT::Point::PointPosition & position,
    const Epetra_SerialDenseVector & nvec,
    const double & curv
)
{
  flvelrel.Scale(0.0);
}



/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
bool XFEM::LevelSetCouplingTwoPhase::RescaleNormal(
    Epetra_SerialDenseVector & normal
)
{
  // compute norm of smoothed normal vector
  const double gradphi_norm = normal.Norm2();

  if (gradphi_norm > 1.0E-12) // the standard case of non-vanishing normal
  {
    normal.Scale(1.0/gradphi_norm); // scale it norm 1
    return true;
  }
  else // 'ngradnorm' == 0.0
  {
    //TODO: this still needs to be checked!!!

    // length of smoothed normal is zero at this node -> node must be on a singularity of the
    // level set function (e.g. "regular level set cone"); all normals add up to zero normal vector
    // -> The fluid convective velocity 'fluidvel' alone constitutes the flame velocity, since the
    //    relative flame velocity 'flvelrel' turns out to be zero due to the zero average normal vector.
    // get the global id for current node
    normal.Scale(0.0);
    return false;
  }
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingTwoPhase::GetInterfaceSlaveMaterial(
  DRT::Element* actele,
  Teuchos::RCP<MAT::Material> & mat
)
{
  XFEM::UTILS::GetVolumeCellMaterial(actele,mat,GEO::CUT::Point::inside);
}



/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingTwoPhase::SetupConfigurationMap()
{
  if (GetAveragingStrategy() == INPAR::XFEM::Harmonic)
  {
    //Configuration of Consistency Terms
    configuration_map_[INPAR::XFEM::F_Con_Row] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::F_Con_Col] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::X_Con_Row] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::X_Con_Col] = std::pair<bool,double>(true,1.0);

    //Configuration of Adjount Consistency Terms
    configuration_map_[INPAR::XFEM::F_Adj_Row] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::F_Adj_Col] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::X_Adj_Row] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::X_Adj_Col] = std::pair<bool,double>(true,1.0);

    //Configuration of Penalty Terms
    configuration_map_[INPAR::XFEM::F_Pen_Row] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::F_Pen_Col] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::X_Pen_Row] = std::pair<bool,double>(true,1.0);
    configuration_map_[INPAR::XFEM::X_Pen_Col] = std::pair<bool,double>(true,1.0);

    //Initialize Traction Jump Terms (Also Remove other version for safety!)
    if(surftensapprox_==INPAR::TWOPHASE::surface_tension_approx_laplacebeltrami)
    {
      configuration_map_[INPAR::XFEM::F_LB_Rhs] = std::pair<bool,double>(true,1.0);
      configuration_map_[INPAR::XFEM::X_LB_Rhs] = std::pair<bool,double>(true,1.0);
      configuration_map_[INPAR::XFEM::F_TJ_Rhs] = std::pair<bool,double>(false,0.0);
      configuration_map_[INPAR::XFEM::X_TJ_Rhs] = std::pair<bool,double>(false,0.0);
    }
    else
    {
      configuration_map_[INPAR::XFEM::F_LB_Rhs] = std::pair<bool,double>(false,0.0);
      configuration_map_[INPAR::XFEM::X_LB_Rhs] = std::pair<bool,double>(false,0.0);
      configuration_map_[INPAR::XFEM::F_TJ_Rhs] = std::pair<bool,double>(true,1.0);
      configuration_map_[INPAR::XFEM::X_TJ_Rhs] = std::pair<bool,double>(true,1.0);
    }
  }
  else if (GetAveragingStrategy() == INPAR::XFEM::invalid)
    dserror("XFEM::LevelSetCouplingTwoPhase: Averaging Strategy not set!");
  else
    dserror("XFEM::LevelSetCouplingTwoPhase: You want to initialize another strategy than harmonic?");
  return;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingTwoPhase::UpdateConfigurationMap_GP(
    double& kappa_m,
    double& visc_m,
    double& visc_s,
    double& visc_stab,
    double& full_stab,
    const LINALG::Matrix<3,1>& x,
    const DRT::Condition* cond)
{
  //Configuration of Penalty Terms
  configuration_map_[INPAR::XFEM::F_Pen_Row].second = full_stab;
  configuration_map_[INPAR::XFEM::X_Pen_Row].second = full_stab;

  if (GetAveragingStrategy() == INPAR::XFEM::Harmonic)
  {
    //Configuration of Consistency Terms
    configuration_map_[INPAR::XFEM::F_Con_Col].second = kappa_m;
    configuration_map_[INPAR::XFEM::X_Con_Col].second = 1.-kappa_m;
    //Configuration of Adjount Consistency Terms
    configuration_map_[INPAR::XFEM::F_Adj_Row].second = kappa_m;
    configuration_map_[INPAR::XFEM::X_Adj_Row].second = 1.-kappa_m;

    //Traction Jump Terms
    if(surftensapprox_==INPAR::TWOPHASE::surface_tension_approx_laplacebeltrami)
    {
      configuration_map_[INPAR::XFEM::F_LB_Rhs].second = 1.-kappa_m;
      configuration_map_[INPAR::XFEM::X_LB_Rhs].second = kappa_m;
    }
    else
    {
      configuration_map_[INPAR::XFEM::F_TJ_Rhs].second = 1.-kappa_m;
      configuration_map_[INPAR::XFEM::X_TJ_Rhs].second = kappa_m;
    }
  }
  else
    dserror("XFEM::LevelSetCouplingTwoPhase: You want to initialize another strategy than harmonic?");
  return;
}


/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingTwoPhase::GetViscositySlave(
    DRT::Element * coup_ele,                   ///< xfluid ele
    double& visc_s)                            ///< viscosity slavesided
{
  Teuchos::RCP<MAT::Material> mat_s;
  XFEM::UTILS::GetVolumeCellMaterial(coup_ele,mat_s,GEO::CUT::Point::inside);
  if (mat_s->MaterialType() == INPAR::MAT::m_fluid)
    visc_s = Teuchos::rcp_dynamic_cast<MAT::NewtonianFluid>(mat_s)->Viscosity();
  else
    dserror("GetCouplingSpecificAverageWeights: Slave Material not a fluid material?");

  return;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingTwoPhase::GetCouplingSpecificAverageWeights(
    DRT::Element * xfele,                      ///< xfluid ele
    DRT::Element * coup_ele,                   ///< coup_ele ele
    double & kappa_m)
{
  if (GetAveragingStrategy() == INPAR::XFEM::Harmonic)
  {
    //Get Materials of master and slave
    double visc_m = 0.0;
    GetViscosityMaster(xfele,visc_m);

    double visc_s = 0.0;
    GetViscositySlave(coup_ele,visc_s);

    kappa_m = visc_s/(visc_m+visc_s);

    if ( kappa_m > 1.0 || kappa_m < 0.0) dserror("Nitsche weights for inverse estimate kappa_m lies not in [0,1]: %d", kappa_m);
  }
  else
    dserror("XFEM::LevelSetCouplingTwoPhase: GetCouplingSpecificAverageWeights not implemented for this averaging strategy!");
  return;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
//! get the smoothed level set gradient at a given node (not necessarily normalized to one)
void XFEM::LevelSetCouplingTwoPhase::GetSmoothedQuantitiesAtNode(
    Epetra_SerialDenseVector & normal,
    double & curvature,
    const DRT::Node* node,
    const Teuchos::RCP<const Epetra_MultiVector> & gradphinp_smoothed_node_col,
    const Teuchos::RCP<const Epetra_MultiVector> & curvaturenp_node_col,
    Teuchos::RCP<DRT::Discretization> & dis,
    const int nds
)
{
  const std::vector<int> lm = dis->Dof(nds, node);


  if(lm.size()!=1) dserror("assume a unique level-set dof in cutterdis-Dofset");

  if(require_nodalcurvature_)
  {
    std::vector<double> local_curvature(1);
    DRT::UTILS::ExtractMyValues(*curvaturenp_node_col, local_curvature, lm);

    if(local_curvature.size() != 1)
      dserror("wrong size of (potentially resized) local matrix!");

    curvature = local_curvature[0];
  }

  if(require_smoothedgradphi_)
  {
    std::vector<double> local_normal(nsd_);
    DRT::UTILS::ExtractMyValues(*gradphinp_smoothed_node_col, local_normal, lm);

    if(local_normal.size() != nsd_)
      dserror("wrong size of (potentially resized) local matrix!");

    // copy local to nvec....
     std::copy(local_normal.begin(), local_normal.begin()+nsd_, normal.A());
  }
}


/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
//! get the smoothed level set gradient at a given node (not necessarily normalized to one)
void XFEM::LevelSetCouplingTwoPhase::GetSmoothedQuantitiesAtElement(
    Epetra_SerialDenseMatrix & normal,
    Epetra_SerialDenseMatrix & curvature,
    const DRT::Element* element,
    const Teuchos::RCP<const Epetra_MultiVector> & gradphinp_smoothed_node_col,
    const Teuchos::RCP<const Epetra_MultiVector> & curvaturenp_node_col,
    Teuchos::RCP<DRT::Discretization> & dis,
    const int nds
)
{
  // get the other nds-set which is connected to the current one via this boundary-cell
  DRT::Element::LocationArray la( 1 );
  element->LocationVector(*dis, la, false );

  const size_t numnode = element->NumNode();

  if(la[0].lm_.size()!=numnode)
    dserror("assume a unique level-set dof in cutterdis-Dofset per node");

  if(require_nodalcurvature_)
  {
    std::vector<double> local_curvature(numnode);
    DRT::UTILS::ExtractMyValues(*curvaturenp_node_col, local_curvature, la[0].lm_);

    if(local_curvature.size() != numnode)
      dserror("wrong size of (potentially resized) local matrix!");

    std::copy(local_curvature.begin(), local_curvature.begin()+numnode, curvature.A());
  }

  if(require_smoothedgradphi_)
  {
    std::vector<double> local_normal(nsd_*numnode);
    DRT::UTILS::ExtractMyValues(*gradphinp_smoothed_node_col, local_normal, la[0].lm_);

    if(local_normal.size() != nsd_*numnode)
      dserror("wrong size of (potentially resized) local matrix!");

    // copy local to normal....
     std::copy(local_normal.begin(), local_normal.begin()+(nsd_*numnode), normal.A());
  }
}


/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingTwoPhase::Output(
    const int step,
    const double time,
    const bool write_restart_data,
    const int lsc_idx
)
{
  LevelSetCoupling::Output(step, time, write_restart_data, lsc_idx);

  std::ostringstream temp;
  temp << lsc_idx;

  // write restart
  if (write_restart_data)
  {
    std::string name_restart;

    if(require_nodalcurvature_)
    {
      name_restart = "cutter_curv_res_"+temp.str();
      cutter_output_->WriteVector(name_restart, curvaturenp_node_);
    }

    if(require_smoothedgradphi_)
    {
      name_restart = "cutter_gradphi_res_"+temp.str();
      cutter_output_->WriteVector(name_restart, gradphinp_smoothed_node_);
    }
  }
}


/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingTwoPhase::ReadRestart(
    const int step,
    const int lsc_idx
)
{
  //-------- boundary discretization
  IO::DiscretizationReader boundaryreader(cutter_dis_, step);

  const double time = boundaryreader.ReadDouble("time");

  if(myrank_ == 0)
  {
    IO::cout << "           RESTART IS PERFORMED FROM STORED VALUES!                            " << IO::endl;
    IO::cout << "ReadRestart for Level Set Cut in Xfluid (time="<< time <<" ; step="<< step <<")" << IO::endl;
  }

  std::ostringstream temp;
  temp << lsc_idx;

  std::string name_restart;

  {
    name_restart= "cutter_phinp_res_"+temp.str();

    boundaryreader.ReadVector(cutter_phinp_,   name_restart);

    if (not (cutter_dis_->DofRowMap(cutter_nds_phi_))->SameAs(cutter_phinp_->Map()))
      dserror("Global dof numbering in maps does not match");
  }

  if(require_nodalcurvature_)
  {
     name_restart = "cutter_curv_res_"+temp.str();

     boundaryreader.ReadVector(curvaturenp_node_,   name_restart);

     if (not (cutter_dis_->DofRowMap(cutter_nds_phi_))->SameAs(curvaturenp_node_->Map()))
       dserror("Global dof numbering in maps does not match");
  }

  if(require_smoothedgradphi_)
  {
    std::string name_restart = "cutter_gradphi_res_"+temp.str();

    boundaryreader.ReadMultiVector(gradphinp_smoothed_node_,   name_restart);

    if (not (cutter_dis_->DofRowMap(cutter_nds_phi_))->SameAs(gradphinp_smoothed_node_->Map()))
      dserror("Global dof numbering in maps does not match");
  }

  ExportGeometricQuantities();
}


/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
XFEM::LevelSetCouplingCombustion::LevelSetCouplingCombustion(
    Teuchos::RCP<DRT::Discretization>&  bg_dis,    ///< background discretization
    const std::string &                 cond_name, ///< name of the condition, by which the derived cutter discretization is identified
    Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< discretization from which the cutter discretization is derived
    const int                           coupling_id,///< id of composite of coupling conditions
    const double                        time,      ///< time
    const int                           step       ///< time step
) : LevelSetCouplingTwoPhase(bg_dis, cond_name, cond_dis, coupling_id, time, step),
    laminar_flamespeed_(0.0),
    mol_diffusivity_(0.0),
    markstein_length_(0.0)
{
}


/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingCombustion::SetParameters_Physical(
    DRT::Condition* cond
)
{
  // get the laminar flame speed (sl)
   laminar_flamespeed_ = cond->GetDouble("laminar_flamespeed");

   // get the molecular diffusivity
   mol_diffusivity_    = cond->GetDouble("mol_diffusivity");

   // get the markstein length
   markstein_length_   = cond->GetDouble("markstein_length");
}


/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingCombustion::SetFlags_InterfaceTransport(
    DRT::Condition* cond
)
{
  // use just the interface-normal part of the convective velocity for computing flame transport velocity?
  transport_direction_ = static_cast<INPAR::TWOPHASE::Transport_Directions>(cond->GetInt("TRANSPORT_DIRECTIONS"));

  // account for curvature in the transport velocity
  transport_curvature_ = (bool)cond->GetInt("transport_curvature");

  if(transport_curvature_)
    dserror("curvature-driven transport not tested so far!");
}


/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingCombustion::SetFlags_GeometricQuantities()
{
  // set the flags as in two-phase flow
  LevelSetCouplingTwoPhase::SetFlags_GeometricQuantities();

  // manipulate it afterwards for combustion two-phase flow
  require_smoothedgradphi_ = true;  // we always use smoothed normals

  // if we do not evaluate it, then we can switch the computation of the curvature off
  // -> otherwise let it as set in the two-phase section
  if(!transport_curvature_)
    require_nodalcurvature_ = false;
}



/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingCombustion::GetPhaseDensities(
    DRT::Element* ele,
    double & rhoplus,
    double & rhominus
    )
{
  // density burnt (+,j) domain (j>i)
  Teuchos::RCP<MAT::Material> matptrplus, matptrminus = Teuchos::null;

  XFEM::UTILS::GetVolumeCellMaterial(ele, matptrplus,GEO::CUT::Point::inside);
  dsassert(matptrplus->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
  const MAT::NewtonianFluid* matplus = static_cast<const MAT::NewtonianFluid*>(matptrplus.get());
  rhoplus = matplus->Density();

  // density unburnt (-,i) domain (j>i) domain
  XFEM::UTILS::GetVolumeCellMaterial(ele, matptrminus,GEO::CUT::Point::outside);
  dsassert(matptrminus->MaterialType() == INPAR::MAT::m_fluid, "material is not of type m_fluid");
  const MAT::NewtonianFluid* matminus = static_cast<const MAT::NewtonianFluid*>(matptrminus.get());
  rhominus = matminus->Density();
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
void XFEM::LevelSetCouplingCombustion::ComputeRelativeTransportVelocity(
    Epetra_SerialDenseVector & flvelrel,
    const GEO::CUT::Point::PointPosition & position,
    const Epetra_SerialDenseVector & nvec,
    const double & curv
)
{
  //------------------------
  // get fluid material parameters
  //------------------------
  // get list of adjacent elements of this node
  if(bg_dis_->NumMyRowElements() == 0)
    dserror("no bg_dis_ row element available to get the fluid material");

  DRT::Element* ele = bg_dis_->lRowElement(0);

  // get material from first (arbitrary!) element adjacent to this node
  double rhoplus  = 0.0; // density burnt   (+,j) domain (j>i) subdomain
  double rhominus = 0.0; // density unburnt (-,i) domain (j>i) subdomain

  GetPhaseDensities(ele, rhoplus, rhominus);

  //---------------------------------------------
  // compute relative flame velocity at the scatra'node position
  //---------------------------------------------

  // compute the speed factor depending on the positioning of the node
  const double wallfac = 1.0; //TODO:

  double speedfac = laminar_flamespeed_*(1.0-markstein_length_*curv);

  if (position == GEO::CUT::Point::inside)
  {
    // interface or burnt domain -> burnt material
    // flame speed factor = laminar flame speed * rho_unburnt / rho_burnt
    speedfac *= rhominus/rhoplus;
  }
  else if (position == GEO::CUT::Point::outside)
  {
    // unburnt domain -> unburnt material
    // flame speed factor = laminar flame speed
  }
  else dserror("what to do now? interface is on a node!");

  //-----------------------------------------------
  // compute the relative flame velocity at this node
  //-----------------------------------------------
  flvelrel=nvec;
  flvelrel.Scale(-wallfac * speedfac);

  if(transport_curvature_)
    dserror("how to account for curvature here?");
}


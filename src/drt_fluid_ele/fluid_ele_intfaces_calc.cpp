/*!----------------------------------------------------------------------
\file fluid_ele_intfaces_calc.cpp
\brief

Integrate internal faces (stabilization) terms on an internal faces element

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

*----------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>

#include "fluid_ele_action.H"
#include "fluid_ele_intfaces_calc.H"
#include "fluid_ele_parameter_std.H"
#include "fluid_ele_parameter_timint.H"
#include "fluid_ele_calc_intfaces_stab.H"

#include "../drt_inpar/inpar_xfem.H"

#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_discret_faces.H"



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidIntFaceImplInterface* DRT::ELEMENTS::FluidIntFaceImplInterface::Impl(const DRT::Element* ele)
{
  switch (ele->Shape())
  {
  case DRT::Element::quad4:
  {
    return FluidIntFaceImpl<DRT::Element::quad4>::Instance();
  }
  case DRT::Element::quad8:
  {
    return FluidIntFaceImpl<DRT::Element::quad8>::Instance();
  }
  case DRT::Element::quad9:
  {
    return FluidIntFaceImpl<DRT::Element::quad9>::Instance();
  }
  case DRT::Element::tri3:
  {
    return FluidIntFaceImpl<DRT::Element::tri3>::Instance();
  }
  case DRT::Element::tri6:
  {
    return FluidIntFaceImpl<DRT::Element::tri6>::Instance();
  }
  case DRT::Element::line2:
  {
    return FluidIntFaceImpl<DRT::Element::line2>::Instance();
  }
  case DRT::Element::line3:
  {
    return FluidIntFaceImpl<DRT::Element::line3>::Instance();
  }
  default:
    dserror("Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode()); break;
  }
  return NULL;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidIntFaceImpl<distype> * DRT::ELEMENTS::FluidIntFaceImpl<distype>::Instance(bool create)
{
  static FluidIntFaceImpl<distype> * instance;
  if (create)
  {
    if (instance==NULL)
      instance = new FluidIntFaceImpl<distype>();
  }
  else
  {
    if (instance!=NULL)
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidIntFaceImpl<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );

}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidIntFaceImpl<distype>::FluidIntFaceImpl()
{
  // pointer to class FluidImplParameterTimInt (access to the time-integration parameter)
  fldparatimint_ = DRT::ELEMENTS::FluidEleParameterTimInt::Instance();
  // pointer to class FluidImplParameter (access to the general parameter)
  fldpara_ = DRT::ELEMENTS::FluidEleParameterStd::Instance();

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::FluidIntFaceImpl<distype>::PrepareAssemble(
    Teuchos::ParameterList &   stabparams,
    Teuchos::ParameterList &   faceparams)
{
  // decide which terms have to be assembled and decide the assembly pattern

  INPAR::XFEM::FaceType face_type = faceparams.get<INPAR::XFEM::FaceType>("facetype");

  // final decision which terms are assembled for the current face
  bool EOS_Pres         = false;
  bool EOS_Conv_Stream  = false;
  bool EOS_Conv_Cross   = false;
  bool EOS_Div_vel_jump = false;
  bool EOS_Div_div_jump = false;
  bool GP_visc          = false;
  bool GP_u_p_2nd       = false;
  bool GP_trans         = false;


  if(face_type == INPAR::XFEM::face_type_std)
  {
    EOS_Pres         = (DRT::ELEMENTS::FluidEleParameterStd::Instance()->EOS_Pres()        == INPAR::FLUID::EOS_PRES_std_eos);
    EOS_Conv_Stream  = (DRT::ELEMENTS::FluidEleParameterStd::Instance()->EOS_Conv_Stream() == INPAR::FLUID::EOS_CONV_STREAM_std_eos);
    EOS_Conv_Cross   = (DRT::ELEMENTS::FluidEleParameterStd::Instance()->EOS_Conv_Cross()  == INPAR::FLUID::EOS_CONV_CROSS_std_eos);
    EOS_Div_vel_jump = (DRT::ELEMENTS::FluidEleParameterStd::Instance()->EOS_Div()         == INPAR::FLUID::EOS_DIV_vel_jump_std_eos);
    EOS_Div_div_jump = (DRT::ELEMENTS::FluidEleParameterStd::Instance()->EOS_Div()         == INPAR::FLUID::EOS_DIV_div_jump_std_eos);

    GP_visc          = false;
    GP_trans         = false;
    GP_u_p_2nd       = false;
  }
  else if(face_type == INPAR::XFEM::face_type_ghost_penalty)
  {
    EOS_Pres         = (DRT::ELEMENTS::FluidEleParameterStd::Instance()->EOS_Pres()        != INPAR::FLUID::EOS_PRES_none);
    EOS_Conv_Stream  = (DRT::ELEMENTS::FluidEleParameterStd::Instance()->EOS_Conv_Stream() != INPAR::FLUID::EOS_CONV_STREAM_none);
    EOS_Conv_Cross   = (DRT::ELEMENTS::FluidEleParameterStd::Instance()->EOS_Conv_Cross()  != INPAR::FLUID::EOS_CONV_CROSS_none);
    EOS_Div_vel_jump = (DRT::ELEMENTS::FluidEleParameterStd::Instance()->EOS_Div()         == INPAR::FLUID::EOS_DIV_vel_jump_std_eos
                     or DRT::ELEMENTS::FluidEleParameterStd::Instance()->EOS_Div()         == INPAR::FLUID::EOS_DIV_vel_jump_xfem_gp);
    EOS_Div_div_jump = (DRT::ELEMENTS::FluidEleParameterStd::Instance()->EOS_Div()         == INPAR::FLUID::EOS_DIV_div_jump_std_eos
                     or DRT::ELEMENTS::FluidEleParameterStd::Instance()->EOS_Div()         == INPAR::FLUID::EOS_DIV_div_jump_xfem_gp);

    GP_visc          = faceparams.get<bool>("visc_ghost_penalty", false);
    GP_trans         = faceparams.get<bool>("trans_ghost_penalty", false);
    GP_u_p_2nd       = faceparams.get<bool>("u_p_ghost_penalty_2nd", false);
  }
  else if(face_type == INPAR::XFEM::face_type_ghost)
  {
    EOS_Pres         = false;
    EOS_Conv_Stream  = false;
    EOS_Conv_Cross   = false;
    EOS_Div_vel_jump = false;
    EOS_Div_div_jump = false;
    GP_visc          = false;
    GP_trans         = false;
    GP_u_p_2nd       = false;
  }
  else dserror("unknown face_type!!!");

  // which pattern has to be activated?
  // TODO: this can be improved if only pressure is assembled and so on!

  if(EOS_Div_div_jump)
  {
    stabparams.set("eos_gp_pattern",INPAR::FLUID::EOS_GP_Pattern_up);
  }
  else
  {
    stabparams.set("eos_gp_pattern",INPAR::FLUID::EOS_GP_Pattern_uvwp);
  }

  stabparams.set<bool>("EOS_Pres",         EOS_Pres);
  stabparams.set<bool>("EOS_Conv_Stream",  EOS_Conv_Stream);
  stabparams.set<bool>("EOS_Conv_Cross",   EOS_Conv_Cross);
  stabparams.set<bool>("EOS_Div_vel_jump", EOS_Div_vel_jump);
  stabparams.set<bool>("EOS_Div_div_jump", EOS_Div_div_jump);
  stabparams.set<bool>("GP_visc",          GP_visc);
  stabparams.set<bool>("GP_trans",         GP_trans);
  stabparams.set<bool>("GP_u_p_2nd",       GP_u_p_2nd);


  stabparams.set("ghost_penalty_reconstruct", faceparams.get<bool>("ghost_penalty_reconstruct", false) );
  stabparams.set("ghost_penalty_fac",         faceparams.get<double>("GHOST_PENALTY_FAC", 0.0));
  stabparams.set("ghost_penalty_trans_fac",   faceparams.get<double>("GHOST_PENALTY_TRANSIENT_FAC", 0.0));


  stabparams.set("action", faceparams.get<int>("action"));

  // return false if no stabilization is required
  if( !EOS_Pres and
      !EOS_Conv_Stream and
      !EOS_Conv_Cross and
      !EOS_Div_vel_jump and
      !EOS_Div_div_jump and
      !GP_visc and
      !GP_trans and
      !GP_u_p_2nd) return false;

  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidIntFaceImpl<distype>::AssembleInternalFacesUsingNeighborData(
    DRT::ELEMENTS::FluidIntFace*         intface,         ///< internal face element
    std::vector<int>&                    nds_master,      ///< nodal dofset w.r.t. master element
    std::vector<int>&                    nds_slave,       ///< nodal dofset w.r.t. slave element
    Teuchos::ParameterList&              params,          ///< parameter list
    DRT::DiscretizationFaces&             discretization,  ///< faces discretization
    Teuchos::RCP<LINALG::SparseMatrix>            systemmatrix,    ///< systemmatrix
    Teuchos::RCP<Epetra_Vector>                   systemvector     ///< systemvector
    )
{
  Teuchos::ParameterList edgebasedparams;

  //decide which terms have to be assembled and decide the assembly pattern, return if no assembly required
  bool stab_required = PrepareAssemble(edgebasedparams, params);

  // do not assemble if no stabilization terms activated for this face
  if(!stab_required) return;

  if (!discretization.Filled()) dserror("FillComplete() was not called");
  if (!discretization.HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  const bool assemblemat = systemmatrix!=Teuchos::null;
  const bool assemblevec = systemvector!=Teuchos::null;

  //--------------------------------------------------------
  /// number of space dimensions of the FluidIntFace element
  static const int facensd = DRT::UTILS::DisTypeToDim<distype>::dim;

  /// number of space dimensions of the parent element
  static const int nsd = facensd+1;

  /// number of dof's per node
  static const int numdofpernode = nsd + 1;


  //----------------------- create patchlm -----------------

  // local maps for patch dofs
  std::vector<int> lm_patch;

  // local maps for master/slave/face dofs
  std::vector<int> lm_master;
  std::vector<int> lm_slave;
  std::vector<int> lm_face;

  // local maps between master/slave dofs and position in patch dofs (lm_patch)
  std::vector<int> lm_masterToPatch;
  std::vector<int> lm_slaveToPatch;
  std::vector<int> lm_faceToPatch;

  // local maps between master/slave nodes and position in patch nodes
  std::vector<int> lm_masterNodeToPatch;
  std::vector<int> lm_slaveNodeToPatch;

  // create patch location vector combining master element, slave element and face element
  intface->PatchLocationVector(   discretization,
                                  nds_master,nds_slave,
                                  lm_patch, lm_master, lm_slave, lm_face,
                                  lm_masterToPatch, lm_slaveToPatch, lm_faceToPatch,
                                  lm_masterNodeToPatch, lm_slaveNodeToPatch
                                  );


  // patch_lm for Velx,Vely,Velz, Pres field
  std::vector<std::vector<int> >patch_components_lm(numdofpernode);
  std::vector<std::vector<int> >patch_components_lmowner(numdofpernode);

  // modify the patch owner to the owner of the internal face element
  std::vector<int> patchlm_owner;
  int owner = intface->Owner();

  for(unsigned i=0; i<lm_patch.size(); i++)
  {
    // i%4 yields the Velx,Vely,Velz,Pres field
    patch_components_lm[(i%numdofpernode)].push_back(lm_patch[i]);
    patch_components_lmowner[(i%numdofpernode)].push_back(owner);

    patchlm_owner.push_back(owner);
  }

  int numnodeinpatch = patch_components_lm[0].size();

#ifdef DEBUG
  for(int isd=0; isd < numdofpernode; isd++)
  if((int)(patch_components_lm[isd].size()) != numnodeinpatch) dserror("patch_components_lm[%d] has wrong size", isd);

#endif

  //------------- create and evaluate block element matrics -----------------

  // define element matrices and vectors
  std::vector<Epetra_SerialDenseMatrix> elemat_blocks;
  std::vector<Epetra_SerialDenseVector> elevec_blocks;

  //------------------------------------------------------------------------------------
  // decide which pattern
  // pattern = "u-v-w-p-diagonal-block matrix pattern";    // assembles each component block separated
  // pattern = "u-p-block matrix pattern";                 // assembles u-block and p-block separated
  // pattern = "full matrix pattern";                      // assembles the whole u-p matrix

  INPAR::FLUID::EOS_GP_Pattern eos_gp_pattern = edgebasedparams.get<INPAR::FLUID::EOS_GP_Pattern>("eos_gp_pattern");


  int numblocks = 0;

  if(eos_gp_pattern == INPAR::FLUID::EOS_GP_Pattern_uvwp)
  {
    // 3D: 4 blocks =  u-u block, v-v block, w-w block and p-p block
    // 2D: 3 blocks =  u-u block, v-v block and p-p block
    numblocks=numdofpernode;
  }
  else if(eos_gp_pattern == INPAR::FLUID::EOS_GP_Pattern_up)
  {
    // 3D: 10 blocks = 3x3 u-u blocks + 1x1 p-p block
    // 3D: 5 blocks  = 2x2 u-u blocks + 1x1 p-p block
    numblocks=nsd*nsd+1; // 10 blocks = 3x3 u-u blocks + 1x1 p-p block
  }
  else if(eos_gp_pattern == INPAR::FLUID::EOS_GP_Pattern_full)
  {
    // 3D: 16 blocks = 4x4 uvwp blocks
    // 2D: 9  blocks = 3x3 uvp blocks
    numblocks=numdofpernode*numdofpernode;
  }

  else dserror("unknown matrix pattern");


  elemat_blocks.resize(numblocks);

  for(int b=0; b<numblocks; b++)
  {
    int err = elemat_blocks[b].Shape(numnodeinpatch,numnodeinpatch); // new shape and init values to zero

    if(err != 0) dserror("element matrix Shape not successful");
  }

  elevec_blocks.resize(numdofpernode); // 3D: 4 vectors for u,v,w,p components, 2D: 3 vectors for u,v,p

  for(int b=0; b<numdofpernode; b++)
  {
    int err = elevec_blocks[b].Size(numnodeinpatch); // new size and init values to zero

    if(err != 0) dserror("element matrix Shape not successful");
  }


  //---------------------------------------------------------------------
  // call the element specific evaluate method

  int err = EvaluateInternalFaces( intface, edgebasedparams, discretization,
                                   lm_patch,
                                   lm_masterToPatch,lm_slaveToPatch,lm_faceToPatch,
                                   lm_masterNodeToPatch, lm_slaveNodeToPatch,
                                   elemat_blocks, elevec_blocks
                                  );
  if (err) dserror("error while evaluating elements");

  //---------------------------------------------------------------------------------------
  // assemble systemmatrix
  {
    TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: assemble FE matrix" );

    // calls the Assemble function for EpetraFECrs matrices including communication of non-row entries
    if (assemblemat)
    {
      if(eos_gp_pattern == INPAR::FLUID::EOS_GP_Pattern_uvwp)
      {
#if(1)
        for(int ij=0; ij<numdofpernode; ij++)
        {
          systemmatrix->FEAssemble(-1, elemat_blocks[ij], patch_components_lm[ij], patch_components_lmowner[ij], patch_components_lm[ij]);
        }
#else // assemble only pressure block
        int ij=nsd;
          systemmatrix->FEAssemble(-1, elemat_blocks[ij], patch_components_lm[ij], patch_components_lmowner[ij], patch_components_lm[ij]);
#endif
      }
      else if(eos_gp_pattern == INPAR::FLUID::EOS_GP_Pattern_up)
      {
        for(int i=0; i<nsd; i++)
        {
          for(int j=0; j<nsd; j++)
          {
            systemmatrix->FEAssemble(-1, elemat_blocks[i*nsd+j], patch_components_lm[i], patch_components_lmowner[i], patch_components_lm[j]);
          }
        }
        systemmatrix->FEAssemble(-1, elemat_blocks[nsd*nsd], patch_components_lm[nsd], patch_components_lmowner[nsd], patch_components_lm[nsd]);
      }
      else if(eos_gp_pattern == INPAR::FLUID::EOS_GP_Pattern_full)
      {
        for(int i=0; i<numdofpernode; i++)
        {
          for(int j=0; j<numdofpernode; j++)
          {
            systemmatrix->FEAssemble(-1, elemat_blocks[i*numdofpernode+j], patch_components_lm[i], patch_components_lmowner[i], patch_components_lm[j]);
          }
        }
      }
      else dserror("unknown matrix pattern");


    }
  }
  //---------------------------------------------------------------------------------------
  //assemble systemvector
  if (assemblevec)
  {
    TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: assemble vector" );

    // REMARK:: call Assemble without lmowner
    // to assemble the residual_col vector on only row elements also column nodes have to be assembled
    // do not exclude non-row nodes (modify the real owner to myowner)
    // after assembly the col vector it has to be exported to the row residual_ vector
    // using the 'Add' flag to get the right value for shared nodes
    for(int i=0; i<numdofpernode; i++)
    {
      LINALG::Assemble(*systemvector,elevec_blocks[i],patch_components_lm[i],patch_components_lmowner[i]);
    }

  }

  return;
}



/*----------------------------------------------------------------------*
 |  Evaluate internal faces (public)                        schott 01/12|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidIntFaceImpl<distype>::EvaluateInternalFaces(   DRT::ELEMENTS::FluidIntFace*       intface,              ///< internal face element
                                                                       Teuchos::ParameterList&            params,               ///< parameter list
                                                                       DRT::Discretization&               discretization,       ///< discretization
                                                                       std::vector<int>&                  patchlm,              ///< patch local map
                                                                       std::vector<int>&                  lm_masterToPatch,     ///< local map between master dofs and patchlm
                                                                       std::vector<int>&                  lm_slaveToPatch,      ///< local map between slave dofs and patchlm
                                                                       std::vector<int>&                  lm_faceToPatch,       ///< local map between face dofs and patchlm
                                                                       std::vector<int>&                  lm_masterNodeToPatch, ///< local map between master nodes and nodes in patch
                                                                       std::vector<int>&                  lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
                                                                       std::vector<Epetra_SerialDenseMatrix>&  elemat_blocks,        ///< element matrix blocks
                                                                       std::vector<Epetra_SerialDenseVector>&  elevec_blocks         ///< element vector blocks
                                                                       )
{
  FLD::IntFaceAction act = FLD::ifa_none;
  act = DRT::INPUT::get<FLD::IntFaceAction>(params,"action");

  switch(act)
  {
  case FLD::EOS_and_GhostPenalty_stabilization:
  {
    return DRT::ELEMENTS::FluidIntFaceStab::Impl(intface)->EvaluateEdgeBasedStabilization(
      intface,
      *fldparatimint_,
      *fldpara_,
      params,
      discretization,
      patchlm,
      lm_masterToPatch,
      lm_slaveToPatch,
      lm_faceToPatch,
      lm_masterNodeToPatch,
      lm_slaveNodeToPatch,
      elemat_blocks, elevec_blocks
    );
    break;
  }
  default:
      dserror("Unknown type of action for FluidIntFace"); break;
  } // end of switch(act)

  return 0;
}


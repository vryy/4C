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


#include "fluid_ele_intfaces_calc.H"
#include "fluid_ele_calc_intfaces_stab.H"
#include "fluid_ele_action.H"
#include "../drt_inpar/inpar_xfem.H"

#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_discret_xfem.H"



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
//  case DRT::Element::nurbs2:    // 1D nurbs boundary element
//  {
//    return FluidIntFaceImpl<DRT::Element::nurbs2>::Instance();
//  }
//  case DRT::Element::nurbs3:    // 1D nurbs boundary element
//  {
//    return FluidIntFaceImpl<DRT::Element::nurbs3>::Instance();
//  }
//  case DRT::Element::nurbs4:    // 2D nurbs boundary element
//  {
//    return FluidIntFaceImpl<DRT::Element::nurbs4>::Instance();
//  }
//  case DRT::Element::nurbs9:    // 2D nurbs boundary element
//  {
//    return FluidIntFaceImpl<DRT::Element::nurbs9>::Instance();
//  }
  default:
    dserror("Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
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
  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidIntFaceImpl<distype>::AssembleInternalFacesUsingNeighborData(
    DRT::ELEMENTS::FluidIntFace*         intface,         ///< internal face element
    std::vector<int>&                    nds_master,      ///< nodal dofset w.r.t. master element
    std::vector<int>&                    nds_slave,       ///< nodal dofset w.r.t. slave element
    ParameterList&                       params,          ///< parameter list
    DRT::DiscretizationXFEM&             discretization,  ///< XFEM discretization
    RCP<LINALG::SparseMatrix>            systemmatrix,    ///< systemmatrix
    RCP<Epetra_Vector>                   systemvector     ///< systemvector
    )
{

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
  vector<vector<int> >patch_components_lm(numdofpernode);
  vector<vector<int> >patch_components_lmowner(numdofpernode);

  // modify the patch owner to the owner of the internal face element
  vector<int> patchlm_owner;
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
  vector<Epetra_SerialDenseMatrix> elemat_blocks;
  vector<Epetra_SerialDenseVector> elevec_blocks;

  //------------------------------------------------------------------------------------
  // decide which pattern
  // pattern = "u-v-w-p-diagonal-block matrix pattern";    // assembles each component block separated
  // pattern = "u-p-block matrix pattern";                 // assembles u-block and p-block separated
  // pattern = "full matrix pattern";                      // assembles the whole u-p matrix

  INPAR::XFEM::EOS_GP_Pattern eos_gp_pattern = params.get<INPAR::XFEM::EOS_GP_Pattern>("eos_gp_pattern");


  int numblocks = 0;

  if(eos_gp_pattern == INPAR::XFEM::EOS_GP_Pattern_uvwp)
  {
    // 3D: 4 blocks =  u-u block, v-v block, w-w block and p-p block
    // 2D: 3 blocks =  u-u block, v-v block and p-p block
    numblocks=numdofpernode;
  }
  else if(eos_gp_pattern == INPAR::XFEM::EOS_GP_Pattern_up)
  {
    // 3D: 10 blocks = 3x3 u-u blocks + 1x1 p-p block
    // 3D: 5 blocks  = 2x2 u-u blocks + 1x1 p-p block
    numblocks=nsd*nsd+1; // 10 blocks = 3x3 u-u blocks + 1x1 p-p block
  }
  else if(eos_gp_pattern == INPAR::XFEM::EOS_GP_Pattern_full)
  {
    // 3D: 16 blocks = 4x4 uvwp blocks
    // 2D: 9  blocks = 3x3 uvp blocks
    numblocks=numdofpernode*numdofpernode;
  }

  else dserror("unknown matrix pattern");


  elemat_blocks.resize(numblocks);

  for(int b=0; b<numblocks; b++)
  {
    elemat_blocks[b].Shape(numnodeinpatch,numnodeinpatch); // new shape and init values to zero
  }

  elevec_blocks.resize(numdofpernode); // 3D: 4 vectors for u,v,w,p components, 2D: 3 vectors for u,v,p

  for(int b=0; b<numdofpernode; b++)
  {
    elevec_blocks[b].Size(numnodeinpatch); // new size and init values to zero
  }


  //---------------------------------------------------------------------
  // call the element specific evaluate method

  int err = EvaluateInternalFaces( intface,params,discretization,
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
      if(eos_gp_pattern == INPAR::XFEM::EOS_GP_Pattern_uvwp)
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
      else if(eos_gp_pattern == INPAR::XFEM::EOS_GP_Pattern_up)
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
      else if(eos_gp_pattern == INPAR::XFEM::EOS_GP_Pattern_full)
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
                                                                       vector<int>&                       patchlm,              ///< patch local map
                                                                       vector<int>&                       lm_masterToPatch,     ///< local map between master dofs and patchlm
                                                                       vector<int>&                       lm_slaveToPatch,      ///< local map between slave dofs and patchlm
                                                                       vector<int>&                       lm_faceToPatch,       ///< local map between face dofs and patchlm
                                                                       std::vector<int>&                  lm_masterNodeToPatch, ///< local map between master nodes and nodes in patch
                                                                       std::vector<int>&                  lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
                                                                       vector<Epetra_SerialDenseMatrix>&  elemat_blocks,        ///< element matrix blocks
                                                                       vector<Epetra_SerialDenseVector>&  elevec_blocks         ///< element vector blocks
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
      dserror("Unknown type of action for FluidIntFace");
  } // end of switch(act)

  return 0;
}


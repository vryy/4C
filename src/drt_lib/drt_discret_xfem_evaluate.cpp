/*!----------------------------------------------------------------------
\file drt_discret_xfem_evaluate.cpp

\brief a class to manage one discretization

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>


#include "drt_discret_xfem.H"
#include "drt_globalproblem.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_combust/combust_defines.H"
#include "../drt_combust/combust3.H"

#include "drt_exporter.H"

#include "drt_utils.H"
#include "../linalg/linalg_utils.H"

#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid_ele/fluid_ele_intfaces_calc.H"
#include "../drt_fluid_ele/fluid_ele_action.H"

#include "../drt_inpar/inpar_xfem.H"

/*----------------------------------------------------------------------*
 |  Evaluate edge-based integrals (public)               rasthofer 02/13|
 *----------------------------------------------------------------------*/
void DRT::DiscretizationXFEM::EvaluateEdgeBasedCombust(
        Teuchos::ParameterList&              params,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
        Teuchos::RCP<Epetra_Vector>          systemvector1
)
{

  std::cout << "->  EvaluateEdgeBased()" << std::endl;

  TEUCHOS_FUNC_TIME_MONITOR( "DRT::DiscretizationXFEM::EdgeBased" );


  Teuchos::RCP<Epetra_Vector> residual_col = LINALG::CreateVector(*(this->DofColMap()),true);

  const Epetra_Map* rmap = NULL;
//  const Epetra_Map* dmap = NULL;

  Teuchos::RCP<Epetra_FECrsMatrix> sysmat_FE;
  if (systemmatrix1 != Teuchos::null)
  {
    rmap = &(systemmatrix1->OperatorRangeMap());
//    dmap = rmap;
    sysmat_FE = Teuchos::rcp(new Epetra_FECrsMatrix(::Copy,*rmap,256,false));
  }
  else dserror("sysmat is NULL!");

  Teuchos::RCP<LINALG::SparseMatrix> sysmat_linalg = Teuchos::rcp(new LINALG::SparseMatrix(sysmat_FE,true,false,LINALG::SparseMatrix::FE_MATRIX));

  const int numrowintfaces = NumMyRowIntFaces();

  for (int i=0; i<numrowintfaces; ++i)
  {
    DRT::Element* actface = lRowIntFace(i);

    DRT::ELEMENTS::Combust3IntFace * ele = dynamic_cast<DRT::ELEMENTS::Combust3IntFace *>(actface);
    if ( ele==NULL ) dserror( "expect Combust3IntFace element" );

    // get the parent combust elements
    DRT::ELEMENTS::Combust3* p_master = ele->ParentMasterElement();
    DRT::ELEMENTS::Combust3* p_slave  = ele->ParentSlaveElement();

    // TODO: num node
//    std::cout << "Hier!" << std::endl;
    size_t p_master_numnode = p_master->NumNode();
    size_t p_slave_numnode  = p_slave->NumNode();

    std::vector<int> nds_master;
    nds_master.reserve(p_master_numnode);

    std::vector<int> nds_slave;
    nds_slave.reserve(p_slave_numnode);

    {
      TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: create nds" );

      for(size_t i=0; i< p_master_numnode; i++)  nds_master.push_back(0);

      for(size_t i=0; i< p_slave_numnode; i++)   nds_slave.push_back(0);
    }

    // call the internal faces stabilization routine for the current side/surface
    TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: AssembleEdgeStabGhostPenalty" );


    // call edge-based stabilization and ghost penalty
    Teuchos::ParameterList edgebasedparams;

    // set action for elements
//    edgebasedparams.set<int>("action",FLD::EOS_and_GhostPenalty_stabilization);
//    edgebasedparams.set("facetype", INPAR::XFEM::face_type_std);

//    // call the egde-based assemble and evaluate routine
    EvaluateEdgeBasedCombust(ele,
                             nds_master,
                             nds_slave,
                             params,
                             *this,
                             sysmat_linalg,
                             residual_col);
  }

  sysmat_linalg->Complete();

  (systemmatrix1)->Add(*sysmat_linalg, false, 1.0, 1.0);

  //------------------------------------------------------------
  // need to export residual_col to systemvector1 (residual_)
  Epetra_Vector res_tmp(systemvector1->Map(),false);
  Epetra_Export exporter(residual_col->Map(),res_tmp.Map());
  int err2 = res_tmp.Export(*residual_col,exporter,Add);
  if (err2) dserror("Export using exporter returned err=%d",err2);
  systemvector1->Update(1.0,res_tmp,1.0);

  return;
}


/*----------------------------------------------------------------------*
 |  Evaluate edge-based integrals (public)               rasthofer 02/13|
 *----------------------------------------------------------------------*/
void DRT::DiscretizationXFEM::EvaluateEdgeBasedCombust(
        DRT::ELEMENTS::Combust3IntFace*      ele,             ///< internal face element
        std::vector<int>&                    nds_master,      ///< nodal dofset w.r.t. master element
        std::vector<int>&                    nds_slave,       ///< nodal dofset w.r.t. slave element
        Teuchos::ParameterList&              params,          ///< parameter list
        DRT::DiscretizationXFEM&             discretization,  ///< XFEM discretization
        RCP<LINALG::SparseMatrix>            systemmatrix,    ///< systemmatrix
        RCP<Epetra_Vector>                   systemvector     ///< systemvector
)
{
    //TODO: denke ich brauche das nicht, weil das bereits alles in der Parameterliste stehen sollte
    //Teuchos::ParameterList edgebasedparams;
    //decide which terms have to be assembled and decide the assembly pattern, return if no assembly required
    //bool stab_required = PrepareAssemble(edgebasedparams, params);

    // do not assemble if no stabilization terms activated for this face
    // TODO: was mach ich damit
    //if(!stab_required) return;

    if (!discretization.Filled()) dserror("FillComplete() was not called");
    if (!discretization.HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

    const bool assemblemat = systemmatrix!=Teuchos::null;
    const bool assemblevec = systemvector!=Teuchos::null;

    //--------------------------------------------------------
    /// number of space dimensions of the Combust3IntFace element
    if (ele->Shape()!=DRT::Element::quad4) dserror("Quad4 expected!");
    static const int facensd = DRT::UTILS::DisTypeToDim<DRT::Element::quad4>::dim;

    /// number of space dimensions of the parent element
    static const int nsd = facensd+1;

    /// number of dof's per node
    // TODO: Vorsicht xfem
    std::cout << "das geht so hier wohl nicht" << std::endl;
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
    ele->PatchLocationVector(discretization,
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
    int owner = ele->Owner();

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

    // the required blocks depends on the chosen edge-based terms
    // pressure stabilization only contributes to the pressure-pressure block
    // convective stabilization only to the velx-velx, vely-vely and velz-velz block
    // conti stabilization contributes to all velocity blocks, but not to the pressure block
    // see also combust3_facemat_standard.H
    INPAR::FLUID::EOS_GP_Pattern eos_gp_pattern = INPAR::FLUID::EOS_GP_Pattern_uvwp;
    INPAR::FLUID::EOS_Div conti_stab = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_Div>(params.sublist("EDGE-BASED STABILIZATION"),"EOS_DIV");
    if (conti_stab == INPAR::FLUID::EOS_DIV_div_jump_std_eos or conti_stab == INPAR::FLUID::EOS_DIV_div_jump_xfem_gp)
        eos_gp_pattern = INPAR::FLUID::EOS_GP_Pattern_up;
    // set pattern in paramterlist
    params.sublist("EDGE-BASED STABILIZATION").set<INPAR::FLUID::EOS_GP_Pattern>("EOS_PATTERN",eos_gp_pattern);

    int numblocks = 0;

    if (eos_gp_pattern == INPAR::FLUID::EOS_GP_Pattern_uvwp)
    {
      // 3D: 4 blocks =  u-u block, v-v block, w-w block and p-p block
      numblocks=numdofpernode;
    }
    else if(eos_gp_pattern == INPAR::FLUID::EOS_GP_Pattern_up)
    {
      // 3D: 10 blocks = 3x3 u-u blocks + 1x1 p-p block
      numblocks=nsd*nsd+1;
    }
    else if(eos_gp_pattern == INPAR::FLUID::EOS_GP_Pattern_full)
    {
      // 3D: 16 blocks = 4x4 uvwp blocks
      dserror("Full pattern does not make sense currently");
      // for the time being no terms that require all blocks are implemented
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

    // TODO: dummy call
//    Epetra_SerialDenseMatrix elemat1;
//    Epetra_SerialDenseMatrix elemat2;
//    Epetra_SerialDenseVector elevec1;
//    Epetra_SerialDenseVector elevec2;
//    Epetra_SerialDenseVector elevec3;
//    int err = ele->Evaluate(params,*this,lm_patch,
//                            elemat1,elemat2,elevec1,elevec2,elevec3);
      int err = ele->Evaluate(params,*this, lm_patch, lm_face, lm_master, lm_slave,
                              lm_masterToPatch,lm_slaveToPatch,lm_faceToPatch,
                              lm_masterNodeToPatch, lm_slaveNodeToPatch,
                              elemat_blocks, elevec_blocks
                              );

//    int err = EvaluateInternalFaces( intface, edgebasedparams, discretization,
//                                     lm_patch,
//                                     lm_masterToPatch,lm_slaveToPatch,lm_faceToPatch,
//                                     lm_masterNodeToPatch, lm_slaveNodeToPatch,
//                                     elemat_blocks, elevec_blocks
//                                    );
    if (err) dserror("error while evaluating elements");

    //---------------------------------------------------------------------------------------
    // assemble systemmatrix
    {
      TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: assemble FE matrix" );

      // calls the Assemble function for EpetraFECrs matrices including communication of non-row entries
      if (assemblemat)
      {
        if (eos_gp_pattern == INPAR::FLUID::EOS_GP_Pattern_uvwp)
        {
          for(int ij=0; ij<numdofpernode; ij++)
          {
            systemmatrix->FEAssemble(-1, elemat_blocks[ij], patch_components_lm[ij], patch_components_lmowner[ij], patch_components_lm[ij]);
          }
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
          dserror("Full pattern does not make sense currently");
          // for the time being no terms that require all blocks are implemented
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

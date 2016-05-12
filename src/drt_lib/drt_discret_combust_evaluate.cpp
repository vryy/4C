/*!----------------------------------------------------------------------
\file drt_discret_combust_evaluate.cpp

\brief Implementation a class to manage one discretization

<pre>
\brief Implementation
\level 2
\maintainer Ursula Rasthofer
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>


#include "drt_discret_combust.H"
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

#include "../drt_xfem/dof_management_element.H"

#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_fluid.H"

/*----------------------------------------------------------------------*
 |  Evaluate edge-based integrals (public)               rasthofer 02/13|
 *----------------------------------------------------------------------*/
void DRT::DiscretizationCombust::EvaluateEdgeBasedCombust(
        Teuchos::ParameterList&              params,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
        Teuchos::RCP<Epetra_Vector>          systemvector1
)
{
  TEUCHOS_FUNC_TIME_MONITOR( "DRT::DiscretizationCombust::EdgeBased" );


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

  Teuchos::RCP<LINALG::SparseMatrix> sysmat_linalg = Teuchos::rcp(new LINALG::SparseMatrix(Teuchos::rcp_static_cast<Epetra_CrsMatrix>(sysmat_FE),LINALG::View,true,false,LINALG::SparseMatrix::FE_MATRIX));

  const int numrowintfaces = NumMyRowFaces();

  // get flags to switch between complete face-based stabilization, i.e., evaluation of all faces,
  // or evaluation of faces belonging to elements intersected by the interface, i.e., ghost penalty
  const bool xfemstab = params.get<bool>("xfemstab");
  bool edge_based_stab = false;
  if (params.sublist("RESIDUAL-BASED STABILIZATION").get<std::string>("STABTYPE")=="edge_based")
    edge_based_stab = true;

  for (int i=0; i<numrowintfaces; ++i)
  {
    DRT::Element* actface = lRowFace(i);

    DRT::ELEMENTS::Combust3IntFace * ele = dynamic_cast<DRT::ELEMENTS::Combust3IntFace *>(actface);
    if ( ele==NULL ) dserror( "expect Combust3IntFace element" );

    // get the parent combust elements
    DRT::ELEMENTS::Combust3* p_master = ele->ParentMasterElement();
    DRT::ELEMENTS::Combust3* p_slave  = ele->ParentSlaveElement();

    // get cut status of master and slave element
    const bool m_intersected = p_master->Bisected();
    const bool s_intersected = p_slave->Bisected();

//    if (p_master->Touched() or p_slave->Touched())
//      dserror("Touched element not yet considered.");

    // - if master is cut and slave is uncut or master is uncut and slave is cut,
    //   the face has to be considered in the xfem stabilization
    // - if master is cut and slave is cut, the face has to be integrated twice,
    //   with respect to the plus and the minus domain
    // - if master is uncut and slave is uncut, the face has only to be considered
    //   for a complete edge-based stabilization

    // check, if we have to evaluate this face
    if ((edge_based_stab == true) or (xfemstab == true and (m_intersected or s_intersected)))
    {
      // do we have to loop this face twice?
      int face_loop = 1;
      if (m_intersected and s_intersected)
        face_loop = 2;

      for (int i = 1; i <= face_loop; i++)
      {
        // call the internal faces stabilization routine for the current side/surface
        TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: AssembleEdgeStabGhostPenalty" );

        // set face type for evaluation with respect to plus and minus domain
        if (face_loop == 2)
        {
          if (i == 1)
            params.set<INPAR::XFEM::FaceType>("facetype", INPAR::XFEM::face_type_double_plus);
          else
            params.set<INPAR::XFEM::FaceType>("facetype", INPAR::XFEM::face_type_double_minus);
        }
        else
          params.set<INPAR::XFEM::FaceType>("facetype", INPAR::XFEM::face_type_std);

        // call the egde-based assemble and evaluate routine
        EvaluateEdgeBasedCombust(ele,
                                 params,
                                 *this,
                                 sysmat_linalg,
                                 residual_col);
      }
    }
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
void DRT::DiscretizationCombust::EvaluateEdgeBasedCombust(
        DRT::ELEMENTS::Combust3IntFace*      ele,             ///< internal face element
        Teuchos::ParameterList&              params,          ///< parameter list
        DRT::DiscretizationCombust&             discretization,  ///< XFEM discretization
        Teuchos::RCP<LINALG::SparseMatrix>            systemmatrix,    ///< systemmatrix
        Teuchos::RCP<Epetra_Vector>                   systemvector     ///< systemvector
)
{
    if (!discretization.Filled()) dserror("FillComplete() was not called");
    if (!discretization.HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

    const bool assemblemat = systemmatrix!=Teuchos::null;
    const bool assemblevec = systemvector!=Teuchos::null;

    /// safety check first
    if (ele->Shape()!=DRT::Element::quad4) dserror("Quad4 expected!");


    // get the parent combust elements
    DRT::ELEMENTS::Combust3* p_master = ele->ParentMasterElement();
    DRT::ELEMENTS::Combust3* p_slave  = ele->ParentSlaveElement();

    // get element dof manager
    Teuchos::RCP<XFEM::ElementDofManager> master_eledofmanager = p_master->GetEleDofManager();
    Teuchos::RCP<XFEM::ElementDofManager> slave_eledofmanager = p_slave->GetEleDofManager();

    // number of fields per node
    static const int numfieldpernode = master_eledofmanager->NumFields();
    if (numfieldpernode != slave_eledofmanager->NumFields())
      dserror("Same number of fields expected");

    //----------------------- create patchlm -----------------

    // local maps for master/slave/face dofs
    std::vector<int> lm_master;
    std::vector<int> lm_slave;

    // local maps between master/slave nodes and position in patch nodes
    std::map<XFEM::PHYSICS::Field,std::vector<int> > lm_masterDofPerFieldToPatch;
    std::map<XFEM::PHYSICS::Field,std::vector<int> > lm_slaveDofPerFieldToPatch;

    // patch_lm for Velx,Vely,Velz, Pres field
    std::map<XFEM::PHYSICS::Field,std::vector<int> > patch_components_lm;
    std::map<XFEM::PHYSICS::Field,std::vector<int> > patch_components_lmowner;

    // set all fields
    const std::vector<XFEM::PHYSICS::Field> fields = master_eledofmanager->GetFields();
    for (std::size_t ifield = 0; ifield < fields.size(); ifield++)
    {
        lm_masterDofPerFieldToPatch[fields[ifield]] = std::vector<int>();
        lm_slaveDofPerFieldToPatch[fields[ifield]] = std::vector<int>();
        patch_components_lm[fields[ifield]] = std::vector<int>();
        patch_components_lmowner[fields[ifield]] = std::vector<int>();
    }

    // create patch location vector combining master element, slave element and face element
    ele->PatchLocationVector(discretization,
                             lm_master, lm_slave,
                             lm_masterDofPerFieldToPatch, lm_slaveDofPerFieldToPatch,
                             patch_components_lm, patch_components_lmowner
                             );

    //------------- create and evaluate block element matrics -----------------

    // define element matrices and vectors
    std::map<std::pair<XFEM::PHYSICS::Field,XFEM::PHYSICS::Field>,Epetra_SerialDenseMatrix> elemat_blocks;
    std::map<XFEM::PHYSICS::Field,Epetra_SerialDenseVector> elevec_blocks;

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

    // define vector of blocks
    std::vector<std::pair<XFEM::PHYSICS::Field,XFEM::PHYSICS::Field> > block_vec;

    if (eos_gp_pattern == INPAR::FLUID::EOS_GP_Pattern_uvwp)
    {
      // 3D: 4 blocks =  u-u block, v-v block, w-w block and p-p block
      for (std::size_t ifield = 0; ifield < fields.size(); ifield++)
        block_vec.push_back(std::make_pair(fields[ifield],fields[ifield]));
    }
    else if(eos_gp_pattern == INPAR::FLUID::EOS_GP_Pattern_up)
    {
      // 3D: 10 blocks = 3x3 u-u blocks + 1x1 p-p block
      for (std::size_t ifield = 0; ifield < fields.size(); ifield++)
      {
        for (std::size_t jfield = 0; jfield < fields.size(); jfield++)
        {
          bool add = false;
          if (fields[ifield] == fields[jfield])
          {
            add = true;
          }
          else
          {
            if (fields[ifield] != XFEM::PHYSICS::Pres and fields[jfield] != XFEM::PHYSICS::Pres)
              add = true;
          }

          if (add)
           block_vec.push_back(std::make_pair(fields[ifield],fields[jfield]));
        }
      }
    }
    else if(eos_gp_pattern == INPAR::FLUID::EOS_GP_Pattern_full)
    {
      // 3D: 16 blocks = 4x4 uvwp blocks
      dserror("Full pattern does not make sense currently");
      // for the time being no terms that require all blocks are implemented
    }
    else dserror("unknown matrix pattern");


    // resize mat end rhs
    for (std::size_t ib=0; ib<block_vec.size(); ib++)
      elemat_blocks[block_vec[ib]] = Epetra_SerialDenseMatrix(patch_components_lm[block_vec[ib].first].size(),patch_components_lm[block_vec[ib].second].size());

    for (std::size_t ib=0; ib<fields.size(); ib++)
      elevec_blocks[fields[ib]] = Epetra_SerialDenseVector(patch_components_lm[fields[ib]].size());


    //---------------------------------------------------------------------
    // call the element specific evaluate method

    int err = ele->Evaluate(params,*this,
                            lm_master, lm_slave,
                            lm_masterDofPerFieldToPatch, lm_slaveDofPerFieldToPatch,
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
          for(std::size_t ib=0; ib<block_vec.size(); ib++)
            // here, the to fields are equal
            systemmatrix->FEAssemble(-1, elemat_blocks[block_vec[ib]],
                                     patch_components_lm[block_vec[ib].first],
                                     patch_components_lmowner[block_vec[ib].first],
                                     patch_components_lm[block_vec[ib].second]);
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
      for(std::size_t ib=0; ib<fields.size(); ib++)
        LINALG::Assemble(*systemvector,elevec_blocks[fields[ib]],patch_components_lm[fields[ib]],patch_components_lmowner[fields[ib]]);
    }


  return;
}

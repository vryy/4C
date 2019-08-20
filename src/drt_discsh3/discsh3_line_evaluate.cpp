/*----------------------------------------------------------------------*/
/*! \file
\brief line evaluate routines for discsh3 element

\level 3

\maintainer Christoph Meier
*/
/*----------------------------------------------------------------------*/
//#ifdef DISCSH3_H
#include <Teuchos_TimeMonitor.hpp>
//#include <algorithm>
#include "discsh3.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"



/*--------------------------------------------------------------------------*
 |  Integrate a Line Neumann boundary condition (public)     mukherjee 08/15|
 *--------------------------------------------------------------------------*/
int DRT::ELEMENTS::DiscSh3Line::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  // do nothing for now

  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Line Neumann boundary condition (public) mukherjee 08/15|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::DiscSh3Line::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  // do nothing for now

  return 0;
}

/*-------------------------------------------------------------------------------------*
 | Evaluate implementation for internal surface stabilization (public)  mukherjee 08/15|
 *-------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::DiscSh3Line::EvaluateEdges(const Teuchos::ParameterList& params,
    DRT::ELEMENTS::DiscSh3Line* edge,        ///< internal face element
    DRT::Discretization& discretization,     ///< discretization
    std::vector<int>& patchlm,               ///< patch local map
    std::vector<int>& lm_masterToPatch,      ///< local map between master dofs and patchlm
    std::vector<int>& lm_slaveToPatch,       ///< local map between slave dofs and patchlm
    std::vector<int>& lm_faceToPatch,        ///< local map between face dofs and patchlm
    std::vector<int>& lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
    std::vector<int>& lm_slaveNodeToPatch,   ///< local map between slave nodes and nodes in patch
    std::vector<Epetra_SerialDenseMatrix>& elemat_blocks,  ///< element matrix blocks
    std::vector<Epetra_SerialDenseVector>& elevec_blocks   ///< element vector blocks
)
{
  // Get master and slave element
  DRT::Element* pele = edge->ParentMasterElement();  // Master element
  DRT::Element* nele = edge->ParentSlaveElement();   // Slave element

  if (pele == NULL) dserror("pele is NULL");
  if (nele == NULL) dserror("nele is NULL");

  DRT::ELEMENTS::DiscSh3* master_ele = dynamic_cast<DRT::ELEMENTS::DiscSh3*>(pele);
  DRT::ELEMENTS::DiscSh3* slave_ele = dynamic_cast<DRT::ELEMENTS::DiscSh3*>(nele);

  // Create a vector with all the degrees of freedom of master element
  // and neighboring element
  std::vector<FAD> x_FAD_ref(NUMDOF_DISCSH3 + 3, 0.0);  // +3 because of additional node of neighbor
  std::vector<FAD> x_FAD_curr(NUMDOF_DISCSH3 + 3, 0.0);  //+3 because of additional node of neighbor
  std::vector<int> connectivity(0, 0.0);                 // Vector storing

  // add Primary dofs for master element
  AddPrimaryDOFsMaster(*master_ele, *slave_ele, connectivity, x_FAD_ref, discretization, true);
  AddPrimaryDOFsMaster(*master_ele, *slave_ele, connectivity, x_FAD_curr, discretization, false);

  // add Primary dofs for neighboring/slave element
  AddPrimaryDOFsSlave(*master_ele, *slave_ele, x_FAD_ref, discretization, true);
  AddPrimaryDOFsSlave(*master_ele, *slave_ele, x_FAD_curr, discretization, false);

  // Sort Primary dofs
  std::vector<FAD> x_FAD_master_ref(NUMDOF_DISCSH3, 0.0);
  std::vector<FAD> x_FAD_neighbour_ref(NUMDOF_DISCSH3, 0.0);
  std::vector<FAD> x_FAD_master_curr(NUMDOF_DISCSH3, 0.0);
  std::vector<FAD> x_FAD_neighbour_curr(NUMDOF_DISCSH3, 0.0);

  SortPrimaryDOFs(*master_ele, *slave_ele, x_FAD_ref, x_FAD_master_ref, x_FAD_neighbour_ref);
  SortPrimaryDOFs(*master_ele, *slave_ele, x_FAD_curr, x_FAD_master_curr, x_FAD_neighbour_curr);

  // Calculate unit normal at reference config. for master/current element
  LINALG::TMatrix<FAD, 1, 3> NormalMasterRef = CalcSurfaceNormalMaster(x_FAD_master_ref);

  // Calculate unit normal at reference config. for neighbor/slave
  LINALG::TMatrix<FAD, 1, 3> NormalSlaveRef = CalcSurfaceNormalSlave(x_FAD_neighbour_ref);

  // Calculate inclusive angle at reference config.
  FAD ThetaRef = CalcTheta(NormalSlaveRef, NormalMasterRef);

  // Calculate unit normal at spatial config. for master/current element
  LINALG::TMatrix<FAD, 1, 3> NormalMasterCurr = CalcSurfaceNormalMaster(x_FAD_master_curr);

  // Calculate unit normal at spatial config. for neighbour
  LINALG::TMatrix<FAD, 1, 3> NormalSlaveCurr = CalcSurfaceNormalSlave(x_FAD_neighbour_curr);

  // Calculate inclusive angle at spatial config.
  FAD ThetaCurr = CalcTheta(NormalSlaveCurr, NormalMasterCurr);

  // Calculate height of element at reference config. for master
  FAD area_ref_master = master_ele->CalcSurfaceArea(discretization, true);
  FAD height_ref_master = 2 * area_ref_master / this->GetRefEdgeLength();

  // Calculate height of element at reference config. for neighbor
  FAD area_ref_neighbour = slave_ele->CalcSurfaceArea(discretization, true);
  FAD height_ref_neighbour = 2 * area_ref_neighbour / this->GetRefEdgeLength();

  // Average span between two triangles coincident with edge
  FAD Span_neighbours = (height_ref_master + height_ref_neighbour) / 6;

  //  Teuchos::ParameterList StatMechParams =
  //  DRT::Problem::Instance()->StatisticalMechanicsParams();
  //  // Penalty parameter
  //  FAD curvature_penalty=StatMechParams.get<double>("CURV_PENALTY",0.0);
  //  if (curvature_penalty==0)
  //    dserror("Please enter a non-zero CURV_PENALTY");

  FAD curvature_penalty = 0.0;

  // Full energy formulation
  FAD BehaviorFunctCurvature =
      (ThetaCurr - ThetaRef) * pow((this->GetRefEdgeLength() / Span_neighbours), 0.5);

  // Total energy at current step
  FAD TotEnergy = 0.5 * curvature_penalty * ThetaCurr * ThetaCurr *
                  (this->GetRefEdgeLength() / Span_neighbours);

  /*************** Begin Calculate forces analytically ****************************/

  // Calculate gradient of normal vector w.r.t location of coordinate
  LINALG::TMatrix<FAD, 3, 3> DnDx1_master(true);
  LINALG::TMatrix<FAD, 3, 3> DnDx3_master(true);
  LINALG::TMatrix<FAD, 3, 3> DnDx4_master(true);


  CalcGradienNormal(x_FAD_master_curr, DnDx1_master, DnDx3_master, DnDx4_master);

  // Calculate gradient of normal vector w.r.t location of coordinate
  LINALG::TMatrix<FAD, 3, 3> DnDx2_neighbour(true);
  LINALG::TMatrix<FAD, 3, 3> DnDx4_neighbour(true);
  LINALG::TMatrix<FAD, 3, 3> DnDx3_neighbour(true);

  CalcGradienNormal(x_FAD_neighbour_curr, DnDx2_neighbour, DnDx4_neighbour, DnDx3_neighbour);

  // \partialCosTheta/\partial x = Dn1/Dx^T.n2 + Dn2/Dx^T.n1
  LINALG::TMatrix<FAD, 1, 12> DcosThetaDx(true);
  LINALG::TMatrix<FAD, 3, 12> DnDx_master(true);     // [DnDx1 DnDx2 DnDx3]
  LINALG::TMatrix<FAD, 3, 12> DnDx_neighbour(true);  // [DnDx1 DnDx2 DnDx3]



  for (int row = 0; row < 3; row++)
    for (int column = 0; column < 3; column++)
    {
      DnDx_master(row, column) = DnDx1_master(row, column);
      DnDx_master(row, column + 3) = 0;
      DnDx_master(row, column + 6) = DnDx3_master(row, column);
      DnDx_master(row, column + 9) = DnDx4_master(row, column);
    }


  for (int row = 0; row < 3; row++)
    for (int column = 0; column < 3; column++)
    {
      DnDx_neighbour(row, column) = 0;
      DnDx_neighbour(row, column + 3) = DnDx2_neighbour(row, column);
      DnDx_neighbour(row, column + 6) = DnDx3_neighbour(row, column);
      DnDx_neighbour(row, column + 9) = DnDx4_neighbour(row, column);
    }



  // Calculate dcosTheta/Dx
  for (int column = 0; column < 12; column++)
    for (int row = 0; row < 3; row++)
    {
      DcosThetaDx(column) += DnDx_master(row, column) * NormalSlaveCurr(row) +
                             DnDx_neighbour(row, column) * NormalMasterCurr(row);
    }


  LINALG::TMatrix<FAD, 1, 3> crossprod_normal(true);
  // Cross Product
  crossprod_normal = CalcCrossProduct(NormalMasterCurr, NormalSlaveCurr);

  // \partialCrossproduct/\partial x



  LINALG::TMatrix<FAD, 3, 3> SkewSymmTensNormalMaster = SkewSymmetricTrans(NormalMasterCurr);
  LINALG::TMatrix<FAD, 3, 3> SkewSymmTensNormalNeighbour = SkewSymmetricTrans(NormalSlaveCurr);

  LINALG::TMatrix<FAD, 3, 12> DCrossProdDx(true);
  for (int row = 0; row < 3; row++)
    for (int column = 0; column < 12; column++)
      for (int k = 0; k < 3; k++)
        DCrossProdDx(row, column) += SkewSymmTensNormalMaster(row, k) * DnDx_neighbour(k, column) -
                                     SkewSymmTensNormalNeighbour(row, k) * DnDx_master(k, column);

  /* Calculate dsinTheta/dx with sinTheta = abs(n1xn2) */
  LINALG::TMatrix<FAD, 1, 12> DsinThetaDx(true);
  for (int column = 0; column < 12; column++)
    for (int row = 0; row < 3; row++)
      DsinThetaDx(column) += DCrossProdDx(row, column) * crossprod_normal(row) /
                             (pow((crossprod_normal.Dot(crossprod_normal)), 0.5));


  FAD CosTheta = NormalMasterCurr.Dot(NormalSlaveCurr);
  FAD SinTheta = pow((crossprod_normal.Dot(crossprod_normal)), 0.5);

  LINALG::TMatrix<FAD, 1, 12> DThetaDx(true);
  for (int i = 0; i < 12; i++)
  {
    DThetaDx(i) = CosTheta * DsinThetaDx(i) - SinTheta * DcosThetaDx(i);
  }

  LINALG::TMatrix<FAD, 1, 12> ana_force_aux(true);
  LINALG::TMatrix<FAD, 1, 12> BehaviorFunctCurvature_dx(true);

  for (int i = 0; i < NUMDOF_DISCSH3 + 3; i++)
  {
    ana_force_aux(i) = curvature_penalty * (ThetaCurr - ThetaRef) * DThetaDx(i) *
                       this->GetRefEdgeLength() / Span_neighbours;

    // Behavior function
    BehaviorFunctCurvature_dx(i) =
        DThetaDx(i) * pow((this->GetRefEdgeLength() / Span_neighbours), 0.5);
  }

  /*************** End Calculate forces analytically ****************************/


  // Analytical force
  LINALG::TMatrix<FAD, 1, 12> FAD_force_aux(true);
  for (int i = 0; i < NUMDOF_DISCSH3 + 3; i++)
  {
    FAD_force_aux(i) = curvature_penalty * BehaviorFunctCurvature.dx(i) * BehaviorFunctCurvature;
  }

  /*
if(master_ele->IfHaveDamping()&& slave_ele->IfHaveDamping())
{
  dserror("stop");
  //get time step size
  FAD dt = params.get<double>("delta time",0.0);

   //damping coefficients for translational and rotational degrees of freedom
  FAD gamma = params.get<double>("ETA",0.0);

  std::vector<FAD> v_FAD_curr(NUMDOF_DISCSH3+3,0.0); // +3 because of additional node of neighbour
  AddCurrVel(*master_ele,*slave_ele,connectivity,v_FAD_curr, discretization);

  LINALG::TMatrix<FAD,1,9> FAD_force_viscous(true);

  for(int i = 0; i < NUMDOF_DISCSH3; i++)
  {
    for(int j = 0; j < NUMDOF_DISCSH3; j++)
    {
      FAD_force_viscous(i) +=
gamma*BehaviorFunctCurvature.dx(i)*BehaviorFunctCurvature.dx(j)*v_FAD_curr[j];
    }
    ana_force_aux(i)+=FAD_force_viscous(i).val();
  }
}
   */
  //

  LINALG::Matrix<1, 12> error(true);
  LINALG::Matrix<1, 12> ana_force_val(true);
  for (int i = 0; i < 12; i++)
  {
    error(i) = FAD_force_aux(i).val() - ana_force_aux(i).val();
    ana_force_val(i) = ana_force_aux(i).val();
  }

  /* Calculate stiffness matrix */
  ////  std::cout<<"ana_force="<<ana_force<<std::endl;
  LINALG::TMatrix<FAD, 12, 12> stiffmatrix_dummy(true);
  LINALG::Matrix<12, 12> stiffmatrix_dummy_val(true);
  for (int i = 0; i < NUMDOF_DISCSH3 + 3; i++)
    for (int j = 0; j < NUMDOF_DISCSH3 + 3; j++)
    {
      stiffmatrix_dummy(i, j) = ana_force_aux(i).dx(j);
      stiffmatrix_dummy_val(i, j) = stiffmatrix_dummy(i, j).val();
    }


  //****Calculate stiffness contribution to keep the length constant (enforce inextensibility)****//
  // LengthConstrtStiffmass(params, edge, x_FAD_curr, discretization,
  // ana_force_aux,stiffmatrix_dummy);


  /*// If Different damping is taken account
 if(master_ele->IfHaveDamping()&& slave_ele->IfHaveDamping())
 {
   //get time step size
   FAD dt = params.get<double>("delta time",0.0);

    //damping coefficients for translational and rotational degrees of freedom
   FAD gamma = params.get<double>("ETA",0.0);


   std::vector<FAD> v_FAD_curr(NUMDOF_DISCSH3+3,0.0); // +3 because of additional node of neighbour
   AddCurrVel(*master_ele,*slave_ele,connectivity,v_FAD_curr, discretization);

   LINALG::TMatrix  <FAD,9,9> stiffmat_viscous(true);

   for (int i=0; i<NUMDOF_DISCSH3; i++)
     for (int j=0; j<NUMDOF_DISCSH3; j++)
     {
       for (int k=0; k<NUMDOF_DISCSH3; k++)
       {
         stiffmat_viscous(i,j)+=  gamma*
 BehaviorFunctCurvature_dx(i).dx(j)*BehaviorFunctCurvature_dx(k)*v_FAD_curr[k]
                                + gamma*
 BehaviorFunctCurvature_dx(i)*BehaviorFunctCurvature_dx(k).dx(j)*v_FAD_curr[k]
                                + gamma*
 BehaviorFunctCurvature_dx(i)*BehaviorFunctCurvature_dx(k)*(j==k)/dt;
       }
       stiffmatrix_dummy(i,j)+=stiffmat_viscous(i,j).val();
     }
 }
   */

  /*%%%%%%%%%%%%%%    DEBUG   %%%%%%%%%%%%%%%*/
  /*LINALG::Matrix<12,12> stiffmatrix_print(true);
 LINALG::Matrix<12,1> force_print(true);
   Epetra_SerialDenseMatrix stiffmatrix_Epetra(12,12,true);
//   stiffmatrix_Epetra.
 for (int i=0; i<NUMDOF_DISCSH3+3; i++)
   for (int j=0; j<NUMDOF_DISCSH3+3; j++)
   {
     stiffmatrix_print(i,j)=stiffmatrix_dummy(i,j).val();
     stiffmatrix_Epetra(i,j)=stiffmatrix_dummy(i,j).val();
     force_print(i)=ana_force_aux(i).val();
   }

 if(edge->Id()==EdgeID)
 {
   std::cout<<"Stiffmatrix_print="<<stiffmatrix_print<<std::endl;
   LINALG::PrintSerialDenseMatrixInMatlabFormat("matrix_element.m",stiffmatrix_Epetra,true);
 }
   */

  LINALG::Matrix<NUMDOF_DISCSH3, 1> elevector_m(true);  // element vector master block
  LINALG::Matrix<NUMDOF_DISCSH3, 1> elevector_s(true);  // element vector slave block


  RemapForceVectoMasterSlave(
      ana_force_aux, elevector_m, elevector_s, connectivity, lm_slaveNodeToPatch);


  // element matrices in block structure master vs. slave
  LINALG::Matrix<NUMDOF_DISCSH3, NUMDOF_DISCSH3> elematrix_mm(
      true);  // element matrix master-master block
  LINALG::Matrix<NUMDOF_DISCSH3, NUMDOF_DISCSH3> elematrix_ms(
      true);  // element matrix master-slave block
  LINALG::Matrix<NUMDOF_DISCSH3, NUMDOF_DISCSH3> elematrix_sm(
      true);  // element matrix slave-master block
  LINALG::Matrix<NUMDOF_DISCSH3, NUMDOF_DISCSH3> elematrix_ss(
      true);  // element matrix slave-slave block

  RemapStiffmattoMasterSlave(stiffmatrix_dummy, elematrix_mm, elematrix_ms, elematrix_sm,
      elematrix_ss, lm_masterNodeToPatch, lm_slaveNodeToPatch, connectivity);


  int numblocks = elemat_blocks.size();

  if (numblocks == 9)
  {
    // 3D: reassemble all uvwp blocks
    for (int idim = 0; idim < NODDOF_DISCSH3; idim++)
    {
      ReassembleRHSBlock(idim, elevec_blocks[idim], elevector_m, elevector_s, lm_masterNodeToPatch,
          lm_slaveNodeToPatch);

      for (int jdim = 0; jdim < NODDOF_DISCSH3; jdim++)
      {
        ReassembleMATBlock(idim, jdim, elemat_blocks[idim * NODDOF_DISCSH3 + jdim], elematrix_mm,
            elematrix_ms, elematrix_sm, elematrix_ss, lm_masterNodeToPatch, lm_slaveNodeToPatch);
      }
    }
  }
  else
    dserror("unknown assembly pattern for given number of epetra block matrices");

  return 0;
}

/*---------------------------------------------------------------------------------------*
 | Assemble internal faces integrals using data from both parent elements mukherjee 08/15|
 *---------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3Line::AssembleInternalFacesUsingNeighborData(
    const Teuchos::ParameterList& params,
    DRT::ELEMENTS::DiscSh3Line* intface,              ///< internal face element
    Teuchos::RCP<MAT::Material>& material,            ///< material for face stabilization
    std::vector<int>& nds_master,                     ///< nodal dofset w.r.t. master element
    std::vector<int>& nds_slave,                      ///< nodal dofset w.r.t. slave element
    DRT::DiscretizationFaces& discretization,         ///< faces discretization
    Teuchos::RCP<LINALG::SparseMatrix> systemmatrix,  ///< systemmatrix
    Teuchos::RCP<Epetra_Vector> systemvector          ///< systemvector
)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::DiscSh3Line::AssembleInternalFacesUsingNeighborData");

  // decide which terms have to be assembled for the current face and decide the assembly pattern,
  // return if no assembly required
  //  bool stab_required = fldpara_intface_->SetFaceSpecificFluidXFEMParameter(face_type, params);

  // do not assemble if no stabilization terms activated for this face
  //  if(!stab_required) return;

  if (!discretization.Filled()) dserror("FillComplete() was not called");
  if (!discretization.HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  const bool assemblemat = systemmatrix != Teuchos::null;
  const bool assemblevec = systemvector != Teuchos::null;

  //----------------------- create patchlm -----------------

  const int numnode_master = intface->ParentMasterElement()->NumNode();
  const int numnode_slave = intface->ParentSlaveElement()->NumNode();
  const int numnode_face = intface->NumNode();

  const int numnodeinpatch = numnode_master + numnode_slave - numnode_face;

  // local maps for patch dofs
  std::vector<int> lm_patch;
  lm_patch.reserve(numnodeinpatch * NODDOF_DISCSH3);

  // local maps between master/slave dofs and position in patch dofs (lm_patch)
  std::vector<int> lm_masterToPatch;
  lm_masterToPatch.reserve(numnode_master * NODDOF_DISCSH3);
  std::vector<int> lm_slaveToPatch;
  lm_slaveToPatch.reserve(numnode_slave * NODDOF_DISCSH3);
  std::vector<int> lm_faceToPatch;
  lm_faceToPatch.reserve(numnode_face * NODDOF_DISCSH3);

  // local maps between master/slave nodes and position in patch nodes
  std::vector<int> lm_masterNodeToPatch;
  lm_masterNodeToPatch.reserve(numnode_master);
  std::vector<int> lm_slaveNodeToPatch;
  lm_slaveNodeToPatch.reserve(numnode_slave);

  // create patch location vector combining master element, slave element and face element
  intface->PatchLocationVector(discretization, nds_master, nds_slave, lm_patch,
      // lm_master, lm_slave, lm_face,
      lm_masterToPatch, lm_slaveToPatch, lm_faceToPatch, lm_masterNodeToPatch, lm_slaveNodeToPatch);

  // patch_lm for Velx,Vely,Velz, Pres field
  std::vector<std::vector<int>> patch_components_lm(NODDOF_DISCSH3);
  std::vector<std::vector<int>> patch_components_lmowner(NODDOF_DISCSH3);

  for (int i = 0; i < NODDOF_DISCSH3; i++)
  {
    patch_components_lm[i].reserve(numnodeinpatch);
    patch_components_lmowner[i].reserve(numnodeinpatch);
  }

  // modify the patch owner to the owner of the internal face element
  const int owner = intface->Owner();
  std::vector<int> patchlm_owner(lm_patch.size(), owner);

  for (unsigned i = 0; i < lm_patch.size(); i++)
  {
    const int field = i % NODDOF_DISCSH3;

    // i%4 yields the Velx,Vely,Velz,Pres field
    patch_components_lm[field].push_back(lm_patch[i]);
    patch_components_lmowner[field].push_back(owner);
  }

#ifdef DEBUG
  for (int isd = 0; isd < NODDOF_DISCSH3; isd++)
    if ((int)(patch_components_lm[isd].size()) != numnodeinpatch)
      dserror("patch_components_lm[%d] has wrong size", isd);

#endif

  int numblocks = NODDOF_DISCSH3 * NODDOF_DISCSH3;

  // define element matrices and vectors
  std::vector<Epetra_SerialDenseMatrix> elemat_blocks(numblocks);
  std::vector<Epetra_SerialDenseVector> elevec_blocks(
      NODDOF_DISCSH3);  // 3D: 4 vectors for u,v,w,p components, 2D: 3 vectors for u,v,p


  for (int b = 0; b < numblocks; b++)
  {
    int err = elemat_blocks[b].Shape(
        numnodeinpatch, numnodeinpatch);  // new shape and init values to zero

    if (err != 0) dserror("element matrix Shape not successful");
  }

  for (int b = 0; b < NODDOF_DISCSH3; b++)
  {
    int err = elevec_blocks[b].Size(numnodeinpatch);  // new size and init values to zero

    if (err != 0) dserror("element matrix Shape not successful");
  }

  //---------------------------------------------------------------------
  // call the element specific evaluate method

  int err =
      EvaluateEdges(params, intface, discretization, lm_patch, lm_masterToPatch, lm_slaveToPatch,
          lm_faceToPatch, lm_masterNodeToPatch, lm_slaveNodeToPatch, elemat_blocks, elevec_blocks);


  if (err) dserror("error while evaluating elements");

  //---------------------------------------------------------------------------------------
  // assemble systemmatrix
  // calls the Assemble function for EpetraFECrs matrices including communication of non-row entries
  if (assemblemat)
  {
    for (int i = 0; i < NODDOF_DISCSH3; i++)
    {
      for (int j = 0; j < NODDOF_DISCSH3; j++)
      {
        systemmatrix->FEAssemble(elemat_blocks[i * NODDOF_DISCSH3 + j], patch_components_lm[i],
            patch_components_lmowner[i], patch_components_lm[j]);
      }
    }
  }

  //---------------------------------------------------------------------------------------
  // assemble systemvector
  if (assemblevec)
  {
    // REMARK:: call Assemble without lmowner
    // to assemble the residual_col vector on only row elements also column nodes have to be
    // assembled do not exclude non-row nodes (modify the real owner to myowner) after assembly the
    // col vector it has to be exported to the row residual_ vector using the 'Add' flag to get the
    // right value for shared nodes
    for (int i = 0; i < NODDOF_DISCSH3; i++)
    {
      LINALG::Assemble(
          *systemvector, elevec_blocks[i], patch_components_lm[i], patch_components_lmowner[i]);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  create the patch location vector (public)           mukherjee 06/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3Line::PatchLocationVector(
    DRT::Discretization& discretization,     ///< discretization
    std::vector<int>& nds_master,            ///< nodal dofset w.r.t master parent element
    std::vector<int>& nds_slave,             ///< nodal dofset w.r.t slave parent element
    std::vector<int>& patchlm,               ///< local map for gdof ids for patch of elements
    std::vector<int>& lm_masterToPatch,      ///< local map between lm_master and lm_patch
    std::vector<int>& lm_slaveToPatch,       ///< local map between lm_slave and lm_patch
    std::vector<int>& lm_faceToPatch,        ///< local map between lm_face and lm_patch
    std::vector<int>& lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
    std::vector<int>& lm_slaveNodeToPatch    ///< local map between slave nodes and nodes in patch
)
{
  //-----------------------------------------------------------------------
  const int m_numnode = ParentMasterElement()->NumNode();
  DRT::Node** m_nodes = ParentMasterElement()->Nodes();


  if (m_numnode != static_cast<int>(nds_master.size()))
  {
    throw std::runtime_error("wrong number of nodes for master element");
  }


  //-----------------------------------------------------------------------
  const int s_numnode = ParentSlaveElement()->NumNode();
  DRT::Node** s_nodes = ParentSlaveElement()->Nodes();

  if (s_numnode != static_cast<int>(nds_slave.size()))
  {
    throw std::runtime_error("wrong number of nodes for slave element");
  }


  //-----------------------------------------------------------------------
  const int f_numnode = NumNode();
  DRT::Node** f_nodes = Nodes();

  // for each master node, the offset for node's dofs in master_lm
  std::map<int, int> m_node_lm_offset;

  //----------------------------------------------------
  int curr_patch_lm_size = 0;  // patch_lm.size() (equal to master_lm.size() during the fill of
                               // patch data with master data)

  // ---------------------------------------------------
  const int dofset = 0;  // assume dofset 0

  int patchnode_count = 0;


  // fill patch lm with master's nodes
  for (int k = 0; k < m_numnode; ++k)
  {
    DRT::Node* node = m_nodes[k];
    std::vector<int> dof;
    discretization.Dof(dof, node, dofset, nds_master[k]);

    const int size = dof.size();

    // insert a pair of node-Id and current length of master_lm ( to get the start offset for node's
    // dofs)
    m_node_lm_offset.insert(std::pair<int, int>(node->Id(), curr_patch_lm_size));

    for (int j = 0; j < size; ++j)
    {
      lm_masterToPatch.push_back(curr_patch_lm_size);

      int actdof = dof[j];
      patchlm.push_back(actdof);
      curr_patch_lm_size++;
    }

    lm_masterNodeToPatch.push_back(patchnode_count);

    patchnode_count++;
  }

  // ---------------------------------------------------
  // fill patch lm with missing slave's nodes and extract slave's lm from patch_lm

  for (int k = 0; k < s_numnode; ++k)
  {
    DRT::Node* node = s_nodes[k];

    // slave node already contained?
    std::map<int, int>::iterator m_offset;
    m_offset = m_node_lm_offset.find(node->Id());


    if (m_offset == m_node_lm_offset.end())  // node not included yet
    {
      std::vector<int> dof;  // = discretization.Dof(dofset,node);
      discretization.Dof(dof, node, dofset, nds_slave[k]);

      const int size = dof.size();

      for (int j = 0; j < size; ++j)
      {
        lm_slaveToPatch.push_back(curr_patch_lm_size);


        int actdof = dof[j];
        patchlm.push_back(actdof);
        curr_patch_lm_size++;
      }

      lm_slaveNodeToPatch.push_back(patchnode_count);


      patchnode_count++;
    }
    else  // node is also a master's node
    {
      // get maximum of numdof per node with the help of master and/or slave element (returns 4 in
      // 3D case, does not return dofset's numdof)
      const int size = NumDofPerNode(*node);

      int offset = m_offset->second;

      for (int j = 0; j < size; ++j)
      {
        // copy from lm_masterToPatch
        lm_slaveToPatch.push_back(lm_masterToPatch[offset + j]);
      }
      if (offset % size != 0)
        dserror("there was at least one node with not %d dofs per node", size);
      int patchnode_index = offset / size;

      lm_slaveNodeToPatch.push_back(patchnode_index);
      // no patchnode_count++; (node already contained)
    }
  }

  // ---------------------------------------------------
  // extract face's lm from patch_lm

  for (int k = 0; k < f_numnode; ++k)
  {
    DRT::Node* node = f_nodes[k];

    // face node must be contained
    std::map<int, int>::iterator m_offset;
    m_offset = m_node_lm_offset.find(node->Id());

    if (m_offset != m_node_lm_offset.end())  // node not included yet
    {
      // get maximum of numdof per node with the help of master and/or slave element (returns 4 in
      // 3D case, does not return dofset's numdof)
      const int size = NumDofPerNode(*node);

      int offset = m_offset->second;

      for (int j = 0; j < size; ++j)
      {
        // copy from lm_masterToPatch
        lm_faceToPatch.push_back(lm_masterToPatch[offset + j]);
      }
    }
    else
      throw std::runtime_error("face's nodes not contained in masternodes_offset map");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  create the patch location vector (public)           mukherjee 06/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3Line::RemapForceVectoMasterSlave(
    LINALG::TMatrix<FAD, 1, 12>& force,  ///< nodal dofset w.r.t master parent element
    LINALG::Matrix<NUMDOF_DISCSH3, 1>& elevector_m, LINALG::Matrix<NUMDOF_DISCSH3, 1>& elevector_s,
    std::vector<int>& connectivity, std::vector<int>& lm_slaveNodeToPatch)
{
  // map force back to the master according to dis
  for (int i = 0; i < NODDOF_DISCSH3; i++)
  {
    // first node of master
    elevector_m(i) = force(3 * connectivity[0] + i).val();
    // second node of master element
    elevector_m(i + 3) = force(3 * connectivity[1] + i).val();
    // third node of master element
    elevector_m(i + 6) = force(3 * connectivity[2] + i).val();
  }

  // slave element
  for (int i = 0; i < NUMNOD_DISCSH3; i++)
  {
    int rpatch = lm_slaveNodeToPatch[i];
    if (rpatch == 3)
    {
      for (int j = 0; j < NODDOF_DISCSH3; j++)
      {
        elevector_s(3 * i + j) = force(j + 3).val();
      }
    }
  }
  return;
}

/*-------------------------------------------------------------------------*
 |  create the patch location vector (public)              mukherjee 06/15 |
 *-------------------------------------------------------------------------*/
void DRT::ELEMENTS::DiscSh3Line::RemapStiffmattoMasterSlave(
    LINALG::TMatrix<FAD, 12, 12>& stiffmatrix,  ///< nodal dofset w.r.t master parent element
    LINALG::Matrix<NUMDOF_DISCSH3, NUMDOF_DISCSH3>& elematrix_mm,
    LINALG::Matrix<NUMDOF_DISCSH3, NUMDOF_DISCSH3>& elematrix_ms,
    LINALG::Matrix<NUMDOF_DISCSH3, NUMDOF_DISCSH3>& elematrix_sm,
    LINALG::Matrix<NUMDOF_DISCSH3, NUMDOF_DISCSH3>& elematrix_ss,
    std::vector<int>& lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
    std::vector<int>& lm_slaveNodeToPatch,   ///< local map between slave nodes and nodes in patch
    std::vector<int>& connectivity)
{
  // master-master block
  for (int i = 0; i < NODDOF_DISCSH3; i++)
    for (int j = 0; j < NODDOF_DISCSH3; j++)
    {
      // rows corresponding to first node of master element
      elematrix_mm(i, j) = stiffmatrix(3 * connectivity[0] + i, 3 * connectivity[0] + j).val();
      elematrix_mm(i, j + 3) = stiffmatrix(3 * connectivity[0] + i, 3 * connectivity[1] + j).val();
      elematrix_mm(i, j + 6) = stiffmatrix(3 * connectivity[0] + i, 3 * connectivity[2] + j).val();
      // rows corresponding to second node of master element
      elematrix_mm(i + 3, j) = stiffmatrix(3 * connectivity[1] + i, 3 * connectivity[0] + j).val();
      elematrix_mm(i + 3, j + 3) =
          stiffmatrix(3 * connectivity[1] + i, 3 * connectivity[1] + j).val();
      elematrix_mm(i + 3, j + 6) =
          stiffmatrix(3 * connectivity[1] + i, 3 * connectivity[2] + j).val();
      // rows corresponding to third node of master element
      elematrix_mm(i + 6, j) = stiffmatrix(3 * connectivity[2] + i, 3 * connectivity[0] + j).val();
      elematrix_mm(i + 6, j + 3) =
          stiffmatrix(3 * connectivity[2] + i, 3 * connectivity[1] + j).val();
      elematrix_mm(i + 6, j + 6) =
          stiffmatrix(3 * connectivity[2] + i, 3 * connectivity[2] + j).val();
    }

  // master-slave block
  for (int row_node = 0; row_node < NUMNOD_DISCSH3; row_node++)
    for (int col_node = 0; col_node < NUMNOD_DISCSH3; col_node++)
    {
      int rpatch = lm_masterNodeToPatch[row_node];
      int cpatch = lm_slaveNodeToPatch[col_node];
      if (rpatch != 3 && cpatch == 3)
      {
        int rindex = connectivity[rpatch];
        for (int i = 0; i < NODDOF_DISCSH3; i++)
          for (int j = 0; j < NODDOF_DISCSH3; j++)
            elematrix_ms(3 * row_node + i, 3 * col_node + j) =
                stiffmatrix(i + 3 * rindex, j + cpatch).val();
      }
    }

  // slave-slave block
  for (int row_node = 0; row_node < NUMNOD_DISCSH3; row_node++)
    for (int col_node = 0; col_node < NUMNOD_DISCSH3; col_node++)
    {
      int rpatch = lm_slaveNodeToPatch[row_node];
      int cpatch = lm_slaveNodeToPatch[col_node];
      if (rpatch == 3 && cpatch == 3)  // 3 because of the slave DOFs
      {
        for (int i = 0; i < NODDOF_DISCSH3; i++)
          for (int j = 0; j < NODDOF_DISCSH3; j++)
            elematrix_ss(3 * row_node + i, 3 * col_node + j) =
                stiffmatrix(i + rpatch, j + cpatch).val();
      }
    }


  // slave-master block
  for (int row_node = 0; row_node < NUMNOD_DISCSH3; row_node++)
    for (int col_node = 0; col_node < NUMNOD_DISCSH3; col_node++)
    {
      int rpatch = lm_slaveNodeToPatch[row_node];
      int cpatch = lm_masterNodeToPatch[col_node];
      if (rpatch == 3 && cpatch != 3)
      {
        int cindex = connectivity[cpatch];

        for (int i = 0; i < NODDOF_DISCSH3; i++)
          for (int j = 0; j < NODDOF_DISCSH3; j++)
            elematrix_sm(3 * row_node + i, 3 * col_node + j) =
                stiffmatrix(i + rpatch, j + 3 * cindex).val();
      }
    }

  return;
}

/*-------------------------------------------------------------------------*
 |  Calculate stiffness contribution to keep the                           |
 |  length constant (enforce inextensibility)              mukherjee 06/15 |
 *-------------------------------------------------------------------------*/

void DRT::ELEMENTS::DiscSh3Line::LengthConstrtStiffmass(Teuchos::ParameterList& params,
    DRT::ELEMENTS::DiscSh3Line* edge,  ///< internal face element
    std::vector<FAD>& x_FAD_curr,
    DRT::Discretization& discretization,  ///< discretization
    LINALG::TMatrix<FAD, 1, 12>& ana_force_aux, LINALG::TMatrix<FAD, 12, 12>& stiff_dummy)
{
  //  Teuchos::ParameterList StatMechParams =
  //  DRT::Problem::Instance()->StatisticalMechanicsParams(); FAD
  //  line_penalty=StatMechParams.get<double>("LINE_PENALTY",0.0); if (line_penalty==0)
  //    dserror("Please enter a non-zero LINE_PENALTY");

  FAD line_penalty = 0.0;

  FAD Edge_length_ref = GetEdgeLengthPrevTimeStep();

  std::vector<double> x_line_curr;
  std::vector<FAD> x_FAD_line_curr(6, 0.0);
  for (int i = 0; i < 6; i++)
  {
    x_line_curr.push_back(x_FAD_curr[i + 6].val());
    x_FAD_line_curr[i] = x_line_curr[i];
    x_FAD_line_curr[i].diff(i, 6);
  }


  FAD Edge_length_curr = GetCurrEdgeLength(x_FAD_line_curr);

  FAD BehaviorFunctLine = std::pow(Edge_length_ref, 0.5) * (1 - Edge_length_curr / Edge_length_ref);

  LINALG::TMatrix<FAD, 1, 6> FAD_force(true);
  // Calculate FAD force
  for (int i = 0; i < 6; i++)
  {
    FAD_force(i) = line_penalty * BehaviorFunctLine.dx(i) * BehaviorFunctLine;
  }

  // Get the expressions analytically
  //  std::vector<FAD> x_FAD(12,0.0);
  // Get Spatial positions 2 nodes, 3 dims
  LINALG::TMatrix<FAD, 1, 6> Ana_force(true);
  LINALG::TMatrix<FAD, 1, 6> grad_edge(true);
  LINALG::Matrix<1, 6> tol(true);

  for (int i = 0; i < 3; i++)
  {
    Ana_force(i) = -line_penalty * BehaviorFunctLine *
                   (x_FAD_line_curr[i] - x_FAD_line_curr[i + 3]) /
                   (Edge_length_curr * std::pow(Edge_length_ref, 0.5));
    Ana_force(i + 3) = line_penalty * BehaviorFunctLine *
                       (x_FAD_line_curr[i] - x_FAD_line_curr[i + 3]) /
                       (Edge_length_curr * std::pow(Edge_length_ref, 0.5));

    //      grad_edge(i)   = (x_FAD_line_curr[i]-x_FAD_line_curr[i+3])/Edge_length_currAlt;
    //      grad_edge(i+3) = -(x_FAD_line_curr[i]-x_FAD_line_curr[i+3])/Edge_length_currAlt;

    //      tol(i)  =Ana_force(i+6).val()-FAD_force(i).val();
    //      tol(i+3)=Ana_force(i+9).val()-FAD_force(i+3).val();
  }

  LINALG::TMatrix<FAD, 6, 6> FAD_stiff(true);
  LINALG::TMatrix<FAD, 6, 6> FAD_stiff_alt(true);
  LINALG::Matrix<6, 6> FAD_stiff_val(true);
  LINALG::TMatrix<FAD, 3, 3> FAD_aux(true);
  LINALG::Matrix<6, 6> FAD_stiff_alt_val(true);

  LINALG::TMatrix<FAD, 6, 6> grad2_edge(true);


  for (int i = 0; i < 6; i++)
  {
    for (int j = 0; j < 6; j++)
    {
      grad2_edge(i, j) = grad_edge(i).dx(j);
      FAD_stiff_alt(i, j) =
          line_penalty * grad_edge(i) * grad_edge(j) / Edge_length_ref -
          line_penalty * (1 - Edge_length_curr / Edge_length_ref) * grad2_edge(i, j);
      FAD_stiff(i, j) = Ana_force(i).dx(j);
      FAD_stiff_val(i, j) = FAD_stiff(i, j).val();
      FAD_stiff_alt_val(i, j) = FAD_stiff_alt(i, j).val();
    }
  }


  for (int i = 0; i < 6; i++)
  {
    ana_force_aux(i + 6) += Ana_force(i).val();
    for (int j = 0; j < 6; j++)
    {
      stiff_dummy(i + 6, j + 6) += FAD_stiff(i, j).val();
    }
  }

  return;
}

//#endif

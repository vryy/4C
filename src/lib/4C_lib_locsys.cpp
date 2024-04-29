/*----------------------------------------------------------------------------*/
/*! \file

\brief Class controlling local coordinate systems on points, lines and
surfaces and supplying all necessary transformation methods for parallel
vectors and matrices.

\level 1


*/
/*----------------------------------------------------------------------------*/

#include "4C_lib_locsys.hpp"

#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_discretization_fem_general_largerotations.hpp"
#include "4C_global_data.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN



/*-------------------------------------------------------------------*
 |  ctor (public)                                         meier 06/13|
 *-------------------------------------------------------------------*/
DRT::UTILS::LocsysManager::LocsysManager(DRT::Discretization& discret)
    : discret_(discret), dim_(-1), numlocsys_(-1), locsysfunct_(false)
{
  // get problem dimension (2D or 3D) and store into dim_
  dim_ = GLOBAL::Problem::Instance()->NDim();

  if (Dim() != 2 && Dim() != 3) FOUR_C_THROW("Locsys problem must be 2D or 3D");

  // get node row layout of discretization
  const Epetra_Map* noderowmap = discret_.NodeRowMap();

  // create locsys vector and initialize to -1
  locsystoggle_ = CORE::LINALG::CreateVector(*noderowmap, false);
  locsystoggle_->PutScalar(-1.0);

  // check for locsys boundary conditions
  Discret().GetCondition("Locsys", locsysconds_);
  numlocsys_ = (int)locsysconds_.size();
  id_.resize(numlocsys_);
  typelocsys_.resize(numlocsys_);

  for (int i = 0; i < NumLocsys(); ++i)
  {
    // Assign internal locsys id (is this still being used?)
    id_[i] = locsysconds_[i]->Id();

    // Check for already existing ConditionID and add it, if not already existing
    const int* locsysId = locsysconds_[i]->GetIf<int>("ConditionID");
    if (locsysId)
    {
      if ((*locsysId) != i) FOUR_C_THROW("Locsys condition has non-matching ID");
    }
    else
    {
      locsysconds_[i]->Add("ConditionID", i);
    }
  }

  // Set boolean that indicates, if a locsys warning has already been thrown, to false
  warning_thrown_ = false;
}

/*-------------------------------------------------------------------*
 *-------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::Update(
    const double time, std::vector<Teuchos::RCP<Epetra_Vector>> nodenormals)
{
  nodenormals_ = std::move(nodenormals);
  // IMPORTANT NOTE:
  // The definition of local coordinate systems only makes sense in
  // combination with Dirichlet boundary conditions. This means that
  // in order to define a boundary condition, both locsys
  // AND Dirichlet condition have to formulated for the same entity
  // (i.e. point, line, surface, volume).

  // LIMITATIONS:
  // - So far locsys only works for 2D and 3D solids
  // - Due to this limitation it's necessary to distinguish between this different element types
  //   by means of there nodal DoFs. If further element types are integrated into locsys
  //   more elaborate criteria might be useful.

  // If we have no functions in the locsys conditions the whole Setup method is only
  // conducted once in the constructor (where time is set to -1.0).
  if (time >= 0.0 and !locsysfunct_) return;

  // get dof row map of discretization
  const Epetra_Map* dofrowmap = discret_.DofRowMap();

  // get node row layout of discretization
  const Epetra_Map* noderowmap = discret_.NodeRowMap();

  // Since also time dependent conditions are possible we clear all local systems in the beginning
  nodalrotvectors_.clear();

  // As for Dirichlet conditions, we keep to a very strict hierarchy
  // for evaluation of the Locsys conditions: Volume locsys conditions
  // are evaluated first, followed by Surface and Line locsys conditions
  // and finally Point locsys conditions. This means that nodes carrying
  // different types of locsys conditions are dominated by the rule "Point
  // above Line above Surface above Volume". When two locsys conditions of
  // the same type are defined for one node, ordering in the input file matters!
  std::vector<CORE::Conditions::ConditionType> geo_hierarchy;
  geo_hierarchy.emplace_back(CORE::Conditions::VolumeLocsys);
  geo_hierarchy.emplace_back(CORE::Conditions::SurfaceLocsys);
  geo_hierarchy.emplace_back(CORE::Conditions::LineLocsys);
  geo_hierarchy.emplace_back(CORE::Conditions::PointLocsys);

  //**********************************************************************
  // read locsys conditions in given hierarchical order
  //**************************+*******************************************
  for (auto geo_level : geo_hierarchy)
  {
    for (int i = 0; i < NumLocsys(); ++i)
    {
      DRT::Condition* currlocsys = locsysconds_[i];

      // safety check
      if (currlocsys->Type() != CORE::Conditions::VolumeLocsys and
          currlocsys->Type() != CORE::Conditions::SurfaceLocsys and
          currlocsys->Type() != CORE::Conditions::LineLocsys and
          currlocsys->Type() != CORE::Conditions::PointLocsys)
        FOUR_C_THROW("Unknown type of locsys condition!");

      if (currlocsys->Type() == geo_level)
      {
        typelocsys_[i] = currlocsys->Type();

        const auto* rotangle = &currlocsys->Get<std::vector<double>>("rotangle");
        const auto* funct = &currlocsys->Get<std::vector<int>>("funct");
        const auto* useUpdatedNodePos = &currlocsys->Get<int>("useupdatednodepos");

        const auto* useConsistentNodeNormal =
            (currlocsys->Type() == CORE::Conditions::SurfaceLocsys or
                currlocsys->Type() == CORE::Conditions::LineLocsys)
                ? currlocsys->GetIf<int>("useconsistentnodenormal")
                : nullptr;

        const std::vector<int>* nodes = currlocsys->GetNodes();

        if (currlocsys->Type() == CORE::Conditions::SurfaceLocsys or
            currlocsys->Type() == CORE::Conditions::LineLocsys)
        {
          // Check, if we have time dependent locsys conditions (through functions)
          if (((*funct)[0] > 0 or (*funct)[1] > 0 or (*funct)[2] > 0) or
              (((*useConsistentNodeNormal) == 1) and ((*useUpdatedNodePos) == 1)))
            locsysfunct_ = true;
        }
        else if (currlocsys->Type() == CORE::Conditions::VolumeLocsys or
                 currlocsys->Type() == CORE::Conditions::PointLocsys)
        {
          // Check, if we have time dependent locsys conditions (through functions)
          if (((*funct)[0] > 0 or (*funct)[1] > 0 or (*funct)[2] > 0)) locsysfunct_ = true;
        }

        // Here we have the convention that 2D problems "live" in the global xy-plane.
        if (Dim() == 2 and ((*rotangle)[0] != 0 or (*rotangle)[1] != 0))
        {
          FOUR_C_THROW(
              "For 2D problems (xy-plane) the vector ROTANGLE has to be parallel to the global "
              "z-axis!");
        }

        if ((currlocsys->Type() == CORE::Conditions::SurfaceLocsys or
                currlocsys->Type() == CORE::Conditions::LineLocsys) and
            (*useConsistentNodeNormal) == 1)
          CalcRotationVectorForNormalSystem(i, time);
        else
        {
          // Check, if the updated node positions shall be used for evaluation of the functions
          // 'funct'
          Teuchos::RCP<const Epetra_Vector> dispnp;
          if (((*useUpdatedNodePos) == 1) && (time >= 0.0))
          {
            dispnp = Discret().GetState("dispnp");
            if (dispnp == Teuchos::null)
            {
              FOUR_C_THROW(
                  "Locsys: Cannot find state 'dispnp'! You need to set the state 'dispnp' before "
                  "calling the locsys setup.");
            }
          }

          // Each component j of the pseudo rotation vector that rotates the global xyz system onto
          // the local system assigned to each node consists of a constant, a time dependent and
          // spatially variable part: currotangle_j(x,t) = rotangle_j * funct_j(t,x)
          CORE::LINALG::Matrix<3, 1> currotangle(true);

          for (int nodeGID : *nodes)
          {
            bool havenode = Discret().HaveGlobalNode(nodeGID);
            if (!havenode) continue;

            // Weights of rotations vector due to temporal and spatial function
            for (int j = 0; j < 3; j++)
            {
              // factor given by spatial function
              double functfac = 1.0;
              if ((*funct)[j] > 0)
              {
                DRT::Node* node = Discret().gNode(nodeGID);

                // Determine node position, which shall be used for evaluating the function, and
                // evaluate it
                if (((*useUpdatedNodePos) == 1) && (time >= 0.0))
                {
                  // Obtain current displacement for node
                  std::vector<int> lm;
                  Discret().Dof(node, lm);

                  std::vector<double> currDisp;
                  currDisp.resize(lm.size());

                  CORE::FE::ExtractMyValues(*dispnp, currDisp, lm);

                  // Calculate current position for node
                  std::vector<double> currPos(Dim());
                  const auto& xp = node->X();

                  for (int dim = 0; dim < Dim(); ++dim) currPos[dim] = xp[dim] + currDisp[dim];

                  // Evaluate function with current node position
                  functfac =
                      (GLOBAL::Problem::Instance()->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(
                           (*funct)[j] - 1))
                          .Evaluate(currPos.data(), time, j);
                }
                else
                {
                  // Evaluate function with reference node position
                  functfac =
                      (GLOBAL::Problem::Instance()->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(
                           (*funct)[j] - 1))
                          .Evaluate(node->X().data(), time, j);
                }
              }
              currotangle(j) = (*rotangle)[j] * functfac;
            }

            nodalrotvectors_[nodeGID] = currotangle;

            int indices = nodeGID;
            double values = i;
            locsystoggle_->ReplaceGlobalValues(1, &values, &indices);
          }
        }
      }
    }
  }

  if (time < 0.0) Print();

  // When building the transformation matrix we apply a node-by-node
  // strategy. The global matrix trafo_ will consist of nodal blocks
  // of dimension (numdof)x(numdof). The following code block is designed
  // for 2D and 3D solid elements. If special fields are constructed with
  // more than dim geometric dofs, i.e. that have to be transformed, then
  // the following code might have to be modified!

  //**********************************************************************
  // Build transformation matrix trafo_
  //**********************************************************************

  // we need to make sure that two nodes sharing the same dofs are not
  // transformed twice. This is a NURBS/periodic boundary feature.
  Teuchos::RCP<Epetra_Vector> already_processed = CORE::LINALG::CreateVector(*dofrowmap, true);
  already_processed->PutScalar(0.0);

  // Perform a check for zero diagonal elements. They will crash the SGS-like preconditioners
  bool sanity_check = false;

  // GIDs of all DoFs subjected to local co-ordinate systems
  std::set<int> locsysdofset;

  trafo_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(*dofrowmap, 3));

  for (int i = 0; i < noderowmap->NumMyElements(); ++i)
  {
    int nodeGID = noderowmap->GID(i);
    DRT::Node* node = Discret().gNode(nodeGID);
    if (!node) FOUR_C_THROW("Cannot find node with gid %", nodeGID);
    std::vector<int> dofs = Discret().Dof(0, node);
    int numdof = (int)dofs.size();
    int locsysindex = (int)(*locsystoggle_)[i];

    // skip nodes whose dofs have already been processed
    for (int rr = 0; rr < numdof; ++rr)
      if ((*already_processed)[dofrowmap->LID(dofs[rr])] > 1e-9) continue;

    // unity matrix for non-locsys node
    if (locsysindex < 0)
    {
      for (int r = 0; r < numdof; ++r) trafo_->SetValue(1.0, dofs[r], dofs[r]);
    }
    // trafo matrix for locsys node
    else
    {
      CORE::LINALG::SerialDenseMatrix nodetrafo(numdof, numdof);
      for (int k = 0; k < numdof; ++k) nodetrafo(k, k) = 1.0;

      CORE::LINALG::Matrix<3, 1> currrotvector = nodalrotvectors_[nodeGID];
      CORE::LINALG::Matrix<3, 3> currrotationmatrix;

      // Compute rotation matrix out of rotation angle
      CORE::LARGEROTATIONS::angletotriad(currrotvector, currrotationmatrix);

      // base vectors of local system
      CORE::LINALG::Matrix<3, 1> vec1;
      CORE::LINALG::Matrix<3, 1> vec2;
      CORE::LINALG::Matrix<3, 1> vec3;

      // The columns of the rotation matrix are the base vectors
      for (int j = 0; j < 3; j++)
      {
        vec1(j) = currrotationmatrix(j, 0);
        vec2(j) = currrotationmatrix(j, 1);
        vec3(j) = currrotationmatrix(j, 2);
      }

      // Check for zero-diagonal elements
      if (fabs(vec1(0)) < 1e-9 || fabs(vec2(1)) < 1e-9 || fabs(vec3(2)) < 1e-9) sanity_check = true;
      // test how big numdofs are
      if (numdof > 3)
      {
        FOUR_C_THROW(
            "The locsys condition is not implemented for elements with numdof "
            "exceeding 3");
      }
      // trafo for 2D case
      if (Dim() == 2)
      {
        for (int dim = 0; dim < 2; dim++)
        {
          nodetrafo(0, dim) = vec1(dim);
          nodetrafo(1, dim) = vec2(dim);
        }
      }
      // trafo for 3D case
      if (Dim() == 3)
      {
        for (int dim = 0; dim < 3; dim++)
        {
          nodetrafo(0, dim) = vec1(dim);
          nodetrafo(1, dim) = vec2(dim);
          nodetrafo(2, dim) = vec3(dim);
        }
      }

      // Assemble the rotation of this dofs ('nodetrafo') into the global matrix
      for (int r = 0; r < numdof; ++r)
        for (int c = 0; c < numdof; ++c) trafo_->SetValue(nodetrafo(r, c), dofs[r], dofs[c]);

      // store the DOF with locsys
      for (int r = 0; r < numdof; ++r) locsysdofset.insert(dofs[r]);

      // node dofs are marked now as already processed
      for (int rr = 0; rr < numdof; ++rr) (*already_processed)[dofrowmap->LID(dofs[rr])] = 1.0;
    }
  }

  // complete transformation matrix
  trafo_->Complete();

  // Throw warning if transformation matrix has zero diagonal elements since
  // they end up on the diagonal of the system matrix. Show this warning only
  // once.
  if ((not warning_thrown_) && sanity_check)
  {
    if (Comm().MyPID() == 0)
    {
      warning_thrown_ = true;

      printf("Locsys warning:\n");
      printf("A zero diagonal element on the transformation matrix occured.\n");
      printf("This will probably cause a crash in the preconditioner.\n");
      printf("Try not to rotate your local coordinate system by 90 degrees \n");
      printf("or more or use the slow version.\n");
      printf("This warning won't be repeated anymore.\n");
    }
  }

  //**********************************************************************
  // Build map holding DOFs linked to nodes with local co-ordinate system
  //**********************************************************************

  // create unique/row map of DOFs subjected to local co-ordinate change
  // transformation matrix for relevent DOFs with local system
  int nummyentries = 0;
  int* myglobalentries = nullptr;
  std::vector<int> locsysdofs;
  if (locsysdofset.size() > 0)
  {
    locsysdofs.reserve(locsysdofset.size());
    locsysdofs.assign(locsysdofset.begin(), locsysdofset.end());
    nummyentries = static_cast<int>(locsysdofs.size());
    myglobalentries = locsysdofs.data();
  }
  locsysdofmap_ = Teuchos::rcp(new Epetra_Map(
      -1, nummyentries, myglobalentries, discret_.DofRowMap()->IndexBase(), discret_.Comm()));
  if (locsysdofmap_ == Teuchos::null) FOUR_C_THROW("Creation failed.");

  // The matrix subtrafo_ is used in order to apply the Dirichlet Conditions in a more efficient
  // manner
  subtrafo_ = trafo_->ExtractDirichletRows(*locsysdofmap_);

  /*
  REMARK:
   The most general approach to apply Dirichlet conditions in a rotated, local system would be:
   1) Transform the system into local coordinates by means of
      K \cdot D = F --> \tilde{K} \cdot \tilde{D} = \tilde{F}
      with \tilde{K} = trafo_ \cdot K \cdot trafo_^T, \tilde{F} = trafo_ \cdot F, \tilde{D} = trafo_
  \cdot D 2) Apply Dirichlet conditions in the rotated system 3) Transform the system back into
  global coordinates, i.e. \tilde{K} \cdot \tilde{D} = \tilde{F} --> K \cdot D = F with K = trafo_^T
  \cdot \tilde{K} \cdot trafo_, F = trafo_^T \cdot \tilde{F}, D = trafo_^T \cdot \tilde{D}

   Nevertheless, we apply a more efficient algorithm which can be shown,to deliver an equivalent
  system of equations: 1) Therefore we only apply one left transformation to our system of equations
  according K \cdot D = F --> trafo_ \cdot K \cdot D = trafo_ \cdot F 2) Afterwards we apply the
  rotated Dirichlet conditions in an appropriate manner, i.e. we zero the corresponding Dirichlet
  line and than insert the corresponding local base vector vec_i of the assigned local system into
  the corresponding 3*3-block, e.g. if the DoFs of the locsys node are represented by fourth, fifth
  and sixth column:
      (*,*,*,*,*,*,*,*,*,*,*,*) --> (0,0,0, vec_i^T, 0,0,0,0,0,0)
      We don't invert the left transformation of our system afterwards. This means, that we don't
  solve the original but a algebraic manipulated system of equations, nevertheless we still solve
  for the original, non-rotated DoFs D. However, this is actually no drawback since e.g.
  zero-diagonal elements resulting from rotated Dirichlet conditions would still exist even if we
  applied the back transformation afterwards.
  */
}


/*----------------------------------------------------------------------*
 |  print manager (public)                                   meier 06/13|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::Print() const
{
  if (Comm().MyPID() == 0)
  {
    IO::cout << "\n-------------------------------------DRT::UTILS::LocsysManager" << IO::endl;
    for (int i = 0; i < NumLocsys(); ++i)
    {
      IO::cout << "*  *  *  *  *  *  *  *  *  *  *  *  *Locsys entity ID: "
               << locsysconds_[i]->Id();
      if (TypeLocsys(i) == CORE::Conditions::PointLocsys)
        IO::cout << " Point   " << IO::endl;
      else if (TypeLocsys(i) == CORE::Conditions::LineLocsys)
        IO::cout << " Line    " << IO::endl;
      else if (TypeLocsys(i) == CORE::Conditions::SurfaceLocsys)
        IO::cout << " Surface " << IO::endl;
      else if (TypeLocsys(i) == CORE::Conditions::VolumeLocsys)
        IO::cout << " Volume  " << IO::endl;
      else
        FOUR_C_THROW("Unknown type of locsys condition!");
    }
    IO::cout << "-------------------------------------------------------------\n\n";
  }
}


/*----------------------------------------------------------------------*
 |  Get the communicator                                                |
 *----------------------------------------------------------------------*/
inline const Epetra_Comm& DRT::UTILS::LocsysManager::Comm() const { return discret_.Comm(); }


/*----------------------------------------------------------------------*
 |  Transform system global -> local (public)                 popp 09/08|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::RotateGlobalToLocal(
    Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat, Teuchos::RCP<Epetra_Vector> rhs) const
{
  // transform rhs vector
  RotateGlobalToLocal(rhs);

  // selective multiplication from left
  Teuchos::RCP<CORE::LINALG::SparseMatrix> temp =
      CORE::LINALG::Multiply(*subtrafo_, false, *sysmat, false, true);
  // put transformed rows back into global matrix
  sysmat->Put(*temp, 1.0, locsysdofmap_);
}


/*----------------------------------------------------------------------*
 |  Transform system matrix global -> local (public)       mueller 05/10|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::RotateGlobalToLocal(
    Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat) const
{
  // selective multiplication from left
  Teuchos::RCP<CORE::LINALG::SparseMatrix> temp = CORE::LINALG::MLMultiply(
      *subtrafo_, *sysmat, sysmat->ExplicitDirichlet(), sysmat->SaveGraph(), true);

  // put transformed rows back into global matrix
  sysmat->Put(*temp, 1.0, locsysdofmap_);
}

/*----------------------------------------------------------------------*
 |  Transform vector global -> local (public)                 popp 09/08|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::RotateGlobalToLocal(
    Teuchos::RCP<Epetra_Vector> vec, bool offset) const
{
  // y = trafo_ . x  with x = vec
  Epetra_Vector tmp(*vec);
  trafo_->Multiply(false, tmp, *vec);
}

/*----------------------------------------------------------------------*
 |  Transform result + system local -> global (public)        popp 09/08|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::RotateLocalToGlobal(Teuchos::RCP<Epetra_Vector> result,
    Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat, Teuchos::RCP<Epetra_Vector> rhs) const
{
  // transform result
  RotateLocalToGlobal(result);

  // transform rhs vector
  RotateLocalToGlobal(rhs);

  // transform system matrix
  Teuchos::RCP<CORE::LINALG::SparseMatrix> temp;
  Teuchos::RCP<CORE::LINALG::SparseMatrix> temp2;

  // We want to keep the SaveGraph() value of sysmat also after transformation.
  // It is not possible to keep ExplicitDirichlet()==true after transformation,
  // so we explicitly set this to false.
  temp = CORE::LINALG::Multiply(*sysmat, false, *trafo_, false, false, sysmat->SaveGraph(), true);
  temp2 = CORE::LINALG::Multiply(*trafo_, true, *temp, false, false, sysmat->SaveGraph(), true);

  // this is a deep copy (expensive!)
  *sysmat = *temp2;
}

/*----------------------------------------------------------------------*
 |  Transform vector local -> global (public)                 popp 09/08|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::RotateLocalToGlobal(
    Teuchos::RCP<Epetra_Vector> vec, bool offset) const
{
  Epetra_Vector tmp(*vec);
  trafo_->Multiply(true, tmp, *vec);
}
/*----------------------------------------------------------------------*
 |  Transform matrix local -> global (public)              mueller 05/10|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::RotateLocalToGlobal(
    Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat) const
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> temp2 =
      CORE::LINALG::Multiply(*trafo_, true, *sysmat, false, false, sysmat->SaveGraph(), true);
  *sysmat = *temp2;
}


/*----------------------------------------------------------------------*
 |  Calculate rotation vector for (mass-consistent) normal              |
 |  system for a given locsys condition                       hahn 07/14|
 *----------------------------------------------------------------------*/
void DRT::UTILS::LocsysManager::CalcRotationVectorForNormalSystem(int numLocsysCond, double time)
{
  // Take care for "negative times", where no information about dispnp_ is available
  if (time < 0.0)
  {
    Teuchos::RCP<Epetra_Vector> zeroVector =
        CORE::LINALG::CreateVector(*Discret().DofRowMap(), true);
    discret_.SetState("dispnp", zeroVector);
  }

  // Calculate (non-normalized) mass-consistent node normals for a given locsys condition
  // *****+++++++++++++++++**************************************************************

  // check if the normals were set
  if (nodenormals_.empty())
  {
    FOUR_C_THROW(
        "The option massConsistentNodeNormals of the Local SYS BC needs the current normals of the "
        "problem. Before calling Setup, please provide the Node Normals.");
  }

  // Declare node normal variable
  Teuchos::RCP<Epetra_Vector> massConsistentNodeNormals = nodenormals_[numLocsysCond];

  // Loop through all nodes in the condition
  // *******************************************************************
  // Obtain desired locsys condition and its corresponding nodes
  DRT::Condition* currLocsysCond = locsysconds_[numLocsysCond];
  const std::vector<int>* nodes = currLocsysCond->GetNodes();

  // Obtain rank of calling processor
  int myrank = discret_.Comm().MyPID();

  // Loop through all nodes in the condition
  for (int nodeGID : *nodes)
  {
    // Don't care about nodes that don't exist
    bool haveNode = discret_.HaveGlobalNode(nodeGID);
    if (!haveNode) continue;

    DRT::Node* node = discret_.gNode(nodeGID);

    // Don't care about nodes that the processor doesn't own
    bool isOwner = (node->Owner() == myrank);
    if (!isOwner) continue;

    // Obtain all gids of all dofs of the node
    std::vector<int> nodeGIDs;
    discret_.Dof(node, nodeGIDs);

    // Obtain node normal for current node and calculate its vector length
    // *******************************************************************
    CORE::LINALG::Matrix<3, 1> nodeNormal;  // massConsistentNodeNormals contains (dim_+1) dofs
                                            // in the fluid case and (dim_) dofs in the ale case,
                                            // but only the first (dim_) are used.
    nodeNormal.Clear();                     // if dim_==2, then the third component is just not used
    double length = 0.0;
    for (int jdim = 0; jdim < dim_; jdim++)
    {
      const int localId = massConsistentNodeNormals->Map().LID(nodeGIDs[jdim]);
      nodeNormal(jdim, 0) = (*massConsistentNodeNormals)[localId];
      length += nodeNormal(jdim, 0) * nodeNormal(jdim, 0);
    }
    length = sqrt(length);
    if (length < 1e-12)
    {
      FOUR_C_THROW(
          "Locsys: CalcRotationVectorForNormalSystem: Node normal length is zero, what shouldn't "
          "happen! Check, if your BC nodeset really contains surface elements!");
    }

    // Normalize the node normal vector
    for (int jdim = 0; jdim < dim_; jdim++)
    {  // massConsistentNodeNormals contains (dim_+1) dofs
       // in the fluid case and (dim_) dofs in the ale case,
       // but only the first (dim_) are used.
      nodeNormal(jdim, 0) = nodeNormal(jdim, 0) / length;
    }

    // Finally, calculate the rotation vector
    // *******************************************************************
    CORE::LINALG::Matrix<3, 1> currNodalRotVector;

    // Determine rotation angle
    const double rotAngle = acos(nodeNormal(0, 0));

    // Calculate the L2-norm of the rotation vector (which is given by (0, -nodeNormal(2,0),
    // nodeNormal(1,0)))
    const double rotVecNorm =
        sqrt(nodeNormal(1, 0) * nodeNormal(1, 0) + nodeNormal(2, 0) * nodeNormal(2, 0));

    // Calculate the requested rotation vector
    if (rotVecNorm > 1e-12)
    {  // normal vector is not (+-1,0,0), thus rotate as planned
      currNodalRotVector(0, 0) = 0.0;
      currNodalRotVector(1, 0) = rotAngle * (-1) * nodeNormal(2, 0) / rotVecNorm;
      currNodalRotVector(2, 0) = rotAngle * nodeNormal(1, 0) / rotVecNorm;
    }
    else if (nodeNormal(0, 0) < 0)
    {  // normal vector is (-1,0,0), thus rotate 180 deg about z-axis, i.e. (0,0,pi)
      currNodalRotVector(0, 0) = 0.0;
      currNodalRotVector(1, 0) = 0.0;
      currNodalRotVector(2, 0) = M_PI;
    }
    else
    {  // normal vector is (+1,0,0), thus no rotation at all
      currNodalRotVector(0, 0) = 0.0;
      currNodalRotVector(1, 0) = 0.0;
      currNodalRotVector(2, 0) = 0.0;
    }

    // Write out nodal rotation vector to global list
    nodalrotvectors_[nodeGID] = currNodalRotVector;

    // Do some locsys voodoo
    int indices = nodeGID;
    double values = numLocsysCond;
    locsystoggle_->ReplaceGlobalValues(1, &values, &indices);
  }
}

FOUR_C_NAMESPACE_CLOSE

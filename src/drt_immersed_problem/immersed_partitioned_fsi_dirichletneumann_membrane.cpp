/*!----------------------------------------------------------------------
\file immersed_partitioned_fsi_dirichletneumann_membrane.cpp

\brief partitioned immersed fsi algorithm for dirichlet-neumann coupling with membrane finite
elements

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*----------------------------------------------------------------------*/
#include "immersed_partitioned_fsi_dirichletneumann_membrane.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

#include "../drt_io/io_control.H"

#include "../drt_fluid/fluid_utils.H"

IMMERSED::ImmersedPartitionedFSIDirichletNeumannMembrane::
    ImmersedPartitionedFSIDirichletNeumannMembrane(const Epetra_Comm& comm)
    : ImmersedPartitionedFSIDirichletNeumann(comm)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 | calc the current artificial velocity                  sfuchs 11/2016 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedFSIDirichletNeumannMembrane::CalcArtificialVelocity()
{
  if (not artificial_velocity_isvalid_)
  {
    // reinitialize the transfer vector
    fluid_artificial_velocity_->Scale(0.0);

    // declare parameter list and set some params
    Teuchos::ParameterList params;
    params.set<Teuchos::RCP<GEO::SearchTree>>("structsearchtree_rcp", structure_SearchTree_);
    params.set<std::map<int, LINALG::Matrix<3, 1>>*>(
        "currpositions_struct", &currpositions_struct_);
    params.set<std::string>("immerseddisname", "structure");

    // pointer to (potentially) immersed fluid element
    DRT::Element* ele;
    DRT::ELEMENTS::FluidImmersedBase* immersedelebase;

    // set the states needed for evaluation
    SetStatesFluidOP();

    // update search trees, etc. ...
    PrepareFluidOp();

    // access fluid discretization
    Teuchos::RCP<DRT::Discretization> fluid_dis = MBFluidField()->Discretization();

    // vector containing nodal dof gids of immersed boundary fluid elements
    std::vector<int> immerseddofs(0);

    // vector containing gids of immersed boundary fluid elements
    std::vector<int> immersedbdrygids(0);


    if (myrank_ == 0)
    {
      std::cout << "###############################################################################"
                   "#################"
                << std::endl;
      std::cout << "###   Least squares approximation of " << fluid_dis->Name()
                << " nodal velocities from immersed elements..." << std::endl;
      std::cout << "###############################################################################"
                   "#################"
                << std::endl;
    }

    /*===============================================================================*
     | search for fluid elements that are immersed by membrane elements              |
     *===============================================================================*/

    // set fluid action type
    params.set<int>("action", FLD::search_immersed_boundary_elements);

    // loop over all potentially immersed fluid elements
    for (std::map<int, std::set<int>>::const_iterator closele = curr_subset_of_fluiddis_.begin();
         closele != curr_subset_of_fluiddis_.end(); closele++)
    {
      for (std::set<int>::const_iterator eleIter = (closele->second).begin();
           eleIter != (closele->second).end(); eleIter++)
      {
        ele = fluid_dis->gElement(*eleIter);

        immersedelebase = dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(ele);
        if (immersedelebase == NULL)
          dserror("dynamic cast from DRT::Element* to DRT::ELEMENTS::FluidImmersedBase* failed");

        // location array containing fluid dofset
        DRT::Element::LocationArray la(1);
        immersedelebase->LocationVector(*fluid_dis, la, false);

        // dummy matrix and vector
        Epetra_SerialDenseMatrix dummyMat;
        Epetra_SerialDenseVector dummyVec;

        // check if current fluid element contains structure nodes or gauss points
        // and set the flags SetIsImmersed(1) and SetHasProjectedDirichlet(1) for the element
        immersedelebase->Evaluate(
            params, *fluid_dis, la[0].lm_, dummyMat, dummyMat, dummyVec, dummyVec, dummyVec);

        // fluid element is immersed by structure
        if (immersedelebase->IsImmersed())
        {
          // append gid of immersed boundary fluid element
          immersedbdrygids.push_back(immersedelebase->Id());

          // loop over nodes of immersed fluid element
          for (int node = 0; node < (ele->NumNode()); node++)
          {
            // access current node
            DRT::Node* currnode = (ele->Nodes()[node]);

            // node is owned by calling proc
            if (currnode->Owner() == myrank_)
            {
              // dof gids of node are already added
              if (static_cast<IMMERSED::ImmersedNode*>(ele->Nodes()[node])->IsMatched()) continue;

              // otherwise mark node as matched and add dof gids
              static_cast<IMMERSED::ImmersedNode*>(ele->Nodes()[node])->SetIsMatched(1);

              // get nodal dof gids
              std::vector<int> dofs = fluid_dis->Dof(0, currnode);

              // append nodal dof gids
              for (int dim = 0; dim < 4; ++dim)
              {
                immerseddofs.push_back(dofs[dim]);
              }
            }
          }  // end loop over nodes of immersed fluid element
        }
      }
    }  // end loop over all potentially immersed fluid elements

    // output of found immersed elements
    std::cout << "CalcArtificialVelocity found " << immersedbdrygids.size()
              << " boundary immersed elements on Proc " << myrank_ << std::endl;

    /*===============================================================================*
     | get DofRowMaps and initialize the map extractor                               |
     *===============================================================================*/

    // build DofRowMap of immersed fluid nodes
    int numimmerseddofs = immerseddofs.size();
    Teuchos::RCP<Epetra_Map> immersed_dofrowmap =
        Teuchos::rcp(new Epetra_Map(-1, numimmerseddofs, &(immerseddofs[0]), 0, fluid_dis->Comm()));

    if (not immersed_dofrowmap->UniqueGIDs())
      dserror("immersed_dofrowmap in CalcArtificialVelocity is not unique!");

    // get DofRowMap of all fluid nodes
    Teuchos::RCP<Epetra_Map> fluid_dofrowmap =
        Teuchos::rcp(new Epetra_Map(*fluid_dis->DofRowMap(0)));

    // initialize (fluid nodal dof <-> immersed nodal dof) map extractor
    Teuchos::RCP<LINALG::MapExtractor> bdry_map_extractor =
        Teuchos::rcp(new LINALG::MapExtractor(*fluid_dofrowmap,
            immersed_dofrowmap  // is condition map
            ));

    // check immersed boundary map extractor
    bdry_map_extractor->CheckForValidMapExtractor();

    /*===============================================================================*
     | initialize least squares system matrix and vectors for immersed boundary      |
     *===============================================================================*/

    // initialize least squares system matrix
    Teuchos::RCP<LINALG::SparseMatrix> ls_mat =
        Teuchos::rcp(new LINALG::SparseMatrix(*immersed_dofrowmap, 32, true, false));

    // initialize least squares right-hand side vector
    Teuchos::RCP<Epetra_Vector> ls_rhs = LINALG::CreateVector(*immersed_dofrowmap, true);

    // initialize least squares solution vector
    Teuchos::RCP<Epetra_Vector> ls_sol = LINALG::CreateVector(*immersed_dofrowmap, true);

    /*===============================================================================*
     | create assemble strategy                                                      |
     *===============================================================================*/

    // create strategy for assembly of least squares system matrix and right-hand side vector
    DRT::AssembleStrategy ls_strategy(
        0,       // row assembly based on number of dofset associated with fluid discretization
        0,       // column assembly based on number of dofset associated with fluid discretization
        ls_mat,  // least squares system matrix
        Teuchos::null,  //
        ls_rhs,         // least squares right-hand side vector
        Teuchos::null,  //
        Teuchos::null   //
    );

    /*===============================================================================*
     | fill least squares system matrix and vector for immersed boundary             |
     *===============================================================================*/

    // set fluid action type
    params.set<int>("action", FLD::least_squares_matrix_rhs_immersed_boundary);

    // loop over all immersed fluid elements
    for (int i = 0; i < (int)immersedbdrygids.size(); ++i)
    {
      int eid = immersedbdrygids[i];

      ele = fluid_dis->gElement(eid);

      immersedelebase = dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(ele);
      if (immersedelebase == NULL)
        dserror("dynamic cast from DRT::Element* to DRT::ELEMENTS::FluidImmersedBase* failed");

      // location array containing fluid dofset
      DRT::Element::LocationArray la(1);
      immersedelebase->LocationVector(*fluid_dis, la, false);

      // get row and col dofset
      int row = ls_strategy.FirstDofSet();
      int col = ls_strategy.SecondDofSet();

      ls_strategy.ClearElementStorage(la[row].Size(), la[col].Size());

      // evaluate the least squares matrix and rhs for the current fluid element
      immersedelebase->Evaluate(params, *fluid_dis, la[row].lm_, ls_strategy.Elematrix1(),
          ls_strategy.Elematrix2(), ls_strategy.Elevector1(), ls_strategy.Elevector2(),
          ls_strategy.Elevector3());

      // assembly
      ls_strategy.AssembleMatrix1(eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_);
      ls_strategy.AssembleVector1(la[row].lm_, la[row].lmowner_);
    }

    // finalize least squares system matrix
    ls_mat->Complete(*immersed_dofrowmap, *immersed_dofrowmap);

    /*===============================================================================*
     | solve resulting least squares system (saddle point system)                    |
     *===============================================================================*/

    // solver setup
    Teuchos::ParameterList param_solve = DRT::Problem::Instance()->UMFPACKSolverParams();
    bool refactor = true;
    bool reset = true;
    Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(new LINALG::Solver(
        param_solve, fluid_dis->Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));
    fluid_dis->ComputeNullSpaceIfNecessary(solver->Params());

    // solve for least squares optimal nodal values
    solver->Solve(ls_mat->EpetraOperator(), ls_sol, ls_rhs, refactor, reset);

    /*===============================================================================*
     | extract the least squares optimal nodal values                                |
     *===============================================================================*/

    bdry_map_extractor->InsertVector(ls_sol, 1, fluid_artificial_velocity_);

    // we just validated the artificial velocity
    artificial_velocity_isvalid_ = true;

    // we just validated the immersed info again.
    // technically this is not entirely true.
    // see remark in constructor of parent class.
    immersed_info_isvalid_ = true;
  }

  return fluid_artificial_velocity_;
}

/*----------------------------------------------------------------------*
 | improve quality of solution near the interface        sfuchs 11/2016 |
 *----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannMembrane::CorrectInterfaceVelocity()
{
  /*
   *  not implemented yet for immersed_partitioned_fsi_dirichletneumann_membrane
   *  (a strategy in order to prevent oscillations of pressure close to the immersed boundary
   * interface)
   */

  return;
}

/*----------------------------------------------------------------------*
 | calc the fluid tractions on the structure             sfuchs 11/2016 |
 *----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFSIDirichletNeumannMembrane::CalcFluidTractionsOnStructure()
{
  /*
   *  the immersed partitioned fsi with membrane finite elements is in a first step implemented
   * solely for structures with dirichlet b.c. prescribed on all nodes so there is no need to
   * determine the fluid traction
   *
   *  remove this overloaded function after implementing the struct_calc_fluid_traction action type
   * in membrane_evaluate
   */

  return;
}

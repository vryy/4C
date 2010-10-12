
#ifdef STKADAPTIVE

#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../stk_lib/stk_discret.H"

#include "str_genalpha.H"
#include "str_problem.H"

#include "../linalg/linalg_solver.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io.H"

#include "../stk_refine/stk_mesh.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
STK::STR::Problem::Problem()
  : meta_( Teuchos::rcp( new MetaMesh ) )
{

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::STR::Problem::Setup()
{
  Teuchos::RCP<DRT::Discretization> actdis = DRT::Problem::Instance()->Dis( genprob.numsf, 0 );
  if ( not actdis->HaveDofs() )
  {
    actdis->FillComplete();
  }

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  Teuchos::RCP<LINALG::Solver> solver =
    rcp(new LINALG::Solver(DRT::Problem::Instance()->StructSolverParams(),
                           actdis->Comm(),
                           DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  dis_ = Teuchos::rcp( new STK::Discretization( actdis->Comm() ) );

  structure_ = Teuchos::rcp( new STK::STR::Structure( *dis_, solver ) );

  // setup mesh part definitions

  dis_->MetaSetup( meta_, *actdis );

  // declare fields

  structure_->declare_fields( meta_->MetaData() );

  // done with meta data

  meta_->Commit();

  // create uniform mesh

  mesh_ = Teuchos::rcp( new Mesh( *meta_, MPI_COMM_WORLD ) );

  // setup mesh

  dis_->MeshSetup( mesh_, *actdis );
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::STR::Problem::Execute()
{
  Setup();

#if 0
  Teuchos::RCP<IO::DiscretizationWriter> output = Teuchos::rcp(new IO::DiscretizationWriter(actdis));
  output->WriteMesh(0,0.0);
  output->NewStep(0,0.0);
  output->WriteElementData();
#endif

  dis_->AdaptMesh( std::vector<stk::mesh::EntityKey>(),
                   std::vector<stk::mesh::EntityKey>() );

  structure_->Integrate();

  Teuchos::RCP<DRT::Discretization> actdis = DRT::Problem::Instance()->Dis( genprob.numsf, 0 );

  DRT::Problem::Instance()->AddFieldTest(structure_->CreateFieldTest());
  DRT::Problem::Instance()->TestAll(actdis->Comm());
}

#endif

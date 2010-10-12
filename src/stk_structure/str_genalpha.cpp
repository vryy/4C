
#ifdef STKADAPTIVE

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_Time.hpp>

#include "str_genalpha.H"
#include "str_resulttest.H"

#include "../drt_io/io_control.H"

#include "../stk_lib/stk_discret.H"
#include "../stk_refine/stk_mesh.H"
#include "../stk_lib/stk_iterator.H"

#include "../stk_lib/stk_algebra.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_element.H"

#include "../linalg/linalg_fixedsparsematrix.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include "../drt_mat/matpar_parameter.H"
#include "../drt_mat/material.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STK::STR::Structure::Structure( STK::Discretization & dis, Teuchos::RCP<LINALG::Solver> solver )
  : dis_( dis ),
    time_( 0.0 ),
    step_( 0 )
{
  const Teuchos::ParameterList& sdyn     = DRT::Problem::Instance()->StructuralDynamicParams();

  stepmax_ = sdyn.get<int>("NUMSTEP");
  maxtime_ = sdyn.get<double>("MAXTIME");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::STR::Structure::Integrate()
{

  for ( int i=step_; i<stepmax_; ++i )
  {
    //if      (predictor==1) ConstantPredictor();
    //else if (predictor==2) ConsistentPredictor();
    FullNewton();
    UpdateandOutput();
    if (time_>=maxtime_)
      break;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> STK::STR::Structure::CreateFieldTest()
{
  return Teuchos::rcp( new STK::STR::StructureResultTest( *this ) );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::STR::Structure::declare_fields( stk::mesh::MetaData & meta )
{

}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::STR::Structure::collect_unknowns( std::vector<stk::mesh::FieldBase*> & fields )
{

}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::STR::Structure::notify_state_changed()
{

}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::STR::Structure::ConstantPredictor()
{

}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::STR::Structure::ConsistentPredictor()
{

}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::STR::Structure::FullNewton()
{

}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::STR::Structure::UpdateandOutput()
{

}

#endif

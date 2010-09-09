#ifdef STKADAPTIVE

#include "fluid_resulttest.H"
#include "fluid_implicit.H"

#include "../stk_lib/stk_discret.H"
#include "../stk_refine/stk_mesh.H"

#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STK::FLD::FluidResultTest::FluidResultTest(Fluid& fluid)
  : mesh_( fluid.Discretization().GetMesh() )
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::FLD::FluidResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  int dis;
  res.ExtractInt("DIS",dis);
  if (dis != 1)
    dserror("fix me: only one fluid discretization supported for testing");

  int node;
  res.ExtractInt("NODE",node);

  stk::mesh::MetaData & meta = mesh_.MetaData();
  stk::mesh::BulkData & bulk = mesh_.BulkData();

  stk::mesh::VectorField * velnp    = meta.get_field<stk::mesh::VectorField>( "velocity" );
  stk::mesh::ScalarField * pressure = meta.get_field<stk::mesh::ScalarField>( "pressure" );

  stk::mesh::Entity * n = bulk.get_entity( stk::mesh::Node, node );
  if ( n!=NULL and n->owner_rank()==bulk.parallel_rank() )
  {
    double result = 0.;

    std::string position;
    res.ExtractString("POSITION",position);
    if (position=="velx")
    {
      double * data = reinterpret_cast<double*>( stk::mesh::field_data( *velnp , *n ) );
      result = data[0];
    }
    else if (position=="vely")
    {
      double * data = reinterpret_cast<double*>( stk::mesh::field_data( *velnp , *n ) );
      result = data[1];
    }
    else if (position=="velz")
    {
      //if (numdim==2)
      //  dserror("Cannot test result for velz in 2D case.");
      double * data = reinterpret_cast<double*>( stk::mesh::field_data( *velnp , *n ) );
      result = data[2];
    }
    else if (position=="pressure")
    {
      double * data = reinterpret_cast<double*>( stk::mesh::field_data( *pressure , *n ) );
      result = data[0];
    }
    else
    {
      dserror("position '%s' not supported in fluid testing", position.c_str());
    }

    nerr += CompareValues(result, res);
    test_count++;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STK::FLD::FluidResultTest::Match(DRT::INPUT::LineDefinition& res)
{
  return res.HaveNamed("FLUID");
}

#endif

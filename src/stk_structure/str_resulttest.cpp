#ifdef STKADAPTIVE

#include "str_resulttest.H"

#include "../stk_lib/stk_discret.H"
#include "../stk_refine/stk_mesh.H"

#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STK::STR::StructureResultTest::StructureResultTest(Structure& structure)
  : mesh_( structure.Discretization().GetMesh() )
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STK::STR::StructureResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  int dis;
  res.ExtractInt("DIS",dis);
  if (dis != 1)
    dserror("fix me: only one structure discretization supported for testing");

  int node;
  res.ExtractInt("NODE",node);

  stk::mesh::MetaData & meta = mesh_.MetaData();
  stk::mesh::BulkData & bulk = mesh_.BulkData();

  //stk::mesh::VectorField * velnp    = meta.get_field<stk::mesh::VectorField>( "velocity" );
  //stk::mesh::ScalarField * pressure = meta.get_field<stk::mesh::ScalarField>( "pressure" );

  stk::mesh::Entity * n = bulk.get_entity( stk::mesh::Node, node );
  if ( n!=NULL and n->owner_rank()==bulk.parallel_rank() )
  {
    double result = 0.;

    std::string position;
    res.ExtractString("POSITION",position);
#if 0
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
#endif
    {
      dserror("position '%s' not supported in structure testing", position.c_str());
    }

    nerr += CompareValues(result, res);
    test_count++;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STK::STR::StructureResultTest::Match(DRT::INPUT::LineDefinition& res)
{
  return res.HaveNamed("STRUCTURE");
}

#endif

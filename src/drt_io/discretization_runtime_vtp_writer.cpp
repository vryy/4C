/*-----------------------------------------------------------------------------------------------*/
/*!
\file discretization_runtime_vtp_writer.cpp

\brief Write visualization output for a discretization in vtk/vtp format at runtime

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

/* headers */
#include "discretization_runtime_vtp_writer.H"

#include "runtime_vtp_writer.H"

#include "io_control.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_lib/drt_dserror.H"


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
DiscretizationRuntimeVtpWriter::DiscretizationRuntimeVtpWriter() :
    runtime_vtpwriter_( Teuchos::rcp( new RuntimeVtpWriter() ) )
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
DiscretizationRuntimeVtpWriter::Initialize(
    Teuchos::RCP<const DRT::Discretization> discretization,
    unsigned int max_number_timesteps_to_be_written,
    bool write_binary_output )
{
  discretization_ = discretization;

  // determine path of output directory
  const std::string outputfilename( DRT::Problem::Instance()->OutputControlFile()->FileName() );

  size_t pos = outputfilename.find_last_of("/");

  if (pos == outputfilename.npos)
    pos = 0ul;
  else
    pos++;

  const std::string output_directory_path( outputfilename.substr(0ul, pos) );


  runtime_vtpwriter_->Initialize(
      discretization_->Comm().MyPID(),
      discretization_->Comm().NumProc(),
      max_number_timesteps_to_be_written,
      output_directory_path,
      DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix(),
      write_binary_output
      );

  SetGeometryFromParticleDiscretization();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
DiscretizationRuntimeVtpWriter::SetGeometryFromParticleDiscretization()
{
  // Note: this will only work for 'particle' discretizations that have nodes but no elements


  // Todo assume 3D for now
  const unsigned int num_spatial_dimensions = 3;

  // count number of nodes and number for each processor; output is completely independent of
  // the number of processors involved
  unsigned int num_row_nodes = (unsigned int) discretization_->NumMyRowNodes();

  // get and prepare storage for point coordinate values
  std::vector<double>& point_coordinates = runtime_vtpwriter_->GetMutablePointCoordinateVector();
  point_coordinates.clear();
  point_coordinates.reserve( num_spatial_dimensions * num_row_nodes );


  // loop over my nodes and collect the geometry/grid data, i.e. reference positions of nodes
  for (unsigned int inode=0; inode < num_row_nodes; ++inode)
  {
    const DRT::Node* node = discretization_->lRowNode(inode);

    for (unsigned int idim=0; idim<num_spatial_dimensions; ++idim)
      point_coordinates.push_back( node->X()[idim] );
  }


  // safety check
  if ( point_coordinates.size() != num_spatial_dimensions * num_row_nodes )
  {
    dserror("DiscretizationRuntimeVtpWriter expected %d coordinate values, but got %d",
        num_spatial_dimensions * num_row_nodes, point_coordinates.size() );
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
DiscretizationRuntimeVtpWriter::ResetTimeAndTimeStep(
    double time,
    unsigned int timestep)
{
  // Todo allow for independent setting of time/timestep and geometry name
  runtime_vtpwriter_->SetupForNewTimeStepAndGeometry( time, timestep, discretization_->Name() );
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
DiscretizationRuntimeVtpWriter::AppendDofBasedResultDataVector(
    const Teuchos::RCP<Epetra_Vector> & result_data_dofbased,
    unsigned int result_num_dofs_per_node,
    const std::string & resultname )
{
  /* the idea is to transform the given data to a 'point data vector' and append it to the
   * collected solution data vectors by calling AppendVisualizationPointDataVector() */

  std::vector<double> vtp_point_result_data;
  vtp_point_result_data.reserve( result_data_dofbased->MyLength() );

  for ( int lid=0; lid<result_data_dofbased->MyLength(); ++lid )
    vtp_point_result_data.push_back( (*result_data_dofbased)[lid] );

  runtime_vtpwriter_->AppendVisualizationPointDataVector(
      vtp_point_result_data, result_num_dofs_per_node, resultname );
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
DiscretizationRuntimeVtpWriter::AppendNodeBasedResultDataVector(
    const Teuchos::RCP<Epetra_MultiVector> & result_data_nodebased,
    unsigned int result_num_components_per_node,
    const std::string & resultname)
{
  /*  the idea is to transform the given data to a 'point data vector' and append it to the
   *  collected solution data vectors by calling AppendVisualizationPointDataVector()
   */

  // count number of nodes for each processor
  const unsigned int num_row_nodes = (unsigned int) result_data_nodebased->Map().NumMyElements();

  // safety check
  if ( (unsigned int) result_data_nodebased->NumVectors() != result_num_components_per_node )
    dserror("DiscretizationRuntimeVtpWriter: expected Epetra_MultiVector with %d columns but got %d",
        result_num_components_per_node, result_data_nodebased->NumVectors() );


  std::vector<double> vtp_point_result_data;
  vtp_point_result_data.reserve( result_num_components_per_node * num_row_nodes );

  for (unsigned int lid=0; lid<num_row_nodes; ++lid)
  {
    for (unsigned int idf=0; idf<result_num_components_per_node; ++idf)
    {
      Epetra_Vector* column = (*result_data_nodebased)(idf);
      vtp_point_result_data.push_back( (*column)[lid] );
    }
  }

  runtime_vtpwriter_->AppendVisualizationPointDataVector(
      vtp_point_result_data, result_num_components_per_node, resultname );
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
DiscretizationRuntimeVtpWriter::WriteFiles()
{
  runtime_vtpwriter_->WriteFiles();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
DiscretizationRuntimeVtpWriter::WriteCollectionFileOfAllWrittenFiles()
{
  runtime_vtpwriter_->WriteCollectionFileOfAllWrittenFiles(
      DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix() +
      "-" + discretization_->Name() );
}

/*-----------------------------------------------------------------------------------------------*/
/*!
\file discretization_runtime_vtu_writer.cpp

\brief Write visualization output for a discretization in vtk/vtu format at runtime

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

/* headers */
#include "discretization_runtime_vtu_writer.H"

#include "runtime_vtu_writer.H"

#include "io_control.H"

#include "../drt_lib/drt_element_vtk_cell_type_register.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_lib/drt_dserror.H"


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
DiscretizationRuntimeVtuWriter::DiscretizationRuntimeVtuWriter() :
    runtime_vtuwriter_( Teuchos::rcp( new RuntimeVtuWriter() ) )
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
DiscretizationRuntimeVtuWriter::Initialize(
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


  runtime_vtuwriter_->Initialize(
      discretization_->Comm().MyPID(),
      discretization_->Comm().NumProc(),
      max_number_timesteps_to_be_written,
      output_directory_path,
      DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix(),
      write_binary_output
      );

  SetGeometryFromDiscretizationStandard();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
DiscretizationRuntimeVtuWriter::SetGeometryFromDiscretizationStandard()
{
  /*  Note:
   *
   *  this will only work for 'standard' elements, whose discretization can directly be
   *  represented by one of the VTK cell types
   *  (see e.g. http://www.vtk.org/doc/nightly/html/vtkCellType_8h.html)
   *
   *  examples where it does not work:
   *    - beam elements using cubic Hermite polynomials for centerline interpolation
   *    - elements using NURBS
   *
   *  see method SetGeometryFromDiscretizationNonStandard() for how to handle these cases
   */

  // Todo assume 3D for now
  const unsigned int num_spatial_dimensions = 3;

  // count number of nodes and number for each processor; output is completely independent of
  // the number of processors involved
  unsigned int num_row_elements = discretization_->NumMyRowElements();
  unsigned int num_nodes = 0;
  for (unsigned int iele=0; iele<num_row_elements; ++iele)
    num_nodes += discretization_->lRowElement(iele)->NumNode();

  // do not need to store connectivity indices here because we create a
  // contiguous array by the order in which we fill the coordinates (otherwise
  // need to adjust order of filling in the coordinates).
  std::vector<double>& point_coordinates = runtime_vtuwriter_->GetMutablePointCoordinateVector();
  point_coordinates.clear();
  point_coordinates.reserve( num_spatial_dimensions * num_nodes );

  std::vector<uint8_t>& cell_types = runtime_vtuwriter_->GetMutableCellTypeVector();
  cell_types.clear();
  cell_types.reserve(num_row_elements);

  std::vector<int32_t>& cell_offsets = runtime_vtuwriter_->GetMutableCellOffsetVector();
  cell_offsets.clear();
  cell_offsets.reserve(num_row_elements);


  // loop over my elements and collect the geometry/grid data
  unsigned int pointcounter = 0;

  for (unsigned int iele=0; iele<num_row_elements; ++iele)
  {
    const DRT::Element* ele = discretization_->lRowElement(iele);

    cell_types.push_back(
        DRT::ELEMENTS::GetVtkCellTypeFromBaciElementShapeType( ele->Shape() ).first );

    const std::vector<int> & numbering =
        DRT::ELEMENTS::GetVtkCellTypeFromBaciElementShapeType( ele->Shape() ).second;

    const DRT::Node* const* nodes = ele->Nodes();

    for (int inode=0; inode<ele->NumNode(); ++inode)
      for (unsigned int idim=0; idim<num_spatial_dimensions; ++idim)
        point_coordinates.push_back( nodes[ numbering[inode] ]->X()[idim] );

    pointcounter += ele->NumNode();

    cell_offsets.push_back(pointcounter);
  }


  // safety checks
  if ( point_coordinates.size() != num_spatial_dimensions * pointcounter )
  {
    dserror("RuntimeVtuWriter expected %d coordinate values, but got %d",
        num_spatial_dimensions * pointcounter, point_coordinates.size() );
  }

  if ( cell_types.size() != num_row_elements )
  {
    dserror("RuntimeVtuWriter expected %d cell type values, but got %d",
        num_row_elements, cell_types.size() );
  }

  if ( cell_offsets.size() != num_row_elements )
  {
    dserror("RuntimeVtuWriter expected %d cell offset values, but got %d",
        num_row_elements, cell_offsets.size() );
  }

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
DiscretizationRuntimeVtuWriter::ResetTimeAndTimeStep(
    double time,
    unsigned int timestep)
{
  // Todo enable independent setting of time/timestep and geometry name in RuntimeVtuWriter
  runtime_vtuwriter_->SetupForNewTimeStepAndGeometry( time, timestep, discretization_->Name() );
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
DiscretizationRuntimeVtuWriter::AppendDofBasedResultDataVector(
    const Teuchos::RCP<Epetra_Vector> & result_data_dofbased,
    unsigned int result_num_dofs_per_node,
    unsigned int read_result_data_from_dofindex,
    const std::string & resultname)
{
  // Todo
  /*  see PostVtuWriter::WriteDofResultStep() for implementation based on PostField
   *
   *  the idea is to transform the given data to a 'point data vector' and append it to the
   *  collected solution data vectors by calling AppendVisualizationPointDataVector()
   */

  dserror("not implemented yet");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
DiscretizationRuntimeVtuWriter::AppendNodeBasedResultDataVector(
    const Teuchos::RCP<Epetra_MultiVector> & result_data_nodebased,
    unsigned int result_num_components_per_node,
    const std::string & resultname)
{
  // Todo
  /*  see PostVtuWriter::WriteDofResultStep() for implementation based on PostField
   *
   *  the idea is to transform the given data to a 'point data vector' and append it to the
   *  collected solution data vectors by calling AppendVisualizationPointDataVector()
   */

  dserror("not implemented yet");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
DiscretizationRuntimeVtuWriter::AppendElementBasedResultDataVector(
    const Teuchos::RCP<Epetra_MultiVector> & result_data_elementbased,
    unsigned int result_num_components_per_element,
    const std::string & resultname)
{
  // Todo
  /*  see PostVtuWriter::WriteDofResultStep() for implementation based on PostField
   *
   *  the idea is to transform the given data to a 'cell data vector' and append it to the
   *  collected solution data vectors by calling AppendVisualizationCellDataVector()
   */

  dserror("not implemented yet");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
DiscretizationRuntimeVtuWriter::WriteFiles()
{
  runtime_vtuwriter_->WriteFiles();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
DiscretizationRuntimeVtuWriter::WriteCollectionFileOfAllWrittenFiles()
{
  runtime_vtuwriter_->WriteCollectionFileOfAllWrittenFiles(
      DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix() +
      "-" + discretization_->Name() );
}

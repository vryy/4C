/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_discretization_runtime_vtu_writer.cpp

\brief Write visualization output for a beam discretization in vtk/vtu format at runtime

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

/* headers */
#include "beam_discretization_runtime_vtu_writer.H"

#include "../drt_io/runtime_vtu_writer.H"
#include "../drt_io/io_control.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_fixedsizematrix.H"

#include "../drt_beam3/beam3_base.H"

#include "../drt_beaminteraction/periodic_boundingbox.H"
#include "../drt_beaminteraction/beaminteraction_calc_utils.H"


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
BeamDiscretizationRuntimeVtuWriter::BeamDiscretizationRuntimeVtuWriter() :
    runtime_vtuwriter_( Teuchos::rcp( new RuntimeVtuWriter() ) ),
    use_absolute_positions_( true )  // Todo set/reset from outside
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
BeamDiscretizationRuntimeVtuWriter::Initialize(
    Teuchos::RCP<DRT::Discretization> discretization,
    Teuchos::RCP<const Epetra_Vector> const& displacement_state_vector,
    bool use_absolute_positions_for_point_coordinates,
    Teuchos::RCP<const GEO::MESHFREE::BoundingBox> const& periodic_boundingbox,
    unsigned int max_number_timesteps_to_be_written,
    bool write_binary_output )
{
  discretization_ = discretization;

  // determine and store local row indices of all beam elements in the given discretization
  for ( unsigned int iele = 0; iele < (unsigned int) discretization_->NumMyRowElements(); ++iele )
  {
    const DRT::Element* ele = discretization_->lRowElement(iele);

    // check for beam element
    const DRT::ELEMENTS::Beam3Base* beamele = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);

    if ( beamele != NULL )
      local_row_indices_beam_elements_.push_back(iele);
  }


  use_absolute_positions_ = use_absolute_positions_for_point_coordinates;

  periodic_boundingbox_ = periodic_boundingbox;
  num_cells_per_element_.resize( local_row_indices_beam_elements_.size() );

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

  SetGeometryFromBeamDiscretization( displacement_state_vector );
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
BeamDiscretizationRuntimeVtuWriter::SetGeometryFromBeamDiscretization(
    Teuchos::RCP<const Epetra_Vector> const& displacement_state_vector )
{
  /*  Note:
   *
   *  The centerline geometry of one element cannot be represented by one simple vtk cell type
   *  because we use cubic Hermite polynomials for the interpolation of the centerline.
   *
   *  Instead, we subdivide each beam element in several linear segments. This corresponds to a
   *  VTK_POLY_LINE (vtk cell type number 4). So one beam element will be visualized as one vtk cell,
   *  but the number of points does not equal the number of FE nodes.
   *
   *  For a list of vtk cell types, see e.g.
   *  http://www.vtk.org/doc/nightly/html/vtkCellType_8h.html
   *
   *  Things get more complicated for 'cut' elements, i.e. when nodes have been 'shifted' due to
   *  periodic boundary conditions. Our approach here is to create two (or more) vtk cells from
   *  one beam element.
   *
   *
   *  Another approach would be to visualize the cubic Hermite polynomials as cubic Lagrange
   *  polynomials and specify four vtk points from e.g. the two FE nodes and two more arbitrary
   *  points along the centerline.
   *  However, the representation of nonlinear geometry in Paraview turned out to not work as
   *  expected (e.g. no visible refinement of subsegments if corresponding parameter is changed).
   */

  // always use 3D for beams
  const unsigned int num_spatial_dimensions = 3;
  num_cells_per_element_.clear();

  // determine number of row BEAM elements for each processor
  // output is completely independent of the number of processors involved
  unsigned int num_beam_row_elements = local_row_indices_beam_elements_.size();
  unsigned int num_vtk_points = num_beam_row_elements * ( BEAMSVTUVISUALSUBSEGMENTS + 1 );

  // do not need to store connectivity indices here because we create a
  // contiguous array by the order in which we fill the coordinates (otherwise
  // need to adjust order of filling in the coordinates).
  std::vector<double>& point_coordinates = runtime_vtuwriter_->GetMutablePointCoordinateVector();
  point_coordinates.clear();
  point_coordinates.reserve( num_spatial_dimensions * num_vtk_points );

  std::vector<uint8_t>& cell_types = runtime_vtuwriter_->GetMutableCellTypeVector();
  cell_types.clear();
  cell_types.reserve(num_beam_row_elements);

  std::vector<int32_t>& cell_offsets = runtime_vtuwriter_->GetMutableCellOffsetVector();
  cell_offsets.clear();
  cell_offsets.reserve(num_beam_row_elements);

  // loop over my elements and collect the geometry/grid data
  unsigned int pointcounter = 0;

  for ( unsigned int ibeamele = 0; ibeamele < num_beam_row_elements; ++ibeamele )
  {
    const DRT::Element* ele = discretization_->lRowElement( local_row_indices_beam_elements_[ibeamele] );

    // cast to beam element
    const DRT::ELEMENTS::Beam3Base* beamele = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);

    // Todo safety check for now, may be removed when better tested
    if ( beamele == NULL )
      dserror("BeamDiscretizationRuntimeVtuWriter expects a beam element here!");


    std::vector<double> beamelement_displacement_vector;

    if ( use_absolute_positions_ )
    {
      BEAMINTERACTION::UTILS::GetCurrentElementDis(
          *discretization_, ele, displacement_state_vector, beamelement_displacement_vector );
    }

    /* loop over the chosen visualization points (equidistant distribution in the element
     * parameter space xi \in [-1,1] ) and determine their interpolated (initial) positions r */
    LINALG::Matrix<3,1> interpolated_position(true);
    LINALG::Matrix<3,1> interpolated_position_priorpoint(true);
    double xi = 0.0;

    for ( unsigned int ipoint = 0; ipoint < BEAMSVTUVISUALSUBSEGMENTS + 1; ++ipoint )
    {
      interpolated_position.Clear();
      xi = -1.0 + ipoint * 2.0 / BEAMSVTUVISUALSUBSEGMENTS;

      if ( use_absolute_positions_ )
        beamele->GetPosAtXi( interpolated_position, xi, beamelement_displacement_vector );
      else
        beamele->GetRefPosAtXi( interpolated_position, xi );

      periodic_boundingbox_->Shift3D( interpolated_position );
      LINALG::Matrix<3,1> unshift_interpolated_position = interpolated_position;

      // if there is a shift between tow consecutive points, double that point and create new cell
      // not for first and last point
      if( ipoint != 0 and periodic_boundingbox_->UnShift3D( unshift_interpolated_position, interpolated_position_priorpoint )
          and ipoint != BEAMSVTUVISUALSUBSEGMENTS)
      {
        for ( unsigned int idim = 0; idim < num_spatial_dimensions; ++idim )
          point_coordinates.push_back( unshift_interpolated_position(idim) );

        ++pointcounter;
        cell_offsets.push_back( pointcounter );
        cell_types.push_back(4);
        ++num_cells_per_element_[ibeamele];
      }

      // in case of last visualization point, we only add the unshifted (compared to former point) configuration
      if ( ipoint == BEAMSVTUVISUALSUBSEGMENTS )
      {
        for ( unsigned int idim = 0; idim < num_spatial_dimensions; ++idim )
          point_coordinates.push_back( unshift_interpolated_position(idim) );
      }
      else
      {
        for ( unsigned int idim = 0; idim < num_spatial_dimensions; ++idim )
          point_coordinates.push_back( interpolated_position(idim) );
      }

      ++pointcounter;

      interpolated_position_priorpoint = interpolated_position;
    }
    // VTK_POLY_LINE (vtk cell type number 4)
    cell_types.push_back(4);
    ++num_cells_per_element_[ibeamele];
    cell_offsets.push_back( pointcounter );

  }

  // safety checks
  if ( cell_types.size() != cell_offsets.size() )
  {
    dserror("RuntimeVtuWriter expected %d cell type values, but got %d",
        num_beam_row_elements, cell_types.size() );
  }

  if ( !periodic_boundingbox_->HavePBC() and
       ( point_coordinates.size() != num_spatial_dimensions * num_vtk_points ) )
  {
    dserror("RuntimeVtuWriter expected %d coordinate values, but got %d",
        num_spatial_dimensions * num_vtk_points, point_coordinates.size() );
  }

  if ( !periodic_boundingbox_->HavePBC() and ( cell_offsets.size() != num_beam_row_elements ) )
  {
    dserror("RuntimeVtuWriter expected %d cell offset values, but got %d",
        num_beam_row_elements, cell_offsets.size() );
  }

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
BeamDiscretizationRuntimeVtuWriter::ResetTimeAndTimeStep(
    double time,
    unsigned int timestep)
{
  // Todo enable independent setting of time/timestep and geometry name in RuntimeVtuWriter
  runtime_vtuwriter_->SetupForNewTimeStepAndGeometry(
      time, timestep, (discretization_->Name() + "-beams") );
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
BeamDiscretizationRuntimeVtuWriter::AppendDisplacementField(
    Teuchos::RCP<const Epetra_Vector> const& displacement_state_vector)
{
  // triads only make sense in 3D
  const unsigned int num_spatial_dimensions = 3;

  // determine number of row BEAM elements for each processor
  // output is completely independent of the number of processors involved
  unsigned int num_beam_row_elements = local_row_indices_beam_elements_.size();
  unsigned int num_vtk_points = num_beam_row_elements * ( BEAMSVTUVISUALSUBSEGMENTS + 1 );
  std::vector<int32_t> const& cell_offsets = runtime_vtuwriter_->GetMutableCellOffsetVector();

  // disp vector
  std::vector<double> displacement_vector;
  displacement_vector.reserve( num_spatial_dimensions * num_vtk_points );

  // number of points so far
  int points_sofar = 0;

  // loop over myrank's beam elements and compute disp for each visualization point
  for ( unsigned int ibeamele = 0; ibeamele < num_beam_row_elements; ++ibeamele )
  {
    const DRT::Element* ele = discretization_->lRowElement( local_row_indices_beam_elements_[ibeamele] );

    // cast to beam element
    const DRT::ELEMENTS::Beam3Base* beamele = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);

    // Todo safety check for now, may be removed when better tested
    if ( beamele == NULL )
      dserror("BeamDiscretizationRuntimeVtuWriter expects a beam element here!");


    // get the displacement state vector for this element
    std::vector<double> beamelement_displacement_vector;

    BEAMINTERACTION::UTILS::GetCurrentElementDis(
        *discretization_, ele, displacement_state_vector, beamelement_displacement_vector );

    /* loop over the chosen visualization points (equidistant distribution in the element
     * parameter space xi \in [-1,1] ) and determine its disp state */
    LINALG::Matrix<3,1> pos_visualization_point;
    LINALG::Matrix<3,1> refpos_visualization_point;
    double xi = 0.0;

    for ( unsigned int ipoint = 0; ipoint < BEAMSVTUVISUALSUBSEGMENTS + 1; ++ipoint )
    {
      xi = -1.0 + ipoint * 2.0 / BEAMSVTUVISUALSUBSEGMENTS;

      pos_visualization_point.Clear();
      refpos_visualization_point.Clear();

      // interpolate
      beamele->GetRefPosAtXi( refpos_visualization_point, xi );
      beamele->GetPosAtXi( pos_visualization_point, xi, beamelement_displacement_vector );

      // in case of periodic boundary conditions, a point (except first and last point of an element) can exists twice,
      // we check this here by looking if current point is in cell offset list and therefore starts of a new cell
      unsigned int num_point_exists = 1;
      if ( ipoint != 0 and ipoint != BEAMSVTUVISUALSUBSEGMENTS )
      {
        unsigned int curr_point_number = points_sofar + 1;
        if ( std::find( cell_offsets.begin(), cell_offsets.end(), curr_point_number ) != cell_offsets.end() )
          num_point_exists = 2;
      }

      // store the information in vectors that can be interpreted by vtu writer (disp = pos - refpos)
      // and update number of point data written
      for( unsigned int i = 0; i < num_point_exists; ++i )
      {
        ++points_sofar;
        for ( unsigned int idim = 0; idim < num_spatial_dimensions; ++idim )
          displacement_vector.push_back( pos_visualization_point( idim, 0 ) - refpos_visualization_point( idim, 0 ) );
      }
    }
  }

  // finally append the solution vectors to the visualization data of the vtu writer object
  runtime_vtuwriter_->AppendVisualizationPointDataVector(
      displacement_vector, num_spatial_dimensions, "displacement" );

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
BeamDiscretizationRuntimeVtuWriter::AppendTriadField(
    Teuchos::RCP<const Epetra_Vector> const& displacement_state_vector)
{
  // triads only make sense in 3D
  const unsigned int num_spatial_dimensions = 3;

  // determine number of row BEAM elements for each processor
  // output is completely independent of the number of processors involved
  unsigned int num_beam_row_elements = local_row_indices_beam_elements_.size();
  unsigned int num_vtk_points = num_beam_row_elements * ( BEAMSVTUVISUALSUBSEGMENTS + 1 );
  std::vector<int32_t> const& cell_offsets = runtime_vtuwriter_->GetMutableCellOffsetVector();

  // we write the triad field as three base vectors at every visualization point
  std::vector<double> base_vector_1;
  base_vector_1.reserve( num_spatial_dimensions * num_vtk_points );

  std::vector<double> base_vector_2;
  base_vector_2.reserve( num_spatial_dimensions * num_vtk_points );

  std::vector<double> base_vector_3;
  base_vector_3.reserve( num_spatial_dimensions * num_vtk_points );

  // number of points so far
  int points_sofar = 0;

  // loop over my elements and collect the data about triads/base vectors
  for ( unsigned int ibeamele = 0; ibeamele < num_beam_row_elements; ++ibeamele )
  {
    const DRT::Element* ele = discretization_->lRowElement( local_row_indices_beam_elements_[ibeamele] );

    // cast to beam element
    const DRT::ELEMENTS::Beam3Base* beamele = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);

    // Todo safety check for now, may be removed when better tested
    if ( beamele == NULL )
      dserror("BeamDiscretizationRuntimeVtuWriter expects a beam element here!");


    // get the displacement state vector for this element
    std::vector<double> beamelement_displacement_vector;

    BEAMINTERACTION::UTILS::GetCurrentElementDis(
        *discretization_, ele, displacement_state_vector, beamelement_displacement_vector );


    /* loop over the chosen visualization points (equidistant distribution in the element
     * parameter space xi \in [-1,1] ) and determine the triad */
    LINALG::Matrix<3,3> triad_visualization_point;
    double xi = 0.0;

    for ( unsigned int ipoint=0; ipoint<BEAMSVTUVISUALSUBSEGMENTS+1; ++ipoint )
    {
      xi = -1.0 + ipoint*2.0/BEAMSVTUVISUALSUBSEGMENTS;

      triad_visualization_point.Clear();

      beamele->GetTriadAtXi( triad_visualization_point, xi, beamelement_displacement_vector );

      // in case of periodic boundary conditions, a point (except first and last point of an element) can exists twice,
      // we check this here by looking if current point is in cell offset list and therefore starts of a new cell
      unsigned int num_point_exists = 1;
      if ( ipoint != 0 and ipoint != BEAMSVTUVISUALSUBSEGMENTS )
      {
        unsigned int curr_point_number = points_sofar + 1;
        if ( std::find( cell_offsets.begin(), cell_offsets.end(), curr_point_number ) != cell_offsets.end() )
          num_point_exists = 2;
      }

      // store the information in vectors that can be interpreted by vtu writer
      // and update number of point data written
      for( unsigned int i = 0; i < num_point_exists; ++i )
      {
        ++points_sofar;
        for ( unsigned int idim = 0; idim < num_spatial_dimensions; ++idim )
        {
          // first column: first base vector
          base_vector_1.push_back( triad_visualization_point(idim,0) );

          // second column: second base vector
          base_vector_2.push_back( triad_visualization_point(idim,1) );

          // third column: third base vector
          base_vector_3.push_back( triad_visualization_point(idim,2) );
        }
      }
    }
  }

  // finally append the solution vectors to the visualization data of the vtu writer object
  runtime_vtuwriter_->AppendVisualizationPointDataVector(
      base_vector_1, num_spatial_dimensions, "base_vector_1" );

  runtime_vtuwriter_->AppendVisualizationPointDataVector(
      base_vector_2, num_spatial_dimensions, "base_vector_2" );

  runtime_vtuwriter_->AppendVisualizationPointDataVector(
      base_vector_3, num_spatial_dimensions, "base_vector_3" );

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
BeamDiscretizationRuntimeVtuWriter::AppendElementOwningProcessor()
{
  // determine number of row BEAM elements for each processor
  // output is completely independent of the number of processors involved
  unsigned int num_beam_row_elements = local_row_indices_beam_elements_.size();

  // processor owning the element
  std::vector<double> owner;
  owner.reserve( num_beam_row_elements );

  // loop over my elements and collect the data about triads/base vectors
  for ( unsigned int ibeamele = 0; ibeamele < num_beam_row_elements; ++ibeamele )
  {
    const DRT::Element* ele = discretization_->lRowElement( local_row_indices_beam_elements_[ibeamele] );

    // cast to beam element
    const DRT::ELEMENTS::Beam3Base* beamele = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);

    // Todo safety check for now, may be removed when better tested
    if ( beamele == NULL )
      dserror("BeamDiscretizationRuntimeVtuWriter expects a beam element here!");

    for( int i = 0; i < num_cells_per_element_[ibeamele]; ++i )
      owner.push_back( ele->Owner() );
  }

  // append the solution vector to the visualization data of the vtu writer object
  runtime_vtuwriter_->AppendVisualizationCellDataVector(
      owner, 1, "element_owner" );

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
BeamDiscretizationRuntimeVtuWriter::AppendElementCircularCrossSectionRadius()
{
  // Todo we assume a circular cross-section shape here; generalize this to other shapes

  // determine number of row BEAM elements for each processor
  // output is completely independent of the number of processors involved
  unsigned int num_beam_row_elements = local_row_indices_beam_elements_.size();

  // we assume a constant cross-section radius over the element length
  std::vector<double> cross_section_radius;
  cross_section_radius.reserve( num_beam_row_elements );

  // loop over my elements and collect the data about triads/base vectors
  for ( unsigned int ibeamele = 0; ibeamele < num_beam_row_elements; ++ibeamele )
  {
    const DRT::Element* ele = discretization_->lRowElement( local_row_indices_beam_elements_[ibeamele] );

    // cast to beam element
    const DRT::ELEMENTS::Beam3Base* beamele = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);

    // Todo safety check for now, may be removed when better tested
    if ( beamele == NULL )
      dserror("BeamDiscretizationRuntimeVtuWriter expects a beam element here!");


    // this needs to be done for all cells that make up a cut element
    for( int i = 0; i < num_cells_per_element_[ibeamele]; ++i )
      cross_section_radius.push_back( beamele->GetCircularCrossSectionRadiusForInteractions() );
  }

  // append the solution vector to the visualization data of the vtu writer object
  runtime_vtuwriter_->AppendVisualizationCellDataVector(
      cross_section_radius, 1, "cross_section_radius" );

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
BeamDiscretizationRuntimeVtuWriter::AppendPointCircularCrossSectionInformationVector(
    Teuchos::RCP<const Epetra_Vector> const& displacement_state_vector)
{
  // assume 3D here
  const unsigned int num_spatial_dimensions = 3;

  // determine number of row BEAM elements for each processor
  // output is completely independent of the number of processors involved
  unsigned int num_beam_row_elements = local_row_indices_beam_elements_.size();
  unsigned int num_vtk_points = num_beam_row_elements * ( BEAMSVTUVISUALSUBSEGMENTS + 1 );


  // a beam with circular cross-section can be visualized as a 'chain' of straight cylinders
  // this is also supported as 'tube' in Paraview
  // to define one cylinder, we use the first base vector as its unit normal vector
  // and scale it with the cross-section radius of the beam
  // Edit: This approach seems not to work with Paraview because the functionality 'Vary Radius'
  //       of the Tube filter is different to what we expected.
  //       However, we keep this method as it could be useful for other visualization approaches
  //       in the future.
  std::vector<double> circular_cross_section_information_vector;
  circular_cross_section_information_vector.reserve( num_spatial_dimensions * num_vtk_points );
  std::vector<int32_t> const& cell_offsets = runtime_vtuwriter_->GetMutableCellOffsetVector();

  // number of points so far
  int points_sofar = 0;

  // loop over my elements and collect the data about triads/base vectors
  for ( unsigned int ibeamele = 0; ibeamele < num_beam_row_elements; ++ibeamele )
  {
    const DRT::Element* ele = discretization_->lRowElement( local_row_indices_beam_elements_[ibeamele] );

    // cast to beam element
    const DRT::ELEMENTS::Beam3Base* beamele = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);

    // Todo safety check for now, may be removed when better tested
    if ( beamele == NULL )
      dserror("BeamDiscretizationRuntimeVtuWriter expects a beam element here!");


    const double circular_cross_section_radius =
        beamele->GetCircularCrossSectionRadiusForInteractions();

    // get the displacement state vector for this element
    std::vector<double> beamelement_displacement_vector;

    BEAMINTERACTION::UTILS::GetCurrentElementDis(
        *discretization_, ele, displacement_state_vector, beamelement_displacement_vector );


    /* loop over the chosen visualization points (equidistant distribution in the element
     * parameter space xi \in [-1,1] ) and determine the triad */
    LINALG::Matrix<3,3> triad_visualization_point;
    double xi = 0.0;
    for ( unsigned int ipoint = 0; ipoint < BEAMSVTUVISUALSUBSEGMENTS + 1; ++ipoint )
    {
      xi = -1.0 + ipoint * 2.0 / BEAMSVTUVISUALSUBSEGMENTS;

      triad_visualization_point.Clear();

      beamele->GetTriadAtXi( triad_visualization_point, xi, beamelement_displacement_vector );

      // in case of periodic boundary conditions, a point (except first and last point of an element) can exists twice,
      // we check this here by looking if current point is in cell offset list and therefore starts of a new cell
      unsigned int num_point_exists = 1;
      if ( ipoint != 0 and ipoint != BEAMSVTUVISUALSUBSEGMENTS )
      {
        unsigned int curr_point_number = points_sofar + 1;
        if ( std::find( cell_offsets.begin(), cell_offsets.end(), curr_point_number ) != cell_offsets.end() )
          num_point_exists = 2;
      }

      // store the information in vectors that can be interpreted by vtu writer (disp = pos - refpos)
      // and update number of point data written
      for( unsigned int i = 0; i < num_point_exists; ++i )
      {
        ++points_sofar;
        for ( unsigned int idim = 0; idim < num_spatial_dimensions; ++idim )
        {
          // first column: first base vector
          circular_cross_section_information_vector.push_back(
              triad_visualization_point(idim,0) * circular_cross_section_radius );
        }
      }
    }
  }

  // finally append the solution vectors to the visualization data of the vtu writer object
  runtime_vtuwriter_->AppendVisualizationPointDataVector(
      circular_cross_section_information_vector, num_spatial_dimensions,
      "circular_cross_section_information_vector" );

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
BeamDiscretizationRuntimeVtuWriter::AppendGaussPointMaterialCrossSectionStrains()
{
  // determine number of row BEAM elements for each processor
  // output is completely independent of the number of processors involved
  unsigned int num_beam_row_elements = local_row_indices_beam_elements_.size();


  // storage for material strain measures at all GPs of all my row elements
  std::vector<double> axial_strain_GPs_all_row_elements;
  std::vector<double> shear_strain_2_GPs_all_row_elements;
  std::vector<double> shear_strain_3_GPs_all_row_elements;

  std::vector<double> twist_GPs_all_row_elements;
  std::vector<double> curvature_2_GPs_all_row_elements;
  std::vector<double> curvature_3_GPs_all_row_elements;


  // storage for material strain measures at all GPs of current element
  std::vector<double> axial_strain_GPs_current_element;
  std::vector<double> shear_strain_2_GPs_current_element;
  std::vector<double> shear_strain_3_GPs_current_element;

  std::vector<double> twist_GPs_current_element;
  std::vector<double> curvature_2_GPs_current_element;
  std::vector<double> curvature_3_GPs_current_element;


  // number of Gauss points must be the same for all elements in the grid
  unsigned int num_GPs_per_element_strains_translational = 0;
  unsigned int num_GPs_per_element_strains_rotational = 0;


  // loop over my elements and collect the data
  for ( unsigned int ibeamele = 0; ibeamele < num_beam_row_elements; ++ibeamele )
  {
    const DRT::Element* ele = discretization_->lRowElement( local_row_indices_beam_elements_[ibeamele] );

    // cast to beam element
    const DRT::ELEMENTS::Beam3Base* beamele = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);

    // Todo safety check for now, may be removed when better tested
    if ( beamele == NULL )
      dserror("BeamDiscretizationRuntimeVtuWriter expects a beam element here!");


    axial_strain_GPs_current_element.clear();
    shear_strain_2_GPs_current_element.clear();
    shear_strain_3_GPs_current_element.clear();

    twist_GPs_current_element.clear();
    curvature_2_GPs_current_element.clear();
    curvature_3_GPs_current_element.clear();


    // get GP strain values from previous element evaluation call
    beamele->GetMaterialCrossSectionStrainsAtAllGPs(
        axial_strain_GPs_current_element,
        shear_strain_2_GPs_current_element,
        shear_strain_3_GPs_current_element,
        twist_GPs_current_element,
        curvature_2_GPs_current_element,
        curvature_3_GPs_current_element);


    // safety check for number of Gauss points per element
    // initialize numbers from first element
    if ( ibeamele == 0 )
    {
      num_GPs_per_element_strains_translational = axial_strain_GPs_current_element.size();
      num_GPs_per_element_strains_rotational = twist_GPs_current_element.size();
    }

    if ( axial_strain_GPs_current_element.size() != num_GPs_per_element_strains_translational or
         shear_strain_2_GPs_current_element.size() != num_GPs_per_element_strains_translational or
         shear_strain_3_GPs_current_element.size() != num_GPs_per_element_strains_translational or
         twist_GPs_current_element.size() != num_GPs_per_element_strains_rotational or
         curvature_2_GPs_current_element.size() != num_GPs_per_element_strains_rotational or
         curvature_3_GPs_current_element.size() != num_GPs_per_element_strains_rotational )
    {
      dserror("number of Gauss points must be the same for all elements in discretization!");
    }

    // store the values of current element in the large vectors collecting data from all elements
    for( int i = 0; i < num_cells_per_element_[ibeamele]; ++i )
    {
      InsertVectorValuesAtBackOfOtherVector(
          axial_strain_GPs_current_element, axial_strain_GPs_all_row_elements);

      InsertVectorValuesAtBackOfOtherVector(
          shear_strain_2_GPs_current_element, shear_strain_2_GPs_all_row_elements);

      InsertVectorValuesAtBackOfOtherVector(
          shear_strain_3_GPs_current_element, shear_strain_3_GPs_all_row_elements);


      InsertVectorValuesAtBackOfOtherVector(
          twist_GPs_current_element, twist_GPs_all_row_elements);

      InsertVectorValuesAtBackOfOtherVector(
          curvature_2_GPs_current_element, curvature_2_GPs_all_row_elements);

      InsertVectorValuesAtBackOfOtherVector(
          curvature_3_GPs_current_element, curvature_3_GPs_all_row_elements);
    }
  }

  // append the solution vectors to the visualization data of the vtu writer object
  runtime_vtuwriter_->AppendVisualizationCellDataVector(
      axial_strain_GPs_all_row_elements,
      num_GPs_per_element_strains_translational,
      "axial_strain_GPs" );

  runtime_vtuwriter_->AppendVisualizationCellDataVector(
      shear_strain_2_GPs_all_row_elements,
      num_GPs_per_element_strains_translational,
      "shear_strain_2_GPs" );

  runtime_vtuwriter_->AppendVisualizationCellDataVector(
      shear_strain_3_GPs_all_row_elements,
      num_GPs_per_element_strains_translational,
      "shear_strain_3_GPs" );


  runtime_vtuwriter_->AppendVisualizationCellDataVector(
      twist_GPs_all_row_elements,
      num_GPs_per_element_strains_rotational,
      "twist_GPs" );

  runtime_vtuwriter_->AppendVisualizationCellDataVector(
      curvature_2_GPs_all_row_elements,
      num_GPs_per_element_strains_rotational,
      "curvature_2_GPs" );

  runtime_vtuwriter_->AppendVisualizationCellDataVector(
      curvature_3_GPs_all_row_elements,
      num_GPs_per_element_strains_rotational,
      "curvature_3_GPs" );

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
BeamDiscretizationRuntimeVtuWriter::AppendGaussPointMaterialCrossSectionStresses()
{
  // determine number of row BEAM elements for each processor
  // output is completely independent of the number of processors involved
  unsigned int num_beam_row_elements = local_row_indices_beam_elements_.size();


  // storage for material strain measures at all GPs of all my row elements
  std::vector<double> axial_stress_GPs_all_row_elements;
  std::vector<double> shear_stress_2_GPs_all_row_elements;
  std::vector<double> shear_stress_3_GPs_all_row_elements;

  std::vector<double> torque_GPs_all_row_elements;
  std::vector<double> bending_moment_2_GPs_all_row_elements;
  std::vector<double> bending_moment_3_GPs_all_row_elements;


  // storage for material stress measures at all GPs of current element
  std::vector<double> axial_stress_GPs_current_element;
  std::vector<double> shear_stress_2_GPs_current_element;
  std::vector<double> shear_stress_3_GPs_current_element;

  std::vector<double> torque_GPs_current_element;
  std::vector<double> bending_moment_2_GPs_current_element;
  std::vector<double> bending_moment_3_GPs_current_element;


  // number of Gauss points must be the same for all elements in the grid
  unsigned int num_GPs_per_element_stresses_translational = 0;
  unsigned int num_GPs_per_element_stresses_rotational = 0;


  // loop over my elements and collect the data
  for ( unsigned int ibeamele = 0; ibeamele < num_beam_row_elements; ++ibeamele )
  {
    const DRT::Element* ele = discretization_->lRowElement( local_row_indices_beam_elements_[ibeamele] );

    // cast to beam element
    const DRT::ELEMENTS::Beam3Base* beamele = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);

    // Todo safety check for now, may be removed when better tested
    if ( beamele == NULL )
      dserror("BeamDiscretizationRuntimeVtuWriter expects a beam element here!");


    axial_stress_GPs_current_element.clear();
    shear_stress_2_GPs_current_element.clear();
    shear_stress_3_GPs_current_element.clear();

    torque_GPs_current_element.clear();
    bending_moment_2_GPs_current_element.clear();
    bending_moment_3_GPs_current_element.clear();


    // get GP stress values from previous element evaluation call
    beamele->GetMaterialCrossSectionStressesAtAllGPs(
        axial_stress_GPs_current_element,
        shear_stress_2_GPs_current_element,
        shear_stress_3_GPs_current_element,
        torque_GPs_current_element,
        bending_moment_2_GPs_current_element,
        bending_moment_3_GPs_current_element);


    // safety check for number of Gauss points per element
    // initialize numbers from first element
    if ( ibeamele == 0 )
    {
      num_GPs_per_element_stresses_translational = axial_stress_GPs_current_element.size();
      num_GPs_per_element_stresses_rotational = torque_GPs_current_element.size();
    }

    if ( axial_stress_GPs_current_element.size() != num_GPs_per_element_stresses_translational or
         shear_stress_2_GPs_current_element.size() != num_GPs_per_element_stresses_translational or
         shear_stress_3_GPs_current_element.size() != num_GPs_per_element_stresses_translational or
         torque_GPs_current_element.size() != num_GPs_per_element_stresses_rotational or
         bending_moment_2_GPs_current_element.size() != num_GPs_per_element_stresses_rotational or
         bending_moment_3_GPs_current_element.size() != num_GPs_per_element_stresses_rotational )
    {
      dserror("number of Gauss points must be the same for all elements in discretization!");
    }

    // store the values of current element in the large vectors collecting data from all elements
    for( int i = 0; i < num_cells_per_element_[ibeamele]; ++i )
    {
      InsertVectorValuesAtBackOfOtherVector(
          axial_stress_GPs_current_element, axial_stress_GPs_all_row_elements);

      InsertVectorValuesAtBackOfOtherVector(
          shear_stress_2_GPs_current_element, shear_stress_2_GPs_all_row_elements);

      InsertVectorValuesAtBackOfOtherVector(
          shear_stress_3_GPs_current_element, shear_stress_3_GPs_all_row_elements);


      InsertVectorValuesAtBackOfOtherVector(
          torque_GPs_current_element, torque_GPs_all_row_elements);

      InsertVectorValuesAtBackOfOtherVector(
          bending_moment_2_GPs_current_element, bending_moment_2_GPs_all_row_elements);

      InsertVectorValuesAtBackOfOtherVector(
          bending_moment_3_GPs_current_element, bending_moment_3_GPs_all_row_elements);
    }
  }

  // append the solution vectors to the visualization data of the vtu writer object
  runtime_vtuwriter_->AppendVisualizationCellDataVector(
      axial_stress_GPs_all_row_elements,
      num_GPs_per_element_stresses_translational,
      "axial_stress_GPs" );

  runtime_vtuwriter_->AppendVisualizationCellDataVector(
      shear_stress_2_GPs_all_row_elements,
      num_GPs_per_element_stresses_translational,
      "shear_stress_2_GPs" );

  runtime_vtuwriter_->AppendVisualizationCellDataVector(
      shear_stress_3_GPs_all_row_elements,
      num_GPs_per_element_stresses_translational,
      "shear_stress_3_GPs" );


  runtime_vtuwriter_->AppendVisualizationCellDataVector(
      torque_GPs_all_row_elements,
      num_GPs_per_element_stresses_rotational,
      "torque_GPs" );

  runtime_vtuwriter_->AppendVisualizationCellDataVector(
      bending_moment_2_GPs_all_row_elements,
      num_GPs_per_element_stresses_rotational,
      "bending_moment_2_GPs" );

  runtime_vtuwriter_->AppendVisualizationCellDataVector(
      bending_moment_3_GPs_all_row_elements,
      num_GPs_per_element_stresses_rotational,
      "bending_moment_3_GPs" );
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
BeamDiscretizationRuntimeVtuWriter::AppendElementElasticEnergy()
{
  dserror("not implemented yet");

//  // count number of nodes and number for each processor; output is completely independent of
//  // the number of processors involved
//  unsigned int num_row_elements = discretization_->NumMyRowElements();
//
//  // processor owning the element
//  std::vector<double> energy_elastic;
//  energy_elastic.reserve( num_row_elements );
//
//
//  // loop over my elements and collect the data about triads/base vectors
//  for (unsigned int iele=0; iele<num_row_elements; ++iele)
//  {
//    const DRT::Element* ele = discretization_->lRowElement(iele);
//
//    // check for beam element
//    const DRT::ELEMENTS::Beam3Base* beamele = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);
//
//    // Todo for now, simply skip all other elements
//    if ( beamele == NULL )
//      continue;
//
//
//    // Todo get Eint_ from previous element evaluation call
//  for( int i = 0; i < num_cells_per_element_[iele]; ++i )
//    energy_elastic.push_back( beamele->GetElasticEnergy() );
//  }
//
//  // append the solution vector to the visualization data of the vtu writer object
//  runtime_vtuwriter_->AppendVisualizationCellDataVector(
//      energy_elastic, 1, "element_elastic_energy" );

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
BeamDiscretizationRuntimeVtuWriter::WriteFiles()
{
  runtime_vtuwriter_->WriteFiles();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
BeamDiscretizationRuntimeVtuWriter::WriteCollectionFileOfAllWrittenFiles()
{
  runtime_vtuwriter_->WriteCollectionFileOfAllWrittenFiles(
      DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix() +
      "-" + discretization_->Name() + "-beams" );
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void
BeamDiscretizationRuntimeVtuWriter::InsertVectorValuesAtBackOfOtherVector(
    const std::vector<double> & vector_input,
    std::vector<double> & vector_output )
{
  vector_output.reserve(
      vector_output.size() + vector_input.size() );

  std::copy( vector_input.begin(), vector_input.end(), std::back_inserter(vector_output) );
}

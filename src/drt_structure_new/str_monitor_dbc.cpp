/*----------------------------------------------------------------------------*/
/*!
\file str_monitor_dbc.cpp

\brief Monitor tagged Dirichlet boundary conditions

\level 3

\maintainer Michael Hiermeier
\date Jan 2, 2018

*/
/*----------------------------------------------------------------------------*/

#include "str_monitor_dbc.H"

#include "str_dbc.H"
#include "str_timint_basedataglobalstate.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/every_iteration_writer.H"

#include "../linalg/linalg_utils.H"
#include "../drt_geometry/element_volume.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MonitorDbc::Init(
    DRT::DiscretizationInterface& discret,
    STR::TIMINT::BaseDataGlobalState& gstate,
    STR::Dbc& dbc )
{
  issetup_ = false;
  isinit_ = false;

  std::vector<const DRT::Condition*> tagged_conds;
  GetTaggedCondition( tagged_conds, "Dirichlet", "monitor_reaction", discret );

  // There are no tagged conditions. This indicates that the reaction forces
  // shall not be monitored thus we can leave.
  isempty_ = ( tagged_conds.size() == 0 );
  if ( isempty_ )
  {
    isinit_ = true;
    return;
  }

  // copy the information of the tagged Dirichlet condition into a new
  // auxiliary "ReactionForce" condition and build the related geometry
  for ( const DRT::Condition* tagged_cond : tagged_conds )
    CreateReactionForceCondition( *tagged_cond, discret );

  // build geometry
  discret.FillComplete(false,false,true);

  discret_ptr_ = &discret;
  gstate_ptr_  = &gstate;
  dbc_ptr_ = &dbc;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MonitorDbc::GetTaggedCondition(
    std::vector<const DRT::Condition*>& tagged_conds,
    const std::string& cond_name,
    const std::string& tag_name,
    const DRT::DiscretizationInterface& discret ) const
{
  tagged_conds.clear();

  std::vector<std::string> cond_names;
  std::vector<Teuchos::RCP<DRT::Condition> > cond_vec;
  discret.GetCondition( cond_name, cond_vec );

  for ( auto& cond_ptr : cond_vec )
  {
    const std::string* cptr = cond_ptr->Get<std::string>( "tag" );

    if ( not cptr )
      continue;

    if ( *cptr == tag_name )
      tagged_conds.push_back( cond_ptr.get() );
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::MonitorDbc::GetUniqueId(
    int tagged_id,
    DRT::Condition::GeometryType gtype ) const
{
  switch ( gtype )
  {
    case DRT::Condition::Point:
      return tagged_id + 100;
    case DRT::Condition::Line:
      return tagged_id + 1000;
    case DRT::Condition::Surface:
      return tagged_id + 10000;
    default:
      dserror("Unsupported geometry type! (enum=%d)", gtype);
      exit( EXIT_FAILURE );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MonitorDbc::CreateReactionForceCondition(
    const DRT::Condition& tagged_cond,
    DRT::DiscretizationInterface& discret ) const
{
  const int new_id = GetUniqueId( tagged_cond.Id(), tagged_cond.GType() );

  Teuchos::RCP<DRT::Condition> rcond_ptr =
      Teuchos::rcp( new DRT::Condition( new_id, DRT::Condition::ElementTag, true,
          tagged_cond.GType() ) );

  rcond_ptr->Add( "onoff", *(tagged_cond.Get<std::vector<int> >("onoff")) );
  rcond_ptr->Add( "Node Ids", *tagged_cond.Nodes() );

  dynamic_cast<DRT::Discretization&>(discret).SetCondition("ReactionForce",rcond_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MonitorDbc::Setup()
{
  ThrowIfNotInit();

  if ( isempty_ )
  {
    issetup_ = true;
    return;
  }

  std::vector<Teuchos::RCP<DRT::Condition> > rconds;
  discret_ptr_->GetCondition("ReactionForce",rconds);
  for ( const Teuchos::RCP<DRT::Condition>& rcond_ptr : rconds )
  {
    DRT::Condition& rcond = *rcond_ptr;
    auto ipair = react_maps_.insert( std::make_pair( rcond.Id(),
        std::vector<Teuchos::RCP<Epetra_Map> >(3,Teuchos::null) ) );

    if ( not ipair.second )
      dserror("The reaction condition id #%d seems to be non-unique!",
          rcond.Id() );

    CreateReactionMaps( *discret_ptr_, rcond, ipair.first->second.data() );
  }

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MonitorDbc::CreateReactionMaps(
    const DRT::DiscretizationInterface& discret,
    const DRT::Condition& rcond,
    Teuchos::RCP<Epetra_Map>* react_maps ) const
{
  const std::vector<int>* onoff = rcond.Get<std::vector<int> >("onoff");
  const std::vector<int>* nids  = rcond.Get<std::vector<int> >("Node Ids");
  std::vector<int> my_dofs[DIM];
  int ndof = 0;
  for ( int i : *onoff )
    ndof += i;

  for ( unsigned i=0; i<DIM; ++i )
   my_dofs[i].reserve( nids->size()*ndof );

  const Epetra_Comm& comm = discret.Comm();
  for ( int nid : *nids )
  {
    const int rlid = discret.NodeRowMap()->LID( nid );
    if ( rlid == -1 )
      continue;

    const DRT::Node* node = discret.lRowNode( rlid );

    for ( unsigned i = 0; i<onoff->size(); ++i )
      if ( (*onoff)[i] == 1 )
        my_dofs[i].push_back( discret.Dof(node,i) );
  }

  for ( unsigned i=0; i<DIM; ++i )
    react_maps[i] = Teuchos::rcp( new Epetra_Map( -1, my_dofs[i].size(),
        my_dofs[i].data(), 0, comm ) );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MonitorDbc::Execute( IO::DiscretizationWriter& writer )
{
  ThrowIfNotInit();
  ThrowIfNotSetup();

  if ( isempty_ )
    return;

  std::vector<Teuchos::RCP<DRT::Condition> > rconds;
  discret_ptr_->GetCondition("ReactionForce",rconds);

  static std::vector<std::string> full_filepaths =
      CreateDirectoryAndFiles( rconds, writer, true );

  double area[2] = {0.0, 0.0};
  double& area_ref = area[0];
  double& area_curr = area[1];
  LINALG::Matrix<DIM,1> rforce_xyz(false);

  auto filepath = full_filepaths.cbegin();
  for ( const Teuchos::RCP<DRT::Condition>& rcond_ptr : rconds )
  {
    std::fill( area, area+2, 0.0 );
    std::fill( rforce_xyz.A(), rforce_xyz.A()+DIM, 0.0 );

    const int rid = rcond_ptr->Id();
    GetArea( area, rcond_ptr.get() );

    GetReactionForce( rforce_xyz, react_maps_[rid].data() );

    WriteResultsToFile( *(filepath++), rforce_xyz, area_ref, area_curr );
    WriteResultsToScreen( rforce_xyz, area_ref, area_curr );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MonitorDbc::WriteResultsToFile(
    const std::string& full_filepath,
    const LINALG::Matrix<DIM,1>& rforce,
    const double& area_ref,
    const double& area_curr ) const
{
  if ( Comm().MyPID() != 0 )
    return;

  std::ofstream of(full_filepath,std::ios_base::out|std::ios_base::app);

  WriteResults( of, OF_WIDTH, OF_PRECISION, gstate_ptr_->GetStepN(), rforce,
      area_ref, area_curr );

  of.close();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MonitorDbc::WriteResultsToScreen(
    const LINALG::Matrix<DIM,1>& rforce,
    const double& area_ref,
    const double& area_curr ) const
{
  IO::cout << "\n\n--- Monitor Dirichlet boundary conditions\n";
  WriteHeader( IO::cout.os(), OS_WIDTH );
  WriteResults( IO::cout.os(), OS_WIDTH, OS_PRECISION, gstate_ptr_->GetStepN(),
      rforce, area_ref, area_curr );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::vector<std::string> STR::MonitorDbc::CreateDirectoryAndFiles(
    const std::vector<Teuchos::RCP<DRT::Condition> >& rconds,
    IO::DiscretizationWriter& writer,
    const bool write_header ) const
{
  std::vector<std::string> full_filepaths( rconds.size() );

  if ( Comm().MyPID() != 0 )
    return full_filepaths;

  const std::string& file_path = writer.Output()->FileName();

  const std::string path( IO::ExtractPath( file_path ) );
  const std::string dir_name( IO::RemoveRestartCountFromFileName(
      IO::ExtractFileName( file_path ) ) + "_monitor_dbc" );

  const std::string full_dirpath( path + dir_name );
  IO::CreateDirectory( full_dirpath, Comm().MyPID() );

  size_t i = 0;
  for ( const Teuchos::RCP<DRT::Condition>& rcond : rconds )
  {
    std::string& full_filepath = full_filepaths[i++];
    full_filepath = full_dirpath + "/" + std::to_string(rcond->Id()) +
        "_monitor_dbc.data";

    if ( write_header )
    {
      // clear old content
      std::ofstream of( full_filepath, std::ios_base::out );
      WriteHeader( of, OF_WIDTH, rcond.get() );
      of.close();
    }
  }

  return full_filepaths;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MonitorDbc::WriteHeader(
    std::ostream& os,
    const int col_width,
    const DRT::Condition* cond ) const
{
  if ( cond )
  {
    cond->Print( os );
    os << "\n\n";
  }

  os << std::setw(col_width) << "step"
     << std::setw(col_width) << "ref_area"
     << std::setw(col_width) << "curr_area"
     << std::setw(col_width) << "f_x"
     << std::setw(col_width) << "f_y"
     << std::setw(col_width) << "f_z\n";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MonitorDbc::WriteResults(
    std::ostream& os,
    const int col_width,
    const int precision,
    const unsigned step,
    const LINALG::Matrix<DIM,1>& rforce,
    const double & area_ref,
    const double & area_curr ) const
{
  os << std::setw(col_width) << step << std::setprecision(precision);
  os << std::setw(col_width) << std::scientific << area_ref
     << std::setw(col_width) << std::scientific << area_curr;

  for ( unsigned i=0; i<DIM; ++i )
    os << std::setw(col_width) << std::scientific << rforce(i,0);

  os << "\n";
  os << std::flush;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Comm& STR::MonitorDbc::Comm() const
{
  return discret_ptr_->Comm();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MonitorDbc::GetArea( double area[],
    const DRT::Condition* rcond ) const
{
  // no area for point DBCs
  if ( rcond->GType() == DRT::Condition::Point )
  {
    std::fill( area, area+2, 0.0 );
    return;
  }

  const DRT::Discretization& discret =
      dynamic_cast<const DRT::Discretization&>( *discret_ptr_ );

  enum AreaType : int { ref=0, curr=1 };
  double larea[2] = { 0.0, 0.0 };
  LINALG::SerialDenseMatrix xyze_ref;
  LINALG::SerialDenseMatrix xyze_curr;

  const std::map<int,Teuchos::RCP<DRT::Element> >& celes = rcond->Geometry();
  Teuchos::RCP<const Epetra_Vector> dispn = gstate_ptr_->GetDisNp();
  Epetra_Vector dispn_col( *discret.DofColMap(), true );
  LINALG::Export( *dispn, dispn_col );

  for ( auto& cele_pair : celes )
  {
    const DRT::Element* cele = cele_pair.second.get();
    const DRT::FaceElement* fele = dynamic_cast<const DRT::FaceElement*>( cele );
    if ( !fele )
      dserror("No face element!");

    if ( !fele->ParentElement() or
         fele->ParentElement()->Owner() != discret.Comm().MyPID() )
      continue;

    const DRT::Node* const * fnodes = fele->Nodes();
    const unsigned num_fnodes = fele->NumNode();
    std::vector<int> fele_dofs;
    fele_dofs.reserve( num_fnodes*DIM );

    for ( unsigned i=0; i<num_fnodes; ++i )
      discret.Dof( fele, fnodes[i], fele_dofs );

    std::vector<double> mydispn;
    DRT::UTILS::ExtractMyValues( dispn_col, mydispn, fele_dofs );

    xyze_ref.Reshape( DIM, num_fnodes );
    xyze_curr.Reshape( DIM, num_fnodes );

    for ( unsigned i=0; i<num_fnodes; ++i )
    {
      const DRT::Node& fnode = *fnodes[i];
      std::copy( fnode.X(), fnode.X()+DIM, &xyze_ref(0,i) );
      std::copy( fnode.X(), fnode.X()+DIM, &xyze_curr(0,i) );

      std::vector<int> ndofs;
      discret.Dof( &fnode, ndofs );

      for ( unsigned d=0; d<ndofs.size(); ++d )
      {
        const int ndof = ndofs[d];

        size_t fedof_count = 0;
        for ( auto cit = fele_dofs.cbegin(); cit!=fele_dofs.cend(); ++cit, ++fedof_count )
        {
          if ( *cit == ndof )
            break;
        }

        if ( fedof_count == fele_dofs.size() )
          dserror("Couln't find the face element dof corresponding to the "
              "current node!" );

        xyze_curr(d,i) += mydispn[fedof_count];
      }
    }

    larea[AreaType::ref] += GEO::ElementArea( fele->Shape(), xyze_ref );
    larea[AreaType::curr] += GEO::ElementArea( fele->Shape(), xyze_curr );
  }

  discret.Comm().SumAll( &larea[0], &area[0], 2 );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::MonitorDbc::GetReactionForce(
    LINALG::Matrix<DIM,1>& rforce_xyz,
    const Teuchos::RCP<Epetra_Map>* react_maps ) const
{
  Epetra_Vector complete_freact( *gstate_ptr_->GetFreactNp() );
  dbc_ptr_->RotateGlobalToLocal( Teuchos::rcpFromRef( complete_freact ) );

  LINALG::Matrix<DIM,1> lrforce_xyz(true);
  for ( unsigned d=0; d<DIM; ++d )
  {
    Teuchos::RCP<Epetra_Vector> partial_freact_ptr =
        LINALG::ExtractMyVector( complete_freact, *(react_maps[d]) );

    double& lrforce_comp = lrforce_xyz(d,0);
    const double* vals = partial_freact_ptr->Values();
    for ( int i=0; i<react_maps[d]->NumMyElements(); ++i )
      lrforce_comp += vals[i];
  }

  discret_ptr_->Comm().SumAll( lrforce_xyz.A(), rforce_xyz.A(), DIM );
  return rforce_xyz.Norm2();
}

/*-----------------------------------------------------------*/
/*! \file

\brief A class handling a (periodic) bounding box as simulation volume

\maintainer Jonas Eichinger

\level 3

*/
/*-----------------------------------------------------------*/

#include "../drt_beaminteraction/periodic_boundingbox.H"
#include "../drt_binstrategy/binning_strategy_utils.H"

#include "../drt_io/discretization_runtime_vtu_writer.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dofset_independent.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_densematrix_communication.H"
#include "../linalg/linalg_utils_densematrix_manipulation.H"

#include "../drt_geometry/intersection_math.H"
#include "../drt_inpar/inpar_binningstrategy.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::MESHFREE::BoundingBox::BoundingBox()
    : isinit_(false),
      issetup_(false),
      boxdiscret_(Teuchos::null),
      disn_row_(Teuchos::null),
      disn_col_(Teuchos::null),
      empty_(true),
      haveperiodicbc_(false),
      havedirichletbc_(false),
      box_(true),
      vtu_writer_ptr_(Teuchos::null)
{
  // initialize arrays
  for (int idim = 0; idim < 3; ++idim)
  {
    pbconoff_[idim] = false;
    edgelength_[idim] = 0.0;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::Init()
{
  issetup_ = false;

  // get bounding box specified in the input file
  // fixme: like this or by eight nodes of element in discret
  box_.PutScalar(1.0e12);
  std::istringstream xaabbstream(Teuchos::getNumericStringParameter(
      DRT::Problem::Instance()->BinningStrategyParams(), "DOMAINBOUNDINGBOX"));
  for (int col = 0; col < 2; ++col)
  {
    for (int row = 0; row < 3; ++row)
    {
      double value = 1.0e12;
      if (xaabbstream >> value)
        box_(row, col) = value;
      else
        dserror(
            " Specify six values for bounding box in three dimensional problem."
            " Fix your input file.");
    }
  }

  // set up boundary conditions
  std::istringstream periodicbc(Teuchos::getNumericStringParameter(
      DRT::Problem::Instance()->BinningStrategyParams(), "PERIODICONOFF"));

  // loop over all spatial directions
  for (int dim = 0; dim < 3; ++dim)
  {
    int val = -1;
    if (periodicbc >> val)
    {
      if (val)
      {
        // set flag
        pbconoff_[dim] = true;

        // set global flag
        haveperiodicbc_ = true;
      }

      // edge length of box
      edgelength_[dim] = box_(dim, 1) - box_(dim, 0);
    }
    else
    {
      dserror(
          "Enter three values to specify each direction as periodic or non periodic."
          "Fix your input file ...");
    }
  }

  // todo
  havedirichletbc_ = false;
  empty_ = false;
  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::Init(
    LINALG::Matrix<3, 2> const& box, std::vector<bool> const& pbconoff)
{
  issetup_ = false;

  // get bounding box specified in the input file
  box_ = box;

  // loop over all spatial directions
  for (unsigned int dim = 0; dim < 3; ++dim)
  {
    pbconoff_[dim] = pbconoff[dim];
    if (pbconoff_[dim]) haveperiodicbc_ = true;

    // edge length of box
    edgelength_[dim] = box_(dim, 1) - box_(dim, 0);
  }

  havedirichletbc_ = false;
  empty_ = false;
  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::Setup()
{
  ThrowIfNotInit();

  // initialize bounding box discretization
  SetupBoundingBoxDiscretization();

  if (boxdiscret_->GetCondition("Dirichlet") != NULL) havedirichletbc_ = true;

  // displacement vector in row and col format
  disn_row_ = LINALG::CreateVector(*boxdiscret_->DofRowMap(), true);
  disn_col_ = LINALG::CreateVector(*boxdiscret_->DofColMap(), true);

  // initialize bounding box runtime output
  if (DRT::Problem::Instance()
          ->IOParams()
          .sublist("RUNTIME VTK OUTPUT")
          .get<int>("INTERVAL_STEPS") != -1)
    InitRuntimeOutput();

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 * (public)                                                                   |
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::SetupBoundingBoxDiscretization()
{
  if (DRT::Problem::Instance()->DoesExistDis("boundingbox"))
  {
    boxdiscret_ = DRT::Problem::Instance()->GetDis("boundingbox");

    if (boxdiscret_->Filled() == false) boxdiscret_->FillComplete(true, false, false);

    // create fully overlapping boundingbox discret
    Teuchos::RCP<Epetra_Map> rednodecolmap = LINALG::AllreduceEMap(*boxdiscret_->NodeRowMap());
    Teuchos::RCP<Epetra_Map> redelecolmap = LINALG::AllreduceEMap(*boxdiscret_->ElementRowMap());

    // do the fully overlapping ghosting of the bounding box element to have everything redundant
    boxdiscret_->ExportColumnNodes(*rednodecolmap);
    boxdiscret_->ExportColumnElements(*redelecolmap);

    boxdiscret_->FillComplete(true, false, false);
  }

  if (not DRT::Problem::Instance()->DoesExistDis("boundingbox") or
      boxdiscret_->NumMyColElements() == 0)
  {
    if (not DRT::Problem::Instance()->DoesExistDis("boundingbox"))
    {
      Teuchos::RCP<Epetra_Comm> com =
          Teuchos::rcp(DRT::Problem::Instance()->GetDis("structure")->Comm().Clone());
      boxdiscret_ = Teuchos::rcp(new DRT::Discretization("boundingbox", com));
    }
    else
    {
      boxdiscret_ = DRT::Problem::Instance()->GetDis("boundingbox");
    }

    // create nodes
    double cornerpos[3];
    int node_ids[8];
    for (int corner_i = 0; corner_i < 8; ++corner_i)
    {
      UndeformedBoxCornerPointPosition(corner_i, cornerpos);
      node_ids[corner_i] = corner_i;

      Teuchos::RCP<DRT::Node> newnode = Teuchos::rcp(new DRT::Node(corner_i, cornerpos, 0));
      boxdiscret_->AddNode(newnode);
    }

    // assign nodes to element
    Teuchos::RCP<DRT::Element> newele = DRT::UTILS::Factory("VELE3", "Polynomial", 0, 0);
    newele->SetNodeIds(8, node_ids);
    boxdiscret_->AddElement(newele);
  }

  // build independent dof set
  Teuchos::RCP<DRT::IndependentDofSet> independentdofset =
      Teuchos::rcp(new DRT::IndependentDofSet(true));
  boxdiscret_->ReplaceDofSet(independentdofset);
  boxdiscret_->FillComplete();
}

/*----------------------------------------------------------------------------*
 * (public)                                                                   |
 *----------------------------------------------------------------------------*/
bool GEO::MESHFREE::BoundingBox::Shift3D(
    LINALG::Matrix<3, 1>& d, LINALG::Matrix<3, 1> const X) const
{
  ThrowIfNotInit();

  bool shifted = false;

  if (not haveperiodicbc_) return shifted;

  // x = X + d
  LINALG::Matrix<3, 1> x(X);
  x.Update(1.0, d, 1.0);

  LINALG::Matrix<3, 1> x_ud(true);
  TransformFromGlobalToUndeformedBoundingBoxSystem(x, x_ud);

  // shift
  for (int dim = 0; dim < 3; ++dim)
    if (Shift1D(dim, x_ud(dim))) shifted = true;

  x.Clear();
  d.Clear();
  TransformFromUndeformedBoundingBoxSystemToGlobal(x_ud, x);

  // d = x - X
  d.Update(1.0, x, -1.0, X);

  return shifted;
}

/*----------------------------------------------------------------------------*
 * (public)                                                                   |
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::GetXiOfIntersection3D(
    LINALG::Matrix<3, 1> const& x1, LINALG::Matrix<3, 1> const& x2, LINALG::Matrix<3, 1>& xi) const
{
  ThrowIfNotInit();
  GetXiOfIntersection3D(x1, x2, xi, box_);
}
/*----------------------------------------------------------------------------*
 * (public)                                                                   |
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::GetXiOfIntersection3D(LINALG::Matrix<3, 1> const& x1,
    LINALG::Matrix<3, 1> const& x2, LINALG::Matrix<3, 1>& xi, LINALG::Matrix<3, 2> const& box) const
{
  ThrowIfNotInit();

  // set default values
  for (unsigned int dim = 0; dim < 3; ++dim) xi(dim) = 2.0;

  LINALG::Matrix<3, 1> x1_ud(true), x2_ud(true);
  TransformFromGlobalToUndeformedBoundingBoxSystem(x1, x1_ud);
  TransformFromGlobalToUndeformedBoundingBoxSystem(x2, x2_ud);

  // check which points resides within in which direction
  std::vector<bool> x1_within_in_dir(3, true), x2_within_in_dir(3, true);
  Within(box, x1_ud, x1_within_in_dir);
  Within(box, x2_ud, x2_within_in_dir);

  for (unsigned int dim = 0; dim < 3; ++dim)
  {
    bool x1_in = x1_within_in_dir[dim];
    bool x2_in = x2_within_in_dir[dim];

    // no intersection in case both are in or both are out
    if ((x1_in or not x2_in) and (x2_in or not x1_in)) continue;

    // x1 out
    if (x2_in)
    {
      if (x1_ud(dim) < box_min(box, dim))
      {
        xi(dim) = (2.0 * box_min(box, dim) - (x2_ud(dim) + x1_ud(dim))) / (x2_ud(dim) - x1_ud(dim));
      }
      else if (x1_ud(dim) > box_max(box, dim))
      {
        xi(dim) = (2.0 * box_max(box, dim) - (x2_ud(dim) + x1_ud(dim))) / (x1_ud(dim) - x2_ud(dim));
      }
      else
      {
        dserror("Your true/false logic is wrong.");
      }
    }
    else if (x1_in)
    {
      if (x2_ud(dim) < box_min(box, dim))
      {
        xi(dim) = (2.0 * box_min(box, dim) - (x2_ud(dim) + x1_ud(dim))) / (x1_ud(dim) - x2_ud(dim));
      }
      else if (x2_ud(dim) > box_max(box, dim))
      {
        xi(dim) = (2.0 * box_max(box, dim) - (x2_ud(dim) + x1_ud(dim))) / (x2_ud(dim) - x1_ud(dim));
      }
      else
      {
        dserror("Your true/false logic is wrong.");
      }
    }
    else
    {
      dserror("Your true/false logic is wrong.");
    }
  }
}

/*----------------------------------------------------------------------------*
 * (public)                                                                   |
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::UnShift3D(
    LINALG::Matrix<3, 1>& d, LINALG::Matrix<3, 1> const& ref, LINALG::Matrix<3, 1> const X) const
{
  ThrowIfNotInit();

  if (not haveperiodicbc_) return;

  // x = X + d
  LINALG::Matrix<3, 1> x(X);
  x.Update(1.0, d, 1.0);

  LINALG::Matrix<3, 1> x_ud(true), ref_ud(true);
  TransformFromGlobalToUndeformedBoundingBoxSystem(x, x_ud);
  TransformFromGlobalToUndeformedBoundingBoxSystem(ref, ref_ud);

  for (int dim = 0; dim < 3; ++dim) UnShift1D(dim, x_ud(dim), ref_ud(dim));

  x.Clear();
  d.Clear();
  TransformFromUndeformedBoundingBoxSystemToGlobal(x_ud, x);

  // d = x - X
  d.Update(1.0, x, -1.0, X);
}

/*----------------------------------------------------------------------------*
 * (public)                                                                   |
 *----------------------------------------------------------------------------*/
bool GEO::MESHFREE::BoundingBox::CheckIfShiftBetweenPoints(LINALG::Matrix<3, 1>& d,
    LINALG::Matrix<3, 1> const& ref, std::vector<bool>& shift_in_dim,
    LINALG::Matrix<3, 1> const X) const
{
  ThrowIfNotInit();

  shift_in_dim.clear();
  shift_in_dim.resize(3);

  if (not haveperiodicbc_)
  {
    shift_in_dim = std::vector<bool>(3, false);
    return false;
  }

  // x = X + d
  LINALG::Matrix<3, 1> x(X);
  x.Update(1.0, d, 1.0);

  LINALG::Matrix<3, 1> x_ud(true), ref_ud(true);
  TransformFromGlobalToUndeformedBoundingBoxSystem(x, x_ud);
  TransformFromGlobalToUndeformedBoundingBoxSystem(ref, ref_ud);

  for (int dim = 0; dim < 3; ++dim) shift_in_dim[dim] = UnShift1D(dim, x_ud(dim), ref_ud(dim));

  x.Clear();
  d.Clear();
  TransformFromUndeformedBoundingBoxSystemToGlobal(x_ud, x);

  // d = x - X
  d.Update(1.0, x, -1.0, X);

  return (shift_in_dim[0] or shift_in_dim[1] or shift_in_dim[2]);
}

/*----------------------------------------------------------------------------*
 * (private)                                                                   |
 *----------------------------------------------------------------------------*/
bool GEO::MESHFREE::BoundingBox::Shift1D(int dim, double& d, double const& X) const
{
  ThrowIfNotInit();

  bool shifted = false;

  if (not pbconoff_[dim]) return shifted;

  double x = d + X;

  if (x < box_min(dim))
  {
    shifted = true;
    d += edgelength_[dim] * std::ceil(std::abs((x - box_(dim, 0)) / edgelength_[dim]));
  }
  else if (x > box_max(dim))
  {
    shifted = true;
    d -= edgelength_[dim] * std::ceil(std::abs((x - box_(dim, 1)) / edgelength_[dim]));
  }

  return shifted;
}

/*----------------------------------------------------------------------------*
 * (private)                                                                   |
 *----------------------------------------------------------------------------*/
bool GEO::MESHFREE::BoundingBox::UnShift1D(
    int dim, double& d, double const& ref, double const& X) const
{
  ThrowIfNotInit();

  bool unshifted = false;

  if (not pbconoff_[dim]) return unshifted;

  double x = d + X;

  if (x - ref < -0.5 * edgelength_[dim])
  {
    unshifted = true;
    d += edgelength_[dim];
  }
  else if (x - ref > 0.5 * edgelength_[dim])
  {
    unshifted = true;
    d -= edgelength_[dim];
  }

  return unshifted;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::MESHFREE::BoundingBox::InBetween(double smin, double smax, double omin, double omax) const
{
  double tol = GEO::TOL7;
  return ((omax > smin - tol) and (smax > omin - tol));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::RandomPosWithin(LINALG::Matrix<3, 1>& randpos) const
{
  ThrowIfNotInit();

  DRT::Problem::Instance()->Random()->SetRandRange(0.0, 1.0);
  std::vector<double> randuni;
  DRT::Problem::Instance()->Random()->Uni(randuni, 3);

  LINALG::Matrix<3, 1> randpos_ud(true);
  for (int dim = 0; dim < 3; ++dim)
    randpos_ud(dim) = box_min(dim) + (edgelength_[dim] * randuni[dim]);

  TransformFromUndeformedBoundingBoxSystemToGlobal(randpos_ud, randpos);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::AddPoint(const double* x)
{
  dserror("Check before use.");
  if (empty_)
  {
    empty_ = false;
    box_(0, 0) = box_(0, 1) = x[0];
    box_(1, 0) = box_(1, 1) = x[1];
    box_(2, 0) = box_(2, 1) = x[2];
  }
  else
  {
    box_(0, 0) = std::min(box_(0, 0), x[0]);
    box_(1, 0) = std::min(box_(1, 0), x[1]);
    box_(2, 0) = std::min(box_(2, 0), x[2]);
    box_(0, 1) = std::max(box_(0, 1), x[0]);
    box_(1, 1) = std::max(box_(1, 1), x[1]);
    box_(2, 1) = std::max(box_(2, 1), x[2]);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::MESHFREE::BoundingBox::Within(BoundingBox const& b) const
{
  dserror("Check before use.");
  if (empty_) return true;
  return (InBetween(box_min(0), box_max(0), b.box_min(0), b.box_max(0)) and
          InBetween(box_min(1), box_max(1), b.box_min(1), b.box_max(1)) and
          InBetween(box_min(2), box_max(2), b.box_min(2), b.box_max(2)));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::MESHFREE::BoundingBox::Within(const double* x, std::vector<bool>& within_in_dir) const
{
  ThrowIfNotInit();

  within_in_dir.resize(3);
  within_in_dir[0] = InBetween(box_min(0), box_max(0), x[0], x[0]);
  within_in_dir[1] = InBetween(box_min(1), box_max(1), x[1], x[1]);
  within_in_dir[2] = InBetween(box_min(2), box_max(2), x[2], x[2]);

  return (within_in_dir[0] and within_in_dir[1] and within_in_dir[2]);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::MESHFREE::BoundingBox::Within(
    LINALG::Matrix<3, 1> const& x, std::vector<bool>& within_in_dir) const
{
  ThrowIfNotInit();

  within_in_dir.resize(3);
  within_in_dir[0] = InBetween(box_min(0), box_max(0), x(0), x(0));
  within_in_dir[1] = InBetween(box_min(1), box_max(1), x(1), x(1));
  within_in_dir[2] = InBetween(box_min(2), box_max(2), x(2), x(2));

  return (within_in_dir[0] and within_in_dir[1] and within_in_dir[2]);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::MESHFREE::BoundingBox::Within(LINALG::Matrix<3, 2> const& box,
    LINALG::Matrix<3, 1> const& x, std::vector<bool>& within_in_dir) const
{
  ThrowIfNotInit();

  within_in_dir.resize(3);
  within_in_dir[0] = InBetween(box_min(box, 0), box_max(box, 0), x(0), x(0));
  within_in_dir[1] = InBetween(box_min(box, 1), box_max(box, 1), x(1), x(1));
  within_in_dir[2] = InBetween(box_min(box, 2), box_max(box, 2), x(2), x(2));

  return (within_in_dir[0] and within_in_dir[1] and within_in_dir[2]);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool GEO::MESHFREE::BoundingBox::Within(const Epetra_SerialDenseMatrix& xyz) const
{
  dserror("Check before use.");
  BoundingBox bb;
  int numnode = xyz.N();
  for (int i = 0; i < numnode; ++i)
  {
    bb.AddPoint(&xyz(0, i));
  }
  return Within(bb);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::Print()
{
  if (empty_)
  {
    std::cout << "  BB: {}\n";
  }
  else
  {
    std::cout << "  BB: {(" << box_(0, 0) << "," << box_(1, 0) << "," << box_(2, 0) << ")-("
              << box_(0, 1) << "," << box_(1, 1) << "," << box_(2, 1) << ")}\n";
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::ApplyDirichlet(double timen)
{
  ThrowIfNotInitOrSetup();

  Teuchos::ParameterList p;
  p.set("total time", timen);

  // disn_ then also holds prescribed new Dirichlet displacements
  boxdiscret_->ClearState();
  boxdiscret_->EvaluateDirichlet(
      p, disn_row_, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  boxdiscret_->ClearState();

  // export to col format
  LINALG::Export(*disn_row_, *disn_col_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::InitRuntimeOutput()
{
  vtu_writer_ptr_ = Teuchos::rcp(new DiscretizationRuntimeVtuWriter());

  // Todo: we need a better upper bound for total number of time steps here
  // however, this 'only' affects the number of leading zeros in the vtk file names
  const unsigned int num_timesteps_in_simulation_upper_bound = 1000000;

  // initialize the writer object
  double time = -1.0;
  if (DRT::Problem::Instance()->Restart())
  {
    IO::DiscretizationReader ioreader(
        DRT::Problem::Instance()->GetDis("structure"), DRT::Problem::Instance()->Restart());
    time = ioreader.ReadDouble("time");
  }
  else
  {
    time = 0.0;
  }

  vtu_writer_ptr_->Initialize(boxdiscret_, num_timesteps_in_simulation_upper_bound, time, true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::RuntimeOutputStepState(double timen, int stepn) const
{
  ThrowIfNotInitOrSetup();

  if (vtu_writer_ptr_ == Teuchos::null) return;

  // reset time and time step of the writer object
  vtu_writer_ptr_->ResetTimeAndTimeStep(timen, stepn);
  vtu_writer_ptr_->AppendDofBasedResultDataVector(disn_col_, 3, 0, "displacement");

  // finalize everything and write all required VTU files to filesystem
  vtu_writer_ptr_->WriteFiles();
  vtu_writer_ptr_->WriteCollectionFileOfAllWrittenFiles();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LINALG::Matrix<3, 1> GEO::MESHFREE::BoundingBox::ReferencePosOfCornerPoint(int i) const
{
  // dof gids of node i (note: each proc just has one element and eight nodes,
  // therefore local numbering from 0 to 7 on each proc)
  DRT::Node* node_i = boxdiscret_->lColNode(i);

  LINALG::Matrix<3, 1> x(true);
  for (int dim = 0; dim < 3; ++dim) x(dim) = node_i->X()[dim];

  return x;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LINALG::Matrix<3, 1> GEO::MESHFREE::BoundingBox::CurrentPositionOfCornerPoint(int i) const
{
  // dof gids of node i (note: each proc just has one element and eight nodes,
  // therefore local numbering from 0 to 7 on each proc)
  LINALG::Matrix<3, 1> x(true);
  if (boxdiscret_ != Teuchos::null)
  {
    DRT::Node* node_i = boxdiscret_->lColNode(i);
    std::vector<int> dofnode = boxdiscret_->Dof(node_i);

    for (int dim = 0; dim < 3; ++dim)
      x(dim) = node_i->X()[dim] + (*disn_col_)[disn_col_->Map().LID(dofnode[dim])];
  }
  else
  {
    x = UndeformedBoxCornerPointPosition(i);
  }

  return x;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::UndeformedBoxCornerPointPosition(int i, double* x) const
{
  // to get numbering according to baci convention of hex eles ( p.122 global report)
  if (i == 2 or i == 6)
    ++i;
  else if (i == 3 or i == 7)
    --i;

  x[0] = ((i & 1) == 1) ? box_max(0) : box_min(0);
  x[1] = ((i & 2) == 2) ? box_max(1) : box_min(1);
  x[2] = ((i & 4) == 4) ? box_max(2) : box_min(2);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LINALG::Matrix<3, 1> GEO::MESHFREE::BoundingBox::UndeformedBoxCornerPointPosition(int i) const
{
  // to get numbering according to baci convention of hex eles ( p.122 global report)
  if (i == 2 or i == 6)
    ++i;
  else if (i == 3 or i == 7)
    --i;

  LINALG::Matrix<3, 1> x(true);
  x(0) = ((i & 1) == 1) ? box_max(0) : box_min(0);
  x(1) = ((i & 2) == 2) ? box_max(1) : box_min(1);
  x(2) = ((i & 4) == 4) ? box_max(2) : box_min(2);

  return x;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::TransformFromUndeformedBoundingBoxSystemToGlobal(
    LINALG::Matrix<3, 1> const& xi, LINALG::Matrix<3, 1>& x) const
{
  ThrowIfNotInit();

  // nothing to do in case periodic bounding is not deforming
  if (not havedirichletbc_)
  {
    x = xi;
    return;
  }

  // reset globcoord variable
  x.Clear();

  // Evaluate lagrangian shape functions at xi
  LINALG::Matrix<8, 1> funct;
  LagrangePolynomialToMapFromUndeformedBoundingBoxSystemToGlobal(funct, xi(0), xi(1), xi(2));

  LINALG::Matrix<3, 8> coord;
  for (unsigned int i = 0; i < 8; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      // use shape function values for interpolation
      coord(j, i) = CurrentPositionOfCornerPoint(i)(j);
      x(j) += funct(i) * coord(j, i);
    }
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::TransformFromUndeformedBoundingBoxSystemToGlobal(
    double const* xi, double* x) const
{
  ThrowIfNotInitOrSetup();

  DRT::Node** mynodes = boxdiscret_->lColElement(0)->Nodes();
  if (!mynodes) dserror("ERROR: LocalToGlobal: Null pointer!");

  // reset globcoord variable
  for (unsigned int dim = 0; dim < 3; ++dim) x[dim] = 0.0;

  // Evaluate lagrangian shape functions at xi
  LINALG::Matrix<8, 1> funct;
  LagrangePolynomialToMapFromUndeformedBoundingBoxSystemToGlobal(funct, xi[0], xi[1], xi[2]);

  LINALG::Matrix<3, 8> coord;
  for (unsigned int i = 0; i < 8; ++i)
  {
    for (unsigned int j = 0; j < 3; ++j)
    {
      // use shape function values for interpolation
      coord(j, i) = CurrentPositionOfCornerPoint(i)(j);
      x[j] += funct(i) * coord(j, i);
    }
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
bool GEO::MESHFREE::BoundingBox::TransformFromGlobalToUndeformedBoundingBoxSystem(
    LINALG::Matrix<3, 1> const& x, LINALG::Matrix<3, 1>& xi) const
{
  ThrowIfNotInit();

  // nothing to do in case periodic bounding is not deforming
  if (not havedirichletbc_)
  {
    xi = x;
    return true;
  }

  // initialize variables
  int const numnode = 8;
  int const ndim = 3;
  double tol = GEO::TOL12;
  bool converged = false;
  LINALG::Matrix<numnode, 1> funct;
  LINALG::Matrix<ndim, numnode> deriv;
  LINALG::Matrix<ndim, numnode> pbbcurrnodepos;
  LINALG::Matrix<ndim, ndim> xjm;
  LINALG::Matrix<ndim, 1> rhs;

  // spatial configuration of this element!
  for (int k = 0; k < numnode; ++k)
    for (int j = 0; j < ndim; ++j) pbbcurrnodepos(j, k) = CurrentPositionOfCornerPoint(k)(j);

  // first estimation for xi
  for (int dim = 0; dim < ndim; ++dim) xi(dim) = x(dim);

  // solve linearized system R(xi) = x(xi) - x = N(xi) x^ - x != 0
  int numiter = 10;
  int j = 0;
  while (not converged and j < numiter)
  {
    // reset matrices
    xjm.Clear();
    deriv.Clear();
    rhs.Clear();

    // evaluate lagrangian shape functions at xi
    LagrangePolynomialToMapFromUndeformedBoundingBoxSystemToGlobal(funct, xi(0), xi(1), xi(2));
    LagrangePolynomialToMapFromUndeformedBoundingBoxSystemToGlobalDeriv1(
        deriv, xi(0), xi(1), xi(2));

    // compute dN/dxi x^
    for (int k = 0; k < numnode; ++k)
      for (int p = 0; p < ndim; ++p)
        for (int l = 0; l < ndim; ++l) xjm(p, l) += deriv(l, k) * pbbcurrnodepos(p, k);

    // rhs of (linearized equation) I
    for (int p = 0; p < ndim; ++p) rhs(p) = -x(p);

    // rhs of (linearized equation) I
    for (int k = 0; k < numnode; ++k)
      for (int p = 0; p < ndim; ++p) rhs(p) += funct(k) * pbbcurrnodepos(p, k);

    // compute norm
    double norm = 0.0;
    for (int p = 0; p < ndim; ++p) norm += rhs(p) * rhs(p);

    // check if we are converged
    if (sqrt(norm) < tol)
    {
      converged = true;
      break;
    }

    // safety check
    if (abs(xjm.Determinant()) < 1e-15)
    {
      std::cout << "WARNING !!! jacobi determinant singular! In DeformedToUndeformed(...)"
                << std::endl;
      std::cout << "JAC= " << xjm.Determinant() << std::endl;
      std::cout << "CONVERGED= " << converged << std::endl;
      dserror("*** WARNING: jacobi singular ***");
      converged = false;
      break;
    }

    // solve equation
    double xjm_invert = xjm.Invert();
    if (abs(xjm_invert) < 1e-15) dserror("ERROR: Singular Jacobian");

    // compute increment
    LINALG::Matrix<ndim, 1> deltaxi(true);
    for (int z = 0; z < ndim; ++z)
      for (int p = 0; p < ndim; ++p) deltaxi(z) -= xjm(z, p) * rhs(p);

    // incremental update
    for (int p = 0; p < ndim; ++p) xi(p) += deltaxi(p);

    ++j;
  }
  return converged;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
bool GEO::MESHFREE::BoundingBox::TransformFromGlobalToUndeformedBoundingBoxSystem(
    double const* x, double* xi) const
{
  ThrowIfNotInit();
  static LINALG::Matrix<3, 1> x_m(true), xi_m(true);
  for (int dim = 0; dim < 3; ++dim) x_m(dim) = x[dim];

  bool converged = TransformFromGlobalToUndeformedBoundingBoxSystem(x_m, xi_m);

  for (int dim = 0; dim < 3; ++dim) xi[dim] = xi_m(dim);

  x_m.Clear();
  xi_m.Clear();
  return converged;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::LagrangePolynomialToMapFromUndeformedBoundingBoxSystemToGlobal(
    LINALG::Matrix<8, 1>& funct,  ///< to be filled with shape function values
    double r, double s, double t) const
{
  ThrowIfNotInit();
  // safety check
#ifdef DEBUG
  for (int dim = 0; dim < 3; ++dim)
    if (abs(edgelength_[dim]) < 1e-12)
      dserror(" you are about to devide by zero, edgelength not correctly initialized.");
#endif

  const double Q = 1.0 / (edgelength_[0] * edgelength_[1] * edgelength_[2]);


  const double rp = r - box_min(0);
  const double rm = box_max(0) - r;
  const double sp = s - box_min(1);
  const double sm = box_max(1) - s;
  const double tp = t - box_min(2);
  const double tm = box_max(2) - t;

  funct(0) = Q * rm * sm * tm;
  funct(1) = Q * rp * sm * tm;
  funct(2) = Q * rp * sp * tm;
  funct(3) = Q * rm * sp * tm;
  funct(4) = Q * rm * sm * tp;
  funct(5) = Q * rp * sm * tp;
  funct(6) = Q * rp * sp * tp;
  funct(7) = Q * rm * sp * tp;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void GEO::MESHFREE::BoundingBox::
    LagrangePolynomialToMapFromUndeformedBoundingBoxSystemToGlobalDeriv1(
        LINALG::Matrix<3, 8>& deriv1,  ///< to be filled with shape function derivative values
        double r, double s, double t) const
{
  ThrowIfNotInit();
  // safety check
#ifdef DEBUG
  for (int dim = 0; dim < 3; ++dim)
    if (abs(edgelength_[dim]) < 1e-12)
      dserror(" you are about to devide by zero, edgelength not correctly initialized.");
#endif

  const double Q = 1.0 / (edgelength_[0] * edgelength_[1] * edgelength_[2]);

  const double rp = r - box_min(0);
  const double rm = box_max(0) - r;
  const double sp = s - box_min(1);
  const double sm = box_max(1) - s;
  const double tp = t - box_min(2);
  const double tm = box_max(2) - t;

  deriv1(0, 0) = -Q * sm * tm;
  deriv1(1, 0) = -Q * tm * rm;
  deriv1(2, 0) = -Q * rm * sm;

  deriv1(0, 1) = Q * sm * tm;
  deriv1(1, 1) = -Q * tm * rp;
  deriv1(2, 1) = -Q * rp * sm;

  deriv1(0, 2) = Q * sp * tm;
  deriv1(1, 2) = Q * tm * rp;
  deriv1(2, 2) = -Q * rp * sp;

  deriv1(0, 3) = -Q * sp * tm;
  deriv1(1, 3) = Q * tm * rm;
  deriv1(2, 3) = -Q * rm * sp;

  deriv1(0, 4) = -Q * sm * tp;
  deriv1(1, 4) = -Q * tp * rm;
  deriv1(2, 4) = Q * rm * sm;

  deriv1(0, 5) = Q * sm * tp;
  deriv1(1, 5) = -Q * tp * rp;
  deriv1(2, 5) = Q * rp * sm;

  deriv1(0, 6) = Q * sp * tp;
  deriv1(1, 6) = Q * tp * rp;
  deriv1(2, 6) = Q * rp * sp;

  deriv1(0, 7) = -Q * sp * tp;
  deriv1(1, 7) = Q * tp * rm;
  deriv1(2, 7) = Q * rm * sp;
}

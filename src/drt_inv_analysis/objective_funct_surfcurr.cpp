#ifdef HAVE_Kokkos
/*----------------------------------------------------------------------*/
/*!
\brief Surface current based objective function

\level 3

\maintainer Sebastian Brandstaeter
*/
/*----------------------------------------------------------------------*/

#include "objective_funct_surfcurr.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_inputreader.H"
#include "../drt_lib/drt_nodereader.H"
#include "../drt_mat/matpar_bundle.H"
#include "../linalg/linalg_utils.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"
#include "../drt_io/io_pstream.H"
#include "../drt_comm/comm_utils.H"

#include "Epetra_SerialSpdDenseSolver.h"
#include "Epetra_SerialSymDenseMatrix.h"

#include <sys/time.h>
#include <random>

/*----------------------------------------------------------------------*/
INVANA::SurfCurrentGroup::SurfCurrentGroup(Teuchos::RCP<DRT::Discretization> discret)
    : sourcedis_(discret), targetdis_(Teuchos::null)
{
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  // set up source discretization
  if (not sourcedis_->Filled() || not sourcedis_->HaveDofs())
    dserror("Such a discretization should not end up here!");

  // set up target discretization
  DRT::Problem* reference_problem = ReadReferenceDiscretization();

  targetdis_ = reference_problem->GetDis("structure");
  if (not targetdis_->Filled() || not targetdis_->HaveDofs()) targetdis_->FillComplete();

  // get the conditions for the current evaluation
  std::vector<DRT::Condition*> scc_source;
  std::vector<DRT::Condition*> scc_target;
  sourcedis_->GetCondition("SurfaceCurrent", scc_source);
  targetdis_->GetCondition("SurfaceCurrent", scc_target);

  // check wether length of both conditions is equal:
  if (scc_source.size() != scc_target.size())
    dserror(
        "Size of conditions in source and target differs. We need the size to be equal for "
        "meaningful computation");

  // initialize parallel environment for the computation of each surface current pair
  Kokkos::initialize();
#if defined(KOKKOS_HAVE_OPENMP)
  std::cout << "Surface Currents in OpenMP-mode" << std::endl;
#else
  std::cout << "Surface Currents in Serial-mode" << std::endl;
#endif

  // find corresponding conditions in source and target discretization
  // and initialize the associated objective functions
  std::vector<Teuchos::RCP<SurfCurrentPair>> dump;
  currents_.push_back(dump);
  for (int i = 0; i < (int)scc_source.size(); i++)
  {
    bool foundit = false;
    int ids = scc_source[i]->GetInt("matching id");
    for (int j = 0; j < (int)scc_target.size(); j++)
    {
      int idt = scc_target[j]->GetInt("matching id");
      if (idt == ids)
      {
        // add to "0" since the loop of multiple measurements over timestill missing
        currents_[0].push_back(Teuchos::rcp(
            new SurfCurrentPair(sourcedis_, targetdis_, scc_source[i], scc_target[j])));
        foundit = true;
        break;
      }
    }

    if (foundit)
      continue;
    else
      dserror("corresponding condition in target not found");
  }

  // currents so far only for the last time step;
  timesteps_.push_back(sdyn.get<double>("MAXTIME"));

  if (currents_[0].size() != scc_source.size())
    dserror(
        "problem in finding corresponding current conditions in source and target discretization");
}

DRT::Problem* INVANA::SurfCurrentGroup::ReadReferenceDiscretization()
{
  // groups stuff in case
  Teuchos::RCP<Epetra_Comm> gcomm = DRT::Problem::Instance()->GetNPGroup()->GlobalComm();
  Teuchos::RCP<Epetra_Comm> lcomm = DRT::Problem::Instance()->GetNPGroup()->LocalComm();
  int ngroups = DRT::Problem::Instance()->GetNPGroup()->NumGroups();
  int groupid = DRT::Problem::Instance()->GetNPGroup()->GroupId();

  // get target discretization from the input file
  const Teuchos::ParameterList& statinvp = DRT::Problem::Instance()->StatInverseAnalysisParams();
  std::string reference_input_file = statinvp.get<std::string>("TARGETDISCRETIZATION");
  if (reference_input_file.compare("none.dat"))
  {
    // check wether absolut path is given and prepend if not
    if (reference_input_file[0] != '/')
    {
      std::string filename = DRT::Problem::Instance()->OutputControlFile()->InputFileName();
      std::string::size_type pos = filename.rfind('/');
      if (pos != std::string::npos)
      {
        std::string path = filename.substr(0, pos + 1);
        reference_input_file.insert(reference_input_file.begin(), path.begin(), path.end());
      }
    }
  }
  else
    dserror("You forgot to specifiy the filename for the target discretization");

  // a new problem instance with npgroup
  int newinstance = DRT::Problem::NumInstances();
  DRT::Problem* reference_problem = DRT::Problem::Instance(newinstance);
  reference_problem->NPGroup(*(DRT::Problem::Instance()->GetNPGroup()));

  // a reader
  DRT::INPUT::DatFileReader refreader(reference_input_file, lcomm, 1);

  // read parameter, materials, functions
  reference_problem->ReadParameter(refreader);

  DRT::Problem::Instance()->Materials()->SetReadFromProblem(1);
  reference_problem->ReadMaterials(refreader);

  reference_problem->ReadTimeFunctionResult(refreader);

  if (groupid == 0)
  {
    // read nodes, elements and conditions
    reference_problem->ReadFields(refreader);

    // read conditions
    reference_problem->ReadConditions(refreader);
  }

  // Broadcast to all groups
  gcomm->Barrier();
  if (ngroups > 1) COMM_UTILS::BroadcastDiscretizations(newinstance);
  gcomm->Barrier();

  // reset were materials are read from
  DRT::Problem::Instance()->Materials()->ResetReadFromProblem();

  return reference_problem;
}

/*----------------------------------------------------------------------*/
int INVANA::SurfCurrentGroup::FindStep(double time)
{
  // find step of the evaluation according to time:
  int step = -1;
  for (int i = 0; i < (int)timesteps_.size(); i++)
  {
    double dt = abs(timesteps_[i] - time);
    if (dt < 1.0e-10) step = i;
  }

  return step;
}

/*----------------------------------------------------------------------*/
void INVANA::SurfCurrentGroup::Evaluate(Teuchos::RCP<Epetra_Vector> state, double time, double& val)
{
  int step = FindStep(time);  // not needed so far
  if (step == -1) dserror("no measurements for the requested time step");

  val = 0.0;

  // so far only the case were one single "measurement" for the final timestep of the simulation is
  // considered
  sourcedis_->SetState("displacements", state);

  struct timeval begin, end;
  gettimeofday(&begin, NULL);

  for (int i = 0; i < (int)currents_[0].size(); ++i)
  {
    // evaluate every single surface combination
    val += currents_[0][i]->WSpaceNorm();
  }
  gettimeofday(&end, NULL);
  double dtime = 1.0 * (end.tv_sec - begin.tv_sec) + 1.0e-6 * (end.tv_usec - begin.tv_usec);

  if (sourcedis_->Comm().MyPID() == 0)
    IO::cout << "SURFACE CURRENT function evaluation took: " << dtime << " seconds" << IO::endl;
}

/*----------------------------------------------------------------------*/
void INVANA::SurfCurrentGroup::EvaluateGradient(
    Teuchos::RCP<Epetra_Vector> state, double time, Teuchos::RCP<Epetra_Vector> gradient)
{
  int step = FindStep(time);  // not needed so far
  if (step == -1) dserror("no measurements for the requested time step");

  gradient->PutScalar(0.0);

  sourcedis_->SetState("displacements", state);

  struct timeval begin, end;
  gettimeofday(&begin, NULL);

  for (int i = 0; i < (int)currents_[0].size(); ++i)
  {
    // evaluate every single surface combination
    currents_[0][i]->GradientWSpaceNorm(gradient);
  }
  gettimeofday(&end, NULL);
  double dtime = 1.0 * (end.tv_sec - begin.tv_sec) + 1.0e-6 * (end.tv_usec - begin.tv_usec);

  if (sourcedis_->Comm().MyPID() == 0)
    IO::cout << "SURFACE CURRENT gradient evaluation took: " << dtime << " seconds" << IO::endl;
}

double INVANA::SurfCurrentGroup::GetScaleFac()
{
  double fac = 0.0;
  for (int i = 0; i < (int)currents_[0].size(); ++i)
  {
    fac += currents_[0][i]->GetScaleFac();
  }

  return fac;
}


/*----------------------------------------------------------------------*/
INVANA::SurfCurrentPair::SurfCurrentPair(Teuchos::RCP<DRT::Discretization> sourcedis,
    Teuchos::RCP<DRT::Discretization> targetdis, DRT::Condition* sourcecond,
    DRT::Condition* targetcond)
    : scalefac_(1.0)
{
  // get the scale of the covariance kernel
  const Teuchos::ParameterList& statinvp = DRT::Problem::Instance()->StatInverseAnalysisParams();
  sigmaW_ = statinvp.get<double>("KERNELSCALE");
  sigmaW_ = 2 * (sigmaW_ * sigmaW_);
  if (sigmaW_ < 1.0e-6) dserror("supposedly the kernelscale is a little too small!");

  // estimation of the variance of the measurement noise
  var_estim_ = statinvp.get<double>("MEASVARESTIM");

  // Setup the triangulations
  tri_target_ = Teuchos::rcp(new INVANA::Triangulation(targetdis, Teuchos::rcp(targetcond, false)));
  tri_source_ = Teuchos::rcp(new INVANA::Triangulation(sourcedis, Teuchos::rcp(sourcecond, false)));

  // want to scale the objective function eventually?
  scaling_ = DRT::INPUT::IntegralValue<bool>(statinvp, "OBJECTIVEFUNCTSCAL");
  if (scaling_) scalefac_ = var_estim_ / (double)tri_target_->NumTris();

  extract_type x_target = tri_target_->Points();
  int xsize = x_target.size();
  trimesh_type x_target_view("X_View_t", xsize);
  trimesh_host_type h_x_target_view = Kokkos::create_mirror_view(x_target_view);

  // Fill View and Copy to device
  ExtractToHView(x_target, h_x_target_view);
  Kokkos::deep_copy(x_target_view, h_x_target_view);

  // normals
  currents_type n_target("N_View_t", xsize);
  Kokkos::parallel_for(xsize, Normal(n_target, x_target_view));
  n_target_ = n_target;

  // centers
  currents_type c_target("C_View_t", xsize);
  Kokkos::parallel_for(xsize, Centers(c_target, x_target_view));
  c_target_ = c_target;

  // TODO: Move this to the objective function base class
  bool noise = DRT::INPUT::IntegralValue<bool>(statinvp, "SYNTHNOISE");
  int seed = statinvp.get<int>("SYNTHNOISESEED");
  if (noise) ApplyNoise(seed);

  double sum = 0.0;
  Kokkos::parallel_reduce(
      xsize, Convolute<currents_type>(c_target, n_target, c_target, n_target, sigmaW_), sum);

  // since the work was done redundantly by all mpi-ranks
  // this value should be the same for all
  preconvtarget_ = sum;
}

/*----------------------------------------------------------------------*/
void INVANA::SurfCurrentPair::ApplyNoise(int seed)
{
  // Bring normals and centers to host
  currents_host_type centers = Kokkos::create_mirror_view(c_target_);
  Kokkos::deep_copy(centers, c_target_);
  currents_host_type normals = Kokkos::create_mirror_view(n_target_);
  Kokkos::deep_copy(normals, n_target_);

  // restore in extract_type format to blurr it
  extract_type extract_centers;
  HViewToExtract(centers, extract_centers);
  extract_type extract_normals;
  HViewToExtract(normals, extract_normals);

  extract_type blurred;
  tri_target_->ApplyNoise(extract_normals, extract_centers, blurred, var_estim_, sigmaW_, seed);

  ExtractToHView(blurred, normals);
  Kokkos::deep_copy(n_target_, normals);

  return;
}

/*----------------------------------------------------------------------*/
double INVANA::SurfCurrentPair::WSpaceNorm()
{
  int myrank = tri_source_->Comm()->MyPID();
  int numrnk = tri_source_->Comm()->NumProc();

  // the structural integrity component of the target is already precomputed
  double val = preconvtarget_;

  // let the triangulation be warped according to the state of whatever, e.g the pde
  // solved by the 'discretization'
  tri_source_->EvaluateWarp();

  //-------------------------------------------------
  // All the points and the data reduced to myrank
  extract_type x_source = tri_source_->Points();
  extract_type disp = tri_source_->Data();
  unsigned xsize = x_source.size();

  // View on all points of the source
  trimesh_type x_source_view("X_View_s", xsize);
  trimesh_host_type h_x_source_view = Kokkos::create_mirror_view(x_source_view);
  ExtractToHView(x_source, h_x_source_view);
  Kokkos::deep_copy(x_source_view, h_x_source_view);

  // View on all the data associated to the points
  trimesh_type d_source_view("D_View", xsize);
  trimesh_host_type h_d_source_view = Kokkos::create_mirror_view(d_source_view);
  ExtractToHView(disp, h_d_source_view);
  Kokkos::deep_copy(d_source_view, h_d_source_view);

  //-------------------------------------------------
  // Normals and Centers
  // push forward points
  Kokkos::parallel_for(xsize, Update(x_source_view, d_source_view));

  // normals
  currents_type n_source("N_View_s", xsize);
  Kokkos::parallel_for(xsize, Normal(n_source, x_source_view));

  // centers
  currents_type c_source("C_View_s", xsize);
  Kokkos::parallel_for(xsize, Centers(c_source, x_source_view));

  //-------------------------------------------------
  // chunks of the data for myrank
  std::pair<currents_type::size_type, currents_type::size_type> mychunk =
      MyIndices(numrnk, myrank, xsize);
  auto my_c_source = Kokkos::subview(c_source, mychunk, Kokkos::ALL());
  auto my_n_source = Kokkos::subview(n_source, mychunk, Kokkos::ALL());
  int my_xsize = my_c_source.dimension_0();

  //-------------------------------------------------
  // Do the convolutions
  double sum = 0.0;
  Kokkos::parallel_reduce(
      my_xsize, Convolute<ViewStride>(my_c_source, my_n_source, c_source, n_source, sigmaW_), sum);
  double gsum = 0.0;
  tri_source_->Comm()->SumAll(&sum, &gsum, 1);
  val += gsum;

  sum = 0.0;
  Kokkos::parallel_reduce(my_xsize,
      Convolute<ViewStride>(my_c_source, my_n_source, c_target_, n_target_, sigmaW_), sum);
  gsum = 0.0;
  tri_source_->Comm()->SumAll(&sum, &gsum, 1);
  val -= 2 * gsum;

  val = val / (2.0 * var_estim_);

  return val;
}

/*----------------------------------------------------------------------*/
void INVANA::SurfCurrentPair::GradientWSpaceNorm(Teuchos::RCP<Epetra_MultiVector> gradient)
{
  int myrank = tri_source_->Comm()->MyPID();
  int numrnk = tri_source_->Comm()->NumProc();

  // let the triangulation be warped according to the state of whatever, e.g the pde
  // solved by the 'discretization'
  tri_source_->EvaluateWarp();

  //-------------------------------------------------
  // All the points and the data reduced to myrank
  extract_type x_source = tri_source_->Points();
  extract_type disp = tri_source_->Data();
  unsigned xsize = x_source.size();

  // View on all points of the source
  trimesh_type x_source_view("X_View_s", xsize);
  trimesh_host_type h_x_source_view = Kokkos::create_mirror_view(x_source_view);
  Teuchos::RCP<Epetra_Map> trimap = Teuchos::rcp(new Epetra_Map(0, 0, *tri_source_->Comm()));
  ExtractToHView(x_source, h_x_source_view, trimap);
  Kokkos::deep_copy(x_source_view, h_x_source_view);

  // View on all the data associated to the points
  trimesh_type d_source_view("D_View", xsize);
  trimesh_host_type h_d_source_view = Kokkos::create_mirror_view(d_source_view);
  ExtractToHView(disp, h_d_source_view);
  Kokkos::deep_copy(d_source_view, h_d_source_view);

  //-------------------------------------------------
  // Normals and Centers
  // push forward points
  Kokkos::parallel_for(xsize, Update(x_source_view, d_source_view));

  // normals
  currents_type n_source("N_View_s", xsize);
  Kokkos::parallel_for(xsize, Normal(n_source, x_source_view));

  // derivative of the normals
  trimesh_vecdata_type dn_source("DN_View_s", xsize);
  Kokkos::parallel_for(xsize, DNormal(x_source_view, dn_source));

  // centers
  currents_type c_source("C_View_s", xsize);
  Kokkos::parallel_for(xsize, Centers(c_source, x_source_view));

  // Initialize the gradient
  trimesh_type grad("G_View", xsize);
  Kokkos::parallel_for(xsize, Init9(grad));

  //-------------------------------------------------
  // chunks of the data for myrank
  std::pair<currents_type::size_type, currents_type::size_type> mychunk =
      MyIndices(numrnk, myrank, xsize);
  auto my_c_source = Kokkos::subview(c_source, mychunk, Kokkos::ALL());
  auto my_n_source = Kokkos::subview(n_source, mychunk, Kokkos::ALL());
  auto my_dn_source = Kokkos::subview(dn_source, mychunk, Kokkos::ALL());
  int my_xsize = my_c_source.dimension_0();

  int off = mychunk.first;
  Kokkos::parallel_for(my_xsize, ConvoluteDN<ViewStride, ViewStride>(my_c_source, my_dn_source,
                                     c_source, n_source, grad, off, sigmaW_, 2.0));
  Kokkos::parallel_for(my_xsize, ConvoluteDN<ViewStride, ViewStride>(my_c_source, my_dn_source,
                                     c_target_, n_target_, grad, off, sigmaW_, -2.0));
  Kokkos::parallel_for(my_xsize, ConvoluteDk<ViewStride>(my_c_source, my_n_source, c_source,
                                     n_source, grad, off, sigmaW_, 1.0, true));
  Kokkos::parallel_for(my_xsize, ConvoluteDk<ViewStride>(my_c_source, my_n_source, c_target_,
                                     n_target_, grad, off, sigmaW_, -2.0, false));

  // Bring back to host
  trimesh_host_type h_grad = Kokkos::create_mirror_view(grad);
  Kokkos::deep_copy(h_grad, grad);

  // restore in extract_type format to pass it to the SetData routine
  extract_type grad_data;
  HViewToExtract(h_grad, *trimap, grad_data);

  // Set gradient data into the global gradient vector
  Teuchos::RCP<Epetra_MultiVector> lgradient =
      Teuchos::rcp(new Epetra_MultiVector(gradient->Map(), true));
  tri_source_->SetDataGlobally(grad_data, *trimap, lgradient);

  gradient->Update(1.0 / (2.0 * var_estim_), *lgradient, 1.0);
}

/*----------------------------------------------------------------------*/
std::pair<int, int> INVANA::SurfCurrentPair::MyIndices(int nrnk, int myrnk, int size)
{
  int l = size % nrnk;
  int s = (int)(size - l) / nrnk;

  int lower = myrnk * s;
  int upper = (myrnk + 1) * s;

  // for the "last" rank
  if (myrnk == (nrnk - 1)) upper = size;

  return std::make_pair(lower, upper);
}

/*----------------------------------------------------------------------*/
template <typename host_data_type>
void INVANA::SurfCurrentPair::ExtractToHView(const extract_type& in, host_data_type& out)
{
  dsassert(in.size() == out.dimension_0(), "dimension mismatch");

  // fill the data in the view
  int i = 0;
  for (auto it = in.begin(); it != in.end(); it++)
  {
    for (int j = 0; j < (int)it->second.size(); j++) out(i, j) = it->second[j];
    i++;
  }
}

/*----------------------------------------------------------------------*/
template <typename host_data_type>
void INVANA::SurfCurrentPair::ExtractToHView(
    const extract_type& in, host_data_type& out, Teuchos::RCP<Epetra_Map> map)
{
  dsassert(in.size() == out.dimension_0(), "dimension mismatch");

  std::vector<int> gids;

  // fill the data in the view
  int i = 0;
  for (auto it = in.begin(); it != in.end(); it++)
  {
    gids.push_back(it->first);
    for (int j = 0; j < (int)it->second.size(); j++) out(i, j) = it->second[j];
    i++;
  }

  *map = Epetra_Map(-1, (int)gids.size(), gids.data(), 0, map->Comm());
}

/*----------------------------------------------------------------------*/
template <typename host_data_type>
void INVANA::SurfCurrentPair::HViewToExtract(
    const host_data_type& in, const Epetra_Map& inmap, extract_type& out)
{
  int size0 = in.dimension_0();
  std::vector<double> dum;
  for (int i = 0; i < size0; i++)
  {
    for (unsigned j = 0; j < in.dimension_1(); j++) dum.push_back(in(i, j));

    out.insert(std::pair<int, std::vector<double>>(inmap.GID(i), dum));
    dum.clear();
  }
}

/*----------------------------------------------------------------------*/
template <typename host_data_type>
void INVANA::SurfCurrentPair::HViewToExtract(const host_data_type& in, extract_type& out)
{
  int size0 = in.dimension_0();
  std::vector<double> dum;
  for (int i = 0; i < size0; i++)
  {
    for (unsigned j = 0; j < in.dimension_1(); j++) dum.push_back(in(i, j));

    out.insert(std::pair<int, std::vector<double>>(i, dum));
    dum.clear();
  }
}

/*----------------------------------------------------------------------*/
void INVANA::Triangulation::EvaluatePoints()
{
  // clear data
  points_.clear();
  std::vector<int> gids;

  int lnumele = surface_->Geometry().size();
  int gnumele;
  surface_->Comm()->SumAll(&lnumele, &gnumele, 1);

  for (auto ele = surface_->Geometry().begin(); ele != surface_->Geometry().end(); ++ele)
  {
    Teuchos::RCP<DRT::Element> element = ele->second;

    // exclude non row map elements to be processed
    Teuchos::RCP<DRT::FaceElement> elef = Teuchos::rcp_dynamic_cast<DRT::FaceElement>(element);
    if (discret_->ElementRowMap()->LID(elef->ParentElementId()) == -1) continue;

    int numnodes = element->NumNode();
    int egid = element->Id();

    // fill maps of normals and centers with data:
    // case tri  -> no special treatment
    // case quad -> split into 2 tris and associate to the second tri a gid as gid=egid+maxele
    std::vector<double> x;  // node coordinates

    if (numnodes == 3)
    {
      for (int i = 0; i < 3; ++i)
      {
        x.push_back(element->Nodes()[i]->X()[0]);
        x.push_back(element->Nodes()[i]->X()[1]);
        x.push_back(element->Nodes()[i]->X()[2]);
      }
      points_.insert(std::pair<int, std::vector<double>>(egid, x));
      gids.push_back(egid);
    }
    else if (numnodes == 4)
    {
      // permutations to split the quad
      int perm0[] = {0, 1, 2};  // nodes of the first tri
      int perm1[] = {0, 2, 3};  // nodes of the second tri

      // first tri
      for (int i = 0; i < 3; ++i)
      {
        x.push_back(element->Nodes()[perm0[i]]->X()[0]);
        x.push_back(element->Nodes()[perm0[i]]->X()[1]);
        x.push_back(element->Nodes()[perm0[i]]->X()[2]);
      }
      points_.insert(std::pair<int, std::vector<double>>(egid, x));
      gids.push_back(egid);

      x.clear();
      // second tri
      for (int i = 0; i < 3; ++i)
      {
        x.push_back(element->Nodes()[perm1[i]]->X()[0]);
        x.push_back(element->Nodes()[perm1[i]]->X()[1]);
        x.push_back(element->Nodes()[perm1[i]]->X()[2]);
      }
      points_.insert(std::pair<int, std::vector<double>>(egid + gnumele, x));
      gids.push_back(egid + gnumele);
    }
    else
      dserror("only linear tris or bilinear quads are processed in here");
  }
  // set up new map of triangles (which is unique)
  trimap_ = Teuchos::rcp(new Epetra_Map(-1, (int)gids.size(), gids.data(), 0, *(surface_->Comm())));

  haspoints_ = true;
}

/*----------------------------------------------------------------------*/
void INVANA::Triangulation::EvaluateDofs()
{
  std::vector<int> gids;

  // get total number of face elements in this condition
  int lnumele = surface_->Geometry().size();
  int gnumele;
  surface_->Comm()->SumAll(&lnumele, &gnumele, 1);

  // loop the surface elements and build dof map for the faces
  facetdofstype facetdofs;
  for (auto ele = surface_->Geometry().begin(); ele != surface_->Geometry().end(); ++ele)
  {
    Teuchos::RCP<DRT::Element> element = ele->second;
    int numnodes = element->NumNode();
    int egid = element->Id();

    // get this element's dofs
    DRT::Element::LocationArray la(discret_->NumDofSets());
    element->LocationVector(*discret_, la, false);

    if (numnodes == 3)
    {
      facetdofs.insert(std::pair<int, std::vector<int>>(egid, la[0].lm_));
      gids.push_back(egid);
    }
    else if (numnodes == 4)
    {
      // permutations to split the quad
      int perm0[] = {0, 1, 2};  // nodes of the first tri
      int perm1[] = {0, 2, 3};  // nodes of the second tri

      std::vector<int> dofs;

      // first tri
      for (int i = 0; i < 3; ++i)
      {
        dofs.push_back(la[0].lm_[perm0[i] * 3 + 0]);
        dofs.push_back(la[0].lm_[perm0[i] * 3 + 1]);
        dofs.push_back(la[0].lm_[perm0[i] * 3 + 2]);
      }
      facetdofs.insert(std::pair<int, std::vector<int>>(egid, dofs));
      gids.push_back(egid);

      dofs.clear();
      // second tri
      for (int i = 0; i < 3; ++i)
      {
        dofs.push_back(la[0].lm_[perm1[i] * 3 + 0]);
        dofs.push_back(la[0].lm_[perm1[i] * 3 + 1]);
        dofs.push_back(la[0].lm_[perm1[i] * 3 + 2]);
      }
      facetdofs.insert(std::pair<int, std::vector<int>>(egid + gnumele, dofs));
      gids.push_back(egid + gnumele);
    }
    else
      dserror("only linear tris or bilinear quads are processed in here");
  }
  // set up new map of triangles (which is not unique)
  Epetra_Map trimap(-1, (int)gids.size(), gids.data(), 0, *(surface_->Comm()));

  // make it unique
  DRT::Exporter ex(trimap, *trimap_, trimap_->Comm());
  ex.Export(facetdofs);

  facetmap_ = facetdofs;
  hasdofs_ = true;
}

/*----------------------------------------------------------------------*/
void INVANA::Triangulation::EvaluateWarp()
{
  // clear data
  data_.clear();
  std::vector<int> gids;

  // get state of displacement
  if (!discret_->HasState("displacements"))
    dserror("state needs to be set in advance to be able to gather data");

  Epetra_Vector disp(*(discret_->GetState("displacements")));

  // get total number of face elements in this condition
  int lnumele = surface_->Geometry().size();
  int gnumele;
  surface_->Comm()->SumAll(&lnumele, &gnumele, 1);

  // triangulation
  std::map<int, Teuchos::RCP<DRT::Element>>& geom = surface_->Geometry();
  std::map<int, Teuchos::RCP<DRT::Element>>::iterator ele;
  for (ele = geom.begin(); ele != geom.end(); ++ele)
  {
    Teuchos::RCP<DRT::Element> element = ele->second;
    int numnodes = element->NumNode();
    int egid = element->Id();

    // get this element's dofs
    DRT::Element::LocationArray la(discret_->NumDofSets());
    element->LocationVector(*discret_, la, false);

    // this element's displacements
    std::vector<double> mydisp(la[0].lm_.size());
    DRT::UTILS::ExtractMyValues(disp, mydisp, la[0].lm_);

    // fill maps of normals and centers with data:
    // case tri  -> no special treatment
    // case quad -> split into 2 tris and associate to the second tri a gid as gid=egid+maxele
    std::vector<double> x;  // node data

    if (numnodes == 3)
    {
      for (int i = 0; i < 3; ++i)
      {
        x.push_back(mydisp[i * 3 + 0]);
        x.push_back(mydisp[i * 3 + 1]);
        x.push_back(mydisp[i * 3 + 2]);
      }
      data_.insert(std::pair<int, std::vector<double>>(egid, x));
      gids.push_back(egid);
    }
    else if (numnodes == 4)
    {
      // permutations to split the quad
      int perm0[] = {0, 1, 2};  // nodes of the first tri
      int perm1[] = {0, 2, 3};  // nodes of the second tri

      // first tri
      for (int i = 0; i < 3; ++i)
      {
        x.push_back(mydisp[perm0[i] * 3 + 0]);
        x.push_back(mydisp[perm0[i] * 3 + 1]);
        x.push_back(mydisp[perm0[i] * 3 + 2]);
      }
      data_.insert(std::pair<int, std::vector<double>>(egid, x));
      gids.push_back(egid);

      x.clear();
      // second tri
      for (int i = 0; i < 3; ++i)
      {
        x.push_back(mydisp[perm1[i] * 3 + 0]);
        x.push_back(mydisp[perm1[i] * 3 + 1]);
        x.push_back(mydisp[perm1[i] * 3 + 2]);
      }
      data_.insert(std::pair<int, std::vector<double>>(egid + gnumele, x));
      gids.push_back(egid + gnumele);
    }
    else
      dserror("only linear tris or bilinear quads are processed in here");
  }
  // set up new map of triangles (which is not unique)
  Epetra_Map trimap(-1, (int)gids.size(), gids.data(), 0, *(surface_->Comm()));

  // make it unique
  DRT::Exporter ex(trimap, *trimap_, trimap_->Comm());
  ex.Export(data_);

  hasdata_ = true;
}

/*----------------------------------------------------------------------*/
void INVANA::Triangulation::SetDataGlobally(
    const extract_type& data, const Epetra_Map& datamap, Teuchos::RCP<Epetra_MultiVector> vector)
{
  // reduce assemble vector to myrank
  Epetra_Map bla(
      -1, vector->Map().NumMyElements(), vector->Map().MyGlobalElements(), 0, vector->Comm());
  Epetra_Map dofmap_red(*(LINALG::AllreduceEMap(bla)));

  // reduce facetmap to myrank
  facetdofstype facetmap_red = facetmap_;
  DRT::Exporter ex1(*trimap_, datamap, trimap_->Comm());
  ex1.Export(facetmap_red);

  // assemble everything on myrank
  Epetra_Vector vector_red(dofmap_red, true);
  std::vector<int> dofs;
  for (auto it = data.begin(); it != data.end(); it++)
  {
    auto facet = facetmap_red.find(it->first);
    if (facet == facetmap_red.end()) dserror("Facet %d could not be found", it->first);

    dofs = facet->second;
    int err = vector_red.SumIntoGlobalValues(9, it->second.data(), dofs.data());
    if (err) dserror("Something went wrong!");
  }

  // bring to the MPI parallel layout
  Epetra_Export ex2(dofmap_red, vector->Map());
  int err2 = vector->Export(vector_red, ex2, Add);
  if (err2) dserror("export into global gradient layout failed with code: %d", err2);
}

/*----------------------------------------------------------------------*/
void INVANA::Triangulation::ApplyNoise(const extract_type& normals, const extract_type& centers,
    extract_type& blurred, double variance, double lengthscale, int seed) const
{
  int size = normals.size();

  // random number generator
  std::mt19937 generator(seed);
  std::normal_distribution<double> distribution(0.0, 1.0);

  Epetra_SerialDenseMatrix n_blur(size, 3);
  int k = 0;
  for (auto it = normals.begin(); it != normals.end(); it++)
  {
    n_blur(k, 0) = distribution(generator);
    n_blur(k, 1) = distribution(generator);
    n_blur(k, 2) = distribution(generator);
    k++;
  }

  // set up Precision Matrix
  Epetra_SerialSymDenseMatrix P;
  P.Shape(size);
  P.SetLower();
  int i = 0;
  int j = 0;
  double x1, x2, x3, y1, y2, y3;
  // loop columns
  for (auto it = centers.begin(); it != centers.end(); it++)
  {
    x1 = it->second[0];
    x2 = it->second[1];
    x3 = it->second[2];
    // loop rows
    j = 0;
    for (auto jt = centers.begin(); jt != centers.end(); jt++)
    {
      if (j > i) break;

      y1 = jt->second[0];
      y2 = jt->second[1];
      y3 = jt->second[2];
      P(i, j) = kernel(x1, x2, x3, y1, y2, y3, lengthscale) / variance;
      j++;
    }
    i++;
  }

  // get the covariance as precision^-1 (in place)
  {
    Epetra_SerialSpdDenseSolver solver;
    solver.SetMatrix(P);
    solver.Invert();
  }

  // the same solver cannot call Factor() on an already
  // inverted matrix; so get another one
  // (P is completely filled now but still symmetric
  // and hasn't lost its UPLO specification!)
  Epetra_SerialSpdDenseSolver solver;
  solver.SetMatrix(P);
  int err = solver.Factor();

  if (err) dserror("Surface current noise covariance factorization failed with err %d", err);

  // get the lower factor
  Epetra_SerialDenseMatrix* L = solver.FactoredMatrix();


  // correlate noise
  Epetra_SerialDenseMatrix n_blurred(size, 3);
  L->Multiply(false, n_blur, n_blurred);

  // add noise to normals and put back to extract type format
  std::vector<double> n(3);
  k = 0;
  for (auto it = normals.begin(); it != normals.end(); it++)
  {
    n[0] = it->second[0] + n_blurred(k, 0);
    n[1] = it->second[1] + n_blurred(k, 1);
    n[2] = it->second[2] + n_blurred(k, 2);
    blurred.insert(std::pair<int, std::vector<double>>(it->first, n));
    k++;
  }

  return;
}

/*----------------------------------------------------------------------*/
INVANA::extract_type INVANA::Triangulation::Points()
{
  if (not haspoints_) dserror("Points were not evaluated so far!");

  // copy since we want to keep the original points
  extract_type p(points_);
  CommunicateData(p);
  return p;
}

/*----------------------------------------------------------------------*/
INVANA::extract_type INVANA::Triangulation::Data()
{
  if (not hasdata_) dserror("Data was not evaluated so far!");

  // copy since we want to keep the original data
  extract_type p(data_);
  CommunicateData(p);
  return p;
}

/*----------------------------------------------------------------------*/
INVANA::extract_type INVANA::Triangulation::MyPoints()
{
  if (not haspoints_) dserror("Points were not evaluated so far!");

  return points_;
}

/*----------------------------------------------------------------------*/
INVANA::extract_type INVANA::Triangulation::MyData()
{
  if (not hasdata_) dserror("Data was not evaluated so far!");

  return data_;
}

/*----------------------------------------------------------------------*/
void INVANA::Triangulation::CommunicateData(extract_type& data)
{
  Epetra_Map map_red(*LINALG::AllreduceEMap(*trimap_));

  // build the exporter
  DRT::Exporter ex(*trimap_, map_red, trimap_->Comm());

  // export
  ex.Export(data);
}

/*----------------------------------------------------------------------*/
const Teuchos::RCP<Epetra_Comm> INVANA::Triangulation::Comm() { return surface_->Comm(); }

/*----------------------------------------------------------------------*/
int INVANA::Triangulation::NumTris()
{
  if (not haspoints_) dserror("Points were not evaluated so far!");

  return trimap_->NumGlobalElements();
}
#endif

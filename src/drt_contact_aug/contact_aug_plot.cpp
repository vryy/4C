/*----------------------------------------------------------------------------*/
/*!

\brief unite all necessary methods to generate the data for external plots in
MATLAB, PGFPlot or other tools.

\level 3

\maintainer Matthias Mayr
\date Aug 24, 2017

*/
/*----------------------------------------------------------------------------*/

#include "contact_aug_plot.H"
#include "contact_augmented_strategy.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_contact.H"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_discret_interface.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/epetra_utils.H"

#include "../solver_nonlin_nox/nox_nln_aux.H"
#include "../solver_nonlin_nox/nox_nln_constraint_group.H"
#include "../solver_nonlin_nox/nox_nln_solver_linesearchbased.H"

#include "../drt_structure_new/str_model_evaluator_contact.H"

#include "../drt_io/every_iteration_writer.H"
#include "../drt_io/io_pstream.H"

#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::Plot::Direction::Direction(const Plot& plot)
    : type_(INPAR::CONTACT::PlotDirection::vague),
      split_(INPAR::CONTACT::PlotDirectionSplit::vague),
      plot_(plot)
{ /* empty */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::Direction::ReadInput(const Teuchos::ParameterList& pp)
{
  type_ = Teuchos::getIntegralValue<INPAR::CONTACT::PlotDirection>(pp, "DIRECTION");

  if (type_ == INPAR::CONTACT::PlotDirection::read_from_file)
  {
    const std::string& input_filepath = pp.get<std::string>("INPUT_FILE_NAME");
    const std::string& dir_file = pp.get<std::string>("DIRECTION_FILE");
    const std::string full_dir_file(GetFullFilePath(input_filepath, dir_file));
    from_file_ = ReadSparseVectorFromMatlab(full_dir_file);
  }

  split_ = Teuchos::getIntegralValue<INPAR::CONTACT::PlotDirectionSplit>(pp, "DIRECTION_SPLIT");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::string CONTACT::AUG::Plot::Direction::GetFullFilePath(
    const std::string& input_file, const std::string& dir_file) const
{
  std::string full_file_path(dir_file);

  // make path relative to input file path if it is not an absolute path
  if (dir_file[0] != '/')
  {
    std::string::size_type pos = input_file.rfind('/');
    if (pos != std::string::npos)
    {
      std::string tmp = input_file.substr(0, pos + 1);
      full_file_path.insert(full_file_path.begin(), tmp.begin(), tmp.end());
    }
  }

  return full_file_path;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> CONTACT::AUG::Plot::Direction::ReadSparseVectorFromMatlab(
    const std::string& dir_file) const
{
  Teuchos::RCP<const Epetra_Map> prbdofs = plot_.strat_->ProblemDofs();
  Teuchos::RCP<const Epetra_Map> lmdofs = plot_.strat_->LMDoFRowMapPtr(false);

  Teuchos::RCP<Epetra_Map> full_map = LINALG::MergeMap(prbdofs, lmdofs, false);
  Teuchos::RCP<Epetra_Vector> direction = Teuchos::rcp(new Epetra_Vector(*full_map, true));


  if (dir_file == "none")
    dserror("No direction file name has been provided! Read input = \"%s\"", dir_file.c_str());

  if (plot_.strat_->Comm().NumProc() != 1)
    dserror(
        "A external direction vector can currently only be considered in "
        "serial mode. This is due to the used input format. -- hiermeier");

  std::string ext_dir_file(dir_file);
  ExtendFileName(ext_dir_file, plot_.filepath_);

  FILE* file_ptr = std::fopen(ext_dir_file.c_str(), "r");
  if (not file_ptr) dserror("The file \"%s\" could not be opened!", ext_dir_file.c_str());

  double* dir_vals = direction->Values();
  const int* mygids = direction->Map().MyGlobalElements();
  char cline[100];
  unsigned count = 0;

  while (fgets(cline, 100, file_ptr))
  {
    int gid = 0;

    // in a first attempt only the global id is extracted
    sscanf(cline, "%d", &gid);
    if (gid != mygids[count]) dserror("Global ID mismatch!");

    // fill the vector
    sscanf(cline, "%d %lf", &gid, &dir_vals[count++]);
  }

  if (count != static_cast<unsigned>(direction->Map().NumMyElements()))
    dserror(
        "Size mismatch! Did you specify the correct DIRECTION_FILE? It seems"
        " as the number of rows in your DIRECTION_FILE is less than the "
        "number of rows in the global DoF map.");

  std::fclose(file_ptr);

  return direction;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::Plot::Direction::ExtendFileName(
    std::string& file_name, const std::string& file_path) const
{
  // check if the file name contains a full path or only a single file name
  if (file_name.find('/') == std::string::npos)
  {
    const std::string path_only = IO::ExtractPath(file_path);
    file_name = path_only + file_name;
    return true;
  }

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::Direction::SplitIntoSlaveMasterBody(const Epetra_Vector& dir,
    Teuchos::RCP<Epetra_Vector>& x_dir_ptr, Teuchos::RCP<Epetra_Vector>& y_dir_ptr) const
{
  if (plot_.strat_->ParRedist()) dserror("Parallel redistribution is not supported!");

  const Epetra_Map& slnodes = plot_.strat_->SlRowNodes();
  const Epetra_Map& manodes = plot_.strat_->MaRowNodes();

  if (slnodes.NumMyElements() > 0)
  {
    const DRT::Node* snode = plot_.discret_->gNode(slnodes.GID(0));
    Teuchos::RCP<Epetra_Map> slbody_dofs = FindConnectedDofs(snode, *plot_.discret_);

    x_dir_ptr = Teuchos::rcp(new Epetra_Vector(*slbody_dofs, true));
    LINALG::ExtractMyVector(dir, *x_dir_ptr);
  }
  else
  {
    Epetra_Map empty_map(0, 0, plot_.discret_->Comm());
    x_dir_ptr = Teuchos::rcp(new Epetra_Vector(empty_map, true));
  }

  if (manodes.NumMyElements() > 0)
  {
    const DRT::Node* mnode = plot_.discret_->gNode(manodes.GID(0));
    Teuchos::RCP<Epetra_Map> mabody_dofs = FindConnectedDofs(mnode, *plot_.discret_);

    y_dir_ptr = Teuchos::rcp(new Epetra_Vector(*mabody_dofs, true));
    LINALG::ExtractMyVector(dir, *y_dir_ptr);
  }
  else
  {
    Epetra_Map empty_map(0, 0, plot_.discret_->Comm());
    y_dir_ptr = Teuchos::rcp(new Epetra_Vector(empty_map, true));
  }

  if (x_dir_ptr->Map().NumGlobalElements() + y_dir_ptr->Map().NumGlobalElements() !=
      plot_.strat_->ProblemDofs()->NumGlobalElements())
    dserror(
        "Split into slave and master dofs failed! This function supports "
        "currently only two distinct bodies. Self-contact, as well as contact"
        "between multiple bodies is not supported!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> CONTACT::AUG::Plot::Direction::FindConnectedDofs(
    const DRT::Node* node, const DRT::DiscretizationInterface& discret) const
{
  std::set<int> done_element_ids;
  std::set<int> connected_node_gids;
  connected_node_gids.insert(node->Id());
  std::vector<const DRT::Node*> connected_nodes(1, node);

  const int mypid = discret.Comm().MyPID();

  unsigned i = 0;
  for (;;)
  {
    const DRT::Node* cnode = connected_nodes[i++];

    const DRT::Element* const* adj_eles = cnode->Elements();
    int num_adj_eles = cnode->NumElement();

    for (int e = 0; e < num_adj_eles; ++e)
    {
      const DRT::Element* ele = adj_eles[e];
      const auto echeck = done_element_ids.insert(ele->Id());
      if (echeck.second == false) continue;

      const DRT::Node* const* nodes = ele->Nodes();
      for (int n = 0; n < ele->NumNode(); ++n)
      {
        const DRT::Node* ele_node = nodes[n];
        if (ele_node->Owner() != mypid) continue;

        const auto ncheck = connected_node_gids.insert(ele_node->Id());
        if (ncheck.second) connected_nodes.push_back(ele_node);
      }
    }

    if (i == connected_nodes.size()) break;
  }

  // use a set to get an ascending order of the GIDs
  std::set<int> dof_set;
  std::vector<int> dof_vec;
  for (const DRT::Node* cnode : connected_nodes)
  {
    dof_vec.reserve(3);
    discret.Dof(cnode, dof_vec);

    dof_set.insert(dof_vec.begin(), dof_vec.end());
    dof_vec.clear();
  }

  dof_vec.resize(dof_set.size());
  std::copy(dof_set.begin(), dof_set.end(), dof_vec.begin());

  Teuchos::RCP<Epetra_Map> connected_dof_map =
      Teuchos::rcp(new Epetra_Map(-1, dof_vec.size(), dof_vec.data(), 0, discret.Comm()));

  return connected_dof_map;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::AUG::Plot::Direction::Get(
    const NOX::Solver::Generic& solver) const
{
  Teuchos::RCP<Epetra_Vector> dir_ptr = Teuchos::null;

  switch (type_)
  {
    case INPAR::CONTACT::PlotDirection::current_search_direction:
    {
      // compute direction
      const NOX::Epetra::Vector& curr_x =
          dynamic_cast<const NOX::Epetra::Vector&>(solver.getSolutionGroup().getX());
      const NOX::Epetra::Vector& old_x =
          dynamic_cast<const NOX::Epetra::Vector&>(solver.getPreviousSolutionGroup().getX());

      dir_ptr = Teuchos::rcp(new Epetra_Vector(curr_x.getEpetraVector()));
      CATCH_EPETRA_ERROR(dir_ptr->Update(-1.0, old_x.getEpetraVector(), 1.0));

      break;
    }
    case INPAR::CONTACT::PlotDirection::read_from_file:
    {
      dir_ptr = from_file_;

      break;
    }
    case INPAR::CONTACT::PlotDirection::zero:
    {
      // compute direction
      const NOX::Epetra::Vector& curr_x =
          dynamic_cast<const NOX::Epetra::Vector&>(solver.getSolutionGroup().getX());
      dir_ptr = Teuchos::rcp(new Epetra_Vector(curr_x.getEpetraVector().Map(), true));

      break;
    }
    default:
    {
      dserror("Unsupported INPAR::CONTACT::PlotDirection.");

      exit(EXIT_FAILURE);
    }
  }

  return dir_ptr.getConst();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::Create(Teuchos::ParameterList& nox_params,
    const Teuchos::ParameterList& plot_params, const CONTACT::CoAbstractStrategy* strat)
{
  if (not Activated(plot_params)) return;

  Teuchos::RCP<Plot> contact_plot = Teuchos::rcp(new Plot);
  contact_plot->Init(plot_params, strat);
  contact_plot->Setup();

  Teuchos::ParameterList& p_sol_opt = nox_params.sublist("Solver Options");

  Teuchos::RCP<NOX::Abstract::PrePostOperator> prepost_solver_ptr =
      Teuchos::rcp(new NOX::NLN::Solver::PrePostOp::CONTACT::Plot(contact_plot));

  NOX::NLN::AUX::AddToPrePostOpVector(p_sol_opt, prepost_solver_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::AUG::Plot::Activated(const Teuchos::ParameterList& plot_params)
{
  bool is_active = false;

  const int step = plot_params.get<int>("STEP");
  const int iter = plot_params.get<int>("ITER");

  const INPAR::CONTACT::PlotMode mode =
      Teuchos::getIntegralValue<INPAR::CONTACT::PlotMode>(plot_params, "MODE");

  switch (mode)
  {
    case INPAR::CONTACT::PlotMode::write_single_iteration_of_step:
    {
      is_active = (step != -1 and iter != -1);

      break;
    }
    case INPAR::CONTACT::PlotMode::write_last_iteration_of_step:
    case INPAR::CONTACT::PlotMode::write_each_iteration_of_step:
    {
      is_active = (step != -1);

      break;
    }
    default:
    {
      // stay inactive
      break;
    }
  }

  return is_active;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::Plot::Plot()
    : dir_(*this),
      file_open_mode_(std::ios_base::out | std::ios_base::trunc),
      mode_(INPAR::CONTACT::PlotMode::off),
      func_type_(INPAR::CONTACT::PlotFuncName::vague),
      type_(INPAR::CONTACT::PlotType::vague),
      reference_type_(INPAR::CONTACT::PlotReferenceType::vague),
      format_(INPAR::CONTACT::PlotFileFormat::vague),
      x_type_(INPAR::CONTACT::PlotSupportType::vague),
      y_type_(INPAR::CONTACT::PlotSupportType::vague)
{ /* empty */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::Init(
    const Teuchos::ParameterList& plot_params, const CONTACT::CoAbstractStrategy* strat)
{
  strat_ = dynamic_cast<const CONTACT::AUG::Strategy*>(strat);
  discret_ = plot_params.get<const DRT::DiscretizationInterface*>("DISCRETIZATION");
  model_ = plot_params.get<STR::MODELEVALUATOR::Contact*>("MODELEVALUATOR");

  const int output_precision = plot_params.get<int>("OUTPUT_PRECISION");
  if (output_precision < 0) dserror("The specified output precision must be positive!");
  opt_.output_precision_ = static_cast<unsigned>(output_precision);

  const int res_x = plot_params.get<int>("RESOLUTION_X");
  if (res_x < 0) dserror("The resolution in x-direction must be positive!");
  opt_.resolution_x_ = static_cast<unsigned>(res_x);

  const int res_y = plot_params.get<int>("RESOLUTION_Y");
  if (res_y < 0) dserror("The resolution in y-direction must be positive!");
  opt_.resolution_y_ = static_cast<unsigned>(res_y);

  opt_.min_x_ = plot_params.get<double>("MIN_X");
  opt_.max_x_ = plot_params.get<double>("MAX_X");

  opt_.min_y_ = plot_params.get<double>("MIN_Y");
  opt_.max_y_ = plot_params.get<double>("MAX_Y");

  filepath_ = plot_params.get<std::string>("OUTPUT_FILE_NAME");

  file_open_mode_ =
      Teuchos::getIntegralValue<std::ios_base::openmode>(plot_params, "FILE_OPEN_MODE");

  mode_ = Teuchos::getIntegralValue<INPAR::CONTACT::PlotMode>(plot_params, "MODE");

  type_ = Teuchos::getIntegralValue<INPAR::CONTACT::PlotType>(plot_params, "TYPE");

  x_type_ = Teuchos::getIntegralValue<INPAR::CONTACT::PlotSupportType>(plot_params, "X_TYPE");
  y_type_ = Teuchos::getIntegralValue<INPAR::CONTACT::PlotSupportType>(plot_params, "Y_TYPE");

  func_type_ = Teuchos::getIntegralValue<INPAR::CONTACT::PlotFuncName>(plot_params, "FUNC_NAME");

  reference_type_ =
      Teuchos::getIntegralValue<INPAR::CONTACT::PlotReferenceType>(plot_params, "REFERENCE_TYPE");

  ReadRefPoints(plot_params);

  format_ = Teuchos::getIntegralValue<INPAR::CONTACT::PlotFileFormat>(plot_params, "FILE_FORMAT");

  wgap_node_gid_ = plot_params.get<int>("WGAP_NODE_GID");

  const int step = plot_params.get<int>("STEP");
  do_plot_.step_ = step;

  const int iter = plot_params.get<int>("ITER");
  do_plot_.iter_ = iter;

  curr_step_np_ = plot_params.get<const int*>("CURRENT_STEP");
  if (not curr_step_np_) dserror("The step pointer is NULL!");

  dir_.ReadInput(plot_params);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::Setup()
{
  if (type_ == INPAR::CONTACT::PlotType::scalar) opt_.resolution_x_ = 1;

  if (type_ == INPAR::CONTACT::PlotType::line or type_ == INPAR::CONTACT::PlotType::scalar)
    opt_.resolution_y_ = 1;

  X_.Reshape(opt_.resolution_x_, opt_.resolution_y_);
  Y_.Reshape(opt_.resolution_x_, opt_.resolution_y_);
  Z_.resize(std::max(static_cast<int>(type_), 0),
      LINALG::SerialDenseMatrix(opt_.resolution_x_, opt_.resolution_y_));

  std::vector<double> x;
  LinSpace(opt_.min_x_, opt_.max_x_, opt_.resolution_x_, x);

  std::vector<double> y;
  LinSpace(opt_.min_y_, opt_.max_y_, opt_.resolution_y_, y);

  for (unsigned i = 0; i < opt_.resolution_x_; ++i)
    for (unsigned j = 0; j < opt_.resolution_y_; ++j)
    {
      X_(i, j) = x[i];
      Y_(i, j) = y[j];
    }

  const std::string path = IO::ExtractPath(filepath_);
  const std::string dir_name(IO::ExtractFileName(filepath_) + "_plot");

  filepath_ = path + dir_name;
  IO::CreateDirectory(filepath_, strat_->Comm().MyPID());

  AddFileNameToPath();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::ReadRefPoints(const Teuchos::ParameterList& plot_params)
{
  ref_points_.resize(2, LINALG::Matrix<3, 1>(true));

  ReadRefPoint(plot_params, "FIRST_REF_POINT", ref_points_[0].A());
  ReadRefPoint(plot_params, "SECOND_REF_POINT", ref_points_[1].A());

  ref_points_[0].Print(std::cout);
  ref_points_[1].Print(std::cout);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::ReadRefPoint(
    const Teuchos::ParameterList& plot_params, const std::string& param_name, double* coords) const
{
  std::istringstream input_stream(Teuchos::getNumericStringParameter(plot_params, param_name));
  std::string word;
  char* input;
  unsigned count = 0;
  while (input_stream >> word)
  {
    coords[count++] = std::strtod(word.c_str(), &input);

    if (count > 3) dserror("Too many coordinates!");
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::AddFileNameToPath()
{
  filepath_ += "/" + INPAR::CONTACT::PlotFuncName2String(func_type_) +
               (func_type_ == INPAR::CONTACT::PlotFuncName::weighted_gap or
                           func_type_ == INPAR::CONTACT::PlotFuncName::weighted_gap_gradient or
                           func_type_ == INPAR::CONTACT::PlotFuncName::weighted_gap_mod_gradient
                       ? "_" + std::to_string(wgap_node_gid_)
                       : "");

  switch (mode_)
  {
    case INPAR::CONTACT::PlotMode::write_single_iteration_of_step:
    {
      filepath_ +=
          "_step_" + std::to_string(do_plot_.step_) + "_iter_" + std::to_string(do_plot_.iter_);
      break;
    }
    case INPAR::CONTACT::PlotMode::write_last_iteration_of_step:
    {
      filepath_ += "_step_" + std::to_string(do_plot_.step_) + "_iter_last";
      break;
    }
    case INPAR::CONTACT::PlotMode::write_each_iteration_of_step:
    {
      dserror("Not yet considered!");
      break;
    }
    default:
      break;
  }

  filepath_ += "." + INPAR::CONTACT::PlotType2String(type_) + "." +
               INPAR::CONTACT::PlotFileFormat2String(format_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::NLN::MeritFunction::MeritFctName CONTACT::AUG::Plot::ConvertPlotFuncName2MeritFuncName(
    const enum INPAR::CONTACT::PlotFuncName pfunc_name) const
{
  switch (pfunc_name)
  {
    case INPAR::CONTACT::PlotFuncName::vague:
      return NOX::NLN::MeritFunction::mrtfct_vague;
    case INPAR::CONTACT::PlotFuncName::lagrangian:
      return NOX::NLN::MeritFunction::mrtfct_lagrangian;
    case INPAR::CONTACT::PlotFuncName::infeasibility:
      return NOX::NLN::MeritFunction::mrtfct_infeasibility_two_norm;
    case INPAR::CONTACT::PlotFuncName::energy:
      return NOX::NLN::MeritFunction::mrtfct_energy;
    default:
      return NOX::NLN::MeritFunction::mrtfct_vague;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum CONTACT::AUG::WGapGradientType CONTACT::AUG::Plot::ConvertPlotFuncName2WGapGradientType(
    const enum INPAR::CONTACT::PlotFuncName pfunc_name) const
{
  switch (pfunc_name)
  {
    case INPAR::CONTACT::PlotFuncName::weighted_gap_mod_gradient:
      return WGapGradientType::force_balance;
    case INPAR::CONTACT::PlotFuncName::weighted_gap_gradient:
      return WGapGradientType::constraint_enforcement;
    default:
      return WGapGradientType::vague;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const CONTACT::AUG::Strategy& CONTACT::AUG::Plot::Strategy() const
{
  if (not strat_) dserror("No augmented strategy has been provided!");

  return *strat_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::LinSpace(
    const double a, const double b, const unsigned n, std::vector<double>& res) const
{
  res.clear();

  if (n == 1)
  {
    res.resize(n, a);
    if (a != b)
      IO::cout << "WARNING: LinSpace(a,b,n,res) has been called with different "
                  "values for a and b, even though n is equal to 1! The result res is "
                  "set to a.\n";

    return;
  }

  res.resize(n, a);
  res.back() = b;
  const double step = (b - a) / static_cast<double>(n - 1);

  for (auto it = res.begin() + 1; it != res.end() - 1; ++it) *it = *(it - 1) + step;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::DoPredictor(const NOX::Solver::Generic& solver)
{
  if (do_plot_.iter_ == 0) Do(solver);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::Do(const NOX::Solver::Generic& solver)
{
  switch (mode_)
  {
    case INPAR::CONTACT::PlotMode::write_single_iteration_of_step:
    {
      if (*curr_step_np_ == do_plot_.step_ and solver.getNumIterations() == do_plot_.iter_)
      {
        Execute(solver);
        return;
      }
      break;
    }
    case INPAR::CONTACT::PlotMode::write_last_iteration_of_step:
    {
      // The cast becomes necessary since the member function of NOX::Solver::Generic
      // class misses the const qualifier.
      const NOX::NLN::Solver::LineSearchBased* nln_solver =
          dynamic_cast<const NOX::NLN::Solver::LineSearchBased*>(&solver);
      if (nln_solver)
      {
        if ((nln_solver->getStatus() == NOX::StatusTest::Converged or
                nln_solver->getStatus() == NOX::StatusTest::Failed) and
            *curr_step_np_ == do_plot_.step_)
        {
          Execute(solver);
        }
      }
      break;
    }
    case INPAR::CONTACT::PlotMode::write_each_iteration_of_step:
    {
      dserror("Currently unsupported!");
      exit(EXIT_FAILURE);
    }
    default:
      /* do nothing */
      break;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::Execute(const NOX::Solver::Generic& solver)
{
  // get the reference group
  const NOX::NLN::CONSTRAINT::Group* ref_grp = GetReferenceGroup(solver);

  // copy the reference solution grp
  NOX::NLN::CONSTRAINT::Group plot_grp = *ref_grp;

  // get direction
  Teuchos::RCP<const Epetra_Vector> dir_ptr = dir_.Get(solver);
  const Epetra_Vector& dir = *dir_ptr;

  switch (type_)
  {
    case INPAR::CONTACT::PlotType::scalar:
    {
      PlotScalar(*ref_grp, dir, plot_grp);
      break;
    }
    case INPAR::CONTACT::PlotType::line:
    {
      PlotLine(*ref_grp, dir, plot_grp);
      break;
    }
    case INPAR::CONTACT::PlotType::surface:
    {
      PlotSurface(*ref_grp, dir, plot_grp);
      break;
    }
    case INPAR::CONTACT::PlotType::vector_field_2d:
    {
      PlotVectorField2D(*ref_grp, dir, plot_grp);
      break;
    }
    default:
      dserror("Unsupported plot type!");
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::GetSupportPoints(
    enum INPAR::CONTACT::PlotSupportType stype, LINALG::SerialDenseMatrix& support_mat)
{
  switch (stype)
  {
    case INPAR::CONTACT::PlotSupportType::step_length:
    {
      // see Setup

      break;
    }
    case INPAR::CONTACT::PlotSupportType::characteristic_element_length:
    {
      std::fill(support_mat.A(), support_mat.A() + (support_mat.M() * support_mat.N()),
          CharacteristicInterfaceElementLength(SideType::slave));

      break;
    }
    case INPAR::CONTACT::PlotSupportType::position_angle:
    {
      ComputeAnglePosition();

      X_.Reshape(position_node_id_map_.size(), X_.N());
      unsigned i = 0;
      for (auto an_cit = position_node_id_map_.begin(); an_cit != position_node_id_map_.end();
           ++an_cit, ++i)
      {
        for (unsigned j = 0; j < static_cast<unsigned>(X_.N()); ++j) X_(i, j) = an_cit->first;
      }

      break;
    }
    case INPAR::CONTACT::PlotSupportType::position_distance:
    {
      ComputeDistancePosition();

      X_.Reshape(position_node_id_map_.size(), X_.N());
      unsigned i = 0;
      for (auto an_cit = position_node_id_map_.begin(); an_cit != position_node_id_map_.end();
           ++an_cit, ++i)
      {
        for (unsigned j = 0; j < static_cast<unsigned>(X_.N()); ++j) X_(i, j) = an_cit->first;
      }

      break;
    }
    default:
      dserror("Unknown PlotSupportType.");
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::ComputeDistancePosition()
{
  const LINALG::Matrix<3, 1> ref_pos(ref_points_[0].A(), true);

  const Epetra_Map& slrownodes = strat_->SlRowNodes();
  const unsigned num_my_nodes = slrownodes.NumMyElements();
  const int* node_gids = slrownodes.MyGlobalElements();

  for (unsigned i = 0; i < num_my_nodes; ++i)
  {
    const int gid = node_gids[i];
    DRT::Node* node = discret_->gNode(gid);

    if (not node) dserror("Couldn't find the node with GID %d!", gid);

    LINALG::Matrix<3, 1> distance(node->X(), false);
    distance.Update(-1.0, ref_pos, 1.0);

    const double d_nrm2 = distance.Norm2();
    position_node_id_map_.insert(std::pair<double, int>(d_nrm2, gid));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::ComputeAnglePosition()
{
  LINALG::Matrix<3, 1> ref12(ref_points_[0], false);
  ref12.Update(1.0, ref_points_[1], -1.0);

  const Epetra_Map& slrownodes = strat_->SlRowNodes();
  const unsigned num_my_nodes = slrownodes.NumMyElements();
  const int* node_gids = slrownodes.MyGlobalElements();

  for (unsigned i = 0; i < num_my_nodes; ++i)
  {
    const int gid = node_gids[i];
    DRT::Node* node = discret_->gNode(gid);

    if (not node) dserror("Couldn't find the node with GID %d!", gid);

    const LINALG::Matrix<3, 1> ref3(node->X(), true);
    LINALG::Matrix<3, 1> ref13(ref_points_[0], false);
    ref13.Update(1.0, ref3, -1.0);

    const double iproduct = ref12.Dot(ref13);
    double angle = std::acos(iproduct / (ref13.Norm2() * ref12.Norm2()));

    position_node_id_map_.insert(std::pair<double, int>(angle, gid));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::PlotScalar(const NOX::NLN::CONSTRAINT::Group& ref_grp,
    const Epetra_Vector& dir, NOX::NLN::CONSTRAINT::Group& plot_grp)
{
  IO::cout << "Start evaluation of the scalar data...\n";

  Epetra_Vector step(dir.Map(), true);
  GetSupportPoints(x_type_, X_);

  ModifyStepLength(x_type_, X_(0, 0), dir, step);
  plot_grp.computeX(ref_grp, step, 1.0);

  plot_grp.computeF();
  Y_(0, 0) = GetValue(func_type_, plot_grp);

  WriteLineDataToFile();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::PlotLine(const NOX::NLN::CONSTRAINT::Group& ref_grp,
    const Epetra_Vector& dir, NOX::NLN::CONSTRAINT::Group& plot_grp)
{
  IO::cout << "Start evaluation of the line data...\n";
  GetSupportPoints(x_type_, X_);
  Y_.Reshape(X_.M(), Y_.N());

  double norm_step = -1.0;
  Epetra_Vector step(dir.Map(), true);

  for (int i = 0; i < X_.M(); ++i)
  {
    IO::cout << "alpha = " << X_(i, 0) << IO::endl;
    ModifyStepLength(x_type_, X_(i, 0), dir, step);

    double curr_norm_step = 0.0;
    step.Norm2(&curr_norm_step);
    if (curr_norm_step != norm_step)
    {
      norm_step = curr_norm_step;
      plot_grp.computeX(ref_grp, step, 1.0);
    }

    plot_grp.computeF();
    Y_(i, 0) = GetValue(func_type_, plot_grp, &X_(i, 0), &dir);
  }

  WriteLineDataToFile();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::PlotSurface(const NOX::NLN::CONSTRAINT::Group& ref_grp,
    const Epetra_Vector& dir, NOX::NLN::CONSTRAINT::Group& plot_grp)
{
  IO::cout << "Start evaluation of the surface data...\n";
  GetSupportPoints(x_type_, X_);
  GetSupportPoints(y_type_, Y_);

  Teuchos::RCP<Epetra_Vector> x_dir_ptr = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> y_dir_ptr = Teuchos::null;

  dir_.SplitIntoSurfaceDirections(dir, x_dir_ptr, y_dir_ptr);

  Epetra_Vector step(dir.Map(), true);

  for (int i = 0; i < X_.M(); ++i)
  {
    for (int j = 0; j < X_.N(); ++j)
    {
      IO::cout << "( alpha, beta ) = ( " << X_(i, j) << ", " << Y_(i, j) << " )\n";

      ModifyStepLength(x_type_, X_(i, j), *x_dir_ptr, step);
      ModifyStepLength(y_type_, Y_(i, j), *y_dir_ptr, step);
      plot_grp.computeX(ref_grp, step, 1.0);

      plot_grp.computeF();
      Z_[0](i, j) = GetValue(func_type_, plot_grp);
    }
  }

  WriteSurfaceDataToFile();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::PlotVectorField2D(const NOX::NLN::CONSTRAINT::Group& ref_grp,
    const Epetra_Vector& dir, NOX::NLN::CONSTRAINT::Group& plot_grp)
{
  if (x_type_ != INPAR::CONTACT::PlotSupportType::step_length or
      y_type_ != INPAR::CONTACT::PlotSupportType::step_length)
    dserror(
        "PlotVectorField2D supports currently only the step_length"
        " PlotSupportType!");

  Teuchos::RCP<Epetra_Vector> x_dir_ptr = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> y_dir_ptr = Teuchos::null;

  dir_.SplitIntoSurfaceDirections(dir, x_dir_ptr, y_dir_ptr);

  Epetra_Vector step(dir.Map(), true);
  std::vector<const Epetra_Vector*> dirs(2, NULL);
  dirs[0] = x_dir_ptr.get();
  dirs[1] = y_dir_ptr.get();

  std::vector<double> vec_vals;

  IO::cout << "Start evaluation of the vector field data...\n";
  for (unsigned i = 0; i < opt_.resolution_x_; ++i)
  {
    for (unsigned j = 0; j < opt_.resolution_y_; ++j)
    {
      IO::cout << "( alpha, beta ) = ( " << X_(i, j) << ", " << Y_(i, j) << " )\n";

      ModifyStepLength(x_type_, X_(i, j), *x_dir_ptr, step);
      ModifyStepLength(y_type_, Y_(i, j), *y_dir_ptr, step);

      plot_grp.computeX(ref_grp, step, 1.0);

      plot_grp.computeF();

      GetVectorValues(func_type_, plot_grp, dirs, vec_vals);

      Z_[0](i, j) = vec_vals[0];
      Z_[1](i, j) = vec_vals[1];
    }
  }

  WriteVectorFieldToFile();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::ModifyStepLength(const INPAR::CONTACT::PlotSupportType stype,
    const double alpha, const Epetra_Vector& full_x_dir, Epetra_Vector& mod_step) const
{
  switch (stype)
  {
    case INPAR::CONTACT::PlotSupportType::step_length:
    {
      LINALG::AssembleMyVector(0.0, mod_step, alpha, full_x_dir);

      break;
    }
    default:
    {
      LINALG::AssembleMyVector(0.0, mod_step, 1.0, full_x_dir);

      break;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::WriteLineDataToFile() const
{
  if (strat_->Comm().MyPID() != 0) return;

  const int nlines = IO::CountLinesInFile(filepath_);
  std::ofstream outputfile(filepath_, file_open_mode_);

  switch (format_)
  {
    case INPAR::CONTACT::PlotFileFormat::matlab:
    {
      outputfile << "X-DATA:\n";
      WriteMatrixToFile(outputfile, X_, opt_.output_precision_);

      outputfile << "\n\nY-DATA:\n";
      WriteMatrixToFile(outputfile, Y_, opt_.output_precision_);

      break;
    }
    case INPAR::CONTACT::PlotFileFormat::pgfplot:
    {
      if (file_open_mode_ != (std::ios_base::out | std::ios_base::app) or nlines < 1)
        outputfile << std::setw(24) << "x" << std::setw(24) << "y\n";

      std::vector<const LINALG::SerialDenseMatrix*> columndata(2, NULL);
      columndata[0] = &X_;
      columndata[1] = &Y_;

      WriteColumnDataToFile(outputfile, columndata, opt_.output_precision_);

      break;
    }
    default:
      dserror("The given format is not supported! (enum=%d)", format_);
      exit(EXIT_FAILURE);
  }

  outputfile.close();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::WriteVectorFieldToFile() const
{
  if (strat_->Comm().MyPID() != 0) return;

  std::ofstream outputfile(filepath_);

  switch (format_)
  {
    case INPAR::CONTACT::PlotFileFormat::matlab:
    {
      outputfile << "X-DATA:\n";
      WriteMatrixToFile(outputfile, X_, opt_.output_precision_);

      outputfile << "\n\nY-DATA:\n";
      WriteMatrixToFile(outputfile, Y_, opt_.output_precision_);

      outputfile << "\n\nU-DATA:\n";
      WriteMatrixToFile(outputfile, Z_[0], opt_.output_precision_);

      outputfile << "\n\nV-DATA:\n";
      WriteMatrixToFile(outputfile, Z_[1], opt_.output_precision_);

      break;
    }
    case INPAR::CONTACT::PlotFileFormat::pgfplot:
    {
      outputfile << std::setw(24) << "x" << std::setw(24) << "y" << std::setw(24) << "u"
                 << std::setw(24) << "v\n";

      std::vector<const LINALG::SerialDenseMatrix*> columndata(4, NULL);
      columndata[0] = &X_;
      columndata[1] = &Y_;
      columndata[2] = &Z_[0];
      columndata[3] = &Z_[1];

      WriteColumnDataToFile(outputfile, columndata, opt_.output_precision_);

      break;
    }
    default:
      dserror("The given format is not supported! (enum=%d)", format_);
      exit(EXIT_FAILURE);
  }

  outputfile.close();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename T>
void CONTACT::AUG::WriteColumnDataToFile(
    std::ofstream& outputfile, const std::vector<const T*>& columndata, const unsigned p)
{
  if (not outputfile.is_open()) dserror("The file must be open!");

  if (columndata.size() < 1 or columndata[0] == NULL) return;

  outputfile << std::setprecision(p);
  for (unsigned i = 0; i < static_cast<unsigned>(columndata[0]->M()); ++i)
  {
    for (unsigned j = 0; j < static_cast<unsigned>(columndata[0]->N()); ++j)
    {
      for (unsigned c = 0; c < columndata.size(); ++c)
      {
        outputfile << std::setw(24) << std::scientific << (*columndata[c])(i, j);
      }
      outputfile << "\n";
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::WriteSurfaceDataToFile() const
{
  if (strat_->Comm().MyPID() != 0) return;

  std::ofstream outputfile(filepath_);

  switch (format_)
  {
    case INPAR::CONTACT::PlotFileFormat::matlab:
    {
      outputfile << "X-DATA:\n";
      WriteMatrixToFile(outputfile, X_, opt_.output_precision_);

      outputfile << "\n\nY-DATA:\n";
      WriteMatrixToFile(outputfile, Y_, opt_.output_precision_);

      outputfile << "\n\nZ-DATA:\n";
      WriteMatrixToFile(outputfile, Z_[0], opt_.output_precision_);

      break;
    }
    default:
      dserror("The given format is currently not supported! (enum=%d)", format_);
      exit(EXIT_FAILURE);
  }

  outputfile.close();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const NOX::NLN::CONSTRAINT::Group* CONTACT::AUG::Plot::GetReferenceGroup(
    const NOX::Solver::Generic& solver) const
{
  const NOX::NLN::CONSTRAINT::Group* ref_grp = NULL;

  switch (reference_type_)
  {
    case INPAR::CONTACT::PlotReferenceType::previous_solution:
    {
      ref_grp =
          dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>(&solver.getPreviousSolutionGroup());
      break;
    }
    case INPAR::CONTACT::PlotReferenceType::current_solution:
    {
      if (dir_.type_ == INPAR::CONTACT::PlotDirection::current_search_direction)
        IO::cout << "WARNING: The reference point is the current solution "
                    "point and the direction the current search direction TO this point. "
                    "Is this really what you want to do?\n";

      ref_grp = dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>(&solver.getSolutionGroup());
      break;
    }
    default:
    {
      dserror("Unsupported INPAR::CONTACT::PlotReferenceType!");
      exit(EXIT_FAILURE);
    }
  }

  if (not ref_grp) dserror("A NOX::NLN::CONSTRAINT::Group object is expected!");

  return ref_grp;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::Direction::SplitIntoSurfaceDirections(const Epetra_Vector& dir,
    Teuchos::RCP<Epetra_Vector>& x_dir_ptr, Teuchos::RCP<Epetra_Vector>& y_dir_ptr) const
{
  switch (split_)
  {
    case INPAR::CONTACT::PlotDirectionSplit::displacement_lagrange_multiplier:
    {
      x_dir_ptr = Teuchos::rcp(new Epetra_Vector(*plot_.strat_->ProblemDofs(), true));
      LINALG::ExtractMyVector(dir, *x_dir_ptr);

      y_dir_ptr = Teuchos::rcp(new Epetra_Vector(plot_.strat_->LMDoFRowMap(false), true));
      LINALG::ExtractMyVector(dir, *y_dir_ptr);

      break;
    }
    case INPAR::CONTACT::PlotDirectionSplit::slave_master_displacements:
    {
      SplitIntoSlaveMasterBody(dir, x_dir_ptr, y_dir_ptr);

      break;
    }
    default:
      dserror("Undefined direction split!");
      break;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::Plot::GetValue(const enum INPAR::CONTACT::PlotFuncName functype,
    NOX::NLN::CONSTRAINT::Group& plot_grp, const double* curr_xy, const Epetra_Vector* dir) const
{
  // try to convert the function type into a merit function type
  const NOX::NLN::MeritFunction::MeritFctName mrt_func_type =
      ConvertPlotFuncName2MeritFuncName(functype);

  if (mrt_func_type != NOX::NLN::MeritFunction::mrtfct_vague)
    return plot_grp.GetModelValue(mrt_func_type);
  else
  {
    switch (functype)
    {
      case INPAR::CONTACT::PlotFuncName::weighted_gap:
      {
        const Epetra_Vector& wgap = strat_->GetWeightedGap(CONTACT::AUG::MapType::all_slave_nodes);
        const int dof_gid = MapSlNodeGID2NDofGID(wgap_node_gid_);

        const int dof_lid = wgap.Map().LID(dof_gid);
        if (dof_lid == -1) dserror("Couldn't find the DoF with GID %d.", dof_gid);

        return wgap[dof_lid];
      }
      case INPAR::CONTACT::PlotFuncName::weighted_gap_gradient:
      case INPAR::CONTACT::PlotFuncName::weighted_gap_mod_gradient:
      {
        if (not dir) dserror("You have to provide a direction vector!");

        std::vector<double> grad_val;
        std::vector<const Epetra_Vector*> dir_vec(1, dir);

        GetVectorValues(functype, plot_grp, dir_vec, grad_val);

        return grad_val[0];
      }
      case INPAR::CONTACT::PlotFuncName::weighted_gap_gradient_error:
      {
        model_->EvaluateWeightedGapGradientError();
        return strat_->GetTotalGradientError();
      }
      case INPAR::CONTACT::PlotFuncName::weighted_gap_gradient_nodal_jacobian_error:
      {
        model_->EvaluateWeightedGapGradientError();
        const std::vector<std::pair<int, double>>& nodal_jac_error =
            strat_->GetNodalGradientError_Jacobian();

        return GetNodalErrorAtPosition(curr_xy, nodal_jac_error);
      }
      case INPAR::CONTACT::PlotFuncName::weighted_gap_gradient_nodal_ma_proj_error:
      {
        model_->EvaluateWeightedGapGradientError();
        const std::vector<std::pair<int, double>>& nodal_ma_error =
            strat_->GetNodalGradientError_MaProj();

        return GetNodalErrorAtPosition(curr_xy, nodal_ma_error);
      }
      default:
      {
        dserror("Not yet supported!");
        exit(EXIT_FAILURE);
      }
    }
  }

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::Plot::GetNodalErrorAtPosition(
    const double* pos, const std::vector<std::pair<int, double>>& nodal_error) const
{
  if (not pos)
    dserror(
        "You have to provide the current x/y support value (a.k.a. "
        "angle/distance or any other scalar position value in this case).");

  const int ngid = position_node_id_map_.at(*pos);
  for (const auto& nje : nodal_error)
    if (nje.first == ngid) return nje.second;

  return -1.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::GetVectorValues(const enum INPAR::CONTACT::PlotFuncName functype,
    NOX::NLN::CONSTRAINT::Group& plot_grp, const std::vector<const Epetra_Vector*>& dirs,
    std::vector<double>& vec_vals) const
{
  vec_vals.clear();
  vec_vals.resize(dirs.size(), 0.0);

  switch (functype)
  {
    case INPAR::CONTACT::PlotFuncName::weighted_gap_gradient:
    case INPAR::CONTACT::PlotFuncName::weighted_gap_mod_gradient:
    {
      plot_grp.computeFandJacobian();

      const enum CONTACT::AUG::WGapGradientType wgap_type =
          ConvertPlotFuncName2WGapGradientType(functype);
      GetWGapDirectionGradients(wgap_type, dirs, vec_vals);

      break;
    }
    case INPAR::CONTACT::PlotFuncName::energy_gradient:
    {
      GetEnergyDirectionGradients(dirs, vec_vals);

      break;
    }
    default:
    {
      dserror(
          "The function \"%s\" is not supported for the vector-field "
          "plot.",
          INPAR::CONTACT::PlotFuncName2String(functype).c_str());
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::GetEnergyDirectionGradients(
    const std::vector<const Epetra_Vector*>& dirs, std::vector<double>& grad_vals) const
{
  const std::vector<INPAR::STR::ModelType> without_contact_model(1, model_->Type());
  Teuchos::RCP<Epetra_Vector> str_gradient =
      model_->AssembleForceOfModels(&without_contact_model, true);

  Epetra_Vector curr_dir(str_gradient->Map());

  for (unsigned i = 0; i < grad_vals.size(); ++i)
  {
    curr_dir.Scale(0.0);
    for (int j = 0; j < curr_dir.Map().NumMyElements(); ++j)
    {
      const int str_gid = curr_dir.Map().GID(j);
      const int dir_lid = dirs[i]->Map().LID(str_gid);
      if (dir_lid == -1) continue;

      curr_dir[j] = (*dirs[i])[dir_lid];
    }

    str_gradient->Dot(curr_dir, &grad_vals[i]);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::Plot::GetWGapDirectionGradients(
    const enum CONTACT::AUG::WGapGradientType wgap_type,
    const std::vector<const Epetra_Vector*>& dirs, std::vector<double>& grad_vals) const
{
  if (dirs.size() != grad_vals.size()) dserror("Size mismatch!");

  const unsigned num_vecs = dirs.size();

  Teuchos::RCP<const LINALG::SparseMatrix> wgap_grad_ptr =
      strat_->GetWeightedGapGradient(wgap_type, MapType::all_slave_nodes);
  const LINALG::SparseMatrix& wgap_grad = *wgap_grad_ptr;

  std::vector<Epetra_Vector> wgap_dir_grads(num_vecs, Epetra_Vector(wgap_grad.RangeMap()));
  Epetra_Vector curr_dir(wgap_grad.DomainMap());

  for (unsigned i = 0; i < wgap_dir_grads.size(); ++i)
  {
    curr_dir.Scale(0.0);
    for (int j = 0; j < curr_dir.Map().NumMyElements(); ++j)
    {
      const int slma_gid = curr_dir.Map().GID(j);
      const int dir_lid = dirs[i]->Map().LID(slma_gid);
      if (dir_lid == -1) continue;

      curr_dir[j] = (*dirs[i])[dir_lid];
    }
    wgap_grad.Multiply(false, curr_dir, wgap_dir_grads[i]);
  }

  const int dof_gid = MapSlNodeGID2NDofGID(wgap_node_gid_);
  const int rlid = wgap_grad.RangeMap().LID(dof_gid);
  if (rlid == -1) dserror("Node to NDof mapping failed! ( %d --> %d )", dof_gid, rlid);

  for (unsigned i = 0; i < num_vecs; ++i) grad_vals[i] = wgap_dir_grads[i][rlid];
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONTACT::AUG::Plot::MapSlNodeGID2NDofGID(int node_gid) const
{
  if (!strat_->SlRowNodes().PointSameAs(strat_->SlNormalDoFRowMap(false)))
    dserror("Mapping is not possible!");

  const int node_lid = strat_->SlRowNodes().LID(node_gid);
  return strat_->SlNormalDoFRowMap(false).GID(node_lid);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AUG::Plot::CharacteristicInterfaceElementLength(
    const enum CONTACT::AUG::SideType stype) const
{
  return strat_->CharacteristicInterfaceElementLength(stype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <typename T>
void CONTACT::AUG::WriteMatrixToFile(std::ofstream& outputfile, const T& mat, const unsigned p)
{
  if (not outputfile.is_open()) dserror("The file must be open!");

  outputfile << std::setprecision(p);
  for (unsigned i = 0; i < static_cast<unsigned>(mat.M()); ++i)
  {
    for (unsigned j = 0; j < static_cast<unsigned>(mat.N()); ++j)
    {
      outputfile << std::setw(24) << std::scientific << mat(i, j);
    }
    outputfile << "\n";
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Solver::PrePostOp::CONTACT::Plot::Plot(const Teuchos::RCP<::CONTACT::AUG::Plot>& plot_ptr)
    : plot_ptr_(plot_ptr), plot_(*plot_ptr_)
{ /* empty */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::PrePostOp::CONTACT::Plot::runPreIterate(const NOX::Solver::Generic& solver)
{
  plot_.DoPredictor(solver);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::PrePostOp::CONTACT::Plot::runPostIterate(const NOX::Solver::Generic& solver)
{
  plot_.Do(solver);
}

template void CONTACT::AUG::WriteMatrixToFile<LINALG::SerialDenseMatrix>(
    std::ofstream& outputfile, const LINALG::SerialDenseMatrix& mat, const unsigned precison);
template void CONTACT::AUG::WriteColumnDataToFile<LINALG::SerialDenseMatrix>(
    std::ofstream& outputfile, const std::vector<const LINALG::SerialDenseMatrix*>& mat,
    const unsigned precision);

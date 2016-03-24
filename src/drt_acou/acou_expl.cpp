/*!----------------------------------------------------------------------
\file acou_expl.cpp
\brief Control routine for acoustic explicit time integration.

<pre>
\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
*----------------------------------------------------------------------*/

#include "acou_expl.H"

#ifdef HAVE_DEAL_II

#include "acou_ele.H"
#include "acou_ele_action.H"

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/function.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/timer.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/lac/parallel_vector.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/operators.h>

#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_face.h>

#include <fstream>
#include <iostream>
#include <iomanip>

#include "time_integrators.h"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/acoustic.H"
#include "../drt_mat/acoustic_sol.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../linalg/linalg_mapextractor.H"

namespace ACOU
{
using namespace dealii;

/*----------------------------------------------------------------------*
 |  Constructor (public)                                 schoeder 03/15 |
 *----------------------------------------------------------------------*/
AcouExplicitTimeInt::AcouExplicitTimeInt(
  const Teuchos::RCP<DRT::DiscretizationHDG>&   actdis,
  const Teuchos::RCP<LINALG::Solver>&           solver,
  const Teuchos::RCP<Teuchos::ParameterList>&   params,
  const Teuchos::RCP<IO::DiscretizationWriter>& output
  ):
  AcouTimeInt(actdis,solver,params,output),
  doubleorfloat_(DRT::INPUT::IntegralValue<bool>(*params_,"DOUBLEORFLOAT"))
{
  dealii::MultithreadInfo::set_thread_limit(1);

  if(numdim_==2 && doubleorfloat_==true)
    wave2dd_ = Teuchos::rcp(new WaveEquationProblem<2,double>(actdis,params,Teuchos::rcp(this,false)));
  else if(numdim_==3 && doubleorfloat_==true)
    wave3dd_ = Teuchos::rcp(new WaveEquationProblem<3,double>(actdis,params,Teuchos::rcp(this,false)));
  else if(numdim_==2 && doubleorfloat_==false)
    wave2df_ = Teuchos::rcp(new WaveEquationProblem<2,float>(actdis,params,Teuchos::rcp(this,false)));
  else if(numdim_==3 && doubleorfloat_==false)
    wave3df_ = Teuchos::rcp(new WaveEquationProblem<3,float>(actdis,params,Teuchos::rcp(this,false)));
  else
    dserror("number of dimensions must be 2 or 3 for explicit time integration of acoustic problems");

  // in case time step is reduced through deal cfl, tell baci
  UpdateTimeStepSize();

  output_->WriteMesh(0,0.0);
}

/*----------------------------------------------------------------------*
 |  Desctructor (public)                                 schoeder 03/15 |
 *----------------------------------------------------------------------*/
AcouExplicitTimeInt::~AcouExplicitTimeInt()
{}

/*----------------------------------------------------------------------*
 |  ReadRestart (public)                                 schoeder 06/15 |
 *----------------------------------------------------------------------*/
void ACOU::AcouExplicitTimeInt::ReadRestart(int step)
{
  ACOU::AcouTimeInt::ReadRestart(step);

  // read initial conditions by loop over elements of discret.
  if(numdim_==2 && doubleorfloat_==true)
  {
    wave2dd_->evaluator->read_initial_conditions(discret_, wave2dd_->solutions);
    wave2dd_->set_time_and_step(time_,step_);
  }
  else if(numdim_==3 && doubleorfloat_==true)
  {
    wave3dd_->evaluator->read_initial_conditions(discret_, wave3dd_->solutions);
    wave3dd_->set_time_and_step(time_,step_);
  }
  else if(numdim_==2 && doubleorfloat_==false)
  {
    wave2df_->evaluator->read_initial_conditions(discret_, wave2df_->solutions);
    wave2df_->set_time_and_step(time_,step_);
  }
  else if(numdim_==3 && doubleorfloat_==false)
  {
    wave3df_->evaluator->read_initial_conditions(discret_, wave3df_->solutions);
    wave3df_->set_time_and_step(time_,step_);
  }
}

/*----------------------------------------------------------------------*
 |  SetInitialZeroField (public)                         schoeder 03/15 |
 *----------------------------------------------------------------------*/
void AcouExplicitTimeInt::SetInitialZeroField()
{
  ACOU::AcouTimeInt::SetInitialZeroField();

  // read initial conditions by loop over elements of discret.
  if(numdim_==2 && doubleorfloat_==true)
    wave2dd_->evaluator->read_initial_conditions(discret_, wave2dd_->solutions);
  else if(numdim_==3 && doubleorfloat_==true)
    wave3dd_->evaluator->read_initial_conditions(discret_, wave3dd_->solutions);
  else if(numdim_==2 && doubleorfloat_==false)
    wave2df_->evaluator->read_initial_conditions(discret_, wave2df_->solutions);
  else if(numdim_==3 && doubleorfloat_==false)
    wave3df_->evaluator->read_initial_conditions(discret_, wave3df_->solutions);

  return;
}

/*----------------------------------------------------------------------*
 |  SetInitial Field (public)                            schoeder 03/15 |
 *----------------------------------------------------------------------*/
void AcouExplicitTimeInt::SetInitialField(int startfuncno)
{
  // set the initial field as usual (later, the element information is needed)
  Epetra_SerialDenseVector elevec1, elevec2, elevec3;
  Epetra_SerialDenseMatrix elemat1, elemat2;

  Teuchos::ParameterList initParams;
  initParams.set<int>("action",ACOU::project_field);
  initParams.set<int>("funct",startfuncno);
  initParams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);
  initParams.set<bool>("padaptivity",false);
  initParams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);

  DRT::Element::LocationArray la(2);
  for (int el=0; el<discret_->NumMyColElements();++el)
  {
    elevec1.Scale(0.0);elevec2.Scale(0.0);
    DRT::Element *ele = discret_->lColElement(el);
    ele->LocationVector(*discret_,la,false);

    if (static_cast<std::size_t>(elevec1.M()) != la[0].lm_.size())
      elevec1.Shape(la[0].lm_.size(), 1);
    if (elevec2.M() != discret_->NumDof(1,ele))
      elevec2.Shape(discret_->NumDof(1,ele), 1);

    ele->Evaluate(initParams,*discret_,la[0].lm_,elemat1,elemat2,elevec1,elevec2,elevec3);
  }

  // read initial conditions by loop over elements of discret.
  if(numdim_==2 && doubleorfloat_==true)
    wave2dd_->evaluator->read_initial_conditions(discret_, wave2dd_->solutions);
  else if(numdim_==3 && doubleorfloat_==true)
    wave3dd_->evaluator->read_initial_conditions(discret_, wave3dd_->solutions);
  else if(numdim_==2 && doubleorfloat_==false)
    wave2df_->evaluator->read_initial_conditions(discret_, wave2df_->solutions);
  else if(numdim_==3 && doubleorfloat_==false)
    wave3df_->evaluator->read_initial_conditions(discret_, wave3df_->solutions);

  return;
}

/*----------------------------------------------------------------------*
 | Initialization by given scatra solution vector (pub)  schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouExplicitTimeInt::SetInitialPhotoAcousticField(
                                                       Teuchos::RCP<Epetra_Vector> light,
                                                       Teuchos::RCP<DRT::Discretization> scatradis,
                                                       bool meshconform)
{
  ACOU::AcouTimeInt::SetInitialPhotoAcousticField(light,scatradis,meshconform);


  // read initial conditions by loop over elements of discret.
  if(numdim_==2 && doubleorfloat_==true)
    wave2dd_->evaluator->read_initial_conditions(discret_, wave2dd_->solutions);
  else if(numdim_==3 && doubleorfloat_==true)
    wave3dd_->evaluator->read_initial_conditions(discret_, wave3dd_->solutions);
  else if(numdim_==2 && doubleorfloat_==false)
    wave2df_->evaluator->read_initial_conditions(discret_, wave2df_->solutions);
  else if(numdim_==3 && doubleorfloat_==false)
    wave3df_->evaluator->read_initial_conditions(discret_, wave3df_->solutions);

  return;
} // SetInitialPhotoAcousticField

/*----------------------------------------------------------------------*
 |  Integrate                                            schoeder 03/15 |
 *----------------------------------------------------------------------*/
void AcouExplicitTimeInt::Integrate(Teuchos::RCP<Epetra_MultiVector> history)
{

  if(numdim_==2 && doubleorfloat_==true)
      wave2dd_->run(history);
  else if(numdim_==3 && doubleorfloat_==true)
      wave3dd_->run(history);
  else if(numdim_==2 && doubleorfloat_==false)
      wave2df_->run(history);
  else if(numdim_==3 && doubleorfloat_==false)
      wave3df_->run(history);

  return;
}

/*----------------------------------------------------------------------*
 | NodalPressureField                                    schoeder 04/15 |
 *----------------------------------------------------------------------*/
void AcouExplicitTimeInt::NodalPressureField(Teuchos::RCP<Epetra_Vector> outvec)
{
  // first: get the deal values and put them in the baci elements (output is probably a bit more expensive than usually)
  if(numdim_==2 && doubleorfloat_==true)
      wave2dd_->write_deal_cell_values();
  else if(numdim_==3 && doubleorfloat_==true)
      wave3dd_->write_deal_cell_values();
  else if(numdim_==2 && doubleorfloat_==false)
    wave2df_->write_deal_cell_values();
  else if(numdim_==3 && doubleorfloat_==false)
    wave3df_->write_deal_cell_values();

  // output of solution
  Teuchos::RCP<Epetra_Vector> pressure;

  // get the node vector
  {
    const Epetra_Map* nodemap = discret_->NodeRowMap();
    pressure.reset(new Epetra_Vector(*nodemap));
  }

  // call element routine for interpolate HDG to elements
  Teuchos::ParameterList params;
  params.set<int>("action",ACOU::interpolate_hdg_to_node);
  params.set<INPAR::ACOU::PhysicalType>("physical type",phys_);
  params.set<bool>("padaptivity",false);
  params.set<int>("useacouoptvecs",-1);
  params.set<bool>("writestress",writestress_);

  DRT::Element::LocationArray la(2);

  Epetra_SerialDenseMatrix dummyMat;
  Epetra_SerialDenseVector dummyVec;
  Epetra_SerialDenseVector interpolVec;
  std::vector<unsigned char> touchCount(pressure->MyLength());
  pressure->PutScalar(0.);

  for (int el=0; el<discret_->NumMyColElements();++el)
  {
    DRT::Element *ele = discret_->lColElement(el);
    ele->LocationVector(*discret_,la,false);
    if (interpolVec.M() == 0)
    {
      if(phys_==INPAR::ACOU::acou_solid)
        interpolVec.Resize(ele->NumNode()*(2*numdim_+2+numdim_*numdim_)+2);
      else
        interpolVec.Resize(ele->NumNode()*(numdim_+2)+1);
    }

    ele->Evaluate(params,*discret_,la[0].lm_,dummyMat,dummyMat,interpolVec,dummyVec,dummyVec);

    // sum values on nodes into vectors and record the touch count (build average of values)
    for (int i=0; i<ele->NumNode(); ++i)
    {
      DRT::Node* node = ele->Nodes()[i];
      const int localIndex = pressure->Map().LID(node->Id());

      if (localIndex < 0)
        continue;

      touchCount[localIndex]++;
      (*pressure)[localIndex] += interpolVec(i+numdim_*ele->NumNode());
    }
  }

  for (int i=0; i<pressure->MyLength(); ++i)
    (*pressure)[i] /= touchCount[i];

  for(int i=0; i<pressure->MyLength(); ++i)
    outvec->ReplaceMyValue(i,0,pressure->operator [](i));

  return;
}


/*----------------------------------------------------------------------*
 |  UpdateTimeStepSize                                   schoeder 03/15 |
 *----------------------------------------------------------------------*/
void AcouExplicitTimeInt::UpdateTimeStepSize()
{
  if(numdim_==2 && doubleorfloat_==true)
    dtp_ = wave2dd_->time_step;
  else if(numdim_==3 && doubleorfloat_==true)
    dtp_ = wave3dd_->time_step;
  else if(numdim_==2 && doubleorfloat_==false)
    dtp_ = wave2df_->time_step;
  else if(numdim_==3 && doubleorfloat_==false)
    dtp_ = wave3df_->time_step;
  return;
}

std::string AcouExplicitTimeInt::Name()
{
  std::string name;
  switch(dyna_)
  {
  case INPAR::ACOU::acou_expleuler:
  {
    name = "ExplEuler";
    break;
  }
  case INPAR::ACOU::acou_classrk4:
  {
    name = "ClassicalRK4";
    break;
  }
  case INPAR::ACOU::acou_lsrk45reg2:
  {
    name = "LowStorageRK45Reg2";
    break;
  }
  case INPAR::ACOU::acou_lsrk33reg2:
  {
    name = "LowStorageRK33Reg2";
    break;
  }
  case INPAR::ACOU::acou_lsrk45reg3:
  {
    name = "LowStorageRK45Reg3";
    break;
  }
  case INPAR::ACOU::acou_ssprk:
  {
    name = "StrongStabilityPreservingRK";
    break;
  }
  case INPAR::ACOU::acou_ader:
  {
    name = "ADER";
    break;
  }
  default:
    name = "Expl";
    break;
  }
  return name;
}

/*----------------------------------------------------------------------*
 |  FROM HERE ON: ADAPTED DEAL STUFF                     schoeder 03/15 |
 *----------------------------------------------------------------------*/

template <int dim>
class ExactSolution : public Function<dim>
{
public:
  ExactSolution (const unsigned int numcomponent,
                 const double time,
                 const unsigned int functno) : Function<dim>(numcomponent, time),
    numcomponent(numcomponent),
    functno (functno)
  {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component=0) const;

private:
  const unsigned int numcomponent;
  const int functno;
};

template <int dim>
double ExactSolution<dim>::value (const Point<dim>   &p,
                                  const unsigned int component) const
{
  double xyz[dim];
  for(unsigned d=0; d<dim; ++d)
    xyz[d] = p(d);
  double t = this->get_time();

  double return_value = 0.0;

  if(functno>=0)
    return_value = DRT::Problem::Instance()->Funct(functno).Evaluate(component,xyz,t,NULL);

  return return_value;
}

template <int dim>
class RightHandSide : public Function<dim>
{
public:
  RightHandSide (const unsigned int n_components,
                 const double time,
                 const unsigned int functno) :
                   Function<dim>(n_components, time),
                   functno(functno)
  {}

  virtual double value (const Point<dim> &p,
                        const unsigned int component=0) const;

  const int functno;
};

template <int dim>
double RightHandSide<dim>::value (const Point<dim>   &p,
                                  const unsigned int component) const
{
  double xyz[dim];
  for(unsigned d=0; d<dim; ++d)
    xyz[d] = p(d);
  double t = this->get_time();

  double return_value = 0.0;

  if(functno>=0)
    return_value = DRT::Problem::Instance()->Funct(functno).Evaluate(component,xyz,t,NULL);

  return return_value;
}

template <int dim>
class DirichletBoundaryFunction : public Function<dim>
{
public:
  DirichletBoundaryFunction (const unsigned int n_components,
                 const double time,
                 Teuchos::RCP<DRT::DiscretizationHDG> discret) :
                   Function<dim>(n_components, time)
  {
    discret->GetCondition("Dirichlet",dirichletBC);
    num_diri = dirichletBC.size();
  }

  virtual double value (const Point<dim> &p,
                        const unsigned int component=0) const;
private:

  std::vector<DRT::Condition*> dirichletBC;
  unsigned num_diri;
};

template <int dim>
double DirichletBoundaryFunction<dim>::value (const Point<dim>   &p,
                                              const unsigned int component) const
{
  double xyz[dim];
  for(unsigned d=0; d<dim; ++d)
    xyz[d] = p(d);
  double t = this->get_time();

  int dimcomp = component % dim;
  int dbcindex = component / dim;

  // warning: only possibilities are "value" and "function" not "curve". if you want to use "curve" in your dbc, implement the evaluation here
  double return_value = 0.0;
  int funct_num = (*(dirichletBC[dbcindex]->template Get<std::vector<int> >("funct")))[0];
  if(funct_num>0) // return value specified by functione evaluation
  {
    return_value = DRT::Problem::Instance()->Funct(funct_num-1).Evaluate(dimcomp,xyz,t,NULL);
  }
  else // return value given by input value "VAL"
  {
    return_value = (*(dirichletBC[dbcindex]->template Get<std::vector<double> >("val")))[dimcomp];
  }

  return return_value;
}

template<int dim, typename Number>
WaveEquationProblem<dim,Number>::WaveEquationProblem(Teuchos::RCP<DRT::DiscretizationHDG>       discretin,
                                              const Teuchos::RCP<Teuchos::ParameterList> params,
                                              Teuchos::RCP<ACOU::AcouExplicitTimeInt>    timeintin)
  :
  pcout (std::cout,
         Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0),
  triangulation(discretin),
  fe(QGaussLobatto<1>(discretin->lRowElement(0)->Degree()+1)),
  fe_post_disp(QGaussLobatto<1>(discretin->lRowElement(0)->Degree()+2)),
  dof_handler(triangulation),
  dof_handler_post_disp(triangulation),
  time(0.0),
  step(0),
  cfl_number (0.3/fe.degree/dim),
  discret(discretin),
  bacitimeint(timeintin),
  invana(params->get<bool>("invana")),
  adjoint(params->get<bool>("adjoint")),
  acouopt(params->get<bool>("acouopt")),
  timereversal(params->get<bool>("timereversal")),
  solid((DRT::INPUT::IntegralValue<INPAR::ACOU::PhysicalType>(*params,"PHYSICAL_TYPE"))==INPAR::ACOU::acou_solid)
{
  make_grid_and_dofs(discret,bacitimeint->AdjointSourceVec());
  const unsigned int fe_degree = discret->lRowElement(0)->Degree();

  // get function numbers for analytic solution and right hand side
  int comps = 1;
  if(solid) comps = dim;
  int sourcetermfuncno = (params->get<int>("SOURCETERMFUNCNO"))-1;
  Teuchos::RCP<Function<dim> > rhs (new RightHandSide<dim>(comps,0.0,sourcetermfuncno));
  Teuchos::RCP<Function<dim> > diriboundary (new DirichletBoundaryFunction<dim>(comps,0.0,discret));

  if (fe_degree==1)
    evaluator.reset(new WaveEquationOperation<dim,1,Number>(dof_handler, discret, diriboundary, rhs, solid, sourcetermfuncno, bacitimeint->AdjointSourceVec()));
  else if (fe_degree==2)
    evaluator.reset(new WaveEquationOperation<dim,2,Number>(dof_handler, discret, diriboundary, rhs, solid, sourcetermfuncno, bacitimeint->AdjointSourceVec()));
  else if (fe_degree==3)
    evaluator.reset(new WaveEquationOperation<dim,3,Number>(dof_handler, discret, diriboundary, rhs, solid, sourcetermfuncno, bacitimeint->AdjointSourceVec()));
  else if (fe_degree==4)
    evaluator.reset(new WaveEquationOperation<dim,4,Number>(dof_handler, discret, diriboundary, rhs, solid, sourcetermfuncno, bacitimeint->AdjointSourceVec()));
  else if (fe_degree==5)
    evaluator.reset(new WaveEquationOperation<dim,5,Number>(dof_handler, discret, diriboundary, rhs, solid, sourcetermfuncno, bacitimeint->AdjointSourceVec()));
  else if (fe_degree==6)
    evaluator.reset(new WaveEquationOperation<dim,6,Number>(dof_handler, discret, diriboundary, rhs, solid, sourcetermfuncno, bacitimeint->AdjointSourceVec()));
  else
    dserror("Only degrees between 1 and 6 are implemented!");

  if(solid)
    solutions.resize(dim*dim+dim+1);
  else
    solutions.resize(dim+1);

  evaluator->initialize_dof_vector(solutions[0]);
  for (unsigned int d=1; d<solutions.size(); ++d)
    solutions[d] = solutions[0];

  if(solid)
    post_quantity.resize(dim);
  else
    post_quantity.resize(1);

  for (unsigned int d=0; d<post_quantity.size(); ++d)
    post_quantity[d].reinit(dof_handler_post_disp.locally_owned_dofs(),MPI_COMM_WORLD);

  // get some parameters:
  final_time = params->get<double>("MAXTIME");
  up_res = params->get<int>("RESULTSEVRY", -1);
  step_max = params->get<int>("NUMSTEP");
  dyna = DRT::INPUT::IntegralValue<INPAR::ACOU::DynamicType>(*params,"TIMEINT");
  exactsolutionfuncno = (params->get<int>("CALCERRORFUNCNO"))-1;

  // init storage vector for optimization of acoustical parameters
  if(invana && adjoint && acouopt)
  {
    stored_forward_solutions.resize(up_res+1);
    for(unsigned int i=0; i<stored_forward_solutions.size(); ++i)
      stored_forward_solutions[i].resize(dim+1);
  }
  // calculation of the time step size by use of the smallest cell
  typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(),endc = triangulation.end();
  double min_cell_diameter = std::numeric_limits<double>::max();
  double diameter = 0.0;
  double given_time_step = params->get<double>("TIMESTEP");

  for (; cell!=endc; ++cell)
  {
    diameter = cell->diameter()/std::sqrt(dim);

    const int element_index = cell->index();
    Teuchos::RCP<MAT::Material> mat = discret->lColElement(element_index)->Material();
    double c = 0.0;
    if(solid)
    {
      MAT::AcousticSolMat* actmat = static_cast<MAT::AcousticSolMat*>(mat.get());
      c = actmat->SpeedofSound();
    }
    else
    {
      MAT::AcousticMat* actmat = static_cast<MAT::AcousticMat*>(mat.get());
      c = actmat->SpeedofSound(element_index);
    }
    diameter /= c;

    if (diameter < min_cell_diameter)
      min_cell_diameter = diameter;
  }
  const double glob_min_cell_diameter = -Utilities::MPI::max(-min_cell_diameter, MPI_COMM_WORLD);

  time_step = cfl_number * glob_min_cell_diameter;
  time_step = (final_time-time)/(1+int((final_time-time)/time_step));

  if( time_step < given_time_step )
  {
    pcout <<"WARNING: Time step size set to "<<time_step<<" to prevent stability problems"<<std::endl;
    pcout <<"   finest cell / speed of sound "<< glob_min_cell_diameter<< std::endl << std::endl;
    if(invana)
      dserror("set your time step size smaller, otherwise mismatch with monitored values");
  }
  else
  {
    time_step = given_time_step;
    pcout << "Time step size: " << time_step << " (smaller than critical timestep) " << std::endl << std::endl;
  }
}

template <int dim, typename Number>
WaveEquationProblem<dim,Number>::~WaveEquationProblem()
{
  evaluator.reset();
  dof_handler.clear();
}



template<int dim, typename Number>
void WaveEquationProblem<dim,Number>::make_grid_and_dofs(Teuchos::RCP<DRT::DiscretizationHDG> discret, Teuchos::RCP<Epetra_MultiVector> adjoint_source)
{
  if(!invana)
  {
    pcout << "Number of global active cells: "
          << triangulation.n_global_active_cells()
          << std::endl;

    pcout << "Number of global active faces: "
          << triangulation.n_active_faces()
          << std::endl;
  }

  assign_dg_dofs(fe, dof_handler);
  assign_dg_dofs(fe_post_disp, dof_handler_post_disp);

  if(!invana)
  {
    pcout << "Number of degrees of freedom \n"
          << "from discontinuous velocitites and fluxes: "
          << (dim+1)*dof_handler.n_dofs()
          << std::endl;
  }
  typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(),endc = triangulation.end();

  // in the follwing, the boundary conditions are brought to deal via boundary ids
  // 0 ------- absorbing boundary
  // 1 ------- monitored boundary (for inverse analysis in the adjoint run)
  // 2 ------- monitored and absorbing boundary
  // 3 ------- free boundary
  // 4 + x --- dirichlet boundary with corresponding id x

  std::vector<DRT::Condition*> absorbingBC;
  discret->GetCondition("Absorbing",absorbingBC);

  std::vector<DRT::Condition*> dirichletBC;
  discret->GetCondition("Dirichlet",dirichletBC);

  std::vector<DRT::Condition*> pressmonBC;
  discret->GetCondition("PressureMonitor",pressmonBC);

  unsigned is_abc = 0;
  unsigned is_diri = 0;
  unsigned is_pmon = 0;

  const Epetra_Map* nodecolmap = discret->NodeColMap();
  for (; cell!=endc; ++cell)
  {
    for(unsigned f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
    {
      if(cell->face(f)->at_boundary())
      {
        // check, if this face is part of bacis boundary conditions!
        std::vector<int> nodeids(GeometryInfo<dim>::vertices_per_face);
        for(unsigned v=0; v<GeometryInfo<dim>::vertices_per_face; ++v)
          nodeids[v] = nodecolmap->GID(cell->face(f)->vertex_index(v));

        is_diri = 0;
        is_abc = 0;
        is_pmon = 0;

        // check for Dirichlet boundaries
        for(unsigned dbc=0; dbc<dirichletBC.size(); ++dbc)
        {
          for(unsigned v=0;  v<GeometryInfo<dim>::vertices_per_face; ++v)
            is_diri += int(dirichletBC[dbc]->ContainsNode(nodeids[v]));

          if(is_diri>=GeometryInfo<dim>::vertices_per_face)
          {
            cell->face(f)->set_boundary_id(5+dbc);
            break;
          }
        }

        if(is_diri < GeometryInfo<dim>::vertices_per_face)
        {
          // check for Absorbing boundaries
          for(unsigned abc=0; abc<absorbingBC.size(); ++abc)
          {
            for(unsigned v=0;  v<GeometryInfo<dim>::vertices_per_face; ++v)
              is_abc += int(absorbingBC[abc]->ContainsNode(nodeids[v]));
          }
          if(is_abc>=GeometryInfo<dim>::vertices_per_face)
            cell->face(f)->set_boundary_id(0);
        }

        if(is_diri < GeometryInfo<dim>::vertices_per_face) // monitored cannot be where Dirichlet is prescribed but it can be absorbing and monitored
        {
          // check for Pressure Monitor conditions
          for(unsigned pmon=0; pmon<pressmonBC.size(); ++pmon)
          {
            for(unsigned v=0;  v<GeometryInfo<dim>::vertices_per_face; ++v)
              is_pmon += int(pressmonBC[pmon]->ContainsNode(nodeids[v]));
          }
          if(is_pmon >= GeometryInfo<dim>::vertices_per_face && timereversal)
          {
            cell->face(f)->set_boundary_id(4);
            break;
          }
          else if(is_pmon >= GeometryInfo<dim>::vertices_per_face)
          {
            if(is_abc >= GeometryInfo<dim>::vertices_per_face)
             cell->face(f)->set_boundary_id(2);
            else
              cell->face(f)->set_boundary_id(1);
            break;
          }
          else if (is_abc >= GeometryInfo<dim>::vertices_per_face)
            break;
        }
        // if you come here (no break previously) then it is free
        if(is_diri < GeometryInfo<dim>::vertices_per_face)
          cell->face(f)->set_boundary_id(3);
      }
    }
  }
}

template <int dim, typename Number>
void
WaveEquationProblem<dim,Number>::compute_post_pressure()
{
  for (unsigned int d=0; d<dim+1; ++d)
    previous_solutions[d] = 0;
  evaluator->apply(solutions,previous_solutions,time);

  QGauss<dim> quadrature_post(dof_handler_post_disp.get_fe().degree+1);
  FEValues<dim> fe_values(dof_handler.get_fe(),quadrature_post,
                          update_values | update_gradients | update_quadrature_points |
                          update_JxW_values);

  FEValues<dim> fe_values_post(dof_handler_post_disp.get_fe(),quadrature_post,
                               update_values | update_gradients | update_quadrature_points |
                               update_JxW_values);

  const unsigned int dofs_per_cell_post = dof_handler_post_disp.get_fe().dofs_per_cell,
                     dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

  const unsigned int n_q_points = quadrature_post.size();

  LAPACKFullMatrix<double> half_matrix(dofs_per_cell_post*dim,dofs_per_cell_post);
  LAPACKFullMatrix<double> cell_matrix(dofs_per_cell_post,dofs_per_cell_post);
  Vector<double> rhs_vector(dofs_per_cell_post),
         replace_vector(dofs_per_cell_post);
  Table<2,double> solution_vector(dim+1, dofs_per_cell);
  Tensor<1,dim> flux_values;

  std::vector<types::global_dof_index> local_dof_indices_post(dofs_per_cell_post),
      local_dof_indices(dofs_per_cell);
  double replace_value;

  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
  typename DoFHandler<dim>::active_cell_iterator cellpost = dof_handler_post_disp.begin_active(),
                                                 endcpost = dof_handler_post_disp.end();

  for (; cellpost!=endcpost; ++cell, ++cellpost)
    {
      const int element_index = cell->index();
      Teuchos::RCP<MAT::Material> mat = discret->lColElement(element_index)->Material();
      MAT::AcousticMat* actmat = static_cast<MAT::AcousticMat*>(mat.get());
      double rho = actmat->Density();

      cell_matrix = 0;
      rhs_vector = 0.;
      replace_vector = 0.;
      fe_values.reinit(cell);
      fe_values_post.reinit(cellpost);
      cellpost->get_dof_indices(local_dof_indices_post);
      cell->get_dof_indices(local_dof_indices);
      replace_value = 0.;
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          for (unsigned int d=0; d<dim; ++d)
            solution_vector[d][i] = previous_solutions[d](local_dof_indices[i]);
          solution_vector[dim][i] = solutions[dim](local_dof_indices[i]);
        }
      for (unsigned int q_index=0; q_index<n_q_points; ++q_index)
        {
          flux_values = 0.;
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              replace_value += fe_values.shape_value(i,q_index)*solution_vector[dim][i]*fe_values.JxW(q_index);
              for (unsigned int d=0; d<dim; ++d)
                flux_values[d] += fe_values.shape_value(i,q_index) * solution_vector[d][i];
            }
          const double weight = std::sqrt(fe_values_post.JxW(q_index));
          for (unsigned int i=0; i<dofs_per_cell_post; ++i)
            {
              for (unsigned int d=0; d<dim; ++d)
                half_matrix(q_index*dim+d,i) = fe_values_post.shape_grad(i,q_index)[d] * weight;

              rhs_vector(i) -= rho*flux_values *
                               fe_values_post.shape_grad(i,q_index) *
                               fe_values_post.JxW(q_index);
              replace_vector(i) += fe_values_post.shape_value(i,q_index) * 1.0 *
                                   fe_values_post.JxW(q_index);
            }
        }
      half_matrix.Tmmult(cell_matrix, half_matrix);
      for (unsigned int i=0; i<dofs_per_cell_post; ++i)
        {
          cell_matrix(0,i) = replace_vector(i);
          //rhs_vector(i) *= -1.;
        }

      cell_matrix.compute_lu_factorization();

      rhs_vector(0) = replace_value;
      cell_matrix.apply_lu_factorization(rhs_vector, false);
      for (unsigned int i=0; i<dofs_per_cell_post; ++i)
        {
          post_quantity[0].local_element(local_dof_indices_post[i]) = rhs_vector(i);
        }
    }
  post_quantity[0].compress(VectorOperation::insert);
}

template <int dim, typename Number>
void
WaveEquationProblem<dim,Number>::compute_post_velocity()
{
  evaluator->set_adjoint_eval(false);
  for (unsigned int d=0; d<dim*dim+dim+1; ++d)
    previous_solutions[d] = 0;
  evaluator->apply(solutions,previous_solutions,time);

  for(unsigned int vcomp = 0; vcomp<dim; ++vcomp)
  {
    QGauss<dim> quadrature_post(dof_handler_post_disp.get_fe().degree+1);
    FEValues<dim> fe_values(dof_handler.get_fe(),quadrature_post,
        update_values | update_gradients | update_quadrature_points |
        update_JxW_values);

    FEValues<dim> fe_values_post(dof_handler_post_disp.get_fe(),quadrature_post,
        update_values | update_gradients | update_quadrature_points |
        update_JxW_values);

    const unsigned int dofs_per_cell_post = dof_handler_post_disp.get_fe().dofs_per_cell,
                       dofs_per_cell      = dof_handler.get_fe().dofs_per_cell;

    const unsigned int n_q_points = quadrature_post.size();

    LAPACKFullMatrix<double> half_matrix(dofs_per_cell_post*dim,dofs_per_cell_post);
    LAPACKFullMatrix<double> cell_matrix(dofs_per_cell_post,dofs_per_cell_post);
    Vector<double> rhs_vector(dofs_per_cell_post),
                   replace_vector(dofs_per_cell_post);
    Table<2,double> solution_vector(dim+1, dofs_per_cell);
    Tensor<1,dim> flux_values;

    std::vector<types::global_dof_index> local_dof_indices_post(dofs_per_cell_post),
        local_dof_indices(dofs_per_cell);
    double replace_value;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
    typename DoFHandler<dim>::active_cell_iterator cellpost = dof_handler_post_disp.begin_active(),
        endcpost = dof_handler_post_disp.end();

    for (; cellpost!=endcpost; ++cell, ++cellpost)
    {
      cell_matrix = 0;
      rhs_vector = 0.;
      replace_vector = 0.;
      fe_values.reinit(cell);
      fe_values_post.reinit(cellpost);
      cellpost->get_dof_indices(local_dof_indices_post);
      cell->get_dof_indices(local_dof_indices);
      replace_value = 0.;
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
        for (unsigned int d=0; d<dim; ++d)
          solution_vector[d][i] = previous_solutions[dim+1+vcomp*dim+d](local_dof_indices[i]);
        solution_vector[dim][i] = solutions[vcomp](local_dof_indices[i]);
      }
      for (unsigned int q_index=0; q_index<n_q_points; ++q_index)
      {
        flux_values = 0.;
        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          replace_value += fe_values.shape_value(i,q_index)*solution_vector[dim][i]*fe_values.JxW(q_index);
          for (unsigned int d=0; d<dim; ++d)
            flux_values[d] += fe_values.shape_value(i,q_index) * solution_vector[d][i];
        }
        const double weight = std::sqrt(fe_values_post.JxW(q_index));
        for (unsigned int i=0; i<dofs_per_cell_post; ++i)
        {
          for (unsigned int d=0; d<dim; ++d)
            half_matrix(q_index*dim+d,i) = fe_values_post.shape_grad(i,q_index)[d] * weight;

          rhs_vector(i) -= flux_values * fe_values_post.shape_grad(i,q_index) * fe_values_post.JxW(q_index);
          replace_vector(i) += fe_values_post.shape_value(i,q_index) * 1.0 * fe_values_post.JxW(q_index);
        }
      }
      half_matrix.Tmmult(cell_matrix, half_matrix);
      for (unsigned int i=0; i<dofs_per_cell_post; ++i)
      {
        cell_matrix(0,i) = replace_vector(i);
        rhs_vector(i) *= -1.;
      }

      cell_matrix.compute_lu_factorization();

      rhs_vector(0) = replace_value;
      cell_matrix.apply_lu_factorization(rhs_vector, false);

      for (unsigned int i=0; i<dofs_per_cell_post; ++i)
      {
        post_quantity[vcomp].local_element(local_dof_indices_post[i]) = rhs_vector(i);
      }
    }
    post_quantity[vcomp].compress(VectorOperation::insert);
  }
}


namespace
{
  template <int dim>
  void get_relevant_set (const DoFHandler<dim> &dof_handler,
                         IndexSet              &relevant_set)
  {
    relevant_set = dof_handler.locally_owned_dofs();

    std::vector<types::global_dof_index> dof_indices;
    std::set<types::global_dof_index> global_dof_indices;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
        endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      dof_indices.resize(cell->get_fe().dofs_per_cell);
      cell->get_dof_indices(dof_indices);

      for (std::vector<types::global_dof_index>::iterator it =
           dof_indices.begin(); it != dof_indices.end(); ++it)
        if (!relevant_set.is_element(*it))
          global_dof_indices.insert(*it);
    }
    relevant_set.add_indices(global_dof_indices.begin(), global_dof_indices.end());
    relevant_set.compress();
  }
}



template <int dim, typename Number>
void
WaveEquationProblem<dim,Number>::output_results (const unsigned int timestep_number, Teuchos::RCP<Epetra_MultiVector> history)
{

  if (exactsolutionfuncno>=0 && this->solid==false)
  {
    for(int d=0; d<dim+1; ++d)
    {
      Vector<double> norm_per_cell_p (triangulation.n_active_cells());

      IndexSet relevant_set;
      get_relevant_set(dof_handler, relevant_set);
      parallel::distributed::Vector<double> ghosted_sol(dof_handler.locally_owned_dofs(),relevant_set,
                                                        solutions[d].get_mpi_communicator());
      ghosted_sol = solutions[d];
      ghosted_sol.update_ghost_values();

      // calculate norm of field
      VectorTools::integrate_difference (dof_handler,
                                         ghosted_sol,
                                         ZeroFunction<dim>(1),
                                         norm_per_cell_p,
                                         QGauss<dim>(fe.degree+1),
                                         VectorTools::L2_norm);
      double solution_mag = (Utilities::MPI::sum (norm_per_cell_p.norm_sqr(), MPI_COMM_WORLD));

      // calculate norm of difference between field and analytic field
      norm_per_cell_p = 0;
      VectorTools::integrate_difference (dof_handler,
          ghosted_sol,
          ExactSolution<dim>(dim,time,exactsolutionfuncno+d),
          norm_per_cell_p,
          QGauss<dim>(fe.degree+2),
          VectorTools::L2_norm);
      double error_mag = (Utilities::MPI::sum (norm_per_cell_p.norm_sqr(), MPI_COMM_WORLD));

      // calculate norm of analytic solution
      norm_per_cell_p =0;
      ghosted_sol*=0.0;
      VectorTools::integrate_difference (dof_handler,
          ghosted_sol,
          ExactSolution<dim>(dim,time,exactsolutionfuncno+d),
          norm_per_cell_p,
          QGauss<dim>(fe.degree+2),
          VectorTools::L2_norm);
      double exact_sol_mag = (Utilities::MPI::sum (norm_per_cell_p.norm_sqr(), MPI_COMM_WORLD));

      pcout<<"d "<<d<<" absolute error "<<std::sqrt(error_mag)<<" exact solution norm "<<std::sqrt(exact_sol_mag)<<" solution norm "<<std::sqrt(solution_mag)<<std::endl;

    }

    // calculate norm of difference betwenn POSTPROCESSED pressure and analytic solution
    /*Vector<double> norm_per_cell_p (triangulation.n_active_cells());

    compute_post_pressure();

    IndexSet relevant_set;
    get_relevant_set(dof_handler_post_disp, relevant_set);
    parallel::distributed::Vector<double> ghosted_post(dof_handler_post_disp.locally_owned_dofs(),relevant_set,
                                                       solutions[dim].get_mpi_communicator());
    ghosted_post = post_quantity[0];
    ghosted_post.update_ghost_values();
    VectorTools::integrate_difference (dof_handler_post_disp,
                                       ghosted_post,
                                       ExactSolution<dim>(dim,time,exactsolutionfuncno+dim),
                                       norm_per_cell_p,
                                       QGauss<dim>(fe.degree+3),
                                       VectorTools::L2_norm);
    double error_post_mag = (Utilities::MPI::sum (norm_per_cell_p.norm_sqr(), MPI_COMM_WORLD));
    pcout<<"post absolute error "<<std::sqrt(error_post_mag)<<std::endl;
  */

  }
  else if (exactsolutionfuncno>=0 && this->solid==true)
  {
    double v_solution_mag = 0.0;
    double v_error_mag = 0.0;

    for(int d=0; d<dim+1+dim*dim; ++d)
    {

      Vector<double> norm_per_cell_v (triangulation.n_active_cells());

      IndexSet relevant_set;
      get_relevant_set(dof_handler, relevant_set);
      parallel::distributed::Vector<double> ghosted_sol(dof_handler.locally_owned_dofs(),relevant_set,
                                                        solutions[d].get_mpi_communicator());
      ghosted_sol = solutions[d];
      ghosted_sol.update_ghost_values();

      // calculate norm of pressure
      VectorTools::integrate_difference (dof_handler,
                                         ghosted_sol,
                                         ZeroFunction<dim>(1),
                                         norm_per_cell_v,
                                         QGauss<dim>(fe.degree+1),
                                         VectorTools::L2_norm);

      v_solution_mag = Utilities::MPI::sum (norm_per_cell_v.norm_sqr(), MPI_COMM_WORLD);

      // calculate norm of difference between velocity and analytic solution
      norm_per_cell_v = 0;
      VectorTools::integrate_difference (dof_handler,
                                         ghosted_sol,
                                         ExactSolution<dim>(dim,time,exactsolutionfuncno+d),
                                         norm_per_cell_v,
                                         QGauss<dim>(fe.degree+2),
                                         VectorTools::L2_norm);
      v_error_mag  = Utilities::MPI::sum (norm_per_cell_v.norm_sqr(), MPI_COMM_WORLD);

      // calculate norm of difference between velocity and analytic solution
      norm_per_cell_v = 0;
      ghosted_sol*=0.0;
      VectorTools::integrate_difference (dof_handler,
                                         ghosted_sol,
                                         ExactSolution<dim>(dim,time,exactsolutionfuncno+d),
                                         norm_per_cell_v,
                                         QGauss<dim>(fe.degree+2),
                                         VectorTools::L2_norm);
      double v_exactsol_mag  = Utilities::MPI::sum (norm_per_cell_v.norm_sqr(), MPI_COMM_WORLD);


      if(d<dim)
        pcout<<"d "<<d<<" absolute error in v "<<std::sqrt(v_error_mag)<<" exakt solution norm "<< std::sqrt(v_exactsol_mag) <<std::endl;
      else if(d==dim)
        pcout<<"d "<<d<<" absolute error in p "<<std::sqrt(v_error_mag)<<" exakt solution norm "<< std::sqrt(v_exactsol_mag) <<std::endl;
      else
        pcout<<"d "<<d<<" absolute error in H "<<std::sqrt(v_error_mag)<<" exakt solution norm "<< std::sqrt(v_exactsol_mag) <<" solution norm "<< std::sqrt(v_solution_mag)<<std::endl;
    }

    // postprocessed solution
    {
      compute_post_velocity();

      for(int d=0; d<dim; ++d)
      {
        Vector<double> norm_per_cell_v (triangulation.n_active_cells());

        IndexSet relevant_set;
        get_relevant_set(dof_handler_post_disp, relevant_set);
        parallel::distributed::Vector<double> ghosted_sol(dof_handler_post_disp.locally_owned_dofs(),relevant_set,
                                                          solutions[0].get_mpi_communicator());
        ghosted_sol = post_quantity[d];
        ghosted_sol.update_ghost_values();

        // calculate norm of difference between velocity and analytic solution
        norm_per_cell_v = 0;
        VectorTools::integrate_difference (dof_handler_post_disp,
                                           ghosted_sol,
                                           ExactSolution<dim>(dim,time,exactsolutionfuncno+d),
                                           norm_per_cell_v,
                                           QGauss<dim>(fe.degree+3),
                                           VectorTools::L2_norm);
        v_error_mag  = Utilities::MPI::sum (norm_per_cell_v.norm_sqr(), MPI_COMM_WORLD);

        pcout<<"d "<<d<<" absolute error in v****** "<<std::sqrt(v_error_mag)<<std::endl;
      }
    }
    // stress sigma_xy
    /*{
      Vector<double> norm_per_cell_v (triangulation.n_active_cells());

      IndexSet relevant_set;
      get_relevant_set(dof_handler, relevant_set);
      parallel::distributed::Vector<double> ghosted_sol(dof_handler.locally_owned_dofs(),relevant_set,
                                                        solutions[dim+2].get_mpi_communicator());
      ghosted_sol = solutions[dim+2];
      ghosted_sol+=solutions[dim+3];
      ghosted_sol.update_ghost_values();

      // calculate norm of pressure
      VectorTools::integrate_difference (dof_handler,
                                         ghosted_sol,
                                         ZeroFunction<dim>(1),
                                         norm_per_cell_v,
                                         QGauss<dim>(fe.degree+1),
                                         VectorTools::L2_norm);

      v_solution_mag = Utilities::MPI::sum (norm_per_cell_v.norm_sqr(), MPI_COMM_WORLD);

      // calculate norm of difference between velocity and analytic solution
      norm_per_cell_v = 0;
      VectorTools::integrate_difference (dof_handler,
                                         ghosted_sol,
                                         ExactSolution<dim>(dim,time,exactsolutionfuncno+7),
                                         norm_per_cell_v,
                                         QGauss<dim>(fe.degree+2),
                                         VectorTools::L2_norm);
      v_error_mag  = Utilities::MPI::sum (norm_per_cell_v.norm_sqr(), MPI_COMM_WORLD);

      // calculate norm of difference between velocity and analytic solution
      norm_per_cell_v = 0;
      ghosted_sol*=0.0;
      VectorTools::integrate_difference (dof_handler,
                                         ghosted_sol,
                                         ExactSolution<dim>(dim,time,exactsolutionfuncno+7),
                                         norm_per_cell_v,
                                         QGauss<dim>(fe.degree+2),
                                         VectorTools::L2_norm);
      double v_exactsol_mag  = Utilities::MPI::sum (norm_per_cell_v.norm_sqr(), MPI_COMM_WORLD);


     pcout<<"stress! absolute error in stressxy "<<std::sqrt(v_error_mag)<<" solution norm "<< std::sqrt(v_solution_mag)<<" exakt solution norm "<< std::sqrt(v_exactsol_mag) <<std::endl;

    }*/


  }

  /* VTU output */
  // if(!invana && step%up_res == 0) // otherwise we  get way too much output files
  // {
  //   DataOut<dim> data_out;
  //
  //   data_out.attach_dof_handler (dof_handler);
  //   data_out.add_data_vector (ghosted_sol, "solution_pressure");
  //   data_out.build_patches (/*2*/);
  //
  //   const std::string filename_pressure = DRT::Problem::Instance()->OutputControlFile()->FileName() +
  //       ".solution-pressure_" +
  //       Utilities::int_to_string(Utilities::MPI::this_mpi_process(solutions[dim].get_mpi_communicator()))
  //       + "-" + Utilities::int_to_string (timestep_number, 3);
  //
  //   std::ofstream output_pressure ((filename_pressure + ".vtu").c_str());
  //   data_out.write_vtu (output_pressure);
  // }

  // write baci output
  write_deal_cell_values();

  //if(!adjoint)
  {
    if(bacitimeint->Step()%bacitimeint->UpRes()==0)
      bacitimeint->Output(history);
    else
    {
      // deal write node output to history please!
      if( history != Teuchos::null )
      {
        IndexSet relevant_set;
        get_relevant_set(dof_handler, relevant_set);
        parallel::distributed::Vector<double> ghosted_sol(dof_handler.locally_owned_dofs(),relevant_set,
                                                          solutions[dim].get_mpi_communicator());
        ghosted_sol = solutions[dim];
        ghosted_sol.update_ghost_values();

        Point<dim> eval_point;

        // monitor boundary condition
        std::string condname = "PressureMonitor";
        std::vector<DRT::Condition*> pressuremon;
        discret->GetCondition(condname,pressuremon);
        const std::vector<int> pressuremonnodes = *(pressuremon[0]->Nodes());
        for(unsigned int i=0; i<pressuremonnodes.size(); ++i)
        {
          if(discret->NodeRowMap()->LID(pressuremonnodes[i])>=0)
          {
            for(unsigned int d=0; d<dim; ++d)
              eval_point[d] = discret->lRowNode(discret->NodeRowMap()->LID(pressuremonnodes[i]))->X()[d];

            double press = VectorTools::point_value(dof_handler,ghosted_sol,eval_point);

            history->ReplaceMyValue(history->Map().LID(pressuremonnodes[i]),bacitimeint->Step(),press);
          }
        }
      }
    }
  }
  return;
}



template<int dim, typename Number>
void WaveEquationProblem<dim,Number>::run(Teuchos::RCP<Epetra_MultiVector> history)
{
  previous_solutions = solutions;

  // if necessary, write a monitor file
  bacitimeint->InitMonitorFile();

  output_results(0,history);

  Teuchos::RCP<ExplicitIntegrator<WaveEquationOperationBase<dim,Number> > > integrator;
  switch(dyna)
  {
  case INPAR::ACOU::acou_expleuler:
  {
    integrator.reset(new ExplicitEuler<WaveEquationOperationBase<dim,Number> >(*evaluator));
    break;
  }
  case INPAR::ACOU::acou_classrk4:
  {
    integrator.reset(new ClassicalRK4<WaveEquationOperationBase<dim,Number> >(*evaluator));
    break;
  }
  case INPAR::ACOU::acou_lsrk45reg2:
  {
    integrator.reset(new LowStorageRK45Reg2<WaveEquationOperationBase<dim,Number> >(*evaluator));
    break;
  }
  case INPAR::ACOU::acou_lsrk33reg2:
  {
    integrator.reset(new LowStorageRK33Reg2<WaveEquationOperationBase<dim,Number> >(*evaluator));
    break;
  }
  case INPAR::ACOU::acou_lsrk45reg3:
  {
    integrator.reset(new LowStorageRK45Reg3<WaveEquationOperationBase<dim,Number> >(*evaluator));
    break;
  }
  case INPAR::ACOU::acou_ssprk:
  {
    integrator.reset(new StrongStabilityPreservingRK<WaveEquationOperationBase<dim,Number> >(*evaluator, 4, 8));
    break;
  }
  case INPAR::ACOU::acou_ader:
  {
    integrator.reset(new ArbitraryHighOrderDG<WaveEquationOperationBase<dim,Number> >(*evaluator));
    break;
  }
  default:
    dserror("unknown explicit time integration scheme");
    break;
  }

  intermediate_integrate(integrator);

  Timer timer;
  double wtime = 0;
  double output_time = 0;
  step=1;
  evaluator->set_adjoint_eval(invana&&adjoint);

  // for percent output in inverse run
  int maxstep = step_max;
  if(int(final_time/time_step)<maxstep)
    maxstep = final_time/time_step;
  int percent_count = 0;

  // the time loop
  for (time+=time_step; time<=final_time && step<step_max; time+=time_step, ++step)
  {
    bacitimeint->IncrementTimeAndStep();
    timer.restart();
    previous_solutions.swap(solutions);

    // in case of adjoint run with acoustic parameter optimization, calculate the gradient contributions
    if(invana&&adjoint&&acouopt)
      evaluator->compute_gradient_contributions(stored_forward_solutions[(step-1)%up_res],stored_forward_solutions[(step-1)%up_res+1],solutions);

    // in case of backward inverse run, we need the driving force (pressure differences)
    if(bacitimeint->AdjointSourceVec()!=Teuchos::null)
      evaluator->set_timestep_source_number(bacitimeint->AdjointSourceVec()->NumVectors()-step-1);

    // do the actual time integration
    integrator->do_time_step(previous_solutions, solutions, time-time_step, time_step);

    // take the time
    wtime += timer.wall_time();

    // output to screen
    if(invana)
    {
      if( int(step*100)/maxstep > 10*percent_count )
      {
        pcout<<"-> "<< int(step*100)/maxstep-1<<"% ";
        percent_count++;
      }
    }
    else
    {
      pcout<< "TIME: "<< std::setw(10) << std::setprecision(5)<< time <<"/"<< std::setw(10) << std::setprecision(5)<<final_time
           << " DT "  << std::setw(10) << std::setprecision(5) <<time_step
           <<" STEP " << std::setw(6)  << step <<"/"<< std::setw(6) << step_max << std::endl;
    }
    pcout<<std::flush;
    timer.restart();

    // output results
    output_results(step / up_res, history);

    output_time += timer.wall_time();

    // intermediate integration for acoutic adjoint with acouopt
    intermediate_integrate(integrator);
  }
  if(invana)
    pcout<<std::endl;

  if(invana&&adjoint&&acouopt)
    evaluator->write_gradient_contributions(discret,time_step);

  pcout << std::endl
        << "   Performed " << step << " time steps."
        << std::endl;

  pcout << "   Average wallclock time per time step: "
        << wtime / step << "s, time per element: "
        << wtime/step/triangulation.n_active_cells()
        << "s" << std::endl;

  pcout << "   Spent " << output_time << "s on output and "
        << wtime << "s on computations." << std::endl;
}

template<int dim, typename Number>
void WaveEquationProblem<dim,Number>::intermediate_integrate(Teuchos::RCP<ExplicitIntegrator<WaveEquationOperationBase<dim,Number> > > integrator)
{
   // only do this all in case of adjoint run and acoustical parameter optimization
  if(acouopt==false || adjoint==false)
    return;

  // only do this if we are at the right step
  if((step)%up_res != 0)
    return;

  // save the adjoint solution and previous solution vector (such that the adjoint run can continue afterwards without recognizing anything)
  std::vector<parallel::distributed::Vector<value_type> > adjoint_solutions(solutions.size());
  std::vector<parallel::distributed::Vector<value_type> > adjoint_previous_solutions(solutions.size());
  for(unsigned i=0; i<solutions.size(); ++i)
  {
    adjoint_solutions[i] = solutions[i];
    adjoint_previous_solutions[i] = previous_solutions[i];
  }

  // read the restart from the forward run
  double restarttime = 0.0;
  int restartstep = -1;
  {
    // output to screen
    if(!bacitimeint->ProcId())
    {
      pcout<< " <";
      pcout<<std::flush;
    }

    restartstep = (step_max<(final_time/time_step)) ? step_max : (final_time/time_step);
    restartstep -= step + up_res;
    if(restartstep < 0) return; // happens in the last adjoint step
    std::stringstream name;
    name<< bacitimeint->Params()->get<std::string>("name");
    name<<"_invforward_acou_run_"<<bacitimeint->Params()->get<int>("outputcount");
    Teuchos::RCP<IO::InputControl> inputfile = Teuchos::rcp(new IO::InputControl(name.str(), discret->Comm()));
    IO::DiscretizationReader reader(discret,inputfile,restartstep);
    restarttime = reader.ReadDouble("time");

    // read the internal field
    Teuchos::RCP<Epetra_Vector> intvelnp = Teuchos::rcp(new Epetra_Vector(*(discret->DofRowMap(1))));
    reader.ReadVector(intvelnp,"intvelnp");

    // bring the internal field to the baci elements
    Teuchos::ParameterList eleparams;
    eleparams.set<int>("action",ACOU::ele_init_from_restart);
    eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna);
    eleparams.set<bool>("padaptivity",false);
    discret->SetState(1,"intvelnp",intvelnp);
    discret->Evaluate(eleparams);
    discret->ClearState(true);

    // bring the internal field from baci to deal values
    evaluator->read_initial_conditions(discret, solutions);

    // store this
    int count = up_res;
    for(unsigned i=0; i<solutions.size(); ++i)
      stored_forward_solutions[count][i] = solutions[i];

    // do the intermediate time loop
    evaluator->set_adjoint_eval(false);
    int maxstep = restartstep + up_res;
    while (restartstep<maxstep)
    {
      // increment time and step
      restartstep++;
      restarttime += time_step;

      // update state vectors
      previous_solutions.swap(solutions);

      // do the actual time step
      integrator->do_time_step(previous_solutions, solutions, restarttime-time_step, time_step);

      // store the forward solutions
      count--;
      for(unsigned i=0; i<solutions.size(); ++i)
        stored_forward_solutions[count][i] = solutions[i];

    }
  }

  // reset the solution vectors and act if nothing ever happened
  for(unsigned i=0; i<solutions.size(); ++i)
  {
    solutions[i] = adjoint_solutions[i];
    previous_solutions[i] = adjoint_previous_solutions[i];
  }
  evaluator->set_adjoint_eval(true);

  return;
}

template<int dim, typename Number>
void WaveEquationProblem<dim,Number>::write_deal_cell_values()
{
  evaluator->write_deal_cell_values(discret,solutions);
}

template<int dim, typename Number>
void WaveEquationProblem<dim,Number>::set_time_and_step(double timein, int stepin)
{
  step = stepin;
  time = timein;
}

// explicit instantiation
template class WaveEquationProblem<2,double>;
template class WaveEquationProblem<3,double>;
template class WaveEquationProblem<2,float>;
template class WaveEquationProblem<3,float>;
}

#endif // HAVE_DEAL_II

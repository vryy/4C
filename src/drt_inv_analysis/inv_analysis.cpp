/*----------------------------------------------------------------------*/
/*!
 * \file inv_analysis.cpp

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/rausch
            089 - 289-15255
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "inv_analysis.H"
#include <ctime>
#include <cstdlib>
#include <iostream>
#include "Epetra_SerialDenseMatrix.h"
#include "../drt_lib/global_inp_control2.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_io/io_hdf.H"
#include "../drt_lib/linalg_ana.H"

using namespace LINALG::ANA;
using namespace std;

#include "../drt_structure/stru_resulttest.H"



// constructor because it is not copied

Inv_analysis::Inv_analysis(ParameterList& params,
                                  DRT::Discretization& dis,
                                  LINALG::Solver& solver,
                                  IO::DiscretizationWriter& output,
                                  bool init)
  : StruGenAlpha(params,dis,solver,output),
    mat_(DRT::Problem::Instance()->Material(0)),
    nonconstmat_(const_cast<MATERIAL&>(mat_))
{
  // error: diference of the measured to the calculated curve
  error_  = 1.0E6;
  error_o_= 1.5E6;

  //  tolerance for the curve fitting
  tol_ = params.get("inv_ana_tol", 1.0);


  // measured points, gives the number how many displacment steps are measured
  mp_   = params_.get<int>   ("nstep",5);

  // displacment vectors
  final_disp_o_.Resize(mp_);                   // calculated displacment
  final_disp_.Resize(mp_);                     // calculated displacment of the previous run
  measured_disp_.Resize(mp_);                  // measured displacment of the experiment
  residual_disp_.Resize(mp_);                  // difference between the measured and the calculated displacment

  get_measured_disp();                         // reads the measured displacment in

  // trainings parameter
  mu_ = 1;                                     // start value
  mu_minus_ = params.get("mu_minus", 1.0);
  mu_plus_  = params.get("mu_plus" , 1.0);


  // material parameters
  p_.Resize(3);
  p_o_.Resize(3);

  p_(0) = sqrt(DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->c);
  p_(1) = sqrt(DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k1);
  p_(2) = sqrt(DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k2);

  numb_run_=0;

   // get the surface neuman nodes
   discret_.GetCondition("SurfaceNeumann",surfneum_ );

   if (surfneum_.size()==1)
     surfneum_nodes_ = *(surfneum_[0]->Nodes());
   else dserror("The inverse analysis only works for a single suface neumann condition!");
}

void Inv_analysis::Integrate()
{
  int max_itter = 10000;

  do {
    int    step    = params_.get<int>   ("step" ,0);
    double maxtime = params_.get<double>("max time",0.0);
    string equil = params_.get<string>("equilibrium iteration","full newton");

    // can have values takes values "constant" consistent"
    string pred  = params_.get<string>("predictor","constant");
    int predictor=-1;

    if      (pred=="constant")   predictor = 1;
    else if (pred=="consistent") predictor = 2;
    else dserror("Unknown type of predictor");

    //full newton equilibrium
    if (equil=="full newton")
    {
      for (int i=step; i<mp_; ++i)
      {
        if      (predictor==1) ConstantPredictor();
        else if (predictor==2) ConsistentPredictor();

        // gets the displacments per timestep
        get_disp_curve(dis_,  i);

        StruGenAlpha::FullNewton();
        StruGenAlpha::UpdateandOutput();
        double time = params_.get<double>("total time",0.0);
        if (time>=maxtime) break;
      }
    }
    else dserror("Inverse analysis is only implemented for full newton");

    evaluate();

   } while (error_>tol_ && numb_run_<max_itter);

   return;
}


void Inv_analysis::evaluate()
{
  discret_.Comm().Barrier();
  for (int proc=0; proc<discret_.Comm().NumProc(); ++proc)
  {
    if (proc==discret_.Comm().MyPID())
    {
      if (proc == 0)
      {
        cout << "********************************************************************************************" << endl;
        cout << "*******************************\t Inverse Analysis \t************************************" << endl;
        cout << "*******************************\t\t run \t\t************************************" << endl;
        cout << "*******************************\t\t  " << numb_run_ <<"\t\t************************************" << endl;
        cout << "********************************************************************************************" << endl;
        calculate_new_parameters();
      }
    }
  }
  discret_.Comm().Barrier();

  reset_parameters();
  discret_.Comm().Broadcast(&p_[0],  3,  0);

  numb_run_++;

  DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->c      = (p_(0)*p_(0));
  DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k1     = (p_(1)*p_(1));
  DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k2     = (p_(2)*p_(2));

  return;
}


void Inv_analysis::calculate_new_parameters()
{

  // initalization of the Jacobi and storage matrix

  Epetra_SerialDenseMatrix J(mp_, 3);
  Epetra_SerialDenseMatrix storage(3, 3);
  Epetra_SerialDenseMatrix storage2(3, mp_);
  Epetra_SerialDenseVector delta_p(3);

  p_0_storage_.push_back(p_[0]);
  p_1_storage_.push_back(p_[1]);
  p_2_storage_.push_back(p_[2]);

  for (int i=0; i<mp_; i++)
  {
    cout << measured_disp_(i) << "\t" << final_disp_(i) << endl;
    residual_disp_(i) = int(final_disp_(i)) - measured_disp_(i);
  }

  // calculating the errors
  if (abs(error_-error_o_)<0.1)
    error_o_=0.0;
  else
    error_o_ = error_;
  error_ = residual_disp_.Norm1();

  // store the error and mu_
  storage_residual_disp_.push_back(error_);
  storage_mu_.push_back(mu_);

  // print the storage
  cout << "error" << "\t" << "mu_" << "\t" << "c" << "\t" << "k1"<< "\t" << "k2" << endl;
  for (unsigned int i=0; i < storage_residual_disp_.size(); i++)
  {
    cout << storage_residual_disp_[i] << "\t" << storage_mu_[i] << "\t" << p_0_storage_[i]<< "\t" << p_1_storage_[i]<< "\t" << p_2_storage_[i] << endl;
  }

  //calculating J(p)
  for (int i=0; i<mp_; i++)
  {
    for (unsigned int j=0; j<3; j++)
    {
      J(i, j) = (final_disp_[i]-final_disp_o_[i])/(p_[j]-p_o_[j]);
    }
  }

  //calculating J.T*J
  storage.Multiply('T',  'N',  1,  J, J,  0);

  //calculating J.T*J+mu*I
  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    {
      if (i==j)
      {
        storage[i][j] = (storage[i][j] + storage[i][j]* mu_);
      }
    }
  }

  //calculating (J.T*J+mu*I).I
  LINALG::NonSymmetricInverse(storage,  3);

  //calculating (J.T*J+mu*I).I*J.T
  storage2.Multiply('N', 'T', 1,  storage, J, 0);

  //calculating (J.T*J+mu*I).I*J.T*residual_disp_
  delta_p.Multiply('N', 'N', 1,  storage2, residual_disp_, 0);
  cout << "Parameters:" << endl;
  //cout << "old:\t" << p_o_[0] << "\t" << p_o_[1] << "\t" << p_o_[2] << endl;
  //cout << "new:\t" << p_[0] << "\t" << p_[1] << "\t" << p_[2] << endl;

  p_o_ = p_;
  for (int i=0; i<3; i++)
  {
    p_(i)   = p_(i) - delta_p(i);
  }

  if (error_o_<error_) mu_ = mu_plus_;
  else  mu_ = mu_minus_;

  final_disp_o_ = final_disp_;

  if (abs(error_-error_o_)<0.1)
  {
    cout << "Reset old parameters to Null!!!" << endl;
    p_=p_o_;
    p_o_.Scale(0.0);
    final_disp_o_.Scale(0.0);
  }


  cout << "********************************************ENDE********************************************" << endl;

  return;
}


void Inv_analysis::get_disp_curve(const RefCountPtr<Epetra_Vector> disp,  int numb)
{

  int disp_direction = 0;
  double nodal_disp_sum = 0;
  double nodal_disp_sum_2 = 0;


  // check in which direction the surface neumann conditions pulls
  if      (numb==0)
    disp_direction = 0;
  else if (abs((*disp)[surfneum_nodes_[0]]) > abs((*disp)[surfneum_nodes_[1]]) && abs((*disp)[surfneum_nodes_[0]]) > abs((*disp)[surfneum_nodes_[2]]))
    disp_direction = 0;
  else if (abs((*disp)[surfneum_nodes_[1]]) > abs((*disp)[surfneum_nodes_[0]]) && abs((*disp)[surfneum_nodes_[1]]) > abs((*disp)[surfneum_nodes_[2]]))
    disp_direction = 1;
  else if (abs((*disp)[surfneum_nodes_[2]]) > abs((*disp)[surfneum_nodes_[0]]) && abs((*disp)[surfneum_nodes_[2]]) > abs((*disp)[surfneum_nodes_[1]]))
    disp_direction = 2;
  else dserror("The displacment direction is not clear!");

  // summing up the displacments at each of the surface neum nodes
  for (unsigned int i=0; i<surfneum_nodes_.size(); i++)
  {
    if (disp->Map().MyGID(surfneum_nodes_[i]*3 + disp_direction))
    {
      nodal_disp_sum = nodal_disp_sum + (*disp)[disp->Map().LID(surfneum_nodes_[i]*3 + disp_direction)];
    }
  }

  discret_.Comm().Barrier();

  for (int proc=0; proc<discret_.Comm().NumProc(); ++proc)
       {
         if (proc==discret_.Comm().MyPID())
         {
           //cout << "Processor: " << proc << " Nodal Sum: \t" << nodal_disp_sum << endl;
         }
       }
  discret_.Comm().Barrier();

  // summing up over all processors
  discret_.Comm().SumAll(&nodal_disp_sum, &nodal_disp_sum_2, 1);


  //cout << "Nodal Sum 2: \t" << nodal_disp_sum_2 << endl;

  final_disp_[numb] =  nodal_disp_sum_2/surfneum_nodes_.size();



  return;
}

void Inv_analysis::get_measured_disp()
{
  for (int i=0; i<mp_; i++)
  {
    measured_disp_[i] = 1084.625603*(1-exp(-pow((0.003123*(500.0/mp_)*i), 2.281512)));
  }
}


void Inv_analysis::reset_parameters()
{
 // reset the values for next run

  params_.set<int>("step", 0);
  params_.set<double>("total time", 0.0);
  params_.set<int> ("numstep",0);

  dis_->PutScalar(0.0);
  vel_->PutScalar(0.0);
  acc_->PutScalar(0.0);
  disn_->PutScalar(0.0);
  veln_->PutScalar(0.0);
  accn_->PutScalar(0.0);
  dism_->PutScalar(0.0);
  velm_->PutScalar(0.0);
  accm_->PutScalar(0.0);
  disi_->PutScalar(0.0);
  fint_->PutScalar(0.0);
  #ifdef STRUGENALPHA_FINTLIKETR
  fintn_->PutScalar(0.0);
  #endif
  finert_->PutScalar(0.0);
  fvisc_->PutScalar(0.0);
  fext_->PutScalar(0.0);
  fextm_->PutScalar(0.0);
  fextn_->PutScalar(0.0);
  fresm_->PutScalar(0.0);
  frobin_->PutScalar(0.0);

  double iotime = params_.get<double>("total time",0.0);
  int iostep = params_.get<int>("step",0);

  output_.NewResultFile((numb_run_));
  output_.WriteMesh(iostep, iotime);

  return;
}


#endif

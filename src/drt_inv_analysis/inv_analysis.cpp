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
#include "../drt_mat/material.H"


using namespace LINALG::ANA;
using namespace std;
using namespace DRT;
using namespace MAT;


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
   // get the surface neuman nodes
   discret_.GetCondition("SurfaceNeumann",surfneum_ );
   discret_.GetCondition("Dirichlet",surfdir_ );

   if (surfneum_.size()==1 && surfdir_.size()==1)
   {
     problem_type_ = 0;
     bc_nodes_ = *(surfneum_[0]->Nodes());
   }
   else if ((surfneum_.size()==0 && surfdir_.size()==2))
   {
     problem_type_ = 1;
     bc_nodes_ = *(surfdir_[1]->Nodes());
   }
   else dserror("The inverse analysis only works for one Neumann or two Dirichlet Conditions!");

  // measured points, gives the number how many displacment steps are measured
  mp_   = params_.get<int>   ("nstep",5);

  // error: diference of the measured to the calculated curve
  error_  = 1.0E6;
  error_o_= 1.5E6;

  //  tolerance for the curve fitting
  tol_ = params.get("inv_ana_tol", 1.0);

  // displacment vectors
  final_value_o_.Resize(mp_);                   // calculated displacment
  final_value_.Resize(mp_);                     // calculated displacment of the previous run
  measured_value_.Resize(mp_);                  // measured displacment of the experiment
  residual_value_.Resize(mp_);                  // difference between the measured and the calculated displacment

  //set load/displacment curve
  double curve_p0 = params.get("measured_curve0", 1.0);
  double curve_p1 = params.get("measured_curve1", 1.0);
  double curve_p2 = params.get("measured_curve2", 1.0);

  for (int i=0; i<mp_; i++)
  {
    measured_value_[i] = curve_p0*(1-exp(-pow((curve_p1*(500.0/mp_)*i), curve_p2)));
  }

  // trainings parameter
  mu_ = 1;                                     // start value
  mu_minus_ = 0.1;
  mu_plus_  = 10;

  // material parameters
  p_.Resize(3);
  p_o_.Resize(3);
  //Which material is used in the input file
  if (dis.lRowElement(0)->Material()->MaterialType() == m_lung_penalty)
  {
    p_(0) = sqrt(DRT::Problem::Instance()->Material(0).m.lung_penalty->c);
    p_(1) = sqrt(DRT::Problem::Instance()->Material(0).m.lung_penalty->k1);
    p_(2) = sqrt(DRT::Problem::Instance()->Material(0).m.lung_penalty->k2);
  }
  else if (dis.lRowElement(0)->Material()->MaterialType() == m_lung_ogden)
  {
    p_(0) = sqrt(DRT::Problem::Instance()->Material(0).m.lung_ogden->c);
    p_(1) = sqrt(DRT::Problem::Instance()->Material(0).m.lung_ogden->k1);
    p_(2) = sqrt(DRT::Problem::Instance()->Material(0).m.lung_ogden->k2);
  }
  else dserror("The inverse analysis is only implemented for the LungOgden and the LungPenalty material");

  numb_run_=0;
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
        if      (predictor==1) {
          ConstantPredictor();

        }
        else if (predictor==2) {
          ConsistentPredictor();
        }

        // gets the displacments per timestep
        if (problem_type_==0)
          get_calculated_curve(dis_,  i);
        else
          get_calculated_curve(fint_,  i);

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

  if (discret_.lRowElement(0)->Material()->MaterialType() == m_lung_penalty)
  {
    DRT::Problem::Instance()->Material(0).m.lung_penalty->c      = (p_(0)*p_(0));
    DRT::Problem::Instance()->Material(0).m.lung_penalty->k1     = (p_(1)*p_(1));
    DRT::Problem::Instance()->Material(0).m.lung_penalty->k2     = (p_(2)*p_(2));
  }
  else if (discret_.lRowElement(0)->Material()->MaterialType() == m_lung_ogden)
  {
    DRT::Problem::Instance()->Material(0).m.lung_ogden->c      = (p_(0)*p_(0));
    DRT::Problem::Instance()->Material(0).m.lung_ogden->k1     = (p_(1)*p_(1));
    DRT::Problem::Instance()->Material(0).m.lung_ogden->k2     = (p_(2)*p_(2));
    //DRT::Problem::Instance()->Material(0).m.lung_ogden->youngs = (p_(0)*p_(0));
  }
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
    cout << measured_value_(i) << "\t" << final_value_(i) << endl;
    residual_value_(i) = int(final_value_(i)) - measured_value_(i);
  }

  // calculating the errors
  if (abs(error_-error_o_)<0.1)
    error_o_=0.0;
  else
    error_o_ = error_;
  error_ = residual_value_.Norm1();

  // store the error and mu_
  storage_residual_value_.push_back(error_);
  storage_mu_.push_back(mu_);

  // print the storage
  cout << "error" << "\t" << "mu_" << "\t" << "c" << "\t" << "k1"<< "\t" << "k2" << endl;
  for (unsigned int i=0; i < storage_residual_value_.size(); i++)
  {
    cout << storage_residual_value_[i] << "\t" << storage_mu_[i] << "\t" << p_0_storage_[i]<< "\t" << p_1_storage_[i]<< "\t" << p_2_storage_[i] << endl;
  }

  //calculating J(p)
  for (int i=0; i<mp_; i++)
  {
    for (unsigned int j=0; j<3; j++)
    {
      J(i, j) = (final_value_[i]-final_value_o_[i])/(p_[j]-p_o_[j]);
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
        storage[i][j] = (storage[i][j] + storage[i][j]* mu_);
    }
  }

  //calculating (J.T*J+mu*I).I
  LINALG::NonSymmetricInverse(storage,  3);

  //calculating (J.T*J+mu*I).I*J.T
  storage2.Multiply('N', 'T', 1,  storage, J, 0);

  //calculating (J.T*J+mu*I).I*J.T*residual_disp_
  delta_p.Multiply('N', 'N', 1,  storage2, residual_value_, 0);

  p_o_ = p_;
  for (int i=0; i<3; i++)
  {
    p_(i)   = p_(i) - delta_p(i);
  }

  if (error_o_<error_) mu_ = mu_plus_;
  else  mu_ = mu_minus_;

  final_value_o_ = final_value_;

  if (abs(error_-error_o_)<0.1)
  {
    cout << "Reset old parameters to Null!!!" << endl;
    p_=p_o_;
    p_o_.Scale(0.0);
    final_value_o_.Scale(0.0);
  }

  cout << "********************************************ENDE********************************************" << endl;
  return;
}


void Inv_analysis::get_calculated_curve(const RefCountPtr<Epetra_Vector> value, int numb)
{
  int direction = 0;
  double nodal_sum = 0;
  double nodal_sum_2 = 0;

  // check in which direction the surface neumann conditions pulls
  if      (numb==0)
    direction = 0;
  else
  {
    if (abs((*value)[bc_nodes_[0]]) > abs((*value)[bc_nodes_[1]]) && abs((*value)[bc_nodes_[0]]) > abs((*value)[bc_nodes_[2]]))
      direction = 0;
    else if (abs((*value)[bc_nodes_[1]]) > abs((*value)[bc_nodes_[0]]) && abs((*value)[bc_nodes_[1]]) > abs((*value)[bc_nodes_[2]]))
      direction = 1;
    else if (abs((*value)[bc_nodes_[2]]) > abs((*value)[bc_nodes_[0]]) && abs((*value)[bc_nodes_[2]]) > abs((*value)[bc_nodes_[1]]))
      direction = 2;
    else dserror("The displacment direction is not clear!");
  }

  // summing up the displacments at each of the surface neum nodes
  for (unsigned int i=0; i<bc_nodes_.size(); i++)
  {
    if (value->Map().MyGID(bc_nodes_[i]*3 + direction))
      nodal_sum = nodal_sum + abs((*value)[value->Map().LID(bc_nodes_[i]*3 + direction)]);
  }

  // summing up over all processors
  discret_.Comm().SumAll(&nodal_sum, &nodal_sum_2, 1);
  final_value_[numb] =  nodal_sum_2/bc_nodes_.size();
  return;
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

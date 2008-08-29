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

#include "../drt_structure/stru_resulttest.H"



using namespace std;

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
  negative_material_parameters_= false;
  error_=1.0E6;
  error_o_=1.5E6;

  mp_   = params_.get<int>   ("nstep",5);

  p_.Resize(3);
  p_o_.Resize(3);

  final_disp_.Resize(mp_);
  final_disp_o_.Resize(mp_);
  measured_disp_.Resize(mp_);
  residual_disp_.Resize(mp_);

  get_measured_disp();

  mu_ =1;                                     // start value
  mu_minus_= params.get("mu_minus", 1.0);
  mu_plus_ = params.get("mu_plus", 1.0);

  //  tolerance for the curve fitting
  tol_ = params.get("inv_ana_tol", 1.0);

  numb_run_=0;

   // get the surface neuman nodes
   discret_.GetCondition("SurfaceNeumann",surfneum_ );

}

void Inv_analysis::get_measured_disp()
{
  //measured_disp_.push_back(params.get("measured_disp0", 0.0));
  //measured_disp_.push_back(params.get("measured_disp1", 0.0));
  //measured_disp_.push_back(params.get("measured_disp2", 0.0));

  Epetra_SerialDenseVector measured_disp(50);

  measured_disp(0) = 0.0;
  measured_disp(1) = 4.0;
  measured_disp(2) = 12.0;
  measured_disp(3) = 23.0;
  measured_disp(4) = 37.0;
  measured_disp(5) = 53.0;
  measured_disp(6) = 72.0;
  measured_disp(7) = 94.0;
  measured_disp(8) = 119.0;
  measured_disp(9) = 146.0;
  measured_disp(10) = 176.0;
  measured_disp(11) = 208.0;
  measured_disp(12) = 242.0;
  measured_disp(13) = 278.0;
  measured_disp(14) = 317.0;
  measured_disp(15) = 356.0;
  measured_disp(16) = 396.0;
  measured_disp(17) = 439.0;
  measured_disp(18) = 483.0;
  measured_disp(19) = 527.0;
  measured_disp(20) = 573.0;
  measured_disp(21) = 620.0;
  measured_disp(22) = 666.0;
  measured_disp(23) = 713.0;
  measured_disp(24) = 760.0;
  measured_disp(25) = 805.0;
  measured_disp(26) = 853.0;
  measured_disp(27) = 900.0;
  measured_disp(28) = 945.0;
  measured_disp(29) = 990.0;
  measured_disp(30) = 1034.0;
  measured_disp(31) = 1077.0;
  measured_disp(32) = 1119.0;
  measured_disp(33) = 1159.0;
  measured_disp(34) = 1198.0;
  measured_disp(35) = 1235.0;
  measured_disp(36) = 1270.0;
  measured_disp(37) = 1302.0;
  measured_disp(38) = 1333.0;
  measured_disp(39) = 1362.0;
  measured_disp(40) = 1388.0;
  measured_disp(41) = 1412.0;
  measured_disp(42) = 1432.0;
  measured_disp(43) = 1451.0;
  measured_disp(44) = 1466.0;
  measured_disp(45) = 1479.0;
  measured_disp(46) = 1489.0;
  measured_disp(47) = 1496.0;
  measured_disp(48) = 1500.0;
  measured_disp(49) = 1510.0;

  if (mp_==50)
  {
    measured_disp_=measured_disp;
  }

  else if (mp_==10)
  {
    for (int i=0;i<mp_;i++)
    {
      measured_disp_(i)=measured_disp(i*5);
    }
  }

  else if (mp_==5)
  {
    for (int i=0;i<mp_;i++)
    {
      measured_disp_(i)=measured_disp(i*12);
    }
  }

  else  dserror("Inverse analysis is only implemented for 50, 10 or 5 measurment values");

}



void Inv_analysis::evaluate()
{

  cout << "********************************************************************************************" << endl;
  cout << "*******************************\t Inverse Analysis \t************************************" << endl;
  cout << "*******************************\t\t run \t\t************************************" << endl;
  cout << "*******************************\t\t  " << numb_run_ <<"\t\t************************************" << endl;
  cout << "********************************************************************************************" << endl;

  calculate_new_parameters();
  numb_run_++;
  reset_parameters();
  return;
}


void Inv_analysis::Integrate()
{
  int max_itter = 100000;

  do {
    int    step    = params_.get<int>   ("step" ,0);
    double maxtime = params_.get<double>("max time",0.0);

    get_surfneum_nodes();

    // can have values "full newton" , "modified newton" , "nonlinear cg"
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

void Inv_analysis::get_surfneum_nodes()
{
  set<int> nodeSet;
  int test = 0;

  // loop for every surfneum condition, for the tension test it should be one time only
  for(unsigned int i=0; i<surfneum_.size(); i++)
  {
    map< int, RefCountPtr<DRT::Element > >           geometryMap = surfneum_[i]->Geometry();
    map< int, RefCountPtr<DRT::Element > >::iterator iterGeo;

    // loop for each geometry should be one as well
    for(iterGeo = geometryMap.begin(); iterGeo != geometryMap.end(); iterGeo++ )
    {
      // loop for each node of the geometry and condition
      for(int inode = 0; inode < iterGeo->second.get()->NumNode(); inode++)
      {
        int nodeId = iterGeo->second.get()->Nodes()[inode]->Id();

    	// if node id is not in the array yet
        if(nodeSet.find(nodeId) == nodeSet.end())
        {
          nodeSet.insert(nodeId);

          //check if the node is already in the array
          for (i=0; i<surfneum_nodes_.size(); i++)
          {
            if (surfneum_nodes_[i]== nodeId) test++;
          }

          //if it is a new node put it in the array
          if (test == 0) surfneum_nodes_.push_back(nodeId);

          // reset test;
          test = 0;
        }
      }
    }
  }


    static int myrank;
    myrank = dis_->Comm().MyPID();
    cout << "Myrank " << myrank << endl;
    cout << "Nodes: " << endl;
    for (int i=0; i<surfneum_nodes_.size();i++) {
      cout << surfneum_nodes_[i] << " ";
    }
    cout << endl;


  return;
}


void Inv_analysis::get_disp_curve(const RefCountPtr<Epetra_Vector> disp,  int numb)
{
  vector<double> final_disp_x;
  vector<double> final_disp_y;
  vector<double> final_disp_z;
  double sum_x = 0.0;
  double sum_y = 0.0;
  double sum_z = 0.0;

  for (unsigned int i=0; i<surfneum_nodes_.size(); i++)
  {
    final_disp_x.push_back((*disp)[surfneum_nodes_[i]*3]);
    sum_x = sum_x+(*disp)[surfneum_nodes_[i]*3];

    final_disp_y.push_back((*disp)[surfneum_nodes_[i]*3+1]);
    sum_y = sum_y+(*disp)[surfneum_nodes_[i]*3+1];

    final_disp_z.push_back((*disp)[surfneum_nodes_[i]*3+2]);
    sum_z = sum_z+(*disp)[surfneum_nodes_[i]*3+2];
  }

  if 	  (abs(sum_x) > abs(sum_y) && abs(sum_x) > abs(sum_z)) 	final_disp_(numb)= sum_x/final_disp_x.size();
  else if (abs(sum_y) > abs(sum_x) && abs(sum_y) > abs(sum_z)) 	final_disp_(numb)= sum_y/final_disp_y.size();
  else 								final_disp_(numb)= sum_z/final_disp_z.size();

  return;
}

void Inv_analysis::calculate_new_parameters()
{

  // initalization of the Jacobi and storage matrix

  Epetra_SerialDenseMatrix J(mp_, 3);
  Epetra_SerialDenseMatrix storage(3, 3);
  Epetra_SerialDenseMatrix storage2(3, mp_);
  Epetra_SerialDenseVector delta_p(3);

  // getting the values for the parameter vector
  p_(0) = DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->c;
  p_(1) = DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k1;
  p_(2) = DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k2;

  p_(0) = sqrt(p_(0));
  p_(1) = sqrt(p_(1));
  p_(2) = sqrt(p_(2));


  p_0_storage_.push_back(p_[0]);
  p_1_storage_.push_back(p_[1]);
  p_2_storage_.push_back(p_[2]);

  // row map for one prosessor
  //const Epetra_Map* final_disp_row_map = final_disp_.Map();

  // calculating delta_p_

  //calculating residual displacment
  for (int i=0; i<mp_; i++)
  {
    residual_disp_(i) = final_disp_(i) - measured_disp_(i);
  }

  // calculating the errors
  error_o_ = error_;
  error_ = residual_disp_.Norm1();

  // store the error and mu_
  storage_residual_disp_.push_back(error_);
  storage_mu_.push_back(mu_);

  // print the storage
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
//  cout << "J"<< endl;
//  cout <<  J << endl;

  //calculating J.T*J
  storage.Multiply('T',  'N',  1,  J, J,  0);
  cout << "J.T*J"<< endl;
  cout <<  storage  << endl;

  //calculating J.T*J+mu*I
  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    {
      if (i==j)
      {
        storage[i][j] = storage[i][j] + mu_;
      }
    }
  }
  cout << "J.T*J+mu*I"<< endl;
  cout <<  storage  << endl;

  //calculating (J.T*J+mu*I).I
  LINALG::NonSymmetricInverse(storage,  3);
//  cout << "(J.T*J+mu*I).I"<< endl;
//  cout <<  storage  << endl;

  //calculating (J.T*J+mu*I).I*J.T
  storage2.Multiply('N', 'T', 1,  storage, J, 0);
//  cout << "(J.T*J+mu*I).I*J.T"<< endl;
//  cout <<  storage2  << endl;

  //calculating (J.T*J+mu*I).I*J.T*residual_disp_
  delta_p.Multiply('N', 'N', 1,  storage2, residual_disp_, 0);
//  cout << "(J.T*J+mu*I).I*J.T*residual_disp_"<< endl;
//  cout <<  storage2  << endl;

  p_o_ = p_;
  for (int i=0; i<3; i++)
  {
    p_(i)   = p_(i) - delta_p(i);
  }

  if (error_o_<error_) mu_ = mu_plus_;
  else  mu_ = mu_minus_;

  final_disp_o_ = final_disp_;

  cout << "Displacment Curve:" << final_disp_ << endl;

  // write back new material properties
  DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->c      = (p_(0)*p_(0));
  DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k1     = (p_(1)*p_(1));
  DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k2     = (p_(2)*p_(2));

  return;
}


#endif

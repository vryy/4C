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

/*#include "../drt_structure/stru_dyn_nln_drt.H"
#include "../drt_structure/stru_genalpha_zienxie_drt.H"
#include "../drt_structure/strugenalpha.H"
#include "../drt_contact/contactstrugenalpha.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"
#include "../drt_structure/stru_resulttest.H"
*/
//#include "../io/io_drt.H"



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

    mu_ = 1.0;
    mu_minus_ = 0.1;
    mu_plus_ = 10.0;
/*    mu_ = 0.0;
    mu_minus_ = 0.0;
    mu_plus_ = 0.0;*/
/*    mu_ = 1.0;
    mu_minus_ = 10.0;
    mu_plus_ = 0.1;*/
/*    mu_ = 1.0;
    mu_minus_ = 1.0;
    mu_plus_ = 1.0;*/
        
    negative_material_parameters_= false;
    error_=1.0E6;
    error_o_=1.5E6;        

    measured_disp_.push_back(params.get("measured_disp0", 0.0));
    measured_disp_.push_back(params.get("measured_disp1", 0.0));
    measured_disp_.push_back(params.get("measured_disp2", 0.0));

    //  tolerance for the curve fitting
    tol_=0.15*measured_disp_[2];

    p_start_.push_back(1100.0);
    p_start_.push_back(3500.0);
    p_start_.push_back(2.5);

    for (unsigned int i=0; i<3; i++) {
    	delta_p_.push_back(0.0);
    	p_o_.push_back(0.0);
    	p_.push_back(0.0);
    	p_neg_.push_back(0.0);
    	final_disp_.push_back(0.0);
    	final_disp_o_.push_back(0.0);
    	residual_disp_.push_back(0.0);
    	measured_disp_.push_back(0.0);
    }

    cout << "initialisierung: neg Parameters:     \t\t\t"<< p_neg_[0] << " \t" << p_neg_[1] << " \t" << p_neg_[2] << endl; 
  }


void Inv_analysis::evaluate()
{
  // get the final displacment at that surface neumans nodes

  cout << "********************************************************************************************" << endl;
  cout << "*******************************\t Inverse Analysis \t************************************" << endl;
  cout << "*******************************\t\t run \t\t************************************" << endl;
  cout << "*******************************\t\t  " << numb_run_ <<"\t\t************************************" << endl;
  cout << "********************************************************************************************" << endl;

  double iotime = params_.get<double>("total time",0.0);
  int iostep = params_.get<int>("step",0);

  // reset the counting values
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

//  // Inverse analyse
//
//  // get the surface neuman nodes
//  vector<DRT::Condition*> surfneum;
//  discret_.GetCondition("SurfaceNeumann",surfneum );
//  get_surfneum_nodes(surfneum);

  // calculate new material parameters
  calculate_new_parameters();
  output_.NewResultFile((numb_run_+1));
  output_.WriteMesh(iostep, iotime);

  numb_run_++;

  return;
}


void Inv_analysis::Integrate()
{
  int max_itter = 100;
  numb_run_= 0;

  do {
    int    step    = params_.get<int>   ("step" ,0);
    int    nstep   = params_.get<int>   ("nstep",5);
    double maxtime = params_.get<double>("max time",0.0);

    // Inverse analyse

    // get the surface neuman nodes
    
    vector<DRT::Condition*> surfneum;

    discret_.GetCondition("SurfaceNeumann",surfneum );
    get_surfneum_nodes(surfneum);


    // can have values "full newton" , "modified newton" , "nonlinear cg"
    string equil = params_.get<string>("equilibrium iteration","full newton");

    // can have values takes values "constant" consistent"
    string pred  = params_.get<string>("predictor","constant");
    int predictor=-1;
    if      (pred=="constant")   predictor = 1;
    else if (pred=="consistent") predictor = 2;
    else dserror("Unknown type of predictor");

    //in case a constraint is defined, use defined algorithm
    if (constrMan_->HaveConstraint())
    {
      string algo = params_.get<string>("uzawa algorithm","newtonlinuzawa");
      for (int i=step; i<nstep; ++i)
      {
        if      (predictor==1) ConstantPredictor();
        else if (predictor==2) ConsistentPredictor();
	        
        //Does predicted displacement satisfy constraint?
	double time = params_.get<double>("total time",0.0);
	//double dt   = params_.get<double>("delta time",0.01);
	// what algorithm is used?
	// - "newtonlinuzawa": 			Potential is linearized wrt displacements and Lagrange multipliers
	//					   			Linear problem is solved with Uzawa algorithm
	// - "augmentedlagrange":		Potential is linearized wrt displacements keeping Lagrange multiplier fixed
	//								Until convergence Lagrange multiplier increased by Uzawa_param*(Vol_err)
	if (algo=="newtonlinuzawa")
	{
	  FullNewtonLinearUzawa();
	}
	else if (algo=="augmentedlagrange")
	{
	  NonLinearUzawaFullNewton(predictor);
	}
	else dserror("Unknown type of algorithm to deal with constraints");
	StruGenAlpha::UpdateandOutput();
	if (time>=maxtime) break;
      }
    }
    else if (equil=="full newton")
    {
      for (int i=step; i<nstep; ++i)
      {
        if      (predictor==1) ConstantPredictor();
	else if (predictor==2) ConsistentPredictor();
        
	if (i==int(nstep/3))
	  get_final_displacment(dis_, 0);
	if (i==2*int(nstep/3))
	  get_final_displacment(dis_, 1);
	if (i==nstep-1)
	  get_final_displacment(dis_, 2);

	StruGenAlpha::FullNewton();
	StruGenAlpha::UpdateandOutput();
	double time = params_.get<double>("total time",0.0);
	if (time>=maxtime) break;
      }
    }
    else if (equil=="modified newton")
    {
      for (int i=step; i<nstep; ++i)
      {
        if      (predictor==1) ConstantPredictor();
	else if (predictor==2) ConsistentPredictor();
	StruGenAlpha::ModifiedNewton();
	StruGenAlpha::UpdateandOutput();
	double time = params_.get<double>("total time",0.0);
	if (time>=maxtime) break;
      }
    }
    else if (equil=="nonlinear cg")
    {
      for (int i=step; i<nstep; ++i)
      {
        if      (predictor==1) ConstantPredictor();
	else if (predictor==2) ConsistentPredictor();
	StruGenAlpha::NonlinearCG();
	StruGenAlpha::UpdateandOutput();
	double time = params_.get<double>("total time",0.0);
	if (time>=maxtime) break;
      }
    }
    else if (equil=="ptc")
    {
      for (int i=step; i<nstep; ++i)
      {
        if      (predictor==1) ConstantPredictor();
	else if (predictor==2) ConsistentPredictor();
	StruGenAlpha::PTC();
	StruGenAlpha::UpdateandOutput();
	double time = params_.get<double>("total time",0.0);
	if (time>=maxtime) break;
      }
    }
    else dserror("Unknown type of equilibrium iteration");

    //Barrier
    // ask drt problem.instance() nach dem epetra comm (communicator) der weiss auf welchem prozessor wir sind myrank und dann nur fuer 0 laufen lassen!
    evaluate();
    //Barrier

  } while (error_>tol_ && numb_run_<max_itter);

  return;
}


void Inv_analysis::get_surfneum_nodes(vector<DRT::Condition*> surfneumConditions) 
{
  set<int> nodeSet;
  int test = 0;

  // loop for every surfneum condition, for the tension test it should be one time only
  for(unsigned int i=0; i<surfneumConditions.size(); i++) 
  {
    map< int, RefCountPtr<DRT::Element > >  geometryMap = surfneumConditions[i]->Geometry();
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
  return;
}


void Inv_analysis::get_final_displacment(const RefCountPtr<Epetra_Vector> disp, int numb)
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

  if 	(abs(sum_x) > abs(sum_y) && abs(sum_x) > abs(sum_z)) 	final_disp_[numb]= sum_x/final_disp_x.size();
  else if (abs(sum_y) > abs(sum_x) && abs(sum_y) > abs(sum_z)) 	final_disp_[numb]= sum_y/final_disp_y.size();
  else 								final_disp_[numb]= sum_z/final_disp_z.size();

  return;
}

void Inv_analysis::calculate_new_parameters()
{
  // initalization of the Jacobi and storage matrix
  double J[3][3];
  double storage[3][3];
  double storage2[3][3];
  for (unsigned int i=0; i<3; i++)
  {
    for (unsigned int j=0; j<3; j++)
    {
      J[i][j] 		= 0.0;
      storage[i][j]	= 0.0;
      storage2[i][j]    = 0.0;
    }
  }
  double storage_det	= 0.0;

  // getting the values for the parameter vector
  p_[0] = DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->c;
  p_[1] = DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k1;
  p_[2] = DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k2;
  
  p_0_storage_.push_back(p_[0]);
  p_1_storage_.push_back(p_[1]);
  p_2_storage_.push_back(p_[2]);

  // calculating delta_p_
  


  //calculating residual displacment
  cout << endl;
  cout << "Old residual disp: " << residual_disp_[0] << "\t" << residual_disp_[1] << "\t" << residual_disp_[2] << endl;
  
  for (unsigned int i=0; i<3; i++)
  {
    residual_disp_[i] = final_disp_[i] - measured_disp_[i];
  }
  
  cout << "New residual disp: " << residual_disp_[0] << "\t" << residual_disp_[1] << "\t" << residual_disp_[2] << endl;
  cout << endl;
  
  error_o_ = error_;
  error_ = abs(residual_disp_[2])+abs(residual_disp_[1])+abs(residual_disp_[0])/3;

  storage_residual_disp_.push_back(error_);
  storage_mu_.push_back(mu_);

  cout << "All old residual displacments:" << endl;
  for (unsigned int i=0; i < storage_residual_disp_.size(); i++)
  {
    cout << storage_residual_disp_[i] << "\t" << storage_mu_[i] << "\t" << p_0_storage_[i]<< "\t" << p_1_storage_[i]<< "\t" << p_2_storage_[i] << endl;
  }
  
  cout << endl;
  cout << endl;
  cout << "final_disp_:        \t\t\t" << final_disp_[0] << " \t" << final_disp_[1]  << " \t" << final_disp_[2]<< endl;
  cout << "final_disp_o_:      \t\t\t" << final_disp_o_[0] << " \t" << final_disp_o_[1]  << " \t" << final_disp_o_[2]<< endl; 
  cout << endl;
  cout << "current Parameters: \t\t\t"<< p_[0] << " \t" << p_[1] << " \t" << p_[2] << endl;
  cout << "old Parameters:     \t\t\t"<< p_o_[0] << " \t" << p_o_[1] << " \t" << p_o_[2] << endl;
  cout << endl;
  cout << endl;
    
  //calculating J(p)	
  for (unsigned int i=0; i<3; i++) 
  {
    for (unsigned int j=0; j<3; j++)
    {
      J[i][j] = (final_disp_[i]-final_disp_o_[i])/(p_[j]-p_o_[j]);
    }
  }
  
  cout << "Neg. Parameter:" << p_neg_[2] << endl;

  cout << "neg Parameters:     \t\t\t"<< p_neg_[0] << " \t" << p_neg_[1] << " \t" << p_neg_[2] << endl; 
  if (abs(p_neg_[2])>1.0E-6) 
  {
    for (unsigned int i=0; i<3; i++) 
    {
      J[i][3] = (final_disp_[i]-final_disp_o_[i])/(p_[2]-p_neg_[2]);
    }
  }

  cout << "J: \t\t\t\t\t" << J[0][0] << " \t" << J[0][1] << " \t" << J[0][2] <<  endl;
  cout << "   \t\t\t\t\t" << J[1][0] << " \t" << J[1][1] << " \t" << J[1][2] <<  endl;
  cout << "   \t\t\t\t\t" << J[2][0] << " \t" << J[2][1] << " \t" << J[2][2] <<  endl;
  cout << endl;

  //calculating (J.T*J+mu*I)
  for (unsigned int a=0; a<3; a++)
  {
    for (unsigned int i=0; i<3; i++)
    {		
      for (unsigned int j=0; j<3; j++)
      {
        storage [a][i] = storage[a][i] + J[j][a]*J[j][i];
      }
      if (i==a) storage[a][i]=storage[a][i]+mu_;
    }
  }

  cout << "(J.T*J+mu*I): \t\t\t\t" << storage[0][0] << " \t" << storage[0][1] << " \t" << storage[0][2] <<  endl;
  cout << "              \t\t\t\t" << storage[1][0] << " \t" << storage[1][1] << " \t" << storage[1][2] <<  endl;
  cout << "              \t\t\t\t" << storage[2][0] << " \t" << storage[2][1] << " \t" << storage[2][2] <<  endl;
  cout << endl;

  // calculating (J.T*J+mu*I).I
  storage_det = storage[0][0]*storage[1][1]*storage[2][2] + storage[0][1]*storage[1][2]*storage[2][0] + storage[0][2]*storage[1][0]*storage[2][1] - storage[0][0]*storage[1][2]*storage[2][1] - storage[0][1]*storage[1][0]*storage[2][2] - storage[0][2]*storage[1][1]*storage[2][0];
  
  storage2[0][0]= 1/storage_det * (storage[1][1]*storage[2][2]-storage[1][2]*storage[2][1]);
  storage2[0][1]= 1/storage_det * (storage[0][2]*storage[2][1]-storage[0][1]*storage[2][2]);
  storage2[0][2]= 1/storage_det * (storage[0][1]*storage[1][2]-storage[0][2]*storage[1][1]);
  
  storage2[1][0]= 1/storage_det * (storage[1][2]*storage[2][0]-storage[1][0]*storage[2][2]);
  storage2[1][1]= 1/storage_det * (storage[0][0]*storage[2][2]-storage[0][2]*storage[2][0]);
  storage2[1][2]= 1/storage_det * (storage[0][2]*storage[1][0]-storage[0][0]*storage[1][2]);
  
  storage2[2][0]= 1/storage_det * (storage[1][0]*storage[2][1]-storage[1][1]*storage[2][0]);
  storage2[2][1]= 1/storage_det * (storage[0][1]*storage[2][0]-storage[0][0]*storage[2][1]);
  storage2[2][2]= 1/storage_det * (storage[0][0]*storage[1][1]-storage[0][1]*storage[1][0]);

  for (unsigned int i=0; i<3; i++)
  {
    for (unsigned int j=0; j<3; j++)
    {
      storage[i][j]=storage2[i][j];
    }
  }
  
  cout << "(J.T*J+mu*I).I: \t\t\t" << storage[0][0] << " \t" << storage[0][1] << " \t" << storage[0][2] <<  endl;
  cout << "                \t\t\t" << storage[1][0] << " \t" << storage[1][1] << " \t" << storage[1][2] <<  endl;
  cout << "                \t\t\t" << storage[2][0] << " " << storage[2][1] << " \t" << storage[2][2] <<  endl;
  cout << endl;
  
  // calculating (J.T*J+mu*I).I*J.T
  for (unsigned int a=0; a<3; a++)
  {
    for (unsigned int i=0; i<3; i++)
    {
      storage2[a][i] = 0.0;
      for (unsigned int j=0; j<3; j++)
      {
        storage2[a][i] = storage2[a][i] + storage[a][j]*J[i][j];
      }
    }
  }
  
  for (unsigned int i=0; i<3; i++)
  {
    for (unsigned int j=0; j<3; j++)
    {
      storage[i][j]=storage2[i][j];
    }
  }

  cout << "(J.T*J+mu*I).I*J.T: \t\t\t" << storage[0][0] << " \t" << storage[0][1] << " \t" << storage[0][2] <<  endl;
  cout << "                    \t\t\t" << storage[1][0] << " \t" << storage[1][1] << " \t" << storage[1][2] <<  endl;
  cout << "                    \t\t\t" << storage[2][0] << " \t" << storage[2][1] << " \t" << storage[2][2] <<  endl;
  cout << endl;
  cout << "residual_disp_:     \t\t\t" << residual_disp_[0] << " \t" << residual_disp_[1] << " \t" << residual_disp_[2] <<  endl;


  // calculating (J.T*J+mu*I).I*J.T*residual_disp_
  for (unsigned int i=0; i<3; i++) 
  {
    for (unsigned int j=0; j<3; j++) 
    {
      delta_p_[i] = delta_p_[i] + storage[i][j]*residual_disp_[j];
    }
  } 

  for (unsigned int i=0; i<3; i++) 
  {
    p_o_[i] = p_[i];
    p_[i]   = p_[i] - delta_p_[i];

  }
 
  // print  out
  cout << "old Parameters:     \t\t\t"<< p_o_[0] << " \t" << p_o_[1] << " \t" << p_o_[2] << endl;
  cout << "current Parameters: \t\t\t"<< p_[0] << " \t" << p_[1] << " \t" << p_[2] << endl;
  cout << "neg Parameters:     \t\t\t"<< p_neg_[0] << " \t" << p_neg_[1] << " \t" << p_neg_[2] << endl; 
  cout << "delta_p_:           \t\t\t" << delta_p_[0] << " \t" << delta_p_[1] << " \t" << delta_p_[2] << endl;
  cout << "residual_disp_:     \t\t\t" << residual_disp_[0] << " \t" << residual_disp_[1]  << " \t" << residual_disp_[2]<< endl;
  cout << "final_disp_:        \t\t\t" << final_disp_[0] << " \t" << final_disp_[1]  << " \t" << final_disp_[2]<< endl;
  cout << "final_disp_o_:      \t\t\t" << final_disp_o_[0] << " \t" << final_disp_o_[1]  << " \t" << final_disp_o_[2]<< endl;
  cout << "measured_disp_:     \t\t\t" << measured_disp_[0] << " \t" << measured_disp_[1]  << " \t" << measured_disp_[2]<< endl;
  cout << "mu                  \t\t\t" << mu_ << endl;
  cout << "Toleranz:           \t\t\t" << tol_ << endl;
  cout << "Error:              \t\t\t" << error_ << endl;
  cout << "Error_old_:         \t\t\t" << error_o_ << endl;
    
  
  if (error_o_<error_) mu_ = mu_plus_;
  else  mu_ = mu_minus_;

  final_disp_o_ = final_disp_;

  // Material Parameter muessen positiv sein
  for (unsigned int i=0; i<3; i++) 
  {
    if (p_[i]<0) 
    {
      p_neg_[i]=p_[i];
      p_[i] = abs(p_o_[i])/2;
    }
    else p_neg_[i]=0.0;
  }

  // write back new material properties
  DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->c      = p_[0];
  DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k1     = p_[1];
  DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k2     = p_[2];

  return;
}


#endif

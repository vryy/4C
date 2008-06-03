#ifdef CCADISCRET

#include "inv_analysis.H"
#include <ctime>
#include <cstdlib>
#include <iostream>
#include "Epetra_SerialDenseMatrix.h"
#include "../drt_lib/global_inp_control2.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
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
	// test
    mu_ = 1;
    mu_minus_ = 0.1;
    mu_plus_ = 10;
    measured_disp_ = 0.02;
    tol_=0.15*measured_disp_;
    for (unsigned int i=0; i<3; i++) {
    	delta_p_.push_back(0.0);
    	p_o_.push_back(0.0);
    	p_.push_back(0.0);
    }

  }


void Inv_analysis::evaluate()
  {

	  // reset the counting values
	     params_.set<int>("step", 0);
	     params_.set<double>("total time", 0.0);
	     params_.set<int> ("numstep",0);
             output_.RevertResultfileChanged();
	 //    if (numb_run_ != 1) output_.SetRestartValues();

	cout << "INV ANALYSIS" << endl;
 // Inverse analyse

  // get the surface neuman nodes
     vector<DRT::Condition*> surfneum;
     discret_.GetCondition("SurfaceNeumann",surfneum );
     get_surfneum_nodes(surfneum);

  // get the final displacment at that surface neuman nodes
     get_final_displacment(dis_);

  // calculate new material parameters
     calculate_new_parameters();

  return;

  }

void Inv_analysis::Integrate()
{
  int max_itter = 100000;
  numb_run_= 0;
  do {

    StruGenAlpha::Integrate();
    cout << "********************************************************************************************" << endl;
    cout << "*******************************\t Inverse Analysis \t************************************" << endl;        
    cout << "*******************************\t\t run \t\t************************************" << endl;      
    cout << "*******************************\t\t  " << numb_run_ <<"\t\t************************************" << endl;   
    cout << "********************************************************************************************" << endl;
    //Barrier
    // ask drt problem.instance() nach dem epetra comm (communicator) der weiss auf welchem prozessor wir sind myrank und dann nur fuer 0 laufen lassen!
    evaluate();
    //Barrier
    numb_run_++;
  } while ((abs(residual_disp_)>tol_) && numb_run_<max_itter);

  return;
}




void Inv_analysis::get_surfneum_nodes(vector<DRT::Condition*> surfneumConditions) {

    set<int> nodeSet;
    int test = 0;

    // loop for every surfneum condition, for the tension test it should be one time only
    for(unsigned int i=0; i<surfneumConditions.size(); i++) {
       map< int, RefCountPtr<DRT::Element > >  geometryMap = surfneumConditions[i]->Geometry();
       map< int, RefCountPtr<DRT::Element > >::iterator iterGeo;
       // loop for each geometry should be one as well
       for(iterGeo = geometryMap.begin(); iterGeo != geometryMap.end(); iterGeo++ ) {
          // loop for each node of the geometry and condition
          for(int inode = 0; inode < iterGeo->second.get()->NumNode(); inode++) {
             int nodeId = iterGeo->second.get()->Nodes()[inode]->Id();
             // if node id is not in the array yet
             if(nodeSet.find(nodeId) == nodeSet.end()) {
                nodeSet.insert(nodeId);
                //check if the node is already in the array
                for (i=0; i<surfneum_nodes_.size(); i++) {
                   if (surfneum_nodes_[i]== nodeId)
                      test++;
                }
                //if it is a new node put it in the array
                if (test == 0)
                   surfneum_nodes_.push_back(nodeId);
                // reset test;
                test = 0;
             }
          }
       }
    }
    return;
}




void Inv_analysis::get_final_displacment(const RefCountPtr<Epetra_Vector> disp)
{
  vector<double> final_disp_x;
  vector<double> final_disp_y;
  vector<double> final_disp_z;
  double sum_x = 0.0;
  double sum_y = 0.0;
  double sum_z = 0.0;

  for (unsigned int i=0; i<surfneum_nodes_.size(); i++) {
    final_disp_x.push_back((*disp)[surfneum_nodes_[i]*3]);
    sum_x = sum_x+(*disp)[surfneum_nodes_[i]*3];

    final_disp_y.push_back((*disp)[surfneum_nodes_[i]*3+1]);
    sum_y = sum_y+(*disp)[surfneum_nodes_[i]*3+1];

    final_disp_z.push_back((*disp)[surfneum_nodes_[i]*3+2]);
    sum_z = sum_z+(*disp)[surfneum_nodes_[i]*3+2];
  }

  if (sum_x > sum_y && sum_x > sum_z)
	  final_disp_= sum_x/final_disp_x.size();
  else if (sum_y > sum_x && sum_y > sum_z)
	  final_disp_= sum_y/final_disp_y.size();
  else
	  final_disp_= sum_z/final_disp_z.size();

  return;
}

void Inv_analysis::calculate_new_parameters(){

  double J[3];											// jacobian
  for (unsigned int i=0; i<3; i++) {
	  J[i]=0;
	  //delta_p_[i]=0;
  }

  double storage[3][3];									// storage
  for (unsigned int i=0; i<3; i++) {
    for (unsigned int j=0; j<3; j++) {
    	storage[i][j] = 0;
    }
  }
  double storage_det = 0;

  // getting the values for the parameter vector
  p_[0] = DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->c;
  p_[1] = DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k1;
  p_[2] = DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k2;

  // calculating delta_p_

  	//calculating residual displa
	  residual_disp_ = final_disp_ - measured_disp_;

	  storage_residual_disp_.push_back(residual_disp_);  
	  cout << "All old residual displacments:" << endl;
	  for (unsigned int i=0; i < storage_residual_disp_.size(); i++)
	  {
		  cout << storage_residual_disp_[i] << endl;  
	  }
	  cout << "---" << endl;	 
	  
	  //calculating J(p)
	  for (unsigned int i=0; i<3; i++) {
		  J[i] = (final_disp_ - final_disp_o_)/(p_[i]-p_o_[i]);
	  }
	  cout << "J: \t\t\t" << J[0] << " " << J[1] << " " << J[2] <<  endl;
	  
	  //calculating (J.T*J+mu*I)
	  for (unsigned int i=0; i<3; i++) {
		  for (unsigned int j=0; j<3; j++) {
			  if (i==j) storage[i][j] = J[i]*J[j] + mu_;
		      else storage[i][j] = J[i]*J[j];
		  }
	  }
	  
	  // calculating (J.T*J+mu*I).I
	  storage_det = storage[0][0]*storage[1][1]*storage[2][2] + storage[0][1]*storage[1][2]*storage[2][0] + storage[0][2]*storage[1][0]*storage[2][1] - storage[0][0]*storage[1][2]*storage[2][1] - storage[0][1]*storage[1][0]*storage[2][2] - storage[0][2]*storage[1][1]*storage[2][0];

	  storage[0][0]= 1/storage_det * (storage[1][1]*storage[2][2]-storage[1][2]*storage[2][1]);
	  storage[0][1]= 1/storage_det * (storage[0][2]*storage[2][1]-storage[0][1]*storage[2][2]);
	  storage[0][2]= 1/storage_det * (storage[0][1]*storage[1][2]-storage[0][2]*storage[1][1]);
	  storage[1][0]= 1/storage_det * (storage[1][2]*storage[2][0]-storage[1][0]*storage[2][2]);
	  storage[1][1]= 1/storage_det * (storage[0][0]*storage[2][2]-storage[0][2]*storage[2][0]);
	  storage[1][2]= 1/storage_det * (storage[0][2]*storage[1][0]-storage[0][0]*storage[1][2]);
	  storage[2][0]= 1/storage_det * (storage[1][0]*storage[2][1]-storage[1][1]*storage[2][0]);
	  storage[2][1]= 1/storage_det * (storage[0][1]*storage[2][0]-storage[0][0]*storage[2][1]);
	  storage[2][2]= 1/storage_det * (storage[0][0]*storage[1][1]-storage[0][1]*storage[1][0]);
	  
	  
	  cout << "----------------------------------------------------------------------------------------" << endl;
	  cout << "delta_p_: \t\t" << delta_p_[0] << " " << delta_p_[1] << " " << delta_p_[2] << endl;	  

	  // calculating (J.T*J+mu*I).I*J.T*residual_disp_
	  for (unsigned int i=0; i<3; i++) {
		  for (unsigned int j=0; j<3; j++) {
			  delta_p_[i] = delta_p_[i]+storage[i][j]*J[j];
			//  delta_p_[i] = storage[i][j]*J[j];
		  }
	  delta_p_[i] = delta_p_[i] * residual_disp_;
	  }
	  
	    cout << "delta_p_: \t\t" << delta_p_[0] << " " << delta_p_[1] << " " << delta_p_[2] << endl;
		  cout << "----------------------------------------------------------------------------------------" << endl;

	  // print  out
	cout << "_________"<< endl;
    cout << "old Parameters: \t\t\t"<< p_o_[0] << " " << p_o_[1] << " " << p_o_[2] << endl;
    cout << "current Parameters:   \t\t\t"<< p_[0] << " " << p_[1] << " " << p_[2] << endl;
    cout << "delta_p_: \t\t" << delta_p_[0] << " " << delta_p_[1] << " " << delta_p_[2] << endl;
    cout << "residual_disp_: \t" << residual_disp_ << endl;    
    cout << "final_disp_: \t\t"<< final_disp_ << endl;
    cout << "final_disp_o_: \t\t" << final_disp_o_ << endl;
    cout << "measured_disp_: \t" << measured_disp_ << endl;
    cout << "mu \t\t\t" << mu_ << endl;

  for (unsigned int i=0; i<3; i++) {
	  p_o_[i]          = p_[i];
	  p_[i]            = p_[i]- delta_p_[i];
	  }


  if (abs(residual_disp_) < abs(measured_disp_ - final_disp_o_))
  	mu_ = mu_plus_;
  else
	mu_ = mu_minus_;

  final_disp_o_ = final_disp_;

  // write back new material properties
  DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->c      = p_[0];
  DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k1     = p_[1];
  DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k2     = p_[2];

  
  cout << "new Parameters:   \t\t\t"<< p_[0] << " " << p_[1] << " " << p_[2] << endl;
  cout << "---------------------ende----------------------------------------------------" << endl;
  cout << "Toleranz: "<< tol_ << endl;
  return;
}


#endif

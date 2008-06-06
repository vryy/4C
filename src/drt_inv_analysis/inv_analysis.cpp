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
    negative_material_parameters_= false;

    measured_disp_.push_back(params.get("measured_disp0", 0.0));
    measured_disp_.push_back(params.get("measured_disp1", 0.0));
    measured_disp_.push_back(params.get("measured_disp2", 0.0));

    //  tolerance for the curve fitting
    tol_=0.02*measured_disp_[2];
    
    p_start_.push_back(1000.0);
    p_start_.push_back(13500.0);
    p_start_.push_back(76.5);
    
    for (unsigned int i=0; i<3; i++) {
    	delta_p_.push_back(0.0);
    	p_o_.push_back(0.0);
    	p_.push_back(0.0);
    	final_disp_.push_back(0.0);
    	final_disp_o_.push_back(0.0);
    	residual_disp_.push_back(0.0);
    	measured_disp_.push_back(0.0);   	    	
    }
    
    cout << "INITIALISIERUNG: final_disp_: \t" << final_disp_[0] << " " << final_disp_[1]  << " " << final_disp_[2]<< endl;  
    
  }


void Inv_analysis::evaluate()
{
	


    // get the final displacment at that surface neumans nodes
    


    
    return;
}

void Inv_analysis::Integrate()
{
  int max_itter = 100000;
  numb_run_= 0;
  do {
	  
	  // reset the counting values
	  params_.set<int>("step", 0);
	  params_.set<double>("total time", 0.0);
	  params_.set<int> ("numstep",0);
	  
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
	  
//	  StruGenAlpha::Integrate();
    
    
	  cout << "********************************************************************************************" << endl;
	  cout << "*******************************\t Inverse Analysis \t************************************" << endl;        
	  cout << "*******************************\t\t run \t\t************************************" << endl;      
	  cout << "*******************************\t\t  " << numb_run_ <<"\t\t************************************" << endl;   
	  cout << "********************************************************************************************" << endl;
	  
	  // calculate new material parameters
	  calculate_new_parameters();
	  
	  
	  //Barrier
	  // ask drt problem.instance() nach dem epetra comm (communicator) der weiss auf welchem prozessor wir sind myrank und dann nur fuer 0 laufen lassen!
	  evaluate();
	  //Barrier
	  
	  //close HDF File
	  //output_.        .Close();

	  output_.NewResultFile((numb_run_+1));

    double time = params_.get<double>("total time",0.0);
    step = params_.get<int>("step",0);

    
	  output_.WriteMesh(step, time);
	  
	  
	  
	  numb_run_++;
  } while ((abs(residual_disp_[2])>tol_) && numb_run_<max_itter);
  


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
  
    for (unsigned int i=0; i<surfneum_nodes_.size(); i++) {
    	
    	final_disp_x.push_back((*disp)[surfneum_nodes_[i]*3]);
        sum_x = sum_x+(*disp)[surfneum_nodes_[i]*3];
        
        final_disp_y.push_back((*disp)[surfneum_nodes_[i]*3+1]);
        sum_y = sum_y+(*disp)[surfneum_nodes_[i]*3+1];

        final_disp_z.push_back((*disp)[surfneum_nodes_[i]*3+2]);
        sum_z = sum_z+(*disp)[surfneum_nodes_[i]*3+2];
    }
       
    if 		(abs(sum_x) > abs(sum_y) && abs(sum_x) > abs(sum_z)) 	final_disp_[numb]= sum_x/final_disp_x.size();
    else if (abs(sum_y) > abs(sum_x) && abs(sum_y) > abs(sum_z)) 	final_disp_[numb]= sum_y/final_disp_y.size();
    else 															final_disp_[numb]= sum_z/final_disp_z.size();
    

  return;
}

void Inv_analysis::calculate_new_parameters(){

  // initalization of the Jacobi and storage matrix
  double J[3][3];
  double storage[3][3];    
  for (unsigned int i=0; i<3; i++){
	  for (unsigned int j=0; j<3; j++){
		  J[i][j] 		= 0.0;
		  storage[i][j]	= 0.0;
	  }
  } 
  double storage_det	= 0.0;

  // getting the values for the parameter vector
  p_[0] = DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->c;
  p_[1] = DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k1;
  p_[2] = DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k2;

  // calculating delta_p_

  cout << "CALCULATION:   final_disp_: \t" << final_disp_[0] << " " << final_disp_[1]  << " " << final_disp_[2]<< endl;  
  cout << "CALCULATION:   final_disp_o_: \t" << final_disp_o_[0] << " " << final_disp_o_[1]  << " " << final_disp_o_[2]<< endl;  
    
  cout << "CALCULATION:   measured_disp_: \t" << measured_disp_[0] << " " << measured_disp_[1]  << " " << measured_disp_[2]<< endl;  
  cout << "CALCULATION:   residual_disp_: \t" << residual_disp_[0] << " " << residual_disp_[1]  << " " << residual_disp_[2]<< endl;  
      
  
  	//calculating residual displacment
  	for (unsigned int i=0; i<3; i++){
  		residual_disp_[i] = final_disp_[i] - measured_disp_[i];
  	}
	storage_residual_disp_.push_back(residual_disp_[2]);    			
  	
	
		cout << "All old residual displacments:" << endl;
		for (unsigned int i=0; i < storage_residual_disp_.size(); i++)
		{
			cout << storage_residual_disp_[i] << endl;  
		}
		cout << "---" << endl;
		if (negative_material_parameters_== true)
		{
			cout << "Negative Material Parameters" << endl;
		}
	
    cout << "CALCULATION:   final_disp_: \t" << final_disp_[0] << " " << final_disp_[1]  << " " << final_disp_[2]<< endl;  
	cout << "CALCULATION:   final_disp_o_: \t" << final_disp_o_[0] << " " << final_disp_o_[1]  << " " << final_disp_o_[2]<< endl; 
    cout << "old Parameters: \t\t\t"<< p_o_[0] << " " << p_o_[1] << " " << p_o_[2] << endl;
    cout << "current Parameters:   \t\t\t"<< p_[0] << " " << p_[1] << " " << p_[2] << endl;
	
	  
	//calculating J(p)
	for (unsigned int i=0; i<3; i++) {
		for (unsigned int j=0; j<3; j++){
			J[i][j] = (final_disp_[i]-final_disp_o_[i])/(p_[j]-p_o_[j]);
		}
	}

	cout << "J: \t\t\t\t\t" << J[0][0] << " " << J[0][1] << " " << J[0][2] <<  endl;
	cout << "   \t\t\t\t\t" << J[1][0] << " " << J[1][1] << " " << J[1][2] <<  endl;
	cout << "   \t\t\t\t\t" << J[2][0] << " " << J[2][1] << " " << J[2][2] <<  endl;
	cout << endl;
		  
	//calculating (J.T*J+mu*I)

/*	for (unsigned int a=0; a<3; a++){
		for (unsigned int i=0; i<3; i++){
			for (unsigned int 	j=0; j<3; j++){
				storage [a][i] = storage[a][i] + J[a][j]*J[i][j];
			}
			if (i==a) storage[a][i]=storage[a][i]+mu_;
		}
	}*/
	
	
    storage[0][0] = J[0][0]*J[0][0] + J[0][1]*J[0][1] + J[0][2]*J[0][2];
    storage[0][1] = J[0][0]*J[1][0] + J[0][1]*J[1][1] + J[0][2]*J[1][2];
    storage[0][2] = J[0][0]*J[2][0] + J[0][1]*J[2][1] + J[0][2]*J[2][2];
                                                                   
    storage[1][0] = J[1][0]*J[0][0] + J[1][1]*J[0][1] + J[1][2]*J[0][2];
    storage[1][1] = J[1][0]*J[1][0] + J[1][1]*J[1][1] + J[1][2]*J[1][2];
    storage[1][2] = J[1][0]*J[2][0] + J[1][1]*J[2][1] + J[1][2]*J[2][2];

    storage[2][0] = J[2][0]*J[0][0] + J[2][1]*J[0][1] + J[2][2]*J[0][2];
    storage[2][1] = J[2][0]*J[1][0] + J[2][1]*J[1][1] + J[2][2]*J[1][2];
    storage[2][2] = J[2][0]*J[2][0] + J[2][1]*J[2][1] + J[2][2]*J[2][2];
    
    cout << "(J.T*J): \t\t\t\t" << storage[0][0] << " " << storage[0][1] << " " << storage[0][2] <<  endl;
    cout << "   \t\t\t\t\t" << storage[1][0] << " " << storage[1][1] << " " << storage[1][2] <<  endl;
    cout << "   \t\t\t\t\t" << storage[2][0] << " " << storage[2][1] << " " << storage[2][2] <<  endl;
    cout << endl;
    	
    cout << "mu:\t\t\t\t\t" << mu_<< endl;
    cout << endl;
    
    storage[0][0] = storage[0][0] +mu_;
    storage[1][1] = storage[1][1] +mu_;       
    storage[2][2] = storage[2][2] +mu_;
        
	
	cout << "(J.T*J+mu*I): \t\t\t\t" << storage[0][0] << " " << storage[0][1] << " " << storage[0][2] <<  endl;
	cout << "   \t\t\t\t\t" << storage[1][0] << " " << storage[1][1] << " " << storage[1][2] <<  endl;
	cout << "   \t\t\t\t\t" << storage[2][0] << " " << storage[2][1] << " " << storage[2][2] <<  endl;
	cout << endl;
	  
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
	
	/*	for (unsigned int a=0; a<3; a++){
			for (unsigned int i=0; i<3; i++){
				for (unsigned int 	j=0; j<3; j++){
					storage2 [a][i] = storage2[a][i] + storage[a][j]*J[i][j];
				}
			}
		}*/
	  
	cout << "(J.T*J+mu*I).I: \t\t\t" << storage[0][0] << " " << storage[0][1] << " " << storage[0][2] <<  endl;
	cout << "   \t\t\t\t\t" << storage[1][0] << " " << storage[1][1] << " " << storage[1][2] <<  endl;
	cout << "   \t\t\t\t\t" << storage[2][0] << " " << storage[2][1] << " " << storage[2][2] <<  endl;  
	cout << "----------------------------------------------------------------------------------------" << endl;
	cout << "delta_p_: \t\t" << delta_p_[0] << " " << delta_p_[1] << " " << delta_p_[2] << endl;	  

	// calculating (J.T*J+mu*I).I*J.T
	storage[0][0] = storage[0][0]*J[0][0] + storage[0][1]*J[0][1] + storage[0][2]*J[0][2];
	storage[0][1] = storage[0][0]*J[1][0] + storage[0][1]*J[1][1] + storage[0][2]*J[1][2];
	storage[0][2] = storage[0][0]*J[2][0] + storage[0][1]*J[2][1] + storage[0][2]*J[2][2];
	                                                                   
	storage[1][0] = storage[1][0]*J[0][0] + storage[1][1]*J[0][1] + storage[1][2]*J[0][2];
	storage[1][1] = storage[1][0]*J[1][0] + storage[1][1]*J[1][1] + storage[1][2]*J[1][2];
	storage[1][2] = storage[1][0]*J[2][0] + storage[1][1]*J[2][1] + storage[1][2]*J[2][2];

	storage[0][0] = storage[2][0]*J[0][0] + storage[2][1]*J[0][1] + storage[2][2]*J[0][2];
	storage[0][1] = storage[2][0]*J[1][0] + storage[2][1]*J[1][1] + storage[2][2]*J[1][2];
	storage[0][2] = storage[2][0]*J[2][0] + storage[2][1]*J[2][1] + storage[2][2]*J[2][2];
	
	cout << "(J.T*J+mu*I).I*J.T: \t\t\t" << storage[0][0] << " " << storage[0][1] << " " << storage[0][2] <<  endl;
	cout << "   \t\t\t\t\t" << storage[1][0] << " " << storage[1][1] << " " << storage[1][2] <<  endl;
	cout << "   \t\t\t\t\t" << storage[2][0] << " " << storage[2][1] << " " << storage[2][2] <<  endl;
	
	// calculating (J.T*J+mu*I).I*J.T*residual_disp_
	for (unsigned int i=0; i<3; i++) {
		for (unsigned int j=0; j<3; j++) {
			delta_p_[i] = delta_p_[i] + storage[i][j]*residual_disp_[j];
		}
	}
	
    // print  out
	cout << "delta_p_: \t\t" << delta_p_[0] << " " << delta_p_[1] << " " << delta_p_[2] << endl;
    cout << "----------------------------------------------------------------------------------------" << endl;
    cout << "_________"<< endl;
    cout << "old Parameters: \t\t\t"<< p_o_[0] << " " << p_o_[1] << " " << p_o_[2] << endl;
    cout << "current Parameters:   \t\t\t"<< p_[0] << " " << p_[1] << " " << p_[2] << endl;
    cout << "delta_p_: \t\t" << delta_p_[0] << " " << delta_p_[1] << " " << delta_p_[2] << endl;
    cout << "residual_disp_: \t" << residual_disp_[0] << " " << residual_disp_[1]  << " " << residual_disp_[2]<< endl;    
    cout << "final_disp_: \t\t" << final_disp_[0] << " " << final_disp_[1]  << " " << final_disp_[2]<< endl;    
    cout << "final_disp_o_: \t\t" << final_disp_o_[0] << " " << final_disp_o_[1]  << " " << final_disp_o_[2]<< endl;   
    cout << "measured_disp_: \t" << measured_disp_[0] << " " << measured_disp_[1]  << " " << measured_disp_[2]<< endl;    
    cout << "mu \t\t\t" << mu_ << endl;

  for (unsigned int i=0; i<3; i++) {
	  p_o_[i]          = p_[i];
	  p_[i]            = p_[i] - delta_p_[i];
	  }


  if (abs(residual_disp_[2]) < abs(measured_disp_[2] - final_disp_o_[2]))
  	mu_ = mu_plus_;
  else
	mu_ = mu_minus_;

  final_disp_o_ = final_disp_;
  
  // Material Parameter muessen positiv sein
  
  for (unsigned int i=0; i<3; i++) {
	  if (p_[i]<0) {
		  cout << "p[i]: " << p_[i] << endl;
		  p_[i] = p_start_[i] + numb_run_/100 * p_start_[i];
		  cout << "p[i]_neu: " << p_[i] << endl;
		  		  
		  negative_material_parameters_=true;
	  }
	  }

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

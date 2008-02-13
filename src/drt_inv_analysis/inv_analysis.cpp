#ifdef CCADISCRET

#include "inv_analysis.H"
#include <ctime>
#include <cstdlib>
#include <iostream>
#include "Epetra_SerialDenseMatrix.h"
#include "../drt_lib/global_inp_control2.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"


extern struct _GENPROB     genprob;
extern struct _IO_FLAGS     ioflags;
extern struct _FILES  allfiles;

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
    
  }


void Inv_analysis::evaluate()
  {
 
 // Inverse analyse
  
  // get the surface neuman nodes
     vector<DRT::Condition*> surfneum;
     discret_.GetCondition("SurfaceNeumann",surfneum );
     get_surfneum_nodes(surfneum);
  
  // get the final displacment at that surface neuman nodes
     get_final_displacment(dis_);
  
  // get measured displacment
     if (numb_run_ == 0)
       get_measured_displacment (0.1);
     
  // compare the measured and the calculated displacment   
     if (compare_displacment ())
       return;
     
  // get old_parameters
     get_old_parameters();
        
  // calculate new material parameters
     calculate_new_parameters();
     

  // reset the counting values
     params_.set<int>("step", 0);
     params_.set<double>("total time", 0.0);       
     params_.set<int> ("numstep",0);
     residual_disp_.clear();
     final_disp_.clear();
     
  
  return;

  }

void Inv_analysis::Integrate()
{
  int max_itter = 10;
  numb_run_= 0;
  do {

    StruGenAlpha::Integrate();
    evaluate();
    numb_run_++;  
  } while (!compare_displacment () && numb_run_<max_itter);

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
};




void Inv_analysis::get_final_displacment(const RefCountPtr<Epetra_Vector> disp)
{
  vector<double> final_disp_x;
  vector<double> final_disp_y;
  vector<double> final_disp_z;
  double sum_x;
  double sum_y;
  double sum_z;

  for (unsigned int i=0; i<surfneum_nodes_.size(); i++) {
    final_disp_x.push_back((*disp)[surfneum_nodes_[i]*3]);
    sum_x = sum_x+(*disp)[surfneum_nodes_[i]*3];
      
    final_disp_y.push_back((*disp)[surfneum_nodes_[i]*3+1]);
    sum_y = sum_y+(*disp)[surfneum_nodes_[i]*3+1];
      
    final_disp_z.push_back((*disp)[surfneum_nodes_[i]*3+2]);
    sum_z = sum_z+(*disp)[surfneum_nodes_[i]*3+2];
  }
   
  final_disp_.push_back(sum_x/final_disp_x.size());
  final_disp_.push_back(sum_y/final_disp_y.size());
  final_disp_.push_back(sum_z/final_disp_z.size());
     
  return;
}

void Inv_analysis::get_measured_displacment (const double measured_displacment_value) {
  
    if (final_disp_[0] > final_disp_[1] && final_disp_[0] > final_disp_[2])
      measured_disp_.push_back(measured_displacment_value);
      measured_disp_.push_back(0);
      measured_disp_.push_back(0);
    if (final_disp_[1] > final_disp_[0] && final_disp_[1] > final_disp_[2])
      measured_disp_.push_back(0);
      measured_disp_.push_back(measured_displacment_value);
      measured_disp_.push_back(0);
    if (final_disp_[2] > final_disp_[0] && final_disp_[2] > final_disp_[0])
      measured_disp_.push_back(0);
      measured_disp_.push_back(0);
      measured_disp_.push_back(measured_displacment_value);
  
  return;
}

bool Inv_analysis::compare_displacment (){
  //L2 norm for the distance between the two displacment vectors
  double comp = 0.0;
  double tolerance = 0.0001;
  
  if (numb_run_ != 0) {
    residual_disp_o_ = residual_disp_;  
  }
  
  for (unsigned int i =0; i<final_disp_.size();i++) {
   residual_disp_.push_back(abs(measured_disp_[i]-final_disp_[i])); 
   comp = comp + (residual_disp_[i]*residual_disp_[i]); 
  }
  
  comp = sqrt(comp)/final_disp_.size();

  return comp < tolerance;
}

void Inv_analysis::get_old_parameters(){
 
  // we expect just one material definition
  // mat_ = DRT::Problem::Instance()->Material(0);
  if (mat_.mattyp!=m_hyper_polyconvex) dserror("The Inv_Analysis is only for hyperpolyconvex material!");
  
  parameters_.clear();
  // built vector
  parameters_.push_back(mat_.m.hyper_polyconvex->c);
  parameters_.push_back(mat_.m.hyper_polyconvex->k1);
  parameters_.push_back(mat_.m.hyper_polyconvex->k2);

  count_++;
}

void Inv_analysis::calculate_new_parameters(){
  
  
  
  vector<double> s_k;                                   // search direction
    for (unsigned int i=0; i<3; i++) {
      s_k.push_back(0);
    }
    
  double J[3][3];                                       // jacobian
  double JTJ [3][3];
    
  double phi = 0;                                       // aim function
  
  vector<double> delta_phi;                             // gradient of the aim function
    for (unsigned int i=0; i<3; i++) {
      delta_phi.push_back(0);
    }
    
  double I[3][3];                                       // identity matrix
    
  double storage[3][3];                                 // storage matrix  
    
  double storage_det = 0;                               // storage determinante for invesion of the matrix
  int alpha = 10;                                        // step width   
  double mu = 0.1;                                        // gibt an ob mehr gradientenverfahren oder mehr gauss newton verfahren verwendet wird                     
  
  // set values for the first run!
  if (numb_run_ == 0) {
    p_o_.push_back(DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->c  + 0.1);
    p_o_.push_back(DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k1 + 20);
    p_o_.push_back(DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k2 + 3);
    p_.push_back(DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->c);
    p_.push_back(DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k1);
    p_.push_back(DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k2);
    residual_disp_o_.push_back(0.01);
    residual_disp_o_.push_back(0.012);
    residual_disp_o_.push_back(0.6);
    }
  
  // getting the weights right
   if ((residual_disp_[0]*residual_disp_[0]+residual_disp_[1]*residual_disp_[1]+residual_disp_[2]*residual_disp_[2])>(residual_disp_o_[0]*
       residual_disp_o_[0]+residual_disp_o_[1]*residual_disp_o_[1]+residual_disp_o_[2]*residual_disp_o_[2]))
     mu = 0.1;
   else
     mu =10;
 

  // getting the values for the parameter vector
  p_[0] = DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->c;
  p_[1] = DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k1;
  p_[2] = DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k2;
  
 // calculating phi
  for (int i=0; i<3; i++) {
    phi = phi + (residual_disp_[i]/2);
  }

  //calculating J(p)

  cout << "residual_disp_o_: " << endl;
  cout << residual_disp_o_[0] << " "<< residual_disp_o_[1] << " "<< residual_disp_o_[2] << endl;
  cout << endl;
  cout << "residual_disp_:   " << endl;
  cout << residual_disp_[0] << " "<< residual_disp_[1] << " "<< residual_disp_[2] << endl;
  cout << endl;
  
  cout << "p_o_:             "<< endl;
  cout << p_o_[0] << " " << p_o_[1] << " " << p_o_[2] << endl;
  cout << endl;
  
  cout << "p_:               "<< endl;
  cout << p_[0] << " " << p_[1] << " " << p_[2] << endl;
  cout << endl;
  
  for (unsigned int i=0; i<3; i++) {
    for (unsigned int j=0; j<3; j++) {
      J[i][j] = (residual_disp_[i] - residual_disp_o_[i])/(p_[j]-p_o_[j]);
    }
  }
  
  cout << "J:" << endl;
  cout << J[0][0] << " " << J[0][1] << " " << J[0][2] <<  endl;
  cout << J[1][0] << " " << J[1][1] << " " << J[1][2] <<  endl;
  cout << J[2][0] << " " << J[2][1] << " " << J[2][2] <<  endl;
  cout  << endl;
 
  //calculating delta_phi(p)
  for (unsigned int i=0; i<3; i++) {
    for (unsigned int j=0; j<3; j++) {
      delta_phi[i] = J[j][i]*residual_disp_[i];
    }
  }
 
  // calculating s_k
  for (unsigned int i=0; i<3; i++) {
    for (unsigned int j=0; j<3; j++) {
      JTJ[i][j] = 0;
     }
  }
  
  for (unsigned int i=0; i<3; i++) {
    for (unsigned int j=0; j<3; j++) {
      for (unsigned int k=0; k<3; k++) {
        JTJ[i][j] = JTJ[i][j] + J[k][i]*J[k][j];
      }
    }
  }
  
  cout << "JTJ:" << endl;
  cout << JTJ[0][0] << " " << JTJ[0][1] << " " << JTJ[0][2] <<  endl;
  cout << JTJ[1][0] << " " << JTJ[1][1] << " " << JTJ[1][2] <<  endl;
  cout << JTJ[2][0] << " " << JTJ[2][1] << " " << JTJ[2][2] <<  endl;  
  cout  << endl;
  
  for (unsigned int i=0; i<3; i++) {
    for (unsigned int j=0; j<3; j++) {
      if (i==j)
        I[i][j]=1;
      else
        I[i][j]=0;
    }
  }
  
  for (unsigned int i=0; i<3; i++) {
    for (unsigned int j=0; j<3; j++) {
      storage[i][j]=JTJ[i][j] + mu * I[i][j];
    }
  }
  
 
    // storage determinante
    storage_det = storage[0][0]*storage[1][1]*storage[2][2] + storage[0][1]*storage[1][2]*storage[2][0] + storage[0][2]*storage[1][0]*storage[2][1] - storage[0][0]*storage[1][2]*storage[2][1] - storage[0][1]*storage[1][0]*storage[2][2] - storage[0][2]*storage[1][1]*storage[2][0];
   
    
    // inverse matrix
  
    storage[0][0]= 1/storage_det * (storage[1][1]*storage[2][2]-storage[1][2]*storage[2][1]);
    storage[0][1]= 1/storage_det * (storage[0][2]*storage[2][1]-storage[0][1]*storage[2][2]);
    storage[0][2]= 1/storage_det * (storage[0][1]*storage[1][2]-storage[0][2]*storage[1][1]);
    storage[1][0]= 1/storage_det * (storage[1][2]*storage[2][0]-storage[1][0]*storage[2][2]);
    storage[1][1]= 1/storage_det * (storage[0][0]*storage[2][2]-storage[0][2]*storage[2][0]);
    storage[1][2]= 1/storage_det * (storage[0][2]*storage[1][0]-storage[0][0]*storage[1][2]);
    storage[2][0]= 1/storage_det * (storage[1][0]*storage[2][1]-storage[1][1]*storage[2][0]);
    storage[2][1]= 1/storage_det * (storage[0][1]*storage[2][0]-storage[0][0]*storage[2][1]);
    storage[2][2]= 1/storage_det * (storage[0][0]*storage[1][1]-storage[0][1]*storage[1][0]);
    
 
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
          s_k[i] = - storage[i][j]*J[j][i]*residual_disp_[i]; //delta_phi[j];
        }
    }
    
    cout << "##########################################################################################################################" << endl;
    
    cout << "p_o_:             "<< endl;
    cout << p_o_[0] << " " << p_o_[1] << " " << p_o_[2] << endl;
    cout << endl;
    
    cout << "p_:               "<< endl;
    cout << p_[0] << " " << p_[1] << " " << p_[2] << endl;
    cout << endl;
 
    p_o_ = p_;
    
  for (unsigned int i=0; i<3; i++) {
    p_[i]=p_[i]+alpha * s_k[i];
  }
  
  cout << "phi:" << endl;
  cout << phi << endl;
  cout  << endl;
  

  cout << "delta_phi:" << endl;
  cout << delta_phi[0] << " " << delta_phi[1] << " " << delta_phi[2] << endl;
  cout  << endl;
  
  
  cout << "J:" << endl;
  cout << J[0][0] << " " << J[0][1] << " " << J[0][2] <<  endl;
  cout << J[1][0] << " " << J[1][1] << " " << J[1][2] <<  endl;
  cout << J[2][0] << " " << J[2][1] << " " << J[2][2] <<  endl;
  cout  << endl;
  
  cout << "JTJ:" << endl;
  cout << JTJ[0][0] << " " << JTJ[0][1] << " " << JTJ[0][2] <<  endl;
  cout << JTJ[1][0] << " " << JTJ[1][1] << " " << JTJ[1][2] <<  endl;
  cout << JTJ[2][0] << " " << JTJ[2][1] << " " << JTJ[2][2] <<  endl;  
  cout  << endl;

  cout << "storage_det" << endl;
  cout << storage_det << endl;
  cout  << endl;
  
  cout << "storage-1:" << endl;
  cout << storage[0][0] << " " << storage[0][1] << " " << storage[0][2] <<  endl;
  cout << storage[1][0] << " " << storage[1][1] << " " << storage[1][2] <<  endl;
  cout << storage[2][0] << " " << storage[2][1] << " " << storage[2][2] <<  endl;
  cout  << endl;
  
  cout << "s_k" << endl;
  cout << s_k[0] << " " << s_k[1] << " " << s_k[2] <<  endl;
  cout  << endl; 
  
  
  cout << "final_disp_:      "<< endl;
  cout << final_disp_[0] <<  " "<< final_disp_[1]<< " " << final_disp_[2] << endl;
  cout << endl;
  
  cout << "residual_disp_o_: " << endl;
  cout << residual_disp_o_[0] << " "<< residual_disp_o_[1] << " "<< residual_disp_o_[2] << endl;
  cout << endl;
  cout << "residual_disp_:   " << endl;
  cout << residual_disp_[0] << " "<< residual_disp_[1] << " "<< residual_disp_[2] << endl;
  cout << endl;
  
  cout << "p_o_:             "<< endl;
  cout << p_o_[0] << " " << p_o_[1] << " " << p_o_[2] << endl;
  cout << endl;
  
  cout << "p_:               "<< endl;
  cout << p_[0] << " " << p_[1] << " " << p_[2] << endl;
  cout << endl;
  
  cout << "mu " << mu << endl;
  
  

  DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->c      = p_[0];
  DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k1     = p_[1];
  DRT::Problem::Instance()->Material(0).m.hyper_polyconvex->k2     = p_[2];
  

  return;
}


#endif

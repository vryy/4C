/*!----------------------------------------------------------------------
\file  mc_utils.cpp
\brief Some utility functions for monte carlo analysis


 <pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
 *!----------------------------------------------------------------------*/

#include "mlmc.H"
#include "randomfield.H"
#include "randomfield_fourier.H"
#include "randomfield_spectral.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/io_control.H"

#include "../drt_inpar/inpar_material.H"
#include "../drt_inpar/inpar_mlmc.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"


void UQ::MLMC::WriteStdVectorToFile(std::vector<double> myvector, std::string FileNameWithPath)
{
  const char * c = FileNameWithPath.c_str();

  if (discret_->Comm().MyPID() == 0)
  {
    std::ofstream File;
    File.open(c,std::ios::out);
    int size = myvector.size();
    for(int i=0;i<size;i++)
    {
      File << myvector[i]<< std::endl;
    }
    File.close();
  }
}

//-----------------------------------------------------------------------------------
double UQ::MLMC::CheckIfNodeInElement(DRT::Node& node, DRT::Element& ele, double* xsi)
{
  // init speudo distance, which is essentially largest values of xsi[i]
  double pseudo_distance = 0.0;
  //local/element coordinates
  xsi[0] = xsi[1] = xsi[2] = 0.0 ;
  // Replace later with problem dimension or element type
  //if (Dim()==3)

  if (ele.Shape()==DRT::Element::hex8)
    {
      // function f (vector-valued)
      double f[3] = {0.0, 0.0, 0.0};
      LINALG::Matrix<3,1> b;
      // gradient of f (df/dxsi[0], df/dxsi[1], df/dxsi[2])
      LINALG::Matrix<3,3> df;
      //convergeence check
      double conv = 0.0;

      for (int k=0;k<num_newton_it_;++k)
        {
          EvaluateF(f,node,ele,xsi);
          conv = sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]);
          //IO::cout << "Iteration " << k << ": -> |f|=" << conv << IO::endl;
          if (conv <= convtol_) break;

          EvaluateGradF(df,node,ele,xsi);

          // solve dxsi = - inv(df) * f
          df.Invert();


          // update xsi
          xsi[0] += -df(0,0)*f[0] - df(1,0)*f[1] - df(2,0)*f[2];
          xsi[1] += -df(0,1)*f[0] - df(1,1)*f[1] - df(2,1)*f[2];
          xsi[2] += -df(0,2)*f[0] - df(1,2)*f[1] - df(2,2)*f[2];
        }

      // Newton iteration unconverged
      if (conv > convtol_)
      {
        dserror("ERROR: CheckIfNodeInElement: Newton unconverged for NodeID %i "
                   "and ElementID %i", node.Id(), ele.Id());
      }
      // Newton iteration converged
      // find largest value of xsi[i] and return as pseudo_distance
      for (int i = 0; i < 3; i++)
      {
        if (fabs(xsi[i]) > pseudo_distance)
        {
          pseudo_distance = fabs(xsi[i]);
        }
      }
      return pseudo_distance;

    }
    else
    {
      dserror("CheckIfNodeInElement only implememted for hex8 Elements");
      return pseudo_distance;
    }
}

//-----------------------------------------------------------------------------------
bool UQ::MLMC::EvaluateF(double* f,DRT::Node& node, DRT::Element& ele,const double* xsi)
{

  LINALG::Matrix<8,1> funct ;
  //DRT::Element::DiscretizationType distype = ele.DiscretizationType;
  const DRT::Element::DiscretizationType distype = DRT::Element::hex8;
  DRT::UTILS::shape_function_3D(funct, xsi[0], xsi[1], xsi[2],distype);

  //LINALG::Matrix<ele.NumNode(),ele->numdim> xrefe;  // material coord. of element
  LINALG::Matrix<8,3> xrefe;  // material coord. of element
    for (int i=0; i<ele.NumNode(); ++i){
      const double* x =ele.Nodes()[i]->X();
      xrefe(i,0) = x[0];
      xrefe(i,1) = x[1];
      xrefe(i,2) = x[2];
    }
    // Calc Difference in Location
    LINALG::Matrix<1,3> point;
    point.MultiplyTN(funct, xrefe);
    f[0]=point(0,0)-node.X()[0];
    f[1]=point(0,1)-node.X()[1];
    f[2]=point(0,2)-node.X()[2];
    return true;
}


//-----------------------------------------------------------------------------------
bool UQ::MLMC::EvaluateGradF(LINALG::Matrix<3,3>& fgrad,DRT::Node& node, DRT::Element& ele,const double* xsi)
{
  //static const int iel = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
  LINALG::Matrix<3,8> deriv1;
  const DRT::Element::DiscretizationType distype = DRT::Element::hex8;

  DRT::UTILS::shape_function_3D_deriv1(deriv1 ,xsi[0], xsi[1], xsi[2],distype);

  LINALG::Matrix<8,3> xrefe;  // material coord. of element
  int NUMNOD_SOH8 = 8;
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    const double* x =ele.Nodes()[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];
  }

  fgrad.MultiplyNN(deriv1,xrefe);

  return true;
}

//-----------------------------------------------------------------------------------
int UQ::MLMC::FindBackgroundElement(DRT::Node node, Teuchos::RCP<DRT::Discretization> background_dis, int* bg_ele_id, double* xsi)
{
  bool outlier = false;
  Teuchos::RCP<Epetra_Vector> element_vector = Teuchos::rcp(new Epetra_Vector(*(background_dis->ElementRowMap()),true));
  // loop over discretization

  double pseudo_distance = 0.0;
  double min_pseudo_distance = 2.0;
  int back_ground_ele_id = 0;
  double background_xsi[3] = {0,0,0};
  //for (int i=0; i< 2 && min_pseudo_distance > InEleRange_ ; i++  )
  for (int i=0; i< element_vector->MyLength() && min_pseudo_distance > InEleRange_ ; i++  )
  {
    int globel_ele_id= background_dis->ElementRowMap()->GID(i);
    DRT::Element* ele = background_dis->gElement(globel_ele_id);
    //inEle = CheckIfNodeInElement(node, *ele);
    pseudo_distance = CheckIfNodeInElement(node, *ele, xsi);
    //IO::cout << " pseudo_distance  " << pseudo_distance  << IO::endl;
    // check if node is in Element
    if (pseudo_distance < min_pseudo_distance)
    {
      min_pseudo_distance = pseudo_distance;
      back_ground_ele_id = globel_ele_id;
      background_xsi[0]=xsi[0];
      background_xsi[1]=xsi[1];
      background_xsi[2]=xsi[2];
    }
  } // end of loop over elements
  // Debug
  if(min_pseudo_distance < InEleRange_)
  {
    //IO::cout << "found background element Ele ID is " <<  back_ground_ele_id << IO::endl;
    //IO::cout << "Local Coordinates are "<< "xsi_0 " << xsi[0] << " xsi_1 " << xsi[1] << " xsi_2 " << xsi[2] << IO::endl;

  }
  else
  {
   //IO::cout << "did not find background element, closest element is: Ele ID: " << back_ground_ele_id << IO::endl;
   //IO::cout << "Local Coordinates are "<< "xsi_0 " << background_xsi[0] << " xsi_1 " << background_xsi[1] << " xsi_2 " << background_xsi[2] << IO::endl;
   // dserror("stop right here");
   // write closest element xsi* into  pointer
    xsi[0]=background_xsi[0];
    xsi[1]=background_xsi[1];
    xsi[2]=background_xsi[2];
    outlier =true;
  }
  *bg_ele_id = back_ground_ele_id;
  return outlier;


}

void UQ::MLMC::WriteStatOutput()
{  //
  std::stringstream name_helper;
  name_helper << filename_ << "_statistics";
  output_->NewResultFile(name_helper.str(),numb_run_);
  output_->WriteMesh(0, 0.01);
  output_->NewStep( 1, 0.01);
  output_->WriteVector("mean_displacements", mean_disp_, output_->dofvector);
  output_->WriteVector("variance_displacements", var_disp_, output_->dofvector);
  output_->WriteVector("mean_gauss_2PK_stresses_xyz", mean_stress_, output_->nodevector);
  output_->WriteVector("variance_gauss_2PK_stresses_xyz", var_stress_, output_->nodevector);
  output_->WriteVector("mean_gauss_GL_strain_xyz", mean_strain_, output_->nodevector);
  output_->WriteVector("variance_gauss_GL_strain_xyz", var_strain_, output_->nodevector);
}

// calculate some statistics
void UQ::MLMC::CalcStatStressDisp(Teuchos::RCP< Epetra_MultiVector> curr_stress,Teuchos::RCP< Epetra_MultiVector> curr_strain,Teuchos::RCP<Epetra_MultiVector> curr_disp)
{
  // in order to avoid saving the stresses and displacements for each step an online
  // algorithm to compute the std deviation is needed .

  // Such an algorithm can be found here:

  // http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance

  // This it how it goes
  //def online_variance(data):
   //   n = 0
   //   mean = 0
   //   M2 = 0
   //
   //   for x in data:
    //      n = n + 1
    //      delta = x - mean
     //     mean = mean + delta/n
      //    M2 = M2 + delta*(x - mean)  # This expression uses the new value of mean
   //
     // variance_n = M2/n
     /// variance = M2/(n - 1)
    //  return variance*/
  // since numb_run_ does not start from zero anymore we need
  int n = numb_run_-start_run_+1;
  // calc mean and variance for displacement
  delta_disp_->Update(1.0,*curr_disp,-1.0,*mean_disp_,0.0);
  mean_disp_->Update(1.0/n,*delta_disp_,1.0);
  m2_helper_var_disp_->Update(1.0,*curr_disp,-1.0,*mean_disp_,0.0);
  m2_var_disp_->Multiply(1.0,*delta_disp_,*m2_helper_var_disp_ ,1.0);
  var_disp_->Update(1.0/(numb_run_),*m2_var_disp_,0.0);
  // do the same for stresses
  delta_stress_->Update(1.0,*curr_stress,-1.0,*mean_stress_,0.0);
  mean_stress_->Update(1.0/n,*delta_stress_,1.0);
  m2_helper_var_stress_->Update(1.0,*curr_stress,-1.0,*mean_stress_,0.0);
  m2_var_stress_->Multiply(1.0,*delta_stress_,*m2_helper_var_stress_ ,1.0);
  var_stress_->Update(1.0/(numb_run_),*m2_var_stress_,0.0);
  // and for strains
  delta_strain_->Update(1.0,*curr_strain,-1.0,*mean_strain_,0.0);
  mean_strain_->Update(1.0/n,*delta_strain_,1.0);
  m2_helper_var_strain_->Update(1.0,*curr_strain,-1.0,*mean_strain_,0.0);
  m2_var_strain_->Multiply(1.0,*delta_strain_,*m2_helper_var_strain_ ,1.0);
  var_strain_->Update(1.0/(numb_run_),*m2_var_strain_,0.0);
}

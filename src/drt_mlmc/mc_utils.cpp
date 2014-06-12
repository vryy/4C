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


void STR::UQ::MLMC::WriteStdVectorToFile(std::vector<double> myvector, std::string FileNameWithPath)
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

// all the functions needed for multilevel sampling


void STR::UQ::MLMC::CalcDifferenceToLowerLevel(Teuchos::RCP< Epetra_MultiVector> stress, Teuchos::RCP< Epetra_MultiVector> strain, Teuchos::RCP<Epetra_MultiVector> disp)
{
  // store difference not in new vectors to save memory
  disp_lower_level_->Update(1.0,*disp,-1.0);
  stress_lower_level_->Update(1.0,*stress,-1.0);
  strain_lower_level_->Update(1.0,*strain,-1.0);

}


//-----------------------------------------------------------------------------------
double STR::UQ::MLMC::CheckIfNodeInElement(DRT::Node& node, DRT::Element& ele, double* xsi)
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
bool STR::UQ::MLMC::EvaluateF(double* f,DRT::Node& node, DRT::Element& ele,const double* xsi)
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
bool STR::UQ::MLMC::EvaluateGradF(LINALG::Matrix<3,3>& fgrad,DRT::Node& node, DRT::Element& ele,const double* xsi)
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
// Read in Stresses and Displacements from corresponding Run on lower Level
void STR::UQ::MLMC::ReadResultsFromLowerLevel()
{
  std::stringstream name_helper;
  // assamble filename and pathfor cluster jobs
  //name_helper << "../../level"<<num_level_-1<< "/START_RUN_"<< start_run_ <<"/"<<filename_lower_level_<<"_prolongated"<< "_run_" << numb_run_ ;

  name_helper << filename_lower_level_ <<"_prolongated"<< "_run_" << numb_run_ ;
  std::string name_ll = name_helper.str();
  // check if bool should be true or false
  Teuchos::RCP<IO::InputControl> inputcontrol = Teuchos::rcp(new IO::InputControl(name_ll, false));
  IO::DiscretizationReader input_coarse(actdis_fine_, inputcontrol,1);

  // read in displacements
  input_coarse.ReadMultiVector(disp_lower_level_, "displacement");
  // read in prolongated GP Stresses
  input_coarse.ReadMultiVector(stress_lower_level_, "prolongated_gauss_2PK_stresses_xyz");
  // read in prolongated GP Strain
   input_coarse.ReadMultiVector(strain_lower_level_, "prolongated_gauss_GL_strains_xyz");
}


//---------------------------------------------------------------------------------------------
void STR::UQ::MLMC::SetupProlongatorParallel()
{
  dserror("no longer maintianed");
/*  // number of outliers
  int num_outliers = 0;
  const Epetra_Map* rmap_disp = NULL;
  const Epetra_Map* dmap_disp = NULL;
  const Epetra_Map* rmap_stress = NULL;
  const Epetra_Map* dmap_stress = NULL;
  rmap_disp= (actdis_fine_->DofRowMap());
  dmap_disp=(actdis_coarse_->DofRowMap());
  //maps for stress prolongator
  rmap_stress = (actdis_fine_->NodeRowMap());
  dmap_stress = (actdis_coarse_->NodeRowMap());
  int bg_ele_id;
  // store location of node
  double xsi[3] = {0.0, 0.0,0.0};
  prolongator_disp_crs_ = Teuchos::rcp(new Epetra_FECrsMatrix(::Copy,*rmap_disp,*dmap_disp,8,false));
  prolongator_stress_crs_ = Teuchos::rcp(new Epetra_FECrsMatrix(::Copy,*rmap_stress,*dmap_stress,8,false));

  //loop over nodes of fine discretization on this proc
  Teuchos::RCP<Epetra_Vector> node_vector = Teuchos::rcp(new Epetra_Vector(*(actdis_fine_->NodeRowMap()),true));
  // loop over dofs
  for (int i=0; i< actdis_fine_->NumMyRowNodes(); i++  )
  //  for (int i=0; i<2; i++  )
  {
    DRT::Node* node = actdis_fine_->lRowNode(i);

    // Get background element and local coordinates
    num_outliers += FindBackgroundElement(*node, actdis_coarse_, &bg_ele_id, xsi);

    // Get element
    DRT::Element* bg_ele = actdis_coarse_->gElement(bg_ele_id);

    Epetra_SerialDenseVector shape_fcts(bg_ele->NumNode());
    DRT::UTILS::shape_function_3D(shape_fcts,xsi[0],xsi[1],xsi[2],bg_ele->Shape());


    int* rows = NULL;
    int* cols = NULL;
    double* values = NULL;
    int numColumns = bg_ele->NumNode();
    // insert rows dof wise
    int numRows   = 1;


    cols = new int[numColumns];
    rows = new int[numRows];
    values = new double[numColumns];

    // fill prolongators
    // DIM = 3
    for (int j = 0; j<3 ; j++)
    {
      int index = 0;
      for (int k=0; k< bg_ele->NumNode(); k ++)
      {
        // store global indices in cols
       // (*rcp_Array)[index]=dmap_disp->GID((bg_ele->Nodes()[k]->Id()*3)+j);
        cols[index]=dmap_disp->GID((bg_ele->Nodes()[k]->Id()*3)+j);
        rows[0]=i*3+j;
        values[index]=shape_fcts[k];
        index++;
      }
      int err=  prolongator_disp_crs_->InsertGlobalValues(1,rows,8,cols,values,Epetra_FECrsMatrix::COLUMN_MAJOR);
      if (err != 0)
      {
        dserror("Could not insert global values");
      }

    } // loop j
    // stress prolongator
    for (int k=0; k< bg_ele->NumNode(); k ++)
    {
      // store global indices in cols
      cols[k]=dmap_stress->GID(bg_ele->Nodes()[k]->Id());
      values[k]=shape_fcts[k];
    }
    rows[0]=i;
    int err=  prolongator_stress_crs_->InsertGlobalValues(1,rows,8,cols,values,Epetra_FECrsMatrix::COLUMN_MAJOR);
    if (err != 0)
    {
      dserror("Could not insert global values");
    }

    delete [] cols;
    delete [] rows;
    delete [] values;

    rows = NULL;
    cols = NULL;
    values = NULL;

  } // End of loop over nodes of fine discretzation

  // Assembly
  prolongator_disp_crs_->GlobalAssemble(*dmap_disp,*rmap_disp,true);
  prolongator_stress_crs_->GlobalAssemble(*dmap_stress,*rmap_stress,true);
  IO::cout << "################################################### " << IO::endl;
  IO::cout << "   SUCCESSFULLY INITIALIZED  PROLONGATOR" << IO::endl;
  IO::cout <<  num_outliers << " Nodes do not lie within a background element " << IO::endl;
  IO::cout << "################################################### " << IO::endl; */
}
//---------------------------------------------------------------------------------------------
void STR::UQ::MLMC::SetupProlongator()
{
  // This functions calculates the prolongtators for the displacement and the nodal stresses

  // 3D Problem
  int num_columns_prolongator_disp = actdis_fine_->NumGlobalNodes()*3;
  int num_columns_prolongator_stress = actdis_fine_->NumGlobalNodes();

  double xsi[3] = {0.0, 0.0,0.0};

  // loop over nodes of fine dis
  int num_nodes;
  int bg_ele_id;
  num_nodes = actdis_fine_->NumGlobalNodes();

  // init prolongators
  IO::cout << "num_columns proongator "  << num_columns_prolongator_disp << IO::endl;
  IO::cout << "num_columns proongator stress  "  << num_columns_prolongator_stress << IO::endl;
  prolongator_disp_ = Teuchos::rcp(new Epetra_MultiVector(*actdis_coarse_->DofRowMap(),num_columns_prolongator_disp,true));
  prolongator_stress_ = Teuchos::rcp(new Epetra_MultiVector(*actdis_coarse_->NodeRowMap(),num_columns_prolongator_stress,true));

  for (int i = 0; i < num_nodes ; i++)
  {

    // Get node
    DRT::Node* node = actdis_fine_->gNode(i);

    // Get background element and local coordinates
    FindBackgroundElement(*node, actdis_coarse_, &bg_ele_id, xsi);
    // Get element
    DRT::Element* bg_ele = actdis_coarse_->gElement(bg_ele_id);

    Epetra_SerialDenseVector shape_fcts(bg_ele->NumNode());
    DRT::UTILS::shape_function_3D(shape_fcts,xsi[0],xsi[1],xsi[2],bg_ele->Shape());

    // fill prolongators add values to prolongator
    for (int j = 0; j<3 ; j++)
    {

      for (int k=0; k< bg_ele->NumNode(); k ++)
      {
        //prolongator_disp_ is dof based
        (*prolongator_disp_)[(i*3)+j][(bg_ele->Nodes()[k]->Id()*3)+j]= shape_fcts[k];
        // prolongator stress_ is node based
        (*prolongator_stress_)[i][(bg_ele->Nodes()[k]->Id())]= shape_fcts[k];
      }
    }
  } // End of loop over nodes of fine discretzation

}
//---------------------------------------------------------------------------------------------
void STR::UQ::MLMC::ProlongateResults()
{
  IO::cout << "Prolongating Resuts " << IO::endl;
  // To avoid messing with the timeintegration we read in the results of the coarse discretization here
  // Get coarse Grid problem instance
  std::stringstream name;
  std::string filename_helper;
  name << filename_ << "_run_"<< numb_run_;
  filename_helper = name.str();

  Teuchos::RCP<IO::InputControl> inputcontrol = Teuchos::rcp(new IO::InputControl(filename_helper, false));

  IO::DiscretizationReader input_coarse(actdis_coarse_, inputcontrol,tsteps_);
  // Vector for displacements
  Teuchos::RCP<Epetra_Vector> dis_coarse = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));

  // read in displacements
  input_coarse.ReadVector(dis_coarse, "displacement");

  Teuchos::RCP<Epetra_MultiVector> dis_fine = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->DofRowMap()),1,true));

  // Try new and shiny prolongator based on crs Matrix
  int error = prolongator_disp_crs_->Multiply(false,*dis_coarse,*dis_fine);
  if(error!=0)
  {
    dserror("stuff went wrong");
  }


  // create new resultfile for prolongated results
  std::stringstream name_prolong;
  std::string filename_helper_prolong;
  name_prolong << filename_ << "_prolongated";
  filename_helper_prolong = name_prolong.str();

  output_fine_->NewResultFile(filename_helper_prolong ,(numb_run_));


  output_fine_->WriteMesh(0, 0.01, meshfilename_);
  // exception for first run
  if(numb_run_==start_run_)
    output_fine_->WriteMesh(0, 0.01);

  output_fine_->NewStep( 1, 0.01);


  // Write interpolated displacement to file
  output_fine_->WriteVector("displacement", dis_fine, output_fine_->dofvector);
  Teuchos::RCP<Epetra_Vector> dis_fine_single = Teuchos::rcp(new Epetra_Vector(*actdis_fine_->DofRowMap(),true));

  // transfer to Epetra_Vector
  for( int i = 0;i< dis_fine->MyLength(); i++)
  {
   (*dis_fine_single)[i]= (*dis_fine)[0][i];
  }

  //#####################################################################################
  //
  //                  prolongate stresses based on interpolated displacement field
  //
  //#####################################################################################

  // set up parameters for Evaluation
  double timen         = 0.9;  // params_.get<double>("total time"             ,0.0);
  double dt            = 0.1; //params_.get<double>("delta time"             ,0.01);
  double alphaf        = 0.459; // params_.get<double>("alpha f"                ,0.459);
  INPAR::STR::StressType iostress =INPAR::STR::stress_2pk; //stress_none;
  INPAR::STR::StrainType iostrain= INPAR::STR::strain_gl; // strain_none;
  Teuchos::RCP<Epetra_Vector>    zeros_ = Teuchos::rcp(new Epetra_Vector(*actdis_fine_->DofRowMap(),true));
  Teuchos::RCP<Epetra_Vector>    dis_ = dis_fine_single;
  Teuchos::RCP<Epetra_Vector>    vel_ = Teuchos::rcp(new Epetra_Vector(*actdis_fine_->DofRowMap(),true));
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // action for elements

  p.set("action","calc_struct_stress");
  // other parameters that might be needed by the elements
  p.set("total time",timen);
  p.set("delta time",dt);
  p.set("alpha f",alphaf);

  Teuchos::RCP<std::vector<char> > stress = Teuchos::rcp(new std::vector<char>());
  Teuchos::RCP<std::vector<char> > strain = Teuchos::rcp(new std::vector<char>());
  // plastic strains need to be init as well
  Teuchos::RCP<std::vector<char> > plstrain = Teuchos::rcp(new std::vector<char>());
  p.set("stress", stress);
  p.set("plstrain",plstrain);
  //

  p.set<int>("iostress", iostress);
  p.set("strain", strain);
  p.set<int>("iostrain", iostrain);
  // set vector values needed by elements
  p.set<double>("random test",5.0);
  actdis_fine_->ClearState();
  actdis_fine_->SetState("residual displacement",zeros_);
  actdis_fine_->SetState("displacement",dis_);
  actdis_fine_->SetState("velocity",vel_);

  // Evaluate Stresses based on interpolated displacements
  //actdis_fine_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  //actdis_fine_->ClearState();
  // Write to file
  // interpolated stresses from disp field
  //output_fine_->WriteVector("gauss_2PK_stresses_xyz",*stress,*(actdis_fine_->ElementRowMap()));




  //#####################################################################################
  //
  //                  prolongate stresses based on interpolated nodal stress field
  //
  //#####################################################################################

  // use same parameter list as above
  Teuchos::RCP<Epetra_Vector>    zeros_coarse = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));
  //Teuchos::RCP<Epetra_Vector>    dis_ = dis_fine_single;
  Teuchos::RCP<Epetra_Vector>    vel_coarse = Teuchos::rcp(new Epetra_Vector(*actdis_coarse_->DofRowMap(),true));
  actdis_coarse_->ClearState();
  actdis_coarse_->SetState("residual displacement",zeros_coarse);
  actdis_coarse_->SetState("displacement",dis_coarse);
  actdis_coarse_->SetState("velocity",vel_coarse);
  // Alrigth lets get the nodal stresses
  p.set("action","calc_global_gpstresses_map");

  const Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > gpstressmap = Teuchos::rcp(new std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> >);
  p.set("gpstressmap", gpstressmap);

  const Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > gpstrainmap = Teuchos::rcp(new std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix> >);
  p.set("gpstrainmap", gpstrainmap);

  actdis_coarse_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  // st action to calc poststresse
  p.set("action","postprocess_stress");
  // Multivector to store poststresses
  Teuchos::RCP<Epetra_MultiVector> poststress =  Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse_->NodeRowMap()),6,true));
  // for fine diskretization as well
  Teuchos::RCP<Epetra_MultiVector> poststress_fine =  Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));

  p.set("poststress", poststress);
  p.set("stresstype","ndxyz");

  actdis_coarse_->ClearState();
  actdis_coarse_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  actdis_coarse_->ClearState();

  // again for strains
  p.set("action","postprocess_stress");
  // Multivector to store poststrains
  Teuchos::RCP<Epetra_MultiVector> poststrain =  Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse_->NodeRowMap()),6,true));
  // for fine diskretization as well
  Teuchos::RCP<Epetra_MultiVector> poststrain_fine =  Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine_->NodeRowMap()),6,true));
  p.set("poststress", poststrain);
  p.set("gpstressmap", gpstrainmap);
  p.set("stresstype","ndxyz");

  actdis_coarse_->ClearState();
  actdis_coarse_->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  actdis_coarse_->ClearState();


  // try new and shiny crs prolongator
  int error2 = prolongator_stress_crs_->Multiply(false,*poststress,*poststress_fine);
  if(error2!=0)
  {
    dserror("stuff went wrong");
  }
  // strains as well
  int error3 = prolongator_stress_crs_->Multiply(false,*poststrain,*poststrain_fine);
  if(error3!=0)
  {
    dserror("stuff went wrong");
  }

  IO::cout << " before prolongated gauss 2pk stresses " << IO::endl;
  output_fine_->WriteVector("prolongated_gauss_2PK_stresses_xyz", poststress_fine, output_fine_->nodevector);
  output_fine_->WriteVector("prolongated_gauss_GL_strains_xyz", poststrain_fine, output_fine_->nodevector);
  // do some statistics
  if(calc_diff_)
  {
    CalcDifferenceToLowerLevel(poststress_fine,poststrain_fine, dis_fine);
    // Write Difference between two Discretizations to File
    output_fine_->WriteVector("diff_to_ll_displacement", disp_lower_level_, output_fine_->dofvector);
    output_fine_->WriteVector("diff_to_ll_prolongated_gauss_2PK_stresses_xyz", stress_lower_level_, output_fine_->nodevector);
    output_fine_->WriteVector("diff_to_ll_prolongated_gauss_GL_strains_xyz", strain_lower_level_, output_fine_->nodevector);
  }
  CalcStatStressDisp(poststress_fine,poststrain_fine, dis_fine);

}
//-----------------------------------------------------------------------------------
int STR::UQ::MLMC::FindBackgroundElement(DRT::Node node, Teuchos::RCP<DRT::Discretization> background_dis, int* bg_ele_id, double* xsi)
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

void STR::UQ::MLMC::WriteStatOutput()
{  //
  std::stringstream name_helper;
  name_helper << filename_ << "_statistics";
  output_fine_->NewResultFile(name_helper.str(),numb_run_);
  output_fine_->WriteMesh(0, 0.01);
  output_fine_->NewStep( 1, 0.01);
  output_fine_->WriteVector("mean_displacements", mean_disp_, output_fine_->dofvector);
  output_fine_->WriteVector("variance_displacements", var_disp_, output_fine_->dofvector);
  output_fine_->WriteVector("mean_gauss_2PK_stresses_xyz", mean_stress_, output_fine_->nodevector);
  output_fine_->WriteVector("variance_gauss_2PK_stresses_xyz", var_stress_, output_fine_->nodevector);
  output_fine_->WriteVector("mean_gauss_GL_strain_xyz", mean_strain_, output_fine_->nodevector);
  output_fine_->WriteVector("variance_gauss_GL_strain_xyz", var_strain_, output_fine_->nodevector);
  // write stats with respect to lower level
  if (calc_diff_)
  {
    output_fine_->WriteVector("diff_mean_displacements", diff_mean_disp_, output_fine_->dofvector);
    output_fine_->WriteVector("diff_variance_displacements", diff_var_disp_, output_fine_->dofvector);
    output_fine_->WriteVector("diff_mean_gauss_2PK_stresses_xyz", diff_mean_stress_, output_fine_->nodevector);
    output_fine_->WriteVector("diff_variance_gauss_2PK_stresses_xyz", diff_var_stress_, output_fine_->nodevector);
    output_fine_->WriteVector("diff_mean_gauss_GL_strain_xyz", diff_mean_strain_, output_fine_->nodevector);
    output_fine_->WriteVector("diff_variance_gauss_GL_strain_xyz", diff_var_strain_, output_fine_->nodevector);
   }
}

// calculate some statistics
void STR::UQ::MLMC::CalcStatStressDisp(Teuchos::RCP< Epetra_MultiVector> curr_stress,Teuchos::RCP< Epetra_MultiVector> curr_strain,Teuchos::RCP<Epetra_MultiVector> curr_disp)
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

  // quick check if we need difference stats
  if (calc_diff_)
  {
  // calc mean and variance for difference between levels
  diff_delta_disp_->Update(1.0,*disp_lower_level_,-1.0,*diff_mean_disp_,0.0);
  diff_mean_disp_->Update(1.0/n,*diff_delta_disp_,1.0);
  diff_m2_helper_var_disp_->Update(1.0,*disp_lower_level_,-1.0,*diff_mean_disp_,0.0);
  diff_m2_var_disp_->Multiply(1.0,*diff_delta_disp_,*diff_m2_helper_var_disp_ ,1.0);
  diff_var_disp_->Update(1.0/(numb_run_),*diff_m2_var_disp_,0.0);
  // do the same for stresses
  diff_delta_stress_->Update(1.0,*stress_lower_level_,-1.0,*diff_mean_stress_,0.0);
  diff_mean_stress_->Update(1.0/n,*diff_delta_stress_,1.0);
  diff_m2_helper_var_stress_->Update(1.0,*stress_lower_level_,-1.0,*diff_mean_stress_,0.0);
  diff_m2_var_stress_->Multiply(1.0,*diff_delta_stress_,*diff_m2_helper_var_stress_ ,1.0);
  diff_var_stress_->Update(1.0/(numb_run_),*diff_m2_var_stress_,0.0);
  // do the same for strains
  diff_delta_strain_->Update(1.0,*strain_lower_level_,-1.0,*diff_mean_strain_,0.0);
  diff_mean_strain_->Update(1.0/n,*diff_delta_strain_,1.0);
  diff_m2_helper_var_strain_->Update(1.0,*strain_lower_level_,-1.0,*diff_mean_strain_,0.0);
  diff_m2_var_strain_->Multiply(1.0,*diff_delta_strain_,*diff_m2_helper_var_strain_ ,1.0);
  diff_var_strain_->Update(1.0/(numb_run_),*diff_m2_var_strain_,0.0);
  }



}

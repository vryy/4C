/*!----------------------------------------------------------------------
\file  mlmc.cpp
\brief Class for performing Multi Level Monte Carlo (MLMC)analysis of structure


 <pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
 *!----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "mlmc.H"
#include <ctime>
#include <cstdlib>
#include <iostream>
#include "Epetra_SerialDenseMatrix.h"
#include "../drt_lib/global_inp_control2.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_io/io_hdf.H"
#include "../drt_mat/material.H"
#include "../drt_mat/aaaneohooke_stopro.H"
#include "../drt_structure/strtimint_create.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

using namespace std;
using namespace DRT;
using namespace MAT;

#include "randomfield.H"
#include "../drt_structure/stru_resulttest.H"


/*----------------------------------------------------------------------*/
/* standard constructor */
STR::MLMC::MLMC(Teuchos::RCP<DRT::Discretization> dis,
                                    Teuchos::RCP<LINALG::Solver> solver,
                                    Teuchos::RCP<IO::DiscretizationWriter> output)
  : discret_(dis),
    solver_(solver),
    output_(output),
    sti_(Teuchos::null)
{
  //int myrank = dis->Comm().MyPID();

  reset_out_count_=0;

  // input parameters structural dynamics
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  // get number of timesteps
  tsteps_ = sdyn.get<int>("NUMSTEP");
  // input parameters multi level monte carlo
  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();

  // Get number of Newton iterations
  num_newton_it_ = mlmcp.get<int>("ITENODEINELE");
  // Get convergence tolerance
  convtol_    = mlmcp.get<double>("CONVTOL");
  // read material parameters from input file

  // In element critirion xsi_i < 1 + eps  eps = MLMCINELETOL
  InEleRange_ = 1.0 + 10e-3;
  //ReadInParameters();

  // controlling parameter
  numb_run_ =  0;     // counter of how many runs were made monte carlo



}


/*----------------------------------------------------------------------*/
/* analyse */
void STR::MLMC::Integrate()
{
  //int myrank = discret_->Comm().MyPID();
  unsigned int random_seed = 1212213;

  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();
  int numruns = mlmcp.get<int>("NUMRUNS");
  do
  {
    // setup stoch mat

    SetupStochMat((random_seed+(unsigned int)numb_run_));
    output_->NewResultFile((numb_run_));


     // get input lists
     const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

     // major switch to different time integrators
     switch (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP"))
     {
     case INPAR::STR::dyna_centr_diff:
       dserror("no central differences in DRT");
       break;
     case INPAR::STR::dyna_gen_alfa:
     case INPAR::STR::dyna_gen_alfa_statics:
     case INPAR::STR::dyna_statics:
     case INPAR::STR::dyna_genalpha:
     case INPAR::STR::dyna_onesteptheta:
     case INPAR::STR::dyna_gemm:
     case INPAR::STR::dyna_ab2:
     case INPAR::STR::dyna_euma :
     case INPAR::STR::dyna_euimsto :
       dyn_nlnstructural_drt();
       break;
     case INPAR::STR::dyna_Gen_EMM:
       dserror("GEMM not supported");
       break;
     default:
       dserror("unknown time integration scheme '%s'", sdyn.get<std::string>("DYNAMICTYP").c_str());
     }
     numb_run_++;
  } while (numb_run_< numruns);

  SetupProlongator();
  ProlongateResults();
  return;
}
//---------------------------------------------------------------------------------------------
void STR::MLMC::SetupProlongator()
{
  // This functions calculates the prolongtators for the displacement and the nodal stresses

  // get the two discretizations
  // Get finest Grid problem instance
  Teuchos::RCP<DRT::Discretization> actdis_fine = Teuchos::null;
  actdis_fine = DRT::Problem::Instance(1)->Dis(genprob.numsf, 0);
  // set degrees of freedom in the discretization
  if (not actdis_fine->Filled()) actdis_fine->FillComplete();

  // Get coarse Grid problem instance
  Teuchos::RCP<DRT::Discretization> actdis_coarse = Teuchos::null;
  actdis_coarse = DRT::Problem::Instance(0)->Dis(genprob.numsf, 0);
  // set degrees of freedom in the discretization
  if (not actdis_coarse->Filled()) actdis_coarse->FillComplete();


  // 3D Problem
  int num_columns_prolongator_disp = actdis_fine->NumGlobalNodes()*3;
  int num_columns_prolongator_stress = actdis_fine->NumGlobalNodes();

  double xsi[3] = {0.0, 0.0,0.0};

  // loop over nodes of fine dis
  int num_nodes;
  int bg_ele_id;
  num_nodes = actdis_fine->NumGlobalNodes();

  // init prolongators
  prolongator_disp_ = rcp(new Epetra_MultiVector(*actdis_coarse->DofRowMap(),num_columns_prolongator_disp,true));
  prolongator_stress_ = rcp(new Epetra_MultiVector(*actdis_coarse->NodeRowMap(),num_columns_prolongator_stress,true));
  for (int i = 0; i < num_nodes ; i++)
  {

    // Get node
    DRT::Node* node = actdis_fine->gNode(i);

    // Get background element and local coordinates
    FindBackgroundElement(*node, actdis_coarse, &bg_ele_id, xsi);
    // Get element
    DRT::Element* bg_ele = actdis_coarse->gElement(bg_ele_id);

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
void STR::MLMC::ProlongateResults()
{
  // To avoid messing with the timeintegration we read in the results of the coarse discretization here
  // Get coarse Grid problem instance
  // access the discretization
  Teuchos::RCP<DRT::Discretization> actdis_coarse = Teuchos::null;
  actdis_coarse = DRT::Problem::Instance(0)->Dis(genprob.numsf, 0);
  //set degrees of freedom in the discretization
   if (not actdis_coarse->Filled()) actdis_coarse->FillComplete();

  string name = "aaa_run_0";
  // context for output and restart
  // check if bool should be true or false
  RCP<IO::InputControl> inputcontrol = rcp(new IO::InputControl(name, false));  IO::DiscretizationReader input_coarse(actdis_coarse, inputcontrol,tsteps_);
  // Vector for displacements
  Teuchos::RCP<Epetra_Vector> dis_coarse = rcp(new Epetra_Vector(*actdis_coarse->DofRowMap(),true));

  // read in displacements
  input_coarse.ReadVector(dis_coarse, "displacement");
  //cout << "chech what we read in " << *dis_coarse << endl;

  // creae multivector
  Teuchos::RCP<Epetra_MultiVector> dis_coarse_mv = rcp(new Epetra_MultiVector(*actdis_coarse->DofRowMap(),1,true));
  // access the fine discretization

  Teuchos::RCP<DRT::Discretization> actdis_fine = Teuchos::null;
  actdis_fine = DRT::Problem::Instance(1)->Dis(genprob.numsf, 0);

  //set degrees of freedom in the discretization
  if (not actdis_fine->Filled()) actdis_fine->FillComplete();


  Teuchos::RCP<Epetra_MultiVector> dis_fine_helper = Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse->DofRowMap()),prolongator_disp_->NumVectors(),true));
  Teuchos::RCP<Epetra_MultiVector> dis_fine = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine->DofRowMap()),1,true));

  // Prolongate results
  dis_fine_helper->Multiply(1,*dis_coarse,*prolongator_disp_,1);
  //write standard output into new file

  // transfer to Multivector
  for (int n = 0; n< dis_fine_helper->NumVectors(); n++)
  {
    for ( int m =0; m < dis_fine_helper->MyLength(); m ++)
    {
         (*dis_fine)[0][n] += (*dis_fine_helper)[n][m];
    }
  }

  // context for output and restart
  Teuchos::RCP<IO::DiscretizationWriter> output_fine= Teuchos::rcp(new IO::DiscretizationWriter(actdis_fine));
  // somehow the order of these following commands is important DONT CHANGE
  output_fine->NewResultFile((23212));
  output_fine->WriteMesh(1, 0.01);
  output_fine->NewStep( 1, 0.01);

  // Write interpolated displacement to file
  output_fine->WriteVector("displacement", dis_fine, output_fine->dofvector);
  Teuchos::RCP<Epetra_Vector> dis_fine_single = rcp(new Epetra_Vector(*actdis_fine->DofRowMap(),true));

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
  RCP<Epetra_Vector>    zeros_ = rcp(new Epetra_Vector(*actdis_fine->DofRowMap(),true));
  RCP<Epetra_Vector>    dis_ = dis_fine_single;
  RCP<Epetra_Vector>    vel_ = rcp(new Epetra_Vector(*actdis_fine->DofRowMap(),true));
  // create the parameters for the discretization
   ParameterList p;
   // action for elements
   p.set("action","calc_struct_stress");
   // other parameters that might be needed by the elements
   p.set("total time",timen);
   p.set("delta time",dt);
   p.set("alpha f",alphaf);

   Teuchos::RCP<std::vector<char> > stress = Teuchos::rcp(new std::vector<char>());
   Teuchos::RCP<std::vector<char> > strain = Teuchos::rcp(new std::vector<char>());
   p.set("stress", stress);

   p.set<int>("iostress", iostress);
   p.set("strain", strain);
   p.set<int>("iostrain", iostrain);
   // set vector values needed by elements
   p.set<double>("random test",5.0);
   actdis_fine->ClearState();
   actdis_fine->SetState("residual displacement",zeros_);
   actdis_fine->SetState("displacement",dis_);
   actdis_fine->SetState("velocity",vel_);
   // Evaluate Stresses based on interpolated displacements
   actdis_fine->Evaluate(p,null,null,null,null,null);
   actdis_fine->ClearState();
   // Write to file
   output_fine->WriteVector("gauss_2PK_stresses_xyz",*stress,*(actdis_fine->ElementRowMap()));


   cout << " Debugging   LINE   " << __LINE__ << endl;

   //#####################################################################################
   //
   //                  prolongate stresses based on interpolated nodal stress field
   //
   //#####################################################################################

   // use same parameter list as above
   RCP<Epetra_Vector>    zeros_coarse = rcp(new Epetra_Vector(*actdis_coarse->DofRowMap(),true));
   //RCP<Epetra_Vector>    dis_ = dis_fine_single;
   RCP<Epetra_Vector>    vel_coarse = rcp(new Epetra_Vector(*actdis_coarse->DofRowMap(),true));
   actdis_coarse->ClearState();
   actdis_coarse->SetState("residual displacement",zeros_coarse);
   actdis_coarse->SetState("displacement",dis_coarse);
   actdis_coarse->SetState("velocity",vel_coarse);
   // Alrigth lets get the nodal stresses

   p.set("action","calc_gobal_gpstresses_map");
  // cout << " Debugging   LINE   " << __LINE__ << endl;
   const RCP<map<int,RCP<Epetra_SerialDenseMatrix> > > gpstressmap = rcp(new std::map<int, RCP<Epetra_SerialDenseMatrix> >);
   p.set("gpstressmap", gpstressmap);

   actdis_coarse->Evaluate(p,null,null,null,null,null);
   //actdis_coarse->ClearState();
   //cout << " Debugging   LINE   " << __LINE__ << endl;
 // const RCP<std::map<int,RCP<Epetra_SerialDenseMatrix> > > data = stress;
   //p.set("gpstressmap", data);
   // st action to calc poststresse
   p.set("action","postprocess_stress");
   // Multivector to store poststresses
   RCP<Epetra_MultiVector> poststress =  Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse->NodeRowMap()),6,true));
   // for fine diskretization as well
   RCP<Epetra_MultiVector> poststress_fine =  Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine->NodeRowMap()),6,true));
   p.set("poststress", poststress);
   p.set("stresstype","ndxyz");
  // cout << " Debugging   LINE   " << __LINE__ << endl;
   actdis_coarse->ClearState();
   actdis_coarse->Evaluate(p,null,null,null,null,null);
   actdis_coarse->ClearState();


    Teuchos::RCP<IO::DiscretizationWriter> output_coarse= Teuchos::rcp(new IO::DiscretizationWriter(actdis_coarse));
    // somehow the order of these following commands is important DONT CHANGE
   //output_coarse->NewResultFile((23213));
   //output_coarse->WriteMesh(1, 0.01);
   //output_coarse->NewStep( 1, 0.01);
   //output_coarse->WriteVector("prolongated_gauss_2PK_stresses_xyz", poststress, output_coarse->nodevector);
   //Teuchos::RCP<Epetra_Vector> dis_fine_single = rcp(new Epetra_Vector(*actdis_fine->DofRowMap(),true));

   Teuchos::RCP<Epetra_MultiVector> poststress_fine_helper = Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse->NodeRowMap()),prolongator_stress_->NumVectors(),true));
  // Teuchos::RCP<Epetra_MultiVector> poststress_fine = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine->NodeRowMap()),6,true));

   // defin poststress coarse helper which stores the ith vector of postress in a new multivector which contains only on vector
   //Teuchos::RCP<Epetra_MultiVector> poststress_helper = Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse->NodeRowMap()),1,true));
   Epetra_DataAccess CV = View;

   //
   // Prolongate results
   //poststress_helper=(*poststress)(0);
   for(int i = 0; i<6; i++)
   {
     Teuchos::RCP<Epetra_MultiVector> poststress_helper = Teuchos::rcp(new Epetra_MultiVector(CV,*poststress,i,1));
     poststress_fine_helper->Multiply(1,*poststress_helper,*prolongator_stress_,0);
    // cout << "*prolongator_stress_ "  << *prolongator_stress_ << endl;
     //cout << "poststress_helper "  << *poststress_helper << endl;
     //cout << "poststress_fine_helper "  << *poststress_fine_helper << endl;
     //write standard output into new file

     // transfer to Multivector
     for (int n = 0; n< poststress_fine_helper->NumVectors(); n++)
     {
       for ( int m =0; m < poststress_fine_helper->MyLength(); m ++)
       {
            (*poststress_fine)[i][n] += (*poststress_fine_helper)[n][m];
       }
     }
   }
   //cout << "posttress  " << *poststress <<endl;
  // cout << "poststress_fine "  << *poststress_fine << endl;
   output_fine->WriteVector("prolongated_gauss_2PK_stresses_xyz", poststress_fine, output_fine->nodevector);


}
//-----------------------------------------------------------------------------------
void STR::MLMC::FindBackgroundElement(DRT::Node node, Teuchos::RCP<DRT::Discretization> background_dis, int* bg_ele_id, double* xsi)
{

  Teuchos::RCP<Epetra_Vector> element_vector = Teuchos::rcp(new Epetra_Vector(*(background_dis->ElementRowMap()),true));
  // loop over discretization

  double pseudo_distance = 0.0;
  double min_pseudo_distance = 2.0;
  int back_ground_ele_id = 0;
  for (int i=0; i< element_vector->MyLength() && min_pseudo_distance > InEleRange_ ; i++  )
  {
    int globel_ele_id= background_dis->ElementRowMap()->GID(i);
    DRT::Element* ele = background_dis->gElement(globel_ele_id);
    //inEle = CheckIfNodeInElement(node, *ele);
    pseudo_distance = CheckIfNodeInElement(node, *ele, xsi);
    // check if node is in Element
    if (pseudo_distance < min_pseudo_distance)
    {
      min_pseudo_distance = pseudo_distance;
      back_ground_ele_id = globel_ele_id;
    }
  } // end of loop over elements
  // Debug
  if(min_pseudo_distance < InEleRange_)
  {
    //cout << "found background element Ele ID is " <<  back_ground_ele_id << endl;
    //cout << "Local Coordinates are "<< "xsi_0 " << xsi[0] << " xsi_1 " << xsi[1] << " xsi_2 " << xsi[2] << endl;
  }
  else
  {
    cout << "did not find background element, closet element is: Ele ID: " << back_ground_ele_id << endl;
    cout << "Local Coordinates are "<< "xsi_0 " << xsi[0] << " xsi_1 " << xsi[1] << " xsi_2 " << xsi[2] << endl;
  }
  *bg_ele_id = back_ground_ele_id;

}

//-----------------------------------------------------------------------------------
double STR::MLMC::CheckIfNodeInElement(DRT::Node& node, DRT::Element& ele, double* xsi)
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
          //cout << "Iteration " << k << ": -> |f|=" << conv << endl;
          if (conv <= convtol_) break;

          EvaluateGradF(df,node,ele,xsi);

          // solve dxsi = - inv(df) * f
          df.Invert();
          //cout << "xsi_0 " << xsi[0] << "xsi_1 " << xsi[1] << "xsi_2 " << xsi[2] << endl;

          // update xsi
          xsi[0] += -df(0,0)*f[0] - df(1,0)*f[1] - df(2,0)*f[2];
          xsi[1] += -df(0,1)*f[0] - df(1,1)*f[1] - df(2,1)*f[2];
          xsi[2] += -df(0,2)*f[0] - df(1,2)*f[1] - df(2,2)*f[2];
          //cout << "iteration " << k<< "xsi: " << xsi[0] <<" " << xsi[1] << " "<<  xsi[2]<<endl ;
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
bool STR::MLMC::EvaluateF(double* f,DRT::Node& node, DRT::Element& ele,const double* xsi)
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
bool STR::MLMC::EvaluateGradF(LINALG::Matrix<3,3>& fgrad,DRT::Node& node, DRT::Element& ele,const double* xsi)
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

// Setup Material Parameters in each element based on Random Field
void STR::MLMC::SetupStochMat(unsigned int random_seed)
{
  // Variables for Random field
  double sigma =0.0 , corrlength = 0.0, beta = 0.0 ,beta_mean = 0.0;

  // Get parameters from stochastic matlaw
  const int myrank = discret_->Comm().MyPID();
  // loop all materials in problem
  const map<int,RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();
  if (myrank == 0) printf("No. material laws considered : %d\n",(int) mats.size());
  map<int,RCP<MAT::PAR::Material> >::const_iterator curr;
  for (curr=mats.begin(); curr != mats.end(); ++curr)
  {
    const RCP<MAT::PAR::Material> actmat = curr->second;
    switch(actmat->Type())
    {
       case INPAR::MAT::m_aaaneohooke_stopro:
       {
         MAT::PAR::AAAneohooke_stopro* params = dynamic_cast<MAT::PAR::AAAneohooke_stopro*>(actmat->Parameter());
         if (!params) dserror("Cannot cast material parameters");
         sigma = params->sigma_0_;
         corrlength = params->corrlength_;
         beta_mean =params->beta_mean_;
       }
       break;
      default:
      {
        // ignore unknown materials ?
        cout << "Mat Type   " << actmat->Type() << endl;
        dserror("Unknown type of material");

      }

    }
  } // EOF loop over mats

  // get elements on proc
  Teuchos::RCP<Epetra_Vector> my_ele = rcp(new Epetra_Vector(*discret_->ElementRowMap(),true));
  // loop over all elements
  for (int i=0; i< my_ele->MyLength(); i++)
    {
    if(discret_->gElement(i)->Material()->MaterialType()==INPAR::MAT::m_aaaneohooke_stopro)
      {
      MAT::AAAneohooke_stopro* aaa_stopro = static_cast <MAT::AAAneohooke_stopro*>(discret_->gElement(i)->Material().get());

      //double sigma = aaa_stopro->Sigma();
      //double corrlength = aaa_stopro->Corrlength();
      RandomField field(random_seed,sigma,corrlength);
      //beta = field.EvalRandomField(0.2,0.2,0.2);
      // get element centercoords
      DRT::Node** nodes = discret_->gElement(i)->Nodes();
      vector<double> ele_center;
      // init to zero
      ele_center.push_back(0.0);
      ele_center.push_back(0.0);
      ele_center.push_back(0.0);

      for (int i = 0; i < 8; i++ )
      {
        ele_center[0] += nodes[i]->X()[0]/8.;
        ele_center[1] += nodes[i]->X()[1]/8.;
        ele_center[2] += nodes[i]->X()[2]/8.;
      }
      // beta = beta_mean = beta_mean * random field value
      beta = beta_mean+beta_mean*field.EvalRandomField(ele_center[0],ele_center[1],ele_center[2]);
      //vector<double> location = dis->gElement(i)->soh8_ElementCenterRefeCoords();

      aaa_stopro->Init(beta);
      }
      cout << "beta  " << beta << endl;
    }
}



#endif

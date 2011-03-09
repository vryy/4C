/*----------------------------------------------------------------------*/
/*!
 * \file mlmc.cpp

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de

</pre>
*/
/*----------------------------------------------------------------------*/
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
#include "../drt_mat/lung_ogden.H"
#include "../drt_mat/lung_penalty.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_mat/neohooke.H"
#include "../drt_structure/strtimint_create.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_matelast/elast_coupanisoexpotwo.H"
#include "../drt_matelast/elast_coupanisoneohooketwo.H"
#include "../drt_matelast/elast_coupblatzko.H"
#include "../drt_matelast/elast_couplogneohooke.H"
#include "../drt_matelast/elast_isoexpo.H"
#include "../drt_matelast/elast_isomooneyrivlin.H"
#include "../drt_matelast/elast_isoneohooke.H"
#include "../drt_matelast/elast_isoyeoh.H"
#include "../drt_matelast/elast_volpenalty.H"
#include "../drt_matelast/elast_vologden.H"
#include "../drt_matelast/elast_volsussmanbathe.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

using namespace std;
using namespace DRT;
using namespace MAT;


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
  //const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  // input parameters multi level monte carlo
  //const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();

  // read material parameters from input file
  //ReadInParameters();
  // controlling parameter
  numb_run_ =  0;     // counter of how many runs were made monte carlo
}


/*----------------------------------------------------------------------*/
/* analyse */
void STR::MLMC::Integrate()
{
  //int myrank = discret_->Comm().MyPID();

  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();
  int numruns = mlmcp.get<int>("NUMRUNS");
  do
  {
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
  //MapResultsToFinestGrid();
  SetupProlongator();
  ProlongateResults();
  return;
}
//---------------------------------------------------------------------------------------------
void STR::MLMC::SetupProlongator()
{
  // This needs to be done only once, on one proc

  // get the two discretizations
  // Get finest Grid problem instance

  Teuchos::RCP<DRT::Discretization> actdis_fine = Teuchos::null;
  actdis_fine = DRT::Problem::Instance(1)->Dis(genprob.numsf, 0);
  // set degrees of freedom in the discretization
  if (not actdis_fine->Filled()) actdis_fine->FillComplete();

  //int num_rows_prolongator = actdis_fine->NumGlobalNodes()*3;

  // Get coarse Grid problem instance
  Teuchos::RCP<DRT::Discretization> actdis_coarse = Teuchos::null;
  actdis_coarse = DRT::Problem::Instance(0)->Dis(genprob.numsf, 0);
  // set degrees of freedom in the discretization
  if (not actdis_coarse->Filled()) actdis_coarse->FillComplete();


  // 3D Problem
  int num_columns_prolongator = actdis_fine->NumGlobalNodes()*3;

  double xsi[3] = {0.0, 0.0,0.0};

  // loop over nodes of fine dis
  int num_nodes;
  int bg_ele_id;
  num_nodes = actdis_fine->NumGlobalNodes();

  // init prolongator
  prolongator = rcp(new Epetra_MultiVector(*actdis_coarse->DofRowMap(),num_columns_prolongator,true));
  for (int i = 0; i < num_nodes ; i++)
  {

    // Get node
    DRT::Node* node = actdis_fine->gNode(i);
    //node = actdis_fine->gNode(72);
    // Get backgron element and local coordinates

    FindBackgroundElement(*node, actdis_coarse, &bg_ele_id, xsi);
    // Get element
    DRT::Element* bg_ele = actdis_coarse->gElement(bg_ele_id);
    //cout << "Background element" << bg_ele_id << endl;

    Epetra_SerialDenseVector shape_fcts(bg_ele->NumNode());
    DRT::UTILS::shape_function_3D(shape_fcts,xsi[0],xsi[1],xsi[2],bg_ele->Shape());

    // add values to prolongator
    for (int j = 0; j<3 ; j++)
    {

      for (int k=0; k< bg_ele->NumNode(); k ++)
      {
        (*prolongator)[(i*3)+j][(bg_ele->Nodes()[k]->Id()*3)+j]= shape_fcts[k];
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

  int last_step = 10;
  string name = "aaa_run_0";
  // context for output and restart
  // check if bool should be true or false
  RCP<IO::InputControl> inputcontrol = rcp(new IO::InputControl(name, false));  IO::DiscretizationReader input_coarse(actdis_coarse, inputcontrol,last_step);
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


  Teuchos::RCP<Epetra_MultiVector> dis_fine_helper = Teuchos::rcp(new Epetra_MultiVector(*(actdis_coarse->DofRowMap()),prolongator->NumVectors(),true));
  Teuchos::RCP<Epetra_MultiVector> dis_fine = Teuchos::rcp(new Epetra_MultiVector(*(actdis_fine->DofRowMap()),1,true));

  // Prolongate results
  dis_fine_helper->Multiply(1,*dis_coarse,*prolongator,1);
  //write standard output into new file


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
  //cout << " What do we save " << *dis_fine << endl;
  output_fine->WriteVector("displacement", dis_fine, output_fine->dofvector);
}
//---------------------------------------------------------------------------------------------
//void STR::MLMC::MapResultsToFinestGrid(Teuchos::RCP<DRT::Discretization> dis_in,
 //                                               Teuchos::RCP<DRT::Discretization> dis_out,
 //                                               int run_id)
void STR::MLMC::MapResultsToFinestGrid()
{
  // Get finest Grid problem instance
  // access the discretization
   Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
   actdis = DRT::Problem::Instance(0)->Dis(genprob.numsf, 0);

   // set degrees of freedom in the discretization
   if (not actdis->Filled()) actdis->FillComplete();

   // context for output and restart
   Teuchos::RCP<IO::DiscretizationWriter> output2
     = Teuchos::rcp(new IO::DiscretizationWriter(actdis));
   // write mesh to file
   output2->NewResultFile((2321));
   output2->WriteMesh(1, 0.01);
   //output2->NewResultFile((2321));
   output2->NewStep( 1, 0.01);

   // Write some vector to file
   // init vector
   Teuchos::RCP<Epetra_MultiVector> test_vector = Teuchos::rcp(new Epetra_MultiVector(*(actdis->DofRowMap()),1,true));
   // fill Vector with arbitrary values for testing
   for (int i = 0; i< test_vector->MyLength();i++)
   {
     (*test_vector)[0][i] = 13.0 ; //->ReplaceMyValue(0,i,13.0);

   }
   // fill with some data
 //  output2->WriteVector("displacement",test_vector,output2->dofvector);

 // cout << "Feine Disketisierung" << actdis<< endl;
 // int num_ele_dis_1 = actdis->NumGlobalElements();
  //cout << "num_ele_dis_1 "<< num_ele_dis_1 <<endl;

  // test some stuff
  Teuchos::RCP<Epetra_Vector> element_vector = Teuchos::rcp(new Epetra_Vector(*(actdis->ElementRowMap()),true));
  // loop over discretization
  bool inEle = false;
  for (int i=0; i< element_vector->MyLength() && !inEle ; i++  )
  {
   // cout <<  "Test what is in element vector "<< (*element_vector)[i] << endl;
    //actdis->ElementRowMap()->Print(std::ostream&);
   //cout<< *(actdis->ElementRowMap()) << "see what camoes here" << endl;
    // DRT::Element* ele = actdis->gElement((*element_vector)[i]);
   //inEle = CheckIfNodeInElement(node, ele);
  }
  // end of test some stuff
  //int globel_ele_id= actdis->NodeRowMap()->GID(12);
  //DRT::Node* node = actdis->gNode(12);
  //FindBackgroundElement(*node, coarse_dis);

}
//-----------------------------------------------------------------------------------
void STR::MLMC::FindBackgroundElement(DRT::Node node, Teuchos::RCP<DRT::Discretization> background_dis, int* bg_ele_id, double* xsi)
{
  double tol_inele = 1e-4;
  Teuchos::RCP<Epetra_Vector> element_vector = Teuchos::rcp(new Epetra_Vector(*(background_dis->ElementRowMap()),true));
  // loop over discretization

  double pseudo_distance = 0.0;
  double min_pseudo_distance = 2.0;
  int back_ground_ele_id = 0;
  for (int i=0; i< element_vector->MyLength() && min_pseudo_distance > 1.0+tol_inele ; i++  )
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

  if(min_pseudo_distance < 1.0+tol_inele)
  {
    cout << "found background element Ele ID is " <<  back_ground_ele_id << endl;
    cout << "Local Coordinates are "<< "xsi_0 " << xsi[0] << " xsi_1 " << xsi[1] << " xsi_2 " << xsi[2] << endl;
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
  // Replace later with problem dimension or element type
  //if (Dim()==3)

  if (true)
    {
      // start in the element center
      DRT::Element::DiscretizationType dt = ele.Shape();
      xsi[0] = xsi[1] = xsi[2] = 0.0 ;
      //if (dt==DRT::Element::tri3 || dt==DRT::Element::tri6)
      if (dt!=DRT::Element::hex8)
      {
        dserror("CheckIfNodeInElement only implememted for HEX8 Elements");
      }


      // function f (vector-valued)
      double f[3] = {0.0, 0.0, 0.0};
      LINALG::Matrix<3,1> b;
      // gradient of f (df/dxsi[0], df/dxsi[1], df/dxsi[2])
      LINALG::Matrix<3,3> df;

      // start iteration
      int k=0;

      // DEFINE SOME CONSTANT THAT LATER NEED TO GO IN THE INPUTFILE
      double MLMCCONVTOL = 10e-5;
      double MLMCNODETOELEMENTITER = 20;
      //convergeence check
      double conv = 0.0;


      for (k=0;k<MLMCNODETOELEMENTITER;++k)
        {
          EvaluateF(f,node,ele,xsi);
          conv = sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]);
          //cout << "Iteration " << k << ": -> |f|=" << conv << endl;
          if (conv <= MLMCCONVTOL) break;

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
      // In element critirion xsi_i < 1 + eps  eps = MLMCINELETOL
           double InEleTol = 1e-3;
           double InEleRange = 1.0 + InEleTol;
      // Newton iteration unconverged
      if (conv > MLMCCONVTOL)
      {
        cout << "ERROR: CheckIfNodeInElement: Newton unconverged for NodeID" << node.Id()<< " ElementID " <<  ele.Id() <<endl;;
        ///return 1000.0; // Node not in element
      }


        //dserror("ERROR: CheckIfNodeInElement: Newton unconverged for NodeID %i "
          //      "and ElementID %i", node.Id(), ele.Id());

      // Newton iteration converged

      //perform check wether node in element

      if(fabs(xsi[0])< InEleRange && fabs(xsi[1])< InEleRange && fabs(xsi[2])< InEleRange)
      {
        //cout << "xsi_0 " << xsi[0] << " xsi_1 " << xsi[1] << " xsi_2 " << xsi[2] << endl;
       // return true;
      }
      // find largest value of xsi[i]
      // init speudo distance
      double pseudo_distance = 0;
      for (int i = 0; i < 3; i++)
      {
        if (fabs(xsi[i]) > pseudo_distance)
        {
          pseudo_distance = fabs(xsi[i]);
        }
      }
      return pseudo_distance;

    }
     else dserror("CheckIfNodeInElement only implememted for hex8 Elements");
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


//-----------------------------------------------------------------------------------
void STR::MLMC::InterpolateResults(DRT::Node& node, DRT::Element& backgroundele,const double* xsi)
{
  // Get nodes for interpolation
  int num_nodes = backgroundele.NumNode();
  int node_id;
  // loop over all nodes
  for(int i = 0 ; i< num_nodes; i++)
  {
    node_id = backgroundele.NodeIds()[i];
  }
  //backgroundele.
  // allocate memory
  //int node_ids[num_nodes];
  //int backgroundnodes;// = {0, 0, 0, 0,0, 0, 0, 0};
 // *node_ids = backgroundele.NodeIds();

  //cout << *backgroundnodes << endl;
  // Get Background Discretization
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numsf, 0);
  // get rcp to current displacement
  Teuchos::RCP<Epetra_MultiVector> test_vector = Teuchos::rcp(new Epetra_MultiVector(*(actdis->DofRowMap()),1,true));
  RCP<const Epetra_Vector> test2 =Teuchos::rcp(new Epetra_Vector(*(actdis->DofRowMap()),true));;
  //check if there are any displacement
  cout << "have we dfound displacement vector" << endl;
  test2 = actdis->GetState("displacement");
  cout << "test2" << test2<< endl;
  if(actdis->HasState("dis"))
  {
    cout << "displacement Vector found" << endl;
  }

}
//-----------------------------------------------------------------------------------
void STR::MLMC::ReadInParameters()
{
  const int myrank = discret_->Comm().MyPID();

  // loop all materials in problem
  const map<int,RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();

  if (myrank == 0) printf("No. material laws considered : %d\n",mats.size());
  map<int,RCP<MAT::PAR::Material> >::const_iterator curr;
  for (curr=mats.begin(); curr != mats.end(); ++curr)
  {
    const RCP<MAT::PAR::Material> actmat = curr->second;
    switch(actmat->Type())
    {
      case INPAR::MAT::m_aaaneohooke:
      {
        MAT::PAR::AAAneohooke* params = dynamic_cast<MAT::PAR::AAAneohooke*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        const int j = p_.Length();
        p_.Resize(j+2);
        p_(j)   = params->youngs_;
        p_(j+1) = params->beta_;
        //p_(j+2) = params->nue_; // need also change resize above to invoke nue
      }
      break;
      case INPAR::MAT::m_neohooke:
      {
        MAT::PAR::NeoHooke* params = dynamic_cast<MAT::PAR::NeoHooke*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        const int j = p_.Length();
        p_.Resize(j+2);
        p_(j)   = params->youngs_;
        p_(j+1) = params->poissonratio_;
      }
      break;
      case INPAR::MAT::m_elasthyper:
      {
        MAT::PAR::ElastHyper* params = dynamic_cast<MAT::PAR::ElastHyper*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        const int nummat               = params->nummat_;
        const std::vector<int>* matids = params->matids_;
        for (int i=0; i<nummat; ++i)
        {
          const int id = (*matids)[i];
          const RCP<MAT::PAR::Material> actelastmat = mats.find(id)->second;
          switch (actelastmat->Type())
          {
            case INPAR::MAT::mes_isoyeoh:
            {
              MAT::ELASTIC::PAR::IsoYeoh* params = dynamic_cast<MAT::ELASTIC::PAR::IsoYeoh*>(actelastmat->Parameter());
              if (!params) dserror("Cannot cast material parameters");
              const int j = p_.Length();
              p_.Resize(j+3);
              p_(j)   = params->c1_;
              p_(j+1) = params->c2_;
              p_(j+2) = params->c3_;
            }
            break;
            case INPAR::MAT::mes_vologden:
            {
              MAT::ELASTIC::PAR::VolOgden* params = dynamic_cast<MAT::ELASTIC::PAR::VolOgden*>(actelastmat->Parameter());
              if (!params) dserror("Cannot cast material parameters");
              const int j = p_.Length();
              p_.Resize(j+1);
              p_(j)   = params->kappa_;
              // p_(j+1) = params->beta_; // need also change resize above to invoke beta_
            }
            break;
            default:
              dserror("Unknown type of elasthyper material");
            break;
          }
        }
      }
      case INPAR::MAT::mes_isoyeoh: // at this level do nothing, its inside the INPAR::MAT::m_elasthyper block
      break;
      case INPAR::MAT::mes_vologden: // at this level do nothing, its inside the INPAR::MAT::m_elasthyper block
      break;
      default:
        // ignore unknown materials ?
        dserror("Unknown type of material");
      break;
    }
  }
  return;
}
//--------------------------------------------------------------------------------------
void STR::MLMC::SetParameters(Epetra_SerialDenseVector p_cur)
{
  const int myrank = discret_->Comm().MyPID();

  // parameters are evaluated on proc 0 only? This should not be neccessary....
  discret_->Comm().Broadcast(&p_cur[0],np_,0);

  // loop all materials in problem
  const map<int,RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();
  int count=0;
  map<int,RCP<MAT::PAR::Material> >::const_iterator curr;
  for (curr=mats.begin(); curr != mats.end(); ++curr)
  {
    const RCP<MAT::PAR::Material> actmat = curr->second;
    switch(actmat->Type())
    {
      case INPAR::MAT::m_aaaneohooke:
      {
        MAT::PAR::AAAneohooke* params = dynamic_cast<MAT::PAR::AAAneohooke*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        // This is a tiny little bit brutal!!!
        const_cast<double&>(params->youngs_) = p_cur[count];
        const_cast<double&>(params->beta_)   = p_cur[count+1];
        //const_cast<double&>(params->nue_)    = p_cur[count+2];
        if (myrank == 0) printf("MAT::PAR::AAAneohooke %20.15e %20.15e\n",p_cur[count],p_cur[count+1]);
        count += 2;
      }
      break;
      case INPAR::MAT::m_neohooke:
      {
        MAT::PAR::NeoHooke* params = dynamic_cast<MAT::PAR::NeoHooke*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        // This is a tiny little bit brutal!!!
        const_cast<double&>(params->youngs_)       = p_cur[count];
        const_cast<double&>(params->poissonratio_) = p_cur[count+1];
        if (myrank == 0) printf("MAT::PAR::NeoHooke %20.15e %20.15e\n",params->youngs_,params->poissonratio_);
        count += 2;
      }
      break;
      case INPAR::MAT::m_elasthyper:
      {
        MAT::PAR::ElastHyper* params = dynamic_cast<MAT::PAR::ElastHyper*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        const int nummat               = params->nummat_;
        const std::vector<int>* matids = params->matids_;
        for (int i=0; i<nummat; ++i)
        {
          const int id = (*matids)[i];
          const RCP<MAT::PAR::Material> actelastmat = mats.find(id)->second;
          switch (actelastmat->Type())
          {
            case INPAR::MAT::mes_isoyeoh:
            {
              MAT::ELASTIC::PAR::IsoYeoh* params = dynamic_cast<MAT::ELASTIC::PAR::IsoYeoh*>(actelastmat->Parameter());
              if (!params) dserror("Cannot cast material parameters");
              const_cast<double&>(params->c1_) = p_cur[count];
              const_cast<double&>(params->c2_) = p_cur[count+1];
              const_cast<double&>(params->c3_) = p_cur[count+2];
              if (myrank == 0) printf("MAT::ELASTIC::PAR::IsoYeoh %20.15e %20.15e %20.15e\n",params->c1_,params->c2_,params->c3_);
              count += 3;
            }
            break;
            case INPAR::MAT::mes_vologden:
            {
              MAT::ELASTIC::PAR::VolOgden* params = dynamic_cast<MAT::ELASTIC::PAR::VolOgden*>(actelastmat->Parameter());
              if (!params) dserror("Cannot cast material parameters");
              const_cast<double&>(params->kappa_) = p_cur[count];
              //const_cast<double&>(params->beta_) = p_cur[count+1];
              if (myrank == 0) printf("MAT::ELASTIC::PAR::VolOgden %20.15e %20.15e\n",params->kappa_,params->beta_);
              count += 1;
            }
            break;
            default:
              dserror("Unknown type of elasthyper material");
            break;
          }
        }
      }
      break;
      case INPAR::MAT::mes_isoyeoh: // at this level do nothing, its inside the INPAR::MAT::m_elasthyper block
      break;
      case INPAR::MAT::mes_vologden: // at this level do nothing, its inside the INPAR::MAT::m_elasthyper block
      break;
      default:
        // ignore unknown materials ?
        dserror("Unknown type of material");
      break;
    }
  }


  return;
}



#endif

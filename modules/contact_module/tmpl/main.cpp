/*---------------------------------------------------------------------------*/
/*!
\file main.cpp

\brief contact module main testing program

\level 2

\maintainer Philipp Farah

*/
/*---------------------------------------------------------------------------*/
//*****************************************************************************
// TEST PROGRAM
//*****************************************************************************
// This test program has no functionality for contact analysis, however it
// mimics the most important steps of the setup of the contact libraries.
//*****************************************************************************

// standard includes
#include <iostream>
#include <vector>

// Epetra (Trilinos)
#include "Epetra_Version.h"
#include "Epetra_Map.h"

// Epetra communicator
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// Teuchos (Trilinos)
#include <Teuchos_Array.hpp>
#include <Teuchos_StrUtils.hpp>
#include <Teuchos_any.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

// contact library includes
#include "drt_contact/friction_node.H"
#include "drt_contact/contact_node.H"
#include "drt_contact/contact_element.H"
#include "drt_contact/contact_interface.H"
#include "drt_contact/contact_penalty_strategy.H"
#include "drt_inpar/inpar_contact.H"
#include "drt_inpar/inpar_wear.H"

#include "validParameters.H"  // provided by contact_parameters module

//-----------------------------------------------------------------------------
// function for creation of contact (friction) nodes
Teuchos::RCP<CONTACT::FriNode> createContactNode(int id,
                                                 double* coords,
                                                 std::vector<int> dofs,
                                                 bool isslave)
{
  Teuchos::RCP<CONTACT::FriNode> n = Teuchos::rcp(new CONTACT::FriNode(id,coords,0,dofs.size(),
                                                                       dofs,isslave,false,false));

  return n;
}


//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
// main function
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

int main(int argc, char *argv[]) {
  std::cout << Epetra_Version() << std::endl << std::endl;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  // create contact nodes (master)
  double coords1[3];  coords1[0] = 0.0;  coords1[1] = 0.0;  coords1[2] = 0.0;
  std::vector<int> dofs1;  dofs1.push_back(0);  dofs1.push_back(1);  dofs1.push_back(2);
  Teuchos::RCP<CONTACT::FriNode> n1 =  createContactNode(0, coords1, dofs1, false);

  double coords2[3];  coords2[0] = 1.0;  coords2[1] = 0.0;  coords2[2] = 0.0;
  std::vector<int> dofs2;  dofs2.push_back(3);  dofs2.push_back(4);  dofs2.push_back(5);
  Teuchos::RCP<CONTACT::FriNode> n2 =  createContactNode(1, coords2, dofs2, false);

  double coords3[3];  coords3[0] = 0.0;  coords3[1] = 1.0;  coords3[2] = 0.0;
  std::vector<int> dofs3;  dofs3.push_back(6);  dofs3.push_back(7);  dofs3.push_back(8);
  Teuchos::RCP<CONTACT::FriNode> n3 =  createContactNode(2, coords3, dofs3, false);

  // create contact element (master)
  std::vector<int> eleids; eleids.push_back(0); eleids.push_back(1); eleids.push_back(2);
  Teuchos::RCP<CONTACT::CoElement> e1 = Teuchos::rcp(new CONTACT::CoElement(0,0,DRT::Element::tri3,eleids.size(),&eleids[0],false));

  // create contact nodes (slave)
  double coords4[3];  coords4[0] = 0.0;  coords4[1] = 0.0;  coords4[2] = 1.0;
  std::vector<int> dofs4;  dofs4.push_back(9);  dofs4.push_back(10);  dofs4.push_back(11);
  Teuchos::RCP<CONTACT::FriNode> n4 =  createContactNode(3, coords4, dofs4, true);

  double coords5[3];  coords5[0] = 0.0;  coords5[1] = 1.0;  coords5[2] = 1.0;
  std::vector<int> dofs5;  dofs5.push_back(12);  dofs5.push_back(13);  dofs5.push_back(14);
  Teuchos::RCP<CONTACT::FriNode> n5 =  createContactNode(4, coords5, dofs5, true);

  double coords6[3];  coords6[0] = 1.0;  coords6[1] = 0.0;  coords6[2] = 1.0;
  std::vector<int> dofs6;  dofs6.push_back(15);  dofs6.push_back(16);  dofs6.push_back(17);
  Teuchos::RCP<CONTACT::FriNode> n6 =  createContactNode(5, coords6, dofs6, true);

   // create contact element (slave)
  std::vector<int> eleids2; eleids2.push_back(3); eleids2.push_back(4); eleids2.push_back(5);
  Teuchos::RCP<CONTACT::CoElement> e2 = Teuchos::rcp(new CONTACT::CoElement(1,0,DRT::Element::tri3,eleids2.size(),&eleids2[0],true));

  // create and fill parameter list
  Teuchos::ParameterList myParams;
  myParams.set<std::string>("FRICTION","Coulomb");
  myParams.set<std::string>("ALGORITHM","Mortar");
  myParams.set<double>("FRCOEFF",0.3);
  myParams.set<double>("PENALTYPARAM",1.0e3);
  myParams.set<double>("PENALTYPARAMTAN",1.0e3);
  myParams.set<std::string>("LM_SHAPEFCN","standard");
  myParams.set<std::string>("STRATEGY","penalty");
  myParams.set<std::string>("LM_NODAL_SCALE","No");
  myParams.set<std::string>("SEARCH_USE_AUX_POS","No");

  // check validity of parameter list
  Teuchos::ParameterList validParams;
  getValidParameters(validParams);
  myParams.validateParametersAndSetDefaults(validParams);
  myParams.set<bool>("NURBS",false); // not an user parameter.

  // add global parameter problem type
  // - unfortunately, we need this little hack due to our code structure, which in the
  //   meantime can also deal with coupled thermo-structure and finite wear problems
  // - this must be added after(!) the validation of the input parameters
  myParams.set<int>("PROBTYPE",INPAR::CONTACT::structure);

  // create a contact interface
  std::cout << "create a Contact interface" << std::endl;
  Teuchos::RCP<CONTACT::CoInterface> interf = Teuchos::rcp(new CONTACT::CoInterface(
                                  0,      /* unique ID */
                                  comm,   /* communicator */
                                  3,      /* 2d only */
                                  myParams, /* parameter list */
                                  false,  /* no self-contact */
                                  INPAR::MORTAR::redundant_master /* not relevant*/
                                  ));

  // add nodes to interface
  interf->AddCoNode(n1);
  interf->AddCoNode(n2);
  interf->AddCoNode(n3);

  interf->AddCoNode(n4);
  interf->AddCoNode(n5);
  interf->AddCoNode(n6);

  // add elements to interface
  interf->AddCoElement(e1);
  interf->AddCoElement(e2);

  // call fill complete on interface
  interf->FillComplete(17);

  // print interface
  interf->Print(std::cout);

  // print friction information
  INPAR::CONTACT::FrictionType fric = DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(myParams,"FRICTION");
  double checkfrcoeff = 0.0;
  if (fric == INPAR::CONTACT::friction_tresca)
  {
    checkfrcoeff = interf->IParams().get<double>("FRBOUND");
    std::cout << "FrBound (Tresca)  " << checkfrcoeff << std::endl;
  }
  else if (fric == INPAR::CONTACT::friction_coulomb)
  {
    checkfrcoeff = interf->IParams().get<double>("FRCOEFF");
    std::cout << "FrCoeff (Coulomb) " << checkfrcoeff << std::endl;
  }

  // create a vector of interfaces
  std::vector<Teuchos::RCP<CONTACT::CoInterface> > myInterfaces;
  myInterfaces.push_back(interf);

  // create map of global degrees of freedom
  Teuchos::RCP<Epetra_Map> myDofMap = Teuchos::rcp(new Epetra_Map(18,0,comm));

   // create map of global nodes
  Teuchos::RCP<Epetra_Map> myNodeMap = Teuchos::rcp(new Epetra_Map(6,0,comm));


  // create a contact penalty strategy
  std::cout << "create a Contact penalty strategy" << std::endl;
  CONTACT::CoPenaltyStrategy penStrat(myDofMap.get(),
                                      myNodeMap.get(),
                                      myParams,
                                      myInterfaces, 3,
                                      Teuchos::rcp(&comm,false), 0.0, 17);

  // create binary search tree
  for (int i=0; i<(int)myInterfaces.size();++i)
    myInterfaces[i]->CreateSearchTree();

  // create zero global displacement vector
  Teuchos::RCP<Epetra_Vector> dis = Teuchos::rcp(new Epetra_Vector(*myDofMap));

  // evaluate contact strategy
  penStrat.SetState(MORTAR::state_new_displacement,*dis);
  penStrat.InitEvalInterface();
  //penStrat.InitEvalMortar();

  // evaluate relative movement for friction
  penStrat.EvaluateRelMov();

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;
}



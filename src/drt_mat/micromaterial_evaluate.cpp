/*----------------------------------------------------------------------*/
/*!
\file micromaterial_evaluate.cpp

\brief class for handling of micro-macro transitions

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>
*/
/*----------------------------------------------------------------------*/


#include "micromaterial.H"
#include "micromaterialgp_static.H"
#include "matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_container.H"
#include "../drt_lib/drt_exporter.H"



// This function has to be separated from the remainder of the
// MicroMaterial class. MicroMaterialGP is NOT a member of
// FILTER_OBJECTS hence the MicroMaterial::Evaluate function that
// builds the connection to MicroMaterialGP is not either. In
// post_drt_evaluation.cpp this function is defined to content the
// compiler. If during postprocessing the MicroMaterial::Evaluate
// function should be called, an error is invoked.
//
// -> see also Makefile.objects and setup-objects.sh
//
// In case of any changes of the function prototype make sure that the
// corresponding prototype in src/filter_common/filter_evaluation.cpp is adapted, too!!

// evaluate for master procs
void MAT::MicroMaterial::Evaluate(LINALG::Matrix<3,3>* defgrd,
                                  LINALG::Matrix<6,6>* cmat,
                                  LINALG::Matrix<6,1>* stress,
                                  double* density,
                                  const int gp,
                                  const int ele_ID)
{
  // activate microscale material

  int microdisnum = MicroDisNum();
  double V0 = InitVol();
  DRT::Problem::Instance()->Materials()->SetReadFromProblem(microdisnum);

  // avoid writing output also for ghosted elements
  const bool eleowner = DRT::Problem::Instance(0)->GetDis("structure")->ElementRowMap()->MyGID(ele_ID);

  // get sub communicator including the supporting procs
  Teuchos::RCP<Epetra_Comm> subcomm = DRT::Problem::Instance(0)->GetNPGroup()->SubComm();

  // tell the supporting procs that the micro material will be evaluated
  int task[2] = {0,ele_ID};
  subcomm->Broadcast(task, 2, 0);

  // container is filled with data for supporting procs
  std::map<int,Teuchos::RCP<DRT::Container> > condnamemap;
  condnamemap[0] = Teuchos::rcp(new DRT::Container());

  condnamemap[0]->Add<3,3>("defgrd",*defgrd);
  condnamemap[0]->Add<6,6>("cmat",*cmat);
  condnamemap[0]->Add<6,1>("stress",*stress);
  condnamemap[0]->Add("density",*density);
  condnamemap[0]->Add("gp",gp);
  condnamemap[0]->Add("microdisnum",microdisnum);
  condnamemap[0]->Add("V0",V0);
  condnamemap[0]->Add("eleowner",eleowner);

	// maps are created and data is broadcast to the supporting procs
  int tag = 0;
  Teuchos::RCP<Epetra_Map> oldmap = Teuchos::rcp(new Epetra_Map(1,1,&tag,0,*subcomm));
  Teuchos::RCP<Epetra_Map> newmap = Teuchos::rcp(new Epetra_Map(1,1,&tag,0,*subcomm));
  DRT::Exporter exporter(*oldmap,*newmap,*subcomm);
  exporter.Export<DRT::Container>(condnamemap);

  // standard evaluation of the micro material
  if (matgp_.find(gp) == matgp_.end())
  {
    matgp_[gp] = Teuchos::rcp(new MicroMaterialGP(gp, ele_ID, eleowner, microdisnum, V0));
  }

  Teuchos::RCP<MicroMaterialGP> actmicromatgp = matgp_[gp];

  // perform microscale simulation and homogenization (if fint and stiff/mass or stress calculation is required)
  actmicromatgp->PerformMicroSimulation(defgrd, stress, cmat, density);

  // reactivate macroscale material
  DRT::Problem::Instance()->Materials()->ResetReadFromProblem();

  return;
}


//evaluate for supporting procs
void MAT::MicroMaterial::Evaluate(LINALG::Matrix<3,3>* defgrd,
                                  LINALG::Matrix<6,6>* cmat,
                                  LINALG::Matrix<6,1>* stress,
                                  double* density,
                                  const int gp,
                                  const int ele_ID,
                                  const int microdisnum,
                                  double V0,
                                  bool eleowner)
{
  DRT::Problem::Instance()->Materials()->SetReadFromProblem(microdisnum);

  if (matgp_.find(gp) == matgp_.end())
  {
    matgp_[gp] = Teuchos::rcp(new MicroMaterialGP(gp, ele_ID, eleowner, microdisnum, V0));
  }

  Teuchos::RCP<MicroMaterialGP> actmicromatgp = matgp_[gp];

  // perform microscale simulation and homogenization (if fint and stiff/mass or stress calculation is required)
  actmicromatgp->PerformMicroSimulation(defgrd, stress, cmat, density);

  // reactivate macroscale material
  DRT::Problem::Instance()->Materials()->ResetReadFromProblem();

  return;
}


// update for all procs
void MAT::MicroMaterial::Update()
{
  // get sub communicator including the supporting procs
  Teuchos::RCP<Epetra_Comm> subcomm = DRT::Problem::Instance(0)->GetNPGroup()->SubComm();
  if(subcomm->MyPID() == 0)
  {
    // tell the supporting procs that the micro material will be evaluated for the element with id eleID
    int eleID = matgp_.begin()->second->eleID();
    int task[2] = {2,eleID};
    subcomm->Broadcast(task, 2, 0);
  }

  std::map<int, Teuchos::RCP<MicroMaterialGP> >::iterator it;
  for (it=matgp_.begin(); it!=matgp_.end(); ++it)
  {
    Teuchos::RCP<MicroMaterialGP> actmicromatgp = (*it).second;
    actmicromatgp->Update();
  }
}


// prepare output for all procs
void MAT::MicroMaterial::PrepareOutput()
{
  // get sub communicator including the supporting procs
  Teuchos::RCP<Epetra_Comm> subcomm = DRT::Problem::Instance(0)->GetNPGroup()->SubComm();
  if(subcomm->MyPID() == 0)
  {
    // tell the supporting procs that the micro material will be prepared for output
    int eleID = matgp_.begin()->second->eleID();
    int task[2] = {1,eleID};
    subcomm->Broadcast(task, 2, 0);
  }

  std::map<int, Teuchos::RCP<MicroMaterialGP> >::iterator it;
  for (it=matgp_.begin(); it!=matgp_.end(); ++it)
  {
    Teuchos::RCP<MicroMaterialGP> actmicromatgp = (*it).second;
    actmicromatgp->PrepareOutput();
  }
}


// output for all procs
void MAT::MicroMaterial::Output()
{
  // get sub communicator including the supporting procs
  Teuchos::RCP<Epetra_Comm> subcomm = DRT::Problem::Instance(0)->GetNPGroup()->SubComm();
  if(subcomm->MyPID() == 0)
  {
    // tell the supporting procs that the micro material will be output
    int eleID = matgp_.begin()->second->eleID();
    int task[2] = {3,eleID};
    subcomm->Broadcast(task, 2, 0);
  }

  std::map<int, Teuchos::RCP<MicroMaterialGP> >::iterator it;
  for (it=matgp_.begin(); it!=matgp_.end(); ++it)
  {
    Teuchos::RCP<MicroMaterialGP> actmicromatgp = (*it).second;
    actmicromatgp->Output();
  }
}


// read restart for master procs
void MAT::MicroMaterial::ReadRestart(const int gp, const int eleID, const bool eleowner)
{
  int microdisnum = MicroDisNum();
  double V0 = InitVol();

  // get sub communicator including the supporting procs
  Teuchos::RCP<Epetra_Comm> subcomm = DRT::Problem::Instance(0)->GetNPGroup()->SubComm();

  // tell the supporting procs that the micro material will restart
  int task[2] = {4,eleID};
  subcomm->Broadcast(task, 2, 0);

  // container is filled with data for supporting procs
  std::map<int,Teuchos::RCP<DRT::Container> > condnamemap;
  condnamemap[0] = Teuchos::rcp(new DRT::Container());

  condnamemap[0]->Add("gp",gp);
  condnamemap[0]->Add("microdisnum",microdisnum);
  condnamemap[0]->Add("V0",V0);
  condnamemap[0]->Add("eleowner",eleowner);

  // maps are created and data is broadcast to the supporting procs
  int tag = 0;
  Teuchos::RCP<Epetra_Map> oldmap = Teuchos::rcp(new Epetra_Map(1,1,&tag,0,*subcomm));
  Teuchos::RCP<Epetra_Map> newmap = Teuchos::rcp(new Epetra_Map(1,1,&tag,0,*subcomm));
  DRT::Exporter exporter(*oldmap,*newmap,*subcomm);
  exporter.Export<DRT::Container>(condnamemap);

  if (matgp_.find(gp) == matgp_.end())
  {
    matgp_[gp] = Teuchos::rcp(new MicroMaterialGP(gp, eleID, eleowner, microdisnum, V0));
  }

  Teuchos::RCP<MicroMaterialGP> actmicromatgp = matgp_[gp];
  actmicromatgp->ReadRestart();
}


// read restart for supporting procs
void MAT::MicroMaterial::ReadRestart(const int gp, const int eleID, const bool eleowner, int microdisnum, double V0)
{
  if (matgp_.find(gp) == matgp_.end())
  {
    matgp_[gp] = Teuchos::rcp(new MicroMaterialGP(gp, eleID, eleowner, microdisnum, V0));
  }

  Teuchos::RCP<MicroMaterialGP> actmicromatgp = matgp_[gp];
  actmicromatgp->ReadRestart();
}


void MAT::MicroMaterial::InvAnaInit(bool eleowner, int eleID)
{
  // get sub communicator including the supporting procs
  Teuchos::RCP<Epetra_Comm> subcomm = DRT::Problem::Instance(0)->GetNPGroup()->SubComm();
  if(subcomm->MyPID() == 0)
  {
    // tell the supporting procs that the micro material initializes inverse analysis
    int task[2] = {5,eleID};
    subcomm->Broadcast(task, 2, 0);
    int owner = eleowner;
    subcomm->Broadcast(&owner, 1, 0);
  }

  std::map<int, Teuchos::RCP<MicroMaterialGP> >::iterator it;
  for (it=matgp_.begin(); it!=matgp_.end(); ++it)
  {
    Teuchos::RCP<MicroMaterialGP> actmicromatgp = (*it).second;
    actmicromatgp->ResetTimeAndStep();
    std::string newfilename;
    actmicromatgp->NewResultFile(eleowner, newfilename);
  }
}

/*---------------------------------------------------------------------*/
/*! \file

\brief Static control for  microstructural problems in case of multiscale
analysis for supporting processors


\level 3

*/
/*---------------------------------------------------------------------*/

#include "stru_multi_microstatic_npsupport.H"
#include "lib_globalproblem.H"
#include "comm_utils.H"
#include "mat_micromaterial.H"
#include "linalg_fixedsizematrix.H"
#include "lib_container.H"
#include "lib_exporter.H"
#include "linalg_utils_densematrix_communication.H"

#include <hdf5.h>


/*----------------------------------------------------------------------*
 | "timeloop" for supporting procs                          ghamm 05/12 |
 *----------------------------------------------------------------------*/
void STRUMULTI::np_support_drt()
{
  // info:
  // Macro processors run on their macro problem (distributed or not)
  // and during run time supporting processors are available.
  // Dependent on the number of processors on the macro scale that need
  // support, the supporting procs are distributed. The distribution is performed
  // in DRT::Problem::ReadMicroFields() where a sub communicator is created
  // that contains one master proc and several supporting procs. The master
  // proc has always MyPID() = 0 in the subcomm. That's why all broadcast commands
  // send from proc 0 in subcomm.
  // Supporting procs always wait in front of "subcomm->Broadcast(task, 2, 0);"
  // until a master proc reaches the same broadcast command. Then all procs
  // in this subcomm go into the same routine. We need one dummy material that
  // is specified in the input file of the supporting procs. It is only used to
  // call the corresponding routines. All necessary information comes from the
  // master proc. Supporting processors always start their work in the routines
  // in this file.

  // one dummy micro material is specified in the supporting input file
  std::map<int, Teuchos::RCP<MAT::MicroMaterial>> dummymaterials;

  // this call is needed in order to increment the unique ids that are distributed
  // by HDF5; the macro procs call output->WriteMesh(0, 0.0) in ADAPTER::Structure
  H5Fcreate("xxxdummyHDF5file", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // get sub communicator including the master proc
  Teuchos::RCP<Epetra_Comm> subcomm = DRT::Problem::Instance(0)->GetCommunicators()->SubComm();

  // bool for checking whether restart has already been called
  bool restart = false;

  // start std::endless loop
  while (true)
  {
    // receive what to do and which element is intended to work
    std::array<int, 2> task = {-1, -1};
    subcomm->Broadcast(task.data(), 2, 0);
    int whattodo = task[0];
    int eleID = task[1];

    // every element needs one micromaterial
    if (dummymaterials[eleID] == Teuchos::null)
      dummymaterials[eleID] =
          Teuchos::rcp_static_cast<MAT::MicroMaterial>(MAT::Material::Factory(1));

    // check what is the next task of the supporting procs
    switch (whattodo)
    {
      case 0:
      {
        // receive data from the master proc
        int tag = 0;
        Teuchos::RCP<Epetra_Map> oldmap = Teuchos::rcp(new Epetra_Map(1, 0, &tag, 0, *subcomm));
        Teuchos::RCP<Epetra_Map> newmap = Teuchos::rcp(new Epetra_Map(1, 1, &tag, 0, *subcomm));
        // create an exporter object that will figure out the communication pattern
        DRT::Exporter exporter(*oldmap, *newmap, *subcomm);
        std::map<int, Teuchos::RCP<DRT::Container>> condnamemap;
        exporter.Export<DRT::Container>(condnamemap);

        // extract received data from the container
        const Epetra_SerialDenseMatrix* defgrdcopy =
            condnamemap[0]->Get<Epetra_SerialDenseMatrix>("defgrd");
        CORE::LINALG::Matrix<3, 3> defgrd(
            *(const_cast<Epetra_SerialDenseMatrix*>(defgrdcopy)), true);
        const Epetra_SerialDenseMatrix* cmatcopy =
            condnamemap[0]->Get<Epetra_SerialDenseMatrix>("cmat");
        CORE::LINALG::Matrix<6, 6> cmat(*(const_cast<Epetra_SerialDenseMatrix*>(cmatcopy)), true);
        const Epetra_SerialDenseMatrix* stresscopy =
            condnamemap[0]->Get<Epetra_SerialDenseMatrix>("stress");
        CORE::LINALG::Matrix<6, 1> stress(
            *(const_cast<Epetra_SerialDenseMatrix*>(stresscopy)), true);
        int gp = condnamemap[0]->GetInt("gp");
        int microdisnum = condnamemap[0]->GetInt("microdisnum");
        double V0 = condnamemap[0]->GetDouble("V0");
        bool eleowner = condnamemap[0]->GetInt("eleowner");

        // dummy material is used to evaluate the micro material
        dummymaterials[eleID]->Evaluate(
            &defgrd, &cmat, &stress, gp, eleID, microdisnum, V0, eleowner);
        break;
      }
      case 1:
      {
        // dummy material is used to prepare the output of the micro material
        dummymaterials[eleID]->PrepareOutput();
        break;
      }
      case 2:
      {
        // dummy material is used to update the micro material
        dummymaterials[eleID]->Update();
        break;
      }
      case 3:
      {
        // dummy material is used to output the micro material
        dummymaterials[eleID]->Output();
        break;
      }
      case 4:
      {
        // receive data from the master proc for restart
        int tag = 0;
        Teuchos::RCP<Epetra_Map> oldmap = Teuchos::rcp(new Epetra_Map(1, 0, &tag, 0, *subcomm));
        Teuchos::RCP<Epetra_Map> newmap = Teuchos::rcp(new Epetra_Map(1, 1, &tag, 0, *subcomm));
        // create an exporter object that will figure out the communication pattern
        DRT::Exporter exporter(*oldmap, *newmap, *subcomm);
        std::map<int, Teuchos::RCP<DRT::Container>> condnamemap;
        exporter.Export<DRT::Container>(condnamemap);

        // extract received data from the container
        int gp = condnamemap[0]->GetInt("gp");
        int microdisnum = condnamemap[0]->GetInt("microdisnum");
        double V0 = condnamemap[0]->GetDouble("V0");
        bool eleowner = condnamemap[0]->GetInt("eleowner");

        // the material is replaced in read mesh within restart on macro scale; this is done here
        // manually one-time
        if (restart == false)
        {
          restart = true;
          dummymaterials.clear();
        }

        // new dummy material is created if necessary
        if (dummymaterials[eleID] == Teuchos::null)
          dummymaterials[eleID] =
              Teuchos::rcp_static_cast<MAT::MicroMaterial>(MAT::Material::Factory(1));

        // dummy material is used to restart the micro material
        dummymaterials[eleID]->ReadRestart(gp, eleID, eleowner, microdisnum, V0);
        break;
      }
      case 5:
      {
        // receive data from the master proc for inverse analysis
        int owner = -1;
        subcomm->Broadcast(&owner, 1, 0);
        const bool eleowner = owner;
        // dummy material is used initialize the inverse analysis on the micro material
        dummymaterials[eleID]->InvAnaInit(eleowner, eleID);
        break;
      }
      case 9:
      {
        // end of simulation after deleting all dummy materials
        dummymaterials.clear();
        return;
      }
      default:
        dserror("Supporting processors do not know what to do (%i)!", whattodo);
        break;
    }

  }  // end while(true) loop

  return;
}

/*----------------------------------------------------------------------*/
/*!
 \file porofluidmultiphase_utils.cpp

 \brief helper function/class for multiphase porous flow problems

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de1
 *----------------------------------------------------------------------*/

#include "porofluidmultiphase_utils.H"

#include "porofluidmultiphase_timint_ost.H"

#include "../drt_mat/material.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"

#include "../drt_adapter/ad_porofluidmultiphase.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::UTILS::SetupMaterial(
    const Epetra_Comm& comm,
    const std::string& struct_disname,
    const std::string& fluid_disname)
{
  // get the fluid discretization
  Teuchos::RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->GetDis(fluid_disname);

  // initialize material map
  std::map<int,int> matmap;
  {
    // get the cloning material map from the .dat file
    std::map<std::pair<std::string,std::string>,std::map<int,int> > clonefieldmatmap =
        DRT::Problem::Instance()->CloningMaterialMap();
    if (clonefieldmatmap.size() < 1)
      dserror("At least one material pairing required in --CLONING MATERIAL MAP.");

    // check if the current discretization is included in the material map
    std::pair<std::string,std::string> key(fluid_disname,struct_disname);
    matmap = clonefieldmatmap[key];
    if (matmap.size() < 1)
      dserror("Key pair '%s/%s' not defined in --CLONING MATERIAL MAP.",
          fluid_disname.c_str(),struct_disname.c_str());
  }


  // number of column elements within fluid discretization
  const int numelements = fluiddis->NumMyColElements();

  // loop over column elements
  for (int i=0; i<numelements; ++i)
  {
    // get current element
    DRT::Element* ele = fluiddis->lColElement(i);

    // find the corresponding material in the matmap
    int src_matid = ele->Material()->Parameter()->Id();
    std::map<int,int>::iterator mat_iter = matmap.find(src_matid);
    if (mat_iter!=matmap.end())
    {
      // get the ID of the secondary material
      const int tar_matid = mat_iter->second;
      // build the material usilng the factory
      Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(tar_matid);

      // add secondary material to poro fluid element
      if(ele->AddMaterial(mat)!=2)
        dserror("unexpected number of materials!");
    }
    else
    {
      // before we stop, print the material id map
      std::cout<<"Material map on PROC "<<comm.MyPID()<<":"<<std::endl;
      for(mat_iter=matmap.begin(); mat_iter != matmap.end(); mat_iter++)
        std::cout<<mat_iter->first<<" -> "<<mat_iter->second<<std::endl;

      dserror("no matching material ID (%d) in map",src_matid);
    }

  } // end loop over column elements

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> POROFLUIDMULTIPHASE::UTILS::ConvertDofVectorToNodeBasedMultiVector(
    const DRT::Discretization& dis,
    const Epetra_Vector& vector,
    const int nds,
    const int numdofpernode)
{
  // initialize multi vector
  Teuchos::RCP<Epetra_MultiVector> multi =
      Teuchos::rcp(new Epetra_MultiVector(*dis.NodeRowMap(),numdofpernode,true));

  // get maps
  const Epetra_BlockMap& vectormap = vector.Map();

  // loop over nodes of the discretization
  for (int inode=0; inode<dis.NumMyRowNodes(); ++inode)
  {
    // get current node
    DRT::Node* node = dis.lRowNode(inode);
    //copy each dof value of node
    for (int idof=0; idof<numdofpernode; ++idof)
      (*multi)[idof][inode] = vector[vectormap.LID(dis.Dof(nds,node,idof))];
  }

  return multi;
}

/*----------------------------------------------------------------------*
 | create algorithm                                                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<ADAPTER::PoroFluidMultiphase> POROFLUIDMULTIPHASE::UTILS::CreateAlgorithm(
    INPAR::POROFLUIDMULTIPHASE::TimeIntegrationScheme timintscheme,
    Teuchos::RCP<DRT::Discretization>                 dis,
    const int                                         linsolvernumber,
    const Teuchos::ParameterList&                     probparams,
    const Teuchos::ParameterList&                     poroparams,
    FILE*                                             errfile,
    Teuchos::RCP<IO::DiscretizationWriter>            output
    )
{
  // Creation of Coupled Problem algortihm.
  Teuchos::RCP<ADAPTER::PoroFluidMultiphase> algo = Teuchos::null;

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------

  switch(timintscheme)
  {
  case INPAR::POROFLUIDMULTIPHASE::timeint_one_step_theta:
  {
    // create algorithm
    algo = Teuchos::rcp(new POROFLUIDMULTIPHASE::TimIntOneStepTheta(
            dis,
            linsolvernumber,
            probparams,
            poroparams,
            errfile,
            output));
    break;
  }
  default:
    dserror("Unknown time-integration scheme for multiphase poro fluid problem");
    break;
  }

  return algo;
}

/*----------------------------------------------------------------------*
 | calculate vector norm                             kremheller 12/17   |
 *----------------------------------------------------------------------*/
double POROFLUIDMULTIPHASE::UTILS::CalculateVectorNorm(
  const enum INPAR::POROFLUIDMULTIPHASE::VectorNorm norm,
  const Teuchos::RCP<const Epetra_Vector> vect
  )
{
  // L1 norm
  // norm = sum_0^i vect[i]
  if (norm == INPAR::POROFLUIDMULTIPHASE::norm_l1)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm;
  }
  // L2/Euclidian norm
  // norm = sqrt{sum_0^i vect[i]^2 }
  else if (norm == INPAR::POROFLUIDMULTIPHASE::norm_l2)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm;
  }
  // RMS norm
  // norm = sqrt{sum_0^i vect[i]^2 }/ sqrt{length_vect}
  else if (norm == INPAR::POROFLUIDMULTIPHASE::norm_rms)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm/sqrt((double) vect->GlobalLength());
  }
  // infinity/maximum norm
  // norm = max( vect[i] )
  else if (norm == INPAR::POROFLUIDMULTIPHASE::norm_inf)
  {
    double vectnorm;
    vect->NormInf(&vectnorm);
    return vectnorm;
  }
  // norm = sum_0^i vect[i]/length_vect
  else if (norm == INPAR::POROFLUIDMULTIPHASE::norm_l1_scaled)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm/((double) vect->GlobalLength());
  }
  else
  {
    dserror("Cannot handle vector norm");
    return 0;
  }
}  // CalculateVectorNorm()

/*----------------------------------------------------------------------*
 |                                                    kremheller 03/17  |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::PrintLogo()
{
 std::cout << "This is a Porous Media problem with multiphase flow" << std::endl;
 std::cout << "       .--..--..--..--..--..--. " << std::endl;
 std::cout << "      .'  \\  (`._   (_)     _   \\ " << std::endl;
 std::cout << "     .'    |  '._)         (_)  | " << std::endl;
 std::cout << "     \\ _.')\\      .----..---.   / " << std::endl;
 std::cout << "     |(_.'  |    /    .-\\-.  \\  | " << std::endl;
 std::cout << "     \\     0|    |   ( O| O) | o| " << std::endl;
 std::cout << "      |  _  |  .--.____.'._.-.  | " << std::endl;
 std::cout << "      \\ (_) | o         -` .-`  | " << std::endl;
 std::cout << "       |    \\   |`-._ _ _ _ _\\ / " << std::endl;
 std::cout << "       \\    |   |  `. |_||_|   | " << std::endl;
 std::cout << "       | o  |    \\_      \\     |     -.   .-. " << std::endl;
 std::cout << "       |.-.  \\     `--..-'   O |     `.`-' .' " << std::endl;
 std::cout << "     _.'  .' |     `-.-'      /-.__   ' .-' " << std::endl;
 std::cout << "   .' `-.` '.|='=.='=.='=.='=|._/_ `-'.' " << std::endl;
 std::cout << "   `-._  `.  |________/\\_____|    `-.' " << std::endl;
 std::cout << "      .'   ).| '=' '='\\/ '=' | " << std::endl;
 std::cout << "      `._.`  '---------------' " << std::endl;
 std::cout << "            //___\\   //___\\ " << std::endl;
 std::cout << "              ||       || " << std::endl;
 std::cout << "              ||_.-.   ||_.-. " << std::endl;
 std::cout << "             (_.--__) (_.--__) " << std::endl;
  return;


}

/*----------------------------------------------------------------------*/
/*!
 \file porofluidmultiphase_utils.cpp

 \brief helper function/class for multiphase porous flow problems

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#include "porofluidmultiphase_utils.H"

#include "../drt_mat/material.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"

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


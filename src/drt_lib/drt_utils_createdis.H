/*----------------------------------------------------------------------*/
/*!
\file drt_utils_createdis.H

\brief Declaration utility functions for automatic creation of a discretization
       from an existing one (e.g. ALE from Fluid)

<pre>
\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef DRT_UTILS_CREATEDIS_H
#define DRT_UTILS_CREATEDIS_H

#include "drt_condition_utils.H"
#include "drt_utils_factory.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_immersed_problem/immersed_node.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matpar_parameter.H"
#include "../drt_io/io_pstream.H"

// required for utility function CloneDiscretization()
#include "drt_globalproblem.H"
#include <Epetra_Time.h>

namespace DRT
{

namespace INPUT
{
  class LineDefinition;
}


namespace UTILS
{

/// class providing basic functionality for cloning discretizations
class DiscretizationCreatorBase
{
public:

  /// constructor
  explicit DiscretizationCreatorBase():
  numeleskips_(0)
  {}

  /// destructor
  virtual ~DiscretizationCreatorBase() {}

  void CopyConditions(
      const Teuchos::RCP<DRT::Discretization> sourcedis,
      Teuchos::RCP<DRT::Discretization> targetdis,
      const std::map<std::string,std::string>& conditions_to_copy);

protected:

  void CreateNodes(
      Teuchos::RCP<DRT::Discretization> sourcedis,
      Teuchos::RCP<DRT::Discretization> targetdis,
      std::set<int>& rownodeset,
      std::set<int>& colnodeset,
      const bool isnurbsdis,
      const bool buildimmersednode);

  Teuchos::RCP<Epetra_Map> CreateNodeRowMap(
      std::set<int>& rownodeset,Teuchos::RCP<DRT::Discretization> targetdis);

  Teuchos::RCP<Epetra_Map> CreateNodeColMap(
      std::set<int>& colnodeset,Teuchos::RCP<DRT::Discretization> targetdis);

  void InitialChecks(
      Teuchos::RCP<DRT::Discretization> sourcedis,
      Teuchos::RCP<DRT::Discretization> targetdis);

  void Finalize(
      const Teuchos::RCP<DRT::Discretization> sourcedis,
      Teuchos::RCP<DRT::Discretization> targetdis);


protected:

  Teuchos::RCP<Epetra_Map> targetnoderowmap_;
  Teuchos::RCP<Epetra_Map> targetnodecolmap_;
  //! local number of skipped elements during cloning
  int numeleskips_;

private:


};

/// class for cloning a new discretization from an existing one
template <class CloneStrategy>
class DiscretizationCreator
 : DiscretizationCreatorBase, CloneStrategy
{
public:

  /// constructor
  explicit DiscretizationCreator() {}
  /// destructor
  virtual ~DiscretizationCreator() {}

  /// method for cloning a new discretization from an existing one
  void CreateMatchingDiscretization(
      Teuchos::RCP<DRT::Discretization> sourcedis,  ///< Teuchos::RCP to source discretization
      Teuchos::RCP<DRT::Discretization> targetdis, ///< Teuchos::RCP to empty target discretization
      const int matid ///< ID of the material which generated elements will get
      )
  {
    // obsolete function call which should not be used anymore!
    // let's use the more general version using an explicit material id mapping
    // as defined in the input file section "--CLONING MATERIAL MAP"

    // check and analyze source discretization (sorcedis must be filled!)
    InitialChecks(sourcedis,targetdis);

    // We have to find out all the material ids of the source discretization.
    // All cloned elements will receive the same material with the provided matid.
    std::map<int,int> matmap;
    int numelements = sourcedis->NumMyColElements();
    if (numelements < 1) dserror("At least one processor has no element");
    for (int i=0; i<numelements; ++i)
    {
      DRT::Element* sourceele = sourcedis->lColElement(i);
      int src_matid = sourceele->Material()->Parameter()->Id();
      // if a new material id is found -> extend the map
      std::map<int,int>::iterator mat_iter = matmap.find(src_matid);
      if (mat_iter==matmap.end())
      {
        std::pair<int,int> matmappair(src_matid,matid);
        matmap.insert(matmappair);
      }
    }

    CreateMatchingDiscretization(sourcedis,targetdis,matmap);

    return;
  };

  /// method for cloning a new discretization from an existing one
  void CreateMatchingDiscretization(
      Teuchos::RCP<DRT::Discretization> sourcedis,  ///< Teuchos::RCP to source discretization
      Teuchos::RCP<DRT::Discretization> targetdis, ///< Teuchos::RCP to empty target discretization
      const std::map<int,int> matmap ///< map of material IDs (source element -> target element)
      )
  {
    // try to cast sourcedis to NurbsDiscretisation
    DRT::NURBS::NurbsDiscretization* nurbsdis =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*(sourcedis)));
    bool isnurbsdis(nurbsdis!=NULL);

    // try to cast source node to immersed node
    IMMERSED::ImmersedNode* inode =
      dynamic_cast<IMMERSED::ImmersedNode*>(sourcedis->gNode(sourcedis->NodeRowMap()->GID(0)));
    bool buildimmersednode(inode!=NULL);

    // check and analyze source discretization
    InitialChecks(sourcedis,targetdis);
    AnalyzeSourceDis(sourcedis,eletype_,egid_,rownodeset_,colnodeset_);

    // do the node business
    CreateNodes(sourcedis,targetdis,rownodeset_,colnodeset_,isnurbsdis,buildimmersednode);
    targetnoderowmap_ = CreateNodeRowMap(rownodeset_,targetdis);
    targetnodecolmap_ = CreateNodeColMap(colnodeset_,targetdis);

    // create elements
    CreateElements(sourcedis,targetdis,matmap,isnurbsdis);

    // copy desired conditions from source to target discretization
    std::map<std::string,std::string> conditions_to_copy = CloneStrategy::ConditionsToCopy();
    CopyConditions(sourcedis,targetdis,conditions_to_copy);

    // call Redistribute, FillComplete etc.
    Finalize(sourcedis,targetdis);
  };

private:

  /// get element type std::strings and global id's and nodes from source discretization
  void AnalyzeSourceDis(
      Teuchos::RCP<DRT::Discretization> sourcedis,
      std::vector<std::string>& eletype,
      std::vector<int>& egid,
      std::set<int>& rownodeset,
      std::set<int>& colnodeset
      )
  {
    const Epetra_Map* noderowmap = sourcedis->NodeRowMap();

    // We need to test for all elements (including ghosted ones) to
    // catch all nodes attached to the elements of the source discretization
    // we will clone only those (-> support for ALE sub-meshes)
    int numelements = sourcedis->NumMyColElements();

    for (int i=0; i<numelements; ++i)
    {
      DRT::Element* actele = sourcedis->lColElement(i);
      bool ismyele = sourcedis->ElementRowMap()->MyGID(actele->Id());

      // we get the element type std::string and a boolean if this element
      // is considered! (see submeshes for Fluid-ALE case!)
      if(CloneStrategy::DetermineEleType(actele, ismyele, eletype_))
      {
        if (ismyele)
          egid.push_back(actele->Id());

        // copy node ids of actele to rownodeset but leave those that do
        // not belong to this processor
        remove_copy_if(actele->NodeIds(), actele->NodeIds()+actele->NumNode(),
                       inserter(rownodeset, rownodeset.begin()),
                       not1(DRT::UTILS::MyGID(noderowmap)));

        copy(actele->NodeIds(), actele->NodeIds()+actele->NumNode(),
             inserter(colnodeset, colnodeset.begin()));
      }
      else
        numeleskips_++;
    } // loop over my elements

    return;
  }

  /// create new elements and add them to the target discretization
  void CreateElements(
      Teuchos::RCP<DRT::Discretization> sourcedis,
      Teuchos::RCP<DRT::Discretization> targetdis,
      std::map<int,int> matmap,
      const bool isnurbsdis)
  {
    // now do the elements
    for (std::map<int,int>::iterator mapit=matmap.begin(); mapit!=matmap.end(); ++mapit)
    {
      int target_id = mapit->second;
      CloneStrategy::CheckMaterialType(target_id);
    }

    // prepare some variables we need
    int myrank = targetdis->Comm().MyPID();

    // construct new elements
    // The order of the elements might be different from that of the
    // source elements. We don't care. There are no dofs to these elements.
    for (std::size_t i=0; i<egid_.size(); ++i)
    {
      DRT::Element* sourceele = sourcedis->gElement(egid_[i]);

      std::string approxtype = "Polynomial";
      if(isnurbsdis)
      {
        if(sourceele->NumNode()==8)
        {
          approxtype="NURBS8";
        }
        else if(sourceele->NumNode()==9)
        {
          approxtype="NURBS9";
        }
        else if(sourceele->NumNode()==4)
        {
          approxtype="NURBS4";
        }
        else if(sourceele->NumNode()==27)
        {
          approxtype="NURBS27";
        }
        else if(sourceele->NumNode()==2)
        {
          approxtype="NURBS2";
        }
        else if(sourceele->NumNode()==3)
        {
          approxtype="NURBS3";
        }
        else
        {
          dserror("unknown type of nurbs element\n");
        }
      }

      // create a new element of desired type with the same global element id
      Teuchos::RCP<DRT::Element> newele = DRT::UTILS::Factory(eletype_[i],approxtype,egid_[i],myrank);

      // get global node ids of fluid element
      std::vector<int> nids;
      nids.reserve(sourceele->NumNode());
      transform(sourceele->Nodes(), sourceele->Nodes()+sourceele->NumNode(),
                back_inserter(nids), std::mem_fun(&DRT::Node::Id));

      // set the same global node ids to the new element
      newele->SetNodeIds(nids.size(), &nids[0]);

      // We need to set material and gauss points to complete element setup.
      // This is again really ugly as we have to extract the actual
      // element type in order to access the material property
      // note: SetMaterial() was reimplemented by the transport element!

      int src_matid = sourceele->Material()->Parameter()->Id();
      std::map<int,int>::iterator mat_iter = matmap.find(src_matid);
      if (mat_iter!=matmap.end())
      {
        int tar_matid = mat_iter->second;
        CloneStrategy::SetElementData(newele,sourceele,tar_matid,isnurbsdis);

        // add new element to discretization
        targetdis->AddElement(newele);
      }
      else
      {
        // before we stop, print the material id map
        std::cout<<"Material map on PROC "<<myrank<<":"<<std::endl;
        for(mat_iter=matmap.begin(); mat_iter != matmap.end(); mat_iter++)
          std::cout<<mat_iter->first<<" -> "<<mat_iter->second<<std::endl;

        dserror("no matching material ID (%d) in map",src_matid);
      }
    }
  }

  //! set of row nodes
  std::set<int> rownodeset_;
  //! set of column nodes
  std::set<int> colnodeset_;
  //! vector for holding global element ids
  std::vector<int> egid_;
  //! vector for holding each (desired) element type std::string
  std::vector<std::string> eletype_;

}; // class DiscretizationCreator


 /// clone target discretization from a given source discretization
 template <class CloneStrategy>
 void CloneDiscretization(
     Teuchos::RCP<DRT::Discretization> sourcedis, ///< source discretization
     Teuchos::RCP<DRT::Discretization> targetdis ///< target discretization
 )
 {
   // access the communicator for time measurement
   const Epetra_Comm& comm = sourcedis->Comm();
   Epetra_Time time(comm);

   // create target discretization using a given clone strategy
   {
     Teuchos::RCP<DRT::UTILS::DiscretizationCreator<CloneStrategy> > clonewizard =
         Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<CloneStrategy>() );

     std::map<std::pair<std::string,std::string>,std::map<int,int> > clonefieldmatmap
     = DRT::Problem::Instance()->CloningMaterialMap();
     if (clonefieldmatmap.size() < 1)
       dserror("At least one material pairing required in --CLONING MATERIAL MAP.");

     std::pair<std::string,std::string> key(sourcedis->Name(),targetdis->Name());
     std::map<int,int> matmap = clonefieldmatmap[key];
     if (matmap.size() < 1)
       dserror("Key pair '%s/%s' not defined in --CLONING MATERIAL MAP.",sourcedis->Name().c_str(),targetdis->Name().c_str());

     clonewizard->CreateMatchingDiscretization(sourcedis,targetdis,matmap);
   }
   if (comm.MyPID()==0)
   {
     IO::cout << "Created discretization " << (targetdis->Name())
              << " as a clone of discretization " << sourcedis->Name()
              << " in...." << time.ElapsedTime() << " secs\n\n";
   }

   return;
 }


 Teuchos::RCP<DRT::INPUT::Lines> ValidCloningMaterialMapLines();

 void PrintCloningMaterialMapDatHeader();


} // namespace UTILS
} // namespace DRT

#endif // #ifndef DRT_UTILS_CREATEDIS_H

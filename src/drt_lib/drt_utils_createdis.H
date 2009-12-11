/*----------------------------------------------------------------------*/
/*!
\file drt_utils_createdis.H

\brief utility functions for automatic creation of a discretization
       from an existing one (e.g. ALE from Fluid)

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef DRT_UTILS_CREATEDIS_H
#define DRT_UTILS_CREATEDIS_H

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_utils.H"

namespace DRT
{
namespace UTILS
{

/// class providing basic functionality for cloning discretizations
class DiscretizationCreatorBase
{
public:

  /// constructor
  explicit DiscretizationCreatorBase(){}

  /// destructor
  virtual ~DiscretizationCreatorBase() {}

protected:

  void CreateNodes(
      Teuchos::RCP<DRT::Discretization> sourcedis,
      Teuchos::RCP<DRT::Discretization> targetdis,
      set<int>& rownodeset,
      set<int>& colnodeset,
      RefCountPtr<Epetra_Map>& targetnoderowmap_,
      RefCountPtr<Epetra_Map>& targetnodecolmap_);

  void InitialChecks(
      Teuchos::RCP<DRT::Discretization> sourcedis,
      Teuchos::RCP<DRT::Discretization> targetdis);

  void CopyConditions(
      const Teuchos::RCP<DRT::Discretization> sourcedis,
      Teuchos::RCP<DRT::Discretization> targetdis,
      const map<string,string>& conditions_to_copy);

  void Finalize(Teuchos::RCP<DRT::Discretization> targetdis);


protected:

  RefCountPtr<Epetra_Map> targetnoderowmap_;
  RefCountPtr<Epetra_Map> targetnodecolmap_;

private:


};

/// class for for cloning a new discretization from an existing one
template <class CloneStrategy>
class DiscretizationCreator
 : DiscretizationCreatorBase, CloneStrategy
{
public:

  /// constructor
  explicit DiscretizationCreator(){}
  /// destructor
  virtual ~DiscretizationCreator() {}

  /// method for cloning a new discretization from an existing one
  void CreateMatchingDiscretization(
      Teuchos::RCP<DRT::Discretization> sourcedis,  ///< RCP to source discretization
      Teuchos::RCP<DRT::Discretization> targetdis, ///< RCP to empty target discretization
      const int matid ///< ID of the material which generated elements will get
      )
  {
    InitialChecks(sourcedis,targetdis);

    AnalyzeSourceDis(sourcedis,eletype_,egid_,rownodeset_,colnodeset_);

    CreateNodes(sourcedis,targetdis,rownodeset_,colnodeset_,targetnoderowmap_,targetnodecolmap_);

    CreateElements(sourcedis,targetdis,matid);

    map<string,string> conditions_to_copy = CloneStrategy::ConditionToCopy();

    CopyConditions(sourcedis,targetdis,conditions_to_copy);

    Finalize(targetdis);
  };

private:

  /// get element type strings and global id's and nodes from source discretization
  void AnalyzeSourceDis(
      Teuchos::RefCountPtr<DRT::Discretization> sourcedis,
      vector<string>& eletype,
      vector<int>& egid,
      set<int>& rownodeset,
      set<int>& colnodeset
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

      // we get the element type string and a boolean if this element
      // is considered! (see submeshes for Fluid-ALE case!)
      if(CloneStrategy::DetermineEleType(eletype_))
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
    } // loop over my elements

    return;
  }

  /// create new elements and add them to the target discretization
  void CreateElements(
      Teuchos::RCP<DRT::Discretization> sourcedis,
      Teuchos::RCP<DRT::Discretization> targetdis,
      const int matid)
  {
    // now do the elements

    CloneStrategy::CheckMaterialType(matid);

    // prepare some variables we need
    int myrank = targetdis->Comm().MyPID();

    // construct scalar transport elements
    // The order of the elements might be different from that of the
    // fluid elements. We don't care. There are no dofs to these elements.
    for (unsigned i=0; i<egid_.size(); ++i)
    {
      DRT::Element* fluidele = sourcedis->gElement(egid_[i]);

      // create a scalar transport element with the same global element id
      RCP<DRT::Element> newele = DRT::UTILS::Factory(eletype_[i],"Polynomial",egid_[i],myrank);

      // get global node ids of fluid element
      vector<int> nids;
      nids.reserve(fluidele->NumNode());
      transform(fluidele->Nodes(), fluidele->Nodes()+fluidele->NumNode(),
                back_inserter(nids), mem_fun(&DRT::Node::Id));

      // set the same global node ids to the new scalar transport element
      newele->SetNodeIds(nids.size(), &nids[0]);

      // We need to set material and gauss points to complete element setup.
      // This is again really ugly as we have to extract the actual
      // element type in order to access the material property
      // note: SetMaterial() was reimplemented by the transport element!
      CloneStrategy::SetElementData(newele,fluidele,matid);

      // add new scalar transport element to discretization
      targetdis->AddElement(newele);
    }
  }

  //! set of row nodes
  set<int> rownodeset_;
  //! set of column nodes
  set<int> colnodeset_;
  //! vector for holding global element ids
  vector<int> egid_;
  //! vector for holding each (desired) element type string
  vector<string> eletype_;

}; // class DiscretizationCreator


} // namespace UTILS
} // namespace DRT

#endif // #ifndef DRT_UTILS_CREATEDIS_H
#endif // #ifdef CCADISCRET

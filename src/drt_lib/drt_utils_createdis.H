/*----------------------------------------------------------------------*/
/*! \file

\brief utility functions for automatic creation of a discretization
       from an existing one (e.g. ALE from Fluid)

\level 1

\maintainer Martin Kronbichler

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
      //! constructor
      explicit DiscretizationCreatorBase() : numeleskips_(0) {}

      //! destructor
      virtual ~DiscretizationCreatorBase() {}

      //! copy given conditions from source to target discretization
      void CopyConditions(const DRT::Discretization& sourcedis, DRT::Discretization& targetdis,
          const std::map<std::string, std::string>& conditions_to_copy) const;


      /**
       * generate identical clone of sourcedis (eles, nodes and conditions, dofs are only cloned
       * with flag clonedofs = true)
       * @param sourcedis: discretization based on which clone will be generated
       * @param targetdisname: name of cloned discretization
       * @param clonedofs: should Dofs be cloned?
       * @param assigndegreesoffreedom: flag for call to fill complete on cloned dis
       * @param initelements: flag for call to fill complete on cloned dis
       * @param doboundaryconditions: flag for call to fill complete on cloned dis
       * @return the cloned discretization
       */
      Teuchos::RCP<DRT::Discretization> CreateMatchingDiscretization(
          const Teuchos::RCP<DRT::Discretization>& sourcedis, const std::string& targetdisname,
          bool clonedofs = true, bool assigndegreesoffreedom = true, bool initelements = true,
          bool doboundaryconditions = true) const;

      //! Base class version for creation of matching discretization without material
      Teuchos::RCP<DRT::Discretization> CreateMatchingDiscretizationFromCondition(
          const DRT::Discretization& sourcedis,  ///< discretization with condition
          const DRT::Condition&
              cond,  ///< condition, from which the derived discretization is derived
          const std::string& discret_name,  ///< name of the new discretization
          const std::string& element_name,  ///< name/type of the elements to be created
          const std::vector<std::string>& conditions_to_copy  ///< list of conditions that will be
                                                              ///< copied to the new discretization
      )
      {
        if (cond.Nodes() == NULL or cond.Nodes()->size() == 0)
          dserror("The condition has no nodes!");

        // make sure connectivity is all set
        // we don't care, whether dofs exist or not
        if (!sourcedis.Filled()) dserror("sourcedis is not filled");

        // get this condition's elements
        std::map<int, Teuchos::RCP<DRT::Element>> sourceelements;
        const std::map<int, Teuchos::RCP<DRT::Element>>& geo = cond.Geometry();
        sourceelements.insert(geo.begin(), geo.end());

        return CreateMatchingDiscretizationFromCondition(
            sourcedis, sourceelements, discret_name, element_name, conditions_to_copy);
      };  // CreateMatchingDiscretizationFromCondition

      //! Base class version for creation of matching discretization without material
      Teuchos::RCP<DRT::Discretization> CreateMatchingDiscretizationFromCondition(
          const DRT::Discretization& sourcedis,  ///< discretization with condition
          const std::string& condname,           ///< name of the condition, by which the derived
                                                 ///< discretization is identified
          const std::string& discret_name,       ///< name of the new discretization
          const std::string&
              element_name,  ///< name/type of the elements to be created, if set to "" the natural
                             ///< element will be created (e.g. Fluid-->FluidBoundary,...)
          const std::vector<std::string>& conditions_to_copy,  ///< list of conditions that will be
                                                               ///< copied to the new discretization
          const int label = -1  ///< consider only conditions with specified label
      )
      {
        // make sure connectivity is all set
        // we don't care, whether dofs exist or not
        if (!sourcedis.Filled()) dserror("sourcedis is not filled");

        // We need to test for all elements (including ghosted ones) to
        // catch all nodes
        std::map<int, Teuchos::RCP<DRT::Element>> sourceelements;
        DRT::UTILS::FindConditionObjects(sourcedis, sourceelements, condname, label);

        return CreateMatchingDiscretizationFromCondition(
            sourcedis, sourceelements, discret_name, element_name, conditions_to_copy);
      };  // CreateDiscretizationFromCondition

      //! method for cloning a new discretization from an existing condition without material
      Teuchos::RCP<DRT::Discretization> CreateMatchingDiscretizationFromCondition(
          const DRT::Discretization& sourcedis,  ///< discretization with condition
          const std::map<int, Teuchos::RCP<DRT::Element>>&
              sourceelements,               ///< element map/geometry of the condition
          const std::string& discret_name,  ///< name of the new discretization
          const std::string&
              element_name,  ///< name/type of the elements to be created, if set to "" the natural
                             ///< element will be created (e.g. Fluid-->FluidBoundary,...)
          const std::vector<std::string>& conditions_to_copy  ///< list of conditions that will be
                                                              ///< copied to the new discretization
      )
      {
        Teuchos::RCP<Epetra_Comm> com = Teuchos::rcp(sourcedis.Comm().Clone());
        const int myrank = com->MyPID();
        const Epetra_Map* sourcenoderowmap = sourcedis.NodeRowMap();

        Teuchos::RCP<DRT::Discretization> targetdis =
            Teuchos::rcp(new DRT::Discretization(discret_name, com));

        // construct new elements
        for (std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator sourceele_iter =
                 sourceelements.begin();
             sourceele_iter != sourceelements.end(); ++sourceele_iter)
        {
          const Teuchos::RCP<DRT::Element> sourceele = sourceele_iter->second;

          // get global node ids
          std::vector<int> nids;
          nids.reserve(sourceele->NumNode());
          transform(sourceele->Nodes(), sourceele->Nodes() + sourceele->NumNode(),
              back_inserter(nids), std::mem_fun(&DRT::Node::Id));

          // check if element has nodes which are not in col map on this proc.
          // this should not be the case since each proc should have all nodes of
          // all owned or ghosted elements in the col map.
          if (std::count_if(nids.begin(), nids.end(), DRT::UTILS::MyGID(sourcedis.NodeColMap())) !=
              (int)(nids.size()))
          {
            dserror("element %d owned by proc %d has remote non-ghost nodes", sourceele->Id(),
                sourceele->Owner());
          }

          copy(nids.begin(), nids.end(), inserter(colnodeset_, colnodeset_.begin()));

          // copy node ids of condition ele to rownodeset but leave those that do
          // not belong to this processor
          remove_copy_if(nids.begin(), nids.end(), inserter(rownodeset_, rownodeset_.begin()),
              not1(DRT::UTILS::MyGID(sourcenoderowmap)));

          // Do not clone ghost elements here! Those will be handled by the
          // discretization itself.
          if (sourceele->Owner() == myrank)
          {
            Teuchos::RCP<DRT::Element> condele;
            if (element_name == "")
            {
              // copy the source ele (created in fill complete of the discretization)
              condele = Teuchos::rcp(sourceele->Clone(), true);
            }
            else
            {
              // create an element with the same global element id
              condele = DRT::UTILS::Factory(element_name, "Polynomial", sourceele->Id(), myrank);
              // set the same global node ids to the new element
              condele->SetNodeIds(nids.size(), &nids[0]);
            }
            // add element
            targetdis->AddElement(condele);
            roweleset_.insert(sourceele->Id());
          }
          coleleset_.insert(sourceele->Id());
        }  // loop over all source elements

        // construct new nodes, which use the same global id as the source nodes
        for (int i = 0; i < sourcenoderowmap->NumMyElements(); ++i)
        {
          const int gid = sourcenoderowmap->GID(i);
          if (rownodeset_.find(gid) != rownodeset_.end())
          {
            const DRT::Node* sourcenode = sourcedis.lRowNode(i);
            targetdis->AddNode(Teuchos::rcp(new DRT::Node(gid, sourcenode->X(), myrank)));
          }
        }

        // we get the element maps almost for free
        targetelerowmap_ = CreateMap(roweleset_, *targetdis);
        targetelecolmap_ = CreateMap(coleleset_, *targetdis);
        // we get the node maps almost for free
        targetnoderowmap_ = CreateMap(rownodeset_, *targetdis);
        targetnodecolmap_ = CreateMap(colnodeset_, *targetdis);

        // copy selected conditions to the new discretization
        for (std::vector<std::string>::const_iterator conditername = conditions_to_copy.begin();
             conditername != conditions_to_copy.end(); ++conditername)
        {
          std::vector<DRT::Condition*> conds;
          sourcedis.GetCondition(*conditername, conds);
          for (unsigned i = 0; i < conds.size(); ++i)
          {
            // We use the same nodal ids and therefore we can just copy the conditions.
            targetdis->SetCondition(*conditername, Teuchos::rcp(new DRT::Condition(*conds[i])));
          }
        }

        // we always skip the safety checks in Finalize()
        // because we create a discretization from a
        // conditioned subset of the source discretiztation.
        numeleskips_++;

        // call Redistribute, FillComplete etc.
        Finalize(sourcedis, *targetdis);

        return targetdis;
      };  // CreateDiscretizationFromCondition without material

     protected:
      //! construct row nodes for cloned target discretization
      void CreateNodes(const DRT::Discretization& sourcedis, DRT::Discretization& targetdis,
          const std::set<int>& rownodeset, const std::set<int>& colnodeset, const bool isnurbsdis,
          const bool buildimmersednode) const;

      //! construct and return Epetra_Map
      Teuchos::RCP<Epetra_Map> CreateMap(
          std::set<int>& gidset, const DRT::Discretization& targetdis) const;

      //! do some checks
      void InitialChecks(
          const DRT::Discretization& sourcedis, const DRT::Discretization& targetdis) const;

      //! export target nodes and elements and perform some checks
      void Finalize(const DRT::Discretization& sourcedis, DRT::Discretization& targetdis) const;

     protected:
      //! set of row nodes
      std::set<int> rownodeset_;
      //! set of column nodes
      std::set<int> colnodeset_;
      //! set of row elements
      std::set<int> roweleset_;
      //! set of column elements
      std::set<int> coleleset_;
      //! vector for holding each (desired) element type std::string
      std::vector<std::string> eletype_;
      //! map containing gids of owned nodes
      Teuchos::RCP<Epetra_Map> targetnoderowmap_;
      //! map containing gids of owned + ghosted nodes
      Teuchos::RCP<Epetra_Map> targetnodecolmap_;
      //! map containing gids of owned elements
      Teuchos::RCP<Epetra_Map> targetelerowmap_;
      //! map containing gids of owned + ghosted elements
      Teuchos::RCP<Epetra_Map> targetelecolmap_;
      //! local number of skipped elements during cloning
      int numeleskips_;

    };  // class DiscretizationCreatorBase


    /// class for cloning a new discretization from an existing one
    template <class CloneStrategy>
    class DiscretizationCreator : DiscretizationCreatorBase, CloneStrategy
    {
     public:
      /// constructor
      explicit DiscretizationCreator(){};

      /// destructor
      virtual ~DiscretizationCreator(){};

      /// Create the clone field material map from the input file
      void CreateCloneFieldMatMap(std::map<int, int>& matmap, const DRT::Discretization& sourcedis,
          const DRT::Discretization& targetdis) const
      {
        if (matmap.size()) dserror("The input material map is supposed to be empty!");

        std::map<std::pair<std::string, std::string>, std::map<int, int>> clonefieldmatmap =
            DRT::Problem::Instance()->CloningMaterialMap();
        if (clonefieldmatmap.size() < 1)
          dserror("At least one material pairing required in --CLONING MATERIAL MAP.");

        std::pair<std::string, std::string> key(sourcedis.Name(), targetdis.Name());
        matmap = clonefieldmatmap[key];
        if (matmap.size() < 1)
          dserror("Key pair '%s/%s' not defined in --CLONING MATERIAL MAP.",
              sourcedis.Name().c_str(), targetdis.Name().c_str());

        return;
      };  // CreateCloneFieldMatMap

      /// method for cloning a new discretization from an existing one
      void CreateMatchingDiscretization(
          Teuchos::RCP<DRT::Discretization> sourcedis,  ///< Teuchos::RCP to source discretization
          Teuchos::RCP<DRT::Discretization>
              targetdis,   ///< Teuchos::RCP to empty target discretization
          const int matid  ///< ID of the material which generated elements will get
      )
      {
        // obsolete function call which should not be used anymore!
        // let's use the more general version using an explicit material id mapping
        // as defined in the input file section "--CLONING MATERIAL MAP"

        // check and analyze source discretization (sorcedis must be filled!)
        InitialChecks(*sourcedis, *targetdis);

        // We have to find out all the material ids of the source discretization.
        // All cloned elements will receive the same material with the provided matid.
        std::map<int, int> matmap;
        int numelements = sourcedis->NumMyColElements();
        if (numelements < 1) dserror("At least one processor has no element");
        for (int i = 0; i < numelements; ++i)
        {
          DRT::Element* sourceele = sourcedis->lColElement(i);
          int src_matid = sourceele->Material()->Parameter()->Id();
          // if a new material id is found -> extend the map
          std::map<int, int>::iterator mat_iter = matmap.find(src_matid);
          if (mat_iter == matmap.end())
          {
            std::pair<int, int> matmappair(src_matid, matid);
            matmap.insert(matmappair);
          }
        }

        CreateMatchingDiscretization(sourcedis, targetdis, matmap);

        return;
      };  // CreateMatchingDiscretization

      /// method for cloning a new discretization from an existing one
      void CreateMatchingDiscretization(
          Teuchos::RCP<DRT::Discretization> sourcedis,  ///< Teuchos::RCP to source discretization
          Teuchos::RCP<DRT::Discretization>
              targetdis,  ///< Teuchos::RCP to empty target discretization
          const std::map<int, int>
              matmap  ///< map of material IDs (source element -> target element)
      )
      {
        // try to cast sourcedis to NurbsDiscretisation
        DRT::NURBS::NurbsDiscretization* nurbsdis =
            dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*(sourcedis)));
        bool isnurbsdis(nurbsdis != NULL);

        // try to cast source node to immersed node
        IMMERSED::ImmersedNode* inode = dynamic_cast<IMMERSED::ImmersedNode*>(
            sourcedis->gNode(sourcedis->NodeRowMap()->GID(0)));
        bool buildimmersednode(inode != NULL);

        // check and analyze source discretization
        InitialChecks(*sourcedis, *targetdis);
        AnalyzeSourceDis(sourcedis, eletype_, rownodeset_, colnodeset_, roweleset_, coleleset_);

        // do the node business
        CreateNodes(
            *sourcedis, *targetdis, rownodeset_, colnodeset_, isnurbsdis, buildimmersednode);
        targetnoderowmap_ = CreateMap(rownodeset_, *targetdis);
        targetnodecolmap_ = CreateMap(colnodeset_, *targetdis);

        // create elements
        CreateElements(sourcedis, targetdis, matmap, isnurbsdis);
        targetelerowmap_ = CreateMap(roweleset_, *targetdis);
        targetelecolmap_ = CreateMap(coleleset_, *targetdis);

        // copy desired conditions from source to target discretization
        std::map<std::string, std::string> conditions_to_copy = CloneStrategy::ConditionsToCopy();
        CopyConditions(*sourcedis, *targetdis, conditions_to_copy);

        // call Redistribute, FillComplete etc.
        Finalize(*sourcedis, *targetdis);
      };  // CreateMatchingDiscretization

      /// method for cloning a new discretization from an existing condition using the actual
      /// condition
      void CreateMatchingDiscretizationFromCondition(
          const DRT::Discretization& sourcedis,  ///< ref. to source discretization
          const std::vector<DRT::Condition*>&
              conds,  ///< vector of conditions containing the elements to clone
          DRT::Discretization& targetdis,  ///< Teuchos::RCP to empty target discretization
          const std::map<int, int>&
              matmap  ///< map of material IDs (source element -> target element)
      )
      {
        // check and analyze source and target discretization
        InitialChecks(sourcedis, targetdis);

        std::vector<DRT::Condition*>::const_iterator cit;
        for (cit = conds.begin(); cit != conds.end(); ++cit)
        {
          // check the source condition
          if ((*cit)->Nodes() == NULL or (*cit)->Nodes()->size() == 0)
            dserror("The condition has no nodes!");
        }

        // get this condition vector's elements
        std::map<int, Teuchos::RCP<DRT::Element>> sourceelements;
        DRT::UTILS::FindConditionObjects(sourceelements, conds);

        CreateMatchingDiscretizationFromCondition(sourcedis, sourceelements, targetdis, matmap);
        return;
      };  // CreateMatchingDiscretizationFromCondition

      /// method for cloning a new discretization from an existing condition using the condition
      /// name
      void CreateMatchingDiscretizationFromCondition(
          const DRT::Discretization& sourcedis,  ///< ref. to source discretization
          const std::string& condname,     ///< string to identify conditioned elements to clone
          DRT::Discretization& targetdis,  ///< Teuchos::RCP to empty target discretization
          const std::map<int, int>&
              matmap  ///< map of material IDs (source element -> target element)
      )
      {
        // check and analyze source discretization
        InitialChecks(sourcedis, targetdis);
        std::map<int, Teuchos::RCP<DRT::Element>> sourceelements;
        DRT::UTILS::FindConditionObjects(sourcedis, sourceelements, condname);

        CreateMatchingDiscretizationFromCondition(sourcedis, sourceelements, targetdis, matmap);
        return;
      };  // CreateMatchingDiscretizationFromCondition

     private:
      /// method for cloning a new discretization from an existing condition with material
      void CreateMatchingDiscretizationFromCondition(
          const DRT::Discretization& sourcedis,  ///< ref. to source discretization
          const std::map<int, Teuchos::RCP<DRT::Element>>&
              sourceelements,              ///< conditioned element map to clone
          DRT::Discretization& targetdis,  ///< Teuchos::RCP to empty target discretization
          const std::map<int, int>&
              matmap  ///< map of material IDs (source element -> target element)
      )
      {
        // try to cast sourcedis to NurbsDiscretisation
        const DRT::NURBS::NurbsDiscretization* nurbsdis_ptr =
            dynamic_cast<const DRT::NURBS::NurbsDiscretization*>(&sourcedis);
        bool isnurbsdis(nurbsdis_ptr != NULL);

        // try to cast source node to immersed node
        IMMERSED::ImmersedNode* inode =
            dynamic_cast<IMMERSED::ImmersedNode*>(sourcedis.gNode(sourcedis.NodeRowMap()->GID(0)));
        bool buildimmersednode(inode != NULL);

        AnalyzeConditionedSourceDis(
            sourcedis, sourceelements, eletype_, rownodeset_, colnodeset_, roweleset_, coleleset_);

        // do the node business
        CreateNodes(sourcedis, targetdis, rownodeset_, colnodeset_, isnurbsdis, buildimmersednode);
        targetnoderowmap_ = CreateMap(rownodeset_, targetdis);
        targetnodecolmap_ = CreateMap(colnodeset_, targetdis);

        // create elements
        CreateElementsFromCondition(sourceelements, targetdis, matmap, isnurbsdis);
        targetelerowmap_ = CreateMap(roweleset_, targetdis);
        targetelecolmap_ = CreateMap(coleleset_, targetdis);

        // copy desired conditions from source to target discretization
        std::map<std::string, std::string> conditions_to_copy = CloneStrategy::ConditionsToCopy();
        CopyConditions(sourcedis, targetdis, conditions_to_copy);

        // call Redistribute, FillComplete etc.
        Finalize(sourcedis, targetdis);
        return;
      };  // CreateMatchingDiscretizationFromCondition with material

      /// get element type std::strings and global id's and nodes from source discretization
      void AnalyzeSourceDis(Teuchos::RCP<DRT::Discretization> sourcedis,
          std::vector<std::string>& eletype, std::set<int>& rownodeset, std::set<int>& colnodeset,
          std::set<int>& roweleset, std::set<int>& coleleset)
      {
        const Epetra_Map* noderowmap = sourcedis->NodeRowMap();

        // We need to test for all elements (including ghosted ones) to
        // catch all nodes attached to the elements of the source discretization
        // we will clone only those (-> support for ALE sub-meshes)
        int numelements = sourcedis->NumMyColElements();

        for (int i = 0; i < numelements; ++i)
        {
          DRT::Element* actele = sourcedis->lColElement(i);
          bool ismyele = sourcedis->ElementRowMap()->MyGID(actele->Id());

          // we get the element type std::string and a boolean if this element
          // is considered! (see submeshes for Fluid-ALE case!)
          if (CloneStrategy::DetermineEleType(actele, ismyele, eletype))
          {
            // we make sure, that the cloned discretization
            // has the same parallel distribution as the
            // source discretization.
            if (ismyele) roweleset.insert(actele->Id());

            coleleset.insert(actele->Id());

            // copy node ids of actele to rownodeset but leave those that do
            // not belong to this processor
            remove_copy_if(actele->NodeIds(), actele->NodeIds() + actele->NumNode(),
                inserter(rownodeset, rownodeset.begin()), not1(DRT::UTILS::MyGID(noderowmap)));

            copy(actele->NodeIds(), actele->NodeIds() + actele->NumNode(),
                inserter(colnodeset, colnodeset.begin()));
          }
          else
            numeleskips_++;
        }  // loop over my elements
        return;
      };  // AnalyzeSourceDis

      /// get element type std::strings and global id's and nodes from conditioned source
      /// discretization
      void AnalyzeConditionedSourceDis(const DRT::Discretization& sourcedis,
          const std::map<int, Teuchos::RCP<DRT::Element>>& sourceelements,
          std::vector<std::string>& eletype, std::set<int>& rownodeset, std::set<int>& colnodeset,
          std::set<int>& roweleset, std::set<int>& coleleset)
      {
        const int myrank = sourcedis.Comm().MyPID();
        const Epetra_Map* sourcenoderowmap = sourcedis.NodeRowMap();
        const Epetra_Map* sourcenodecolmap = sourcedis.NodeColMap();

        // construct new elements
        std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator sourceele_iter;
        for (sourceele_iter = sourceelements.begin(); sourceele_iter != sourceelements.end();
             ++sourceele_iter)
        {
          const Teuchos::RCP<DRT::Element> actele = sourceele_iter->second;
          const bool ismyele = (actele->Owner() == myrank);

          // we get the element type std::string and a boolean if this element
          // is considered! (see submeshes for Fluid-ALE case!)
          if (CloneStrategy::DetermineEleType(&(*actele), ismyele, eletype))
          {
            // we make sure, that the cloned discretization
            // has the same parallel distribution as the
            // source discretization.
            // as prerequisite it is required that the parallel
            // distribution of the condition geometry matches
            // the distribution of the underlying parent elements.
            // this is made sure in drt_discret_condition.cpp
            // (rauch 10/16).
            if (ismyele) roweleset.insert(actele->Id());

            coleleset.insert(actele->Id());

            // get global node ids
            std::vector<int> nids;
            nids.reserve(actele->NumNode());
            transform(actele->Nodes(), actele->Nodes() + actele->NumNode(), back_inserter(nids),
                std::mem_fun(&DRT::Node::Id));

            // check if element has nodes, which are not in col map on this proc.
            // this should not be, since each proc should have all nodes of all
            // owned, or ghosted elements in the col map.
            if (std::count_if(nids.begin(), nids.end(), DRT::UTILS::MyGID(sourcenodecolmap)) !=
                (int)(nids.size()))
            {
              dserror("element %d owned by proc %d has remote non-ghost nodes", actele->Id(),
                  actele->Owner());
            }

            // copy node ids of condition ele to set of column nodes
            copy(nids.begin(), nids.end(), inserter(colnodeset, colnodeset.begin()));

            // copy node ids of condition ele to rownodeset except for those which do
            // not belong to this processor
            remove_copy_if(nids.begin(), nids.end(), inserter(rownodeset, rownodeset.begin()),
                std::not1(DRT::UTILS::MyGID(sourcenoderowmap)));
          }
        }

        // we always skip the safety checks in Finalize()
        // because we create a discretization from a
        // conditioned subset of the source discretiztation.
        numeleskips_++;
        return;
      };  // AnalyzeConditionedSourceDis

      /// create new elements and add them to the target discretization
      void CreateElements(Teuchos::RCP<DRT::Discretization> sourcedis,
          Teuchos::RCP<DRT::Discretization> targetdis, std::map<int, int> matmap,
          const bool isnurbsdis)
      {
        // now do the elements
        for (std::map<int, int>::iterator mapit = matmap.begin(); mapit != matmap.end(); ++mapit)
        {
          int target_id = mapit->second;
          CloneStrategy::CheckMaterialType(target_id);
        }

        // prepare some variables we need
        int myrank = targetdis->Comm().MyPID();

        // construct new elements
        // The order of the elements might be different from that of the
        // source elements. We don't care. There are no dofs to these elements.
        std::set<int>::iterator it = roweleset_.begin();
        for (std::size_t i = 0; i < roweleset_.size(); ++i)
        {
          DRT::Element* sourceele = sourcedis->gElement(*it);

          std::string approxtype = "Polynomial";
          if (isnurbsdis)
          {
            if (sourceele->NumNode() == 8)
            {
              approxtype = "NURBS8";
            }
            else if (sourceele->NumNode() == 9)
            {
              approxtype = "NURBS9";
            }
            else if (sourceele->NumNode() == 4)
            {
              approxtype = "NURBS4";
            }
            else if (sourceele->NumNode() == 27)
            {
              approxtype = "NURBS27";
            }
            else if (sourceele->NumNode() == 2)
            {
              approxtype = "NURBS2";
            }
            else if (sourceele->NumNode() == 3)
            {
              approxtype = "NURBS3";
            }
            else
            {
              dserror("unknown type of nurbs element\n");
            }
          }

          // create a new element of desired type with the same global element id
          Teuchos::RCP<DRT::Element> newele =
              DRT::UTILS::Factory(eletype_[i], approxtype, *it, myrank);

          // get global node ids of source element
          std::vector<int> nids;
          nids.reserve(sourceele->NumNode());
          transform(sourceele->Nodes(), sourceele->Nodes() + sourceele->NumNode(),
              back_inserter(nids), std::mem_fun(&DRT::Node::Id));

          // set the same global node ids to the new element
          newele->SetNodeIds(nids.size(), &nids[0]);

          // We need to set material and gauss points to complete element setup.
          // This is again really ugly as we have to extract the actual
          // element type in order to access the material property
          // note: SetMaterial() was reimplemented by the transport element!

          int src_matid = sourceele->Material()->Parameter()->Id();
          std::map<int, int>::iterator mat_iter = matmap.find(src_matid);
          if (mat_iter != matmap.end())
          {
            int tar_matid = mat_iter->second;
            CloneStrategy::SetElementData(newele, sourceele, tar_matid, isnurbsdis);

            // add new element to discretization
            targetdis->AddElement(newele);
          }
          else
          {
            // before we stop, print the material id map
            std::cout << "Material map on PROC " << myrank << ":" << std::endl;
            for (mat_iter = matmap.begin(); mat_iter != matmap.end(); mat_iter++)
              std::cout << mat_iter->first << " -> " << mat_iter->second << std::endl;

            dserror("no matching material ID (%d) in map", src_matid);
          }
          it++;
        }
        return;
      };  // CreateElements

      /// create new elements from the condition and add them to the target discretization
      void CreateElementsFromCondition(
          const std::map<int, Teuchos::RCP<DRT::Element>>& sourceelements,
          DRT::Discretization& targetdis, const std::map<int, int>& matmap, const bool& isnurbsdis)
      {
        // now do the elements
        for (std::map<int, int>::const_iterator mapit = matmap.begin(); mapit != matmap.end();
             ++mapit)
        {
          int target_id = mapit->second;
          CloneStrategy::CheckMaterialType(target_id);
        }

        // prepare some variables we need
        int myrank = targetdis.Comm().MyPID();

        // construct new elements
        // The order of the elements might be different from that of the
        // source elements. We don't care. There are no dofs to these elements.
        std::set<int>::iterator it = roweleset_.begin();
        for (std::size_t i = 0; i < roweleset_.size(); ++i)
        {
          std::map<int, Teuchos::RCP<DRT::Element>>::const_iterator src_ele_citer =
              sourceelements.find(*it);
          if (src_ele_citer == sourceelements.end())
            dserror(
                "The source element %d could not be found in the source "
                "condition element map!",
                *it);

          DRT::Element* sourceele = src_ele_citer->second.get();
          if (sourceele == NULL) dserror("The sourceele pointer is NULL!");

          std::string approxtype = "Polynomial";
          if (isnurbsdis)
          {
            if (sourceele->NumNode() == 8)
            {
              approxtype = "NURBS8";
            }
            else if (sourceele->NumNode() == 9)
            {
              approxtype = "NURBS9";
            }
            else if (sourceele->NumNode() == 4)
            {
              approxtype = "NURBS4";
            }
            else if (sourceele->NumNode() == 27)
            {
              approxtype = "NURBS27";
            }
            else if (sourceele->NumNode() == 2)
            {
              approxtype = "NURBS2";
            }
            else if (sourceele->NumNode() == 3)
            {
              approxtype = "NURBS3";
            }
            else
            {
              dserror("unknown type of nurbs element\n");
            }
          }

          // get owner of source element
          const int sourceeleowner = sourceele->Owner();
          if (myrank != sourceeleowner) dserror("roweleset_ should only contain my element gids!");

          // create a new element of desired type with the same global element id and same owner as
          // source element
          Teuchos::RCP<DRT::Element> newele =
              DRT::UTILS::Factory(eletype_[i], approxtype, *it, sourceeleowner);

          // get global node ids of fluid element
          std::vector<int> nids;
          nids.reserve(sourceele->NumNode());
          transform(sourceele->Nodes(), sourceele->Nodes() + sourceele->NumNode(),
              back_inserter(nids), std::mem_fun(&DRT::Node::Id));

          // set the same global node ids to the new element
          newele->SetNodeIds(nids.size(), &nids[0]);

          // We need to set material and gauss points to complete element setup.
          // This is again really ugly as we have to extract the actual
          // element type in order to access the material property
          // note: SetMaterial() was reimplemented by the transport element!
          Teuchos::RCP<MAT::Material> mat_ptr = sourceele->Material();
          /* Check if the material pointer is null. If necessary, try to cast
           * the condition element to a FaceElement and ask the parent element for
           * the material.                                                      */
          if (mat_ptr.is_null())
          {
            DRT::FaceElement* src_face_element = dynamic_cast<DRT::FaceElement*>(sourceele);
            if (src_face_element != NULL) mat_ptr = src_face_element->ParentElement()->Material();
          }
          // It is no FaceElement or the material pointer of the parent element is NULL.
          if (mat_ptr.is_null()) dserror("The condition element has no material!");

          int src_matid = mat_ptr->Parameter()->Id();
          std::map<int, int>::const_iterator mat_iter = matmap.find(src_matid);
          if (mat_iter != matmap.end())
          {
            int tar_matid = mat_iter->second;
            CloneStrategy::SetElementData(newele, sourceele, tar_matid, isnurbsdis);

            // add new element to discretization
            targetdis.AddElement(newele);
          }
          else
          {
            // before we stop, print the material id map
            std::cout << "Material map on PROC " << myrank << ":" << std::endl;
            for (mat_iter = matmap.begin(); mat_iter != matmap.end(); mat_iter++)
              std::cout << mat_iter->first << " -> " << mat_iter->second << std::endl;

            dserror("no matching material ID (%d) in map", src_matid);
          }
          it++;
        }
        return;
      }  // CreateElementsFromCondition

    };  // class DiscretizationCreator


    /// clone target discretization from a given source discretization
    template <class CloneStrategy>
    void CloneDiscretization(
        Teuchos::RCP<DRT::Discretization> sourcedis,  ///< source discretization
        Teuchos::RCP<DRT::Discretization> targetdis   ///< target discretization
    )
    {
      // access the communicator for time measurement
      const Epetra_Comm& comm = sourcedis->Comm();
      Epetra_Time time(comm);

      // create target discretization using a given clone strategy
      {
        Teuchos::RCP<DRT::UTILS::DiscretizationCreator<CloneStrategy>> clonewizard =
            Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<CloneStrategy>());

        std::map<int, int> matmap;
        clonewizard->CreateCloneFieldMatMap(matmap, *sourcedis, *targetdis);

        clonewizard->CreateMatchingDiscretization(sourcedis, targetdis, matmap);
      }
      if (comm.MyPID() == 0)
      {
        IO::cout << "Created discretization " << targetdis->Name()
                 << " as a clone of discretization " << sourcedis->Name() << " in...."
                 << time.ElapsedTime() << " secs\n\n";
      }
      return;
    };  // CloneDiscretization

    /// clone target discretization from a given condition of the source discretization
    template <class CloneStrategy>
    void CloneDiscretizationFromCondition(
        const DRT::DiscretizationInterface& sourcedis,  ///< source discretization
        DRT::DiscretizationInterface& targetdis,        ///< target discretization
        const std::vector<DRT::Condition*>& conds       ///< source conditions to clone from
    )
    {
      const DRT::Discretization* sourcedis_ptr =
          dynamic_cast<const DRT::Discretization*>(&sourcedis);
      if (sourcedis_ptr == NULL) dserror("Cast of the source discretization failed!");
      DRT::Discretization* targetdis_ptr = dynamic_cast<DRT::Discretization*>(&targetdis);
      if (targetdis_ptr == NULL) dserror("Cast of the target discretization failed!");

      // access the communicator for time measurement
      const Epetra_Comm& comm = sourcedis_ptr->Comm();
      Epetra_Time time(comm);

      // create target discretization using a given clone strategy
      {
        Teuchos::RCP<DRT::UTILS::DiscretizationCreator<CloneStrategy>> clonewizard =
            Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<CloneStrategy>());

        std::map<int, int> matmap;
        clonewizard->CreateCloneFieldMatMap(matmap, *sourcedis_ptr, *targetdis_ptr);

        clonewizard->CreateMatchingDiscretizationFromCondition(
            *sourcedis_ptr, conds, *targetdis_ptr, matmap);
      }
      if (comm.MyPID() == 0)
      {
        IO::cout << "Created discretization " << targetdis_ptr->Name()
                 << " as a clone from the condition(s) with ID(s)=";
        for (unsigned int i = 0; i < conds.size(); ++i) IO::cout << " " << conds[i]->Id();
        IO::cout << " of the discretization " << sourcedis_ptr->Name() << " in...."
                 << time.ElapsedTime() << " secs\n\n";
      }
      return;
    };  // CloneDiscretizationFromCondition

    /// clone target discretization from a given condition of the source discretization
    template <class CloneStrategy>
    void CloneDiscretizationFromCondition(
        const DRT::Discretization& sourcedis,  ///< source discretization
        DRT::Discretization& targetdis,        ///< target discretization
        const std::string& condname            ///< source condition name to clone from
    )
    {
      // access the communicator for time measurement
      const Epetra_Comm& comm = sourcedis.Comm();
      Epetra_Time time(comm);

      // create target discretization using a given clone strategy
      {
        Teuchos::RCP<DRT::UTILS::DiscretizationCreator<CloneStrategy>> clonewizard =
            Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<CloneStrategy>());

        std::map<int, int> matmap;
        clonewizard->CreateCloneFieldMatMap(matmap, sourcedis, targetdis);

        clonewizard->CreateMatchingDiscretizationFromCondition(
            sourcedis, condname, targetdis, matmap);
      }
      if (comm.MyPID() == 0)
      {
        IO::cout << "Created discretization " << targetdis.Name()
                 << " as a clone from the condition \"" << condname.c_str()
                 << "\" of the discretization " << sourcedis.Name() << " in...."
                 << time.ElapsedTime() << " secs\n\n";
      }
      return;
    };  // CloneDiscretizationFromCondition

    //! construct and return cloning material map input lines
    Teuchos::RCP<DRT::INPUT::Lines> ValidCloningMaterialMapLines();

    //! helper method for printout
    void PrintCloningMaterialMapDatHeader();

  }  // namespace UTILS
}  // namespace DRT

#endif  // #ifndef DRT_UTILS_CREATEDIS_H

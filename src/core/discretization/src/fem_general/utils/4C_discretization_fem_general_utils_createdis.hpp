/*----------------------------------------------------------------------*/
/*! \file

\brief utility functions for automatic creation of a discretization
       from an existing one (e.g. ALE from Fluid)

\level 1


*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_DISCRETIZATION_FEM_GENERAL_UTILS_CREATEDIS_HPP
#define FOUR_C_DISCRETIZATION_FEM_GENERAL_UTILS_CREATEDIS_HPP

#include "4C_config.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_discretization_condition_utils.hpp"
#include "4C_discretization_fem_general_immersed_node.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_io_pstream.hpp"
#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_nurbs_discret.hpp"

#include <Teuchos_Time.hpp>

FOUR_C_NAMESPACE_OPEN

namespace INPUT
{
  class LineDefinition;
}

namespace CORE::FE
{
  /// class providing basic functionality for cloning discretizations
  class DiscretizationCreatorBase
  {
   public:
    //! constructor
    explicit DiscretizationCreatorBase() : numeleskips_(0) {}

    //! destructor
    virtual ~DiscretizationCreatorBase() = default;
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
    Teuchos::RCP<DRT::Discretization> create_matching_discretization(
        const Teuchos::RCP<DRT::Discretization>& sourcedis, const std::string& targetdisname,
        bool clonedofs = true, bool assigndegreesoffreedom = true, bool initelements = true,
        bool doboundaryconditions = true) const;

    //! Base class version for creation of matching discretization without material
    Teuchos::RCP<DRT::Discretization> create_matching_discretization_from_condition(
        const DRT::Discretization& sourcedis,  ///< discretization with condition
        const CORE::Conditions::Condition&
            cond,  ///< condition, from which the derived discretization is derived
        const std::string& discret_name,  ///< name of the new discretization
        const std::string& element_name,  ///< name/type of the elements to be created
        const std::vector<std::string>& conditions_to_copy  ///< list of conditions that will be
                                                            ///< copied to the new discretization
    )
    {
      if (cond.GetNodes() == nullptr or cond.GetNodes()->size() == 0)
        FOUR_C_THROW("The condition has no nodes!");

      // make sure connectivity is all set
      // we don't care, whether dofs exist or not
      if (!sourcedis.Filled()) FOUR_C_THROW("sourcedis is not filled");

      // get this condition's elements
      std::map<int, Teuchos::RCP<CORE::Elements::Element>> sourceelements;
      const std::map<int, Teuchos::RCP<CORE::Elements::Element>>& geo = cond.Geometry();
      sourceelements.insert(geo.begin(), geo.end());

      return create_matching_discretization_from_condition(
          sourcedis, sourceelements, discret_name, element_name, conditions_to_copy);
    };  // create_matching_discretization_from_condition

    //! Base class version for creation of matching discretization without material
    Teuchos::RCP<DRT::Discretization> create_matching_discretization_from_condition(
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
      if (!sourcedis.Filled()) FOUR_C_THROW("sourcedis is not filled");

      // We need to test for all elements (including ghosted ones) to
      // catch all nodes
      std::map<int, Teuchos::RCP<CORE::Elements::Element>> sourceelements;
      CORE::Conditions::FindConditionObjects(sourcedis, sourceelements, condname, label);

      return create_matching_discretization_from_condition(
          sourcedis, sourceelements, discret_name, element_name, conditions_to_copy);
    };  // create_discretization_from_condition

    //! method for cloning a new discretization from an existing condition without material
    Teuchos::RCP<DRT::Discretization> create_matching_discretization_from_condition(
        const DRT::Discretization& sourcedis,  ///< discretization with condition
        const std::map<int, Teuchos::RCP<CORE::Elements::Element>>&
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
          Teuchos::rcp(new DRT::Discretization(discret_name, com, sourcedis.n_dim()));

      // construct new elements
      for (std::map<int, Teuchos::RCP<CORE::Elements::Element>>::const_iterator sourceele_iter =
               sourceelements.begin();
           sourceele_iter != sourceelements.end(); ++sourceele_iter)
      {
        const Teuchos::RCP<CORE::Elements::Element> sourceele = sourceele_iter->second;

        // get global node ids
        std::vector<int> nids;
        nids.reserve(sourceele->num_node());
        transform(sourceele->Nodes(), sourceele->Nodes() + sourceele->num_node(),
            back_inserter(nids), std::mem_fn(&CORE::Nodes::Node::Id));

        // check if element has nodes which are not in col map on this proc.
        // this should not be the case since each proc should have all nodes of
        // all owned or ghosted elements in the col map.
        if (std::count_if(nids.begin(), nids.end(),
                CORE::Conditions::MyGID(sourcedis.NodeColMap())) != (int)(nids.size()))
        {
          FOUR_C_THROW("element %d owned by proc %d has remote non-ghost nodes", sourceele->Id(),
              sourceele->Owner());
        }

        copy(nids.begin(), nids.end(), inserter(colnodeset_, colnodeset_.begin()));

        // copy node ids of condition ele to rownodeset but leave those that do
        // not belong to this processor
        remove_copy_if(nids.begin(), nids.end(), inserter(rownodeset_, rownodeset_.begin()),
            std::not_fn(CORE::Conditions::MyGID(sourcenoderowmap)));

        // Do not clone ghost elements here! Those will be handled by the
        // discretization itself.
        if (sourceele->Owner() == myrank)
        {
          Teuchos::RCP<CORE::Elements::Element> condele;
          if (element_name == "")
          {
            // copy the source ele (created in fill complete of the discretization)
            condele = Teuchos::rcp(sourceele->Clone(), true);
          }
          else
          {
            // create an element with the same global element id
            condele = CORE::COMM::Factory(element_name, "Polynomial", sourceele->Id(), myrank);
            // set the same global node ids to the new element
            condele->SetNodeIds(nids.size(), nids.data());
          }
          // add element
          targetdis->add_element(condele);
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
          const CORE::Nodes::Node* sourcenode = sourcedis.lRowNode(i);
          targetdis->AddNode(Teuchos::rcp(new CORE::Nodes::Node(gid, sourcenode->X(), myrank)));
        }
      }

      // we get the element maps almost for free
      targetelerowmap_ = create_map(roweleset_, *targetdis);
      targetelecolmap_ = create_map(coleleset_, *targetdis);
      // we get the node maps almost for free
      targetnoderowmap_ = create_map(rownodeset_, *targetdis);
      targetnodecolmap_ = create_map(colnodeset_, *targetdis);

      // copy selected conditions to the new discretization
      for (const auto& cond_name : conditions_to_copy)
      {
        std::vector<CORE::Conditions::Condition*> conds;
        sourcedis.GetCondition(cond_name, conds);
        for (const auto& cond : conds)
        {
          // We use the same nodal ids and therefore we can just copy the conditions.
          targetdis->SetCondition(cond_name, cond->copy_without_geometry());
        }
      }

      // we always skip the safety checks in Finalize()
      // because we create a discretization from a
      // conditioned subset of the source discretiztation.
      numeleskips_++;

      // call Redistribute, fill_complete etc.
      finalize(sourcedis, *targetdis);

      return targetdis;
    };  // create_discretization_from_condition without material

   protected:
    //! construct row nodes for cloned target discretization
    void create_nodes(const DRT::Discretization& sourcedis, DRT::Discretization& targetdis,
        const std::set<int>& rownodeset, const std::set<int>& colnodeset, const bool isnurbsdis,
        const bool buildimmersednode) const;

    //! construct and return Epetra_Map
    Teuchos::RCP<Epetra_Map> create_map(
        std::set<int>& gidset, const DRT::Discretization& targetdis) const;

    //! do some checks
    void initial_checks(
        const DRT::Discretization& sourcedis, const DRT::Discretization& targetdis) const;

    //! export target nodes and elements and perform some checks
    void finalize(const DRT::Discretization& sourcedis, DRT::Discretization& targetdis) const;

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

    /// Create the clone field material map from the input file
    void create_clone_field_mat_map(std::map<int, int>& matmap,
        const DRT::Discretization& sourcedis, const DRT::Discretization& targetdis,
        const std::map<std::pair<std::string, std::string>, std::map<int, int>>& clonefieldmatmap)
        const
    {
      if (matmap.size()) FOUR_C_THROW("The input material map is supposed to be empty!");

      if (clonefieldmatmap.size() < 1)
        FOUR_C_THROW("At least one material pairing required in --CLONING MATERIAL MAP.");

      std::pair<std::string, std::string> key(sourcedis.Name(), targetdis.Name());
      matmap = clonefieldmatmap.at(key);
      if (matmap.size() < 1)
        FOUR_C_THROW("Key pair '%s/%s' not defined in --CLONING MATERIAL MAP.",
            sourcedis.Name().c_str(), targetdis.Name().c_str());

      return;
    };  // create_clone_field_mat_map

    /// method for cloning a new discretization from an existing one
    void create_matching_discretization(
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
      initial_checks(*sourcedis, *targetdis);

      // We have to find out all the material ids of the source discretization.
      // All cloned elements will receive the same material with the provided matid.
      std::map<int, int> matmap;
      int numelements = sourcedis->NumMyColElements();
      if (numelements < 1) FOUR_C_THROW("At least one processor has no element");
      for (int i = 0; i < numelements; ++i)
      {
        CORE::Elements::Element* sourceele = sourcedis->lColElement(i);
        int src_matid = sourceele->Material()->Parameter()->Id();
        // if a new material id is found -> extend the map
        std::map<int, int>::iterator mat_iter = matmap.find(src_matid);
        if (mat_iter == matmap.end())
        {
          std::pair<int, int> matmappair(src_matid, matid);
          matmap.insert(matmappair);
        }
      }

      create_matching_discretization(sourcedis, targetdis, matmap);

      return;
    };  // create_matching_discretization

    /// method for cloning a new discretization from an existing one
    void create_matching_discretization(
        Teuchos::RCP<DRT::Discretization> sourcedis,  ///< Teuchos::RCP to source discretization
        Teuchos::RCP<DRT::Discretization>
            targetdis,                   ///< Teuchos::RCP to empty target discretization
        const std::map<int, int> matmap  ///< map of material IDs (source element -> target element)
    )
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (!(sourcedis->HaveGlobalNode(sourcedis->NodeRowMap()->GID(0))))
        FOUR_C_THROW("Cloning not possible since node with GID %d is not stored on this proc!",
            sourcedis->NodeRowMap()->GID(0));
#endif
      // try to cast sourcedis to NurbsDiscretisation
      DRT::NURBS::NurbsDiscretization* nurbsdis =
          dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*(sourcedis)));
      bool isnurbsdis(nurbsdis != nullptr);

      // try to cast source node to immersed node
      CORE::Nodes::ImmersedNode* inode = dynamic_cast<CORE::Nodes::ImmersedNode*>(
          sourcedis->gNode(sourcedis->NodeRowMap()->GID(0)));
      bool buildimmersednode(inode != nullptr);

      // check and analyze source discretization
      initial_checks(*sourcedis, *targetdis);
      analyze_source_dis(sourcedis, eletype_, rownodeset_, colnodeset_, roweleset_, coleleset_);

      // do the node business
      create_nodes(*sourcedis, *targetdis, rownodeset_, colnodeset_, isnurbsdis, buildimmersednode);
      targetnoderowmap_ = create_map(rownodeset_, *targetdis);
      targetnodecolmap_ = create_map(colnodeset_, *targetdis);

      // create elements
      create_elements(sourcedis, targetdis, matmap, isnurbsdis);
      targetelerowmap_ = create_map(roweleset_, *targetdis);
      targetelecolmap_ = create_map(coleleset_, *targetdis);

      // copy desired conditions from source to target discretization
      const auto conditions_to_copy = CloneStrategy::conditions_to_copy();
      CopyConditions(*sourcedis, *targetdis, conditions_to_copy);

      // call Redistribute, fill_complete etc.
      finalize(*sourcedis, *targetdis);
    };  // create_matching_discretization

    /// method for cloning a new discretization from an existing condition using the actual
    /// condition
    void create_matching_discretization_from_condition(
        const DRT::Discretization& sourcedis,  ///< ref. to source discretization
        const std::vector<CORE::Conditions::Condition*>&
            conds,                       ///< vector of conditions containing the elements to clone
        DRT::Discretization& targetdis,  ///< Teuchos::RCP to empty target discretization
        const std::map<int, int>&
            matmap  ///< map of material IDs (source element -> target element)
    )
    {
      // check and analyze source and target discretization
      initial_checks(sourcedis, targetdis);

      std::vector<CORE::Conditions::Condition*>::const_iterator cit;
      for (cit = conds.begin(); cit != conds.end(); ++cit)
      {
        // check the source condition
        if ((*cit)->GetNodes() == nullptr or (*cit)->GetNodes()->size() == 0)
          FOUR_C_THROW("The condition has no nodes!");
      }

      // get this condition vector's elements
      std::map<int, Teuchos::RCP<CORE::Elements::Element>> sourceelements;
      CORE::Conditions::FindConditionObjects(sourceelements, conds);

      create_matching_discretization_from_condition(sourcedis, sourceelements, targetdis, matmap);
      return;
    };  // create_matching_discretization_from_condition

    /// method for cloning a new discretization from an existing condition using the condition
    /// name
    void create_matching_discretization_from_condition(
        const DRT::Discretization& sourcedis,  ///< ref. to source discretization
        const std::string& condname,           ///< string to identify conditioned elements to clone
        DRT::Discretization& targetdis,        ///< Teuchos::RCP to empty target discretization
        const std::map<int, int>&
            matmap  ///< map of material IDs (source element -> target element)
    )
    {
      // check and analyze source discretization
      initial_checks(sourcedis, targetdis);
      std::map<int, Teuchos::RCP<CORE::Elements::Element>> sourceelements;
      CORE::Conditions::FindConditionObjects(sourcedis, sourceelements, condname);

      create_matching_discretization_from_condition(sourcedis, sourceelements, targetdis, matmap);
      return;
    };  // create_matching_discretization_from_condition

   private:
    /// method for cloning a new discretization from an existing condition with material
    void create_matching_discretization_from_condition(
        const DRT::Discretization& sourcedis,  ///< ref. to source discretization
        const std::map<int, Teuchos::RCP<CORE::Elements::Element>>&
            sourceelements,              ///< conditioned element map to clone
        DRT::Discretization& targetdis,  ///< Teuchos::RCP to empty target discretization
        const std::map<int, int>&
            matmap  ///< map of material IDs (source element -> target element)
    )
    {
      // try to cast sourcedis to NurbsDiscretisation
      const DRT::NURBS::NurbsDiscretization* nurbsdis_ptr =
          dynamic_cast<const DRT::NURBS::NurbsDiscretization*>(&sourcedis);
      bool isnurbsdis(nurbsdis_ptr != nullptr);

      // try to cast source node to immersed node
      CORE::Nodes::ImmersedNode* inode =
          dynamic_cast<CORE::Nodes::ImmersedNode*>(sourcedis.gNode(sourcedis.NodeRowMap()->GID(0)));
      bool buildimmersednode(inode != nullptr);

      analyze_conditioned_source_dis(
          sourcedis, sourceelements, eletype_, rownodeset_, colnodeset_, roweleset_, coleleset_);

      // do the node business
      create_nodes(sourcedis, targetdis, rownodeset_, colnodeset_, isnurbsdis, buildimmersednode);
      targetnoderowmap_ = create_map(rownodeset_, targetdis);
      targetnodecolmap_ = create_map(colnodeset_, targetdis);

      // create elements
      create_elements_from_condition(sourceelements, targetdis, matmap, isnurbsdis);
      targetelerowmap_ = create_map(roweleset_, targetdis);
      targetelecolmap_ = create_map(coleleset_, targetdis);

      // copy desired conditions from source to target discretization
      const auto conditions_to_copy = CloneStrategy::conditions_to_copy();
      CopyConditions(sourcedis, targetdis, conditions_to_copy);

      // call Redistribute, fill_complete etc.
      finalize(sourcedis, targetdis);
    };  // create_matching_discretization_from_condition with material

    /// get element type std::strings and global id's and nodes from source discretization
    void analyze_source_dis(Teuchos::RCP<DRT::Discretization> sourcedis,
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
        CORE::Elements::Element* actele = sourcedis->lColElement(i);
        bool ismyele = sourcedis->ElementRowMap()->MyGID(actele->Id());

        // we get the element type std::string and a boolean if this element
        // is considered! (see submeshes for Fluid-ALE case!)
        if (CloneStrategy::determine_ele_type(actele, ismyele, eletype))
        {
          // we make sure, that the cloned discretization
          // has the same parallel distribution as the
          // source discretization.
          if (ismyele) roweleset.insert(actele->Id());

          coleleset.insert(actele->Id());

          // copy node ids of actele to rownodeset but leave those that do
          // not belong to this processor
          remove_copy_if(actele->NodeIds(), actele->NodeIds() + actele->num_node(),
              inserter(rownodeset, rownodeset.begin()),
              std::not_fn(CORE::Conditions::MyGID(noderowmap)));

          copy(actele->NodeIds(), actele->NodeIds() + actele->num_node(),
              inserter(colnodeset, colnodeset.begin()));
        }
        else
          numeleskips_++;
      }  // loop over my elements
      return;
    };  // analyze_source_dis

    /// get element type std::strings and global id's and nodes from conditioned source
    /// discretization
    void analyze_conditioned_source_dis(const DRT::Discretization& sourcedis,
        const std::map<int, Teuchos::RCP<CORE::Elements::Element>>& sourceelements,
        std::vector<std::string>& eletype, std::set<int>& rownodeset, std::set<int>& colnodeset,
        std::set<int>& roweleset, std::set<int>& coleleset)
    {
      const int myrank = sourcedis.Comm().MyPID();
      const Epetra_Map* sourcenoderowmap = sourcedis.NodeRowMap();
      const Epetra_Map* sourcenodecolmap = sourcedis.NodeColMap();

      // construct new elements
      std::map<int, Teuchos::RCP<CORE::Elements::Element>>::const_iterator sourceele_iter;
      for (sourceele_iter = sourceelements.begin(); sourceele_iter != sourceelements.end();
           ++sourceele_iter)
      {
        const Teuchos::RCP<CORE::Elements::Element> actele = sourceele_iter->second;
        const bool ismyele = (actele->Owner() == myrank);

        // we get the element type std::string and a boolean if this element
        // is considered! (see submeshes for Fluid-ALE case!)
        if (CloneStrategy::determine_ele_type(&(*actele), ismyele, eletype))
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
          nids.reserve(actele->num_node());
          transform(actele->Nodes(), actele->Nodes() + actele->num_node(), back_inserter(nids),
              std::mem_fn(&CORE::Nodes::Node::Id));

          // check if element has nodes, which are not in col map on this proc.
          // this should not be, since each proc should have all nodes of all
          // owned, or ghosted elements in the col map.
          if (std::count_if(nids.begin(), nids.end(), CORE::Conditions::MyGID(sourcenodecolmap)) !=
              (int)(nids.size()))
          {
            FOUR_C_THROW("element %d owned by proc %d has remote non-ghost nodes", actele->Id(),
                actele->Owner());
          }

          // copy node ids of condition ele to set of column nodes
          copy(nids.begin(), nids.end(), inserter(colnodeset, colnodeset.begin()));

          // copy node ids of condition ele to rownodeset except for those which do
          // not belong to this processor
          remove_copy_if(nids.begin(), nids.end(), inserter(rownodeset, rownodeset.begin()),
              std::not_fn(CORE::Conditions::MyGID(sourcenoderowmap)));
        }
      }

      // we always skip the safety checks in Finalize()
      // because we create a discretization from a
      // conditioned subset of the source discretiztation.
      numeleskips_++;
      return;
    };  // analyze_conditioned_source_dis

    /// create new elements and add them to the target discretization
    void create_elements(Teuchos::RCP<DRT::Discretization> sourcedis,
        Teuchos::RCP<DRT::Discretization> targetdis, std::map<int, int> matmap,
        const bool isnurbsdis)
    {
      // now do the elements
      for (std::map<int, int>::iterator mapit = matmap.begin(); mapit != matmap.end(); ++mapit)
      {
        int target_id = mapit->second;
        CloneStrategy::check_material_type(target_id);
      }

      // prepare some variables we need
      int myrank = targetdis->Comm().MyPID();

      // construct new elements
      // The order of the elements might be different from that of the
      // source elements. We don't care. There are no dofs to these elements.
      std::set<int>::iterator it = roweleset_.begin();
      for (std::size_t i = 0; i < roweleset_.size(); ++i)
      {
        CORE::Elements::Element* sourceele = sourcedis->gElement(*it);

        std::string approxtype = "Polynomial";
        if (isnurbsdis)
        {
          if (sourceele->num_node() == 8)
          {
            approxtype = "NURBS8";
          }
          else if (sourceele->num_node() == 9)
          {
            approxtype = "NURBS9";
          }
          else if (sourceele->num_node() == 4)
          {
            approxtype = "NURBS4";
          }
          else if (sourceele->num_node() == 27)
          {
            approxtype = "NURBS27";
          }
          else if (sourceele->num_node() == 2)
          {
            approxtype = "NURBS2";
          }
          else if (sourceele->num_node() == 3)
          {
            approxtype = "NURBS3";
          }
          else
          {
            FOUR_C_THROW("unknown type of nurbs element\n");
          }
        }

        // create a new element of desired type with the same global element id
        Teuchos::RCP<CORE::Elements::Element> newele =
            CORE::COMM::Factory(eletype_[i], approxtype, *it, myrank);

        // get global node ids of source element
        std::vector<int> nids;
        nids.reserve(sourceele->num_node());
        transform(sourceele->Nodes(), sourceele->Nodes() + sourceele->num_node(),
            back_inserter(nids), std::mem_fn(&CORE::Nodes::Node::Id));

        // set the same global node ids to the new element
        newele->SetNodeIds(nids.size(), nids.data());

        // We need to set material and gauss points to complete element setup.
        // This is again really ugly as we have to extract the actual
        // element type in order to access the material property
        // note: SetMaterial() was reimplemented by the transport element!

        int src_matid = sourceele->Material()->Parameter()->Id();
        std::map<int, int>::iterator mat_iter = matmap.find(src_matid);
        if (mat_iter != matmap.end())
        {
          int tar_matid = mat_iter->second;
          CloneStrategy::set_element_data(newele, sourceele, tar_matid, isnurbsdis);

          // add new element to discretization
          targetdis->add_element(newele);
        }
        else
        {
          // before we stop, print the material id map
          std::cout << "Material map on PROC " << myrank << ":" << std::endl;
          for (mat_iter = matmap.begin(); mat_iter != matmap.end(); mat_iter++)
            std::cout << mat_iter->first << " -> " << mat_iter->second << std::endl;

          FOUR_C_THROW("no matching material ID (%d) in map", src_matid);
        }
        it++;
      }
      return;
    };  // create_elements

    /// create new elements from the condition and add them to the target discretization
    void create_elements_from_condition(
        const std::map<int, Teuchos::RCP<CORE::Elements::Element>>& sourceelements,
        DRT::Discretization& targetdis, const std::map<int, int>& matmap, const bool& isnurbsdis)
    {
      // now do the elements
      for (std::map<int, int>::const_iterator mapit = matmap.begin(); mapit != matmap.end();
           ++mapit)
      {
        int target_id = mapit->second;
        CloneStrategy::check_material_type(target_id);
      }

      // prepare some variables we need
      int myrank = targetdis.Comm().MyPID();

      // construct new elements
      // The order of the elements might be different from that of the
      // source elements. We don't care. There are no dofs to these elements.
      std::set<int>::iterator it = roweleset_.begin();
      for (std::size_t i = 0; i < roweleset_.size(); ++i)
      {
        std::map<int, Teuchos::RCP<CORE::Elements::Element>>::const_iterator src_ele_citer =
            sourceelements.find(*it);
        if (src_ele_citer == sourceelements.end())
          FOUR_C_THROW(
              "The source element %d could not be found in the source "
              "condition element map!",
              *it);

        CORE::Elements::Element* sourceele = src_ele_citer->second.get();
        if (sourceele == nullptr) FOUR_C_THROW("The sourceele pointer is nullptr!");

        std::string approxtype = "Polynomial";
        if (isnurbsdis)
        {
          if (sourceele->num_node() == 8)
          {
            approxtype = "NURBS8";
          }
          else if (sourceele->num_node() == 9)
          {
            approxtype = "NURBS9";
          }
          else if (sourceele->num_node() == 4)
          {
            approxtype = "NURBS4";
          }
          else if (sourceele->num_node() == 27)
          {
            approxtype = "NURBS27";
          }
          else if (sourceele->num_node() == 2)
          {
            approxtype = "NURBS2";
          }
          else if (sourceele->num_node() == 3)
          {
            approxtype = "NURBS3";
          }
          else
          {
            FOUR_C_THROW("unknown type of nurbs element\n");
          }
        }

        // get owner of source element
        const int sourceeleowner = sourceele->Owner();
        if (myrank != sourceeleowner)
          FOUR_C_THROW("roweleset_ should only contain my element gids!");

        // create a new element of desired type with the same global element id and same owner as
        // source element
        Teuchos::RCP<CORE::Elements::Element> newele =
            CORE::COMM::Factory(eletype_[i], approxtype, *it, sourceeleowner);

        // get global node ids of fluid element
        std::vector<int> nids;
        nids.reserve(sourceele->num_node());
        transform(sourceele->Nodes(), sourceele->Nodes() + sourceele->num_node(),
            back_inserter(nids), std::mem_fn(&CORE::Nodes::Node::Id));

        // set the same global node ids to the new element
        newele->SetNodeIds(nids.size(), nids.data());

        // We need to set material and gauss points to complete element setup.
        // This is again really ugly as we have to extract the actual
        // element type in order to access the material property
        // note: SetMaterial() was reimplemented by the transport element!
        Teuchos::RCP<CORE::MAT::Material> mat_ptr = sourceele->Material();
        /* Check if the material pointer is null. If necessary, try to cast
         * the condition element to a FaceElement and ask the parent element for
         * the material.                                                      */
        if (mat_ptr.is_null())
        {
          CORE::Elements::FaceElement* src_face_element =
              dynamic_cast<CORE::Elements::FaceElement*>(sourceele);
          if (src_face_element != nullptr) mat_ptr = src_face_element->parent_element()->Material();
        }
        // It is no FaceElement or the material pointer of the parent element is nullptr.
        if (mat_ptr.is_null()) FOUR_C_THROW("The condition element has no material!");

        int src_matid = mat_ptr->Parameter()->Id();
        std::map<int, int>::const_iterator mat_iter = matmap.find(src_matid);
        if (mat_iter != matmap.end())
        {
          int tar_matid = mat_iter->second;
          CloneStrategy::set_element_data(newele, sourceele, tar_matid, isnurbsdis);

          // add new element to discretization
          targetdis.add_element(newele);
        }
        else
        {
          // before we stop, print the material id map
          std::cout << "Material map on PROC " << myrank << ":" << std::endl;
          for (mat_iter = matmap.begin(); mat_iter != matmap.end(); mat_iter++)
            std::cout << mat_iter->first << " -> " << mat_iter->second << std::endl;

          FOUR_C_THROW("no matching material ID (%d) in map", src_matid);
        }
        it++;
      }
      return;
    }  // create_elements_from_condition

  };  // class DiscretizationCreator


  /// clone target discretization @p targetdis from a given source discretization @p sourcedis.
  /// The @p clonefieldmatmap is required from the global CloningMaterialMap.
  template <class CloneStrategy>
  void CloneDiscretization(Teuchos::RCP<DRT::Discretization> sourcedis,
      Teuchos::RCP<DRT::Discretization> targetdis,
      const std::map<std::pair<std::string, std::string>, std::map<int, int>>& clonefieldmatmap)
  {
    // access the communicator for time measurement
    const Epetra_Comm& comm = sourcedis->Comm();
    Teuchos::Time time("", true);

    // create target discretization using a given clone strategy
    {
      Teuchos::RCP<CORE::FE::DiscretizationCreator<CloneStrategy>> clonewizard =
          Teuchos::rcp(new CORE::FE::DiscretizationCreator<CloneStrategy>());

      std::map<int, int> matmap;
      clonewizard->create_clone_field_mat_map(matmap, *sourcedis, *targetdis, clonefieldmatmap);

      clonewizard->create_matching_discretization(sourcedis, targetdis, matmap);
    }
    if (comm.MyPID() == 0)
    {
      CORE::IO::cout << "Created discretization " << targetdis->Name()
                     << " as a clone of discretization " << sourcedis->Name() << " in...."
                     << time.totalElapsedTime(true) << " secs\n\n";
    }
    return;
  };  // CloneDiscretization

  /// clone target discretization @p targetdis from a given source discretization @p sourcedis
  /// based on conditions @p conds. The @p clonefieldmatmap is required from the global
  /// CloningMaterialMap.
  template <class CloneStrategy>
  void CloneDiscretizationFromCondition(const DRT::Discretization& sourcedis,
      DRT::Discretization& targetdis, const std::vector<CORE::Conditions::Condition*>& conds,
      const std::map<std::pair<std::string, std::string>, std::map<int, int>>& clonefieldmatmap)
  {
    const DRT::Discretization* sourcedis_ptr = dynamic_cast<const DRT::Discretization*>(&sourcedis);
    if (sourcedis_ptr == nullptr) FOUR_C_THROW("Cast of the source discretization failed!");
    DRT::Discretization* targetdis_ptr = dynamic_cast<DRT::Discretization*>(&targetdis);
    if (targetdis_ptr == nullptr) FOUR_C_THROW("Cast of the target discretization failed!");

    // access the communicator for time measurement
    const Epetra_Comm& comm = sourcedis_ptr->Comm();
    Teuchos::Time time("", true);

    // create target discretization using a given clone strategy
    {
      Teuchos::RCP<CORE::FE::DiscretizationCreator<CloneStrategy>> clonewizard =
          Teuchos::rcp(new CORE::FE::DiscretizationCreator<CloneStrategy>());

      std::map<int, int> matmap;
      clonewizard->create_clone_field_mat_map(
          matmap, *sourcedis_ptr, *targetdis_ptr, clonefieldmatmap);

      clonewizard->create_matching_discretization_from_condition(
          *sourcedis_ptr, conds, *targetdis_ptr, matmap);
    }
    if (comm.MyPID() == 0)
    {
      CORE::IO::cout << "Created discretization " << targetdis_ptr->Name()
                     << " as a clone from the condition(s) with ID(s)=";
      for (unsigned int i = 0; i < conds.size(); ++i) CORE::IO::cout << " " << conds[i]->Id();
      CORE::IO::cout << " of the discretization " << sourcedis_ptr->Name() << " in...."
                     << time.totalElapsedTime(true) << " secs\n\n";
    }
    return;
  };  // CloneDiscretizationFromCondition

  /// clone target discretization @p targetdis from a given source discretization @p sourcedis
  /// based on the name of a condition @p condname. The @p clonefieldmatmap is required from the
  /// global CloningMaterialMap.
  template <class CloneStrategy>
  void CloneDiscretizationFromCondition(const DRT::Discretization& sourcedis,
      DRT::Discretization& targetdis, const std::string& condname,
      const std::map<std::pair<std::string, std::string>, std::map<int, int>>& clonefieldmatmap)
  {
    // access the communicator for time measurement
    const Epetra_Comm& comm = sourcedis.Comm();
    Teuchos::Time time("", true);

    // create target discretization using a given clone strategy
    {
      Teuchos::RCP<CORE::FE::DiscretizationCreator<CloneStrategy>> clonewizard =
          Teuchos::rcp(new CORE::FE::DiscretizationCreator<CloneStrategy>());

      std::map<int, int> matmap;
      clonewizard->create_clone_field_mat_map(matmap, sourcedis, targetdis, clonefieldmatmap);

      clonewizard->create_matching_discretization_from_condition(
          sourcedis, condname, targetdis, matmap);
    }
    if (comm.MyPID() == 0)
    {
      CORE::IO::cout << "Created discretization " << targetdis.Name()
                     << " as a clone from the condition \"" << condname.c_str()
                     << "\" of the discretization " << sourcedis.Name() << " in...."
                     << time.totalElapsedTime(true) << " secs\n\n";
    }
    return;
  };  // CloneDiscretizationFromCondition

  //! Return valid cloning material map input lines.
  std::vector<INPUT::LineDefinition> valid_cloning_material_map_lines();

}  // namespace CORE::FE

FOUR_C_NAMESPACE_CLOSE

#endif

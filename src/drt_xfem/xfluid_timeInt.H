/*----------------------------------------------------------------------*/
/*! \file

\brief provides the xfluid timeIntegration,
       maps vectors from old interface position to vectors at new interface position,
       determines the reconstruction method for missing and unreasonable ghost and standard values

\level 2

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249

*/
/*----------------------------------------------------------------------*/


#ifndef XFLUID_TIMEINTEGRATION_H_
#define XFLUID_TIMEINTEGRATION_H_


#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_cut/cut_utils.H"
#include "../drt_cut/cut_node.H"


class Epetra_Map;
class Epetra_Vector;



namespace DRT
{
  class Discretization;
  class Element;
  class PackBuffer;
  class Exporter;
  class Node;
}  // namespace DRT

namespace GEO
{
  class CutWizard;

  namespace CUT
  {
    class SideHandle;
  }
}  // namespace GEO

namespace LINALG
{
  class MapExtractor;
}

namespace Teuchos
{
  class ParameterList;
}

namespace XFEM
{
  class XFEMDofSet;
  class ConditionManager;

  /*!
  \brief this class is the basic TIMEINT class for the projection, adaption or
         something else in XFEM-problems between consecutive time steps
   */
  class XFluidTimeInt
  {
   public:
    //! constructor
    explicit XFluidTimeInt(
        const bool is_newton_increment_transfer,  /// monolithic newton increment transfer or time
                                                  /// step transfer?
        const Teuchos::RCP<DRT::Discretization>& dis,                   /// discretization
        const Teuchos::RCP<XFEM::ConditionManager>& condition_manager,  /// condition manager
        const Teuchos::RCP<GEO::CutWizard>& wizard_old,                 /// cut wizard at t^n
        const Teuchos::RCP<GEO::CutWizard>& wizard_new,                 /// cut wizard at t^(n+1)
        const Teuchos::RCP<XFEM::XFEMDofSet>& dofset_old,               /// XFEM dofset at t^n
        const Teuchos::RCP<XFEM::XFEMDofSet>& dofset_new,               /// XFEM dofset at t^(n+1)
        const INPAR::XFEM::XFluidTimeIntScheme xfluid_timintapproach,   /// xfluid_timintapproch
        std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt>>&
            node_to_reconstr_method,  /// reconstruction map for nodes and its dofsets
        std::map<INPAR::XFEM::XFluidTimeInt, std::map<int, std::set<int>>>&
            reconstr_method_to_node,  /// inverse reconstruction map for nodes and its dofsets
        const int step,               /// global time step
        const bool xfluid_timint_check_interfacetips  /// check interfacetips?
    );

    //! destructor
    ~XFluidTimeInt() { return; }

    /// set and print reconstruction status for nodes
    void SetAndPrintStatus(const bool screenout);


    /// transfer standard and ghost dofs to new map as far as possible and mark dofs for
    /// reconstruction
    void TransferDofsToNewMap(
        const std::vector<Teuchos::RCP<const Epetra_Vector>>&
            oldRowStateVectors,  /// row map based vectors w.r.t old interface position
        const std::vector<Teuchos::RCP<Epetra_Vector>>&
            newRowStateVectors,  /// row map based vectors w.r.t new interface position
        const Teuchos::RCP<std::set<int>>
            dbcgids  /// set of dof gids that must not be changed by ghost penalty reconstruction
    );

    /// transfer standard and ghost dofs to new map as far as possible and mark dofs for
    /// reconstruction for given vector of node gids
    void TransferDofsToNewMap(
        const std::vector<Teuchos::RCP<const Epetra_Vector>>&
            oldRowStateVectors,  /// row map based vectors w.r.t old interface position
        const std::vector<Teuchos::RCP<Epetra_Vector>>&
            newRowStateVectors,  /// row map based vectors w.r.t new interface position
        const Teuchos::RCP<std::set<int>>
            dbcgids,  /// set of dof gids that must not be changed by ghost penalty reconstruction
        const std::vector<int>& node_gids  /// vector of node gids
    );

    /// get for each type of reconstruction method the number of nodes for that this method has to
    /// be applied on this proc
    std::map<INPAR::XFEM::XFluidTimeInt, int>& Get_Reconstr_Counts() { return reconstr_counts_; };

    /// get for each type of reconstruction method the node ids with corresponding dof its for that
    /// this method has to be applied on this proc
    std::map<int, std::set<int>>& Get_NodeToDofMap_For_Reconstr(
        INPAR::XFEM::XFluidTimeInt reconstr);

    /// get permutation map for ghost dofs
    Teuchos::RCP<std::map<int, int>> GetPermutationMap() { return permutation_map_; };

    /// timint output for reconstruction methods
    void Output();

   private:
    /// transfer standard and ghost dofs to new map as far as possible and mark dofs for
    /// reconstruction for a given node gid
    void TransferNodalDofsToNewMap(
        const std::vector<Teuchos::RCP<const Epetra_Vector>>&
            oldRowStateVectors,  /// row map based vectors w.r.t old interface position
        const std::vector<Teuchos::RCP<Epetra_Vector>>&
            newRowStateVectors,  /// row map based vectors w.r.t new interface position
        const Teuchos::RCP<std::set<int>>
            dbcgids,  /// set of dof gids that must not be changed by ghost penalty reconstruction
        int gid       /// nodal gid
    );

    /// returns matching std::string for each reconstruction method
    std::string MapMethodEnumToString(const enum INPAR::XFEM::XFluidTimeInt term);

    /// all surrounding elements non-intersected?
    bool NonIntersectedElements(DRT::Node* n, const Teuchos::RCP<GEO::CutWizard> wizard);

    /// find all ghost dofsets around this node and its std-dofset
    void FindSurroundingGhostDofsets(
        std::map<int, std::set<int>>&
            ghostDofsets,       /// surrounding ghost dofsets to be filled, map of ghost nodes and
                                /// correponding ghost nds index w.r.t given std nodal dofset
        const DRT::Node* node,  /// node
        const int nds_new       /// dofset of node used for finding the surrounding ghost dofs
    );

    /// copy dofs from old vectors to new vector for all row vectors
    void CopyDofs(const DRT::Node* node,          /// drt node
        const int nds_new,                        /// nodal dofset at t^(n+1)
        const int nds_old,                        /// nodal dofset at t^n
        const INPAR::XFEM::XFluidTimeInt method,  /// reconstruction method
        const std::vector<Teuchos::RCP<Epetra_Vector>>&
            newRowStateVectors,  /// row map based state vectors at t^(n+1)
        const std::vector<Teuchos::RCP<const Epetra_Vector>>&
            oldRowStateVectors,                    /// row map based state vectors at t^n
        const Teuchos::RCP<std::set<int>> dbcgids  /// set of DBC global ids
    );

    /// mark nodal dofs of vector w.r.t new interface position for reconstruction
    void MarkDofs(const DRT::Node* node,  /// drt node
        const int nds_new,                /// nodal dofset at t^(n+1)
        const std::vector<Teuchos::RCP<Epetra_Vector>>&
            newRowStateVectors,                   /// row map based state vectors at t^(n+1)
        const INPAR::XFEM::XFluidTimeInt method,  /// reconstruction method
        const Teuchos::RCP<std::set<int>>
            dbcgids  /// set of dof gids that must not be changed by ghost penalty reconstruction
    );

    /// mark one specific nodal dofset with used for export
    void MarkDofsForExport(const int nid,  /// node id
        const int dofset,                  /// ghost dofset number
        const INPAR::XFEM::XFluidTimeInt
            method  /// reconstruction method used for marking the nodal dofset
    );

    /// set the reconstruction method for current nodal dofset, return if set
    bool SetReconstrMethod(const DRT::Node* node,  /// drt node
        const int nds_new,                         /// nodal dofset at t^(n+1)
        const INPAR::XFEM::XFluidTimeInt method    /// which type of reconstruction method
    );

    /// identify cellsets at time t^n with cellsets at time t^(n+1)
    int IdentifyOldSets(const GEO::CUT::Node* n_old,  /// node w.r.t to old wizard
        const GEO::CUT::Node* n_new,                  /// node w.r.t to new wizard
        const std::vector<Teuchos::RCP<GEO::CUT::NodalDofSet>>&
            dof_cellsets_old,  /// all dofcellsets at t^n
        const GEO::CUT::NodalDofSet*
            cell_set_new  /// dofcellset at t^(n+1) which has to be identified
    );

    /// special check if the node slides along the cut surface
    bool SpecialCheck_SlidingOnSurface(bool& changed_side,
        const GEO::CUT::Node* n_old,  /// node w.r.t to old wizard
        const GEO::CUT::Node* n_new   /// node w.r.t to new wizard
    );

    /// check if the node has changed the side w.r.t identified sides at t^n and t^(n+1), return if
    /// check was successful
    bool SpecialCheck_InterfaceTips(bool& changed_side,  /// did the node change the side ?
        std::vector<int>& identified_sides,              /// side Id of identified side
        const GEO::CUT::Node* n_old,                     /// node w.r.t to old wizard
        const GEO::CUT::Node* n_new                      /// node w.r.t to new wizard
    );

    /// special check for level-set based interface tips, note: currently empty, can be filled if
    /// necessary
    bool SpecialCheck_InterfaceTips_Levelset(bool& changed_side  /// did the node change the side ?
    );

    /// special check for mesh based interface tips
    bool SpecialCheck_InterfaceTips_SpaceTime(
        bool& changed_side,  /// did the node change the side ?
        DRT::Element* side, const int coup_sid,
        const LINALG::Matrix<3, 1>& n_coord  /// node coodinates
    );

    /// check if the node is within the space time side
    template <DRT::Element::DiscretizationType side_distype,
        DRT::Element::DiscretizationType space_time_distype>
    bool WithinSpaceTimeSide(bool& within_space_time_side,  /// within the space time side
        DRT::Element* side, const int coup_sid,
        const LINALG::Matrix<3, 1>& n_coord  /// node coodinates
    );

    /// check the volume of the space time side, distorted space-time side ?
    template <DRT::Element::DiscretizationType space_time_distype, const int numnode_space_time>
    bool CheckSTSideVolume(const LINALG::Matrix<3, numnode_space_time>& xyze_st);

    /// export data about reconstruction method to neighbor proc and receive data from previous proc
    void ExportMethods(
        const std::vector<Teuchos::RCP<Epetra_Vector>>&
            newRowStateVectors,  /// row map based vectors w.r.t new interface position
        const Teuchos::RCP<std::set<int>>
            dbcgids  /// set of dof gids that must not be changed by ghost penalty reconstruction
    );

    /*!
    \brief Basic function sending data to destination and receiving data from source
     */
    void sendData(DRT::PackBuffer& dataSend,  //!< pack buffer
        int& dest,                            //!< destination proc
        int& source,                          //!< source proc
        std::vector<char>& dataRecv           //!< received data
        ) const;

    const bool is_newton_increment_transfer_;  /// monolithic newton increment transfer or time step
                                               /// transfer?

    Teuchos::RCP<DRT::Discretization> dis_;  /// background  discretization

    Teuchos::RCP<XFEM::ConditionManager> condition_manager_;  /// condition manager

    Teuchos::RCP<GEO::CutWizard> wizard_old_;  /// old cut wizard w.r.t old interface position
    Teuchos::RCP<GEO::CutWizard> wizard_new_;  /// new cut wizard w.r.t new interface position

    Teuchos::RCP<XFEM::XFEMDofSet> dofset_old_;  /// old XFEM dofset w.r.t old interface position
    Teuchos::RCP<XFEM::XFEMDofSet> dofset_new_;  /// new XFEM dofset w.r.t new interface position

    // current processor id and number of procs
    int myrank_;
    int numproc_;

    INPAR::XFEM::XFluidTimeIntScheme timeint_scheme_;  /// which type of time integration scheme
                                                       /// shall be used for reconstructing values

    std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt>>&
        node_to_reconstr_method_;  /// map of nodeIds and reconstruction method for all nodal
                                   /// dofsets

    std::map<INPAR::XFEM::XFluidTimeInt, std::map<int, std::set<int>>>&
        reconstr_method_to_node_;  /// map of reconstruction method to map of node-id and nodal
                                   /// dofset-id (inverse to node_to_reconstr_method)

    std::map<INPAR::XFEM::XFluidTimeInt, int>
        reconstr_counts_;  /// counts the number of dofsets that have to be reconstructed using each
                           /// method

    std::map<int, std::map<int, int>>
        dofset_marker_export_;  /// std::map<nid, std::map<dofset_number,method> > that contains
                                /// marker for reconstruction method of single nodal dofsets used
                                /// for export in parallel

    // global timestep
    int step_;

    // name check interfacetips in timeintegration
    bool xfluid_timint_check_interfacetips_;

    Teuchos::RCP<std::map<int, int>> permutation_map_;
  };

}  // namespace XFEM

#endif /*XFLUID_TIMEINTEGRATION_H_*/

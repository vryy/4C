/*!-----------------------------------------------------------------------------------------------*
\file xfluid_timeInt.H

\brief provides the xfluid timeIntegration,
       maps vectors from old interface position to vectors at new interface position,
       determines the reconstruction method for missing and unreasonable ghost and standard values

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*------------------------------------------------------------------------------------------------*/


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
}
// namespace CUT
// {
//   class Cmp;
//   class Node;
// }
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
  class FluidWizard;

/*!
\brief this class is the basic TIMEINT class for the projection, adaption or
       something else in XFEM-problems between consecutive time steps
 */
class XFluidTimeInt
{
public:

  //! constructor
  explicit XFluidTimeInt(
      const Teuchos::RCP<DRT::Discretization> dis,                               /// discretization
      const Teuchos::RCP<XFEM::FluidWizard>   wizard_old,                        /// fluid wizard at t^n
      const Teuchos::RCP<XFEM::FluidWizard>   wizard_new,                        /// fluid wizard at t^(n+1)
      const Teuchos::RCP<XFEM::FluidDofSet>   dofset_old,                        /// dofset at t^n
      const Teuchos::RCP<XFEM::FluidDofSet>   dofset_new,                        /// dofset at t^(n+1)
      const Teuchos::ParameterList&           params,                            /// parameter list
      INPAR::XFEM::XFluidTimeIntScheme        xfluid_timintapproach,             /// xfluid_timintapproch
      std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >& reconstr_method,  /// reconstruction map for nodes and its dofsets
      int                                     step                               /// global time step
  );

  //! destructor
  ~XFluidTimeInt()
  {
    return;
  }

  /// set and print reconstruction status for nodes
  void SetAndPrintStatus(bool screenout);


  /// transfer standard and ghost dofs to new map as far as possible and mark dofs for reconstruction
  void TransferDofsToNewMap(
      const Epetra_Map&                       olddofrowmap,                     /// dof row map w.r.t old interface position
      const Epetra_Map&                       olddofcolmap,                     /// dof col map w.r.t old interface position
      std::vector<RCP<const Epetra_Vector> >& oldRowStateVectors,               /// row map based vectors w.r.t old interface position
      std::vector<RCP<Epetra_Vector> >&       newRowStateVectors,               /// row map based vectors w.r.t new interface position
      std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >& reconstr_method, /// reconstruction map for nodes and its dofsets
      Teuchos::RCP<std::set<int> >            dbcgids                           /// set of dof gids that must not be changed by ghost penalty reconstruction
      );


  /// get for each type of reconstruction method the number of nodes for that this method has to be applied on this proc
  std::map<INPAR::XFEM::XFluidTimeInt,int>& Get_Reconstr_Counts(){return reconstr_counts_;};

  /// timint output for reconstruction methods
  void Output();

  private:

  /// returns matching std::string for each reconstruction method
  std::string MapMethodEnumToString
  (
     const enum INPAR::XFEM::XFluidTimeInt term
  );

  /// all surrounding elements uncut?
  bool UncutEles(GEO::CUT::Node* n);

  /// find all ghost dofsets around this node and its std-dofset
  void FindSurroundingGhostDofsets(
      std::map<int,std::vector<int> >& ghostDofsets,   /// surrounding ghost dofsets to be filled
      DRT::Node*                       node,           /// node
      int                              nds_new         /// dofset of node used for finding the surrounding ghost dofs
  );

  /// copy dofs from old vectors to new vector for all row vectors
  void CopyDofs(
      DRT::Node*                                 node,               /// drt node
      const int                                  nds_new,            /// nodal dofset at t^(n+1)
      const int                                  nds_old,            /// nodal dofset at t^n
      INPAR::XFEM::XFluidTimeInt                 method,               /// reconstruction method
      std::vector<RCP<Epetra_Vector> >&          newRowStateVectors, /// row map based state vectors at t^(n+1)
      std::vector<RCP<const Epetra_Vector> >&    oldRowStateVectors, /// row map based state vectors at t^n
      Teuchos::RCP<std::set<int> >               dbcgids             /// set of DBC global ids
  );

  /// mark nodal dofs of vector w.r.t new interface position for reconstruction
  void MarkDofs(
      DRT::Node*                          node,                 /// drt node
      const int                           nds_new,              /// nodal dofset at t^(n+1)
      std::vector<RCP<Epetra_Vector> >&   newRowStateVectors,   /// row map based state vectors at t^(n+1)
      INPAR::XFEM::XFluidTimeInt          method,               /// reconstruction method
      Teuchos::RCP<std::set<int> >        dbcgids               /// set of dof gids that must not be changed by ghost penalty reconstruction
  );

  /// mark one specific nodal dofset with used for export
  void MarkDofsForExport(
      int                        nid,         /// node id
      int                        dofset,      /// ghost dofset number
      INPAR::XFEM::XFluidTimeInt method       /// reconstruction method used for marking the nodal dofset
      );

  /// set the reconstruction method for current nodal dofset, return if set
  bool SetReconstrMethod(
      DRT::Node*                     node,                 /// drt node
      const int                      nds_new,              /// nodal dofset at t^(n+1)
      INPAR::XFEM::XFluidTimeInt     method                /// which type of reconstruction method
  );

  /// find the standard dofset, return the dofset number of std dofset
  int FindStdDofSet(
      GEO::CUT::Node*                                                               node,          /// cut node
      const std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp > >& dof_cell_sets  /// dofcellsets of node
  );

  /// is this node a standard or ghost node w.r.t current set
  bool Is_Std_CellSet(
      GEO::CUT::Node*                                  node,     /// cut node
      const std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp>&  cell_set  /// set of volumecells
  );

  /// identify cellsets at time t^n with cellsets at time t^(n+1)
  void IdentifyOldSets(
      int &                                                          nds_old,            /// set identified nodal dofset at t^n
      std::vector<int> &                                             identified_sides,   /// set identified using sides (side-Ids)
      const std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp> >&  dof_cellsets_old,   /// all dofcellsets at t^n
      const std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp>& cell_set_new        /// dofcellset at t^(n+1) which has to be identified
  );

  /// check if the node has changed the side w.r.t identified sides at t^n and t^(n+1), return if check was successful
  bool CheckChangingSide(
      bool&                                                          changed_side,        /// did the node change the side ?
      std::vector<int> &                                             identified_sides,    /// side Id of identified side
      GEO::CUT::Node *                                               n_old,               /// node w.r.t to old wizard
      GEO::CUT::Node *                                               n_new,               /// node w.r.t to new wizard
      const std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp>& cell_set_old,        /// dofcellset at t^n
      const std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp>& cell_set_new         /// dofcellset at t^(n+1)
  );

  /// check if the node is within the space time side
  template<DRT::Element::DiscretizationType side_distype, DRT::Element::DiscretizationType space_time_distype>
  bool WithinSpaceTimeSide(
      bool&                   within_space_time_side,    /// within the space time side
      GEO::CUT::SideHandle*   side_old,                  /// side w.r.t old interface
      GEO::CUT::SideHandle*   side_new,                  /// side w.r.t new interface
      LINALG::Matrix<3,1>&    n_coord                    /// node coodinates
  );

  /// check the volume of the space time side, distorted space-time side ?
  template<DRT::Element::DiscretizationType space_time_distype, const int numnode_space_time>
  bool CheckSTSideVolume( LINALG::Matrix<3,numnode_space_time>&  xyze_st );

  /// export data about reconstruction method to neighbor proc and receive data from previous proc
  void ExportMethods(
      std::vector<RCP<Epetra_Vector> >&       newRowStateVectors,               /// row map based vectors w.r.t new interface position
      Teuchos::RCP<std::set<int> >            dbcgids                           /// set of dof gids that must not be changed by ghost penalty reconstruction
  );

  /*!
  \brief Basic function sending data to destination and receiving data from source
   */
  void sendData(
      DRT::PackBuffer& dataSend,                        //!< pack buffer
      int& dest,                                        //!< destination proc
      int& source,                                      //!< source proc
      std::vector<char>& dataRecv                       //!< received data
  ) const;

  Teuchos::RCP<DRT::Discretization> dis_;         /// background  discretization

  Teuchos::RCP<XFEM::FluidWizard> wizard_old_;    /// old fluid wizard w.r.t old interface position
  Teuchos::RCP<XFEM::FluidWizard> wizard_new_;    /// new fluid wizard w.r.t new interface position

  Teuchos::RCP<XFEM::FluidDofSet> dofset_old_;    /// old dofset w.r.t old interface position
  Teuchos::RCP<XFEM::FluidDofSet> dofset_new_;    /// new dofset w.r.t new interface position

  Teuchos::ParameterList params_;                 /// parameter list

  // current processor id and number of procs
  int myrank_;
  int numproc_;

  INPAR::XFEM::XFluidTimeIntScheme timeint_scheme_; /// which type of time integration scheme shall be used for reconstructing values

  std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >& reconstr_method_; /// map of nodeIds and reconstruction method for all nodal dofsets

  std::map<INPAR::XFEM::XFluidTimeInt,int> reconstr_counts_; /// counts the number of dofsets that have to be reconstructed using each method

  std::map<int, std::map<int,int> > dofset_marker_export_; /// std::map<nid, std::map<dofset_number,method> > that contains marker for reconstruction method of single nodal dofsets used for export in parallel

  // global timestep
  int step_;

};

} // namespace XFEM

#endif /*XFLUID_TIMEINTEGRATION_H_*/


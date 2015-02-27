/*!-----------------------------------------------------------------------------------------------*
\file xfluidfluid_timeInt.H

\brief class provides the following functionalities for
       xfluid-fluid time integration:
       - mapping of state vectors from old intersection states (interface positions) to new ones
       - reconstruction of missing / unreasonable ghost values
       - subsequent enforcement of incompressibility on xfluid element patch
         after projection step

<pre>
Maintainer:
             Raffaela Kruse
             kruse@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef XFLUIDFLUID_TIMEINTEGRATION_H_
#define XFLUIDFLUID_TIMEINTEGRATION_H_


#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_lib/drt_elementtype.H"


class Epetra_Map;
class Epetra_Vector;

namespace GEO
{
  class CutWizard;
  class SearchTree;
}

namespace DRT
{
  class Discretization;
  class Element;
  class PackBuffer;
  class Exporter;
  class Node;
}

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

/*!
\brief basic time integration class for xfluid-fluid problems
 */
class XFluidFluidTimeIntegration
{
public:

  //! constructor
  /*!
   * \param bgdis xfem fluid discretization
   * \param embdis embedded fluid discretization
   * \param cut-wizard
   * \param step time step
   * \param xfem_timeintapproach type of xfluid-fluid time integration approach (projection,etc.)
   * \param params parameter list
   */
  explicit XFluidFluidTimeIntegration(
    const Teuchos::RCP<DRT::Discretization> &  bgdis,
    const Teuchos::RCP<DRT::Discretization> &  embdis,
    Teuchos::RCP<GEO::CutWizard>               wizard,
    int                                        step,
    enum INPAR::XFEM::XFluidFluidTimeInt       xfem_timeintapproach,
    const Teuchos::ParameterList&              params
    );

  //! destructor
  ~XFluidFluidTimeIntegration()
  {
  }


  /*!
   * \brief Save the old state (maps of node-ids to dof-ids) and create a new one
   * \return true, if the maps have changed compared to the old state
   */
  bool SaveBgNodeMapsAndCreateNew(Teuchos::RCP<GEO::CutWizard> wizard);

  //! Create the background node maps
  void CreateBgNodeMapsForRestart(Teuchos::RCP<GEO::CutWizard> wizard);


  /*!
   * \brief Determine state vectors for background fluid nodes (either by projection from emb.
   * fluid discretization or by keeping the ghost values)
   */
  void SetNewBgStateVectors(
    const Teuchos::RCP<Epetra_Vector>        bgstate_velnp,
    const Teuchos::RCP<Epetra_Vector>        bgstate_veln,
    const Teuchos::RCP<Epetra_Vector>        bgstate_accn,
    const Teuchos::RCP<const Epetra_Vector>  bgstaten_velnp,
    const Teuchos::RCP<const Epetra_Vector>  bgstaten_veln,
    const Teuchos::RCP<const Epetra_Vector>  bgstaten_accn,
    const Teuchos::RCP<const Epetra_Vector>  embstate_velnp,
    const Teuchos::RCP<const Epetra_Vector>  embstate_veln,
    const Teuchos::RCP<const Epetra_Vector>  embstate_accn,
    const Teuchos::RCP<const Epetra_Vector>  aledispn);

  /*!
   * \brief Determine embedded state vectors after ALE-relaxation in interpolation-based fluid-fluid FSI
   */
  void SetNewEmbStateVectors(
      const Teuchos::RCP<Epetra_Vector>        embstate_velnp,
      const Teuchos::RCP<Epetra_Vector>        embstate_veln,
      const Teuchos::RCP<Epetra_Vector>        embstate_accn,
      const Teuchos::RCP<const Epetra_Vector>  embstaten_velnp,
      const Teuchos::RCP<const Epetra_Vector>  embstaten_veln,
      const Teuchos::RCP<const Epetra_Vector>  embstaten_accn,
      const Teuchos::RCP<const Epetra_Vector>  bgstaten_velnp,
      const Teuchos::RCP<const Epetra_Vector>  bgstaten_veln,
      const Teuchos::RCP<const Epetra_Vector>  bgstaten_accn,
      const Teuchos::RCP<const Epetra_Vector>  aledispnp,
      const Teuchos::RCP<const Epetra_Vector>  aledispn);

  /*!
   * \brief Subsequent enforcement of incompressibility on background element patch after projection
   * by solution of an optimization problem
   */
  void EnforceIncompAfterProjection(
      Teuchos::RCP<GEO::CutWizard>                      wizard,         ///< cut wizard
      Teuchos::RCP<GEO::CutWizard>                      wizard_n,       ///< cut wizard from timestep n
      Teuchos::RCP<Epetra_Vector> &                     bgstate_velnp,  ///< background fluid velocity at timestep n+1
      Teuchos::RCP<Epetra_Vector> &                     bgstate_veln,   ///< background fluid velocity at timestep n
      const Teuchos::RCP<const LINALG::MapExtractor> &  dbcmaps         ///< background fluid dirichlet map extractor
      );

  //! @name Gmsh-output
  //@{
  void GmshOutput();

  void GmshOutputForInterpolateFSI(const Teuchos::RCP<const Epetra_Vector> & aledispnp,
                                   const Teuchos::RCP<const Epetra_Vector> & aledispn);
  //@}

private:

  /*!
   * \brief Container class for source state vectors, accessed by encapsulating class
   */
  class SourceState
  {

  public:

    //! ctor
    SourceState(
      Teuchos::RCP<const Epetra_Vector> velnp,
      Teuchos::RCP<const Epetra_Vector> veln  = Teuchos::null,
      Teuchos::RCP<const Epetra_Vector> accn  = Teuchos::null
    ) :
      velnp_(velnp),
      veln_(veln),
      accn_(accn)
    {
    }

    //! @name State vectors
    //@{
    Teuchos::RCP<const Epetra_Vector> velnp_;
    Teuchos::RCP<const Epetra_Vector> veln_;
    Teuchos::RCP<const Epetra_Vector> accn_;
    //@}

  };

  /*!
   * \brief Container class for target state vectors, accessed by encapsulating class
   */
  class TargetState
  {

  public:

    //! ctor
    TargetState(
      Teuchos::RCP<Epetra_Vector> velnp,
      Teuchos::RCP<Epetra_Vector> veln  = Teuchos::null,
      Teuchos::RCP<Epetra_Vector> accn  = Teuchos::null
    ) :
      velnp_(velnp),
      veln_(veln),
      accn_(accn)
    {
    }

    //! @name State vectors
    //@{
    Teuchos::RCP<Epetra_Vector> velnp_;
    Teuchos::RCP<Epetra_Vector> veln_;
    Teuchos::RCP<Epetra_Vector> accn_;
    //@}

  };

  /// Create maps of bgfluid's nodes needed for time integration
  void CreateBgNodeMaps(Teuchos::RCP<GEO::CutWizard> wizard);

  /// Save the state (maps of node-ids to dof-ids)
  void SaveBgNodeMaps();

  /*!
   * \brief Determine state vectors for background fluid nodes (either by projection from emb.
   * fluid discretization or by keeping the ghost values)
   */
  void SetNewBgStateVectors(const TargetState &                       bgfluid_state,
                            const SourceState &                       bgfluid_state_n,
                            const SourceState &                       embfluid_state,
                            const Teuchos::RCP<const Epetra_Vector> & aledispn);

  /// set values of background state vector (do not do the projection if ghost-values are available)
  void SetNewBgStateVectorKeepGhostValues(const TargetState &                       bgfluid_state,
                                          const SourceState &                       bgfluid_state_n,
                                          const SourceState &                       embfluid_state,
                                          const Teuchos::RCP<const Epetra_Vector> & aledispn);

  /// set values of background state vector (always project values from embedded discretization if no history values are available)
  void SetNewBgStateVectorFullProjection( const TargetState &                       bgfluid_state,
                                          const SourceState &                       bgfluid_state_n,
                                          const SourceState &                       embfluid_state,
                                          const Teuchos::RCP<const Epetra_Vector> & aledispn);

  /*!
   * \brief Determine state vectors for embedded fluid nodes (either by projection from emb.
   * mesh at previous location or from background fluid nodes)
   */
  void SetNewEmbStateVectors( const TargetState &                       embfluid_state,
                              const SourceState &                       embfluid_state_n,
                              const SourceState &                       bgfluid_state_n,
                              const Teuchos::RCP<const Epetra_Vector> & aledispnp,
                              const Teuchos::RCP<const Epetra_Vector> & aledispn);

  /// Transfer data in the round robin communication pattern
  void CommunicateNodes(std::vector<LINALG::Matrix<3,1> >  &      bgnodes_xyz,
                        std::vector<LINALG::Matrix<12,1> > &      interpolated_vecs,
                        std::vector<int>  &                       bgnodes_nohistory,
                        const Teuchos::RCP<const Epetra_Vector> & aledispn,
                        const TargetState &                       bgfluid_state,
                        const SourceState &                       bgfluid_state_n,
                        const SourceState &                       embfluid_state,
                        const Teuchos::RCP<DRT::Discretization>&  discret);

  /// receive a block in the round robin communication pattern
  void ReceiveBlock(std::vector<char>  & rblock,
                    DRT::Exporter      & exporter,
                    MPI_Request        & request);

  /// send a block in the round robin communication pattern
  void SendBlock(std::vector<char>  & sblock,
                 DRT::Exporter      & exporter,
                 MPI_Request        & request );

  /// pack values in the round robin communication pattern
  void PackValues(std::vector<LINALG::Matrix<3,1> >  & bgnodes_xyz,
                  std::vector<LINALG::Matrix<12,1> > & interpolatedvec,
                  std::vector<int>                   & bgnodes_nohistory,
                  std::vector<int>                   & have_values,
                  std::vector<char>                  & sblock);

  /// determine the embedded element, that previously covered the nodes and interpolate values
  void FindEmbEleAndInterpolateValues(
    std::vector<LINALG::Matrix<3,1> >        & nodes_xyz,
    std::vector<LINALG::Matrix<12,1> >       & interpolated_vecs,
    std::vector<int>                         & have_values,
    const Teuchos::RCP<const Epetra_Vector>  & aledispn,
    const SourceState                        & embfluid_state);


  /// setup a search tree for the embedded (ALE-) discretization
  void SetupSearchTree(const Teuchos::RCP<const Epetra_Vector>  & aledispn);

  /// determine the search radius for the search tree
  template<DRT::Element::DiscretizationType distype>
  void FindSearchRadius();

  /// Write the values of node from old to new background fluid state vector
  void WriteValuesToBgStateVector(const DRT::Node*                           bgnode,
                                  std::vector<int>                           gdofs_n,
                                  const Teuchos::RCP<Epetra_Vector> &        bgstatevec_np,
                                  const Teuchos::RCP<const Epetra_Vector> &  bgstatevec_n);

  /*!
   * \brief check whether the unknown background fluid node is covered by the embedded element -
   * if yes, interpolate the state vectors and return true
   *
   * \param (in) src_ele the projection source (embedded) element, the node is possibly covered by
   * \param (in) src_state src_state container with the state vectors from projection source
   * \param (in) src_disp displacements of the projection source (ALE-)element
   * \param src_dis source discretization (embedded ALE-fluid, this is checked!)
   * \param bgnode_xyz the coordinates of the xfem fluid node without values
   * \param interpolatedvec the vector with the values after state vector interpolation
   * \return true, if the node is covered by the embedded element, false otherwise
   */
  template<DRT::Element::DiscretizationType distype>
  bool ComputeSpatialToElementCoordAndProject(
    const DRT::Element *                      src_ele,
    const SourceState &                       src_state,
    const Teuchos::RCP<const Epetra_Vector> & src_disp,
    const DRT::Discretization &               src_dis,
    const LINALG::Matrix<3,1> &               bgnode_xyz,
    LINALG::Matrix<12,1> &                    interpolatedvec
  );

  //! @name Methods for subsequent enforcement of incompressibility on background element patch after projection
  //@{
  void BuildElementPatchForIncompOpt( Teuchos::RCP<GEO::CutWizard>                      wizard_n,
                                      Teuchos::RCP<GEO::CutWizard>                      wizard_np,
                                      const Teuchos::RCP<const LINALG::MapExtractor>  & dbcmaps);


  /// build a discretization from the element patch, on which incompressibility shall be enforced
  void PrepareIncompOptDiscret(Teuchos::RCP<GEO::CutWizard> wizard_np);

  /// evaluate incompressibility on element level
  void EvaluateIncompOpt(Teuchos::RCP<GEO::CutWizard> wizard);

  /// solve optimization problem
  void SolveIncompOptProblem(Teuchos::RCP<Epetra_Vector> initialvel);
  //@}

  const Teuchos::RCP<DRT::Discretization> & bgdis_;   ///< xfem fluid discretization
  const Teuchos::RCP<DRT::Discretization> & embdis_;  ///< embedded fluid discretization
  Teuchos::RCP<DRT::Discretization> incompdis_;       ///< auxiliary discretization to enforce incompressibility after projection

  //! @name map of node-gid to dof-gid for standard nodes
  //@{
  //! map for current time
  std::map<int, std::vector<int> > nodeToDof_std_np_;
  //! map for previous time
  std::map<int, std::vector<int> > nodeToDof_std_n_;
  //@}

  //! @name map of node-gid to dof-gid for enriched nodes
  //@{
  //! map for current time
  std::map<int, std::vector<int> > nodeToDof_enriched_np_;
  //! map for previous time
  std::map<int, std::vector<int> > nodeToDof_enriched_n_;
  //@}

  //! map of every node id to the ids of it's parent elements
  std::map<int, std::set<int> > nodeToParentEle_;

  //! map of the parent elements (from the map above) to a vector of nodal dofsets
  std::map<int, std::vector< std::vector< int > > > parentEleToNds_;

  //! set of IDs of projected nodes
  std::set<int> projectedNodes_;

  //! current time step
  int step_;

  //! current proc id
  int myrank_;

  //! number of procs
  int numproc_;

  //! current xfluid dof-rowmap
  const Epetra_Map* currentbgdofmap_;

  //! old xfluid dof-rowmap
  const Epetra_Map* oldbgdofmap_;

  //! map of node Ids from embedded discretization to previous node positions
  std::map<int,LINALG::Matrix<3,1> > emb_nodepositions_n_;

  //! 3D seach tree for embedded discretization
  Teuchos::RCP<GEO::SearchTree> searchTree_;

  //! if the status of nodes remains the same don't do any interpolation from the embedded fluid
  bool samemaps_;

  //! the type of time integration approach
  enum INPAR::XFEM::XFluidFluidTimeInt timeintapproach_;

  //! incompressibility vector
  Teuchos::RCP<Epetra_Vector> C_;

  //! set of node ids for incompressibility problem
  std::set<int> incompnodeids_set_;
  //! set of element ids for incompressibility problem
  std::set<int> incompelementids_set_;

  //! parameter list
  Teuchos::ParameterList params_;

  //! flag for gmsh-output
  bool gmsh_debug_out_;

  //! min. radius needed for the search tree
  double minradius_;

  //! search radius factor
  double searchradius_fac_;
};

} // namespace XFEM

#endif /*TIMEINTEGRATION_H_*/


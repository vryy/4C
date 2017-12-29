/*----------------------------------------------------------------------------*/
/*!
\file adapter_coupling.H

\brief Coupling of two discretizations (surface- or volume-coupling)

\level 2

\maintainer Christoph Ager

*/
/*----------------------------------------------------------------------------*/

#ifndef ADAPTER_COUPLING_H
#define ADAPTER_COUPLING_H

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace LINALG
{
  class SparseMatrix;
} // namespace LINALG

namespace DRT
{
  class Discretization;
  class DiscretizationInterface;
} // namespace DRT

/*----------------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_FEVector.h>
#include <Epetra_IntVector.h>

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace ADAPTER
{
  /*! \class CouplingBase
   *  \brief Abstract interface to coupling managers
   *  This class is an interface to the coupling adapters.
   *
   *  The standard version ADAPTER::Coupling (see below) can couple two
   *  (matching) sets of nodes from two discretizations at an surface-interface
   *  (e.g. fluid-structure) or within a volume (e.g. fluid-ale).
   *
   *  For non-matching sets of nodes within a volume, the ADAPTER::MortarVolCoupl
   *  can be used.
   *
   *  For non-matching sets of nodes at an interface, the ADAPTER::CouplingMortar
   *  can be used.
   */
  class CouplingBase
  {
  public:
    /// empty constructor
    CouplingBase(){};

    /// virtual destructor
    virtual ~CouplingBase(){};

    /// @name Conversion between master and slave
    //@{

    /// transfer a dof vector from master to slave
    virtual Teuchos::RCP<Epetra_Vector> MasterToSlave(
        Teuchos::RCP<Epetra_Vector> mv ///< master vector (to be transferred)
        )const = 0;

    /// transfer a dof vector from slave to master
    virtual Teuchos::RCP<Epetra_Vector> SlaveToMaster(
        Teuchos::RCP<Epetra_Vector> sv ///< slave vector (to be transferred)
        ) const = 0;

    /// transfer a dof vector from master to slave
    virtual Teuchos::RCP<Epetra_MultiVector> MasterToSlave(
        Teuchos::RCP<Epetra_MultiVector> mv ///< master vector (to be transferred)
        ) const = 0;

    /// transfer a dof vector from slave to master
    virtual Teuchos::RCP<Epetra_MultiVector> SlaveToMaster(
        Teuchos::RCP<Epetra_MultiVector> sv ///< slave vector (to be transferred)
        ) const = 0;

    /// transfer a dof vector from master to slave
    virtual Teuchos::RCP<Epetra_Vector> MasterToSlave(
        Teuchos::RCP<const Epetra_Vector> mv ///< master vector (to be transferred)
        ) const = 0;

    /// transfer a dof vector from slave to master
    virtual Teuchos::RCP<Epetra_Vector> SlaveToMaster(
        Teuchos::RCP<const Epetra_Vector> sv ///< slave vector (to be transferred)
        ) const = 0;

    /// transfer a dof vector from master to slave
    virtual Teuchos::RCP<Epetra_MultiVector> MasterToSlave(
        Teuchos::RCP<const Epetra_MultiVector> mv ///< master vector (to be transferred)
        ) const = 0;

    /// transfer a dof vector from slave to master
    virtual Teuchos::RCP<Epetra_MultiVector> SlaveToMaster(
        Teuchos::RCP<const Epetra_MultiVector> sv ///< slave vector (to be transferred)
        ) const = 0;

    /// transfer a dof vector from master to slave
    virtual void MasterToSlave(
        Teuchos::RCP<const Epetra_MultiVector> mv, ///< master vector (to be transferred)
        Teuchos::RCP<Epetra_MultiVector> sv ///< slave vector (containing result)
        ) const = 0;

    /// transfer a dof vector from slave to master
    virtual void SlaveToMaster(
        Teuchos::RCP<const Epetra_MultiVector> sv, ///< slave vector (to be transferred)
        Teuchos::RCP<Epetra_MultiVector> mv ///< master vector (containing result)
        ) const = 0;

    //@}

    /** \name Coupled maps */
    //@{

    /// the interface dof map of the master side
    virtual Teuchos::RCP<const Epetra_Map> MasterDofMap() const = 0;

    /// the interface dof map of the slave side
    virtual Teuchos::RCP<const Epetra_Map> SlaveDofMap()  const = 0;

    //@}
  };

  /*! \class Coupling
   *  \brief Management of coupling between two (matching) sets of nodes from
   *  two discretizations
   *
   *  Interface coupled problems with matching meshes need to transfer
   *  results between their interfaces. That is values which belong to
   *  the dofs of one side have to be accessed by the other side as
   *  well. In a parallel setting, of course, these sides do not in
   *  general reside on the same processors.
   *
   *  This class handles the transfer of dof values at the interface
   *  (Epetra_Vectors build on the interface dof map of either side) to
   *  the other side. To distinguish both sides lets speak of master and
   *  slave, even though no side really dominates the other. On the
   *  contrary we provide perfect symmetry once the setup is done.
   *
   *  The idea is simple: We have a master dof map that describes the
   *  distribution of the master's interface dofs. And we have a slave
   *  dof map that describes the distribution of the slave's interface
   *  dofs. Both maps, however, live in their respective
   *  communicator. We cannot transfer data between them. Furthermore
   *  the dof gids are totally independent of each other.
   *
   *  So we build a mapping between the nodes from both meshes during
   *  the setup phase. This establishes the connection. From that
   *  connection we deduce (in a somewhat painful way) a permuted
   *  master dof map and a permuted slave dof map. These permuted
   *  maps are bound to have the same layout as the normal maps from the
   *  other side. So we can exchange dof values between fields by simply
   *  copying from a normal Epetra_Vector to the permuted
   *  Epetra_Vector from the other side without actually looking at the
   *  respective maps. Afterwards the communication happens within one
   *  field in the usual fashion. So the transfer functions
   *  MasterToSlave() and SlaveToMaster() are quite simple. The hard
   *  work happens (once) during setup.
   *
   *  \author u.kue
   *  \date 06/07
   */
  class Coupling: public CouplingBase
  {
  public:

    /// empty constructor
    Coupling();

    /** \name Setup */
    //@{

    /// setup coupling of nodes marked with condition
    /*!
      Setup the whole thing. Find matching nodes via octtree and build
      appropriate dof maps.

     */
    void SetupConditionCoupling(const DRT::Discretization& masterdis,
                                Teuchos::RCP<const Epetra_Map> mastercondmap,
                                const DRT::Discretization& slavedis,
                                Teuchos::RCP<const Epetra_Map> slavecondmap,
                                const std::string& condname,
                                const int numdof,
                                bool matchall=true,
                                const int nds_master=0,
                                const int nds_slave=0);

    /// setup coupling of nodes marked with condition1 not belonging
    /// also to condition2
    /*!
      Setup the whole thing. Find matching nodes via octtree, check
      affiliation to other condition and build appropriate dof maps.

     */
    void SetupConstrainedConditionCoupling(
        const DRT::Discretization& masterdis, ///< discretization of master side
        Teuchos::RCP<const Epetra_Map> mastercondmap, ///< map with condition DOFs of master side
        const DRT::Discretization& slavedis, ///< discretization of slave side
        Teuchos::RCP<const Epetra_Map> slavecondmap, ///< map with condition DOFs of slave side
        const std::string& condname1, ///< condition name
        const std::string& condname2, ///< condition name
        const int numdof, ///< number of DOFs to be coupled at each node
        bool matchall = true ///< Do all nodes need to match exactly?
        );

    /*! \brief Setup coupling of given nodes
     *
     *  Setup the whole thing. Find matching nodes via octtree and build
     *  appropriate dof maps.
     *
     *  \note All slave nodes are required to find a match. In contrast
     *  master nodes do not need to have a match. So it is legal to hand
     *  in more master nodes that slave nodes as long as there is one
     *  master node for each slave node.
     *
     *  We need some way to guess the tolerance for the octree. It must not be
     *  too small, otherwise we won't find matching nodes. Too large a tolerance
     *  will not hurt that much. It just means we will have to test more nodes.
     *  (But it could hurt, if the geometry is very small. It should be
     *   around the size of the smallest element)
     *
     *  \param masterdis   (i) master side mesh
     *  \param slavedis    (i) slave side mesh
     *  \param masternodes (i) gids of nodes on master side to be coupled
     *  \param slavenodes  (i) gids of nodes on slave side to be coupled
     *  \param numdof      (i) number of dofs per node to be coupled
     *  \param matchall    (i) flag indicating matching slave and master nodes
     *  \param tolerance   (i) tolerance for octree for node matching
     */
    void SetupCoupling(const DRT::Discretization& masterdis,
        const DRT::Discretization& slavedis,
        const std::vector<int>& masternodes,
        const std::vector<int>& slavenodes,
        const int numdof,
        const bool matchall=true,
        const double tolerance=1.e-3,
        const int nds_master=0,
        const int nds_slave=0);

    /*! \brief Setup coupling of given nodes
     *
     *  Setup the whole thing. Find matching nodes via octtree and build
     *  appropriate dof maps.
     *
     *  \note All slave nodes are required to find a match. In contrast
     *  master nodes do not need to have a match. So it is legal to hand
     *  in more master nodes that slave nodes as long as there is one
     *  master node for each slave node.
     *
     *  We need some way to guess the tolerance for the octree. It must not be
     *  too small, otherwise we won't find matching nodes. Too large a tolerance
     *  will not hurt that much. It just means we will have to test more nodes.
     *
     *  \param masterdis   (i) master side mesh
     *  \param slavedis    (i) slave side mesh
     *  \param masternodes (i) gids of nodes on master side to be coupled
     *  \param slavenodes  (i) gids of nodes on slave side to be coupled
     *  \param numdof      (i) number of dofs per node to be coupled
     *  \param matchall    (i) flag indicating matching slave and master nodes
     *  \param tolerance   (i) tolerance for octree for node matching
     */
    void SetupCoupling(const DRT::Discretization& masterdis,
        const DRT::Discretization& slavedis,
        const Epetra_Map& masternodes,
        const Epetra_Map& slavenodes,
        const int numdof,
        const bool matchall=true,
        const double tolerance=1.e-3,
        const int nds_master=0,
        const int nds_slave=0);

    /*! \brief Setup coupling of given nodes and node maps
     *
     *  Setup the whole thing. Build appropriate dof maps (for given, matching nodes on both sides).
     *
     *  \note All master and slave nodes are required to match.
     *  The permuted slave node map is given as input, thus no search tree to find matching nodes is used.
     *  Only the corresponding dof maps are build here (+ some checks).
     *  This setup method can be used, if one discretization was cloned from the others discretization,
     *  and thus matching node GIDs are assured (for Fluid-ALE or poro, for instance).
     *  In this case, this setup gives the same results as the standard version above, but
     *  without the need to perform the matching node search.
     *
     *  \param masterdis         (i) master side mesh
     *  \param slavedis          (i) slave side mesh
     *  \param masternodemap     (i) gids of nodes on master side to be coupled
     *  \param slavenodemap      (i) gids of nodes on slave side to be coupled
     *  \param permslavenodemap  (i) permuted node map of the slave side
     *  \param numdof            (i) number of dofs per node to be coupled
     */
    void SetupCoupling(const DRT::DiscretizationInterface& masterdis,
        const DRT::DiscretizationInterface& slavedis,
        const Epetra_Map& masternodemap,
        const Epetra_Map& slavenodemap,
        const Epetra_Map& permslavenodemap,
        const int numdof);

    /*! \brief Setup coupling of given master and slave discretization
     *
     *  Setup the whole thing for matching dofs on both sides
     *
     *  \note All master and slave dofs are required to match.
     *  This setup method can be used, if one discretization was cloned from the other discretization,
     *  and thus matching dof GIDs are assured. All dofs of all nodes are coupled this way.
     *  In this case, this setup gives the same results as the standard version above, but
     *  neither a matching node search nor building dof maps is necessary.
     *
     *  \param masterdis         (i) master side mesh
     *  \param slavedis          (i) slave side mesh
     */
    void SetupCoupling(const DRT::Discretization& masterdis,
        const DRT::Discretization& slavedis);

    //@}

    /** \name Conversion between master and slave */
    //@{
    /// There are different versions to satisfy all needs. The basic
    /// idea is the same for all of them.

    /// transfer a dof vector from master to slave
    Teuchos::RCP<Epetra_Vector> MasterToSlave(
        Teuchos::RCP<Epetra_Vector> mv ///< master vector (to be transferred)
        ) const
    { return MasterToSlave(mv.getConst()); }

    /// transfer a dof vector from slave to master
    Teuchos::RCP<Epetra_Vector> SlaveToMaster(
        Teuchos::RCP<Epetra_Vector> sv ///< slave vector (to be transferred)
        ) const
    { return SlaveToMaster(sv.getConst()); }

    /// transfer a dof vector from master to slave
    Teuchos::RCP<Epetra_FEVector> MasterToSlave(
        Teuchos::RCP<Epetra_FEVector> mv ///< master vector (to be transferred)
        ) const
    { return MasterToSlave(mv.getConst()); }

    /// transfer a dof vector from slave to master
    Teuchos::RCP<Epetra_FEVector> SlaveToMaster(
        Teuchos::RCP<Epetra_FEVector> sv ///< slave vector (to be transferred)
        ) const
    { return SlaveToMaster(sv.getConst()); }

    /// transfer a dof vector from master to slave
    Teuchos::RCP<Epetra_MultiVector> MasterToSlave(
        Teuchos::RCP<Epetra_MultiVector> mv ///< master vector (to be transferred)
        ) const
    { return MasterToSlave(mv.getConst()); }

    /// transfer a dof vector from slave to master
    Teuchos::RCP<Epetra_MultiVector> SlaveToMaster(
        Teuchos::RCP<Epetra_MultiVector> sv ///< slave vector (to be transferred)
        ) const
    { return SlaveToMaster(sv.getConst()); }

    /// transfer a dof vector from master to slave
    Teuchos::RCP<Epetra_Vector> MasterToSlave(
        Teuchos::RCP<const Epetra_Vector> mv ///< master vector (to be transferred)
        ) const;

    /// transfer a dof vector from slave to master
    Teuchos::RCP<Epetra_Vector> SlaveToMaster(
        Teuchos::RCP<const Epetra_Vector> sv ///< slave vector (to be transferred)
        ) const;

    /// transfer a dof vector from master to slave
    Teuchos::RCP<Epetra_FEVector> MasterToSlave(
        Teuchos::RCP<const Epetra_FEVector> mv ///< master vector (to be transferred)
        ) const;

    /// transfer a dof vector from slave to master
    Teuchos::RCP<Epetra_FEVector> SlaveToMaster(
        Teuchos::RCP<const Epetra_FEVector> sv ///< slave vector (to be transferred)
        ) const;

    /// transfer a dof vector from master to slave
    Teuchos::RCP<Epetra_MultiVector> MasterToSlave(
        Teuchos::RCP<const Epetra_MultiVector> mv  ///< master vector (to be transferred)
        ) const;

    /// transfer a dof vector from slave to master
    Teuchos::RCP<Epetra_MultiVector> SlaveToMaster(
        Teuchos::RCP<const Epetra_MultiVector> sv ///< slave vector (to be transferred)
        ) const;

    /// transfer a dof vector from master to slave
    void MasterToSlave(
        Teuchos::RCP<const Epetra_MultiVector> mv, ///< master vector (to be transferred)
        Teuchos::RCP<Epetra_MultiVector> sv ///< slave vector (containing result)
        ) const;

    /// transfer a dof vector from slave to master
    void SlaveToMaster(
        Teuchos::RCP<const Epetra_MultiVector> sv, ///< slave vector (to be transferred)
        Teuchos::RCP<Epetra_MultiVector> mv ///< master vector (containing result)
        ) const;

    /// transfer a dof vector from master to slave
    void MasterToSlave(
        const Epetra_IntVector& mv, ///< master vector (to be transferred)
        Epetra_IntVector& sv ///< slave vector (containing result)
        ) const;

    /// transfer a dof vector from slave to master
    void SlaveToMaster(
        const Epetra_IntVector& sv, ///< slave vector (to be transferred)
        Epetra_IntVector& mv ///< master vector (containing result)
        ) const;

    //@}

    //! @name Access to coupled maps
    //@{

    /// the interface dof map of the master side
    Teuchos::RCP<const Epetra_Map> MasterDofMap() const { return masterdofmap_; }

    /// the interface dof map of the slave side
    Teuchos::RCP<const Epetra_Map> SlaveDofMap()  const { return slavedofmap_; }

    /// the permuted interface dof map of the master side
    Teuchos::RCP<const Epetra_Map> PermMasterDofMap() const { return permmasterdofmap_; }

    /// the permuted interface dof map of the slave side
    Teuchos::RCP<const Epetra_Map> PermSlaveDofMap()  const { return permslavedofmap_; }

    //@}

    //! @name Matrix transform
    //@{

    /// fill rowmap with master -> slave pairs
    void FillMasterToSlaveMap(std::map<int,int>& rowmap) const;

    /// fill rowmap with slave -> master pairs
    void FillSlaveToMasterMap(std::map<int,int>& rowmap) const;

    /// fill partial mastermap with gid of partial slavemap
    Teuchos::RCP<Epetra_Map> SlaveToMasterMap(Teuchos::RCP<Epetra_Map> slave);

    /// fill partial slavemap with gid of partial mastermap
    Teuchos::RCP<Epetra_Map> MasterToSlaveMap(Teuchos::RCP<Epetra_Map> master);

    /// redistribute crsmatrix from master row map to permuted master row map
    Teuchos::RCP<LINALG::SparseMatrix> MasterToPermMaster(const LINALG::SparseMatrix& sm) const;

    /// redistribute crsmatrix from slave row map to permuted slave row map
    Teuchos::RCP<LINALG::SparseMatrix> SlaveToPermSlave(const LINALG::SparseMatrix& sm) const;

    //@}

    /// \name Lagrangian coupling helpers

    /// create coupling matrices for Lagrangian coupling conditions
    void SetupCouplingMatrices(const Epetra_Map& shiftedmastermap,
        const Epetra_Map& masterdomainmap,
        const Epetra_Map& slavedomainmap);

    Teuchos::RCP<Epetra_CrsMatrix> MasterToMasterMat() const { return matmm_; }
    Teuchos::RCP<Epetra_CrsMatrix> SlaveToMasterMat() const { return matsm_; }
    Teuchos::RCP<Epetra_CrsMatrix> MasterToMasterMatTrans() const { return matmm_trans_; }
    Teuchos::RCP<Epetra_CrsMatrix> SlaveToMasterMatTrans() const { return matsm_trans_; }

    //@}

  protected:
    virtual void BuildDofMaps(
        const DRT::DiscretizationInterface& masterdis,
        const DRT::DiscretizationInterface& slavedis,
        const Teuchos::RCP<const Epetra_Map>& masternodemap,
        const Teuchos::RCP<const Epetra_Map>& slavenodemap,
        const Teuchos::RCP<const Epetra_Map>& permmasternodemap,
        const Teuchos::RCP<const Epetra_Map>& permslavenodemap,
        const int& numdof,
        const int nds_master=0,
        const int nds_slave=0);

  private:

    /*! \brief Do the actual matching of the master <-> slave pairs
     *
     *  Here the octree is used. Afterwards all not paired master nodes
     *  are removed.
     *
     *  \param masterdis   (i) master side mesh
     *  \param slavedis    (i) slave side mesh
     *  \param masternodes (i/o) all master node gids. on output those that have a match
     *  \param permslavenodes (o) slave node gids permuted to match master node gids
     *  \param slavenodes (i) slave node gids
     *  \param matchall (i) bool indicating match of all slave and master nodes
     *  \param tolerance   (i) tolerance for octree for node matching
     */
    void MatchNodes(const DRT::Discretization& masterdis,
        const DRT::Discretization& slavedis,
        std::vector<int>& masternodes,
        std::vector<int>& permslavenodes,
        const std::vector<int>& slavenodes,
        const bool matchall,
        const double tolerance);

    /// build slave to master permutation and dof all maps
    void FinishCoupling(const DRT::DiscretizationInterface& masterdis,
        const DRT::DiscretizationInterface& slavedis,
        Teuchos::RCP<Epetra_Map> masternodemap,
        Teuchos::RCP<Epetra_Map> slavenodemap,
        Teuchos::RCP<Epetra_Map> permslavenodemap,
        const int numdof,
        const int nds_master=0,
        const int nds_slave=0);

    /// build dof maps from node maps
    /*!
      \note It is assumed that the first ndim dofs of each
      node are of interest.
     */
    void BuildDofMaps(const DRT::DiscretizationInterface& dis,
        Teuchos::RCP<const Epetra_Map> nodemap,
        Teuchos::RCP<const Epetra_Map> permnodemap,
        Teuchos::RCP<const Epetra_Map>& dofmap,
        Teuchos::RCP<const Epetra_Map>& permdofmap,
        Teuchos::RCP<Epetra_Export>& exporter,
        const int numdof,
        const int nds=0) const;

  protected:
    /// @name accessors to the private class members for derived classes
    /// @{

    /// access the interface DoF map of the master side
    Teuchos::RCP<const Epetra_Map>& MaDofMapPtr();
    const Epetra_Map& MaDofMap() const;

    /// access the permuted interface DoF map of the master side
    Teuchos::RCP<const Epetra_Map>& PermutedMaDofMapPtr();
    const Epetra_Map& PermutedMaDofMap() const;

    /// access the interface DoF map of the slave side
    Teuchos::RCP<const Epetra_Map>& SlDofMapPtr();
    const Epetra_Map& SlDofMap() const;

    /// access the permuted interface DoF map of the slave side
    Teuchos::RCP<const Epetra_Map>& PermutedSlDofMapPtr();
    const Epetra_Map& PermutedSlDofMap() const;

    /// access the permuted master dof map to master dof map exporter
    Teuchos::RCP<Epetra_Export>& MaExporterPtr();
    const Epetra_Export& MaExporter() const;

    /// access the permuted slave dof map to slave dof map exporter
    Teuchos::RCP<Epetra_Export>& SlExporterPtr();
    const Epetra_Export& SlExporter() const;

    /// @}


  private:
    /*! @name Fundamental dof maps
     *
     *  We keep the master and slave dof map as well as permuted versions
     *  that match the respective other side.
     */
    //@{

    //! the interface dof map of the master side
    Teuchos::RCP<const Epetra_Map> masterdofmap_;

    //! the permuted interface dof map of the master side
    Teuchos::RCP<const Epetra_Map> permmasterdofmap_;

    //! the interface dof map of the slave side
    Teuchos::RCP<const Epetra_Map> slavedofmap_;

    //! the permuted interface dof map of the slave side
    Teuchos::RCP<const Epetra_Map> permslavedofmap_;

    //@}

    //! @name Communication object
    //@{

    //! permuted master dof map to master dof map exporter
    Teuchos::RCP<Epetra_Export> masterexport_;

    //! permuted slave dof map to slave dof map exporter
    Teuchos::RCP<Epetra_Export> slaveexport_;

    //@}

    //! @name coupling matrices for Lagrangian multiplier coupling
    //@{

    Teuchos::RCP<Epetra_CrsMatrix> matmm_;
    Teuchos::RCP<Epetra_CrsMatrix> matsm_;
    Teuchos::RCP<Epetra_CrsMatrix> matmm_trans_;
    Teuchos::RCP<Epetra_CrsMatrix> matsm_trans_;

    //@}
  };


  /*! \class CouplingConverter
   *  \brief abstract converter base
   *
   *  The point is that many generic coupling algorithms that transfer data
   *  between master and slave might be used in both directions. These
   *  algorithms can utilize a Converter to enable use in both directions.
   *
   *  \author u.kue
   *  \date 05/08
   */
  class CouplingConverter
  {
  public:
    virtual ~CouplingConverter() {}

    virtual Teuchos::RCP<Epetra_Vector> SrcToDst(Teuchos::RCP<const Epetra_Vector> s) const = 0;

    virtual Teuchos::RCP<Epetra_Vector> DstToSrc(Teuchos::RCP<const Epetra_Vector> d) const = 0;

    virtual Teuchos::RCP<const Epetra_Map> SrcMap() const = 0;

    virtual Teuchos::RCP<const Epetra_Map> DstMap() const = 0;

    virtual Teuchos::RCP<const Epetra_Map> PermSrcMap() const = 0;

    virtual Teuchos::RCP<const Epetra_Map> PermDstMap() const = 0;

    virtual void FillSrcToDstMap(std::map<int,int>& rowmap) const = 0;

  };

  /// master to slave converter
  class CouplingMasterConverter : public CouplingConverter
  {
  public:
    explicit CouplingMasterConverter(const Coupling& coup) : coup_(coup) {}

    virtual Teuchos::RCP<Epetra_Vector> SrcToDst(Teuchos::RCP<const Epetra_Vector> s) const
    { return coup_.MasterToSlave(s); }

    virtual Teuchos::RCP<Epetra_Vector> DstToSrc(Teuchos::RCP<const Epetra_Vector> d) const
    { return coup_.SlaveToMaster(d); }

    virtual Teuchos::RCP<const Epetra_Map> SrcMap() const
    { return coup_.MasterDofMap(); }

    virtual Teuchos::RCP<const Epetra_Map> DstMap() const
    { return coup_.SlaveDofMap(); }

    virtual Teuchos::RCP<const Epetra_Map> PermSrcMap() const
    { return coup_.PermMasterDofMap(); }

    virtual Teuchos::RCP<const Epetra_Map> PermDstMap() const
    { return coup_.PermSlaveDofMap(); }

    virtual void FillSrcToDstMap(std::map<int,int>& rowmap) const
    { coup_.FillMasterToSlaveMap(rowmap); }

  private:
    const Coupling& coup_;
  };

  /// slave to master converter
  class CouplingSlaveConverter : public CouplingConverter
  {
  public:
    explicit CouplingSlaveConverter(const Coupling& coup) : coup_(coup) {}

    virtual Teuchos::RCP<Epetra_Vector> SrcToDst(Teuchos::RCP<const Epetra_Vector> s) const
    { return coup_.SlaveToMaster(s); }

    virtual Teuchos::RCP<Epetra_Vector> DstToSrc(Teuchos::RCP<const Epetra_Vector> d) const
    { return coup_.MasterToSlave(d); }

    virtual Teuchos::RCP<const Epetra_Map> SrcMap() const
    { return coup_.SlaveDofMap(); }

    virtual Teuchos::RCP<const Epetra_Map> DstMap() const
    { return coup_.MasterDofMap(); }

    virtual Teuchos::RCP<const Epetra_Map> PermSrcMap() const
    { return coup_.PermSlaveDofMap(); }

    virtual Teuchos::RCP<const Epetra_Map> PermDstMap() const
    { return coup_.PermMasterDofMap(); }

    virtual void FillSrcToDstMap(std::map<int,int>& rowmap) const
    { coup_.FillSlaveToMasterMap(rowmap); }

  private:
    const Coupling& coup_;
  };

}

#endif /* ADAPTER_COUPLING_H */

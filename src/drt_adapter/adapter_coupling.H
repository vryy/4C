/*----------------------------------------------------------------------*/
/*!
\file adapter_coupling.H

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef FSI_COUPLING_H
#define FSI_COUPLING_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

namespace LINALG
{
  class SparseMatrix;
}

namespace DRT
{
  class Discretization;
}

namespace ADAPTER
{
  /// Management of coupling between two (matching) sets of nodes from two discretizations
  /*!
    Interface coupled problems with matching meshes need to transfer
    results between their interfaces. That is values which belong to
    the dofs of one side have to be accessed by the other side as
    well. In a parallel setting, of course, these sides do not in
    general reside on the same processors.

    This class handles the transfer of dof values at the interface
    (Epetra_Vectors build on the interface dof map of either side) to
    the other side. To distinguish both sides lets speak of master and
    slave, even though no side really dominates the other. On the
    contrary we provide perfect symmetry once the setup is done.

    The idea is simple: We have a master dof map that describes the
    distribution of the master's interface dofs. And we have a slave
    dof map that describes the distribution of the slave's interface
    dofs. Both maps, however, live in there respective
    communicator. We cannot transfer data between them. Furthermore
    the dof gids are totally independent of each other.

    So we build a mapping between the nodes from both meshes during
    the setup phase. This establishes the connection. From that
    connection we deduce (in a somewhat painful way) a permutated
    master dof map and a permutated slave dof map. These permutated
    maps are bound to have the same layout as the normal maps from the
    other side. So we can exchange dof values between fields by simply
    copying from a normal Epetra_Vector to the permutated
    Epetra_Vector from the other side without actually looking at the
    respective maps. Afterwards the communication happens within one
    field in the usual fashion. So the transfer functions
    MasterToSlave() and SlaveToMaster() are quite simple. The hard
    work happens (once) during setup.

    \author u.kue
    \date 06/07
   */
  class Coupling
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
                                bool matchall=true);

    /// setup coupling of nodes marked with condition1 not belonging
    /// also to condition2
    /*!
      Setup the whole thing. Find matching nodes via octtree, check
      affiliation to other condition and build appropriate dof maps.

     */
    void SetupConstrainedConditionCoupling(const DRT::Discretization& masterdis,
                                           Teuchos::RCP<const Epetra_Map> mastercondmap,
                                           const DRT::Discretization& slavedis,
                                           Teuchos::RCP<const Epetra_Map> slavecondmap,
                                           const std::string& condname1,
                                           const std::string& condname2,
                                           const int numdof,
                                           bool matchall=true);

    /// setup coupling of given nodes
    /*!
      Setup the whole thing. Find matching nodes via octtree and build
      appropriate dof maps.

      \note All slave nodes are required to find a match. In contrast
      master nodes do not need to have a match. So it is legal to hand
      in more master nodes that slave nodes as long as there is one
      master node for each slave node.

      \param masterdis   (i) master side mesh
      \param slavedis    (i) slave side mesh
      \param masternodes (i) gids of nodes on master side to be coupled
      \param slavenodes  (i) gids of nodes on slave side to be coupled
      \param numdof      (i) number of dofs per node to be coupled
      \param matchall    (i) flag indicating matching slave and master nodes
     */
    void SetupCoupling(const DRT::Discretization& masterdis,
                       const DRT::Discretization& slavedis,
                       const std::vector<int>& masternodes,
                       const std::vector<int>& slavenodes,
                       const int numdof,
                       bool matchall=true);

    /// setup coupling of given nodes
    /*!
      Setup the whole thing. Find matching nodes via octtree and build
      appropriate dof maps.

      \note All slave nodes are required to find a match. In contrast
      master nodes do not need to have a match. So it is legal to hand
      in more master nodes that slave nodes as long as there is one
      master node for each slave node.

      \param masterdis   (i) master side mesh
      \param slavedis    (i) slave side mesh
      \param masternodes (i) gids of nodes on master side to be coupled
      \param slavenodes  (i) gids of nodes on slave side to be coupled
      \param numdof      (i) number of dofs per node to be coupled
      \param matchall    (i) flag indicating matching slave and master nodes
     */
    void SetupCoupling(const DRT::Discretization& masterdis,
                       const DRT::Discretization& slavedis,
                       const Epetra_Map& masternodes,
                       const Epetra_Map& slavenodes,
                       const int numdof,
                       bool matchall=true);

    //@}

    /** \name Conversion between master and slave */
    //@{
    /// There are different versions to satisfy all needs. The basic
    /// idea is the same for all of them.

    /// transfer a dof vector from master to slave
    Teuchos::RCP<Epetra_Vector> MasterToSlave(Teuchos::RCP<Epetra_Vector> mv) const
      { return MasterToSlave(Teuchos::rcp_static_cast<const Epetra_Vector>(mv)); }

    /// transfer a dof vector from slave to master
    Teuchos::RCP<Epetra_Vector> SlaveToMaster(Teuchos::RCP<Epetra_Vector> sv) const
      { return SlaveToMaster(Teuchos::rcp_static_cast<const Epetra_Vector>(sv)); }

    /// transfer a dof vector from master to slave
    Teuchos::RCP<Epetra_MultiVector> MasterToSlave(Teuchos::RCP<Epetra_MultiVector> mv) const
      { return MasterToSlave(Teuchos::rcp_static_cast<const Epetra_MultiVector>(mv)); }

    /// transfer a dof vector from slave to master
    Teuchos::RCP<Epetra_MultiVector> SlaveToMaster(Teuchos::RCP<Epetra_MultiVector> sv) const
      { return SlaveToMaster(Teuchos::rcp_static_cast<const Epetra_MultiVector>(sv)); }

    /// transfer a dof vector from master to slave
    Teuchos::RCP<Epetra_Vector> MasterToSlave(Teuchos::RCP<const Epetra_Vector> mv) const;

    /// transfer a dof vector from slave to master
    Teuchos::RCP<Epetra_Vector> SlaveToMaster(Teuchos::RCP<const Epetra_Vector> sv) const;

    /// transfer a dof vector from master to slave
    Teuchos::RCP<Epetra_MultiVector> MasterToSlave(Teuchos::RCP<const Epetra_MultiVector> mv) const;

    /// transfer a dof vector from slave to master
    Teuchos::RCP<Epetra_MultiVector> SlaveToMaster(Teuchos::RCP<const Epetra_MultiVector> sv) const;

    /// transfer a dof vector from master to slave
    void MasterToSlave(Teuchos::RCP<const Epetra_MultiVector> mv, Teuchos::RCP<Epetra_MultiVector> sv) const;

    /// transfer a dof vector from slave to master
    void SlaveToMaster(Teuchos::RCP<const Epetra_MultiVector> sv, Teuchos::RCP<Epetra_MultiVector> mv) const;

    //@}

    /** \name Coupled maps */
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

    /** \name Matrix transform */
    //@{

    /// fill rowmap with master -> slave pairs
    void FillMasterToSlaveMap(std::map<int,int>& rowmap) const;

    /// fill rowmap with slave -> master pairs
    void FillSlaveToMasterMap(std::map<int,int>& rowmap) const;

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

  private:

    /// do the actual matching of the master <-> slave pairs
    /*!
      Here the octree is used. Afterwards all not paired master nodes
      are removed.

      \param masterdis   (i) master side mesh
      \param slavedis    (i) slave side mesh
      \param masternodes (i/o) all master node gids. on output those that have a match
      \param permslavenodes (o) slave node gids permutated to match master node gids
      \param slavenodes (i) slave node gids
      \param matchall (i) bool indicating match of all slave and master nodes
     */
    void MatchNodes(const DRT::Discretization& masterdis,
                    const DRT::Discretization& slavedis,
                    std::vector<int>& masternodes,
                    std::vector<int>& permslavenodes,
                    const std::vector<int>& slavenodes,
                    bool matchall=true);

    /// build slave to master permutation and dof all maps
    void FinishCoupling(const DRT::Discretization& masterdis,
                        const DRT::Discretization& slavedis,
                        Teuchos::RCP<Epetra_Map> masternodemap,
                        Teuchos::RCP<Epetra_Map> slavenodemap,
                        Teuchos::RCP<Epetra_Map> permslavenodemap,
                        const int numdof);

    /// build dof maps from node maps
    /*!
      \note It is assumed that the first ndim dofs of each
      node are of interest.
     */
    void BuildDofMaps(const DRT::Discretization& dis,
                      Teuchos::RCP<const Epetra_Map> nodemap,
                      Teuchos::RCP<const Epetra_Map> permnodemap,
                      Teuchos::RCP<const Epetra_Map>& dofmap,
                      Teuchos::RCP<const Epetra_Map>& permdofmap,
                      Teuchos::RCP<Epetra_Export>& exporter,
                      const int numdof);

    /** \name Fundamental dof maps */
    //@{
    /// We keep the master and slave dof map as well as permutated versions
    /// that match the respective other side.

    Teuchos::RCP<const Epetra_Map> masterdofmap_;
    Teuchos::RCP<const Epetra_Map> permmasterdofmap_;
    Teuchos::RCP<const Epetra_Map> slavedofmap_;
    Teuchos::RCP<const Epetra_Map> permslavedofmap_;

    //@}

    /** \name Communication object */
    //@{

    /// permuted master dof map to master dof map exporter
    Teuchos::RCP<Epetra_Export> masterexport_;

    /// permuted slave dof map to slave dof map exporter
    Teuchos::RCP<Epetra_Export> slaveexport_;

    //@}

    /** \name coupling matrices for Lagrangian multiplier coupling */
    //@{

    Teuchos::RCP<Epetra_CrsMatrix> matmm_;
    Teuchos::RCP<Epetra_CrsMatrix> matsm_;
    Teuchos::RCP<Epetra_CrsMatrix> matmm_trans_;
    Teuchos::RCP<Epetra_CrsMatrix> matsm_trans_;

    //@}
  };


  /// abstract converter base
  /*!
    The point is that many generic coupling algorithms that transfer data
    between master and slave might be used in both directions. These
    algorithms can utilize a Converter to enable use in both directions.

    \author u.kue
    \date 05/08
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

#endif

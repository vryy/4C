
#ifdef CCADISCRET

#ifndef FSI_COUPLING_H
#define FSI_COUPLING_H

#include <string>
#include <vector>
#include <set>

#include <Teuchos_RCP.hpp>
#include <Epetra_Export.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

#include "../drt_lib/drt_discret.H"


namespace FSI
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

    /// setup coupling of nodes marked with condition
    /*!
      Setup the whole thing. Find matching nodes via octtree and build
      appropiate dof maps.

      \param masterdis  (i) master side mesh
      \param slavedis   (i) slave side mesh
      \param coupname   (i) name of condition on coupling nodes
     */
    void SetupConditionCoupling(const DRT::Discretization& masterdis,
                                const DRT::Discretization& slavedis,
                                std::string coupname);

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
     */
    void SetupCoupling(const DRT::Discretization& masterdis,
                       const DRT::Discretization& slavedis,
                       const std::vector<int>& masternodes,
                       const std::vector<int>& slavenodes);

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
     */
    void SetupCoupling(const DRT::Discretization& masterdis,
                       const DRT::Discretization& slavedis,
                       const Epetra_Map& masternodes,
                       const Epetra_Map& slavenodes);

    /// transfer a dof vector from master to slave
    Teuchos::RCP<Epetra_Vector> MasterToSlave(Teuchos::RCP<Epetra_Vector> mv) const;

    /// transfer a dof vector from slave to master
    Teuchos::RCP<Epetra_Vector> SlaveToMaster(Teuchos::RCP<Epetra_Vector> sv) const;

    /// the interface dof map of the master side
    Teuchos::RCP<Epetra_Map> MasterDofMap() const { return masterdofmap_; }

    /// the interface dof map of the slave side
    Teuchos::RCP<Epetra_Map> SlaveDofMap()  const { return slavedofmap_; }

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
     */
    void MatchNodes(const DRT::Discretization& masterdis,
                    const DRT::Discretization& slavedis,
                    std::vector<int>& masternodes,
                    std::vector<int>& permslavenodes,
                    const std::vector<int>& slavenodes);

    /// build slave to master permutation and dof all maps
    void FinishCoupling(const DRT::Discretization& masterdis,
                        const DRT::Discretization& slavedis,
                        Teuchos::RCP<Epetra_Map> masternodemap,
                        Teuchos::RCP<Epetra_Map> slavenodemap,
                        Teuchos::RCP<Epetra_Map> permslavenodemap);

    /// build dof maps from node maps
    /*!
      \note It is assumed that the first genprob.ndim dofs of each
      node are of interest.
     */
    void BuildDofMaps(const DRT::Discretization& dis,
                      RCP<Epetra_Map> nodemap,
                      RCP<Epetra_Map> permnodemap,
                      RCP<Epetra_Map>& dofmap,
                      RCP<Epetra_Map>& permdofmap,
                      RCP<Epetra_Export>& exporter);

    /// find all local nodes from discretization marked with condition
    void FindCondNodes(const DRT::Discretization& dis,
                       std::string condname,
                       std::set<int>& nodes);

    Teuchos::RCP<Epetra_Map> masterdofmap_;
    Teuchos::RCP<Epetra_Map> permmasterdofmap_;
    Teuchos::RCP<Epetra_Map> slavedofmap_;
    Teuchos::RCP<Epetra_Map> permslavedofmap_;

    Teuchos::RCP<Epetra_Export> masterexport_;
    Teuchos::RCP<Epetra_Export> slaveexport_;
  };

}

#endif
#endif

/*----------------------------------------------------------------------------*/
/**
\file xfield_field_coupling.H

\brief Special coupling routines to handle the possible changing number of DoFs
       from node to node during XFEM simulations

\maintainer Michael Hiermeier

\date Sep 28, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_XFEM_XFIELD_FIELD_COUPLING_H_
#define SRC_DRT_XFEM_XFIELD_FIELD_COUPLING_H_


#include "../drt_adapter/adapter_coupling.H"
#include "../drt_lib/drt_dserror.H"
// due to enum use ...
#include "xfem_multi_field_mapextractor.H"

namespace XFEM {
namespace XFieldField {
  class Coupling : public ADAPTER::Coupling
  {
  public:
    /** \brief enumerator specifying which interface discretization (slave
     *  or master) has the minimum number of DoF's at each node
     *
     *  min_dof_slave  : The slave discretization has the minimum number
     *                   of DoF's at each interface node.
     *
     *  min_dof_master : The master discretization has the minimum number
     *                   of DoF's at each interface node.
     *
     *  min_dof_unknown: The discretization with the minimum number of DoF's
     *                   per node is unknown or cannot be identified, since it
     *                   changes from node to node. This case needs extra
     *                   communication effort and is currently unsupported.
     *
     * \author hiermeier
     * \date 09/16 */
    enum MinDofDiscretization
    {
      min_dof_slave,
      min_dof_master,
      min_dof_unknown
    };

  public:
    /// constructor
    Coupling();

    /// destructor
    virtual ~Coupling() {};

    /// initialize class member variables
    void Init(const enum MinDofDiscretization& min_dof_dis);

    /** \name Conversion between master and slave
     *
     *  In contradiction to the base class versions, nodal distributed as well
     *  as dof distributed vectors can be transfered. */
    //@{
    /// There are different versions to satisfy all needs. The basic
    /// idea is the same for all of them.

    /** \brief transfer a nodel/dof vector from master to slave
     *
     *  \param mv       (in) : master vector (to be transferred)
     *  \param map_type (in) : map type of the master vector */
    inline Teuchos::RCP<Epetra_Vector> MasterToSlave(
        const Teuchos::RCP<Epetra_Vector>& mv,
        const enum XFEM::MapType & map_type) const
    { return MasterToSlave(mv.getConst(),map_type); }

    /** \brief transfer a nodal/dof vector from slave to master
     *
     *  \param sv       (in) : slave vector (to be transferred)
     *  \param map_type (in) : map type of the slave vector */
    inline Teuchos::RCP<Epetra_Vector> SlaveToMaster(
        Teuchos::RCP<Epetra_Vector> sv,
        const enum XFEM::MapType & map_type) const
    { return SlaveToMaster(sv.getConst(),map_type); }

    /** \brief transfer a nodal/dof multi vector from master to slave
     *
     *  \param mv       (in) : master multi vector (to be transferred)
     *  \param map_type (in) : map type of the master vector */
    inline Teuchos::RCP<Epetra_MultiVector> MasterToSlave(
        Teuchos::RCP<Epetra_MultiVector> mv,
        const enum XFEM::MapType & map_type) const
    { return MasterToSlave(mv.getConst(),map_type); }

     /** \brief transfer a nodal/dof multi vector from slave to master
      *
      *  \param sv       (in) : slave multi vector (to be transferred)
      *  \param map_type (in) : map type of the slave vector */
    inline Teuchos::RCP<Epetra_MultiVector> SlaveToMaster(
        Teuchos::RCP<Epetra_MultiVector> sv,
        const enum XFEM::MapType & map_type) const
    { return SlaveToMaster(sv.getConst(),map_type); }

    /** \brief transfer a nodel/dof vector from master to slave
     *
     *  \param mv       (in) : master vector (to be transferred)
     *  \param map_type (in) : map type of the master vector */
    Teuchos::RCP<Epetra_Vector> MasterToSlave(
        const Teuchos::RCP<const Epetra_Vector>& mv,
        const enum XFEM::MapType & map_type) const;

    /** \brief transfer a nodal/dof vector from slave to master
     *
     *  \param sv       (in) : slave vector (to be transferred)
     *  \param map_type (in) : map type of the slave vector */
    Teuchos::RCP<Epetra_Vector> SlaveToMaster(
        const Teuchos::RCP<const Epetra_Vector>& sv,
        const enum XFEM::MapType & map_type) const;

    /** \brief transfer a nodel/dof vector from master to slave
     *
     *  \param mv       (in) : master vector (to be transferred)
     *  \param map_type (in) : map type of the master vector */
    Teuchos::RCP<Epetra_MultiVector> MasterToSlave(
        const Teuchos::RCP<const Epetra_MultiVector>& mv,
        const enum XFEM::MapType & map_type) const;

    /** \brief transfer a nodal/dof multi vector from slave to master
     *
     *  \param sv       (in) : slave multi vector (to be transferred)
     *  \param map_type (in) : map type of the slave vector */
    Teuchos::RCP<Epetra_MultiVector> SlaveToMaster(
        const Teuchos::RCP<const Epetra_MultiVector>& sv,
        const enum XFEM::MapType & map_type) const;

    /** \brief transfer a nodel/dof multi vector from master to slave
     *
     *  \param mv       (in) : master multi vector (to be transferred/source)
     *  \param map_type (in) : map type of the master vector
     *  \param sv       (out): slave multi vector (target)*/
    void MasterToSlave(
        const Teuchos::RCP<const Epetra_MultiVector>& mv,
        const enum XFEM::MapType & map_type,
        Teuchos::RCP<Epetra_MultiVector> sv) const;

    /** \brief transfer a nodal/dof multi vector from slave to master
     *
     *  \param sv       (in) : slave multi vector (to be transferred)
     *  \param map_type (in) : map type of the slave vector
     *  \param mv       (out): master multi vector (target)*/
    void SlaveToMaster(
        const Teuchos::RCP<const Epetra_MultiVector>& sv,
        const enum XFEM::MapType & map_type,
        Teuchos::RCP<Epetra_MultiVector> mv) const;

    //@}

  protected:
    /// check the isinit_ flag
    inline void CheckInit() const
    {
      if (not isinit_)
        dserror("Call Init first!");
    }

    /** \brief build dof maps from node maps
     *
     *  \note It is assumed that the first numdof DoF's of each
     *  node are of interest. If numdof is equal -1 the maximum
     *  number of shared DoF's per node will be coupled, i.e.
     *  numdof is equivalent to the dof number of MinDofDiscretization
     *  at each node.
     *  Otherwise the base class variant is called.
     *
     *  \author  hiermeier
     *  \date 09/16  */
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
    /** \brief store the nodal maps and create the corresponding exporter objects
     *
     *  This becomes necessary, since we want to transfer dof's as well as nodal
     *  distributed vectors
     *
     *  \author hiermeier
     *  \date 10/16 */
    void SaveNodeMaps(
        const Teuchos::RCP<const Epetra_Map>& masternodemap,
        const Teuchos::RCP<const Epetra_Map>& slavenodemap,
        const Teuchos::RCP<const Epetra_Map>& permmasternodemap,
        const Teuchos::RCP<const Epetra_Map>& permslavenodemap);

    /** \brief Identify the number of DoF's of the min-dof discretization as well
     *  as the actual dof GID's
     *
     *  This method takes all current DoF's of each node in the min-dof-discretization
     *  and sets the necessary maps for the coupling.
     */
    void BuildMinDofMaps(
        const DRT::DiscretizationInterface& min_dis,
        const Epetra_Map& min_nodemap,
        const Epetra_Map& min_permnodemap,
        Teuchos::RCP<const Epetra_Map>& min_dofmap,
        Teuchos::RCP<const Epetra_Map>& min_permdofmap,
        Teuchos::RCP<Epetra_Export>& min_exporter,
        const Epetra_Map& max_nodemap,
        std::map<int,unsigned>& my_mindofpernode) const;

    void BuildMaxDofMaps(
        const DRT::DiscretizationInterface& max_dis,
        const Epetra_Map& max_nodemap,
        const Epetra_Map& max_permnodemap,
        Teuchos::RCP<const Epetra_Map>& max_dofmap,
        Teuchos::RCP<const Epetra_Map>& max_permdofmap,
        Teuchos::RCP<Epetra_Export>& max_exporter,
        const std::map<int,unsigned>& my_mindofpernode) const;

    inline const enum MinDofDiscretization & MinDofDis() const
    {
      return min_dof_dis_;
    }

  private:
    bool isinit_;

    enum MinDofDiscretization min_dof_dis_;

    Teuchos::RCP<const Epetra_Map> masternodemap_;

    Teuchos::RCP<const Epetra_Map> slavenodemap_;

    Teuchos::RCP<const Epetra_Map> permmasternodemap_;

    Teuchos::RCP<const Epetra_Map> permslavenodemap_;

    //! @name Nodal communication objects
    //@{

    //! permuted master dof map to master dof map exporter
    Teuchos::RCP<Epetra_Export> nodal_masterexport_;

    //! permuted slave dof map to slave dof map exporter
    Teuchos::RCP<Epetra_Export> nodal_slaveexport_;

    //@}
  };  // class Coupling
} // namespace XFieldField
} // namespace XFEM



#endif /* SRC_DRT_XFEM_XFIELD_FIELD_COUPLING_H_ */

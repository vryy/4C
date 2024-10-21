// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_CONDITION_PERIODIC_HPP
#define FOUR_C_FEM_CONDITION_PERIODIC_HPP

#include "4C_config.hpp"

#include "4C_fem_condition.hpp"

#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::DOFSets
{
  class PBCDofSet;
}  // namespace Core::DOFSets


namespace Core::Conditions
{
  /*!
  \brief  update dofrowmap and dofsets for periodic boundary conditions

          o effects distribution of master and slave nodes. Master and
            slave are owned by one proc afterwards.

          o effects ghosting of nodes. Master and slave nodes are ghosted
            in pairs.

          o passes list of coupled nodes to the dofset


  \author gammi

  */
  class PeriodicBoundaryConditions
  {
   public:
    /*!
    \brief Standard Constructor

    This method gets the periodic boundary conditions from the discretisation
    and creates some empty maps which will be filled by the node matching.


    \param  dis (i) Discretisation which contains the nodes of the
                    periodic boundary condition


    \return void

    */
    PeriodicBoundaryConditions(Teuchos::RCP<Core::FE::Discretization> dis, bool verbose = true);

    /*!
    \brief Destructor

    */
    virtual ~PeriodicBoundaryConditions() = default;

    //! @name Framework

    /*!
      \brief This method gathers all slave nodes on the master node's
      proc, eliminates their dofs and redistributes (nodes and
      master/slave couples) afterwards to get a better distribution
      of dofs among processors.

      After this function, the nodes and elements are redistributed among
      the processors according to the generated connectivity. During this
      procedure, the dof set is changed (coupled dofs are eliminated)
      and the dofs are renumbered.

    \return void

    */
    void update_dofs_for_periodic_boundary_conditions();
    //@}

    //! @name methods to provide constructed node to node coupling
    Teuchos::RCP<std::map<int, std::vector<int>>> return_all_coupled_col_nodes()
    {
      return (allcoupledcolnodes_);
    }
    //@}

    //! @name methods to provide constructed node to node coupling
    Teuchos::RCP<std::map<int, std::vector<int>>> return_all_coupled_row_nodes()
    {
      return (allcoupledrownodes_);
    }
    //@}


    bool has_pbc()
    {
      if (numpbcpairs_ > 0)
      {
        return (true);
      }
      else
      {
        return (false);
      }
    }


    //! @name method that returns a pointer to the vector of the conditions
    std::vector<Core::Conditions::Condition *> *return_surface_pb_cs() { return &mysurfpbcs_; }
    //@}


    //! @name handling of all pbc pairs

    /*!
    \brief generate maps of master->slave couples and send the slaves
    to the master's proc. Update the dofset.

     o generate master->slave connectivities
        * allcoupledrownodes_
        * allcoupledcolnodes_ (including ghosted master/slave nodes)

        This method contains a loop over all existing pairs of
        periodic boundary conditions (always two surfaces (3d) or lines
        (2d) belong to the same condition. They are usually referred to
        as corresponding "master" and "slave" conditons).

        For every pair a node coupling is generated by a call to the
        method create_node_coupling_for_single_pbc. This means, that we have
        the information which masternode id is coupled to which
        slavenode id (map midtosid).

        This node coupling is added to the connectivity map of all
        previous periodic boundary conditions (this allows pbcs in
        several spatial directions which are nescessary for example
        for 3d channel flows).

     o send slave nodes to master proc.

     o Generate a new dofset in which slaves do not have their own dofs
       anymore; they just point to the master's dofs

       After this function call, the nodes and elements are distributed
       among the processors according to the generated connectivity.
       During this procedure, the dof set is changed (coupled dofs are
       eliminated) and the dofs are renumbered.

       \return void

    */
    void put_all_slaves_to_masters_proc();

   protected:
    /*!
    \brief Couple nodes for specific pair of periodic boundary conditions.

    On input, this method needs two lists of nodes (master- and
    slavenodes). For every masternode, the method will look for the
    "closest" slavenode and insert this connectivity in the midtosid map.
    In this context, "closest" means, that one coordinate of the slavenode
    is replaced by the coordinate of the masternode plane for the search
    of the closest point. This is controlled by the parameter dofsforpbcplane,
    an int vector of length 3:


                       (1,1,0)    xy-plane

                       (0,1,1)    yz-plane

                       (1,0,1)    xz-plane

    The search algorithm will substitute the coordinate with 0 of the
    slavenodes by the corresponding coordinate of the masternode plane
    to calculate the closest node.

    Keep in mind that the masternodes and the slavenodes are distributed
    among several processors. The connectivity information will be only
    available on the proc which owns the specific masternode --- but the
    corresponding closest slavenode will be searched on all procs!

    The node matching is created in two steps. Step one generates a
    processor local octree that contains all masternodes available on
    this proc (NodeMatchingOctree nodematchingoctree).
    The second part searches for the closest slave node on
    all processors (nodematchingoctree.CreateGlobalNodeMatching).


    \param  midtosid        (o) map from master to slavenodes
    \param  masternodeids   (i) all master node ids
    \param  slavenodeids    (i) all slave node ids
    \param  dofsforpbcplane (i) periodic boundary is in this plane
    \param  rotangle        (i) angle (RAD) for rotation of slave plane
    \return void

    */
    void create_node_coupling_for_single_pbc(std::map<int, std::vector<int>> &midtosid,
        const std::vector<int> masternodeids, const std::vector<int> slavenodeids,
        const std::vector<int> dofsforpbcplane, const double rotangle, const double abstol);

    /*!
    \brief Add the connectivity from this condition (midtosid, on input)
    to the connectivity of all previously processed periodic boundary
    conditions (allcoupledrownodes_).

    We are interested in a map of all coupled nodes. For nodes, which have
    several periodic boundary conditions (for example edge nodes for
    3d channel flow, see figure below), only one node (the pure master node)
    is allowed to have its own degrees of freedom. The mixed master/slave
    nodes must be eliminated and all connectivity information must be
    transferred to the pure master.

          MM----------MS

           |xxxxxxxxxx|

           |xxxxxxxxxx|

           |xxxxxxxxxx|

          SM----------SS

    Since not all information is available on every proc, the completion
    of this coupling requires some additional communication of the coupling
    behaviour of multiple coupled nodes.


    \param  midtosid        (i) map from master to slavenodes of current condition
    \param  pbcid           (i) the id of the current periodic boundary condition

    \return void

    */
    void add_connectivity(std::map<int, std::vector<int>> &midtosid, const int pbcid);

    //@}


    //! @name Methods to finish the coupling after the list of coupled nodes is known.

    /*!
    \brief Redistribute the nodes and assign the dofs to the
    current distribution of nodes

    o Create a new rownodemap using the map allcoupledrownodes_

    o Build graph and temporary nodecolmap

    o Fix ghosting (if slaves are ghosted, masternodes have to
      be ghosted, too)
      This leads to the newnodecolmap

    o Create a new dofset specialisation for periodic boundary
      conditions

    o Redistribute nodes and call FillComplete including a call
      to the modified AssignDegreeesOfFreedom

    \return void

    */
    void redistribute_and_create_dof_coupling();

    //@}

    /*!
    \brief generate a better distribution of nodes for more efficiency
    in a parallel setting

     o adjust weights of slavenodes
       they need a small weight since they do not contribute any dofs
       to the linear system

     o compute connectivity
       iterate all elements on this proc including ghosted ones. Include
       additional connections between master and slave nodes

     o set weights of edges between master/slave pairs to a high value
       in order to keep both on the same proc when redistributing

     o gather all data to proc 1, do partitioning using METIS with
       both, weights for edges and nodes

     o redistribute nodes without assigning dofs

     o repair master/slave distribution, finally assign dofs

     \return void

    */
    void balance_load();

    //@}

    //!\brief the discretisation
    Teuchos::RCP<Core::FE::Discretization> discret_;

    //!\brief a flag controlling the verbosity, i.e. the amount of output
    // displayed on the screen
    bool verbose_;

    //!\brief the pbc-dofset created by this object
    Teuchos::RCP<Core::DOFSets::PBCDofSet> pbcdofset_;

    //!\brief number of pairs of periodic boundary conditions
    int numpbcpairs_;

    //!\brief vector of periodic surface boundary conditions
    std::vector<Core::Conditions::Condition *> mysurfpbcs_;

    //!\brief map connecting master to slave nodes owned by this proc
    //       master node -> list of his slave node(s)
    Teuchos::RCP<std::map<int, std::vector<int>>> allcoupledrownodes_;

    //!\brief map connecting master to slave nodes owned or ghosted by this proc
    //       master node -> list of his slave node(s)
    Teuchos::RCP<std::map<int, std::vector<int>>> allcoupledcolnodes_;

    //!\brief time measurement (total)
    Teuchos::RCP<Teuchos::Time> timepbctot_;
    //!\brief time measurement (create master slave matching for pairs)
    Teuchos::RCP<Teuchos::Time> timepbcmidtosid_;
    //!\brief time measurement (create octree)
    Teuchos::RCP<Teuchos::Time> timepbcmidoct_;
    //!\brief time measurement (search in octree)
    Teuchos::RCP<Teuchos::Time> timepbcmidmatch_;
    //!\brief time measurement (redistribute nodes)
    Teuchos::RCP<Teuchos::Time> timepbcreddis_;
    //!\brief time measurement (add connectivity to previous pbcs)
    Teuchos::RCP<Teuchos::Time> timepbcaddcon_;
    //!\brief time measurement (repair ghosting)
    Teuchos::RCP<Teuchos::Time> timepbcghost_;
    //!\brief time measurement (make colmap for ghosting)
    Teuchos::RCP<Teuchos::Time> timepbcmakeghostmap_;
    //!\brief time measurement (discret->redistribute)
    Teuchos::RCP<Teuchos::Time> timepbcrenumdofs_;


    //!\brief time measurement (total)
    Teuchos::RCP<Teuchos::TimeMonitor> tm0_ref_;
    //!\brief time measurement (create master slave matching for pairs)
    Teuchos::RCP<Teuchos::TimeMonitor> tm1_ref_;
    //!\brief time measurement (create octree)
    Teuchos::RCP<Teuchos::TimeMonitor> tm2_ref_;
    //!\brief time measurement (search in octree)
    Teuchos::RCP<Teuchos::TimeMonitor> tm3_ref_;
    //!\brief time measurement (add connectivity to previous pbcs)
    Teuchos::RCP<Teuchos::TimeMonitor> tm4_ref_;
    //!\brief time measurement (redistribute nodes)
    Teuchos::RCP<Teuchos::TimeMonitor> tm5_ref_;
    //!\brief time measurement (make row and colmap for ghosting)
    Teuchos::RCP<Teuchos::TimeMonitor> tm6_ref_;
    //!\brief time measurement (repair ghosting)
    Teuchos::RCP<Teuchos::TimeMonitor> tm7_ref_;
    //!\brief time measurement (discret->redistribute)
    Teuchos::RCP<Teuchos::TimeMonitor> tm8_ref_;
  };
}  // namespace Core::Conditions

FOUR_C_NAMESPACE_CLOSE

#endif

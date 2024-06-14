/*----------------------------------------------------------------------*/
/*! \file

\brief utils class for use of binning strategy

\level 2

*----------------------------------------------------------------------*/


#include "4C_binstrategy_utils.hpp"

#include "4C_comm_exporter.hpp"
#include "4C_io_pstream.hpp"
#include "4C_rebalance_binning_based.hpp"
#include "4C_rebalance_print.hpp"

FOUR_C_NAMESPACE_OPEN

namespace BINSTRATEGY
{
  namespace UTILS
  {
    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void ExtendDiscretizationGhosting(Teuchos::RCP<Core::FE::Discretization> discret,
        Teuchos::RCP<Epetra_Map> const& extendedelecolmap, bool assigndegreesoffreedom,
        bool initelements, bool doboundaryconditions)
    {
      // make sure that all procs are either filled or unfilled
      // oldmap in export_column_elements must be reset() on every proc or nowhere
      discret->CheckFilledGlobally();

      // adapt layout to extended ghosting in discret
      // first export the elements according to the processor local element column maps
      discret->export_column_elements(*extendedelecolmap);

      // get the node ids of the elements that are to be ghosted
      // and create a proper node column map for their export
      std::set<int> nodes;
      for (int lid = 0; lid < extendedelecolmap->NumMyElements(); ++lid)
      {
        Core::Elements::Element* ele = discret->gElement(extendedelecolmap->GID(lid));
        const int* nodeids = ele->NodeIds();
        for (int inode = 0; inode < ele->num_node(); ++inode) nodes.insert(nodeids[inode]);
      }

      std::vector<int> colnodes(nodes.begin(), nodes.end());
      Teuchos::RCP<Epetra_Map> nodecolmap = Teuchos::rcp(
          new Epetra_Map(-1, (int)colnodes.size(), colnodes.data(), 0, discret->Comm()));

      // now ghost the nodes
      discret->ExportColumnNodes(*nodecolmap);

      // fillcomplete discret with extended ghosting
      discret->fill_complete(assigndegreesoffreedom, initelements, doboundaryconditions);

#ifdef FOUR_C_ENABLE_ASSERTIONS
      // print distribution after standard ghosting
      if (discret->Comm().MyPID() == 0)
        std::cout << "parallel distribution with extended ghosting" << std::endl;
      Core::Rebalance::UTILS::print_parallel_distribution(*discret);
#endif

      return;
    }

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void CommunicateElements(Teuchos::RCP<Core::FE::Discretization>& discret,
        std::map<int, std::vector<Core::Elements::Element*>> const& toranktosendeles)
    {
      // build exporter
      Core::Communication::Exporter exporter(discret->Comm());
      int const numproc = discret->Comm().NumProc();

      // -----------------------------------------------------------------------
      // send
      // -----------------------------------------------------------------------
      // ---- pack data for sending -----
      std::map<int, std::vector<char>> sdata;
      std::vector<int> targetprocs(numproc, 0);
      std::map<int, std::vector<Core::Elements::Element*>>::const_iterator p;
      for (p = toranktosendeles.begin(); p != toranktosendeles.end(); ++p)
      {
        std::vector<Core::Elements::Element*>::const_iterator iter;
        for (iter = p->second.begin(); iter != p->second.end(); ++iter)
        {
          Core::Communication::PackBuffer data;
          (*iter)->pack(data);
          sdata[p->first].insert(sdata[p->first].end(), data().begin(), data().end());
        }
        targetprocs[p->first] = 1;
      }

      // ---- send ----
      const int length = sdata.size();
      std::vector<MPI_Request> request(length);
      int tag = 0;
      for (std::map<int, std::vector<char>>::const_iterator p = sdata.begin(); p != sdata.end();
           ++p)
      {
        exporter.i_send(discret->Comm().MyPID(), p->first, (p->second).data(),
            (int)(p->second).size(), 1234, request[tag]);
        ++tag;
      }
      if (tag != length) FOUR_C_THROW("Number of messages is mixed up");

      // -----------------------------------------------------------------------
      // receive
      // -----------------------------------------------------------------------
      // ---- prepare receiving procs -----
      std::vector<int> summedtargets(numproc, 0);
      discret->Comm().SumAll(targetprocs.data(), summedtargets.data(), numproc);

      // ---- receive ----
      for (int rec = 0; rec < summedtargets[discret->Comm().MyPID()]; ++rec)
      {
        std::vector<char> rdata;
        int length = 0;
        int tag = -1;
        int from = -1;
        exporter.ReceiveAny(from, tag, rdata, length);
        if (tag != 1234)
          FOUR_C_THROW("Received on proc %i data with wrong tag from proc %i",
              discret->Comm().MyPID(), from);

        // ---- unpack ----
        {
          // Put received nodes into discretization
          std::vector<char>::size_type index = 0;
          while (index < rdata.size())
          {
            std::vector<char> data;
            Core::Communication::ParObject::extract_from_pack(index, rdata, data);
            // this Teuchos::rcp holds the memory of the node
            Teuchos::RCP<Core::Communication::ParObject> object =
                Teuchos::rcp(Core::Communication::Factory(data), true);
            Teuchos::RCP<Core::Elements::Element> element =
                Teuchos::rcp_dynamic_cast<Core::Elements::Element>(object);
            if (element == Teuchos::null) FOUR_C_THROW("Received object is not a element");

            // safety check
            if (discret->HaveGlobalElement(element->Id()) != true)
              FOUR_C_THROW(
                  "%i is getting owner of element %i without having it ghosted before, "
                  "this is not intended.",
                  discret->Comm().MyPID(), element->Id());

            // delete already existing element (as it has wrong internal variables)
            discret->DeleteElement(element->Id());
            // add node (ownership already adapted on sending proc)
            discret->add_element(element);
          }
          if (index != rdata.size())
            FOUR_C_THROW(
                "Mismatch in size of data %d <-> %d", static_cast<int>(rdata.size()), index);
        }
      }

      // wait for all communications to finish
      for (int i = 0; i < length; ++i) exporter.Wait(request[i]);
      // safety, should be a no time operation if everything works fine before
      discret->Comm().Barrier();
    }

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void CommunicateDistributionOfTransferredElementsToBins(
        Teuchos::RCP<Core::FE::Discretization>& discret,
        std::map<int, std::vector<std::pair<int, std::vector<int>>>> const& toranktosendbinids,
        std::map<int, std::set<int>>& bintorowelemap)
    {
      // build exporter
      Core::Communication::Exporter exporter(discret->Comm());
      int const numproc = discret->Comm().NumProc();

      // -----------------------------------------------------------------------
      // send
      // -----------------------------------------------------------------------
      // ---- pack data for sending -----
      std::map<int, std::vector<char>> sdata;
      std::vector<int> targetprocs(numproc, 0);
      std::map<int, std::vector<std::pair<int, std::vector<int>>>>::const_iterator p;
      for (p = toranktosendbinids.begin(); p != toranktosendbinids.end(); ++p)
      {
        std::vector<std::pair<int, std::vector<int>>>::const_iterator iter;
        for (iter = p->second.begin(); iter != p->second.end(); ++iter)
        {
          Core::Communication::PackBuffer data;
          Core::Communication::ParObject::add_to_pack(data, *iter);
          sdata[p->first].insert(sdata[p->first].end(), data().begin(), data().end());
        }
        targetprocs[p->first] = 1;
      }

      // ---- send ----
      const int length = sdata.size();
      std::vector<MPI_Request> request(length);
      int tag = 0;
      for (std::map<int, std::vector<char>>::const_iterator p = sdata.begin(); p != sdata.end();
           ++p)
      {
        exporter.i_send(discret->Comm().MyPID(), p->first, (p->second).data(),
            (int)(p->second).size(), 1234, request[tag]);
        ++tag;
      }
      if (tag != length) FOUR_C_THROW("Number of messages is mixed up");

      // -----------------------------------------------------------------------
      // receive
      // -----------------------------------------------------------------------
      // ---- prepare receiving procs -----
      std::vector<int> summedtargets(numproc, 0);
      discret->Comm().SumAll(targetprocs.data(), summedtargets.data(), numproc);

      // ---- receive ----
      for (int rec = 0; rec < summedtargets[discret->Comm().MyPID()]; ++rec)
      {
        std::vector<char> rdata;
        int length = 0;
        int tag = -1;
        int from = -1;
        exporter.ReceiveAny(from, tag, rdata, length);
        if (tag != 1234)
          FOUR_C_THROW("Received on proc %i data with wrong tag from proc %i",
              discret->Comm().MyPID(), from);

        // ---- unpack ----
        {
          // Put received nodes into discretization
          std::vector<char>::size_type index = 0;
          while (index < rdata.size())
          {
            std::pair<int, std::vector<int>> pair;
            Core::Communication::ParObject::extract_from_pack(index, rdata, pair);
            std::vector<int>::const_iterator j;
            for (j = pair.second.begin(); j != pair.second.end(); ++j)
              bintorowelemap[*j].insert(pair.first);
          }
          if (index != rdata.size())
            FOUR_C_THROW(
                "Mismatch in size of data %d <-> %d", static_cast<int>(rdata.size()), index);
        }
      }

      // wait for all communications to finish
      for (int i = 0; i < length; ++i) exporter.Wait(request[i]);
      // safety, should be a no time operation if everything works fine before
      discret->Comm().Barrier();
    }

    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    void GetCurrentNodePos(Teuchos::RCP<const Core::FE::Discretization> const discret,
        Core::Nodes::Node const* node,
        std::function<Core::Nodes::Node const*(Core::Nodes::Node const* node)>
            correct_beam_center_node,
        Teuchos::RCP<const Epetra_Vector> const disnp, double* currpos)
    {
      // Todo make this nicer

      // the problem is that we might have nodes without position DoFs
      // (e.g. for beam elements with 'interior' nodes that are only used for
      // triad interpolation)
      // instead of the node position itself, we return the position of the
      // first node of the  element here (for the sake of binning)

      // standard case
      Core::Nodes::Node const* node_with_position_Dofs = correct_beam_center_node(node);

      if (disnp != Teuchos::null)
      {
        const int gid = discret->Dof(node_with_position_Dofs, 0);
        const int lid = disnp->Map().LID(gid);
        if (lid < 0)
          FOUR_C_THROW(
              "Your displacement is incomplete (need to be based on a column map"
              " as this function is also called from a loop over elements and "
              "each proc does (usually) not own all nodes of his row elements ");
        for (int dim = 0; dim < 3; ++dim)
        {
          currpos[dim] = node_with_position_Dofs->X()[dim] + (*disnp)[lid + dim];
        }
      }
      else
      {
        for (int dim = 0; dim < 3; ++dim) currpos[dim] = node_with_position_Dofs->X()[dim];
      }
    }
  }  // namespace UTILS
}  // namespace BINSTRATEGY

FOUR_C_NAMESPACE_CLOSE

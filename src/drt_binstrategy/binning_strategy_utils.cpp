/*----------------------------------------------------------------------*/
/*! \file

\brief utils class for use of binning strategy

\level 2

\maintainer Jonas Eichinger
*----------------------------------------------------------------------*/


#include "binning_strategy_utils.H"

#include "../drt_scatra_ele/scatra_ele.H"
#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_bele3/bele3.H"
#include "../drt_beam3/beam3_base.H"
#include "../drt_rigidsphere/rigidsphere.H"
#include "../drt_so3/so_base.H"

#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_exporter.H"

namespace BINSTRATEGY
{
  namespace UTILS
  {
    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void ExtendDiscretizationGhosting(Teuchos::RCP<DRT::Discretization> discret,
        Teuchos::RCP<Epetra_Map> const& extendedelecolmap, bool assigndegreesoffreedom,
        bool initelements, bool doboundaryconditions)
    {
      // make sure that all procs are either filled or unfilled
      // oldmap in ExportColumnElements must be Reset() on every proc or nowhere
      discret->CheckFilledGlobally();

      // adapt layout to extended ghosting in discret
      // first export the elements according to the processor local element column maps
      discret->ExportColumnElements(*extendedelecolmap);

      // get the node ids of the elements that are to be ghosted
      // and create a proper node column map for their export
      std::set<int> nodes;
      for (int lid = 0; lid < extendedelecolmap->NumMyElements(); ++lid)
      {
        DRT::Element* ele = discret->gElement(extendedelecolmap->GID(lid));
        const int* nodeids = ele->NodeIds();
        for (int inode = 0; inode < ele->NumNode(); ++inode) nodes.insert(nodeids[inode]);
      }

      std::vector<int> colnodes(nodes.begin(), nodes.end());
      Teuchos::RCP<Epetra_Map> nodecolmap =
          Teuchos::rcp(new Epetra_Map(-1, (int)colnodes.size(), &colnodes[0], 0, discret->Comm()));

      // now ghost the nodes
      discret->ExportColumnNodes(*nodecolmap);

      // fillcomplete discret with extended ghosting
      discret->FillComplete(assigndegreesoffreedom, initelements, doboundaryconditions);

#ifdef DEBUG
      // print distribution after standard ghosting
      if (discret->Comm().MyPID() == 0)
        std::cout << "parallel distribution with extended ghosting" << std::endl;
      DRT::UTILS::PrintParallelDistribution(*discret);
#endif

      return;
    }

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    BINSTRATEGY::UTILS::BinContentType ConvertElementToBinContentType(
        DRT::Element const* const eleptr)
    {
      // (Todo make this nicer and cheaper)

      if (dynamic_cast<DRT::ELEMENTS::Transport const*>(eleptr) != NULL)
      {
        return BINSTRATEGY::UTILS::Scatra;
      }
      else if (dynamic_cast<DRT::ELEMENTS::Fluid const*>(eleptr) != NULL)
      {
        return BINSTRATEGY::UTILS::Fluid;
      }
      else if (dynamic_cast<DRT::ELEMENTS::Bele3 const*>(eleptr) != NULL)
      {
        return BINSTRATEGY::UTILS::BELE3;
      }
      else if (dynamic_cast<DRT::ELEMENTS::Beam3Base const*>(eleptr) != NULL)
      {
        return BINSTRATEGY::UTILS::Beam;
      }
      else if (dynamic_cast<DRT::ELEMENTS::Rigidsphere const*>(eleptr) != NULL)
      {
        return BINSTRATEGY::UTILS::RigidSphere;
      }
      else if (dynamic_cast<DRT::ELEMENTS::So_base const*>(eleptr) != NULL)
      {
        return BINSTRATEGY::UTILS::Solid;
      }
      else
      {
        dserror(
            " Element you are about to assign to a bin could not be converted"
            " to a valid bin content type. ");
        exit(EXIT_FAILURE);
      }
    }

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void CommunicateElements(Teuchos::RCP<DRT::Discretization>& discret,
        std::map<int, std::vector<DRT::Element*>> const& toranktosendeles)
    {
      // build exporter
      DRT::Exporter exporter(discret->Comm());
      int const numproc = discret->Comm().NumProc();

      // -----------------------------------------------------------------------
      // send
      // -----------------------------------------------------------------------
      // ---- pack data for sending -----
      std::map<int, std::vector<char>> sdata;
      std::vector<int> targetprocs(numproc, 0);
      std::map<int, std::vector<DRT::Element*>>::const_iterator p;
      for (p = toranktosendeles.begin(); p != toranktosendeles.end(); ++p)
      {
        std::vector<DRT::Element*>::const_iterator iter;
        for (iter = p->second.begin(); iter != p->second.end(); ++iter)
        {
          DRT::PackBuffer data;
          (*iter)->Pack(data);
          data.StartPacking();
          (*iter)->Pack(data);
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
        exporter.ISend(discret->Comm().MyPID(), p->first, &((p->second)[0]),
            (int)(p->second).size(), 1234, request[tag]);
        ++tag;
      }
      if (tag != length) dserror("Number of messages is mixed up");

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
          dserror("Received on proc %i data with wrong tag from proc %i", discret->Comm().MyPID(),
              from);

        // ---- unpack ----
        {
          // Put received nodes into discretization
          std::vector<char>::size_type index = 0;
          while (index < rdata.size())
          {
            std::vector<char> data;
            DRT::ParObject::ExtractfromPack(index, rdata, data);
            // this Teuchos::rcp holds the memory of the node
            Teuchos::RCP<DRT::ParObject> object = Teuchos::rcp(DRT::UTILS::Factory(data), true);
            Teuchos::RCP<DRT::Element> element = Teuchos::rcp_dynamic_cast<DRT::Element>(object);
            if (element == Teuchos::null) dserror("Received object is not a element");

            // safety check
            if (discret->HaveGlobalElement(element->Id()) != true)
              dserror(
                  "%i is getting owner of element %i without having it ghosted before, "
                  "this is not intended.",
                  discret->Comm().MyPID(), element->Id());

            // delete already existing element (as it has wrong internal variables)
            discret->DeleteElement(element->Id());
            // add node (ownership already adapted on sending proc)
            discret->AddElement(element);
          }
          if (index != rdata.size())
            dserror("Mismatch in size of data %d <-> %d", static_cast<int>(rdata.size()), index);
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
        Teuchos::RCP<DRT::Discretization>& discret,
        std::map<int, std::vector<std::pair<int, std::vector<int>>>> const& toranktosendbinids,
        std::map<int, std::set<int>>& bintorowelemap)
    {
      // build exporter
      DRT::Exporter exporter(discret->Comm());
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
          DRT::PackBuffer data;
          DRT::ParObject::AddtoPack(data, *iter);
          data.StartPacking();
          DRT::ParObject::AddtoPack(data, *iter);
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
        exporter.ISend(discret->Comm().MyPID(), p->first, &((p->second)[0]),
            (int)(p->second).size(), 1234, request[tag]);
        ++tag;
      }
      if (tag != length) dserror("Number of messages is mixed up");

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
          dserror("Received on proc %i data with wrong tag from proc %i", discret->Comm().MyPID(),
              from);

        // ---- unpack ----
        {
          // Put received nodes into discretization
          std::vector<char>::size_type index = 0;
          while (index < rdata.size())
          {
            std::pair<int, std::vector<int>> pair;
            DRT::ParObject::ExtractfromPack(index, rdata, pair);
            std::vector<int>::const_iterator j;
            for (j = pair.second.begin(); j != pair.second.end(); ++j)
              bintorowelemap[*j].insert(pair.first);
          }
          if (index != rdata.size())
            dserror("Mismatch in size of data %d <-> %d", static_cast<int>(rdata.size()), index);
        }
      }

      // wait for all communications to finish
      for (int i = 0; i < length; ++i) exporter.Wait(request[i]);
      // safety, should be a no time operation if everything works fine before
      discret->Comm().Barrier();
    }

    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    void GetCurrentNodePos(Teuchos::RCP<const DRT::Discretization> const discret,
        DRT::Node const* node, Teuchos::RCP<const Epetra_Vector> const disnp, double* currpos)
    {
      // Todo make this nicer

      // the problem is that we might have nodes without position DoFs
      // (e.g. for beam elements with 'interior' nodes that are only used for
      // triad interpolation)
      // instead of the node position itself, we return the position of the
      // first node of the  element here (for the sake of binning)

      // standard case
      DRT::Node const* node_with_position_Dofs = node;

      const DRT::Element* element = node->Elements()[0];
      const DRT::ELEMENTS::Beam3Base* beamelement =
          dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(element);

      // fixme: should be do get position at xi with xi = 0.0?
      // if the node does not have position DoFs, we return the position of the first
      // node of the corresponding element
      if (beamelement != NULL and not beamelement->IsCenterlineNode(*node))
      {
        node_with_position_Dofs = beamelement->Nodes()[0];
      }

      if (disnp != Teuchos::null)
      {
        const int gid = discret->Dof(node_with_position_Dofs, 0);
        const int lid = disnp->Map().LID(gid);
        if (lid < 0)
          dserror(
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

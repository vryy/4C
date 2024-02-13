/*-----------------------------------------------------------*/
/*! \file
\brief calc utils for beam interaction framework


\level 3
*/
/*-----------------------------------------------------------*/

#include "baci_beaminteraction_calc_utils.hpp"

#include "baci_beam3_base.hpp"
#include "baci_beaminteraction_crosslinking_params.hpp"
#include "baci_beaminteraction_link.hpp"
#include "baci_beaminteraction_periodic_boundingbox.hpp"
#include "baci_beaminteraction_spherebeamlinking_params.hpp"
#include "baci_binstrategy.hpp"
#include "baci_discretization_geometry_intersection_math.hpp"
#include "baci_global_data.hpp"
#include "baci_io_pstream.hpp"
#include "baci_linalg_serialdensematrix.hpp"
#include "baci_linalg_serialdensevector.hpp"
#include "baci_linalg_utils_densematrix_communication.hpp"
#include "baci_linalg_utils_densematrix_multiply.hpp"
#include "baci_linalg_utils_sparse_algebra_create.hpp"
#include "baci_linalg_utils_sparse_algebra_manipulation.hpp"
#include "baci_rigidsphere.hpp"
#include "baci_so3_base.hpp"
#include "baci_solid_ele.hpp"

#include <Epetra_FEVector.h>

BACI_NAMESPACE_OPEN

namespace BEAMINTERACTION
{
  namespace UTILS
  {
    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    bool IsBeamElement(DRT::Element const& element)
    {
      return (dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(&element) != nullptr) ? true : false;
    }

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    bool IsRigidSphereElement(DRT::Element const& element)
    {
      return (dynamic_cast<const DRT::ELEMENTS::Rigidsphere*>(&element) != nullptr) ? true : false;
    }

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    bool IsBeamNode(DRT::Node const& node)
    {
      bool beameles = false;
      bool othereles = false;

      for (int i = 0; i < static_cast<int>(node.NumElement()); ++i)
      {
        if (IsBeamElement(*(node.Elements())[i]))
          beameles = true;
        else
          othereles = true;
      }

      if (beameles and othereles)
        dserror(
            "Beam elements and other (solid, rigid sphere) elements sharing the same node is "
            "currently not allowed in BACI!");

      return beameles;
    }

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    bool IsBeamCenterlineNode(DRT::Node const& node)
    {
      bool beamclnode = false;

      // TODO: actually we would have to check all elements of all processors!!! Gather?
      for (int i = 0; i < static_cast<int>(node.NumElement()); ++i)
      {
        const DRT::ELEMENTS::Beam3Base* beamele =
            dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(node.Elements()[i]);

        if (beamele != nullptr and beamele->IsCenterlineNode(node)) beamclnode = true;
      }
      return beamclnode;
    }

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    bool IsRigidSphereNode(DRT::Node const& node)
    {
      bool sphereele = false;
      bool othereles = false;

      // TODO: actually we would have to check all elements of all processors!!! Gather?
      for (int i = 0; i < node.NumElement(); ++i)
      {
        if (IsRigidSphereElement(*(node.Elements())[i]))
          sphereele = true;
        else
          othereles = true;
      }

      if (sphereele and othereles)
        dserror(
            "Rigid sphere elements and other (solid, beam) elements sharing "
            "the same node is currently not allowed in BACI!");

      return sphereele;
    }

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void PeriodicBoundaryConsistentDisVector(Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<const CORE::GEO::MESHFREE::BoundingBox> const& pbb,
        Teuchos::RCP<const DRT::Discretization> const& discret)
    {
      CORE::LINALG::Matrix<3, 1> d;
      CORE::LINALG::Matrix<3, 1> X;
      int doflid[3];

      for (int i = 0; i < discret->NumMyRowNodes(); ++i)
      {
        d.Clear();
        X.Clear();

        // get a pointer at i-th row node
        DRT::Node* node = discret->lRowNode(i);

        /* Hermite Interpolation: Check whether node is a beam node which is NOT
         * used for centerline interpolation if so, we simply skip it because
         * it does not have position DoFs */
        if (IsBeamNode(*node) and not IsBeamCenterlineNode(*node)) continue;

        // get GIDs of this node's degrees of freedom
        std::vector<int> dofnode = discret->Dof(node);

        for (int dim = 0; dim < 3; ++dim)
        {
          doflid[dim] = dis->Map().LID(dofnode[dim]);
          d(dim) = (*dis)[doflid[dim]];
          X(dim) = node->X()[dim];
        }
        // shift
        pbb->Shift3D(d, X);

        for (int dim = 0; dim < 3; ++dim)
        {
          (*dis)[doflid[dim]] = d(dim);
        }
      }
    }

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    std::vector<int> Permutation(int number)
    {
      // auxiliary variable
      int j = 0;

      GLOBAL::Problem::Instance()->Random()->SetRandRange(0.0, 1.0);

      // result vector initialized with ordered numbers from 0 to N-1
      std::vector<int> randorder(number, 0);
      for (int i = 0; i < (int)randorder.size(); i++) randorder[i] = i;

      for (int i = 0; i < number; ++i)
      {
        // generate random number between 0 and i
        j = (int)floor((i + 1.0) * GLOBAL::Problem::Instance()->Random()->Uni());

        /*exchange values at positions i and j (note: value at position i is i due to above
         *initialization and because so far only positions <=i have been changed*/
        randorder[i] = randorder[j];
        randorder[j] = i;
      }

      return randorder;
    }

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void GetCurrentElementDis(DRT::Discretization const& discret, DRT::Element const* ele,
        Teuchos::RCP<const Epetra_Vector> const& ia_discolnp, std::vector<double>& eledisp)
    {
      // clear
      eledisp.clear();

      std::vector<int> lm, lmowner, lmstride;

      ele->LocationVector(discret, lm, lmowner, lmstride);
      DRT::UTILS::ExtractMyValues(*ia_discolnp, eledisp, lm);
    }

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void GetCurrentUnshiftedElementDis(DRT::Discretization const& discret, DRT::Element const* ele,
        Teuchos::RCP<const Epetra_Vector> const& ia_discolnp,
        CORE::GEO::MESHFREE::BoundingBox const& pbb, std::vector<double>& eledisp)
    {
      GetCurrentElementDis(discret, ele, ia_discolnp, eledisp);

      // cast to beambase element
      DRT::ELEMENTS::Beam3Base const* beamele = dynamic_cast<DRT::ELEMENTS::Beam3Base const*>(ele);

      // so far, only beam elements can be cut by a periodic boundary
      if (beamele == nullptr) return;

      beamele->UnShiftNodePosition(eledisp, pbb);
    }

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    template <typename T>
    void SetFilamentBindingSpotPositions(
        Teuchos::RCP<DRT::Discretization> discret, Teuchos::RCP<T> params)
    {
      // todo: set somewhere else
      double const tol = CORE::GEO::TOL7;

      // temporarily extend ghosting, should be reverted somewhere later on as it is needed
      // in this fashion only here
      std::set<int> relevantfilaments;
      ExtendGhostingForFilamentBspotSetup(relevantfilaments, discret);

      // get pointers to all filament number conditions set
      std::vector<DRT::Condition*> filamentconditions(0);
      discret->GetCondition("BeamLineFilamentCondition", filamentconditions);

      // compute number of linker types
      std::vector<INPAR::BEAMINTERACTION::CrosslinkerType> linkertypes = params->LinkerTypes();

      // loop over all relevant (on myrank) filaments
      for (auto const& filiter : relevantfilaments)
      {
        // loop over all nodes of current filament, sort elements and calculate total filament
        // length
        std::vector<int> const* nodeids = filamentconditions[filiter]->GetNodes();
        std::vector<DRT::Element*> sortedfilamenteles(0);
        double filreflength = 0.0;
        ComputeFilamentLengthAndSortItsElements(sortedfilamenteles, nodeids, filreflength, discret);

        // loop over all linking types
        for (unsigned int linkertype_i = 0; linkertype_i < linkertypes.size(); ++linkertype_i)
        {
          // set start and end arc length parameter for filament binding spot
          double start = 0.0;
          double end = filreflength;
          double filamentbspotinterval =
              params->FilamentBspotIntervalGlobal(linkertypes[linkertype_i]);
          std::pair<double, double> filamentbspotrangelocal =
              params->FilamentBspotRangeLocal(linkertypes[linkertype_i]);
          std::pair<double, double> filamentbspotrangeglobal =
              params->FilamentBspotRangeGlobal(linkertypes[linkertype_i]);

          // in case certain range of filament was specified
          if (filamentbspotrangelocal.first > 0.0 and filamentbspotrangelocal.first < 1.0)
            start = filamentbspotrangelocal.first * filreflength;
          if (filamentbspotrangelocal.second > 0.0 and filamentbspotrangelocal.second < 1.0)
            end *= filamentbspotrangelocal.second;

          // in case binding spot interval is constant on every filament
          if (filamentbspotrangeglobal.first > 0.0)
          {
            // get arc parameter range for binding spot positions for current filament
            start = (filreflength < filamentbspotrangeglobal.first)
                        ? (filreflength + 1.0)
                        : (filamentbspotrangeglobal.first);
            end = (filreflength < filamentbspotrangeglobal.second or
                      filamentbspotrangeglobal.second < 0.0)
                      ? filreflength
                      : filamentbspotrangeglobal.second;
          }

          // set different filament binding spot interval for each filament dependent on its
          // reference length
          if (params->FilamentBspotIntervalLocal(linkertypes[linkertype_i]) > 0.0)
            filamentbspotinterval =
                (end - start) * params->FilamentBspotIntervalLocal(linkertypes[linkertype_i]);

          // get number of binding spots for current filament in bonds
          int numbspot = (start < filreflength and ((end - start) > 0.0))
                             ? (std::floor((end + tol - start) / filamentbspotinterval) + 1)
                             : 0;

          // in case current filament has no binding spots, we go to next filament
          if (numbspot == 0) continue;

          // print number of binding spots for current filament
          IO::cout(IO::debug) << "\n---------------------------------------------------------------"
                              << IO::endl;
          IO::cout(IO::debug) << numbspot << " binding spots of type "
                              << INPAR::BEAMINTERACTION::CrosslinkerType2String(
                                     linkertypes[linkertype_i])
                              << " on filament " << filiter << " (consists of "
                              << static_cast<int>(sortedfilamenteles.size()) << " elements)"
                              << " with: " << IO::endl;

          // set xis on element level
          SetBindingSpotsPositionsOnFilament(sortedfilamenteles, start, linkertypes[linkertype_i],
              numbspot, filamentbspotinterval, tol);

          IO::cout(IO::debug) << "---------------------------------------------------------------\n"
                              << IO::endl;
        }
      }
    }

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void ExtendGhostingForFilamentBspotSetup(
        std::set<int>& relevantfilaments, Teuchos::RCP<DRT::Discretization> discret)
    {
      std::set<int> setofnodegidswithrequiredelecloud;
      DetermineOffMyRankNodesWithRelevantEleCloudForFilamentBspotSetup(
          relevantfilaments, setofnodegidswithrequiredelecloud, discret);

      // collect element gids that need to be ghosted so that each proc ghosts all
      // elements of a filament containing at least one node of myrank
      // do communication to gather all elements to temporarily extend ghosting
      std::set<int> coleleset;
      for (int iproc = 0; iproc < discret->Comm().NumProc(); ++iproc)
      {
        // myrank == iproc: copy set to vector in order to broadcast data
        std::vector<int> requirednodes(0);
        if (iproc == discret->Comm().MyPID())
          requirednodes.insert(requirednodes.begin(), setofnodegidswithrequiredelecloud.begin(),
              setofnodegidswithrequiredelecloud.end());

        // proc i tells all procs how many nodegids it has
        int numnodes = requirednodes.size();
        discret->Comm().Broadcast(&numnodes, 1, iproc);

        // proc i actually sends nodegids
        requirednodes.resize(numnodes);
        discret->Comm().Broadcast(requirednodes.data(), numnodes, iproc);

        std::set<int> sdata;
        std::set<int> rdata;

        // each proc looks in node row map and adds element clouds of owned nodes to sdata
        for (int i = 0; i < numnodes; ++i)
        {
          // only if myrank is owner of requested node
          if (discret->NodeRowMap()->LID(requirednodes[i]) < 0) continue;

          // insert element cloud of current node
          DRT::Node* currnode = discret->gNode(requirednodes[i]);
          for (int j = 0; j < currnode->NumElement(); ++j)
            sdata.insert(currnode->Elements()[j]->Id());
        }

        // gather and store information on iproc
        CORE::LINALG::Gather<int>(sdata, rdata, 1, &iproc, discret->Comm());
        if (iproc == discret->Comm().MyPID()) coleleset = rdata;
      }

      // insert previous ghosting
      for (int lid = 0; lid < discret->NumMyColElements(); ++lid)
        coleleset.insert(discret->ElementColMap()->GID(lid));
      std::vector<int> colgids(coleleset.begin(), coleleset.end());

      // create new ele col map
      Teuchos::RCP<Epetra_Map> newelecolmap = Teuchos::rcp(
          new Epetra_Map(-1, static_cast<int>(colgids.size()), colgids.data(), 0, discret->Comm()));

      // temporarily extend ghosting
      BINSTRATEGY::UTILS::ExtendDiscretizationGhosting(discret, newelecolmap, true, false, true);
    }

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void DetermineOffMyRankNodesWithRelevantEleCloudForFilamentBspotSetup(
        std::set<int>& relevantfilaments, std::set<int>& setofrequirednodes,
        Teuchos::RCP<DRT::Discretization> discret)
    {
      // loop over all row nodes
      for (int rown = 0; rown < discret->NumMyRowNodes(); ++rown)
      {
        // get filament number of current node ( requirement: node belongs to only one filament)
        DRT::Condition* cond = discret->lRowNode(rown)->GetCondition("BeamLineFilamentCondition");

        // in case node (e.g. node of rigid sphere element) does not belong to a filament, go to
        // next node
        if (cond == nullptr) continue;

        // get filament number
        int const currfilnum = cond->GetInt("FilamentId");

        // if a filament has already been examined --> continue with next node
        if (relevantfilaments.find(currfilnum) != relevantfilaments.end()) continue;
        // filament is examined for the first time --> new entry for relevant filaments on myrank
        else
          relevantfilaments.insert(currfilnum);

        // loop over all nodes of current filament and store gids of nodes not owned by myrank
        std::vector<int> const* nodeids = cond->GetNodes();
        for (int i = 0; i < static_cast<int>(nodeids->size()); ++i)
          if (discret->NodeRowMap()->LID((*nodeids)[i]) < 0)
            setofrequirednodes.insert((*nodeids)[i]);
      }
    }

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void ComputeFilamentLengthAndSortItsElements(std::vector<DRT::Element*>& sortedfilamenteles,
        std::vector<int> const* nodeids, double& filreflength,
        Teuchos::RCP<DRT::Discretization> discret)
    {
      // loop over all nodes associated with current filament
      for (int nodei = 0; nodei < static_cast<int>(nodeids->size()); ++nodei)
      {
        // insert element cloud of current node
        DRT::Node* node = discret->gNode((*nodeids)[nodei]);
        for (int j = 0; j < node->NumElement(); ++j)
        {
          // only if element has not yet been added to the filaments elements
          if (std::find(sortedfilamenteles.begin(), sortedfilamenteles.end(),
                  node->Elements()[j]) != sortedfilamenteles.end())
            continue;

          DRT::ELEMENTS::Beam3Base* currbeamele =
              dynamic_cast<DRT::ELEMENTS::Beam3Base*>(node->Elements()[j]);

#ifdef BACI_DEBUG
          if (currbeamele == nullptr)
            dserror("DESIGN LINE BEAM FILAMENT CONDITIONS only applicable to beam elements.");
#endif

          // add element reference length of new element to filament reference length
          filreflength += currbeamele->RefLength();
          // add element
          sortedfilamenteles.push_back(node->Elements()[j]);
        }
      }
    }

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void SetBindingSpotsPositionsOnFilament(std::vector<DRT::Element*>& sortedfilamenteles,
        double start, INPAR::BEAMINTERACTION::CrosslinkerType linkertype, int numbspot,
        double filamentbspotinterval, double tol)
    {
      // default no occupied binding spots
      std::pair<double, double> elearcinterval = std::make_pair(0.0, 0.0);
      int bspotcounter = 0;
      double currbspotarcparam = start;
      // loop over all elements in sorted order
      for (unsigned int se_iter = 0; se_iter < sortedfilamenteles.size(); ++se_iter)
      {
        DRT::ELEMENTS::Beam3Base* beamele =
            dynamic_cast<DRT::ELEMENTS::Beam3Base*>(sortedfilamenteles[se_iter]);

        // init variables set in beam eles
        std::vector<double> bspotposxi;

        // get arc parameter range of current element
        elearcinterval =
            std::make_pair(elearcinterval.second, elearcinterval.second + beamele->RefLength());

        while ((currbspotarcparam > (elearcinterval.first - tol) and
                   currbspotarcparam <= (elearcinterval.second + tol)) and
               bspotcounter < numbspot)
        {
          // linear mapping from arc pos to xi (local element system)
          double xi = (2.0 * currbspotarcparam - elearcinterval.first - elearcinterval.second) /
                      beamele->RefLength();
          xi = (abs(round(xi) - xi) < tol) ? round(xi) : xi;
          bspotposxi.push_back(xi);

          // print to screen
          IO::cout(IO::debug) << bspotcounter + 1 << ". binding spot: "
                              << "xi = " << xi << " on element " << se_iter << " (gid "
                              << beamele->Id() << ")" << IO::endl;

          ++bspotcounter;
          currbspotarcparam += filamentbspotinterval;
        }

        // set element variables for binding spot position and status
        beamele->SetPositionsOfBindingSpotType(linkertype, bspotposxi);
      }
    }

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void GetPosAndTriadOfBindingSpot(DRT::Element* ele,
        Teuchos::RCP<CORE::GEO::MESHFREE::BoundingBox> const& pbb,
        INPAR::BEAMINTERACTION::CrosslinkerType linkertype, int locbspotnum,
        CORE::LINALG::Matrix<3, 1>& bspotpos, CORE::LINALG::Matrix<3, 3>& bspottriad,
        std::vector<double>& eledisp)
    {
      // cast to beambase element
      DRT::ELEMENTS::Beam3Base* beamele = dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ele);

#ifdef BACI_DEBUG
      if (beamele == nullptr) dserror("Dynamic cast to beam3base failed");
#endif

      // get current position at binding spot xi
      beamele->GetPosOfBindingSpot(bspotpos, eledisp, linkertype, locbspotnum, *pbb);

      // get current triad at binding spot xi
      beamele->GetTriadOfBindingSpot(bspottriad, eledisp, linkertype, locbspotnum);
    }

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void GetPosAndTriadOfBindingSpot(DRT::Discretization const& discret, DRT::Element* ele,
        Teuchos::RCP<Epetra_Vector> const& ia_discolnp,
        Teuchos::RCP<CORE::GEO::MESHFREE::BoundingBox> const& pbb,
        INPAR::BEAMINTERACTION::CrosslinkerType linkertype, int locbspotnum,
        CORE::LINALG::Matrix<3, 1>& bspotpos, CORE::LINALG::Matrix<3, 3>& bspottriad)
    {
      std::vector<double> eledisp;
      GetCurrentUnshiftedElementDis(discret, ele, ia_discolnp, *pbb, eledisp);

      GetPosAndTriadOfBindingSpot(ele, pbb, linkertype, locbspotnum, bspotpos, bspottriad, eledisp);
    }

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    bool IsDistanceOutOfRange(CORE::LINALG::Matrix<3, 1> const& pos1,
        CORE::LINALG::Matrix<3, 1> const& pos2, double const lowerbound, double const upperbound)
    {
      CORE::LINALG::Matrix<3, 1> dist_vec(true);
      dist_vec.Update(1.0, pos1, -1.0, pos2);

      const double distance = dist_vec.Norm2();

      return (distance < lowerbound or distance > upperbound) ? true : false;
    }

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    bool IsEnclosedAngleOutOfRange(CORE::LINALG::Matrix<3, 1> const& direction1,
        CORE::LINALG::Matrix<3, 1> const& direction2, double const lowerbound,
        double const upperbound)
    {
      // cosine of angle is scalar product of vectors divided by their norms
      // direction vectors should be unit vectors since they come from triads, but anyway ...
      double cos_angle = direction1.Dot(direction2) / direction1.Norm2() / direction2.Norm2();

      // to enable parallel directions that are numerically slightly > 1 ( < -1)
      // ( would lead to NaN in std::acos )
      if (abs(abs(cos_angle) - 1.0) < CORE::GEO::TOL12) cos_angle = std::round(cos_angle);

      double angle = std::acos(cos_angle);

      // acos returns angle \in [0,\pi] but we always want the acute angle here
      if (angle > 0.5 * M_PI) angle = M_PI - angle;

      return (angle < lowerbound or angle > upperbound) ? true : false;
    }

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    bool DoBeamElementsShareNodes(DRT::Element const* const beam, DRT::Element const* const nbbeam)
    {
      // check if two considered eles share nodes
      for (unsigned int i = 0; i < 2; ++i)
      {
        // node 0 and 1 are always first and last node, respectively
        if (beam->NodeIds()[i] == nbbeam->NodeIds()[0] ||
            beam->NodeIds()[i] == nbbeam->NodeIds()[1])
          return true;
      }

      return false;
    }

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void FEAssembleEleForceStiffIntoSystemVectorMatrix(const DRT::Discretization& discret,
        std::vector<int> const& elegid, std::vector<CORE::LINALG::SerialDenseVector> const& elevec,
        std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> const& elemat,
        Teuchos::RCP<Epetra_FEVector> fe_sysvec, Teuchos::RCP<CORE::LINALG::SparseMatrix> fe_sysmat)
    {
      // the entries of elevecX  belong to the Dofs of the element with GID elegidX
      // the rows    of elematXY belong to the Dofs of the element with GID elegidX
      // the columns of elematXY belong to the Dofs of the element with GID elegidY
      const DRT::Element* ele1 = discret.gElement(elegid[0]);
      const DRT::Element* ele2 = discret.gElement(elegid[1]);

      // get element location vector and ownerships
      std::vector<int> lmrow1;
      std::vector<int> lmrow2;
      std::vector<int> lmrowowner1;
      std::vector<int> lmrowowner2;
      std::vector<int> lmstride;

      ele1->LocationVector(discret, lmrow1, lmrowowner1, lmstride);
      ele2->LocationVector(discret, lmrow2, lmrowowner2, lmstride);

      // assemble both element vectors into global system vector
      if (fe_sysvec != Teuchos::null)
      {
        fe_sysvec->SumIntoGlobalValues(elevec[0].length(), lmrow1.data(), elevec[0].values());
        fe_sysvec->SumIntoGlobalValues(elevec[1].length(), lmrow2.data(), elevec[1].values());
      }

      // and finally also assemble stiffness contributions
      if (fe_sysmat != Teuchos::null)
      {
        fe_sysmat->FEAssemble(elemat[0][0], lmrow1, lmrow1);
        fe_sysmat->FEAssemble(elemat[0][1], lmrow1, lmrow2);
        fe_sysmat->FEAssemble(elemat[1][0], lmrow2, lmrow1);
        fe_sysmat->FEAssemble(elemat[1][1], lmrow2, lmrow2);
      }
    }
    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/

    /**
     *
     */
    unsigned int GetNumberOfElementCenterlineDof(const DRT::Element* ele)
    {
      std::vector<unsigned int> local_centerline_dof_indices;
      auto beam_element = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);
      if (beam_element != nullptr)
      {
        beam_element->CenterlineDofIndicesOfElement(local_centerline_dof_indices);
        return local_centerline_dof_indices.size();
      }
      else
        dserror("GetNumberOfElementCenterlineDof: The given element has to be a beam element.");
    }


    /**
     *
     */
    template <unsigned int n_centerline_dof>
    void GetElementCenterlineGIDIndices(DRT::Discretization const& discret, const DRT::Element* ele,
        CORE::LINALG::Matrix<n_centerline_dof, 1, int>& centerline_gid)
    {
      // First we need the local IDs of the centerline DOFs.
      std::vector<unsigned int> local_centerline_dof_indices;
      auto beam_element = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);
      if (beam_element != nullptr)
      {
        beam_element->CenterlineDofIndicesOfElement(local_centerline_dof_indices);
        if (local_centerline_dof_indices.size() != n_centerline_dof)
          dserror(
              "GetElementCenterlineGIDIndices: The size of the local_centerline_dof_indices (%d) "
              "does not match the number of beam dofs (%d).",
              local_centerline_dof_indices.size(), n_centerline_dof);
      }
      else
        dserror("GetElementCenterlineGIDIndices: The given element has to be a beam element.");

      // Get the GIDs of the whole beam element (including rotational DOFs).
      std::vector<int> lmrow;
      std::vector<int> lmrowowner;
      std::vector<int> lmstride;
      ele->LocationVector(discret, lmrow, lmrowowner, lmstride);

      // Get the GIDs of only the centerline DOFs.
      for (unsigned int i_dof = 0; i_dof < n_centerline_dof; i_dof++)
        centerline_gid(i_dof) = lmrow[local_centerline_dof_indices[i_dof]];
    }

    /**
     *
     */
    void GetElementCenterlineDOFIndices(DRT::Discretization const& discret, const DRT::Element* ele,
        std::vector<unsigned int>& ele_centerline_dof_indices, unsigned int& num_dof)
    {
      // Todo implement method in DRT::Element or find alternative way of doing this
      // find out the elements' number of Dofs (=dimension of element vector/matrices)
      std::vector<int> lmrow;
      std::vector<int> dummy1, dummy2;

      ele->LocationVector(discret, lmrow, dummy1, dummy2);
      num_dof = lmrow.size();

      const DRT::ELEMENTS::Beam3Base* beamele = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);

      if (beamele != nullptr)
      {
        beamele->CenterlineDofIndicesOfElement(ele_centerline_dof_indices);
      }
      else
      {
        ele_centerline_dof_indices.resize(num_dof);
        for (unsigned int i = 0; i < num_dof; ++i) ele_centerline_dof_indices[i] = i;
      }
    }

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void AssembleCenterlineDofForceStiffIntoElementForceStiff(DRT::Discretization const& discret,
        std::vector<int> const& elegid,
        std::vector<CORE::LINALG::SerialDenseVector> const& eleforce_centerlineDOFs,
        std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> const& elestiff_centerlineDOFs,
        std::vector<CORE::LINALG::SerialDenseVector>* eleforce,
        std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>>* elestiff)
    {
      std::vector<unsigned int> numdof_ele(2);
      std::vector<std::vector<unsigned int>> ele_centerlinedofindices(2);

      for (unsigned int iele = 0; iele < 2; ++iele)
      {
        DRT::Element* ele = discret.gElement(elegid[iele]);
        GetElementCenterlineDOFIndices(
            discret, ele, ele_centerlinedofindices[iele], numdof_ele[iele]);
      }


      // assemble centerline DOF values correctly into element DOFvec vectors/matrices
      if (eleforce != nullptr)
      {
        for (unsigned int iele = 0; iele < 2; ++iele)
        {
          // resize and clear variable
          ((*eleforce)[iele]).size(numdof_ele[iele]);

          // safety check: dimensions
          if ((unsigned int)eleforce_centerlineDOFs[iele].numRows() !=
              ele_centerlinedofindices[iele].size())
            dserror(
                "size mismatch! need to assemble %d values of centerline-Dof based "
                "force vector into element vector but only got %d element-local Dof indices",
                eleforce_centerlineDOFs[iele].numRows(), ele_centerlinedofindices[iele].size());

          // Todo maybe use a more general 'SerialDenseAssemble' method here
          for (unsigned int idof = 0; idof < ele_centerlinedofindices[iele].size(); ++idof)
            ((*eleforce)[iele])(ele_centerlinedofindices[iele][idof]) =
                eleforce_centerlineDOFs[iele](idof);
        }
      }

      if (elestiff != nullptr)
      {
        for (unsigned int iele = 0; iele < 2; ++iele)
        {
          for (unsigned int jele = 0; jele < 2; ++jele)
          {
            // resize and clear variable
            ((*elestiff)[iele][jele]).shape(numdof_ele[iele], numdof_ele[jele]);

            // safety check: dimensions
            if ((unsigned int)elestiff_centerlineDOFs[iele][jele].numRows() !=
                ele_centerlinedofindices[iele].size())
              dserror(
                  "size mismatch! need to assemble %d row values of centerline-Dof based "
                  "stiffness matrix into element matrix but only got %d element-local Dof indices",
                  elestiff_centerlineDOFs[iele][jele].numRows(),
                  ele_centerlinedofindices[iele].size());

            if ((unsigned int)elestiff_centerlineDOFs[iele][jele].numCols() !=
                ele_centerlinedofindices[jele].size())
              dserror(
                  "size mismatch! need to assemble %d column values of centerline-Dof based "
                  "stiffness matrix into element matrix but only got %d element-local Dof indices",
                  elestiff_centerlineDOFs[iele][jele].numCols(),
                  ele_centerlinedofindices[jele].size());

            for (unsigned int idof = 0; idof < ele_centerlinedofindices[iele].size(); ++idof)
              for (unsigned int jdof = 0; jdof < ele_centerlinedofindices[jele].size(); ++jdof)
                ((*elestiff)[iele][jele])(
                    ele_centerlinedofindices[iele][idof], ele_centerlinedofindices[jele][jdof]) =
                    elestiff_centerlineDOFs[iele][jele](idof, jdof);
          }
        }
      }
    }

    /**
     *
     */
    void AssembleCenterlineDofColMatrixIntoElementColMatrix(DRT::Discretization const& discret,
        const DRT::Element* element,
        CORE::LINALG::SerialDenseMatrix const& row_matrix_centerlineDOFs,
        CORE::LINALG::SerialDenseMatrix& row_matrix_elementDOFs)
    {
      // Get the centerline DOFs of the element.
      unsigned int numdof_ele;
      std::vector<unsigned int> ele_centerlinedofindices;
      GetElementCenterlineDOFIndices(discret, element, ele_centerlinedofindices, numdof_ele);

      // Fill in the centerline matrix into the full element matrix.
      // Resize and clear output matrix variable.
      row_matrix_elementDOFs.shape(row_matrix_centerlineDOFs.numRows(), numdof_ele);

      // Safety check: dimensions.
      if ((unsigned int)row_matrix_centerlineDOFs.numCols() != ele_centerlinedofindices.size())
        dserror(
            "Size mismatch! Need to assemble %d col values of centerline-Dof based "
            "stiffness matrix into element matrix but only got %d element-local Dof indices",
            row_matrix_centerlineDOFs.numCols(), ele_centerlinedofindices.size());

      // Fill in the values.
      for (unsigned int i_row = 0; i_row < (unsigned int)row_matrix_elementDOFs.numRows(); ++i_row)
        for (unsigned int i_col = 0; i_col < ele_centerlinedofindices.size(); ++i_col)
          row_matrix_elementDOFs(i_row, ele_centerlinedofindices[i_col]) =
              row_matrix_centerlineDOFs(i_row, i_col);
    }

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void ExtractPosDofVecAbsoluteValues(DRT::Discretization const& discret, DRT::Element const* ele,
        Teuchos::RCP<const Epetra_Vector> const& ia_discolnp,
        std::vector<double>& element_posdofvec_absolutevalues)
    {
      std::vector<double> eledispvec;

      // extract the Dof values of this element from displacement vector
      GetCurrentElementDis(discret, ele, ia_discolnp, eledispvec);

      DRT::ELEMENTS::Beam3Base const* beam_element_ptr =
          dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);

      if (beam_element_ptr != nullptr)
      {
        // get the current absolute values for those Dofs relevant for centerline interpolation
        // initial values are added by element itself
        beam_element_ptr->ExtractCenterlineDofValuesFromElementStateVector(
            eledispvec, element_posdofvec_absolutevalues, true);
      }
      else
      {
        element_posdofvec_absolutevalues = eledispvec;
        for (unsigned int dim = 0; dim < 3; ++dim)
          for (int node = 0; node < ele->NumNode(); ++node)
            element_posdofvec_absolutevalues[3 * node + dim] += ele->Nodes()[node]->X()[dim];
      }
    }

    /*-----------------------------------------------------------------------------*
     *-----------------------------------------------------------------------------*/
    void ExtractPosDofVecValues(DRT::Discretization const& discret, DRT::Element const* ele,
        Teuchos::RCP<const Epetra_Vector> const& ia_discolnp,
        std::vector<double>& element_posdofvec_values)
    {
      std::vector<double> eledispvec;

      // extract the Dof values of this element from displacement vector
      GetCurrentElementDis(discret, ele, ia_discolnp, eledispvec);

      DRT::ELEMENTS::Beam3Base const* beam_element_ptr =
          dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);

      if (beam_element_ptr != nullptr)
      {
        // get the current  values for those Dofs relevant for centerline interpolation
        beam_element_ptr->ExtractCenterlineDofValuesFromElementStateVector(
            eledispvec, element_posdofvec_values, false);
      }
      else
      {
        // for non beam elements, just return the element Dof vector
        element_posdofvec_values = eledispvec;
      }
    }

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    template <class T1, class T2>
    void ApplyBindingSpotForceToParentElements(DRT::Discretization const& discret,
        Teuchos::RCP<CORE::GEO::MESHFREE::BoundingBox> const& pbb,
        const Teuchos::RCP<Epetra_Vector> disp_np_col,
        const Teuchos::RCP<BEAMINTERACTION::BeamLink> elepairptr,
        std::vector<CORE::LINALG::SerialDenseVector> const& bspotforce,
        std::vector<CORE::LINALG::SerialDenseVector>& eleforce)
    {
      // auxiliary transformation matrix, will be resized and reused
      CORE::LINALG::SerialDenseMatrix trafomatrix;

      T1* cast_ele1 = dynamic_cast<T1*>(discret.gElement(elepairptr->GetEleGid(0)));
      T2* cast_ele2 = dynamic_cast<T2*>(discret.gElement(elepairptr->GetEleGid(1)));

      for (unsigned int elei = 0; elei < 2; ++elei)
      {
        // get current element displacements
        std::vector<double> eledisp;
        GetCurrentUnshiftedElementDis(
            discret, discret.gElement(elepairptr->GetEleGid(elei)), disp_np_col, *pbb, eledisp);
        const int numdof_ele = eledisp.size();

        // zero out and set correct size of transformation matrix
        trafomatrix.shape(6, numdof_ele);

        // I_variations
        if (elei == 0)
          cast_ele1->GetGeneralizedInterpolationMatrixVariationsAtXi(trafomatrix,
              cast_ele1->GetBindingSpotXi(
                  elepairptr->GetLinkerType(), elepairptr->GetLocBSpotNum(elei)),
              eledisp);
        else
          cast_ele2->GetGeneralizedInterpolationMatrixVariationsAtXi(trafomatrix,
              cast_ele2->GetBindingSpotXi(
                  elepairptr->GetLinkerType(), elepairptr->GetLocBSpotNum(elei)),
              eledisp);

        eleforce[elei].size(numdof_ele);
        CORE::LINALG::multiplyTN(eleforce[elei], trafomatrix, bspotforce[elei]);
      }
    }

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    template <class T1, class T2>
    void ApplyBindingSpotStiffToParentElements(DRT::Discretization const& discret,
        Teuchos::RCP<CORE::GEO::MESHFREE::BoundingBox> const& pbb,
        const Teuchos::RCP<Epetra_Vector> disp_np_col,
        const Teuchos::RCP<BEAMINTERACTION::BeamLink> elepairptr,
        std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> const& bspotstiff,
        std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>>& elestiff)
    {
      // Todo grill 02/17
      dserror(
          "we can not evaluate the tangent stiffness matrix without evaluating "
          "the residual in case of Kirchhoff beam Beam3k. This is because we have a coupled "
          "term and the residual vector is required to compute the full linearization. "
          "Use the combined evaluation method instead!");

      // todo: put this in a loop
      DRT::Element* ele1 = discret.gElement(elepairptr->GetEleGid(0));
      DRT::Element* ele2 = discret.gElement(elepairptr->GetEleGid(1));

      // get current element displacements
      std::vector<double> ele1disp;
      GetCurrentUnshiftedElementDis(discret, ele1, disp_np_col, *pbb, ele1disp);
      const int numdof_ele1 = ele1disp.size();

      std::vector<double> ele2disp;
      GetCurrentUnshiftedElementDis(discret, ele2, disp_np_col, *pbb, ele2disp);
      const int numdof_ele2 = ele2disp.size();

      // transformation matrix, will be resized and reused for various needs
      CORE::LINALG::SerialDenseMatrix trafomatrix;

      // auxiliary matrices required to store intermediate results after first of two
      // consecutive matrix-matrix products for each of four stiffness matrices
      std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> auxmat(
          2, std::vector<CORE::LINALG::SerialDenseMatrix>(2));

      T1* cast_ele1 = dynamic_cast<T1*>(ele1);
      T2* cast_ele2 = dynamic_cast<T2*>(ele2);

      // element 1:
      {
        // zero out and set correct size of transformation matrix
        trafomatrix.shape(6, numdof_ele1);

        // i) I_variations
        cast_ele1->GetGeneralizedInterpolationMatrixVariationsAtXi(trafomatrix,
            cast_ele1->GetBindingSpotXi(elepairptr->GetLinkerType(), elepairptr->GetLocBSpotNum(0)),
            ele1disp);

        auxmat[0][0].shape(numdof_ele1, 6);
        CORE::LINALG::multiplyTN(auxmat[0][0], trafomatrix, bspotstiff[0][0]);

        auxmat[0][1].shape(numdof_ele1, 6);
        CORE::LINALG::multiplyTN(auxmat[0][1], trafomatrix, bspotstiff[0][1]);

        // ii) I_increments
        trafomatrix.shape(6, numdof_ele1);

        cast_ele1->GetGeneralizedInterpolationMatrixIncrementsAtXi(trafomatrix,
            cast_ele1->GetBindingSpotXi(elepairptr->GetLinkerType(), elepairptr->GetLocBSpotNum(0)),
            ele1disp);

        elestiff[0][0].shape(numdof_ele1, numdof_ele1);
        CORE::LINALG::multiply(elestiff[0][0], auxmat[0][0], trafomatrix);

        auxmat[1][0].shape(6, numdof_ele1);
        CORE::LINALG::multiply(auxmat[1][0], bspotstiff[1][0], trafomatrix);
      }

      // element 2:
      {
        // i) I_variations
        trafomatrix.shape(6, numdof_ele2);

        cast_ele2->GetGeneralizedInterpolationMatrixVariationsAtXi(trafomatrix,
            cast_ele2->GetBindingSpotXi(elepairptr->GetLinkerType(), elepairptr->GetLocBSpotNum(1)),
            ele2disp);

        elestiff[1][0].shape(numdof_ele2, numdof_ele1);
        CORE::LINALG::multiplyTN(elestiff[1][0], trafomatrix, auxmat[1][0]);

        auxmat[1][1].shape(numdof_ele2, 6);
        CORE::LINALG::multiplyTN(auxmat[1][1], trafomatrix, bspotstiff[1][1]);

        // ii) I_increments
        trafomatrix.shape(6, numdof_ele2);

        cast_ele2->GetGeneralizedInterpolationMatrixIncrementsAtXi(trafomatrix,
            cast_ele2->GetBindingSpotXi(elepairptr->GetLinkerType(), elepairptr->GetLocBSpotNum(1)),
            ele2disp);

        elestiff[0][1].shape(numdof_ele1, numdof_ele2);
        CORE::LINALG::multiply(elestiff[0][1], auxmat[1][0], trafomatrix);

        elestiff[1][1].shape(numdof_ele2, numdof_ele2);
        CORE::LINALG::multiply(elestiff[1][1], auxmat[1][1], trafomatrix);
      }
    }

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    template <class T1, class T2>
    void ApplyBindingSpotForceStiffToParentElements(DRT::Discretization const& discret,
        Teuchos::RCP<CORE::GEO::MESHFREE::BoundingBox> const& pbb,
        const Teuchos::RCP<Epetra_Vector> disp_np_col,
        const Teuchos::RCP<BEAMINTERACTION::BeamLink> elepairptr,
        std::vector<CORE::LINALG::SerialDenseVector> const& bspotforce,
        std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> const& bspotstiff,
        std::vector<CORE::LINALG::SerialDenseVector>& eleforce,
        std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>>& elestiff)
    {
      // apply force on binding spots and corresponding linearizations to parent elements
      // and get discrete element force vectors and stiffness matrices

      /* Of course, the underlying idea is much more general:
       * Consider some kind of abstract interaction between two points \xi_1,\xi_2 \in [-1,1]
       * on two beam elements (or the same beam element).
       * We computed a discrete force/moment exerted on each of the two points; in addition,
       * we computed the linearizations with respect to the position/rotation vector of
       * this point and also the linearizations with respect to the position/rotation vector
       * of the other point (in general, the forces/moments depend on position/rotation
       * vector of both points).
       * Starting from here, we calculate the generalized interpolation matrices for
       * variations and increments of the position/rotation vectors at the two points
       * by expressing their dependency on the primary nodal DoFs of the corresponding
       * element (this is element specific, of course). Finally, we can 'transform' the
       * force vectors and stiffness matrices to discrete element force vectors and
       * stiffness matrices */

      // todo @grill: put this in a loop
      DRT::Element* ele1 = discret.gElement(elepairptr->GetEleGid(0));
      DRT::Element* ele2 = discret.gElement(elepairptr->GetEleGid(1));

      // get current element displacements
      std::vector<double> ele1disp;
      GetCurrentUnshiftedElementDis(discret, ele1, disp_np_col, *pbb, ele1disp);
      const int numdof_ele1 = ele1disp.size();

      std::vector<double> ele2disp;
      GetCurrentUnshiftedElementDis(discret, ele2, disp_np_col, *pbb, ele2disp);
      const int numdof_ele2 = ele2disp.size();

      // transformation matrix, will be resized and reused for various needs
      CORE::LINALG::SerialDenseMatrix trafomatrix;

      // auxiliary matrices required to store intermediate results after first of two
      // consecutive matrix-matrix products for each of four stiffness matrices
      std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> auxmat(
          2, std::vector<CORE::LINALG::SerialDenseMatrix>(2));

      // contribution to stiffmat from linearization of generalized interpolation matrix for
      // variations
      CORE::LINALG::SerialDenseMatrix stiffmat_lin_Ivar;

      // zero out and set correct size of transformation matrix
      trafomatrix.shape(6, numdof_ele1);

      T1* cast_ele1 = dynamic_cast<T1*>(ele1);
      T2* cast_ele2 = dynamic_cast<T2*>(ele2);

      // todo:
      // element 1:
      {
        // i) I_variations
        cast_ele1->GetGeneralizedInterpolationMatrixVariationsAtXi(trafomatrix,
            cast_ele1->GetBindingSpotXi(elepairptr->GetLinkerType(), elepairptr->GetLocBSpotNum(0)),
            ele1disp);

        eleforce[0].size(numdof_ele1);
        CORE::LINALG::multiplyTN(eleforce[0], trafomatrix, bspotforce[0]);

        auxmat[0][0].shape(numdof_ele1, 6);
        CORE::LINALG::multiplyTN(auxmat[0][0], trafomatrix, bspotstiff[0][0]);

        auxmat[0][1].shape(numdof_ele1, 6);
        CORE::LINALG::multiplyTN(auxmat[0][1], trafomatrix, bspotstiff[0][1]);

        // ii) I_increments
        trafomatrix.shape(6, numdof_ele1);

        cast_ele1->GetGeneralizedInterpolationMatrixIncrementsAtXi(trafomatrix,
            cast_ele1->GetBindingSpotXi(elepairptr->GetLinkerType(), elepairptr->GetLocBSpotNum(0)),
            ele1disp);

        elestiff[0][0].shape(numdof_ele1, numdof_ele1);
        CORE::LINALG::multiply(elestiff[0][0], auxmat[0][0], trafomatrix);

        auxmat[1][0].shape(6, numdof_ele1);
        CORE::LINALG::multiply(auxmat[1][0], bspotstiff[1][0], trafomatrix);


        // additional contribution from linearization of generalized interpolation matrix for
        // variations
        stiffmat_lin_Ivar.shape(numdof_ele1, numdof_ele1);

        cast_ele1->GetStiffmatResultingFromGeneralizedInterpolationMatrixAtXi(stiffmat_lin_Ivar,
            cast_ele1->GetBindingSpotXi(elepairptr->GetLinkerType(), elepairptr->GetLocBSpotNum(0)),
            ele1disp, bspotforce[0]);

        CORE::LINALG::Update(1.0, stiffmat_lin_Ivar, 1.0, elestiff[0][0]);
      }

      // element 2
      {
        // i) I_variations
        trafomatrix.shape(6, numdof_ele2);

        cast_ele2->GetGeneralizedInterpolationMatrixVariationsAtXi(trafomatrix,
            cast_ele2->GetBindingSpotXi(elepairptr->GetLinkerType(), elepairptr->GetLocBSpotNum(1)),
            ele2disp);

        eleforce[1].size(numdof_ele2);
        CORE::LINALG::multiplyTN(eleforce[1], trafomatrix, bspotforce[1]);

        elestiff[1][0].shape(numdof_ele2, numdof_ele1);
        CORE::LINALG::multiplyTN(elestiff[1][0], trafomatrix, auxmat[1][0]);

        auxmat[1][1].shape(numdof_ele2, 6);
        CORE::LINALG::multiplyTN(auxmat[1][1], trafomatrix, bspotstiff[1][1]);

        // ii) I_increments
        trafomatrix.shape(6, numdof_ele2);

        cast_ele2->GetGeneralizedInterpolationMatrixIncrementsAtXi(trafomatrix,
            cast_ele2->GetBindingSpotXi(elepairptr->GetLinkerType(), elepairptr->GetLocBSpotNum(1)),
            ele2disp);

        elestiff[0][1].shape(numdof_ele1, numdof_ele2);
        CORE::LINALG::multiply(elestiff[0][1], auxmat[0][1], trafomatrix);

        elestiff[1][1].shape(numdof_ele2, numdof_ele2);
        CORE::LINALG::multiply(elestiff[1][1], auxmat[1][1], trafomatrix);


        // additional contribution from linearization of generalized interpolation matrix for
        // variations
        stiffmat_lin_Ivar.shape(numdof_ele2, numdof_ele2);

        cast_ele2->GetStiffmatResultingFromGeneralizedInterpolationMatrixAtXi(stiffmat_lin_Ivar,
            cast_ele2->GetBindingSpotXi(elepairptr->GetLinkerType(), elepairptr->GetLocBSpotNum(1)),
            ele2disp, bspotforce[1]);

        CORE::LINALG::Update(1.0, stiffmat_lin_Ivar, 1.0, elestiff[1][1]);
      }
    }

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void SetupEleTypeMapExtractor(Teuchos::RCP<const DRT::Discretization> const& discret,
        Teuchos::RCP<CORE::LINALG::MultiMapExtractor>& eletypeextractor)
    {
      std::vector<std::set<int>> eletypeset(3);

      for (int i = 0; i < discret->NumMyRowElements(); ++i)
      {
        // get ele pointer
        DRT::Element* eleptr = discret->lRowElement(i);

        if (dynamic_cast<DRT::ELEMENTS::Beam3Base const*>(eleptr) != nullptr)
        {
          eletypeset[0].insert(eleptr->Id());
        }
        else if (dynamic_cast<DRT::ELEMENTS::Rigidsphere const*>(eleptr) != nullptr)
        {
          eletypeset[1].insert(eleptr->Id());
        }
        else if (dynamic_cast<DRT::ELEMENTS::So_base const*>(eleptr) != nullptr ||
                 dynamic_cast<DRT::ELEMENTS::Solid const*>(eleptr) != nullptr)
        {
          eletypeset[2].insert(eleptr->Id());
        }
        else
        {
          dserror("eletype multi map extractor cannot yet handle current element type.");
        }
      }

      std::vector<Teuchos::RCP<const Epetra_Map>> maps(eletypeset.size());
      for (int i = 0; i < static_cast<int>(eletypeset.size()); ++i)
      {
        std::vector<int> mapvec(eletypeset[i].begin(), eletypeset[i].end());
        eletypeset[i].clear();
        maps[i] =
            Teuchos::rcp(new Epetra_Map(-1, mapvec.size(), mapvec.data(), 0, discret->Comm()));
      }

      eletypeextractor->Setup(*discret()->ElementRowMap(), maps);
    }

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    void UpdateDofMapOfVector(Teuchos::RCP<DRT::Discretization> discret,
        Teuchos::RCP<Epetra_Vector>& dofmapvec, Teuchos::RCP<Epetra_Vector> old)
    {
      if (dofmapvec != Teuchos::null)
      {
        if (old == Teuchos::null) old = dofmapvec;
        dofmapvec = CORE::LINALG::CreateVector(*discret->DofRowMap(), true);
        CORE::LINALG::Export(*old, *dofmapvec);
      }
    }

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    long long CantorPairing(std::pair<int, int> const& pair)
    {
      long long z = 0.5 * (pair.first + pair.second) * (pair.first + pair.second + 1) + pair.second;

#ifdef BACI_DEBUG
      if (z > std::numeric_limits<long long>::max())
        dserror(" Your cantor paired value exceeds limit of data type int.");
      if (pair != CantorDePairing(z))
        dserror(
            " %i and %i cannot be paired using Cantor pairing function", pair.first, pair.second);
#endif

      return z;
    }

    /*----------------------------------------------------------------------------*
     *----------------------------------------------------------------------------*/
    std::pair<int, int> CantorDePairing(long long z)
    {
#ifdef BACI_DEBUG
      if (8.0 * z > std::numeric_limits<long long>::max())
        dserror(" Your cantor paired value exceeds limit of data type int.");
#endif

      long long w = std::floor((std::sqrt(8.0 * z + 1.0) - 1.0) * 0.5);
      long long t = (w + 1) * w * 0.5;

      std::pair<int, int> pair;
      pair.second = z - t;
      pair.first = w - pair.second;

      return pair;
    }


    //-----------------------------------------------------------------------------
    // explicit template instantiation (to please every compiler)
    //-----------------------------------------------------------------------------
    template void SetFilamentBindingSpotPositions(
        Teuchos::RCP<DRT::Discretization>, Teuchos::RCP<BEAMINTERACTION::CrosslinkingParams>);
    template void SetFilamentBindingSpotPositions(
        Teuchos::RCP<DRT::Discretization>, Teuchos::RCP<BEAMINTERACTION::SphereBeamLinkingParams>);

    template void
    ApplyBindingSpotForceToParentElements<DRT::ELEMENTS::Beam3Base, DRT::ELEMENTS::Beam3Base>(
        DRT::Discretization const&, Teuchos::RCP<CORE::GEO::MESHFREE::BoundingBox> const&,
        const Teuchos::RCP<Epetra_Vector>, const Teuchos::RCP<BEAMINTERACTION::BeamLink>,
        std::vector<CORE::LINALG::SerialDenseVector> const&,
        std::vector<CORE::LINALG::SerialDenseVector>&);
    template void
    ApplyBindingSpotForceToParentElements<DRT::ELEMENTS::Rigidsphere, DRT::ELEMENTS::Beam3Base>(
        DRT::Discretization const&, Teuchos::RCP<CORE::GEO::MESHFREE::BoundingBox> const&,
        const Teuchos::RCP<Epetra_Vector>, const Teuchos::RCP<BEAMINTERACTION::BeamLink>,
        std::vector<CORE::LINALG::SerialDenseVector> const&,
        std::vector<CORE::LINALG::SerialDenseVector>&);

    template void
    ApplyBindingSpotStiffToParentElements<DRT::ELEMENTS::Beam3Base, DRT::ELEMENTS::Beam3Base>(
        DRT::Discretization const&, Teuchos::RCP<CORE::GEO::MESHFREE::BoundingBox> const&,
        const Teuchos::RCP<Epetra_Vector>, const Teuchos::RCP<BEAMINTERACTION::BeamLink>,
        std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> const&,
        std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>>&);
    template void
    ApplyBindingSpotStiffToParentElements<DRT::ELEMENTS::Rigidsphere, DRT::ELEMENTS::Beam3Base>(
        DRT::Discretization const&, Teuchos::RCP<CORE::GEO::MESHFREE::BoundingBox> const&,
        const Teuchos::RCP<Epetra_Vector>, const Teuchos::RCP<BEAMINTERACTION::BeamLink>,
        std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> const&,
        std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>>&);

    template void
    ApplyBindingSpotForceStiffToParentElements<DRT::ELEMENTS::Beam3Base, DRT::ELEMENTS::Beam3Base>(
        DRT::Discretization const&, Teuchos::RCP<CORE::GEO::MESHFREE::BoundingBox> const&,
        const Teuchos::RCP<Epetra_Vector>, const Teuchos::RCP<BEAMINTERACTION::BeamLink>,
        std::vector<CORE::LINALG::SerialDenseVector> const&,
        std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> const&,
        std::vector<CORE::LINALG::SerialDenseVector>&,
        std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>>&);
    template void ApplyBindingSpotForceStiffToParentElements<DRT::ELEMENTS::Rigidsphere,
        DRT::ELEMENTS::Beam3Base>(DRT::Discretization const&,
        Teuchos::RCP<CORE::GEO::MESHFREE::BoundingBox> const&, const Teuchos::RCP<Epetra_Vector>,
        const Teuchos::RCP<BEAMINTERACTION::BeamLink>,
        std::vector<CORE::LINALG::SerialDenseVector> const&,
        std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> const&,
        std::vector<CORE::LINALG::SerialDenseVector>&,
        std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>>&);

    template void GetElementCenterlineGIDIndices<6>(
        DRT::Discretization const&, const DRT::Element*, CORE::LINALG::Matrix<6, 1, int>&);
    template void GetElementCenterlineGIDIndices<12>(
        DRT::Discretization const&, const DRT::Element*, CORE::LINALG::Matrix<12, 1, int>&);

  }  // namespace UTILS
}  // namespace BEAMINTERACTION

BACI_NAMESPACE_CLOSE

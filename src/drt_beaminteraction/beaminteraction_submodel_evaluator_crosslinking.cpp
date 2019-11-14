/*-----------------------------------------------------------*/
/*! \file
\brief class for submodel crosslinking

\maintainer Jonas Eichinger

\level 3
*/
/*-----------------------------------------------------------*/

#include "../drt_beaminteraction/beaminteraction_submodel_evaluator_crosslinking.H"
#include "../drt_beaminteraction/beaminteraction_calc_utils.H"
#include "../drt_beaminteraction/beaminteraction_data.H"
#include "../drt_beaminteraction/drt_utils_parallel_proctoproc.H"
#include "../drt_beaminteraction/crosslinking_params.H"
#include "../drt_beaminteraction/str_model_evaluator_beaminteraction_datastate.H"
#include "../drt_beaminteraction/crosslinker_node.H"
#include "../drt_beaminteraction/periodic_boundingbox.H"
#include "../drt_beaminteraction/beam_link.H"
#include "../drt_beaminteraction/beam_link_beam3r_lin2_rigidjointed.H"
#include "../drt_beaminteraction/beam_link_beam3r_lin2_pinjointed.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"
#include "../drt_io/runtime_vtp_writer.H"
#include "../drt_io/discretization_runtime_vtp_writer.H"

#include "../drt_structure_new/str_timint_basedataglobalstate.H"
#include "../drt_structure_new/str_timint_basedataio.H"
#include "../drt_structure_new/str_timint_basedataio_runtime_vtp_output.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_geometry/intersection_math.H"

#include "../drt_binstrategy/drt_meshfree_multibin.H"

#include "../drt_mat/crosslinkermat.H"

#include "../drt_beam3/beam3_base.H"

#include <Teuchos_TimeMonitor.hpp>

#include <unordered_set>
#include "beam_crosslinker_handler.H"


/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::Crosslinking()
    : crosslinking_params_ptr_(Teuchos::null),
      cl_exporter_(Teuchos::null),
      beam_exporter_(Teuchos::null),
      vtp_writer_ptr_(Teuchos::null),
      linker_disnp_(Teuchos::null),
      dis_at_last_redistr_(Teuchos::null),
      half_interaction_distance_(0.0),
      cl_noderowmap_prior_redistr_(Teuchos::null),
      cl_nodecolmap_prior_redistr_(Teuchos::null),
      beam_elerowmap_prior_redistr_(Teuchos::null),
      beam_elecolmap_prior_redistr_(Teuchos::null)
{
  crosslinker_data_.clear();
  beam_data_.clear();
  doublebondcl_.clear();
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::Setup()
{
  CheckInit();

  // construct, init and setup data container for crosslinking
  crosslinking_params_ptr_ = Teuchos::rcp(new BEAMINTERACTION::CrosslinkingParams());
  crosslinking_params_ptr_->Init(GState());
  crosslinking_params_ptr_->Setup();

  // set binding spot positions on filament elements according input file specifications
  SetFilamentTypes();
  // this includes temporary change in ghosting
  BEAMINTERACTION::UTILS::SetFilamentBindingSpotPositions(DiscretPtr(), crosslinking_params_ptr_);

  // add crosslinker to bin discretization
  AddCrosslinkerToBinDiscretization();

  // build runtime vtp writer
  if (GInOutput().GetRuntimeVtpOutputParams() != Teuchos::null) InitOutputRuntimeVtpStructure();

  // store old maps prior to redistribution
  cl_noderowmap_prior_redistr_ = Teuchos::rcp(new Epetra_Map(*BinDiscret().NodeRowMap()));
  cl_nodecolmap_prior_redistr_ = Teuchos::rcp(new Epetra_Map(*BinDiscret().NodeColMap()));
  beam_elerowmap_prior_redistr_ = Teuchos::rcp(new Epetra_Map(*EleTypeMapExtractor().BeamMap()));
  beam_elecolmap_prior_redistr_ = Teuchos::rcp(new Epetra_Map(*Discret().ElementColMap()));

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PostPartitionProblem()
{
  CheckInitSetup();

  // init beam data container
  unsigned int numcolele = Discret().NumMyColElements();
  beam_data_.clear();
  beam_data_.resize(Discret().NumMyColElements());
  for (unsigned int i = 0; i < numcolele; ++i)
  {
    // beam element i for which data will be collected
    DRT::ELEMENTS::Beam3Base* beamele_i =
        dynamic_cast<DRT::ELEMENTS::Beam3Base*>(DiscretPtr()->lColElement(i));

    // go to next element in case the current one is not a beam element
    if (beamele_i == NULL) continue;

    std::vector<double> eledisp;
    BEAMINTERACTION::UTILS::GetCurrentUnshiftedElementDis(Discret(), beamele_i,
        BeamInteractionDataStatePtr()->GetMutableDisColNp(), PeriodicBoundingBox(), eledisp);

    beam_data_[i] = Teuchos::rcp(new BEAMINTERACTION::DATA::BeamData());
    beam_data_[i]->SetId(beamele_i->Id());

    // loop over all binding spots of current element
    // loop over binding spot types of current element
    for (auto const& iter : beamele_i->GetBindingSpots())
    {
      // loop over all binding spots of current type j of current element
      int unsigned const numbbspot = beamele_i->GetNumberOfBindingSpots(iter.first);
      LINALG::Matrix<3, 1> pos(true);
      LINALG::Matrix<3, 3> triad(true);
      for (int unsigned k = 0; k < numbbspot; ++k)
      {
        BEAMINTERACTION::UTILS::GetPosAndTriadOfBindingSpot(
            beamele_i, PeriodicBoundingBoxPtr(), iter.first, k, pos, triad, eledisp);

        beam_data_[i]->SetBSpotPosition(iter.first, k, pos);
        beam_data_[i]->SetBSpotTriad(iter.first, k, triad);
        beam_data_[i]->SetBSpotStatus(iter.first, k, std::set<int>());
      }
    }
  }

  std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::CrosslinkerData>> newlinker;
  // map key is crosslinker gid to be able to uniquely address one entry over all procs
  std::map<int, NewDoubleBonds> mynewdbondcl;
  SetAllPossibleInitialDoubleBondedCrosslinker(newlinker, mynewdbondcl);

  // set row map of newly created linker discretization
  BinDiscretPtr()->FillComplete(false, false, false);

  // init crosslinker data container
  crosslinker_data_.clear();
  unsigned int numrowcl = BinDiscret().NumMyRowNodes();
  crosslinker_data_.resize(BinDiscret().NumMyColNodes());
  for (unsigned int i = 0; i < numrowcl; ++i)
  {
    DRT::Node* cl = BinDiscret().lRowNode(i);
    crosslinker_data_[cl->LID()] = Teuchos::rcp(new BEAMINTERACTION::DATA::CrosslinkerData());

    crosslinker_data_[cl->LID()]->SetId(cl->Id());
  }

  // set initially set crosslinker
  for (auto& iter : newlinker)
  {
    int cl_lid = BinDiscret().NodeColMap()->LID(iter->GetId());
    crosslinker_data_[cl_lid]->SetBSpots(iter->GetBSpots());
    crosslinker_data_[cl_lid]->SetNumberOfBonds(iter->GetNumberOfBonds());
    crosslinker_data_[cl_lid]->SetPosition(iter->GetPosition());
  }

  // setup new double bonds and insert them in doublebondcl_
  CreateNewDoubleBondedCrosslinkerElementPairs(mynewdbondcl);

  // store maps
  StoreMapsPriorRedistribution();
  UpdateAndExportCrosslinkerData();
  UpdateAndExportBeamData(false);

  // local flag if one proc has new linker
  int loc_newlinks = static_cast<int>(newlinker.size() > 0);
  // global flag
  int g_newlinks = 0;
  BinDiscret().Comm().MaxAll(&loc_newlinks, &g_newlinks, 1);

  newlinker.clear();

  // we need to partition our problem again
  return static_cast<bool>(g_newlinks);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PostSetup()
{
  CheckInitSetup();

  if (not DRT::Problem::Instance()->Restart())
  {
    // in case of initially set crosslinker
    if (crosslinking_params_ptr_->TotalNumInitCrosslinker() > 0) UpdateMyDoubleBondsRemoteIdList();

    // store displacement of restart step as displacement state of last redistribution
    dis_at_last_redistr_ = Teuchos::rcp(new Epetra_Vector(*BinDiscret().DofRowMap(), true));
    for (int i = 0; i < BinDiscret().NumMyRowNodes(); ++i)
    {
      CROSSLINKING::CrosslinkerNode* crosslinker_i =
          dynamic_cast<CROSSLINKING::CrosslinkerNode*>(BinDiscret().lRowNode(i));

      // std::vector holding gids of dofs
      std::vector<int> dofnode = BinDiscret().Dof(crosslinker_i);

      // loop over all dofs
      for (unsigned int dim = 0; dim < 3; ++dim)
      {
        int doflid = dis_at_last_redistr_->Map().LID(dofnode[dim]);
        (*dis_at_last_redistr_)[doflid] = crosslinker_i->X()[dim];
      }
    }
    // build up column linker information
    UpdateAndExportCrosslinkerData();
    UpdateAndExportBeamData(false);
  }

  // store maps
  StoreMapsPriorRedistribution();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::InitSubmodelDependencies(
    Teuchos::RCP<STR::MODELEVALUATOR::BeamInteraction::Map> const submodelmap)
{
  CheckInitSetup();
  // no active influence on other submodels
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::SetFilamentTypes()
{
  CheckInit();

  std::set<int> examined_fils;

  // loop over all col nodes
  for (int coln = 0; coln < Discret().NumMyColNodes(); ++coln)
  {
    DRT::Node* currnode = Discret().lColNode(coln);
    // get filament number of current node ( requirement: node belongs to only one filament)
    DRT::Condition* cond = currnode->GetCondition("BeamLineFilamentCondition");

    // in case node (e.g. node of rigid sphere element) does not belong to a filament, go to next
    // node
    if (cond == NULL) continue;

    // get filament type
    INPAR::BEAMINTERACTION::FilamentType filtype =
        INPAR::BEAMINTERACTION::String2FilamentType(*(cond->Get<std::string>("Type")));

    for (int i = 0; i < currnode->NumElement(); ++i)
    {
      DRT::ELEMENTS::Beam3Base* beamele =
          dynamic_cast<DRT::ELEMENTS::Beam3Base*>(currnode->Elements()[i]);

#ifdef DEBUG
      if (beamele == NULL)
        dserror(" DESIGN LINE BEAM FILAMENT CONDITIONS only suitable for beam elements.");
#endif

      beamele->SetFilamentType(filtype);
    }
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::SetAllPossibleInitialDoubleBondedCrosslinker(
    std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::CrosslinkerData>>& newlinker,
    std::map<int, NewDoubleBonds>& mynewdbondcl)
{
  CheckInit();

  // in case no initial crosslinker are supposed to be set, return, no repartitioning
  // of discretization required in this case
  if (crosslinking_params_ptr_->TotalNumInitCrosslinker() == 0) return;

  // get all possible bspot partners
  std::vector<BEAMINTERACTION::DATA::BspotLinkerData> my_bspot_linker;
  std::map<int, std::vector<BEAMINTERACTION::DATA::BspotLinkerData>> global_bspot_linker;

  // loop over all binding spots and find matching ones
  GetAllPossibleBspotLinks(my_bspot_linker);

  // create same list on procs
  CommunicateInitialLinker(my_bspot_linker, global_bspot_linker);

  // all procs make decisions for a valid link based the same rule
  std::vector<int> newlinkermatid;
  UnambiguousDecisionsOnAllProcs(newlinker, global_bspot_linker, newlinkermatid);

  // setupt new double bonded linker
  SetupMyInitialDoubleBondedLinker(newlinker, mynewdbondcl, newlinkermatid);
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::GetAllPossibleBspotLinks(
    std::vector<BEAMINTERACTION::DATA::BspotLinkerData>& my_bspot_linker)
{
  CheckInit();

  my_bspot_linker.clear();

  // loop over all row beam elements
  int unsigned const numbeams = EleTypeMapExtractorPtr()->BeamMap()->NumMyElements();
  my_bspot_linker.reserve(numbeams);
  for (unsigned int rowbeam_i = 0; rowbeam_i < numbeams; ++rowbeam_i)
  {
    DRT::ELEMENTS::Beam3Base* beamele = dynamic_cast<DRT::ELEMENTS::Beam3Base*>(
        Discret().gElement(EleTypeMapExtractorPtr()->BeamMap()->GID(rowbeam_i)));

    BEAMINTERACTION::DATA::BeamData const* beamdata_i = beam_data_[beamele->LID()].get();

    // loop over all binding spot types of current filament
    for (auto const& iter : beamdata_i->GetBSpotStatus())
    {
      // loop over all binding spots of current binding spot type of current element
      unsigned int numbspots = beamele->GetNumberOfBindingSpots(iter.first);
      for (unsigned int locbspot_i = 0; locbspot_i < numbspots; ++locbspot_i)
      {
        // get bin of current bspot
        int const bingid =
            BinStrategy().ConvertPosToGid(beamdata_i->GetBSpotPosition(iter.first, locbspot_i));

        // get neighboring bins
        // note: interaction distance cl to beam needs to be smaller than the bin size
        std::vector<int> neighboring_binIds;
        neighboring_binIds.reserve(27);
        // do not check on existence here -> shifted to GetBinContent
        BinStrategyPtr()->GetNeighborAndOwnBinIds(bingid, neighboring_binIds);

        // get set of neighboring beam elements (i.e. elements that somehow touch nb bins)
        // we also need col elements (flag = false) here (in contrast to "normal" crosslinking)
        std::set<DRT::Element*> neighboring_beams;
        std::vector<BINSTRATEGY::UTILS::BinContentType> bc(1, BINSTRATEGY::UTILS::Beam);
        BinStrategyPtr()->GetBinContent(neighboring_beams, bc, neighboring_binIds, false);

        // in case there are no neighbors, go to next binding spot
        if (neighboring_beams.empty()) continue;

        // loop over all neighboring beam elements
        for (auto const& nb_ele_i : neighboring_beams)
        {
          DRT::ELEMENTS::Beam3Base* nb_beamele = dynamic_cast<DRT::ELEMENTS::Beam3Base*>(nb_ele_i);

          BEAMINTERACTION::DATA::BeamData const* nb_beamdata_i =
              beam_data_[nb_beamele->LID()].get();

          // exclude linking of touching elements
          if (BEAMINTERACTION::UTILS::DoBeamElementsShareNodes(beamele, nb_beamele)) continue;

          // loop over binding spots of neighboring element
          for (unsigned int nb_locbspot_i = 0;
               nb_locbspot_i < nb_beamele->GetNumberOfBindingSpots(iter.first); ++nb_locbspot_i)
          {
            // loop over different linker types and check for feasibility of link
            std::vector<int> const& matcrosslinkerpertype =
                crosslinking_params_ptr_->MatCrosslinkerPerType();
            for (unsigned int type_i = 0; type_i < matcrosslinkerpertype.size(); ++type_i)
            {
              // check criteria for feasible bond
              Teuchos::RCP<MAT::CrosslinkerMat> mat =
                  Teuchos::rcp_dynamic_cast<MAT::CrosslinkerMat>(
                      MAT::Material::Factory(matcrosslinkerpertype[type_i]));

              // linker type needs to match binding spot type
              if (mat->LinkerType() != iter.first) continue;

              // distance
              // minimum and maximum distance at which a double-bond crosslink can be established
              double const linkdistmin = mat->LinkingLength() - mat->LinkingLengthTolerance();
              double const linkdistmax = mat->LinkingLength() + mat->LinkingLengthTolerance();

#ifdef DEBUG
              if (linkdistmax > BinStrategy().CutoffRadius())
                dserror(
                    "The allowed binding distance of your linker material is greater than the "
                    "cutoff radius, this can lead to missed binding events");
#endif

              if (BEAMINTERACTION::UTILS::IsDistanceOutOfRange(
                      beamdata_i->GetBSpotPosition(iter.first, locbspot_i),
                      nb_beamdata_i->GetBSpotPosition(iter.first, nb_locbspot_i), linkdistmin,
                      linkdistmax))
                continue;

              // orientation of centerline tangent vectors at binding spots:
              // a crosslink (double-bonded crosslinker) will only be established if the
              // enclosed angle is in the specified range
              double const linkanglemin = mat->LinkingAngle() - mat->LinkingAngleTolerance();
              double const linkanglemax = mat->LinkingAngle() + mat->LinkingAngleTolerance();

              // get tangent of binding spot on beamele
              LINALG::Matrix<3, 1> bindingspot_beam_tangent(true);
              for (unsigned int idim = 0; idim < 3; ++idim)
                bindingspot_beam_tangent(idim) =
                    beamdata_i->GetBSpotTriad(iter.first, locbspot_i)(idim, 0);

              // get tangent of binding spot on nb_beamele
              LINALG::Matrix<3, 1> nb_bindingspot_beam_tangent(true);
              for (unsigned int idim = 0; idim < 3; ++idim)
                nb_bindingspot_beam_tangent(idim) =
                    nb_beamdata_i->GetBSpotTriad(iter.first, nb_locbspot_i)(idim, 0);

              if (BEAMINTERACTION::UTILS::IsEnclosedAngleOutOfRange(bindingspot_beam_tangent,
                      nb_bindingspot_beam_tangent, linkanglemin, linkanglemax))
                continue;

              // if all criteria were met for a certain linker, add this bond to linker specific
              // bond list
              BEAMINTERACTION::DATA::BspotLinkerData bspotlinkerdata;
              bspotlinkerdata.SetEleGid1(beamele->Id());
              bspotlinkerdata.SetLocBspotId1(locbspot_i);
              bspotlinkerdata.SetEleGid2(nb_beamele->Id());
              bspotlinkerdata.SetLocBspotId2(nb_locbspot_i);
              bspotlinkerdata.SetType(static_cast<int>(iter.first));
              bspotlinkerdata.SetMatId(matcrosslinkerpertype[type_i]);

              // append to all my linker
              my_bspot_linker.push_back(bspotlinkerdata);

              // go to next binding spot, if current one will possibly bind
              break;
            }
          }
        }
      }
    }
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::CommunicateInitialLinker(
    std::vector<BEAMINTERACTION::DATA::BspotLinkerData> const& my_bspot_linker,
    std::map<int, std::vector<BEAMINTERACTION::DATA::BspotLinkerData>>& global_bspot_linker)
{
  CheckInit();

  global_bspot_linker.clear();

  // gather all data over all procs
  Epetra_Comm const& com = Discret().Comm();
  const int numproc = com.NumProc();
  int numpairs = static_cast<int>(my_bspot_linker.size());
  std::vector<int> elegid_1, elegid_2, locbspot_1, locbspot_2, type, mat;
  for (int iproc = 0; iproc < numproc; ++iproc)
  {
    numpairs = my_bspot_linker.size();
    std::vector<int> elegid_1_i(numpairs, 0), elegid_2_i(numpairs, 0), locbspot_1_i(numpairs, 0),
        locbspot_2_i(numpairs, 0), type_i(numpairs, 0), mat_i(numpairs, 0);
    if (iproc == GState().GetMyRank())
    {
      for (unsigned int i = 0; i < my_bspot_linker.size(); ++i)
      {
        elegid_1_i[i] = my_bspot_linker[i].GetEleGid1();
        elegid_2_i[i] = my_bspot_linker[i].GetEleGid2();
        locbspot_1_i[i] = my_bspot_linker[i].GetLocBspotId1();
        locbspot_2_i[i] = my_bspot_linker[i].GetLocBspotId2();
        type_i[i] = my_bspot_linker[i].GetType();
        mat_i[i] = my_bspot_linker[i].GetMatId();
      }
    }

    // first: proc i tells all procs how many pairs it has
    com.Broadcast(&numpairs, 1, iproc);
    // second: proc i tells all procs which pairs it has
    elegid_1_i.resize(numpairs);
    elegid_2_i.resize(numpairs);
    locbspot_1_i.resize(numpairs);
    locbspot_2_i.resize(numpairs);
    type_i.resize(numpairs);
    mat_i.resize(numpairs);

    com.Broadcast(&elegid_1_i[0], numpairs, iproc);
    com.Broadcast(&elegid_2_i[0], numpairs, iproc);
    com.Broadcast(&locbspot_1_i[0], numpairs, iproc);
    com.Broadcast(&locbspot_2_i[0], numpairs, iproc);
    com.Broadcast(&type_i[0], numpairs, iproc);
    com.Broadcast(&mat_i[0], numpairs, iproc);

    elegid_1.insert(elegid_1.end(), elegid_1_i.begin(), elegid_1_i.end());
    elegid_2.insert(elegid_2.end(), elegid_2_i.begin(), elegid_2_i.end());
    locbspot_1.insert(locbspot_1.end(), locbspot_1_i.begin(), locbspot_1_i.end());
    locbspot_2.insert(locbspot_2.end(), locbspot_2_i.begin(), locbspot_2_i.end());
    type.insert(type.end(), type_i.begin(), type_i.end());
    mat.insert(mat.end(), mat_i.begin(), mat_i.end());
  }

  // rebuild bspot_linker_data
  int numglobalpairs = elegid_1.size();
  std::map<int, std::map<long long, BEAMINTERACTION::DATA::BspotLinkerData>> sort_data;
  long long bspotgid = -1;
  for (int i = 0; i < numglobalpairs; ++i)
  {
    BEAMINTERACTION::DATA::BspotLinkerData bspotlinkerdata;
    bspotlinkerdata.SetEleGid1(elegid_1[i]);
    bspotlinkerdata.SetLocBspotId1(locbspot_1[i]);
    bspotlinkerdata.SetEleGid2(elegid_2[i]);
    bspotlinkerdata.SetLocBspotId2(locbspot_2[i]);
    bspotlinkerdata.SetType(type[i]);
    bspotlinkerdata.SetMatId(mat[i]);

    bspotgid = BEAMINTERACTION::UTILS::CantorPairing(std::make_pair(elegid_2[i], locbspot_2[i]));

    sort_data[elegid_1[i]][bspotgid] = bspotlinkerdata;
  }

  for (auto const& iter_sort_i : sort_data)
    for (auto const& iter_sort_j : iter_sort_i.second)
      global_bspot_linker[iter_sort_i.first].push_back(iter_sort_j.second);
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UnambiguousDecisionsOnAllProcs(
    std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::CrosslinkerData>>& newlinker,
    std::map<int, std::vector<BEAMINTERACTION::DATA::BspotLinkerData>> const& global_bspot_linker,
    std::vector<int>& newlinkermatid)
{
  // initialize a box within linker are spawned
  std::vector<bool> dummy(3, false);
  Teuchos::RCP<GEO::MESHFREE::BoundingBox> linker_init_box =
      Teuchos::rcp(new GEO::MESHFREE::BoundingBox());
  linker_init_box->Init(
      crosslinking_params_ptr_->LinkerInitializationBox(), dummy);  // no Setup() call needed here

  // loop over bspotpairs and make decision
  newlinker.reserve(global_bspot_linker.size());
  std::map<long long, int> bondsperbindingspot;
  std::set<std::pair<long long, long long>, BEAMINTERACTION::UTILS::StdPairComparatorOrderCounts>
      doublebonds;
  std::map<int, int> numpertype;
  LINALG::Matrix<3, 1> clpos(true);
  for (auto const& iter_sort : global_bspot_linker)
  {
    for (auto const& iter : iter_sort.second)
    {
      long long bspotgid = BEAMINTERACTION::UTILS::CantorPairing(
          std::make_pair(iter.GetEleGid1(), iter.GetLocBspotId1()));
      long long nb_bspotgid = BEAMINTERACTION::UTILS::CantorPairing(
          std::make_pair(iter.GetEleGid2(), iter.GetLocBspotId2()));

      INPAR::BEAMINTERACTION::CrosslinkerType linkertype =
          static_cast<INPAR::BEAMINTERACTION::CrosslinkerType>(iter.GetType());

      // check if binding spot has reached its maximum number of bonds
      if (bondsperbindingspot.find(bspotgid) != bondsperbindingspot.end() and
          bondsperbindingspot[bspotgid] >=
              crosslinking_params_ptr_->MaxNumberOfBondsPerFilamentBspot(linkertype))
        continue;

      // check if maximal number of crosslinker for current type is already exceeded
      if ((numpertype[iter.GetMatId()] + 1) >
          crosslinking_params_ptr_->NumInitCrosslinkerPerCrosslinkerMatId(iter.GetMatId()))
        continue;

      // to ensure that double bonds are recognized as such independent of the order of the
      // bspotgids
      std::pair<long long, long long> currdoublebond = (bspotgid < nb_bspotgid)
                                                           ? std::make_pair(bspotgid, nb_bspotgid)
                                                           : std::make_pair(nb_bspotgid, bspotgid);
      // to check:
      // i)  check if binding spot has reached its maximum number of bonds
      // ii) check if identical bond already exists
      if ((bondsperbindingspot.find(nb_bspotgid) != bondsperbindingspot.end() and
              bondsperbindingspot[nb_bspotgid] >=
                  crosslinking_params_ptr_->MaxNumberOfBondsPerFilamentBspot(linkertype)) or
          doublebonds.find(currdoublebond) != doublebonds.end())
        continue;

      // update variables
      doublebonds.insert(currdoublebond);
      ++bondsperbindingspot[nb_bspotgid];
      ++bondsperbindingspot[bspotgid];
      ++numpertype[iter.GetMatId()];

      // check ownership
      if (Discret().ElementRowMap()->LID(iter.GetEleGid1()) == -1) continue;

      // store data of new crosslinker
      Teuchos::RCP<BEAMINTERACTION::DATA::CrosslinkerData> cldata =
          Teuchos::rcp(new BEAMINTERACTION::DATA::CrosslinkerData());

      // set positions
      clpos.Clear();
      SetPositionOfDoubleBondedCrosslinkerPBCconsistent(clpos,
          beam_data_[DiscretPtr()->gElement(iter.GetEleGid1())->LID()]->GetBSpotPosition(
              linkertype, iter.GetLocBspotId1()),
          beam_data_[DiscretPtr()->gElement(iter.GetEleGid2())->LID()]->GetBSpotPosition(
              linkertype, iter.GetLocBspotId2()));

      // check if linker is within linker initialization box
      if (not linker_init_box->Within(clpos, dummy)) continue;

      // set current binding spot status of crosslinker
      cldata->SetPosition(clpos);
      cldata->SetBspot(0, std::make_pair(iter.GetEleGid1(), iter.GetLocBspotId1()));
      cldata->SetBspot(1, std::make_pair(iter.GetEleGid2(), iter.GetLocBspotId2()));
      // set number of bonds
      cldata->SetNumberOfBonds(2);

      newlinker.push_back(cldata);
      newlinkermatid.push_back(iter.GetMatId());
    }
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::SetupMyInitialDoubleBondedLinker(
    std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::CrosslinkerData>>& newlinker,
    std::map<int, NewDoubleBonds>& mynewdbondcl, std::vector<int> const& newlinkermatid)
{
  // determine unique gids on each proc (ascending over all procs)
  // gather numbers of new linker on each proc
  Epetra_Comm const& com = Discret().Comm();
  std::vector<int> nummynewlinks(1);
  nummynewlinks[0] = static_cast<int>(newlinker.size());
  // initialize std::vector for communication
  std::vector<int> numnewlinks(com.NumProc(), 0);
  // communicate
  com.GatherAll(&nummynewlinks[0], &numnewlinks[0], nummynewlinks.size());
  com.Barrier();

  // calculate starting index on myrank
  int mystartgid = 0;
  for (int i = 0; i < GState().GetMyRank(); ++i) mystartgid += numnewlinks[i];

  // loop over new linker on myrank
  int gid = BinDiscretPtr()->NodeRowMap()->MaxAllGID() + 1 + mystartgid;
  std::vector<double> X(3);
  for (unsigned int i = 0; i < newlinker.size(); ++i)
  {
    for (unsigned int dim = 0; dim < 3; ++dim) X[dim] = newlinker[i]->GetPosition()(dim);

    Teuchos::RCP<CROSSLINKING::CrosslinkerNode> newcrosslinker =
        Teuchos::rcp(new CROSSLINKING::CrosslinkerNode(gid, &X[0], GState().GetMyRank()));
    newcrosslinker->SetMaterial(Teuchos::rcp_dynamic_cast<MAT::CrosslinkerMat>(
        MAT::Material::Factory(newlinkermatid[i])));  // HACK HACK HACK
    BinDiscretPtr()->AddNode(newcrosslinker);

    NewDoubleBonds dbondcl;
    std::vector<std::pair<int, int>> bspots = newlinker[i]->GetBSpots();
    dbondcl.id = gid;
    dbondcl.eleids.push_back(bspots[0]);
    dbondcl.eleids.push_back(bspots[1]);

    int colelelid = DiscretPtr()->ElementColMap()->LID(bspots[0].first);
    int nb_colelelid = DiscretPtr()->ElementColMap()->LID(bspots[1].first);
    dbondcl.bspotposs.push_back(beam_data_[colelelid]->GetBSpotPosition(
        newcrosslinker->GetMaterial()->LinkerType(), bspots[0].second));
    dbondcl.bspotposs.push_back(beam_data_[nb_colelelid]->GetBSpotPosition(
        newcrosslinker->GetMaterial()->LinkerType(), bspots[1].second));
    dbondcl.bspottriads.push_back(beam_data_[colelelid]->GetBSpotTriad(
        newcrosslinker->GetMaterial()->LinkerType(), bspots[0].second));
    dbondcl.bspottriads.push_back(beam_data_[nb_colelelid]->GetBSpotTriad(
        newcrosslinker->GetMaterial()->LinkerType(), bspots[1].second));

    mynewdbondcl[dbondcl.id] = dbondcl;

    // set correct states for linker
    newlinker[i]->SetId(gid);

    std::vector<double> newpos(3, 0.0);
    for (int dim = 0; dim < 3; ++dim) newpos[dim] = newlinker[i]->GetPosition()(dim);
    newcrosslinker->SetPos(newpos);

    beam_data_[DiscretPtr()->ElementColMap()->LID(newlinker[i]->GetBSpots()[0].first)]
        ->AddBondToBindingSpot(
            newcrosslinker->GetMaterial()->LinkerType(), newlinker[i]->GetBSpots()[0].second, gid);

    if (Discret().HaveGlobalElement(newlinker[i]->GetBSpots()[1].first))
      beam_data_[DiscretPtr()->ElementColMap()->LID(newlinker[i]->GetBSpots()[1].first)]
          ->AddBondToBindingSpot(newcrosslinker->GetMaterial()->LinkerType(),
              newlinker[i]->GetBSpots()[1].second, gid);

    // update gid
    ++gid;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::AddCrosslinkerToBinDiscretization()
{
  CheckInit();

  // initialize a box within linker are spawned
  std::vector<bool> dummy(3, false);
  Teuchos::RCP<GEO::MESHFREE::BoundingBox> linker_init_box =
      Teuchos::rcp(new GEO::MESHFREE::BoundingBox());
  linker_init_box->Init(
      crosslinking_params_ptr_->LinkerInitializationBox(), dummy);  // no Setup() call needed here

  // loop over all linker types that should be added to simulation volume
  // (only proc 0 is doing this (as the number of crosslinker is manageable)
  std::vector<int> const& numcrosslinkerpertype = crosslinking_params_ptr_->NumCrosslinkerPerType();
  std::vector<int> const& matcrosslinkerpertype = crosslinking_params_ptr_->MatCrosslinkerPerType();
  int gid = 0;
  if (GState().GetMyRank() == 0)
  {
    for (int cltype_i = 0; cltype_i < static_cast<int>(numcrosslinkerpertype.size()); ++cltype_i)
    {
      for (int cltype_i_cl_j = 0; cltype_i_cl_j < numcrosslinkerpertype[cltype_i]; ++cltype_i_cl_j)
      {
        // random reference position of crosslinker in bounding box
        LINALG::Matrix<3, 1> Xmat;
        linker_init_box->RandomPosWithin(Xmat);

        std::vector<double> X(3);
        for (int dim = 0; dim < 3; ++dim) X[dim] = Xmat(dim);

        // construct node, init data container, set material and add to bin discret
        Teuchos::RCP<CROSSLINKING::CrosslinkerNode> newcrosslinker =
            Teuchos::rcp(new CROSSLINKING::CrosslinkerNode(gid++, &X[0], GState().GetMyRank()));
        newcrosslinker->SetMaterial(matcrosslinkerpertype[cltype_i]);
        BinDiscretPtr()->AddNode(newcrosslinker);
      }
    }
  }

  // set row map of newly created linker discretization
  BinDiscretPtr()->FillComplete(false, false, false);
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::Reset()
{
  CheckInitSetup();

  // reset crosslinker pairs
  for (auto const& iter : doublebondcl_)
  {
    Teuchos::RCP<BEAMINTERACTION::BeamLink> elepairptr = iter.second;

#ifdef DEBUG

    CROSSLINKING::CrosslinkerNode* cl_i =
        dynamic_cast<CROSSLINKING::CrosslinkerNode*>(BinDiscretPtr()->gNode(elepairptr->Id()));
    // safety check
    BEAMINTERACTION::DATA::CrosslinkerData const* cldata_i = crosslinker_data_[cl_i->LID()].get();

    if (cldata_i->GetNumberOfBonds() != 2)
      dserror("Cl with gid %i Owner %i on myrank %i and numbonds %i", elepairptr->Id(),
          cl_i->Owner(), GStatePtr()->GetMyRank(), cldata_i->GetNumberOfBonds());
#endif

    // init positions and triads
    std::vector<LINALG::Matrix<3, 1>> pos(2);
    std::vector<LINALG::Matrix<3, 3>> triad(2);

    for (int i = 0; i < 2; ++i)
    {
      int elegid = elepairptr->GetEleGid(i);

      // safety check
      if (DiscretPtr()->ElementColMap()->LID(elegid) < 0)
      {
        elepairptr->Print(std::cout);
        dserror("Reset(): elegid %i not there on proc %i ", elegid, GState().GetMyRank());
      }

      int locbspotnum = elepairptr->GetLocBSpotNum(i);
      DRT::Element* ele = DiscretPtr()->gElement(elegid);

      BEAMINTERACTION::UTILS::GetPosAndTriadOfBindingSpot(Discret(), ele,
          BeamInteractionDataStatePtr()->GetMutableDisColNp(), PeriodicBoundingBoxPtr(),
          elepairptr->GetLinkerType(), locbspotnum, pos[i], triad[i]);
    }

    // unshift one of the positions if both are separated by a periodic boundary
    // condition, i.e. have been shifted before
    PeriodicBoundingBoxPtr()->UnShift3D(pos[1], pos[0]);

    // safety check until code is better tested for potential problems with periodic boundary
    // conditions
    // **************************** DEBUG ****************************************
    LINALG::Matrix<3, 1> dist(true);
    dist.Update(1.0, pos[0], -1.0, pos[1]);
    for (unsigned int i = 0; i < 3; ++i)
    {
      if (std::abs(dist(i)) > 0.5 * PeriodicBoundingBoxPtr()->EdgeLength(i))
      {
        dserror(
            "You are trying to set the binding spot positions of this crosslinker "
            "in at least one direction\n at a distance larger than %f, which is "
            " half of the period length in the respective direction",
            0.5 * PeriodicBoundingBoxPtr()->EdgeLength(i));
      }
    }
    // ********************** END DEBUG ****************************************

    // finally reset state
    elepairptr->ResetState(pos, triad);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::EvaluateForce()
{
  CheckInitSetup();

  // force and moment exerted on the two connection sites due to the mechanical connection
  std::vector<LINALG::SerialDenseVector> bspotforce(2, LINALG::SerialDenseVector(6));

  // resulting discrete element force vectors of the two parent elements
  std::vector<LINALG::SerialDenseVector> eleforce(2);

  std::vector<std::vector<LINALG::SerialDenseMatrix>> dummystiff;

  // element gids of interacting elements
  std::vector<int> elegids(2);

  for (auto const& iter : doublebondcl_)
  {
    Teuchos::RCP<BEAMINTERACTION::BeamLink> elepairptr = iter.second;

    for (int i = 0; i < 2; ++i)
    {
      elegids[i] = elepairptr->GetEleGid(i);
      bspotforce[i].Zero();
    }

    // evaluate beam linkage object to get forces and moments on binding spots
    elepairptr->EvaluateForce(bspotforce[0], bspotforce[1]);

    // apply forces on binding spots to parent elements
    // and get their discrete element force vectors
    BEAMINTERACTION::UTILS::ApplyBindingSpotForceToParentElements<DRT::ELEMENTS::Beam3Base,
        DRT::ELEMENTS::Beam3Base>(Discret(), PeriodicBoundingBoxPtr(),
        BeamInteractionDataStatePtr()->GetMutableDisColNp(), elepairptr, bspotforce, eleforce);

    // assemble the contributions into force vector class variable
    // f_crosslink_np_ptr_, i.e. in the DOFs of the connected nodes
    BEAMINTERACTION::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrix(Discret(), elegids,
        eleforce, dummystiff, BeamInteractionDataStatePtr()->GetMutableForceNp(), Teuchos::null);
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::EvaluateStiff()
{
  CheckInitSetup();

  /* linearizations, i.e. stiffness contributions due to forces on the two
   * connection sites due to the mechanical connection */
  std::vector<std::vector<LINALG::SerialDenseMatrix>> bspotstiff(
      2, std::vector<LINALG::SerialDenseMatrix>(2, LINALG::SerialDenseMatrix(6, 6)));

  // linearizations, i.e. discrete stiffness contributions to the two parent elements
  // we can't handle this separately for both elements because there are entries which
  // couple the two element stiffness blocks
  std::vector<std::vector<LINALG::SerialDenseMatrix>> elestiff(
      2, std::vector<LINALG::SerialDenseMatrix>(2));

  std::vector<LINALG::SerialDenseVector> dummyforce;

  // element gids of interacting elements
  std::vector<int> elegids(2);

  for (auto const& iter : doublebondcl_)
  {
    Teuchos::RCP<BEAMINTERACTION::BeamLink> elepairptr = iter.second;

    for (int i = 0; i < 2; ++i)
    {
      elegids[i] = elepairptr->GetEleGid(i);

      for (int j = 0; j < 2; ++j) bspotstiff[i][j].Zero();
    }

    // evaluate beam linkage object to get linearizations of forces and moments on binding spots
    elepairptr->EvaluateStiff(
        bspotstiff[0][0], bspotstiff[0][1], bspotstiff[1][0], bspotstiff[1][1]);

    // apply linearizations to parent elements and get their discrete element stiffness matrices
    BEAMINTERACTION::UTILS::ApplyBindingSpotStiffToParentElements<DRT::ELEMENTS::Beam3Base,
        DRT::ELEMENTS::Beam3Base>(Discret(), PeriodicBoundingBoxPtr(),
        BeamInteractionDataStatePtr()->GetMutableDisColNp(), elepairptr, bspotstiff, elestiff);

    // assemble the contributions into stiffness matrix class variable
    // stiff_crosslink_ptr_, i.e. in the DOFs of the connected nodes
    BEAMINTERACTION::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrix(Discret(), elegids,
        dummyforce, elestiff, Teuchos::null, BeamInteractionDataStatePtr()->GetMutableStiff());
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::EvaluateForceStiff()
{
  CheckInitSetup();

  // force and moment exerted on the two connection sites due to the mechanical connection
  std::vector<LINALG::SerialDenseVector> bspotforce(2, LINALG::SerialDenseVector(6));

  /* linearizations, i.e. stiffness contributions due to forces on the two
   * connection sites due to the mechanical connection */
  std::vector<std::vector<LINALG::SerialDenseMatrix>> bspotstiff(
      2, std::vector<LINALG::SerialDenseMatrix>(2, LINALG::SerialDenseMatrix(6, 6)));

  // resulting discrete element force vectors of the two parent elements
  std::vector<LINALG::SerialDenseVector> eleforce(2);

  // linearizations, i.e. discrete stiffness contributions to the two parent elements
  // we can't handle this separately for both elements because there are entries which
  // couple the two element stiffness blocks
  std::vector<std::vector<LINALG::SerialDenseMatrix>> elestiff(
      2, std::vector<LINALG::SerialDenseMatrix>(2));

  // element gids of interacting elements
  std::vector<int> elegids(2);

  for (auto const& iter : doublebondcl_)
  {
    Teuchos::RCP<BEAMINTERACTION::BeamLink> elepairptr = iter.second;
    for (int i = 0; i < 2; ++i)
    {
      elegids[i] = elepairptr->GetEleGid(i);
      bspotforce[i].Zero();

      for (int j = 0; j < 2; ++j) bspotstiff[i][j].Zero();
    }

    // evaluate beam linkage object to get forces and moments on binding spots
    elepairptr->EvaluateForceStiff(bspotforce[0], bspotforce[1], bspotstiff[0][0], bspotstiff[0][1],
        bspotstiff[1][0], bspotstiff[1][1]);

    // apply forces on binding spots and corresponding linearizations to parent elements
    // and get their discrete element force vectors and stiffness matrices
    BEAMINTERACTION::UTILS::ApplyBindingSpotForceStiffToParentElements<DRT::ELEMENTS::Beam3Base,
        DRT::ELEMENTS::Beam3Base>(Discret(), PeriodicBoundingBoxPtr(),
        BeamInteractionDataStatePtr()->GetMutableDisColNp(), elepairptr, bspotforce, bspotstiff,
        eleforce, elestiff);

    // assemble the contributions into force and stiffness class variables
    // f_crosslink_np_ptr_, stiff_crosslink_ptr_, i.e. in the DOFs of the connected nodes
    BEAMINTERACTION::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrix(Discret(), elegids,
        eleforce, elestiff, BeamInteractionDataStatePtr()->GetMutableForceNp(),
        BeamInteractionDataStatePtr()->GetMutableStiff());
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UpdateStepState(const double& timefac_n)
{
  CheckInitSetup();

  // crosslinker diffusion: - according to browninan dyn for free cl
  //                        - according to beams for single and double bonded
  DiffuseCrosslinker();
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PreUpdateStepElement(bool beam_redist)
{
  CheckInitSetup();

#ifdef DEBUG
  // safety check
  if (not dis_at_last_redistr_->Map().SameAs(*BinDiscret().DofRowMap()))
    dserror(
        "current linker dof map and map of disp vector after last redistribution are\n "
        "are not the same. Something went wrong");
#endif

  linker_disnp_ = Teuchos::rcp(new Epetra_Vector(*BinDiscret().DofRowMap(), true));
  Teuchos::RCP<Epetra_Vector> dis_increment =
      Teuchos::rcp(new Epetra_Vector(*BinDiscret().DofRowMap(), true));

  LINALG::Matrix<3, 1> d;
  LINALG::Matrix<3, 1> ref;
  int doflid[3];

  for (int i = 0; i < BinDiscret().NumMyRowNodes(); ++i)
  {
    d.Clear();
    ref.Clear();

    // get a pointer at i-th row node
    DRT::Node* node = BinDiscret().lRowNode(i);

    // get GIDs of this node's degrees of freedom
    std::vector<int> dofnode = BinDiscret().Dof(node);

    for (int dim = 0; dim < 3; ++dim)
    {
      doflid[dim] = dis_at_last_redistr_->Map().LID(dofnode[dim]);
      d(dim) = (*dis_at_last_redistr_)[doflid[dim]];
      (*linker_disnp_)[doflid[dim]] = ref(dim) = node->X()[dim];
    }
    // unshift
    PeriodicBoundingBox().UnShift3D(d, ref);

    for (int dim = 0; dim < 3; ++dim) (*dis_increment)[doflid[dim]] = d(dim) - ref(dim);
  }

  // get maximal displacement increment since last redistribution over all procs
  double extrema[2] = {0.0, 0.0};
  dis_increment->MinValue(&extrema[0]);
  dis_increment->MaxValue(&extrema[1]);
  const double gmaxdisincr = std::max(-extrema[0], extrema[1]);

  // some screen output
  if (GState().GetMyRank() == 0)
    IO::cout(IO::debug) << " max linker movement " << gmaxdisincr << IO::endl;

  bool linker_redist =
      ((half_interaction_distance_ + gmaxdisincr) > (0.5 * BinStrategy().BinSizeLowerBound()));

  // store old maps prior to potential redistribution
  // this needs to be stored even no redistribution takes place later one
  StoreMapsPriorRedistribution();

  if (linker_redist or beam_redist)
  {
    // current displacement state gets new reference state
    dis_at_last_redistr_ = Teuchos::rcp(new Epetra_Vector(*linker_disnp_));
    // transfer crosslinker to new bins
    Teuchos::RCP<std::list<int>> lostcl = BeamCrosslinkerHandlerPtr()->TransferLinker(true);
    if (not lostcl->empty()) dserror("Crosslinker got lost during transfer, something went wrong");
  }

  return linker_redist;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UpdateStepElement(bool repartition_was_done)
{
  CheckInitSetup();

  if (repartition_was_done)
  {
    // adapt map of vector to map after redistribution
    BEAMINTERACTION::UTILS::UpdateDofMapOfVector(BinDiscretPtr(), dis_at_last_redistr_);

    // update double bonded linker
    UpdateMyDoubleBondsAfterRedistribution();
  }

  // gather data for all column crosslinker and column beams initially
  // this needs to be done every time since e.g. beam positions and triads
  // change every time
  UpdateAndExportCrosslinkerData();
  UpdateAndExportBeamData();

  // manage binding and unbinding events of crosslinker
  BindAndUnbindCrosslinker();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PostUpdateStepElement()
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::map<STR::EnergyType, double> BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::GetEnergy()
    const
{
  CheckInitSetup();
  std::map<STR::EnergyType, double> cl_energies;

  for (auto db_iter : doublebondcl_)
  {
    cl_energies[STR::beam_to_beam_link_internal_energy] += db_iter.second->GetInternalEnergy();
    cl_energies[STR::beam_to_beam_link_kinetic_energy] += db_iter.second->GetKineticEnergy();
  }

  return cl_energies;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::OutputStepState(
    IO::DiscretizationWriter& iowriter) const
{
  CheckInitSetup();
  // not used, we are writing output during runtime
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::RuntimeOutputStepState() const
{
  CheckInitSetup();

  TEUCHOS_FUNC_TIME_MONITOR(
      "BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::"
      "RuntimeOutputStepState");

  if (vtp_writer_ptr_ != Teuchos::null) WriteOutputRuntimeVtpStructure();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::InitOutputRuntimeVtpStructure()
{
  CheckInit();

  vtp_writer_ptr_ = Teuchos::rcp(new DiscretizationRuntimeVtpWriter());

  // Todo: we need a better upper bound for total number of time steps here
  // however, this 'only' affects the number of leading zeros in the vtk file names
  const unsigned int num_timesteps_in_simulation_upper_bound = 1000000;

  // initialize the writer object
  vtp_writer_ptr_->Initialize(BinDiscretPtr(), BinDiscretPtr()->Name(),
      num_timesteps_in_simulation_upper_bound, GState().GetTimeN(),
      GInOutput().GetRuntimeVtpOutputParams()->WriteBinaryOutput());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::WriteOutputRuntimeVtpStructure() const
{
  CheckInitSetup();

  // ************** BEGIN RUNTIME VTP OUTPUT *** OPTION 1: DISCRETIZATION *********
  // this section compiles and seems to do the job correctly :-)

  // initialize the writer object
  vtp_writer_ptr_->SetGeometryFromParticleDiscretization();

  // reset time and time step of the writer object
  vtp_writer_ptr_->ResetTimeAndTimeStep(GState().GetTimeN(), GState().GetStepN());

  // append all desired node and dof output data to the writer object's storage
  DRT::Discretization const& bindis = BinDiscret();
  // node
  Teuchos::RCP<Epetra_Vector> numbond = LINALG::CreateVector(*bindis.NodeRowMap(), true);
  Teuchos::RCP<Epetra_Vector> owner = LINALG::CreateVector(*bindis.NodeRowMap(), true);

  // dof
  Teuchos::RCP<Epetra_Vector> dis = LINALG::CreateVector(*bindis.DofRowMap(), true);
  Teuchos::RCP<Epetra_Vector> orientation = LINALG::CreateVector(*bindis.DofRowMap(), true);
  Teuchos::RCP<Epetra_Vector> force = LINALG::CreateVector(*bindis.DofRowMap(), true);

  FillStateDataVectorsForOutput(dis, orientation, numbond, owner, force);

  // append displacement vector if desired
  // append displacement if desired
  //   if ( GInOutput().GetRuntimeVtpOutputParams()->OutputDisplacementState() )
  //     vtp_writer_ptr_->AppendDofBasedResultDataVector( dis, 3, "displacement" );

  // append owner if desired
  if (GInOutput().GetRuntimeVtpOutputParams()->OutputOwner())
    vtp_writer_ptr_->AppendNodeBasedResultDataVector(owner, 1, "owner");

  // append orientation vector if desired
  if (GInOutput().GetRuntimeVtpOutputParams()->OutputOrientationAndLength())
    vtp_writer_ptr_->AppendDofBasedResultDataVector(orientation, 3, "orientation");

  // append number of bonds if desired
  if (GInOutput().GetRuntimeVtpOutputParams()->OutputNumberOfBonds())
    vtp_writer_ptr_->AppendNodeBasedResultDataVector(numbond, 1, "numberofbonds");

  // append number of bonds if desired
  if (GInOutput().GetRuntimeVtpOutputParams()->OutputLinkingForce())
    vtp_writer_ptr_->AppendDofBasedResultDataVector(force, 3, "force");

  // finalize everything and write all required VTU files to file system
  vtp_writer_ptr_->WriteFiles();

  // write collection files
  vtp_writer_ptr_->WriteCollectionFileOfAllWrittenFiles(BinDiscretPtr()->Name());



  // ************** BEGIN RUNTIME VTP OUTPUT *** OPTION 2: DIRECTLY *********
  // this section is just to get the idea and needs some minor modifications (indicated by Fixme)
  // this may serve as a template for visualization of stochastic&viscous forces, contact forces,
  // ...

  //  Teuchos::RCP<RuntimeVtpWriter> vtp_writer_ptr =
  //      Teuchos::rcp( new RuntimeVtpWriter() );
  //
  //  // Todo: we need a better upper bound for total number of time steps here
  //  // however, this 'only' affects the number of leading zeros in the vtk file names
  //  const unsigned int num_timesteps_in_simulation_upper_bound = 1000000;
  //
  //
  //  // initialize the writer object
  //  vtp_writer_ptr->Initialize();   // Fixme
  //
  //  // set geometry manually
  //  const unsigned int num_spatial_dimensions = 3;
  //  unsigned int num_row_points = 2000;
  //
  //  // get and prepare storage for point coordinate values
  //  std::vector< double >& point_coordinates = vtp_writer_ptr->GetMutablePointCoordinateVector();
  //  point_coordinates.clear();
  //  point_coordinates.reserve( num_spatial_dimensions * num_row_points );
  //
  //
  //  // loop over my points and collect the geometry/grid data, i.e. reference positions of nodes
  //  for (unsigned int inode=0; inode < num_row_points; ++inode)
  //  {
  //    const DRT::Node* node = BinDiscretPtr()->lRowNode(inode);
  //
  //    for (unsigned int idim=0; idim<num_spatial_dimensions; ++idim)
  //      point_coordinates.push_back( node->X()[idim] );
  //  }
  //
  //
  //  // reset time and time step and geometry name in the writer object
  //  vtp_writer_ptr->SetupForNewTimeStepAndGeometry(
  //      GState().GetTimeN(), GState().GetStepN(), BinDiscretPtr()->Name() );
  //
  //
  //
  //  // append all desired output data to the writer object's storage
  //
  //  // number of bonds: collect data and append to visualization results if desired
  //  std::vector< double > num_bonds( num_row_points );
  //
  //  for ( unsigned int i = 0; i < num_row_points; ++i )
  //  {
  //    CROSSLINKING::CrosslinkerNode *crosslinker_i =
  //        dynamic_cast< CROSSLINKING::CrosslinkerNode* >( bindis.lRowNode(i) );
  //
  //    num_bonds[i] = crosslinker_i->ClData()->GetNumberOfBonds();
  //  }
  //
  //  vtp_writer_ptr->AppendVisualizationPointDataVector( num_bonds, 1, "num_bonds" );
  //
  //
  //
  //  // finalize everything and write all required VTU files to filesystem
  //  vtp_writer_ptr->WriteFiles();
  //
  //
  //  // write a collection file summarizing all previously written files
  //  vtp_writer_ptr->WriteCollectionFileOfAllWrittenFiles();   // Fixme
  // ************** END RUNTIME VTP OUTPUT ***************************************
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::FillStateDataVectorsForOutput(
    Teuchos::RCP<Epetra_Vector> displacement, Teuchos::RCP<Epetra_Vector> orientation,
    Teuchos::RCP<Epetra_Vector> numberofbonds, Teuchos::RCP<Epetra_Vector> owner,
    Teuchos::RCP<Epetra_Vector> force) const
{
  CheckInitSetup();
  DRT::Discretization const& bindis = BinDiscret();

  LINALG::SerialDenseVector bspotforce(true);

  // todo: this is of course not nice, this needs to be done somewhere else
  for (int i = 0; i < bindis.NumMyRowNodes(); ++i)
  {
    DRT::Node* crosslinker_i = bindis.lRowNode(i);
    // std::vector holding gids of dofs
    std::vector<int> dofnode = bindis.Dof(crosslinker_i);
    int numbonds = crosslinker_data_[crosslinker_i->LID()]->GetNumberOfBonds();

    Teuchos::RCP<BEAMINTERACTION::BeamLink> beamlink = Teuchos::null;

    if (numbonds == 2)
    {
      beamlink = doublebondcl_.at(crosslinker_i->Id());
      beamlink->GetBindingSpotForce(0, bspotforce);
    }

    // loop over all dofs
    for (unsigned int dim = 0; dim < 3; ++dim)
    {
      int doflid = displacement->Map().LID(dofnode[dim]);
      (*displacement)[doflid] = crosslinker_i->X()[dim];

      if (numbonds == 2)
      {
        (*orientation)[doflid] =
            beamlink->GetBindSpotPos1()(dim) - beamlink->GetBindSpotPos2()(dim);
        (*force)[doflid] = bspotforce(dim);
      }
      else
      {
        (*orientation)[doflid] = 0.0;
        (*force)[doflid] = 0.0;
      }
    }

    (*numberofbonds)[i] = numbonds;
    (*owner)[i] = crosslinker_i->Owner();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ResetStepState() { CheckInitSetup(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::WriteRestart(
    IO::DiscretizationWriter& ia_writer, IO::DiscretizationWriter& bin_writer) const
{
  CheckInitSetup();

  // -------------------------------------------------------------------------
  // 1) write list of double bonded crosslinker on each proc
  // -------------------------------------------------------------------------
  DRT::PackBuffer linker_buffer;
  for (auto const& iter : doublebondcl_)
  {
    Teuchos::RCP<BEAMINTERACTION::BeamLink> btbl = iter.second;
    btbl->Pack(linker_buffer);
  }
  linker_buffer.StartPacking();
  for (auto const& iter : doublebondcl_)
  {
    Teuchos::RCP<BEAMINTERACTION::BeamLink> btbl = iter.second;
    btbl->Pack(linker_buffer);
  }

  Teuchos::RCP<std::vector<char>> db_linker = Teuchos::rcp(new std::vector<char>);
  std::swap(*db_linker, linker_buffer());

  // -------------------------------------------------------------------------
  // 2) write crosslinker data
  // -------------------------------------------------------------------------
  DRT::PackBuffer cldata_buffer;
  unsigned int numrowcl = BinDiscret().NumMyRowNodes();
  for (unsigned int i = 0; i < numrowcl; ++i)
  {
    int const clgid = BinDiscret().NodeRowMap()->GID(i);
    Teuchos::RCP<BEAMINTERACTION::DATA::CrosslinkerData> cl_data_i =
        crosslinker_data_[BinDiscret().NodeColMap()->LID(clgid)];

    cl_data_i->Pack(cldata_buffer);
  }
  cldata_buffer.StartPacking();
  for (unsigned int i = 0; i < numrowcl; ++i)
  {
    int const clgid = BinDiscret().NodeRowMap()->GID(i);
    Teuchos::RCP<BEAMINTERACTION::DATA::CrosslinkerData> cl_data_i =
        crosslinker_data_[BinDiscret().NodeColMap()->LID(clgid)];

    cl_data_i->Pack(cldata_buffer);
  }

  Teuchos::RCP<std::vector<char>> cldata = Teuchos::rcp(new std::vector<char>);
  std::swap(*cldata, cldata_buffer());

  // -------------------------------------------------------------------------
  // 3) beam data
  // -------------------------------------------------------------------------
  DRT::PackBuffer beamdata_buffer;
  unsigned int numrowbeam = EleTypeMapExtractor().BeamMap()->NumMyElements();
  for (unsigned int i = 0; i < numrowbeam; ++i)
  {
    int const beamgid = EleTypeMapExtractor().BeamMap()->GID(i);
    Teuchos::RCP<BEAMINTERACTION::DATA::BeamData> beam_data_i =
        beam_data_[Discret().ElementColMap()->LID(beamgid)];

#ifdef DEBUG
    if (beam_data_i == Teuchos::null)
      dserror("beam data of row beam with gid %i not there", beamgid);
#endif

    beam_data_i->Pack(beamdata_buffer);
  }
  beamdata_buffer.StartPacking();
  for (unsigned int i = 0; i < numrowbeam; ++i)
  {
    int const beamgid = EleTypeMapExtractor().BeamMap()->GID(i);
    Teuchos::RCP<BEAMINTERACTION::DATA::BeamData> beam_data_i =
        beam_data_[Discret().ElementColMap()->LID(beamgid)];

    beam_data_i->Pack(beamdata_buffer);
  }

  Teuchos::RCP<std::vector<char>> beamdata = Teuchos::rcp(new std::vector<char>);
  std::swap(*beamdata, beamdata_buffer());

  // -------------------------------------------------------------------------
  // write data
  // -------------------------------------------------------------------------
  bin_writer.WriteCharVector("Linker", db_linker);
  bin_writer.WriteCharVector("ClData", cldata);
  bin_writer.WriteCharVector("BeamData", beamdata);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PreReadRestart()
{
  CheckInitSetup();
  BeamCrosslinkerHandlerPtr()->RemoveAllLinker();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ReadRestart(
    IO::DiscretizationReader& ia_reader, IO::DiscretizationReader& bin_reader)
{
  CheckInitSetup();

  // -------------------------------------------------------------------------
  // 1) read list of double bonded crosslinker on each proc
  // -------------------------------------------------------------------------
  Teuchos::RCP<std::vector<char>> linkercharvec;
  bin_reader.ReadCharVector(linkercharvec, "Linker");

  std::vector<char>::size_type index = 0;
  while (index < linkercharvec->size())
  {
    std::vector<char> data;
    DRT::ParObject::ExtractfromPack(index, *linkercharvec, data);
    Teuchos::RCP<DRT::ParObject> object = Teuchos::rcp(DRT::UTILS::Factory(data), true);
    Teuchos::RCP<BEAMINTERACTION::BeamLink> beamtobeamlink =
        Teuchos::rcp_dynamic_cast<BEAMINTERACTION::BeamLink>(object);
    if (beamtobeamlink == Teuchos::null) dserror("Failed to build a node from the node data");

    // insert in my list of double bonded crosslinker
    doublebondcl_[beamtobeamlink->Id()] = beamtobeamlink;
  }

  // -------------------------------------------------------------------------
  // 2) read crosslinker data
  // -------------------------------------------------------------------------
  Teuchos::RCP<std::vector<char>> cldata_charvec;
  bin_reader.ReadCharVector(cldata_charvec, "ClData");
  crosslinker_data_.resize(BinDiscret().NumMyColNodes());

  std::map<int, Teuchos::RCP<BEAMINTERACTION::DATA::CrosslinkerData>> cl_not_owned;
  std::set<int> not_owned_gids;
  index = 0;
  std::map<int, std::vector<char>> cl_datapacks;
  std::vector<int> read_node_ids;
  while (index < cldata_charvec->size())
  {
    // unpack
    std::vector<char> recv_singlecontainer_data;
    DRT::ParObject::ExtractfromPack(index, *cldata_charvec, recv_singlecontainer_data);

    Teuchos::RCP<BEAMINTERACTION::DATA::CrosslinkerData> cl_data = Teuchos::rcp(
        BEAMINTERACTION::DATA::CreateDataContainer<BEAMINTERACTION::DATA::CrosslinkerData>(
            recv_singlecontainer_data),
        true);

    int const cl_gid = cl_data->GetId();
    read_node_ids.push_back(cl_gid);

    // repack it for communication to find owner
    DRT::PackBuffer data;
    cl_data->Pack(data);
    data.StartPacking();
    cl_data->Pack(data);
    cl_datapacks[cl_gid].insert(cl_datapacks[cl_gid].end(), data().begin(), data().end());
  }

  // build dummy map according to read data on myrank
  Teuchos::RCP<Epetra_Map> dummy_cl_map = Teuchos::rcp(
      new Epetra_Map(-1, read_node_ids.size(), &read_node_ids[0], 0, BinDiscret().Comm()));

  // build exporter object
  Teuchos::RCP<DRT::Exporter> exporter = Teuchos::rcp(
      new DRT::Exporter(*dummy_cl_map, *BinDiscret().NodeColMap(), BinDiscret().Comm()));

  // export
  exporter->Export(cl_datapacks);

  // rebuild data container
  crosslinker_data_.resize(BinDiscret().NumMyColNodes());
  for (auto& iter : cl_datapacks)
  {
    // this needs to done
    // Fixme
    std::vector<char>::size_type position = 0;
    std::vector<char> data;
    DRT::ParObject::ExtractfromPack(position, iter.second, data);

    Teuchos::RCP<BEAMINTERACTION::DATA::CrosslinkerData> cl_data = Teuchos::rcp(
        BEAMINTERACTION::DATA::CreateDataContainer<BEAMINTERACTION::DATA::CrosslinkerData>(data),
        true);
    crosslinker_data_[BinDiscret().NodeColMap()->LID(cl_data->GetId())] = cl_data;
  }

  // -------------------------------------------------------------------------
  // 3) read beam data
  // -------------------------------------------------------------------------
  Teuchos::RCP<std::vector<char>> beamdata_charvec;
  bin_reader.ReadCharVector(beamdata_charvec, "BeamData");
  beam_data_.resize(Discret().NumMyColElements());

  std::map<int, Teuchos::RCP<BEAMINTERACTION::DATA::BeamData>> beams_not_owned;
  not_owned_gids.clear();
  index = 0;
  std::map<int, std::vector<char>> beam_datapacks;
  std::vector<int> read_ele_ids;
  while (index < beamdata_charvec->size())
  {
    // unpack
    std::vector<char> recv_singlecontainer_data;
    DRT::ParObject::ExtractfromPack(index, *beamdata_charvec, recv_singlecontainer_data);

    Teuchos::RCP<BEAMINTERACTION::DATA::BeamData> beam_data =
        Teuchos::rcp(BEAMINTERACTION::DATA::CreateDataContainer<BEAMINTERACTION::DATA::BeamData>(
                         recv_singlecontainer_data),
            true);

    int const beam_gid = beam_data->GetId();
    read_ele_ids.push_back(beam_gid);

    // repack it for communication to find owner
    DRT::PackBuffer data;
    beam_data->Pack(data);
    data.StartPacking();
    beam_data->Pack(data);
    beam_datapacks[beam_gid].insert(beam_datapacks[beam_gid].end(), data().begin(), data().end());
  }

  // build dummy map according to read data on myrank
  Teuchos::RCP<Epetra_Map> dummy_beam_map =
      Teuchos::rcp(new Epetra_Map(-1, read_ele_ids.size(), &read_ele_ids[0], 0, Discret().Comm()));

  // build exporter object
  exporter = Teuchos::rcp(
      new DRT::Exporter(*dummy_beam_map, *Discret().ElementColMap(), Discret().Comm()));

  // export
  exporter->Export(beam_datapacks);

  // rebuild data container
  beam_data_.resize(Discret().NumMyColElements());
  for (auto& iter : beam_datapacks)
  {
    // this needs to be done
    // Fixme
    std::vector<char>::size_type position = 0;
    std::vector<char> data;
    DRT::ParObject::ExtractfromPack(position, iter.second, data);

    Teuchos::RCP<BEAMINTERACTION::DATA::BeamData> beam_data = Teuchos::rcp(
        BEAMINTERACTION::DATA::CreateDataContainer<BEAMINTERACTION::DATA::BeamData>(data), true);
    beam_data_[Discret().ElementColMap()->LID(beam_data->GetId())] = beam_data;
  }

  // init maps
  StoreMapsPriorRedistribution();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PostReadRestart()
{
  CheckInitSetup();

  // bring each object in doublebondcl_ map to its correct owner
  UpdateMyDoubleBondsRemoteIdList();

  // store displacement of restart step as displacement state of last redistribution
  dis_at_last_redistr_ = Teuchos::rcp(new Epetra_Vector(*BinDiscret().DofRowMap(), true));
  for (int i = 0; i < BinDiscret().NumMyRowNodes(); ++i)
  {
    DRT::Node* crosslinker_i = BinDiscret().lRowNode(i);

    // std::vector holding gids of dofs
    std::vector<int> dofnode = BinDiscret().Dof(crosslinker_i);

    // loop over all dofs
    for (unsigned int dim = 0; dim < 3; ++dim)
    {
      int doflid = dis_at_last_redistr_->Map().LID(dofnode[dim]);
      (*dis_at_last_redistr_)[doflid] = crosslinker_i->X()[dim];
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::AddBinsToBinColMap(std::set<int>& colbins)
{
  // nothing to do
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::AddBinsWithRelevantContentForIaDiscretColMap(
    std::set<int>& colbins) const
{
  CheckInitSetup();

  bool map_changed = false;
  if (not(cl_noderowmap_prior_redistr_->SameAs(*BinDiscret().NodeRowMap()) and
          cl_nodecolmap_prior_redistr_->SameAs(*BinDiscret().NodeColMap())))
    map_changed = true;

  std::set<int> colbinsext;
  for (auto const& iter : colbins)
  {
    // get current bin
    DRT::Element* currbin = BinDiscret().gElement(iter);
    // get all crosslinker in current bin
    DRT::Node** clincurrentbin = currbin->Nodes();

    // loop over all crosslinker in CurrentBin
    for (int i = 0; i < currbin->NumNode(); ++i)
    {
      if (map_changed or (crosslinker_data_.size() == 0) or
          (crosslinker_data_[clincurrentbin[i]->LID()]->GetNumberOfBonds() > 0))
      {
        std::vector<int> binvec;
        BinStrategy().GetNeighborBinIds(iter, binvec);
        colbinsext.insert(binvec.begin(), binvec.end());
        // go to next bin
        break;
      }
    }
  }

  colbins.insert(colbinsext.begin(), colbinsext.end());

  // the following would be a default additional one ghost layer

  //
  //  std::set< int > colbinsext(colbins.begin(),colbins.end());
  //  std::set< int >::const_iterator biniter;
  //  for ( biniter = colbinsext.begin(); biniter != colbinsext.end() ; ++biniter )
  //  {
  //    std::vector< int > binvec;
  //    BinStrategy().GetNeighborAndOwnBinIds( *biniter, binvec );
  //    colbinsext.insert( binvec.begin(), binvec.end() );
  //  }
  //
  //  colbins.insert( colbinsext.begin(), colbinsext.end() );
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::GetHalfInteractionDistance(
    double& half_interaction_distance)
{
  CheckInitSetup();

  // loop over all linker of all linker types and get largest linker (also
  // considering tolerance)
  double curr_ia_dist = -1.0;
  double local_half_interaction_distance = 0.0;
  const int numrowcl = BinDiscretPtr()->NumMyRowNodes();
  for (int rowcli = 0; rowcli < numrowcl; ++rowcli)
  {
    // get current linker
    CROSSLINKING::CrosslinkerNode* crosslinker_i =
        dynamic_cast<CROSSLINKING::CrosslinkerNode*>(BinDiscretPtr()->lRowNode(rowcli));

    curr_ia_dist = 0.5 * (crosslinker_i->GetMaterial()->LinkingLength() +
                             crosslinker_i->GetMaterial()->LinkingLengthTolerance());

    local_half_interaction_distance = (curr_ia_dist > local_half_interaction_distance)
                                          ? curr_ia_dist
                                          : local_half_interaction_distance;
  }

  // get global maximum
  double global_half_interaction_distance = 0.0;
  // build sum over all procs
  MPI_Allreduce(&local_half_interaction_distance, &global_half_interaction_distance, 1, MPI_DOUBLE,
      MPI_MAX, dynamic_cast<const Epetra_MpiComm*>(&(Discret().Comm()))->Comm());
  half_interaction_distance_ = global_half_interaction_distance;

  // some screen output
  if (GState().GetMyRank() == 0)
    IO::cout(IO::verbose) << " beam to beam crosslinking half interaction distance "
                          << global_half_interaction_distance << IO::endl;

  half_interaction_distance = (half_interaction_distance_ > half_interaction_distance)
                                  ? half_interaction_distance_
                                  : half_interaction_distance;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::DiffuseCrosslinker()
{
  CheckInitSetup();

  // loop over all row crosslinker (beam binding status not touched here)
  const int numrowcl = BinDiscretPtr()->NumMyRowNodes();
  for (int rowcli = 0; rowcli < numrowcl; ++rowcli)
  {
    // get current linker
    CROSSLINKING::CrosslinkerNode* crosslinker_i =
        dynamic_cast<CROSSLINKING::CrosslinkerNode*>(BinDiscretPtr()->lRowNode(rowcli));

#ifdef DEBUG
    if (crosslinker_i->NumElement() != 1) dserror("More than one element for this crosslinker");
#endif

    BEAMINTERACTION::DATA::CrosslinkerData* cldata_i =
        crosslinker_data_[crosslinker_i->LID()].get();

    // different treatment according to number of bonds a crosslinker has
    switch (cldata_i->GetNumberOfBonds())
    {
      case 0:
      {
        // crosslinker has zero bonds, i.e. is free to diffuse according to
        // brownian dynamics
        DiffuseUnboundCrosslinker(crosslinker_i, cldata_i);
        break;
      }
      case 1:
      {
        // get clbspot that is currently bonded
        int occbspotid = GetSingleOccupiedClBspot(cldata_i->GetBSpots());

        // get current position of binding spot of filament partner
        // note: we can not use our beam data container, as bspot position is not current position
        // (as this is the result of a sum, you can not have a reference to that)
        const int elegid = cldata_i->GetBSpots()[occbspotid].first;

        DRT::ELEMENTS::Beam3Base* ele =
            dynamic_cast<DRT::ELEMENTS::Beam3Base*>(DiscretPtr()->gElement(elegid));

#ifdef DEBUG
        // safety check
        const int colelelid = DiscretPtr()->ElementColMap()->LID(elegid);
        if (colelelid < 0)
          dserror(
              "Crosslinker has %i bonds but his binding partner with gid %i "
              "is \nnot ghosted/owned on proc %i (owner of crosslinker)",
              cldata_i->GetNumberOfBonds(), elegid, GState().GetMyRank());
        // safety check
        if (ele == NULL)
          dserror("Dynamic cast of ele with gid %i failed on proc ", elegid, GState().GetMyRank());
#endif

        // get current position of filament binding spot
        LINALG::Matrix<3, 1> bbspotpos;
        std::vector<double> eledisp;
        BEAMINTERACTION::UTILS::GetCurrentUnshiftedElementDis(Discret(), ele,
            BeamInteractionDataStatePtr()->GetMutableDisColNp(), PeriodicBoundingBox(), eledisp);
        ele->GetPosOfBindingSpot(bbspotpos, eledisp, crosslinker_i->GetMaterial()->LinkerType(),
            cldata_i->GetBSpots()[occbspotid].second, PeriodicBoundingBox());

        // note: a crosslinker can not leave the computational domain here, as no beam binding
        // spot can be outside the periodic box at this point
        std::vector<double> newpos(3, 0.0);
        for (int dim = 0; dim < 3; ++dim) newpos[dim] = bbspotpos(dim);

        crosslinker_i->SetPos(newpos);
        cldata_i->SetPosition(bbspotpos);

        break;
      }
      case 2:
      {
        // crosslinker has two bonds (cl gets current mid position between the filament
        // binding spot it is attached to)
        // -----------------------------------------------------------------
        // partner one
        // -----------------------------------------------------------------
        int elegid = cldata_i->GetBSpots()[0].first;

#ifdef DEBUG
        if (elegid < 0 or cldata_i->GetBSpots()[0].second < 0)
          dserror(
              " double bonded crosslinker has stored beam partner gid or loc bsponum of -1, "
              " something went wrong");
        // safety check
        int colelelid = DiscretPtr()->ElementColMap()->LID(elegid);
        if (colelelid < 0)
          dserror(
              "Crosslinker has %i bonds but his binding partner with gid %i "
              "is not \nghosted/owned on proc %i (owner of crosslinker)",
              cldata_i->GetNumberOfBonds(), elegid, GState().GetMyRank());
#endif

        DRT::ELEMENTS::Beam3Base* ele =
            dynamic_cast<DRT::ELEMENTS::Beam3Base*>(DiscretPtr()->gElement(elegid));

#ifdef DEBUG
        // safety check
        if (ele == NULL)
          dserror("Dynamic cast of ele with gid %i failed on proc ", elegid, GState().GetMyRank());
#endif

        // get current position of filament binding spot
        LINALG::Matrix<3, 1> bbspotposone;
        std::vector<double> eledisp;
        BEAMINTERACTION::UTILS::GetCurrentUnshiftedElementDis(Discret(), ele,
            BeamInteractionDataStatePtr()->GetMutableDisColNp(), PeriodicBoundingBox(), eledisp);
        ele->GetPosOfBindingSpot(bbspotposone, eledisp, crosslinker_i->GetMaterial()->LinkerType(),
            cldata_i->GetBSpots()[0].second, PeriodicBoundingBox());

        // -----------------------------------------------------------------
        // partner two
        // -----------------------------------------------------------------
        elegid = cldata_i->GetBSpots()[1].first;

#ifdef DEBUG
        // safety check
        if (elegid < 0 or cldata_i->GetBSpots()[1].second < 0)
          dserror(
              " double bonded crosslinker has stored beam partner gid or loc bsponum of -1, "
              " something went wrong");
        colelelid = DiscretPtr()->ElementColMap()->LID(elegid);
        if (colelelid < 0)
          dserror(
              "Crosslinker has %i bonds but his binding partner with gid %i "
              "is \nnot ghosted/owned on proc %i (owner of crosslinker)",
              cldata_i->GetNumberOfBonds(), elegid, GState().GetMyRank());
#endif

        ele = dynamic_cast<DRT::ELEMENTS::Beam3Base*>(DiscretPtr()->gElement(elegid));

#ifdef DEBUG
        // safety check
        if (ele == NULL)
          dserror("Dynamic cast of ele with gid %i failed on proc ", elegid, GState().GetMyRank());
#endif

        // get current position of filament binding spot
        LINALG::Matrix<3, 1> bbspotpostwo;
        BEAMINTERACTION::UTILS::GetCurrentUnshiftedElementDis(Discret(), ele,
            BeamInteractionDataStatePtr()->GetMutableDisColNp(), PeriodicBoundingBox(), eledisp);
        ele->GetPosOfBindingSpot(bbspotpostwo, eledisp, crosslinker_i->GetMaterial()->LinkerType(),
            cldata_i->GetBSpots()[1].second, PeriodicBoundingBox());

        LINALG::Matrix<3, 1> clpos(true);
        SetPositionOfDoubleBondedCrosslinkerPBCconsistent(clpos, bbspotposone, bbspotpostwo);

        std::vector<double> newpos(3, 0.0);
        for (int dim = 0; dim < 3; ++dim) newpos[dim] = clpos(dim);

        crosslinker_i->SetPos(newpos);
        cldata_i->SetPosition(clpos);

        break;
      }
      default:
      {
        dserror("Unrealistic number %i of bonds for a crosslinker.", cldata_i->GetNumberOfBonds());
        exit(EXIT_FAILURE);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::DiffuseUnboundCrosslinker(
    DRT::Node* crosslinker_i, BEAMINTERACTION::DATA::CrosslinkerData* cldata_i)
{
  CheckInit();

  CROSSLINKING::CrosslinkerNode* crosslinker =
      dynamic_cast<CROSSLINKING::CrosslinkerNode*>(crosslinker_i);

#ifdef DEBUG
  if (crosslinker == NULL) dserror("Dynamic cast to CrosslinkerNode failed");
#endif

  // get standard deviation and mean value for crosslinker that are free to
  // diffuse
  double standarddev = std::sqrt(2.0 * crosslinking_params_ptr_->KT() /
                                 (3.0 * M_PI * crosslinking_params_ptr_->Viscosity() *
                                     crosslinker->GetMaterial()->LinkingLength()) *
                                 crosslinking_params_ptr_->DeltaTime());
  double meanvalue = 0.0;
  // Set mean value and standard deviation of normal distribution
  // FixMe standard deviation = sqrt(variance) check this for potential error !!!
  DRT::Problem::Instance()->Random()->SetMeanVariance(meanvalue, standarddev);

  // diffuse crosslinker according to brownian dynamics
  LINALG::Matrix<3, 1> newclpos(true);
  std::vector<double> randvec;
  int count = 3;
  // maximal diffusion given by cutoff radius (sqrt(3) = 1.73..)
  double const maxmov = BinStrategy().BinSizeLowerBound() / 1.74;
  DRT::Problem::Instance()->Random()->Normal(randvec, count);
  for (int dim = 0; dim < 3; ++dim)
  {
    if (abs(randvec[dim]) > maxmov)
    {
      double old = randvec[dim];
      randvec[dim] = (abs(randvec[dim]) / randvec[dim]) * maxmov;
      IO::cout(IO::verbose) << "Movement of free crosslinker " << crosslinker->Id()
                            << " was restricted by cutoff radius"
                               " in "
                            << dim << " direction. " << old << " to " << randvec[dim]
                            << "\nThis should not happen to often "
                               "to stay physical. Increase cutoff or reduce movement"
                            << IO::endl;
    }
    newclpos(dim) = crosslinker->X()[dim] + randvec[dim];
  }

  // check compliance with periodic boundary conditions
  PeriodicBoundingBox().Shift3D(newclpos);
  std::vector<double> newpos(3, 0.0);
  for (int dim = 0; dim < 3; ++dim) newpos[dim] = newclpos(dim);
  crosslinker->SetPos(newpos);

  cldata_i->SetPosition(newclpos);
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
int BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::GetSingleOccupiedClBspot(
    std::vector<std::pair<int, int>> const& clbspots) const
{
  CheckInit();

  if (clbspots[0].first > -1)
    return 0;
  else if (clbspots[1].first > -1)
    return 1;
  else
    dserror("numbond = 1 but both binding spots store invalid element GIDs!");

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::
    SetPositionOfDoubleBondedCrosslinkerPBCconsistent(LINALG::Matrix<3, 1>& clpos,
        LINALG::Matrix<3, 1> const& bspot1pos, LINALG::Matrix<3, 1> const& bspot2pos) const
{
  /* the position of (the center) of a double-bonded crosslinker is defined as
   * midpoint between the two given binding spot positions. (imagine a linker
   * being a slender body with a binding domain at each of both ends) */

  /* if the two binding spots are separated by a periodic boundary, we need to
   * shift one position back to get the interpolation right */
  clpos = bspot2pos;
  PeriodicBoundingBox().UnShift3D(clpos, bspot1pos);

  // fixme: to avoid senseless dserror in debug mode
  LINALG::Matrix<3, 1> dummy(clpos);
  clpos.Update(0.5, bspot1pos, 0.5, dummy);

  // shift the interpolated position back in the periodic box if necessary
  PeriodicBoundingBox().Shift3D(clpos);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::SetPositionOfNewlyFreeCrosslinker(
    CROSSLINKING::CrosslinkerNode* crosslinker, BEAMINTERACTION::DATA::CrosslinkerData* cldata)
{
  CheckInit();

  // generate vector in random direction
  // of length half the linking length to "reset" crosslink molecule position: it may now
  // reenter or leave the bonding proximity
  // todo: does this make sense?
  LINALG::Matrix<3, 1> clpos(cldata->GetPosition());
  LINALG::Matrix<3, 1> cldeltapos_i;
  std::vector<double> randunivec(3);
  int count = 3;
  DRT::Problem::Instance()->Random()->Uni(randunivec, count);
  for (unsigned int dim = 0; dim < 3; ++dim) cldeltapos_i(dim) = randunivec[dim];

  cldeltapos_i.Scale(crosslinker->GetMaterial()->LinkingLength() / cldeltapos_i.Norm2());

  clpos.Update(1.0, cldeltapos_i, 1.0);

  PeriodicBoundingBox().Shift3D(clpos);

  std::vector<double> newpos(3, 0.0);
  for (int dim = 0; dim < 3; ++dim) newpos[dim] = clpos(dim);
  crosslinker->SetPos(newpos);

  cldata->SetPosition(clpos);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::SetPositionOfNewlySingleBondedCrosslinker(
    CROSSLINKING::CrosslinkerNode* crosslinker, BEAMINTERACTION::DATA::CrosslinkerData* cldata,
    int stayoccpotid)
{
  CheckInit();

  // update postion
  const int collidoccbeam =
      DiscretPtr()->ElementColMap()->LID(cldata->GetBSpots()[stayoccpotid].first);

#ifdef DEBUG
  // safety check
  if (collidoccbeam < 0)
    dserror("element with gid %i not ghosted on proc %i", cldata->GetBSpots()[stayoccpotid].first,
        GState().GetMyRank());
#endif

  BEAMINTERACTION::DATA::BeamData const* beamdata_i = beam_data_[collidoccbeam].get();
  LINALG::Matrix<3, 1> clpos(beamdata_i->GetBSpotPosition(
      crosslinker->GetMaterial()->LinkerType(), cldata->GetBSpots()[stayoccpotid].second));

  std::vector<double> newpos(3, 0.0);
  for (int dim = 0; dim < 3; ++dim) newpos[dim] = clpos(dim);
  crosslinker->SetPos(newpos);

  cldata->SetPosition(clpos);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::StoreMapsPriorRedistribution()
{
  CheckInit();

  *cl_noderowmap_prior_redistr_ = *BinDiscret().NodeRowMap();
  *cl_nodecolmap_prior_redistr_ = *BinDiscret().NodeColMap();
  *beam_elerowmap_prior_redistr_ = *EleTypeMapExtractor().BeamMap();
  *beam_elecolmap_prior_redistr_ = *Discret().ElementColMap();
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UpdateAndExportCrosslinkerData()
{
  CheckInit();

  TEUCHOS_FUNC_TIME_MONITOR(
      "BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::"
      "UpdateAndExportCrosslinkerData");

  //  // if one proc has changed maps, all procs have to update exporter
  //  int loc_changed_maps = static_cast< int >( (cl_exporter_ == Teuchos::null) or not
  //      ( cl_exporter_->SourceMap().SameAs( *cl_noderowmap_prior_redistr_ ) and
  //        cl_nodecolmap_prior_redistr_->SameAs(*BinDiscret().NodeColMap() ) ) );
  //
  //  //global filled flag (is true / one if and only if loc_changed_maps == true on each processor
  //  int g_changed_maps = 0;
  //
  //  /*the global flag is set to the maximal value of any local flag
  //   * i.e. if on any processor loc_changed_maps == true, the flag g_changed_maps is set to
  //   * one*/
  //  BinDiscret().Comm().MaxAll( &loc_changed_maps, &g_changed_maps, 1 );
  //
  //  if ( g_changed_maps )
  cl_exporter_ = Teuchos::rcp(new DRT::Exporter(
      *cl_noderowmap_prior_redistr_, *BinDiscret().NodeColMap(), BinDiscret().Comm()));

  // we first need to pack our stuff into and std::vector< char > for communication
  std::map<int, std::vector<char>> allpacks;
  unsigned int numrowcl = cl_noderowmap_prior_redistr_->NumMyElements();
  for (unsigned int i = 0; i < numrowcl; ++i)
  {
    int const clgid = cl_noderowmap_prior_redistr_->GID(i);
    Teuchos::RCP<BEAMINTERACTION::DATA::CrosslinkerData> cl_data_i =
        crosslinker_data_[cl_nodecolmap_prior_redistr_->LID(clgid)];

    DRT::PackBuffer data;
    cl_data_i->Pack(data);
    data.StartPacking();
    cl_data_i->Pack(data);
    allpacks[clgid].insert(allpacks[clgid].end(), data().begin(), data().end());
  }

  // export
  cl_exporter_->Export(allpacks);

  // rebuild data container
  crosslinker_data_.resize(BinDiscret().NumMyColNodes());
  for (auto& iter : allpacks)
  {
    // this needs to done
    // Fixme
    std::vector<char>::size_type position = 0;
    std::vector<char> data;
    DRT::ParObject::ExtractfromPack(position, iter.second, data);

    Teuchos::RCP<BEAMINTERACTION::DATA::CrosslinkerData> cl_data = Teuchos::rcp(
        BEAMINTERACTION::DATA::CreateDataContainer<BEAMINTERACTION::DATA::CrosslinkerData>(data),
        true);
    crosslinker_data_[BinDiscret().NodeColMap()->LID(cl_data->GetId())] = cl_data;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UpdateAndExportBeamData(bool update_states)
{
  CheckInit();

  //  TEUCHOS_FUNC_TIME_MONITOR("BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::"
  //      "UpdateAndExportBeamData");
  //
  //  // if maps have change on one proc, all procs have to rebuild the exporter object
  //  int loc_changed_maps = static_cast< int >( (beam_exporter_ == Teuchos::null) or not
  //      ( beam_exporter_->SourceMap().SameAs( *beam_elerowmap_prior_redistr_ ) and
  //        beam_elecolmap_prior_redistr_->SameAs( *Discret().ElementColMap() ) ) );
  //
  //  //global filled flag (is true / one if and only if loc_changed_maps == true on each processor
  //  int g_changed_maps = 0;
  //
  //  /*the global flag is set to the maximal value of any local flag
  //   * i.e. if on any processor loc_changed_maps == true, the flag g_changed_maps is set to
  //   * one*/
  //  Discret().Comm().MaxAll( &loc_changed_maps, &g_changed_maps, 1 );
  //
  //  if ( g_changed_maps )
  beam_exporter_ = Teuchos::rcp(new DRT::Exporter(
      *beam_elerowmap_prior_redistr_, *Discret().ElementColMap(), Discret().Comm()));

  // we first need to pack our row stuff into and std::vector< char > for communication
  std::map<int, std::vector<char>> allpacks;
  unsigned int numrowele = beam_elerowmap_prior_redistr_->NumMyElements();
  for (unsigned int i = 0; i < numrowele; ++i)
  {
    int const elegid = beam_elerowmap_prior_redistr_->GID(i);

    Teuchos::RCP<BEAMINTERACTION::DATA::BeamData> beam_data_i =
        beam_data_[beam_elecolmap_prior_redistr_->LID(elegid)];

    // safety check
    if (beam_data_i == Teuchos::null)
      dserror("beam data container for beam with gid %i not there on rank %i ", elegid,
          GState().GetMyRank());

    if (update_states)
    {
      // safety check
      if (Discret().ElementColMap()->LID(elegid) < 0)
        dserror(" Element %i has moved too far between two redistributions.", elegid);

      // beam element i for which data will be collected
      DRT::ELEMENTS::Beam3Base* beamele_i =
          dynamic_cast<DRT::ELEMENTS::Beam3Base*>(Discret().gElement(elegid));

      // go to next element in case the current one is not a beam element
#ifdef DEBUG
      if (beamele_i == NULL) dserror("cast did not work");
#endif

      std::vector<double> eledisp;
      BEAMINTERACTION::UTILS::GetCurrentUnshiftedElementDis(Discret(), beamele_i,
          BeamInteractionDataStatePtr()->GetMutableDisColNp(), PeriodicBoundingBox(), eledisp);

      // loop over binding spot types of current element
      for (auto const& iter : beamele_i->GetBindingSpots())
      {
        // loop over all binding spots of current type j of current element
        int unsigned const numbbspot = beamele_i->GetNumberOfBindingSpots(iter.first);
        LINALG::Matrix<3, 1> pos(true);
        LINALG::Matrix<3, 3> triad(true);
        for (int unsigned k = 0; k < numbbspot; ++k)
        {
          BEAMINTERACTION::UTILS::GetPosAndTriadOfBindingSpot(
              beamele_i, PeriodicBoundingBoxPtr(), iter.first, k, pos, triad, eledisp);

          beam_data_i->SetBSpotPosition(iter.first, k, pos);
          beam_data_i->SetBSpotTriad(iter.first, k, triad);
        }
      }
    }

    DRT::PackBuffer data;
    beam_data_i->Pack(data);
    data.StartPacking();
    beam_data_i->Pack(data);
    allpacks[elegid].insert(allpacks[elegid].end(), data().begin(), data().end());
    //    allpacks[elegid] = data();
  }

  // export
  beam_exporter_->Export(allpacks);

  // rebuild data container
  beam_data_.resize(Discret().NumMyColElements());
  for (auto& iter : allpacks)
  {
    // this needs to be done
    // Fixme
    std::vector<char>::size_type position = 0;
    std::vector<char> data;
    DRT::ParObject::ExtractfromPack(position, iter.second, data);

    Teuchos::RCP<BEAMINTERACTION::DATA::BeamData> beam_data = Teuchos::rcp(
        BEAMINTERACTION::DATA::CreateDataContainer<BEAMINTERACTION::DATA::BeamData>(data), true);
    beam_data_[Discret().ElementColMap()->LID(beam_data->GetId())] = beam_data;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::BindAndUnbindCrosslinker()
{
  CheckInitSetup();

  TEUCHOS_FUNC_TIME_MONITOR(
      "BEAMINTERACTION::SUBMODELEVALUATOR::"
      "Crosslinking::BindAndUnbindCrosslinker");

  // manage binding events
  int num_new_linker = BindCrosslinker();

  // manage unbinding events
  int num_dissolved_linker = UnBindCrosslinker();

  // write some information to screen
  std::vector<int> num_local(3, 0);
  std::vector<int> num_global(3, 0);
  num_local[0] = static_cast<int>(doublebondcl_.size());
  num_local[1] = num_new_linker;
  num_local[2] = num_dissolved_linker;
  MPI_Reduce(&num_local[0], &num_global[0], 3, MPI_INT, MPI_SUM, 0,
      dynamic_cast<const Epetra_MpiComm*>(&(Discret().Comm()))->Comm());
  if (GState().GetMyRank() == 0)
  {
    IO::cout(IO::standard) << "\n************************************************" << IO::endl;
    IO::cout(IO::standard) << "Beam to Beam Links: " << num_global[0];
    IO::cout(IO::standard) << " (New: " << num_global[1];
    IO::cout(IO::standard) << " Dissolved: " << num_global[2];
    IO::cout(IO::standard) << ")\n************************************************\n" << IO::endl;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::BindCrosslinker()
{
  CheckInit();

  DRT::Problem::Instance()->Random()->SetRandRange(0.0, 1.0);

  // intended bonds of row crosslinker on myrank (key is clgid)
  std::map<int, Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>> mybonds;
  mybonds.clear();
  // intended bond col crosslinker to row element (key is owner of crosslinker != myrank)
  std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>> undecidedbonds;
  undecidedbonds.clear();

  // fill binding event maps
  FindPotentialBindingEvents(mybonds, undecidedbonds);

  // bind events where myrank only owns the elements, cl are taken care
  // of by their owner (key is clgid)
  std::map<int, Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>> myelebonds;

  // now each row owner of a linker gets requests, makes a random decision and
  // informs back its requesters
  ManageBindingInParallel(mybonds, undecidedbonds, myelebonds);

  // actual update of binding states is done here
  return UpdateMyCrosslinkerAndElementBindingStates(mybonds, myelebonds);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::FindPotentialBindingEvents(
    std::map<int, Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>& mybonds,
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>>& undecidedbonds)
{
  CheckInit();

  TEUCHOS_FUNC_TIME_MONITOR(
      "BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::FindPotentialBindingEvents");

  // this variable is used to check if a beam binding spot is linked twice on
  // myrank during a time step
  // ( first key is linkertype, second key is locbspotid, set holds gids of bonded crosslinker)
  std::map<int, std::vector<std::map<int, std::set<int>>>> intendedbeambonds;
  for (unsigned int i = 0; i < crosslinking_params_ptr_->LinkerTypes().size(); ++i)
    intendedbeambonds[crosslinking_params_ptr_->LinkerTypes()[i]].resize(
        DiscretPtr()->NumMyRowElements());

  // store bins that have already been examined
  std::vector<int> examinedbins(BinDiscretPtr()->NumMyColElements(), 0);
  // loop over all column crosslinker in random order
  // create random order of indices
  std::vector<int> rordercolcl =
      BEAMINTERACTION::UTILS::Permutation(BinDiscretPtr()->NumMyColNodes());

  for (auto const& icl : rordercolcl)
  {
    DRT::Node* currcrosslinker = BinDiscretPtr()->lColNode(icl);

#ifdef DEBUG
    if (currcrosslinker == NULL) dserror("Node not there");
    if (currcrosslinker->NumElement() != 1) dserror("More than one element for this crosslinker");
#endif

    // get bin that contains this crosslinker (can only be one)
    DRT::Element* currentbin = currcrosslinker->Elements()[0];

#ifdef DEBUG
    if (currentbin->Id() < 0) dserror(" negative bin id number %i ", currentbin->Id());
#endif

    // if a bin has already been examined --> continue with next crosslinker
    if (examinedbins[currentbin->LID()]) continue;
    // else: bin is examined for the first time --> new entry in examinedbins_
    else
      examinedbins[currentbin->LID()] = 1;

    FindPotentialBindingEventsInBinAndNeighborhood(
        currentbin, mybonds, undecidedbonds, intendedbeambonds, true);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::
    FindPotentialBindingEventsInBinAndNeighborhood(DRT::Element* bin,
        std::map<int, Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>& mybonds,
        std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>>&
            undecidedbonds,
        std::map<int, std::vector<std::map<int, std::set<int>>>>& intendedbeambonds,
        bool checklinkingprop)
{
  CheckInit();

  // get neighboring bins
  // note: interaction distance cl to beam needs to be smaller than the bin size
  std::vector<int> neighboring_binIds;
  neighboring_binIds.reserve(27);
  // do not check on existence here -> shifted to GetBinContent
  BinStrategyPtr()->GetNeighborAndOwnBinIds(bin->Id(), neighboring_binIds);

  // get set of neighboring beam elements (i.e. elements that somehow touch nb bins)
  // as explained above, we only need row elements (true flag in GetBinContent())
  std::set<DRT::Element*> neighboring_row_beams;
  std::vector<BINSTRATEGY::UTILS::BinContentType> bc_beam(1, BINSTRATEGY::UTILS::Beam);
  BinStrategyPtr()->GetBinContent(neighboring_row_beams, bc_beam, neighboring_binIds, true);
  std::set<DRT::Element*> neighboring_col_spheres;
  std::vector<BINSTRATEGY::UTILS::BinContentType> bc_sphere(1, BINSTRATEGY::UTILS::RigidSphere);
  BinStrategyPtr()->GetBinContent(neighboring_col_spheres, bc_sphere, neighboring_binIds, false);


  // in case there are no neighbors, go to next crosslinker (an therefore bin)
  if (neighboring_row_beams.empty()) return;

  // get all crosslinker in current bin
  DRT::Node** clincurrentbin = bin->Nodes();
  const int numcrosslinker = bin->NumNode();

  // obtain random order in which crosslinker are addressed
  std::vector<int> randorder = BEAMINTERACTION::UTILS::Permutation(numcrosslinker);

  // loop over all crosslinker in CurrentBin in random order
  for (auto const& randcliter : randorder)
  {
    // get random crosslinker in current bin
    DRT::Node* crosslinker_i = clincurrentbin[randcliter];

    // todo: this can be done more efficiently
    if (CheckIfSphereProhibitsBinding(neighboring_col_spheres, crosslinker_i)) continue;

    // get all potential binding events on myrank
    PrepareBinding(crosslinker_i, neighboring_row_beams, mybonds, undecidedbonds, intendedbeambonds,
        checklinkingprop);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::CheckIfSphereProhibitsBinding(
    std::set<DRT::Element*> const& neighboring_col_spheres, DRT::Node* node_i) const
{
  CheckInit();

  CROSSLINKING::CrosslinkerNode* crosslinker_i =
      dynamic_cast<CROSSLINKING::CrosslinkerNode*>(node_i);

  if (std::abs(crosslinker_i->GetMaterial()->NoBondDistSphere()) < 1.0e-8) return false;

  BEAMINTERACTION::DATA::CrosslinkerData* cldata_i = crosslinker_data_[crosslinker_i->LID()].get();

  for (auto const& sphere_iter : neighboring_col_spheres)
  {
    // init position of linker nodes
    LINALG::Matrix<3, 1> sphere_pos(true);

    // sphere current position
    std::vector<double> sphereeledisp;
    BEAMINTERACTION::UTILS::GetCurrentElementDis(
        Discret(), sphere_iter, BeamInteractionDataState().GetDisColNp(), sphereeledisp);

    // note: sphere has just one node (with three translational dofs)
    for (unsigned int dim = 0; dim < 3; ++dim)
      sphere_pos(dim) = sphere_iter->Nodes()[0]->X()[dim] + sphereeledisp[dim];

    LINALG::Matrix<3, 1> dist_vec(true);
    dist_vec.Update(1.0, sphere_pos, -1.0, cldata_i->GetPosition());
    const double distance = dist_vec.Norm2();

    if (distance < crosslinker_i->GetMaterial()->NoBondDistSphere()) return true;
  }

  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PrepareBinding(DRT::Node* node_i,
    std::set<DRT::Element*> const& neighboring_beams,
    std::map<int, Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>& mybonds,
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>>& undecidedbonds,
    std::map<int, std::vector<std::map<int, std::set<int>>>>& intendedbeambonds,
    bool checklinkingprop)
{
  CheckInit();

  CROSSLINKING::CrosslinkerNode* crosslinker_i =
      dynamic_cast<CROSSLINKING::CrosslinkerNode*>(node_i);

  // get precomputed data of crosslinker i
  BEAMINTERACTION::DATA::CrosslinkerData* cldata_i = crosslinker_data_[crosslinker_i->LID()].get();

  // -------------------------------------------------------------------------
  // We now check all criteria that need to be passed for a binding event one
  // after the other
  // -------------------------------------------------------------------------
  // 1. criterion: in case crosslinker is double bonded, we can leave here
  if (cldata_i->GetNumberOfBonds() == 2) return;

  // loop over all neighboring beam elements in random order (keep in mind
  // we are only looping over row elements)
  std::vector<DRT::Element*> beamvec(neighboring_beams.begin(), neighboring_beams.end());

  // -------------------------------------------------------------------------
  // NOTE: This is crucial for reproducibility to ensure that computation does
  // not depend on pointer addresses (see also comment of class Less)
  // -------------------------------------------------------------------------
  std::sort(beamvec.begin(), beamvec.end(), BEAMINTERACTION::UTILS::Less());

  std::vector<int> randorder =
      BEAMINTERACTION::UTILS::Permutation(static_cast<int>(beamvec.size()));
  for (auto const& randiter : randorder)
  {
    // get neighboring (nb) beam element
    DRT::Element* nbbeam = beamvec[randiter];

    // get pre computed data of current nbbeam
    BEAMINTERACTION::DATA::BeamData* beamdata_i = beam_data_[nbbeam->LID()].get();

    if (cldata_i->GetNumberOfBonds() == 1)
    {
      int cl_bondedtogid =
          cldata_i->GetBSpots()[GetSingleOccupiedClBspot(cldata_i->GetBSpots())].first;

      // safety check
      if (Discret().ElementColMap()->LID(cl_bondedtogid) < 0)
        dserror("Element %i not ghosted on rank %i", cl_bondedtogid, GState().GetMyRank());

      // 2. criterion:
      // exclude binding of a single bonded crosslinker in close proximity on the
      // same filament (i.e. element cloud of old element binding partner is excluded)
      if (BEAMINTERACTION::UTILS::DoBeamElementsShareNodes(
              Discret().gElement(cl_bondedtogid), nbbeam))
        continue;
    }

    // loop over all binding spots of current element in random order
    std::vector<int> randbspot = BEAMINTERACTION::UTILS::Permutation(
        beamdata_i->GetNumberOfBindingSpotsOfType(crosslinker_i->GetMaterial()->LinkerType()));

    for (auto const& rbspotiter : randbspot)
    {
      // get local number of binding spot in element
      const int locnbspot = rbspotiter;

      // we are now doing some additional checks if a binding event is feasible
      if (not CheckBindEventCriteria(crosslinker_i, nbbeam, cldata_i, beamdata_i, locnbspot,
              intendedbeambonds, checklinkingprop))
        continue;

      // ---------------------------------------------------------------------
      // if we made it this far, we can add this potential binding event to its
      // corresponding map
      // ---------------------------------------------------------------------
      Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData> bindeventdata =
          Teuchos::rcp(new BEAMINTERACTION::DATA::BindEventData());
      // default permission is true, is changed if owner of cl has something against it
      bindeventdata->Init(crosslinker_i->Id(), nbbeam->Id(), locnbspot, GState().GetMyRank(), 1);

      // in case myrank is owner, we add it to the mybonds map
      if (crosslinker_i->Owner() == GState().GetMyRank())
      {
        mybonds[bindeventdata->GetClId()] = bindeventdata;
      }
      else
      {
        // myrank is not owner, we add it to the map of events that need to be
        // communicated to make a decision
        undecidedbonds[crosslinker_i->Owner()].push_back(bindeventdata);
      }

      // as we allow only one binding event for each cl in one time step,
      // we are done here, if we made it so far (i.e met criteria 1. - 7.)
      return;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::CheckBindEventCriteria(
    CROSSLINKING::CrosslinkerNode const* const crosslinker_i,
    DRT::Element const* const potbeampartner, BEAMINTERACTION::DATA::CrosslinkerData* cldata_i,
    BEAMINTERACTION::DATA::BeamData const* beamdata_i, int locnbspot,
    std::map<int, std::vector<std::map<int, std::set<int>>>>& intendedbeambonds,
    bool checklinkingprop) const
{
  CheckInit();

  int const potbeampartnerrowlid = Discret().ElementRowMap()->LID(potbeampartner->Id());
  INPAR::BEAMINTERACTION::CrosslinkerType linkertype = crosslinker_i->GetMaterial()->LinkerType();

  // check compatibility of crosslinker type and filament type (some linker can only
  // bind to certain filament types)
  if (not CheckLinkerAndFilamentTypeCompatibility(linkertype,
          dynamic_cast<DRT::ELEMENTS::Beam3Base const* const>(potbeampartner)->GetFilamentType()))
    return false;

  // a crosslink is set if and only if it passes the probability check
  // for a binding event to happen
  double plink = 1.0 - exp((-1.0) * crosslinking_params_ptr_->DeltaTime() *
                           crosslinker_i->GetMaterial()->KOn());

  if (checklinkingprop and (DRT::Problem::Instance()->Random()->Uni() > plink)) return false;

  // criterion:
  // first check if binding spot has free bonds left
  if (static_cast<int>(beamdata_i->GetBSpotStatusAt(linkertype, locnbspot).size()) >=
      crosslinking_params_ptr_->MaxNumberOfBondsPerFilamentBspot(linkertype))
    return false;

  // exclude multiple identical crosslinks
  if (not ReturnFalseIfIdenticalBondAlreadyExists(
          crosslinker_i, cldata_i, intendedbeambonds, beamdata_i, locnbspot, potbeampartnerrowlid))
    return false;

  /* check RELEVANT distance criterion
   * if free:
   *   distance between crosslinker center and current beam binding spot
   * if singly bound:
   *   distance between already bound bspot of crosslinker and current beam binding spot
   * note: as we set the crosslinker position to coincide with beam bspot position if singly bound,
   *       we can also use cldata_i.clpos in the second case*/

  // get current position and tangent vector of filament axis at free binding spot
  LINALG::Matrix<3, 1> const& currbbspos = beamdata_i->GetBSpotPosition(linkertype, locnbspot);

  // minimum and maximum distance at which a double-bond crosslink can be established
  double const linkdistmin = crosslinker_i->GetMaterial()->LinkingLength() -
                             crosslinker_i->GetMaterial()->LinkingLengthTolerance();
  double const linkdistmax = crosslinker_i->GetMaterial()->LinkingLength() +
                             crosslinker_i->GetMaterial()->LinkingLengthTolerance();

#ifdef DEBUG
  // safety check
  if (linkdistmax > BinStrategy().CutoffRadius())
    dserror(
        "The allowed binding distance of linker %i (in case it is single bonded) is"
        "\ngreater than the cutoff radius, this could lead to missing a binding event",
        crosslinker_i->Id());
#endif

  if ((cldata_i->GetNumberOfBonds() == 0 and
          BEAMINTERACTION::UTILS::IsDistanceOutOfRange(
              cldata_i->GetPosition(), currbbspos, 0.5 * linkdistmin, 0.5 * linkdistmax)) or
      (cldata_i->GetNumberOfBonds() == 1 and
          BEAMINTERACTION::UTILS::IsDistanceOutOfRange(
              cldata_i->GetPosition(), currbbspos, linkdistmin, linkdistmax)))
    return false;

  // orientation of centerline tangent vectors at binding spots
  // a crosslink (double-bonded crosslinker) will only be established if the
  // enclosed angle is in the specified range
  double const linkanglemin = crosslinker_i->GetMaterial()->LinkingAngle() -
                              crosslinker_i->GetMaterial()->LinkingAngleTolerance();
  double const linkanglemax = crosslinker_i->GetMaterial()->LinkingAngle() +
                              crosslinker_i->GetMaterial()->LinkingAngleTolerance();

  // if crosslinker is singly bound, we fetch the orientation vector
  LINALG::Matrix<3, 1> occ_bindingspot_beam_tangent(true);
  if (cldata_i->GetNumberOfBonds() == 1)
    GetOccupiedClBSpotBeamTangent(
        crosslinker_i, cldata_i, occ_bindingspot_beam_tangent, crosslinker_i->Id());

  // note: we use first base vector instead of tangent vector here
  LINALG::Matrix<3, 1> curr_bindingspot_beam_tangent(true);
  for (unsigned int idim = 0; idim < 3; ++idim)
    curr_bindingspot_beam_tangent(idim) = beamdata_i->GetBSpotTriad(linkertype, locnbspot)(idim, 0);

  if (cldata_i->GetNumberOfBonds() == 1 and
      BEAMINTERACTION::UTILS::IsEnclosedAngleOutOfRange(
          occ_bindingspot_beam_tangent, curr_bindingspot_beam_tangent, linkanglemin, linkanglemax))
    return false;

  // check if current beam binding spot yet intended to bind this timestep
  // by a crosslinker that came before in this random order
  if (static_cast<int>(intendedbeambonds.at(linkertype)[potbeampartnerrowlid][locnbspot].size() +
                       beamdata_i->GetBSpotStatusAt(linkertype, locnbspot).size()) >=
      crosslinking_params_ptr_->MaxNumberOfBondsPerFilamentBspot(
          crosslinker_i->GetMaterial()->LinkerType()))
  {
    /* note: it is possible that the binding event that rejects the current one is rejected itself
     * later during communication with other procs and therefore the current one could be
     * valid. Just neglecting this here is a slight inconsistency, but should be ok as such an
     * coincidence is extremely rare in a simulation with realistic proportion of crosslinker
     * to beam binding spots. Additionally missing one event would not change any physics.
     * (Could be cured with additional communication)
     */
    if (Discret().Comm().NumProc() > 1)
      IO::cout(IO::verbose)
          << " Warning: There is a minimal chance of missing a regular binding event on "
             "rank "
          << GState().GetMyRank() << IO::endl;
    return false;
  }
  else
  {
    intendedbeambonds.at(linkertype)[potbeampartnerrowlid][locnbspot].insert(crosslinker_i->Id());
  }

  // bind event can happen
  // (rejection afterwards only possible in case of parallel simulation)
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ReturnFalseIfIdenticalBondAlreadyExists(
    CROSSLINKING::CrosslinkerNode const* const crosslinker_i,
    BEAMINTERACTION::DATA::CrosslinkerData* cldata_i,
    std::map<int, std::vector<std::map<int, std::set<int>>>>& intendedbeambonds,
    BEAMINTERACTION::DATA::BeamData const* beamdata_i, int locnbspot,
    int potbeampartnerrowlid) const
{
  INPAR::BEAMINTERACTION::CrosslinkerType linkertype = crosslinker_i->GetMaterial()->LinkerType();

  if (not(cldata_i->GetNumberOfBonds() == 1 and
          crosslinking_params_ptr_->MaxNumberOfBondsPerFilamentBspot(linkertype) > 1))
    return true;

  // get element and gid of beam element to which current linker is already bonded to
  int occbspotid = GetSingleOccupiedClBspot(cldata_i->GetBSpots());
  int const elegid = cldata_i->GetBSpots()[occbspotid].first;
  int const locbspotnum = cldata_i->GetBSpots()[occbspotid].second;

  int const elecollid = Discret().ElementColMap()->LID(elegid);

#ifdef DEBUG
  // safety check
  if (elecollid < 0) dserror("element with gid %i not on proc %i", elegid, GState().GetMyRank());
#endif

  // loop over crosslinker that are already bonded to binding spot to which current linker is bonded
  // to
  for (auto const& iter : beam_data_[elecollid]->GetBSpotStatusAt(linkertype, locbspotnum))
  {
    // this is needed in case a binding event was allowed in this time step in opposite direction
    if (intendedbeambonds.at(linkertype)[potbeampartnerrowlid][locnbspot].find(iter) !=
        intendedbeambonds.at(linkertype)[potbeampartnerrowlid][locnbspot].end())
      return false;
  }

  // loop over crosslinker that are already bonded to potential new binding spot
  for (auto const iter : beamdata_i->GetBSpotStatusAt(linkertype, locnbspot))
  {
    DRT::Node* bonded_crosslinker_i = BinDiscret().gNode(iter);

#ifdef DEBUG
    // safety check
    if (bonded_crosslinker_i == NULL)
      dserror(" Linker with gid %i not on rank %i", iter, GState().GetMyRank());
#endif

    BEAMINTERACTION::DATA::CrosslinkerData const* bondedcl_data_i =
        crosslinker_data_[bonded_crosslinker_i->LID()].get();

#ifdef DEBUG
    // safety check
    if (bonded_crosslinker_i == NULL)
      dserror("Data for crosslinker %i not there on rank %i", bonded_crosslinker_i->Id(),
          GState().GetMyRank());
#endif

    if (bondedcl_data_i->GetNumberOfBonds() == 1)
    {
      // this is needed in case a binding event was allowed in this time step in opposite direction
      int const elelid = Discret().ElementRowMap()->LID(elegid);
      if (elelid != -1 and
          intendedbeambonds.at(linkertype)[elelid][locbspotnum].find(bonded_crosslinker_i->Id()) !=
              intendedbeambonds.at(linkertype)[elelid][locbspotnum].end())
        return false;
    }
    else if (bondedcl_data_i->GetNumberOfBonds() == 2)
    {
      // if intended bond between two filament binding spots already exists, reject current intended
      // bond
      if ((bondedcl_data_i->GetBSpots()[0].first == elegid and
              bondedcl_data_i->GetBSpots()[0].second == locbspotnum) or
          (bondedcl_data_i->GetBSpots()[1].first == elegid and
              bondedcl_data_i->GetBSpots()[1].second == locbspotnum))
        return false;
    }
    else
    {
      dserror(
          " unrealistic number of bonds (%i) for crosslinker (gid %i) at this point. Beam %i local "
          "%i ",
          bondedcl_data_i->GetNumberOfBonds(), bonded_crosslinker_i->Owner(), elegid, locbspotnum);
    }
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::CheckLinkerAndFilamentTypeCompatibility(
    INPAR::BEAMINTERACTION::CrosslinkerType linkertype,
    INPAR::BEAMINTERACTION::FilamentType filamenttype) const
{
  switch (linkertype)
  {
    case INPAR::BEAMINTERACTION::linkertype_arbitrary:
    {
      // no check of filament type necessary
      return true;
      break;
    }
    case INPAR::BEAMINTERACTION::linkertype_actin:
    {
      if (filamenttype == INPAR::BEAMINTERACTION::filtype_actin)
        return true;
      else
        return false;
      break;
    }
    case INPAR::BEAMINTERACTION::linkertype_collagen:
    {
      if (filamenttype == INPAR::BEAMINTERACTION::filtype_collagen)
        return true;
      else
        return false;
      break;
    }
    default:
    {
      dserror("Unknown linker type.");
      exit(EXIT_FAILURE);
    }
  }

  // default false
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::GetOccupiedClBSpotBeamTangent(
    CROSSLINKING::CrosslinkerNode const* const crosslinker_i,
    BEAMINTERACTION::DATA::CrosslinkerData* cldata_i,
    LINALG::Matrix<3, 1>& occ_bindingspot_beam_tangent, int clgid) const
{
  CheckInitSetup();

  int occbspotid = GetSingleOccupiedClBspot(cldata_i->GetBSpots());

  const int locbspotnum = cldata_i->GetBSpots()[occbspotid].second;
  const int elegid = cldata_i->GetBSpots()[occbspotid].first;
  const int elecollid = Discret().ElementColMap()->LID(elegid);

#ifdef DEBUG
  if (elecollid < 0)
    dserror(" Element with gid %i bonded to cl %i on rank %i not even ghosted", elegid, clgid,
        GState().GetMyRank());
#endif

  // note: we use first base vector instead of tangent vector here
  for (unsigned int idim = 0; idim < 3; ++idim)
    occ_bindingspot_beam_tangent(idim) = beam_data_[elecollid]->GetBSpotTriad(
        crosslinker_i->GetMaterial()->LinkerType(), locbspotnum)(idim, 0);
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ManageBindingInParallel(
    std::map<int, Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>& mybonds,
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>>& undecidedbonds,
    std::map<int, Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>& myelebonds) const
{
  CheckInit();

  TEUCHOS_FUNC_TIME_MONITOR(
      "BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ManageBindingInParallel");

  // -------------------------------------------------------------------------
  // 1) each procs makes his requests and receives the request of other procs
  // -------------------------------------------------------------------------
  // store requested cl and its data
  std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>> requestedcl;
  CommunicateUndecidedBonds(undecidedbonds, requestedcl);

  // -------------------------------------------------------------------------
  // 2) now myrank needs to decide which proc is allowed to set the requested
  //    link
  // -------------------------------------------------------------------------
  std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>> decidedbonds;
  DecideBindingInParallel(requestedcl, mybonds, decidedbonds);

  // -------------------------------------------------------------------------
  // 3) communicate the binding decisions made on myrank, receive decisions
  //    made for its own requests and create colbondmap accordingly
  // -------------------------------------------------------------------------
  CommunicateDecidedBonds(decidedbonds, myelebonds);
  decidedbonds.clear();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UpdateMyCrosslinkerAndElementBindingStates(
    std::map<int, Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>& mybonds,
    std::map<int, Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>& myelebonds)
{
  CheckInit();

  // map key is crosslinker gid to be able to uniquely address one entry over all procs
  std::map<int, NewDoubleBonds> mynewdbondcl;

  // myrank owner of crosslinker and most elements
  UpdateMyCrosslinkerBindingStates(mybonds, mynewdbondcl);

  // myrank only owner of current binding partner ele
  UpdateMyElementBindingStates(myelebonds);

  // setup new double bonds and insert them in doublebondcl_
  CreateNewDoubleBondedCrosslinkerElementPairs(mynewdbondcl);

  return static_cast<int>(mynewdbondcl.size());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UpdateMyCrosslinkerBindingStates(
    std::map<int, Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>> const& mybonds,
    std::map<int, NewDoubleBonds>& mynewdbondcl)
{
  CheckInit();

  for (auto const& cliter : mybonds)
  {
    // get binding event data
    BEAMINTERACTION::DATA::BindEventData* binevdata = cliter.second.get();

#ifdef DEBUG
    if (binevdata->GetPermission() != 1)
      dserror(
          " Rank %i wants to bind crosslinker %i without permission, "
          " something went wrong",
          GState().GetMyRank(), cliter.first);
#endif

    // get current linker and beam data
    const int clcollid = BinDiscretPtr()->NodeColMap()->LID(cliter.first);
    const int colelelid = DiscretPtr()->ElementColMap()->LID(binevdata->GetEleId());

#ifdef DEBUG
    // safety checks
    if (clcollid < 0) dserror("Crosslinker not even ghosted, should be owned here.");
    if (colelelid < 0) dserror("Element with gid %i not ghosted.", binevdata->GetEleId());
#endif
    CROSSLINKING::CrosslinkerNode* crosslinker_i =
        dynamic_cast<CROSSLINKING::CrosslinkerNode*>(BinDiscretPtr()->lColNode(clcollid));

    // get crosslinker data
    BEAMINTERACTION::DATA::CrosslinkerData* cldata_i = crosslinker_data_[clcollid].get();

    DRT::Element* beamele_i = DiscretPtr()->lColElement(colelelid);
    BEAMINTERACTION::DATA::BeamData* beamdata_i = beam_data_[colelelid].get();

#ifdef DEBUG
    // safety checks
    if (cliter.first != binevdata->GetClId())
      dserror("Map key does not match crosslinker gid of current binding event.");

    if (crosslinker_i->Owner() != GState().GetMyRank())
      dserror("Only row owner of crosslinker is changing its status");

    if (colelelid < 0)
      dserror(
          "Binding element partner of current row crosslinker is not ghosted, "
          "this must be the case though.");
#endif

    // -------------------------------------------------------------------------
    // different treatment according to number of bonds crosslinker had before
    // this binding event
    // -------------------------------------------------------------------------
    switch (cldata_i->GetNumberOfBonds())
    {
      case 0:
      {
        // -----------------------------------------------------------------
        // update crosslinker status
        // -----------------------------------------------------------------
        // store gid and bspot local number of this element, first binding spot
        // always bonded first
        cldata_i->SetBspot(0, std::make_pair(binevdata->GetEleId(), binevdata->GetBSpotLocN()));

        // update number of bonds
        cldata_i->SetNumberOfBonds(1);

        // update position
        std::vector<double> newpos(3, 0.0);
        for (int dim = 0; dim < 3; ++dim)
          newpos[dim] = beamdata_i->GetBSpotPosition(
              crosslinker_i->GetMaterial()->LinkerType(), binevdata->GetBSpotLocN())(dim);
        crosslinker_i->SetPos(newpos);
        cldata_i->SetPosition(beamdata_i->GetBSpotPosition(
            crosslinker_i->GetMaterial()->LinkerType(), binevdata->GetBSpotLocN()));

        // -----------------------------------------------------------------
        // update beam status
        // -----------------------------------------------------------------
        // store crosslinker gid in status of beam binding spot if myrank
        // is owner of beam
        if (beamele_i->Owner() == GState().GetMyRank())
          beamdata_i->AddBondToBindingSpot(crosslinker_i->GetMaterial()->LinkerType(),
              binevdata->GetBSpotLocN(), binevdata->GetClId());

#ifdef DEBUG
        // safety check
        if (not(cldata_i->GetBSpots()[1].first < 0))
          dserror("Numbond does not fit to clbspot vector.");
#endif

        break;
      }
      case 1:
      {
        // get clbspot that is currently bonded
        int occbspotid = GetSingleOccupiedClBspot(cldata_i->GetBSpots());
        int freebspotid = 1;
        if (occbspotid == 1) freebspotid = 0;

        // -----------------------------------------------------------------
        // update crosslinker status
        // -----------------------------------------------------------------
        // store gid and bspot local number of this element
        cldata_i->SetBspot(
            freebspotid, std::make_pair(binevdata->GetEleId(), binevdata->GetBSpotLocN()));

        // update number of bonds
        cldata_i->SetNumberOfBonds(2);

        // update position
        LINALG::Matrix<3, 1> clpos(cldata_i->GetPosition());
        SetPositionOfDoubleBondedCrosslinkerPBCconsistent(clpos,
            beamdata_i->GetBSpotPosition(crosslinker_i->GetMaterial()->LinkerType(),
                cldata_i->GetBSpots()[freebspotid].second),
            cldata_i->GetPosition());

        std::vector<double> newpos(3, 0.0);
        for (int dim = 0; dim < 3; ++dim) newpos[dim] = clpos(dim);
        crosslinker_i->SetPos(newpos);

        cldata_i->SetPosition(clpos);

        // create double bond cl data
        int occ_colelelid =
            DiscretPtr()->ElementColMap()->LID(cldata_i->GetBSpots()[occbspotid].first);
        NewDoubleBonds dbondcl;
        dbondcl.id = binevdata->GetClId();
        if (cldata_i->GetBSpots()[freebspotid].first > cldata_i->GetBSpots()[occbspotid].first)
        {
          dbondcl.eleids.push_back(cldata_i->GetBSpots()[freebspotid]);
          dbondcl.eleids.push_back(cldata_i->GetBSpots()[occbspotid]);
          dbondcl.bspotposs.push_back(
              beam_data_[colelelid]->GetBSpotPosition(crosslinker_i->GetMaterial()->LinkerType(),
                  cldata_i->GetBSpots()[freebspotid].second));
          dbondcl.bspotposs.push_back(beam_data_[occ_colelelid]->GetBSpotPosition(
              crosslinker_i->GetMaterial()->LinkerType(),
              cldata_i->GetBSpots()[occbspotid].second));
          dbondcl.bspottriads.push_back(
              beam_data_[colelelid]->GetBSpotTriad(crosslinker_i->GetMaterial()->LinkerType(),
                  cldata_i->GetBSpots()[freebspotid].second));
          dbondcl.bspottriads.push_back(
              beam_data_[occ_colelelid]->GetBSpotTriad(crosslinker_i->GetMaterial()->LinkerType(),
                  cldata_i->GetBSpots()[occbspotid].second));
        }
        else
        {
          dbondcl.eleids.push_back(cldata_i->GetBSpots()[occbspotid]);
          dbondcl.eleids.push_back(cldata_i->GetBSpots()[freebspotid]);
          dbondcl.bspotposs.push_back(beam_data_[occ_colelelid]->GetBSpotPosition(
              crosslinker_i->GetMaterial()->LinkerType(),
              cldata_i->GetBSpots()[occbspotid].second));
          dbondcl.bspotposs.push_back(
              beam_data_[colelelid]->GetBSpotPosition(crosslinker_i->GetMaterial()->LinkerType(),
                  cldata_i->GetBSpots()[freebspotid].second));
          dbondcl.bspottriads.push_back(
              beam_data_[occ_colelelid]->GetBSpotTriad(crosslinker_i->GetMaterial()->LinkerType(),
                  cldata_i->GetBSpots()[occbspotid].second));
          dbondcl.bspottriads.push_back(
              beam_data_[colelelid]->GetBSpotTriad(crosslinker_i->GetMaterial()->LinkerType(),
                  cldata_i->GetBSpots()[freebspotid].second));
        }

        // insert pair in mypairs
        mynewdbondcl[dbondcl.id] = dbondcl;

        // first check if myrank is owner of element of current binding event
        // (additionally to being owner of cl)
        if (beamele_i->Owner() == GState().GetMyRank())
          beamdata_i->AddBondToBindingSpot(crosslinker_i->GetMaterial()->LinkerType(),
              binevdata->GetBSpotLocN(), binevdata->GetClId());

        break;
      }
      default:
      {
        dserror(
            "You should not be here, crosslinker has unrealistic number "
            "%i of bonds.",
            cldata_i->GetNumberOfBonds());
        exit(EXIT_FAILURE);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UpdateMyElementBindingStates(
    std::map<int, Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>> const& myelebonds)
{
  CheckInit();

  /*
   * 1| 2__|2  or  1| 3__|2  or  1| 2__|1
   * 1|    |2      1|    |2      1|    |1
   * legend: | = beam; __= cl; 2,3 = owner; 1 = myrank
   */
  // loop through all binding events
  for (auto const& cliter : myelebonds)
  {
    // get binding event data
    BEAMINTERACTION::DATA::BindEventData* binevdata = cliter.second.get();

    // get linker data and beam data
    CROSSLINKING::CrosslinkerNode* linker =
        dynamic_cast<CROSSLINKING::CrosslinkerNode*>(BinDiscretPtr()->gNode(cliter.first));
    int const colelelid = DiscretPtr()->ElementColMap()->LID(binevdata->GetEleId());

#ifdef DEBUG
    // safety checks
    if (linker == NULL) dserror("Crosslinker needs to be ghosted, but this isn't the case.");
    if (colelelid < 0)
      dserror("element with gid %i not ghosted on proc %i", binevdata->GetEleId(),
          GState().GetMyRank());
#endif

    // linker
    int const clcollid = linker->LID();
    BEAMINTERACTION::DATA::CrosslinkerData* cldata_i = crosslinker_data_[clcollid].get();

    BEAMINTERACTION::DATA::BeamData* beamdata_i = beam_data_[colelelid].get();

#ifdef DEBUG
    // safety checks
    if (DiscretPtr()->lColElement(colelelid)->Owner() != GState().GetMyRank())
      dserror("Only row owner of element is allowed to change its status");
    if (linker->Owner() == GState().GetMyRank())
      dserror("myrank should not be owner of this crosslinker");
#endif

    // different treatment according to number of bonds crosslinker had before
    // this binding event
    switch (cldata_i->GetNumberOfBonds())
    {
      case 0:
      {
        beamdata_i->AddBondToBindingSpot(
            linker->GetMaterial()->LinkerType(), binevdata->GetBSpotLocN(), binevdata->GetClId());
        break;
      }
      case 1:
      {
        beamdata_i->AddBondToBindingSpot(
            linker->GetMaterial()->LinkerType(), binevdata->GetBSpotLocN(), binevdata->GetClId());
        break;
      }
      default:
      {
        dserror(
            "You should not be here, crosslinker has unrealistic number "
            "%i of bonds.",
            cldata_i->GetNumberOfBonds());
        exit(EXIT_FAILURE);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::CreateNewDoubleBondedCrosslinkerElementPairs(
    std::map<int, NewDoubleBonds> const& mynewdbondcl)
{
  CheckInit();

  for (auto const& iter : mynewdbondcl)
  {
    NewDoubleBonds const& newdoublebond_i = iter.second;
    CROSSLINKING::CrosslinkerNode* cl_node =
        dynamic_cast<CROSSLINKING::CrosslinkerNode*>(BinDiscretPtr()->gNode(iter.first));

    // create and initialize objects of beam-to-beam connections
    // Todo move this inside the create routines (or one create routine in BeamLink class)
    Teuchos::RCP<BEAMINTERACTION::BeamLink> linkelepairptr;
    if (cl_node->GetMaterial()->JointType() == INPAR::BEAMINTERACTION::beam3r_lin2_rigid)
      linkelepairptr = BEAMINTERACTION::BeamLinkRigidJointed::Create();
    else if (cl_node->GetMaterial()->JointType() == INPAR::BEAMINTERACTION::beam3r_lin2_pin or
             cl_node->GetMaterial()->JointType() == INPAR::BEAMINTERACTION::truss)
      linkelepairptr =
          BEAMINTERACTION::BeamLinkPinJointed::Create(cl_node->GetMaterial()->JointType());

    // finally initialize and setup object
    linkelepairptr->Init(iter.first, newdoublebond_i.eleids, newdoublebond_i.bspotposs,
        newdoublebond_i.bspottriads, cl_node->GetMaterial()->LinkerType(), GState().GetTimeNp());
    linkelepairptr->Setup(cl_node->GetMaterial()->BeamElastHyperMatNum());

    // add to my double bonds
    doublebondcl_[linkelepairptr->Id()] = linkelepairptr;

#ifdef DEBUG
    // safety check
    DRT::Node* crosslinker_i = BinDiscretPtr()->gNode(linkelepairptr->Id());

    // safety check
    BEAMINTERACTION::DATA::CrosslinkerData const* cldata_i =
        crosslinker_data_[crosslinker_i->LID()].get();

    if (cldata_i->GetNumberOfBonds() != 2)
      dserror("Setup: Cl with gid %i Owner %i on myrank %i and numbonds %i", linkelepairptr->Id(),
          crosslinker_i->Owner(), GStatePtr()->GetMyRank(), cldata_i->GetNumberOfBonds());
#endif
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::DoubleBindCrosslinkerInBinsAndNeighborhood(
    std::set<int> const& bingids)
{
  CheckInit();

  unsigned int maxbonds = 2;

  for (unsigned int bond_i = 0; bond_i < maxbonds; ++bond_i)
  {
    // intended bonds of row crosslinker on myrank (key is clgid)
    std::map<int, Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>> mybonds;
    // intended bond col crosslinker to row element (key is owner of crosslinker != myrank)
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>> undecidedbonds;

    // store bins that have already been examined
    std::vector<int> examinedbins(BinDiscretPtr()->NumMyColElements(), 0);
    // this variable is used to check if a beam binding spot is linked twice on
    // myrank during a time step
    // ( first key is linkertype, second key is locbspotid, set holds gids of bonded crosslinker)
    std::map<int, std::vector<std::map<int, std::set<int>>>> intendedbeambonds;
    for (unsigned int i = 0; i < crosslinking_params_ptr_->LinkerTypes().size(); ++i)
      intendedbeambonds[crosslinking_params_ptr_->LinkerTypes()[i]].resize(
          DiscretPtr()->NumMyRowElements());

    for (auto const& b_iter : bingids)
    {
      // get neighboring bins
      std::vector<int> nb_binIds;
      nb_binIds.reserve(27);
      // do not check on existence here -> shifted to GetBinContent
      BinStrategyPtr()->GetNeighborAndOwnBinIds(b_iter, nb_binIds);

      for (auto const& nb_iter : nb_binIds)
      {
        // check on existence of bin on this proc
        if (not BinDiscretPtr()->HaveGlobalElement(nb_iter)) continue;

        DRT::Element* currentbin = BinDiscretPtr()->gElement(nb_iter);

        // if a bin has already been examined --> continue with next crosslinker
        if (examinedbins[currentbin->LID()]) continue;
        // else: bin is examined for the first time --> new entry in examinedbins_
        else
          examinedbins[currentbin->LID()] = 1;

        FindPotentialBindingEventsInBinAndNeighborhood(
            currentbin, mybonds, undecidedbonds, intendedbeambonds, false);
      }
    }

    // bind events where myrank only owns the elements, cl are taken care
    // of by their owner (key is clgid)
    std::map<int, Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>> myelebonds;

    // now each row owner of a linker gets requests, makes a random decision and
    // informs back its requesters
    ManageBindingInParallel(mybonds, undecidedbonds, myelebonds);

    // actual update of binding states is done here
    UpdateMyCrosslinkerAndElementBindingStates(mybonds, myelebonds);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UnBindCrosslinker()
{
  CheckInit();

  DRT::Problem::Instance()->Random()->SetRandRange(0.0, 1.0);

  // data containing information about elements that need to be updated on
  // procs != myrank
  std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::UnBindEventData>>> sendunbindevents;
  sendunbindevents.clear();
  // elements that need to be updated on myrank
  std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::UnBindEventData>> myrankunbindevents;
  myrankunbindevents.clear();
  int num_db_dissolved = 0;

  // loop over all row linker (in random order) and dissolve bond if probability
  // criterion is met
  /* note: we loop over all row crosslinker, i.e. myrank needs to update all
   * crosslinker information. As it possible that a row crosslinker is linked
   * to col element, we potentially need to communicate if such an element
   * needs to be updated*/
  const int numrowcl = BinDiscretPtr()->NumMyRowNodes();
  std::vector<int> rorderrowcl = BEAMINTERACTION::UTILS::Permutation(numrowcl);
  for (auto const& rowcli : rorderrowcl)
  {
    CROSSLINKING::CrosslinkerNode* linker =
        dynamic_cast<CROSSLINKING::CrosslinkerNode*>(BinDiscretPtr()->lRowNode(rowcli));

    // only consider unbinding in case off rate is unequal zero
    if (linker->GetMaterial()->KOff() < 1e-08) continue;

    const int clcollid = linker->LID();
    BEAMINTERACTION::DATA::CrosslinkerData* cldata_i = crosslinker_data_[clcollid].get();

    // probability with which a crosslink breaks up in the current time step
    double p_unlink =
        1.0 - exp((-1.0) * crosslinking_params_ptr_->DeltaTime() * linker->GetMaterial()->KOff());

    // different treatment according to number of bonds of a crosslinker
    switch (cldata_i->GetNumberOfBonds())
    {
      case 0:
      {
        // nothing to do here
        break;
      }
      case 1:
      {
        // if probability criterion is not met, we are done here
        if (DRT::Problem::Instance()->Random()->Uni() > p_unlink) break;

        // dissolve bond and update states
        DissolveBond(linker, GetSingleOccupiedClBspot(cldata_i->GetBSpots()), 1, sendunbindevents,
            myrankunbindevents);

        break;
      }
      case 2:
      {
        // calc unbind probability in case of force dependent off rate
        std::vector<double> p_unlink_db(2, 0.0);
        if (abs(linker->GetMaterial()->DeltaBellEq()) > 1.0e-8)
          CalcBellsForceDependentUnbindProbability(
              linker, doublebondcl_[linker->Id()], p_unlink_db);
        else
          p_unlink_db[0] = p_unlink_db[1] = p_unlink;

        // loop through crosslinker bonds in random order
        std::vector<int> ro = BEAMINTERACTION::UTILS::Permutation(cldata_i->GetNumberOfBonds());
        for (auto const& clbspotiter : ro)
        {
          // if probability criterion isn't met, go to next spot
          if (DRT::Problem::Instance()->Random()->Uni() > p_unlink_db[clbspotiter]) continue;

          // dissolve bond and update states
          DissolveBond(linker, clbspotiter, 2, sendunbindevents, myrankunbindevents);
          ++num_db_dissolved;

          // we only want to dissolve one bond per timestep, therefore we go to
          // next crosslinker if we made it so far (i.e. a bond got dissolved)
          break;
        }

        break;
      }
      default:
      {
        dserror("Unrealistic number %i of bonds for a crosslinker.", cldata_i->GetNumberOfBonds());
        exit(EXIT_FAILURE);
      }
    }
  }

  // communicate which elements need to be updated on rank != myrank
  CommunicateCrosslinkerUnbinding(sendunbindevents, myrankunbindevents);

  // update binding status of beam binding partners on myrank
  UpdateBeamBindingStatusAfterUnbinding(myrankunbindevents);

  return num_db_dissolved;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::CalcBellsForceDependentUnbindProbability(
    CROSSLINKING::CrosslinkerNode* linker,
    Teuchos::RCP<BEAMINTERACTION::BeamLink> const& elepairptr,
    std::vector<double>& punlinkforcedependent) const
{
  CheckInitSetup();

  /* characteristic vond length (=nodereldis)
   * from B. Gui and W. Guilford: Mechanics of actomyosin bonds in different nucleotide states are
   * tuned to muscle contraction Fig 2: slip pathway, ADP, delta = 0.0004; Note: delta < 0 -> catch
   * bond, delta > 0 -> bond-weakening see Kai Müller Dis p. 67/68 */
  double const delta = linker->GetMaterial()->DeltaBellEq();
  double const kt = crosslinking_params_ptr_->KT();
  double const koff = linker->GetMaterial()->KOff();
  double const dt = crosslinking_params_ptr_->DeltaTime();

  // safety check
  if (kt < 1e-08)
    dserror(" Thermal energy (KT) set to zero, although you are about to divide by it. ");

  // force and moment exerted on the two binding sites of crosslinker with clgid
  std::vector<LINALG::SerialDenseVector> bspotforce(2, LINALG::SerialDenseVector(6));
  elepairptr->EvaluateForce(bspotforce[0], bspotforce[1]);

  // check if linker is stretched -> sgn+ or compressed -> sgn- by checking orientation of force
  // vector note: this works only if there are no other forces (like inertia, stochastic, damping)
  // acting on the cl
  LINALG::Matrix<3, 1> dist_vec(true);
  LINALG::Matrix<3, 1> bspotforceone(true);
  dist_vec.Update(1.0, elepairptr->GetBindSpotPos1(), -1.0, elepairptr->GetBindSpotPos2());
  for (unsigned int j = 0; j < 3; ++j) bspotforceone(j) = bspotforce[0](j);
  double sgn = (dist_vec.Dot(bspotforceone) < 0.0) ? -1.0 : 1.0;

  /* note: you have a sign criterion that is dependent on the axial strain, but the force you are
   * using also contains shear parts. This means in cases of axial strains near 0 and (large) shear
   * forces you are doing something strange (e.g. jumping between behaviour). Think about a two or
   * three dimensional calculation of the new koff, considering shear and axial forces independently
   */
  //  if ( GState().GetMyRank() == 0 )
  //    std::cout << "Warning: in cases of high shear forces and low axial strain you might "
  //                 "be doing something strange using the force dependent off rate ..." <<
  //                 std::endl;


  // NOTE if you want to add a different force dependent unbinding law: The forces of your linker
  // that were newly set this time step are zero as they are set stress free. Dissolving them here
  // because of a zero force does not make sense. This needs to be considered.


  // calculate new off rate
  std::vector<double> clbspotforcenorm(2, 0.0);
  std::vector<double> forcedependentkoff(2, 0.0);
  for (unsigned int i = 0; i < 2; ++i)
  {
    // currently, only forces (not moments) considered
    bspotforce[i].Reshape(3, 1);
    clbspotforcenorm[i] = bspotforce[i].Norm2();

    // adjusted off-rate according to Bell's equation (Howard, eq 5.10, p.89)
    forcedependentkoff[i] = koff * exp(sgn * clbspotforcenorm[i] * delta / kt);

    // get respective force dependent unbind probability for each cl binding spot
    punlinkforcedependent[i] = 1.0 - exp((-1.0) * dt * forcedependentkoff[i]);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UpdateBeamBindingStatusAfterUnbinding(
    std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::UnBindEventData>> const& unbindevent)
{
  CheckInit();

  // loop through all unbinding events on myrank
  for (auto const& iter : unbindevent)
  {
    // get data
    const int elegidtoupdate = iter->GetEleToUpdate().first;
    const int bspotlocn = iter->GetEleToUpdate().second;
    const int colelelid = Discret().ElementColMap()->LID(elegidtoupdate);

#ifdef DEBUG
    // safety check
    if (Discret().ElementRowMap()->LID(elegidtoupdate) < 0)
      dserror("element with gid %i not owned by proc %i", elegidtoupdate, GState().GetMyRank());
#endif

    // erase current bond
    beam_data_[colelelid]->EraseBondFromBindingSpot(
        iter->GetLinkerType(), bspotlocn, iter->GetClId());
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UpdateMyDoubleBondsAfterRedistribution()
{
  CheckInit();

  std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::BeamLink>>> dbcltosend;
  std::set<int> dbtoerase;

  // loop over all double bonds on myrank
  std::map<int, Teuchos::RCP<BEAMINTERACTION::BeamLink>>::iterator iter;
  for (iter = doublebondcl_.begin(); iter != doublebondcl_.end(); ++iter)
  {
    const int clgid = iter->first;

    // safety check
    if (BinDiscretPtr()->NodeColMap()->LID(clgid) < 0)
      dserror(
          "Crosslinker %i moved further than the bin length in one time step on rank %i, "
          "this is not allowed (maybe increase cutoff radius). ",
          clgid, GState().GetMyRank());

    DRT::Node* doublebondedcl_i = BinDiscretPtr()->gNode(clgid);

    // check ownership
    int owner = doublebondedcl_i->Owner();
    if (owner != GState().GetMyRank())
    {
#ifdef DEBUG
      if (not doublebondcl_.count(clgid))
        dserror("willing to delete double bond %i which is not existing", clgid);
#endif
      dbcltosend[owner].push_back(iter->second);
      dbtoerase.insert(clgid);
    }
  }

  for (auto const& i : dbtoerase) doublebondcl_.erase(i);

  // add new double bonds
  CommunicateBeamLinkAfterRedistribution(dbcltosend);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UpdateMyDoubleBondsRemoteIdList()
{
  CheckInit();

  // loop over all double bonded crosslinker that were read during restart and
  // get the ones that do not belong to myrank
  std::set<int> notonmyrank;
  std::map<int, Teuchos::RCP<BEAMINTERACTION::BeamLink>>::iterator iter;
  for (iter = doublebondcl_.begin(); iter != doublebondcl_.end(); ++iter)
  {
    // double bonded crosslinker gid
    const int clgid = iter->first;

    // not owned
    if (BinDiscretPtr()->NodeRowMap()->LID(clgid) < 0) notonmyrank.insert(clgid);
  }

  int const size = static_cast<int>(notonmyrank.size());
  std::vector<int> unique_clgidlist(notonmyrank.begin(), notonmyrank.end());
  std::vector<int> unique_pidlist(size);

  // find new host procs for double bonded crosslinker by communication
  int err = BinDiscretPtr()->NodeRowMap()->RemoteIDList(
      size, unique_clgidlist.data(), unique_pidlist.data(), NULL);
  if (err < 0) dserror("Epetra_BlockMap::RemoteIDList returned err=%d", err);

  std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::BeamLink>>> dbcltosend;
  for (unsigned int i = 0; i < static_cast<unsigned int>(unique_clgidlist.size()); ++i)
    dbcltosend[unique_pidlist[i]].push_back(doublebondcl_[unique_clgidlist[i]]);

  // update myrank's map
  for (auto const& i : notonmyrank) doublebondcl_.erase(i);

  // send and receive double bonds
  CommunicateBeamLinkAfterRedistribution(dbcltosend);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UnbindCrosslinkerInBinsAndNeighborhood(
    std::set<int> const& bingids)
{
  CheckInit();

  UnbindCrosslinkerInBinsAndNeighborhood(bingids, false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::UnbindCrosslinkerInBinsAndNeighborhood(
    std::set<int> const& bingids, bool doubleunbind)
{
  CheckInitSetup();

  std::set<int> binsonmyrank;
  DetermineResponsilbeProcsForForcedCrosslinkerUnbinding(bingids, binsonmyrank);

  // data containing information about elements that need to be updated on
  // procs != myrank
  std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::UnBindEventData>>> sendunbindevents;
  // elements that need to be updated on myrank
  std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::UnBindEventData>> myrankunbindevents;

  for (auto const& nb_iter : binsonmyrank)
  {
    // get all crosslinker in current bin
    DRT::Element* currentbin = BinDiscretPtr()->gElement(nb_iter);
    DRT::Node** clincurrentbin = currentbin->Nodes();
    const int numcrosslinker = currentbin->NumNode();

    // loop over all crosslinker in current bin
    for (int i = 0; i < numcrosslinker; ++i)
    {
      // get crosslinker in current bin
      DRT::Node* crosslinker_i = clincurrentbin[i];
      BEAMINTERACTION::DATA::CrosslinkerData* cldata_i =
          crosslinker_data_[crosslinker_i->LID()].get();

#ifdef DEBUG
      // safety checks
      if (crosslinker_i->Owner() != GState().GetMyRank())
        dserror(
            " Only row owner of crosslinker changes its state, rank %i is not owner "
            "of linker with gid %i, but rank %i",
            GState().GetMyRank(), crosslinker_i->Id(), crosslinker_i->Owner());
#endif

      switch (cldata_i->GetNumberOfBonds())
      {
        case 0:
        {
          break;
        }
        case 1:
        {
          // dissolve bond and update states
          DissolveBond(crosslinker_i, GetSingleOccupiedClBspot(cldata_i->GetBSpots()),
              cldata_i->GetNumberOfBonds(), sendunbindevents, myrankunbindevents);
          break;
        }
        case 2:
        {
          // dissolve random bond and update states
          DissolveBond(crosslinker_i,
              BEAMINTERACTION::UTILS::Permutation(cldata_i->GetNumberOfBonds())[0],
              cldata_i->GetNumberOfBonds(), sendunbindevents, myrankunbindevents);

          // in case we want to allow transition from double bonded to free, take same linker
          // again
          if (doubleunbind) --i;
          break;
        }
        default:
        {
          dserror(
              " Unrealistic number %i of bonds for a crosslinker.", cldata_i->GetNumberOfBonds());
          exit(EXIT_FAILURE);
        }
      }
    }
  }

  // communicate which elements need to be updated on rank != myrank
  CommunicateCrosslinkerUnbinding(sendunbindevents, myrankunbindevents);

  // update binding status of beam binding partners on myrank
  UpdateBeamBindingStatusAfterUnbinding(myrankunbindevents);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::
    DetermineResponsilbeProcsForForcedCrosslinkerUnbinding(
        std::set<int> const& bingids, std::set<int>& binsonmyrank) const
{
  CheckInitSetup();

  std::set<int> checkedbins;
  std::map<int, std::vector<int>> binstosend;

  // determine all bins that need to dissolve all bonds on myrank and
  // on other procs (these need to be informed by myrank)
  for (auto const& b_iter : bingids)
  {
    std::vector<int> nb_binIds;
    nb_binIds.reserve(27);
    // do not check on existence here
    BinStrategy().GetNeighborAndOwnBinIds(b_iter, nb_binIds);

    for (auto const& nb_iter : nb_binIds)
    {
      // safety check
      if (not BinDiscret().HaveGlobalElement(nb_iter))
        dserror("Not entire neighborhood ghosted, this is a problem in the following ");

      // if a bin has already been examined --> continue with next bin
      // like this we get a unique vector that myrank sends
      if (checkedbins.find(nb_iter) != checkedbins.end()) continue;
      // else: bin is examined for the first time --> new entry in examinedbins_
      else
        checkedbins.insert(nb_iter);

      // decide who needs to dissolve bonds
      const int owner = BinDiscret().gElement(nb_iter)->Owner();
      if (owner == GState().GetMyRank())
        binsonmyrank.insert(nb_iter);
      else
        binstosend[owner].push_back(nb_iter);
    }
  }

  CommunicateBinIds(binstosend, binsonmyrank);
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::CommunicateBinIds(
    std::map<int, std::vector<int>> const& binstosend, std::set<int>& binsonmyrank) const
{
  CheckInitSetup();

  // build exporter
  DRT::Exporter exporter(Discret().Comm());
  int const numproc = Discret().Comm().NumProc();
  int const myrank = GState().GetMyRank();

  // ---- send ---- ( we do not need to pack anything)
  int const length = binstosend.size();
  std::vector<MPI_Request> request(length);
  int tag = 0;
  std::map<int, std::vector<int>>::const_iterator p;
  std::vector<int> targetprocs(numproc, 0);
  for (p = binstosend.begin(); p != binstosend.end(); ++p)
  {
    targetprocs[p->first] = 1;
    exporter.ISend(myrank, p->first, &((p->second)[0]), static_cast<int>((p->second).size()), 1234,
        request[tag]);
    ++tag;
  }
  if (tag != length) dserror("Number of messages is mixed up");

  // ---- prepare receiving procs -----
  std::vector<int> summedtargets(numproc, 0);
  Discret().Comm().SumAll(targetprocs.data(), summedtargets.data(), numproc);

  // ---- receive ----- (we do not need to unpack anything)
  for (int rec = 0; rec < summedtargets[myrank]; ++rec)
  {
    std::vector<int> rdata;
    int length = 0;
    int tag = -1;
    int from = -1;
    exporter.ReceiveAny(from, tag, rdata, length);
    if (tag != 1234) dserror("Received on proc %i data with wrong tag from proc %i", myrank, from);

    // insert in binsonmyrank
    binsonmyrank.insert(rdata.begin(), rdata.end());
  }

  // wait for all communication to finish
  Wait(exporter, request, static_cast<int>(binstosend.size()));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::DissolveBond(DRT::Node* linker,
    int freedbspotid, int numbondsold,
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::UnBindEventData>>>&
        sendunbindevents,
    std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::UnBindEventData>>& myrankunbindevents)
{
  CheckInitSetup();

#ifdef DEBUG
  // safety check
  if (numbondsold < 1) dserror("dissolution of free crosslinker does not make any sense");
#endif

  CROSSLINKING::CrosslinkerNode* crosslinker = dynamic_cast<CROSSLINKING::CrosslinkerNode*>(linker);

  // get linker data
  const int clcollid = crosslinker->LID();
  BEAMINTERACTION::DATA::CrosslinkerData* cldata = crosslinker_data_[clcollid].get();

  // store unbinding event data
  Teuchos::RCP<BEAMINTERACTION::DATA::UnBindEventData> unbindevent =
      Teuchos::rcp(new BEAMINTERACTION::DATA::UnBindEventData());
  unbindevent->SetClId(linker->Id());
  unbindevent->SetEleToUpdate(cldata->GetBSpots()[freedbspotid]);
  unbindevent->SetLinkerType(crosslinker->GetMaterial()->LinkerType());

  // owner of beam
  const int beamowner = DiscretPtr()->gElement(unbindevent->GetEleToUpdate().first)->Owner();

  // check who needs to update the element status
  if (beamowner == GState().GetMyRank())
    myrankunbindevents.push_back(unbindevent);
  else
    sendunbindevents[beamowner].push_back(unbindevent);

  // -----------------------------------------------------------------
  // update crosslinker status
  // -----------------------------------------------------------------
  // update binding status of linker
  cldata->SetBspot(freedbspotid, std::make_pair(-1, -1));

  // update number of bonds
  cldata->SetNumberOfBonds(numbondsold - 1);

  if (numbondsold == 1)
  {
    SetPositionOfNewlyFreeCrosslinker(crosslinker, cldata);
  }
  else if (numbondsold == 2)
  {
    int stayoccpotid = 0;
    if (freedbspotid == 0) stayoccpotid = 1;

    SetPositionOfNewlySingleBondedCrosslinker(crosslinker, cldata, stayoccpotid);

#ifdef DEBUG
    // safety check
    if (not doublebondcl_.count(linker->Id()))
      dserror("crosslinker %i with %i bonds is not in double bonded map of rank %i", linker->Id(),
          cldata->GetNumberOfBonds() + 1, GStatePtr()->GetMyRank());
#endif

    // erase crosslinker from double bonded crosslinker list
    doublebondcl_.erase(linker->Id());
  }
  else
  {
    dserror("dissolution of free linker does not make any sense, something went wrong.");
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::CommunicateUndecidedBonds(
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>>& undecidedbonds,
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>>& requestedcl)
    const
{
  CheckInit();

  // do communication
  std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>> recvbindevent;
  ISendRecvAny(undecidedbonds, recvbindevent);

  for (unsigned int i = 0; i < recvbindevent.size(); ++i)
    requestedcl[recvbindevent[i]->GetClId()].push_back(recvbindevent[i]);

  recvbindevent.clear();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::CommunicateDecidedBonds(
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>>& decidedbonds,
    std::map<int, Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>& myelebonds) const
{
  CheckInit();

  // communicate decided bonds
  std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>> recvbindevent;
  ISendRecvAny(decidedbonds, recvbindevent);

  // loop over received binding events
  for (unsigned int i = 0; i < recvbindevent.size(); ++i)
  {
    // add binding events to new colbond map
    if (recvbindevent[i]->GetPermission())
      myelebonds[recvbindevent[i]->GetClId()] = recvbindevent[i];
  }
  recvbindevent.clear();
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::DecideBindingInParallel(
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>>& requestedcl,
    std::map<int, Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>& mybonds,
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>>& decidedbonds)
    const
{
  CheckInit();

  std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>>::iterator cliter;
  // loop over all requested cl (note myrank is owner of these)
  for (cliter = requestedcl.begin(); cliter != requestedcl.end(); ++cliter)
  {
    // check if myrank wants to bind this crosslinker
    bool myrankbond = false;
    if (mybonds.find(cliter->first) != mybonds.end()) myrankbond = true;

    // ---------------------------------------------------------------------
    // if only one request and myrank does not want to bind this cl,
    // requesting proc gets the permission to do so
    // ---------------------------------------------------------------------
    if (static_cast<int>(cliter->second.size()) == 1 and not myrankbond)
    {
      // we send back the permission to the relevant proc, because myrank as row
      // owner of bspot needs to set the respective stuff for the element of this
      // binding event
      // note: permission = true was send as default, so this can be sent back
      // without changes
      decidedbonds[cliter->second[0]->GetRequestProc()].push_back(cliter->second[0]);

#ifdef DEBUG
      if (cliter->second[0]->GetPermission() != 1)
        dserror(
            " something during communication went wrong, default true permission "
            " not received");
#endif

      // insert this new binding event in map of myrank, because as row owner of
      // this cl he is responsible to set the respective stuff for the crosslinker
      // of this binding event
      mybonds[cliter->first] = cliter->second[0];

      // go to next crosslinker
      continue;
    }

    // ---------------------------------------------------------------------
    // in case number of requesting procs >1 for this cl or myrank wants to
    // set it itself
    // ---------------------------------------------------------------------
    int numrequprocs = static_cast<int>(cliter->second.size());
    if (myrankbond) numrequprocs += 1;

    // get random proc out of affected ones
    DRT::Problem::Instance()->Random()->SetRandRange(0.0, 1.0);
    // fixme: what if random number exactly = 1?
    int rankwithpermission = std::floor(numrequprocs * DRT::Problem::Instance()->Random()->Uni());

    // myrank is allowed to set link
    if (myrankbond and rankwithpermission == (numrequprocs - 1))
    {
      // note: this means link is set between row cl and row ele on myrank,
      // all relevant information for myrank is stored in mybonds
      // loop over all requesters and store their veto
      std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>::iterator iter;
      for (iter = cliter->second.begin(); iter != cliter->second.end(); ++iter)
      {
        (*iter)->SetPermission(0);
        decidedbonds[(*iter)->GetRequestProc()].push_back(*iter);
      }
    }
    // certain requester is allowed to set the link
    else
    {
      // loop over all requesters and store veto for all requester except for one
      std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>::iterator iter;

      int counter = 0;
      for (iter = cliter->second.begin(); iter != cliter->second.end(); ++iter)
      {
        if (rankwithpermission == counter)
        {
          // permission for this random proc
          decidedbonds[(*iter)->GetRequestProc()].push_back(*iter);

#ifdef DEBUG
          if ((*iter)->GetPermission() != 1)
            dserror(
                " something during communication went wrong, default true permission "
                " not received");
#endif

          // erase old binding event
          if (myrankbond) mybonds.erase(cliter->first);

          // insert new binding event
          mybonds[cliter->first] = *iter;
        }
        else
        {
          (*iter)->SetPermission(0);
          decidedbonds[(*iter)->GetRequestProc()].push_back(*iter);
        }
        counter++;
      }
    }
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::CommunicateBeamLinkAfterRedistribution(
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::BeamLink>>>& dbondcltosend)
{
  CheckInit();

  // build exporter
  DRT::Exporter exporter(DiscretPtr()->Comm());
  int const numproc = DiscretPtr()->Comm().NumProc();

  // -----------------------------------------------------------------------
  // send
  // -----------------------------------------------------------------------
  // ---- pack data for sending -----
  std::map<int, std::vector<char>> sdata;
  std::vector<int> targetprocs(numproc, 0);
  std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::BeamLink>>>::const_iterator p;
  for (p = dbondcltosend.begin(); p != dbondcltosend.end(); ++p)
  {
    std::vector<Teuchos::RCP<BEAMINTERACTION::BeamLink>>::const_iterator iter;
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
  for (std::map<int, std::vector<char>>::const_iterator p = sdata.begin(); p != sdata.end(); ++p)
  {
    exporter.ISend(GState().GetMyRank(), p->first, &((p->second)[0]), (int)(p->second).size(), 1234,
        request[tag]);
    ++tag;
  }
  if (tag != length) dserror("Number of messages is mixed up");

  // -----------------------------------------------------------------------
  // receive
  // -----------------------------------------------------------------------
  // ---- prepare receiving procs -----
  std::vector<int> summedtargets(numproc, 0);
  DiscretPtr()->Comm().SumAll(targetprocs.data(), summedtargets.data(), numproc);

  // myrank receive all packs that are sent to him
  for (int rec = 0; rec < summedtargets[GState().GetMyRank()]; ++rec)
  {
    std::vector<char> rdata;
    int length = 0;
    int tag = -1;
    int from = -1;
    exporter.ReceiveAny(from, tag, rdata, length);
    if (tag != 1234)
      dserror("Received on proc %i data with wrong tag from proc %i", GState().GetMyRank(), from);

    // store received data
    std::vector<char>::size_type position = 0;
    while (position < rdata.size())
    {
      std::vector<char> data;
      DRT::ParObject::ExtractfromPack(position, rdata, data);
      // this Teuchos::rcp holds the memory
      Teuchos::RCP<DRT::ParObject> object = Teuchos::rcp(DRT::UTILS::Factory(data), true);
      Teuchos::RCP<BEAMINTERACTION::BeamLink> beamtobeamlink =
          Teuchos::rcp_dynamic_cast<BEAMINTERACTION::BeamLink>(object);
      if (beamtobeamlink == Teuchos::null) dserror("Received object is not a beam to beam linkage");

#ifdef DEBUG
      // some safety checks
      if (BinDiscretPtr()->gNode(beamtobeamlink->Id())->Owner() != GStatePtr()->GetMyRank())
        dserror(
            " A double bond was sent to rank %i, although it is not the owner of "
            "the cl with gid %i ",
            GStatePtr()->GetMyRank(), beamtobeamlink->Id());
      if (doublebondcl_.count(beamtobeamlink->Id()))
        dserror(" Rank %i got sent double bonded crosslinker %i which it already has ",
            GStatePtr()->GetMyRank(), beamtobeamlink->Id());
#endif

      // insert new double bonds in my list
      doublebondcl_[beamtobeamlink->Id()] = beamtobeamlink;
    }

    if (position != rdata.size())
      dserror("Mismatch in size of data %d <-> %d", static_cast<int>(rdata.size()), position);
  }

  // wait for all communication to finish
  Wait(exporter, request, static_cast<int>(sdata.size()));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::CommunicateCrosslinkerUnbinding(
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::UnBindEventData>>>&
        sendunbindevent,
    std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::UnBindEventData>>& myrankunbindevent) const
{
  CheckInit();

  ISendRecvAny(sendunbindevent, myrankunbindevent);
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
template <typename T>
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ISend(DRT::Exporter& exporter,
    std::vector<MPI_Request>& request,
    std::map<int, std::vector<Teuchos::RCP<T>>> const& send) const
{
  CheckInit();

  // ---- pack data for sending -----
  std::map<int, std::vector<char>> sdata;
  typename std::map<int, std::vector<Teuchos::RCP<T>>>::const_iterator p;
  for (p = send.begin(); p != send.end(); ++p)
  {
    typename std::vector<Teuchos::RCP<T>>::const_iterator iter;

    for (iter = p->second.begin(); iter != p->second.end(); ++iter)
    {
      DRT::PackBuffer data;
      (*iter)->Pack(data);
      data.StartPacking();
      (*iter)->Pack(data);
      sdata[p->first].insert(sdata[p->first].end(), data().begin(), data().end());
    }
  }

  // ---- send ----
  const int length = sdata.size();
  request.resize(length);
  int tag = 0;
  for (std::map<int, std::vector<char>>::const_iterator p = sdata.begin(); p != sdata.end(); ++p)
  {
    exporter.ISend(GState().GetMyRank(), p->first, &((p->second)[0]),
        static_cast<int>((p->second).size()), 1234, request[tag]);
    ++tag;
  }
  if (tag != length) dserror("Number of messages is mixed up");
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
template <typename T>
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PrepareReceivingProcs(
    std::map<int, std::vector<Teuchos::RCP<T>>> const& datasenttorank,
    std::vector<int>& summedtargets) const
{
  CheckInit();

  const int numproc = Discret().Comm().NumProc();

  // get number of procs from which myrank receives data
  std::vector<int> targetprocs(numproc, 0);
  typename std::map<int, std::vector<Teuchos::RCP<T>>>::const_iterator prociter;
  for (prociter = datasenttorank.begin(); prociter != datasenttorank.end(); ++prociter)
    targetprocs[prociter->first] = 1;
  // store number of messages myrank receives
  summedtargets.resize(numproc, 0);
  Discret().Comm().SumAll(targetprocs.data(), summedtargets.data(), numproc);
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
template <typename T>
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::RecvAny(
    DRT::Exporter& exporter, int receivesize, std::vector<Teuchos::RCP<T>>& recv) const
{
  CheckInit();

  // myrank receive all packs that are sent to him
  for (int rec = 0; rec < receivesize; ++rec)
  {
    std::vector<char> rdata;
    int length = 0;
    int tag = -1;
    int from = -1;
    exporter.ReceiveAny(from, tag, rdata, length);
    if (tag != 1234)
      dserror("Received on proc %i data with wrong tag from proc %i", GState().GetMyRank(), from);

    // store received data
    std::vector<char>::size_type position = 0;
    while (position < rdata.size())
    {
      std::vector<char> data;

      DRT::ParObject::ExtractfromPack(position, rdata, data);

      Teuchos::RCP<T> data_container =
          Teuchos::rcp(BEAMINTERACTION::DATA::CreateDataContainer<T>(data), true);

      // add received data to list
      recv.push_back(data_container);
    }

    if (position != rdata.size())
      dserror("Mismatch in size of data %d <-> %d", static_cast<int>(rdata.size()), position);
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
template <typename T>
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ISendRecvAny(
    std::map<int, std::vector<Teuchos::RCP<T>>> const& send,
    std::vector<Teuchos::RCP<T>>& recv) const
{
  CheckInit();

  // build exporter
  DRT::Exporter exporter(BinDiscret().Comm());

  // -----------------------------------------------------------------------
  // send
  // -----------------------------------------------------------------------
  // unblocking send
  std::vector<MPI_Request> request;
  ISend(exporter, request, send);

  // -----------------------------------------------------------------------
  // prepare receive
  // -----------------------------------------------------------------------
  std::vector<int> summedtargets;
  PrepareReceivingProcs(send, summedtargets);

  // -----------------------------------------------------------------------
  // receive
  // -----------------------------------------------------------------------
  int receivesize = summedtargets[GState().GetMyRank()];
  RecvAny(exporter, receivesize, recv);

  // wait for all communication to finish
  Wait(exporter, request, static_cast<int>(send.size()));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::Wait(
    DRT::Exporter& exporter, std::vector<MPI_Request>& request, int length) const
{
  CheckInit();

  // wait for all communication to finish
  for (int i = 0; i < length; ++i) exporter.Wait(request[i]);

  // note: if we have done everything correct, this should be a no time operation
  BinDiscret().Comm().Barrier();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PrintAndCheckBindEventData(
    Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData> bindeventdata) const
{
  CheckInit();

  // extract data
  std::cout << "\n Rank: " << GState().GetMyRank() << std::endl;
  std::cout << " crosslinker gid " << bindeventdata->GetClId() << std::endl;
  std::cout << " element gid " << bindeventdata->GetEleId() << std::endl;
  std::cout << " bspot local number " << bindeventdata->GetBSpotLocN() << std::endl;
  std::cout << " requesting proc " << bindeventdata->GetRequestProc() << std::endl;
  std::cout << " permission " << bindeventdata->GetPermission() << std::endl;

  if (bindeventdata->GetClId() < 0 or bindeventdata->GetEleId() < 0 or
      bindeventdata->GetBSpotLocN() < 0 or bindeventdata->GetRequestProc() < 0 or
      not(bindeventdata->GetPermission() == 0 or bindeventdata->GetPermission() == 1))
    dserror(" your bindevent does not make sense.");
}


//-----------------------------------------------------------------------------
// explicit template instantiation (to please every compiler)
//-----------------------------------------------------------------------------
template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ISend(DRT::Exporter&,
    std::vector<MPI_Request>&,
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::CrosslinkerData>>> const&) const;
template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ISend(DRT::Exporter&,
    std::vector<MPI_Request>&,
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BeamData>>> const&) const;
template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ISend(DRT::Exporter&,
    std::vector<MPI_Request>&,
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>> const&) const;
template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ISend(DRT::Exporter&,
    std::vector<MPI_Request>&,
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::UnBindEventData>>> const&) const;


template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PrepareReceivingProcs(
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::CrosslinkerData>>> const&,
    std::vector<int>&) const;
template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PrepareReceivingProcs(
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BeamData>>> const&,
    std::vector<int>&) const;
template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PrepareReceivingProcs(
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>> const&,
    std::vector<int>&) const;
template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::PrepareReceivingProcs(
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::UnBindEventData>>> const&,
    std::vector<int>&) const;


template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::RecvAny(DRT::Exporter&, int const,
    std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::CrosslinkerData>>&) const;
template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::RecvAny(
    DRT::Exporter&, int const, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BeamData>>&) const;
template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::RecvAny(DRT::Exporter&, int const,
    std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>&) const;
template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::RecvAny(DRT::Exporter&, int const,
    std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::UnBindEventData>>&) const;


template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ISendRecvAny(
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::CrosslinkerData>>> const&,
    std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::CrosslinkerData>>&) const;
template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ISendRecvAny(
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BeamData>>> const&,
    std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BeamData>>&) const;
template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ISendRecvAny(
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>> const&,
    std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::BindEventData>>&) const;
template void BEAMINTERACTION::SUBMODELEVALUATOR::Crosslinking::ISendRecvAny(
    std::map<int, std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::UnBindEventData>>> const&,
    std::vector<Teuchos::RCP<BEAMINTERACTION::DATA::UnBindEventData>>&) const;

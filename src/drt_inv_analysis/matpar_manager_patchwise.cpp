/*----------------------------------------------------------------------*/
/*! \file
\brief Creating patches from an elementwise layout

\level 3

\maintainer Sebastian Brandstaeter
!*/

/*----------------------------------------------------------------------*/
/* headers */
#include "matpar_manager_patchwise.H"

#include "invana_utils.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../linalg/linalg_utils_densematrix_communication.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_element.H"

#include "../drt_mat/material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"

#include <queue>

typedef std::map<int, std::vector<int>> PATCHES;

/*----------------------------------------------------------------------*/
INVANA::MatParManagerPerPatch::MatParManagerPerPatch(Teuchos::RCP<DRT::Discretization> discret)
    : MatParManagerPerElement(discret),
      qthresh_(0.1),
      map_restart_file_("none"),
      map_restart_step_(-1),
      max_num_levels_(1)
{
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerPatch::Setup()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "-----------------------------" << std::endl;
    std::cout << "MatParManagerPatch Setup:" << std::endl;
  }

  const Teuchos::ParameterList& invp = Inpar();
  map_restart_step_ = invp.get<int>("MAP_REDUCT_RESTART");
  map_restart_file_ = Teuchos::getNumericStringParameter(invp, "MAP_REDUCT_RESTARTFILE");

  max_num_levels_ = invp.get<int>("NUM_REDUCT_LEVELS");

  if (max_num_levels_ < 1) dserror("Choose at least NUM_LEVELS = 1 for the patch creation!");

  // call setup of the Base class to have all the
  // layout of the elementwise distribution
  MatParManagerPerElement::Setup();
  optparams_elewise_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_, 1, true));
  elewise_map_ = Teuchos::rcp(new Epetra_Map(*paramlayoutmap_));

  // read map approximation to perform the reduction on
  ReadMAPApproximation();

  // sort elementwise solution into histogram
  MakeHistogram();

  // create sparse approximation of the MAP solution
  CreateProjection();

  // initialize parameters
  InitParameters();

  // Some user information
  if (Comm().MyPID() == 0) std::cout << std::endl;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerPatch::FillParameters(Teuchos::RCP<Epetra_MultiVector> params)
{
  params->PutScalar(0.0);

  // Inject into the elementwise solution space
  int err = projector_->Multiply(true, *optparams_, *optparams_elewise_);
  if (err != 0) dserror("Application of prolongator failed.");

  // loop the parameter blocks
  for (int k = 0; k < paramapextractor_->NumMaps(); k++)
  {
    Teuchos::RCP<Epetra_Vector> tmp =
        paramapextractor_->ExtractVector(*(*optparams_elewise_)(0), k);
    for (int i = 0; i < tmp->MyLength(); i++)
    {
      int pgid =
          tmp->Map().GID(i);  // !! the local id of the partial map is not the local parameter id!!
      int plid = paramapextractor_->FullMap()->LID(pgid);
      params->ReplaceGlobalValue(ParamsLIDtoeleGID()[plid], k, (*tmp)[i]);
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerPatch::ApplyParametrization(
    DcsMatrix& matrix, Teuchos::RCP<Epetra_MultiVector> diagonals)
{
  // this is ok here since we have a sparse approximation
  Teuchos::RCP<Epetra_CrsMatrix> fullmatrix = matrix.FillMatrix();

  // todo: this is not ok! loop over the single columns of matrix
  // and extract only the diagonal component.
  // matrix * projector_
  Teuchos::RCP<Epetra_CrsMatrix> mr = LINALG::Multiply(fullmatrix, false, projector_, false);
  // projector'*matrix*projector
  Teuchos::RCP<Epetra_CrsMatrix> pmr = LINALG::Multiply(projector_, true, mr, false);

  Epetra_Vector diagonal(pmr->RowMap(), true);
  pmr->ExtractDiagonalCopy(diagonal);

  // loop the parameter blocks
  for (int k = 0; k < paramapextractor_->NumMaps(); k++)
  {
    Teuchos::RCP<Epetra_Vector> tmp = paramapextractor_->ExtractVector(diagonal, k);
    for (int i = 0; i < tmp->MyLength(); i++)
    {
      int pgid =
          tmp->Map().GID(i);  // !! the local id of the partial map is not the local parameter id!!
      int plid = paramapextractor_->FullMap()->LID(pgid);
      diagonals->ReplaceGlobalValue(ParamsLIDtoeleGID()[plid], k, (*tmp)[i]);
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerPatch::InitParameters()
{
  // sanity checks
  if (not projector_->DomainMap().PointSameAs(optparams_elewise_->Map()))
    dserror("Restrictor->DomainMap error.");

  if (not projector_->RangeMap().PointSameAs(optparams_->Map()))
    dserror("Restrictor->RangeMap error");

  // parameters are not initialized from input but
  // from the elementwise layout
  int err = projector_->Multiply(false, *optparams_elewise_, *optparams_);
  if (err != 0) dserror("Application of restrictor failed.");

  // set initial values
  optparams_initial_->Scale(1.0, *optparams_);
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerPatch::ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint,
    double val, int elepos, int paraposglobal, int paraposlocal)
{
  if (EleGIDtoparamsLID().find(elepos) == EleGIDtoparamsLID().end())
    dserror("proc %d, ele %d not in this map", Discret()->Comm().MyPID(), elepos);

  // parameter in the 'elementwise' (there can be multiple
  // per element) parameter layout
  int plid = EleGIDtoparamsLID()[elepos].at(paraposlocal);

  // check to which patch it belongs
  int patchid = pidtopatch_[plid];

  int success = dfint->SumIntoMyValue(patchid, 0, val);
  if (success != 0) dserror("gid %d is not on this processor", plid);
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerPatch::Finalize(
    Teuchos::RCP<Epetra_MultiVector> source, Teuchos::RCP<Epetra_MultiVector> target)
{
  // sum across processor
  std::vector<double> val(source->MyLength(), 0.0);
  Discret()->Comm().SumAll((*source)(0)->Values(), &val[0], source->MyLength());

  for (int i = 0; i < target->MyLength(); i++) target->SumIntoGlobalValue(i, 0, val[i]);

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerPatch::ReadMAPApproximation()
{
  // Create the input control file object
  Teuchos::RCP<IO::InputControl> input =
      Teuchos::rcp(new IO::InputControl(map_restart_file_, Discret()->Comm()));

  // and the discretization reader to read from the input file
  IO::DiscretizationReader reader(Discret(), input, map_restart_step_);

  if (not Discret()->Comm().MyPID())
  {
    std::cout << "  Reading MAP approximation: ";
    std::cout << "  step " << map_restart_step_ << " (from: " << input->FileName() << ")"
              << std::endl;
  }

  // MAP solution
  reader.ReadMultiVector(optparams_elewise_, "solution");

  // Read lbfgs matrix storage
  Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector>> sstore =
      Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, elewise_map_.get(), true));
  Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector>> ystore =
      Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, elewise_map_.get(), true));

  int actsize = reader.ReadInt("storage_size");

  // initialize storage
  sstore->Resize(-actsize + 1, 0, elewise_map_.get(), true);
  ystore->Resize(-actsize + 1, 0, elewise_map_.get(), true);

  Teuchos::RCP<Epetra_MultiVector> storage =
      Teuchos::rcp(new Epetra_MultiVector(*elewise_map_, actsize, false));

  reader.ReadMultiVector(storage, "sstore");
  for (int i = 0; i < actsize; i++) sstore->UpdateSteps(*(*storage)(i));

  storage->Scale(0.0);
  reader.ReadMultiVector(storage, "ystore");
  for (int i = 0; i < actsize; i++) ystore->UpdateSteps(*(*storage)(i));

  // ---- create covariance matrix
  double scalefac = 1.0;
  bool objfuncscal = DRT::INPUT::IntegralValue<bool>(Inpar(), "OBJECTIVEFUNCTSCAL");
  if (objfuncscal) scalefac = Objfunct().GetScaleFac();

  bool initscal = DRT::INPUT::IntegralValue<bool>(Inpar(), "LBFGSINITSCAL");

  fullcovariance_ = Teuchos::rcp(new DcsMatrix(sstore, ystore, initscal, objfuncscal, scalefac));

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerPatch::CreateProjection()
{
  // create the orthogonal dictionary for each level
  // and check the approximation power
  double quality = 100.0;
  int level = 1;
  while (quality > qthresh_ and level <= max_num_levels_)
  {
    CreateLevelDictionary(level);
    quality = CheckApproximation();
    level += 1;
  }

  // Some user information
  if (Comm().MyPID() == 0)
  {
    std::cout << "  Reached approximation quality of " << quality << std::endl;
    std::cout << "  using the first " << level - 1 << " maxima of the solution histogram"
              << std::endl;
  }

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerPatch::CreateLevelDictionary(int patchlevel)
{
  // for convenience
  Epetra_Vector* optparams = (*optparams_elewise_)(0);

  const Epetra_Map& elemap = graph_->RowMap();

  if (not elemap.SameAs(optparams->Map())) dserror("Graph and solution map don't match. Fatal!");

  // -------- define levels for patching according to the histogram
  // the values of the patches
  std::vector<double> patchvalues;
  // copy the histogram for manipulation
  std::map<int, int> histbins = histbins_;
  for (int i = 0; i < patchlevel; i++)
  {
    // find bin with maximum entries
    int imax;
    FindHistogramMax(histbins, imax);

    // get mid point of this bin as patchvalue
    double val = (histvalues_[imax] - histvalues_[imax + 1]) / 2.0 + histvalues_[imax];
    patchvalues.push_back(val);

    // remove this bin from the histogram
    histbins.erase(imax);
  }
  // --------

  // -------- clustering the map solution according to the patchlevels
  PATCHES levelgids;
  for (int i = 0; i < optparams->MyLength(); i++)
  {
    int patchval = 0;
    double currval = (*optparams)[i];
    double dist = 1.0e10;
    for (int j = 0; j < (int)patchvalues.size(); j++)
    {
      if (abs(currval - patchvalues[j]) < dist)
      {
        patchval = j;
        dist = abs(currval - patchvalues[j]);
      }
    }
    levelgids[patchval].push_back(optparams->Map().GID(i));
  }

  // keep the different levels separated in a mapextractor
  std::vector<Teuchos::RCP<const Epetra_Map>> partials;
  for (int i = 0; i < patchlevel; i++)
  {
    // this will create empty maps for levels not having any elements
    partials.push_back(Teuchos::rcp(
        new Epetra_Map(-1, levelgids[i].size(), levelgids[i].data(), 0, optparams->Comm())));
  }
  Teuchos::RCP<LINALG::MultiMapExtractor> levelmap =
      Teuchos::rcp(new LINALG::MultiMapExtractor(elemap, partials));
  // --------

  // -------- check connectivity within levels
  // loop all the different levels and check for their connectivity
  PATCHES patches;
  int patchcount = 0;
  for (int i = 0; i < levelmap->NumMaps(); i++)
  {
    // map for this level
    const Teuchos::RCP<const Epetra_Map> level = levelmap->Map(i);

    if (level->NumGlobalElements() == 0) continue;

    // allreduced version to check existence of
    // neighbours across processors in this level
    Teuchos::RCP<Epetra_Map> levelallred = LINALG::AllreduceEMap(*level);

    // collect neighbouring information for this level
    PATCHES neighbours;
    for (int i = 0; i < level->NumMyElements(); i++)
    {
      int egid = level->GID(i);
      int lid = optparams->Map().LID(egid);

      int numcols;
      double* values;
      int* indices;
      graph_->ExtractMyRowView(lid, numcols, values, indices);

      std::vector<int> myneighbours;
      for (int j = 0; j < numcols; j++)
      {
        int ngid = graph_->ColMap().GID(indices[j]);
        int nlid = levelallred->LID(ngid);
        // dont make itself and off-level elements a neighbour
        if (ngid != egid and nlid != -1) myneighbours.push_back(ngid);
      }
      neighbours[egid] = myneighbours;
    }

    // bring neighbour information to proc zero
    Teuchos::RCP<Epetra_Map> tomap = LINALG::AllreduceEMap(*level, 0);
    DRT::Exporter expo(*level, *tomap, level->Comm());
    expo.Export(neighbours);


    // let only proc 0 find the relevant connectivities
    if (level->Comm().MyPID() == 0)
    {
      PATCHES lpatches;
      FindLevelConnectivity(neighbours, lpatches);

      PATCHES::iterator it;
      for (it = lpatches.begin(); it != lpatches.end(); it++)
      {
        patches[patchcount] = it->second;
        patchcount++;
      }
    }
    level->Comm().Barrier();
  }
  // --------

  // build paramlayout maps anew
  paramlayoutmapunique_ = Teuchos::rcp(new Epetra_Map(-1, patchcount, 0, optparams->Map().Comm()));
  paramlayoutmap_ = LINALG::AllreduceEMap(*paramlayoutmapunique_);

  // export patch information to all procs
  DRT::Exporter expo(*paramlayoutmapunique_, *paramlayoutmap_, optparams->Map().Comm());
  expo.Export(patches);

  // fill patchmap
  std::vector<Teuchos::RCP<const Epetra_Map>> patchmaps;

  PATCHES::iterator it;
  for (it = patches.begin(); it != patches.end(); it++)
  {
    std::vector<int> mygids;
    std::vector<int>::iterator jt;
    for (jt = it->second.begin(); jt != it->second.end(); jt++)
    {
      int lid = optparams->Map().LID(*jt);

      if (lid == -1) continue;

      mygids.push_back(*jt);
      pidtopatch_[lid] = it->first;
    }

    patchmaps.push_back(
        Teuchos::rcp(new Epetra_Map(-1, mygids.size(), mygids.data(), 0, optparams->Comm())));
  }
  patchmap_ = Teuchos::rcp(new LINALG::MultiMapExtractor(elemap, patchmaps));

  // build restrictor
  int maxbw = 0;
  for (it = patches.begin(); it != patches.end(); it++)
  {
    if (maxbw < (int)it->second.size()) maxbw = it->second.size();
  }
  Teuchos::RCP<Epetra_Map> colmap = LINALG::AllreduceEMap(*patchmap_->FullMap(), 0);
  projector_ =
      Teuchos::rcp(new Epetra_CrsMatrix(Copy, *paramlayoutmapunique_, *colmap, maxbw, false));

  for (int i = 0; i < projector_->NumMyRows(); i++)
  {
    int numentries = patches[i].size();
    std::vector<double> values(numentries, 1.0 / sqrt(numentries));

    int err = projector_->InsertGlobalValues(i, numentries, &values[0], patches[i].data());
    if (err < 0) dserror("Restrictor/Prolongator insertion failed.");
  }
  int err = projector_->FillComplete(*patchmap_->FullMap(), *paramlayoutmapunique_, true);
  if (err != 0) dserror("Restrictor/Prolongator FillComplete failed.");

  optparams_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmapunique_, 1, true));
  optparams_initial_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmapunique_, 1, true));

  return;
}

/*----------------------------------------------------------------------*/
double INVANA::MatParManagerPerPatch::CheckApproximation()
{
  // compute 'optimal' optimization parameters
  int err = projector_->Multiply(false, *optparams_elewise_, *optparams_);
  if (err != 0) dserror("Application of restrictor failed.");

  // project to elementwise solution space
  Teuchos::RCP<Epetra_MultiVector> projection =
      Teuchos::rcp(new Epetra_MultiVector(*elewise_map_, 1, false));
  err = projector_->Multiply(true, *optparams_, *projection);
  if (err != 0) dserror("Application of restrictor failed.");

  // compute metric
  projection->Update(-1.0, *optparams_elewise_, 1.0);
  // projection->Print(std::cout);
  double metric;
  projection->Norm2(&metric);

  double base;
  optparams_elewise_->Norm2(&base);

  // make metric relative
  metric /= base;

  return metric;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerPatch::FindLevelConnectivity(PATCHES& neighbourhood, PATCHES& patches)
{
  // initialize distances to -1
  std::map<int, int> distance;
  PATCHES::iterator it;
  for (it = neighbourhood.begin(); it != neighbourhood.end(); it++) distance[it->first] = -1;

  int distsum = -1;
  int iter = 0;
  while (distsum != (int)neighbourhood.size())
  {
    std::queue<int> nodes;
    std::vector<int> patchvec;

    // init queue with the first node that has not a distance yet
    std::map<int, int>::iterator kt;
    for (kt = distance.begin(); kt != distance.end(); kt++)
    {
      if (kt->second == -1)
      {
        nodes.push(kt->first);
        patchvec.push_back(kt->first);
        kt->second = 1;
        break;
      }
    }

    while (not nodes.empty())
    {
      int current = nodes.front();
      nodes.pop();
      std::vector<int>::iterator it;
      for (it = neighbourhood[current].begin(); it != neighbourhood[current].end(); it++)
      {
        if (distance[*it] == -1)
        {
          distance[*it] = 1;
          nodes.push(*it);
          patchvec.push_back(*it);
        }
      }
    }
    patches[iter] = patchvec;

    iter++;
    distsum = 0;
    std::map<int, int>::iterator it;
    for (it = distance.begin(); it != distance.end(); it++) distsum += it->second;
  }

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerPatch::MakeHistogram()
{
  // communicate solution to all procs
  Teuchos::RCP<Epetra_Map> alllocal = LINALG::AllreduceEMap(*elewise_map_);
  Epetra_Vector data(*alllocal, true);

  // bring to every proc the same data
  Epetra_Import importer(*alllocal, *elewise_map_);
  int err = data.Import(*(*optparams_elewise_)(0), importer, Insert);
  if (err) dserror("import of data to every processor failed with code: %d", err);

  // compute histogram bin values
  histvalues_.clear();
  int nbins = 10;
  double min = 0.0;
  data.MinValue(&min);
  double max = 0.0;
  data.MaxValue(&max);
  histvalues_.push_back(min);
  double dd = (max - min) / nbins;
  for (int i = 1; i <= nbins; i++) histvalues_.push_back(min + i * dd);

  // init bins
  histbins_.clear();
  for (int i = 0; i < nbins; i++) histbins_[i] = 0;

  // sort data
  for (int j = 0; j < data.MyLength(); j++)
  {
    for (int i = 0; i < nbins; i++)
    {
      if ((data[j] >= histvalues_[i]) && (data[j] < histvalues_[i + 1])) histbins_[i] += 1;
    }
  }

#if 0
  // print histogram
  std::cout << "histogram: ";
  for (int i=0; i<nbins; i++)
    std::cout << histbins_[i] << ", ";
  std::cout << std::endl;
  for (int i=0; i<nbins; i++)
    std::cout << histvalues_[i] << " " << histvalues_[i+1] << ", ";
  std::cout << std::endl;
#endif

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerPatch::FindHistogramMax(std::map<int, int>& histbins, int& index)
{
  int max = -1;

  std::map<int, int>::iterator it;
  for (it = histbins.begin(); it != histbins.end(); it++)
  {
    if (it->second > max)
    {
      max = it->second;
      index = it->first;
    }
  }


  return;
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix> INVANA::MatParManagerPerPatch::InitialCovariance()
{
  // projector * covariance
  Teuchos::RCP<Epetra_CrsMatrix> interm = Teuchos::rcp(new Epetra_CrsMatrix(
      Copy, projector_->RowMap(), projector_->ColMap(), projector_->MaxNumEntries(), false));

  Teuchos::RCP<Epetra_Vector> column =
      Teuchos::rcp(new Epetra_Vector(fullcovariance_->RowMap(), false));
  Teuchos::RCP<Epetra_Vector> col_interm =
      Teuchos::rcp(new Epetra_Vector(projector_->RowMap(), false));
  int maxnumentries = fullcovariance_->MaxNumEntries();
  for (int i = 0; i < maxnumentries; i++)
  {
    // get column
    fullcovariance_->ExtractGlobalColumnCopy(i, *column);

    // project this column
    int err = projector_->Multiply(false, *column, *col_interm);
    if (err != 0) dserror("Projection failed.");

    // put column into intermediate matrix
    double* vals;
    int colind = i;
    col_interm->ExtractView(&vals);
    for (int j = 0; j < col_interm->MyLength(); j++)
    {
      int gid = projector_->RowMap().GID(j);
      interm->InsertGlobalValues(gid, 1, &vals[j], &colind);
    }
  }
  interm->FillComplete(projector_->ColMap(), projector_->RangeMap());

  // interm * projector'
  Teuchos::RCP<Epetra_CrsMatrix> cov = LINALG::Multiply(*interm, false, *projector_, true);
  return cov;
}

/*----------------------------------------------------------------------*/
/*!
\file matpar_manager_patchwise.cpp
\brief Creating patches from an elementwise layout
<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
!*/

/*----------------------------------------------------------------------*/
/* headers */
#include "matpar_manager_patchwise.H"

#include "invana_utils.H"
#include "../linalg/linalg_utils.H"
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

typedef std::map<int, std::vector<int> > PATCHES;

/*----------------------------------------------------------------------*/
INVANA::MatParManagerPerPatch::MatParManagerPerPatch(Teuchos::RCP<DRT::Discretization> discret)
   :MatParManagerPerElement(discret),
map_restart_file_("none"),
map_restart_step_(-1),
num_levels_(1)
{}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerPatch::Setup()
{

  const Teuchos::ParameterList& invp = DRT::Problem::Instance()->StatInverseAnalysisParams();
  map_restart_step_ = invp.get<int>("MAP_RESTART");
  map_restart_file_ = Teuchos::getNumericStringParameter(invp,"MAP_RESTARTFILE");

  num_levels_ = invp.get<int>("NUM_PATCH_LEVELS");

  if (num_levels_<=1)
    dserror("Choose at least NUM_PATCH_LEVELS = 2 for the patch creation!");

  // call setup of the Base class to have all the
  // layout of the elementwise distribution
  MatParManagerPerElement::Setup();
  optparams_elewise_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_,1,true));
  elewise_map_ = Teuchos::rcp(new Epetra_Map(*paramlayoutmap_));

  // read map approximation to perform the reduction on
  ReadMAPApproximation();

  // reduce set of parameters
  ReduceBasis();

  // initialize parameters
  InitParameters();
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerPatch::FillParameters(Teuchos::RCP<Epetra_MultiVector> params)
{
  params->PutScalar(0.0);

  // loop the parameter blocks
  for (int k=0; k<paramapextractor_->NumMaps(); k++)
  {
    Teuchos::RCP<Epetra_Vector> optparams = paramapextractor_->ExtractVector(*(*optparams_)(0),k);

    // make optparams redundant on every proc
    Teuchos::RCP<Epetra_Vector> optparams_allred = Teuchos::rcp(new Epetra_Vector(*paramlayoutmap_,false));
    Epetra_Import importer(optparams_allred->Map(), optparams_->Map());
    int err = optparams_allred->Import(*optparams_, importer, Insert);
    if (err)
      dserror("Export using exporter returned err=%d", err);

    for (int i=0; i<patchmap_->NumMaps(); i++)
    {
      double val = (*optparams_allred)[i];

      const Teuchos::RCP<const Epetra_Map> patch = patchmap_->Map(i);
      for (int j=0; j<patch->NumMyElements(); j++)
      {
        int pgid = patch->GID(j);
        int plid = patchmap_->FullMap()->LID(pgid);
        params->ReplaceGlobalValue(ParamsLIDtoeleGID()[plid],k,val);
      }
    }
  }
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerPatch::InitParameters()
{
  // sanity checks
  if ( not restrictor_->DomainMap().PointSameAs(optparams_elewise_->Map()))
    dserror("Restrictor->DomainMap error.");

  if ( not restrictor_->RangeMap().PointSameAs(optparams_->Map()))
    dserror("Restrictor->RangeMap error");

  // parameters are not initialized from input but
  // from the elementwise layout
  int err = restrictor_->Multiply(false,*optparams_elewise_,*optparams_);
  if (err!=0)
    dserror("Application of restrictor failed.");

  // set initial values
  optparams_initial_->Scale(1.0,*optparams_);

}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerPatch::ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint,
                                                            double val,
                                                            int elepos,
                                                            int paraposglobal,
                                                            int paraposlocal)
{
  if (EleGIDtoparamsLID().find(elepos) == EleGIDtoparamsLID().end())
    dserror("proc %d, ele %d not in this map", Discret()->Comm().MyPID(), elepos);

  // parameter in the 'elementwise' (there can be multiple
  // per element) parameter layout
  int plid = EleGIDtoparamsLID()[elepos].at(paraposlocal);

  // check to which patch it belongs
  int patchid = pidtopatch_[plid];

  int success = dfint->SumIntoMyValue(patchid,0,val);
  if (success!=0) dserror("gid %d is not on this processor", plid);
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerPatch::Finalize(Teuchos::RCP<Epetra_MultiVector> source,
    Teuchos::RCP<Epetra_MultiVector> target)
{
  // sum across processor
  std::vector<double> val(source->MyLength(),0.0);
  Discret()->Comm().SumAll((*source)(0)->Values(),&val[0],source->MyLength());

  for (int i=0; i<target->MyLength(); i++)
    target->SumIntoGlobalValue(i,0,val[i]);

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerPatch::ReadMAPApproximation()
{
  // Create the input control file object
  Teuchos::RCP<IO::InputControl> input = Teuchos::rcp(new
      IO::InputControl(map_restart_file_, Discret()->Comm()));

  // and the discretization reader to read from the input file
  IO::DiscretizationReader reader(Discret(),input,map_restart_step_);

  if (not Discret()->Comm().MyPID())
  {
    std::cout << "Reading MAP approximation for parameter basis reduction ";
    std::cout << "step " << map_restart_step_ << " from file: " << input->FileName() << std::endl;
  }

  reader.ReadMultiVector(optparams_elewise_,"solution");


  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerPatch::ReduceBasis()
{
  // for convenience
  Teuchos::RCP<Epetra_CrsMatrix> adjacency = GetConnectivityData()->AdjacencyMatrix();
  Epetra_Vector* optparams = (*optparams_elewise_)(0);

  const Epetra_Map& elemap = adjacency->RowMap();

  if (not elemap.SameAs(optparams->Map()))
    dserror("Graph and solution map don't match. Fatal!");

  // read from inputfile
  int patchlevels = num_levels_;

  // -------- define levels for patching
  double patchmin;
  double patchmax;
  optparams->MinValue(&patchmin);
  optparams->MaxValue(&patchmax);

  std::vector<double> patchvalues;
  double incr = (patchmax-patchmin)/(patchlevels-1);
  for (int i=0; i<patchlevels; i++)
    patchvalues.push_back(patchmin+(i*incr));
  // --------

  // -------- clustering the map solution according to the patchlevels
  PATCHES levelgids;
  for (int i=0; i<optparams->MyLength(); i++)
  {
    int patchval = 0;
    double currval = (*optparams)[i];
    for (int j=1; j<(int)patchvalues.size(); j++)
    {
      if (abs(currval-patchvalues[j]) < abs(currval-patchvalues[patchval]))
        patchval = j;
    }
    levelgids[patchval].push_back(optparams->Map().GID(i));
  }

  // keep the different levels separated in a mapextractor
  std::vector< Teuchos::RCP<const Epetra_Map> > partials;
  for (int i=0; i<patchlevels; i++)
  {
    // this will create empty maps for levels not having any elements
    partials.push_back(Teuchos::rcp(new
        Epetra_Map(-1,levelgids[i].size(),levelgids[i].data(),0,optparams->Comm())));
  }
  Teuchos::RCP<LINALG::MultiMapExtractor> levelmap =
      Teuchos::rcp(new LINALG::MultiMapExtractor(elemap,partials));
  // --------

  // -------- check connectivity within levels
  // loop all the different levels and check for their connectivity
  PATCHES patches;
  int patchcount = 0;
  for (int i=0; i<levelmap->NumMaps(); i++)
  {
    // map for this level
    const Teuchos::RCP<const Epetra_Map> level = levelmap->Map(i);

    if (level->NumGlobalElements()==0)
      continue;

    // allreduced version to check existence of
    // neighbours across processors in this level
    Teuchos::RCP<Epetra_Map> levelallred = LINALG::AllreduceEMap(*level);

    // collect neighbouring information for this level
    PATCHES neighbours;
    for (int i=0; i<level->NumMyElements(); i++)
    {
      int egid = level->GID(i);
      int lid = optparams->Map().LID(egid);

      int numcols;
      double* values;
      int* indices;
      adjacency->ExtractMyRowView(lid,numcols,values,indices);

      std::vector<int> myneighbours;
      for (int j=0; j<numcols; j++)
      {
        int ngid = adjacency->ColMap().GID(indices[j]);
        int nlid = levelallred->LID(ngid);
        // dont make itself and off-level elements a neighbour
        if (ngid != egid and nlid != -1)
          myneighbours.push_back(ngid);
      }
      neighbours[egid] = myneighbours;
    }

    // bring neighbour information to proc zero
    Teuchos::RCP<Epetra_Map> tomap = LINALG::AllreduceEMap(*level,0);
    DRT::Exporter expo(*level,*tomap,level->Comm());
    expo.Export(neighbours);


    // let only proc 0 find the relevant connectivities
    if (level->Comm().MyPID() == 0)
    {
      PATCHES lpatches;
      FindLevelConnectivity(neighbours,lpatches);

      PATCHES::iterator it;
      for (it=lpatches.begin(); it!=lpatches.end(); it++)
      {
        patches[patchcount] = it->second;
        patchcount++;
      }
    }
    level->Comm().Barrier();
  }
  // --------

  // build paramlayout maps anew
  paramlayoutmapunique_ = Teuchos::rcp(new Epetra_Map(-1,patchcount,0,optparams->Map().Comm()));
  paramlayoutmap_ = LINALG::AllreduceEMap(*paramlayoutmapunique_);

  // export patch information to all procs
  DRT::Exporter expo(*paramlayoutmapunique_,*paramlayoutmap_,optparams->Map().Comm());
  expo.Export(patches);

  // fill patchmap
  std::vector< Teuchos::RCP<const Epetra_Map> > patchmaps;

  PATCHES::iterator it;
  for (it=patches.begin(); it!=patches.end(); it++)
  {
    std::vector<int> mygids;
    std::vector<int>::iterator jt;
    for (jt=it->second.begin(); jt!=it->second.end(); jt++)
    {
      int lid = optparams->Map().LID(*jt);

      if (lid == -1)
        continue;

      mygids.push_back(*jt);
      pidtopatch_[lid] = it->first;
    }

    patchmaps.push_back(Teuchos::rcp(new
            Epetra_Map(-1,mygids.size(),mygids.data(),0,optparams->Comm())));
  }
  patchmap_ = Teuchos::rcp(new LINALG::MultiMapExtractor(elemap,patchmaps));

  // build restrictor
  int maxbw=0;
  for (it=patches.begin(); it!=patches.end(); it++)
  {
    if (maxbw < (int)it->second.size())
      maxbw = it->second.size();
  }
  Teuchos::RCP<Epetra_Map> colmap = LINALG::AllreduceEMap(*patchmap_->FullMap(),0);
  restrictor_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*paramlayoutmapunique_,*colmap,maxbw,false));

  for (int i=0; i<restrictor_->NumMyRows(); i++)
  {
    int numentries = patches[i].size();
    std::vector<double> values(numentries, 1.0/numentries);

    int err = restrictor_->InsertGlobalValues(i,numentries,&values[0],patches[i].data());
    if (err < 0)
      dserror("Restrictor insertion failed.");
  }
  int err = restrictor_->FillComplete(*patchmap_->FullMap(), *paramlayoutmapunique_,true);
  if (err != 0)
    dserror("Restrictor FillComplete failed.");


  optparams_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmapunique_,1,true));
  optparams_initial_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmapunique_,1,true));

  // todo: extend all the graph based stuff to multiple physical parameters
  std::vector< Teuchos::RCP<const Epetra_Map> > tmp;
  tmp.push_back(paramlayoutmapunique_);
  paramapextractor_ = Teuchos::rcp(new
      LINALG::MultiMapExtractor(*paramlayoutmapunique_,tmp));


  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerPatch::FindLevelConnectivity(
    PATCHES& neighbourhood, PATCHES& patches)
{
  // initialize distances to -1
  std::map<int, int> distance;
  PATCHES::iterator it;
  for (it=neighbourhood.begin(); it!=neighbourhood.end(); it++)
    distance[it->first] = -1;

  int distsum = -1;
  int iter = 0;
  while (distsum != (int)neighbourhood.size())
  {

    std::queue<int> nodes;
    std::vector<int> patchvec;

    // init queue with the first node that has not a distance yet
    std::map<int, int>::iterator kt;
    for (kt=distance.begin(); kt!=distance.end(); kt++)
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
      for (it=neighbourhood[current].begin(); it!=neighbourhood[current].end(); it++)
      {
        if (distance[*it]== -1)
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
    for (it=distance.begin(); it!=distance.end(); it++)
      distsum += it->second;

  }

  return;
}

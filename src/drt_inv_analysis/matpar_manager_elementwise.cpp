/*----------------------------------------------------------------------*/
/*!
\file matpar_manager_elementwise.cpp
\brief Class to handle calls to material parameters from an optimization routine

<pre>
\level 3
\maintainer Sebastian Brandstaeter
            brandstaeter@lnm.mw.tum.de
            089 - 289-15276
</pre>

!*/

/*----------------------------------------------------------------------*/
/* headers */
#include "matpar_manager_elementwise.H"

#include "invana_utils.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_element.H"

#include "../drt_mat/material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"

INVANA::MatParManagerPerElement::MatParManagerPerElement(Teuchos::RCP<DRT::Discretization> discret)
    : MatParManager(discret), has_graph_(false)
{
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerElement::Setup()
{
  // temp map to keep correspondence of parameter block position and eleids
  // used to build the mapextractor and the various maps to keep track of parameters and elements
  std::map<int, std::vector<int>> elemap;
  int nummyparams = 0;

  // fill it
  for (int i = 0; i < Discret()->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = Discret()->lRowElement(i);
    int elematid = ElementOptMat(actele);

    if (elematid == -1) continue;

    std::vector<int> actparapos = ParaPos().at(elematid);
    std::vector<int>::const_iterator it;
    for (it = actparapos.begin(); it != actparapos.end(); it++)
    {
      elemap[*it].push_back(actele->Id());
      nummyparams++;
    }
  }

  // generate global ids plus build map paramsLIDtoeleGID_
  std::map<int, std::vector<int>> gids;
  int count = 0;
  for (int i = 0; i < Discret()->Comm().NumProc(); i++)
  {
    if (Discret()->Comm().MyPID() == i)
    {
      for (int j = 0; j < NumParams(); j++)
      {
        for (int k = 0; k < (int)elemap[j].size(); k++)
        {
          gids[j].push_back(count);
          paramsLIDtoeleGID_.push_back(elemap[j].at(k));
          count++;
        }
      }
    }
    Discret()->Comm().Broadcast(&count, 1, i);
  }

  // build map eleGIDtoparamsLID_
  for (int i = 0; i < (int)paramsLIDtoeleGID_.size(); i++)
  {
    // the blocks are ordered so this is simply
    eleGIDtoparamsLID_[paramsLIDtoeleGID_[i]].push_back(i);
  }

  // the full map of the vector layout
  paramlayoutmap_ = Teuchos::rcp(
      new Epetra_Map(-1, nummyparams, 0, *(DRT::Problem::Instance()->GetNPGroup()->LocalComm())));
  paramlayoutmapunique_ = Teuchos::rcp(new Epetra_Map(*paramlayoutmap_));

  // the partial maps:
  std::vector<Teuchos::RCP<const Epetra_Map>> partials;
  for (int i = 0; i < NumParams(); i++)
  {
    partials.push_back(
        Teuchos::rcp(new Epetra_Map(-1, gids[i].size(), gids[i].data(), 0, Discret()->Comm())));
  }

  // finally build the MapExtractor
  paramapextractor_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*paramlayoutmap_, partials));

  optparams_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_, 1, true));
  optparams_initial_ = Teuchos::rcp(new Epetra_MultiVector(*paramlayoutmap_, 1, true));

  // set up restrictor and prolongator (identity matrices)
  projector_ = Teuchos::rcp(
      new Epetra_CrsMatrix(Copy, *paramlayoutmapunique_, *paramlayoutmapunique_, 1, false));

  // insert ones onto the diagonal
  double values = 1.0;
  for (int i = 0; i < projector_->NumMyRows(); i++)
  {
    int gid = paramlayoutmapunique_->GID(i);
    projector_->InsertGlobalValues(gid, 1, &values, &gid);
  }
  projector_->FillComplete();

  // initialize parameter vector from material parameters given in the input file
  InitParams();

  // Create elemetwise graph structure, before some specialized
  // class might overwrite some maps
  CreateGraph();
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerElement::FillParameters(Teuchos::RCP<Epetra_MultiVector> params)
{
  params->PutScalar(0.0);

  // loop the parameter blocks
  for (int k = 0; k < paramapextractor_->NumMaps(); k++)
  {
    Teuchos::RCP<Epetra_Vector> tmp = paramapextractor_->ExtractVector(*(*optparams_)(0), k);
    for (int i = 0; i < tmp->MyLength(); i++)
    {
      int pgid =
          tmp->Map().GID(i);  // !! the local id of the partial map is not the local parameter id!!
      int plid = paramlayoutmap_->LID(pgid);
      params->ReplaceGlobalValue(paramsLIDtoeleGID_[plid], k, (*tmp)[i]);
    }
  }
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerElement::InitParameters(int parapos, double val)
{
  Teuchos::RCP<Epetra_Vector> tmp =
      Teuchos::rcp(new Epetra_Vector(*paramapextractor_->Map(parapos), false));
  tmp->PutScalar(val);

  paramapextractor_->InsertVector(tmp, parapos, Teuchos::rcp((*optparams_)(0), false));
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerElement::ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint,
    double val, int elepos, int paraposglobal, int paraposlocal)
{
  if (eleGIDtoparamsLID_.find(elepos) == eleGIDtoparamsLID_.end())
    dserror("proc %d, ele %d not in this map", Discret()->Comm().MyPID(), elepos);

  int plid = eleGIDtoparamsLID_[elepos].at(paraposlocal);
  int success = dfint->SumIntoMyValue(plid, 0, val);
  if (success != 0) dserror("gid %d is not on this processor", plid);
}

void INVANA::MatParManagerPerElement::Finalize(
    Teuchos::RCP<Epetra_MultiVector> source, Teuchos::RCP<Epetra_MultiVector> target)
{
  // nothing to be summed across processors. since both maps are
  // the same by construction they can be just added up
  target->Update(1.0, *source, 1.0);

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerElement::ApplyParametrization(
    DcsMatrix& matrix, Teuchos::RCP<Epetra_MultiVector> diagonals)
{
  Teuchos::RCP<Epetra_Vector> diagonal =
      Teuchos::rcp(new Epetra_Vector(*paramlayoutmapunique_, true));
  matrix.ExtractDiagonalCopy(*diagonal);

  // loop the parameter blocks
  for (int k = 0; k < paramapextractor_->NumMaps(); k++)
  {
    Teuchos::RCP<Epetra_Vector> tmp = paramapextractor_->ExtractVector(*diagonal, k);
    for (int i = 0; i < tmp->MyLength(); i++)
    {
      int pgid =
          tmp->Map().GID(i);  // !! the local id of the partial map is not the local parameter id!!
      int plid = paramlayoutmap_->LID(pgid);
      diagonals->ReplaceGlobalValue(paramsLIDtoeleGID_[plid], k, (*tmp)[i]);
    }
  }

  return;
}

void INVANA::MatParManagerPerElement::CreateGraph()
{
  if (has_graph_)
    dserror("the graph should only be computed once. Maps might not match anylonger!");

  int maxbw = 6;  // based on connectivity for hex8 elements
  graph_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *(paramapextractor_->FullMap()), maxbw, false));

  for (int i = 0; i < paramapextractor_->NumMaps(); i++)
    FillAdjacencyMatrix(*(paramapextractor_->Map(i)), graph_);

  // Finalize the graph ...
  graph_->FillComplete();

  // put zeros one the diagonal; the diagonal is the "self weight" and it should never
  // be used somewhere since its meaningless but its better to have 0.0 than some
  // random value resulting from redundant inserting during FillAdjacencyMatrix
  Teuchos::RCP<Epetra_Vector> diagonal =
      Teuchos::rcp(new Epetra_Vector(*(paramapextractor_->FullMap()), true));
  diagonal->PutScalar(0.0);
  graph_->ReplaceDiagonalValues(*diagonal);

  has_graph_ = true;

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerElement::FillAdjacencyMatrix(
    const Epetra_Map& paramrowmap, Teuchos::RCP<Epetra_CrsMatrix> graph)
{
  /*-------------------------------------------------------------------
   * STEP 1: loop elements in elerowmap and store map of faces(vector
   * of NodeIds in sorted order) with corresponding gids of the elements
   */
  std::map<std::vector<int>, std::vector<int>>
      facemap;                                    // map faces to corresponding elements/parameters
  std::map<std::vector<int>, double> faceweight;  // map of faces to its area
  for (int i = 0; i < paramrowmap.NumMyElements(); i++)
  {
    // the current element
    int pgid =
        paramrowmap.GID(i);  // !! the local id of the partial map is not the local parameter id!!
    int plid = paramapextractor_->FullMap()->LID(pgid);
    int elegid = paramsLIDtoeleGID_[plid];
    DRT::Element* ele = Discret()->gElement(elegid);
    if (ele == NULL) dserror("element not found here");

    // decide whether we deal with a 2D or 3D discretization
    unsigned int nele = 0;
    const DRT::Element::DiscretizationType distype = ele->Shape();
    std::vector<std::vector<int>> faces;
    if (ele->NumSurface() > 1)  // 2D boundary element and 3D parent element
    {
      nele = ele->NumSurface();
      faces = DRT::UTILS::getEleNodeNumberingSurfaces(distype);
    }
    else if (ele->NumSurface() == 1)  // 1D boundary element and 2D parent element
    {
      nele = ele->NumLine();
      faces = DRT::UTILS::getEleNodeNumberingLines(distype);
    }
    else
      dserror("creating internal faces for 1D elements (would be points) not implemented yet");

    if (nele != faces.size()) dserror("number of surfaces or lines does not match!");

    // get the surface/line elements for area computation
    std::vector<Teuchos::RCP<DRT::Element>> surfs;
    if (ele->NumSurface() > 1)
      surfs = ele->Surfaces();
    else if (ele->NumSurface() == 1)
      surfs = ele->Lines();
    else
      dserror("creating internal faces for 1D elements (would be points) not implemented yet");

    // get nodes of each of this element's face
    for (unsigned int iele = 0; iele < nele; iele++)
    {
      // allocate node vectors
      unsigned int nnode = faces[iele].size();
      std::vector<int> nodeids(nnode);

      // get connectivity info
      for (unsigned int inode = 0; inode < nnode; inode++)
      {
        nodeids[inode] = ele->NodeIds()[faces[iele][inode]];
      }

      // get the area of this face
      Teuchos::ParameterList p;
      SetAction(p);
      p.set("area", 0.0);
      DRT::Element::LocationArray la(Discret()->NumDofSets());
      surfs[iele]->LocationVector(*Discret(), la, false);
      // initialize element vectors
      int ndof = ele->NumNode() * 3;
      Epetra_SerialDenseMatrix elematrix1(ndof, ndof, false);
      Epetra_SerialDenseMatrix elematrix2(ndof, ndof, false);
      Epetra_SerialDenseVector elevector1(ndof);
      Epetra_SerialDenseVector elevector2(ndof);
      Epetra_SerialDenseVector elevector3(ndof);
      surfs[iele]->Evaluate(
          p, *Discret(), la, elematrix1, elematrix2, elevector1, elevector2, elevector3);
      double area = p.get("area", -1.0);
      if (area < 0.0) dserror("area computation of surface failed");

      // sort the nodes in faces to make them to be used for keys in the facemap.
      // they are not unique (in a parallel layout sense) though since a face can
      // be on multiple processors
      std::sort(nodeids.begin(), nodeids.end());

      // store parameter global id corresponding to this face
      facemap[nodeids].push_back(pgid);
      faceweight[nodeids] = area;
    }
  }  // loop local elements


  /*-------------------------------------------------------------------
   *  STEP 2: Compute weight with respect to the average area of all
   *  faces. Make faceweight redundant all an procs to be able to
   *  compute the correct average face area on all procs; alternatively
   *  just put weight of one between all elements
   */
  INPAR::INVANA::StatInvGraphWeight weighttype;
  weighttype =
      DRT::INPUT::IntegralValue<INPAR::INVANA::StatInvGraphWeight>(Inpar(), "GRAPHWEIGHTS");

  switch (weighttype)
  {
    case INPAR::INVANA::stat_inv_graph_area:
    {
      LINALG::GatherAll<std::vector<int>, double>(faceweight, paramrowmap.Comm());
      std::map<std::vector<int>, double>::iterator weightsit;

      double avgarea = 0.0;
      for (weightsit = faceweight.begin(); weightsit != faceweight.end(); weightsit++)
        avgarea += weightsit->second;

      avgarea = avgarea / faceweight.size();

      for (weightsit = faceweight.begin(); weightsit != faceweight.end(); weightsit++)
        weightsit->second = avgarea / weightsit->second;
    }
    break;
    case INPAR::INVANA::stat_inv_graph_unity:
    {
      std::map<std::vector<int>, double>::iterator weightsit;
      for (weightsit = faceweight.begin(); weightsit != faceweight.end(); weightsit++)
        weightsit->second = 1.0;
    }
    break;
    default:
      dserror("no proper method to specify adjacency matrix specidied");
      break;
  }

  /*-------------------------------------------------------------------
   * STEP 3: gather for each face on each proc the same set of corresponding
   * elements; ntargetprocs is equal to the total number of processors to make
   * data redundant on all procs the face weights don't need to be communicated
   * since they are locally available anyways and have furthemore already been
   * made redundant.
   */
  const int numprocs = Discret()->Comm().NumProc();
  int allproc[numprocs];
  for (int i = 0; i < numprocs; ++i) allproc[i] = i;

  for (int i = 0; i < Discret()->Comm().NumProc(); i++)
  {
    // by looping the procs one by one we need to make sure that every proc
    // does the same number of loops regardless of its own lenght of keys
    // in the facemap. so the actual number of loops is defined by the current
    // proc and distributed in localmapsize
    int localmapsize;
    if (Discret()->Comm().MyPID() == i) localmapsize = facemap.size();
    Discret()->Comm().Broadcast(&localmapsize, 1, i);

    // now iterate as often as there are faces on proc i
    std::map<std::vector<int>, std::vector<int>>::iterator face_it(facemap.begin());
    for (int j = 0; j < localmapsize; j++)
    {
      // get length of current face-key on all procs
      int keylength = 0;
      if (Discret()->Comm().MyPID() == i) keylength = face_it->first.size();
      Discret()->Comm().Broadcast(&(keylength), 1, i);

      // send current face-key to all procs
      std::vector<int> facekey(keylength, 0);
      if (Discret()->Comm().MyPID() == i) facekey = face_it->first;
      Discret()->Comm().Broadcast(&(facekey[0]), keylength, i);

      // check whether one of the other procs also has this key and write IDs of
      // procs who own this face in "sowningprocs" and distribute this knowledge
      // among all procs to rowningprocs
      std::map<std::vector<int>, std::vector<int>>::iterator face_abroad(facemap.find(facekey));
      std::vector<int> sowningprocs;
      std::vector<int> rowningprocs;
      if (face_abroad != facemap.end() && Discret()->Comm().MyPID() != i)
        sowningprocs.push_back(Discret()->Comm().MyPID());
      LINALG::Gather(sowningprocs, rowningprocs, numprocs, allproc, Discret()->Comm());

      // now bring parameters corresponding to this face on the other procs to proc i
      // (they are send to all procs but only proc i stores them in the map with
      // the current face-key)
      std::vector<int> sparams;
      std::vector<int> rparams;
      if (std::find(rowningprocs.begin(), rowningprocs.end(), Discret()->Comm().MyPID()) !=
          rowningprocs.end())
        sparams = facemap[facekey];
      LINALG::Gather(sparams, rparams, numprocs, allproc, Discret()->Comm());

      // store additional elements on proc i
      if (Discret()->Comm().MyPID() == i)
      {
        for (int iele = 0; iele < (int)rparams.size(); iele++)
        {
          // depending on which proc comes first it might be that parameters
          // are added redundantly to a face which leads to undesired summation
          // upon inserting weights into the matrix. So only add if not existent yet
          if (std::find(facemap[facekey].begin(), facemap[facekey].end(), rparams[iele]) ==
              facemap[facekey].end())
            facemap[facekey].push_back(rparams[iele]);
        }
      }

      // increase face pointer on this proc
      if (Discret()->Comm().MyPID() == i) face_it++;

    }  // faces on each proc
  }    // procs

  /*-------------------------------------------------------------------
   * STEP 4: sort elements into graph according to neighbour information
   * in facemap
   */
  std::map<std::vector<int>, std::vector<int>>::iterator faces;
  for (faces = facemap.begin(); faces != facemap.end(); faces++)
  {
    // all parameter ids connected to this face
    std::vector<int> parameters = faces->second;

    // weight corresponding to these parameters
    std::vector<double> weights(parameters.size(), faceweight[faces->first]);

    for (int iele = 0; iele < (int)parameters.size(); iele++)
    {
      int globalrow = parameters[iele];
      if (paramrowmap.MyGID(globalrow))
      {
        // like this the diagonal entries are inserted redundantly and summed up
        // after FillComplete() is called; they are more or less useless anyways
        graph->InsertGlobalValues(globalrow, parameters.size(), &weights[0], &parameters[0]);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix> INVANA::MatParManagerPerElement::InitialCovariance()
{
  // the best inital guess to get from datfile input is a
  // unit diagonal. This is already available via the projector_
  return projector_;
}


/*----------------------------------------------------------------------*/
void INVANA::MatParManagerPerElement::SetAction(Teuchos::ParameterList& p)
{
  p.set("action", "calc_struct_area");
  return;
}

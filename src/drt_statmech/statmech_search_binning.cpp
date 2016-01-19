/*!----------------------------------------------------------------------
\file statmech_search_binning.cpp
\brief volume partitioning/binning search
<pre>
Maintainer: Kei Müller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
*----------------------------------------------------------------------*/
#include "statmech_search.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 |  constructor (public)                                   mueller 11/13|
 *----------------------------------------------------------------------*/
STATMECH::SEARCH::BinSearch::BinSearch(Teuchos::RCP<DRT::Discretization> discret,
                                       const std::vector<double>&        rootboxdim,
                                       const std::vector<int>&           resolution):
discret_(discret),
resolution_(resolution)
{
  // some sanity checks
  if((int)resolution_.size()!=3)
    dserror("Incorrect size %i of resolution vector!", (int)resolution_.size());
  if(resolution[0]<=0 || resolution[1]<=0 || resolution[2]<=0)
    dserror("Set resolution to non-zero values for all dimensions!");
  if(resolution[0]<3 || resolution[1]<3 || resolution[2]<3)
    dserror("Search resolution needs to be at least 3 in each direction!");
  if(rootboxdim[0]<=0.0 || rootboxdim[1]<=0.0 || rootboxdim[2]<=0.0)
    dserror("Incorrect root box dimensions: %d, %d, %d !", rootboxdim[0], rootboxdim[1],rootboxdim[2]);

  SetupRootBin(rootboxdim);

  BuildBinMaps();

  return;
}

STATMECH::SEARCH::BinSearch::BinSearch::Bin::Bin(const int& binid):
binid_(binid)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Assign positions to bins (public)                      mueller 11/13|
 *----------------------------------------------------------------------*/
void STATMECH::SEARCH::BinSearch::AssignPositionsToBins(const Teuchos::RCP<Epetra_MultiVector>& fullpositions)
{
  // note: positions have to be fully overlapping!
  // step 1: assign positions to bins (processor-wise unique)
  std::vector<std::vector<int> > binhaspositionsrowstd(binrowmap_->NumMyElements(),std::vector<int>());

  for(int i=0; i<fullpositions->MyLength(); i++)
  {
    std::vector<int> indices(fullpositions->NumVectors(),-1);
    for(int j=0; j<(int)indices.size(); j++)
    {
      indices[j] = (int)std::floor(((*fullpositions)[j][i]-rootlimits_(2*j))/(rootlimits_(2*j+1)-rootlimits_(2*j))*(double)resolution_[j]);
      if(indices[j]>=resolution_[j] || indices[j]<0)
        dserror("Bin index indices[%i] = %i is outside of the rootlimits! Search resolution: %i", j, indices[j], resolution_[j]);
    }

    const int binid = CalculateBinId(indices);

    if(binrowmap_->LID(binid)!=-1)
    {
      int posgid = fullpositions->Map().GID(i);
      binhaspositionsrowstd[binrowmap_->LID(binid)].push_back(posgid);
    }
  }

  // step 2: make full information accessible to all Procs
  // find largest number of positions in bins for creation of Epetra_MultiVector
  int maxlocalentries = -1;
  for(int i=0; i<(int)binhaspositionsrowstd.size(); i++)
    maxlocalentries = std::max(maxlocalentries,(int)binhaspositionsrowstd[i].size());
  int maxglobalentries = 0;
  discret_->Comm().MaxAll(&maxlocalentries,&maxglobalentries,1);

  Teuchos::RCP<Epetra_MultiVector> binhaspositionsrow = Teuchos::rcp(new Epetra_MultiVector(*binrowmap_,maxglobalentries));
  Teuchos::RCP<Epetra_MultiVector> binhaspositions = Teuchos::rcp(new Epetra_MultiVector(*binfullmap_,maxglobalentries));
  // -2.0 signals empty/unused entry
  binhaspositionsrow->PutScalar(-2.0);
  if(maxglobalentries==0)
    maxglobalentries = 1;

  // translate to Epetra_MultiVector
  for(int i=0; i<binhaspositionsrow->MyLength(); i++)
    for(int j=0; j<(int)binhaspositionsrowstd[i].size(); j++)
      (*binhaspositionsrow)[j][i] = (double)binhaspositionsrowstd[i][j];

  // communication: make information redundant on all procs
  CommunicateMultiVector(binhaspositionsrow,binhaspositions,false,true);

  // Sort into Bins
  int checksum = 0;
  for(int i=0; i<binhaspositions->MyLength(); i++)
    for(int j=0; j<binhaspositions->NumVectors(); j++)
    {
      if((*binhaspositions)[j][i]<0.0)
        break;
      bins_[i].AddMember((int)(*binhaspositions)[j][i]);
      checksum++;
    }

  return;
}

int STATMECH::SEARCH::BinSearch::CalculateBinId(const std::vector<int>& indices)
{
  if((int)indices.size()<2 || (int)indices.size()>3)
    dserror("Check the size (== %i) of your index vector", indices.size());
  const int bin = indices[0]*resolution_[0]*resolution_[0] + indices[1]*resolution_[1] + indices[2];
  return bin;
}

const std::vector<int> STATMECH::SEARCH::BinSearch::GetSurroundingBins(const std::vector<int> indices)
{
  if((int)indices.size()<2 || (int)indices.size()>3)
    dserror("Check the size (== %i) of your index vector", indices.size());
  std::vector<int> bins;
  bins.clear();

  for(int i=indices[0]-1; i<indices[0]+2; i++)
    for(int j=indices[1]-1; j<indices[1]+2; j++)
      for(int k=indices[2]-1; k<indices[2]+2; k++)
      {
        int currbin = ((i+resolution_[0])%(resolution_[0]))*resolution_[0]*resolution_[0]+((j+resolution_[1])%(resolution_[1]))*resolution_[1]+((k+resolution_[2])%resolution_[2]);
        bins.push_back(currbin);
      }

  return bins;
}


void STATMECH::SEARCH::BinSearch::SetupRootBin(const std::vector<double>& rootboxdim)
{
  for(int i=0; i<(int)rootlimits_.M(); i++)
  {
    if(i%2==0)
      rootlimits_(i) = 0.0;
    else
      rootlimits_(i) = rootboxdim[(i-1)/2];
  }
  return;
}

void STATMECH::SEARCH::BinSearch::BuildBinMaps()
{
  int numbins = resolution_[0]*resolution_[1]*resolution_[2];

  // build bin maps
  std::vector<int> gids;
  for (int i=0 ; i<numbins; i++ )
    gids.push_back(i);
  // build bin  column and row map
  binrowmap_ = Teuchos::rcp(new Epetra_Map((int)gids.size(), 0, discret_->Comm()));
  binfullmap_ = Teuchos::rcp(new Epetra_Map(-1, (int)gids.size(), &gids[0], 0, discret_->Comm()));

  // build fully overlapping bin vectors
  bins_.clear();
  for(int i=0; i<binfullmap_->NumMyElements(); i++)
    bins_.push_back(binfullmap_->GID(i));

  return;
}

bool STATMECH::SEARCH::BinSearch::HaveBin(const int& binid)
{
  if(binrowmap_->LID(binid)!=-1)
    return true;
  else
    return false;
}

/*-----------------------------------------------------------------------*
 *  | communicate MultiVector to all Processors               mueller 11/11 |
 *   *-----------------------------------------------------------------------*/
void STATMECH::SEARCH::BinSearch::CommunicateMultiVector(Teuchos::RCP<Epetra_MultiVector> InVec,
                                                         Teuchos::RCP<Epetra_MultiVector> OutVec,
                                                         bool                             doexport,
                                                         bool                             doimport,
                                                         bool                             zerofy,
                                                         bool                             exportinsert)
{
  Epetra_Export exporter(OutVec->Map(), InVec->Map());
  Epetra_Import importer(OutVec->Map(), InVec->Map());
  if(doexport)
  {
    if(discret_->Comm().MyPID()!=0 && zerofy)
      OutVec->PutScalar(0.0);
    if(exportinsert)
      InVec->Export(*OutVec, exporter, Insert);
    else
      InVec->Export(*OutVec, exporter, Add);
  }
  if(doimport)
    OutVec->Import(*InVec,importer,Insert);
  return;
}

void STATMECH::SEARCH::BinSearch::Bin::AddMember(const int& memberid)
{
  binmembers_.push_back(memberid);
  return;
}

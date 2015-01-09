/*!----------------------------------------------------------------------
\file xfluidfluid_state.cpp
\brief State class for (in)stationary XFEM fluid problems involving embedded
fluid meshes

<pre>
Maintainer:  Raffaela Kruse
             kruse@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15240
</pre>

*----------------------------------------------------------------------*/

#include "xfluidfluid_state.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_utils.H"

#include "../drt_fluid/fluid_utils.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

/*----------------------------------------------------------------------*
 |  Constructor for XFluidFluidState                         kruse 01/15 |
 *----------------------------------------------------------------------*/
FLD::XFluidFluidState::XFluidFluidState(
  Teuchos::RCP<GEO::CutWizard>& wizard,
  Teuchos::RCP<XFEM::XFEMDofSet>& dofset,
  Teuchos::RCP<const Epetra_Map> & xfluiddofrowmap,
  Teuchos::RCP<const Epetra_Map> & embfluiddofrowmap) :
  XFluidState(wizard,dofset,xfluiddofrowmap),
  xffluiddofrowmap_(LINALG::MergeMap(xfluiddofrowmap,embfluiddofrowmap,false)),
  xffluidsplitter_(Teuchos::rcp(new FLD::UTILS::FluidXFluidMapExtractor())),
  xffluidvelpressplitter_(Teuchos::rcp(new LINALG::MapExtractor())),
  embfluiddofrowmap_(embfluiddofrowmap)
{
  xffluidsplitter_->Setup(*xffluiddofrowmap_,xfluiddofrowmap,embfluiddofrowmap);
  InitSystemMatrix();
  InitStateVectors();
}

/*----------------------------------------------------------------------*
 |  Initialize (coupling) matrices & rhs-vectors
 |                                                         kruse 01/15
 *----------------------------------------------------------------------*/
void FLD::XFluidFluidState::InitSystemMatrix()
{
  // the combined fluid system matrix is not of FECrs-type - it is solely composed out of
  // fully assembled submatrices
  xffluidsysmat_   = Teuchos::rcp(new LINALG::SparseMatrix(*xffluiddofrowmap_,108,false,true));
}

/*----------------------------------------------------------------------*
 |  Initialize state vectors                                kruse 01/15 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluidState::InitStateVectors()
{
  // matrices & vectors for merged background & embedded fluid
  xffluidresidual_ = LINALG::CreateVector(*xffluiddofrowmap_,true);
  xffluidincvel_   = LINALG::CreateVector(*xffluiddofrowmap_,true);
  xffluidvelnp_    = LINALG::CreateVector(*xffluiddofrowmap_,true);
  xffluidveln_     = LINALG::CreateVector(*xffluiddofrowmap_,true);
  xffluidzeros_    = LINALG::CreateVector(*xffluiddofrowmap_,true);
}

/*----------------------------------------------------------------------*
 |  Create merged DBC map extractor                         kruse 01/15 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluidState::CreateMergedDBCMapExtractor(Teuchos::RCP<const LINALG::MapExtractor> embfluiddbcmaps)
{
  // create merged dbc map from both fluids
  std::vector<Teuchos::RCP<const Epetra_Map> > dbcmaps;
  dbcmaps.push_back(XFluidState::dbcmaps_->CondMap());
  dbcmaps.push_back(embfluiddbcmaps->CondMap());

  Teuchos::RCP<const Epetra_Map> xffluiddbcmap = LINALG::MultiMapExtractor::MergeMaps(dbcmaps);

  std::vector<Teuchos::RCP<const Epetra_Map> > othermaps;
  othermaps.push_back(XFluidState::dbcmaps_->OtherMap());
  othermaps.push_back(embfluiddbcmaps->OtherMap());
  Teuchos::RCP<const Epetra_Map> xffluidothermap = LINALG::MultiMapExtractor::MergeMaps(othermaps);

  xffluiddbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor(*xffluiddofrowmap_,xffluiddbcmap,xffluidothermap));
}

/*----------------------------------------------------------------------*
 |  Set dirichlet- and velocity/pressure-map extractor      kruse 01/15 |
 *----------------------------------------------------------------------*/
void FLD::XFluidFluidState::SetupMapExtractors( const Teuchos::RCP<DRT::Discretization> & xfluiddiscret,
                                                const Teuchos::RCP<DRT::Discretization> & embfluiddiscret,
                                                const double & time)
{
  // create merged dirichlet map extractor
  XFluidState::SetupMapExtractors(xfluiddiscret,time);
  xffluidsplitter_->Setup(*xffluiddofrowmap_,embfluiddofrowmap_,XFluidState::xfluiddofrowmap_);

  FLD::UTILS::SetupFluidFluidVelPresSplit(*xfluiddiscret,DRT::Problem::Instance()->NDim(),
      *embfluiddiscret,*xffluidvelpressplitter_,xffluiddofrowmap_);
}

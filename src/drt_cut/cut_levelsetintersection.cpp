/*!-----------------------------------------------------------------------------------------------*

\brief provides the basic functionality for cutting a mesh with a level set function

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249

\level 2
 *------------------------------------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>

#include "cut_side.H"
#include "cut_levelsetside.H"
#include "cut_levelsetintersection.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::LevelSetIntersection::LevelSetIntersection(const Epetra_Comm& comm, bool create_side)
    : ParentIntersection(comm.MyPID()), side_(Teuchos::null), comm_(&comm)
{
  if (create_side) AddCutSide(1);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::LevelSetIntersection::LevelSetIntersection(int myrank, bool create_side)
    : ParentIntersection(myrank), side_(Teuchos::null), comm_(NULL)
{
  if (create_side) AddCutSide(1);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::LevelSetIntersection::AddCutSide(int levelset_sid)
{
  if (!side_.is_null()) dserror("currently only one levelset-side is supported");

  // create the levelset-side
  side_ = Teuchos::rcp(Side::CreateLevelSetSide(levelset_sid));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::ElementHandle* GEO::CUT::LevelSetIntersection::AddElement(int eid,
    const std::vector<int>& nids, const Epetra_SerialDenseMatrix& xyz,
    DRT::Element::DiscretizationType distype, const double* lsv, const bool lsv_only_plus_domain,
    const bool& check_lsv)
{
  int numnode = nids.size();
  if (numnode != xyz.N()) dserror("node coordinate number mismatch");

  bool ltz = false;
  bool gtz = false;

  if (check_lsv)
  {
    // make sure element has LSV +ve and -ve in one of its nodes
    // ensures this is a cut element
    for (int i = 0; i < numnode; ++i)
    {
      if (lsv[i] <= REFERENCETOL) ltz = true;
      if (lsv[i] >= -REFERENCETOL) gtz = true;
    }
  }

  /* add all cut elements (different signs of levelset values) OR
   * if only plus domain is a physical field we have to add also
   * elements with only negative values (as they are not allowed to
   * carry DOFS at the end) */
  if ((not check_lsv) or (ltz and gtz) or (lsv_only_plus_domain and ltz))
  {
    // add all nodes to mesh
    for (int i = 0; i < numnode; ++i)
    {
      NormalMesh().GetNode(nids[i], &xyz(0, i), lsv[i]);
    }

    // create element
    return mesh_.CreateElement(eid, nids, distype);
  }

  return NULL;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::LevelSetIntersection::Cut_Mesh(bool screenoutput)
{
  TEUCHOS_FUNC_TIME_MONITOR("GEO::CUT --- 1/3 --- Cut");

  Mesh& m = NormalMesh();

  double t_diff = 0.0;
  // m.Status();
  {
    const double t_start = Teuchos::Time::wallTime();
    if (myrank_ == 0 and screenoutput)
      std::cout << "\n\t * 1/6 Cut ......................" << std::flush;

    m.Cut(*side_);

    if (myrank_ == 0 and screenoutput)
    {
      if (comm_) comm_->Barrier();
      t_diff = Teuchos::Time::wallTime() - t_start;

      printf("success! ( %10.4e secs )", t_diff);
    }
  }

  {
    const double t_start = Teuchos::Time::wallTime();
    if (myrank_ == 0 and screenoutput)
      std::cout << "\n\t * 2/6 MakeCutLines ............." << std::flush;

    m.MakeCutLines();

    if (myrank_ == 0 and screenoutput)
    {
      if (comm_) comm_->Barrier();
      t_diff = Teuchos::Time::wallTime() - t_start;

      printf("success! ( %10.4e secs )", t_diff);
    }
  }

  {
    const double t_start = Teuchos::Time::wallTime();
    if (myrank_ == 0 and screenoutput)
      std::cout << "\n\t * 3/6 MakeFacets ..............." << std::flush;

    m.MakeFacets();

    if (myrank_ == 0 and screenoutput)
    {
      if (comm_) comm_->Barrier();
      t_diff = Teuchos::Time::wallTime() - t_start;

      printf("success! ( %10.4e secs )", t_diff);
    }
  }

  {
    const double t_start = Teuchos::Time::wallTime();
    if (myrank_ == 0 and screenoutput)
      std::cout << "\n\t * 4/6 MakeVolumeCells .........." << std::flush;

    m.MakeVolumeCells();

    if (myrank_ == 0 and screenoutput)
    {
      if (comm_) comm_->Barrier();
      t_diff = Teuchos::Time::wallTime() - t_start;

      printf("success! ( %10.4e secs )", t_diff);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void GEO::CUT::LevelSetIntersection::Cut(
    bool include_inner, bool screenoutput, INPAR::CUT::VCellGaussPts VCellGP)
{
  //###########################################################################
  // STEP 1/3 CUT THE MESH
  //###########################################################################

  // m.Status();
  Cut_Mesh(screenoutput);


  //###########################################################################
  // STEP 2/3 ASSIGN DOFS
  //###########################################################################

  Mesh& m = NormalMesh();

  if (options_.FindPositions())
  {
    m.FindLSNodePositions();

    //=====================================================================

    m.FindNodalDOFSets(include_inner);

    //=====================================================================
  }
  //#############################################################################
  // STEP 3/3 FINALIZE, ASSIGN INTEGRATIONRULES
  //#############################################################################

  //  VCellGP=INPAR::CUT::VCellGaussPts_DirectDivergence;
  //  std::cout << "VCellGP: " << VCellGP << std::endl;
  double t_diff = 0.0;
  {
    const double t_start = Teuchos::Time::wallTime();
    if (myrank_ == 0 and screenoutput)
      std::cout << "\n\t * 5/6 CreateIntegrationCells ..." << std::flush;

    if (VCellGP == INPAR::CUT::VCellGaussPts_Tessellation)
      m.CreateIntegrationCells(0);
    else
      m.DirectDivergenceGaussRule(true, INPAR::CUT::BCellGaussPts_Tessellation);

    if (myrank_ == 0 and screenoutput)
    {
      if (comm_) comm_->Barrier();
      t_diff = Teuchos::Time::wallTime() - t_start;

      printf("success! ( %10.4e secs )", t_diff);
    }
  }

  {
    const double t_start = Teuchos::Time::wallTime();
    if (myrank_ == 0 and screenoutput)
      std::cout << "\n\t * 6/6 TestElementVolume ........" << std::flush;

    m.TestElementVolume(true, VCellGP);

    if (myrank_ == 0 and screenoutput)
    {
      if (comm_) comm_->Barrier();
      t_diff = Teuchos::Time::wallTime() - t_start;

      printf("success! ( %10.4e secs )\n\n", t_diff);
    }
  }
  // m.RemoveEmptyVolumeCells();

#ifdef DEBUGCUTLIBRARY
  m.DumpGmsh("mesh.pos");
  m.DumpGmshVolumeCells("volumecells.pos", true);
  m.DumpGmshIntegrationCells("integrationcells.pos");
#endif

  // m.TestElementVolume( true, VCellGP );

#ifdef DEBUGCUTLIBRARY
  // m.TestVolumeSurface(); //Broken test, needs to be fixed for proper usage.
  // If set sharply (i.e. 10^-11) this will cause some test cases to fail for Combust.
  //  thus this makes it likely that the cut is sensitive to 10^-11? Why?
  m.TestFacetArea();
#endif

  // Test if the volume-cell integration of DD and Tes are producing the same volumes
#ifdef DEBUGCUTLIBRARY
  {
    if (VCellGP == INPAR::CUT::VCellGaussPts_Tessellation)
    {
      DebugCut(m);
    }
  }
#endif


  //######################################################################################

}  // GEO::CUT::LevelSetIntersection::Cut

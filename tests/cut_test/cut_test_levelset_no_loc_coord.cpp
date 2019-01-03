/*!----------------------------------------------------------------------
\brief Test for the CUT Library
\file cut_test_levelset_no_loc_coord.cpp

\level 1

\maintainer Ager Christoph
*----------------------------------------------------------------------*/
// Issues with non-local cut!
// Tried to be resolved here (Only for small deformation created from cubit).
// Sheared elements needs to be tested separately
// m.w. 05.2015

#include "../../src/drt_cut/cut_options.H"
#include "../../src/drt_cut/cut_mesh.H"
#include "../../src/drt_cut/cut_element.H"
#include "../../src/drt_cut/cut_levelsetintersection.H"

// Added
#include "../../src/drt_cut/cut_volumecell.H"
#include "../../src/drt_cut/cut_facet.H"
#include "../../src/drt_cut/cut_integrationcell.H"
#include "../../src/drt_cut/cut_utils.H"
#include "../../src/drt_cut/sorted_vector.H"
#include "../../src/drt_cut/cut_combintersection.H"



#include "cut_test_utils.H"

#include <iterator>

#define PRECISION24
//#define GMSH_OUTPUT_LSNOLOC_CUT_TEST

Epetra_SerialDenseMatrix GetLocalHex8Coords()
{
  Epetra_SerialDenseMatrix xyze(3, 8);

  xyze(0, 0) = -1;
  xyze(1, 0) = -1;
  xyze(2, 0) = -1;

  xyze(0, 1) = 1;
  xyze(1, 1) = -1;
  xyze(2, 1) = -1;

  xyze(0, 2) = 1;
  xyze(1, 2) = 1;
  xyze(2, 2) = -1;

  xyze(0, 3) = -1;
  xyze(1, 3) = 1;
  xyze(2, 3) = -1;

  xyze(0, 4) = -1;
  xyze(1, 4) = -1;
  xyze(2, 4) = 1;

  xyze(0, 5) = 1;
  xyze(1, 5) = -1;
  xyze(2, 5) = 1;

  xyze(0, 6) = 1;
  xyze(1, 6) = 1;
  xyze(2, 6) = 1;

  xyze(0, 7) = -1;
  xyze(1, 7) = 1;
  xyze(2, 7) = 1;

  return xyze;
}


GEO::CUT::CombIntersection Cut_With_Tesselation(std::vector<int> nids, std::vector<double> lsvs,
    Epetra_SerialDenseMatrix xyze, std::string testname)
{
  // non-planar cut surface
  GEO::CUT::CombIntersection ci(-1);
  ci.GetOptions().Init_for_Cuttests();
  ci.AddLevelSetSide(1);

  ci.AddElement(1, nids, xyze, DRT::Element::hex8, &lsvs[0], false);

  ci.Cut(true);

  ci.NormalMesh().FindLSNodePositions();
  ci.NormalMesh().FindNodalDOFSets(true);

  ci.Cut_Finalize(true, INPAR::CUT::VCellGaussPts_Tessellation,
      INPAR::CUT::BCellGaussPts_Tessellation, false, true);
  std::cout << "TESSELATION, Cut without error." << std::endl;

#ifdef GMSH_OUTPUT_LSNOLOC_CUT_TEST
  // Gmsh-output Tesselation
  ci.NormalMesh().DumpGmsh(testname + "_tes.CUT.pos");
  ci.NormalMesh().DumpGmshVolumeCells(testname + "_tes.CUT_volumecells.pos", true);
  ci.DumpGmshIntegrationCells(testname + "_tes.CUT_integrationcells.pos");
#endif

  return ci;
}

GEO::CUT::CombIntersection Cut_With_DirectDivergence(std::vector<int> nids,
    std::vector<double> lsvs, Epetra_SerialDenseMatrix xyze, std::string testname)
{
  // non-planar cut surface
  GEO::CUT::CombIntersection ci(-1);
  ci.GetOptions().Init_for_Cuttests();
  ci.AddLevelSetSide(1);

  ci.AddElement(1, nids, xyze, DRT::Element::hex8, &lsvs[0], false);

  ci.Cut(true);

  ci.NormalMesh().FindLSNodePositions();
  ci.NormalMesh().FindNodalDOFSets(true);

  ci.Cut_Finalize(true, INPAR::CUT::VCellGaussPts_DirectDivergence,
      INPAR::CUT::BCellGaussPts_Tessellation, false, true);
  std::cout << "DIRECT DIVERGENCE, Cut without error." << std::endl;

#ifdef GMSH_OUTPUT_LSNOLOC_CUT_TEST
  // Gmsh-output DD
  ci.NormalMesh().DumpGmsh(testname + "_dd.CUT.pos");
  ci.NormalMesh().DumpGmshVolumeCells(testname + "_dd.CUT_volumecells.pos", true);
  // ci.DumpGmshIntegrationCells(testname+"_dd.CUT_integrationcells.pos");
#endif

  return ci;
}


void Test_LevelSetCut_Tesselation_and_DD(std::vector<int> nids, std::vector<double> lsvs,
    Epetra_SerialDenseMatrix xyze, std::string testname)
{
  // non-planar cut surface
  GEO::CUT::CombIntersection ci(-1);
  ci.GetOptions().Init_for_Cuttests();
  GEO::CUT::CombIntersection cidd(-1);
  cidd.GetOptions().Init_for_Cuttests();
  ci.AddLevelSetSide(1);
  cidd.AddLevelSetSide(1);

  ci.AddElement(1, nids, xyze, DRT::Element::hex8, &lsvs[0], false);
  cidd.AddElement(1, nids, xyze, DRT::Element::hex8, &lsvs[0], false);

  ci.Cut(true);
  cidd.Cut(true);

  ci.NormalMesh().FindLSNodePositions();
  ci.NormalMesh().FindNodalDOFSets(true);
  cidd.NormalMesh().FindLSNodePositions();
  cidd.NormalMesh().FindNodalDOFSets(true);

  ci.Cut_Finalize(true, INPAR::CUT::VCellGaussPts_Tessellation,
      INPAR::CUT::BCellGaussPts_Tessellation, false, true);
  std::cout << "TESSELATION, Cut without error." << std::endl;
  cidd.Cut_Finalize(true, INPAR::CUT::VCellGaussPts_DirectDivergence,
      INPAR::CUT::BCellGaussPts_Tessellation, false, true);
  std::cout << "DIRECT DIVERGENCE, Cut without error." << std::endl;

#ifdef GMSH_OUTPUT_LSNOLOC_CUT_TEST
  // Gmsh-output Tesselation
  ci.NormalMesh().DumpGmsh(testname + "_tes.CUT.pos");
  ci.NormalMesh().DumpGmshVolumeCells(testname + "_tes.CUT_volumecells.pos", true);
  ci.DumpGmshIntegrationCells(testname + "_tes.CUT_integrationcells.pos");
  // Gmsh-output DD
  cidd.NormalMesh().DumpGmsh(testname + "_dd.CUT.pos");
  cidd.NormalMesh().DumpGmshVolumeCells(testname + "_dd.CUT_volumecells.pos", true);
  // ci.DumpGmshIntegrationCells(testname+"_dd.CUT_integrationcells.pos");
#endif


  // TEST IF VOLUMES PREDICTED OF CELLS ARE SAME:
  //###########################################################
  std::vector<double> tessVol, dirDivVol;

  // Sum Tessellation volume of test
  const std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>& other_cells = ci.NormalMesh().VolumeCells();
  std::cout << "# Volume Cells Tesselation: " << other_cells.size() << std::endl;
  int iteration_VC = 0;
  for (std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>::const_iterator i = other_cells.begin();
       i != other_cells.end(); ++i)
  {
    GEO::CUT::VolumeCell* vc = &**i;
    tessVol.push_back(vc->Volume());

    // TEST THAT TESSELATION DOES NOT HAVE EMPTY VolumeCells.
    //++++++++++++++++++++
    iteration_VC++;
    std::cout << "VC(" << iteration_VC << "):" << std::endl;

    const GEO::CUT::plain_integrationcell_set& integrationcells = vc->IntegrationCells();
    std::cout << "Has #IC=" << integrationcells.size() << std::endl;

    if (integrationcells.size() == 0)
    {
      throw std::runtime_error("VolumeCell contains 0 integration cells.");
    }
    //++++++++++++++++++++
  }
  //------------------------------------------------------

  // Sum Direct Divergence of test
  const std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>& other_cellsdd =
      cidd.NormalMesh().VolumeCells();
  std::cout << "# Volume Cells Direct Divergence: " << other_cellsdd.size() << std::endl;
  for (std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>::const_iterator idd = other_cellsdd.begin();
       idd != other_cellsdd.end(); ++idd)
  {
    GEO::CUT::VolumeCell* vc = &**idd;
    dirDivVol.push_back(vc->Volume());
  }
  //------------------------------------------------------


  // Test Direct Divergence against Tesselation
  //------------------------------------------------------
  bool error = false;
  std::cout << "the volumes predicted by\n tessellation \t DirectDivergence\n";
  for (unsigned i = 0; i < tessVol.size(); i++)
  {
    std::cout << tessVol[i] << "\t" << dirDivVol[i] << "\n";
    if (fabs(tessVol[i] - dirDivVol[i]) > 1e-9) error = true;
  }
  if (error)
  {
    throw std::runtime_error("Volume predicted by one of the methods is wrong.");
  }
  //------------------------------------------------------
  //###########################################################
}


void Test_LevelSetCut_Tesselation_and_DD_same_VC(std::vector<int> nids, std::vector<double> lsvs,
    Epetra_SerialDenseMatrix xyze, std::string testname)
{
  GEO::CUT::CombIntersection ci = Cut_With_Tesselation(nids, lsvs, xyze, testname);

  // TEST IF VOLUMES PREDICTED OF CELLS ARE SAME:
  //###########################################################
  std::vector<double> tessVol, dirDivVol;
  double diff_tol = BASICTOL;  // 1e-20; //How sharp to test for.

  // Sum Tessellation volume of test
  const std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>& other_cells = ci.NormalMesh().VolumeCells();
  std::cout << "# Volume Cells Tesselation: " << other_cells.size() << std::endl;
  int iteration_VC = 0;
  for (std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>::const_iterator i = other_cells.begin();
       i != other_cells.end(); ++i)
  {
    GEO::CUT::VolumeCell* vc = &**i;
    tessVol.push_back(vc->Volume());

    // TEST THAT TESSELATION DOES NOT HAVE EMPTY VolumeCells.
    //++++++++++++++++++++
    iteration_VC++;
    std::cout << "VC(" << iteration_VC << "):" << std::endl;

    const GEO::CUT::plain_integrationcell_set& integrationcells = vc->IntegrationCells();
    std::cout << "Has #IC=" << integrationcells.size() << std::endl;

    if (integrationcells.size() == 0)
    {
      throw std::runtime_error("VolumeCell contains 0 integration cells.");
    }
    //++++++++++++++++++++
  }
  //------------------------------------------------------

  // Cut with DirectDivergence as well
  int counter = 1;
  for (std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>::const_iterator i = other_cells.begin();
       i != other_cells.end(); ++i)
  {
    GEO::CUT::VolumeCell* vc = &**i;
    //    std::cout << "DIRCETDIVERGENCE CUTTING VOLUMECELL #" << counter <<  ".....!" << std::endl;
    //    std::cout << "Volumecell Position: " << vc->Position() << std::endl;
    vc->DirectDivergenceGaussRule(
        vc->ParentElement(), ci.NormalMesh(), true, INPAR::CUT::BCellGaussPts_Tessellation);
    //    std::cout << "DIRCETDIVERGENCE CUT VOLUMECELL #" << counter <<  " WITHOUT ERROR!" <<
    //    std::endl << std::endl;
    counter++;
    dirDivVol.push_back(vc->Volume());
  }

  // Test Direct Divergence against Tesselation
  //------------------------------------------------------
  bool error = false;
  double sumtessvol = 0.0;
  for (unsigned i = 0; i < tessVol.size(); i++) sumtessvol += tessVol[i];

  std::cout << "Tolerance: " << diff_tol * sumtessvol << std::endl;

  std::cout << "the volumes predicted by\n Tessellation \t DirectDivergence\n";
  for (unsigned i = 0; i < tessVol.size(); i++)
  {
    std::cout << std::setprecision(17) << std::scientific << tessVol[i] << "\t"
              << dirDivVol[i];  //<<"\n";
    if (fabs(tessVol[i] - dirDivVol[i]) >
        diff_tol * (sumtessvol))  // Is the relative difference greater than the Tolerance?
    {
      std::cout << "\t difference: " << fabs(tessVol[i] - dirDivVol[i]);  // << "\n";
      error = true;
    }
    else if (isnan(tessVol[i]) or isnan(dirDivVol[i]))
      error = true;
    std::cout << "\n";
  }
  if (error)
  {
    throw std::runtime_error("Volume predicted by one of the methods is wrong.");
  }
  //------------------------------------------------------
  //###########################################################
}

void test_ls_hex8_magnus1()
{
  // simple hex8 element
  std::vector<int> nids(8);
  std::vector<double> lsvs(8);
  Epetra_SerialDenseMatrix xyze(3, 8);
  // Epetra_SerialDenseMatrix xyze_local( 3, 8 );


  for (int i = 0; i < 8; ++i)
  {
    nids[i] = i;
  }

#ifdef PRECISION24
  std::cout << "Precision 24" << std::endl;

  // eleID=47018 (precision 24)
  xyze(0, 0) = 0.449999999999999955591079;
  xyze(1, 0) = -0.250000000000000055511151;
  xyze(2, 0) = 0.0999999999999999500399639;

  xyze(0, 1) = 0.449999999999999955591079;
  xyze(1, 1) = -0.25;
  xyze(2, 1) = 0.0499999999999998639976795;

  xyze(0, 2) = 0.450000000000000066613381;
  xyze(1, 2) = -0.200000000000000066613381;
  xyze(2, 2) = 0.0499999999999998778754673;

  xyze(0, 3) = 0.449999999999999955591079;
  xyze(1, 3) = -0.200000000000000066613381;
  xyze(2, 3) = 0.0999999999999999222843883;

  xyze(0, 4) = 0.5;
  xyze(1, 4) = -0.25;
  xyze(2, 4) = 0.100000000000000005551115;

  xyze(0, 5) = 0.5;
  xyze(1, 5) = -0.25;
  xyze(2, 5) = 0.0499999999999999195088307;

  xyze(0, 6) = 0.500000000000000111022302;
  xyze(1, 6) = -0.199999999999999983346655;
  xyze(2, 6) = 0.0499999999999999333866185;

  xyze(0, 7) = 0.5;
  xyze(1, 7) = -0.200000000000000066613381;
  xyze(2, 7) = 0.0999999999999999777955395;

  lsvs[0] = 0.024404424085075815398227;
  lsvs[1] = 0.0172040216394300227165104;
  lsvs[2] = -0.00502525316941665467496136;
  lsvs[3] = 0.00249378105604447508625299;
  lsvs[4] = 0.0678908345800273149706072;
  lsvs[5] = 0.0612486080160912216285851;
  lsvs[6] = 0.0408326913195984353421863;
  lsvs[7] = 0.0477225575051661854431018;
#else
  std::cout << "Precision 16" << std::endl;

  // qhull with QdB fails this one?!
  // should prolly not be qhull at fault... Something else is fundamentaly wrong..

  // eleID=47018 (precision 16)
  xyze(0, 0) = 0.45;
  xyze(1, 0) = -0.2500000000000001;
  xyze(2, 0) = 0.09999999999999995;

  xyze(0, 1) = 0.45;
  xyze(1, 1) = -0.25;
  xyze(2, 1) = 0.04999999999999986;

  xyze(0, 2) = 0.4500000000000001;
  xyze(1, 2) = -0.2000000000000001;
  xyze(2, 2) = 0.04999999999999988;

  xyze(0, 3) = 0.45;
  xyze(1, 3) = -0.2000000000000001;
  xyze(2, 3) = 0.09999999999999992;

  xyze(0, 4) = 0.5;
  xyze(1, 4) = -0.25;
  xyze(2, 4) = 0.1;

  xyze(0, 5) = 0.5;
  xyze(1, 5) = -0.25;
  xyze(2, 5) = 0.04999999999999992;

  xyze(0, 6) = 0.5000000000000001;
  xyze(1, 6) = -0.2;
  xyze(2, 6) = 0.04999999999999993;

  xyze(0, 7) = 0.5;
  xyze(1, 7) = -0.2000000000000001;
  xyze(2, 7) = 0.09999999999999998;

  lsvs[0] = 0.02440442408507582;
  lsvs[1] = 0.01720402163943002;
  lsvs[2] = -0.005025253169416655;
  lsvs[3] = 0.002493781056044475;
  lsvs[4] = 0.06789083458002731;
  lsvs[5] = 0.06124860801609122;
  lsvs[6] = 0.04083269131959844;
  lsvs[7] = 0.04772255750516619;
#endif


  std::string testname("xxx_cut_test_ls_hex8_magnus1");

#ifdef LOCAL_COORDS
  Test_LevelSetCut_Tesselation_and_DD_same_VC(nids, lsvs, GetLocalHex8Coords(), testname);
#else
  Test_LevelSetCut_Tesselation_and_DD_same_VC(nids, lsvs, xyze, testname);
#endif
}

void test_ls_hex8_magnus2()
{
  // simple hex8 element
  std::vector<int> nids(8);
  std::vector<double> lsvs(8);
  Epetra_SerialDenseMatrix xyze(3, 8);

  for (int i = 0; i < 8; ++i)
  {
    nids[i] = i;
  }

  //  //eleID=43985  ve=0.000125  vc=0.000124881  vd= 1.18574e-07  err=0.00094859
#ifdef PRECISION24
  std::cout << "Precision 24" << std::endl;

  xyze(0, 0) = 0.350000000000000088817842;
  xyze(1, 0) = -0.0499999999999999472644063;
  xyze(2, 0) = -0.25;

  xyze(0, 1) = 0.350000000000000088817842;
  xyze(1, 1) = -0.0499999999999999333866185;
  xyze(2, 1) = -0.300000000000000044408921;

  xyze(0, 2) = 0.350000000000000088817842;
  xyze(1, 2) = 0;
  xyze(2, 2) = -0.300000000000000099920072;

  xyze(0, 3) = 0.350000000000000088817842;
  xyze(1, 3) = 0;
  xyze(2, 3) = -0.25;

  xyze(0, 4) = 0.400000000000000133226763;
  xyze(1, 4) = -0.0499999999999999611421941;
  xyze(2, 4) = -0.25;

  xyze(0, 5) = 0.40000000000000002220446;
  xyze(1, 5) = -0.0499999999999999472644063;
  xyze(2, 5) = -0.300000000000000044408921;

  xyze(0, 6) = 0.400000000000000133226763;
  xyze(1, 6) = -4.1633363423443376428862e-18;
  xyze(2, 6) = -0.300000000000000155431223;

  xyze(0, 7) = 0.400000000000000133226763;
  xyze(1, 7) = -3.46944695195361418882385e-18;
  xyze(2, 7) = -0.25;

  lsvs[0] = -0.0669872981077805906835465;
  lsvs[1] = -0.0363190752252147142087324;
  lsvs[2] = -0.0390227771353555130673385;
  lsvs[3] = -0.069883736647868621716384;
  lsvs[4] = -0.025658350974743004968559;
  lsvs[5] = 0.00249378105604458610855545;
  lsvs[6] = 2.22044604925031308084726e-16;
  lsvs[7] = -0.028300943397169708859451;
#else
  std::cout << "Precision 16" << std::endl;

  //  EleID=43985
  xyze(0, 0) = 0.3500000000000001;
  xyze(1, 0) = -0.04999999999999995;
  xyze(2, 0) = -0.25;

  xyze(0, 1) = 0.3500000000000001;
  xyze(1, 1) = -0.04999999999999993;
  xyze(2, 1) = -0.3;

  xyze(0, 2) = 0.3500000000000001;
  xyze(1, 2) = 0;
  xyze(2, 2) = -0.3000000000000001;

  xyze(0, 3) = 0.3500000000000001;
  xyze(1, 3) = 0;
  xyze(2, 3) = -0.25;

  xyze(0, 4) = 0.4000000000000001;
  xyze(1, 4) = -0.04999999999999996;
  xyze(2, 4) = -0.25;

  xyze(0, 5) = 0.4;
  xyze(1, 5) = -0.04999999999999995;
  xyze(2, 5) = -0.3;

  xyze(0, 6) = 0.4000000000000001;
  xyze(1, 6) = -4.163336342344338e-18;
  xyze(2, 6) = -0.3000000000000002;

  xyze(0, 7) = 0.4000000000000001;
  xyze(1, 7) = -3.469446951953614e-18;
  xyze(2, 7) = -0.25;

  lsvs[0] = -0.06698729810778059;
  lsvs[1] = -0.03631907522521471;
  lsvs[2] = -0.03902277713535551;
  lsvs[3] = -0.06988373664786862;
  lsvs[4] = -0.025658350974743;
  lsvs[5] = 0.002493781056044586;
  lsvs[6] = 2.220446049250313e-16;
  lsvs[7] = -0.02830094339716971;
#endif

  std::string testname("xxx_cut_test_ls_hex8_magnus2");

#ifdef LOCAL_COORDS
  Test_LevelSetCut_Tesselation_and_DD_same_VC(nids, lsvs, GetLocalHex8Coords(), testname);
#else
  Test_LevelSetCut_Tesselation_and_DD_same_VC(nids, lsvs, xyze, testname);
#endif
}


void test_ls_hex8_magnus3()
{
  // simple hex8 element
  std::vector<int> nids(8);
  std::vector<double> lsvs(8);
  Epetra_SerialDenseMatrix xyze(3, 8);

  for (int i = 0; i < 8; ++i)
  {
    nids[i] = i;
  }

  // eleID = 45458 prec:24
#ifdef PRECISION24
  std::cout << "Precision 24" << std::endl;

  xyze(0, 0) = 0.40000000000000002220446;
  xyze(1, 0) = -0.200000000000000066613381;
  xyze(2, 0) = 0.0999999999999999777955395;

  xyze(0, 1) = 0.40000000000000002220446;
  xyze(1, 1) = -0.199999999999999955591079;
  xyze(2, 1) = 0.0499999999999999333866185;

  xyze(0, 2) = 0.399999999999999911182158;
  xyze(1, 2) = -0.149999999999999966693309;
  xyze(2, 2) = 0.0499999999999999056310429;

  xyze(0, 3) = 0.40000000000000002220446;
  xyze(1, 3) = -0.149999999999999994448885;
  xyze(2, 3) = 0.0999999999999999777955395;

  xyze(0, 4) = 0.449999999999999955591079;
  xyze(1, 4) = -0.200000000000000066613381;
  xyze(2, 4) = 0.0999999999999999222843883;

  xyze(0, 5) = 0.450000000000000066613381;
  xyze(1, 5) = -0.200000000000000066613381;
  xyze(2, 5) = 0.0499999999999998778754673;

  xyze(0, 6) = 0.449999999999999955591079;
  xyze(1, 6) = -0.149999999999999994448885;
  xyze(2, 6) = 0.0499999999999998639976795;

  xyze(0, 7) = 0.449999999999999955591079;
  xyze(1, 7) = -0.150000000000000049960036;
  xyze(2, 7) = 0.0999999999999999777955395;

  lsvs[0] = -0.0417424305044158949762334;
  lsvs[1] = -0.0499999999999999888977698;
  lsvs[2] = -0.0698837366478687882498377;
  lsvs[3] = -0.0612517806303938816547827;
  lsvs[4] = 0.00249378105604447508625299;
  lsvs[5] = -0.00502525316941665467496136;
  lsvs[6] = -0.0230303992915272814911987;
  lsvs[7] = -0.0152320142583671214175922;
#else
  std::cout << "Precision 16" << std::endl;

  // eleID = 45458 prec:16
  xyze(0, 0) = 0.4;
  xyze(1, 0) = -0.2000000000000001;
  xyze(2, 0) = 0.09999999999999998;

  xyze(0, 1) = 0.4;
  xyze(1, 1) = -0.2;
  xyze(2, 1) = 0.04999999999999993;

  xyze(0, 2) = 0.3999999999999999;
  xyze(1, 2) = -0.15;
  xyze(2, 2) = 0.04999999999999991;

  xyze(0, 3) = 0.4;
  xyze(1, 3) = -0.15;
  xyze(2, 3) = 0.09999999999999998;

  xyze(0, 4) = 0.45;
  xyze(1, 4) = -0.2000000000000001;
  xyze(2, 4) = 0.09999999999999992;

  xyze(0, 5) = 0.4500000000000001;
  xyze(1, 5) = -0.2000000000000001;
  xyze(2, 5) = 0.04999999999999988;

  xyze(0, 6) = 0.45;
  xyze(1, 6) = -0.15;
  xyze(2, 6) = 0.04999999999999986;

  xyze(0, 7) = 0.45;
  xyze(1, 7) = -0.15;
  xyze(2, 7) = 0.09999999999999998;

  lsvs[0] = -0.04174243050441589;
  lsvs[1] = -0.04999999999999999;
  lsvs[2] = -0.06988373664786879;
  lsvs[3] = -0.06125178063039388;
  lsvs[4] = 0.002493781056044475;
  lsvs[5] = -0.005025253169416655;
  lsvs[6] = -0.02303039929152728;
  lsvs[7] = -0.01523201425836712;
#endif

  std::string testname("xxx_cut_test_ls_hex8_magnus3");

#ifdef LOCAL_COORDS
  Test_LevelSetCut_Tesselation_and_DD_same_VC(nids, lsvs, GetLocalHex8Coords(), testname);
#else
  Test_LevelSetCut_Tesselation_and_DD_same_VC(nids, lsvs, xyze, testname);
#endif
}

void test_ls_hex8_magnus4()
{
  // simple hex8 element
  std::vector<int> nids(8);
  std::vector<double> lsvs(8);
  Epetra_SerialDenseMatrix xyze(3, 8);

  for (int i = 0; i < 8; ++i)
  {
    nids[i] = i;
  }

#ifdef PRECISION24
  std::cout << "Precision 24" << std::endl;
  // ELEID: 37941 PREC:24
  xyze(0, 0) = 0.15000000000000002220446;
  xyze(1, 0) = 0.399999999999999911182158;
  xyze(2, 0) = -0.050000000000000044408921;

  xyze(0, 1) = 0.15000000000000002220446;
  xyze(1, 1) = 0.399999999999999911182158;
  xyze(2, 1) = -0.100000000000000116573418;

  xyze(0, 2) = 0.150000000000000049960036;
  xyze(1, 2) = 0.450000000000000177635684;
  xyze(2, 2) = -0.100000000000000116573418;

  xyze(0, 3) = 0.15000000000000002220446;
  xyze(1, 3) = 0.450000000000000177635684;
  xyze(2, 3) = -0.0500000000000000721644966;

  xyze(0, 4) = 0.200000000000000066613381;
  xyze(1, 4) = 0.399999999999999911182158;
  xyze(2, 4) = -0.050000000000000044408921;

  xyze(0, 5) = 0.200000000000000122124533;
  xyze(1, 5) = 0.399999999999999911182158;
  xyze(2, 5) = -0.100000000000000088817842;

  xyze(0, 6) = 0.200000000000000066613381;
  xyze(1, 6) = 0.450000000000000066613381;
  xyze(2, 6) = -0.100000000000000088817842;

  xyze(0, 7) = 0.200000000000000066613381;
  xyze(1, 7) = 0.449999999999999955591079;
  xyze(2, 7) = -0.0500000000000000721644966;

  lsvs[0] = -0.0698837366478687327386865;
  lsvs[1] = -0.0612517806303939371659339;
  lsvs[2] = -0.0152320142583668993729873;
  lsvs[3] = -0.0230303992915270039354425;
  lsvs[4] = -0.050000000000000044408921;
  lsvs[5] = -0.0417424305044160615096871;
  lsvs[6] = 0.00249378105604458610855545;
  lsvs[7] = -0.0050252531694167101861126;
#else
  std::cout << "Precision 16" << std::endl;

  // ELEID: 37941 PREC: 16
  xyze(0, 0) = 0.15;
  xyze(1, 0) = 0.3999999999999999;
  xyze(2, 0) = -0.05000000000000004;

  xyze(0, 1) = 0.15;
  xyze(1, 1) = 0.3999999999999999;
  xyze(2, 1) = -0.1000000000000001;

  xyze(0, 2) = 0.15;
  xyze(1, 2) = 0.4500000000000002;
  xyze(2, 2) = -0.1000000000000001;

  xyze(0, 3) = 0.15;
  xyze(1, 3) = 0.4500000000000002;
  xyze(2, 3) = -0.05000000000000007;

  xyze(0, 4) = 0.2000000000000001;
  xyze(1, 4) = 0.3999999999999999;
  xyze(2, 4) = -0.05000000000000004;

  xyze(0, 5) = 0.2000000000000001;
  xyze(1, 5) = 0.3999999999999999;
  xyze(2, 5) = -0.1000000000000001;

  xyze(0, 6) = 0.2000000000000001;
  xyze(1, 6) = 0.4500000000000001;
  xyze(2, 6) = -0.1000000000000001;

  xyze(0, 7) = 0.2000000000000001;
  xyze(1, 7) = 0.45;
  xyze(2, 7) = -0.05000000000000007;

  lsvs[0] = -0.06988373664786873;
  lsvs[1] = -0.06125178063039394;
  lsvs[2] = -0.0152320142583669;
  lsvs[3] = -0.023030399291527;
  lsvs[4] = -0.05000000000000004;
  lsvs[5] = -0.04174243050441606;
  lsvs[6] = 0.002493781056044586;
  lsvs[7] = -0.00502525316941671;
#endif



  std::string testname("xxx_cut_test_ls_hex8_magnus4");

#ifdef LOCAL_COORDS
  Test_LevelSetCut_Tesselation_and_DD_same_VC(nids, lsvs, GetLocalHex8Coords(), testname);
#else
  Test_LevelSetCut_Tesselation_and_DD_same_VC(nids, lsvs, xyze, testname);
#endif
}

void test_ls_hex8_magnus5()
{
  // simple hex8 element
  std::vector<int> nids(8);
  std::vector<double> lsvs(8);
  Epetra_SerialDenseMatrix xyze(3, 8);

  for (int i = 0; i < 8; ++i)
  {
    nids[i] = i;
  }

#ifdef PRECISION24
  std::cout << "Precision 24" << std::endl;
  // ELEID: 37941 PREC:24
  //  WARNING: One Volume-cell is empty!!
  //   !!!!!!!!!!! volume test failed: !!!!!!!!!!!!!!!eleID=22974  ve=0.000125  vc=0.000117163
  //   vd= 7.83715e-06  err=0.0626972
  xyze(0, 0) = -0.299999999999999933386619;
  xyze(1, 0) = -0.300000000000000044408921;
  xyze(2, 0) = 0.300000000000000044408921;

  xyze(0, 1) = -0.300000000000000044408921;
  xyze(1, 1) = -0.29999999999999998889777;
  xyze(2, 1) = 0.25;

  xyze(0, 2) = -0.300000000000000044408921;
  xyze(1, 2) = -0.25;
  xyze(2, 2) = 0.25;

  xyze(0, 3) = -0.300000000000000044408921;
  xyze(1, 3) = -0.25;
  xyze(2, 3) = 0.29999999999999998889777;

  xyze(0, 4) = -0.25;
  xyze(1, 4) = -0.299999999999999933386619;
  xyze(2, 4) = 0.29999999999999998889777;

  xyze(0, 5) = -0.25;
  xyze(1, 5) = -0.299999999999999933386619;
  xyze(2, 5) = 0.25;

  xyze(0, 6) = -0.25;
  xyze(1, 6) = -0.25;
  xyze(2, 6) = 0.25;

  xyze(0, 7) = -0.250000000000000055511151;
  xyze(1, 7) = -0.249999999999999972244424;
  xyze(2, 7) = 0.300000000000000044408921;

  lsvs[0] = 0.0196152422706632467708232;
  lsvs[1] = -0.0075571099101947591947237;
  lsvs[2] = -0.0363190752252147697198836;
  lsvs[3] = -0.0075571099101947591947237;
  lsvs[4] = -0.00755710991019481470587493;
  lsvs[5] = -0.0363190752252148252310349;
  lsvs[6] = -0.066987298107780701705849;
  lsvs[7] = -0.0363190752252147697198836;
#else
  std::cout << "Precision 16" << std::endl;
  xyze(0, 0) = -0.2999999999999999;
  xyze(1, 0) = -0.3;
  xyze(2, 0) = 0.3;

  xyze(0, 1) = -0.3;
  xyze(1, 1) = -0.3;
  xyze(2, 1) = 0.25;

  xyze(0, 2) = -0.3;
  xyze(1, 2) = -0.25;
  xyze(2, 2) = 0.25;

  xyze(0, 3) = -0.3;
  xyze(1, 3) = -0.25;
  xyze(2, 3) = 0.3;

  xyze(0, 4) = -0.25;
  xyze(1, 4) = -0.2999999999999999;
  xyze(2, 4) = 0.3;

  xyze(0, 5) = -0.25;
  xyze(1, 5) = -0.2999999999999999;
  xyze(2, 5) = 0.25;

  xyze(0, 6) = -0.25;
  xyze(1, 6) = -0.25;
  xyze(2, 6) = 0.25;

  xyze(0, 7) = -0.2500000000000001;
  xyze(1, 7) = -0.25;
  xyze(2, 7) = 0.3;

  lsvs[0] = 0.01961524227066325;
  lsvs[1] = -0.007557109910194759;
  lsvs[2] = -0.03631907522521477;
  lsvs[3] = -0.007557109910194759;
  lsvs[4] = -0.007557109910194815;
  lsvs[5] = -0.03631907522521483;
  lsvs[6] = -0.0669872981077807;
  lsvs[7] = -0.03631907522521477;
#endif

  std::string testname("xxx_cut_test_ls_hex8_magnus5");

#ifdef LOCAL_COORDS
  Test_LevelSetCut_Tesselation_and_DD_same_VC(nids, lsvs, GetLocalHex8Coords(), testname);
#else
  Test_LevelSetCut_Tesselation_and_DD_same_VC(nids, lsvs, xyze, testname);
#endif
}



void test_ls_hex8_magnus12()
{
  // simple hex8 element
  std::vector<int> nids(8);
  std::vector<double> lsvs(8);
  Epetra_SerialDenseMatrix xyze(3, 8);
  Epetra_SerialDenseMatrix xyze_local(3, 8);


  for (int i = 0; i < 8; ++i)
  {
    nids[i] = i;
  }

  xyze_local(0, 0) = -1;
  xyze_local(1, 0) = -1;
  xyze_local(2, 0) = -1;

  xyze_local(0, 1) = 1;
  xyze_local(1, 1) = -1;
  xyze_local(2, 1) = -1;

  xyze_local(0, 2) = 1;
  xyze_local(1, 2) = 1;
  xyze_local(2, 2) = -1;

  xyze_local(0, 3) = -1;
  xyze_local(1, 3) = 1;
  xyze_local(2, 3) = -1;

  xyze_local(0, 4) = -1;
  xyze_local(1, 4) = -1;
  xyze_local(2, 4) = 1;

  xyze_local(0, 5) = 1;
  xyze_local(1, 5) = -1;
  xyze_local(2, 5) = 1;

  xyze_local(0, 6) = 1;
  xyze_local(1, 6) = 0.5;
  xyze_local(2, 6) = 1;

  xyze_local(0, 7) = -1;
  xyze_local(1, 7) = 1;
  xyze_local(2, 7) = 1;


  // Test with local coord.
  //  //eleID=47018 (precision 16)
  //  xyze(  0,0 ) = 0.45;
  //  xyze(  1,0 ) = -0.2500000000000001;
  //  xyze(  2,0 ) = 0.09999999999999995;
  //
  //  xyze(  0,1 ) = 0.45;
  //  xyze(  1,1 ) = -0.25;
  //  xyze(  2,1 ) = 0.04999999999999986;
  //
  //  xyze(  0,2 ) = 0.4500000000000001;
  //  xyze(  1,2 ) = -0.2000000000000001;
  //  xyze(  2,2 ) = 0.04999999999999988;
  //
  //  xyze(  0,3 ) = 0.45;
  //  xyze(  1,3 ) = -0.2000000000000001;
  //  xyze(  2,3 ) = 0.09999999999999992;
  //
  //  xyze(  0,4 ) = 0.5;
  //  xyze(  1,4 ) = -0.25;
  //  xyze(  2,4 ) = 0.1;
  //
  //  xyze(  0,5 ) = 0.5;
  //  xyze(  1,5 ) = -0.25;
  //  xyze(  2,5 ) = 0.04999999999999992;
  //
  //  xyze(  0,6 ) = 0.5000000000000001;
  //  xyze(  1,6 ) = -0.2;
  //  xyze(  2,6 ) = 0.04999999999999993;
  //
  //  xyze(  0,7 ) = 0.5;
  //  xyze(  1,7 ) = -0.2000000000000001;
  //  xyze(  2,7 ) = 0.09999999999999998;

  //  lsvs[0] =0.02440442408507582;
  //  lsvs[1] =0.01720402163943002;
  //  lsvs[2] =-0.005025253169416655;
  //  lsvs[3] =0.002493781056044475;
  //  lsvs[4] =0.06789083458002731;
  //  lsvs[5] =0.06124860801609122;
  //  lsvs[6] =0.04083269131959844;
  //  lsvs[7] =0.04772255750516619;


  //  //eleID=47018 (precision 24)
  //  xyze(  0,0 ) = 0.449999999999999955591079;
  //  xyze(  1,0 ) = -0.250000000000000055511151;
  //  xyze(  2,0 ) = 0.0999999999999999500399639;
  //
  //  xyze(  0,1 ) = 0.449999999999999955591079;
  //  xyze(  1,1 ) = -0.25;
  //  xyze(  2,1 ) = 0.0499999999999998639976795;
  //
  //  xyze(  0,2 ) = 0.450000000000000066613381;
  //  xyze(  1,2 ) = -0.200000000000000066613381;
  //  xyze(  2,2 ) = 0.0499999999999998778754673;
  //
  //  xyze(  0,3 ) = 0.449999999999999955591079;
  //  xyze(  1,3 ) = -0.200000000000000066613381;
  //  xyze(  2,3 ) = 0.0999999999999999222843883;
  //
  //  xyze(  0,4 ) = 0.5;
  //  xyze(  1,4 ) = -0.25;
  //  xyze(  2,4 ) = 0.100000000000000005551115;
  //
  //  xyze(  0,5 ) = 0.5;
  //  xyze(  1,5 ) = -0.25;
  //  xyze(  2,5 ) = 0.0499999999999999195088307;
  //
  //  xyze(  0,6 ) = 0.500000000000000111022302;
  //  xyze(  1,6 ) = -0.199999999999999983346655;
  //  xyze(  2,6 ) = 0.0499999999999999333866185;
  //
  //  xyze(  0,7 ) = 0.5;
  //  xyze(  1,7 ) = -0.200000000000000066613381;
  //  xyze(  2,7 ) = 0.0999999999999999777955395;

  //  lsvs[0] =0.024404424085075815398227;
  //  lsvs[1] =0.0172040216394300227165104;
  //  lsvs[2] =-0.00502525316941665467496136;
  //  lsvs[3] =0.00249378105604447508625299;
  //  lsvs[4] =0.0678908345800273149706072;
  //  lsvs[5] =0.0612486080160912216285851;
  //  lsvs[6] =0.0408326913195984353421863;
  //  lsvs[7] =0.0477225575051661854431018;


  // 7: (0.5, -0.2, 0.1)
  // 0: (0.45, -0.25, 0.1)
  // 1: (0.45, -0.25, 0.05)

  double x_mid = 0.475;
  double y_mid = -0.225;
  double z_mid = 0.075;

  double delh = 0.05;

  // Stupid-projection to "local coordinates".
  Epetra_SerialDenseMatrix xyze_scaled(3, 8);

  for (int i = 0; i < 8; ++i)
  {
    xyze_scaled(0, i) = 2.0 * (xyze(0, i) - x_mid) / delh;
    xyze_scaled(1, i) = 2.0 * (xyze(1, i) - y_mid) / delh;
    xyze_scaled(2, i) = 2.0 * (xyze(2, i) - z_mid) / delh;
  }

  std::string testname("xxx_cut_test_ls_hex8_magnus12");

  // Only works for Tesselation as in this config the VC is undefined.
  // Cut_With_Tesselation(nids, lsvs, xyze_local, testname);

  //#ifdef LOCAL_COORDS
  //  Test_LevelSetCut_Tesselation_and_DD_same_VC(nids, lsvs, xyze_local, testname);
  //#else
  //  Test_LevelSetCut_Tesselation_and_DD_same_VC(nids, lsvs, xyze_local, testname);
  //#endif
}

// From variable surface tension problem in Combust!
void test_ls_hex8_magnus6()
{
  // simple hex8 element
  std::vector<int> nids(8);
  std::vector<double> lsvs(8);
  Epetra_SerialDenseMatrix xyze(3, 8);

  for (int i = 0; i < 8; ++i)
  {
    nids[i] = i;
  }

#ifdef PRECISION24
  std::cout << "Precision 24" << std::endl;
  //!!!!!!!!!!! volume test failed: !!!!!!!!!!!!!!!eleID=1  ve=9.86254e-13  vc=9.28731e-13
  //! vd= 5.75238e-14  err=0.0583255
  xyze(0, 0) = 0.00129103448275862077096465;
  xyze(1, 0) = -0.00069517241379310349205789;
  xyze(2, 0) = 5.00000000000000023960868e-05;

  xyze(0, 1) = 0.00129103448275862077096465;
  xyze(1, 1) = -0.00069517241379310349205789;
  xyze(2, 1) = -5.00000000000000023960868e-05;

  xyze(0, 2) = 0.00129103448275862077096465;
  xyze(1, 2) = -0.000595862068965517170486546;
  xyze(2, 2) = -5.00000000000000023960868e-05;

  xyze(0, 3) = 0.00129103448275862098780509;
  xyze(1, 3) = -0.000595862068965517062066328;
  xyze(2, 3) = 5.00000000000000023960868e-05;

  xyze(0, 4) = 0.00139034482758620720095621;
  xyze(1, 4) = -0.000695172413793103275217455;
  xyze(2, 4) = 5.00000000000000023960868e-05;

  xyze(0, 5) = 0.00139034482758620720095621;
  xyze(1, 5) = -0.000695172413793103383637673;
  xyze(2, 5) = -5.00000000000000023960868e-05;

  xyze(0, 6) = 0.00139034482758620741779665;
  xyze(1, 6) = -0.000595862068965517495747197;
  xyze(2, 6) = -5.00000000000000023960868e-05;

  xyze(0, 7) = 0.00139034482758620698411578;
  xyze(1, 7) = -0.000595862068965517170486546;
  xyze(2, 7) = 5.00000000000000023960868e-05;
  lsvs[0] = 2.68300282587983526902831e-05;
  lsvs[1] = 2.68300282587983323614924e-05;
  lsvs[2] = -1.77862543325451184349911e-05;
  lsvs[3] = -1.77862543325450438960918e-05;
  lsvs[4] = 0.000115213400540953834778196;
  lsvs[5] = 0.000115213400540954024513576;
  lsvs[6] = 7.35608797896675984783615e-05;
  lsvs[7] = 7.3560879789667029272221e-05;
#else
  std::cout << "Precision 16" << std::endl;
  xyze(0, 0) = -0.2999999999999999;
  xyze(1, 0) = -0.3;
  xyze(2, 0) = 0.3;

  xyze(0, 1) = -0.3;
  xyze(1, 1) = -0.3;
  xyze(2, 1) = 0.25;

  xyze(0, 2) = -0.3;
  xyze(1, 2) = -0.25;
  xyze(2, 2) = 0.25;

  xyze(0, 3) = -0.3;
  xyze(1, 3) = -0.25;
  xyze(2, 3) = 0.3;

  xyze(0, 4) = -0.25;
  xyze(1, 4) = -0.2999999999999999;
  xyze(2, 4) = 0.3;

  xyze(0, 5) = -0.25;
  xyze(1, 5) = -0.2999999999999999;
  xyze(2, 5) = 0.25;

  xyze(0, 6) = -0.25;
  xyze(1, 6) = -0.25;
  xyze(2, 6) = 0.25;

  xyze(0, 7) = -0.2500000000000001;
  xyze(1, 7) = -0.25;
  xyze(2, 7) = 0.3;

  lsvs[0] = 0.01961524227066325;
  lsvs[1] = -0.007557109910194759;
  lsvs[2] = -0.03631907522521477;
  lsvs[3] = -0.007557109910194759;
  lsvs[4] = -0.007557109910194815;
  lsvs[5] = -0.03631907522521483;
  lsvs[6] = -0.0669872981077807;
  lsvs[7] = -0.03631907522521477;
#endif

  std::string testname("xxx_cut_test_ls_hex8_magnus6");

#ifdef LOCAL_COORDS
  Test_LevelSetCut_Tesselation_and_DD_same_VC(nids, lsvs, GetLocalHex8Coords(), testname);
#else
  Test_LevelSetCut_Tesselation_and_DD(nids, lsvs, xyze, testname);
#endif
}

void test_ls_hex8_tes_dd_simple()
{
  // simple hex8 element
  std::vector<int> nids(8);
  std::vector<double> lsvs(8, -1);
  Epetra_SerialDenseMatrix xyze(3, 8);

  for (int i = 0; i < 8; ++i)
  {
    nids[i] = i;
  }

  // lsvs[2] = -.99981;

  //  lsvs[0] = 0.5;
  //  lsvs[1] = 0.43;
  //  lsvs[4] = -0.4123;
  //  lsvs[5] = -0.900091;

  //    lsvs[0] = 1.0;
  //    lsvs[3] = 0.0;
  //    lsvs[4] = 1.0;
  //    lsvs[7] = 0.0;
  //    lsvs[1] = 0.0;
  //    lsvs[5] = 0.0;


  lsvs[0] = 1.0;

  //  std::fill( &lsvs[0], &lsvs[4],   1. );
  //  std::fill( &lsvs[4], &lsvs[8],  -1. );

  xyze(0, 0) = 0;
  xyze(1, 0) = 0;
  xyze(2, 0) = 0;

  xyze(0, 1) = 1;
  xyze(1, 1) = 0;
  xyze(2, 1) = 0;

  xyze(0, 2) = 1;
  xyze(1, 2) = 1;
  xyze(2, 2) = 0;

  xyze(0, 3) = 0;
  xyze(1, 3) = 1;
  xyze(2, 3) = 0;

  xyze(0, 4) = 0;
  xyze(1, 4) = 0;
  xyze(2, 4) = 1;

  xyze(0, 5) = 1;
  xyze(1, 5) = 0;
  xyze(2, 5) = 1;

  xyze(0, 6) = 1;
  xyze(1, 6) = 1;
  xyze(2, 6) = 1;

  xyze(0, 7) = 0;
  xyze(1, 7) = 1;
  xyze(2, 7) = 1;

  std::string testname("xxx_cut_test_ls_hex8_magnus_simple");

#ifdef LOCAL_COORDS
  Test_LevelSetCut_Tesselation_and_DD_same_VC(nids, lsvs, GetLocalHex8Coords(), testname);
#else
  Test_LevelSetCut_Tesselation_and_DD_same_VC(nids, lsvs, xyze, testname);
#endif
}


void test_ls_mesh_hex8_simple()
{
  // simple hex8 element
  std::vector<int> nids(8);
  std::vector<double> lsvs(8, -1);
  Epetra_SerialDenseMatrix xyze(3, 8);


  GEO::CUT::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids_mesh;

  int sidecount = 0;
  {
    Epetra_SerialDenseMatrix tri3_xyze(3, 3);

    tri3_xyze(0, 2) = -1.0;
    tri3_xyze(1, 2) = 0.0;
    tri3_xyze(2, 2) = -1.0;
    tri3_xyze(0, 1) = -1.0;
    tri3_xyze(1, 1) = -1.0;
    tri3_xyze(2, 1) = 0.0;
    tri3_xyze(0, 0) = 0.0;
    tri3_xyze(1, 0) = -1.0;
    tri3_xyze(2, 0) = -1.0;
    nids_mesh.clear();
    nids_mesh.push_back(761);
    nids_mesh.push_back(857);
    nids_mesh.push_back(859);
    intersection.AddCutSide(++sidecount, nids_mesh, tri3_xyze, DRT::Element::tri3);
  }

  for (int i = 0; i < 8; ++i)
  {
    nids[i] = i;
  }

  // CUT WITH MESH
  intersection.AddElement(1, nids, GetLocalHex8Coords(), DRT::Element::hex8);

  intersection.Status();
  intersection.CutTest_Cut(true, INPAR::CUT::VCellGaussPts_Tessellation);

  std::vector<double> tessVol, dirDivVol;

  GEO::CUT::Mesh mesh = intersection.NormalMesh();
  const std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>& other_cells = mesh.VolumeCells();
  for (std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>::const_iterator i = other_cells.begin();
       i != other_cells.end(); ++i)
  {
    GEO::CUT::VolumeCell* vc = &**i;
    tessVol.push_back(vc->Volume());
  }

  int counter = 1;
  for (std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>::const_iterator i = other_cells.begin();
       i != other_cells.end(); ++i)
  {
    GEO::CUT::VolumeCell* vc = &**i;
    std::cout << "DIRCETDIVERGENCE CUTTING VOLUMECELL USING MESH #" << counter << ".....!"
              << std::endl;
    std::cout << "Volumecell Position: " << vc->Position() << std::endl;
    vc->DirectDivergenceGaussRule(
        vc->ParentElement(), mesh, true, INPAR::CUT::BCellGaussPts_Tessellation);
    std::cout << "DIRCETDIVERGENCE CUT VOLUMECELL USING MESH #" << counter << " WITHOUT ERROR!"
              << std::endl
              << std::endl;
    counter++;
    dirDivVol.push_back(vc->Volume());
  }

  // Test Direct Divergence against Tesselation
  //------------------------------------------------------
  bool error = false;
  std::cout << "the volumes predicted by\n Tessellation \t DirectDivergence\n";
  for (unsigned i = 0; i < tessVol.size(); i++)
  {
    std::cout << tessVol[i] << "\t" << dirDivVol[i] << "\n";
    if (fabs(tessVol[i] - dirDivVol[i]) > 1e-9) error = true;
  }
  if (error)
  {
    throw std::runtime_error("Volume predicted by one of the methods is wrong.");
    // dserror("volume predicted by either one of the method is wrong");
  }
  //------------------------------------------------------
  //###########################################################

  // CUT WITH LEVELSET in similar configuration.

  lsvs[0] = 1.0;

  std::string testname("xxx_cut_test_ls_hex8_magnus_ls_and_mesh");

#ifdef LOCAL_COORDS
  Test_LevelSetCut_Tesselation_and_DD_same_VC(nids, lsvs, GetLocalHex8Coords(), testname);
#else  // So far no Global Coord test.
  Test_LevelSetCut_Tesselation_and_DD_same_VC(nids, lsvs, GetLocalHex8Coords(), testname);
#endif
}

// See what happens
void test_ls_hex8_experiment_magnus()
{
  GEO::CUT::CombIntersection ci(-1);
  ci.GetOptions().Init_for_Cuttests();
  ci.AddLevelSetSide(1);

  // simple hex8 element
  std::vector<int> nids(8);
  std::vector<double> lsvs(8);
  Epetra_SerialDenseMatrix xyze(3, 8);

  for (int i = 0; i < 8; ++i)
  {
    nids[i] = i;
  }

  //  lsvs[0] =  1;
  //  lsvs[1] =  1;
  //  lsvs[2] =  -1;
  //  lsvs[3] =  -1;
  //  lsvs[4] =  1;
  //  lsvs[5] =  -1;
  //  lsvs[6] =  1;
  //  lsvs[7] =  -1;

  // TEST:
  lsvs[0] = 1;
  lsvs[1] = -1;
  lsvs[2] = 1;
  lsvs[3] = 1;
  lsvs[4] = 1;
  lsvs[5] = 1;
  lsvs[6] = 1;
  lsvs[7] = -1;

  xyze(0, 0) = -1;
  xyze(1, 0) = -1;
  xyze(2, 0) = -1;

  xyze(0, 1) = 1;
  xyze(1, 1) = -1;
  xyze(2, 1) = -1;

  xyze(0, 2) = 1;
  xyze(1, 2) = 1;
  xyze(2, 2) = -1;

  xyze(0, 3) = -1;
  xyze(1, 3) = 1;
  xyze(2, 3) = -1;

  xyze(0, 4) = -1;
  xyze(1, 4) = -1;
  xyze(2, 4) = 1;

  xyze(0, 5) = 1;
  xyze(1, 5) = -1;
  xyze(2, 5) = 1;

  xyze(0, 6) = 1;
  xyze(1, 6) = 1;
  xyze(2, 6) = 1;

  xyze(0, 7) = -1;
  xyze(1, 7) = 1;
  xyze(2, 7) = 1;

  std::string testname("xxx_cut_test_ls_hex8_experiment");

  // Test_LevelSetCut_Tesselation_and_DD(nids, lsvs, xyze, testname);


  ci.AddElement(1, nids, xyze, DRT::Element::hex8, &lsvs[0], false);

  ci.Cut(true);
  // ci.PrintCellStats();

  //  ci.NormalMesh().FindLSNodePositions();
  //  ci.NormalMesh().FindNodalDOFSets( true );

  ci.FindNodePositions();

  // Changed to DD
  ci.Cut_Finalize(true, INPAR::CUT::VCellGaussPts_Tessellation,
      INPAR::CUT::BCellGaussPts_Tessellation, false, true);
  // ci.Cut_Finalize( true, INPAR::CUT::VCellGaussPts_DirectDivergence,
  // INPAR::CUT::BCellGaussPts_Tessellation, false, true );


  //#ifdef GMSH_OUTPUT_LSNOLOC_CUT_TEST
  // Gmsh-output
  ci.NormalMesh().DumpGmsh("xxx_cut_test_ls_hex8_magnus6.CUT.pos");
  ci.NormalMesh().DumpGmshVolumeCells("xxx_cut_test_ls_hex8_magnus6.CUT_volumecells.pos", true);
  ci.DumpGmshIntegrationCells("xxx_cut_test_ls_hex8_magnus6.CUT_integrationcells.pos");
  //#endif

  const std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>& other_cells = ci.NormalMesh().VolumeCells();
  std::cout << "# Volume Cells: " << other_cells.size() << std::endl;

  int iteration_VC = 0;
  for (std::list<Teuchos::RCP<GEO::CUT::VolumeCell>>::const_iterator i = other_cells.begin();
       i != other_cells.end(); ++i)
  {
    iteration_VC++;
    std::cout << "VC(" << iteration_VC << "):" << std::endl;
    GEO::CUT::VolumeCell* vc = &**i;

    const GEO::CUT::plain_integrationcell_set& integrationcells = vc->IntegrationCells();
    std::cout << "Has #IC=" << integrationcells.size() << std::endl;

    if (integrationcells.size() == 0)
    {
      throw std::runtime_error("VolumeCell contains 0 integration cells.");
    }
  }
}

// This one will fail in GLOBAL DirectDiv configuration.
void test_ls_hex8_magnus7()
{
  std::vector<int> nids;

  std::vector<double> lsvs(8);
  lsvs[0] = 0.0280187;
  lsvs[1] = 0.0280187;
  lsvs[2] = 0.016173;
  lsvs[3] = 0.016173;
  lsvs[4] = 0.0123106;
  lsvs[5] = 0.0123106;
  lsvs[6] = -8.6e-09;
  lsvs[7] = -8.6e-09;

  Epetra_SerialDenseMatrix hex8_xyze(3, 8);

  nids.clear();
  hex8_xyze(0, 0) = -0.34;
  hex8_xyze(1, 0) = -0.26;
  hex8_xyze(2, 0) = 0.025;
  nids.push_back(0);
  hex8_xyze(0, 1) = -0.34;
  hex8_xyze(1, 1) = -0.26;
  hex8_xyze(2, 1) = -0.025;
  nids.push_back(1);
  hex8_xyze(0, 2) = -0.34;
  hex8_xyze(1, 2) = -0.24;
  hex8_xyze(2, 2) = -0.025;
  nids.push_back(2);
  hex8_xyze(0, 3) = -0.34;
  hex8_xyze(1, 3) = -0.24;
  hex8_xyze(2, 3) = 0.025;
  nids.push_back(3);
  hex8_xyze(0, 4) = -0.32;
  hex8_xyze(1, 4) = -0.26;
  hex8_xyze(2, 4) = 0.025;
  nids.push_back(4);
  hex8_xyze(0, 5) = -0.32;
  hex8_xyze(1, 5) = -0.26;
  hex8_xyze(2, 5) = -0.025;
  nids.push_back(5);
  hex8_xyze(0, 6) = -0.32;
  hex8_xyze(1, 6) = -0.24;
  hex8_xyze(2, 6) = -0.025;
  nids.push_back(6);
  hex8_xyze(0, 7) = -0.32;
  hex8_xyze(1, 7) = -0.24;
  hex8_xyze(2, 7) = 0.025;
  nids.push_back(7);

  std::string testname("xxx_cut_test_ls_hex8_magnus7");

  // Test_LevelSetCut_Tesselation_and_DD(nids, lsvs, hex8_xyze, testname);

  Cut_With_DirectDivergence(nids, lsvs, hex8_xyze, testname);
}

// void test_ls_hex8_magnus4()
//{
//
//  // FIX THIS TEST CASE!!!-> Cut normal way.
//
//  GEO::CUT::LevelSetIntersection lsi;
//
//  GEO::CUT::Options options;
//  lsi.GetOptions(options);
//  options.SetSimpleShapes( false );
//
//  // simple hex8 element
//  std::vector<int> nids( 8 );
//  std::vector<double> lsvs( 8 );
//  Epetra_SerialDenseMatrix xyze( 3, 8 );
//
//  for ( int i=0; i<8; ++i )
//  {
//    nids[i] = i;
//  }
//
////  WARNING: In eleId=37941 area mismatch: a1=0 a2=6.36112e-05 diff=-6.36112e-05
////   !!!!!!!!!!! volume test failed: !!!!!!!!!!!!!!!eleID=37941  ve=0.000125  vc=0.000124945
/// vd= 5.48009e-08  err=0.000438408
//  // ELEID: 37941 PREC: 16
//  //  xyze(  0,0 ) = 0.15;
//  //  xyze(  1,0 ) = 0.3999999999999999;
//  //  xyze(  2,0 ) = -0.05000000000000004;
//  //
//  //  xyze(  0,1 ) = 0.15;
//  //  xyze(  1,1 ) = 0.3999999999999999;
//  //  xyze(  2,1 ) = -0.1000000000000001;
//  //
//  //  xyze(  0,2 ) = 0.15;
//  //  xyze(  1,2 ) = 0.4500000000000002;
//  //  xyze(  2,2 ) = -0.1000000000000001;
//  //
//  //  xyze(  0,3 ) = 0.15;
//  //  xyze(  1,3 ) = 0.4500000000000002;
//  //  xyze(  2,3 ) = -0.05000000000000007;
//  //
//  //  xyze(  0,4 ) = 0.2000000000000001;
//  //  xyze(  1,4 ) = 0.3999999999999999;
//  //  xyze(  2,4 ) = -0.05000000000000004;
//  //
//  //  xyze(  0,5 ) = 0.2000000000000001;
//  //  xyze(  1,5 ) = 0.3999999999999999;
//  //  xyze(  2,5 ) = -0.1000000000000001;
//  //
//  //  xyze(  0,6 ) = 0.2000000000000001;
//  //  xyze(  1,6 ) = 0.4500000000000001;
//  //  xyze(  2,6 ) = -0.1000000000000001;
//  //
//  //  xyze(  0,7 ) = 0.2000000000000001;
//  //  xyze(  1,7 ) = 0.45;
//  //  xyze(  2,7 ) = -0.05000000000000007;
//  //
//  //  lsvs[0] =-0.06988373664786873;
//  //  lsvs[1] =-0.06125178063039394;
//  //  lsvs[2] =-0.0152320142583669;
//  //  lsvs[3] =-0.023030399291527;
//  //  lsvs[4] =-0.05000000000000004;
//  //  lsvs[5] =-0.04174243050441606;
//  //  lsvs[6] =0.002493781056044586;
//  //  lsvs[7] =-0.00502525316941671;
//
//  // ELEID: 37941 PREC:24
//  xyze(  0,0 ) = 0.15000000000000002220446;
//  xyze(  1,0 ) = 0.399999999999999911182158;
//  xyze(  2,0 ) = -0.050000000000000044408921;
//
//  xyze(  0,1 ) = 0.15000000000000002220446;
//  xyze(  1,1 ) = 0.399999999999999911182158;
//  xyze(  2,1 ) = -0.100000000000000116573418;
//
//  xyze(  0,2 ) = 0.150000000000000049960036;
//  xyze(  1,2 ) = 0.450000000000000177635684;
//  xyze(  2,2 ) = -0.100000000000000116573418;
//
//  xyze(  0,3 ) = 0.15000000000000002220446;
//  xyze(  1,3 ) = 0.450000000000000177635684;
//  xyze(  2,3 ) = -0.0500000000000000721644966;
//
//  xyze(  0,4 ) = 0.200000000000000066613381;
//  xyze(  1,4 ) = 0.399999999999999911182158;
//  xyze(  2,4 ) = -0.050000000000000044408921;
//
//  xyze(  0,5 ) = 0.200000000000000122124533;
//  xyze(  1,5 ) = 0.399999999999999911182158;
//  xyze(  2,5 ) = -0.100000000000000088817842;
//
//  xyze(  0,6 ) = 0.200000000000000066613381;
//  xyze(  1,6 ) = 0.450000000000000066613381;
//  xyze(  2,6 ) = -0.100000000000000088817842;
//
//  xyze(  0,7 ) = 0.200000000000000066613381;
//  xyze(  1,7 ) = 0.449999999999999955591079;
//  xyze(  2,7 ) = -0.0500000000000000721644966;
//
//  lsvs[0] =-0.0698837366478687327386865;
//  lsvs[1] =-0.0612517806303939371659339;
//  lsvs[2] =-0.0152320142583668993729873;
//  lsvs[3] =-0.0230303992915270039354425;
//  lsvs[4] =-0.050000000000000044408921;
//  lsvs[5] =-0.0417424305044160615096871;
//  lsvs[6] =0.00249378105604458610855545;
//  lsvs[7] =-0.0050252531694167101861126;
//
//  lsi.AddElement( 1, nids, xyze, DRT::Element::hex8, &lsvs[0]  );
//  lsi.Cutest_Cut(true,true);
//
//  lsi.PrintCellStats();
//  GEO::CUT::Mesh mesh = lsi.NormalMesh();
//
////  ci.Cut_Finalize( true, INPAR::CUT::VCellGaussPts_Tessellation,
/// INPAR::CUT::BCellGaussPts_Tessellation, false, true );
//
//  //Gmsh-output
//  mesh.DumpGmsh("xxx_cut_test_ls_hex8_magnus3.CUT.pos");
//  lsi.DumpGmshIntegrationCells("xxx_cut_test_ls_hex8_magnus3.CUT_integrationcells.pos");
//  mesh.DumpGmshVolumeCells("xxx_cut_test_ls_hex8_magnus3.CUT_volcells(mesh).pos",true);
//
//  const std::list<Teuchos::RCP<GEO::CUT::VolumeCell> > & other_cells = mesh.VolumeCells();
//  std::cout << "# Volume Cells: " << other_cells.size() << std::endl;
//  int iteration_VC = 0;
//  for ( std::list<Teuchos::RCP<GEO::CUT::VolumeCell> >::const_iterator i=other_cells.begin();
//      i!=other_cells.end();
//      ++i )
//  {
//    iteration_VC++;
//    std::cout << "#VC: " << iteration_VC << std::endl;
//    GEO::CUT::VolumeCell * vc = &**i;
//
//    const GEO::CUT::plain_integrationcell_set & integrationcells = vc->IntegrationCells();
//    std::cout << "Has #IC=" << integrationcells.size() << std::endl;
//    //vc->Print(std::cout);
//    const GEO::CUT::plain_facet_set & facete = vc->Facets();
//    // check whether all facets of this vc are oriented in a plane
//    // if not then some sides are warped
//    // we need to generate quadrature rules in global coordinates
//    for(GEO::CUT::plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
//    {
//      GEO::CUT::Facet *fe = *i;
//      std::vector<GEO::CUT::Point*> corn = fe->CornerPoints();
//
//      fe->Print();
//
//      bool isPlanar = fe->IsPlanar( mesh, corn );
//      if ( isPlanar == false )
//      {
//        std::cout << "WARNING: Facet is NOT planar!!!" << std::endl;
//        //dserror("VolumeCell is not planar!!");
//      }
//    }
//
//    std::cout<<std::endl;
//    //tessVol.push_back(vc->Volume());
//  }
//}

/*----------------------------------------------------------------------------*/
/** \file
\brief compute the intersection of the zero level-set iso-contour on the
  SCATRA surface discretization and evaluate important quantities


\level 3
*/
/*----------------------------------------------------------------------------*/


#include "xcontact_levelset_intersection.H"

#include "../drt_cut/cut_boundarycell.H"
#include "../drt_cut/cut_position.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_io/io_pstream.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XCONTACT::LEVELSET::Intersection::Intersection() : SCATRA::LEVELSET::Intersection()
{
  check_lsv_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::Intersection::CaptureZeroLevelSet(const Epetra_Vector& phi,
    const DRT::Discretization& scatradis,
    std::map<int, GEO::BoundaryIntCellPtrs>& elementBoundaryIntCells)
{
  // reset class members
  Reset();

  // herein the actual capturing happens
  GetZeroLevelSet(phi, scatradis, elementBoundaryIntCells, false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::Intersection::CaptureZeroLevelSet(const Epetra_Vector& phi,
    const DRT::Discretization& scatradis, double& volumedomainminus, double& volumedomainplus,
    double& zerosurface, std::map<int, GEO::BoundaryIntCellPtrs>& elementBoundaryIntCells)
{
  CaptureZeroLevelSet(phi, scatradis, elementBoundaryIntCells);

  // collect contributions from all procs and store in respective variables
  scatradis.Comm().SumAll(&VolumePlus(), &volumedomainplus, 1);
  scatradis.Comm().SumAll(&VolumeMinus(), &volumedomainminus, 1);
  scatradis.Comm().SumAll(&Surface(), &zerosurface, 1);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::Intersection::CheckBoundaryCellType(
    DRT::Element::DiscretizationType distype_bc) const
{
  if (distype_bc != DRT::Element::point1 and distype_bc != DRT::Element::line2)
  {
    dserror("unexpected type of boundary integration cell: %s",
        DRT::DistypeToString(distype_bc).c_str());
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::Intersection::AddToBoundaryIntCellsPerEle(
    const LINALG::SerialDenseMatrix& xyze, const GEO::CUT::BoundaryCell& bcell,
    DRT::Element::DiscretizationType distype_ele)
{
  DRT::Element::DiscretizationType distype_bc = bcell.Shape();
  CheckBoundaryCellType(distype_bc);

  const int numnodebc = DRT::UTILS::getNumberOfElementNodes(distype_bc);

  // get physical coordinates of this cell
  LINALG::SerialDenseMatrix coord = bcell.Coordinates();

  // transfer to element coordinates
  LINALG::SerialDenseMatrix localcoord(3, numnodebc, true);

  for (int ivertex = 0; ivertex < numnodebc; ivertex++)
  {
    LINALG::Matrix<3, 1> pcoord(&coord(0, ivertex), true);
    if (pcoord.M() != static_cast<unsigned>(coord.M())) dserror("row dimension mismatch!");

    Teuchos::RCP<GEO::CUT::Position> pos = GEO::CUT::Position::Create(xyze, pcoord, distype_ele);
    if (pos->Compute())
    {
      LINALG::Matrix<3, 1> lcoord(&localcoord(0, ivertex), true);
      pos->LocalCoordinates(lcoord);
    }
    else
      dserror("The local position calculation failed!");
  }


  std::cout << "--- local coordinates ---\n";
  std::cout << localcoord << std::endl;

  // store boundary element and sum area into surface
  // be careful, we only set physical coordinates
  BoundaryIntCellsPerEle<GEO::BoundaryIntCellPtrs>().push_back(
      Teuchos::rcp(GEO::BoundaryIntCell::Create(distype_bc, -1, localcoord, NULL, coord, true)));
}

/*----------------------------------------------------------------------------*/
/** \file

  \brief

  allow for computing intersection of zero level-set iso-contour with discretization
  and related quantities, i.g., volume of subdomains, interface discretization, ...


  \level 2


*/
/*----------------------------------------------------------------------------*/


#ifndef BACI_LEVELSET_INTERSECTION_UTILS_H
#define BACI_LEVELSET_INTERSECTION_UTILS_H

#include "baci_cut_point.H"
#include "baci_discretization_geometry_geo_utils.H"

#include <Epetra_MpiComm.h>
#include <Epetra_MultiVector.h>
#include <Teuchos_RCP.hpp>

namespace DRT
{
  class Discretization;
  class PackBuffer;
}  // namespace DRT

namespace CORE::GEO
{
  namespace CUT
  {
    class ElementHandle;
    class LevelSetIntersection;
  }  // namespace CUT
}  // namespace CORE::GEO

namespace SCATRA
{
  namespace LEVELSET
  {
    /** \brief level-set intersection utils class
     *
     *  Level-Set intersection functions wrapped in a class, thus inheritance
     *  becomes possible.
     *
     *  \author hiermeier \date 11/16 */
    class Intersection
    {
     public:
      /// constructor
      Intersection();

      /// destructor
      virtual ~Intersection() = default;

      /** \brief construct zero iso-contour of level-set field
       *
       *  \author rasthofer \date 09/13 */
      void CaptureZeroLevelSet(const Teuchos::RCP<const Epetra_Vector>& phi,
          const Teuchos::RCP<const DRT::Discretization>& scatradis, double& volumedomainminus,
          double& volumedomainplus, double& zerosurface,
          std::map<int, CORE::GEO::BoundaryIntCells>& elementBoundaryIntCells);

      /** \brief Set desired positions
       *
       *  \param desired_pos (in) : vector containing desired domain positions
       *                            ( i.e. inside, outside )
       *
       *  We will extract the boundary cells from the volume cells corresponding
       *  to the here defined positions. If no position vector is given, the
       *  outside domain will be considered.
       *
       *  \author hiermeier \date 11/16 */
      void SetDesiredPositions(
          const std::vector<CORE::GEO::CUT::Point::PointPosition>& desired_pos);

     protected:
      /// reset class member variables
      void Reset();

      template <typename T>
      void GetZeroLevelSet(const Epetra_Vector& phi, const DRT::Discretization& scatradis,
          std::map<int, T>& elementBoundaryIntCells, bool cut_screenoutput = false);

      /** \brief export boundary integration cells from this proc to parallel distribution
       *
       * \author henke \date 12/09 */
      void ExportInterface(
          std::map<int, CORE::GEO::BoundaryIntCells>& myinterface, const Epetra_Comm& comm);

      /** \brief pack boundary integration cells from set into char array
       *
       *  \author henke \date 12/09 */
      void packBoundaryIntCells(
          const std::map<int, CORE::GEO::BoundaryIntCells>& intcellmap, DRT::PackBuffer& dataSend);

      /** brief unpack boundary integration cells from char array
       *
       * \author henke \date 12/09 */
      void unpackBoundaryIntCells(const std::vector<char>& dataRecv,
          std::map<int, CORE::GEO::BoundaryIntCells>& intcellmap);

      /// return the volume of the plus domain
      inline double& VolumePlus() { return volumeplus_; };

      /// return the volume of the minus domain
      inline double& VolumeMinus() { return volumeminus_; };

      /** \brief add volume corresponding to the given PointPosition
       *
       *  Small inconsistency in the name convention:
       *
       *      outside --> plus domain
       *      inside  --> minus domain
       *
       *  \author hiermeier \date 11/16 */
      void AddToVolume(CORE::GEO::CUT::Point::PointPosition pos, double vol);

      /// access the boundary cell surface value
      inline double& Surface() { return surface_; };

      /** \brief prepare the cut algorithm
       *
       *  \author hiermeier \date 11/16 */
      void PrepareCut(const DRT::Element* ele, const DRT::Discretization& scatradis,
          const Epetra_Vector& phicol, CORE::LINALG::SerialDenseMatrix& xyze,
          std::vector<double>& phi_nodes, std::vector<int>& node_ids) const;

      /// perform the cut operation
      CORE::GEO::CUT::ElementHandle* Cut(CORE::GEO::CUT::LevelSetIntersection& levelset,
          const CORE::LINALG::SerialDenseMatrix& xyze, const std::vector<double>& phi_nodes,
          bool cut_screenoutput) const;

      /// collect the cut elements after a successful cut operation
      void CollectCutEles(CORE::GEO::CUT::ElementHandle& ehandle,
          CORE::GEO::CUT::plain_element_set& cuteles,
          DRT::Element::DiscretizationType distype) const;

      /** \brief check the point position (OR-combination)
       *
       *  \param curr_pos (in) : current position of the volume cell
       *
       *  \author hiermeier \date 11/16 */
      bool IsPointPosition(const CORE::GEO::CUT::Point::PointPosition& curr_pos)
      {
        return IsPointPosition(curr_pos, DesiredPositions());
      }
      bool IsPointPosition(const CORE::GEO::CUT::Point::PointPosition& curr_pos,
          const std::vector<CORE::GEO::CUT::Point::PointPosition>& desired_pos) const;

      /** \brief get the zero level-set
       *
       *  \author rasthofer \date 09/13 */
      void GetZeroLevelSetContour(const CORE::GEO::CUT::plain_element_set& cuteles,
          const CORE::LINALG::SerialDenseMatrix& xyze, DRT::Element::DiscretizationType distype);

      /// check for supported boundary cell discretization types
      virtual void CheckBoundaryCellType(DRT::Element::DiscretizationType distype_bc) const;

      virtual void AddToBoundaryIntCellsPerEle(const CORE::LINALG::SerialDenseMatrix& xyze,
          const CORE::GEO::CUT::BoundaryCell& bcell, DRT::Element::DiscretizationType distype_ele);

      /// access the private boundary cell vector
      template <typename T>
      T& BoundaryIntCellsPerEle();

      const std::vector<CORE::GEO::CUT::Point::PointPosition>& DesiredPositions();

     protected:
      /** check the level set values before we add a new element to the
       *  CORE::GEO::CUT::LevelSetIntersection object */
      bool check_lsv_;

      /// vector containing the desired positions ( default: outside )
      std::vector<CORE::GEO::CUT::Point::PointPosition> desired_positions_;

     private:
      /// boundary cell vector
      CORE::GEO::BoundaryIntCells listBoundaryIntCellsperEle_;

      // boundary cell pointer vector
      CORE::GEO::BoundaryIntCellPtrs boundary_cells_per_ele_;

      /// accumulated value of the plus domain volumes (POSITION == outside)
      double volumeplus_;

      /// accumulated value of the minus domain volumes (POSITION == inside)
      double volumeminus_;

      /// accumulated value of the boundary cell surfaces
      double surface_;
    };  // class Intersection

    /*----------------------------------------------------------------------------*/
    template <>
    inline CORE::GEO::BoundaryIntCells&
    Intersection::BoundaryIntCellsPerEle<CORE::GEO::BoundaryIntCells>()
    {
      return listBoundaryIntCellsperEle_;
    }

    /*----------------------------------------------------------------------------*/
    template <>
    inline CORE::GEO::BoundaryIntCellPtrs&
    Intersection::BoundaryIntCellsPerEle<CORE::GEO::BoundaryIntCellPtrs>()
    {
      return boundary_cells_per_ele_;
    }

  }  // namespace LEVELSET
}  // namespace SCATRA


#endif  // LEVELSET_INTERSECTION_UTILS_H

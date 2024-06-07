/*----------------------------------------------------------------------------*/
/** \file

  \brief

  allow for computing intersection of zero level-set iso-contour with discretization
  and related quantities, i.g., volume of subdomains, interface discretization, ...


  \level 2


*/
/*----------------------------------------------------------------------------*/


#ifndef FOUR_C_LEVELSET_INTERSECTION_UTILS_HPP
#define FOUR_C_LEVELSET_INTERSECTION_UTILS_HPP

#include "4C_config.hpp"

#include "4C_cut_point.hpp"
#include "4C_fem_geometry_geo_utils.hpp"

#include <Epetra_MpiComm.h>
#include <Epetra_MultiVector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::Communication
{
  class PackBuffer;
}

namespace Discret
{
  class Discretization;
}  // namespace Discret

namespace Core::Geo
{
  namespace Cut
  {
    class ElementHandle;
    class LevelSetIntersection;
  }  // namespace Cut
}  // namespace Core::Geo

namespace ScaTra
{
  namespace LevelSet
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
          const Teuchos::RCP<const Discret::Discretization>& scatradis, double& volumedomainminus,
          double& volumedomainplus, double& zerosurface,
          std::map<int, Core::Geo::BoundaryIntCells>& elementBoundaryIntCells);

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
          const std::vector<Core::Geo::Cut::Point::PointPosition>& desired_pos);

     protected:
      /// reset class member variables
      void reset();

      template <typename T>
      void get_zero_level_set(const Epetra_Vector& phi, const Discret::Discretization& scatradis,
          std::map<int, T>& elementBoundaryIntCells, bool cut_screenoutput = false);

      /** \brief export boundary integration cells from this proc to parallel distribution
       *
       * \author henke \date 12/09 */
      void export_interface(
          std::map<int, Core::Geo::BoundaryIntCells>& myinterface, const Epetra_Comm& comm);

      /** \brief pack boundary integration cells from set into char array
       *
       *  \author henke \date 12/09 */
      void pack_boundary_int_cells(const std::map<int, Core::Geo::BoundaryIntCells>& intcellmap,
          Core::Communication::PackBuffer& dataSend);

      /** brief unpack boundary integration cells from char array
       *
       * \author henke \date 12/09 */
      void unpack_boundary_int_cells(const std::vector<char>& dataRecv,
          std::map<int, Core::Geo::BoundaryIntCells>& intcellmap);

      /// return the volume of the plus domain
      inline double& volume_plus() { return volumeplus_; };

      /// return the volume of the minus domain
      inline double& volume_minus() { return volumeminus_; };

      /** \brief add volume corresponding to the given PointPosition
       *
       *  Small inconsistency in the name convention:
       *
       *      outside --> plus domain
       *      inside  --> minus domain
       *
       *  \author hiermeier \date 11/16 */
      void add_to_volume(Core::Geo::Cut::Point::PointPosition pos, double vol);

      /// access the boundary cell surface value
      inline double& surface() { return surface_; };

      /** \brief prepare the cut algorithm
       *
       *  \author hiermeier \date 11/16 */
      void prepare_cut(const Core::Elements::Element* ele, const Discret::Discretization& scatradis,
          const Epetra_Vector& phicol, Core::LinAlg::SerialDenseMatrix& xyze,
          std::vector<double>& phi_nodes, std::vector<int>& node_ids) const;

      /// perform the cut operation
      Core::Geo::Cut::ElementHandle* cut(Core::Geo::Cut::LevelSetIntersection& levelset,
          const Core::LinAlg::SerialDenseMatrix& xyze, const std::vector<double>& phi_nodes,
          bool cut_screenoutput) const;

      /// collect the cut elements after a successful cut operation
      void collect_cut_eles(Core::Geo::Cut::ElementHandle& ehandle,
          Core::Geo::Cut::plain_element_set& cuteles, Core::FE::CellType distype) const;

      /** \brief check the point position (OR-combination)
       *
       *  \param curr_pos (in) : current position of the volume cell
       *
       *  \author hiermeier \date 11/16 */
      bool is_point_position(const Core::Geo::Cut::Point::PointPosition& curr_pos)
      {
        return is_point_position(curr_pos, desired_positions());
      }
      bool is_point_position(const Core::Geo::Cut::Point::PointPosition& curr_pos,
          const std::vector<Core::Geo::Cut::Point::PointPosition>& desired_pos) const;

      /** \brief get the zero level-set
       *
       *  \author rasthofer \date 09/13 */
      void get_zero_level_set_contour(const Core::Geo::Cut::plain_element_set& cuteles,
          const Core::LinAlg::SerialDenseMatrix& xyze, Core::FE::CellType distype);

      /// check for supported boundary cell discretization types
      virtual void check_boundary_cell_type(Core::FE::CellType distype_bc) const;

      virtual void add_to_boundary_int_cells_per_ele(const Core::LinAlg::SerialDenseMatrix& xyze,
          const Core::Geo::Cut::BoundaryCell& bcell, Core::FE::CellType distype_ele);

      /// access the private boundary cell vector
      template <typename T>
      T& boundary_int_cells_per_ele();

      const std::vector<Core::Geo::Cut::Point::PointPosition>& desired_positions();

     protected:
      /** check the level set values before we add a new element to the
       *  Core::Geo::Cut::LevelSetIntersection object */
      bool check_lsv_;

      /// vector containing the desired positions ( default: outside )
      std::vector<Core::Geo::Cut::Point::PointPosition> desired_positions_;

     private:
      /// boundary cell vector
      Core::Geo::BoundaryIntCells list_boundary_int_cellsper_ele_;

      // boundary cell pointer vector
      Core::Geo::BoundaryIntCellPtrs boundary_cells_per_ele_;

      /// accumulated value of the plus domain volumes (POSITION == outside)
      double volumeplus_;

      /// accumulated value of the minus domain volumes (POSITION == inside)
      double volumeminus_;

      /// accumulated value of the boundary cell surfaces
      double surface_;
    };  // class intersection

    /*----------------------------------------------------------------------------*/
    template <>
    inline Core::Geo::BoundaryIntCells&
    Intersection::boundary_int_cells_per_ele<Core::Geo::BoundaryIntCells>()
    {
      return list_boundary_int_cellsper_ele_;
    }

    /*----------------------------------------------------------------------------*/
    template <>
    inline Core::Geo::BoundaryIntCellPtrs&
    Intersection::boundary_int_cells_per_ele<Core::Geo::BoundaryIntCellPtrs>()
    {
      return boundary_cells_per_ele_;
    }

  }  // namespace LevelSet
}  // namespace ScaTra


FOUR_C_NAMESPACE_CLOSE

#endif

/*---------------------------------------------------------------------*/
/*! \file

\brief Create and handle integrationcells for the Tessellation routine

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_INTEGRATIONCELLCREATOR_HPP
#define FOUR_C_CUT_INTEGRATIONCELLCREATOR_HPP

#include "4C_config.hpp"

#include "4C_cut_volumecell.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  namespace Cut
  {
    class Mesh;

    /*!
    \brief Special cases for integration cell creation: a library of cuts.

    Some cuts lead to known shapes those meshes are known. Here we collect
    these cases. Feel free to add your own.
   */
    class IntegrationCellCreator
    {
     public:
      static bool create_cells(Mesh& mesh, Element* element, const plain_volumecell_set& cells);

      static bool create_cell(Mesh& mesh, Core::FE::CellType shape, VolumeCell* cell);

     private:
      /** \brief loop over all volume cells and initiate the volume and boundary integration cell
       *  creation process
       *
       *  This is done after a successful initiation of the volume cells.*/
      void execute(Mesh& mesh);

      /// fill the point1 boundary cell
      bool create_point1_cell(Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets);

      /// fill the line2 volume cell
      bool create_line2_cell(Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets);

      /** \brief fill the 2-D volume cells
       *
       *  \author hiermeier \date 01/17 */
      template <Core::FE::CellType celltype,
          Core::FE::CellType facetype = Core::FE::DisTypeToFaceShapeType<celltype>::shape,
          unsigned numfaces = Core::FE::num_faces<celltype>>
      bool create_2d_cell(Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets);

      /// fill the tet4 volume cell
      bool create_tet4_cell(Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets);

      /// fill the hex8 volume cell
      bool create_hex8_cell(Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets);

      /// fill the wedge6 volume cell
      bool create_wedge6_cell(Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets);

      /// fill the pyramid5 volume cell
      bool create_pyramid5_cell(Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets);

      bool create_special_cases(Mesh& mesh, VolumeCell* cell, const plain_facet_set& facets);

      bool hex8_horizontal_cut(Mesh& mesh, Element* element, VolumeCell* cell,
          const plain_facet_set& facets, int axis, double r);

      /// add the volume cell information during the Create<Shape>Cell calls
      void add(VolumeCell* vc, Core::FE::CellType shape, const std::vector<Point*>& points)
      {
        Volume& v = cells_[vc];
        std::vector<Ic>& cells = v.domain_;
        cells.push_back(Ic());
        Ic& cell = cells.back();
        cell.shape_ = shape;
        cell.points_.reserve(points.size());
        cell.points_.assign(points.begin(), points.end());
      }

      /** \brief add the side (boundary cell) information during Create<Shape>Cell calls
       *         for the desired boundary cell positions
       *
       *  \author hiermeier \date 01/17   */
      void add_side(Core::Geo::Cut::BoundaryCellPosition bcell_position, VolumeCell* vc,
          Facet* facet, Core::FE::CellType shape, const std::vector<Point*>& side);

      /// add the side (boundary cell) information during the Create<Shape>Cell calls
      void add_side(
          VolumeCell* vc, Facet* facet, Core::FE::CellType shape, const std::vector<Point*>& side)
      {
        Volume& cell = cells_[vc];
        std::vector<Bc>& bcells = cell.boundary_;
        bcells.push_back(Bc());
        Bc& bcell = bcells.back();
        bcell.shape_ = shape;
        bcell.facet_ = facet;
        bcell.side_.reserve(side.size());
        bcell.side_.assign(side.begin(), side.end());
      }

      /*-------------------------------------------------------------------------*/
      /// construction of a boundary integration cell
      struct Bc
      {
        Core::FE::CellType shape_;
        std::vector<Point*> side_;
        Facet* facet_;

        /// actual creation of the desired boundary integration cell
        void execute(Mesh& mesh, VolumeCell* vc)
        {
          vc->new_boundary_cell(mesh, shape_, facet_, side_);
        }
      };

      /*-------------------------------------------------------------------------*/
      /// construction of a volume integration cell
      struct Ic
      {
        Core::FE::CellType shape_;
        std::vector<Point*> points_;

        /// actual creation of the desired volume integration cell
        void execute(Mesh& mesh, VolumeCell* vc)
        {
          vc->new_integration_cell(mesh, shape_, points_);
        }
      };

      /*-------------------------------------------------------------------------*/
      /** \brief Create all volume and boundary integration cells                */
      struct Volume
      {
        std::vector<Ic> domain_;
        std::vector<Bc> boundary_;

        /** \brief This routine initiates the volume and boundary cell creation.
         *
         *  At this point the necessary information for the creation must have been
         *  already added by the Add() and add_side() methods. */
        void execute(Mesh& mesh, VolumeCell* vc)
        {
          for (std::vector<Ic>::iterator i = domain_.begin(); i != domain_.end(); ++i)
          {
            Ic& cell = *i;
            cell.execute(mesh, vc);
          }
          for (std::vector<Bc>::iterator i = boundary_.begin(); i != boundary_.end(); ++i)
          {
            Bc& bcell = *i;
            bcell.execute(mesh, vc);
          }
        }
      };  // struct volume

      std::map<VolumeCell*, Volume> cells_;
    };

  }  // namespace Cut
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif

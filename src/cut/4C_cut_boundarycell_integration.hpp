/*---------------------------------------------------------------------*/
/*! \file

\brief used in boundary cell integration

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_BOUNDARYCELL_INTEGRATION_HPP
#define FOUR_C_CUT_BOUNDARYCELL_INTEGRATION_HPP

#include "4C_config.hpp"

#include "4C_cut_facet_integration.hpp"
#include "4C_linalg_serialdensematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CORE::GEO
{
  namespace CUT
  {
    class Element;
    class Facet;

    /*
    \brief This class integrates the base functions used in the moment fitting equations over the
    cut facets of the volumecell
     */
    class BoundarycellIntegration
    {
     public:
      /*!
       \brief Constructor
       */
      BoundarycellIntegration(Element* elem, Facet* bcell,
          const CORE::GEO::CUT::Point::PointPosition posi, int num_func)
          : elem1_(elem), bcell_(bcell), position_(posi), num_func_(num_func)
      {
      }

      /*!
      \brief Generate integration rule for the considered boundarycell.
      Unlike facet integration facets, whose x-component of normal is zero, cannot be eliminated
      from the integration.
       */
      CORE::LINALG::SerialDenseVector generate_boundary_cell_integration_rule();

      /*!
      \brief Returns the location of Gauss points over the boundarycell
       */
      std::vector<std::vector<double>> get_bcell_gauss_point_location() { return bcellgaus_pts_; }

     private:
      Element* elem1_;
      Facet* bcell_;
      const CORE::GEO::CUT::Point::PointPosition position_;
      int num_func_;
      std::vector<std::vector<double>> bcellgaus_pts_;

      /*!
      \brief Distribute the Gauss points over the boundarycell
       */
      void distribute_boundary_cell_gauss_points(std::vector<double> eqn,
          std::vector<std::vector<double>> corners, std::vector<std::vector<double>>& bcGausspts,
          int ptNos);

      /*!
      \brief Moment fitting matrix is formed for the boundarycell integration
      */
      void moment_fitting_matrix(
          std::vector<std::vector<double>>& mom, std::vector<std::vector<double>> gauspts);

      /*!
      \brief The geometry of boundarycell and the location of Gaussian points are written in GMSH
      output file for visualization and for debugging
      */
      void bcell_gauss_point_gmsh(const std::vector<std::vector<double>> bcGausspts,
          const std::vector<std::vector<double>> corners);
    };
  }  // namespace CUT
}  // namespace CORE::GEO

FOUR_C_NAMESPACE_CLOSE

#endif

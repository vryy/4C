/*----------------------------------------------------------------------*/
/*! \file

\brief service methods for accessing knot vector and weights for a given
       element of a nurbs discretisations

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_FEM_NURBS_DISCRETIZATION_UTILS_HPP
#define FOUR_C_FEM_NURBS_DISCRETIZATION_UTILS_HPP

#include "4C_config.hpp"

#include "4C_fem_nurbs_discretization.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  namespace Nurbs
  {
    /*----------------------------------------------------------------------*/
    /*!
    \brief A service method for accessing knotvector and weights for
           an isogeometric element


    \param dis         (i) the discretisation
    \param ele         (i) a pointer to the element
    \param myknots     (o) knot vector (to be filled)
    \param weights     (o) weight vector (to be filled)

    \date 12/10
    */
    template <class WG>
    bool GetMyNurbsKnotsAndWeights(const Core::FE::Discretization& dis,
        const Core::Elements::Element* ele, std::vector<Core::LinAlg::SerialDenseVector>& myknots,
        WG& weights)
    {
      // try to cast dis to a nurbs discretisation
      const Core::FE::Nurbs::NurbsDiscretization* nurbsdis =
          dynamic_cast<const Core::FE::Nurbs::NurbsDiscretization*>(&(dis));
      if (nurbsdis == nullptr) FOUR_C_THROW("Received discretization which is not Nurbs!");

      // get local knot vector entries and check for zero sized elements
      const bool zero_size = (*((*nurbsdis).get_knot_vector())).get_ele_knots(myknots, ele->id());

      // if we have a zero sized element due to a interpolated
      // point --- exit here and tell the outside world about that
      if (zero_size)
      {
        return (zero_size);
      }
      // you are still here? So get the node weights for the nurbs element as well
      const Core::Nodes::Node* const* nodes = ele->nodes();
      const int nen = ele->num_node();
      for (int inode = 0; inode < nen; inode++)
      {
        const Core::FE::Nurbs::ControlPoint* cp =
            dynamic_cast<const Core::FE::Nurbs::ControlPoint*>(nodes[inode]);
        weights(inode) = cp->w();
      }

      // goodbye
      return zero_size;
    }  // GetMyNurbsKnotsAndWeights()


    //! determine whether a given element is a NURBS element or not
    inline bool IsNurbs(Core::FE::CellType distype)
    {
      switch (distype)
      {
        case Core::FE::CellType::nurbs8:
        case Core::FE::CellType::nurbs27:
        case Core::FE::CellType::nurbs4:
        case Core::FE::CellType::nurbs9:
        case Core::FE::CellType::nurbs2:
        case Core::FE::CellType::nurbs3:
        {
          return true;
          break;
        }
        default:
          return false;
      }
    };


    /*!
    \brief A service method for accessing knotvector and weights for
           an isogeometric boundary element


    \param boundaryele      (i) a pointer to the boundary element
    \param localsurfaceid   (i) local id of this boundary element
    \param parenteleid      (i) global id of parent element
    \param dis              (i) the discretisation
    \param myknots          (o) parent knot vector (to be filled)
    \param myknots          (o) knot vector for boundary element (to be filled)
    \param weights          (o) weight vector for boundary element (to be filled)
    \param normalfac        (o) normalfac (to be filled)

    \date 12/10
    */
    template <class WG>
    bool GetKnotVectorAndWeightsForNurbsBoundary(const Core::Elements::Element* boundaryele,
        const int localsurfaceid, const int parenteleid,
        const Core::FE::Discretization& discretization,
        std::vector<Core::LinAlg::SerialDenseVector>& mypknots,
        std::vector<Core::LinAlg::SerialDenseVector>& myknots, WG& weights, double& normalfac)
    {
      // get knotvector(s)
      const Core::FE::Nurbs::NurbsDiscretization* nurbsdis =
          dynamic_cast<const Core::FE::Nurbs::NurbsDiscretization*>(&(discretization));

      Teuchos::RCP<const Core::FE::Nurbs::Knotvector> knots = (*nurbsdis).get_knot_vector();

      bool zero_size = knots->get_boundary_ele_and_parent_knots(
          mypknots, myknots, normalfac, parenteleid, localsurfaceid);

      // if we have a zero sized element due to a interpolated
      // point --- exit here and tell the outside world about that
      if (zero_size)
      {
        return (zero_size);
      }
      // you are still here? So get the node weights as well
      const Core::Nodes::Node* const* nodes = boundaryele->nodes();
      const int boundarynen = boundaryele->num_node();
      for (int inode = 0; inode < boundarynen; inode++)
      {
        const Core::FE::Nurbs::ControlPoint* cp =
            dynamic_cast<const Core::FE::Nurbs::ControlPoint*>(nodes[inode]);
        weights(inode) = cp->w();
      }

      // goodbye
      return zero_size;
    }  // GetKnotVectorAndWeightsForNurbsBoundary()

    /*!
    \brief A service method for accessing knotvector and weights for
           an isogeometric boundary element and parent element


    \param boundaryele      (i) a pointer to the boundary element
    \param localsurfaceid   (i) local id of this boundary element
    \param parenteleid      (i) global id of parent element
    \param dis              (i) the discretisation
    \param pmyknots         (o) parent knot vector (to be filled)
    \param myknots          (o) knot vector for boundary element (to be filled)
    \param pweights         (o) weight vector for parent element (to be filled)
    \param weights          (o) weight vector for boundary element (to be filled)
    \param normalfac        (o) normalfac (to be filled)

    \date 01/14
    */
    template <class WG>
    bool GetKnotVectorAndWeightsForNurbsBoundaryAndParent(Core::Elements::Element* parentele,
        Core::Elements::Element* boundaryele, const int localsurfaceid,
        const Core::FE::Discretization& discretization,
        std::vector<Core::LinAlg::SerialDenseVector>& mypknots,
        std::vector<Core::LinAlg::SerialDenseVector>& myknots, WG& pweights, WG& weights,
        double& normalfac)
    {
      // get knotvector(s)
      const Core::FE::Nurbs::NurbsDiscretization* nurbsdis =
          dynamic_cast<const Core::FE::Nurbs::NurbsDiscretization*>(&(discretization));

      Teuchos::RCP<const Core::FE::Nurbs::Knotvector> knots = (*nurbsdis).get_knot_vector();

      bool zero_size = knots->get_boundary_ele_and_parent_knots(
          mypknots, myknots, normalfac, parentele->id(), localsurfaceid);

      // if we have a zero sized element due to a interpolated
      // point --- exit here and tell the outside world about that
      if (zero_size)
      {
        return (zero_size);
      }
      // you are still here? So get the node weights as well
      Core::Nodes::Node** nodes = boundaryele->nodes();
      const int boundarynen = boundaryele->num_node();
      for (int inode = 0; inode < boundarynen; inode++)
      {
        Core::FE::Nurbs::ControlPoint* cp =
            dynamic_cast<Core::FE::Nurbs::ControlPoint*>(nodes[inode]);
        weights(inode) = cp->w();
      }

      Core::Nodes::Node** pnodes = parentele->nodes();
      const int pnen = parentele->num_node();
      for (int inode = 0; inode < pnen; inode++)
      {
        Core::FE::Nurbs::ControlPoint* cp =
            dynamic_cast<Core::FE::Nurbs::ControlPoint*>(pnodes[inode]);
        pweights(inode) = cp->w();
      }

      // goodbye
      return zero_size;
    }  // GetKnotVectorAndWeightsForNurbsBoundary()

  }  // namespace Nurbs

}  // namespace Core::FE

FOUR_C_NAMESPACE_CLOSE

#endif

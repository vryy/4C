/*----------------------------------------------------------------------*/
/*! \file

\brief service methods for accessing knot vector and weights for a given
       element of a nurbs discretisations

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_NURBS_DISCRET_NURBS_UTILS_HPP
#define FOUR_C_NURBS_DISCRET_NURBS_UTILS_HPP

#include "4C_config.hpp"

#include "4C_nurbs_discret.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace NURBS
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
    bool GetMyNurbsKnotsAndWeights(const DRT::Discretization& dis,
        const CORE::Elements::Element* ele, std::vector<CORE::LINALG::SerialDenseVector>& myknots,
        WG& weights)
    {
      // try to cast dis to a nurbs discretisation
      const DRT::NURBS::NurbsDiscretization* nurbsdis =
          dynamic_cast<const DRT::NURBS::NurbsDiscretization*>(&(dis));
      if (nurbsdis == nullptr) FOUR_C_THROW("Received discretization which is not Nurbs!");

      // get local knot vector entries and check for zero sized elements
      const bool zero_size = (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots, ele->Id());

      // if we have a zero sized element due to a interpolated
      // point --- exit here and tell the outside world about that
      if (zero_size)
      {
        return (zero_size);
      }
      // you are still here? So get the node weights for the nurbs element as well
      const DRT::Node* const* nodes = ele->Nodes();
      const int nen = ele->num_node();
      for (int inode = 0; inode < nen; inode++)
      {
        const DRT::NURBS::ControlPoint* cp =
            dynamic_cast<const DRT::NURBS::ControlPoint*>(nodes[inode]);
        weights(inode) = cp->W();
      }

      // goodbye
      return zero_size;
    }  // GetMyNurbsKnotsAndWeights()


    //! determine whether a given element is a NURBS element or not
    inline bool IsNurbs(CORE::FE::CellType distype)
    {
      switch (distype)
      {
        case CORE::FE::CellType::nurbs8:
        case CORE::FE::CellType::nurbs27:
        case CORE::FE::CellType::nurbs4:
        case CORE::FE::CellType::nurbs9:
        case CORE::FE::CellType::nurbs2:
        case CORE::FE::CellType::nurbs3:
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
    bool GetKnotVectorAndWeightsForNurbsBoundary(const CORE::Elements::Element* boundaryele,
        const int localsurfaceid, const int parenteleid, const DRT::Discretization& discretization,
        std::vector<CORE::LINALG::SerialDenseVector>& mypknots,
        std::vector<CORE::LINALG::SerialDenseVector>& myknots, WG& weights, double& normalfac)
    {
      // get knotvector(s)
      const DRT::NURBS::NurbsDiscretization* nurbsdis =
          dynamic_cast<const DRT::NURBS::NurbsDiscretization*>(&(discretization));

      Teuchos::RCP<const DRT::NURBS::Knotvector> knots = (*nurbsdis).GetKnotVector();

      bool zero_size = knots->get_boundary_ele_and_parent_knots(
          mypknots, myknots, normalfac, parenteleid, localsurfaceid);

      // if we have a zero sized element due to a interpolated
      // point --- exit here and tell the outside world about that
      if (zero_size)
      {
        return (zero_size);
      }
      // you are still here? So get the node weights as well
      const DRT::Node* const* nodes = boundaryele->Nodes();
      const int boundarynen = boundaryele->num_node();
      for (int inode = 0; inode < boundarynen; inode++)
      {
        const DRT::NURBS::ControlPoint* cp =
            dynamic_cast<const DRT::NURBS::ControlPoint*>(nodes[inode]);
        weights(inode) = cp->W();
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
    bool GetKnotVectorAndWeightsForNurbsBoundaryAndParent(CORE::Elements::Element* parentele,
        CORE::Elements::Element* boundaryele, const int localsurfaceid,
        const DRT::Discretization& discretization,
        std::vector<CORE::LINALG::SerialDenseVector>& mypknots,
        std::vector<CORE::LINALG::SerialDenseVector>& myknots, WG& pweights, WG& weights,
        double& normalfac)
    {
      // get knotvector(s)
      const DRT::NURBS::NurbsDiscretization* nurbsdis =
          dynamic_cast<const DRT::NURBS::NurbsDiscretization*>(&(discretization));

      Teuchos::RCP<const DRT::NURBS::Knotvector> knots = (*nurbsdis).GetKnotVector();

      bool zero_size = knots->get_boundary_ele_and_parent_knots(
          mypknots, myknots, normalfac, parentele->Id(), localsurfaceid);

      // if we have a zero sized element due to a interpolated
      // point --- exit here and tell the outside world about that
      if (zero_size)
      {
        return (zero_size);
      }
      // you are still here? So get the node weights as well
      DRT::Node** nodes = boundaryele->Nodes();
      const int boundarynen = boundaryele->num_node();
      for (int inode = 0; inode < boundarynen; inode++)
      {
        DRT::NURBS::ControlPoint* cp = dynamic_cast<DRT::NURBS::ControlPoint*>(nodes[inode]);
        weights(inode) = cp->W();
      }

      DRT::Node** pnodes = parentele->Nodes();
      const int pnen = parentele->num_node();
      for (int inode = 0; inode < pnen; inode++)
      {
        DRT::NURBS::ControlPoint* cp = dynamic_cast<DRT::NURBS::ControlPoint*>(pnodes[inode]);
        pweights(inode) = cp->W();
      }

      // goodbye
      return zero_size;
    }  // GetKnotVectorAndWeightsForNurbsBoundary()

  }  // namespace NURBS

}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif

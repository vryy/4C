/*-----------------------------------------------------------*/
/*! \file

\brief A class handling a (periodic) bounding box as simulation volume


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_BEAMINTERACTION_PERIODIC_BOUNDINGBOX_H
#define BACI_BEAMINTERACTION_PERIODIC_BOUNDINGBOX_H

#include "baci_config.H"

#include "baci_linalg_fixedsizematrix.H"

#include <Teuchos_RCP.hpp>

class Epetra_Vector;

BACI_NAMESPACE_OPEN

// forward declaration
namespace IO
{
  class DiscretizationVisualizationWriterMesh;
}
namespace DRT
{
  class Discretization;
}

namespace CORE::GEO
{
  namespace MESHFREE
  {
    /*!
     \brief Construction of periodic boundingbox over entire considered simulation volume
     */
    class BoundingBox
    {
     public:
      BoundingBox();

      virtual ~BoundingBox() = default;

      /// initialize bounding box object
      void Init();

      /// initialize bounding box object
      void Init(CORE::LINALG::Matrix<3, 2> const& box, std::vector<bool> const& pbconoff);

      /// setup bounding box object, setup call is needed in case of box dirichlet
      void Setup();

      /// get edge length
      double EdgeLength(int dim) const { return edgelength_[dim]; }

      /// get box
      CORE::LINALG::Matrix<3, 2> const& Box() const { return box_; }

      /// get flag indicating if periodic boundary conditions are active
      bool HavePBC() const { return haveperiodicbc_; }

      /// get const bounding box discretization
      DRT::Discretization const& BoundingBoxDiscret() const { return *boxdiscret_; }

      /// get corner points
      double operator()(int i, int j) const { return box_(i, j); }

      /// initialize bounding box discretization
      void SetupBoundingBoxDiscretization();

      /*!
      \brief shift node (if outside) back in box if periodic boundary conditions
      */
      bool Shift3D(CORE::LINALG::Matrix<3, 1>& d,
          CORE::LINALG::Matrix<3, 1> const X = CORE::LINALG::Matrix<3, 1>(true)) const;

      /*!
      \brief get xi of intersection between two points
      */
      void GetXiOfIntersection3D(CORE::LINALG::Matrix<3, 1> const& x1,
          CORE::LINALG::Matrix<3, 1> const& x2, CORE::LINALG::Matrix<3, 1>& xi) const;
      void GetXiOfIntersection3D(CORE::LINALG::Matrix<3, 1> const& x1,
          CORE::LINALG::Matrix<3, 1> const& x2, CORE::LINALG::Matrix<3, 1>& xi,
          CORE::LINALG::Matrix<3, 2> const& box) const;

      /*! Check the distance to a reference point position (e.g. node of the
       * same element). If the distance is larger than half of the period
       * length, the point position has been shifted before.
       *
       * Warning: This assumes that the distance between point and reference
       *          point is not larger than half of the period length unless we
       *          shift it. For beam elements, this restricts the element length
       *          to be smaller than this value throughout the entire simulation.
       *          So far, we only check this once in the beginning.
       *
       * Note: this should be equivalent to the previously applied criterion
       *       that the distance between given point and reference point
       *       decreases by either adding or subtracting the period length. */
      void UnShift3D(CORE::LINALG::Matrix<3, 1>& d, CORE::LINALG::Matrix<3, 1> const& ref,
          CORE::LINALG::Matrix<3, 1> const X = CORE::LINALG::Matrix<3, 1>(true)) const;

      bool CheckIfShiftBetweenPoints(CORE::LINALG::Matrix<3, 1>& d,
          CORE::LINALG::Matrix<3, 1> const& ref, std::vector<bool>& shift_in_dim,
          CORE::LINALG::Matrix<3, 1> const X = CORE::LINALG::Matrix<3, 1>(true)) const;

      /*!
      \brief get random position inside box
      */
      void RandomPosWithin(CORE::LINALG::Matrix<3, 1>& pos) const;

      /*!
       \brief If necessary make the boundingbox larger to include this point as one of the corners
       of the box
       */
      void AddPoint(double const* x);

      /*!
       \brief Check whether "b" is within this boundingbox
       */
      bool Within(const BoundingBox& b) const;

      /*!
       \brief Check if the point is within this boundingbox
       */
      bool Within(const double* x, std::vector<bool>& within_in_dir) const;
      bool Within(CORE::LINALG::Matrix<3, 1> const& x, std::vector<bool>& within_in_dir) const;
      bool Within(CORE::LINALG::Matrix<3, 2> const& box, CORE::LINALG::Matrix<3, 1> const& x,
          std::vector<bool>& within_in_dir) const;

      /*!
       \brief Check these points are within this boundingbox
       */
      bool Within(const CORE::LINALG::SerialDenseMatrix& xyz) const;

      /*!
       \brief Print the corner points of boundingbox on the screen
       */
      void Print();

      /*!
       \brief get min of box in certain dim
       */
      double box_min(int dim) const { return box_min(box_, dim); }
      double box_min(CORE::LINALG::Matrix<3, 2> const& box, int dim) const { return box(dim, 0); }

      /*!
       \brief get max of box in certain dim
       */
      double box_max(int dim) const { return box_max(box_, dim); }
      double box_max(CORE::LINALG::Matrix<3, 2> const& box, int dim) const { return box(dim, 1); }

      /*!
       \brief Get the outmost point of the boundingbox
       */
      void UndeformedBoxCornerPointPosition(int i, std::vector<double>& x) const;
      CORE::LINALG::Matrix<3, 1> UndeformedBoxCornerPointPosition(int i) const;
      /*!
       \brief get reference position of corner point i
       */
      CORE::LINALG::Matrix<3, 1> ReferencePosOfCornerPoint(int i) const;

      /*!
       \brief get current position of corner point i
       */
      CORE::LINALG::Matrix<3, 1> CurrentPositionOfCornerPoint(int i) const;

      /*!
       \brief print box
      */
      void Print(std::ostream& out) const;

      /*!
       \brief Write output
      */
      void RuntimeOutputStepState(double timen, int stepn) const;

      /*!
       \brief Apply dirichlet condition according to input file
      */
      void ApplyDirichlet(double timen);

      /*!
       \brief init runtime output object for bounding box discretization
      */
      void InitRuntimeOutput();

      //! @name public function dealing with mapping of positions in case of a deforming bounding
      //! box
      //! @{

      //! transform from undeformed to global
      void TransformFromUndeformedBoundingBoxSystemToGlobal(
          CORE::LINALG::Matrix<3, 1> const& xi, CORE::LINALG::Matrix<3, 1>& x) const;

      void TransformFromUndeformedBoundingBoxSystemToGlobal(double const* xi, double* x) const;

      //! transform from global to undeformed
      bool TransformFromGlobalToUndeformedBoundingBoxSystem(
          CORE::LINALG::Matrix<3, 1> const& x,  ///< input  -> global position
          CORE::LINALG::Matrix<3, 1>& xi  ///< output -> position in undeformed bounding box system
      ) const;
      bool TransformFromGlobalToUndeformedBoundingBoxSystem(
          double const* x,  ///< input  -> global position
          double* xi        ///< output -> position in undeformed bounding box system
      ) const;

      //! @}

     protected:
      //! returns init state
      inline bool IsInit() const { return isinit_; };

      //! returns setup state
      inline bool IsSetup() const { return issetup_; };

      //! Check the init state
      inline void ThrowIfNotInit() const
      {
        if (not IsInit()) dserror("Call Init() first!");
      }

      //! Check the init and setup state
      inline void ThrowIfNotInitOrSetup() const
      {
        if (not IsInit() or not IsSetup()) dserror("Call Init() and Setup() first!");
      }

     private:
      /*!
      \brief shift node (if outside) back in box if periodic boundary conditions
      */
      bool Shift1D(int dim, double& d, double const& X = 0.0) const;

      /*!
      \brief shift node out of box if it was shifted in previously
      */
      bool UnShift1D(int dim, double& d, double const& ref, double const& X = 0.0) const;

      bool InBetween(double smin, double smax, double omin, double omax) const;

      //! @name private function dealing with mapping of positions in case of a deforming bounding
      //! box
      //! @{

      //! evaluate lagrange polynomial that maps from undeformed to global at xi
      void LagrangePolynomialToMapFromUndeformedBoundingBoxSystemToGlobal(
          CORE::LINALG::Matrix<8, 1>& funct,  ///< to be filled with shape function values
          double r, double s, double t) const;

      //! evaluate first derivative of lagrange polynomial that maps from undeformed to global at xi
      void LagrangePolynomialToMapFromUndeformedBoundingBoxSystemToGlobalDeriv1(
          CORE::LINALG::Matrix<3, 8>&
              deriv1,  ///< to be filled with shape function derivative values
          double r, double s, double t) const;

      //! @}

     protected:
      //! @name member variables

      //! indicates if the Init() function has been called
      bool isinit_;

      //! indicates if the Setup() function has been called
      bool issetup_;

     private:
      /// discretization with one volume element representing the box ( used e.g. for output)
      Teuchos::RCP<DRT::Discretization> boxdiscret_;
      /// box displacement vector
      Teuchos::RCP<Epetra_Vector> disn_row_;
      Teuchos::RCP<Epetra_Vector> disn_col_;

      bool empty_;
      /// set global pbc flag
      bool haveperiodicbc_;
      /// set global dbc flag
      bool havedirichletbc_;
      /// box corners
      CORE::LINALG::Matrix<3, 2> box_;
      /// flags for existence of periodic boundary conditions in x, y, z direction
      bool pbconoff_[3];
      ///< box edge lengths in x, y, z direction
      double edgelength_[3];

      //! bounding box discretization runtime visualization writer
      Teuchos::RCP<IO::DiscretizationVisualizationWriterMesh> visualization_output_writer_ptr_;
    };

  }  // namespace MESHFREE
}  // namespace CORE::GEO

BACI_NAMESPACE_CLOSE

#endif

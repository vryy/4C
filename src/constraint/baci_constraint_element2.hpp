/*----------------------------------------------------------------------*/
/*! \file
\brief A 2D constraint element with no physics attached
\level 2


*----------------------------------------------------------------------*/

#ifndef BACI_CONSTRAINT_ELEMENT2_HPP
#define BACI_CONSTRAINT_ELEMENT2_HPP


#include "baci_config.hpp"

#include "baci_comm_parobjectfactory.hpp"
#include "baci_lib_element.hpp"
#include "baci_lib_elementtype.hpp"
#include "baci_lib_node.hpp"
#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_linalg_serialdensematrix.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    // forward declarations
    // class ConstraintElementLine;

    class ConstraintElement2Type : public DRT::ElementType
    {
     public:
      std::string Name() const override { return "ConstraintElement2Type"; }

      static ConstraintElement2Type& Instance();

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      CORE::LINALG::SerialDenseMatrix ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp) override;

     private:
      static ConstraintElement2Type instance_;
    };

    /*!
     */
    class ConstraintElement2 : public DRT::Element
    {
     public:
      //! @name Friends

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      */
      ConstraintElement2(int id,  ///< A unique global id
          int owner               ///< element owner
      );

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      ConstraintElement2(const ConstraintElement2& old);

      /*!
      \brief Deep copy this instance of ConstraintElement2 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const override;

      /*!
      \brief Get shape type of element
      */
      CORE::FE::CellType Shape() const override
      {
        dserror("ConstraintElement2 has no shape!");
        return CORE::FE::CellType::dis_none;
      };

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override
      {
        return ConstraintElement2Type::Instance().UniqueParObjectId();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Pack(CORE::COMM::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Unpack(const std::vector<char>& data) override;


      //@}

      //! @name Access methods


      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual DRT::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      int NumDofPerNode(const DRT::Node& node) const override { return 2; }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual DRT::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      int NumDofPerElement() const override { return 0; }

      DRT::ElementType& ElementType() const override { return ConstraintElement2Type::Instance(); }

      //@}

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      \return 0 if successful, negative otherwise
      */
      int Evaluate(Teuchos::ParameterList& params,   ///< ParameterList for communication
          DRT::Discretization& discretization,       ///< discretization
          std::vector<int>& lm,                      ///< location vector
          CORE::LINALG::SerialDenseMatrix& elemat1,  ///< first matrix to be filled by element
          CORE::LINALG::SerialDenseMatrix& elemat2,  ///< second matrix to be filled by element
          CORE::LINALG::SerialDenseVector& elevec1,  ///< third matrix to be filled by element
          CORE::LINALG::SerialDenseVector& elevec2,  ///< first vector to be filled by element
          CORE::LINALG::SerialDenseVector& elevec3   ///< second vector to be filled by element
          ) override;


      /*!
      \brief Evaluate a Neumann boundary condition

      Since the element has no physics attached this method will give a dserror

      \return 0 if successful, negative otherwise
      */
      int EvaluateNeumann(Teuchos::ParameterList& params,  ///< ParameterList for communication
          DRT::Discretization& discretization,             ///< discretization
          DRT::Condition& condition,                       ///< Neumann condition to evaluate
          std::vector<int>& lm,                            ///< location vector
          CORE::LINALG::SerialDenseVector& elevec1,        ///< vector to be filled by element
          CORE::LINALG::SerialDenseMatrix* elemat1 = nullptr) override;


     private:
      //! action parameters recognized by constraint element
      enum ActionType
      {
        none,
        calc_MPC_dist_stiff,
        calc_MPC_angle_stiff
      };

      //! vector of surfaces of this element (length 1)
      std::vector<DRT::Element*> surface_;

      // don't want = operator
      ConstraintElement2& operator=(const ConstraintElement2& old);


      /*!
      \brief Create matrix with material configuration for 2 dimensions and 3 nodes
      */
      inline void MaterialConfiguration(
          CORE::LINALG::Matrix<3, 2>& x  ///< nodal coords in material frame
      ) const
      {
        const int numnode = 3;
        const int numdim = 2;
        for (int i = 0; i < numnode; ++i)
        {
          for (int j = 0; j < numdim; ++j)
          {
            x(i, j) = Nodes()[i]->X()[j];
          }
        }
        return;
      }

      /*!
      \brief Create matrix with spatial configuration for 2 dimensions and 3 nodes
      */
      inline void SpatialConfiguration(
          CORE::LINALG::Matrix<3, 2>& x,  ///< nodal coords in spatial frame
          const std::vector<double> disp  ///< displacements
      ) const
      {
        const int numnode = 3;
        const int numdim = 2;
        for (int i = 0; i < numnode; ++i)
        {
          for (int j = 0; j < numdim; ++j)
          {
            x(i, j) = Nodes()[i]->X()[j] + disp[i * numdim + j];
          }
        }
        return;
      }


      /// compute normal for 2D case using the first two nodes to specify a line
      void ComputeNormal(const CORE::LINALG::Matrix<3, 2>& xc,  ///< nodal coords in spatial frame
          CORE::LINALG::Matrix<2, 1>& elenormal                 ///< resulting element normal
      );

      /// Compute normal distance between line and third node
      double ComputeNormalDist(
          const CORE::LINALG::Matrix<3, 2>& xc,        ///< nodal coords in spatial frame
          const CORE::LINALG::Matrix<2, 1>& elenormal  ///< element normal
      );

      /// Compute first derivatives of normal distance with respect to the nodal displacements
      void ComputeFirstDerivDist(
          const CORE::LINALG::Matrix<3, 2>& xc,        ///< nodal coords in spatial frame
          CORE::LINALG::SerialDenseVector& elevector,  ///< vector to store results into
          const CORE::LINALG::Matrix<2, 1>& elenormal  ///< element normal
      );

      /// Compute second derivatives of normal distance with respect to the nodal displacements
      void ComputeSecondDerivDist(
          const CORE::LINALG::Matrix<3, 2>& xc,        ///< nodal coords in spatial frame
          CORE::LINALG::SerialDenseMatrix& elematrix,  ///< matrix to store results into
          const CORE::LINALG::Matrix<2, 1>& elenormal  ///< element normal
      );

      /// Compute angle at second node
      double ComputeAngle(const CORE::LINALG::Matrix<3, 2>& xc  ///< nodal coords in spatial frame
      );

      /// Compute first derivatives of angle at second node with respect to the nodal displacements
      void ComputeFirstDerivAngle(
          const CORE::LINALG::Matrix<3, 2>& xc,       ///< nodal coords in spatial frame
          CORE::LINALG::SerialDenseVector& elevector  ///< vector to store results into
      );

      /// Compute second derivatives of angle at second node with respect to the nodal displacements
      void ComputeSecondDerivAngle(
          const CORE::LINALG::Matrix<3, 2>& xc,       ///< nodal coords in spatial frame
          CORE::LINALG::SerialDenseMatrix& elematrix  ///< matrix to store results into
      );

    };  // class ConstraintElement2


    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================


  }  // namespace ELEMENTS
}  // namespace DRT


BACI_NAMESPACE_CLOSE

#endif

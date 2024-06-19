/*----------------------------------------------------------------------*/
/*! \file

\brief prestress functionality in solid elements

\level 2

*----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_PRESTRESS_HPP
#define FOUR_C_SO3_PRESTRESS_HPP

#include "4C_config.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_comm_parobjectfactory.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace ELEMENTS
  {
    class PreStressType : public Core::Communication::ParObjectType
    {
     public:
      std::string Name() const override { return "PreStressType"; }

      static PreStressType& Instance() { return instance_; };

     private:
      static PreStressType instance_;
    };

    /*!
    \brief A class for handling the prestressing in finite deformations

    */
    class PreStress : public Core::Communication::ParObject
    {
     public:
      /*!
      \brief Standard Constructor
      */
      PreStress(const int numnode, const int ngp, const bool istet4 = false);

      /*!
      \brief Copy Constructor
      */
      PreStress(const Discret::ELEMENTS::PreStress& old);

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override;
      /*!
      \brief Pack this class so it can be communicated

      \ref pack and \ref unpack are used to communicate this node

      */
      void pack(Core::Communication::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref pack and \ref unpack are used to communicate this node

      */
      void unpack(const std::vector<char>& data) override;

      /// get history of deformation gradient
      inline Core::LinAlg::SerialDenseMatrix& FHistory() const { return *fhist_; }

      /// get history of of reference configuration (inverse of Jacobian)
      inline Core::LinAlg::SerialDenseMatrix& JHistory() const { return *inv_jhist_; }

      /// put a matrix to storage
      inline void MatrixtoStorage(const int gp, const Core::LinAlg::Matrix<3, 3>& Mat,
          Core::LinAlg::SerialDenseMatrix& gpMat) const
      {
        for (int i = 0; i < gpMat.numCols(); ++i) gpMat(gp, i) = Mat.A()[i];
        return;
      }

      /// put a matrix to storage
      inline void MatrixtoStorage(const int gp, const Core::LinAlg::Matrix<4, 3>& Mat,
          Core::LinAlg::SerialDenseMatrix& gpMat) const
      {
        for (int i = 0; i < gpMat.numCols(); ++i) gpMat(gp, i) = Mat.A()[i];
        return;
      }

      /// get matrix from storage
      inline void StoragetoMatrix(const int gp, Core::LinAlg::Matrix<3, 3>& Mat,
          const Core::LinAlg::SerialDenseMatrix& gpMat) const
      {
        for (int i = 0; i < gpMat.numCols(); ++i) Mat.A()[i] = gpMat(gp, i);
        return;
      }

      /// get matrix from storage
      inline void StoragetoMatrix(const int gp, Core::LinAlg::Matrix<4, 3>& Mat,
          const Core::LinAlg::SerialDenseMatrix& gpMat) const
      {
        for (int i = 0; i < gpMat.numCols(); ++i) Mat.A()[i] = gpMat(gp, i);
        return;
      }

      /// get indication whether class is initialized (important for restarts)
      bool& is_init() { return isinit_; }

     private:
      /// flagindicating whether material configuration has been initialized
      bool isinit_;

      /// no. nodal points of element
      int numnode_;

      /// history of deformation gradient
      Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> fhist_;

      /// updated Lagrange inverse of Jacobian
      Teuchos::RCP<Core::LinAlg::SerialDenseMatrix> inv_jhist_;

      /// get number of gaussian points considered
      inline int num_gp() const { return fhist_->numRows(); }

      /// get no. of nodal points
      inline int num_node() const { return numnode_; }

    };  // class PreStress
  }     // namespace ELEMENTS
}  // namespace Discret


FOUR_C_NAMESPACE_CLOSE

#endif

/*----------------------------------------------------------------------*/
/*! \file

\brief prestress functionality in solid elements

\level 2

\maintainer Christoph Meier
*----------------------------------------------------------------------*/
#ifndef PRESTRESS_H
#define PRESTRESS_H

#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace DRT
{
  namespace ELEMENTS
  {
    class PreStressType : public DRT::ParObjectType
    {
     public:
      std::string Name() const { return "PreStressType"; }

      static PreStressType& Instance() { return instance_; };

     private:
      static PreStressType instance_;
    };

    /*!
    \brief A class for handling the prestressing in finite deformations

    \author gee (gee@lnm.mw.tum.de)
    */
    class PreStress : public DRT::ParObject
    {
     public:
      /*!
      \brief Standard Constructor
      */
      PreStress(const int numnode, const int ngp, const bool istet4 = false);

      /*!
      \brief Copy Constructor
      */
      PreStress(const DRT::ELEMENTS::PreStress& old);

      /*!
      \brief Destructor

      */
      virtual ~PreStress() {}

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      virtual int UniqueParObjectId() const;
      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this node

      */
      virtual void Pack(DRT::PackBuffer& data) const;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this node

      */
      virtual void Unpack(const std::vector<char>& data);

      /// get history of deformation gradient
      inline Epetra_SerialDenseMatrix& FHistory() const { return *Fhist_; }

      /// get history of of reference configuration (inverse of Jacobian)
      inline Epetra_SerialDenseMatrix& JHistory() const { return *invJhist_; }

      /// put a matrix to storage
      inline void MatrixtoStorage(
          const int gp, const LINALG::Matrix<3, 3>& Mat, Epetra_SerialDenseMatrix& gpMat) const
      {
        for (int i = 0; i < gpMat.N(); ++i) gpMat(gp, i) = Mat.A()[i];
        return;
      }

      /// put a matrix to storage
      inline void MatrixtoStorage(
          const int gp, const LINALG::Matrix<4, 3>& Mat, Epetra_SerialDenseMatrix& gpMat) const
      {
        for (int i = 0; i < gpMat.N(); ++i) gpMat(gp, i) = Mat.A()[i];
        return;
      }

      /// get matrix from storage
      inline void StoragetoMatrix(
          const int gp, LINALG::Matrix<3, 3>& Mat, const Epetra_SerialDenseMatrix& gpMat) const
      {
        for (int i = 0; i < gpMat.N(); ++i) Mat.A()[i] = gpMat(gp, i);
        return;
      }

      /// get matrix from storage
      inline void StoragetoMatrix(
          const int gp, LINALG::Matrix<4, 3>& Mat, const Epetra_SerialDenseMatrix& gpMat) const
      {
        for (int i = 0; i < gpMat.N(); ++i) Mat.A()[i] = gpMat(gp, i);
        return;
      }

      /// get indication whether class is initialized (important for restarts)
      bool& IsInit() { return isinit_; }

     private:
      /// flagindicating whether material configuration has been initialized
      bool isinit_;

      /// no. nodal points of element
      int numnode_;

      /// history of deformation gradient
      Teuchos::RCP<Epetra_SerialDenseMatrix> Fhist_;

      /// updated Lagrange inverse of Jacobian
      Teuchos::RCP<Epetra_SerialDenseMatrix> invJhist_;

      /// get number of gaussian points considered
      inline int NGP() const { return Fhist_->M(); }

      /// get no. of nodal points
      inline int NumNode() const { return numnode_; }

    };  // class PreStress
  }     // namespace ELEMENTS
}  // namespace DRT


#endif  // #ifndef PreStress

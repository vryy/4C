/*!----------------------------------------------------------------------
\file hyperpolyconvex.H
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifndef HYPERPOLYCONVEX_H_
#define HYPERPOLYCONVEX_H_
#ifdef CCADISCRET

#include "material.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{

  /// Wrapper for hyperelastic polyconvex anisotropic material
  class HyperPolyconvex : public Material
  {
  public:

    /// construct empty material object
    HyperPolyconvex();

    /// construct the material object given material parameters
    explicit HyperPolyconvex(MATERIAL* matdata);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ParObject_HyperPolyconvex; }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    virtual void Pack(vector<char>& data) const;

    /*!
      \brief Unpack data from a char vector into this class

      The vector data contains all information to rebuild the
      exact copy of an instance of a class on a different processor.
      The first entry in data has to be an integer which is the unique
      parobject id defined at the top of this file and delivered by
      UniqueParObjectId().

      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    virtual void Unpack(const vector<char>& data);

    //@}

    /// material type
    virtual MATERIAL_TYP MaterialType() const { return m_hyper_polyconvex; }

    // wrapper for less intensive error output information
    void Evaluate(const Epetra_SerialDenseVector* glstrain,
                  Epetra_SerialDenseMatrix* cmat,
                  Epetra_SerialDenseVector* stress)
    {
      Epetra_SerialDenseMatrix defgrd(3,3);
      Evaluate(glstrain,&defgrd,-1,-1,-1.0,cmat,stress);
      return;
    }

    // wrapper for less intensive error output information
    void Evaluate(LINALG::FixedSizeSerialDenseMatrix<6,1>* glstrain,
                  LINALG::FixedSizeSerialDenseMatrix<6,6>* cmat,
                  LINALG::FixedSizeSerialDenseMatrix<6,1>* stress)
    {
      Epetra_SerialDenseMatrix defgrd(3,3);
      Epetra_SerialDenseVector glstrain_e(Copy,glstrain->A(),glstrain->Rows());
      Epetra_SerialDenseVector stress_e(View,stress->A(),stress->Rows());
      Epetra_SerialDenseMatrix cmat_e(View,cmat->A(),cmat->Rows(),cmat->Rows(),cmat->Columns());

      Evaluate(&glstrain_e,&defgrd,-1,-1,-1.0,&cmat_e,&stress_e);
      return;
    }

    // there are a lot of parameters just for error output,
    void Evaluate(const Epetra_SerialDenseVector* glstrain,
                  const Epetra_SerialDenseMatrix* defgrd,
                  const int gp,  const int ele_ID, const double time,
                  Epetra_SerialDenseMatrix* cmat,
                  Epetra_SerialDenseVector* stress);

    // there are a lot of parameters just for error output,
    void Evaluate(LINALG::FixedSizeSerialDenseMatrix<6,1>* glstrain,
                  LINALG::FixedSizeSerialDenseMatrix<3,3>* defgrd,
                  const int gp,  const int ele_ID, const double time,
                  LINALG::FixedSizeSerialDenseMatrix<6,6>* cmat,
                  LINALG::FixedSizeSerialDenseMatrix<6,1>* stress)
    {
      Epetra_SerialDenseVector glstrain_e(Copy,glstrain->A(),glstrain->Rows());
      Epetra_SerialDenseVector stress_e(View,stress->A(),stress->Rows());
      Epetra_SerialDenseMatrix defgrd_e(Copy,defgrd->A(),defgrd->Rows(),defgrd->Rows(),defgrd->Columns());
      Epetra_SerialDenseMatrix cmat_e(View,cmat->A(),cmat->Rows(),cmat->Rows(),cmat->Columns());

      Evaluate(&glstrain_e, &defgrd_e, gp, ele_ID, time, &cmat_e, &stress_e);
    }

    double Density();


  private:
    void hyper_ElastSymTensorMultiply(Epetra_SerialDenseMatrix& C,
                                     const double ScalarAB,
                                     const Epetra_SerialDenseMatrix& A,
                                     const Epetra_SerialDenseMatrix& B,
                                     const double ScalarThis);
    void hyper_ElastSymTensorMultiplyAddSym(Epetra_SerialDenseMatrix& C,
                                     const double ScalarAB,
                                     const Epetra_SerialDenseMatrix& A,
                                     const Epetra_SerialDenseMatrix& B,
                                     const double ScalarThis);
    void hyper_ElastSymTensor_o_Multiply(Epetra_SerialDenseMatrix& C,
                                     const double ScalarAB,
                                     const Epetra_SerialDenseMatrix& A,
                                     const Epetra_SerialDenseMatrix& B,
                                     const double ScalarThis);
//    Epetra_SerialDenseMatrix tensorproduct(const Epetra_SerialDenseMatrix A,
//                                           const Epetra_SerialDenseMatrix B,
//                                           const double scalarA,
//                                           const double scalarB);


    /// my material parameters
    MATERIAL* matdata_;

  };
}

#endif /*CCADISCRET*/
#endif /*HYPERPOLYCONVEX_H_*/

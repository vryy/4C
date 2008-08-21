/*!----------------------------------------------------------------------
\file inversedesign.H

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)
#ifndef INVERSEDESIGN_H
#define INVERSEDESIGN_H

#include "../drt_lib/drt_parobject.H"

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

namespace DRT
{
namespace ELEMENTS
{
  using namespace Teuchos;
  using namespace std;
  class So_hex8;

/*!
\brief A class for handling the prestressing in finite deformations

\author gee (gee@lnm.mw.tum.de)
*/
class InvDesign : public DRT::ParObject
{
public:

  /*!
  \brief Standard Constructor
  */
  InvDesign(const int numnod, const int ngp);

  /*!
  \brief Copy Constructor
  */
  InvDesign(const DRT::ELEMENTS::InvDesign& old);

  /*!
  \brief Destructor

  */
  virtual ~InvDesign() {}

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return ParObject_InvDesign;}

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this node

  */
  virtual void Pack(vector<char>& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this node

  */
  virtual void Unpack(const vector<char>& data);

  /// put a matrix to storage
  inline void MatrixtoStorage(const int gp, 
                              const Epetra_SerialDenseMatrix& Mat,
                              Epetra_SerialDenseMatrix& gpMat) const
  {
    for (int i=0; i<gpMat.N(); ++i) gpMat(gp,i) = Mat.A()[i];
    return;
  }
  
  /// get matrix from storage
  inline void StoragetoMatrix(const int gp, 
                              Epetra_SerialDenseMatrix& Mat,
                              const Epetra_SerialDenseMatrix& gpMat) const
  {
    for (int i=0; i<gpMat.N(); ++i) Mat.A()[i] = gpMat(gp,i);
    return;
  }
  
  void soh8_nlnstiffmass(
      DRT::ELEMENTS::So_hex8*   ele,            ///< this element
      vector<int>&              lm,             ///< location matrix
      vector<double>&           disp,           ///< current displacements
      vector<double>&           residual,       ///< current residuum
      Epetra_SerialDenseMatrix* stiffmatrix,    ///< element stiffness matrix
      Epetra_SerialDenseMatrix* massmatrix,     ///< element mass matrix
      Epetra_SerialDenseVector* force,          ///< element internal force vector
      Epetra_SerialDenseMatrix* elestress,      ///< stresses at GP
      Epetra_SerialDenseMatrix* elestrain,      ///< strains at GP
      ParameterList&            params,         ///< algorithmic parameters e.g. time
      const bool                cauchy=false,   ///< stress output option
      const bool                euler_almansi=false);  ///< strain output option

private:

  int                      numnod_;
  int                      ngp_;

  void BuildLambda(double Lambda4[][3][3][3], LINALG::SerialDenseMatrix& L, const LINALG::SerialDenseMatrix& F) const;
  void BuildLambdaT(LINALG::SerialDenseMatrix& L, const LINALG::SerialDenseMatrix& F) const;
  void BuildIF(LINALG::SerialDenseMatrix& IF, const LINALG::SerialDenseMatrix& F) const;
  void BuildTheta(double Theta4[][3][3][3],LINALG::SerialDenseMatrix& Theta, const LINALG::SerialDenseMatrix& F) const;
  void BuildYpsilon(double Y4[][3][3][3],
                    LINALG::SerialDenseMatrix& Y, 
                    const LINALG::SerialDenseMatrix& F,
                    const LINALG::SerialDenseMatrix& S) const;

  // these were for debugging
  void FDLambda(LINALG::SerialDenseMatrix& L, const LINALG::SerialDenseMatrix& f) const;
  void FDLambdaT(LINALG::SerialDenseMatrix& L, const LINALG::SerialDenseMatrix& f) const;
  void FDTheta(double Theta4[][3][3][3],LINALG::SerialDenseMatrix& Theta, const LINALG::SerialDenseMatrix& F) const;
  void FDYpsilon(double Y4[][3][3][3],LINALG::SerialDenseMatrix& Y, 
                 const LINALG::SerialDenseMatrix& F,
                 const LINALG::SerialDenseMatrix& S) const;
  void FDstiffmatrix(Epetra_SerialDenseMatrix& stiff,
                     const vector<double>& disp,
                     const int gp,
                     DRT::ELEMENTS::So_hex8* ele,
                     ParameterList& params) const;
  void FD_djdX(Epetra_SerialDenseMatrix& djdX,
                     const vector<double>& disp,
                     const int gp,
                     DRT::ELEMENTS::So_hex8* ele) const;
  void FD_dISdX(Epetra_SerialDenseMatrix& stiff,
                     const vector<double>& disp,
                     const int gp,
                     DRT::ELEMENTS::So_hex8* ele,
                     ParameterList& params) const;
  void FD_dISdf(Epetra_SerialDenseMatrix& dISdf,
                const Epetra_SerialDenseVector& stress,
                const LINALG::SerialDenseMatrix& f) const;
  void TensorMultiply(Epetra_SerialDenseMatrix& sum,double fac,
                      double Y4[][3][3][3], double L4[][3][3][3]);
 
}; // class InvDesign
} // namespace ELEMENTS
} // namespace DRT











#endif  // #ifndef INVERSEDESIGN_H
#endif  // #if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)
#endif  // #ifdef CCADISCRET

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
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_fixedsizematrix.H"

#include "so_hex8.H"
#include "so_weg6.H"
#include "so_tet4.H"


// forward declarations
struct _MATERIAL;

namespace DRT
{

namespace ELEMENTS
{

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
  InvDesign(const int numnod, const int ngp,  const bool istet4=false);

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

  // M could be a FixedSize<3,3> or a FixedSize<4,3>, for tet4.
  /// put a matrix to storage
  template <typename M>
  inline void MatrixtoStorage(const int gp,
                              const M& Mat,
                              Epetra_SerialDenseMatrix& gpMat) const
  {
    for (int i=0; i<gpMat.N(); ++i) gpMat(gp,i) = Mat.A()[i];
    return;
  }

  // M could be a FixedSize<3,3> or a FixedSize<4,3>, for tet4.
  /// get matrix from storage
  template <typename M>
  inline void StoragetoMatrix(const int gp,
                              M& Mat,
                              const Epetra_SerialDenseMatrix& gpMat) const
  {
    for (int i=0; i<gpMat.N(); ++i) Mat.A()[i] = gpMat(gp,i);
    return;
  }

  /// get history of deformation gradient
  inline Epetra_SerialDenseMatrix& FHistory() const { return *Fhist_; }

  /// get history of reference configuration (inverse of Jacobian)
  inline Epetra_SerialDenseMatrix& JHistory() const { return *invJhist_; }

  /// get history of determinant of J in reference configuration
  inline vector<double>& DetJHistory() { return detJ_; }

  /// access the init flag
  inline bool& IsInit() { return isinit_; }

  /// id problem for hex8 element
  void soh8_nlnstiffmass(
      DRT::ELEMENTS::So_hex8*   ele,            ///< this element
      vector<int>&              lm,             ///< location matrix
      vector<double>&           disp,           ///< current displacements
      vector<double>&           residual,       ///< current residuum
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOH8,NUMDOF_SOH8>* stiffmatrix,    ///< element stiffness matrix
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOH8,NUMDOF_SOH8>* massmatrix,     ///< element mass matrix
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOH8,          1>* force,          ///< element internal force vector
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestress,      ///< stresses at GP
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOH8,NUMSTR_SOH8>* elestrain,      ///< strains at GP
      ParameterList&            params,         ///< algorithmic parameters e.g. time
      const bool                cauchy=false,   ///< stress output option
      const bool                euler_almansi=false);  ///< strain output option

  /// store material config for hex8 element
  void soh8_StoreMaterialConfiguration(DRT::ELEMENTS::So_hex8* ele,
                                       const vector<double>& disp);


  /// id problem for wedge6 element
void sow6_nlnstiffmass(
      DRT::ELEMENTS::So_weg6*   ele,            ///< this element
      vector<int>&              lm,             ///< location matrix
      vector<double>&           disp,           ///< current displacements
      vector<double>&           residual,       ///< current residuum
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_WEG6,NUMDOF_WEG6>* stiffmatrix,    ///< element stiffness matrix
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_WEG6,NUMDOF_WEG6>* massmatrix,     ///< element mass matrix
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_WEG6,          1>* force,          ///< element internal force vector
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_WEG6,NUMSTR_WEG6>* elestress,      ///< stresses at GP
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_WEG6,NUMSTR_WEG6>* elestrain,      ///< strains at GP
      ParameterList&            params,         ///< algorithmic parameters e.g. time
      const bool                cauchy=false,         ///< stress output option
      const bool                euler_almansi=false); ///< strain output option

  /// store material config for wedge6 element
  void sow6_StoreMaterialConfiguration(DRT::ELEMENTS::So_weg6* ele,
                                       const vector<double>& disp);


  /// id problem for tet4 element
void so_tet4_nlnstiffmass(
      DRT::ELEMENTS::So_tet4*   ele,            ///< this element
      vector<int>&              lm,             ///< location matrix
      vector<double>&           disp,           ///< current displacements
      vector<double>&           residual,       ///< current residuum
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOTET4,NUMDOF_SOTET4>* stiffmatrix,    ///< element stiffness matrix
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOTET4,NUMDOF_SOTET4>* massmatrix,     ///< element mass matrix
      LINALG::FixedSizeSerialDenseMatrix<NUMDOF_SOTET4,            1>* force,          ///< element internal force vector
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOTET4,NUMSTR_SOTET4>* elestress,      ///< stresses at GP
      LINALG::FixedSizeSerialDenseMatrix<NUMGPT_SOTET4,NUMSTR_SOTET4>* elestrain,      ///< strains at GP
      struct _MATERIAL*         material,       ///< element material data
      const bool                cauchy=false,   ///< stress output options
      const bool                ea=false);

  /// store material config for wedge6 element
  void sot4_StoreMaterialConfiguration(DRT::ELEMENTS::So_tet4* ele,
                                       const vector<double>& disp);

private:

  int                                     numnod_;
  int                                     ngp_;
  bool                                    isinit_;
  Teuchos::RCP<Epetra_SerialDenseMatrix>  Fhist_;
  Teuchos::RCP<Epetra_SerialDenseMatrix>  invJhist_;
  vector<double>                          detJ_;

  /// get number of gaussian points considered
  inline int NGP() { return Fhist_->M(); }

  /// get no. of nodal points
  inline int NumNode() { return numnod_; }

  void BuildLambda(LINALG::FixedSizeSerialDenseMatrix<9,9>& L, const LINALG::FixedSizeSerialDenseMatrix<3,3>& F) const;
  void BuildIF(LINALG::FixedSizeSerialDenseMatrix<6,6>& IF, const LINALG::FixedSizeSerialDenseMatrix<3,3>& F) const;
  void BuildTheta(LINALG::FixedSizeSerialDenseMatrix<6,9>& Theta, const LINALG::FixedSizeSerialDenseMatrix<3,3>& F) const;
  void BuildYpsilon(LINALG::FixedSizeSerialDenseMatrix<6,9>& Y,
                    const LINALG::FixedSizeSerialDenseMatrix<3,3>& F,
                    const LINALG::FixedSizeSerialDenseMatrix<3,3>& S) const;


  // these were for debugging
#if 0
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
#endif


}; // class InvDesign
} // namespace ELEMENTS
} // namespace DRT











#endif  // #ifndef INVERSEDESIGN_H
#endif  // #if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)
#endif  // #ifdef CCADISCRET

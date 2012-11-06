/*!----------------------------------------------------------------------
\file inversedesign.H

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifndef INVERSEDESIGN_H
#define INVERSEDESIGN_H

#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_mat/material.H"


namespace LINALG
{
  class SerialDenseVector;
}
namespace DRT
{

namespace ELEMENTS
{
class So_tet4;
class So_weg6;
class So_hex8;

class InvDesignType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "InvDesignType"; }

  static InvDesignType & Instance() { return instance_; };

private:

  static InvDesignType instance_;
};

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
  virtual int UniqueParObjectId() const { return InvDesignType::Instance().UniqueParObjectId(); }

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this node

  */
  virtual void Pack(DRT::PackBuffer& data) const;

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
      LINALG::Matrix<24,24>* stiffmatrix,    ///< element stiffness matrix
      LINALG::Matrix<24,24>* massmatrix,     ///< element mass matrix
      LINALG::Matrix<24,1>* force,          ///< element internal force vector
      LINALG::Matrix<8,6>* elestress,      ///< stresses at GP
      LINALG::Matrix<8,6>* elestrain,      ///< strains at GP
      Teuchos::ParameterList&            params,         ///< algorithmic parameters e.g. time
      const INPAR::STR::StressType iostress,    ///< stress output option
      const INPAR::STR::StrainType iostrain);   ///< strain output option

  /// store material config for hex8 element
  void soh8_StoreMaterialConfiguration(DRT::ELEMENTS::So_hex8* ele,
                                       const vector<double>& disp);


  /// id problem for wedge6 element
void sow6_nlnstiffmass(
      DRT::ELEMENTS::So_weg6*   ele,            ///< this element
      vector<int>&              lm,             ///< location matrix
      vector<double>&           disp,           ///< current displacements
      vector<double>&           residual,       ///< current residuum
      LINALG::Matrix<18,18>* stiffmatrix,    ///< element stiffness matrix
      LINALG::Matrix<18,18>* massmatrix,     ///< element mass matrix
      LINALG::Matrix<18,1>* force,          ///< element internal force vector
      LINALG::Matrix<6,6>* elestress,      ///< stresses at GP
      LINALG::Matrix<6,6>* elestrain,      ///< strains at GP
      Teuchos::ParameterList&            params,         ///< algorithmic parameters e.g. time
      const INPAR::STR::StressType iostress,    ///< stress output option
      const INPAR::STR::StrainType iostrain);   ///< strain output option

  /// store material config for wedge6 element
  void sow6_StoreMaterialConfiguration(DRT::ELEMENTS::So_weg6* ele,
                                       const vector<double>& disp);


  /// id problem for tet4 element
void so_tet4_nlnstiffmass(
      Teuchos::ParameterList&           params,
      DRT::ELEMENTS::So_tet4*   ele,            ///< this element
      vector<int>&              lm,             ///< location matrix
      vector<double>&           disp,           ///< current displacements
      vector<double>&           residual,       ///< current residuum
      LINALG::Matrix<12,12>* stiffmatrix,    ///< element stiffness matrix
      LINALG::Matrix<12,12>* massmatrix,     ///< element mass matrix
      LINALG::Matrix<12,1>* force,          ///< element internal force vector
      LINALG::Matrix<1,6>* elestress,      ///< stresses at GP
      LINALG::Matrix<1,6>* elestrain,      ///< strains at GP
      const INPAR::STR::StressType iostress,    ///< stress output option
      const INPAR::STR::StrainType iostrain);   ///< strain output option

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

  void BuildLambda(LINALG::Matrix<9,9>& L, const LINALG::Matrix<3,3>& F) const;
  void BuildIF(LINALG::Matrix<6,6>& IF, const LINALG::Matrix<3,3>& F) const;
  void BuildTheta(LINALG::Matrix<6,9>& Theta, const LINALG::Matrix<3,3>& F) const;
  void BuildYpsilon(LINALG::Matrix<6,9>& Y,
                    const LINALG::Matrix<3,3>& F,
                    const LINALG::Matrix<3,3>& S) const;


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

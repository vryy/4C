/*!----------------------------------------------------------------------
\file contact_interpolator.H

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Prof. Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*
 | definitions                                             farah 09/14 |
 *---------------------------------------------------------------------*/
#ifndef CONTACT_INTERPOLATOR_H_
#define CONTACT_INTERPOLATOR_H_

/*---------------------------------------------------------------------*
 | headers                                                 farah 09/14 |
 *---------------------------------------------------------------------*/
#include "../headers/pairedvector.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_inpar/inpar_wear.H"

/*---------------------------------------------------------------------*
 | forward declarations                                    farah 09/14 |
 *---------------------------------------------------------------------*/
namespace MORTAR
{
class MortarElement;
}

namespace LINALG
{
class SerialDenseVector;
class SerialDenseMatrix;
}

/*!
 \brief CONTACT: namespace of the BACI contact problems

 */
namespace CONTACT
{
class CoNode;

class CoInterpolator
{
public:

  /*!
   \brief Constructor

   */
  CoInterpolator(Teuchos::ParameterList& params);

  /*!
   \brief Destructor

   */
  virtual ~CoInterpolator()
  {
  };

  /*!
   \brief Interpolate for 2D problems

   */
  void Interpolate2D(MORTAR::MortarElement& sele,
                     std::vector<MORTAR::MortarElement*> meles);

  /*!
   \brief Interpolate for 3D problems

   */
  void Interpolate3D(MORTAR::MortarElement& sele,
                     std::vector<MORTAR::MortarElement*> meles);

  /*!
   \brief lin 2D projection

   */
  void DerivXiGP2D(MORTAR::MortarElement& sele,
                   MORTAR::MortarElement& mele,
                   double& sxigp, double& mxigp,
                   const GEN::pairedvector<int,double>& derivsxi,
                   GEN::pairedvector<int,double>& derivmxi,
                   int& linsize);

  /*!
   \brief lin 3D projection

   */
  void DerivXiGP3D(MORTAR::MortarElement& sele,
                   MORTAR::MortarElement& mele,
                   double* sxigp, double* mxigp,
                   const std::vector<GEN::pairedvector<int,double> >& derivsxi,
                   std::vector<GEN::pairedvector<int,double> >& derivmxi,
                   double& alpha);

  /*!
   \brief node-wise D/M calculation

   */
  void nwDM2D(CoNode& mynode,
            MORTAR::MortarElement& mele,
            LINALG::SerialDenseVector& mval,
            LINALG::SerialDenseMatrix& mderiv,
            GEN::pairedvector<int,double>& dmxi);

  /*!
   \brief node-wise D/M calculation for 3D problems

   */
  void nwDM3D(CoNode& mynode,
             MORTAR::MortarElement& mele,
             LINALG::SerialDenseVector& mval,
             LINALG::SerialDenseMatrix& mderiv,
             std::vector<GEN::pairedvector<int,double> >& dmxi);

  /*!
   \brief node-wise gap calculation

   */
  void nwGap2D(CoNode& mynode,
               MORTAR::MortarElement& mele,
               LINALG::SerialDenseVector& mval,
               LINALG::SerialDenseMatrix& mderiv,
               GEN::pairedvector<int,double>& dmxi,
               double* gpn);

  /*!
   \brief node-wise gap calculation for 3D problems

   */
  void nwGap3D(CoNode& mynode,
               MORTAR::MortarElement& mele,
               LINALG::SerialDenseVector& mval,
               LINALG::SerialDenseMatrix& mderiv,
               std::vector<GEN::pairedvector<int,double> >& dmxi,
               double* gpn);

  /*!
   \brief node-wise slip calculation

   */
  void nwSlip2D(CoNode& mynode,
                MORTAR::MortarElement& mele,
                LINALG::SerialDenseVector& mval,
                LINALG::SerialDenseMatrix& mderiv,
                LINALG::SerialDenseMatrix& scoord,
                LINALG::SerialDenseMatrix& mcoord,
                Teuchos::RCP<LINALG::SerialDenseMatrix> scoordold,
                Teuchos::RCP<LINALG::SerialDenseMatrix> mcoordold,
                int& snodes,
                int& linsize,
                GEN::pairedvector<int,double>& dmxi);

  /*!
   \brief node-wise wear calculation (internal state var.)

   */
  void nwWear2D(CoNode& mynode,
                MORTAR::MortarElement& mele,
                LINALG::SerialDenseVector& mval,
                LINALG::SerialDenseMatrix& mderiv,
                LINALG::SerialDenseMatrix& scoord,
                LINALG::SerialDenseMatrix& mcoord,
                Teuchos::RCP<LINALG::SerialDenseMatrix> scoordold,
                Teuchos::RCP<LINALG::SerialDenseMatrix> mcoordold,
                Teuchos::RCP<LINALG::SerialDenseMatrix> lagmult,
                int& snodes,
                int& linsize,
                double& jumpval,
                double& area,
                double* gpn,
                GEN::pairedvector<int,double>& dmxi,
                GEN::pairedvector<int,double>& dslipmatrix,
                GEN::pairedvector<int,double>& dwear);

  /*!
   \brief node-wise wear calculation (primary variable)

   */
  void nwTE2D(CoNode& mynode,
              double& area,
              double& jumpval,
              GEN::pairedvector<int,double>& dslipmatrix);

protected:

  Teuchos::ParameterList&                          iparams_;     // containing contact input parameters
  bool                                             pwslip_;      // point-wise evaluated slip increment

  // wear inputs from parameter list
  INPAR::CONTACT::WearLaw                          wearlaw_;     // type of wear law
  bool                                             wearimpl_;    // flag for implicit wear algorithm
  INPAR::CONTACT::WearSide                         wearside_;    // definition of wear surface
  INPAR::CONTACT::WearType                         weartype_;    // definition of contact wear algorithm
  INPAR::CONTACT::WearShape                        wearshapefcn_;// type of wear shape function
  double                                           wearcoeff_;   // wear coefficient
  double                                           wearcoeffm_;  // wear coefficient master
  bool                                             sswear_;      // flag for steady state wear
  double                                           ssslip_;      // fixed slip for steady state wear

};


}// namespace CONTACT

#endif /* CONTACT_INTERPOLATOR_H_ */

/*!----------------------------------------------------------------------
\file contact_integrator.H

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

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/
#ifndef CONTACT_INTEGRATOR_H
#define CONTACT_INTEGRATOR_H

#include "../drt_mortar/mortar_integrator.H"
#include "../drt_inpar/inpar_wear.H"
#include "../headers/pairedvector.H"

// forward declarations
namespace LINALG
{
class SerialDenseVector;
}
/*!
 \brief CONTACT: namespace of the BACI contact and meshtying algorithms

 */
namespace CONTACT
{
/*!
 \brief A class to perform Gaussian integration and assembly of Mortar
 matrices on the overlap of two MortarElements (1 Slave, 1 Master)
 in 1D (which is equivalent to a 2D coupling problem) and in 2D
 (which is equivalent to a 3D coupling problem).

 This is a derived class from MORTAR::MortarIntegrator which does
 the contact-specific stuff for 3d mortar coupling.

 \author popp (popp@lnm.mw.tum.de)
 */

class CoIntegrator
{
public:

  /*!
   \brief Constructor  with shape function specification

   Constructs an instance of this class using a specific type of shape functions.<br>
   Note that this is \b not a collective call as overlaps are
   integrated in parallel by individual processes.<br>
   Note also that this constructor relies heavily on the
   DRT::UTILS::IntegrationPoints structs to get Gauss points
   and corresponding weights.

   */
  CoIntegrator(Teuchos::ParameterList& params,
      DRT::Element::DiscretizationType eletype);

  CoIntegrator(Teuchos::ParameterList& params,
      DRT::Element::DiscretizationType eletype, const Epetra_Comm& comm);

  /*!
   \brief Destructor

   */
  virtual ~CoIntegrator()
  {
  }

  /*!
   \brief Initialize Gauss rule (points, weights) for this MortarIntegrator

   */
  void InitializeGP(DRT::Element::DiscretizationType eletype);

  //! @name 2D and 3D integration methods

  /*!
   \brief check for boundary segmentation in 2D

   */
  bool BoundarySegmCheck2D(MORTAR::MortarElement& sele,
      std::vector<MORTAR::MortarElement*> meles);

  /*!
   \brief check for boundary segmentation in 2D

   */
  bool BoundarySegmCheck3D(MORTAR::MortarElement& sele,
      std::vector<MORTAR::MortarElement*> meles);

  /*!
   \brief evaluate D2-matrix entries at GP

   */
  void inline GP_D2(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& lm2val, LINALG::SerialDenseVector& m2val,
      double& jac, double& wgt, const Epetra_Comm& comm);

  /*!
   \brief evaluate D/M-matrix entries at GP

   */
  void inline GP_DM(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& lmval, LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& mval, double& jac, double& wgt, int& nrow,
      int& ncol, int& ndof, bool& bound);

  /*!
   \brief evaluate D/M-matrix entries at GP (3D quadratic)

   */
  void inline GP_3D_DM_Quad(MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele, MORTAR::IntElement& sintele,
      LINALG::SerialDenseVector& lmval, LINALG::SerialDenseVector& lmintval,
      LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& mval,
      const double& jac, double& wgt, const int& nrow, const int& nintrow,
      const int& ncol, const int& ndof, bool& bound);

  /*!
   \brief lin D/M-matrix entries at GP for bound case

   */
  void inline GP_DM_Lin_bound(
      int& iter,
      bool& duallin,
      MORTAR::MortarElement& sele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& lmderiv,
      double& dsxideta,
      double& dxdsxi,
      double& dxdsxidsxi,
      double& wgt,
      double& sxi,
      const GEN::pairedvector<int, double>& dsxigp
      ,
      const GEN::pairedvector<int, double>& derivjac
      ,
      const std::vector<GEN::pairedvector<int, double> >& ximaps
      ,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap);

  /*!
   \brief lin D/M-matrix entries at GP for bound case

   */
  void inline GP_2D_DM_Lin(
      int& iter,
      bool& bound,
      bool& linlm,
      MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& mderiv,
      LINALG::SerialDenseMatrix& lmderiv,
      double& dsxideta,
      double& dxdsxi,
      double& dxdsxidsxi,
      double& wgt,
      const GEN::pairedvector<int, double>& dsxigp,
      const GEN::pairedvector<int, double>& dmxigp,
      const GEN::pairedvector<int, double>& derivjac,
      const std::vector<GEN::pairedvector<int, double> >& ximaps,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap);

  /*!
   \brief lin D/M-matrix entries at GP for elebased integration

   */
  void inline GP_2D_DM_Ele_Lin(
      int& iter,
      bool& bound,
      MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& mderiv,
      double& dxdsxi,
      double& wgt,
      const GEN::pairedvector<int, double>& dmxigp,
      const GEN::pairedvector<int, double>& derivjac,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap);

  /*!
   \brief lin D/M-matrix entries at GP for bound case

   */
  void inline GP_3D_DM_Lin(
      int& iter,
      bool& duallin,
      MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& mderiv,
      LINALG::SerialDenseMatrix& lmderiv,
      double& wgt,
      double& jac,
      std::vector<GEN::pairedvector<int, double> >& dsxigp,
      std::vector<GEN::pairedvector<int, double> >& dmxigp,
      GEN::pairedvector<int, double>& jacintcellmap,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap);

  /*!
   \brief lin D/M-matrix entries at GP for elebased integration

   */
  void inline GP_3D_DM_Ele_Lin(
      int& iter,
      bool& duallin,
      bool& dualquad3d,
      MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& svalmod,
      LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& mderiv,
      double& wgt,
      double& jacslave,
      std::vector<GEN::pairedvector<int, double> >& dmxigp,
      GEN::pairedvector<int, double>& jacslavemap,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap);
  /*!
   \brief lin D/M-matrix entries at GP for bound case (3D quad)

   */
  void inline GP_3D_DM_Quad_Lin(int& iter, bool& duallin,
      MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& svalmod,
      LINALG::SerialDenseVector& mval, LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseMatrix& mderiv,
      LINALG::SerialDenseMatrix& lmderiv, double& wgt, double& jac,
      const std::vector<GEN::pairedvector<int, double> >& dpsxigp,
      const std::vector<GEN::pairedvector<int, double> >& dpmxigp,
      const GEN::pairedvector<int, double>& jacintcellmap,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap,
      bool dualquad3d);

  void inline GP_3D_DM_Quad_pwlin_Lin(int& iter, MORTAR::MortarElement& sele,
      MORTAR::MortarElement& sintele, MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseVector& lmintval, LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& mderiv, LINALG::SerialDenseMatrix& lmintderiv,
      double& wgt, double& jac,
      const std::vector<GEN::pairedvector<int, double> >& dsxigp,
      const std::vector<GEN::pairedvector<int, double> >& dpsxigp,
      const std::vector<GEN::pairedvector<int, double> >& dpmxigp,
      const GEN::pairedvector<int, double>& jacintcellmap);

  /*!
   \brief evaluate weighted Gap entries at GP

   */
  void inline GP_2D_G(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseVector& lmval, LINALG::SerialDenseMatrix& scoord,
      LINALG::SerialDenseMatrix& mcoord, LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& mderiv, double* gap, double* gpn,
      double* lengthn, double& dsxideta, double& dxdsxi, double& wgt,
      const GEN::pairedvector<int, double> & dsxigp,
      const GEN::pairedvector<int, double> & dmxigp,
      GEN::pairedvector<int, double> & dgapgp,
      std::vector<GEN::pairedvector<int, double> >& dnmap_unit,
      int& linsize);

  /*!
   \brief evaluate weighted Gap entries at GP

   */
  void inline GP_3D_G(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseVector& lmval, LINALG::SerialDenseMatrix& scoord,
      LINALG::SerialDenseMatrix& mcoord, LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& mderiv, double* gap, double* gpn,
      double* lengthn, double& jac, double& wgt,
      std::vector<GEN::pairedvector<int, double> >& dsxigp,
      std::vector<GEN::pairedvector<int, double> >& dmxigp,
      GEN::pairedvector<int, double> & dgapgp,
      std::vector<GEN::pairedvector<int, double> >& dnmap_unit,
      bool quadratic, int nintrow = 0);

  /*!
   \brief evaluate weighted Gap entries at GP (quad-pwlin)

   */
  void inline GP_3D_G_Quad_pwlin(MORTAR::MortarElement& sele,
      MORTAR::IntElement& sintele, MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseVector& lmintval, LINALG::SerialDenseMatrix& scoord,
      LINALG::SerialDenseMatrix& mcoord, LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& mderiv, double* gap, double* gpn,
      double* lengthn, double& jac, double& wgt,
      const std::vector<GEN::pairedvector<int, double> >& dsxigp,
      const std::vector<GEN::pairedvector<int, double> >& dmxigp,
      GEN::pairedvector<int, double> & dgapgp,
      std::vector<GEN::pairedvector<int, double> >& dnmap_unit);

  /*!
   \brief evaluate weighted Gap entries at GP

   */
  void inline GP_2D_G_Lin(
      int& iter,
      MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& lmderiv,
      double& gap,
      double* gpn,
      double& dsxideta,
      double& dxdsxi,
      double& dxdsxidsxi,
      double& wgt,
      const GEN::pairedvector<int, double>& dgapgp,
      const GEN::pairedvector<int, double>& dsxigp,
      const GEN::pairedvector<int, double>& dmxigp,
      const GEN::pairedvector<int, double>& derivjac,
      const std::vector<GEN::pairedvector<int, double> >& ximaps,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap);

  /*!
   \brief evaluate weighted Gap entries at GP

   */
  void inline GP_2D_G_Ele_Lin(
      int& iter,
      MORTAR::MortarElement& sele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseVector& lmval,
      double& gap,
      double& dxdsxi,
      double& wgt,
      const GEN::pairedvector<int, double>& dgapgp,
      const GEN::pairedvector<int, double>& derivjac,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap);

  /*!
   \brief evaluate weighted Gap entries at GP

   */
  void inline GP_3D_G_Lin(
      int& iter,
      MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& lmderiv,
      double& gap,
      double *gpn,
      double& jac,
      double& wgt,
      bool& duallin,
      GEN::pairedvector<int, double>& dgapgp,
      GEN::pairedvector<int, double>& jacintcellmap,
      std::vector<GEN::pairedvector<int, double> >& dsxigp,
      std::vector<GEN::pairedvector<int, double> >& dmxigp,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap);

  /*!
   \brief evaluate weighted Gap entries at GP (quad)

   */
  void inline GP_3D_G_Quad_Lin(int& iter, MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& svalmod, LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseMatrix& lmderiv,
      double& gap, double *gpn, double& jac, double& wgt, bool& duallin,
      const GEN::pairedvector<int, double>& dgapgp,
      const GEN::pairedvector<int, double>& jacintcellmap,
      const std::vector<GEN::pairedvector<int, double> >& dpsxigp,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap,
      bool dualquad3d);

  /*!
   \brief evaluate lin. weighted Gap entries at GP for ele-based
   integration

   */
  void inline GP_3D_G_Ele_Lin(
      int& iter,
      MORTAR::MortarElement& sele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& svalmod,
      LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& lmderiv,
      double& gap,
      double& jacslave,
      double& wgt,
      bool& duallin,
      bool& dualquad3d,
      GEN::pairedvector<int, double>& dgapgp,
      GEN::pairedvector<int, double>& jacslavemap,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap);

  /*!
   \brief evaluate weighted Gap entries at GP (quad)

   */
  void inline GP_3D_G_Quad_pwlin_Lin(int& iter, MORTAR::IntElement& sintele,
      LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& lmintval,
      LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseMatrix& lmintderiv,
      double& gap, double *gpn, double& jac, double& wgt,
      const GEN::pairedvector<int, double>& dgapgp,
      const GEN::pairedvector<int, double>& jacintcellmap,
      const std::vector<GEN::pairedvector<int, double> >& dsxigp);

  /*!
   \brief evaluate scaling at GP

   */
  void inline GP_2D_Scaling(MORTAR::MortarElement& sele,
      LINALG::SerialDenseVector& sval, double& dsxideta, double& wgt);

  /*!
   \brief evaluate and lin scaling at GP

   */
  void inline GP_2D_Scaling_Lin(int& iter, MORTAR::MortarElement& sele,
      LINALG::SerialDenseVector& sval, LINALG::SerialDenseMatrix& sderiv,
      double& dsxideta, double& wgt,
      const GEN::pairedvector<int, double>& dsxigp,
      const std::vector<GEN::pairedvector<int, double> >& ximaps);

  /*!
   \brief evaluate scaling at GP

   */
  void inline GP_3D_Scaling(MORTAR::MortarElement& sele,
      LINALG::SerialDenseVector& sval, double& jac, double& wgt, double* sxi);

  /*!
   \brief evaluate and lin scaling at GP

   */
  void inline GP_3D_Scaling_Lin(int& iter, MORTAR::MortarElement& sele,
      LINALG::SerialDenseVector& sval, LINALG::SerialDenseMatrix& sderiv,
      double& jac, double& wgt, double& jacsele,
      const GEN::pairedvector<int, double>& derivjacsele,
      const GEN::pairedvector<int, double>& jacintcellmap,
      const std::vector<GEN::pairedvector<int, double> >& dsxigp,
      double* derivjacselexi);

  /*!
   \brief evaluate and lin slipincr at GP

   */
  void inline GP_2D_SlipIncr(MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& mval, LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& scoord, LINALG::SerialDenseMatrix& mcoord,
      Teuchos::RCP<LINALG::SerialDenseMatrix> scoordold,
      Teuchos::RCP<LINALG::SerialDenseMatrix> mcoordold,
      LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseMatrix& mderiv,
      double& dsxideta, double& dxdsxi, double& wgt, double* jumpvalv,
      const GEN::pairedvector<int, double>& dsxigp,
      const GEN::pairedvector<int, double>& dmxigp,
      GEN::pairedvector<int, double>& dslipgp,
      int& linsize);

  /*!
   \brief evaluate and lin slipincr at GP

   */
  void inline GP_3D_SlipIncr(MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& mval, LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& scoord, LINALG::SerialDenseMatrix& mcoord,
      Teuchos::RCP<LINALG::SerialDenseMatrix> scoordold,
      Teuchos::RCP<LINALG::SerialDenseMatrix> mcoordold,
      LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseMatrix& mderiv,
      double& jac, double& wgt, double* jumpvalv,
      const std::vector<GEN::pairedvector<int, double> >& dsxigp,
      const std::vector<GEN::pairedvector<int, double> >& dmxigp,
      std::vector<GEN::pairedvector<int, double> >& dslipgp);

  /*!
   \brief evaluate and lin slipincr at GP at node

   */
  void inline GP_2D_SlipIncr_Lin(
      int& iter,
      MORTAR::MortarElement& sele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& lmderiv,
      double& dsxideta,
      double& dxdsxi,
      double& dxdsxidsxi,
      double& wgt,
      double* jumpvalv,
      const GEN::pairedvector<int, double>& dsxigp,
      const GEN::pairedvector<int, double>& dslipgp,
      const std::vector<GEN::pairedvector<int, double> >& ximaps,
      const GEN::pairedvector<int, double>& derivjac,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap);

  void inline GP_3D_SlipIncr_Lin(
      int& iter,
      MORTAR::MortarElement& sele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& lmderiv,
      double& jac,
      double& wgt,
      double* jumpvalv,
      const GEN::pairedvector<int, double>& jacintcellmap,
      const std::vector<GEN::pairedvector<int, double> >& dslipgp,
      const std::vector<GEN::pairedvector<int, double> >& dsxigp,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap);
  /*!
   \brief evaluate  T and E matrix

   */
  void inline GP_TE(MORTAR::MortarElement& sele,
      LINALG::SerialDenseVector& lmval, LINALG::SerialDenseVector& sval,
      double& jac, double& wgt, double* jumpval);

  /*!
   \brief evaluate  T and E matrix

   */
  void inline GP_TE_Master(MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele, LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseVector& lm2val, LINALG::SerialDenseVector& sval,
      double& jac, double& wgt, double* jumpval, const Epetra_Comm& comm);

  /*!
   \brief evaluate Lin T and E matrix

   */
  void inline GP_2D_TE_Lin(
      int& iter, //like k
      MORTAR::MortarElement& sele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& lmderiv,
      double& dsxideta,
      double& dxdsxi,
      double& dxdsxidsxi,
      double& wgt,
      double* jumpval,
      const GEN::pairedvector<int, double>& dsxigp,
      const GEN::pairedvector<int, double>& derivjac,
      const GEN::pairedvector<int, double>& dsliptmatrixgp,
      const std::vector<GEN::pairedvector<int, double> >& ximaps,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap);

  /*!
   \brief evaluate Lin T and E matrix

   */
  void inline GP_2D_TE_Master_Lin(
      int& iter, //like k
      MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseVector& lmval, LINALG::SerialDenseMatrix& mderiv,
      LINALG::SerialDenseMatrix& lmderiv, double& dsxideta, double& dxdsxi,
      double& dxdsxidsxi, double& wgt, double* jumpval,
      const GEN::pairedvector<int, double>& dsxigp,
      const GEN::pairedvector<int, double>& dmxigp,
      const GEN::pairedvector<int, double>& derivjac,
      const GEN::pairedvector<int, double>& dsliptmatrixgp,
      const std::vector<GEN::pairedvector<int, double> >& ximaps,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap,
      const Epetra_Comm& comm);

  /*!
   \brief evaluate Lin T and E matrix

   */
  void inline GP_3D_TE_Lin(
      int& iter,
      bool& duallin, //like k
      MORTAR::MortarElement& sele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& lmderiv,
      double& jac,
      double& wgt,
      double* jumpval,
      const std::vector<GEN::pairedvector<int, double> >& dsxigp,
      const GEN::pairedvector<int, double>& jacintcellmap,
      const GEN::pairedvector<int, double>& dsliptmatrixgp,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap);

  /*!
   \brief evaluate Lin T and E matrix (Master)

   */
  void inline GP_3D_TE_Master_Lin(
      int& iter,
      bool& duallin, //like k
      MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& sval, LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseVector& lmval, LINALG::SerialDenseVector& lm2val,
      LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseMatrix& mderiv,
      LINALG::SerialDenseMatrix& lmderiv, LINALG::SerialDenseMatrix& lm2deriv,
      double& jac, double& wgt, double* jumpval,
      const std::vector<GEN::pairedvector<int, double> >& dsxigp,
      const std::vector<GEN::pairedvector<int, double> >& dmxigp,
      const GEN::pairedvector<int, double>& jacintcellmap,
      const GEN::pairedvector<int, double>& dsliptmatrixgp,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dual2map,
      const Epetra_Comm& comm);
  /*!
   \brief evaluate matrix A for tsi-contact

   */
  void inline GP_TSI_A(MORTAR::MortarElement& sele,
      LINALG::SerialDenseVector& lmval, double& jac, double& wgt, int& nrow,
      int& ncol, int& ndof);

  /*!
   \brief evaluate matrix B for tsi-contact

   */
  void inline GP_TSI_B(MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& mval, double& jac, double& wgt, int& ncol,
      int& ndof);

  /*!
   \brief evaluate mechanical dissipation(tsi-contact)

   */
  void inline GP_TSI_MechDiss(MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& mval, LINALG::SerialDenseVector& lmval,
      double& jac, double& mechdiss, double& wgt, int& nrow, int& ncol,
      int& ndof, bool& thermolagmult);

  /*!
   \brief evaluate wear + lin at GP

   */
  void inline GP_2D_Wear(MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseMatrix& mderiv, LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& lmderiv, LINALG::SerialDenseMatrix& scoord,
      Teuchos::RCP<LINALG::SerialDenseMatrix> scoordold,
      LINALG::SerialDenseMatrix& mcoord,
      Teuchos::RCP<LINALG::SerialDenseMatrix> mcoordold,
      Teuchos::RCP<LINALG::SerialDenseMatrix> lagmult, double* gpn,
      double& dsxideta, double& dxdsxi, double& dxdsxidsxi, double& wgt,
      double* jumpval, double* wearval,
      const GEN::pairedvector<int, double>& dsxigp,
      const GEN::pairedvector<int, double>& dmxigp,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap,
      const std::vector<GEN::pairedvector<int, double> >& ximaps,
      const std::vector<GEN::pairedvector<int, double> >& dnmap_unit,
      GEN::pairedvector<int, double> & dsliptmatrixgp,
      GEN::pairedvector<int, double> & dweargp,
      int& linsize);

  /*!
   \brief evaluate wear + lin at GP

   */
  void inline GP_3D_Wear(MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele, LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseMatrix& sderiv, LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseMatrix& mderiv, LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& lmderiv, LINALG::SerialDenseMatrix& scoord,
      Teuchos::RCP<LINALG::SerialDenseMatrix> scoordold,
      LINALG::SerialDenseMatrix& mcoord,
      Teuchos::RCP<LINALG::SerialDenseMatrix> mcoordold,
      Teuchos::RCP<LINALG::SerialDenseMatrix> lagmult, double* gpn, double& jac,
      double& wgt, double* jumpval, double* wearval,
      GEN::pairedvector<int, double> & dsliptmatrixgp,
      GEN::pairedvector<int, double> & dweargp,
      const std::vector<GEN::pairedvector<int, double> >& dsxigp,
      const std::vector<GEN::pairedvector<int, double> >& dmxigp,
      const std::vector<GEN::pairedvector<int, double> >& dnmap_unit,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap,
      double& mechdiss);

  /*!
   \brief lin weighted wear at GP

   */
  void inline GP_2D_Wear_Lin(
      int& iter,
      MORTAR::MortarElement& sele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& lmderiv,
      double& dsxideta,
      double& dxdsxi,
      double& dxdsxidsxi,
      double* gpn,
      double& wgt,
      double& wearval,
      double* jumpval,
      const GEN::pairedvector<int, double>& dsxigp,
      const GEN::pairedvector<int, double>& dweargp,
      const std::vector<GEN::pairedvector<int, double> >& ximaps,
      const GEN::pairedvector<int, double>& derivjac,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap);

  /*!
   \brief lin weighted wear at GP

   */
  void inline GP_3D_Wear_Lin(
      int& iter,
      MORTAR::MortarElement& sele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& lmderiv,
      double& jac,
      double* gpn,
      double& wgt,
      double& wearval,
      double* jumpval,
      const GEN::pairedvector<int, double>& dweargp,
      const GEN::pairedvector<int, double>& jacintcellmap,
      const std::vector<GEN::pairedvector<int, double> >& dsxigp,
      const std::vector<std::vector<GEN::pairedvector<int, double> > >& dualmap);


  /*!
  \brief evaluate scalar normal coupling condition for poro no penetration entries at GP (poro-contact)

  */
  void inline GP_3D_NCOUP_DERIV(
       MORTAR::MortarElement& sele,
       MORTAR::MortarElement& mele,
       LINALG::SerialDenseVector& sval,
       LINALG::SerialDenseVector& mval,
       LINALG::SerialDenseVector& lmval,
       LINALG::SerialDenseMatrix& sderiv,
       LINALG::SerialDenseMatrix& mderiv,
       double* ncoup, double* gpn, double* lengthn,
       double& jac,
       double& wgt,
       const std::vector<GEN::pairedvector<int,double> >& dsxigp,
       const std::vector<GEN::pairedvector<int,double> >& dmxigp,
       std::map<int,double> & dncoupgp,
       std::map<int,double> & dvelncoupgp,
       std::vector<GEN::pairedvector<int,double> >& dnmap_unit, bool quadratic,
       int nintrow=0);

  /*!
  \brief evaluate weighted normal coupling entries at GP

  */
  void inline GP_3D_NCOUP_LIN(
       int& iter,
       MORTAR::MortarElement& sele,
       MORTAR::MortarElement& mele,
       LINALG::SerialDenseVector& sval,
       LINALG::SerialDenseVector& mval,
       LINALG::SerialDenseVector& lmval,
       LINALG::SerialDenseMatrix& sderiv,
       LINALG::SerialDenseMatrix& lmderiv,
       double& gap, double *gpn,double& jac,
       double& wgt, bool& duallin,
       const std::map<int,double>& dncoupgp,
       const std::map<int,double>& dvelncoupgp,
       const GEN::pairedvector<int,double>& jacintcellmap,
       const std::vector<GEN::pairedvector<int,double> >& dsxigp,
       const std::vector<GEN::pairedvector<int,double> >& dmxigp,
       const std::vector<std::vector<GEN::pairedvector<int,double> > >& dualmap);

  /*!
   \brief Build all integrals and linearizations without segmentation -- 2D
   (i.e. M, g, LinM, Ling and possibly D, LinD)

   */
  void IntegrateDerivEle2D(MORTAR::MortarElement& sele,
      std::vector<MORTAR::MortarElement*> meles, bool *boundary_ele);
  /*!
   \brief Build all integrals and linearizations on a 1D slave /
   master overlap (i.e. M, g, LinM, Ling and possibly D, LinD and
   wear)

   */
  void IntegrateDerivSegment2D(MORTAR::MortarElement& sele, double& sxia,
      double& sxib, MORTAR::MortarElement& mele, double& mxia, double& mxib,
      const Epetra_Comm& comm);

  /*!
   \brief Build all integrals and linearizations without segmentation -- 3D
   (i.e. M, g, LinM, Ling and possibly D, LinD)

   */
  void IntegrateDerivEle3D(MORTAR::MortarElement& sele,
      std::vector<MORTAR::MortarElement*> meles, bool *boundary_ele,
      bool *proj_, const Epetra_Comm& comm);

  /*!
   \brief Build all integrals and linearizations on a 2D slave /
   master integration cell (i.e. M, g, LinM, Ling and possibly D, LinD)
   for the auxiliary plane coupling case

   */
  void IntegrateDerivCell3DAuxPlane(MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele, Teuchos::RCP<MORTAR::IntCell> cell,
      double* auxn, const Epetra_Comm& comm);

  /*!
   \brief Build all integrals and linearizations on a 2D slave /
   master integration cell (i.e. M, g, LinM, Ling) for
   the auxiliary plane coupling case with quadratic interpolation

   */
  void IntegrateDerivCell3DAuxPlaneQuad(MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele, MORTAR::IntElement& sintele,
      MORTAR::IntElement& mintele, Teuchos::RCP<MORTAR::IntCell> cell,
      double* auxn);

  /*!
   \brief Compute penalty scaling factor kappa on slave element

   */
  void IntegrateKappaPenalty(MORTAR::MortarElement& sele, double* sxia,
      double* sxib, Teuchos::RCP<Epetra_SerialDenseVector> gseg);

  /*!
   \brief Compute penalty scaling factor kappa on slave integration element
   (special version for the 3D quadratic case)

   */
  void IntegrateKappaPenalty(MORTAR::MortarElement& sele,
      MORTAR::IntElement& sintele, double* sxia, double* sxib,
      Teuchos::RCP<Epetra_SerialDenseVector> gseg);

  //@}

  //! @name 2D and 3D linearization methods

  /*!
   \brief Compute directional derivative of segment end coordinates
   Xi on a 1D slave / master overlap

   */
  void DerivXiAB2D(MORTAR::MortarElement& sele, double& sxia, double& sxib,
      MORTAR::MortarElement& mele, double& mxia, double& mxib,
      std::vector<GEN::pairedvector<int, double> >& derivxi,
      bool& startslave, bool& endslave,
      int& linsize);

  /*!
   \brief Compute directional derivative of master Gauss point
   coordinates XiGP on a 1D slave / master overlap

   */
  void DerivXiGP2D(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
      double& sxigp, double& mxigp,
      const GEN::pairedvector<int, double>& derivsxi,
      GEN::pairedvector<int, double>& derivmxi,
      int& linsize);

  /*!
   \brief Compute directional derivative of master Gauss point
   coordinates XiGP on a 2D slave / master integration cell

   */
  void DerivXiGP3D(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
      double* sxigp, double* mxigp,
      const std::vector<GEN::pairedvector<int, double> >& derivsxi,
      std::vector<GEN::pairedvector<int, double> >& derivmxi,
      double& alpha);

  /*!
   \brief Compute directional derivative of slave / master Gauss point
   coordinates XiGP on a 2D slave / master integration cell
   (This is the AuxPlane version, thus master and slave are projected)

   */
  void DerivXiGP3DAuxPlane(MORTAR::MortarElement& ele, double* xigp,
      double* auxn, std::vector<GEN::pairedvector<int, double> >& derivxi
      , double& alpha,
      std::vector<GEN::pairedvector<int, double> >& derivauxn,
      std::vector<GEN::pairedvector<int, double> >& derivgp);

  /*!
   \brief Assemble D contribution of current overlap into slave nodes

   */
  bool AssembleD(const Epetra_Comm& comm, MORTAR::MortarElement& sele,
      Epetra_SerialDenseMatrix& D_seg);

  /*!
   \brief Assemble scale factor contribution of current overlap into slave nodes

   */
  bool AssembleScale(const Epetra_Comm& comm, MORTAR::MortarElement& sele,
      Epetra_SerialDenseVector& scseg);

  /*!
   \brief Assemble T contribution of current overlap into slave nodes

   */
  bool AssembleT(const Epetra_Comm& comm, MORTAR::MortarElement& sele,
      Epetra_SerialDenseMatrix& T_seg);

  /*!
   \brief Assemble E contribution of current overlap into slave nodes

   */
  bool AssembleE(const Epetra_Comm& comm, MORTAR::MortarElement& sele,
      Epetra_SerialDenseMatrix& E_seg);

  /*!
   \brief Assemble Ws contribution of current overlap into slave nodes

   */
  bool AssembleWS(const Epetra_Comm& comm, MORTAR::MortarElement& sele,
      Epetra_SerialDenseMatrix& Ws_seg);
  /*!
   \brief Assemble D2 contribution of current overlap into master nodes

   */
  bool AssembleD2(const Epetra_Comm& comm, MORTAR::MortarElement& mele,
      Epetra_SerialDenseMatrix& D2_seg);

  /*!
   \brief Assemble D contribution of current overlap into slave nodes
   (special version for 3D quadratic mortar with piecewise linear LM interpolation)

   */
  virtual bool AssembleD(const Epetra_Comm& comm, MORTAR::MortarElement& sele,
      MORTAR::IntElement& sintele, Epetra_SerialDenseMatrix& D_seg);

  /*!
   \brief Assemble M contribution of current overlap into slave nodes

   */
  bool AssembleM(const Epetra_Comm& comm, MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele, Epetra_SerialDenseMatrix& M_seg);

  /*!
   \brief Assemble M contribution of current overlap into slave nodes
   (special version for 3D quadratic mortar with piecewise linear LM interpolation)

   */
  virtual bool AssembleM(const Epetra_Comm& comm, MORTAR::IntElement& sintele,
      MORTAR::MortarElement& mele, Epetra_SerialDenseMatrix& M_seg);

  /*!
   \brief Assemble M contribution of current slave-master pairing for EleBased int.

   */
  virtual bool AssembleM_EleBased(const Epetra_Comm& comm,
      MORTAR::MortarElement& sele, std::vector<MORTAR::MortarElement*> meles,
      Epetra_SerialDenseMatrix& M_seg);
  /*!
   \brief Assemble g~ contribution of current overlap into slave nodes

   */
  bool AssembleG(const Epetra_Comm& comm, MORTAR::MortarElement& sele,
      Epetra_SerialDenseVector& g_seg);

  /*!
   \brief Assemble u~ contribution of current overlap into slave nodes
   --> this is the relative tangential weighted slip integrated at the gp

   */
  bool AssembleU(const Epetra_Comm& comm, MORTAR::MortarElement& sele,
      Epetra_SerialDenseMatrix& u_seg);

  /*!
   \brief Assemble g~ contribution of current overlap into slave nodes
   (special version for 3D quadratic mortar with piecewise linear LM interpolation)

   */
  bool AssembleG(const Epetra_Comm& comm, MORTAR::IntElement& sintele,
      Epetra_SerialDenseVector& g_seg);

  /*!
   \brief Assemble mechanical dissipation contribution of current overlap
   into slave nodes

   */
  bool AssembleMechDissSlave(const Epetra_Comm& comm,
      MORTAR::MortarElement& sele, Epetra_SerialDenseVector& mdissseg);
  /*!
   \brief Assemble mechanical dissipation contribution of current overlap
   into master nodes

   */
  bool AssembleMechDissMaster(const Epetra_Comm& comm,
      MORTAR::MortarElement& sele, Epetra_SerialDenseVector& mdissseg);

  /*!
   \brief Assemble A contribution of current overlap into slave nodes

   */
  bool AssembleA(const Epetra_Comm& comm, MORTAR::MortarElement& sele,
      Epetra_SerialDenseMatrix& aseg);

  /*!
   \brief Assemble B contribution of current overlap into master nodes

   */
  bool AssembleB(const Epetra_Comm& comm, MORTAR::MortarElement& mele,
      Epetra_SerialDenseMatrix& bseg);

  /*!
   \brief Assemble wear contribution of current overlap into slave nodes

   */
  bool AssembleWear(const Epetra_Comm& comm, MORTAR::MortarElement& sele,
      Epetra_SerialDenseVector& wseg);

  /*!
   \brief Return the Wear shape fcn type (wear weighting...)

   */
  INPAR::CONTACT::WearShape WearShapeFcn()
  {
    return wearshapefcn_;
  }

  /*!
   \brief Return type of wear surface definition

   */
  INPAR::CONTACT::WearSide WearSide()
  {
    return wearside_;
  }

  /*!
   \brief Return type of wear algorithm

   */
  INPAR::CONTACT::WearType WearType()
  {
    return weartype_;
  }

  /*!
   \brief Return the LM shape fcn type

   */
  INPAR::MORTAR::ShapeFcn ShapeFcn()
  {
    return shapefcn_;
  }

  /*!
   \brief Return the LM interpolation / testing type for quadratic FE

   */
  INPAR::MORTAR::LagMultQuad LagMultQuad()
  {
    return lagmultquad_;
  }
  //@}

  // GP calls
  /*!
   \brief Return number of Gauss points for this instance

   */
  int& nGP()
  {
    return ngp_;
  }

  /*!
   \brief Return coordinates of a specific GP in 1D/2D CElement

   */
  double& Coordinate(int& gp, int dir)
  {
    return coords_(gp, dir);
  }

  /*!
   \brief Return weight of a specific GP in 1D/2D CElement

   */
  double& Weight(int& gp)
  {
    return weights_[gp];
  }

  /*!
   \brief Get problem dimension

   Note that only 2D and 3D are possible here as this refers to the global
   problem dimension. On integration level this corresponds to 1D integration
   (dim_==2) and 2D integration (dim_==3) on the interface!

   */
  virtual const int& Dim()
  {
    return dim_;
  }
  ;

protected:

    // don't want = operator and cctor
    CoIntegrator operator = (const CoIntegrator& old);
    CoIntegrator(const CoIntegrator& old);

    Teuchos::ParameterList&                          imortar_;     // containing contact input parameters
    DRT::Node**                                      mynodes_;     // slave element nodes
    LINALG::SerialDenseMatrix                        scoord_;      // slave element nodal coordinates
    LINALG::SerialDenseMatrix                        mcoord_;      // master element nodal coordinates
    LINALG::SerialDenseMatrix                        ssecderiv_;   // storage for 2nd derivatives of standard shape fct.
    bool                                             duallin_;     // flag indicating linearization of dual shape fct.
    std::vector<std::vector<std::map<int,double> > > dualmap_;     // storage for dual shape fct. derivatives
    std::vector<std::map<int,double> >               ximaps_;      // storage for 2d segment endpoints derivatives
    const Epetra_Comm&                               Comm_;

    int                                              ngp_;         // number of Gauss points
    Epetra_SerialDenseMatrix                         coords_;      // Gauss point coordinates
    std::vector<double>                              weights_;     // Gauss point weights
    int                                              dim_;         // dimension of problem (2D or 3D)

    // inputs from parameter list
    INPAR::MORTAR::ShapeFcn                          shapefcn_;    // lm shape function type
    INPAR::MORTAR::LagMultQuad                       lagmultquad_; // type of lm interpolation for quadr. FE
    bool                                             nodalscale_;  // flag for nodal scaling
    bool                                             gpslip_;      // gp-wise evaluated slip increment

    // wear inputs from parameter list
    INPAR::CONTACT::WearLaw                          wearlaw_;     // type of wear law
    bool                                             wearimpl_;    // flag for implicit wear algorithm
    INPAR::CONTACT::WearSide                         wearside_;    // definition of wear surface
    INPAR::CONTACT::WearType                         weartype_;    // definition of contact wear algorithm
    INPAR::CONTACT::WearShape                        wearshapefcn_;// type of wear shape function
    double                                           wearcoeff_;   // wear coefficient
    double                                           wearcoeffm_;  // wear coefficient master
}; // class CoIntegrator
}  // namespace CONTACT


#endif  // #ifndef CONTACT_INTEGRATOR_H

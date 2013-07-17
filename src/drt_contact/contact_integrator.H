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


/*!
\brief CONTACT: namespace of the BACI contact and meshtying algorithms

*/
namespace CONTACT
{

// forward declarations

/*!
\brief A class to perform Gaussian integration and assembly of Mortar
       matrices on the overlap of two MortarElements (1 Slave, 1 Master)
       in 1D (which is equivalent to a 2D coupling problem) and in 2D
       (which is equivalent to a 3D coupling problem).
       
       This is a derived class from MORTAR::MortarIntegrator which does
       the contact-specific stuff for 3d mortar coupling.

\author popp (popp@lnm.mw.tum.de)
*/

class CoIntegrator : public MORTAR::MortarIntegrator
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

  /*!
  \brief Destructor

  */
  virtual ~CoIntegrator() {}

  //! @name 2D and 3D integration methods

  //! @name 2D and 3D integration methods

  /*!
  \brief Build all integrals and linearizations on a 1D or 2D
         slave element (i.e. onl D, LinD possible).

  This method is outdated!!! We now obtain the D-matrix by
  integrating the mortar segments (just like the M-matrix entries).

  */
  void IntegrateDerivSlave2D3D(
      MORTAR::MortarElement& sele, double* sxia, double* sxib,
      Teuchos::RCP<Epetra_SerialDenseMatrix> dseg);

  /*!
  \brief Build all integrals and linearizations without segmentation -- 2D
         (i.e. M, g, LinM, Ling and possibly D, LinD)

  */
  void EleBased_Integration(
      MORTAR::MortarElement& sele,
      std::vector<MORTAR::MortarElement*> meles,
      Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
      Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
      Teuchos::RCP<Epetra_SerialDenseVector> gseg,
      Teuchos::RCP<Epetra_SerialDenseVector> wseg,
      bool *boundary_ele);
  /*!
  \brief Build all integrals and linearizations on a 1D slave /
         master overlap (i.e. M, g, LinM, Ling and possibly D, LinD and
         wear)

  */
  void IntegrateDerivSegment2D(
      MORTAR::MortarElement& sele, double& sxia, double& sxib,
      MORTAR::MortarElement& mele, double& mxia, double& mxib,
      Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
      Teuchos::RCP<Epetra_SerialDenseMatrix> d2seg,
      Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
      Teuchos::RCP<Epetra_SerialDenseVector> gseg,
      Teuchos::RCP<Epetra_SerialDenseVector> wseg);

  /*!
  \brief Build all integrals and linearizations on a 2D slave /
         master integration cell (i.e. M, g, LinM, Ling and possibly D, LinD)

  */
  void IntegrateDerivCell3D(
      MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
      Teuchos::RCP<MORTAR::IntCell> cell,
      Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
      Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
      Teuchos::RCP<Epetra_SerialDenseVector> gseg);
    
  /*!
  \brief Build all integrals and linearizations without segmentation -- 3D
         (i.e. M, g, LinM, Ling and possibly D, LinD)

  */
  void IntegrateDerivCell3D_EleBased(
      MORTAR::MortarElement& sele, std::vector<MORTAR::MortarElement*> meles,
      Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
      Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
      Teuchos::RCP<Epetra_SerialDenseVector> gseg,
      INPAR::MORTAR::LagMultQuad lmtype,
      bool *boundary_ele,
      bool *proj_);

  /*!
  \brief Build all integrals and linearizations on a 2D slave /
         master integration cell (i.e. M, g, LinM, Ling and possibly D, LinD)
         for the auxiliary plane coupling case

  */
  void IntegrateDerivCell3DAuxPlane(
     MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
     Teuchos::RCP<MORTAR::IntCell> cell, double* auxn,
     Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
     Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
     Teuchos::RCP<Epetra_SerialDenseVector> gseg,
     Teuchos::RCP<Epetra_SerialDenseVector> mdisssegs,
     Teuchos::RCP<Epetra_SerialDenseVector> mdisssegm,
     Teuchos::RCP<Epetra_SerialDenseMatrix> aseg,
     Teuchos::RCP<Epetra_SerialDenseMatrix> bseg,
     Teuchos::RCP<Epetra_SerialDenseVector> wseg);

  /*!
  \brief Build all integrals and linearizations on a 2D slave /
         master integration cell (i.e. M, g, LinM, Ling) for
         the auxiliary plane coupling case with quadratic interpolation

  */
  void IntegrateDerivCell3DAuxPlaneQuad(
     MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
     MORTAR::IntElement& sintele, MORTAR::IntElement& mintele,
     Teuchos::RCP<MORTAR::IntCell> cell, double* auxn,
     Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
     Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
     Teuchos::RCP<Epetra_SerialDenseVector> gseg);
  
  /*!
  \brief Compute penalty scaling factor kappa on slave element

  */
  void IntegrateKappaPenalty(MORTAR::MortarElement& sele,
                             double* sxia, double* sxib,
                             Teuchos::RCP<Epetra_SerialDenseVector> gseg);
    
  /*!
  \brief Compute penalty scaling factor kappa on slave integration element
  (special version for the 3D quadratic case)

  */
  void IntegrateKappaPenalty(MORTAR::MortarElement& sele,
                             MORTAR::IntElement& sintele,
                             double* sxia, double* sxib,
                             Teuchos::RCP<Epetra_SerialDenseVector> gseg);
  
  //@}

  //! @name 2D and 3D linearization methods

  /*!
  \brief Compute directional derivative of segment end coordinates
         Xi on a 1D slave / master overlap

  */
  void DerivXiAB2D(MORTAR::MortarElement& sele, double& sxia, double& sxib,
                   MORTAR::MortarElement& mele, double& mxia, double& mxib,
                   std::vector<std::map<int,double> >& derivxi,
                   bool& startslave, bool& endslave);

  /*!
  \brief Compute directional derivative of master Gauss point
         coordinates XiGP on a 1D slave / master overlap

  */
  void DerivXiGP2D(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
                   double& sxigp, double& mxigp, const std::map<int,double>& derivsxi,
                   std::map<int,double>& derivmxi);
  
  /*!
  \brief Compute directional derivative of master Gauss point
         coordinates XiGP on a 2D slave / master integration cell

  */
  void DerivXiGP3D(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
                   double* sxigp, double* mxigp, const std::vector<std::map<int,double> >& derivsxi,
                   std::vector<std::map<int,double> >& derivmxi, double& alpha);

  /*!
  \brief Compute directional derivative of slave / master Gauss point
         coordinates XiGP on a 2D slave / master integration cell
         (This is the AuxPlane version, thus master and slave are projected)

  */
  void DerivXiGP3DAuxPlane(MORTAR::MortarElement& ele, double* xigp, double* auxn,
                   std::vector<std::map<int,double> >& derivxi, double& alpha,
                   const std::vector<std::map<int,double> >& derivauxn,
                   const std::vector<std::map<int,double> >& derivgp);
  
  /*!
  \brief Assemble D contribution of current overlap into slave nodes

  */
  bool AssembleD(const Epetra_Comm& comm,
                 MORTAR::MortarElement& sele,
                 Epetra_SerialDenseMatrix& D_seg);

  bool AssembleD2(const Epetra_Comm& comm,
                 MORTAR::MortarElement& mele,
                 Epetra_SerialDenseMatrix& D2_seg);

  /*!
  \brief Assemble D contribution of current overlap into slave nodes
         (special version for 3D quadratic mortar with piecewise linear LM interpolation)

  */
  virtual bool AssembleD(const Epetra_Comm& comm,
                 MORTAR::MortarElement& sele,
                 MORTAR::IntElement& sintele,
                 Epetra_SerialDenseMatrix& D_seg);
    
  /*!
  \brief Assemble M contribution of current overlap into slave nodes

  */
  bool AssembleM(const Epetra_Comm& comm,
                 MORTAR::MortarElement& sele,
                 MORTAR::MortarElement& mele,
                 Epetra_SerialDenseMatrix& M_seg);
  
  /*!
  \brief Assemble M contribution of current overlap into slave nodes
         (special version for 3D quadratic mortar with piecewise linear LM interpolation)

  */
  virtual bool AssembleM(const Epetra_Comm& comm,
                 MORTAR::IntElement& sintele,
                 MORTAR::MortarElement& mele,
                 Epetra_SerialDenseMatrix& M_seg);
  

  /*!
  \brief Assemble M contribution of current slave-master pairing for EleBased int.

  */
  virtual bool AssembleM_EleBased(const Epetra_Comm& comm,
                 MORTAR::MortarElement& sele,
                 std::vector<MORTAR::MortarElement*> meles,
                 Epetra_SerialDenseMatrix& M_seg);
  /*!
  \brief Assemble g~ contribution of current overlap into slave nodes

  */
  bool AssembleG(const Epetra_Comm& comm,
                 MORTAR::MortarElement& sele,
                 Epetra_SerialDenseVector& g_seg);
  
  /*!
  \brief Assemble g~ contribution of current overlap into slave nodes
         (special version for 3D quadratic mortar with piecewise linear LM interpolation)

  */
  bool AssembleG(const Epetra_Comm& comm,
                 MORTAR::IntElement& sintele,
                 Epetra_SerialDenseVector& g_seg);
  
  /*!
  \brief Assemble mechanical dissipation contribution of current overlap 
  into slave nodes 

  */
  bool AssembleMechDissSlave(const Epetra_Comm& comm,
                             MORTAR::MortarElement& sele,
                             Epetra_SerialDenseVector& mdissseg);
  /*!
  \brief Assemble mechanical dissipation contribution of current overlap 
  into master nodes 
  
  */
  bool AssembleMechDissMaster(const Epetra_Comm& comm,
                              MORTAR::MortarElement& sele,
                              Epetra_SerialDenseVector& mdissseg);

  /*!
  \brief Assemble A contribution of current overlap into slave nodes

  */
  bool AssembleA(const Epetra_Comm& comm,
                 MORTAR::MortarElement& sele,
                 Epetra_SerialDenseMatrix& aseg);

  /*!
  \brief Assemble B contribution of current overlap into master nodes

  */
  bool AssembleB(const Epetra_Comm& comm,
                 MORTAR::MortarElement& mele,
                 Epetra_SerialDenseMatrix& bseg);
  
  /*!
  \brief Assemble wear contribution of current overlap into slave nodes

  */
  bool AssembleWear(const Epetra_Comm& comm,
                    MORTAR::MortarElement& sele,
                    Epetra_SerialDenseVector& wseg);    

  //@}
   
protected:

    // don't want = operator and cctor
    CoIntegrator operator = (const CoIntegrator& old);
    CoIntegrator(const CoIntegrator& old);

    DRT::Node**                                      mynodes_;   // slave element nodes
    LINALG::SerialDenseMatrix                        scoord_;    // slave element nodal coordinates
    LINALG::SerialDenseMatrix                        mcoord_;    // master element nodal coordinates
    LINALG::SerialDenseMatrix                        ssecderiv_; // storage for 2nd derivatives of standard shape fct.
    bool                                             duallin_;   // flag indicating linearization of dual shape fct.
    bool                                             wear_both_;   // flag indicating linearization of dual shape fct.
    std::vector<std::vector<std::map<int,double> > > dualmap_;   // storage for dual shape fct. derivatives
    std::vector<std::map<int,double> >               ximaps_;    // storage for 2d segment endpoints derivatives

}; // class CoIntegrator
}  // namespace CONTACT


#endif  // #ifndef CONTACT_INTEGRATOR_H

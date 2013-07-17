/*!----------------------------------------------------------------------
\file mortar_integrator.H

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
#ifndef MORTAR_INTEGRATOR_H
#define MORTAR_INTEGRATOR_H

#include "Epetra_Comm.h"
#include "../drt_lib/drt_element.H"
#include "../drt_inpar/inpar_mortar.H"
#include "../drt_inpar/inpar_contact.H"


/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms

*/
namespace MORTAR
{

// forward declarations
class MortarElement;
class IntElement;
class IntCell;

/*!
\brief A class to perform Gaussian integration and assembly of Mortar
       matrices on the overlap of two MortarElements (1 Slave, 1 Master)
       in 1D (which is equivalent to a 2D coupling problem) and in 2D
       (which is equivalent to a 3D coupling problem)

\author popp (popp@lnm.mw.tum.de)
*/

class MortarIntegrator
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
  MortarIntegrator(Teuchos::ParameterList& params,
                   DRT::Element::DiscretizationType eletype);

  /*!
  \brief Destructor

  */
  virtual ~MortarIntegrator() {}

  /*!
  \brief Initialize Gauss rule (points, weights) for this MortarIntegrator

  */
  void InitializeGP(DRT::Element::DiscretizationType eletype);
    
  //! @name Access methods

  /*!
  \brief Get problem dimension

  Note that only 2D and 3D are possible here as this refers to the global
  problem dimension. On integration level this corresponds to 1D integration
  (dim_==2) and 2D integration (dim_==3) on the interface!

  */
  virtual const int& Dim() {return dim_;};

  /*!
  \brief Return the LM shape fcn type

  */
  INPAR::MORTAR::ShapeFcn ShapeFcn()
      { return DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(imortar_,"SHAPEFCN"); }

  /*!
  \brief Return type of wear surface definition

  */
  INPAR::CONTACT::WearSide WearSide()
      { return DRT::INPUT::IntegralValue<INPAR::CONTACT::WearSide>(imortar_,"BOTH_SIDED_WEAR"); }

  /*!
  \brief Return the LM interpolation / testing type for quadratic FE

  */
  INPAR::MORTAR::LagMultQuad LagMultQuad() { return DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(imortar_,"LAGMULT_QUAD") ; }

  /*!
  \brief Return number of Gauss points for this instance

  */
  int nGP() {return ngp_;}

  /*!
  \brief Return coordinates of a specific GP in 1D/2D CElement

  */
  double Coordinate(int& gp, int dir) {return coords_(gp,dir);}

  /*!
  \brief Return weight of a specific GP in 1D/2D CElement

  */
  double Weight(int& gp) {return weights_[gp];}

  //@}

  //! @name 2D and 3D integration methods

  /*!
  \brief Build all integrals and linearizations on a 1D or 2D
         slave element (i.e. only D, LinD possible).

  This method is outdated!!! We now obtain the D-matrix by
  integrating the mortar segments (just like the M-matrix entries).

  */
  virtual void IntegrateDerivSlave2D3D(
      MORTAR::MortarElement& sele, double* sxia, double* sxib,
      Teuchos::RCP<Epetra_SerialDenseMatrix> dseg);

  /*!
  \brief Perform mortar-integration without previous segmentation -- 2D

  */
  virtual void EleBased_Integration(
         Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
         Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
         MORTAR::MortarElement& sele,
         std::vector<MORTAR::MortarElement*> meles,
         bool *boundary_ele);

  /*!
  \brief Build all integrals and linearizations on a 1D slave /
         master overlap (i.e. D, M, g, LindD, LinM, Ling)

  */
  virtual void IntegrateDerivSegment2D(
      MORTAR::MortarElement& sele, double& sxia, double& sxib,
      MORTAR::MortarElement& mele, double& mxia, double& mxib,
      Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
      Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
      Teuchos::RCP<Epetra_SerialDenseVector> gseg);

  /*!
  \brief Integrate modification Mmod on a 1D slave / master overlap

  This modification is based on a paper by Puso/Wohlmuth, 2005.<br>
  It is necessary in the case of linear slave side elements
  and dual shape functions for the Lagrange multipliers, when
  the interface is curved (but only for mesh tying)!

  */
  Teuchos::RCP<Epetra_SerialDenseMatrix> IntegrateMmod2D(
       MORTAR::MortarElement& sele,
       double& sxia, double& sxib,
       MORTAR::MortarElement& mele,
       double& mxia, double& mxib);

  /*!
  \brief Build all integrals and linearizations on a 2D slave /
         master integration cell (i.e. M, g, LinM, Ling and possibly D, LinD)

  */
  virtual void IntegrateDerivCell3D(
      MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
      Teuchos::RCP<MORTAR::IntCell> cell,
      Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
      Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
      Teuchos::RCP<Epetra_SerialDenseVector> gseg);
    
  /*!
  \brief Build all integrals and linearizations without segmentation
         (i.e. M, g, LinM, Ling and possibly D, LinD)

  */
  virtual void IntegrateDerivCell3D_EleBased(
      MORTAR::MortarElement& sele, std::vector<MORTAR::MortarElement*> meles,
      Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
      Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
      Teuchos::RCP<Epetra_SerialDenseVector> gseg,
      bool *boundary_ele);

  /*!
  \brief Build all integrals and linearizations on a 2D slave /
         master integration cell (i.e. M, g, LinM, Ling and possibly D, LinD)
         for the auxiliary plane coupling case

  */
  virtual void IntegrateDerivCell3DAuxPlane(
     MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
     Teuchos::RCP<MORTAR::IntCell> cell, double* auxn,
     Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
     Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
     Teuchos::RCP<Epetra_SerialDenseVector> gseg);

  /*!
  \brief Build all integrals and linearizations on a 2D slave /
         master integration cell (i.e. M, g, LinM, Ling and possibly D, LinD)
         with quadratic interpolation

  */
  virtual void IntegrateDerivCell3DQuad(
     MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
     MORTAR::IntElement& sintele, MORTAR::IntElement& mintele,
     Teuchos::RCP<MORTAR::IntCell> cell,
     Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
     Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
     Teuchos::RCP<Epetra_SerialDenseVector> gseg);

  /*!
  \brief Build all integrals and linearizations on a 2D slave /
         master integration cell (i.e. M, g, LinM, Ling) for
         the auxiliary plane coupling case with quadratic interpolation

  */
  virtual void IntegrateDerivCell3DAuxPlaneQuad(
     MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
     MORTAR::IntElement& sintele, MORTAR::IntElement& mintele,
     Teuchos::RCP<MORTAR::IntCell> cell, double* auxn,
     Teuchos::RCP<Epetra_SerialDenseMatrix> dseg,
     Teuchos::RCP<Epetra_SerialDenseMatrix> mseg,
     Teuchos::RCP<Epetra_SerialDenseVector> gseg);
  
  //@}
  
  //! @name 2D and 3D assembly methods
  
  /*!
  \brief Assemble D contribution of current overlap into slave nodes

  */
  virtual bool AssembleD(const Epetra_Comm& comm,
                 MORTAR::MortarElement& sele,
                 Epetra_SerialDenseMatrix& D_seg);
  
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
  virtual bool AssembleM(const Epetra_Comm& comm,
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
  \brief Assemble M contribution of current slave-master pairing for EleBased_Integration

  */
  virtual bool AssembleM_EleBased(const Epetra_Comm& comm,
                 MORTAR::MortarElement& sele,
                 std::vector<MORTAR::MortarElement*> meles,
                 Epetra_SerialDenseMatrix& M_seg);

  /*!
  \brief Assemble Mmod contribution of current overlap into slave nodes

  */
  bool AssembleMmod(const Epetra_Comm& comm,
                    MORTAR::MortarElement& sele,
                    MORTAR::MortarElement& mele,
                    Epetra_SerialDenseMatrix& Mmod_seg);

  //@}
  
protected:

    // don't want = operator and cctor
    MortarIntegrator operator = (const MortarIntegrator& old);
    MortarIntegrator(const MortarIntegrator& old);

    Teuchos::ParameterList&       imortar_;    // containing contact input parameters
    int                           dim_;         // dimension of problem (2D or 3D)
    int                           ngp_;         // number of Gauss points
    Epetra_SerialDenseMatrix      coords_;      // Gauss point coordinates
    std::vector<double>           weights_;     // Gauss point weights

}; // class MortarIntegrator
}  // namespace MORTAR


#endif  // #ifndef MORTAR_INTEGRATOR_H

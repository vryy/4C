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
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef MORTAR_INTEGRATOR_H
#define MORTAR_INTEGRATOR_H

#include "Epetra_Comm.h"
#include "mortar_element.H"
#include "mortar_interface.H"
#include "mortar_coupling3d.H"
#include "mortar_interface.H"

using namespace std;

/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms

*/
namespace MORTAR
{

// forward declarations

/*!
\brief A class to perform Gaussian integration and assembly of Mortar
       matrices on the overlap of two MortarElements (1 Slave, 1 Master)
       in 1D (which is equivalent to a 2D coupling problem) and in 2D
       (which is equivalent to a 3D coupling problem)

\author popp (popp@lnm.mw.tum.de)
*/

class Integrator
{
public:

  /*!
  \brief Standard constructor

  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as overlaps are
  integrated in parallel by individual processes.<br>
  Note also that this constructor relies heavily on the
  DRT::UTILS::IntegrationPoints structs to get Gauss points
  and corresponding weights.
  
  Note: This version of the constructor creates an Integrator instance with undefined type of shape functions.
  As a result, no calls to functions relying on the evaluation of shape functions is allowed. 
  To be able to evaluate them, the Integrator have to be created with the alternative constructor (see below). 

  */
  explicit Integrator(DRT::Element::DiscretizationType eletype);

  /*!
  \brief Constructor  with shape function specification
  
  Constructs an instance of this class using a specific type of shape functions.<br>
  Note that this is \b not a collective call as overlaps are
  integrated in parallel by individual processes.<br>
  Note also that this constructor relies heavily on the
  DRT::UTILS::IntegrationPoints structs to get Gauss points
  and corresponding weights.
  
  */
  explicit Integrator(const MortarInterface::ShapeFcnType shapefcn,
                      DRT::Element::DiscretizationType eletype);

  /*!
  \brief Destructor

  */
  virtual ~Integrator() {}

  //! @name Access methods

  /*!
  \brief Get problem dimension

  Note that only 2D and 3D are possible here as this refers to the global
  problem dimension. On integration level this corresponds to 1D integration
  (dim_==2) and 2D integration (dim_==3) on the interface!

  */
  virtual inline const int& Dim() {return dim_;};

  /*!
  \brief Return number of Gauss points for this instance

  */
  inline int nGP() {return ngp_;}

  /*!
  \brief Return coordinates of a specific GP in 1D/2D CElement

  */
  inline double Coordinate(int& gp, int dir) {return coords_(gp,dir);}

  /*!
  \brief Return weight of a specific GP in 1D/2D CElement

  */
  inline double Weight(int& gp) {return weights_[gp];}

  //@}

  //! @name 2D and 3D integration methods

  /*!
  \brief Build all integrals and linearizations on a 1D or 2D
         slave element (i.e. D, LinD).

  This method is only called if the flag MORTARONELOOP is
  not set. Otherwise we obtain the D-matrix entries via a fundamental
  Mortar equation from the M-matrix entries.

  */
  void IntegrateDerivSlave2D3D(
      MORTAR::MortarElement& sele, double* sxia, double* sxib,
      RCP<Epetra_SerialDenseMatrix> dseg);

  /*!
  \brief Build all integrals and linearizations on a 1D slave /
         master overlap (i.e. M, g, LinM, Ling and possibly D, LinD)

  */
  void IntegrateDerivSegment2D(
      MORTAR::MortarElement& sele, double& sxia, double& sxib,
      MORTAR::MortarElement& mele, double& mxia, double& mxib,
      RCP<Epetra_SerialDenseMatrix> dseg,
      RCP<Epetra_SerialDenseMatrix> mseg,
      RCP<Epetra_SerialDenseVector> gseg);

  /*!
  \brief Integrate modification Mmod on a 1D slave / master overlap

  This modification is based on a paper by Puso/Wohlmuth, 2005.<br>
  It is necessary in the case of linear slave side elements
  and dual shape functions for the Lagrange multipliers, when
  the interface is curved (but only for mesh tying)!

  */
  RCP<Epetra_SerialDenseMatrix> IntegrateMmod2D(MORTAR::MortarElement& sele,
                                                double& sxia, double& sxib,
                                                MORTAR::MortarElement& mele,
                                                double& mxia, double& mxib);

  /*!
  \brief Build all integrals and linearizations on a 2D slave /
         master integration cell (i.e. M, g, LinM, Ling and possibly D, LinD)

  */
  void IntegrateDerivCell3D(
     MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
     RCP<MORTAR::Intcell> cell,
     RCP<Epetra_SerialDenseMatrix> dseg,
     RCP<Epetra_SerialDenseMatrix> mseg,
     RCP<Epetra_SerialDenseVector> gseg);

  /*!
  \brief Build all integrals and linearizations on a 2D slave /
         master integration cell (i.e. M, g, LinM, Ling and possibly D, LinD)
         for the auxiliary plane coupling case

  */
  void IntegrateDerivCell3DAuxPlane(
     MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
     RCP<MORTAR::Intcell> cell, double* auxn,
     RCP<Epetra_SerialDenseMatrix> dseg,
     RCP<Epetra_SerialDenseMatrix> mseg,
     RCP<Epetra_SerialDenseVector> gseg);

  /*!
  \brief Build all integrals and linearizations on a 2D slave /
         master integration cell (i.e. M, g, LinM, Ling) for
         the auxiliary plane coupling case with quadratic interpolation

  */
  void IntegrateDerivCell3DAuxPlaneQuad(
     MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
     MORTAR::IntElement& sintele, MORTAR::IntElement& mintele,
     RCP<MORTAR::Intcell> cell, double* auxn,
     RCP<Epetra_SerialDenseMatrix> dseg,
     RCP<Epetra_SerialDenseMatrix> mseg,
     RCP<Epetra_SerialDenseVector> gseg);
    
  /*!
  \brief Assemble D contribution of current overlap into slave nodes

  */
  bool AssembleD(const Epetra_Comm& comm,
                 MORTAR::MortarElement& sele,
                 Epetra_SerialDenseMatrix& D_seg);

  /*!
  \brief Assemble M contribution of current overlap into slave nodes

  */
  bool AssembleM(const Epetra_Comm& comm,
                 MORTAR::MortarElement& sele,
                 MORTAR::MortarElement& mele,
                 Epetra_SerialDenseMatrix& M_seg);

  /*!
  \brief Assemble Mmod contribution of current overlap into slave nodes

  */
  bool AssembleMmod(const Epetra_Comm& comm,
                    MORTAR::MortarElement& sele,
                    MORTAR::MortarElement& mele,
                    Epetra_SerialDenseMatrix& Mmod_seg);

  //@}

private:
  
  /*!
  \brief Initialize Gauss rule (points, weights) for this Integrator

  */
  void InitializeGP(DRT::Element::DiscretizationType eletype);
  
protected:

    // don't want = operator and cctor
    Integrator operator = (const Integrator& old);
    Integrator(const Integrator& old);

    const MortarInterface::ShapeFcnType shapefcn_;    // employed type of shape function set
    
    int                           dim_;         // dimension of problem (2D or 3D)
    int                           ngp_;         // number of Gauss points
    Epetra_SerialDenseMatrix      coords_;      // Gauss point coordinates
    vector<double>                weights_;     // Gauss point weights

}; // class Integrator
}  // namespace MORTAR


#endif  // #ifndef MORTAR_INTEGRATOR_H
#endif  // #ifdef CCADISCRET

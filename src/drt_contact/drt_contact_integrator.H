/*!----------------------------------------------------------------------
\file drt_contact_integrator.H

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
#ifndef DRT_CONTACT_INTEGRATOR_H
#define DRT_CONTACT_INTEGRATOR_H

#include "Epetra_Comm.h"
#include "drt_celement.H"
#include "drt_contact_interface.H"
#include "drt_contact_coupling.H"

using namespace std;

/*!
\brief CONTACT: namespace of the baci contact algorithms

*/
namespace CONTACT
{

// forward declarations

/*!
\brief A class to perform Gaussian integration and assembly of Mortar
       matrices on the overlap of two CElements (1 Slave, 1 Master)
       in 1D (which is equivalent to a 2D contact problem) and in 2D
       (which is equivalent to a 3D contact problem)

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
  
  */
  explicit Integrator(DRT::Element::DiscretizationType eletype);

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
  \brief Integrate Mortar matrix D on a 1D/2D slave element
  
  This method is only called if the flag CONTACTONEMORTARLOOP
  is not set. Otherwise we obtain the D-matrix entries via a
  fundamental Mortar equation from the M-matrix entries. Therefore,
  when AssembleM is called, Assemble D is done too in that case.
  
  */
  RCP<Epetra_SerialDenseMatrix> IntegrateD(CONTACT::CElement& sele,
                                           double* sxia, double* sxib);
  
  /*!
  \brief Compute directional derivative of D on a 1D/2D slave element
  
  This method is only called if the flag CONTACTONEMORTARLOOP
  is not set. Linearization of D is not yet implemented for the
  CONTACTONEMORTARLOOP case and will throw a dserror.
  
  */
  void DerivD(CONTACT::CElement& sele, double* sxia, double* sxib);
    
  /*!
  \brief Integrate Mortar matrix M on a 1D slave / master overlap
      
  */
  RCP<Epetra_SerialDenseMatrix> IntegrateM(CONTACT::CElement& sele,
                                           double& sxia, double& sxib,
                                           CONTACT::CElement& mele,
                                           double& mxia, double& mxib);
  
  /*!
  \brief Compute directional derivative of M on a 1D slave / master overlap
        
  */
  void DerivM(CONTACT::CElement& sele, double& sxia, double& sxib,
              CONTACT::CElement& mele, double& mxia, double& mxib);
  
  /*!
  \brief Compute directional derivative of segment end coordinates
         Xi on a 1D slave / master overlap
          
  */
  void DerivXiAB(CONTACT::CElement& sele, double& sxia, double& sxib,
                 CONTACT::CElement& mele, double& mxia, double& mxib,
                 vector<map<int,double> >& derivxi,
                 bool& startslave, bool& endslave);
  
  /*!
  \brief Compute directional derivative of master Gauss point
         coordinates XiGP on a 1D slave / master overlap 
            
  */
  void DerivXiGP(CONTACT::CElement& sele, CONTACT::CElement& mele,
                 double& sxigp, double& mxigp, const map<int,double>& derivsxi,
                 map<int,double>& derivmxi);
    
  /*!
  \brief Integrate modification Mmod on a 1D slave / master overlap
   
  This modification is based on a paper by Puso/Wohlmuth, 2005.<br>
  It is necessary in the case of linear slave side elements
  and dual shape functions for the Lagrange multipliers, when
  the interface is curved (but only for mesh tying)!
  
  */
  RCP<Epetra_SerialDenseMatrix> IntegrateMmod(CONTACT::CElement& sele,
                                              double& sxia, double& sxib,
                                              CONTACT::CElement& mele,
                                              double& mxia, double& mxib);
  
  /*!
  \brief Integrate weighted gap g~ on a 1D slave / master overlap
          
  */
  RCP<Epetra_SerialDenseVector> IntegrateG(CONTACT::CElement& sele,
                                           double& sxia, double& sxib,
                                           CONTACT::CElement& mele,
                                           double& mxia, double& mxib);
  
  /*!
  \brief Compute directional derivative of weighted gap g~ on
         a 1D slave / master overlap
        
  */
  void DerivG(CONTACT::CElement& sele, double& sxia, double& sxib,
              CONTACT::CElement& mele, double& mxia, double& mxib);
    
  /*!
  \brief Integrate Mortar matrix M on a 2D slave / master integration cell
        
  */
  RCP<Epetra_SerialDenseMatrix> IntegrateM3D(CONTACT::CElement& sele,
                                             CONTACT::CElement& mele,
                                             RCP<CONTACT::Intcell> cell);
  
  // THIS IS A FINITE DIFFERENCE VERSION ONLY!!!
  RCP<Epetra_SerialDenseMatrix> IntegrateM3D(CONTACT::CElement& sele,
                                             CONTACT::CElement& mele,
                                             RCP<CONTACT::Intcell> cell,
                                             vector<vector<double> >& testgps,
                                             vector<vector<double> >& testgpm,
                                             vector<vector<double> >& testjs,
                                             vector<vector<double> >& testji);
  
  /*!
  \brief Compute directional deriv. of M on a 2D slave / master integration cell
        
  */
  void DerivM3D(CONTACT::CElement& sele, CONTACT::CElement& mele,
                RCP<CONTACT::Intcell> cell, bool printderiv=false);
  
  /*!
  \brief Compute directional derivative of master Gauss point
         coordinates XiGP on a 2D slave / master integration cell 
            
  */
  void DerivXiGP3D(CONTACT::CElement& sele, CONTACT::CElement& mele,
                   double* sxigp, double* mxigp, const vector<map<int,double> >& derivsxi,
                   vector<map<int,double> >& derivmxi, double& alpha);
    
  /*!
  \brief Integrate Mortar matrix M on a 2D slave / master integration cell
         NOTE: This version works with the Auxiliary Plane of 3D Coupling!
        
  */
  RCP<Epetra_SerialDenseMatrix> IntegrateMAuxPlane3D(CONTACT::CElement& sele,
                                                     CONTACT::CElement& mele,
                                                     RCP<CONTACT::Intcell> cell,
                                                     double* auxn);
  
  /*!
  \brief Integrate weighted gap g~ on a 2D slave / master integration cell
        
  */
  RCP<Epetra_SerialDenseVector> IntegrateG3D(CONTACT::CElement& sele,
                                             CONTACT::CElement& mele,
                                             RCP<CONTACT::Intcell> cell);
  
  /*!
  \brief Compute directional deriv. of weighted gap g~ on a 2D slave /
         master integration cell
        
  */
  void DerivG3D(CONTACT::CElement& sele, CONTACT::CElement& mele,
                RCP<CONTACT::Intcell> cell);
    
    
  /*!
  \brief Integrate weighted gap g~ on a 2D slave / master integration cell
         NOTE: This version works with the Auxiliary Plane of 3D Coupling!
        
  */
  RCP<Epetra_SerialDenseVector> IntegrateGAuxPlane3D(CONTACT::CElement& sele,
                                                     CONTACT::CElement& mele,
                                                     RCP<CONTACT::Intcell> cell,
                                                     double* auxn);
    
  /*!
  \brief Assemble D contribution of current overlap into slave nodes
  
  This method is only called if the flag CONTACTONEMORTARLOOP
  is not set. Otherwise, when AssembleM is called, Assemble D
  is done there too.
  
  */
  bool AssembleD(const Epetra_Comm& comm,
                 CONTACT::CElement& sele,
                 Epetra_SerialDenseMatrix& D_seg);
  
  /*!
  \brief Assemble M contribution of current overlap into slave nodes
    
  */
  bool AssembleM(const Epetra_Comm& comm,
                 CONTACT::CElement& sele,
                 CONTACT::CElement& mele,
                 Epetra_SerialDenseMatrix& M_seg);
  
  /*!
  \brief Assemble Mmod contribution of current overlap into slave nodes
      
  */
  bool AssembleMmod(const Epetra_Comm& comm,
                    CONTACT::CElement& sele,
                    CONTACT::CElement& mele,
                    Epetra_SerialDenseMatrix& Mmod_seg);
  
  /*!
  \brief Assemble g~ contribution of current overlap into slave nodes
      
  */
  bool AssembleG(const Epetra_Comm& comm,
                 CONTACT::CElement& sele,
                 Epetra_SerialDenseVector& g_seg);
       
  //@}
      
protected:
    
    // don't want = operator and cctor
    Integrator operator = (const Integrator& old);
    Integrator(const Integrator& old);
  
    int                        dim_;         // dimension of problem (2D or 3D)
    int                        ngp_;         // number of Gauss points
    Epetra_SerialDenseMatrix   coords_;      // Gauss point coordinates
    vector<double>             weights_;     // Gauss point weights

}; // class Integrator
}  // namespace CONTACT


#endif  // #ifndef DRT_CONTACT_INTEGRATOR_H
#endif  // #ifdef CCADISCRET

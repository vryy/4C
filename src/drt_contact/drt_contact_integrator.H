/*!----------------------------------------------------------------------
\file drt_contact_integrator.H

\class CONTACT::Integrator

\brief A class to perform intgerations of Mortar matrices on the overlap
        of two CElements in 1D and 2D

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

#include "drt_celement.H"
#include "drt_contact_interface.H"

using namespace std;

/*!
\brief CONTACT: namespace of the baci contact algorithms

*/
namespace CONTACT
{

// forward declarations

/*!
\class Integrator

\brief A class to perform Gaussian integration and assembly of Mortar
       matrices on the overlap of two CElements (1 Slave, 1 Master)
       in 1D (which is equivalent to a 2D contact problem) and in 2D
       (which is equivalent to a 3D contact problem)


\author popp (popp@lnm.mw.tum.de)
*/

class Integrator
{
public:
  
  //! @name Constructors and destructors
  
  /*!
  \brief Standard constructor
  
  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as overlaps are
  integrated in parallel by individual processes.<br>
  Note also that this constructor relies heavily on the
  DRT::UTILS::IntegrationPoints structs to get Gauss points
  and corresponding weights.
  
  \param twoD (in): true if problem is 2D, false if 3D
  */
  explicit Integrator(int ngp, bool oned);

  /*!
  \brief Destructor

  */
  virtual ~Integrator() {}
  
  //@}
  
  //! @name Access methods
  
  /*!
  \brief Return whether this instance was created for 1D or 2D integrations
  
  */
  inline bool IsOneDimensional() {return oned_;}
  
  /*!
  \brief Return number of Gauss points for this instance
    
  */
  inline int nGP() {return ngp_;}
  
  /*!
  \brief Return coordinate of a specific GP in 1D CElement
    
  */
  inline double Coordinate(int gp) {return coords_[gp];}
  
  /*!
  \brief Return weight of a specific GP in 1D/2D CElement
      
  */
  inline double Weight(int gp) {return weights_[gp];}
  
  //@}
  
  //! @name 2D and 3D integration methods
  
  /*!
  \brief Integrate Mortar matrix D on a 1D slave element overlap
    
  */
  RCP<Epetra_SerialDenseMatrix> IntegrateD(CONTACT::CElement& sele,
                                           double sxia, double sxib);
  
  /*!
  \brief Integrate Mortar matrix M on a 1D slave / master overlap
      
  */
  RCP<Epetra_SerialDenseMatrix> IntegrateM(CONTACT::CElement& sele,
                                           double sxia, double sxib,
                                           CONTACT::CElement& mele,
                                           double mxia, double mxib);
  
  /*!
  \brief Integrate modification Mmod on a 1D slave / master overlap
   
   This modification is based on a paper by Puso/Wohlmuth, 2005.<br>
   It is necessary in the case of linear slave side elements
   and dual shape functions for the Lagrange multipliers, when
   the interface is curved!

  */
  RCP<Epetra_SerialDenseMatrix> IntegrateMmod(CONTACT::CElement& sele,
                                              double sxia, double sxib,
                                              CONTACT::CElement& mele,
                                              double mxia, double mxib);
  
  /*!
  \brief Integrate weighted gap g~ on a 1D slave / master overlap
          
  */
  RCP<Epetra_SerialDenseVector> IntegrateG(CONTACT::CElement& sele,
                                           double sxia, double sxib,
                                           CONTACT::CElement& mele,
                                           double mxia, double mxib);
  
  
  /*!
  \brief Assemble D contribution of current overlap into slave nodes
  
  */
  bool AssembleD(CONTACT::Interface& inter, CONTACT::CElement& sele,
                 Epetra_SerialDenseMatrix& D_seg);
  
  /*!
  \brief Assemble M contribution of current overlap into slave nodes
    
  */
  bool AssembleM(CONTACT::Interface& inter, CONTACT::CElement& sele,
                 CONTACT::CElement& mele,
                 Epetra_SerialDenseMatrix& M_seg);
  
  /*!
  \brief Assemble Mmod contribution of current overlap into slave nodes
      
  */
  bool AssembleMmod(CONTACT::Interface& inter, CONTACT::CElement& sele,
                    CONTACT::CElement& mele,
                    Epetra_SerialDenseMatrix& Mmod_seg);
  
  /*!
  \brief Assemble g~ contribution of current overlap into slave nodes
      
  */
  bool AssembleG(CONTACT::Interface& inter, CONTACT::CElement& sele,
                 Epetra_SerialDenseVector& g_seg);
    
    
  //@}
      
protected:
    
    // don't want = operator and cctor
    Integrator operator = (const Integrator& old);
    Integrator(const Integrator& old);
  
    bool             oned_;           // intgeration dimension, true if 1D
    int              ngp_;            // number of Gauss points
    vector<double>   coords_;         // Gauss point coordinates
    vector<double>   weights_;        // Gauss point weights

}; // class Integrator
}  // namespace CONTACT


#endif  // #ifndef DRT_CONTACT_INTEGRATOR_H
#endif  // #ifdef CCADISCRET

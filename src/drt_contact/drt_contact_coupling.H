/*!----------------------------------------------------------------------
\file drt_contact_coupling.H

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
<\pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef DRT_CONTACT_COUPLING_H
#define DRT_CONTACT_COUPLING_H

#include "Epetra_Comm.h"
#include "drt_celement.H"
#include "../drt_lib/drt_discret.H"

using namespace std;

/*!
\brief CONTACT: namespace of the baci contact algorithms

*/
namespace CONTACT
{

// forward declarations

/*!
\brief A class representing the framework for mortar coupling of ONE
       slave element and the master surface of a contact interface in
       2D and 3D. Concretely, this class controls projection, overlap
       detection and finally integration of the mortar coupling matrix M.

\author popp (popp@lnm.mw.tum.de)
*/

class Coupling
{
public:
  
  /*!
  \brief Standard constructor
  
  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.
  
  */
  explicit Coupling(DRT::Discretization& idiscret,
                    CONTACT::CElement& sele, int dim,
                    Epetra_SerialDenseMatrix& csegs);

  /*!
  \brief Destructor

  */
  virtual ~Coupling() {}
  
  //! @name Access methods
  
  /*!
  \brief Get interface discretization

  */
  virtual inline DRT::Discretization& Discret() const { return idiscret_; }
  
  /*!
  \brief Get communicator

  */
  virtual const Epetra_Comm& Comm() const { return idiscret_.Comm();}
    
  /*!
  \brief Get coupling slave element

  */
  virtual inline CONTACT::CElement& SlaveElement() const { return sele_; }
    
  /*!
  \brief Get problem dimension
  
  Note that only 2D and 3D are possible here as this refers to the global
  problem dimension. On interface level this corresponds to 1D interfaces
  (dim_==2) and 2D interfaces (dim_==3)!
      
  */
  virtual inline const int& Dim() {return dim_;};
  
  /*!
  \brief Return matrix of contact segment coordinates
  
  */
  virtual Epetra_SerialDenseMatrix& CSegs() {return contactsegs_;};
    
  //@}
    
  //! @name Evlauation methods
  
  /*!
  \brief Projection of slave / master pair (2D)
  
  This method projects the nodes of the slave CElement sele_ onto
  the master CElement mele and vice versa. The parameter hasproj
  stores a boolean variable for each of the 4 end nodes, indicating
  whether a feasible projection was found or not. The local element
  coordinates of the 4 projection points are stored in xiproj.
  
  */
  virtual bool Project2D(CONTACT::CElement& mele,
                         vector<bool>& hasproj,
                         vector<double>& xiproj);
  
  /*!
  \brief Detect overlap of slave / master pair (2D)
  
  This method evaluates the overlap of the current CElement pair
  sele_ / mele based on the projection status of the 4 end nodes
  (hasproj) and the coordinates of the projection points (xiproj).
  According to the detected overlap case, the integration limits
  are determined and written into xiproj and the overlap status
  is returned by a boolean variable.
  
  */
  virtual bool DetectOverlap2D(CONTACT::CElement& mele,
                               vector<bool>& hasproj,
                               vector<double>& xiproj);
  
  /*!
  \brief Integrate overlap of slave / master pair (2D)
  
  This method integrates the overlap of the current CElement
  pair sele_ / mele based on the integration limits (xiproj). The
  integration always includes the Mortar matrix M and the gap g.
  In some cases it also includes the Mortar matrix D (if the
  flag CONTACTONEMORTARLOOP is chosen) and the Mortar matrix
  Mmod (if the respective conditions are fulfilled).
        
  */
  virtual bool IntegrateOverlap2D(CONTACT::CElement& mele,
                                  vector<double>& xiproj);
  
  //@}
    
      
protected:
    
    // don't want = operator and cctor
    Coupling operator = (const Coupling& old);
    Coupling(const Coupling& old);
  
    
    DRT::Discretization&       idiscret_;    // discretization of the contact interface
    CONTACT::CElement&         sele_;        // slave element to perform coupling for
    int                        dim_;         // dimension of problem (2D or 3D)
    Epetra_SerialDenseMatrix&  contactsegs_; // contains coordinates of contact segments (gmsh)
  

}; // class Coupling
}  // namespace CONTACT


#endif  // #ifndef DRT_CONTACT_COUPLING_H
#endif  // #ifdef CCADISCRET

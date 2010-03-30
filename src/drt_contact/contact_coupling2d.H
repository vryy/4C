/*!----------------------------------------------------------------------
\file contact_coupling2d.H

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
#ifndef CONTACT_COUPLING2D_H
#define CONTACT_COUPLING2D_H

#include "../drt_mortar/mortar_coupling2d.H"

using namespace std;

/*!
\brief CONTACT: namespace of the BACI contact and meshtying algorithms

*/
namespace CONTACT
{

/*!
\brief A class representing the framework for mortar coupling of ONE
       slave element and ONE master element of a contact interface in
       2D. This is a derived class from MORTAR::Coupling2d which does
       the contact-specific stuff for 2d mortar coupling.

\author popp (popp@lnm.mw.tum.de)
*/

class CoCoupling2d : public MORTAR::Coupling2d
{
public:

  /*!
  \brief Standard constructor

  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.
  
  Note: This version of the constructor creates an CoCoupling2d instance with undefined type of shape functions.
  As a result, no calls to functions relying on the evaluation of shape functions is allowed. 
  To be able to evaluate them, the CoCoupling2d have to be created with the alternative constructor (see below).

  */
  CoCoupling2d(DRT::Discretization& idiscret, int dim,
               MORTAR::MortarElement& sele, MORTAR::MortarElement& mele);
  
  /*!
  \brief Constructor with shape function specification
  
  Constructs an instance of this class and enables custom shape function types.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.
  
  */  
  CoCoupling2d(const INPAR::MORTAR::ShapeFcn shapefcn, DRT::Discretization& idiscret, int dim,
               MORTAR::MortarElement& sele, MORTAR::MortarElement& mele);
  
  /*!
  \brief Destructor

  */
  virtual ~CoCoupling2d() {}

  //! @name Evlauation methods

  /*!
  \brief Integrate overlap of slave / master pair (2D)

  Derived version! Most importantly, in this derived version
  a CONTACT::CoIntegrator instance is created, which also
  does integration of the mortar quantity linearizations
  
  This method integrates the overlap of the current MortarElement
  pair sele_ / mele_ based on the integration limits (xiproj). The
  integration always includes the Mortar matrix M and the gap g.
  In some cases it also includes the Mortar matrix D (if the
  flag MORTARONELOOP is chosen) and the Mortar matrix
  Mmod (if the respective conditions are fulfilled). 

  */
  bool IntegrateOverlap(vector<double>& xiproj);

  //@}


protected:

    // don't want = operator and cctor
    CoCoupling2d operator = (const CoCoupling2d& old);
    CoCoupling2d(const CoCoupling2d& old);
    
    // no new variables as compared to base class

}; // class CoCoupling2d
}  // namespace CONTACT


#endif  // #ifndef CONTACT_COUPLING2D_H
#endif  // #ifdef CCADISCRET

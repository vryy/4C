/*!----------------------------------------------------------------------
\file mortar_coupling2d.H

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
#ifndef MORTAR_COUPLING2D_H
#define MORTAR_COUPLING2D_H

#include "Epetra_Comm.h"
#include "mortar_element.H"
#include "mortar_interface.H"
#include "../drt_lib/drt_discret.H"

using namespace std;

/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms

*/
namespace MORTAR
{

/*!
\brief A class representing the framework for mortar coupling of ONE
       slave element and ONE master element of a mortar interface in
       2D. Concretely, this class controls projection, overlap detection
       and finally integration of the mortar coupling matrices D and M
       and possibly of the weighted gap vector g~.

\author popp (popp@lnm.mw.tum.de)
*/

class Coupling2d{

public:

  /*!
  \brief Standard constructor

  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.
  
  Note: This version of the constructor creates an Coupling2d instance with undefined type of shape functions.
  As a result, no calls to functions relying on the evaluation of shape functions is allowed. 
  To be able to evaluate them, the Coupling2d have to be created with the alternative constructor (see below).

  */
  explicit Coupling2d(DRT::Discretization& idiscret, int dim,
                      MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
                      bool nonlinear = false);
  
  /*!
  \brief Constructor with shape function specification
  
  Constructs an instance of this class and enables custom shape function types.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.
  
  */  
  explicit Coupling2d(const MortarInterface::ShapeFcnType shapefcn, DRT::Discretization& idiscret, int dim,
                      MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
                      bool nonlinear = false);

  /*!
  \brief Destructor

  */
  virtual ~Coupling2d() {}

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
  virtual inline MORTAR::MortarElement& SlaveElement() const { return sele_; }

  /*!
  \brief Get coupling master element

  */
  virtual inline MORTAR::MortarElement& MasterElement() const { return mele_; }

  /*!
  \brief Get problem dimension (here: 2D)

  */
  virtual inline const int& Dim() {return dim_;};

  //@}

  //! @name Evlauation methods

  /*!
  \brief Projection of slave / master pair (2D)

  This method projects the nodes of the slave MortarElement sele_ onto
  the master MortarElement mele_ and vice versa. The parameter hasproj
  stores a boolean variable for each of the 4 end nodes, indicating
  whether a feasible projection was found or not. The local element
  coordinates of the 4 projection points are stored in xiproj.

  */
  virtual bool Project(vector<bool>& hasproj,
                       vector<double>& xiproj);

  /*!
  \brief Detect overlap of slave / master pair (2D)

  This method evaluates the overlap of the current MortarElement pair
  sele_ / mele_ based on the projection status of the 4 end nodes
  (hasproj) and the coordinates of the projection points (xiproj).
  According to the detected overlap case, the integration limits
  are determined and written into xiproj and the overlap status
  is returned by a boolean variable.

  */
  virtual bool DetectOverlap(vector<bool>& hasproj,
                             vector<double>& xiproj);

  /*!
  \brief Integrate overlap of slave / master pair (2D)

  This method integrates the overlap of the current MortarElement
  pair sele_ / mele_ based on the integration limits (xiproj). The
  integration always includes the Mortar matrix M and the gap g.
  In some cases it also includes the Mortar matrix D (if the
  flag MORTARONELOOP is chosen) and the Mortar matrix
  Mmod (if the respective conditions are fulfilled).

  */
  virtual bool IntegrateOverlap(vector<double>& xiproj);

  //@}


protected:

    // don't want = operator and cctor
    Coupling2d operator = (const Coupling2d& old);
    Coupling2d(const Coupling2d& old);
    
    const MortarInterface::ShapeFcnType shapefcn_;  // employed type of shape function set
    bool                                nonlinear_; // nonlinear version (incl. linearizations)

    DRT::Discretization&       idiscret_;    // discretization of the contact interface
    int                        dim_;         // problem dimension (here: 2D)
    MORTAR::MortarElement&     sele_;        // slave element to perform coupling for
    MORTAR::MortarElement&     mele_;        // master element to perform coupling for

}; // class Coupling2d
}  // namespace MORTAR


#endif  // #ifndef MORTAR_COUPLING2D_H
#endif  // #ifdef CCADISCRET

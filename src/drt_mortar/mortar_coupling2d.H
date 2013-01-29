/*!----------------------------------------------------------------------
\file mortar_coupling2d.H

\brief Classes for mortar coupling in 2D.

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
#ifndef MORTAR_COUPLING2D_H
#define MORTAR_COUPLING2D_H

#include <vector>
#include "Epetra_Comm.h"
#include "../drt_inpar/inpar_mortar.H"


// forward declarations
namespace DRT
{
  class Discretization;
}

/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms

*/
namespace MORTAR
{

// forward declarations
class MortarElement;

/*!
\brief A class representing the framework for mortar coupling of ONE
       slave element and ONE master element of a mortar interface in
       2D. Concretely, this class controls projection, overlap detection
       and finally integration of the mortar coupling matrices D and M
       and possibly of the weighted gap vector g~.

\author popp (popp@lnm.mw.tum.de)
*/

class Coupling2d
{
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
  Coupling2d(DRT::Discretization& idiscret,
             int dim, bool quad, INPAR::MORTAR::LagMultQuad lmtype,
             MORTAR::MortarElement& sele, MORTAR::MortarElement& mele);
  
  /*!
  \brief Constructor with shape function specification
  
  Constructs an instance of this class and enables custom shape function types.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.
  
  */  
  Coupling2d(const INPAR::MORTAR::ShapeFcn shapefcn, DRT::Discretization& idiscret,
             int dim, bool quad, INPAR::MORTAR::LagMultQuad lmtype,
             MORTAR::MortarElement& sele, MORTAR::MortarElement& mele);

  /*!
  \brief Destructor

  */
  virtual ~Coupling2d() {}

  //! @name Access methods

  /*!
  \brief Get interface discretization

  */
  virtual DRT::Discretization& Discret() const { return idiscret_; }

  /*!
  \brief Get communicator

  */
  virtual const Epetra_Comm& Comm() const;

  /*!
  \brief Get coupling slave element

  */
  virtual MORTAR::MortarElement& SlaveElement() const { return sele_; }

  /*!
  \brief Get coupling master element

  */
  virtual MORTAR::MortarElement& MasterElement() const { return mele_; }

  /*!
  \brief Get problem dimension (here: 2D)

  */
  virtual const int& Dim() {return dim_;};

  /*!
  \brief Get coupling / FE ansatz type (true = quadratic)

  */
  virtual const bool& Quad() {return quad_;};

  /*!
  \brief Return the LM interpolation / testing type for quadratic FE

  */
  INPAR::MORTAR::LagMultQuad LagMultQuad() { return lmtype_; }

  /*!
  \brief Get projection status of the four end nodes

  */
  virtual const std::vector<bool>& HasProj() {return hasproj_;}

  /*!
  \brief Get overlap regions in parameter spaces

  */
  virtual const std::vector<double>& XiProj() {return xiproj_;}

  /*!
  \brief Get overlap status

  */
  virtual const bool& Overlap() {return overlap_;};

  //@}

  //! @name Evaluation methods
  
  /*!
  \brief Checks roughly whether the two elements are near

  This methods checks the orientation of slave and master element.
  Projection and overlap detection only make sense if the scalar
  product of the two center normals is negative (i.e. if the two
  elements form an angle smaller than 90°).

  */
  virtual bool RoughCheckOrient();

  /*!
  \brief Projection of slave / master pair

  This method projects the nodes of the slave MortarElement sele_ onto
  the master MortarElement mele_ and vice versa. The variable hasproj
  stores a boolean variable for each of the 4 end nodes, indicating
  whether a feasible projection was found or not. The local element
  coordinates of the 4 projection points are stored in xiproj.

  */
  virtual bool Project();

  /*!
  \brief Detect overlap of slave / master pair

  This method evaluates the overlap of the current MortarElement pair
  sele_ / mele_ based on the projection status of the 4 end nodes
  (hasproj) and the coordinates of the projection points (xiproj).
  According to the detected overlap case, the integration limits
  are determined and written into xiproj and the overlap status
  is stored in the boolean variable overlap.

  */
  virtual bool DetectOverlap();

  /*!
  \brief Integrate overlap of slave / master pair

  This method integrates the overlap of the current MortarElement
  pair sele_ / mele_ based on the integration limits (xiproj). The
  integration always includes the Mortar matrices D/M and the gap g.

  */
  virtual bool IntegrateOverlap();

  //@}


protected:

    // don't want = operator and cctor
    Coupling2d operator = (const Coupling2d& old);
    Coupling2d(const Coupling2d& old);
    
    const INPAR::MORTAR::ShapeFcn shapefcn_;  // employed type of shape function set
    DRT::Discretization&          idiscret_;  // discretization of the contact interface
    int                           dim_;       // problem dimension (here: 2D)
    bool                          quad_;      // flag indicating coupling type (true = quadratic)
    INPAR::MORTAR::LagMultQuad    lmtype_;    // type of quadratic LM interpolation and testing
    MORTAR::MortarElement&        sele_;      // slave element to perform coupling for
    MORTAR::MortarElement&        mele_;      // master element to perform coupling for

    std::vector<bool>             hasproj_;   // projection status of the four end nodes
    std::vector<double>           xiproj_;    // overlap regions in parameter spaces
    bool                          overlap_;   // overlap status

}; // class Coupling2d


/*!
\brief A class representing the framework for mortar coupling of ONE
       slave element and SEVERAL master elements of a mortar interface in
       2D. This class simply stores and manages several Coupling2d objects.

\author popp (popp@lnm.mw.tum.de)
*/

class Coupling2dManager{

public:

  /*!
  \brief Standard constructor

  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.

  Note: This version of the constructor creates an Coupling2dManager instance with undefined type of
  shape functions. As a result, no calls to functions relying on the evaluation of shape functions is
  allowed. To be able to evaluate them, the Coupling2dManager have to be created with the alternative
  constructor (see below).

  */
  Coupling2dManager(DRT::Discretization& idiscret,
                    int dim, bool quad, INPAR::MORTAR::LagMultQuad lmtype,
                    INPAR::MORTAR::IntType inttype,
                    MORTAR::MortarElement* sele, std::vector<MORTAR::MortarElement*> mele);

  /*!
  \brief Constructor with shape function specification

  Constructs an instance of this class and enables custom shape function types.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.

  */
  Coupling2dManager(const INPAR::MORTAR::ShapeFcn shapefcn, DRT::Discretization& idiscret,
                    int dim, bool quad, INPAR::MORTAR::LagMultQuad lmtype,
                    INPAR::MORTAR::IntType inttype,
                    MORTAR::MortarElement* sele, std::vector<MORTAR::MortarElement*> mele);

  /*!
  \brief Destructor

  */
  virtual ~Coupling2dManager() {}

  /*!
  \brief Get coupling slave element

  */
  virtual MORTAR::MortarElement& SlaveElement() const { return *sele_; }

  /*!
  \brief Get one specific coupling master element

  */
  virtual MORTAR::MortarElement& MasterElement(int k) const { return *(mele_[k]); }

  /*!
  \brief Get all coupling master elements

  */
  virtual std::vector<MORTAR::MortarElement*> MasterElements() const { return mele_; }

  /*!
  \brief Get coupling pairs

  */
  virtual std::vector<Teuchos::RCP<MORTAR::Coupling2d> >& Coupling() { return coup_; }

  /*!
  \brief Get type of quadratic LM interpolation

  */
  virtual INPAR::MORTAR::LagMultQuad LagMultQuad() { return lmtype_; }

  /*!
  \brief Get integration type

  */
  virtual INPAR::MORTAR::IntType IntType() { return inttype_; }

  /*!
  \brief Evaluate coupling pairs

  */
  virtual bool EvaluateCoupling();

  /*!
  \brief Get coupling type

  */
  virtual const bool& Quad() {return quad_;};

  //@}


protected:

    // don't want = operator and cctor
    Coupling2dManager operator = (const Coupling2dManager& old);
    Coupling2dManager(const Coupling2dManager& old);

    const INPAR::MORTAR::ShapeFcn           shapefcn_;  // employed type of shape function set
    DRT::Discretization&                    idiscret_;  // discretization of the contact interface
    int                                     dim_;       // problem dimension (here: 2D)
    bool                                    quad_;      // flag indicating coupling type (true = quadratic)
    INPAR::MORTAR::LagMultQuad              lmtype_;    // type of quadratic LM interpolation and testing
    INPAR::MORTAR::IntType                  inttype_;   // type of numerical integration scheme
    MORTAR::MortarElement*                  sele_;      // slave element
    std::vector<MORTAR::MortarElement*>     mele_;      // master elements
    std::vector<Teuchos::RCP<Coupling2d> >  coup_;      // coupling pairs

}; // class Coupling2dManager
}  // namespace MORTAR

#endif  // #ifndef MORTAR_COUPLING2D_H

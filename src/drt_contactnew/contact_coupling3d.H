/*!----------------------------------------------------------------------
\file contact_coupling3d.H

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
#ifndef CONTACT_COUPLING3D_H
#define CONTACT_COUPLING3D_H

#include "../drt_mortar/mortar_coupling3d.H"

using namespace std;

/*!
\brief CONTACT: namespace of the BACI contact and meshtying algorithms

*/
namespace CONTACT
{

// forward declarations


/*!
\brief A class representing the framework for mortar coupling of ONE
       slave element and ONE master element of a mortar interface in
       3D. Concretely, this class controls projection, overlap detection
       and finally integration of the mortar coupling matrices D and M
       and possibly the weighted gap vector g~.
       Note that 3D Coupling can EITHER be done in physical space (this is
       the case when an auxiliary plane is used) or in the slave element
       parameter space (this is the case when everything is done directly
       on the slave surface without any auxiliary plane). The boolean class
       variable auxplane_ decides about this (true = auxiliary plane).
       
       This is a derived class from MORTAR::Coupling3d which does the
       contact-specific stuff for 3d mortar coupling.

\author popp (popp@lnm.mw.tum.de)
*/

class CoCoupling3d : public MORTAR::Coupling3d
{
public:

  /*!
  \brief Standard constructor

  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.

  Note: This version of the constructor creates an CCooupling3d instance with undefined type of shape functions.
  As a result, no calls to functions relying on the evaluation of shape functions is allowed. 
  To be able to evaluate them, the CoCoupling3d have to be created with the alternative constructor (see below).
  
  */
  explicit CoCoupling3d(DRT::Discretization& idiscret, int dim, bool quad,
                        bool auxplane,
                        MORTAR::MortarElement& sele, MORTAR::MortarElement& mele);

  /*!
  \brief Constructor with shape function specification
  
  Constructs an instance of this class and enables custom shape function types.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.
  
  */
  explicit CoCoupling3d(const MORTAR::MortarInterface::ShapeFcnType shapefcn, 
                        DRT::Discretization& idiscret, int dim, bool quad,
                        bool auxplane,
                        MORTAR::MortarElement& sele, MORTAR::MortarElement& mele);

  /*!
  \brief Destructor

  */
  virtual ~CoCoupling3d() {}

  //! @name Evlauation methods

  /*!
  \brief Build auxiliary plane from slave element (3D)

  Derived version, also doing normal linearization.
  
  This method builds an auxiliary plane based on the possibly
  warped slave element of this coupling class. This plane is
  defined by the slave normal at the slave element center.

  */
  bool AuxiliaryPlane();

  /*!
  \brief Integrate the integration cells (3D)

  Derived version! Most importantly, in this derived version
  a CONTACT::CoIntegrator instance is created, which also
  does integration of the mortar quantity linearizations

  This method creates an integrator object for the cell triangles,
  then projects the Gauss points back onto slave and master elements
  (1st case, aux. plane) or only back onto the master element (2nd case)
  in order to evaluate the respective shape function there. Then
  entries of the mortar matrix M and the weighted gap g are integrated
  and assembled into the slave element nodes.

  */
  virtual bool IntegrateCells();
  
  //@}

  //! @name Linearization methods
  
  /*!
  \brief Linearization of clip vertex coordinates (3D)

  This method computes and returns full linearizations of all
  clip polygon vertices. We distinguish three possible cases here,
  namely the vertex being a slave node, a projected master node in
  slave element parameter space or a line-clipping intersection in
  slave element paramater space. NOT implemented for AuxPlane case!

  */
  virtual bool VertexLinearization(vector<vector<map<int,double> > >& linvertex,
                           map<int,double>& projpar, bool printderiv=false);

  /*!
  \brief Linearization of clip vertex coordinates (3D)

  Sub-method of VertexLinearization for slave linearization.
  ONLY necessary for for AuxPlane case!

  */
  virtual bool SlaveVertexLinearization(vector<map<int,double> >& currlin,
                                int sid);

  /*!
  \brief Linearization of clip vertex coordinates (3D)

  Sub-method of VertexLinearization for master linearization.
  Two versions for (a) slave parameter space (b) AuxPlane case.

  */
  virtual bool MasterVertexLinearization(MORTAR::Vertex& currv,
                                 vector<map<int,double> >& currlin,
                                 int mid, double alpha);
  virtual bool MasterVertexLinearization(vector<map<int,double> >& currlin,
                                 int mid);

  /*!
  \brief Linearization of clip vertex coordinates (3D)

  Sub-method of VertexLinearization for lineclip linearization.
  Two versions for (a) slave parameter space (b) AuxPlane case.
  For the AuxPlane case version, we just combine the correct
  slave and master vertex linearizations here, which were already
  computed earlier in VertexLinearization3D!

  */
  virtual bool LineclipVertexLinearization(MORTAR::Vertex& currv,
                                   vector<map<int,double> >& currlin,
                                   MORTAR::Vertex* sv1, MORTAR::Vertex* sv2,
                                   MORTAR::Vertex* mv1, MORTAR::Vertex* mv2,
                                   map<int,double>& projpar);
  virtual bool LineclipVertexLinearization(MORTAR::Vertex& currv,
                                   vector<map<int,double> >& currlin,
                                   MORTAR::Vertex* sv1, MORTAR::Vertex* sv2,
                                   MORTAR::Vertex* mv1, MORTAR::Vertex* mv2,
                                   vector<vector<map<int,double> > >& linsnodes,
                                           vector<vector<map<int,double> > >& linmnodes);

  /*!
  \brief Linearization of clip vertex coordinates (3D)

  This method computes and returns the full linearization of
  the clip polygon center, which itself is obtained from the
  clip polygon vertices by centroid formulas. NOT implemented
  for AuxPlane case!

  */
  virtual bool CenterLinearization(const vector<vector<map<int,double> > >& linvertex,
                           vector<map<int,double> >& lincenter);
  
  //@}

protected:

    // don't want = operator and cctor
    CoCoupling3d operator = (const CoCoupling3d& old);
    CoCoupling3d(const CoCoupling3d& old);

    // no new variables as compared to base class
    
}; // class CoCoupling3d

/*!
\brief A class representing the framework for mortar coupling of ONE
       slave element and ONE master element of a mortar interface in
       3D. Concretely, this class controls projection, overlap
       detection and finally integration of the mortar coupling matrices
       D and M and possibly the weighted gap vector g~.
       
       This is a special derived class for 3D quadratic mortar coupling
       with the use of auxiliary planes. This approach is based on
       "Puso, M.A., Laursen, T.A., Solberg, J., A segment-to-segment
       mortar contact method for quadratic elements and large deformations,
       CMAME, 197, 2008, pp. 555-566". For this type of formulation, a
       quadratic MortarElement is split into several linear IntElements,
       on which the geometrical coupling is performed. Thus, we additionally
       hand in in two IntElements to CoCoupling3dQuad.
       
       This is a derived class from MORTAR::CoCoupling3d which does the
       contact-specific stuff for 3d quadratic mortar coupling.

\author popp (popp@lnm.mw.tum.de)
*/

class CoCoupling3dQuad : public CoCoupling3d{

public:

  /*!
  \brief Standard constructor

  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.

  Note: This version of the constructor creates an CoCoupling3dQuad instance with undefined type of shape functions.
  As a result, no calls to functions relying on the evaluation of shape functions is allowed. 
  To be able to evaluate them, the CoCoupling3dQuad have to be created with the alternative constructor (see below).
  
  */
  explicit CoCoupling3dQuad(DRT::Discretization& idiscret,
                            int dim, bool quad, bool auxplane,
                            MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
                            MORTAR::IntElement& sintele, MORTAR::IntElement& mintele);

  /*!
  \brief Constructor with shape function specification
  
  Constructs an instance of this class and enables custom shape function types.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.
  
  */
  explicit CoCoupling3dQuad(const MORTAR::MortarInterface::ShapeFcnType shapefcn,
                            DRT::Discretization& idiscret,
                            int dim, bool quad, bool auxplane,
                            MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
                            MORTAR::IntElement& sintele, MORTAR::IntElement& mintele);
  
  /*!
  \brief Destructor

  */
  virtual ~CoCoupling3dQuad() {}
  
  //! @name Access methods

  /*!
  \brief Get coupling slave integration element

  */
  MORTAR::IntElement& SlaveIntElement() const { return sintele_; }

  /*!
  \brief Get coupling master integration element

  */
  MORTAR::IntElement& MasterIntElement() const { return mintele_; }

  //@}
  
protected:

    // don't want = operator and cctor
    CoCoupling3dQuad operator = (const CoCoupling3dQuad& old);
    CoCoupling3dQuad(const CoCoupling3dQuad& old);

    MORTAR::IntElement&   sintele_;     // slave sub-integration element
    MORTAR::IntElement&   mintele_;     // slave sub-integration element

}; // class CoCoupling3dQuad
}  // namespace CONTACT


#endif  // #ifndef CONTACT_COUPLING3D_H
#endif  // #ifdef CCADISCRET

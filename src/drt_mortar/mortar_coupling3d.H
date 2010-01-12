/*!----------------------------------------------------------------------
\file mortar_coupling3d.H

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
#ifndef MORTAR_COUPLING3D_H
#define MORTAR_COUPLING3D_H

#include "Epetra_Comm.h"
#include "mortar_element.H"
#include "mortar_interface.H"
#include "mortar_coupling3d_classes.H"
#include "../drt_lib/drt_discret.H"

using namespace std;

/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms

*/
namespace MORTAR
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

\author popp (popp@lnm.mw.tum.de)
*/

class Coupling3d{

public:

  /*!
  \brief Standard constructor

  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.

  Note: This version of the constructor creates an Coupling3d instance with undefined type of shape functions.
  As a result, no calls to functions relying on the evaluation of shape functions is allowed. 
  To be able to evaluate them, the Coupling3d have to be created with the alternative constructor (see below).
  
  */
  explicit Coupling3d(DRT::Discretization& idiscret, int dim, bool quad,
                      bool auxplane,
                      MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
                      bool nonlinear = false);

  /*!
  \brief Constructor with shape function specification
  
  Constructs an instance of this class and enables custom shape function types.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.
  
  */
  explicit Coupling3d(const MortarInterface::ShapeFcnType shapefcn, 
                      DRT::Discretization& idiscret, int dim, bool quad,
                      bool auxplane,
                      MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
                      bool nonlinear = false);

  /*!
  \brief Destructor

  */
  virtual ~Coupling3d() {}

  //! @name Access methods

  /*!
  \brief Get interface discretization

  */
  virtual DRT::Discretization& Discret() const { return idiscret_; }

  /*!
  \brief Get communicator

  */
  virtual const Epetra_Comm& Comm() const { return idiscret_.Comm();}

  /*!
  \brief Get problem dimension (here: 3D)

  */
  virtual const int& Dim() {return dim_;};

  /*!
  \brief Get coupling / FE ansatz type (true = quadratic)

  */
  virtual const bool& Quad() {return quad_;};

  /*!
  \brief Get coupling slave element

  */
  virtual MORTAR::MortarElement& SlaveElement() const { return sele_; }

  /*!
  \brief Get coupling master element

  */
  virtual MORTAR::MortarElement& MasterElement() const { return mele_; }

  /*!
  \brief Get coupling slave integration element

  Note that (here) for linear ansatz functions in 3D, this is IDENTICAL to
  the SlaveElement() as no splitting of the MortarElements is performed.
  For the 3D quadratic case with the use of auxiliary planes, this
  methd is overloaded by the derived class MORTAR::Coupling3dQuad!

  */
  virtual MORTAR::MortarElement& SlaveIntElement() const { return sele_; }

  /*!
  \brief Get coupling master integration element

  Note that (here) for linear ansatz functions in 3D, this is IDENTICAL to
  the SlaveElement() as no splitting of the MortarElements is performed.
  For the 3D quadratic case with the use of auxiliary planes, this
  methd is overloaded by the derived class MORTAR::Coupling3dQuad!

  */
  virtual MORTAR::MortarElement& MasterIntElement() const{ return mele_; }

  /*!
  \brief Get 3D coupling strategy (true = auxiliary plane)

  Note that 3D Coupling can EITHER be done in physical space (this is
  the case when an auxiliary plane is used) or in the slave element
  parameter space (this is the case when everything is done directly
  on the slave surface without any auxiliary plane).

  */
  virtual const bool& CouplingInAuxPlane() {return auxplane_;};

  /*!
  \brief Return center of auxiliary plane

  */
  virtual double* Auxc() { return auxc_; }

  /*!
  \brief Return normal of auxiliary plane

  */
  virtual double* Auxn() { return auxn_; }

  /*!
  \brief Return length of Auxn() before normalization

  */
  virtual double& Lauxn() {return lauxn_; }

  /*!
  \brief Return vector of (projected) slave node vertex objects

  */
  virtual vector<Vertex>& SlaveVertices() {return svertices_; }

  /*!
  \brief Return vector of projected master node vertex objects

  */
  virtual vector<Vertex>& MasterVertices() {return mvertices_; }

  /*!
  \brief Return vector of clip polygon vertex objects

  */
  virtual vector<Vertex>& Clip() {return clip_; }

  /*!
  \brief Return vector of integration cells

  */
  virtual vector<RCP<Intcell> >& Cells() {return cells_; }

  /*!
  \brief Return the 'DerivAuxn' map (vector) of this coupling pair

  */
  virtual vector<map<int,double> >& GetDerivAuxn() { return derivauxn_; }

  //@}

  //! @name Evlauation methods

  /*!
  \brief Evaluate coupling (3D)

  */
  virtual bool EvaluateCoupling();

  /*!
  \brief Checks roughly whether the two elements are near (3D)

  This methods computes the distance of the two element centers.
  If they are not close, then coupling is stopped for the pair.

  */
  virtual bool RoughCheckCenters();
  
  /*!
  \brief Checks roughly whether the two elements are near (3D)

  This methods computes the distance of all master nodes to the
  slave element (auxiliary plane). If they are not close, then
  coupling is stopped for the pair.

  */
  virtual bool RoughCheckNodes();

  /*!
  \brief Build auxiliary plane from slave element (3D)

  This method builds an auxiliary plane based on the possibly
  warped slave element of this coupling class. This plane is
  defined by the slave normal at the slave element center.

  */
  virtual bool AuxiliaryPlane();

  /*!
  \brief Projection of slave element onto aux. plane (3D)

  This method projects the nodes of the given slave CElement
  onto the auxiliary plane derived before.

  */
  virtual bool ProjectSlave();

  /*!
  \brief Projection of master element onto aux. plane (3D)

  This method projects the nodes of the current master CElement
  onto the auxiliary plane derived from the slave CElement before.

  */
  virtual bool ProjectMaster();

  /*!
  \brief Clipping of slave and master element (3D)

  Note that PolygonClipping3D can EITHER be done in physical space (this
  is the case when an auxiliary plane is used for 3D coupling) or in the
  slave element parameter space (this is the case when 3D coupling is
  performed on the slave surface without any auxiliary plane).

  This method applies a polygon clipping algorithm to find the
  polygon defined by the projection of the slave and master element
  in the auxiliary plane (1st case) or defined by the slave element
  itself and the projection of the master element into the slave
  parameter space (2nd case). Of course, in the 2nd case the clipping
  plane has the normal (0,0,1). As input variables the clipping algorithm
  requires two sets of vertices, which define sele_ / mele_ respectively.
  The clipping algorithm is based on doubly linked lists as data structure
  (Greiner, G., Hormann, K. (1998), Efficient Clipping of Arbitrary
  Polygons, ACM Transactions on Graphics, vol. 17, pp. 71-83).
  Returns a vector of vertices defining the clipped polygon.

  \param poly1list (in): vector of vertex objects for input polygon 1 (slave)
  \param poly2list (in): vector of vertex objects for input polygon 2 (master)
  \param respoly (out): vector of vertex objects for result polygon
  \param tol (in): clipping tolerance for close vertices detection

  */
  virtual void PolygonClipping(vector<Vertex>& poly1list,
                               vector<Vertex>& poly2list,
                               vector<Vertex>& respoly,
                               double& tol);
  
  /*!
  \brief Clipping of slave and master element (3D) (NEW VERSION based on convex hull)

  Note that PolygonClipping3D can EITHER be done in physical space (this
  is the case when an auxiliary plane is used for 3D coupling) or in the
  slave element parameter space (this is the case when 3D coupling is
  performed on the slave surface without any auxiliary plane).

  This method applies a polygon clipping algorithm to find the
  polygon defined by the projection of the slave and master element
  in the auxiliary plane (1st case) or defined by the slave element
  itself and the projection of the master element into the slave
  parameter space (2nd case). Of course, in the 2nd case the clipping
  plane has the normal (0,0,1). As input variables the clipping algorithm
  requires two sets of vertices, which define sele_ / mele_ respectively.
  The clipping algorithm is based on finding the convex hull of three
  sets of vertices: sele_ vertices, mele_ vertices, line intersections,
  using an algorithm quite similar to the well-known "Graham Scan".
  This makes the algorithm FAR MORE ROBUST than the old approach above.
  Returns a vector of vertices defining the clipped polygon.

  In some cases, projection of the master element may lead to a non-convex
  input polygon poly2list, which leads to a dserror. Before the error is
  thrown, this method (similarly to RoughCheck()) again checks if the two
  elements to tbe clipped are close to each other at all and returns a
  boolean flag. If they are not close, then coupling is stopped for the pair.
  
  \param poly1list (in): vector of vertex objects for input polygon 1 (slave)
  \param poly2list (in): vector of vertex objects for input polygon 2 (master)
  \param respoly (out): vector of vertex objects for result polygon
  \param tol (in): clipping tolerance for close vertices detection

  */
  virtual bool PolygonClippingConvexHull(vector<Vertex>& poly1list,
                                         vector<Vertex>& poly2list,
                                         vector<Vertex>& respoly,
                                         double& tol);

  /*!
  \brief Check / set projection status of slave nodes (3D)

  This method checks for all slave nodes if they are part of the clip
  polygon (equal to any vertex). If so the HasProj status is set true!

  */
  virtual bool HasProjStatus();

  /*!
  \brief Triangulation of clip polygon (3D)

  This method first finds the center of the clip polygon by computing
  the centroid / geometric center (2nd case) or by simple averaging of
  vertex coordinates (arithmetic center, 1st case, aux. plane!). In a
  second step, linearizations of all clip polygon vertices and of the
  clip center are computed and stored into local maps. Then the clip
  polygon is triangulated into as many integration cells (tri3) as
  necessary. The method stores these cells in the variable cells_.

  */
  virtual bool Triangulation(map<int,double>& projpar);

  /*!
  \brief Integrate the integration cells (3D)

  This method creates an integrator object for the cell triangles,
  then projects the Gauss points back onto slave and master elements
  (1st case, aux. plane) or only back onto the master element (2nd case)
  in order to evaluate the respective shape function there. Then
  entries of the mortar matrix M and the weighted gap g are integrated
  and assembled into the slave element nodes.

  */
  virtual bool IntegrateCells();
  
  //@}

  //! @name Linearization methods (for mortar contact only!)
  
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
  virtual bool MasterVertexLinearization(Vertex& currv,
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
  virtual bool LineclipVertexLinearization(Vertex& currv,
                                           vector<map<int,double> >& currlin,
                                           Vertex* sv1, Vertex* sv2, Vertex* mv1, Vertex* mv2,
                                           map<int,double>& projpar);
  virtual bool LineclipVertexLinearization(Vertex& currv,
                                           vector<map<int,double> >& currlin,
                                           Vertex* sv1, Vertex* sv2, Vertex* mv1, Vertex* mv2,
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
    Coupling3d operator = (const Coupling3d& old);
    Coupling3d(const Coupling3d& old);

    const MortarInterface::ShapeFcnType shapefcn_;  // employed type of shape function set
    bool                                nonlinear_; // nonlinear version (incl. linearizations)

    DRT::Discretization&       idiscret_;    // discretization of the contact interface
    int                        dim_;         // problem dimension (here: 3D)
    bool                       quad_;        // flag indicating coupling type (true = quadratic)
    bool                       auxplane_;    // flag indicating coupling strategy (true = aux. plane)
    MORTAR::MortarElement&     sele_;        // slave element to perform coupling for
    MORTAR::MortarElement&     mele_;        // master element to perform coupling for

    double                     auxc_[3];     // center of auxiliary plane
    double                     auxn_[3];     // normal of auxiliary plane
    double                     lauxn_;       // length of interpolated Auxn() before normalization
    vector<Vertex>             svertices_;   // slave node vertex objects
    vector<Vertex>             mvertices_;   // master node vertex objects
    vector<Vertex>             clip_;        // clipped polygon vertex objects

    vector<RCP<Intcell> >      cells_;       // vector of integration cells
    vector<map<int,double> >   derivauxn_;   // derivatives of auxiliary plane normal

}; // class Coupling3d

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
       hand in in two IntElements to Coupling3dQuad.

\author popp (popp@lnm.mw.tum.de)
*/

class Coupling3dQuad : public Coupling3d{

public:

  /*!
  \brief Standard constructor

  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.

  Note: This version of the constructor creates an Coupling3dQuad instance with undefined type of shape functions.
  As a result, no calls to functions relying on the evaluation of shape functions is allowed. 
  To be able to evaluate them, the Coupling3dQuad have to be created with the alternative constructor (see below).
  
  */
  explicit Coupling3dQuad(DRT::Discretization& idiscret,
                          int dim, bool quad, bool auxplane,
                          MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
                          MORTAR::IntElement& sintele, MORTAR::IntElement& mintele,
                          bool nonlinear = false);

  /*!
  \brief Constructor with shape function specification
  
  Constructs an instance of this class and enables custom shape function types.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.
  
  */
  explicit Coupling3dQuad(const MortarInterface::ShapeFcnType shapefcn,
                          DRT::Discretization& idiscret,
                          int dim, bool quad, bool auxplane,
                          MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
                          MORTAR::IntElement& sintele, MORTAR::IntElement& mintele,
                          bool nonlinear = false);
  
  /*!
  \brief Destructor

  */
  virtual ~Coupling3dQuad() {}

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
    Coupling3dQuad operator = (const Coupling3dQuad& old);
    Coupling3dQuad(const Coupling3dQuad& old);

    MORTAR::IntElement&   sintele_;     // slave sub-integration element
    MORTAR::IntElement&   mintele_;     // slave sub-integration element

}; // class Coupling3dQuad
}  // namespace MORTAR


#endif  // #ifndef MORTAR_COUPLING3D_H
#endif  // #ifdef CCADISCRET

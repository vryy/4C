/*!----------------------------------------------------------------------
\file drt_contact_coupling3d_classes.H

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
#ifndef DRT_CONTACT_COUPLING3D_CLASSES_H
#define DRT_CONTACT_COUPLING3D_CLASSES_H

#include "Epetra_Comm.h"
#include "drt_celement.H"
#include "../drt_lib/drt_discret.H"
#include "drt_contact_interface.H"

using namespace std;

/*!
\brief CONTACT: namespace of the baci contact algorithms

*/
namespace CONTACT
{

// forward declarations

/*!
\brief A class representing a special contact element for 3D quadratic
       mortar coupling with the use auf auxiliary planes. This approach
       is based on "Puso, M.A., Laursen, T.A., Solberg, J., A segment-to-
       segment mortar contact method for quadratic elements and large
       deformations, CMAME, 197, 2008, pp. 555-566". For this type of
       quadratic formulation, a quadratic CElement is split into several
       linear IntElements, on which the geometrical coupling is performed.
       Note that this is a derived class from CONTACT::CElement, with the
       only difference that we explicitly set the node pointers here!

\author popp (popp@lnm.mw.tum.de)
*/

class IntElement : public CElement{

public:

  /*!
  \brief Standard constructor

  \param lid   (in): A locally unique element id
  \param id    (in): Parent element global id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  \param parshape (in): shape of parent element
  \param shape (in): shape of this element
  \param numnode (in): Number of nodes to this element
  \param nodeids (in): ids of nodes adjacent to this element
  \param nodes (in): pointers to nodes adjacent to this element
  \param isslave (in): flag indicating whether element is slave or master side
  */
  explicit IntElement(int lid, int id, ElementType etype, int owner,
                      const DRT::Element::DiscretizationType& parshape,
                      const DRT::Element::DiscretizationType& shape,
                      const int numnode,
                      const int* nodeids,
                      vector<DRT::Node*> nodes,
                      const bool isslave);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an IntElement

  */
  IntElement(const IntElement& old);

  /*!
  \brief Destructor

  */
  virtual ~IntElement() {}

  /*!
  \brief Return local id of this IntElement

  */
  inline int Lid() const { return lid_; }

  /*!
  \brief Get vector of ptrs to nodes

  The pointers to the nodes are usually built in
  DRT::Discretization::FillComplete. A standalone
  element that has not been added to a discretization
  (or the discretization has not been called FillComplete),
  such as is the case here, does not have pointers to nodes,
  unless we explicitly define them. Thus, we hand them in
  to the constructur here and overload the base class method
  Nodes() such that the user-provided nodes are used first!

  */
  virtual inline DRT::Node** Nodes()
  { if (quadnode_.size()) return &quadnode_[0]; else return NULL;}

  /*!
  \brief Get shape type of parent element

  */
  virtual inline DiscretizationType ParShape() const { return parshape_;}

  /*!
  \brief Affine map of IntElement coordinates to parent element

  */
  virtual bool MapToParent(const double* xi, double* parxi);

  /*!
  \brief Affine map of IntElement coordinate derivatives to parent element

  */
  virtual bool MapToParent(const vector<map<int,double> >& dxi,
                             vector<map<int,double> >& dparxi);

protected:

  int                  lid_;        // local IntElement id
  vector<DRT::Node*>   quadnode_;   // explicit node pointers
  DRT::Element::DiscretizationType   parshape_;   // shape of parent element

}; //class IntElement


/*!
\brief A class representing one Integration Cell after triangulation
       of the clip polygon of slave and master element from the Coupling
       class. This class provides some basic funcitonality a CElement would
       also provide (coords, shape functions and derivatives, Jacobian, ...).
       Note that an Intcell can EITHER live in physical space (this is the
       case when an auxiliary plane is used for 3D coupling) or in the slave
       element parameter space (this is the case when 3D coupling is performed
       on the slave surface without any auxiliary plane). Of course, in the
       2nd case the third coordinate of all intcell point is zero!

\author popp (popp@lnm.mw.tum.de)
*/

class Intcell{

public:

  /*!
  \brief Standard constructor

  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.

  */
  explicit Intcell(int id, int nvertices, Epetra_SerialDenseMatrix& coords,
                   double* auxn, const DRT::Element::DiscretizationType& shape,
                   bool auxplane,
                   vector<map<int,double> >& linv1,
                   vector<map<int,double> >& linv2,
                   vector<map<int,double> >& linv3,
                   vector<map<int,double> >& linauxn);

  /*!
  \brief Copy Constructor

  Makes a deep copy of an Intcell

  */
  Intcell(const Intcell& old);

  /*!
  \brief Destructor

  */
  virtual ~Intcell() {}

  //! @name Access methods

  /*!
  \brief Return ID this intcell

  */
  inline int Id() const { return id_; }

  /*!
  \brief Return number of vertices of this intcell

  */
  inline int NumVertices() const { return nvertices_; }

  /*!
  \brief Return current area

  */
  virtual inline double& Area() { return area_; }

  /*!
  \brief Return coordinates of intcell vertices

  */
  inline const Epetra_SerialDenseMatrix& Coords() { return coords_; }

  /*!
  \brief Return normal of auxiliary plane of this intcell

  */
  virtual inline double* Auxn() { return auxn_; }

  /*!
  \brief Get shape type of element

  */
  virtual DRT::Element::DiscretizationType Shape() const { return shape_;}

  /*!
  \brief Get 3D coupling strategy (true = auxiliary plane)

  Note that 3D Coupling can EITHER be done in physical space (this is
  the case when an auxiliary plane is used) or in the slave element
  parameter space (this is the case when everything is done directly
  on the slave surface without any auxiliary plane).

  */
  virtual inline const bool& CouplingInAuxPlane() {return auxplane_;};

  /*!
  \brief Return one of the three 'DerivVertex' maps (vectors) of this node

  These maps contain the directional derivatives of the intcell's vertex
  coordinates with respect to the slave and master displacements. A vector
  is used because the coordinates themselves are a vector (2 components for
  the 'coupling in slave parameter space' case). The respective vertex (1,2,3)
  is addressed by an int-variable and checked internally.

  */
  virtual inline vector<map<int,double> >& GetDerivVertex(int i)
  {
    if (i<0 || i>2) dserror("ERROR: Intcell has 3 vertex linearizations only!");
    return linvertex_[i];
  }

  /*!
  \brief Return the 'DerivAuxn' map (vector) of this intcell

  */
  virtual inline vector<map<int,double> >& GetDerivAuxn() { return linauxn_; }

  //@}

  //! @name Evaluation methods

  /*!
  \brief Interpolate global coordinates for given local intcell coordinates

  This method interpolates global coordinates for a given local intcell
  coordinate variable using the intcell vertex coordinates. For interpolation
  one can choose between shape functions or shape function derivatives!

  \param xi (in)        : local intcell coordinates
  \param inttype (in)   : set to 0 for shape function usage,
                          set to 1 for derivative xi usage
                          set to 2 for derivative eta usage
  \param globccord (out): interpolated global coordinates
  */
  virtual bool LocalToGlobal(const double* xi, double* globcoord, int inttype);

  /*!
  \brief Evaluate shape functions and derivatives

  */
  virtual bool EvaluateShape(const double* xi, LINALG::SerialDenseVector& val,
                             LINALG::SerialDenseMatrix& deriv);

  /*!
  \brief Evaluate Jacobian determinant for parameter space integration

  */
  virtual double Jacobian(double* xi);

  /*!
  \brief Compute Jacobian determinant derivative
         Note that this is a linearization with respect to the intcell
         vertices, which themselves have to be linearized later, of course!
         Two versions for (a) slave parameter space (b) AuxPlane case.
  */
  virtual void DerivJacobian(double* xi, vector<double>& derivjac);
  virtual void DerivJacobian(double* xi, map<int,double>& derivjac);

  //@}

protected:

  int                        id_;              // local ID of this intcell
  int                        nvertices_;       // number of vertices (always 3)
  double                     area_;            // integration cell area
  Epetra_SerialDenseMatrix   coords_;          // coords of cell vertices
  double                     auxn_[3];         // normal of auxiliary plane (3D)

  DRT::Element::DiscretizationType   shape_;     // shape of this element (always tri3)
  bool                               auxplane_;  // flag indicating coupling strategy (true = aux. plane)
  vector<vector<map<int,double> > >  linvertex_; // derivatives of the 3 vertices
  vector<map<int,double> >           linauxn_;   // derivatives of auxiliary plane normal

}; // class Intcell


/*!
\brief A class representing one Vertex during the polygon clipping of
       slave and master element from the Coupling class. Besides the
       vertex coordinates this class provides different pointers to
       build up doubly-linked list structures.
       Note that a Vertex can EITHER live in physical space (this is the
       case when an auxiliary plane is used for 3D coupling) or in the slave
       element parameter space (this is the case when 3D coupling is performed
       on the slave surface without any auxiliary plane). Of course, in the
       2nd case the third coordinate of all vertices is zero!

\author popp (popp@lnm.mw.tum.de)
*/

class Vertex{

public:

  //! @name Enums and Friends
  enum vType // vertex types recognized by Vertex
    {
      slave,        // slave node
      projmaster,   // projected master node
      lineclip      // clipping point of two lines
    };

  //@}

  /*!
  \brief Standard constructor

  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.

  */
  explicit Vertex(vector<double> coord, Vertex::vType type,
                  vector<int> nodeids, Vertex* next, Vertex* prev,
                  bool intersect, bool entryexit, Vertex* neighbor,
                  double alpha);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Vertex

  */
  Vertex(const Vertex& old);

  /*!
  \brief Destructor

  */
  virtual ~Vertex() {}

  //! @name Access methods

  /*!
  \brief Return vector of vertex coordinates (length 3)

  */
  virtual inline vector<double>& Coord() {return coord_; }

  /*!
  \brief Return vertex type (slave, projmaster or lineclip)

  */
  virtual inline Vertex::vType& VType() {return type_; }

  /*!
  \brief Return pointer to next vertex on polygon

  */
  virtual inline Vertex* Next() {return next_; }

  /*!
  \brief Assign pointer to next vertex on polygon

  */
  virtual inline void AssignNext(Vertex* assign) {next_=assign; return; }

  /*!
  \brief Return pointer to previous vertex on polygon

  */
  virtual inline Vertex* Prev() {return prev_; }

  /*!
  \brief Assign pointer to previous vertex on polygon

  */
  virtual inline void AssignPrev(Vertex* assign) {prev_=assign; return; }

  /*!
  \brief Return intersection status of this vertex
  True if vertex is an intersection point of the polygons.

  */
  virtual inline bool& Intersect() {return intersect_; }

  /*!
  \brief Return entry / exit status of this vertex
  True if vertex is an entry intersection point, false if
  vertex is an exit intesection point with respect to the
  respective other polygon. Irrelevant if intersect_==false.

  */
  virtual inline bool& EntryExit() {if (!intersect_) dserror("EntryExit only for intersections");
                                    return entryexit_; }

  /*!
  \brief Return pointer to neighbor on other polygon
  This pointer can only be set for an intersrection vertex,
  i.e. if intersect_==true. It then points to the identical
  vertex on the other polygon.

  */
  virtual inline Vertex* Neighbor() {if (!intersect_) dserror("Neighbor only for intersections");
                                     return neighbor_; }

  /*!
  \brief Assign pointer to neighbor on other polygon

  */
  virtual inline void AssignNeighbor(Vertex* assign) {if (!intersect_) dserror("Neighbor only for intersections");
                                                      neighbor_=assign; return; }

  /*!
  \brief Return intersection parameter alpha. Note that
  valid intersections yield an alpha in the range [0,1].

  */
  virtual inline double& Alpha() {return alpha_; }

  /*!
  \brief Return vector of relevant node ids (length 1 or 4)
  Note that for a slave or projmaster type vertex only the
  respective node id is relevant, therefore length 1. For a
  lineclip type vertex we need the 4 node ids of both the
  corresponding slave and master lines which intersect!

  */
  virtual inline vector<int>& Nodeids()
  { //if(type_==Vertex::slave && nodeids_.size()!=1) dserror("Error: Vertex Ids");
    //if(type_==Vertex::projmaster && nodeids_.size()!=1) dserror("Error: Vertex Ids");
    //if(type_==Vertex::lineclip && nodeids_.size()!=4) dserror("Error: Vertex Ids");
    return nodeids_; }

  //@}


protected:

  vector<double>   coord_;        // vertex coordinates (length 3)
  Vertex::vType    type_;         // vertex type (slave,projmaster,lineclip)
  vector<int>      nodeids_;      // relevant ids (1 if slave or master, 4 if lineclip)
  Vertex*          next_;         // pointer to next vertex on polygon
  Vertex*          prev_;         // pointer to previous vertex on polygon
  bool             intersect_;    // if true, this is an intersection vertex
  bool             entryexit_;    // if true, this is an entry vertex
  Vertex*          neighbor_;     // pointer to neighbor vertex on other polygon
  double           alpha_;        // intersection parameter

}; // class Vertex

}  // namespace CONTACT

#endif  // #ifndef DRT_CONTACT_COUPLING3D_CLASSES_H
#endif  // #ifdef CCADISCRET

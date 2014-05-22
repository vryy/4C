/*!----------------------------------------------------------------------
\file mortar_coupling3d_classes.H

\brief A class for mortar coupling of ONE slave element and ONE master
       element of a mortar interface in 3D (definition of sub-classes).

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
#ifndef MORTAR_COUPLING3D_CLASSES_H
#define MORTAR_COUPLING3D_CLASSES_H

#include "mortar_element.H"


/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms

*/
namespace MORTAR
{

// forward declarations

/*!
\brief A class representing a special mortar element for 3D quadratic
       mortar coupling with the use auf auxiliary planes. This approach
       is based on "Puso, M.A., Laursen, T.A., Solberg, J., A segment-to-
       segment mortar contact method for quadratic elements and large
       deformations, CMAME, 197, 2008, pp. 555-566". For this type of
       quadratic formulation, a quadratic MortarElement is split into several
       linear IntElements, on which the geometrical coupling is performed.
       Note that this is a derived class from MORTAR::MortarElement, with the
       only difference that we explicitly set the node pointers here!

\author popp (popp@lnm.mw.tum.de)
*/

class IntElement : public MortarElement{

public:

  /*!
  \brief Standard constructor

  \param lid   (in):    A locally unique element id
  \param id    (in):    Parent element global id
  \param owner (in):    owner processor of the element
  \param parshape (in): shape of parent element
  \param shape (in):    shape of this element
  \param numnode (in):  Number of nodes to this element
  \param nodeids (in):  ids of nodes adjacent to this element
  \param nodes (in):    pointers to nodes adjacent to this element
  \param isslave (in):  flag indicating whether element is slave or master side
  */
  IntElement(int lid, int id, int owner,
             const DRT::Element::DiscretizationType& parshape,
             const DRT::Element::DiscretizationType& shape,
             const int numnode,
             const int* nodeids,
             std::vector<DRT::Node*> nodes,
             const bool isslave);

  /*!
  \brief Destructor

  */
  virtual ~IntElement() {}

  /*!
  \brief Return local id of this IntElement

  */
  int Lid() const { return lid_; }

  /*!
  \brief Get shape type of parent element

  */
  virtual DiscretizationType ParShape() const { return parshape_;}

  /*!
  \brief Affine map of IntElement coordinates to parent element

  */
  virtual bool MapToParent(const double* xi, double* parxi);

  /*!
  \brief Affine map of IntElement coordinate derivatives to parent element

  */
  virtual bool MapToParent(const std::vector<std::map<int,double> >& dxi,
                           std::vector<std::map<int,double> >& dparxi);

protected:

  // don't want = operator and cctor
  IntElement operator = (const IntElement& old);
  IntElement(const IntElement& old);

  int                               lid_;        // local IntElement id
  DRT::Element::DiscretizationType  parshape_;   // shape of parent element

}; //class IntElement


/*!
\brief A class representing one Integration Cell after triangulation
       of the clip polygon of slave and master element from the Coupling
       class. This class provides some basic funcitonality a MortarElement would
       also provide (coords, shape functions and derivatives, Jacobian, ...).
       Note that an IntCell can EITHER live in physical space (this is the
       case when an auxiliary plane is used for 3D coupling) or in the slave
       element parameter space (this is the case when 3D coupling is performed
       on the slave surface without any auxiliary plane). Of course, in the
       2nd case the third coordinate of all intcell points is zero!

\author popp (popp@lnm.mw.tum.de)
*/

class IntCell{

public:

  /*!
  \brief Standard constructor

  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.

  */
  IntCell(int id, int nvertices, Epetra_SerialDenseMatrix& coords,
          double* auxn, const DRT::Element::DiscretizationType& shape,
          std::vector<std::map<int,double> >& linv1,
          std::vector<std::map<int,double> >& linv2,
          std::vector<std::map<int,double> >& linv3,
          std::vector<std::map<int,double> >& linauxn);

  /*!
  \brief Destructor

  */
  virtual ~IntCell() {}

  //! @name Access methods

  /*!
  \brief Return ID this intcell

  */
  int Id() const { return id_; }

  /*!
  \brief Return number of vertices of this intcell

  */
  int NumVertices() const { return nvertices_; }

  /*!
  \brief Return current area

  */
  virtual double& Area() { return area_; }

  /*!
  \brief Return coordinates of intcell vertices

  */
  const Epetra_SerialDenseMatrix& Coords() { return coords_; }

  /*!
  \brief Return normal of auxiliary plane of this intcell

  */
  virtual double* Auxn() { return auxn_; }

  /*!
  \brief Get shape type of element

  */
  virtual DRT::Element::DiscretizationType Shape() const { return shape_;}

  /*!
  \brief Return one of the three 'DerivVertex' maps (vectors) of this node

  These maps contain the directional derivatives of the intcell's vertex
  coordinates with respect to the slave and master displacements. A vector
  is used because the coordinates themselves are a vector (2 components for
  the 'coupling in slave parameter space' case). The respective vertex (1,2,3)
  is addressed by an int-variable and checked internally.

  */
  virtual std::vector<std::map<int,double> >& GetDerivVertex(int i)
  {
    if (i<0 || i>2) dserror("ERROR: IntCell has 3 vertex linearizations only!");
    return linvertex_[i];
  }

  /*!
  \brief Return the 'DerivAuxn' map (vector) of this intcell

  */
  virtual std::vector<std::map<int,double> >& GetDerivAuxn() { return linauxn_; }

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
  */
  virtual void DerivJacobian(double* xi, std::map<int,double>& derivjac);

  //@}

protected:

  // don't want = operator and cctor
  IntCell operator = (const IntCell& old);
  IntCell(const IntCell& old);

  int                        id_;              // local ID of this intcell
  int                        nvertices_;       // number of vertices (always 3)
  double                     area_;            // integration cell area
  Epetra_SerialDenseMatrix   coords_;          // coords of cell vertices
  double                     auxn_[3];         // normal of auxiliary plane (3D)

  bool                                              auxplane_;  // flag indicating coupling strategy (true = auxplane)
  DRT::Element::DiscretizationType                  shape_;     // shape of this element (always tri3)
  std::vector<std::vector<std::map<int,double> > >  linvertex_; // derivatives of the 3 vertices
  std::vector<std::map<int,double> >                linauxn_;   // derivatives of auxiliary plane normal

}; // class IntCell


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
  Vertex(std::vector<double> coord, Vertex::vType type,
         std::vector<int> nodeids, Vertex* next, Vertex* prev,
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
  virtual std::vector<double>& Coord() {return coord_; }

  /*!
  \brief Return vertex type (slave, projmaster or lineclip)

  */
  virtual Vertex::vType& VType() {return type_; }

  /*!
  \brief Return pointer to next vertex on polygon

  */
  virtual Vertex* Next() {return next_; }

  /*!
  \brief Assign pointer to next vertex on polygon

  */
  virtual void AssignNext(Vertex* assign) {next_=assign; return; }

  /*!
  \brief Return pointer to previous vertex on polygon

  */
  virtual Vertex* Prev() {return prev_; }

  /*!
  \brief Assign pointer to previous vertex on polygon

  */
  virtual void AssignPrev(Vertex* assign) {prev_=assign; return; }

  /*!
  \brief Return intersection status of this vertex
  True if vertex is an intersection point of the polygons.

  */
  virtual bool& Intersect() {return intersect_; }

  /*!
  \brief Return entry / exit status of this vertex
  True if vertex is an entry intersection point, false if
  vertex is an exit intesection point with respect to the
  respective other polygon. Irrelevant if intersect_==false.

  */
  virtual bool& EntryExit() {if (!intersect_) dserror("EntryExit only for intersections");
                                    return entryexit_; }

  /*!
  \brief Return pointer to neighbor on other polygon
  This pointer can only be set for an intersrection vertex,
  i.e. if intersect_==true. It then points to the identical
  vertex on the other polygon.

  */
  virtual Vertex* Neighbor() {if (!intersect_) dserror("Neighbor only for intersections");
                                     return neighbor_; }

  /*!
  \brief Assign pointer to neighbor on other polygon

  */
  virtual void AssignNeighbor(Vertex* assign) {if (!intersect_) dserror("Neighbor only for intersections");
                                                      neighbor_=assign; return; }

  /*!
  \brief Return intersection parameter alpha. Note that
  valid intersections yield an alpha in the range [0,1].

  */
  virtual double& Alpha() {return alpha_; }

  /*!
  \brief Return vector of relevant node ids (length 1 or 4)
  Note that for a slave or projmaster type vertex only the
  respective node id is relevant, therefore length 1. For a
  lineclip type vertex we need the 4 node ids of both the
  corresponding slave and master lines which intersect!

  */
  virtual std::vector<int>& Nodeids()
  { //if(type_==Vertex::slave && nodeids_.size()!=1) dserror("Error: Vertex Ids");
    //if(type_==Vertex::projmaster && nodeids_.size()!=1) dserror("Error: Vertex Ids");
    //if(type_==Vertex::lineclip && nodeids_.size()!=4) dserror("Error: Vertex Ids");
    return nodeids_; }

  //@}


protected:

  std::vector<double>  coord_;        // vertex coordinates (length 3)
  Vertex::vType        type_;         // vertex type (slave,projmaster,lineclip)
  std::vector<int>     nodeids_;      // relevant ids (1 if slave or master, 4 if lineclip)
  Vertex*              next_;         // pointer to next vertex on polygon
  Vertex*              prev_;         // pointer to previous vertex on polygon
  bool                 intersect_;    // if true, this is an intersection vertex
  bool                 entryexit_;    // if true, this is an entry vertex
  Vertex*              neighbor_;     // pointer to neighbor vertex on other polygon
  double               alpha_;        // intersection parameter

}; // class Vertex

}  // namespace MORTAR

#endif  // #ifndef MORTAR_COUPLING3D_CLASSES_H

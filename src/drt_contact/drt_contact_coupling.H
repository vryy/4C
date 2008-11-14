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
</pre>

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
       slave element and ONE master element of a contact interface in
       2D and 3D. Concretely, this class controls projection, overlap
       detection and finally integration of the mortar coupling matrix M.

\author popp (popp@lnm.mw.tum.de)
*/

class Coupling{
  
private:
  
/*!
\brief A class representing one Vertex during the polygon clipping of
       slave and master element from the Coupling class. Besides the
       vertex coordinates this class provides different pointers to
       build up doubly-linked list structures.

\author popp (popp@lnm.mw.tum.de)
*/
  
class Vertex{
  
public:
  
  /*!
  \brief Standard constructor
  
  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.
  
  */
  explicit Vertex(vector<double> coord, Vertex* next, Vertex* prev, bool intersect,
                  bool entryexit, Vertex* neighbor, double alpha);
  
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
  virtual inline vector<double> Coord() {return coord_; }
  
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
    
  //@}
  
  
protected:
  
  vector<double>   coord_;        // vertex coordinates (length 3)
  Vertex*          next_;         // pointer to next vertex on polygon
  Vertex*          prev_;         // pointer to previous vertex on polygon
  bool             intersect_;    // if true, this is an intersection vertex
  bool             entryexit_;    // if true, this is an entry vertex
  Vertex*          neighbor_;     // pointer to neighbor vertex on other polygon
  double           alpha_;        // intersection parameter
  
  
}; // class Vertex

// class Coupling
public:
  
  /*!
  \brief Standard constructor
  
  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as coupling is
  performed in parallel by individual processes.
  
  */
  explicit Coupling(DRT::Discretization& idiscret,
                    CONTACT::CElement& sele, CONTACT::CElement& mele,
                    int dim, Epetra_SerialDenseMatrix& csegs);

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
  \brief Get coupling master element

  */
  virtual inline CONTACT::CElement& MasterElement() const { return mele_; }
    
  /*!
  \brief Get problem dimension
  
  Note that only 2D and 3D are possible here as this refers to the global
  problem dimension. On interface level this corresponds to 1D interfaces
  (dim_==2) and 2D interfaces (dim_==3)!
      
  */
  virtual inline const int& Dim() {return dim_;};
  
  /*!
  \brief Return center of auxiliary plane (3D only)
  
  */
  virtual inline double* Auxc() { if (dim_!=3) dserror("Aux. only for 3D!");
                                  return auxc_; }
  
  /*!
  \brief Return normal of auxiliary plane (3D only)
  
  */
  virtual inline double* Auxn() { if (dim_!=3) dserror("Aux. only for 3D!");
                                  return auxn_; }
  
  /*!
  \brief Return vector of slave node coordinates
  
  */
  virtual inline vector<vector<double> >& SlaveVertices() {return svertices_; }
  
  /*!
  \brief Return vector of master node coordinates
  
  */
  virtual inline vector<vector<double> >& MasterVertices() {return mvertices_; }
    
  /*!
  \brief Return matrix of contact segment coordinates
  
  */
  virtual Epetra_SerialDenseMatrix& CSegs() {return contactsegs_;};
    
  //@}
    
  //! @name Evlauation methods
  
  /*!
  \brief Projection of slave / master pair (2D)
  
  This method projects the nodes of the slave CElement sele_ onto
  the master CElement mele_ and vice versa. The parameter hasproj
  stores a boolean variable for each of the 4 end nodes, indicating
  whether a feasible projection was found or not. The local element
  coordinates of the 4 projection points are stored in xiproj.
  
  */
  virtual bool Project2D(vector<bool>& hasproj,
                         vector<double>& xiproj);
  
  /*!
  \brief Detect overlap of slave / master pair (2D)
  
  This method evaluates the overlap of the current CElement pair
  sele_ / mele_ based on the projection status of the 4 end nodes
  (hasproj) and the coordinates of the projection points (xiproj).
  According to the detected overlap case, the integration limits
  are determined and written into xiproj and the overlap status
  is returned by a boolean variable.
  
  */
  virtual bool DetectOverlap2D(vector<bool>& hasproj,
                               vector<double>& xiproj);
  
  /*!
  \brief Integrate overlap of slave / master pair (2D)
  
  This method integrates the overlap of the current CElement
  pair sele_ / mele_ based on the integration limits (xiproj). The
  integration always includes the Mortar matrix M and the gap g.
  In some cases it also includes the Mortar matrix D (if the
  flag CONTACTONEMORTARLOOP is chosen) and the Mortar matrix
  Mmod (if the respective conditions are fulfilled).
        
  */
  virtual bool IntegrateOverlap2D(vector<double>& xiproj);
  
  /*!
  \brief Build auxiliary plane from slave element (3D)
  
  This method builds an auxiliary plane based on the possibly
  warped slave element of this coupling class. This plane is
  defined by the slave normal at the slave element center.
        
  */
  virtual bool AuxiliaryPlane3D();
  
  /*!
  \brief Projection of slave element onto aux. plane (3D)
  
  This method projects the nodes of the given slave CElement
  onto the auxiliary plane derived before.
  
  */
  virtual bool ProjectSlave3D();
    
  /*!
  \brief Projection of master element onto aux. plane (3D)
  
  This method projects the nodes of the current master CElement
  onto the auxiliary plane derived from the slave CElement before.
  
  */
  virtual bool ProjectMaster3D();
  
  /*!
  \brief Clipping of slave and master element in aux. plane (3D)
  
  This method applies a polygon clipping algorithm to find the
  polygon defined by the projection of the slave and master element
  in the auxiliary plane. As input variables the clipping algorithm
  requires two sets of vertices, which define sele_ / mele_ respectively.
  The clipping algorithm is based on doubly linked lists as data structure
  (Greiner, G., Hormann, K. (1998), Efficient Clipping of Arbitrary
  Polygons, ACM Transactions on Graphics, vol. 17, pp. 71-83).
  Returns a vector of vertices defining the clipped polygon.
  
  \param polygon1 (in): vector of vertices for input polygon 1
  \param polygon2 (in): vector of vertices for input polygon 2
  
  */
  virtual vector<vector<double> > PolygonClipping(vector<vector<double> > polygon1,
                                                  vector<vector<double> > polygon2);
  //@}
    
      
protected:
    
    // don't want = operator and cctor
    Coupling operator = (const Coupling& old);
    Coupling(const Coupling& old);
  
    
    DRT::Discretization&       idiscret_;    // discretization of the contact interface
    CONTACT::CElement&         sele_;        // slave element to perform coupling for
    CONTACT::CElement&         mele_;        // master element to perform coupling for
    int                        dim_;         // dimension of problem (2D or 3D)
    
    double                     auxc_[3];     // center of auxiliary plane (3D)
    double                     auxn_[3];     // normal of auxiliary plane (3D)
    vector<vector<double> >    svertices_;   // slave node coordinates
    vector<vector<double> >    mvertices_;   // master node coordinates
    vector<vector<double> >    clip_;        // clipped polygon coordinates
    
    Epetra_SerialDenseMatrix&  contactsegs_; // contains coordinates of contact segments (gmsh)
  

}; // class Coupling
}  // namespace CONTACT


#endif  // #ifndef DRT_CONTACT_COUPLING_H
#endif  // #ifdef CCADISCRET

/*-----------------------------------------------------------------------*/
/*! \file
\brief Classes for mortar coupling in 3D.

\level 2

*/
/*-----------------------------------------------------------------------*/
#ifndef MORTAR_COUPLING3D_H
#define MORTAR_COUPLING3D_H

#include <vector>
#include "Epetra_Comm.h"
#include "mortar_coupling3d_classes.H"
#include "../drt_inpar/inpar_mortar.H"
#include "../headers/pairedvector.H"

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace MORTAR
{
  // forward declarations
  class MortarElement;
  class ParamsInterface;

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

   */

  class Coupling3d
  {
   public:
    /*!
     \brief Constructor with shape function specification

     Constructs an instance of this class and enables custom shape function types.<br>
     Note that this is \b not a collective call as coupling is
     performed in parallel by individual processes.

     */
    Coupling3d(DRT::Discretization& idiscret, int dim, bool quad, Teuchos::ParameterList& params,
        MORTAR::MortarElement& sele, MORTAR::MortarElement& mele);

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
    virtual const Epetra_Comm& Comm() const;

    /*!
     \brief Get problem dimension (here: 3D)

     */
    virtual const int& Dim() { return dim_; };

    /*!
     \brief Get coupling / FE ansatz type (true = quadratic)

     */
    virtual const bool& Quad() { return quad_; };

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
    virtual MORTAR::MortarElement& MasterIntElement() const { return mele_; }

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
    virtual double& Lauxn() { return lauxn_; }

    /*!
     \brief Return vector of (projected) slave node vertex objects

     */
    virtual std::vector<Vertex>& SlaveVertices() { return svertices_; }

    /*!
     \brief Return vector of projected master node vertex objects

     */
    virtual std::vector<Vertex>& MasterVertices() { return mvertices_; }

    /*!
     \brief Return vector of clip polygon vertex objects

     */
    virtual std::vector<Vertex>& Clip() { return clip_; }

    /*!
     \brief Return vector of integration cells

     */
    virtual std::vector<Teuchos::RCP<IntCell>>& Cells() { return cells_; }

    /*!
     \brief Return the 'DerivAuxn' map (vector) of this coupling pair

     */
    virtual std::vector<GEN::pairedvector<int, double>>& GetDerivAuxn() { return derivauxn_; }

    /*!
     \brief Return the LM interpolation / testing type for quadratic FE

     */
    virtual INPAR::MORTAR::LagMultQuad LagMultQuad() { return lmquadtype_; };

    /*!
     \brief Get interface contact parameter list

     */
    virtual Teuchos::ParameterList& InterfaceParams() { return imortar_; };

    /*!
     \brief Return the LM shape fcn type

     */
    INPAR::MORTAR::ShapeFcn ShapeFcn() { return shapefcn_; };

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

     This methods checks the orientation of slave and master element.
     Projection and overlap detection only make sense if the scalar
     product of the two center normals is negative (i.e. if the two
     elements form an angle smaller than 90°).

     */
    virtual bool RoughCheckOrient();

    /*!
     \brief Integrate the integration cells (3D)

     This method creates an integrator object for the cell triangles,
     then projects the Gauss points back onto slave and master elements
     (1st case, aux. plane) or only back onto the master element (2nd case)
     in order to evaluate the respective shape function there. Then
     entries of the mortar matrix M and the weighted gap g are integrated
     and assembled into the slave element nodes.

     */
    virtual bool IntegrateCells(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

    //@}

    //! @name Visualization and Debugging methods

    /*!
     \brief Visualize integration cells with gmsh

     */
    virtual void GmshOutputCells(int lid);

    //@}

    //! @name Linearization methods (for mortar contact only!)

    /*!
     \brief Linearization of clip vertex coordinates (3D)

     */
    virtual bool VertexLinearization(
        std::vector<std::vector<GEN::pairedvector<int, double>>>& linvertex,
        std::map<int, double>& projpar, bool printderiv = false)
    {
      return true;
    }

    /*!
     \brief Linearization of clip vertex coordinates (3D)

     */
    virtual bool CenterLinearization(
        const std::vector<std::vector<GEN::pairedvector<int, double>>>& linvertex,
        std::vector<GEN::pairedvector<int, double>>& lincenter)
    {
      return true;
    }

    //@}
   private:
    /*!
     \brief Build auxiliary plane from slave element (3D)

     This method builds an auxiliary plane based on the possibly
     warped slave element of this coupling class. This plane is
     defined by the slave normal at the slave element center.

     */
    virtual bool AuxiliaryPlane();

    /*!
     \brief Triangulation of clip polygon (3D)

     In a first step, this method computes linearizations of all clip
     polygon vertices and stores them into local maps. Then, one of the
     two available triangulation algorithms (Delaunay, Center) is called
     to do the actual triangulation of the convex clip polygon.

     DelaunayTriangulation() is more efficient, as it always generates N-2 triangles
     for a clip polygon with N vertices. Moreover, the triangle shape is optimized
     according to the Delaunay criterion, which maximizes the smallest internal
     angle. Thus, DelaunayTriangulation() is used BY DEFAULT here.

     CenterTriangulation() is an old version that was based on first finding the
     center of the clip polygon by computing the centroid / geometric center
     (see M. Puso, A 3D mortar method for solid mechanics, IJNME, 2004).
     The clip polygon is triangulated into integration cells all containing two
     vertices and the center point. Thus, we always generate N triangles here
     for a clip polygon with N vertices (except for the special case N=3).
     Thus, CenterTriangulation() is NOT used anymore.

     */
    virtual bool Triangulation(std::map<int, double>& projpar, double tol);
    virtual bool DelaunayTriangulation(
        std::vector<std::vector<GEN::pairedvector<int, double>>>& linvertex, double tol);
    virtual bool CenterTriangulation(
        std::vector<std::vector<GEN::pairedvector<int, double>>>& linvertex, double tol);

    /*!
     \brief Check / set projection status of slave nodes (3D)

     This method checks for all slave nodes if they are part of the clip
     polygon (equal to any vertex). If so the HasProj status is set true!

     */
    virtual bool HasProjStatus();

    /*!
     \brief Compute and return area of clipping polygon (3D)

     This method automatically returns the relevant(!) area of the
     current clip polygon, depending on whether coupling is performed
     in the auxiliary plane or in the slave parameter space.

     */
    virtual double PolygonArea();

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
    virtual void PolygonClipping(std::vector<Vertex>& poly1list, std::vector<Vertex>& poly2list,
        std::vector<Vertex>& respoly, double& tol);

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
    virtual bool PolygonClippingConvexHull(std::vector<Vertex>& poly1list,
        std::vector<Vertex>& poly2list, std::vector<Vertex>& respoly, double& tol);

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
     \brief Checks roughly whether the two elements are near (3D)

     This methods computes the distance of all master nodes to the
     slave element (auxiliary plane). If they are not close, then
     coupling is stopped for the pair.

     */
    virtual bool RoughCheckNodes();

   protected:
    /*!
     \brief Compute and return area of slave element (3D)

     This method automatically returns the relevant(!) area of the
     current slave element, depending on whether coupling is performed
     in the auxiliary plane or in the slave parameter space.

     */
    virtual double SlaveElementArea();

    // don't want = operator and cctor
    Coupling3d operator=(const Coupling3d& old);
    Coupling3d(const Coupling3d& old);

    DRT::Discretization& idiscret_;          // discretization of the contact interface
    int dim_;                                // problem dimension (here: 3D)
    INPAR::MORTAR::ShapeFcn shapefcn_;       // lm shape function type
    bool quad_;                              // flag indicating coupling type (true = quadratic)
    INPAR::MORTAR::LagMultQuad lmquadtype_;  // type of quadratic lm interpolation

    MORTAR::MortarElement& sele_;      // slave element to perform coupling for
    MORTAR::MortarElement& mele_;      // master element to perform coupling for
    Teuchos::ParameterList& imortar_;  // containing contact input parameters
    double auxc_[3];                   // center of auxiliary plane
    double auxn_[3];                   // normal of auxiliary plane
    double lauxn_;                     // length of interpolated Auxn() before normalization
    std::vector<Vertex> svertices_;    // slave node vertex objects
    std::vector<Vertex> mvertices_;    // master node vertex objects
    std::vector<Vertex> clip_;         // clipped polygon vertex objects

    std::vector<Teuchos::RCP<IntCell>> cells_;  // vector of integration cells
    std::vector<GEN::pairedvector<int, double>>
        derivauxn_;  // derivatives of auxiliary plane normal
  };
  // class Coupling3d

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

   */

  class Coupling3dQuad : public Coupling3d
  {
   public:
    /*!
     \brief Constructor with shape function specification

     Constructs an instance of this class and enables custom shape function types.<br>
     Note that this is \b not a collective call as coupling is
     performed in parallel by individual processes.

     */
    Coupling3dQuad(DRT::Discretization& idiscret, int dim, bool quad,
        Teuchos::ParameterList& params, MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        MORTAR::IntElement& sintele, MORTAR::IntElement& mintele);

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

   protected:
    // don't want = operator and cctor
    Coupling3dQuad operator=(const Coupling3dQuad& old);
    Coupling3dQuad(const Coupling3dQuad& old);

    MORTAR::IntElement& sintele_;  // slave sub-integration element
    MORTAR::IntElement& mintele_;  // slave sub-integration element
  };
  // class Coupling3dQuad

  /*!
   \brief A class representing the framework for mortar coupling of ONE
   slave element and SEVERAL master elements of a mortar interface in
   3D. Concretely, this class simply stores several Coupling3d objects.

   */

  class Coupling3dManager
  {
   public:
    /*!
     \brief Constructor with shape function specification

     Constructs an instance of this class and enables custom shape function types.<br>
     Note that this is \b not a collective call as coupling is
     performed in parallel by individual processes.

     */
    Coupling3dManager(DRT::Discretization& idiscret, int dim, bool quad,
        Teuchos::ParameterList& params, MORTAR::MortarElement* sele,
        std::vector<MORTAR::MortarElement*> mele);

    /*!
     \brief Destructor

     */
    virtual ~Coupling3dManager() {}

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
    virtual std::vector<Teuchos::RCP<MORTAR::Coupling3d>>& Coupling() { return coup_; }

    /*!
     \brief Get type of integration scheme

     */
    virtual INPAR::MORTAR::IntType IntType() { return integrationtype_; }

    /*!
     \brief Get coupling type

     */
    virtual const bool& Quad() { return quad_; };

    /*!
     \brief Get flag indicating dual consistent lm

     */
    virtual bool LMDualConsistent() { return lmdualconsistent_; }

    /*!
     \brief Get communicator

     */
    virtual const Epetra_Comm& Comm() const;

    /*!
     \brief Evaluate coupling pairs

     */
    virtual bool EvaluateCoupling(Teuchos::RCP<MORTAR::ParamsInterface> mparams_ptr);

    /*!
     \brief Return the LM interpolation / testing type for quadratic FE

     */
    virtual INPAR::MORTAR::LagMultQuad LagMultQuad() { return lmquadtype_; };

    /*!
     \brief Return the LM shape fcn type

     */
    virtual INPAR::MORTAR::ShapeFcn ShapeFcn() { return shapefcn_; }
    //@}

   protected:
    /*!
     \brief Evaluate mortar-coupling pairs

     */
    virtual void IntegrateCoupling(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

    /*!
     \brief Calculate consistent dual shape functions in boundary elements

     */
    virtual void ConsistDualShape();

    // don't want = operator and cctor
    Coupling3dManager operator=(const Coupling3dManager& old);
    Coupling3dManager(const Coupling3dManager& old);

    DRT::Discretization& idiscret_;                       // discretization of the contact interface
    int dim_;                                             // problem dimension (here: 3D)
    INPAR::MORTAR::IntType integrationtype_;              // integration type
    INPAR::MORTAR::ShapeFcn shapefcn_;                    // lm shape function type
    INPAR::MORTAR::LagMultQuad lmquadtype_;               // type of quadratic lm interpolation
    INPAR::MORTAR::ConsistentDualType lmdualconsistent_;  // flag indicating dual consistent lm
    bool quad_;                                 // flag indicating coupling type (true = quadratic)
    Teuchos::ParameterList& imortar_;           // containing contact input parameters
    MORTAR::MortarElement* sele_;               // slave element
    std::vector<MORTAR::MortarElement*> mele_;  // master elements
    std::vector<Teuchos::RCP<MORTAR::Coupling3d>> coup_;  // coupling pairs
  };
  // class Coupling3dManager

  /*!
   \brief A class representing the framework for mortar coupling of ONE
   slave element and SEVERAL master elements of a mortar interface in
   3D with QUADRATIC ELEMENTS. Concretely, this class simply sub-divides
   elements and then stores several Coupling3d objects.

  */

  class Coupling3dQuadManager : public Coupling3dManager
  {
   public:
    /*!
     \brief Constructor
     */
    Coupling3dQuadManager(DRT::Discretization& idiscret, int dim, bool quad,
        Teuchos::ParameterList& params, MORTAR::MortarElement* sele,
        std::vector<MORTAR::MortarElement*> mele);

    /*!
     \brief Destructor

     */
    virtual ~Coupling3dQuadManager() {}

    INPAR::MORTAR::LagMultQuad LagMultQuad()
    {
      return DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(imortar_, "LM_QUAD");
    }

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
     \brief Get integration type

     */
    INPAR::MORTAR::IntType IntType() { return integrationtype_; }

    /*!
     \brief Get coupling type

     */
    virtual const bool& Quad() { return quad_; }

    // @

   protected:
    /*!
     \brief Evaluate mortar-coupling pairs

     */
    virtual void IntegrateCoupling(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

    /*!
     \brief Split MortarElements into IntElements for 3D quadratic coupling

     */
    virtual bool SplitIntElements(
        MORTAR::MortarElement& ele, std::vector<Teuchos::RCP<MORTAR::IntElement>>& auxele);

    // don't want = operator and cctor
    Coupling3dQuadManager operator=(const Coupling3dQuadManager& old);
    Coupling3dQuadManager(const Coupling3dQuadManager& old);
  };

}  // namespace MORTAR

#endif  // #ifndef MORTAR_COUPLING3D_H

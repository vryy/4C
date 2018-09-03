/*---------------------------------------------------------------------*/
/*!
\file cut_element.H

\brief cut element

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#ifndef CUT_ELEMENT_H
#define CUT_ELEMENT_H

#include "cut_integrationcell.H"
#include "cut_utils.H"
#include "cut_node.H"

#include "../drt_inpar/inpar_cut.H"

namespace GEO
{
  namespace CUT
  {
    class VolumeCell;
    class IntegrationCell;
    class BoundaryCell;
    class BoundingBox;

    /*--------------------------------------------------------------------------*/
    /*! \brief Base class of all elements.
     *
     * Contains all the information about volumecells, facets, sides and nodes
     * associated with this element. This class handles only linear elements. When
     * a quadratic element is passed to cut library, it is split into a number of
     * linear elements through elementhandle.
     *
     * \note The discretization specific stuff is moved to the ConcreteElement
     * class. All methods of the concrete class have (or should have) an accessor
     * in the base class. In this way we can circumvent an explicit knowledge
     * about the discretization type, problem dimension etc., while using the
     * different methods.                               note by hiermeier 08/16 */
    class Element
    {
     public:
      /// \brief create an element with the given type
      static Teuchos::RCP<GEO::CUT::Element> Create(
          const DRT::Element::DiscretizationType& elementtype, const int& eid,
          const std::vector<Side*>& sides, const std::vector<Node*>& nodes, const bool& active);

      /// create an element with the given shards-key (coming from trilinos library)
      static Teuchos::RCP<GEO::CUT::Element> Create(const unsigned& shardskey, const int& eid,
          const std::vector<Side*>& sides, const std::vector<Node*>& nodes, const bool& active);

      //! constructor
      Element(
          int eid, const std::vector<Side*>& sides, const std::vector<Node*>& nodes, bool active);

      //! destructor
      virtual ~Element() {}

      /*! \brief Returns the ID of the element */
      int Id() const { return eid_; }

      /*! \brief Returns the shape of the element */
      virtual DRT::Element::DiscretizationType Shape() const = 0;

      //! element dimension
      virtual unsigned Dim() const = 0;

      //! element number of nodes
      virtual unsigned NumNodes() const = 0;

      //! problem dimension
      virtual unsigned ProbDim() const = 0;

      /*! \brief Returns true if the point lies inside the element */
      virtual bool PointInside(Point* p) = 0;

      /*! \brief Returns element local coordinates "rst" of a point from its global coordinates
       * "xyz"
       *
       *  \remark This variant uses the LINALG::Matrices as input and does a dimension check.
       *  Call this function, to be on the safe side.
       *
       *  \param xyz (in) : Global spatial coordinate
       *  \param rst (out): Local parameter space coordinate */
      template <class T1, class T2>
      bool LocalCoordinates(const T1& xyz, T2& rst)
      {
        if (static_cast<unsigned>(xyz.M()) < ProbDim())
          dserror("The dimension of xyz is wrong! (probdim = %d)", ProbDim());
        if (static_cast<unsigned>(rst.M()) < Dim())
          dserror("The dimension of rst is wrong! (dim = %d)", Dim());

        bool success = LocalCoordinates(xyz.A(), rst.A());

        std::fill(rst.A() + Dim(), rst.A() + rst.M(), 0.0);

        return success;
      }

      /*! Returns element global coordinates "xyz" of a point from its local coordinates "rst"
       *
       *  \remark This variant uses the LINALG::Matrices as input and does a dimensional check.
       *  Call this function to be on the safe side.
       *
       *  \param rst (in)  : Local parameter space coordinate
       *  \param xyz (out) : Global spatial coordinate */
      template <class T1, class T2>
      void GlobalCoordinates(const T1& rst, T2& xyz)
      {
        if (static_cast<unsigned>(xyz.M()) < ProbDim())
          dserror("The dimension of xyz is wrong! (probdim = %d)", ProbDim());
        if (static_cast<unsigned>(rst.M()) < Dim())
          dserror("The dimension of rst is wrong! (dim = %d)", Dim());

        GlobalCoordinates(rst.A(), xyz.A());

        std::fill(xyz.A() + ProbDim(), xyz.A() + xyz.M(), 0.0);
      }

      /*! \brief Returns the element center global coordinates
       *
       *  \remark This variant uses the LINALG::Matrix as input and does a dimensional check.
       *  Call this function to be on the safe side.
       *
       *  \param midpoint (out): The midpoint of the element in global spatial coordinates */
      template <class T>
      void ElementCenter(T& midpoint)
      {
        if (static_cast<unsigned>(midpoint.M()) != ProbDim())
          dserror("The dimension of midpoint is wrong! (probdim = %d)", ProbDim());

        return ElementCenter(midpoint.A());
      }

      /*! \brief Find the scalar value at a particular point inside the element
       *  specified by its local coordinates rst.
       *
       *  \remark This variant uses the LINALG::Matrix as input and does a dimensional check.
       *  Call this function to be on the safe side.
       *
       *  \param ns (in)   : contains the value of the scalar at the corner nodes of the element
       *  \param rst (in)  : local parameter space coordinate inside the element */
      template <class T>
      double Scalar(const std::vector<double>& ns, const T& rst)
      {
        if (static_cast<unsigned>(rst.M()) != Dim())
          dserror("The dimension of rst is wrong! (dim = %d)", Dim());
        return Scalar(ns, rst.A());
      }

      //! \brief Cutting the element with the given cut side
      bool Cut(Mesh& mesh, Side& cut_side, int recursion);

      /*! \brief cut this element with its cut faces */
      void FindCutPoints(Mesh& mesh, int recursion = 0);

      /*! \brief cut this element with given cut_side */
      bool FindCutPoints(Mesh& mesh, Side& cut_side, int recursion);

      /*! \brief Create cut lines over this element by connecting appropriate cut points */
      void MakeCutLines(Mesh& mesh);

      /*! \brief Create facets */
      void MakeFacets(Mesh& mesh);

      /*! \brief Create volumecells */
      void MakeVolumeCells(Mesh& mesh);

      /*! \brief Create integrationcells by performing tessellation over the volumecells
       *  of the element. This uses QHULL */
      void CreateIntegrationCells(Mesh& mesh, int count, bool tetcellsonly = false);

      /*! \brief Try to create simpled shaped integrationcells */
      bool CreateSimpleShapedIntegrationCells(Mesh& mesh);

      /*! \brief Construct quadrature rules for the volumecells by solving the moment
       *  fitting equations */
      void MomentFitGaussWeights(
          Mesh& mesh, bool include_inner, INPAR::CUT::BCellGaussPts Bcellgausstype);

      /*! \brief Construct quadrature rules for the volumecells by triangulating the
       *  facets and applying divergence theorem */
      void DirectDivergenceGaussRule(
          Mesh& mesh, bool include_inner, INPAR::CUT::BCellGaussPts Bcellgausstype);

      /*! \brief Return the level set value at the given global coordinate
       *  which has to be INSIDE the element. */
      template <class T>
      double GetLevelSetValue(const T& xyz)
      {
        if (static_cast<unsigned>(xyz.M()) < ProbDim())
          dserror("The dimension of xyz is wrong! (probdim = %d)", ProbDim());

        return GetLevelSetValue(xyz.A());
      }

      /*! \brief Return the level set value at the given local coordinate
       *  which has to be INSIDE the element. */
      template <class T>
      double GetLevelSetValueAtLocalCoords(const T& rst)
      {
        if (static_cast<unsigned>(rst.M()) < Dim())
          dserror("The dimension of rst is wrong! (dim = %d)", Dim());

        return GetLevelSetValueAtLocalCoords(rst.A());
      }

      /*! \brief Return the level set gradient at the given global coordinate
       *  which has to be INSIDE the element.
       *
       *  \comment Might need to add a similar function for cut_side as it uses
       *  the LevelSetGradient information to determine cut cases with
       *  4 nodes on a side. Right now this function is called.
       *
       *  \remark This function is necessary because the orientation of a facet is not
       *  considered in its creation, this could be solved by introducing this
       *  information earlier in the facet creation.
       *
       *  For example:
       *  + Get cut_node and its two edges on the cut_side. Calculate the
       *    cross-product
       *                                   -> Get the orientation of the node.
       *  + Compare to LS info from its two edges which are not on a cut_side.
       *  + Make sure, the facet is created according to the calculated
       *    orientation.
       *                                           originally     winter 07/15
       *                                           generalized hiermeier 07/16 */
      template <class T>
      std::vector<double> GetLevelSetGradient(const T& xyz)
      {
        if (static_cast<unsigned>(xyz.M()) < ProbDim())
          dserror("The dimension of xyz is wrong! (probdim = %d)", ProbDim());

        return GetLevelSetGradient(xyz.A());
      }
      template <class T>
      std::vector<double> GetLevelSetGradientAtLocalCoords(const T& rst)
      {
        if (static_cast<unsigned>(rst.M()) < Dim())
          dserror("The dimension of rst is wrong! (dim = %d)", Dim());

        return GetLevelSetGradientAtLocalCoords(rst.A());
      }
      template <class T>
      std::vector<double> GetLevelSetGradientInLocalCoords(const T& xyz)
      {
        if (static_cast<unsigned>(xyz.M()) < ProbDim())
          dserror("The dimension of xyz is wrong! (probdim = %d)", ProbDim());

        return GetLevelSetGradientInLocalCoords(xyz.A());
      }
      template <class T>
      std::vector<double> GetLevelSetGradientAtLocalCoordsInLocalCoords(const T& rst)
      {
        if (static_cast<unsigned>(rst.M()) < Dim())
          dserror("The dimension of rst is wrong! (probdim = %d)", ProbDim());

        return GetLevelSetGradientAtLocalCoordsInLocalCoords(rst.A());
      }
      /*! \brief Does the element contain a level set side. */
      bool HasLevelSetSide();

      /** \brief Artificial cells with zero volumes are removed
       *
       *  After the integration and boundary cells have been created, there might
       *  still be some empty volume cells. Those are artificial cells with zero
       *  volume and need to be removed from their elements, before we count nodal
       *  dof sets. */
      void RemoveEmptyVolumeCells();

      /*! \brief Determine the inside/outside/on cut-surface position for the element's nodes */
      void FindNodePositions();

      /** \brief main routine to compute the position based on the angle between the
       *  line-vec (p-c) and an appropriate cut-side
       *
       *  \remark The following inside/outside position is based on comparisons of the
       *  line vec between point and the cut-point and the normal vector w.r.t cut-side
       *  "angle-comparison". In case the cut-side is not unique we have to determine at
       *  least one cut-side which can be used for the "angle-criterion" such that the right
       *  position is guaranteed.
       *  In case that the cut-point lies on an edge between two different cut-sides or even
       *  on a node between several cut-sides, we have to find a side (maybe not unique,
       *  but at least one of these) that defines the right position based on the "angle-criterion".
       *
       * \param p        : the point for that the position has to be computed
       * \param cutpoint : the point on cut side which is connected to p via a facets line)
       * \param f        : the facet via which p and cutpoint are connected
       * \param s        : the current cut side, the cutpoint lies on */
      bool ComputePosition(Point* p, Point* cutpoint, Facet* f, Side* s);

      /** \brief determine the position of point p based on the angle between the line (p-c)
       * and the side's normal vector, return \TRUE if successful */
      bool PositionByAngle(Point* p, Point* cutpoint, Side* s);

      /** \brief returns true in case that any cut-side cut with the element produces cut points,
       *  i.e. also for touched cases (at points, edges or sides), or when an element side has more
       *  than one facet or is touched by fully/partially by the cut side */
      bool IsCut();

      /** \brief return true if the element has more than one volume-cell and therefore
       *  is intersected by a cut-side */
      bool IsIntersected() { return (this->NumVolumeCells() > 1); }

      bool OnSide(Facet* f);

      bool OnSide(const std::vector<Point*>& facet_points);

      /*! \brief Get the integrationcells created from this element */
      INPAR::CUT::ElementIntegrationType GetElementIntegrationType() { return eleinttype_; }

      /*! \brief Get the integrationcells created from this element */
      void GetIntegrationCells(plain_integrationcell_set& cells);

      /*! \brief Get the boundarycells created from cut facets of this element */
      void GetBoundaryCells(plain_boundarycell_set& bcells);

      /*! \brief Returns all the sides of this element */
      const std::vector<Side*>& Sides() const { return sides_; }

      /*! \brief Returns true if the side is one of the element sides */
      bool OwnedSide(Side* side)
      {
        return std::find(sides_.begin(), sides_.end(), side) != sides_.end();
      }

      /*! \brief Get the coordinates of all the nodes of this linear (shadow) element */
      template <class T>
      void Coordinates(T& xyze) const
      {
        if (static_cast<unsigned>(xyze.M()) != ProbDim())
          dserror("xyze has the wrong number of rows! (probdim = %d)", ProbDim());

        if (static_cast<unsigned>(xyze.N()) != NumNodes())
          dserror("xyze has the wrong number of columns! (numNodes = %d)", NumNodes());

        Coordinates(xyze.A());
      }

      /*! \brief Get the coordinates of all the nodes of this linear shadow element */
      virtual void Coordinates(double* xyze) const = 0;

      /*! \brief Get the coordinates of all node of parent Quad element */
      virtual void CoordinatesQuad(double* xyze) = 0;

      /*! \brief Get all the nodes of this element */
      const std::vector<Node*>& Nodes() const { return nodes_; }

      const std::vector<Node*>& QuadCorners() const { return quadCorners_; }

      /*! \brief Points that are associated with nodes of this shadow element */
      const std::vector<Point*>& Points() const { return points_; }

      /*! \brief Get the cutsides of this element */
      const plain_side_set& CutSides() const { return cut_faces_; }

      /*! \brief Get cutpoints of this element */
      void GetCutPoints(PointSet& cut_points);

      /*! \brief Get the volumecells of this element */
      const plain_volumecell_set& VolumeCells() const { return cells_; }

      /*! \brief Get the number of volume-cells */
      int NumVolumeCells() const { return cells_.size(); }

      /*! \brief Get the facets of this element */
      const plain_facet_set& Facets() const { return facets_; }

      /*! \brief Are there any facets? */
      bool IsFacet(Facet* f) { return facets_.count(f) > 0; }

      /*! \brief check if the side's normal vector is orthogonal to the line
       *  between p and the cutpoint */
      bool IsOrthogonalSide(Side* s, Point* p, Point* cutpoint);

      /*!
      \brief Get total number of Gaussinan points generated over all the volumecells of this element
       */
      int NumGaussPoints(DRT::Element::DiscretizationType shape);

      void DebugDump();

      void GnuplotDump();

      /*! \brief When the cut library breaks down, this writes the element
       * geometry and all the cut sides to inspect visually whether the cut
       * configuration is appropriate */
      void GmshFailureElementDump();

      void DumpFacets();

      /*!
      \brief Inset this volumecell to this element
       */
      void AssignOtherVolumeCell(VolumeCell* vc) { cells_.insert(vc); }

      /*!
      \brief Calculate the volume of all the volumecells associated with this element only when
      tessellation is used
       */
      void CalculateVolumeOfCellsTessellation();

      /*!
      \brief Transform the scalar variable either from local coordinates to the global coordinates
      or the other way specified by TransformType which can be "LocalToGlobal" or "GlobalToLocal".
      If set shadow = true, then the mapping is done w.r to the parent Quad element from which this
      shadow element is derived
       */
      template <int probdim, DRT::Element::DiscretizationType distype>
      double ScalarFromLocalToGlobal(double scalar, std::string transformType, bool shadow = false)
      {
        const int nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
        const int dim = DRT::UTILS::DisTypeToDim<distype>::dim;
        LINALG::Matrix<probdim, nen> xyze;
        LINALG::Matrix<probdim, 1> xei;
        xei = 0.1;  // the determinant of the Jacobian is independent of the considered point (true
                    // only for linear elements???)

        // get coordiantes of parent Quad element
        if (shadow) CoordinatesQuad(xyze.A());
        // get coordinates of linear shadow element
        else
          Coordinates(xyze.A());

        LINALG::Matrix<dim, nen> deriv;
        LINALG::Matrix<probdim, probdim> xjm;
        DRT::UTILS::shape_function_deriv1<distype>(xei, deriv);

        double det = 0;
        if (dim == probdim)
        {
          xjm.MultiplyNT(deriv, xyze);
          det = xjm.Determinant();
        }
        // element dimension is smaller than the problem dimension (manifold)
        else if (dim < probdim)
        {
          LINALG::Matrix<dim, dim> metrictensor;
          DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype, probdim>(
              xyze, deriv, metrictensor, det, NULL);
        }
        else
          dserror("The element dimension is larger than the problem dimension! (%d > %d)", dim,
              probdim);
        double VarReturn = 0.0;
        if (transformType == "LocalToGlobal")
          VarReturn = scalar * det;
        else if (transformType == "GlobalToLocal")
          VarReturn = scalar / det;
        else
          dserror("Undefined transformation option");
        return VarReturn;
      }

      /*! \brief Integrate pre-defined functions over each volumecell
       * created from this element when using Tessellation */
      void integrateSpecificFunctionsTessellation();

      /*! \brief Assign the subelement his parent id (equal to the subelement
       * id eid_ for linear elements) */
      void ParentId(int id) { parent_id_ = id; }

      /*! \brief Gets the parent id of the subelement (equal to the sub-element
       * id eid_ for linear elements) */
      int GetParentId() { return parent_id_; }

      /*!
       \brief set this element as shadow element
       */
      void setAsShadowElem() { isShadow_ = true; }

      /*!
      \brief Return true if this is a shadow element
       */
      bool isShadow() { return isShadow_; }

      /*!
      \brief Store the corners of parent Quad element from which this shadow element is derived
       */
      void setQuadCorners(Mesh& mesh, const std::vector<int>& nodeids);

      /*!
      \brief Get corner nodes of parent Quad element from which this shadow element is derived
       */
      std::vector<Node*> getQuadCorners();

      /*!
      \brief Set the discretization type of parent quad element
       */
      void setQuadShape(DRT::Element::DiscretizationType dis) { quadshape_ = dis; }

      /*!
      \brief Get shape of parent Quad element from which this shadow element is derived
       */
      DRT::Element::DiscretizationType getQuadShape() { return quadshape_; }

      /*!
      \brief calculate the local coordinates of "xyz" with respect to the parent Quad element from
      which this shadow element is derived
       */
      void LocalCoordinatesQuad(const LINALG::Matrix<3, 1>& xyz, LINALG::Matrix<3, 1>& rst);

      /*!
      \brief Add cut face
       */
      void AddCutFace(Side* cutface) { cut_faces_.insert(cutface); }

      /*!
      \brief get the bounding volume
       */
      const BoundingBox& GetBoundingVolume() const { return *boundingvolume_; }

     protected:
      // const std::vector<Side*> & Sides() { return sides_; }

      /*! @name All these functions have to be implemented in the dervied classes!
       *
       *  \remark Please note, that none of these functions has any inherent
       *          dimension checks! Be careful if you access them directly. Each of
       *          these functions has a PUBLIC alternative which checks the dimensions
       *          and is therefore much safer.                      hiermeier 07/16 */
      //! @{
      /*! Returns element local coordinates "rst" of a point from its global coordinates "xyz"
       *
       *  \param xyz (in) : Global spatial coordinate
       *  \param rst (out): Local parameter space coordinate */
      virtual bool LocalCoordinates(const double* xyz, double* rst) = 0;

      /*! Returns element global coordinates "xyz" of a point from its local coordinates "rst"
       *
       *  \param rst (in): Local parameter space coordinate
       *  \param xyz (out) : Global spatial coordinate */
      virtual void GlobalCoordinates(const double* rst, double* xyz) = 0;

      /*! \brief Returns the element center global coordinates
       *
       *  \param midpoint (out): The midpoint of the element in global spatial coordinates */
      virtual void ElementCenter(double* midpoint) = 0;

      /*! \brief Find the scalar value at a particular point inside the element
       *  specified by its local coordinates rst.
       *
       *  \param ns (in)   : contains the value of the scalar at the corner nodes of the element
       *  \param rst (in)  : local parameter space coordinate inside the element */
      virtual double Scalar(const std::vector<double>& ns, const double* rst) = 0;

      /*! \brief Return the level set value at the given global coordinate which
       *   has to be INSIDE the element.
       *
       *  \param xyz (in) : Global spatial coordinates */
      virtual double GetLevelSetValue(const double* xyz) = 0;

      /*! \brief Return the level set value at the given local coordinate which
       *   has to be INSIDE the element.
       *
       *  \param rst (in) : local parameter space coordinates */
      virtual double GetLevelSetValueAtLocalCoords(const double* rst) = 0;

      /*! \brief Return the level set gradient in global coordinates
       *
       *  \param xyz (in) : global spatial coordinates */
      virtual std::vector<double> GetLevelSetGradient(const double* xyz) = 0;

      /*! \brief Return the level set gradient in global coordinates
       *
       *  \param rst (in) : local parameter space coordinates */
      virtual std::vector<double> GetLevelSetGradientAtLocalCoords(const double* rst) = 0;

      /*! \brief Return the level set gradient in local (parameter space) coordinates
       *
       *  \param xyz (in) : global spatial coordinates */
      virtual std::vector<double> GetLevelSetGradientInLocalCoords(const double* xyz) = 0;

      /*! \brief Return the level set gradient in local (parameter space) coordinates
       *
       *  \param rst (in) : local parameter space coordinates */
      virtual std::vector<double> GetLevelSetGradientAtLocalCoordsInLocalCoords(
          const double* rst) = 0;
      //! @}
     private:
      /*! \brief Returns true if there is at least one cut point between the background
       *  element side and the cut side ("other") */
      bool FindCutPoints(Mesh& mesh, Side& side, Side& other, int recursion);

      bool FindCutLines(Mesh& mesh, Side& side, Side& other);

      /// element Id of the linear element (sub-element for quadratic elements)
      int eid_;

      /*!\brief element Id of the original quadratic/linear element
       *
       * \remark for time integration and parallel cut it is necessary that
       * parent_id is set also for linear elements (then eid_=parent_id_) */
      int parent_id_;

      bool active_;

      std::vector<Side*> sides_;

      std::vector<Node*> nodes_;

      std::vector<Point*> points_;

      /** (sub-)sides that cut this element, also sides that just touch the element
       * at a single point, edge or the whole side
       * all sides for that the intersection with the element finds points */
      plain_side_set cut_faces_;

      plain_facet_set facets_;

      plain_volumecell_set cells_;

      /** For originally linear elements, this just stores the corner points. However,
       * if this linear element is derived from a Quadratic element, then this contains
       * the corner points of quad element. */
      std::vector<Node*> quadCorners_;

      /// True if this is a linear (shadow) element derived from a Quadratic element
      bool isShadow_;

      /// Discretization shape of the parent quad element
      DRT::Element::DiscretizationType quadshape_;

      /// the bounding volume of the element
      Teuchos::RCP<BoundingBox> boundingvolume_;

      /// type of integration-rule for this element
      INPAR::CUT::ElementIntegrationType eleinttype_;

    };  // class Element

    /*--------------------------------------------------------------------------*/
    /*! \brief Implementation of the concrete cut element
     *
     *  \author hiermeier */
    template <unsigned probdim, DRT::Element::DiscretizationType elementtype,
        unsigned numNodesElement =
            DRT::UTILS::DisTypeToNumNodePerEle<elementtype>::numNodePerElement,
        unsigned dim = DRT::UTILS::DisTypeToDim<elementtype>::dim>
    class ConcreteElement : public Element
    {
     public:
      /// constructor
      ConcreteElement(
          int eid, const std::vector<Side*>& sides, const std::vector<Node*>& nodes, bool active)
          : Element(eid, sides, nodes, active)
      {
      }

      //! destructor
      virtual ~ConcreteElement(){};

      //! return element shape
      virtual DRT::Element::DiscretizationType Shape() const { return elementtype; }

      //! get problem dimension
      virtual unsigned ProbDim() const { return probdim; }

      //! get element dimension
      virtual unsigned Dim() const { return dim; }

      //! get the number of nodes
      virtual unsigned NumNodes() const { return numNodesElement; }

      bool PointInside(Point* p);

      void Coordinates(LINALG::Matrix<probdim, numNodesElement>& xyze) const
      {
        Coordinates(xyze.A());
      }

      /*! \brief Get the coordinates of all the nodes of this linear shadow element */
      virtual void Coordinates(double* xyze) const
      {
        double* x = xyze;
        for (std::vector<Node*>::const_iterator i = Nodes().begin(); i != Nodes().end(); ++i)
        {
          Node& n = **i;
          n.Coordinates(x);
          x += probdim;
        }
      }

      /*! \brief Get the coordinates of all node of parent Quad element */
      virtual void CoordinatesQuad(double* xyze)
      {
        double* x = xyze;
        for (std::vector<Node*>::const_iterator i = QuadCorners().begin(); i != QuadCorners().end();
             ++i)
        {
          Node& n = **i;
          n.Coordinates(x);
          x += probdim;
        }
      }

      /*! \brief Returns element local coordinates "rst" of a point from its global coordinates
       * "xyz"
       *
       *  This variant uses the LINALG::Matrices as input and does a dimension check.
       *  Call this function, to be on the safe side.
       *
       *  \param xyz (in)        : Global spatial coordinate
       *  \param rst (out)       : Local parameter space coordinate */
      bool LocalCoordinates(const LINALG::Matrix<probdim, 1>& xyz, LINALG::Matrix<dim, 1>& rst);

      /*! Returns element global coordinates "xyz" of a point from its local coordinates "rst"
       *
       *  This variant uses the LINALG::Matrices as input and does a dimensional check.
       *  Call this function to be on the safe side.
       *
       *  \param rst (in)  : Local parameter space coordinate
       *  \param xyz (out) : Global spatial coordinate */
      void GlobalCoordinates(const LINALG::Matrix<dim, 1>& rst, LINALG::Matrix<probdim, 1>& xyz)
      {
        LINALG::Matrix<numNodesElement, 1> funct;
        DRT::UTILS::shape_function<elementtype>(rst, funct);

        xyz = 0;

        const std::vector<Node*>& nodes = Nodes();
        for (unsigned i = 0; i < numNodesElement; ++i)
        {
          LINALG::Matrix<probdim, 1> x(nodes[i]->point()->X());
          xyz.Update(funct(i), x, 1);
        }
      }

      /** \brief get the global coordinates in the element at given local coordinates
       *
       *  \param rst (in)  : local coordinates
       *  \param xyz (out) : corresponding global coordinates
       *
       *  \author hiermeier
       *  \date 08/16 */
      void PointAt(const LINALG::Matrix<dim, 1>& rst, LINALG::Matrix<probdim, 1>& xyz)
      {
        LINALG::Matrix<numNodesElement, 1> funct(true);
        LINALG::Matrix<probdim, numNodesElement> xyze(true);
        this->Coordinates(xyze);

        DRT::UTILS::shape_function<elementtype>(rst, funct);
        xyz.Multiply(xyze, funct);
      }

      /*! \brief Returns the element center global coordinates
       *
       *  This variant uses the LINALG::Matrix as input and does a dimensional check.
       *  Call this function to be on the safe side.
       *
       *  \param midpoint (out): The midpoint of the element in global spatial coordinates */
      void ElementCenter(LINALG::Matrix<probdim, 1>& midpoint)
      {
        // get the element center in the parameter coordinates (dim)
        LINALG::Matrix<dim, 1> center_rst(DRT::UTILS::getLocalCenterPosition<dim>(elementtype));
        PointAt(center_rst, midpoint);
      }

      /*! \brief Find the scalar value at a particular point inside the element
       *  specified by its local coordinates rst.
       *
       *  This variant uses the LINALG::Matrix as input and does a dimensional check.
       *  Call this function to be on the safe side.
       *
       *  \param ns (in)   : contains the value of the scalar at the corner nodes of the element
       *  \param rst (in)  : local parameter space coordinate inside the element */
      double Scalar(const std::vector<double>& ns, const LINALG::Matrix<dim, 1>& rst)
      {
        LINALG::Matrix<numNodesElement, 1> funct;
        DRT::UTILS::shape_function<elementtype>(rst, funct);

        LINALG::Matrix<numNodesElement, 1> scalar(&ns[0]);
        LINALG::Matrix<1, 1> res;
        res.MultiplyTN(funct, scalar);
        return res(0);
      }

      /*! \brief Return the level set value at the given global coordinate which
       *   has to be INSIDE the element.
       *
       *  \param xyz (in) : Global spatial coordinates */
      double GetLevelSetValue(const LINALG::Matrix<probdim, 1>& xyz)
      {
        LINALG::Matrix<dim, 1> rst;
        LocalCoordinates(xyz, rst);
        return GetLevelSetValueAtLocalCoords(rst);
      }

      /*! \brief Return the level set value at the given local coordinate which
       *   has to be INSIDE the element.
       *
       *  \param rst (in) : local parameter space coordinates */
      double GetLevelSetValueAtLocalCoords(const LINALG::Matrix<dim, 1>& rst)
      {
        LINALG::Matrix<numNodesElement, 1> funct;

        DRT::UTILS::shape_function<elementtype>(rst, funct);

        const std::vector<Node*> ele_node = this->Nodes();

        // Extract Level Set values from element.
        LINALG::Matrix<numNodesElement, 1> escaa;
        int mm = 0;
        for (std::vector<Node*>::const_iterator i = ele_node.begin(); i != ele_node.end(); i++)
        {
          Node* nod = *i;
          escaa(mm, 0) = nod->LSV();
          mm++;
        }

        return funct.Dot(escaa);
      }

      /*! \brief Return the level set gradient in global coordinates
       *
       *  \param xyz (in) : global spatial coordinates */
      std::vector<double> GetLevelSetGradient(const LINALG::Matrix<probdim, 1>& xyz)
      {
        LINALG::Matrix<dim, 1> rst;
        LocalCoordinates(xyz, rst);
        return GetLevelSetGradientAtLocalCoords(rst);
      }

      /*! \brief Return the level set gradient in global coordinates
       *
       *  \param rst (in) : local parameter space coordinates */
      std::vector<double> GetLevelSetGradientAtLocalCoords(const LINALG::Matrix<dim, 1>& rst)
      {
        // Calculate global derivatives
        //----------------------------------
        LINALG::Matrix<probdim, numNodesElement> deriv1;
        LINALG::Matrix<probdim, numNodesElement> xyze;
        Coordinates(xyze);
        // transposed jacobian dxyz/drst
        LINALG::Matrix<probdim, probdim> xjm;
        // inverse of transposed jacobian drst/dxyz
        LINALG::Matrix<probdim, probdim> xij;
        // nodal spatial derivatives dN_i/dr * dr/dx + dN_i/ds * ds/dx + ...
        LINALG::Matrix<probdim, numNodesElement> derxy;
        // only filled for manifolds
        LINALG::Matrix<dim, dim> metrictensor;
        LINALG::Matrix<probdim, 1> normalvec1;
        LINALG::Matrix<probdim, 1> normalvec2;

        EvalDerivsInParameterSpace<probdim, elementtype>(
            xyze, rst, deriv1, metrictensor, xjm, &xij, &normalvec1, &normalvec2, true);

        // compute global first derivates
        derxy.Multiply(xij, deriv1);
        //----------------------------------

        const std::vector<Node*> ele_node = this->Nodes();

        // Extract Level Set values from element.
        LINALG::Matrix<1, numNodesElement> escaa;
        int mm = 0;
        for (std::vector<Node*>::const_iterator i = ele_node.begin(); i != ele_node.end(); i++)
        {
          Node* nod = *i;
          escaa(0, mm) = nod->LSV();
          mm++;
        }
        LINALG::Matrix<probdim, 1> phi_deriv1;
        phi_deriv1.MultiplyNT(derxy, escaa);

        std::vector<double> normal_facet(phi_deriv1.A(), phi_deriv1.A() + probdim);
        if (normal_facet.size() != probdim) dserror("Something went wrong!");

        return normal_facet;
      }

      /*! \brief Return the level set gradient in local (parameter space) coordinates
       *
       *  \param xyz (in) : global spatial coordinates */
      std::vector<double> GetLevelSetGradientInLocalCoords(const LINALG::Matrix<probdim, 1>& xyz)
      {
        LINALG::Matrix<dim, 1> rst;
        LocalCoordinates(xyz, rst);
        return GetLevelSetGradientAtLocalCoordsInLocalCoords(rst);
      }

      /*! \brief Return the level set gradient in local (parameter space) coordinates
       *
       *  \param rst (in) : local parameter space coordinates */
      std::vector<double> GetLevelSetGradientAtLocalCoordsInLocalCoords(
          const LINALG::Matrix<dim, 1>& rst)
      {
        LINALG::Matrix<dim, numNodesElement> deriv1;

        DRT::UTILS::shape_function_deriv1<elementtype>(rst, deriv1);

        const std::vector<Node*> ele_node = this->Nodes();

        // Extract Level Set values from element.
        LINALG::Matrix<1, numNodesElement> escaa;
        int mm = 0;
        for (std::vector<Node*>::const_iterator i = ele_node.begin(); i != ele_node.end(); i++)
        {
          Node* nod = *i;
          escaa(0, mm) = nod->LSV();
          mm++;
        }
        LINALG::Matrix<dim, 1> phi_deriv1;
        phi_deriv1.MultiplyNT(deriv1, escaa);

        std::vector<double> normal_facet(phi_deriv1.A(), phi_deriv1.A() + dim);
        if (normal_facet.size() != dim) dserror("Something went wrong! (dim=%d)", dim);

        return normal_facet;
      }

     protected:
      //! derived
      virtual bool LocalCoordinates(const double* xyz, double* rst)
      {
        const LINALG::Matrix<probdim, 1> xyz_mat(xyz, true);  // create view
        LINALG::Matrix<dim, 1> rst_mat(rst, true);            // create view
        return LocalCoordinates(xyz_mat, rst_mat);
      }

      //! derived
      virtual void GlobalCoordinates(const double* rst, double* xyz)
      {
        const LINALG::Matrix<dim, 1> rst_mat(rst, true);  // create view
        LINALG::Matrix<probdim, 1> xyz_mat(xyz, true);    // create view
        GlobalCoordinates(rst_mat, xyz_mat);
      }

      //! derived
      virtual void ElementCenter(double* midpoint)
      {
        LINALG::Matrix<probdim, 1> midpoint_mat(midpoint, true);  // create view
        ElementCenter(midpoint_mat);
      }

      //! derived
      virtual double Scalar(const std::vector<double>& ns, const double* rst)
      {
        const LINALG::Matrix<dim, 1> rst_mat(rst, true);  // create view
        return Scalar(ns, rst_mat);
      }

      //! derived
      virtual double GetLevelSetValue(const double* xyz)
      {
        const LINALG::Matrix<probdim, 1> xyz_mat(xyz, true);  // create view
        return GetLevelSetValue(xyz_mat);
      }

      virtual double GetLevelSetValueAtLocalCoords(const double* rst)
      {
        const LINALG::Matrix<dim, 1> rst_mat(rst, true);  // create view
        return GetLevelSetValueAtLocalCoords(rst_mat);
      }

      //! derived
      virtual std::vector<double> GetLevelSetGradient(const double* xyz)
      {
        const LINALG::Matrix<probdim, 1> xyz_mat(xyz, true);  // create view
        return GetLevelSetGradient(xyz_mat);
      }

      //! derived
      virtual std::vector<double> GetLevelSetGradientAtLocalCoords(const double* rst)
      {
        const LINALG::Matrix<dim, 1> rst_mat(rst, true);  // create view
        return GetLevelSetGradientAtLocalCoords(rst_mat);
      }

      //! derived
      virtual std::vector<double> GetLevelSetGradientInLocalCoords(const double* xyz)
      {
        const LINALG::Matrix<probdim, 1> xyz_mat(xyz, true);  // create view
        return GetLevelSetGradientInLocalCoords(xyz_mat);
      }

      //! derived
      virtual std::vector<double> GetLevelSetGradientAtLocalCoordsInLocalCoords(const double* rst)
      {
        const LINALG::Matrix<dim, 1> rst_mat(rst, true);  // create view
        return GetLevelSetGradientAtLocalCoordsInLocalCoords(rst_mat);
      }
    };  // class ConcreteElement

    /*--------------------------------------------------------------------------*/
    class ElementFactory
    {
     public:
      /// constructor
      ElementFactory(){};

      Teuchos::RCP<Element> CreateElement(DRT::Element::DiscretizationType elementtype, int eid,
          const std::vector<Side*>& sides, const std::vector<Node*>& nodes, bool active) const;

     private:
      template <DRT::Element::DiscretizationType elementtype>
      Element* CreateConcreteElement(int eid, const std::vector<Side*>& sides,
          const std::vector<Node*>& nodes, bool active, int probdim) const
      {
        Element* e = NULL;
        // sanity check
        if (probdim < DRT::UTILS::DisTypeToDim<elementtype>::dim)
          dserror("Problem dimension is smaller than the element dimension!");

        switch (probdim)
        {
          case 2:
            e = new ConcreteElement<2, elementtype>(eid, sides, nodes, active);
            break;
          case 3:
            e = new ConcreteElement<3, elementtype>(eid, sides, nodes, active);
            break;
          default:
            dserror("Unsupported problem dimension! ( probdim = %d )", probdim);
            break;
        }
        return e;
      }
    };

    // typedef EntityIdLess<Element> ElementIdLess;
    //
    // inline int EntityId( const Element & e ) { return e.Id(); }

  }  // namespace CUT
}  // namespace GEO

#endif

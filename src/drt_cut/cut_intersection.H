/*----------------------------------------------------------------------*/
/*!
\file cut_intersection.H

\brief here the intersection of a (plane) surface with a line is performed

\level 2

\maintainer Christoph Ager
*/
/*----------------------------------------------------------------------*/

#ifndef CUT_INTERSECTION_H
#define CUT_INTERSECTION_H

#include "cut_mesh.H"
#include "cut_kernel.H"
#include "cut_output.H"
#include "cut_node.H"
#include "cut_edge.H"
#include "cut_side.H"
#include "cut_options.H"

#include <Teuchos_TimeMonitor.hpp>


namespace GEO
{
  namespace CUT
  {
    class Side;
    class Edge;

    enum IntersectionStatus
    {
      intersect_newton_failed = -2,      // if newton failed
      intersect_unevaluated = -1,        // before ComputeEdgeSideIntersection has been called
      intersect_no_cut_point = 0,        // no cut point was found
      intersect_single_cut_point = 1,    // one single cut point was found
      intersect_multiple_cut_points = 2  // parallel cases
    };


    enum ParallelIntersectionStatus
    {
      intersection_not_possible = -1,
      intersection_not_found = 0,
      intersection_found = 1,
    };

    //! Map IntersectionStatus to std::string
    static inline std::string IntersectionStatus2String(IntersectionStatus istatus)
    {
      switch (istatus)
      {
        case intersect_unevaluated:
          return "intersect_unevaluated";
        case intersect_no_cut_point:
          return "intersect_no_cut_point";
        case intersect_single_cut_point:
          return "intersect_single_cut_point";
        case intersect_multiple_cut_points:
          return "intersect_multiple_cut_points";
        case intersect_newton_failed:
          return "intersect_newton_failed";
        default:
          return "Unknown IntersectionStatus";
      }
      exit(EXIT_FAILURE);
    };

    inline IntersectionStatus IntersectionStatus2Enum(unsigned num_cut_points)
    {
      switch (num_cut_points)
      {
        case 0:
          return intersect_no_cut_point;
        case 1:
          return intersect_single_cut_point;
        default:
          return intersect_multiple_cut_points;
      }
      exit(EXIT_FAILURE);
    }

    /*--------------------------------------------------------------------------*/
    /** \brief Base class to calculate the intersection of an edge with a side.
     *
     *  \author ager, hiermeier*/
    class IntersectionBase
    {
     public:
      static Teuchos::RCP<IntersectionBase> Create(
          const DRT::Element::DiscretizationType& edge_type,
          const DRT::Element::DiscretizationType& side_type);

     public:
      /// constructor
      IntersectionBase()
          : isinit_(false),
            isscaled_(false),
            isshifted_(false),
            useboundingbox_(false),
            mesh_ptr_(NULL),
            edge_ptr_(NULL),
            side_ptr_(NULL),
            options_ptr_(NULL)
      {
      }

      /** Lean Init() routine w/o mesh, edge or side objects
       *
       *  \remark If you use this Init() routine, you won't be able to call
       *  the Intersect() routine. Simply due to the fact, that you haven't
       *  passed the necessary input objects. Use the 2-nd (standard) Init()
       *  routine, instead. Anyhow, this Init() routine is the right one,
       *  if you want to intersect two edges. Just use the routine
       *  ComputeEdgeSideIntersection() afterwards.
       *
       *  \param xyze_lineElement    (in) : global nodal coordinates of the edge element
       *  \param xyze_surfaceElement (in) : global nodal coordinates of the side element
       *  \param usescaling          (in) : switch scaling on/off
       *  \param useshifting         (in) : switch shifting on/off
       *  \param useboundingbox      (in) : switch the bounding box checks on/off
       *
       *  \author hiermeier \date 08/16 */
      template <class T1, class T2>
      void Init(T1& xyze_lineElement, T2& xyze_surfaceElement, bool usescaling, bool useshifting,
          bool useboundingbox, Options* options)
      {
        isscaled_ = usescaling;
        isshifted_ = useshifting;
        useboundingbox_ = useboundingbox;

        mesh_ptr_ = NULL;
        edge_ptr_ = NULL;
        side_ptr_ = NULL;
        options_ptr_ = options;

        if (static_cast<unsigned>(xyze_lineElement.M()) != ProbDim() or
            static_cast<unsigned>(xyze_lineElement.N()) != NumNodesEdge())
          dserror(
              "Dimension mismatch of xyze_lineElement! \n"
              "expected input: %d x %d (rows x cols)\n"
              "current input : %d x %d (rows x cols)",
              ProbDim(), NumNodesEdge(), xyze_lineElement.M(), xyze_lineElement.N());

        if (static_cast<unsigned>(xyze_surfaceElement.M()) != ProbDim() or
            static_cast<unsigned>(xyze_surfaceElement.N()) != NumNodesSide())
          dserror(
              "Dimension mismatch of xyze_surfaceElement! \n"
              "expected input: %d x %d (rows x cols)\n"
              "current input : %d x %d (rows x cols)",
              ProbDim(), NumNodesSide(), xyze_surfaceElement.M(), xyze_surfaceElement.N());

        SetCoordinates(xyze_surfaceElement.A(), xyze_lineElement.A());
        ScaleAndShift();

        isinit_ = true;
      }

      /** \brief Standard Init() routine
       *
       *  \param mesh_ptr       (in) : pointer to the underlying mesh
       *  \param edge_ptr       (in) : pointer to the intersecting edge object
       *  \param side_ptr       (in) : pointer to the side which will be intersected
       *  \param usescaling     (in) : switch scaling on/off
       *  \param useshifting    (in) : switch shifting on/off
       *  \param useboundingbox (in) : switch the bounding box checks on/off
       *
       *  \author hiermeier \date 08/16 */
      void Init(Mesh* mesh_ptr, Edge* edge_ptr, Side* side_ptr, bool usescaling, bool useshifting,
          bool useboundingbox)
      {
        isscaled_ = usescaling;
        isshifted_ = useshifting;
        useboundingbox_ = useboundingbox;

        mesh_ptr_ = mesh_ptr;
        edge_ptr_ = edge_ptr;
        side_ptr_ = side_ptr;
        options_ptr_ = &(mesh_ptr->CreateOptions());

        SetCoordinates();
        ScaleAndShift();

        isinit_ = true;
      }

      /// destructor
      virtual ~IntersectionBase(){};

      /** \brief Calculate the actual intersection of an edge and a side ( or 2-nd edge )
       *
       *  See derived class for more information.
       *
       *  \author hiermeier \date 08/16 */
      virtual IntersectionStatus ComputeEdgeSideIntersection(
          double& tolerance, bool check_inside = true, std::vector<int>* touched_edges = NULL) = 0;

      /** \brief Computes the intersection points of the edge with the specified side
       *  and stores the points in cuts
       *
       *  See derived class for more information.
       *
       *  \author hiermeier \date 08/16 */
      virtual bool Intersect(PointSet& cuts) = 0;

      virtual ParallelIntersectionStatus HandleParallelIntersection(
          PointSet& cuts, int id = -1, bool output = false) = 0;

      virtual bool TriangulatedIntersection(PointSet& cuts) = 0;

      virtual bool HandleSpecialCases() = 0;

      /** \brief Get the final cut point global coordinates
       *
       *  Only allowed if there was only one cut point!
       *
       *  \author hiermeier \date 08/16 */
      virtual double* FinalPoint() = 0;

      virtual double* FinalPoint(unsigned cp_id) = 0;

      /** Access the cut point local coordinates on the side element
       * ( also working for multiple cut points )
       *
       *  \author hiermeier \date 01/17 */
      template <unsigned dimside>
      void LocalSideCoordinates(std::vector<LINALG::Matrix<dimside, 1>>& side_rs_cuts)
      {
        if (GetIntersectionStatus() < intersect_single_cut_point)
          dserror("INVALID IntersectionStatus! ( istatus = \"%s\" )",
              IntersectionStatus2String(GetIntersectionStatus()).c_str());

        side_rs_cuts.clear();
        side_rs_cuts.reserve(NumCutPoints());

        for (unsigned i = 0; i < NumCutPoints(); ++i)
          side_rs_cuts.push_back(LINALG::Matrix<dimside, 1>(LocalSideCoordinates(i), true));
      }

      /** Access the final cut point global coordinates
       * ( also working for multiple cut points )
       *
       *  \author hiermeier \date 01/17 */
      template <unsigned probdim>
      void FinalPoints(std::vector<LINALG::Matrix<probdim, 1>>& xyz_cuts)
      {
        if (GetIntersectionStatus() < intersect_single_cut_point)
          dserror("INVALID IntersectionStatus! ( istatus = \"%s\" )",
              IntersectionStatus2String(GetIntersectionStatus()).c_str());

        xyz_cuts.clear();
        xyz_cuts.reserve(NumCutPoints());

        for (unsigned i = 0; i < NumCutPoints(); ++i)
          xyz_cuts.push_back(LINALG::Matrix<probdim, 1>(FinalPoint(i), false));
      }

      virtual double* LocalCoordinates() = 0;

      virtual double* LocalSideCoordinates(unsigned cp_id) = 0;

      virtual bool SurfaceWithinLimits(double tol = REFERENCETOL) const = 0;

      virtual bool LineWithinLimits(double tol = REFERENCETOL) const = 0;

     protected:
      virtual unsigned NumCutPoints() const = 0;

      virtual const IntersectionStatus& GetIntersectionStatus() const = 0;

      inline void CheckInit() const
      {
        if (not isinit_) dserror("The Intersection object is not initialized! Call Init() first.");
      }

      virtual unsigned ProbDim() const = 0;
      virtual unsigned NumNodesSide() const = 0;
      virtual unsigned NumNodesEdge() const = 0;

      virtual void SetCoordinates() = 0;
      virtual void SetCoordinates(double* xyze_surfaceElement, double* xyze_lineElement) = 0;

      virtual void ScaleAndShift() = 0;

      /// Are the global coordinates scaled?
      const bool& IsScaled() const { return isscaled_; }

      /// Are the global coordinates shifted?
      const bool& IsShifted() const { return isshifted_; }

      /// Shall we use the bounding box information?
      const bool& UseBoundingBox() const { return useboundingbox_; }

      /// get a reference to the mesh object
      Mesh& GetMesh()
      {
        if (mesh_ptr_ != NULL) return *mesh_ptr_;
        dserror("The mesh pointer is not yet initialized!");
        exit(EXIT_FAILURE);
      }

      /// get a pointer to the mesh object
      Mesh* GetMeshPtr()
      {
        if (mesh_ptr_ != NULL) return mesh_ptr_;
        dserror("The mesh pointer is not yet initialized!");
        exit(EXIT_FAILURE);
      }

      /// get a reference to the edge object
      Edge& GetEdge()
      {
        if (edge_ptr_ != NULL) return *edge_ptr_;
        dserror("The edge pointer is not yet initialized!");
        exit(EXIT_FAILURE);
      }

      /// get a pointer to the edge object
      Edge* GetEdgePtr()
      {
        if (edge_ptr_ != NULL) return edge_ptr_;
        dserror("The edge pointer is not yet initialized!");
        exit(EXIT_FAILURE);
      }

      /// get a reference to the side object
      Side& GetSide()
      {
        if (side_ptr_ != NULL) return *side_ptr_;
        dserror("The side pointer is not yet initialized!");
        exit(EXIT_FAILURE);
      }

      /// get a pointer to the side object
      Side* GetSidePtr()
      {
        if (side_ptr_ != NULL) return side_ptr_;
        dserror("The side pointer is not yet initialized!");
        exit(EXIT_FAILURE);
      }

      /// get a pointer to the cut options object
      Options* GetOptionsPtr()
      {
        if (options_ptr_ != NULL) return options_ptr_;
        dserror("The option pointer is not yet initialized!");
        exit(EXIT_FAILURE);
      }

     private:
      /// flag which indicates whether the Init() has been called or not.
      bool isinit_;

      /// Did we scale the position vectors?
      bool isscaled_;

      /// Did we shift the position vectors?
      bool isshifted_;

      /// Shall we use bounding boxes to speed things up?
      bool useboundingbox_;

      /// mesh pointer
      Mesh* mesh_ptr_;

      /// edge pointer
      Edge* edge_ptr_;

      /// side pointer
      Side* side_ptr_;

      /// cut options pointer
      Options* options_ptr_;

    };  // class IntersectionBase


    /*--------------------------------------------------------------------------*/
    /** \brief Concrete class to calculate the intersection of an edge with a side.
     *
     *  The core class where all the cut points are actually calculated. It is
     *  also meaningful to use this class to calculate the intersection of two
     *  edges, if the related Init() routine is used.
     *
     *  \author ager, hiermeier */
    template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
        DRT::Element::DiscretizationType sidetype, bool debug = false,
        unsigned dimedge = DRT::UTILS::DisTypeToDim<edgetype>::dim,
        unsigned dimside = DRT::UTILS::DisTypeToDim<sidetype>::dim,
        unsigned numNodesEdge = DRT::UTILS::DisTypeToNumNodePerEle<edgetype>::numNodePerElement,
        unsigned numNodesSide = DRT::UTILS::DisTypeToNumNodePerEle<sidetype>::numNodePerElement>
    class Intersection : public IntersectionBase
    {
     public:
      /// constructor
      Intersection()
          : IntersectionBase(),
            xsi_side_(xsi_.A(), true),
            xsi_edge_(xsi_.A() + dimside, true),
            multiple_xsi_side_(0),
            multiple_xsi_edge_(0),
            num_cut_points_(0),
            istatus_(intersect_unevaluated),
            scale_(1.0),
            shift_(0.0)
      {
        /* intentionally left blank */
      }

      /// destructor
      virtual ~Intersection(){};

      // No public access to these methods! Use the base class accessors, instead.
     protected:
      /// get the number of detected feasible cut points
      virtual unsigned NumCutPoints() const
      {
        if (num_cut_points_ > 1 and (multiple_xsi_edge_.size() != num_cut_points_ or
                                        multiple_xsi_side_.size() != num_cut_points_))
          dserror("Size mismatch!");

        return num_cut_points_;
      }

      /// get the intersection status
      virtual const IntersectionStatus& GetIntersectionStatus() const { return istatus_; }

      /// get the local cut coordinates
      virtual double* LocalCoordinates() { return xsi_.A(); }

      /** \brief access the local coordinates of the cut point corresponding to the
       *  cut point ID \c cp_id on the side element
       *
       *  \param cp_id (in) : cut point id
       *
       *  \author hiermeier \date 01/17 */
      virtual double* LocalSideCoordinates(unsigned cp_id)
      {
        if (NumCutPoints() < 2) return xsi_side_.A();

        return multiple_xsi_side_[cp_id].A();
      }

      /** \brief access the local coordinates of the cut point corresponding to the
       *  cut point ID \c cp_id on the edge element
       *
       *  \param cp_id (in) : cut point id
       *
       *  \author hiermeier \date 01/17 */
      const LINALG::Matrix<dimedge, 1>& LocalEdgeCoordinates(const unsigned& cp_id) const
      {
        if (NumCutPoints() < 2) return xsi_edge_;

        return multiple_xsi_edge_[cp_id];
      }

      /// get the final cut point global coordinates
      void FinalPoint(const LINALG::Matrix<dimedge, 1>& xsi_edge, LINALG::Matrix<probdim, 1>& x)
      {
        CheckInit();

        // get final point
        x = 0;
        LINALG::Matrix<numNodesEdge, 1> lineFunct;
        DRT::UTILS::shape_function<edgetype>(xsi_edge, lineFunct);
        for (unsigned inode = 0; inode < numNodesEdge; ++inode)
          for (unsigned isd = 0; isd < probdim; ++isd)
            x(isd) += xyze_lineElement_(isd, inode) * lineFunct(inode);

        // first un-do the shifting
        x.Update(1.0, shift_, 1.0);
        // second un-do the scaling
        x.Scale(scale_);
      }
      virtual double* FinalPoint()
      {
        if (istatus_ != intersect_single_cut_point)
          dserror(
              "INVALID IntersectionStatus: This routine is restricted to one single "
              "cut point only! ( istatus_ = \"%s\" )",
              IntersectionStatus2String(istatus_).c_str());

        FinalPoint(xsi_edge_, x_);
        return x_.A();
      }
      virtual double* FinalPoint(unsigned cp_id)
      {
        FinalPoint(LocalEdgeCoordinates(cp_id), x_);
        return x_.A();
      }

      // Remove all the edges from list of touching edges that are further away than 1e-14 (
      // TOPOLOGICAL_TOLERANCE ) distance from the point
      template <unsigned int dim, INPAR::CUT::CUT_Floattype floattype>
      void FixDistantTouchingEdges(
          LINALG::Matrix<dim, 1>& p_coord, std::vector<int>& touching_edges)
      {
        bool signeddistance = false;
        double distance = 0;
        LINALG::Matrix<dim, 1> xsi;
        LINALG::Matrix<dim, numNodesEdge> xyze_edge;
        const std::vector<Edge*>& side_edges = GetSide().Edges();

        for (std::vector<int>::iterator it = touching_edges.begin(); it != touching_edges.end();)
        {
          side_edges[*it]->Coordinates(xyze_edge.A());

          KERNEL::ComputeDistance<dim, edgetype, floattype> cd(xsi);
          bool conv = cd(xyze_edge, p_coord, distance, signeddistance);
          KERNEL::PointOnSurfaceLoc loc = cd.GetSideLocation();

          if (conv)
          {
            if (not loc.OnSide())
            {
              // safety check if it is larger then some (rather arbitrary distance)
              if (distance > 1e-10)
              {
                std::ofstream file("far_touched_edges.pos");
                Edge* e = side_edges[*it];
                GEO::CUT::OUTPUT::GmshEdgeDump(file, e, std::string("FarEdge"));
                GEO::CUT::OUTPUT::GmshNewSection(file, "Point", false);
                LINALG::Matrix<3, 1> p(p_coord.A());
                GEO::CUT::OUTPUT::GmshCoordDump(file, p, 0);
                GEO::CUT::OUTPUT::GmshEndSection(file);
                file.close();
                GenerateGmshDump();
                dserror("Distance between point touching edge is too high! Check this case!");
              }
              it = touching_edges.erase(it);
            }
            else
              ++it;
          }
          else
            dserror("Newton did not converge for simple ComputeDistance between point and a line");
        }
      }

      /** \brief Calculate the actual intersection of an edge and a side ( or 2-nd edge )
       *
       *  \param tolerance (out) : tolerance specified by the CUT::KERNEL
       *                           intersection method
       *
       *  This function returns GEO::CUT::IntersectionStatus. There are three different
       *  outcomes:
       *
       *  ( 1 ) Multiple cut points are detected during the CheckParallelism call.
       *
       *  ( 2 ) A single cut point is detected during the CheckParallelism or the
       *        intersection call.
       *
       *  ( 3 ) There is no feasible cut point.
       *
       *  All feasible cut points are within the given element bounds.
       *
       *  \author hiermeier \date 08/16 */
      virtual IntersectionStatus ComputeEdgeSideIntersection(
          double& tolerance, bool check_inside = true, std::vector<int>* touched_edges = NULL)
      {
        switch (GetOptionsPtr()->GeomIntersect_Floattype())
        {
          case INPAR::CUT::floattype_cln:
          {
            return ComputeEdgeSideIntersectionT<INPAR::CUT::floattype_cln>(
                tolerance, check_inside, touched_edges);
          }
          case INPAR::CUT::floattype_double:
          {
            return ComputeEdgeSideIntersectionT<INPAR::CUT::floattype_double>(
                tolerance, check_inside, touched_edges);
          }
          default:
            throw std::runtime_error("Unexpected floattype for ComputeEdgeSideIntersectionT!");
        }
      }

      template <INPAR::CUT::CUT_Floattype floattype>
      IntersectionStatus ComputeEdgeSideIntersectionT(
          double& tolerance, bool check_inside = true, std::vector<int>* touched_edges = NULL)
      {
        CheckInit();
        TEUCHOS_FUNC_TIME_MONITOR("ComputeEdgeSideIntersection");

        const bool success = CheckParallelism(multiple_xsi_side_, multiple_xsi_edge_, tolerance);

        /* The parallelism check was successful and we are done. At this point
         * it is possible, that we find more than one cut point. A special
         * treatment becomes necessary for multiple cut points. */
        if (success)
        {
          istatus_ = GEO::CUT::IntersectionStatus2Enum(num_cut_points_);
          return (istatus_);
        }

        KERNEL::ComputeIntersection<probdim, edgetype, sidetype,
            (floattype == INPAR::CUT::floattype_cln)>
            ci(xsi_);
        // KERNEL::DebugComputeIntersection<probdim,edgetype,sidetype,(floattype ==
        // INPAR::CUT::floattype_cln)> ci( xsi_ );

        bool conv = ci(xyze_surfaceElement_, xyze_lineElement_);
        tolerance = ci.GetTolerance();


        if (probdim > dimside + dimedge)
        {
          // edge-edge intersection, we might create point even if newton did not converge
          //
          double line_distance = ci.DistanceBetween();
          if ((line_distance < SIDE_DETECTION_TOLERANCE) and (ci.GetEdgeLocation().WithinSide()) and
              (ci.GetSideLocation().WithinSide()))
          {
            istatus_ = intersect_single_cut_point;
          }
          else
          {
            istatus_ = intersect_no_cut_point;
          }
        }


        // normal intersection
        else
        {
          /* Check if the found point is within the limits of the side and edge element,
           * if the Newton scheme did converge */
          if (check_inside)
          {
            if (conv)
            {
              if ((ci.GetEdgeLocation().WithinSide()) and (ci.GetSideLocation().WithinSide()))
                istatus_ = intersect_single_cut_point;
              // converged but is outside
              else
                istatus_ = intersect_no_cut_point;  // limits will be checked later
            }
            else
              istatus_ = intersect_newton_failed;
          }
          else
          {
            num_cut_points_ = conv;
            istatus_ = GEO::CUT::IntersectionStatus2Enum(num_cut_points_);
          }
          // if we want to get edges back
          if (touched_edges)
          {
            std::vector<int>& touched = *touched_edges;
            ci.GetTouchedSideEdges(touched);

            // this should not happen, as all the touching edges must be identified by edge-edge
            // intersections
            if ((istatus_ == intersect_no_cut_point) && (touched.size() > 0))
            {
              std::stringstream err_msg;
              err_msg << "Touching " << touched.size()
                      << " edges, but no intersection! This should not happen! ";
              GenerateGmshDump();
              dserror(err_msg.str());
            }
          }
        }

        return istatus_;
      }

      /** \brief Computes the intersection points of the edge with the specified side
       *   and stores the points in cuts
       *
       *  WARNING: Intersection just works for planes ( TRI3, QUAD4 unwarped! ) with lines!!!
       *
       *  (1) try to find not overlapping geometries with bounding-boxes to avoid a big load
       *  of work ... has to be done. This is here just for performance ... intersection
       *  should also be robust without that!!!
       *
       *  (2) first we start to calculate the distance with both end points of a line,
       *  go get rid of parallel cases (where intersection wouldn't converge) and also get
       *  rid of cases, where the line is just on one side of the surface --> definitely
       *  no intersection!!!
       *
       *  \remark As for QUAD4 where the projected end points of the line are outside the
       *  element, we do always get reliable results (normal can flip outside the element),
       *  generally the distance is computed to the two triangles (be aware of the fact,
       *  that this is just possible because we limit this function to plane ( unwarped )
       *  QUAD4 sides!!!
       *
       *  (3) Perform edge-edge intersection of line surface edges, bounding box is also
       *  applied here to sppeed up calculations
       *
       *  (4) try to calculate the intersection point directly with the Newton
       *  this will basically fail if the system is conditioned badly --> means that line
       *  and plane are parallel ( which shouldn't be the case anymore as it was already
       *  captured in point (2) ) or element is distorted or it's a QUAD4 and the
       *  intersection point is outside the element and is not part of the interpolation
       *  space! These cases should be treated separately later!!!
       *  If TRIANGULATED_INTERSECTION flag is enabled, intersection of the quad4 with
       *  the line is splitted into intersection of line with two tri3, obtained from the
       *  quad4. In the case it should always converge
       *
       *  (5) throw dserror in case this intersection wasn't treated right --> this means
       *  there is still handling of some special cases missing in the code & it does not
       *  mean that there is no intersection point.
       *
       *  \author ager */
      virtual bool Intersect(PointSet& cuts);

      // Try to find possible intersection points, if this intersection is between parallell size
      // and edge without using real ComputeIntersection
      virtual ParallelIntersectionStatus HandleParallelIntersection(
          PointSet& cuts, int id = -1, bool output = false);

      virtual void GenerateGmshDump();

      // Handle cases for which normal intersection procedure did not work
      virtual bool HandleSpecialCases();

      // compute intersection by splitting quad4 into two triangles
      virtual bool TriangulatedIntersection(PointSet& cuts);

      /** \brief Will return TRUE, if local side coordinates are within the side
       *  element parameter space bounds */
      virtual bool SurfaceWithinLimits(double tol = REFERENCETOL) const
      {
        return GEO::CUT::KERNEL::WithinLimits<sidetype>(xsi_side_, tol);
      }

      /** \brief Will return TRUE, if local side coordinates are within the TRI3
       *  side element parameter space bounds */
      bool Tri3WithinLimits(double tol = REFERENCETOL) const
      {
        return GEO::CUT::KERNEL::WithinLimits<DRT::Element::tri3>(xsi_side_, tol);
      }

      /** \brief Will return TRUE, if local edge coordinate is within the line
       *  element parameter space bounds */
      virtual bool LineWithinLimits(double tol = REFERENCETOL) const
      {
        return GEO::CUT::KERNEL::WithinLimits<edgetype>(xsi_edge_, tol);
      }

      /// access the problem dimension
      virtual unsigned ProbDim() const { return probdim; }

      /// access the number of nodes of the side ( or 2-nd edge ) element
      virtual unsigned NumNodesSide() const { return numNodesSide; }

      /// access the number of nodes of the edge
      virtual unsigned NumNodesEdge() const { return numNodesEdge; }

      /// set the edge and side coordinates
      virtual void SetCoordinates();

      /// set the edge and side coordinates
      virtual void SetCoordinates(double* xyze_surfaceElement, double* xyze_lineElement)
      {
        xyze_lineElement_.SetCopy(xyze_lineElement);
        xyze_surfaceElement_.SetCopy(xyze_surfaceElement);
      }

      /** \brief scale and shift the nodal positions of the given line and surface element
       *
       *  This can help to get a better conditioned system of equations and
       *  makes the used tolerances more reliable. The same procedure is used for
       *  the position calculation.
       *
       *  \author hiermeier
       *  \date 08/16 */
      virtual void ScaleAndShift()
      {
        // ---------------------------------------------------------------------
        // scale the input elements if desired
        // ---------------------------------------------------------------------
        if (not IsScaled())
          scale_ = 1.0;
        else
        {
          GetElementScale<probdim>(xyze_surfaceElement_, scale_);

          xyze_lineElement_.Scale(1. / scale_);
          xyze_surfaceElement_.Scale(1. / scale_);
        }
        // ---------------------------------------------------------------------
        // shift the input elements if desired
        // ---------------------------------------------------------------------
        if (not IsShifted())
          shift_ = 0;
        else
        {
          GetElementShift<probdim>(xyze_surfaceElement_, shift_);

          for (unsigned i = 0; i < numNodesSide; ++i)
          {
            LINALG::Matrix<probdim, 1> x1(&xyze_surfaceElement_(0, i), true);
            x1.Update(-1, shift_, 1);
          }
          for (unsigned i = 0; i < numNodesEdge; ++i)
          {
            LINALG::Matrix<probdim, 1> x1(&xyze_lineElement_(0, i), true);
            x1.Update(-1, shift_, 1);
          }
        }
      }

      /** check if the given local coordinates are at one of the edges of the side element,
       *  i.e. at the boundaries of the side element. */
      template <class T>
      bool AtEdge(const T& xsi)
      {
        return GEO::CUT::KERNEL::AtEdge<sidetype>(xsi);
      }

     private:
      /// Do the bounding box overlap check for the class internal edge and side variables
      bool CheckBoundingBoxOverlap();

      /// Check if the edge \c ebb and surface \c sbb bounding boxes overlap
      bool CheckBoundingBoxOverlap(BoundingBox& ebb, BoundingBox& sbb) const;

      /** \brief Checks the side dimension and calls the corresponding method
       *
       *  Currently surface and line elements are supported.
       *
       *  \author hiermeier \date 12/16 */
      bool CheckParallelism(std::vector<LINALG::Matrix<dimside, 1>>& side_rs_intersect,
          std::vector<LINALG::Matrix<dimedge, 1>>& edge_r_intersect, double& tolerance);

      /** \brief Check if the two lines are collinear, end points are on the line, or
       *  the distance values imply that no intersection is possible
       *
       *  \author hiermeier \date 12/16 */
      bool CheckCollinearity(std::vector<LINALG::Matrix<dimside, 1>>& side_rs_corner_intersect,
          std::vector<LINALG::Matrix<dimedge, 1>>& edge_r_corner_intersect, double& tolerance);

      /** \brief Check the angle between two edge lines.
       *
       *  This is a quick check to skip cases which are definitely not parallel.
       *
       *  \author hiermeier \date 12/16 */
      bool CheckAngleCriterionBetweenTwoEdges();

      /** ToDo This method is currently unused, since this case should be treated by
       *  the Intersect() method.
       *
       *  \author hiermeier \date 12/16 */
      bool CheckParallelismBetweenSideAndEdge(
          std::vector<LINALG::Matrix<dimside, 1>>& side_rs_intersect,
          std::vector<LINALG::Matrix<dimedge, 1>>& edge_r_intersect, double& tolerance);

      /** \brief Check the angle between a edge and a surface normal.
       *
       *  This is a quick check to skip cases which are definitely not parallel.
       *
       *  \author hiermeier \date 12/16 */
      bool CheckAngleCriterionBetweenSideNormalAndEdge();

      /// find the local coordinate of a given edge end point ( i.e. -1 or +1 )
      bool FindLocalCoordinateOfEdgeEndPoint(
          double& pos, const LINALG::Matrix<probdim, 1>& xyz, const double& tol) const;

      /** \brief Compute the intersection of an edge with a TRI3 surface element,
       *         which is created by splitting a QUAD4 element into two TRI3 elements
       *
       *  \param tolerance (out) : location status of the compute distance routine
       *  \param triangleid (in) : ID of the desired triangle ( 0 or 1 )
       *
       *  Returns TRUE if the calculation was successful. This does not imply, that the
       *  calculated intersection point is feasible!
       *
       *  \author hiermeier \date 08/16 */
      bool ComputeEdgeTri3Intersection(int triangleid, KERNEL::PointOnSurfaceLoc& location)
      {
        switch (GetOptionsPtr()->GeomIntersect_Floattype())
        {
          case INPAR::CUT::floattype_cln:
          {
            return ComputeEdgeTri3IntersectionT<INPAR::CUT::floattype_cln>(triangleid, location);
          }
          case INPAR::CUT::floattype_double:
          {
            return ComputeEdgeTri3IntersectionT<INPAR::CUT::floattype_double>(triangleid, location);
          }
          default:
            throw std::runtime_error("Unexpected floattype for ComputeEdgeTri3IntersectionT!");
        }
      }

      template <INPAR::CUT::CUT_Floattype floattype>
      bool ComputeEdgeTri3IntersectionT(int triangleid, KERNEL::PointOnSurfaceLoc& location)
      {
        if (triangleid < 0) dserror("The triangle id has to be positive!");

        TEUCHOS_FUNC_TIME_MONITOR("ComputeEdgeTri3Intersection");
        LINALG::Matrix<3, 1> xsi;
        if (xsi_.M() != 3)
          dserror("xsi_ has the wrong dimension! (dimedge + 2 = %d + 2)", dimedge);
        else
          xsi.SetView(xsi_.A());

        KERNEL::ComputeIntersection<3, edgetype, DRT::Element::tri3, floattype> ci(xsi);
        // KERNEL::DebugComputeIntersection<probdim,edgetype,DRT::Element::tri3,floattype> ci( xsi
        // );

        LINALG::Matrix<3, 3> xyze_triElement;
        GetTriangle(xyze_triElement, triangleid);
        LINALG::Matrix<3, numNodesEdge> xyze_lineElement(xyze_lineElement_.A(), true);

        bool conv = ci(xyze_triElement, xyze_lineElement);
        location = ci.GetSideLocation();

        return conv;
      }

      // Computes tri3 edge intersection used for quad4 -> 2 tri3 splits
      IntersectionStatus ComputeEdgeTri3IntersectionQuad4Split(
          int triangleid, bool* close_to_shared_edge = NULL)
      {
        switch (GetOptionsPtr()->GeomIntersect_Floattype())
        {
          case INPAR::CUT::floattype_cln:
          {
            return ComputeEdgeTri3IntersectionQuad4SplitT<INPAR::CUT::floattype_cln>(
                triangleid, close_to_shared_edge);
          }
          case INPAR::CUT::floattype_double:
          {
            return ComputeEdgeTri3IntersectionQuad4SplitT<INPAR::CUT::floattype_double>(
                triangleid, close_to_shared_edge);
          }
          default:
            throw std::runtime_error(
                "Unexpected floattype for ComputeEdgeTri3IntersectionQuad4SplitT!");
        }
      }

      template <INPAR::CUT::CUT_Floattype floattype>
      IntersectionStatus ComputeEdgeTri3IntersectionQuad4SplitT(
          int triangleid, bool* close_to_shared_edge = NULL)
      {
        if (triangleid < 0) dserror("The triangle id has to be positive!");

        TEUCHOS_FUNC_TIME_MONITOR("ComputeEdgeTri3Intersection");
        LINALG::Matrix<3, 1> xsi;
        if (xsi_.M() != 3)
          dserror("xsi_ has the wrong dimension! (dimedge + 2 = %d + 2)", dimedge);
        else
          xsi.SetView(xsi_.A());

        KERNEL::ComputeIntersection<3, edgetype, DRT::Element::tri3, floattype> ci(xsi);
        // KERNEL::DebugComputeIntersection<probdim,edgetype,DRT::Element::tri3,floattype> ci( xsi
        // );

        LINALG::Matrix<3, 3> xyze_triElement;
        GetTriangle(xyze_triElement, triangleid);
        LINALG::Matrix<3, numNodesEdge> xyze_lineElement(xyze_lineElement_.A(), true);

        bool conv = ci(xyze_triElement, xyze_lineElement);

        if (conv)
        {
          if (ci.GetEdgeLocation().WithinSide() and ci.GetSideLocation().WithinSide())
            istatus_ = intersect_single_cut_point;
          else
            istatus_ = intersect_no_cut_point;
        }
        else
          istatus_ = intersect_newton_failed;
        // if is done during triangulation
        if (close_to_shared_edge)
          *close_to_shared_edge = (ci.GetSideLocationTriangleSplit().WithinSide());

        return istatus_;
      }
      /** \brief get one of the two triangles with id 0 or 1
       *  of a QUAD4 element
       *
       *  tri3_id=0 ---> Quad4 nodes = {0 1 2}
       *  tri3_id=1 ---> Quad4 nodes = {2 3 0} */
      void GetTriangle(LINALG::Matrix<3, 3>& xyze_triElement, const unsigned& tri3_id)
      {
        if (sidetype == DRT::Element::quad4)
        {
          /* here it is important that the triangle is created in the same rotation
           * as the QUAD4 is, to get normal in the same direction and therefore the
           * same signed distance!!! */
          KERNEL::SplitQuad4IntoTri3(xyze_surfaceElement_, tri3_id, xyze_triElement);
        }
        else
        {
          dserror("Cut::Intersection::GetTriangle: For Triangulation a QUAD4 is expected!");
        }
      };

      /** \brief Update ComputeDistance routine to get information about location from the
       * cut_kernel This function is used for normal ComputeDistance (without triangulation) */
      bool ComputeDistance(LINALG::Matrix<probdim, 1> point, double& distance, double& tolerance,
          bool& zeroarea, KERNEL::PointOnSurfaceLoc& loc, std::vector<int>& touched_edges,
          bool signeddistance = false)
      {
        switch (GetOptionsPtr()->GeomDistance_Floattype())
        {
          case INPAR::CUT::floattype_cln:
          {
            return ComputeDistanceT<INPAR::CUT::floattype_cln>(
                point, distance, tolerance, zeroarea, loc, touched_edges, signeddistance);
          }
          case INPAR::CUT::floattype_double:
          {
            return ComputeDistanceT<INPAR::CUT::floattype_double>(
                point, distance, tolerance, zeroarea, loc, touched_edges, signeddistance);
          }
          default:
            throw std::runtime_error("Unexpected floattype for ComputeDistanceT!");
        }
      }

      template <INPAR::CUT::CUT_Floattype floattype>
      bool ComputeDistanceT(LINALG::Matrix<probdim, 1> point, double& distance, double& tolerance,
          bool& zeroarea, KERNEL::PointOnSurfaceLoc& loc, std::vector<int>& touched_edges,
          bool signeddistance = false)
      {
        TEUCHOS_FUNC_TIME_MONITOR("ComputeDistance");

        if (dimside + dimedge != probdim)
          dserror(
              "This ComputeDistance variant won't work! Think about using "
              "a GEO::CUT::Position object instead!");
        LINALG::Matrix<probdim, 1> xsi(xsi_.A(), true);

        KERNEL::ComputeDistance<probdim, sidetype, floattype> cd(xsi);

        bool conv = cd(xyze_surfaceElement_, point, distance, signeddistance);
        tolerance = cd.GetTolerance();
        zeroarea = cd.ZeroArea();
        loc = cd.GetSideLocation();
        cd.GetTouchedSideEdges(touched_edges);
        if (not(loc.WithinSide()))
        {
          touched_edges.clear();
        }
        FixDistantTouchingEdges<probdim, floattype>(point, touched_edges);

        return conv;
      }

      bool ComputeDistance(Point* p, double& distance, double& tolerance, bool& zeroarea,
          KERNEL::PointOnSurfaceLoc& loc, std::vector<int>& touched_edges,
          bool signeddistance = false)
      {
        LINALG::Matrix<probdim, 1> point(p->X());
        return ComputeDistance(
            point, distance, tolerance, zeroarea, loc, touched_edges, signeddistance);
      }


      // Transform IDs of the edges in the one of triangle in the splitted quad4 into the ids of
      // quad4 edges
      void GetQuadEdgeIdsFromTri(std::vector<int>& quad4_touched_edges,
          const std::vector<int>& tri_touched_edges_ids, int tri_id)
      {
        // NOTE: Transformation follow from the transformation function in the cut_kernel
        // SplitQuad4IntoTri3, see notes about ids there first transform normal edges
        int triangle = 2 * tri_id;
        std::vector<int> allowed_ids;  // (0, 2) ( 1 is diagonal and is ignored)
        int count_id;
        for (std::vector<int>::const_iterator e_it = tri_touched_edges_ids.begin();
             e_it != tri_touched_edges_ids.end(); ++e_it)
        {
          if ((*e_it) == 0)
            count_id = 0;
          else if ((*e_it) == 1)
            count_id = 1;
          // else it is diagnoal
          else
            continue;
          int quad4_id = triangle + count_id;
          quad4_touched_edges.push_back(quad4_id);
        }

        if (quad4_touched_edges.size() > 4) dserror("this should not be possible");
      }


      // Calculates if all nodal points of this quad4 belong to the same plane
      // (if any nodal point lie on the plane created by other 3).
      // Splitting into triangles is done in the same way, as in the other routines
      // Calculation is done both with same tolerance as in cut_kernel, as well, as
      // more tight tolerance of 1e-30
      std::pair<bool, bool> IsQuad4Distorted();

      /* Used during splitting of quad4 into tri3 and computing distance to them */
      bool ComputeDistance(LINALG::Matrix<3, 1> point, double& distance, double& tolerance,
          bool& zeroarea, KERNEL::PointOnSurfaceLoc& loc, std::vector<int>& touched_edges,
          bool signeddistance, int tri3_id, bool& extended_tri_tolerance_loc_triangle_split)
      {
        switch (GetOptionsPtr()->GeomDistance_Floattype())
        {
          case INPAR::CUT::floattype_cln:
          {
            return ComputeDistanceT<INPAR::CUT::floattype_cln>(point, distance, tolerance, zeroarea,
                loc, touched_edges, signeddistance, tri3_id,
                extended_tri_tolerance_loc_triangle_split);
          }
          case INPAR::CUT::floattype_double:
          {
            return ComputeDistanceT<INPAR::CUT::floattype_double>(point, distance, tolerance,
                zeroarea, loc, touched_edges, signeddistance, tri3_id,
                extended_tri_tolerance_loc_triangle_split);
          }
          default:
            throw std::runtime_error("Unexpected floattype for ComputeDistanceT!");
        }
      }

      template <INPAR::CUT::CUT_Floattype floattype>
      bool ComputeDistanceT(LINALG::Matrix<3, 1> point, double& distance, double& tolerance,
          bool& zeroarea, KERNEL::PointOnSurfaceLoc& loc, std::vector<int>& touched_edges,
          bool signeddistance, int tri3_id, bool& extended_tri_tolerance_loc_triangle_split)
      {
        if (sidetype != DRT::Element::quad4)
          dserror(
              "This ComputeDistance routine is only meaningful for "
              "QUAD4 side elements! But you passed in a side element "
              "of type %i | %s.",
              sidetype, DRT::DistypeToString(sidetype).c_str());

        TEUCHOS_FUNC_TIME_MONITOR("ComputeDistance");

        // dimension of xsi: element dimension of 2 + 1 entry for the distance
        if (xsi_.M() != 3) dserror("xsi_ has the wrong dimension! (dimedge + 2 = %d + 2)", dimedge);
        LINALG::Matrix<3, 1> xsi(xsi_.A(), true);

        // KERNEL::DebugComputeDistance<probdim,DRT::Element::tri3, floattype> cd( xsi );
        KERNEL::ComputeDistance<3, DRT::Element::tri3, floattype> cd(xsi);

        LINALG::Matrix<3, 3> xyze_triElement;
        GetTriangle(xyze_triElement, tri3_id);

        bool conv = cd(xyze_triElement, point, distance, signeddistance);
        tolerance = cd.GetTolerance();
        loc = cd.GetSideLocation();
        zeroarea = cd.ZeroArea();
        std::vector<int> tri_touched_edges;
        cd.GetTouchedSideEdges(tri_touched_edges);

        GetQuadEdgeIdsFromTri(touched_edges, tri_touched_edges, tri3_id);
        extended_tri_tolerance_loc_triangle_split =
            (cd.GetSideLocationTriangleSplit().WithinSide());
        FixDistantTouchingEdges<3, floattype>(point, touched_edges);

        return conv;
      }

      /// get the coordinates of the point and call the related ComputeDistance routine
      bool ComputeDistance(Point* p, double& distance, double& tolerance, bool& zeroarea,
          KERNEL::PointOnSurfaceLoc& loc, std::vector<int>& touched_edges, bool signeddistance,
          int tri3_id, bool& extended_tri_tolerance_loc_triangle_split)
      {
        LINALG::Matrix<3, 1> point(p->X());
        return ComputeDistance(point, distance, tolerance, zeroarea, loc, touched_edges,
            signeddistance, tri3_id, extended_tri_tolerance_loc_triangle_split);
      }

      /** \brief check if the two given edges \c e1 and \c e2 intersect
       *
       *  This routine computes the intersection point between edge \c sedge and
       *  edge \c eedge, and returns TRUE if the computation was successful AND
       *  the local coordinate of the cut point is within the element limits
       *  of the two edges.
       *
       *  \param sedge         (in)  : one edge of the side element
       *  \param eedge         (in)  : edge which is supposed to cut the side
       *  \param side          (in)  : side
       *  \param ee_cut_points (out) : cut points between the two given edges
       *  \param tolerance     (out) : used internal adaptive tolerance
       *                           ( specified by the CUT::KERNEL )
       *
       *  \author hiermeier \date 08/16 */
      bool ComputeCut(
          Edge* sedge, Edge* eedge, Side* side, PointSet& ee_cut_points, double& tolerance);

      /// add cut point that is a node to all edges and sides it touches
      void InsertCut(Node* n, PointSet& cuts)
      {
        cuts.insert(Point::InsertCut(GetEdgePtr(), GetSidePtr(), n));
      }

      void AddConnectivityInfo(Point* p, const LINALG::Matrix<probdim, 1>& xreal,
          const std::vector<int>& touched_vertices_ids, const std::vector<int>& touched_edges_ids);

      void AddConnectivityInfo(Point* p, const LINALG::Matrix<probdim, 1>& xreal,
          const std::vector<int>& touched_edges_ids,
          const std::set<std::pair<Side*, Edge*>>& touched_cut_pairs);

      void GetConnectivityInfo(const LINALG::Matrix<probdim, 1>& xreal,
          const std::vector<int>& touched_edges_ids, std::set<std::pair<Side*, Edge*>>& out);

      bool RefinedBBOverlapCheck(int maxstep = 10);

     protected:
      static LINALG::Matrix<probdim, numNodesEdge> xyze_lineElement_;
      static LINALG::Matrix<probdim, numNodesSide> xyze_surfaceElement_;

      static LINALG::Matrix<dimedge + dimside, 1> xsi_;
      LINALG::Matrix<dimside, 1> xsi_side_;
      LINALG::Matrix<dimedge, 1> xsi_edge_;
      static LINALG::Matrix<probdim, 1> x_;

      std::vector<LINALG::Matrix<dimside, 1>> multiple_xsi_side_;
      std::vector<LINALG::Matrix<dimedge, 1>> multiple_xsi_edge_;

      unsigned num_cut_points_;

      /// intersection status
      IntersectionStatus istatus_;

      /// scaling calculated based on the input element
      double scale_;

      /// shifting calculated based on the input element
      LINALG::Matrix<probdim, 1> shift_;
    };  // class Intersection

    /*--------------------------------------------------------------------------*/
    /** \brief Create a Intersection object
     *
     *  \author hiermeier \date 12/16 */
    class IntersectionFactory
    {
     public:
      IntersectionFactory(){};

      Teuchos::RCP<IntersectionBase> CreateIntersection(DRT::Element::DiscretizationType edge_type,
          DRT::Element::DiscretizationType side_type) const;

     private:
      template <DRT::Element::DiscretizationType edgeType>
      IntersectionBase* CreateIntersection(
          DRT::Element::DiscretizationType side_type, int probdim) const
      {
        switch (side_type)
        {
          case DRT::Element::quad4:
            return CreateConcreteIntersection<edgeType, DRT::Element::quad4>(probdim);
          case DRT::Element::quad8:
            return CreateConcreteIntersection<edgeType, DRT::Element::quad8>(probdim);
          case DRT::Element::quad9:
            return CreateConcreteIntersection<edgeType, DRT::Element::quad9>(probdim);
          case DRT::Element::tri3:
            return CreateConcreteIntersection<edgeType, DRT::Element::tri3>(probdim);
          case DRT::Element::line2:
            return CreateConcreteIntersection<edgeType, DRT::Element::line2>(probdim);
          default:
            dserror(
                "Unsupported SideType! If meaningful, add your sideType here. \n"
                "Given SideType = %s",
                DRT::DistypeToString(side_type).c_str());
            exit(EXIT_FAILURE);
        }
        exit(EXIT_FAILURE);
      }

      template <DRT::Element::DiscretizationType edgeType,
          DRT::Element::DiscretizationType sideType>
      IntersectionBase* CreateConcreteIntersection(const int& probdim) const
      {
        GEO::CUT::IntersectionBase* inter_ptr = NULL;
        switch (probdim)
        {
          case 2:
            inter_ptr = new GEO::CUT::Intersection<2, edgeType, sideType>();
            break;
          case 3:
            inter_ptr = new GEO::CUT::Intersection<3, edgeType, sideType>();
            break;
          default:
            dserror("Unsupported ProbDim! ( probdim = %d )", probdim);
            exit(EXIT_FAILURE);
        }
        return inter_ptr;
      };

    };  // class IntersectionFactory

  }  // namespace CUT
}  // namespace GEO

// static members of Intersection base class
template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
LINALG::Matrix<probdim, numNodesEdge> GEO::CUT::Intersection<probdim, edgetype, sidetype, debug,
    dimedge, dimside, numNodesEdge, numNodesSide>::xyze_lineElement_;
template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
LINALG::Matrix<probdim, numNodesSide> GEO::CUT::Intersection<probdim, edgetype, sidetype, debug,
    dimedge, dimside, numNodesEdge, numNodesSide>::xyze_surfaceElement_;
template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
LINALG::Matrix<dimedge + dimside, 1> GEO::CUT::Intersection<probdim, edgetype, sidetype, debug,
    dimedge, dimside, numNodesEdge, numNodesSide>::xsi_;
template <unsigned probdim, DRT::Element::DiscretizationType edgetype,
    DRT::Element::DiscretizationType sidetype, bool debug, unsigned dimedge, unsigned dimside,
    unsigned numNodesEdge, unsigned numNodesSide>
LINALG::Matrix<probdim, 1> GEO::CUT::Intersection<probdim, edgetype, sidetype, debug, dimedge,
    dimside, numNodesEdge, numNodesSide>::x_;

#endif

/*----------------------------------------------------------------------*/
/*! \file

\brief here the intersection of a (plane) surface with a line is performed

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_CUT_INTERSECTION_H
#define BACI_CUT_INTERSECTION_H

#include "baci_cut_edge.H"
#include "baci_cut_kernel.H"
#include "baci_cut_mesh.H"
#include "baci_cut_node.H"
#include "baci_cut_options.H"
#include "baci_cut_output.H"
#include "baci_cut_side.H"

#include <Teuchos_TimeMonitor.hpp>


namespace CORE::GEO
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
          const CORE::FE::CellType& edge_type, const CORE::FE::CellType& side_type);

     public:
      /// constructor
      IntersectionBase()
          : isinit_(false),
            isscaled_(false),
            isshifted_(false),
            useboundingbox_(false),
            mesh_ptr_(nullptr),
            edge_ptr_(nullptr),
            side_ptr_(nullptr),
            options_ptr_(nullptr)
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

        mesh_ptr_ = nullptr;
        edge_ptr_ = nullptr;
        side_ptr_ = nullptr;
        options_ptr_ = options;

        if (static_cast<unsigned>(xyze_lineElement.numRows()) != ProbDim() or
            static_cast<unsigned>(xyze_lineElement.numCols()) != NumNodesEdge())
          dserror(
              "Dimension mismatch of xyze_lineElement! \n"
              "expected input: %d x %d (rows x cols)\n"
              "current input : %d x %d (rows x cols)",
              ProbDim(), NumNodesEdge(), xyze_lineElement.numRows(), xyze_lineElement.numCols());

        if (static_cast<unsigned>(xyze_surfaceElement.numRows()) != ProbDim() or
            static_cast<unsigned>(xyze_surfaceElement.numCols()) != NumNodesSide())
          dserror(
              "Dimension mismatch of xyze_surfaceElement! \n"
              "expected input: %d x %d (rows x cols)\n"
              "current input : %d x %d (rows x cols)",
              ProbDim(), NumNodesSide(), xyze_surfaceElement.numRows(),
              xyze_surfaceElement.numCols());

        SetCoordinates(xyze_surfaceElement.values(), xyze_lineElement.values());
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
      virtual ~IntersectionBase() = default;

      /** \brief Calculate the actual intersection of an edge and a side ( or 2-nd edge )
       *
       *  See derived class for more information.
       *
       *  \author hiermeier \date 08/16 */
      virtual IntersectionStatus ComputeEdgeSideIntersection(double& tolerance,
          bool check_inside = true, std::vector<int>* touched_edges = nullptr) = 0;

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

      /// Get the coordinates of the computed point from Edge-Edge Intersection
      virtual double* FinalPointEdgeEdge() = 0;

      /** Access the cut point local coordinates on the side element
       * ( also working for multiple cut points )
       *
       *  \author hiermeier \date 01/17 */
      template <unsigned dimside>
      void LocalSideCoordinates(std::vector<CORE::LINALG::Matrix<dimside, 1>>& side_rs_cuts)
      {
        if (GetIntersectionStatus() < intersect_single_cut_point)
          dserror("INVALID IntersectionStatus! ( istatus = \"%s\" )",
              IntersectionStatus2String(GetIntersectionStatus()).c_str());

        side_rs_cuts.clear();
        side_rs_cuts.reserve(NumCutPoints());

        for (unsigned i = 0; i < NumCutPoints(); ++i)
          side_rs_cuts.push_back(CORE::LINALG::Matrix<dimside, 1>(LocalSideCoordinates(i), true));
      }

      /** Access the final cut point global coordinates
       * ( also working for multiple cut points )
       *
       *  \author hiermeier \date 01/17 */
      template <unsigned probdim>
      void FinalPoints(std::vector<CORE::LINALG::Matrix<probdim, 1>>& xyz_cuts)
      {
        if (GetIntersectionStatus() < intersect_single_cut_point)
          dserror("INVALID IntersectionStatus! ( istatus = \"%s\" )",
              IntersectionStatus2String(GetIntersectionStatus()).c_str());

        xyz_cuts.clear();
        xyz_cuts.reserve(NumCutPoints());

        for (unsigned i = 0; i < NumCutPoints(); ++i)
          xyz_cuts.push_back(CORE::LINALG::Matrix<probdim, 1>(FinalPoint(i), false));
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
        if (mesh_ptr_ != nullptr) return *mesh_ptr_;
        dserror("The mesh pointer is not yet initialized!");
        exit(EXIT_FAILURE);
      }

      /// get a pointer to the mesh object
      Mesh* GetMeshPtr()
      {
        if (mesh_ptr_ != nullptr) return mesh_ptr_;
        dserror("The mesh pointer is not yet initialized!");
        exit(EXIT_FAILURE);
      }

      /// get a reference to the edge object
      Edge& GetEdge()
      {
        if (edge_ptr_ != nullptr) return *edge_ptr_;
        dserror("The edge pointer is not yet initialized!");
        exit(EXIT_FAILURE);
      }

      /// get a pointer to the edge object
      Edge* GetEdgePtr()
      {
        if (edge_ptr_ != nullptr) return edge_ptr_;
        dserror("The edge pointer is not yet initialized!");
        exit(EXIT_FAILURE);
      }

      /// get a reference to the side object
      Side& GetSide()
      {
        if (side_ptr_ != nullptr) return *side_ptr_;
        dserror("The side pointer is not yet initialized!");
        exit(EXIT_FAILURE);
      }

      /// get a pointer to the side object
      Side* GetSidePtr()
      {
        if (side_ptr_ != nullptr) return side_ptr_;
        dserror("The side pointer is not yet initialized!");
        exit(EXIT_FAILURE);
      }

      /// get a pointer to the cut options object
      Options* GetOptionsPtr()
      {
        if (options_ptr_ != nullptr) return options_ptr_;
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
    template <unsigned probdim, CORE::FE::CellType edgetype, CORE::FE::CellType sidetype,
        bool debug = false, unsigned dimedge = CORE::DRT::UTILS::DisTypeToDim<edgetype>::dim,
        unsigned dimside = CORE::DRT::UTILS::DisTypeToDim<sidetype>::dim,
        unsigned numNodesEdge = CORE::FE::num_nodes<edgetype>,
        unsigned numNodesSide = CORE::FE::num_nodes<sidetype>>
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

      // No public access to these methods! Use the base class accessors, instead.
     protected:
      /// get the number of detected feasible cut points
      unsigned NumCutPoints() const override
      {
        if (num_cut_points_ > 1 and (multiple_xsi_edge_.size() != num_cut_points_ or
                                        multiple_xsi_side_.size() != num_cut_points_))
          dserror("Size mismatch!");

        return num_cut_points_;
      }

      /// get the intersection status
      const IntersectionStatus& GetIntersectionStatus() const override { return istatus_; }

      /// get the local cut coordinates
      double* LocalCoordinates() override { return xsi_.A(); }

      /** \brief access the local coordinates of the cut point corresponding to the
       *  cut point ID \c cp_id on the side element
       *
       *  \param cp_id (in) : cut point id
       *
       *  \author hiermeier \date 01/17 */
      double* LocalSideCoordinates(unsigned cp_id) override
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
      const CORE::LINALG::Matrix<dimedge, 1>& LocalEdgeCoordinates(const unsigned& cp_id) const
      {
        if (NumCutPoints() < 2) return xsi_edge_;

        return multiple_xsi_edge_[cp_id];
      }

      /// We need choose the edge, based on which we compute the global coordinates in a smart
      /// way -> If we choose so that the cut point will be close to endpoint, we essentially
      /// extend its tolerance and it therefore would lead to problems in the cut
      double* FinalPointEdgeEdge() override
      {
        CheckInit();
        if (dimedge != dimside) dserror("This method only works for edge-edge intersection!");

        CORE::LINALG::Matrix<probdim, 1> x_edge_1;
        x_edge_1 = 0;
        CORE::LINALG::Matrix<numNodesEdge, 1> edge_funct_1;
        CORE::DRT::UTILS::shape_function<edgetype>(xsi_edge_, edge_funct_1);
        for (unsigned inode = 0; inode < numNodesEdge; ++inode)
          for (unsigned isd = 0; isd < probdim; ++isd)
            x_edge_1(isd) += xyze_lineElement_(isd, inode) * edge_funct_1(inode);

        // first un-do the shifting
        x_edge_1.Update(1.0, shift_, 1.0);
        x_edge_1.Scale(scale_);

        CORE::LINALG::Matrix<probdim, 1> x_edge_2;
        x_edge_2 = 0;
        CORE::LINALG::Matrix<numNodesSide, 1> edge_funct_2;
        CORE::DRT::UTILS::shape_function<sidetype>(xsi_side_, edge_funct_2);
        for (unsigned inode = 0; inode < numNodesSide; ++inode)
          for (unsigned isd = 0; isd < probdim; ++isd)
            x_edge_2(isd) += xyze_surfaceElement_(isd, inode) * edge_funct_2(inode);

        x_edge_2.Update(1.0, shift_, 1.0);
        x_edge_2.Scale(scale_);

        bool will_be_merged[2];
        will_be_merged[0] =
            IsCloseToEndpoints(xyze_lineElement_, x_edge_1, SIDE_DETECTION_TOLERANCE) or
            IsCloseToEndpoints(xyze_surfaceElement_, x_edge_1, SIDE_DETECTION_TOLERANCE);

        will_be_merged[1] =
            IsCloseToEndpoints(xyze_lineElement_, x_edge_2, SIDE_DETECTION_TOLERANCE) or
            IsCloseToEndpoints(xyze_surfaceElement_, x_edge_2, SIDE_DETECTION_TOLERANCE);

        if (will_be_merged[0] and will_be_merged[1])
        {
          // NOTE: If such case ever occur, one might consider to create intersection on the
          // line between x_edge_2 and x_edge_1, in a way so that it will not be merged in any of
          // the participating edge endpoint
          dserror(
              "Cannot decide which edge should serve as the basis for global coordinates in "
              "edge-edge intersection");
        }
        else if (will_be_merged[0])
        {
          x_ = x_edge_2;
        }
        else if (will_be_merged[1])
        {
          x_ = x_edge_1;
        }
        else
        {
          // Both edges could serve as a good basis, we take the edge where the intersection point
          // is further away from the end point.
          if (std::abs(xsi_edge_(0, 0)) < std::abs(xsi_side_(0, 0)))
            x_ = x_edge_1;
          else
            x_ = x_edge_2;
        }
        return x_.A();
      }

      /// get the final cut point global coordinates
      void FinalPoint(
          const CORE::LINALG::Matrix<dimedge, 1>& xsi_edge, CORE::LINALG::Matrix<probdim, 1>& x)
      {
        CheckInit();

        // get final point
        x = 0;
        CORE::LINALG::Matrix<numNodesEdge, 1> lineFunct;
        CORE::DRT::UTILS::shape_function<edgetype>(xsi_edge, lineFunct);
        for (unsigned inode = 0; inode < numNodesEdge; ++inode)
          for (unsigned isd = 0; isd < probdim; ++isd)
            x(isd) += xyze_lineElement_(isd, inode) * lineFunct(inode);

        // first un-do the shifting
        x.Update(1.0, shift_, 1.0);
        // second un-do the scaling
        x.Scale(scale_);
      }
      double* FinalPoint() override
      {
        if (istatus_ != intersect_single_cut_point)
          dserror(
              "INVALID IntersectionStatus: This routine is restricted to one single "
              "cut point only! ( istatus_ = \"%s\" )",
              IntersectionStatus2String(istatus_).c_str());

        FinalPoint(xsi_edge_, x_);
        return x_.A();
      }
      double* FinalPoint(unsigned cp_id) override
      {
        FinalPoint(LocalEdgeCoordinates(cp_id), x_);
        return x_.A();
      }

      // Remove all the edges from list of touching edges that are further away than 1e-14 (
      // TOPOLOGICAL_TOLERANCE ) distance from the point
      template <unsigned int dim, INPAR::CUT::CUT_Floattype floattype>
      void FixDistantTouchingEdges(
          CORE::LINALG::Matrix<dim, 1>& p_coord, std::vector<int>& touching_edges)
      {
        bool signeddistance = false;
        double distance = 0;
        CORE::LINALG::Matrix<dim, 1> xsi;
        CORE::LINALG::Matrix<dim, numNodesEdge> xyze_edge;
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
                CORE::GEO::CUT::OUTPUT::GmshEdgeDump(file, e, std::string("FarEdge"));
                CORE::GEO::CUT::OUTPUT::GmshNewSection(file, "Point", false);
                CORE::LINALG::Matrix<3, 1> p(p_coord.A());
                CORE::GEO::CUT::OUTPUT::GmshCoordDump(file, p, 0);
                CORE::GEO::CUT::OUTPUT::GmshEndSection(file);
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
       *  This function returns CORE::GEO::CUT::IntersectionStatus. There are three different
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
      IntersectionStatus ComputeEdgeSideIntersection(double& tolerance, bool check_inside = true,
          std::vector<int>* touched_edges = nullptr) override
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
          double& tolerance, bool check_inside = true, std::vector<int>* touched_edges = nullptr)
      {
        CheckInit();
        TEUCHOS_FUNC_TIME_MONITOR("ComputeEdgeSideIntersection");

        const bool success = CheckParallelism(multiple_xsi_side_, multiple_xsi_edge_, tolerance);

        /* The parallelism check was successful and we are done. At this point
         * it is possible, that we find more than one cut point. A special
         * treatment becomes necessary for multiple cut points. */
        if (success)
        {
          istatus_ = CORE::GEO::CUT::IntersectionStatus2Enum(num_cut_points_);
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
            istatus_ = CORE::GEO::CUT::IntersectionStatus2Enum(num_cut_points_);
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
      bool Intersect(PointSet& cuts) override;

      // Try to find possible intersection points, if this intersection is between parallell size
      // and edge without using real ComputeIntersection
      ParallelIntersectionStatus HandleParallelIntersection(
          PointSet& cuts, int id = -1, bool output = false) override;

      virtual void GenerateGmshDump();

      // Handle cases for which normal intersection procedure did not work
      bool HandleSpecialCases() override;

      // compute intersection by splitting quad4 into two triangles
      bool TriangulatedIntersection(PointSet& cuts) override;

      /** \brief Will return TRUE, if local side coordinates are within the side
       *  element parameter space bounds */
      bool SurfaceWithinLimits(double tol = REFERENCETOL) const override
      {
        return CORE::GEO::CUT::KERNEL::WithinLimits<sidetype>(xsi_side_, tol);
      }

      /** \brief Will return TRUE, if local side coordinates are within the TRI3
       *  side element parameter space bounds */
      bool Tri3WithinLimits(double tol = REFERENCETOL) const
      {
        return CORE::GEO::CUT::KERNEL::WithinLimits<CORE::FE::CellType::tri3>(xsi_side_, tol);
      }

      /** \brief Will return TRUE, if local edge coordinate is within the line
       *  element parameter space bounds */
      bool LineWithinLimits(double tol = REFERENCETOL) const override
      {
        return CORE::GEO::CUT::KERNEL::WithinLimits<edgetype>(xsi_edge_, tol);
      }

      /// access the problem dimension
      unsigned ProbDim() const override { return probdim; }

      /// access the number of nodes of the side ( or 2-nd edge ) element
      unsigned NumNodesSide() const override { return numNodesSide; }

      /// access the number of nodes of the edge
      unsigned NumNodesEdge() const override { return numNodesEdge; }

      /// set the edge and side coordinates
      void SetCoordinates() override;

      /// set the edge and side coordinates
      void SetCoordinates(double* xyze_surfaceElement, double* xyze_lineElement) override
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
      void ScaleAndShift() override
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
            CORE::LINALG::Matrix<probdim, 1> x1(&xyze_surfaceElement_(0, i), true);
            x1.Update(-1, shift_, 1);
          }
          for (unsigned i = 0; i < numNodesEdge; ++i)
          {
            CORE::LINALG::Matrix<probdim, 1> x1(&xyze_lineElement_(0, i), true);
            x1.Update(-1, shift_, 1);
          }
        }
      }

      /** check if the given local coordinates are at one of the edges of the side element,
       *  i.e. at the boundaries of the side element. */
      template <class T>
      bool AtEdge(const T& xsi)
      {
        return CORE::GEO::CUT::KERNEL::AtEdge<sidetype>(xsi);
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
      bool CheckParallelism(std::vector<CORE::LINALG::Matrix<dimside, 1>>& side_rs_intersect,
          std::vector<CORE::LINALG::Matrix<dimedge, 1>>& edge_r_intersect, double& tolerance);

      /** \brief Check if the two lines are collinear, end points are on the line, or
       *  the distance values imply that no intersection is possible
       *
       *  \author hiermeier \date 12/16 */
      bool CheckCollinearity(
          std::vector<CORE::LINALG::Matrix<dimside, 1>>& side_rs_corner_intersect,
          std::vector<CORE::LINALG::Matrix<dimedge, 1>>& edge_r_corner_intersect,
          double& tolerance);

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
          std::vector<CORE::LINALG::Matrix<dimside, 1>>& side_rs_intersect,
          std::vector<CORE::LINALG::Matrix<dimedge, 1>>& edge_r_intersect, double& tolerance);

      /** \brief Check the angle between a edge and a surface normal.
       *
       *  This is a quick check to skip cases which are definitely not parallel.
       *
       *  \author hiermeier \date 12/16 */
      bool CheckAngleCriterionBetweenSideNormalAndEdge();

      /// find the local coordinate of a given edge end point ( i.e. -1 or +1 )
      bool FindLocalCoordinateOfEdgeEndPoint(
          double& pos, const CORE::LINALG::Matrix<probdim, 1>& xyz, const double& tol) const;

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
        CORE::LINALG::Matrix<3, 1> xsi;
        if (xsi_.M() != 3)
          dserror("xsi_ has the wrong dimension! (dimedge + 2 = %d + 2)", dimedge);
        else
          xsi.SetView(xsi_.A());

        KERNEL::ComputeIntersection<3, edgetype, CORE::FE::CellType::tri3, floattype> ci(xsi);
        // KERNEL::DebugComputeIntersection<probdim,edgetype,CORE::FE::CellType::tri3,floattype>
        // ci( xsi
        // );

        CORE::LINALG::Matrix<3, 3> xyze_triElement;
        GetTriangle(xyze_triElement, triangleid);
        CORE::LINALG::Matrix<3, numNodesEdge> xyze_lineElement(xyze_lineElement_.A(), true);

        bool conv = ci(xyze_triElement, xyze_lineElement);
        location = ci.GetSideLocation();

        return conv;
      }

      // Computes tri3 edge intersection used for quad4 -> 2 tri3 splits
      IntersectionStatus ComputeEdgeTri3IntersectionQuad4Split(
          int triangleid, bool* close_to_shared_edge = nullptr)
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
          int triangleid, bool* close_to_shared_edge = nullptr)
      {
        if (triangleid < 0) dserror("The triangle id has to be positive!");

        TEUCHOS_FUNC_TIME_MONITOR("ComputeEdgeTri3Intersection");
        CORE::LINALG::Matrix<3, 1> xsi;
        if (xsi_.M() != 3)
          dserror("xsi_ has the wrong dimension! (dimedge + 2 = %d + 2)", dimedge);
        else
          xsi.SetView(xsi_.A());

        KERNEL::ComputeIntersection<3, edgetype, CORE::FE::CellType::tri3, floattype> ci(xsi);
        // KERNEL::DebugComputeIntersection<probdim,edgetype,CORE::FE::CellType::tri3,floattype>
        // ci( xsi
        // );

        CORE::LINALG::Matrix<3, 3> xyze_triElement;
        GetTriangle(xyze_triElement, triangleid);
        CORE::LINALG::Matrix<3, numNodesEdge> xyze_lineElement(xyze_lineElement_.A(), true);

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
      void GetTriangle(CORE::LINALG::Matrix<3, 3>& xyze_triElement, const unsigned& tri3_id)
      {
        if (sidetype == CORE::FE::CellType::quad4)
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
      bool ComputeDistance(CORE::LINALG::Matrix<probdim, 1> point, double& distance,
          double& tolerance, bool& zeroarea, KERNEL::PointOnSurfaceLoc& loc,
          std::vector<int>& touched_edges, bool signeddistance = false)
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
      bool ComputeDistanceT(CORE::LINALG::Matrix<probdim, 1> point, double& distance,
          double& tolerance, bool& zeroarea, KERNEL::PointOnSurfaceLoc& loc,
          std::vector<int>& touched_edges, bool signeddistance = false)
      {
        TEUCHOS_FUNC_TIME_MONITOR("ComputeDistance");

        if (dimside + dimedge != probdim)
          dserror(
              "This ComputeDistance variant won't work! Think about using "
              "a CORE::GEO::CUT::Position object instead!");
        CORE::LINALG::Matrix<probdim, 1> xsi(xsi_.A(), true);

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
        CORE::LINALG::Matrix<probdim, 1> point(p->X());
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

      /// Detects in the point is close to an endpoint of the edge
      template <unsigned int numNodes, unsigned int probDim>
      bool IsCloseToEndpoints(const CORE::LINALG::Matrix<probDim, numNodes>& surf,
          const CORE::LINALG::Matrix<probDim, 1>& p, double tol = TOPOLOGICAL_TOLERANCE)
      {
        for (unsigned int node_id = 0; node_id < numNodes; ++node_id)
        {
          const CORE::LINALG::Matrix<probDim, 1> edge_point(surf.A() + node_id * probDim, true);
          if (CORE::GEO::CUT::DistanceBetweenPoints(edge_point, p) <= tol) return true;
        }
        return false;
      }

      // Calculates if all nodal points of this quad4 belong to the same plane
      // (if any nodal point lie on the plane created by other 3).
      // Splitting into triangles is done in the same way, as in the other routines
      // Calculation is done both with same tolerance as in cut_kernel, as well, as
      // more tight tolerance of 1e-30
      std::pair<bool, bool> IsQuad4Distorted();

      /* Used during splitting of quad4 into tri3 and computing distance to them */
      bool ComputeDistance(CORE::LINALG::Matrix<3, 1> point, double& distance, double& tolerance,
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
      bool ComputeDistanceT(CORE::LINALG::Matrix<3, 1> point, double& distance, double& tolerance,
          bool& zeroarea, KERNEL::PointOnSurfaceLoc& loc, std::vector<int>& touched_edges,
          bool signeddistance, int tri3_id, bool& extended_tri_tolerance_loc_triangle_split)
      {
        if (sidetype != CORE::FE::CellType::quad4)
          dserror(
              "This ComputeDistance routine is only meaningful for "
              "QUAD4 side elements! But you passed in a side element "
              "of type %i | %s.",
              sidetype, CORE::FE::CellTypeToString(sidetype).c_str());

        TEUCHOS_FUNC_TIME_MONITOR("ComputeDistance");

        // dimension of xsi: element dimension of 2 + 1 entry for the distance
        if (xsi_.M() != 3) dserror("xsi_ has the wrong dimension! (dimedge + 2 = %d + 2)", dimedge);
        CORE::LINALG::Matrix<3, 1> xsi(xsi_.A(), true);

        // KERNEL::DebugComputeDistance<probdim,CORE::FE::CellType::tri3, floattype>
        // cd( xsi );
        KERNEL::ComputeDistance<3, CORE::FE::CellType::tri3, floattype> cd(xsi);

        CORE::LINALG::Matrix<3, 3> xyze_triElement;
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
        CORE::LINALG::Matrix<3, 1> point(p->X());
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

      void AddConnectivityInfo(Point* p, const CORE::LINALG::Matrix<probdim, 1>& xreal,
          const std::vector<int>& touched_vertices_ids, const std::vector<int>& touched_edges_ids);

      void AddConnectivityInfo(Point* p, const CORE::LINALG::Matrix<probdim, 1>& xreal,
          const std::vector<int>& touched_edges_ids,
          const std::set<std::pair<Side*, Edge*>>& touched_cut_pairs);

      void GetConnectivityInfo(const CORE::LINALG::Matrix<probdim, 1>& xreal,
          const std::vector<int>& touched_edges_ids, std::set<std::pair<Side*, Edge*>>& out);

      bool RefinedBBOverlapCheck(int maxstep = 10);

     protected:
      static CORE::LINALG::Matrix<probdim, numNodesEdge> xyze_lineElement_;
      static CORE::LINALG::Matrix<probdim, numNodesSide> xyze_surfaceElement_;

      static CORE::LINALG::Matrix<dimedge + dimside, 1> xsi_;
      CORE::LINALG::Matrix<dimside, 1> xsi_side_;
      CORE::LINALG::Matrix<dimedge, 1> xsi_edge_;
      static CORE::LINALG::Matrix<probdim, 1> x_;

      std::vector<CORE::LINALG::Matrix<dimside, 1>> multiple_xsi_side_;
      std::vector<CORE::LINALG::Matrix<dimedge, 1>> multiple_xsi_edge_;

      unsigned num_cut_points_;

      /// intersection status
      IntersectionStatus istatus_;

      /// scaling calculated based on the input element
      double scale_;

      /// shifting calculated based on the input element
      CORE::LINALG::Matrix<probdim, 1> shift_;
    };  // class Intersection

    /*--------------------------------------------------------------------------*/
    /** \brief Create a Intersection object
     *
     *  \author hiermeier \date 12/16 */
    class IntersectionFactory
    {
     public:
      IntersectionFactory(){};

      Teuchos::RCP<IntersectionBase> CreateIntersection(
          CORE::FE::CellType edge_type, CORE::FE::CellType side_type) const;

     private:
      template <CORE::FE::CellType edgeType>
      IntersectionBase* CreateIntersection(CORE::FE::CellType side_type, int probdim) const
      {
        switch (side_type)
        {
          case CORE::FE::CellType::quad4:
            return CreateConcreteIntersection<edgeType, CORE::FE::CellType::quad4>(probdim);
          case CORE::FE::CellType::quad8:
            return CreateConcreteIntersection<edgeType, CORE::FE::CellType::quad8>(probdim);
          case CORE::FE::CellType::quad9:
            return CreateConcreteIntersection<edgeType, CORE::FE::CellType::quad9>(probdim);
          case CORE::FE::CellType::tri3:
            return CreateConcreteIntersection<edgeType, CORE::FE::CellType::tri3>(probdim);
          case CORE::FE::CellType::line2:
            return CreateConcreteIntersection<edgeType, CORE::FE::CellType::line2>(probdim);
          default:
            dserror(
                "Unsupported SideType! If meaningful, add your sideType here. \n"
                "Given SideType = %s",
                CORE::FE::CellTypeToString(side_type).c_str());
            exit(EXIT_FAILURE);
        }
        exit(EXIT_FAILURE);
      }

      template <CORE::FE::CellType edgeType, CORE::FE::CellType sideType>
      IntersectionBase* CreateConcreteIntersection(const int& probdim) const
      {
        CORE::GEO::CUT::IntersectionBase* inter_ptr = nullptr;
        switch (probdim)
        {
          case 2:
            inter_ptr = new CORE::GEO::CUT::Intersection<2, edgeType, sideType>();
            break;
          case 3:
            inter_ptr = new CORE::GEO::CUT::Intersection<3, edgeType, sideType>();
            break;
          default:
            dserror("Unsupported ProbDim! ( probdim = %d )", probdim);
            exit(EXIT_FAILURE);
        }
        return inter_ptr;
      };

    };  // class IntersectionFactory

  }  // namespace CUT
}  // namespace CORE::GEO

// static members of Intersection base class
template <unsigned probdim, CORE::FE::CellType edgetype, CORE::FE::CellType sidetype, bool debug,
    unsigned dimedge, unsigned dimside, unsigned numNodesEdge, unsigned numNodesSide>
CORE::LINALG::Matrix<probdim, numNodesEdge> CORE::GEO::CUT::Intersection<probdim, edgetype,
    sidetype, debug, dimedge, dimside, numNodesEdge, numNodesSide>::xyze_lineElement_;
template <unsigned probdim, CORE::FE::CellType edgetype, CORE::FE::CellType sidetype, bool debug,
    unsigned dimedge, unsigned dimside, unsigned numNodesEdge, unsigned numNodesSide>
CORE::LINALG::Matrix<probdim, numNodesSide> CORE::GEO::CUT::Intersection<probdim, edgetype,
    sidetype, debug, dimedge, dimside, numNodesEdge, numNodesSide>::xyze_surfaceElement_;
template <unsigned probdim, CORE::FE::CellType edgetype, CORE::FE::CellType sidetype, bool debug,
    unsigned dimedge, unsigned dimside, unsigned numNodesEdge, unsigned numNodesSide>
CORE::LINALG::Matrix<dimedge + dimside, 1> CORE::GEO::CUT::Intersection<probdim, edgetype, sidetype,
    debug, dimedge, dimside, numNodesEdge, numNodesSide>::xsi_;
template <unsigned probdim, CORE::FE::CellType edgetype, CORE::FE::CellType sidetype, bool debug,
    unsigned dimedge, unsigned dimside, unsigned numNodesEdge, unsigned numNodesSide>
CORE::LINALG::Matrix<probdim, 1> CORE::GEO::CUT::Intersection<probdim, edgetype, sidetype, debug,
    dimedge, dimside, numNodesEdge, numNodesSide>::x_;

#endif

/*----------------------------------------------------------------------*/
/*! \file

\brief here the intersection of a (plane) surface with a line is performed

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_INTERSECTION_HPP
#define FOUR_C_CUT_INTERSECTION_HPP

#include "4C_config.hpp"

#include "4C_cut_edge.hpp"
#include "4C_cut_kernel.hpp"
#include "4C_cut_mesh.hpp"
#include "4C_cut_node.hpp"
#include "4C_cut_options.hpp"
#include "4C_cut_output.hpp"
#include "4C_cut_side.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

namespace CORE::GEO
{
  namespace CUT
  {
    class Side;
    class Edge;

    enum IntersectionStatus
    {
      intersect_newton_failed = -2,      // if newton failed
      intersect_unevaluated = -1,        // before compute_edge_side_intersection has been called
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
       *  compute_edge_side_intersection() afterwards.
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

        if (static_cast<unsigned>(xyze_lineElement.numRows()) != prob_dim() or
            static_cast<unsigned>(xyze_lineElement.numCols()) != num_nodes_edge())
          FOUR_C_THROW(
              "Dimension mismatch of xyze_lineElement! \n"
              "expected input: %d x %d (rows x cols)\n"
              "current input : %d x %d (rows x cols)",
              prob_dim(), num_nodes_edge(), xyze_lineElement.numRows(), xyze_lineElement.numCols());

        if (static_cast<unsigned>(xyze_surfaceElement.numRows()) != prob_dim() or
            static_cast<unsigned>(xyze_surfaceElement.numCols()) != num_nodes_side())
          FOUR_C_THROW(
              "Dimension mismatch of xyze_surfaceElement! \n"
              "expected input: %d x %d (rows x cols)\n"
              "current input : %d x %d (rows x cols)",
              prob_dim(), num_nodes_side(), xyze_surfaceElement.numRows(),
              xyze_surfaceElement.numCols());

        set_coordinates(xyze_surfaceElement.values(), xyze_lineElement.values());
        scale_and_shift();

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

        set_coordinates();
        scale_and_shift();

        isinit_ = true;
      }

      /// destructor
      virtual ~IntersectionBase() = default;

      /** \brief Calculate the actual intersection of an edge and a side ( or 2-nd edge )
       *
       *  See derived class for more information.
       *
       *  \author hiermeier \date 08/16 */
      virtual IntersectionStatus compute_edge_side_intersection(double& tolerance,
          bool check_inside = true, std::vector<int>* touched_edges = nullptr) = 0;

      /** \brief Computes the intersection points of the edge with the specified side
       *  and stores the points in cuts
       *
       *  See derived class for more information.
       *
       *  \author hiermeier \date 08/16 */
      virtual bool Intersect(PointSet& cuts) = 0;

      virtual ParallelIntersectionStatus handle_parallel_intersection(
          PointSet& cuts, int id = -1, bool output = false) = 0;

      virtual bool triangulated_intersection(PointSet& cuts) = 0;

      virtual bool HandleSpecialCases() = 0;

      /** \brief Get the final cut point global coordinates
       *
       *  Only allowed if there was only one cut point!
       *
       *  \author hiermeier \date 08/16 */
      virtual double* FinalPoint() = 0;

      virtual double* FinalPoint(unsigned cp_id) = 0;

      /// Get the coordinates of the computed point from Edge-Edge intersection
      virtual double* FinalPointEdgeEdge() = 0;

      /** Access the cut point local coordinates on the side element
       * ( also working for multiple cut points )
       *
       *  \author hiermeier \date 01/17 */
      template <unsigned dimside>
      void local_side_coordinates(std::vector<CORE::LINALG::Matrix<dimside, 1>>& side_rs_cuts)
      {
        if (get_intersection_status() < intersect_single_cut_point)
          FOUR_C_THROW("INVALID IntersectionStatus! ( istatus = \"%s\" )",
              IntersectionStatus2String(get_intersection_status()).c_str());

        side_rs_cuts.clear();
        side_rs_cuts.reserve(num_cut_points());

        for (unsigned i = 0; i < num_cut_points(); ++i)
          side_rs_cuts.push_back(CORE::LINALG::Matrix<dimside, 1>(local_side_coordinates(i), true));
      }

      /** Access the final cut point global coordinates
       * ( also working for multiple cut points )
       *
       *  \author hiermeier \date 01/17 */
      template <unsigned probdim>
      void FinalPoints(std::vector<CORE::LINALG::Matrix<probdim, 1>>& xyz_cuts)
      {
        if (get_intersection_status() < intersect_single_cut_point)
          FOUR_C_THROW("INVALID IntersectionStatus! ( istatus = \"%s\" )",
              IntersectionStatus2String(get_intersection_status()).c_str());

        xyz_cuts.clear();
        xyz_cuts.reserve(num_cut_points());

        for (unsigned i = 0; i < num_cut_points(); ++i)
          xyz_cuts.push_back(CORE::LINALG::Matrix<probdim, 1>(FinalPoint(i), false));
      }

      virtual double* local_coordinates() = 0;

      virtual double* local_side_coordinates(unsigned cp_id) = 0;

      virtual bool SurfaceWithinLimits(double tol = REFERENCETOL) const = 0;

      virtual bool LineWithinLimits(double tol = REFERENCETOL) const = 0;

     protected:
      virtual unsigned num_cut_points() const = 0;

      virtual const IntersectionStatus& get_intersection_status() const = 0;

      inline void check_init() const
      {
        if (not isinit_)
          FOUR_C_THROW("The intersection object is not initialized! Call Init() first.");
      }

      virtual unsigned prob_dim() const = 0;
      virtual unsigned num_nodes_side() const = 0;
      virtual unsigned num_nodes_edge() const = 0;

      virtual void set_coordinates() = 0;
      virtual void set_coordinates(double* xyze_surfaceElement, double* xyze_lineElement) = 0;

      virtual void scale_and_shift() = 0;

      /// Are the global coordinates scaled?
      const bool& is_scaled() const { return isscaled_; }

      /// Are the global coordinates shifted?
      const bool& is_shifted() const { return isshifted_; }

      /// Shall we use the bounding box information?
      const bool& use_bounding_box() const { return useboundingbox_; }

      /// get a reference to the mesh object
      Mesh& get_mesh()
      {
        if (mesh_ptr_ != nullptr) return *mesh_ptr_;
        FOUR_C_THROW("The mesh pointer is not yet initialized!");
        exit(EXIT_FAILURE);
      }

      /// get a pointer to the mesh object
      Mesh* get_mesh_ptr()
      {
        if (mesh_ptr_ != nullptr) return mesh_ptr_;
        FOUR_C_THROW("The mesh pointer is not yet initialized!");
        exit(EXIT_FAILURE);
      }

      /// get a reference to the edge object
      Edge& get_edge()
      {
        if (edge_ptr_ != nullptr) return *edge_ptr_;
        FOUR_C_THROW("The edge pointer is not yet initialized!");
        exit(EXIT_FAILURE);
      }

      /// get a pointer to the edge object
      Edge* get_edge_ptr()
      {
        if (edge_ptr_ != nullptr) return edge_ptr_;
        FOUR_C_THROW("The edge pointer is not yet initialized!");
        exit(EXIT_FAILURE);
      }

      /// get a reference to the side object
      Side& get_side()
      {
        if (side_ptr_ != nullptr) return *side_ptr_;
        FOUR_C_THROW("The side pointer is not yet initialized!");
        exit(EXIT_FAILURE);
      }

      /// get a pointer to the side object
      Side* get_side_ptr()
      {
        if (side_ptr_ != nullptr) return side_ptr_;
        FOUR_C_THROW("The side pointer is not yet initialized!");
        exit(EXIT_FAILURE);
      }

      /// get a pointer to the cut options object
      Options* get_options_ptr()
      {
        if (options_ptr_ != nullptr) return options_ptr_;
        FOUR_C_THROW("The option pointer is not yet initialized!");
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
        bool debug = false, unsigned dimedge = CORE::FE::dim<edgetype>,
        unsigned dimside = CORE::FE::dim<sidetype>,
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
      unsigned num_cut_points() const override
      {
        if (num_cut_points_ > 1 and (multiple_xsi_edge_.size() != num_cut_points_ or
                                        multiple_xsi_side_.size() != num_cut_points_))
          FOUR_C_THROW("Size mismatch!");

        return num_cut_points_;
      }

      /// get the intersection status
      const IntersectionStatus& get_intersection_status() const override { return istatus_; }

      /// get the local cut coordinates
      double* local_coordinates() override { return xsi_.A(); }

      /** \brief access the local coordinates of the cut point corresponding to the
       *  cut point ID \c cp_id on the side element
       *
       *  \param cp_id (in) : cut point id
       *
       *  \author hiermeier \date 01/17 */
      double* local_side_coordinates(unsigned cp_id) override
      {
        if (num_cut_points() < 2) return xsi_side_.A();

        return multiple_xsi_side_[cp_id].A();
      }

      /** \brief access the local coordinates of the cut point corresponding to the
       *  cut point ID \c cp_id on the edge element
       *
       *  \param cp_id (in) : cut point id
       *
       *  \author hiermeier \date 01/17 */
      const CORE::LINALG::Matrix<dimedge, 1>& local_edge_coordinates(const unsigned& cp_id) const
      {
        if (num_cut_points() < 2) return xsi_edge_;

        return multiple_xsi_edge_[cp_id];
      }

      /// We need choose the edge, based on which we compute the global coordinates in a smart
      /// way -> If we choose so that the cut point will be close to endpoint, we essentially
      /// extend its tolerance and it therefore would lead to problems in the cut
      double* FinalPointEdgeEdge() override
      {
        check_init();
        if (dimedge != dimside) FOUR_C_THROW("This method only works for edge-edge intersection!");

        CORE::LINALG::Matrix<probdim, 1> x_edge_1;
        x_edge_1 = 0;
        CORE::LINALG::Matrix<numNodesEdge, 1> edge_funct_1;
        CORE::FE::shape_function<edgetype>(xsi_edge_, edge_funct_1);
        for (unsigned inode = 0; inode < numNodesEdge; ++inode)
          for (unsigned isd = 0; isd < probdim; ++isd)
            x_edge_1(isd) += xyze_lineElement_(isd, inode) * edge_funct_1(inode);

        // first un-do the shifting
        x_edge_1.Update(1.0, shift_, 1.0);
        x_edge_1.Scale(scale_);

        CORE::LINALG::Matrix<probdim, 1> x_edge_2;
        x_edge_2 = 0;
        CORE::LINALG::Matrix<numNodesSide, 1> edge_funct_2;
        CORE::FE::shape_function<sidetype>(xsi_side_, edge_funct_2);
        for (unsigned inode = 0; inode < numNodesSide; ++inode)
          for (unsigned isd = 0; isd < probdim; ++isd)
            x_edge_2(isd) += xyze_surfaceElement_(isd, inode) * edge_funct_2(inode);

        x_edge_2.Update(1.0, shift_, 1.0);
        x_edge_2.Scale(scale_);

        bool will_be_merged[2];
        will_be_merged[0] =
            is_close_to_endpoints(xyze_lineElement_, x_edge_1, SIDE_DETECTION_TOLERANCE) or
            is_close_to_endpoints(xyze_surfaceElement_, x_edge_1, SIDE_DETECTION_TOLERANCE);

        will_be_merged[1] =
            is_close_to_endpoints(xyze_lineElement_, x_edge_2, SIDE_DETECTION_TOLERANCE) or
            is_close_to_endpoints(xyze_surfaceElement_, x_edge_2, SIDE_DETECTION_TOLERANCE);

        if (will_be_merged[0] and will_be_merged[1])
        {
          // NOTE: If such case ever occur, one might consider to create intersection on the
          // line between x_edge_2 and x_edge_1, in a way so that it will not be merged in any of
          // the participating edge endpoint
          FOUR_C_THROW(
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
        check_init();

        // get final point
        x = 0;
        CORE::LINALG::Matrix<numNodesEdge, 1> lineFunct;
        CORE::FE::shape_function<edgetype>(xsi_edge, lineFunct);
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
          FOUR_C_THROW(
              "INVALID IntersectionStatus: This routine is restricted to one single "
              "cut point only! ( istatus_ = \"%s\" )",
              IntersectionStatus2String(istatus_).c_str());

        FinalPoint(xsi_edge_, x_);
        return x_.A();
      }
      double* FinalPoint(unsigned cp_id) override
      {
        FinalPoint(local_edge_coordinates(cp_id), x_);
        return x_.A();
      }

      // Remove all the edges from list of touching edges that are further away than 1e-14 (
      // TOPOLOGICAL_TOLERANCE ) distance from the point
      template <unsigned int dim, INPAR::CUT::CutFloattype floattype>
      void fix_distant_touching_edges(
          CORE::LINALG::Matrix<dim, 1>& p_coord, std::vector<int>& touching_edges)
      {
        bool signeddistance = false;
        double distance = 0;
        CORE::LINALG::Matrix<dim, 1> xsi;
        CORE::LINALG::Matrix<dim, numNodesEdge> xyze_edge;
        const std::vector<Edge*>& side_edges = get_side().Edges();

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
                generate_gmsh_dump();
                FOUR_C_THROW("Distance between point touching edge is too high! Check this case!");
              }
              it = touching_edges.erase(it);
            }
            else
              ++it;
          }
          else
            FOUR_C_THROW(
                "Newton did not converge for simple compute_distance between point and a line");
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
       *  ( 1 ) Multiple cut points are detected during the check_parallelism call.
       *
       *  ( 2 ) A single cut point is detected during the check_parallelism or the
       *        intersection call.
       *
       *  ( 3 ) There is no feasible cut point.
       *
       *  All feasible cut points are within the given element bounds.
       *
       *  \author hiermeier \date 08/16 */
      IntersectionStatus compute_edge_side_intersection(double& tolerance, bool check_inside = true,
          std::vector<int>* touched_edges = nullptr) override
      {
        switch (get_options_ptr()->geom_intersect_floattype())
        {
          case INPAR::CUT::floattype_cln:
          {
            return compute_edge_side_intersection_t<INPAR::CUT::floattype_cln>(
                tolerance, check_inside, touched_edges);
          }
          case INPAR::CUT::floattype_double:
          {
            return compute_edge_side_intersection_t<INPAR::CUT::floattype_double>(
                tolerance, check_inside, touched_edges);
          }
          default:
            FOUR_C_THROW("Unexpected floattype for compute_edge_side_intersection_t!");
        }
      }

      template <INPAR::CUT::CutFloattype floattype>
      IntersectionStatus compute_edge_side_intersection_t(
          double& tolerance, bool check_inside = true, std::vector<int>* touched_edges = nullptr)
      {
        check_init();
        TEUCHOS_FUNC_TIME_MONITOR("compute_edge_side_intersection");

        const bool success = check_parallelism(multiple_xsi_side_, multiple_xsi_edge_, tolerance);

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
              generate_gmsh_dump();
              FOUR_C_THROW(err_msg.str());
            }
          }
        }

        return istatus_;
      }

      /** \brief Computes the intersection points of the edge with the specified side
       *   and stores the points in cuts
       *
       *  WARNING: intersection just works for planes ( TRI3, QUAD4 unwarped! ) with lines!!!
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
       *  (5) throw FOUR_C_THROW in case this intersection wasn't treated right --> this means
       *  there is still handling of some special cases missing in the code & it does not
       *  mean that there is no intersection point.
       *
       *  \author ager */
      bool Intersect(PointSet& cuts) override;

      // Try to find possible intersection points, if this intersection is between parallell size
      // and edge without using real ComputeIntersection
      ParallelIntersectionStatus handle_parallel_intersection(
          PointSet& cuts, int id = -1, bool output = false) override;

      virtual void generate_gmsh_dump();

      // Handle cases for which normal intersection procedure did not work
      bool HandleSpecialCases() override;

      // compute intersection by splitting quad4 into two triangles
      bool triangulated_intersection(PointSet& cuts) override;

      /** \brief Will return TRUE, if local side coordinates are within the side
       *  element parameter space bounds */
      bool SurfaceWithinLimits(double tol = REFERENCETOL) const override
      {
        return CORE::GEO::CUT::KERNEL::within_limits<sidetype>(xsi_side_, tol);
      }

      /** \brief Will return TRUE, if local side coordinates are within the TRI3
       *  side element parameter space bounds */
      bool tri3_within_limits(double tol = REFERENCETOL) const
      {
        return CORE::GEO::CUT::KERNEL::within_limits<CORE::FE::CellType::tri3>(xsi_side_, tol);
      }

      /** \brief Will return TRUE, if local edge coordinate is within the line
       *  element parameter space bounds */
      bool LineWithinLimits(double tol = REFERENCETOL) const override
      {
        return CORE::GEO::CUT::KERNEL::within_limits<edgetype>(xsi_edge_, tol);
      }

      /// access the problem dimension
      unsigned prob_dim() const override { return probdim; }

      /// access the number of nodes of the side ( or 2-nd edge ) element
      unsigned num_nodes_side() const override { return numNodesSide; }

      /// access the number of nodes of the edge
      unsigned num_nodes_edge() const override { return numNodesEdge; }

      /// set the edge and side coordinates
      void set_coordinates() override;

      /// set the edge and side coordinates
      void set_coordinates(double* xyze_surfaceElement, double* xyze_lineElement) override
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
      void scale_and_shift() override
      {
        // ---------------------------------------------------------------------
        // scale the input elements if desired
        // ---------------------------------------------------------------------
        if (not is_scaled())
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
        if (not is_shifted())
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
      bool at_edge(const T& xsi)
      {
        return CORE::GEO::CUT::KERNEL::AtEdge<sidetype>(xsi);
      }

     private:
      /// Do the bounding box overlap check for the class internal edge and side variables
      bool check_bounding_box_overlap();

      /// Check if the edge \c ebb and surface \c sbb bounding boxes overlap
      bool check_bounding_box_overlap(BoundingBox& ebb, BoundingBox& sbb) const;

      /** \brief Checks the side dimension and calls the corresponding method
       *
       *  Currently surface and line elements are supported.
       *
       *  \author hiermeier \date 12/16 */
      bool check_parallelism(std::vector<CORE::LINALG::Matrix<dimside, 1>>& side_rs_intersect,
          std::vector<CORE::LINALG::Matrix<dimedge, 1>>& edge_r_intersect, double& tolerance);

      /** \brief Check if the two lines are collinear, end points are on the line, or
       *  the distance values imply that no intersection is possible
       *
       *  \author hiermeier \date 12/16 */
      bool check_collinearity(
          std::vector<CORE::LINALG::Matrix<dimside, 1>>& side_rs_corner_intersect,
          std::vector<CORE::LINALG::Matrix<dimedge, 1>>& edge_r_corner_intersect,
          double& tolerance);

      /** \brief Check the angle between two edge lines.
       *
       *  This is a quick check to skip cases which are definitely not parallel.
       *
       *  \author hiermeier \date 12/16 */
      bool check_angle_criterion_between_two_edges();

      /** ToDo This method is currently unused, since this case should be treated by
       *  the Intersect() method.
       *
       *  \author hiermeier \date 12/16 */
      bool check_parallelism_between_side_and_edge(
          std::vector<CORE::LINALG::Matrix<dimside, 1>>& side_rs_intersect,
          std::vector<CORE::LINALG::Matrix<dimedge, 1>>& edge_r_intersect, double& tolerance);

      /** \brief Check the angle between a edge and a surface normal.
       *
       *  This is a quick check to skip cases which are definitely not parallel.
       *
       *  \author hiermeier \date 12/16 */
      bool check_angle_criterion_between_side_normal_and_edge();

      /// find the local coordinate of a given edge end point ( i.e. -1 or +1 )
      bool find_local_coordinate_of_edge_end_point(
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
      bool compute_edge_tri3_intersection(int triangleid, KERNEL::PointOnSurfaceLoc& location)
      {
        switch (get_options_ptr()->geom_intersect_floattype())
        {
          case INPAR::CUT::floattype_cln:
          {
            return compute_edge_tri3_intersection_t<INPAR::CUT::floattype_cln>(
                triangleid, location);
          }
          case INPAR::CUT::floattype_double:
          {
            return compute_edge_tri3_intersection_t<INPAR::CUT::floattype_double>(
                triangleid, location);
          }
          default:
            FOUR_C_THROW("Unexpected floattype for compute_edge_tri3_intersection_t!");
        }
      }

      template <INPAR::CUT::CutFloattype floattype>
      bool compute_edge_tri3_intersection_t(int triangleid, KERNEL::PointOnSurfaceLoc& location)
      {
        if (triangleid < 0) FOUR_C_THROW("The triangle id has to be positive!");

        TEUCHOS_FUNC_TIME_MONITOR("compute_edge_tri3_intersection");
        CORE::LINALG::Matrix<3, 1> xsi;
        if (xsi_.M() != 3)
          FOUR_C_THROW("xsi_ has the wrong dimension! (dimedge + 2 = %d + 2)", dimedge);
        else
          xsi.SetView(xsi_.A());

        KERNEL::ComputeIntersection<3, edgetype, CORE::FE::CellType::tri3, floattype> ci(xsi);
        // KERNEL::DebugComputeIntersection<probdim,edgetype,CORE::FE::CellType::tri3,floattype>
        // ci( xsi
        // );

        CORE::LINALG::Matrix<3, 3> xyze_triElement;
        get_triangle(xyze_triElement, triangleid);
        CORE::LINALG::Matrix<3, numNodesEdge> xyze_lineElement(xyze_lineElement_.A(), true);

        bool conv = ci(xyze_triElement, xyze_lineElement);
        location = ci.GetSideLocation();

        return conv;
      }

      // Computes tri3 edge intersection used for quad4 -> 2 tri3 splits
      IntersectionStatus compute_edge_tri3_intersection_quad4_split(
          int triangleid, bool* close_to_shared_edge = nullptr)
      {
        switch (get_options_ptr()->geom_intersect_floattype())
        {
          case INPAR::CUT::floattype_cln:
          {
            return compute_edge_tri3_intersection_quad4_split_t<INPAR::CUT::floattype_cln>(
                triangleid, close_to_shared_edge);
          }
          case INPAR::CUT::floattype_double:
          {
            return compute_edge_tri3_intersection_quad4_split_t<INPAR::CUT::floattype_double>(
                triangleid, close_to_shared_edge);
          }
          default:
            FOUR_C_THROW("Unexpected floattype for compute_edge_tri3_intersection_quad4_split_t!");
        }
      }

      template <INPAR::CUT::CutFloattype floattype>
      IntersectionStatus compute_edge_tri3_intersection_quad4_split_t(
          int triangleid, bool* close_to_shared_edge = nullptr)
      {
        if (triangleid < 0) FOUR_C_THROW("The triangle id has to be positive!");

        TEUCHOS_FUNC_TIME_MONITOR("compute_edge_tri3_intersection");
        CORE::LINALG::Matrix<3, 1> xsi;
        if (xsi_.M() != 3)
          FOUR_C_THROW("xsi_ has the wrong dimension! (dimedge + 2 = %d + 2)", dimedge);
        else
          xsi.SetView(xsi_.A());

        KERNEL::ComputeIntersection<3, edgetype, CORE::FE::CellType::tri3, floattype> ci(xsi);
        // KERNEL::DebugComputeIntersection<probdim,edgetype,CORE::FE::CellType::tri3,floattype>
        // ci( xsi
        // );

        CORE::LINALG::Matrix<3, 3> xyze_triElement;
        get_triangle(xyze_triElement, triangleid);
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
          *close_to_shared_edge = (ci.get_side_location_triangle_split().WithinSide());

        return istatus_;
      }
      /** \brief get one of the two triangles with id 0 or 1
       *  of a QUAD4 element
       *
       *  tri3_id=0 ---> Quad4 nodes = {0 1 2}
       *  tri3_id=1 ---> Quad4 nodes = {2 3 0} */
      void get_triangle(CORE::LINALG::Matrix<3, 3>& xyze_triElement, const unsigned& tri3_id)
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
          FOUR_C_THROW("Cut::intersection::get_triangle: For Triangulation a QUAD4 is expected!");
        }
      };

      /** \brief Update compute_distance routine to get information about location from the
       * cut_kernel This function is used for normal compute_distance (without triangulation) */
      bool compute_distance(CORE::LINALG::Matrix<probdim, 1> point, double& distance,
          double& tolerance, bool& zeroarea, KERNEL::PointOnSurfaceLoc& loc,
          std::vector<int>& touched_edges, bool signeddistance = false)
      {
        switch (get_options_ptr()->geom_distance_floattype())
        {
          case INPAR::CUT::floattype_cln:
          {
            return compute_distance_t<INPAR::CUT::floattype_cln>(
                point, distance, tolerance, zeroarea, loc, touched_edges, signeddistance);
          }
          case INPAR::CUT::floattype_double:
          {
            return compute_distance_t<INPAR::CUT::floattype_double>(
                point, distance, tolerance, zeroarea, loc, touched_edges, signeddistance);
          }
          default:
            FOUR_C_THROW("Unexpected floattype for compute_distance_t!");
        }
      }

      template <INPAR::CUT::CutFloattype floattype>
      bool compute_distance_t(CORE::LINALG::Matrix<probdim, 1> point, double& distance,
          double& tolerance, bool& zeroarea, KERNEL::PointOnSurfaceLoc& loc,
          std::vector<int>& touched_edges, bool signeddistance = false)
      {
        TEUCHOS_FUNC_TIME_MONITOR("compute_distance");

        if (dimside + dimedge != probdim)
          FOUR_C_THROW(
              "This compute_distance variant won't work! Think about using "
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
        fix_distant_touching_edges<probdim, floattype>(point, touched_edges);

        return conv;
      }

      bool compute_distance(Point* p, double& distance, double& tolerance, bool& zeroarea,
          KERNEL::PointOnSurfaceLoc& loc, std::vector<int>& touched_edges,
          bool signeddistance = false)
      {
        CORE::LINALG::Matrix<probdim, 1> point(p->X());
        return compute_distance(
            point, distance, tolerance, zeroarea, loc, touched_edges, signeddistance);
      }


      // Transform IDs of the edges in the one of triangle in the splitted quad4 into the ids of
      // quad4 edges
      void get_quad_edge_ids_from_tri(std::vector<int>& quad4_touched_edges,
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

        if (quad4_touched_edges.size() > 4) FOUR_C_THROW("this should not be possible");
      }

      /// Detects in the point is close to an endpoint of the edge
      template <unsigned int numNodes, unsigned int probDim>
      bool is_close_to_endpoints(const CORE::LINALG::Matrix<probDim, numNodes>& surf,
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
      std::pair<bool, bool> is_quad4_distorted();

      /* Used during splitting of quad4 into tri3 and computing distance to them */
      bool compute_distance(CORE::LINALG::Matrix<3, 1> point, double& distance, double& tolerance,
          bool& zeroarea, KERNEL::PointOnSurfaceLoc& loc, std::vector<int>& touched_edges,
          bool signeddistance, int tri3_id, bool& extended_tri_tolerance_loc_triangle_split)
      {
        switch (get_options_ptr()->geom_distance_floattype())
        {
          case INPAR::CUT::floattype_cln:
          {
            return compute_distance_t<INPAR::CUT::floattype_cln>(point, distance, tolerance,
                zeroarea, loc, touched_edges, signeddistance, tri3_id,
                extended_tri_tolerance_loc_triangle_split);
          }
          case INPAR::CUT::floattype_double:
          {
            return compute_distance_t<INPAR::CUT::floattype_double>(point, distance, tolerance,
                zeroarea, loc, touched_edges, signeddistance, tri3_id,
                extended_tri_tolerance_loc_triangle_split);
          }
          default:
            FOUR_C_THROW("Unexpected floattype for compute_distance_t!");
        }
      }

      template <INPAR::CUT::CutFloattype floattype>
      bool compute_distance_t(CORE::LINALG::Matrix<3, 1> point, double& distance, double& tolerance,
          bool& zeroarea, KERNEL::PointOnSurfaceLoc& loc, std::vector<int>& touched_edges,
          bool signeddistance, int tri3_id, bool& extended_tri_tolerance_loc_triangle_split)
      {
        if (sidetype != CORE::FE::CellType::quad4)
          FOUR_C_THROW(
              "This compute_distance routine is only meaningful for "
              "QUAD4 side elements! But you passed in a side element "
              "of type %i | %s.",
              sidetype, CORE::FE::CellTypeToString(sidetype).c_str());

        TEUCHOS_FUNC_TIME_MONITOR("compute_distance");

        // dimension of xsi: element dimension of 2 + 1 entry for the distance
        if (xsi_.M() != 3)
          FOUR_C_THROW("xsi_ has the wrong dimension! (dimedge + 2 = %d + 2)", dimedge);
        CORE::LINALG::Matrix<3, 1> xsi(xsi_.A(), true);

        // KERNEL::DebugComputeDistance<probdim,CORE::FE::CellType::tri3, floattype>
        // cd( xsi );
        KERNEL::ComputeDistance<3, CORE::FE::CellType::tri3, floattype> cd(xsi);

        CORE::LINALG::Matrix<3, 3> xyze_triElement;
        get_triangle(xyze_triElement, tri3_id);

        bool conv = cd(xyze_triElement, point, distance, signeddistance);
        tolerance = cd.GetTolerance();
        loc = cd.GetSideLocation();
        zeroarea = cd.ZeroArea();
        std::vector<int> tri_touched_edges;
        cd.GetTouchedSideEdges(tri_touched_edges);

        get_quad_edge_ids_from_tri(touched_edges, tri_touched_edges, tri3_id);
        extended_tri_tolerance_loc_triangle_split =
            (cd.get_side_location_triangle_split().WithinSide());
        fix_distant_touching_edges<3, floattype>(point, touched_edges);

        return conv;
      }

      /// get the coordinates of the point and call the related compute_distance routine
      bool compute_distance(Point* p, double& distance, double& tolerance, bool& zeroarea,
          KERNEL::PointOnSurfaceLoc& loc, std::vector<int>& touched_edges, bool signeddistance,
          int tri3_id, bool& extended_tri_tolerance_loc_triangle_split)
      {
        CORE::LINALG::Matrix<3, 1> point(p->X());
        return compute_distance(point, distance, tolerance, zeroarea, loc, touched_edges,
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
      bool compute_cut(
          Edge* sedge, Edge* eedge, Side* side, PointSet& ee_cut_points, double& tolerance);

      /// add cut point that is a node to all edges and sides it touches
      void insert_cut(Node* n, PointSet& cuts)
      {
        cuts.insert(Point::insert_cut(get_edge_ptr(), get_side_ptr(), n));
      }

      void add_connectivity_info(Point* p, const CORE::LINALG::Matrix<probdim, 1>& xreal,
          const std::vector<int>& touched_vertices_ids, const std::vector<int>& touched_edges_ids);

      void add_connectivity_info(Point* p, const CORE::LINALG::Matrix<probdim, 1>& xreal,
          const std::vector<int>& touched_edges_ids,
          const std::set<std::pair<Side*, Edge*>>& touched_cut_pairs);

      void get_connectivity_info(const CORE::LINALG::Matrix<probdim, 1>& xreal,
          const std::vector<int>& touched_edges_ids, std::set<std::pair<Side*, Edge*>>& out);

      bool refined_bb_overlap_check(int maxstep = 10);

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
    };  // class intersection

    /*--------------------------------------------------------------------------*/
    /** \brief Create a intersection object
     *
     *  \author hiermeier \date 12/16 */
    class IntersectionFactory
    {
     public:
      IntersectionFactory(){};

      Teuchos::RCP<IntersectionBase> create_intersection(
          CORE::FE::CellType edge_type, CORE::FE::CellType side_type) const;

     private:
      template <CORE::FE::CellType edgeType>
      IntersectionBase* create_intersection(CORE::FE::CellType side_type, int probdim) const
      {
        switch (side_type)
        {
          case CORE::FE::CellType::quad4:
            return create_concrete_intersection<edgeType, CORE::FE::CellType::quad4>(probdim);
          case CORE::FE::CellType::quad8:
            return create_concrete_intersection<edgeType, CORE::FE::CellType::quad8>(probdim);
          case CORE::FE::CellType::quad9:
            return create_concrete_intersection<edgeType, CORE::FE::CellType::quad9>(probdim);
          case CORE::FE::CellType::tri3:
            return create_concrete_intersection<edgeType, CORE::FE::CellType::tri3>(probdim);
          case CORE::FE::CellType::line2:
            return create_concrete_intersection<edgeType, CORE::FE::CellType::line2>(probdim);
          default:
            FOUR_C_THROW(
                "Unsupported SideType! If meaningful, add your sideType here. \n"
                "Given SideType = %s",
                CORE::FE::CellTypeToString(side_type).c_str());
            exit(EXIT_FAILURE);
        }
        exit(EXIT_FAILURE);
      }

      template <CORE::FE::CellType edgeType, CORE::FE::CellType sideType>
      IntersectionBase* create_concrete_intersection(const int& probdim) const
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
            FOUR_C_THROW("Unsupported ProbDim! ( probdim = %d )", probdim);
            exit(EXIT_FAILURE);
        }
        return inter_ptr;
      };

    };  // class IntersectionFactory

  }  // namespace CUT
}  // namespace CORE::GEO

// static members of intersection base class
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

FOUR_C_NAMESPACE_CLOSE

#endif

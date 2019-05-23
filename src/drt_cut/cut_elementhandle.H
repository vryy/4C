/*----------------------------------------------------------------------*/
/*!
\brief Outside world interface to element. Converts quadratic to linear element. This provides the
  Gaussian rules generated from the cut

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249

\level 2
*/
/*----------------------------------------------------------------------*/

#ifndef CUT_ELEMENTHANDLE_H
#define CUT_ELEMENTHANDLE_H


#include "cut_element.H"

#include "../../src/drt_fem_general/drt_utils_gausspoints.H"


namespace GEO
{
  namespace CUT
  {
    class Point;
    class Node;
    class Element;
    class Mesh;


    /*!
    \brief Outside world interface to element. Converts quadratic to linear element. This provides
    the Gaussian rules generated from the cut
     */
    //--------------------------------------------------------------------------//
    // the ElementHandle base class
    //--------------------------------------------------------------------------//
    class ElementHandle
    {
     public:
      virtual ~ElementHandle() {}

      //--------------------------------------------------------------------------//
      //! @name access to cut/intersection status of the element

      /*!
     \brief returns true in case that any cut-side cut with the element produces cut points,
            i.e. also for touched cases (at points, edges or sides),
            or when an element side has more than one facet or is touched by fully/partially
            by the cut side
       */
      virtual bool IsCut() = 0;

      /// return true if the element has more than one volume-cell and therefore is intersected by a
      /// cut-side
      virtual bool IsIntersected() = 0;

      //@}

      //--------------------------------------------------------------------------//
      //! @name access to basic element data

      //! get the shape of the element
      virtual DRT::Element::DiscretizationType Shape() = 0;

      //! get the nodes of the element
      virtual const std::vector<Node*>& Nodes() = 0;

      //! compute local coordinates of the element for a given point
      const LINALG::Matrix<3, 1>& LocalCoordinates(Point* p)
      {
        std::map<Point*, LINALG::Matrix<3, 1>>::iterator i = local_coordinates_.find(p);
        if (i != local_coordinates_.end())
        {
          return i->second;
        }
        LINALG::Matrix<3, 1>& rst = local_coordinates_[p];
        LINALG::Matrix<3, 1> xyz;
        p->Coordinates(xyz.A());
        LocalCoordinates(xyz, rst);
        return rst;
      }

      //! compute local coordinates of the element for given global coordinates
      virtual void LocalCoordinates(const LINALG::Matrix<3, 1>& xyz, LINALG::Matrix<3, 1>& rst) = 0;


      //@}

      //--------------------------------------------------------------------------//
      //! @name access to the element's GEO::CUT::Elements (sub-elements)

      virtual void CollectElements(plain_element_set& elements) = 0;

      //@}

      //--------------------------------------------------------------------------//
      //! @name access to the element's volume-cells

      //! Collect all volume-cells belonging to this elements
      virtual void GetVolumeCells(plain_volumecell_set& cells) = 0;

      //! Collect all volume-cells belonging to this element ordered by position
      virtual void CollectVolumeCells(
          plain_volumecell_set& cells_inside, plain_volumecell_set& cells_outside) = 0;

      //@}

      //--------------------------------------------------------------------------//
      //! @name access to the element's volume integration specific quantities

      //! get the element's volumetric integration cells (just for Tessellation)
      virtual void GetIntegrationCells(plain_integrationcell_set& cells) = 0;


      /*!
       \brief  Collect the Gaussian points of all volume-cells belonging to this element in such a
       way that Gaussian rule for every volume-cell can be separated
        */
      void VolumeCellGaussPoints(
          plain_volumecell_set& cells, std::vector<DRT::UTILS::GaussIntegration>& intpoints);

      void AppendVolumeCellGaussPoints_Tessellation(
          Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc, GEO::CUT::VolumeCell* vc);

      void AppendVolumeCellGaussPoints_MomentFitting(
          Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc, GEO::CUT::VolumeCell* vc);

      void AppendVolumeCellGaussPoints_DirectDivergence(
          Teuchos::RCP<DRT::UTILS::GaussPointsComposite> gpc, GEO::CUT::VolumeCell* vc);


      /*!
       \brief Collect the Gaussian points of all the volume-cells belonging to this element.
              The integration rules over all the volume-cells are connected.
       */
      Teuchos::RCP<DRT::UTILS::GaussPointsComposite> GaussPointsConnected(
          plain_volumecell_set& cells, INPAR::CUT::VCellGaussPts gausstype);


      //! ...
      virtual int NumVolumeCellGaussPoints() = 0;

      //@}

      //--------------------------------------------------------------------------//
      //! @name access to the element's interface integration specific quantities

      //! get the element's boundary integration cells
      // TODO note which ones???
      virtual void GetBoundaryCells(plain_boundarycell_set& bcells) = 0;

#if (0)
      /*!
      \brief Unused
       */
      void BoundaryCellGaussPoints(
          const std::map<int, std::vector<GEO::CUT::BoundaryCell*>>& bcells,
          std::map<int, std::vector<DRT::UTILS::GaussIntegration>>& intpoints);

      /*!
      \brief unused: old implementation of boundarycell Gauss points collection
       */
      void BoundaryCellGaussPoints(MeshIntersection& mesh, int mi,
          const std::map<int, std::vector<GEO::CUT::BoundaryCell*>>& bcells,
          std::map<int, std::vector<DRT::UTILS::GaussIntegration>>& intpoints);
#endif

      /*!
      \brief Collect the Gauss points of all the boundary-cells belong to this element.
      This is the method used now in the new implementation
       */
      void BoundaryCellGaussPointsLin(
          const std::map<int, std::vector<GEO::CUT::BoundaryCell*>>& bcells,
          std::map<int, std::vector<DRT::UTILS::GaussIntegration>>& intpoints);

      //@}


      //--------------------------------------------------------------------------//
      //! @name access to the element's sets of volume-cells, nds-vectors and integration points

      //! get all the element' sets of volume-cells, nds-vectors and integration points
      virtual bool GetCellSets_DofSets_GaussPoints(std::vector<plain_volumecell_set>& cell_sets,
          std::vector<std::vector<int>>& nds_sets,
          std::vector<std::vector<DRT::UTILS::GaussIntegration>>& intpoints_sets,
          bool include_inner);


      //! get the element's sets of volume-cells ordered by inside/outside position
      virtual void VolumeCellSets() = 0;

      //! get all the element' sets of volume-cells and nds-vectors
      void GetVolumeCellsDofSets(std::vector<plain_volumecell_set>& cellsets,
          std::vector<std::vector<int>>& nds_sets, bool include_inner);

      void GetBoundaryCellSets(const std::vector<GEO::CUT::Point::PointPosition>& desired_positions,
          std::vector<plain_boundarycell_set>& bcellsets);

      void GetBoundaryCellSets(GEO::CUT::Point::PointPosition desired_position,
          std::vector<plain_boundarycell_set>& bcellsets)
      {
        const std::vector<GEO::CUT::Point::PointPosition> desired_positions(1, desired_position);
        GetBoundaryCellSets(desired_positions, bcellsets);
      }

     protected:
      /// fill the derived class internal boundary cell sets
      virtual void BoundaryCellSet(Point::PointPosition position) = 0;

      /// access the class internal boundary cells sets
      virtual const std::vector<plain_boundarycell_set>& GetBoundaryCellSet(
          Point::PointPosition position) = 0;

     public:
      //! get the element's sets of volume-cells with inside position
      virtual const std::vector<plain_volumecell_set>& GetVcSetsInside() = 0;

      //! get the element's sets of volume-cells with outside position
      virtual const std::vector<plain_volumecell_set>& GetVcSetsOutside() = 0;

      //! get the element's sets of nodal dofset vectors (nds-vectors) with inside position
      std::vector<std::vector<int>>& GetNodalDofSet_VcSets_Inside()
      {
        return nodaldofset_vc_sets_inside_;
      };

      //! get the element's sets of nodal dofset vectors (nds-vectors) with outside position
      std::vector<std::vector<int>>& GetNodalDofSet_VcSets_Outside()
      {
        return nodaldofset_vc_sets_outside_;
      };

      //! ...
      std::vector<std::map<int, int>>& Get_NodeDofsetMap_VcSets_Inside_forCommunication()
      {
        return vcsets_nid_dofsetnumber_toComm_inside_;
      };

      //! ...
      std::vector<std::map<int, int>>& Get_NodeDofsetMap_VcSets_Outside_forCommunication()
      {
        return vcsets_nid_dofsetnumber_toComm_outside_;
      };

      //@}


     private:
      /*!
       \brief Project the integration rule available in the local coordinates of the
       integation-cells to the local coordinates of background element
       */
      template <DRT::Element::DiscretizationType distype>
      Teuchos::RCP<DRT::UTILS::GaussPoints> CreateProjected(
          const std::vector<GEO::CUT::Point*>& cpoints,
          Teuchos::RCP<DRT::UTILS::GaussPoints> gp_ic);

#if 0
  template <DRT::Element::DiscretizationType distype>
  Teuchos::RCP<DRT::UTILS::GaussPoints> TransformVolumeCellGaussPoints( GEO::CUT::IntegrationCell * ic );
#endif

      std::map<Point*, LINALG::Matrix<3, 1>> local_coordinates_;

     protected:
      /// dof set number of all element nodes, contains the dofset numbers for all nodes of the
      /// superior element (i.e. 20 for hex20 superior element)
      std::vector<std::vector<int>> nodaldofset_vc_sets_inside_;
      std::vector<std::vector<int>> nodaldofset_vc_sets_outside_;

      /// for each set of volumecells a map containing nids and dofsetnumbers that has to be
      /// communicated for the superior element (i.e. 20 for hex20 superior element)
      std::vector<std::map<int, int>> vcsets_nid_dofsetnumber_toComm_inside_;
      std::vector<std::map<int, int>> vcsets_nid_dofsetnumber_toComm_outside_;
    };



    //--------------------------------------------------------------------------//
    // the linear ElementHandle for linear elements
    //--------------------------------------------------------------------------//
    /// linear element handle
    class LinearElementHandle : public ElementHandle
    {
     public:
      LinearElementHandle() : element_(NULL), cells_set_(false), bcell_sets_(0) {}

      explicit LinearElementHandle(Element* e) : element_(e), cells_set_(false), bcell_sets_(0)
      {
        // set also the parent element Id which is trivial the same as the element Id
        element_->ParentId(e->Id());
      }

      //--------------------------------------------------------------------------//
      //! @name access to cut/intersection status of the element

      /*!
       \brief returns true in case that any cut-side cut with the element produces cut points,
              i.e. also for touched cases (at points, edges or sides),
              or when an element side has more than one facet or is touched by fully/partially by
       the cut side
       */
      virtual bool IsCut() { return element_->IsCut(); }

      /// return true if the element has more than one volume-cell and therefore is intersected by a
      /// cut-side
      virtual bool IsIntersected() { return element_->IsIntersected(); }

      //@}


      //--------------------------------------------------------------------------//
      //! @name access to basic element data

      //! get the shape of the element
      virtual DRT::Element::DiscretizationType Shape() { return element_->Shape(); }

      //! get the nodes of the element
      virtual const std::vector<Node*>& Nodes() { return element_->Nodes(); }

      //! compute local coordinates of the element for given global coordinates
      virtual void LocalCoordinates(const LINALG::Matrix<3, 1>& xyz, LINALG::Matrix<3, 1>& rst)
      {
        element_->LocalCoordinates(xyz, rst);
      }

      //@}


      //--------------------------------------------------------------------------//
      //! @name access to the element's GEO::CUT::Elements (subelements)

      //! collect all sub-elements
      virtual void CollectElements(plain_element_set& elements) { elements.insert(element_); }
      //@}

      //--------------------------------------------------------------------------//
      //! @name access to the element's volume-cells

      //! Collect all volume-cells belonging to this elements
      virtual void GetVolumeCells(plain_volumecell_set& cells);

      //! Collect all volume-cells belonging to this element ordered by position
      virtual void CollectVolumeCells(
          plain_volumecell_set& cells_inside, plain_volumecell_set& cells_outside);

      //@}

      //--------------------------------------------------------------------------//
      //! @name access to the element's volume integration specific quantities

      //! get the element's volumetric integration cells (just for Tessellation)
      virtual void GetIntegrationCells(plain_integrationcell_set& cells)
      {
        element_->GetIntegrationCells(cells);
      }

      //! ...
      virtual int NumVolumeCellGaussPoints() { return element_->NumGaussPoints(Shape()); }

      //@}



      //--------------------------------------------------------------------------//
      //! @name access to the element's interface integration specific quantities

      //! get the element's boundary integration cells
      // TODO note which ones???
      virtual void GetBoundaryCells(plain_boundarycell_set& bcells)
      {
        dserror(
            "Deprecated function! Use the GetBoundaryCellSet( Point::PointPosition ), "
            "instead!");
        element_->GetBoundaryCells(bcells);
      }

     protected:
      /** \brief Get the boundary cell sets with given position
       *
       *  \param position (in) : desired position of the boundary cell
       *
       *  \author hiermeier \date 01/17 */
      virtual const std::vector<plain_boundarycell_set>& GetBoundaryCellSet(
          Point::PointPosition position)
      {
        BoundaryCellSet(position);
        return bcell_sets_.at(position);
      }

      /** \brief Fill the class internal boundary sets for the given position
       *
       *  This function does a direct return, if the boundary cells have been already
       *  extracted, otherwise we extract the boundary cells from the volume cells
       *  and add them.
       *
       *  \param position (in) : desired position of the boundary cell
       *
       *  \author hiermeier \date 01/17 */
      virtual void BoundaryCellSet(Point::PointPosition position);

      //@}

     public:
      //--------------------------------------------------------------------------//
      //! @name access to the element's sets of volume-cells, nds-vectors and integration points

      //! get the element's sets of volume-cells ordered by inside/outside position
      virtual void VolumeCellSets();

      //! get the element's sets of volume-cells with inside position
      const std::vector<plain_volumecell_set>& GetVcSetsInside()
      {
        if (!cells_set_) VolumeCellSets();  // build the volume cell sets
        return vc_sets_inside_;
      };

      //! get the element's sets of volume-cells with outside position
      const std::vector<plain_volumecell_set>& GetVcSetsOutside()
      {
        if (!cells_set_) VolumeCellSets();  // build the volume cell sets
        return vc_sets_outside_;
      };

      //@}


     private:
      Element* element_;  ///< the unique element in the cut library


     protected:
      bool cells_set_;  ///< have the cells already set?
      std::vector<plain_volumecell_set>
          vc_sets_inside_;  ///< connected sets of volume-cells with inside position
      std::vector<plain_volumecell_set>
          vc_sets_outside_;  ///< connected sets of volume-cells with outside position


      GEN::pairedvector<Point::PointPosition, std::vector<plain_boundarycell_set>> bcell_sets_;
    };



    //--------------------------------------------------------------------------//
    // the quadratic ElementHandle base class for quadratic elements
    //--------------------------------------------------------------------------//
    /// quadratic element handle
    class QuadraticElementHandle : public ElementHandle
    {
     public:
      QuadraticElementHandle() : cells_connected_(false), connected_bcell_sets_(0) {}

      //--------------------------------------------------------------------------//
      //! @name access to cut/intersection status of the element

      /*!
        \brief returns true in case that any cut-side cut with the element produces cut points,
               i.e. also for touched cases (at points, edges or sides),
               or when an element side has more than one facet or is touched by fully/partially by
        the cut side
       */
      virtual bool IsCut();

      /// return true if the element has more than one volume-cell and therefore is intersected by a
      /// cut-side
      virtual bool IsIntersected();

      //@}


      //--------------------------------------------------------------------------//
      //! @name access to basic element data

      //! get the nodes of the element
      virtual const std::vector<Node*>& Nodes() { return nodes_; }

      //@}


      //--------------------------------------------------------------------------//
      //! @name access to the element's GEO::CUT::Elements (subelements)

      //! collect all sub-elements
      virtual void CollectElements(plain_element_set& elements)
      {
        std::copy(
            subelements_.begin(), subelements_.end(), std::inserter(elements, elements.begin()));
      }
      //@}

      //--------------------------------------------------------------------------//
      //! @name access to the element's volume-cells

      //! Collect all volume-cells belonging to this elements
      virtual void GetVolumeCells(plain_volumecell_set& cells);

      //! Collect all volume-cells belonging to this element ordered by position
      virtual void CollectVolumeCells(
          plain_volumecell_set& cells_inside, plain_volumecell_set& cells_outside);

      //! Collect volume-cells belonging to this element with the given position
      void CollectVolumeCells(
          Point::PointPosition position, plain_volumecell_set& evolcells_position) const;


      //@}

      //--------------------------------------------------------------------------//
      //! @name access to the element's volume integration specific quantities

      //! get the quadratic element's volumetric integration cells (just for Tessellation)
      virtual void GetIntegrationCells(plain_integrationcell_set& cells);

      //! ...
      virtual int NumVolumeCellGaussPoints()
      {
        int numgp = 0;
        for (std::vector<Element*>::iterator i = subelements_.begin(); i != subelements_.end(); ++i)
        {
          Element* e = *i;
          numgp += e->NumGaussPoints(Shape());
        }
        return numgp;
      }

      //@}



      //--------------------------------------------------------------------------//
      //! @name access to the element's interface integration specific quantities

      /// @{

      //! get all the quadratic element's boundary integration cells
      // TODO note which ones???
      virtual void GetBoundaryCells(plain_boundarycell_set& bcells);

     protected:
      /** \brief Fill the class internal connected boundary sets for the given position
       *
       *  This function does a direct return, if the connected boundary cells have been
       *  already created, otherwise we extract the boundary cells from the connected
       *  volume cells and add them.
       *
       *  \param position (in) : desired position of the boundary cell
       *
       *  \author hiermeier \date 01/17 */
      virtual void BoundaryCellSet(Point::PointPosition position);

      /** \brief Get the connected boundary cell sets with given position
       *  ( read access )
       *
       *  \param position (in) : desired position of the boundary cell
       *
       *  \author hiermeier \date 01/17 */
      virtual const std::vector<plain_boundarycell_set>& GetBoundaryCellSet(
          Point::PointPosition position)
      {
        BoundaryCellSet(position);
        return connected_bcell_sets_.at(position);
      }

     private:
      /** \brief Connect volume cells first and start to create the corresponding
       *  boundary cells afterwards.
       *
       *  \param position (in) : desired position of the boundary cell
       *
       *  \author hiermeier \date 01/17 */
      void ConnectBoundaryCells(Point::PointPosition position);

      /** \brief After the volume cells have been connected, the connected boundary
       *  cells are build
       *
       *  \param connected_vcell_set (in)  : connected volume cell sets
       *  \param connected_bcell_set (out) : connected boundary cell sets
       *
       *  \author hiermeier \date 01/17 */
      void BuildBoundaryCellSets(const std::vector<plain_volumecell_set>& connected_vcell_set,
          std::vector<plain_boundarycell_set>& connected_bcell_set) const;

      /// @}

     public:
      //--------------------------------------------------------------------------//
      //! @name access to the element's sets of volume-cells, nds-vectors and integration points
      /// @{
      //! get the element's sets of volume-cells ordered by inside/outside position
      virtual void VolumeCellSets();

      //! get connections/sets of volume-cells between sub-elements ordered by inside position
      const std::vector<plain_volumecell_set>& GetVcSetsInside()
      {
        if (!cells_connected_) VolumeCellSets();
        return connected_vc_sets_inside_;
      };

      //! get connections/sets of volume-cells between sub-elements ordered by outside position
      const std::vector<plain_volumecell_set>& GetVcSetsOutside()
      {
        if (!cells_connected_) VolumeCellSets();
        return connected_vc_sets_outside_;
      };

      /// @}


      //--------------------------------------------------------------------------//
      //! @name specific routines to connect volume-cells between sub-elements which belong to one
      //! connection

      //! connect volume-cells to sets of volume-cells
      virtual void ConnectVolumeCells();

      //! build sets
      void BuildCellSets(plain_volumecell_set& cells_to_connect,
          std::vector<plain_volumecell_set>& connected_sets);

      //@}


     protected:
      bool
          cells_connected_;  ///< have the volume-cells been already connected between sub-elements?
      std::vector<plain_volumecell_set>
          connected_vc_sets_inside_;  ///< connected volume-cells with inside position
      std::vector<plain_volumecell_set>
          connected_vc_sets_outside_;  ///< connected volume-cells with outside position

      GEN::pairedvector<Point::PointPosition, std::vector<plain_boundarycell_set>>
          connected_bcell_sets_;

      std::vector<Element*> subelements_;  ///< the quadratic element's linear sub-elements
      std::vector<Node*> nodes_;           ///< the quadratic element's nodes
    };

    /// hex20 element handle
    class Hex20ElementHandle : public QuadraticElementHandle
    {
     public:
      //! constructor
      Hex20ElementHandle(Mesh& mesh, int eid, const std::vector<int>& nids);

      //! get the shape of the element
      virtual DRT::Element::DiscretizationType Shape() { return DRT::Element::hex20; }

      //! compute local coordinates of the element for given global coordinates
      virtual void LocalCoordinates(const LINALG::Matrix<3, 1>& xyz, LINALG::Matrix<3, 1>& rst);
    };

    /// hex27 element handle
    class Hex27ElementHandle : public QuadraticElementHandle
    {
     public:
      //! constructor
      Hex27ElementHandle(Mesh& mesh, int eid, const std::vector<int>& nids);

      //! get the shape of the element
      virtual DRT::Element::DiscretizationType Shape() { return DRT::Element::hex27; }

      //! compute local coordinates of the element for given global coordinates
      virtual void LocalCoordinates(const LINALG::Matrix<3, 1>& xyz, LINALG::Matrix<3, 1>& rst);
    };

    /// tet10 element handle
    class Tet10ElementHandle : public QuadraticElementHandle
    {
     public:
      //! constructor
      Tet10ElementHandle(Mesh& mesh, int eid, const std::vector<int>& nids);

      //! get the shape of the element
      virtual DRT::Element::DiscretizationType Shape() { return DRT::Element::tet10; }

      //! compute local coordinates of the element for given global coordinates
      virtual void LocalCoordinates(const LINALG::Matrix<3, 1>& xyz, LINALG::Matrix<3, 1>& rst);
    };

    /// wedge15 element handle
    class Wedge15ElementHandle : public QuadraticElementHandle
    {
     public:
      //! constructor
      Wedge15ElementHandle(Mesh& mesh, int eid, const std::vector<int>& nids);

      //! get the shape of the element
      virtual DRT::Element::DiscretizationType Shape() { return DRT::Element::wedge15; }

      //! compute local coordinates of the element for given global coordinates
      virtual void LocalCoordinates(const LINALG::Matrix<3, 1>& xyz, LINALG::Matrix<3, 1>& rst);
    };

  }  // namespace CUT
}  // namespace GEO

#endif

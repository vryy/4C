/*---------------------------------------------------------------------*/
/*!
\file facet_integration.H

\brief Integrates base functions over the facet for both volumecell facets and for boundarycells
equations

\level 2

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifndef FACET_INTEGRATION_H
#define FACET_INTEGRATION_H

#include "cut_element.H"
#include "cut_mesh.H"
#include "cut_enum.H"

//#define DIRECTDIV_EXTENDED_DEBUG_OUTPUT
//#define TRIANGULATE_ALL_FACETS_FOR_DIVERGENCECELLS

namespace GEO
{
  namespace CUT
  {
    class Element;
    class Facet;

    /*!
    \brief This class performs the integration of base functions over the facet. The points of the facet
     should be arranged in anti-clockwise manner when looking the facet away from the volume
    this ensures outward normal vector when divergence theorem is used
    */
    class FacetIntegration
    {
      public:
        FacetIntegration(Facet *face1,Element *element1,const GEO::CUT::Point::PointPosition posi,bool bcellInt,bool global)
        :face1_(face1),               //facet under consideration
         elem1_(element1),            //the element for which the facet is a part of
         position_(posi),             //position
         bcellInt_(bcellInt),         //"true" if it is boundarycell integration
         global_(global),             //"true" for performing integration in global coordinates
         orderingComputed_(false)     //whether cw or acw ordering of vertices computed
        {
        }

        ~FacetIntegration(){};

        /*!
        \brief Select the base function to be integrated
        */
        void set_integ_number(int inte_num){inte_num_ = inte_num;}

        /*!
        \brief Performs the integration of a function over the facet
        */
        double integrate_facet();

        /*!
        \brief Computes the equation of the plane that contains this facet
        */
        std::vector<double> equation_plane(const std::vector<std::vector<double> > & cornersLocal);

        /*!
        \brief Returns the equation of plane that contains this facet
        */
        std::vector<double> get_equation() {return eqn_plane_;}

        /*!
        \brief Return whether the vertices numbering of the facet is clockwise
        */
        bool IsClockwiseOrdering();

        /*!
        \brief Generate Gaussian points over the considered facet by triangulating it. This is used when
        DirectDivergence option is used for Gauss point generation
        */
        void DivergenceIntegrationRule( Mesh &mesh,
            Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> & cgp );

        /*!
        \brief Generate Gaussian points over the considered facet by triangulating it. This is used when
        DirectDivergence option is used for Gauss point generation
        */
        void DivergenceIntegrationRuleNew( Mesh &mesh,
            Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> & cgp );

      private:
        /*!
        \brief Check whether the vertices numbering of the facet is clockwise
        */
        void IsClockwise(const std::vector<double> & eqn_plane,const std::vector<std::vector<double> > & cornersLocal);

        /*
        \brief Compute the function which replaces "x" when projecting the facet into coordinate plane
        */
        std::vector<double> compute_alpha(std::vector<double> & eqn_plane,GEO::CUT::ProjectionDirection intType);

        /*!
        \brief Get normal of the considered facet in a particular coordinate direction defined by intType
        */
        double getNormal(GEO::CUT::ProjectionDirection intType);

        /*!
        \brief Perform integration of base functions over boundarycells
        */
        void BoundaryFacetIntegration(const std::vector<std::vector<double> > & cornersLocal,double& facet_integ,
                                      GEO::CUT::ProjectionDirection intType);

        /*!
        \brief Generate boundary cells for the considered facet. May need to perform triangulatio
        */
        void GenerateDivergenceCells(bool divergenceRule,
                                     Mesh &mesh,
                                     std::list<Teuchos::RCP<BoundaryCell> > & divCells);

        /*!
        \brief Generate boundary cells for the considered facet. May need to perform triangulatio
        */
        void GenerateDivergenceCellsNew(bool                                   divergenceRule,
                                     Mesh                                   &mesh,
                                     std::list<Teuchos::RCP<BoundaryCell> > &divCells,
                                     const std::vector<Point*>              &cornersGlobal);


        /*!
        \brief Temporarily create Tri3 cell. This is not stored in Mesh
        */
        void TemporaryTri3(const std::vector<Point*>& corners, std::list<Teuchos::RCP<BoundaryCell> >& divCells);

        /*!
        \brief Temporarily create Quad4 cell. This is not stored in Mesh
        */
        void TemporaryQuad4(const std::vector<Point*>& corners, std::list<Teuchos::RCP<BoundaryCell> >& divCells);

        /*!
        \brief Check whether the facet area calculated by splitting and triangulation are the same
         */
        void DebugAreaCheck( std::list<Teuchos::RCP<BoundaryCell> > & divCells,
                             std::string alreadyDone,
                             Mesh &mesh );

        //! considered facet
        Facet *face1_;

        //! background element which was cut to produce this facet
        Element *elem1_;

        //! position of the facet
        const GEO::CUT::Point::PointPosition position_;

        //! True for boundarycell integration
        bool bcellInt_;

        //! True if performing integration in global coordinates
        bool global_;

        //! set the base function to be integrated
        int inte_num_;

        //! True if nodes of facet are arranged to give inward normal
        bool clockwise_;

        //! whether the clockwise or ACW ordering of facet computed already
        bool orderingComputed_;

        //! equation of plane that contains the facet
        std::vector<double> eqn_plane_;

        std::list<Teuchos::RCP<BoundaryCell> > boundarycells_;

    };
  }
}

#endif

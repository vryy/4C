/*---------------------------------------------------------------------*/
/*!
\file boundarycell_integration.H

\brief used in boundary cell integration

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#ifndef BOUNDARYCELL_INTEGRATION_H
#define BOUNDARYCELL_INTEGRATION_H

#include "facet_integration.H"

#include "../linalg/linalg_serialdensematrix.H"

namespace GEO
{
  namespace CUT
  {
    class Element;
    class Facet;

    /*
    \brief This class integrates the base functions used in the moment fitting equations over the cut facets of the volumecell
     */
    class BoundarycellIntegration
    {
      public:
      /*!
       \brief Constructor
       */
        BoundarycellIntegration(Element* elem, Facet* bcell, const GEO::CUT::Point::PointPosition posi, int num_func)
          :elem1_(elem),
          bcell_(bcell),
          position_(posi),
          num_func_(num_func)
        {
        }

        /*!
         \brief Destructor
         */
        ~BoundarycellIntegration(){};

        /*!
        \brief Generate integration rule for the considered boundarycell.
        Unlike facet integration facets, whose x-component of normal is zero, cannot be eliminated from the integration.
         */
        Epetra_SerialDenseVector GenerateBoundaryCellIntegrationRule();

        /*!
        \brief Returns the location of Gauss points over the boundarycell
         */
        std::vector<std::vector<double> > getBcellGaussPointLocation(){return BcellgausPts_;}

      private:
        Element* elem1_;
        Facet* bcell_;
        const GEO::CUT::Point::PointPosition position_;
        int num_func_;
        std::vector<std::vector<double> > BcellgausPts_;

        /*!
        \brief Distribute the Gauss points over the boundarycell
         */
        void DistributeBoundaryCellGaussPoints(std::vector<double> eqn, std::vector<std::vector<double> > corners,
                            std::vector<std::vector<double> >& bcGausspts,int ptNos);

        /*!
        \brief Moment fitting matrix is formed for the boundarycell integration
        */
        void momentFittingMatrix(std::vector<std::vector<double> >&mom,std::vector<std::vector<double> >gauspts);

        /*!
        \brief The geometry of boundarycell and the location of Gaussian points are written in GMSH output file for visualization
        and for debugging
        */
        void BcellGaussPointGmsh(const std::vector<std::vector<double> > bcGausspts, const std::vector<std::vector<double> > corners);
    };
  }
}

#endif

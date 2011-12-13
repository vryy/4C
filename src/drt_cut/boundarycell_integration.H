#ifndef BOUNDARYCELL_INTEGRATION_H
#define BOUNDARYCELL_INTEGRATION_H

#include<vector>
#include "facet_integration.H"

#include "../linalg/linalg_serialdensematrix.H"

namespace GEO
{
  namespace CUT
  {
    class Element;
    class Facet;
    class BoundarycellIntegration
    {
      public:
        BoundarycellIntegration(Element* elem, Facet* bcell, const GEO::CUT::Point::PointPosition posi, int num_func)
          :elem1_(elem),
          bcell_(bcell),
          position_(posi),
          num_func_(num_func)
        {
        }

        ~BoundarycellIntegration(){};
        Epetra_SerialDenseVector GenerateBoundaryCellIntegrationRule();
        std::vector<std::vector<double> > getBcellGaussPointLocation(){return BcellgausPts_;}

      private:
        Element* elem1_;
        Facet* bcell_;
        const GEO::CUT::Point::PointPosition position_;
        int num_func_;
        std::vector<std::vector<double> > BcellgausPts_;
        void DistributeBoundaryCellGaussPoints(std::vector<double> eqn, std::vector<std::vector<double> > corners,
                            std::vector<std::vector<double> >& bcGausspts,int ptNos);
        void momentFittingMatrix(std::vector<std::vector<double> >&mom,std::vector<std::vector<double> >gauspts);
        void BcellGaussPointGmsh(const std::vector<std::vector<double> > bcGausspts, const std::vector<std::vector<double> > corners);
    };
  }
}

#endif

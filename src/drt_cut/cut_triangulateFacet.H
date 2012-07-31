#ifndef CUT_TRIANGULATEFACET_H
#define CUT_TRIANGULATEFACET_H

//#include "cut_facet.H"


namespace GEO
{
namespace CUT
{
class Point;
class Facet;

/*!
\brief A class to split a facet into tri and quad cells
 */
class TriangulateFacet
{
public:
  /*!
  \brief Constructor
   */
  TriangulateFacet( Facet* face, std::vector<Point*> ptlist )
    :facet1_(face),
     ptlist_(ptlist)
    {
    }

  /*!
  \brief Destructor
   */
  ~TriangulateFacet(){};

  /*!
  \brief Split the facet into appropriate number of tri and quad
   */
  void SplitFacet();

  /*!
  \brief A general facet is triangulated with ear clipping method.
  When triOnly=true calls conventional Earclipping method. Otherwise it creates both Tri and Quad cells
  to reduce the number of Gaussian points
   */
  void EarClipping( std::vector<int> ptConcavity,
                    bool triOnly=false  );          //create triangles only?

  /*!
  \brief Returns Tri and Quad cells that are created by facet splitting
   */
  std::vector<std::vector<Point*> > GetSplitCells(){return split_;}

private:

  /*!
  \brief Split a concave 4 noded facet into a 2 tri
  */
  void Split4nodeFacet( std::vector<Point*> &poly, bool callFromSplitAnyFacet=false );

  /*!
  \brief Split a convex facet or a facet with only one concave point into 1 Tri and few Quad cells
   */
  void SplitConvex_1ptConcave_Facet( std::vector<int> ptConcavity );

  /*!
  \brief A concave facet which has more than 2 concavity points are split into appropriate cells
  */
  void SplitGeneralFacet( std::vector<int> ptConcavity );

  /*!
  \brief check whether any two adjacent polygonal points are concave
   */
  bool HasTwoContinuousConcavePts( std::vector<int> ptConcavity );

  Facet* facet1_;
  std::vector<Point*> ptlist_;
  std::vector<std::vector<Point*> > split_;
};
}
}

#endif

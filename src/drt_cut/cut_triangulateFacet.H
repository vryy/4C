#ifndef CUT_TRIANGULATEFACET_H
#define CUT_TRIANGULATEFACET_H

#include "cut_facet.H"


namespace GEO
{
namespace CUT
{
class Point;

/*!
\brief A class to split a facet into tri and quad cells
 */
class TriangulateFacet
{
public:
  /*!
  \brief Constructor
   */
  TriangulateFacet( Facet* face, Mesh& mesh, std::vector<Point*> ptlist )
    :facet1_(face),
     mesh_(mesh),
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

  std::vector<std::vector<Point*> > GetSplitCells(){return split_;}

  std::vector<std::vector<Point*> > split_;

private:

  /*!
  \brief Check whether more than two points fall on a line, and delete all points other than the end points
   */
  void DeleteInlinePts();

  /*!
  \brief Split a 5 noded facet into a tri and a quad
  */
  void Split5nodeFacet();

  /*!
  \brief Split a 6 noded facet into 2 quad
   */
  void Split6nodeFacet();

  /*!
  \brief Split a 7 noded facet into 2 quad and a tri
   */
  void Split7nodeFacet();

  /*!
  \brief Split a 8 noded facet into 2 quad and a tri
   */
  void Split8nodeFacet();

  /*!
  \brief Split a 9 noded facet into 2 quad and a tri
   */
  void Split9nodeFacet();

  /*!
  \brief Split a 10 noded facet into 4 quad
   */
  void Split10nodeFacet();

  /*!
  \brief A concave facet which has more than 2 concavity points are split into appropriate cells. This also
  splits facets those has more than 10 corners
  */
  void SplitAnyFacet( std::vector<int> ptConcavity );

  Facet* facet1_;
  Mesh& mesh_;
  std::vector<Point*> ptlist_;
};
}
}

#endif

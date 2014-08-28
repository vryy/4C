/*----------------------------------------------------------------------*/
/*!
\file SplitHexIntoTwoWedges.H

\brief Split one Hex element in a discretization into two Wedge elements
and rebuild the discretization

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef SPLITHEXINTOTWOWEDGES_H_
#define SPLITHEXINTOTWOWEDGES_H_

#include "propagateCrack.H"

namespace DRT
{
namespace CRACK
{

class SplitHexIntoTwoWedges
{
public:

  /*!
  \brief Constructor
   */
  SplitHexIntoTwoWedges( Teuchos::RCP<DRT::Discretization> dis, int material_id )
    :discret_(dis),
     material_id_( material_id )
  {
  }

  /*!
  \brief Destructor
   */
  virtual ~SplitHexIntoTwoWedges(){};

  /*!
  \brief Main routine performing all the related operations
   */
  void DoAllSplittingOperations( DRT::Element * ele, std::vector<int> tipnodes, std::vector<int> splitnodes, int neweleid1, int neweleid2 );

private:

  /*!
  \brief split the given Quad into two Triangles along the diagonal given by splitnodes
   */
  void SplitThisQuad( std::vector<int> & face,
                      std::vector<int> & tipnodes,
                      std::vector<int> & splitnodes,
                      std::vector<int> & tri1,
                      std::vector<int> & tri2 );

  /*!
  \brief Add a WEDGE element with given nodal data into the discretization
   */
  void AddThisWedge( int owner, int eleid, std::vector<int> elenodes );

  //! underlying discretization
  Teuchos::RCP<DRT::Discretization> discret_;

  //! material id of elements that are in structural discretization
  int material_id_;
};

} /* namespace CRACK */
} /* namespace DRT */
#endif /* SPLITHEXINTOTWOWEDGES_H_ */

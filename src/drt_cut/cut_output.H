/*!-----------------------------------------------------------------------------------------------*
\file cut_output.H

\brief Handles file writing of all cut related stuff

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef CUTOUTPUT_H_
#define CUTOUTPUT_H_

namespace GEO {
namespace CUT {

class Element;
class Side;
class Node;

namespace OUTPUT {


/*!
\brief Write output of background element geometry
 */
void GmshElementDump( std::ofstream & file, Element * ele );

/*!
\brief Write output of background element geometry
 */
void GmshElementDump( std::ofstream & file, const std::vector<Node*> & nodes, char elementtype );

/*!
\brief Write output of a side
 */
void GmshSideDump( std::ofstream & file, const Side* s );

/*!
\brief Write output of the background element and all the cut sides corresponding to this element
 */
void GmshCompleteCutElement( std::ofstream & file, Element * ele );

}

} /* namespace CUT */
} /* namespace GEO */
#endif /* CUTOUTPUT_H_ */

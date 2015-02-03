/*!----------------------------------------------------------------------
\file xfem_utils.H
\brief Basic tools used in XFEM routines

<pre>
Maintainer:  Raffaela Kruse
             kruse@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
*----------------------------------------------------------------------*/

#ifndef XFEM_UTILS_H_
#define XFEM_UTILS_H_

#include "../drt_lib/drt_discret.H"

namespace XFEM
{
  namespace UTILS
  {
    void PrintDiscretizationToStream(
      Teuchos::RCP<DRT::Discretization> dis,
      const std::string& disname,
      bool elements,
      bool elecol,
      bool nodes,
      bool nodecol,
      bool faces,
      bool facecol,
      std::ostream& s,
      std::map<int, LINALG::Matrix<3,1> >* curr_pos = NULL
    );

    //! extract the nodal vectors and store them in node-vector-map
    //! \author schott \date 01/13
    void ExtractNodeVectors(
      Teuchos::RCP<DRT::Discretization>    dis,
      std::map<int, LINALG::Matrix<3,1> >& nodevecmap);
  }
}


#endif /* XFEM_UTILS_H_ */

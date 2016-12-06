/*----------------------------------------------------------------------*/
/*!
\file xfem_discretization_utils.H

\brief Basic discretization-related tools used in XFEM routines

\level 1

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef XFEM_DISCRETIZATION_UTILS_H_
#define XFEM_DISCRETIZATION_UTILS_H_

#include <vector>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../linalg/linalg_fixedsizematrix.H"

namespace DRT
{
  class Discretization;
  class Condition;
}

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

    class XFEMDiscretizationBuilder
    {

    public:

      void SetupXFEMDiscretization(
        const Teuchos::ParameterList&      xgen_params,
        Teuchos::RCP<DRT::Discretization>  dis,
        int numdof = 4);

      //! setup xfem discretization and embedded discretization
      void SetupXFEMDiscretization(
        const Teuchos::ParameterList&      xgen_params,
        Teuchos::RCP<DRT::Discretization>  dis,
        Teuchos::RCP<DRT::Discretization>  embedded_dis,
        int numdof = 4);

    private:

      //! split a discretization into two by removing conditioned nodes
      //! in source and adding to target
      void SplitDiscretizationByCondition(
        Teuchos::RCP<DRT::Discretization>  sourcedis,           //< initially contains all
        Teuchos::RCP<DRT::Discretization>  targetdis,           //< initially empty
        std::vector<DRT::Condition*>&      conditions,          //< conditioned nodes to be shifted to target
        const std::string&                 element_name,        //< element name (defining phys. field)
        const std::vector<std::string>&    conditions_to_copy   //< conditions to copy to target
      );

      //! re-partitioning of newly created discretizations (e.g. split by condition)
      void Redistribute(
        Teuchos::RCP<DRT::Discretization> dis,
        std::vector<int>& noderowvec,
        std::vector<int>& nodecolvec
      );

    };
  }
}


#endif /* XFEM_DISCRETIZATION_UTILS_H_ */

/*!----------------------------------------------------------------------
\file patspec.H

\brief A collection of methods to modify patient specific geometries

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

#ifndef PATSPEC_H
#define PATSPEC_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

class Epetra_Vector;

namespace DRT
{
  class Discretization;
}

namespace Teuchos
{
  class ParameterList;
}

namespace LINALG
{
  class SparseOperator;
}

namespace IO
{
  class DiscretizationWriter;
}

namespace PATSPEC
{

  void PatientSpecificGeometry(Teuchos::RCP<DRT::Discretization> dis,
	                       Teuchos::RCP<Teuchos::ParameterList> params);

  void ComputeEleNormalizedLumenDistance(Teuchos::RCP<DRT::Discretization> dis,
	                                 Teuchos::RCP<Teuchos::ParameterList> params);

  void ComputeEleLocalRadius(Teuchos::RCP<DRT::Discretization> dis);

  std::vector<double> GetCenterline(std::string filename);

  void ComputeEleStrength(Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<Teuchos::ParameterList> params);

  void GetILTDistance(const int eleid, 
                      Teuchos::ParameterList& params, 
                      DRT::Discretization& dis);

  void GetLocalRadius(const int eleid, 
                      Teuchos::ParameterList& params, 
                      DRT::Discretization& dis);

  void PatspecOutput(Teuchos::RCP<IO::DiscretizationWriter> output_,
		     Teuchos::RCP<DRT::Discretization> discret_,
		     Teuchos::RCP<Teuchos::ParameterList> params);
                      
  void CheckEmbeddingTissue(Teuchos::RCP<DRT::Discretization> discret,
                            Teuchos::RCP<LINALG::SparseOperator> stiff,
                            Teuchos::RCP<Epetra_Vector> fint,
                            Teuchos::RCP<Epetra_Vector> disn);




}
#endif /* PATSPEC_H */

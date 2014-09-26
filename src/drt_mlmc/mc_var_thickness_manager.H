/*----------------------------------------------------------------------*/
/*!
\file mc_var_thickness_manager.H
\brief class to modify wall thickness of cardiovascular structures on the fly

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            089 - 28915276
</pre>

!*/


#ifdef HAVE_FFTW
#ifndef MC_VAR_THICKNESS_MANAGER_H_
#define MC_VAR_THICKNESS_MANAGER_H_

/*----------------------------------------------------------------------*/
/* headers */

#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_dserror.H"
#include "../drt_ale/ale.H"
// forward declarations
namespace DRT
{
  class Discretization;
}

namespace LINALG
{
  class MultiMapExtractor;
}


namespace UQ
{

class RandomField;
/*----------------------------------------------------------------------*/
class MCVarThicknessManager
{
public:

  //! Constructor
  MCVarThicknessManager(Teuchos::RCP<DRT::Discretization> discret,const int my_thickness_field_id);

  //! Destructor
  ~MCVarThicknessManager(){};

  //brief create new realization of all random fields and push parameters down to material
  void SetUpThickness(unsigned int myseed, double para_cont_parameter, bool reuse_rf);

  //brief compute new set of random Phase angles/ new realizations of all random fields
  void CreateNewRealizationOfRandomField(unsigned int myseed);

  //brief write all random variables that are used for the simulation of stochastic processes to file
  void WriteRandomVariablesToFile(std::string filename, int numrun);

  //! brief reset geometry
  void ResetGeometry();

  //! brief evaluate underlying randomfield at a specific location
  double EvalThicknessAtLocation(std::vector<double> myloc, double para_cont_parameter);

  double InitialThickness(){return initial_wall_thickness_;};


private:

  //! brief set ALE displacement conditions to move the uncertain surface
  void ModifyConditions( const std::map<int, std::vector<double> >& ale_bc_nodes);

  //! brief add the given node ids and Displacement conditions to discretization
  void AddConditions( Teuchos::RCP<DRT::Discretization> discret,
                                           const std::map<int, std::vector<double> >& ale_bc_nodes );

  //! brief Delete all Dirichlet conditions that has only one associated node
  void DeleteConditions( Teuchos::RCP<DRT::Discretization> discret );

  //! brief Build ALE system matrix and solve the system
  void ALESolve();

  //! brief perform all necessary ALE steps
  void ALEStep( const std::map<int, std::vector<double> >& ale_bc_nodes);

  //! change position of uncertain surface nodes and perform ALE step to
  // adjust position of all nodes in the wall
  void SetThickness(double para_cont_parameter);

  //! brief create new random field based on input data
  Teuchos::RCP<RandomField> CreateRandomField(int random_field_id, unsigned int myseed);

  //! brief compute surface normal for all nodes in uncertain surface condition
  //! to accomplish this task this function creates a mortar interface using
  //  a special setup routine from the mortar adapter and then extracts
  // the surface normals computed by the mortar interface.
  // To create a mortar interface just to compute surface normals is probably a bit of an
  // overkill, however it works well with little implementational effort. Furthermore,
  // we only need to do this once in a UQ simulation and thus efficiency is not an issue here
  void ComputeNormals();

  //! structure discretization for convenience
  Teuchos::RCP<DRT::Discretization> discret_;

  //! underlying ALE
  Teuchos::RCP<ALE::Ale>   ale_;

  //! random field to model wall thickness
  Teuchos::RCP<RandomField>  randomfield_;

  const std::vector<int>* my_uncert_nodeids_;

  //! total number of nodes in uncertain surface condition
  int my_uncert_nnodes_;

  //! original node positions of all nodes in uncertain surface condition (redundant on all procs)
  std::map<int,std::vector<double > >  my_uncertain_nodes_org_pos_;

  //! delta in node positions of all nodes in uncertain surface condition based on random field (redundant on all procs)
  std::map<int,std::vector<double > >  my_uncertain_nodes_pos_;

  //! averaged nodal normals of all nodes in uncertain surface condition (redundant on all procs)
  std::map<int,std::vector<double> >  my_uncertain_nodes_normals_;

  // store initial geometry of the whole discretization, needed for reset
  Teuchos::RCP<Epetra_Vector> org_geom_;

  // brief initial wall thickness (asssumed uniform on geometry)
  double initial_wall_thickness_;

};
} // namespace MLMC



#endif /* MC_VAR_THICKNESS_MANAGER_H_ */
#endif // FFTW

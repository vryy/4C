/*----------------------------------------------------------------------*/
/*!
\file mc_var_thickness_manager.H
\brief class to modify wall thickness of cardiovascular structures on the fly

\maintainer Jonas Biehler
\level 3
*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#ifdef HAVE_FFTW
#ifndef MC_VAR_THICKNESS_MANAGER_H_
#define MC_VAR_THICKNESS_MANAGER_H_

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace DRT
{
  class Discretization;
}

namespace LINALG
{
  class MultiMapExtractor;
}

namespace ADAPTER
{
  class Ale;
  class Coupling;
}
namespace UQ
{
class RandomField;

/*----------------------------------------------------------------------------*/
/// Class to handle uncertain wall thickness in an uq routine
class MCVarThicknessManager
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor
  */
  MCVarThicknessManager(
      Teuchos::RCP<DRT::Discretization> discret, ///< discretization
      const int my_thickness_field_id);          ///< id of random field

  //! Destructor
  ~MCVarThicknessManager(){};

  //@}

  //!brief create new realization of all random fields and push parameters down to material
  void SetUpThickness(
      unsigned int myseed,
      double para_cont_parameter,
      bool reuse_rf
      );

  //brief compute new set of random Phase angles/ new realizations of all random fields
  void CreateNewRealizationOfRandomField(
      unsigned int myseed         ///< random seed for generation of random phase angles
      );

  //! \brief reset discretization to its original state
  void ResetGeometry();


  /*! \brief Evaluate random field for variable wall thickness
  *   Thickness will be evaluated at closest point to element center
  */
  double EvalThicknessAtLocation(
      int ele_GID,  ///< spatial location is ecoded by element GID
      double para_cont_parameter  ///< parameter continuation parameter for blending over thickness
      );

  //! brief return constant part of global, initial wall thickness
  double InitialThickness(){return initial_wall_thickness_;};

private:

  /*! brief add point dirichlet conditions to ale discretization
  *
  * Add the given node ids and displacement conditions to ale discretization
  * Type of the condition added is PointDirichlet
  * For each node, a separate condition is added. This is helpful to identify
  * these conditions later and delete them at appropriate time
  */
  void AddPointDBCToAleDiscretization(
      const std::map<int, std::vector<double> >& ale_bc_nodes_val,  ///< dbc values
      const std::map<int, std::vector<int> >& ale_bc_nodes_onoff,   ///< dbc on off toggle
      const std::map<int, std::vector<int> >& ale_bc_nodes_funct    ///< dbc functions
      );


  /*! brief Delete all point dirichlet conditions from ale discretization
   *
   * Delete all DBC point condition that have only one associated node.
   * These are the conditions that we added to move the uncertain surface around
  */
  void DeletePointDBCConditionsFromAleDiscretization();

  //! brief create new random field based on input data
  Teuchos::RCP<RandomField> CreateRandomField(
      int random_field_id,
      unsigned int myseed
      );

  void ComputeNewAleDBCFromRF();

  /*! brief Change position of uncertain surface nodes using ALE
   *  Change the geometry of the structure discretization based on random field.
   *  Therefor, an ALE discretization is used to which point DBCs are applied.
   *  After solving the ALE field, the resulting displacement field is used to
   *  update the reference configuarion of the structure
   */
  void ModifyGeometryBasedOnRF(const unsigned int myseed);


  /*! \brief Compute surface normal for all nodes in uncertain surface condition.
  *
  * To accomplish this task this function creates a mortar interface using
  * a special setup routine from the mortar adapter and then extracts
  * the surface normals computed by the mortar interface.
  * To create a mortar interface just to compute surface normals is probably a bit of an
  * overkill, however it works well with little implementational effort. Furthermore,
  * we only need to do this once in a UQ simulation and thus efficiency is not an issue here
  */
  void ComputeNormals();


  /*! \brief sort uncertain nodes into bins to speed up closest point search
  *
  */
  void SortUncertainNodesIntoBins();


  void ComputeMapWallElementsClosestSurfacePoint();


  //! brief structure discretization for convenience
  Teuchos::RCP<DRT::Discretization> discret_;

  //! bief ale discretization for convenience
  Teuchos::RCP<DRT::Discretization> aledis_;

  //! brief random field to model wall thickness
  Teuchos::RCP<RandomField> randomfield_;

  //! brief vector with all nodes subject to uncertainty
  std::vector<int> my_uncert_nodeids_;

  //! brief map with all nodes subject to uncertainty sorted into bins with unique bin id
  std::map<int, std::vector<int> > my_uncert_nodeids_in_bins_;

  //! brief map with giving closest uncertain nodes for each wall element
  // ele_GID , node_ID
  std::map<int, int > closest_uncertain_node_;

  //! brief total number of nodes in uncertain surface condition
  int my_uncert_nnodes_;

  //! brief original node positions of all nodes in uncertain surface condition (redundant on all procs)
  std::map<int, std::vector<double> > my_uncertain_nodes_org_pos_;

  //! brief onoff toggle vector for point dbc
  std::map<int, std::vector<int> > my_uncertain_nodes_dbc_onoff_;

  //! brief value vector for point dbc which contains local mean value of RF
  std::map<int, std::vector<double> > my_uncertain_nodes_dbc_val_;

  //! brief funct vector for point dbc
  std::map<int, std::vector<int> > my_uncertain_nodes_dbc_funct_;

  //! brief  delta in node positions of all nodes in uncertain surface condition based on random field
  std::map<int, std::vector<double> > my_uncertain_nodes_delta_pos_;

  //! brief store initial geometry of the whole discretization, needed for reset
  Teuchos::RCP<Epetra_Vector> org_geom_;

  //! brief initial wall thickness (asssumed uniform on geometry)
  double initial_wall_thickness_;

  //! brief use median local wall thickness based on condition input
  bool use_median_field_from_cond_;

  //! brief start random geometry above bifurcation of AAA only
  bool start_rand_geo_above_bif_;

  //!brief location of bifurcation of AAA plus an offset
  double z_pos_start_rf_;

  //!brief width of transition domain in which random thickness is ramped up
  double transition_width_;

  //! coupling of structure and ale (whole field)
  Teuchos::RCP<ADAPTER::Coupling> coupsa_;

  /*! brief averaged nodal normals of all nodes in uncertain surface condition
  *
  *  These nodal normals are currently computed using the structure
  *  discretization. It is assumed that node IDs of structure and
  *  ALE discretization are the same
  */
  std::map<int, std::vector<double> > my_uncertain_nodes_normals_;

  };
} // namespace MLMC



#endif /* MC_VAR_THICKNESS_MANAGER_H_ */
#endif // FFTW

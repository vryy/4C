/*---------------------------------------------------------------------*/
/*! \file
\brief Abstract data container class for contact solution strategies

\level 2


*/
/*---------------------------------------------------------------------*/


#include "4C_contact_abstract_data_container.hpp"

FOUR_C_NAMESPACE_OPEN

CONTACT::AbstractStratDataContainer::AbstractStratDataContainer()
    : glmdofrowmap_(Teuchos::null),
      gsnoderowmap_(Teuchos::null),
      gmnoderowmap_(Teuchos::null),
      gsdofrowmap_(Teuchos::null),
      gmdofrowmap_(Teuchos::null),
      gndofrowmap_(Teuchos::null),
      gsmdofrowmap_(Teuchos::null),
      gdisprowmap_(Teuchos::null),
      gactivenodes_(Teuchos::null),
      gactivedofs_(Teuchos::null),
      ginactivenodes_(Teuchos::null),
      ginactivedofs_(Teuchos::null),
      gactiven_(Teuchos::null),
      gactivet_(Teuchos::null),
      gslipnodes_(Teuchos::null),
      gslipdofs_(Teuchos::null),
      gslipt_(Teuchos::null),
      gsdof_vertex_(Teuchos::null),
      gsdof_edge_(Teuchos::null),
      gsdof_surf_(Teuchos::null),
      unbalance_evaluation_time_(0),
      unbalance_num_slave_elements_(0),
      pglmdofrowmap_(Teuchos::null),
      pgsdofrowmap_(Teuchos::null),
      pgmdofrowmap_(Teuchos::null),
      pgsmdofrowmap_(Teuchos::null),
      pgsdirichtoggle_(Teuchos::null),
      partype_(Inpar::Mortar::ParallelRedist::redist_none),
      initial_elecolmap_(Teuchos::null),
      dmatrix_(Teuchos::null),
      mmatrix_(Teuchos::null),
      wgap_(Teuchos::null),
      tangrhs_(Teuchos::null),
      inactiverhs_(Teuchos::null),
      str_contact_rhs_ptr_(Teuchos::null),
      constrrhs_(Teuchos::null),
      lindmatrix_(Teuchos::null),
      linmmatrix_(Teuchos::null),
      kteffnew_(Teuchos::null),
      dold_(Teuchos::null),
      mold_(Teuchos::null),
      z_(Teuchos::null),
      zold_(Teuchos::null),
      zincr_(Teuchos::null),
      zuzawa_(Teuchos::null),
      stressnormal_(Teuchos::null),
      stresstangential_(Teuchos::null),
      forcenormal_(Teuchos::null),
      forcetangential_(Teuchos::null),
      stepnp_(-1),
      iter_(-1),
      isincontact_(false),
      wasincontact_(false),
      wasincontactlts_(false),
      isselfcontact_(false),
      friction_(false),
      non_smooth_contact_(false),
      regularized_(false),
      dualquadslavetrafo_(false),
      trafo_(Teuchos::null),
      invtrafo_(Teuchos::null),
      dmatrixmod_(Teuchos::null),
      doldmod_(Teuchos::null),
      inttime_(0.0),
      ivel_(0),
      stype_(Inpar::CONTACT::solution_vague),
      constr_direction_(Inpar::CONTACT::constr_vague)
{
  return;
}

FOUR_C_NAMESPACE_CLOSE

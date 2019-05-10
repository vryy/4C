/*!----------------------------------------------------------------------

\brief manage access to and provide data globally in immersed problems

\level 3

\maintainer Jonas Eichinger

*----------------------------------------------------------------------*/
#include "immersed_field_exchange_manager.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

///----------------------------------------------------------------------*/
/// the instance
///----------------------------------------------------------------------*/
DRT::ImmersedFieldExchangeManager* DRT::ImmersedFieldExchangeManager::instance_;

//----------------------------------------------------------------------*/
//    definition of the instance
//----------------------------------------------------------------------*/
DRT::ImmersedFieldExchangeManager* DRT::ImmersedFieldExchangeManager::Instance(bool create)
{
  if (create)
  {
    if (instance_ == NULL)
    {
      instance_ = new ImmersedFieldExchangeManager();
    }
  }
  else
  {
    if (instance_ != NULL) delete instance_;
    instance_ = NULL;
  }
  return instance_;
}

//----------------------------------------------------------------------*/
//    destruction method
//----------------------------------------------------------------------*/
void DRT::ImmersedFieldExchangeManager::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}

//----------------------------------------------------------------------*/
//    constructor
//----------------------------------------------------------------------*/
DRT::ImmersedFieldExchangeManager::ImmersedFieldExchangeManager()
{
  // initialize members
  num_isimmersedbry_ = 0;
  gap_max_ = 0.0;
  gap_min_ = 0.0;
  void_min_ = 0.0;
  void_max_ = 0.0;
  delta_porosity_max_ = 0.0;
  points_to_ecm_penalty_traction_ = Teuchos::null;
  isinitialized_ = false;
  currpositions_immerseddis_ = NULL;
  simple_ecm_interaction_constant = 0.0;
  numnlniter_ = 0;

  points_to_ecm_adhesion_force_ = Teuchos::null;
  points_to_curr_subset_of_backgrounddis_ = NULL;
  points_to_cell_adhesion_force_ = Teuchos::null;
  points_to_cell_adhesion_fixpoints_ = Teuchos::null;

  isfluidinteraction_ = false;
  isprotrusionformation = false;
  isPureAdhesionSimulation_ = false;
  isPureConfinementSimulation_ = false;
  isPureProtrusionFormationSimulation_ = false;

  points_to_phin_ = Teuchos::null;
  points_to_phinp_ = Teuchos::null;
  points_to_phidtn_ = Teuchos::null;
  points_to_phidtnp_ = Teuchos::null;
  points_to_Rates_ = Teuchos::null;
  points_to_Rates_Actin_ = Teuchos::null;

  points_to_auxdis_ = Teuchos::null;

  not_assemble_to_pseudoboundary_switch = false;
}


//----------------------------------------------------------------------*/
//----------------------------------------------------------------------*/
void DRT::ImmersedFieldExchangeManager::InitializePorosityAtGPMap()
{
  DRT::Condition* mycondition =
      DRT::Problem::Instance()->GetDis("cell")->GetCondition("IMMERSEDCoupling");
  if (mycondition == NULL) dserror("Failed to get condition IMMERSEDCoupling");
  std::map<int, Teuchos::RCP<DRT::Element>> mygeometry = mycondition->Geometry();
  int mygeometrysize = mygeometry.size();
  if (mygeometrysize == 0) dserror("My geometry is empty");

  //  // DEBUG output
  //  for (std::map<int,Teuchos::RCP<DRT::Element> >::iterator it=mygeometry.begin();
  //  it!=mygeometry.end(); ++it)
  //     std::cout <<"PROC "<<DRT::Problem::Instance()->GetDis("cell")->Comm().MyPID()<<": "<<
  //     it->first << " => " << it->second->Id() << '\n';

  for (std::map<int, Teuchos::RCP<DRT::Element>>::iterator it = mygeometry.begin();
       it != mygeometry.end(); ++it)
  {
    // vector o the 4 gp values of the porosity for element with key id
    std::vector<double> vectortoinsert(4);
    std::vector<LINALG::Matrix<3, 1>> vectorofmatricestoinsert(4);
    for (int gp = 0; gp < 4; ++gp)
    {
      LINALG::Matrix<3, 1> matrixtoinsert(true);
      vectorofmatricestoinsert[gp] = matrixtoinsert;
    }
    vectortoinsert[0] = 0.8;
    vectortoinsert[1] = 0.8;
    vectortoinsert[2] = 0.8;
    vectortoinsert[3] = 0.8;

    int id = it->second->Id();

    porosity_at_gp_.insert(std::pair<int, std::vector<double>>(id, vectortoinsert));
    porosity_at_gp_old_timestep_.insert(std::pair<int, std::vector<double>>(id, vectortoinsert));
    penalty_traction_at_gp_.insert(
        std::pair<int, std::vector<LINALG::Matrix<3, 1>>>(id, vectorofmatricestoinsert));
  }

  return;
}

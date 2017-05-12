/*!----------------------------------------------------------------------
\file immersed_field_exchange_manager.H

\brief manage access to and provide data globally in immersed problems

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

\level 3

*----------------------------------------------------------------------*/

#ifndef IMMERSED_FIELD_EXCHANGE_MANAGER_H
#define IMMERSED_FIELD_EXCHANGE_MANAGER_H

#include <map>
#include <cstdlib>
#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"

#include "../drt_lib/drt_singletondestruction.H"
#include "../linalg/linalg_fixedsizematrix.H"

namespace GEO
{
class SearchTree;
}

namespace ADAPTER
{
class FSIStructureWrapperImmersed;
}

namespace DRT
{

class ImmersedFieldExchangeManager : public DRT::SingletonDestruction
{

public:

  /// virtual empty destructor (dtor)
  virtual ~ImmersedFieldExchangeManager() {}

  /// singleton access method
  static ImmersedFieldExchangeManager* Instance( bool create=true );

  /// called upon destruction
  virtual void Done();

  //! @name Initialization Flag
  //@{
  bool isinitialized_; //!< true if this data manager has been properly initialized
  //@}

  //! @name selective assembly specific (e.g. for treatment of artificial ECM)
  //@{
  std::vector<int> eles_;  //!< vector of elements
  std::vector<int> beles_; //!< vector of elements cut by the boundary
  bool not_assemble_to_pseudoboundary_switch; //!< global switch to turn selective assembly on/off
  //@}

  //! @name Simulation Type Flags
  //@{
  bool isfluidinteraction_;
  bool isprotrusionformation;
  bool isPureAdhesionSimulation_;
  bool isPureConfinementSimulation_;
  bool isPureProtrusionFormationSimulation_;
  //@}

  //! @name Misc. Counters
  //@{
  int num_isimmersedbry_;
  int numnlniter_;
  //@}

  //! @name Confinement Simulation Variables
  //@{
  double gap_max_;
  double gap_min_;
  double void_min_;
  double void_max_;
  double delta_porosity_max_;
  double simple_ecm_interaction_constant;
  //@}

  //! @name Misc. std maps
  //@{
  std::map<int,std::vector<double> > porosity_at_gp_;
  std::map<int,std::vector<double> > porosity_at_gp_old_timestep_;
  std::map<int,std::vector<double> > porositiy_old_;
  std::map<int,LINALG::Matrix<3,1> >* currpositions_immerseddis_;
  std::map<int,std::vector<LINALG::Matrix<3,1> > > penalty_traction_at_gp_;
  std::map<int,std::set<int> >* points_to_curr_subset_of_backgrounddis_;
  //@}

  //! @name Misc. epetra vectors
  //@{
  Teuchos::RCP<Epetra_Vector> points_to_ecm_penalty_traction_;
  Teuchos::RCP<Epetra_Vector> points_to_ecm_penalty_gap_;
  Teuchos::RCP<Epetra_Vector> cell_current_nodal_normals_;
  Teuchos::RCP<Epetra_Vector> points_to_ecm_adhesion_force_;
  Teuchos::RCP<Epetra_Vector> points_to_surface_traction_;
  //@}

  //! @name Misc. epetra multivectors
  //@{
  Teuchos::RCP<Epetra_MultiVector> points_to_phin_;
  Teuchos::RCP<Epetra_MultiVector> points_to_phinp_;
  Teuchos::RCP<Epetra_MultiVector> points_to_phidtn_;
  Teuchos::RCP<Epetra_MultiVector> points_to_phidtnp_;
  Teuchos::RCP<Epetra_MultiVector> points_to_Rates_;
  Teuchos::RCP<Epetra_MultiVector> points_to_Rates_Actin_;
  //@}

  //! @name Misc. linalg matrices
  //@{
  LINALG::Matrix<3,1> max_gap_searchdirection_;
  LINALG::Matrix<3,1> max_gap_spacepoint_;
  LINALG::Matrix<3,1> min_gap_searchdirection_;
  LINALG::Matrix<3,1> min_gap_spacepoint_;
  //@}

  /// pointer to adapter FSIStructureWrapperImmersed
  Teuchos::RCP<ADAPTER::FSIStructureWrapperImmersed> FSIStructureWrapperImmersed_;

  /// pointer to immersed search tree
  Teuchos::RCP<GEO::SearchTree> immersed_SearchTree_;

  /*!
    \brief initialize the map containing the porosity at integration points (gauss points)

    \note only needed during confinement simulation

   */
  void InitializePorosityAtGPMap();

  /*!
    \brief update the map containing the porosity at integration points (gauss points)

    \note only needed during confinement simulation

   */
  void UpdatePorosityAtGPOldTimestep(){ porosity_at_gp_old_timestep_ = porosity_at_gp_; return;};

  //! @name Set and get routines for all members
  //@{
  void SetIsInitialized(bool trueorfalse){isinitialized_=trueorfalse; return;};

  bool IsInitialized(){return isinitialized_;};

  void SetEles(std::vector<int> eles){eles_=eles; return;};

  std::vector<int>* GetEles(){return &eles_;};

  void SetBdryEles(std::vector<int> eles){beles_=eles; return;};

  std::vector<int>* GetBdryEles(){return &beles_;};

  void SetPseudoBoundarySwitch(bool trueorfalse){not_assemble_to_pseudoboundary_switch=trueorfalse;};
  bool GetPseudoBoundarySwitch(){return not_assemble_to_pseudoboundary_switch;};

  void SetIsFluidInteraction(bool trueorfalse){isfluidinteraction_=trueorfalse; return;};

  bool IsFluidInteraction(){return isfluidinteraction_;};

  void SetIsProtrusionFormation(bool trueorfalse){isprotrusionformation=trueorfalse; return;};

  bool IsProtrusionFormation(){return isprotrusionformation;};

  void SetImmersedSearchTree(Teuchos::RCP<GEO::SearchTree> immersed_SearchTree)
  {
    immersed_SearchTree_ = immersed_SearchTree;
    return;
  };

  Teuchos::RCP<GEO::SearchTree>  GetImmersedSearchTree()
  {
    return immersed_SearchTree_;
  };

  void SetCurrentPositionsImmersedDis(std::map<int,LINALG::Matrix<3,1> >* currpositions_immerseddis)
  {
    currpositions_immerseddis_ = currpositions_immerseddis;
    return;
  };

  std::map<int,LINALG::Matrix<3,1> >* GetCurrentPositionsImmersedDis()
  {
    return currpositions_immerseddis_;
  };

  void IncrementNumIsImmersedBoundary(){num_isimmersedbry_++; return;};

   int GetNumIsImmersedBoundary(){return num_isimmersedbry_;};

   void SetNumIsImmersedBoundary(int value){num_isimmersedbry_=value; return;};

   void SetOldPorosity(int eid, int gpnum, double porosity){porositiy_old_.at(eid); return;};

   void SetGapMax(double value){gap_max_=value; return;};

   void SetGapMin(double value){gap_min_=value; return;};

   double GetGapMax(){return gap_max_;};

   double GetGapMin(){return gap_min_;};

   void SetVoidMax(double value){void_max_=value; return;};

   void SetVoidMin(double value){void_min_=value; return;};

   double GetVoidMax(){return void_max_;};

   double GetVoidMin(){return void_min_;};

   void SetDeltaPorosityMax(double value){delta_porosity_max_=value; return;};

   double GetDeltaPorosityMax(){return delta_porosity_max_;};

   std::map<int,std::vector<double> >* GetPointerToPorosityAtGP(){return &porosity_at_gp_;};

   std::map<int,std::vector<double> >* GetPointerToPorosityAtGPOldTimestep(){return &porosity_at_gp_old_timestep_;};

   void SetPointerToECMPenaltyTraction(const Teuchos::RCP<Epetra_Vector> vector){points_to_ecm_penalty_traction_=vector; return;};

   Teuchos::RCP<Epetra_Vector> GetPointerToECMPenaltyTraction(){return points_to_ecm_penalty_traction_;};

   std::map<int,std::vector<LINALG::Matrix<3,1> > >* GetPointerToPenaltyTractionAtGPMap(){return &penalty_traction_at_gp_;};

   void SetPointerToGap(const Teuchos::RCP<Epetra_Vector> vector){points_to_ecm_penalty_gap_=vector; return;};

   Teuchos::RCP<Epetra_Vector> GetPointerToGap(){return points_to_ecm_penalty_gap_;};

   void SetPointerToCurrentNodalNormals(const Teuchos::RCP<Epetra_Vector> vector){cell_current_nodal_normals_=vector; return;};

   Teuchos::RCP<Epetra_Vector> GetPointerToCurrentNodalNormals(){return cell_current_nodal_normals_;};

   void SetMaxGapSearchDirection(const LINALG::Matrix<3,1> vector){max_gap_searchdirection_=vector; return;};

   void SetMinGapSearchDirection(const LINALG::Matrix<3,1> vector){min_gap_searchdirection_=vector; return;};

   void SetMaxGapSpacePoint(const LINALG::Matrix<3,1> vector){max_gap_spacepoint_=vector; return;};

   void SetMinGapSpacePoint(const LINALG::Matrix<3,1> vector){min_gap_spacepoint_=vector; return;};

   const LINALG::Matrix<3,1>* GetMaxGapSearchDirection(){return &max_gap_searchdirection_;};

   const LINALG::Matrix<3,1>* GetMinGapSearchDirection(){return &min_gap_searchdirection_;};

   const LINALG::Matrix<3,1>* GetMaxGapSpacePoint(){return &max_gap_spacepoint_;};

   const LINALG::Matrix<3,1>* GetMinGapSpacePoint(){return &min_gap_spacepoint_;};

   void SetSimpleECMInteractionConstant(double value){simple_ecm_interaction_constant = value; return;};

   double GetSimpleECMInteractionConstant(){return simple_ecm_interaction_constant;};

   void SetAdapter(Teuchos::RCP<ADAPTER::FSIStructureWrapperImmersed> adapter){FSIStructureWrapperImmersed_ = adapter; return;};

   Teuchos::RCP<ADAPTER::FSIStructureWrapperImmersed> GetAdapter(){return FSIStructureWrapperImmersed_;};

   void SetPointerECMAdhesionForce(const Teuchos::RCP<Epetra_Vector> vector){points_to_ecm_adhesion_force_=vector; return;};

   Teuchos::RCP<Epetra_Vector> GetPointerECMAdhesionForce(){return points_to_ecm_adhesion_force_;};

   void SetPointerSurfaceTraction(Teuchos::RCP<Epetra_Vector> vector){points_to_surface_traction_=vector; return;};

   Teuchos::RCP<Epetra_Vector>& GetPointerSurfaceTraction(){return points_to_surface_traction_;};

   void SetPointerToCurrentSubsetOfBackgrdDis(std::map<int,std::set<int> >* pointer){points_to_curr_subset_of_backgrounddis_ = pointer;};

   std::map<int,std::set<int> >* GetPointerToCurrentSubsetOfBackgrdDis(){return points_to_curr_subset_of_backgrounddis_;};

   void SetIsPureAdhesionSimulation(bool trueorfalse){isPureAdhesionSimulation_=trueorfalse; return;};

   bool IsPureAdhesionSimulation(){return isPureAdhesionSimulation_;};

   void SetIsPureConfinementSimulation(bool trueorfalse){isPureConfinementSimulation_=trueorfalse; return;};

   bool IsPureConfinementSimulation(){return isPureConfinementSimulation_;};

   void SetIsPureProtrusionFormationSimulation(bool trueorfalse){isPureProtrusionFormationSimulation_=trueorfalse; return;};

   bool IsPureProtrusionFormationSimulation(){return isPureProtrusionFormationSimulation_;};

   void SetPointerToPhins(const Teuchos::RCP<Epetra_MultiVector> vector){points_to_phin_=vector; return;};

   Teuchos::RCP<Epetra_MultiVector> GetPointerToPhins(){return points_to_phin_;};

   void SetPointerToPhinps(const Teuchos::RCP<Epetra_MultiVector> vector){points_to_phinp_=vector; return;};

   Teuchos::RCP<Epetra_MultiVector> GetPointerToPhinps(){return points_to_phinp_;};

   void SetPointerToPhidtns(const Teuchos::RCP<Epetra_MultiVector> vector){points_to_phidtn_=vector; return;};

   Teuchos::RCP<Epetra_MultiVector> GetPointerToPhidtns(){return points_to_phidtn_;};

   void SetPointerToPhidtnps(const Teuchos::RCP<Epetra_MultiVector> vector){points_to_phidtnp_=vector; return;};

   Teuchos::RCP<Epetra_MultiVector> GetPointerToPhidtnps(){return points_to_phidtnp_;};

   void SetPointerToRates(const Teuchos::RCP<Epetra_MultiVector> vector){points_to_Rates_=vector; return;};

   Teuchos::RCP<Epetra_MultiVector> GetPointerToRates(){return points_to_Rates_;};

   void SetPointerToRatesActin(const Teuchos::RCP<Epetra_MultiVector> vector){points_to_Rates_Actin_=vector; return;};

   Teuchos::RCP<Epetra_MultiVector> GetPointerToRatesActin(){return points_to_Rates_Actin_;};
   //@}

protected:

  /// protected constructor since we are a Singleton.
  ImmersedFieldExchangeManager();

  /// declaration of the static instance
  static DRT::ImmersedFieldExchangeManager* instance_;


}; // class ImmersedFieldExchangeManager

} // namespace DRT


#endif /*! IMMERSED_FIELD_EXCHANGE_MANAGER_H */

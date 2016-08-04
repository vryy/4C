/*!-----------------------------------------------------------------------------------------------*
\file xfem_coupling_base.H

\brief is the base for the different types of mesh and level-set based coupling conditions and thereby builds the bridge between the
xfluid class and the cut-library

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef XFEM_COUPLING_BASE_H_
#define XFEM_COUPLING_BASE_H_


#include <vector>
#include <Teuchos_RCP.hpp>
//
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

#include <Epetra_Vector.h>

namespace XFEM
{

typedef std::pair<INPAR::XFEM::EleCouplingCondType, DRT::Condition*> EleCoupCond;

INPAR::XFEM::EleCouplingCondType CondType_stringToEnum(const std::string& condname);

class CouplingBase
{



public:

  //! which boolean set operator used to combine current field with previous one
  enum LevelSetBooleanType
  {
    ls_none = 0,           // used for first Boundary condition level-setcoupling
    ls_cut = 1,            // latex: \cap:         Omega 1 \cap \Omega 2
    ls_union = 2,          // latex: \cup          Omega 1 \cup \Omega 2
    ls_difference = 3,     // latex: \backslash    Omega 1 - Omega 2
    ls_sym_difference =4   // latex: \triangle     (Omega 1 - Omega 2) \cup (Omega 2 - \Omega 1)
  };

  //! constructor
  explicit CouplingBase(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,    ///< background discretization
      const std::string &                 cond_name, ///< name of the condition, by which the derived cutter discretization is identified
      Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< full discretization from which the cutter discretization is derived
      const int                           coupling_id,///< id of composite of coupling conditions
      const double                        time,      ///< time
      const int                           step       ///< time step
  );

  //! destructor
  virtual ~CouplingBase()
  {
    return;
  }

  void SetTimeAndStep(
      const double time,
      const int step
  )
  {
    time_ = time;
    step_ = step;
  }

  void IncrementTimeAndStep(
      const double dt
  )
  {
    dt_    = dt;
    time_ += dt;
    step_ += 1;
  }

  void GetConditionByCouplingId(
      const std::vector<DRT::Condition*> & mycond,
      const int coupling_id,
      std::vector<DRT::Condition*> & mynewcond
  );

  void Status(
      const int coupling_idx,
      const int side_start_gid
  );


  const std::string DisNameToString(
      Teuchos::RCP<DRT::Discretization> dis
  )
  {
    if(dis == Teuchos::null)
      return "---";

    return dis->Name();
  }

  const std::string TypeToStringForPrint(const INPAR::XFEM::EleCouplingCondType& type)
  {
    if     (type == INPAR::XFEM::CouplingCond_SURF_FSI_PART)             return "XFSI Partitioned";
    else if(type == INPAR::XFEM::CouplingCond_SURF_FSI_MONO)             return "XFSI Monolithic";
    else if(type == INPAR::XFEM::CouplingCond_SURF_CRACK_FSI_PART)       return "XFSI Partitioned with Crack";
    else if(type == INPAR::XFEM::CouplingCond_SURF_FLUIDFLUID)           return "FLUID-FLUID Coupling";
    else if(type == INPAR::XFEM::CouplingCond_LEVELSET_WEAK_DIRICHLET)   return "WEAK DIRICHLET BC / LS";
    else if(type == INPAR::XFEM::CouplingCond_LEVELSET_NEUMANN)          return "NEUMANN BC        / LS";
    else if(type == INPAR::XFEM::CouplingCond_LEVELSET_NAVIER_SLIP)      return "NAVIER SLIP BC    / LS";
    else if(type == INPAR::XFEM::CouplingCond_LEVELSET_TWOPHASE)         return "TWO-PHASE Coupling";
    else if(type == INPAR::XFEM::CouplingCond_LEVELSET_COMBUSTION)       return "COMBUSTION Coupling";
    else if(type == INPAR::XFEM::CouplingCond_SURF_WEAK_DIRICHLET)       return "WEAK DIRICHLET BC / MESH";
    else if(type == INPAR::XFEM::CouplingCond_SURF_NEUMANN)              return "NEUMANN BC        / MESH";
    else if(type == INPAR::XFEM::CouplingCond_SURF_NAVIER_SLIP)          return "NAVIER SLIP BC    / MESH";
    else dserror("unsupported coupling condition type %i", type);

    return "UNKNOWN";
  }

  const std::string AveragingToStringForPrint(const INPAR::XFEM::AveragingStrategy& strategy)
  {
    if     (strategy == INPAR::XFEM::Xfluid_Sided)      return "XFLUID-sided averaging";
    else if(strategy == INPAR::XFEM::Embedded_Sided)    return "EMBEDDED-sided averaging";
    else if(strategy == INPAR::XFEM::Mean)              return "MEAN averaging";
    else if(strategy == INPAR::XFEM::Harmonic)          return "HARMONIC averaging";
    else if(strategy == INPAR::XFEM::invalid)           return "INVALID";
    else dserror("unsupported averaging strategy %i", strategy);

    return "UNKNOWN";
  }


  const EleCoupCond & GetCouplingCondition(
      const int gid ///< global element element id w.r.t cutter discretization (bgele->Id for LevelsetCoupling cut and side-Id for MeshCoupling)
  )
  {
    int lid = cutter_dis_->ElementColMap()->LID(gid);
    return cutterele_conds_[lid];
  }

  const std::string & GetName(){return cond_name_;}

  Teuchos::RCP<DRT::Discretization> GetCutterDis(){return cutter_dis_;}
  Teuchos::RCP<DRT::Discretization> GetCouplingDis(){return coupl_dis_;}
  Teuchos::RCP<DRT::Discretization> GetCondDis(){return cond_dis_;}

  INPAR::XFEM::AveragingStrategy GetAveragingStrategy(){return averaging_strategy_;}

  virtual void PrepareSolve(){};

  virtual bool HasMovingInterface() = 0;

  virtual void EvaluateCouplingConditions(
      LINALG::Matrix<3,1>& ivel,
      LINALG::Matrix<3,1>& itraction,
      const LINALG::Matrix<3,1>& x,
      const DRT::Condition* cond
  ){dserror("EvaluateCouplingConditions should be implemented by derived class");};

  virtual void EvaluateCouplingConditionsOldState(
      LINALG::Matrix<3,1>& ivel,
      LINALG::Matrix<3,1>& itraction,
      const LINALG::Matrix<3,1>& x,
      const DRT::Condition* cond
  ){dserror("EvaluateCouplingConditionsOldState should be implemented by derived class");};

  /// set material pointer for coupling slave side
  virtual void GetInterfaceSlaveMaterial(
    DRT::Element* actele,
    Teuchos::RCP<MAT::Material> & mat
  )
  {
    mat = Teuchos::null;
  }

  /// get the sliplength for the specific coupling condition
  virtual void GetSlipCoefficient(
      double& slipcoeff,
      const LINALG::Matrix<3,1>& x,
      const DRT::Condition* cond)
  {
    slipcoeff = 0.0;
  }

  std::map<INPAR::XFEM::CoupTerm, std::pair<bool,double> >& GetConfigurationmap(
      double& kappa_m,    //< fluid sided weighting
      double& visc_m,     //< fluid sided weighting
      double& visc_s,     //< slave sided dynamic viscosity
      double& visc_stab,  //< viscous NIT Penalty scaling
      double& full_stab,  //< full NIT Penalty scaling
      const LINALG::Matrix<3,1>& x, //<Position x
      const DRT::Condition* cond //< Condition
      )
  {
    UpdateConfigurationMap_GP(kappa_m,visc_m,visc_s, visc_stab, full_stab,x,cond);
    #ifdef DEBUG
    //Do some safety checks, every combination which is not handled correct by the element level should be caught here ...

    // At the moment you cannot use different consistent constraints between Adjoint and Penalty terms
    //Check if we need a more general implementation (If constraints between Adjoint and Penalty are not the same!)
    if(fabs(configuration_map_.at(INPAR::XFEM::X_Adj_Col).second - configuration_map_.at(INPAR::XFEM::X_Pen_Col).second) > 1e-16 ||
       fabs(configuration_map_.at(INPAR::XFEM::X_Adj_n_Col).second - configuration_map_.at(INPAR::XFEM::X_Pen_n_Col).second) > 1e-16 ||
       fabs(configuration_map_.at(INPAR::XFEM::X_Adj_t_Col).second - configuration_map_.at(INPAR::XFEM::X_Pen_t_Col).second) > 1e-16)
      dserror("Your consistent constraint for Penalty and Adjoint term is not equal, go to element level and split up velint_diff_ for penalty and adjoint!");

    //At the moment you cannot use different consistent constraints for normal or full direction!
    if(fabs(configuration_map_.at(INPAR::XFEM::X_Adj_Col).second - configuration_map_.at(INPAR::XFEM::X_Adj_n_Col).second) > 1e-16 ||
        fabs(configuration_map_.at(INPAR::XFEM::F_Adj_Col).second - configuration_map_.at(INPAR::XFEM::F_Adj_n_Col).second) > 1e-16)
      dserror("Your consistent constraint for Adjoint in normal and full direction are different, this is an issure for the pressure adjoint term, fix velint_diff_ used for this term!; ");
    #endif
    return configuration_map_;
  }

protected:

  //! Initializes configurationmap
  virtual void InitConfigurationMap();

  //! Updates configurationmap for specific Gausspoint
  virtual void UpdateConfigurationMap_GP(
      double& kappa_m,    //< fluid sided weighting
      double& visc_m,     //< fluid sided weighting
      double& visc_s,     //< slave sided dynamic viscosity
      double& visc_stab,  //< viscous NIT Penalty scaling
      double& full_stab,  //< full NIT Penalty scaling
      const LINALG::Matrix<3,1>& x, //<Position x
      const DRT::Condition* cond //< Condition
      ){return;}

  void SetElementConditions();

  void SetAveragingStrategy();

  void SetCouplingDiscretization();

  void EvaluateDirichletFunction(
      LINALG::Matrix<3,1>& ivel,
      const LINALG::Matrix<3,1>& x,
      const DRT::Condition* cond,
      double time
  );

  void EvaluateNeumannFunction(
      LINALG::Matrix<3,1>& itraction,
      const LINALG::Matrix<3,1>& x,
      const DRT::Condition* cond,
      double time
  );

  void EvaluateFunction(
      std::vector<double>& final_values,
      const double* x,
      const DRT::Condition* cond,
      const double time
  );

  void EvaluateScalarFunction(
      double & final_value,
      const double* x,
      const double & val,
      const DRT::Condition* cond,
      const double time
  );

  ///< background discretization
  Teuchos::RCP<DRT::Discretization> bg_dis_;

  ///------------------------
  // CUTTER-DISCRETIZATION specific member
  ///------------------------

  ///< name of the condition, by which the derived cutter discretization is identified
  const std::string cond_name_;

  ///< discretization from which the cutter discretization is derived
  Teuchos::RCP<DRT::Discretization> cond_dis_;

  ///< id of composite of coupling conditions
  const int coupling_id_;

  ///< discretization w.r.t which the interface is described and w.r.t which the state vectors describing the interface position are defined (bgdis for LevelSetCoupling and boundary dis for MeshCoupling)
  Teuchos::RCP<DRT::Discretization> cutter_dis_;

  ///< pairs of condition type and pointer to DRT::Condition for all column elements of the cutter discretization (bgdis for LevelSetCoupling and boundary dis for MeshCoupling)
  std::vector<EleCoupCond> cutterele_conds_;

  std::vector<std::string> conditions_to_copy_;   ///< list of conditions that will be copied to the new discretization and used to set for each cutter element


  ///------------------------
  // Coupling-DISCRETIZATION specific member
  ///------------------------

  ///< discretization with which the background discretization is coupled (structural dis, fluid dis, poro dis, scatra dis, boundary dis),
  ///< Teuchos::null in case that no coupling terms but only boundary terms are evaluated
  Teuchos::RCP<DRT::Discretization> coupl_dis_;

  ///< averaging strategy, type of weighting
  INPAR::XFEM::AveragingStrategy averaging_strategy_;

  int myrank_;

  double dt_; ///< current time step size

  double time_;

  int step_;

  ///< map which configures element level (which terms are evaluated & scaled with which value)
  std::map<INPAR::XFEM::CoupTerm, std::pair<bool,double> > configuration_map_;

};

} // namespace XFEM

#endif /*XFEM_COUPLING_BASE_H_*/

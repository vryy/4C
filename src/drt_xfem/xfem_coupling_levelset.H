/*----------------------------------------------------------------------*/
/*!
\file xfem_coupling_levelset.H

\brief manages the different types of level-set based coupling conditions and thereby builds the bridge between the
xfluid class and the cut-library

\level 2

<pre>
\maintainer Benedikt Schott, Magnus Winter
            {schott,winter}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef XFEM_COUPLING_LEVELSET_H_
#define XFEM_COUPLING_LEVELSET_H_

#include "xfem_coupling_base.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_inpar/inpar_twophase.H"
#include "../drt_cut/cut_point.H"

namespace GEO
{
  class CutWizard;
}

namespace XFEM
{

/*!
\brief
 */
class LevelSetCoupling : public CouplingBase
{
public:

  //! constructor
  explicit LevelSetCoupling(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
      const std::string &                 cond_name, ///< name of the condition, by which the derived cutter discretization is identified
      Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< full discretization from which the cutter discretization is derived
      const int                           coupling_id,///< id of composite of coupling conditions
      const double                        time,      ///< time
      const int                           step       ///< time step
  );


  //! destructor
  virtual ~LevelSetCoupling()
  {
    return;
  }

  void SetCouplingDofsets();

  bool HaveMatchingNodes(
      const Teuchos::RCP<DRT::Discretization> & dis_A,
      const Teuchos::RCP<DRT::Discretization> & dis_B
  );

  void MapCutterToBgVector(
      const Teuchos::RCP<DRT::Discretization> & source_dis,
      const Teuchos::RCP<Epetra_Vector> &       source_vec_dofbased,
      const int                                 source_nds,
      const Teuchos::RCP<DRT::Discretization> & target_dis,
      const Teuchos::RCP<Epetra_Vector> &       target_vec_nodebased,
      const int                                 target_nds
  );

  //TODO: sort the functions...

  virtual void SetCutterDiscretization();

  virtual void SetConditionSpecificParameters(){};

  virtual void PrepareCutterOutput();

  virtual void DoConditionSpecificSetup();


  /// set levelset field by function, return if interface moved compared to last time step
  bool SetLevelSetField(const double time);

  /// initialize level set based state vectors
  virtual void InitStateVectors();

  virtual void InitStateVectors_Bg();

  virtual void InitStateVectors_Cutter();

  /// set level-boolean type
  virtual void SetLevelSetBooleanType();

  virtual bool ApplyComplementaryOperator();

  virtual void Output(
      const int step,
      const double time,
      const bool write_restart_data,
      const int lsc_idx = 0
  );

  virtual void GmshOutput(
      const std::string & filename_base,
      const int step,
      const int gmsh_step_diff,
      const bool gmsh_debug_out_screen
  );

  Teuchos::RCP<Epetra_Vector> & GetLevelSetField(){ return cutter_phinp_; }

  virtual void ReadRestart(
      const int step,
      const int lsc_idx = 0
  );

  virtual bool HasMovingInterface(){ return true; }

  virtual void GetInterfaceSlaveMaterial(
    DRT::Element* actele,
    Teuchos::RCP<MAT::Material> & mat
  )
  {
    mat = Teuchos::null;
  }

  XFEM::CouplingBase::LevelSetBooleanType GetBooleanCombination(){ return ls_boolean_type_; }

  //! export row vectors storing geometric quantities to col vectors
  virtual void ExportGeometricQuantities(){};


private:

  void SetConditionsToCopy();

  /// set level-set field implemented in this routine
  double FunctImplementation(
      const int      func_no,
      const double * coords,
      const double t
  );

protected:

  //! Output specific
  Teuchos::RCP<IO::DiscretizationWriter> bg_output_;

  //! @name fluid discretization related state vectors

  //! fluid-dis (bgdis) state vectors for levelset applications
  Teuchos::RCP<Epetra_Vector>    phinp_;


  //@}

  //! @name scatra discretization related state vectors

  //! scatra-dis (cutterdis) state vectors for levelset applications, prepares nonmatching discretizations between scatra and fluid
  Teuchos::RCP<Epetra_Vector>    cutter_phinp_;

  //! The nodal curvature and smoothed gradient of the levelset field. (Stored w.r.t to the scatra-dis = cutter-dis)
  Teuchos::RCP<Epetra_Vector>        curvaturenp_node_;
  Teuchos::RCP<Epetra_MultiVector>   gradphinp_smoothed_node_;
  //Teuchos::RCP<Epetra_MultiVector>   gradphi2np_smoothed_node_;

  //! and column versions
  Teuchos::RCP<Epetra_Vector>        curvaturenp_node_col_;
  Teuchos::RCP<Epetra_MultiVector>   gradphinp_smoothed_node_col_;

  //! boolean operation type on level-set for current ls-field and previous combination of level-set fields
  XFEM::CouplingBase::LevelSetBooleanType ls_boolean_type_;


  // Specify way of creating the projection matrix
  INPAR::XFEM::ProjToSurface projtosurf_;

  //@}

  int bg_nds_phi_;     ///<
  int cutter_nds_phi_; ///<

  double normal_orientation_; ///< correction factor between normal of phi-gradient and normal in xfluid

  bool have_nodematching_dis_; ///< are bgdis and cutterdis node-matching?


};

class LevelSetCouplingBC : public LevelSetCoupling
{
public:

  //! constructor
  explicit LevelSetCouplingBC(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
      const std::string &                 cond_name,///< name of the condition, by which the derived cutter discretization is identified
      Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< full discretization from which the cutter discretization is derived
      const int                           coupling_id,///< id of composite of coupling conditions
      const double                        time,      ///< time
      const int                           step       ///< time step
  );

  //! destructor
  virtual ~LevelSetCouplingBC()
  {
    return;
  }


  virtual void PrepareSolve();

  virtual bool HasMovingInterface();

protected:

  bool has_interface_moved_; ///< did interface move compared to the last time step?

};


/*!
\brief
 */
class LevelSetCouplingWeakDirichlet : public LevelSetCouplingBC
{
public:

  //! constructor
  explicit LevelSetCouplingWeakDirichlet(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
      const std::string &                 cond_name, ///< name of the condition, by which the derived cutter discretization is identified
      Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< full discretization from which the cutter discretization is derived
      const int                           coupling_id,///< id of composite of coupling conditions
      const double                        time,      ///< time
      const int                           step       ///< time step
  ) : LevelSetCouplingBC(bg_dis, cond_name, cond_dis, coupling_id, time, step)
  {
  }


  //! destructor
  virtual ~LevelSetCouplingWeakDirichlet()
  {
    return;
  }

public:

  void EvaluateCouplingConditions(
      LINALG::Matrix<3,1>& ivel,
      LINALG::Matrix<3,1>& itraction,
      const LINALG::Matrix<3,1>& x,
      const DRT::Condition* cond
  );

  void EvaluateCouplingConditionsOldState(
      LINALG::Matrix<3,1>& ivel,
      LINALG::Matrix<3,1>& itraction,
      const LINALG::Matrix<3,1>& x,
      const DRT::Condition* cond
  );

protected:

  //! Initializes configurationmap
  virtual void SetupConfigurationMap();

  //! Updates configurationmap for specific Gausspoint
  virtual void UpdateConfigurationMap_GP( double& kappa_m, double& visc_m, double& visc_s,
      double& visc_stab, double& full_stab, const LINALG::Matrix<3,1>& x, const DRT::Condition* cond);
};


/*!
\brief
 */
class LevelSetCouplingNeumann : public LevelSetCouplingBC
{
public:

  //! constructor
  explicit LevelSetCouplingNeumann(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
      const std::string &                 cond_name, ///< name of the condition, by which the derived cutter discretization is identified
      Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< full discretization from which the cutter discretization is derived
      const int                           coupling_id,///< id of composite of coupling conditions
      const double                        time,      ///< time
      const int                           step       ///< time step
  ) : LevelSetCouplingBC(bg_dis, cond_name, cond_dis, coupling_id, time, step)
  {}


  //! destructor
  virtual ~LevelSetCouplingNeumann()
  {
    return;
  }

public:

  void EvaluateCouplingConditions(
      LINALG::Matrix<3,1>& ivel,
      LINALG::Matrix<3,1>& itraction,
      const LINALG::Matrix<3,1>& x,
      const DRT::Condition* cond
  );

  void EvaluateCouplingConditionsOldState(
      LINALG::Matrix<3,1>& ivel,
      LINALG::Matrix<3,1>& itraction,
      const LINALG::Matrix<3,1>& x,
      const DRT::Condition* cond
  );

};

/*!
\brief
 */
class LevelSetCouplingNavierSlip : public LevelSetCouplingBC
{
public:

  //! constructor
  explicit LevelSetCouplingNavierSlip(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
      const std::string &                 cond_name, ///< name of the condition, by which the derived cutter discretization is identified
      Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< full discretization from which the cutter discretization is derived
      const int                           coupling_id,///< id of composite of coupling conditions
      const double                        time,      ///< time
      const int                           step       ///< time step
  );


  //! destructor
  virtual ~LevelSetCouplingNavierSlip()
  {
    return;
  }

public:

  void EvaluateCouplingConditions(
      LINALG::Matrix<3,1>& ivel,
      LINALG::Matrix<3,1>& itraction,
      const LINALG::Matrix<3,1>& x,
      const DRT::Condition* cond
  );

  void EvaluateCouplingConditionsOldState(
      LINALG::Matrix<3,1>& ivel,
      LINALG::Matrix<3,1>& itraction,
      const LINALG::Matrix<3,1>& x,
      const DRT::Condition* cond
  );

  void GetSlipCoefficient(
      double& slipcoeff,
      const LINALG::Matrix<3,1>& x,
      const DRT::Condition* cond);


  /*!
   Return prescribed velocities and traction vectors for a GNBC boundary condition.
   Also returns the projection matrix (to the plane of the surface) needed for the GNBC condition.
   */
  template <DRT::Element::DiscretizationType DISTYPE,
  class V1, class V2, class X1, class T1, class M1, class M2, class M3>
  void EvaluateCouplingConditions(
      V1&       ivel,                          ///< prescribed velocity at interface
      V2&       itraction,                     ///< prescribed traction at interface
      X1&       x,                             ///< coordinates of gauss point
      const DRT::Condition* cond,              ///< condition prescribed to this surface
      T1&       itraction_jump_matrix,         ///< Laplace-Beltrami matrix for surface tension calculations
      int       eid,                           ///< element ID
      M1&       funct,                         ///< local shape function for Gauss Point (from fluid element)
      M2&       derxy,                         ///< local derivatives of shape function for Gauss Point (from fluid element)
      M3&       normal,                        ///< surface normal of cut element
      double& kappa_m,                         ///< fluid sided weighting
      double& visc_m,                          ///< fluid sided weighting
      double& visc_s                           ///< slave sided dynamic viscosity
  )
{
    EvalProjectionMatrix<DISTYPE>(itraction_jump_matrix,eid,funct,derxy,normal);
    EvaluateCouplingConditions(ivel,itraction,x,cond);

    if(has_neumann_jump_)
    {
      // This is maybe not the most efficient implementation as we evaluate dynvisc as well as the sliplenght twice
      // evaluate interface traction (given by Neumann condition)
      // Add this to the veljump!
      double sliplength = 0.0;
      GetSlipCoefficient(sliplength,x,cond);

      if(sliplength < 0.0)
        dserror("The slip length can not be negative.");

      if ( sliplength != 0.0)
      {
        double sl_visc_fac  = sliplength/(kappa_m*visc_m + (1.0-kappa_m)*visc_s);
        V2 tmp_itraction(true);
        tmp_itraction.MultiplyTN(itraction_jump_matrix,itraction);
         //Project this into tangential direction!!!
        ivel.Update(sl_visc_fac,tmp_itraction,1.0);
      }
      itraction.Clear();
    }

    /*Here one could do a projection of ivel to only point in the tangential direction.
      This would enforce that no spurious velocities occur in the normal direction (i.e. no-penetration always enforced).
      However, these will occur in an XFSI. Thus a solution has to be found which can handle this the best.
    */

    if(forcetangvel_)
    {

#if(1) //If 0 project in normal direction
      LINALG::Matrix<3,1> tmp_ivel(true);
      tmp_ivel.MultiplyTN(itraction_jump_matrix,ivel); //apply Projection matrix from the right. (u_0 * P^t)
      ivel.Update(1.0,tmp_ivel,0.0);
#else
      //If we want only normal direction for testing it can be "hacked" in here.
//      //  itraction_jump_matrix.Print(std::cout);
      LINALG::Matrix<3,3> proj_normal(true);
//      proj_normal.Scale(0.0);
      for(unsigned i=0; i<3; i++) proj_normal(i,i) = 1.0;
      proj_normal.Update(-1.0,itraction_jump_matrix,1.0);

      LINALG::Matrix<3,1> tmp_ivel(true);
      tmp_ivel.MultiplyTN(proj_normal,ivel); //apply Projection matrix from the right. (u_0 * P^n)
      ivel.Update(1.0,tmp_ivel,0.0);
#endif

    }

};

  /*!
   Return a smoothed/non-smoothed tangiential projection of the level set surface.
   */
  template <DRT::Element::DiscretizationType DISTYPE,
  class T1, class M1, class M2, class M3>
  void EvalProjectionMatrix(
      T1&       itraction_jump_matrix,         ///< Projection matrix
      int       eid,                           ///< element ID
      M1&       funct,                         ///< local shape function for Gauss Point (from fluid element)
      M2&       derxy,                         ///< local derivatives of shape function for Gauss Point (from fluid element)
      M3&       normal                         ///< surface normal of cut element
  )
  {

    //Properties of a projection matrix:
    //-------------------------------------------------------------------------
    // 1) P is singular (i.e. not of full rank, no inverse exists).
    // 2) P*P = P
    // 3) P^T = P
    // 4) a*P*a \geq 0 \forall a
    //-------------------------------------------------------------------------

    // number space dimensions for element
    const size_t nsd = DRT::UTILS::DisTypeToDim<DISTYPE>::dim;

    // number of nodes of element
    const size_t nen = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    //Should this be provided as well by the input?
    DRT::Element* actele = cutter_dis_->gElement(eid);

    //Create the identity matrix (Probably not the fastest way...) Might make it global?
    LINALG::Matrix<nsd,nsd> eye(true);
    for(unsigned int i =0; i<nsd; ++i)
      eye(i,i)=1;

    //   Non-smoothed projection matrix
    LINALG::Matrix<nsd,nsd> p_matrix(true);
    if(projtosurf_!=INPAR::XFEM::Proj_smoothed)
    {
      for(unsigned int i =0; i<nsd; ++i)
      {
        for(unsigned int j =0; j<nsd; ++j)
        {
          p_matrix(i,j)          = eye(i,j) - normal(i,0) * normal(j,0);
        }
      }
    }

    //   Smoothed projection matrix
    LINALG::Matrix<nsd,nsd> p_smoothed_matrix(true);
    // Does the chosen Laplace-Beltrami configuration necessitate a smoothed projection matrix?
    if(projtosurf_!=INPAR::XFEM::Proj_normal)
    {
      if(gradphinp_smoothed_node_==Teuchos::null)
        dserror("Smoothed projection only available if a smoothed level set gradient field is existing.");

      LINALG::Matrix<nsd,nen> esmoothedgradphi(true);

      //Make this one global in class? Thus removing unnecessary calls to ExtractMyNodeBasedValues?
      std::vector<double> mygradphi;
      DRT::UTILS::ExtractMyNodeBasedValues(actele,mygradphi,*gradphinp_smoothed_node_);

      for (unsigned int inode=0; inode<nen; ++inode)  // number of nodes
      {
        for(unsigned int idim=0; idim<nsd; ++idim)   // number of dimensions
        {
          (esmoothedgradphi)(idim,inode) = mygradphi[idim+(inode*nsd)];
        }  // end for(idim)
      }
      //Gradients @ GaussPoints
      LINALG::Matrix<nsd,1> gradphi;
      gradphi.Multiply(esmoothedgradphi,funct);

      const double normgradphi=gradphi.Norm2();
      if(normgradphi>1e-9) //1e-9 is set to create a reasonable scaling.
        gradphi.Scale(1.0/normgradphi);
      else
        gradphi.Scale(0.0); //This to catch the cases when gradphi \approx 0

      for(unsigned int i =0; i<nsd; ++i)
      {
        for(unsigned int j =0; j<nsd; ++j)
        {
          p_smoothed_matrix(i,j) = eye(i,j) - gradphi(i,0) * gradphi(j,0);
        }
      }
    }

    // If the smoothed and boundary normal are to be combined this option exists:
    if(projtosurf_==INPAR::XFEM::Proj_normal_smoothed_comb)
    {
      if(gradphinp_smoothed_node_==Teuchos::null)
        dserror("Smoothed projection only available if a smoothed level set gradient field is existing.");

      LINALG::Matrix<nsd,nen> esmoothedgradphi(true);

      //Make this one global in class? Thus removing unnecessary calls to ExtractMyNodeBasedValues?
      std::vector<double> mygradphi;
      DRT::UTILS::ExtractMyNodeBasedValues(actele,mygradphi,*gradphinp_smoothed_node_);

      for (unsigned int inode=0; inode<nen; ++inode)  // number of nodes
      {
        for(unsigned int idim=0; idim<nsd; ++idim)   // number of dimensions
        {
          (esmoothedgradphi)(idim,inode) = mygradphi[idim+(inode*nsd)];
        }  // end for(idim)
      }
      //Gradients @ GaussPoints
      LINALG::Matrix<nsd,1> gradphi;
      gradphi.Multiply(esmoothedgradphi,funct);

      const double normgradphi=gradphi.Norm2();
      if(normgradphi>1e-9) //1e-9 is set to create a reasonable scaling.
        gradphi.Scale(1.0/normgradphi);
      else
        gradphi.Scale(0.0); //This to catch the cases when gradphi \approx 0


      // normal_comb = alpha_n * normal + (1-alpha)*gradphi
      LINALG::Matrix<nsd,1> normal_comb(true);
      double alpha_n = 0.3;
      normal_comb.Update(alpha_n,normal,-(1.0-alpha_n),gradphi);

      const double normnormal_comb=normal_comb.Norm2();
      if(normnormal_comb>1e-9) //1e-9 is set to create a reasonable scaling.
        normal_comb.Scale(1.0/normnormal_comb);
      else
        normal_comb.Scale(0.0); //This to catch the cases when normal_comb \approx 0

      for(unsigned int i =0; i<nsd; ++i)
      {
        for(unsigned int j =0; j<nsd; ++j)
        {
          p_smoothed_matrix(i,j) = eye(i,j) - normal_comb(i,0) * normal_comb(j,0);
        }
      }
    }

    // normal vector extracted from the levelset field without smoothing
    if(projtosurf_==INPAR::XFEM::Proj_normal_phi)
    {

      LINALG::Matrix<nen,1> ephi(true);

      //Make this one global in class? Thus removing unnecessary calls to ExtractMyNodeBasedValues?
      std::vector<double> myphi;
      DRT::UTILS::ExtractMyNodeBasedValues(actele,myphi,*cutter_phinp_);

      for (unsigned int inode=0; inode<nen; ++inode)  // number of nodes
      {
          (ephi)(inode,0) = myphi[inode];
      }

      //Gradients @ GaussPoints
      LINALG::Matrix<nsd,1> gradphi;
      gradphi.Multiply(derxy,ephi);

      const double normgradphi=gradphi.Norm2();
      if(normgradphi>1e-9) //1e-9 is set to create a reasonable scaling.
        gradphi.Scale(1.0/normgradphi);
      else
        gradphi.Scale(0.0); //This to catch the cases when gradphi \approx 0

      for(unsigned int i =0; i<nsd; ++i)
      {
        for(unsigned int j =0; j<nsd; ++j)
        {
          p_smoothed_matrix(i,j) = eye(i,j) - gradphi(i,0) * gradphi(j,0);
        }
      }
    }

    switch(projtosurf_)
    {
    case INPAR::XFEM::Proj_normal:
    {
      itraction_jump_matrix.Update(p_matrix);
      break;
    }
    case INPAR::XFEM::Proj_smoothed:
    case INPAR::XFEM::Proj_normal_smoothed_comb:
    case INPAR::XFEM::Proj_normal_phi:
    {
      itraction_jump_matrix.Update(p_smoothed_matrix);
      break;
    }
    default:
    {
      dserror("This option for a projection matrix %d does not exist. \n",projtosurf_);
      break;
    }
    }

    return;
  }

protected:


  virtual void SetElementConditions();

  void SetElementSpecificConditions(
      std::vector<DRT::Condition* > &  cutterele_cond,
      const std::string &              cond_name,
      const int &                      robin_id
  );

  virtual void SetConditionSpecificParameters();

  void GetConditionByRobinId(
      const std::vector<DRT::Condition*> & mycond,
      const int coupling_id,
      std::vector<DRT::Condition*> & mynewcond
  );

  //! Initializes configurationmap
  virtual void SetupConfigurationMap();

  //! Updates configurationmap for specific Gausspoint
  virtual void UpdateConfigurationMap_GP( double& kappa_m, double& visc_m, double& visc_s,
      double& visc_stab, double& full_stab, const LINALG::Matrix<3,1>& x, const DRT::Condition* cond);

private:

  bool forcetangvel_;
  bool is_constant_sliplength_;
  bool has_neumann_jump_;
  double sliplength_;

  // ID given the Robin Dirichlet/-Neumann conditions (need to match the one given in the "MAIN condition")
  int robin_dirichlet_id_;
  int robin_neumann_id_;

  // Get the condition for the dirichlet and neumann condition associated with the Robin-condition
  std::vector<DRT::Condition* > cutterele_cond_robin_dirichlet_;
  std::vector<DRT::Condition* > cutterele_cond_robin_neumann_;


}; //End LevelSetCouplingNavierSlip


class LevelSetCouplingTwoPhase : public LevelSetCoupling
{
public:

  //! constructor
  explicit LevelSetCouplingTwoPhase(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,    ///< background discretization
      const std::string &                 cond_name, ///< name of the condition, by which the derived cutter discretization is identified
      Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< discretization from which the cutter discretization is derived
      const int                           coupling_id,///< id of composite of coupling conditions
      const double                        time,      ///< time
      const int                           step       ///< time step
  );

  virtual ~LevelSetCouplingTwoPhase()
  {
    return;
  }

  virtual void SetCutterDiscretization();

  virtual void InitStateVectors();

  /// set levelset field by given vector
  void WriteAccess_GeometricQuantities(
      Teuchos::RCP<Epetra_Vector> &      scalaraf,
      Teuchos::RCP<Epetra_MultiVector> & smoothed_gradphiaf,
      Teuchos::RCP<Epetra_Vector> &      curvatureaf
  );

  //! export row vectors storing geometric quantities to col vectors
  virtual void ExportGeometricQuantities();

  //! set twophase flow specific parameters from the input file
  virtual void SetConditionSpecificParameters();


  /// set material pointer for coupling slave side
  void GetInterfaceSlaveMaterial(
    DRT::Element* actele,
    Teuchos::RCP<MAT::Material> & mat
  );

  template <DRT::Element::DiscretizationType DISTYPE,
  class M1, class M2>
  void EvaluateCurvature(
      double&   icurvature,         ///< curvature to be computed
      int       eid,                ///< element ID
      M1&       funct,              ///< local shape function for Gauss Point (from fluid element)
      M2&       derxy               ///< local derivatives of shape function for Gauss Point (from fluid element)
  )
  {
    CheckForValidVectors();

    if(surftensapprox_==INPAR::TWOPHASE::surface_tension_approx_none)
    {
      icurvature=0.0;
      return;
    }

    if(surftensapprox_==INPAR::TWOPHASE::surface_tension_approx_laplacebeltrami)
      dserror("wrong routine for laplace beltrami surface tension approximation!");


    // number space dimensions for element
    const size_t nsd = DRT::UTILS::DisTypeToDim<DISTYPE>::dim;

    // number of nodes of element
    const size_t nen = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    //Should this be provided as well by the input?
    DRT::Element* cutter_ele = cutter_dis_->gElement(eid);

    // smoothed normal at cutter element nodes, the Gaussian point lies in
    Epetra_SerialDenseMatrix esmoothedgradphi(nsd, nen);
    Epetra_SerialDenseMatrix esmoothedcurvature(nen, 1);

    LINALG::Matrix<nsd,nen> esmoothedgradphi_T(esmoothedgradphi, View);
    LINALG::Matrix<nen,1>   esmoothedcurvature_T(esmoothedcurvature, View);


    GetSmoothedQuantitiesAtElement(
        esmoothedgradphi,
        esmoothedcurvature,
        cutter_ele,
        gradphinp_smoothed_node_col_,
        curvaturenp_node_col_,
        cutter_dis_,
        cutter_nds_phi_);

    icurvature = InterpolateCurvature<DISTYPE>(
        funct,
        derxy,
        esmoothedgradphi_T,
        esmoothedcurvature_T,
        nen);

    return;
  }

  inline void CheckForValidVectors()
  {
    if(!col_vectors_valid_)
      dserror("you need to call ExportGeometricQuantities() first");
  }

  template <DRT::Element::DiscretizationType DISTYPE, class M1, class M2, class M3, class M4>
  double InterpolateCurvature(
      const M1 & funct,
      const M2 & derxy,
      const M3 & esmoothedgradphi_T,
      const M4 & esmoothedcurvature_T,
      const int numnode
  )
  {
    // number space dimensions for element
    const size_t nsd = DRT::UTILS::DisTypeToDim<DISTYPE>::dim;

    double icurvature = 0.0;

    if(surftensapprox_==INPAR::TWOPHASE::surface_tension_approx_divgrad_normal)
    {
      LINALG::Matrix<nsd,1> grad_phi(false);
      grad_phi.Multiply(esmoothedgradphi_T,funct);

      LINALG::Matrix<nsd*nsd,1> grad_phi2(true);

      // compute 2nd phi derivatives based on smoothed first derivatives
      DeriveSmoothedGradPhi(grad_phi2, derxy, esmoothedgradphi_T, numnode);


      //=========================================== calculate curvature =========================
      double grad_phi_norm = grad_phi.Norm2();

      // check norm of normal gradient
      if (fabs(grad_phi_norm) < 1.0E-5)// 'ngradnorm' == 0.0
      {
        std::cout << "grad phi is small -> set to 1.0E12" << grad_phi_norm << std::endl;
        // phi gradient too small -> there must be a local max or min in the level-set field
        // set curvature to a large value (it will be cut off based on the element size)
        icurvature = 1.0E12;
      }
      else
      {
        icurvature = EvaluateCurvatureFormula(grad_phi_norm, grad_phi, grad_phi2 );
      }
    }
    else if(surftensapprox_==INPAR::TWOPHASE::surface_tension_approx_nodal_curvature)
    {
      //++++++++++++++++++++++
      //       WARNING
      //++++++++++++++++++++++
      //This curvature calculation needs to be checked!!!!
      // Should there really be a minus here?!?!
      icurvature= (-1.0) * funct.Dot(esmoothedcurvature_T); //No idea why there is a minus here?!
    }
    else
      dserror("unsupported type of curvature computation!");

    return icurvature;
  }

  template <class M1, class M2, class M3>
  void DeriveSmoothedGradPhi(
      M1 & grad_phi2,
      const M2 & derxy,
      const M3 & esmoothedgradphi,
      const int numnode
  )
  {
    // loop over nodes
    for(int i = 0; i< numnode; ++i)
    {
      if(nsd_==3)
      {
        grad_phi2(0) += derxy(0,i)*esmoothedgradphi(0,i); // ,xx
        grad_phi2(1) += derxy(1,i)*esmoothedgradphi(1,i); // ,yy
        grad_phi2(2) += derxy(2,i)*esmoothedgradphi(2,i); // ,zz
        grad_phi2(3) += derxy(1,i)*esmoothedgradphi(0,i); // ,xy
        grad_phi2(4) += derxy(2,i)*esmoothedgradphi(0,i); // ,xz
        grad_phi2(5) += derxy(2,i)*esmoothedgradphi(1,i); // ,yz
        grad_phi2(6) += derxy(0,i)*esmoothedgradphi(1,i); // ,yx
        grad_phi2(7) += derxy(0,i)*esmoothedgradphi(2,i); // ,zx
        grad_phi2(8) += derxy(1,i)*esmoothedgradphi(2,i); // ,zy
      }
      else
        dserror("Only 3D-elements are supported as of now.");
    }
  }

  template <class M1, class M2>
  double EvaluateCurvatureFormula(
      double & grad_phi_norm,
      M1 & grad_phi,
      M2 & grad_phi2
  )
  {
    double val = grad_phi_norm*grad_phi_norm*grad_phi_norm;
    double invval = 1.0 / val;
    double curvature =
        -invval*(
            grad_phi(0)*grad_phi(0)*grad_phi2(0)
          + grad_phi(1)*grad_phi(1)*grad_phi2(1)
          + grad_phi(2)*grad_phi(2)*grad_phi2(2)
        )
        -invval*(
            grad_phi(0)*grad_phi(1)*( grad_phi2(3) + grad_phi2(6) )
          + grad_phi(0)*grad_phi(2)*( grad_phi2(4) + grad_phi2(7) )
          + grad_phi(1)*grad_phi(2)*( grad_phi2(5) + grad_phi2(8) )
        )
        +1.0/grad_phi_norm*(
            grad_phi2(0) + grad_phi2(1) + grad_phi2(2)
        );

    return curvature;
  }

  template <DRT::Element::DiscretizationType DISTYPE,
  class T1, class M1, class M2>
  void EvalProjectionMatrix(
      T1&       itraction_jump_matrix,         ///< Laplace-Beltrami matrix for surface tension calculations
      int       eid,                           ///< element ID
      M1&       funct,                         ///< local shape function for Gauss Point (from fluid element)
      M2&       normal                         ///< surface normal of cut element
  )
  {
    CheckForValidVectors();

    // number space dimensions for element
    const size_t nsd = DRT::UTILS::DisTypeToDim<DISTYPE>::dim;

    // Should this be provided as well by the input?
    DRT::Element* cutter_ele = cutter_dis_->gElement(eid);

    LINALG::Matrix<nsd,nsd> p_matrix(false);
    LINALG::Matrix<nsd,nsd> p_smoothed_matrix(false);

    //----------------------------------------------------
    // Non-smoothed projection matrix
    if(laplacebeltrami_!=INPAR::TWOPHASE::matrix_smoothed)
      SetupProjectionMatrix( p_matrix, normal);

    //----------------------------------------------------
    // Smoothed projection matrix
    if(laplacebeltrami_!=INPAR::TWOPHASE::matrix_non_smoothed)
    {
      CheckForValidVectors();

      LINALG::Matrix<nsd,1> gradphi(true);
      double curvature = 0.0;

      GetSmoothedQuantitiesAtPoint<DISTYPE>(
          gradphi,
          curvature,
          cutter_ele,
          funct
      );

      LINALG::Matrix<nsd,1> normal_smoothed(true);
      GetNormalFromLevelsetGradient(normal_smoothed, gradphi);

      SetupProjectionMatrix( p_smoothed_matrix, normal_smoothed);
    }

    //  -----------------------------------------------------------
    //
    // HERE DEPENDING ON HOW THE PROJECTION SHOULD BE CALCULATED
    // different ways of the projection matrix should be tested.
    //
    // Might, only want to use matrix_mixed_smoothed option as this is the most promising!
    //   Check for stabilized Laplace-Beltrami option!
    //    Burman, Erik and Hansbo, Peter and Larson, Mats G
    //    A stabilized cut finite element method for partial differential equations on surfaces: The Laplace--Beltrami operator
    //    Computer Methods in Applied Mechanics and Engineering
    //       2015
    //
    //  -----------------------------------------------------------

    SetupConcatenatedProjectionMatrix(
        itraction_jump_matrix,
        p_matrix,
        p_smoothed_matrix
    );

    return;
  }

  template <class M1, class M2>
  inline void SetupProjectionMatrix(
      M1 & proj_matrix,
      const M2 & normal
  )
  {
    double n_j = 0.0;

    for(unsigned int j=0; j<nsd_; ++j)
    {
      n_j = normal(j,0);

      for(unsigned int i=0; i<nsd_; ++i)
      {
        proj_matrix(i,j) = ((i==j) ? 1.0 : 0.0) - normal(i,0) * n_j;
      }
    }

    return;
  }


  template <class M1>
  inline void SetupConcatenatedProjectionMatrix(
      M1 & itraction_jump_matrix,
      const M1 & p_matrix,
      const M1 & p_smoothed_matrix
  )
  {
    switch(laplacebeltrami_)
    {
    case INPAR::TWOPHASE::matrix_non_smoothed:
    {
      itraction_jump_matrix.Update(1.0, p_matrix, 0.0);
      break;
    }
    case INPAR::TWOPHASE::matrix_mixed_smoothed:
    {
      // As suggested by:
      // S. Groß and A. Reusken. Finite element discretization error analysis of a surface
      // tension force in two-phase incompressible flows. SIAM J. Numer. Anal., 45(4):1679–1700, 2007.
      itraction_jump_matrix.Multiply(p_matrix,p_smoothed_matrix);
      break;
    }
    case INPAR::TWOPHASE::matrix_smoothed:
    {
      itraction_jump_matrix.Update(1.0, p_smoothed_matrix, 0.0);
      break;
    }
    default:
    {
      dserror("The Laplace-Beltrami option does not exist. Please review your input.\n");
      break;
    }
    }
  }

  template <class M1>
  void GetNormalFromLevelsetGradient(
      M1 & normal_smoothed,
      const M1 & gradphi
  )
  {
    const double normgradphi=gradphi.Norm2();
    double scale = 1.0;

    if(normgradphi>1e-9) //1e-9 is set to create a reasonable scaling.
      scale=1.0/normgradphi;
    else
      scale=0.0; //This to catch the cases when gradphi \approx 0

    scale*=normal_orientation_;

    normal_smoothed.Update(scale, gradphi, 0.0);
  }

  template <DRT::Element::DiscretizationType DISTYPE,
  class T1, class T2, class M1, class M2, class M3>
  void EvaluateTractionDiscontinuity(
      T1&       itraction_jump,         ///< Traction jump vector
      T2&       itraction_jump_matrix,  ///< Laplace-Beltrami matrix for surface tension calculations
      int       eid,                    ///< element ID
      M1&       funct,                  ///< local shape function for Gauss Point (from fluid element)
      M2&       derxy,                  ///< local derivatives of shape function for Gauss Point (from fluid element)
      M3&       normal,                 ///< surface normal of cut element
      double&   surf_coeff              ///< surface tension coefficient gamma
  )
  {

    //Maybe instead of ifs -> template the functions with surftensapprox_.

    if(not (surftensapprox_==INPAR::TWOPHASE::surface_tension_approx_laplacebeltrami))
    {

      // where [*] = (*)^m - (*)^s = (*)^+ - (*)^-
      // n = n^m = n^+
      // [sigma*n] = gamma * curv * n   with curv = div(grad(phi)/||grad(phi)||)

      if(surftensapprox_==INPAR::TWOPHASE::surface_tension_approx_fixed_curvature)
      {
        itraction_jump.Update(surf_coeff,normal,0.0);
        return;
      }

      // Calculations for evaluation of curvature
      double curvature_int = 0.0;
      EvaluateCurvature<DISTYPE>(curvature_int, eid, funct, derxy);

      itraction_jump.Update(curvature_int * surf_coeff, normal, 0.0);
    }
    else
    {
      //itraction_jump_matrix.Scale(0.0);

      // P = (\delta_{ij} - n_{i}n{j})
      // Calculations for evaluation with Laplace-Beltrami operator
      EvalProjectionMatrix<DISTYPE>(itraction_jump_matrix, eid, funct, normal);
      itraction_jump_matrix.Scale(surf_coeff);
    }

    return;
  }


  virtual const Teuchos::RCP<const Epetra_Vector> ComputeTransportVelocity(
      const Teuchos::RCP<GEO::CutWizard> & wizard,                   ///< the cut wizard
      const Teuchos::RCP<const Epetra_Vector>& convective_velocity   ///< the convective fluid velocity based on initial dofmap (just velocity)
  );

  /// get viscosity of the slave fluid
  void GetViscositySlave(
      DRT::Element * coup_ele,                   ///< xfluid ele
      double& visc_s                             ///< viscosity slavesided
      );

  /// get scaling of the master side for penalty (viscosity, E-modulus for solids)
  virtual void GetPenaltyScalingSlave(
      DRT::Element * coup_ele,                      ///< xfluid ele
      double& penscaling_s)                         ///< penalty scaling slavesided
  {
    GetViscositySlave(coup_ele,penscaling_s);
  }

  /// get coupling specific weighting paramters
  virtual void  GetCouplingSpecificAverageWeights(
      DRT::Element * xfele,                      ///< xfluid ele
      DRT::Element * coup_ele,                   ///< coup_ele ele
      double & kappa_m
      );

protected:

  /// Check ifSetParameters_SurfaceTension() has been called
  inline void CheckInit_SurfaceTension() const
  {
    if (not surfacetension_init_)
      dserror("Call SetParameters_SurfaceTension() first!");
  }

  virtual void SetParameters_Physical( DRT::Condition* cond ){};

  virtual void SetFlags_InterfaceTransport( DRT::Condition* cond );

  virtual void SetFlags_GeometricQuantities();

  void SetParameters_SurfaceTension();

  //! get the smoothed level set gradient at a given node (not necessarily normalized to one)
  void GetSmoothedQuantitiesAtNode(
      Epetra_SerialDenseVector & normal,
      double & curvature,
      const DRT::Node* node,
      const Teuchos::RCP<const Epetra_MultiVector> & gradphinp_smoothed_node_col,
      const Teuchos::RCP<const Epetra_MultiVector> & curvaturenp_node_col,
      Teuchos::RCP<DRT::Discretization> & dis,
      const int nds
  );


  void GetSmoothedQuantitiesAtElement(
      Epetra_SerialDenseMatrix & normal,
      Epetra_SerialDenseMatrix & curvature,
      const DRT::Element* element,
      const Teuchos::RCP<const Epetra_MultiVector> & gradphinp_smoothed_node_col,
      const Teuchos::RCP<const Epetra_MultiVector> & curvaturenp_node_col,
      Teuchos::RCP<DRT::Discretization> & dis,
      const int nds
  );

  template<DRT::Element::DiscretizationType DISTYPE,
    class M1, class M2>
    void GetSmoothedQuantitiesAtPoint(
        M1 & gradphi,
        double & curvature,
        const DRT::Element * cutter_ele,
        const M2 & funct
    )
  {
    // number space dimensions for element
    const size_t nsd = DRT::UTILS::DisTypeToDim<DISTYPE>::dim;

    // number of nodes of element
    const size_t nen = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    Epetra_SerialDenseMatrix esmoothedgradphi(nsd, nen);
    Epetra_SerialDenseMatrix esmoothedcurvature(nen, 1);

    LINALG::Matrix<nsd,nen> esmoothedgradphi_T(esmoothedgradphi, View);
    LINALG::Matrix<nen,1>   esmoothedcurvature_T(esmoothedcurvature, View);

    GetSmoothedQuantitiesAtElement(
        esmoothedgradphi,
        esmoothedcurvature,
        cutter_ele,
        gradphinp_smoothed_node_col_,
        curvaturenp_node_col_,
        cutter_dis_,
        cutter_nds_phi_);

    //Gradients @ GaussPoints
    if(require_smoothedgradphi_)
      gradphi.Multiply(esmoothedgradphi_T,funct);

    if(require_nodalcurvature_)
      curvature = esmoothedcurvature_T.Dot(funct);

  }

  //! Initializes configurationmap
  virtual void SetupConfigurationMap();

  //! Updates configurationmap for specific Gausspoint
  virtual void UpdateConfigurationMap_GP( double& kappa_m, double& visc_m, double& visc_s,
      double& visc_stab, double& full_stab, const LINALG::Matrix<3,1>& x, const DRT::Condition* cond);

  //! do normalization or a manipulation of the smoothed normal vector
  bool RescaleNormal( Epetra_SerialDenseVector & normal );

  virtual void ComputeRelativeTransportVelocity(
      Epetra_SerialDenseVector & flvelrel,
      const GEO::CUT::Point::PointPosition & position,
      const Epetra_SerialDenseVector & nvec,
      const double & curv
  );

  virtual void Output(
      const int step,
      const double time,
      const bool write_restart_data,
      const int lsc_idx = 0
  );

  /// Read Restart data for ScaTra coupled level set
  virtual void ReadRestart(
      const int step,
      const int lsc_idx = 0
  );


  Teuchos::RCP<Epetra_Vector>    cutter_transport_vel_;

  const Teuchos::ParameterList& twophasedyn_;

  INPAR::TWOPHASE::SurfaceTensionApprox surftensapprox_;  ///< type of surface tension approximation
  INPAR::TWOPHASE::LaplaceBeltramiCalc  laplacebeltrami_; ///< type of laplace beltrami surface tension operator

  bool surfacetension_init_;     ///< have surface tension specific quantities been initialized?

  bool require_smoothedgradphi_; ///< do we need to compute a smoothed normal field?
  bool require_nodalcurvature_;  ///< do we require to compute nodal curvatures?

  INPAR::TWOPHASE::Transport_Directions transport_direction_; ///< use just normal part of convective velocity for transporting the interface


  int cutter_nds_vel_; // dofset id for transport velocity in cutter discretization

  bool col_vectors_valid_; ///< are column vectors uptodate?

};

class LevelSetCouplingCombustion : public LevelSetCouplingTwoPhase
{
public:

  //! constructor
  explicit LevelSetCouplingCombustion(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,    ///< background discretization
      const std::string &                 cond_name, ///< name of the condition, by which the derived cutter discretization is identified
      Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< discretization from which the cutter discretization is derived
      const int                           coupling_id,///< id of composite of coupling conditions
      const double                        time,      ///< time
      const int                           step       ///< time step
  );

  //! destructor
  virtual ~LevelSetCouplingCombustion()
  {
    return;
  }


  //! evaluate the element's coupling condition for the combustion interface
  template <DRT::Element::DiscretizationType DISTYPE,
  class T1, class M1, class M2, class M3>
  void EvaluateCouplingConditions(
      T1&       ivelint_jump,           ///< velocity jump vector
      T1&       itraction_jump,         ///< Traction jump vector
      const DRT::Condition* cond,       ///< coupling condition
      int       eid,                    ///< element ID
      M1&       funct,                  ///< local shape function for Gauss Point (from fluid element)
      M2&       derxy,                  ///< local derivatives of shape function for Gauss Point (from fluid element)
      M3&       normal,                 ///< surface normal of cut element
      const double densaf_master,       ///< master density
      const double densaf_slave         ///< slave density
  )
  {
    CheckForValidVectors();

    if(!have_nodematching_dis_) // TODO: use the information from the cutter dis and not from the bg_dis!!! template over cutterdis!
      dserror("you need to evaluate at the Gaussian point and so interpolate the smoothed geometric quantities from the cutterdis");

//    const int bg_eid     = eid;
    const int cutter_eid = eid;

    // number space dimensions for element
    const size_t nsd = DRT::UTILS::DisTypeToDim<DISTYPE>::dim;

    DRT::Element* cutter_ele = cutter_dis_->gElement(cutter_eid);
    if(cutter_ele == NULL)
      dserror("element not available");

    LINALG::Matrix<nsd,1> gradphi(true);
    double curvature = 0.0;

    GetSmoothedQuantitiesAtPoint<DISTYPE>(
        gradphi,
        curvature,
        cutter_ele,
        funct
    );

    // smoothed normal interpolated at Gaussian point
    LINALG::Matrix<nsd,1> normal_smoothed(true);
    GetNormalFromLevelsetGradient(normal_smoothed, gradphi);

    if(transport_curvature_)
      dserror("how to account for curvature here?");

    // convention from thesis by B. Schott
    // where [*] = (*)^i - (*)^j = (*)^- - (*)^+ (i<j)
    //           = (*)^m - (*)^s
    // n = n^{ij} = n^i = -n^j = n^m = -n^s = n^+ = -n^-
    // note: n^+ = -grad(phi)

    // M = -(rho^i * sl) mass flow rate
    // [u]       =  -M   * [rho^(-1)] * n_smoothed = -M   * [rho^(-1)] (-grad(phi))
    // [sigma*n] =   M*M * [rho^(-1)] * n

    // jump of inverse densities
    const double dens_times_speed  =  densaf_master*laminar_flamespeed_;
    const double mass_flow_rate    = -dens_times_speed;
    const double inv_dens_jump     =  1.0/densaf_master - 1.0/densaf_slave;
    const double vel_jump_fac      = -1.0*mass_flow_rate*inv_dens_jump;
    const double traction_jump_fac =  mass_flow_rate * mass_flow_rate * inv_dens_jump;

    ivelint_jump.Update  ( vel_jump_fac,      normal_smoothed, 0.0);
    itraction_jump.Update( traction_jump_fac, normal         , 0.0);

    //TODO: add dependency on interface curvature!

    //std::cout << "ivelint_jump" << ivelint_jump << std::endl;
    //std::cout << "itraction_jump" << itraction_jump << std::endl;
    return;
  }



  //! get the material densities of the two adjacent phases
  void GetPhaseDensities(
      DRT::Element* ele,
      double & rhoplus,
      double & rhominus
  );

protected:

  virtual void SetParameters_Physical( DRT::Condition* cond );

  virtual void SetFlags_InterfaceTransport( DRT::Condition* cond );

  virtual void SetFlags_GeometricQuantities();

  virtual void ComputeRelativeTransportVelocity(
      Epetra_SerialDenseVector & flvelrel,
      const GEO::CUT::Point::PointPosition & position,
      const Epetra_SerialDenseVector & nvec,
      const double & curv
  );


private:

  double laminar_flamespeed_;  /// the laminar flame speed
  double mol_diffusivity_;     /// the molecular diffusivity
  double markstein_length_;    /// the markstein length

  bool transport_curvature_; ///< account for curvature in computing the transport velocity

};

} // namespace XFEM

#endif /*XFEM_COUPLING_LEVELSET_H_*/

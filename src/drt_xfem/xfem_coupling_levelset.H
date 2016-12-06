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
      const int                           coupling_id,///< id of composite of coupling conditions
      const double                        time,      ///< time
      const int                           step       ///< time step
  );


  //! destructor
  virtual ~LevelSetCoupling()
  {
    return;
  }


  /// set levelset field by function, return if interface moved compared to last time step
  bool SetLevelSetField(const double time);

  virtual void SetLevelSetBooleanType();

  virtual bool ApplyComplementaryOperator();

  virtual void Output(
      const int step,
      const double time,
      const bool write_restart_data,
      const int lsc_idx = 0
  );

  Teuchos::RCP<Epetra_Vector> & GetLevelSetField(){ return phinp_; }

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


private:

  void SetConditionsToCopy();

  /// set level-set field implemented in this routine
  double FunctImplementation(
      const int      func_no,
      const double * coords,
      const double t
  );

protected:


  //! @name state vectors

  //! scatra-dis (cutterdis) state vectors for levelset applications, prepares nonmatching discretizations between scatra and fluid
  Teuchos::RCP<Epetra_Vector>    phinp_;

  //! The nodal curvature and smoothed gradient of the levelset field.
  Teuchos::RCP<Epetra_Vector>        curvaturenp_node_;
  Teuchos::RCP<Epetra_MultiVector>   gradphinp_smoothed_node_;
  Teuchos::RCP<Epetra_MultiVector>   gradphi2np_smoothed_node_;

  //! boolean operation type on level-set for current ls-field and previous combination of level-set fields
  XFEM::CouplingBase::LevelSetBooleanType ls_boolean_type_;

  //! Output specific
  Teuchos::RCP<IO::DiscretizationWriter> ls_output_;

  // Specify way of creating the projection matrix
  INPAR::XFEM::ProjToSurface projtosurf_;

  //@}

};

class LevelSetCouplingBC : public LevelSetCoupling
{
public:

  //! constructor
  explicit LevelSetCouplingBC(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
      const std::string &                 cond_name,///< name of the condition, by which the derived cutter discretization is identified
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

private:

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
      const int                           coupling_id,///< id of composite of coupling conditions
      const double                        time,      ///< time
      const int                           step       ///< time step
  ) : LevelSetCouplingBC(bg_dis, cond_name, coupling_id, time, step)
  {
    // initialize the configuration map
    InitConfigurationMap();
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
  void InitConfigurationMap();

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
      const int                           coupling_id,///< id of composite of coupling conditions
      const double                        time,      ///< time
      const int                           step       ///< time step
  ) : LevelSetCouplingBC(bg_dis, cond_name, coupling_id, time, step)
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
      const int                           coupling_id,///< id of composite of coupling conditions
      const double                        time,      ///< time
      const int                           step       ///< time step
  ) : LevelSetCouplingBC(bg_dis, cond_name, coupling_id, time, step)
  {
    SetConditionSpecificParameters();

    //Might save these differently....
    //  For instance in a std::map< robin_id, DRT::Condition* > as is done in NavSlip-Mesh
    SetElementSpecificConditions(cutterele_cond_robin_dirichlet_,"XFEMRobinDirichletVol",robin_dirichlet_id_);
    if(has_neumann_jump_)
      SetElementSpecificConditions(cutterele_cond_robin_neumann_,"XFEMRobinNeumannVol",robin_neumann_id_);

    // initialize the configuration map
    InitConfigurationMap();
  }


  //! destructor
  virtual ~LevelSetCouplingNavierSlip()
  {
    return;
  }

public:

  void SetElementSpecificConditions(
      std::vector<DRT::Condition* > &  cutterle_cond,
      const std::string &              cond_name,
      const int &                      robin_id
  );

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
      DRT::UTILS::ExtractMyNodeBasedValues(actele,myphi,*phinp_);

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
    {
      itraction_jump_matrix.Update(p_smoothed_matrix);
      break;
    }
    case INPAR::XFEM::Proj_normal_smoothed_comb:
    {
      itraction_jump_matrix.Update(p_smoothed_matrix);
      break;
    }
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

  void SetConditionSpecificParameters();

  void GetConditionByRobinId(
      const std::vector<DRT::Condition*> & mycond,
      const int coupling_id,
      std::vector<DRT::Condition*> & mynewcond
  );

  //! Initializes configurationmap
  void InitConfigurationMap();

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
      const int                           coupling_id,///< id of composite of coupling conditions
      const double                        time,      ///< time
      const int                           step       ///< time step
  ) : LevelSetCoupling(bg_dis, cond_name, coupling_id, time, step),
      surfacetension_init_(false),
      surftensapprox_(INPAR::TWOPHASE::surface_tension_approx_none),
      laplacebeltrami_(INPAR::TWOPHASE::matrix_mixed_smoothed)
  {
    // initialize the configuration map
    InitConfigurationMap();
  }

  virtual ~LevelSetCouplingTwoPhase()
  {}

  /// set levelset field by given vector
  void SetLevelSetField(
      Teuchos::RCP<const Epetra_Vector> scalaraf,
      Teuchos::RCP<const Epetra_Vector> curvatureaf,
      Teuchos::RCP<Epetra_MultiVector>  smoothed_gradphiaf,
      Teuchos::RCP<DRT::Discretization> scatradis
  );

  /// initialize surface tension specific parameters
  void SetSurfaceTensionSpecifcParameters(
      INPAR::TWOPHASE::SurfaceTensionApprox surftensapprox,
      INPAR::TWOPHASE::LaplaceBeltramiCalc  laplacebeltrami
      );

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
//      M1&       posXiDomain,        ///< position of gaussian point for evaluating of surface tension terms
      M1&       funct,              ///< local shape function for Gauss Point (from fluid element)
      M2&       derxy               ///< local derivatives of shape function for Gauss Point (from fluid element)
  )
  {

    // number space dimensions for element
    const size_t nsd = DRT::UTILS::DisTypeToDim<DISTYPE>::dim;

    // number of nodes of element
    const size_t nen = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // get element-wise velocity/pressure field for current time step
    LINALG::Matrix<nen,1> ecurvature(true);
    LINALG::Matrix<nsd,nen> esmoothedgradphi(true);
    LINALG::Matrix<nsd*nsd,nen> esmoothedgradphi2(true);

    //Should this be provided as well by the input?
    DRT::Element* actele = cutter_dis_->gElement(eid);

    if(curvaturenp_node_!=Teuchos::null)
    {
      //Make this one global in class? Thus removing unnecessary calls to ExtractMyNodeBasedValues?
      {
        std::vector<double> mymatrix_curv;
        DRT::UTILS::ExtractMyNodeBasedValues(actele,mymatrix_curv,*curvaturenp_node_);

        for (unsigned int inode=0; inode<nen; ++inode)  // number of nodes
        {
          (ecurvature)(inode,0) = mymatrix_curv[(inode*1)];
        }
      }
      //++++++++++++++++++++++
      //       WARNING
      //++++++++++++++++++++++
      //This curvature calculation needs to be checked!!!!
      // Should there really be a minus here?!?!
      icurvature= (-1.0) * funct.Dot(ecurvature); //No idea why there is a minus here?!
    }
    else if(gradphinp_smoothed_node_!=Teuchos::null)
    {
      //Make this one global in class? Thus removing unnecessary calls to ExtractMyNodeBasedValues?
      {
        std::vector<double> mygradphi;
        DRT::UTILS::ExtractMyNodeBasedValues(actele,mygradphi,*gradphinp_smoothed_node_);

        for (unsigned int inode=0; inode<nen; ++inode)  // number of nodes
        {
          for(unsigned int idim=0; idim<nsd; ++idim)   // number of dimensions
          {
            (esmoothedgradphi)(idim,inode) = mygradphi[idim+(inode*nsd)];
          }  // end for(idim)
        }
      }
      //Gradients @ GaussPoints
      LINALG::Matrix<nsd,1> grad_phi;
      LINALG::Matrix<nsd*nsd,1> grad_phi2;

      grad_phi.Multiply(esmoothedgradphi,funct);
      if(gradphi2np_smoothed_node_==Teuchos::null)
      {
        // loop over nodes
        for(size_t i = 0; i< nen; ++i)
        {
          if(nsd==3)
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
      else
        dserror("Smoothed second derivatives are not supported as of yet.");

      //=========================================== calculate curvature =========================

      double grad_phi_norm = grad_phi.Norm2();

      // check norm of normal gradient
      if (fabs(grad_phi_norm) < 1.0E-5)// 'ngradnorm' == 0.0
      {
        //std::cout << "grad phi is small -> set to 1.0E12" << grad_phi_norm << std::endl;
        // phi gradient too small -> there must be a local max or min in the level-set field
        // set curvature to a large value (it will be cut off based on the element size)
        icurvature = 1.0E12;
      }
      else
      {
        double val = grad_phi_norm*grad_phi_norm*grad_phi_norm;
        double invval = 1.0 / val;
        icurvature = -invval*(  grad_phi(0)*grad_phi(0)*grad_phi2(0)
            + grad_phi(1)*grad_phi(1)*grad_phi2(1)
            + grad_phi(2)*grad_phi(2)*grad_phi2(2)  )
            -invval*(  grad_phi(0)*grad_phi(1)*( grad_phi2(3) + grad_phi2(6) )
                + grad_phi(0)*grad_phi(2)*( grad_phi2(4) + grad_phi2(7) )
                + grad_phi(1)*grad_phi(2)*( grad_phi2(5) + grad_phi2(8)) )
                +1.0/grad_phi_norm * ( grad_phi2(0) + grad_phi2(1) + grad_phi2(2) );
      }
    }
    else
      dserror("Why are you evaluating the surface tension when you've specified to have none? Check your GAMMA in Material input.");

    return;
  }

  template <DRT::Element::DiscretizationType DISTYPE,
  class T1, class M1, class M2, class M3>
  void EvalProjectionMatrix(
      T1&       itraction_jump_matrix,         ///< Laplace-Beltrami matrix for surface tension calculations
      int       eid,                           ///< element ID
      M1&       funct,                         ///< local shape function for Gauss Point (from fluid element)
      M2&       derxy,                         ///< local derivatives of shape function for Gauss Point (from fluid element)
      M3&       normal                         ///< surface normal of cut element
  )
  {

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
    if(laplacebeltrami_!=INPAR::TWOPHASE::matrix_smoothed)
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
    if(laplacebeltrami_!=INPAR::TWOPHASE::matrix_non_smoothed)
    {
      if(gradphinp_smoothed_node_==Teuchos::null)
        dserror("Laplace-Beltrami with smoothed gradient chosen but no smoothed normal vector provided!!");

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

    switch(laplacebeltrami_)
    {
     case INPAR::TWOPHASE::matrix_non_smoothed:
     {
       itraction_jump_matrix.Update(p_matrix);
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
       itraction_jump_matrix.Update(p_smoothed_matrix);
       break;
     }
     default:
     {
       dserror("The Laplace-Beltrami option does not exist. Please review your input.\n");
       break;
     }
     }

    return;
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
      EvaluateCurvature<DISTYPE>(curvature_int,eid,funct,derxy);

      itraction_jump.Update(curvature_int * surf_coeff, normal, 0.0);
    }
    else
    {
      //itraction_jump_matrix.Scale(0.0);

      // P = (\delta_{ij} - n_{i}n{j})
      // Calculations for evaluation with Laplace-Beltrami operator
      EvalProjectionMatrix<DISTYPE>(itraction_jump_matrix,eid, funct, derxy, normal);
      itraction_jump_matrix.Scale(surf_coeff);
    }

    return;
  }

  void ReadRestart(
      const int step,
      const int lsc_idx = 0
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

  //! Initializes configurationmap
  void InitConfigurationMap();

  //! Updates configurationmap for specific Gausspoint
  virtual void UpdateConfigurationMap_GP( double& kappa_m, double& visc_m, double& visc_s,
      double& visc_stab, double& full_stab, const LINALG::Matrix<3,1>& x, const DRT::Condition* cond);

private:

  bool surfacetension_init_;
  INPAR::TWOPHASE::SurfaceTensionApprox surftensapprox_;
  INPAR::TWOPHASE::LaplaceBeltramiCalc  laplacebeltrami_;

};

} // namespace XFEM

#endif /*XFEM_COUPLING_LEVELSET_H_*/

/*!-----------------------------------------------------------------------------------------------*
\file xfem_coupling_levelset.H

\brief

<pre>
Maintainer: Benedikt Schott, Magnus Winter
            {schott,winter}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


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
      const double                        time,      ///< time
      const int                           step       ///< time step
  );


  //! destructor
  virtual ~LevelSetCoupling()
  {
    return;
  }


  /// set levelset field by function
  void SetLevelSetField(const double time);

  virtual void Output(
      const int step,
      const double time,
      const bool write_restart_data
  );

  Teuchos::RCP<Epetra_Vector> GetLevelSetField(){ return phinp_; }

  virtual void ReadRestart(
      const int step
  );

  virtual bool HasMovingInterface(){ return true; }

  virtual void GetInterfaceSlaveMaterial(
    DRT::Element* actele,
    Teuchos::RCP<MAT::Material> & mat
  )
  {
    mat = Teuchos::null;
  }

private:

  void SetConditionsToCopy();

protected:


  //! @name state vectors

  //! scatra-dis (cutterdis) state vectors for levelset applications, prepares nonmatching discretizations between scatra and fluid
  Teuchos::RCP<Epetra_Vector>    phinp_;

  //! The nodal curvature and smoothed gradient of the levelset field.
  Teuchos::RCP<Epetra_Vector>        curvaturenp_node_;
  Teuchos::RCP<Epetra_MultiVector>   gradphinp_smoothed_node_;
  Teuchos::RCP<Epetra_MultiVector>   gradphi2np_smoothed_node_;

  //! Output specific
  Teuchos::RCP<IO::DiscretizationWriter> ls_output_;

  //@}

};

class LevelSetCouplingBC : public LevelSetCoupling
{
public:

  //! constructor
  explicit LevelSetCouplingBC(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
      const std::string &                 cond_name,///< name of the condition, by which the derived cutter discretization is identified
      const double                        time,      ///< time
      const int                           step       ///< time step
  ) : LevelSetCoupling(bg_dis, cond_name, time, step)
  {}

  //! destructor
  virtual ~LevelSetCouplingBC()
  {
    return;
  }

protected:

  virtual void PrepareSolve();

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
      const double                        time,      ///< time
      const int                           step       ///< time step
  ) : LevelSetCouplingBC(bg_dis, cond_name, time, step)
  {}


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
      const double                        time,      ///< time
      const int                           step       ///< time step
  ) : LevelSetCouplingBC(bg_dis, cond_name, time, step)
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


class LevelSetCouplingTwoPhase : public LevelSetCoupling
{
public:

  //! constructor
  explicit LevelSetCouplingTwoPhase(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,    ///< background discretization
      const std::string &                 cond_name, ///< name of the condition, by which the derived cutter discretization is identified
      const double                        time,      ///< time
      const int                           step       ///< time step
  ) : LevelSetCoupling(bg_dis, cond_name, time, step),
      surfacetension_init_(false)
  {}

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
      INPAR::TWOPHASE::LaplaceBeltramiCalc  laplacebeltrami)
  {
    surftensapprox_ = surftensapprox;
    laplacebeltrami_ = laplacebeltrami;

    surfacetension_init_ = true;

  }

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
        for(size_t i = 0; i< nen; i++)
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
    for(unsigned int i =0; i<nsd; i++)
      for(unsigned int j =0; j<nsd; j++)
        if(i==j)
          eye(i,j)=1;

    //   Non-smoothed projection matrix
    LINALG::Matrix<nsd,nsd> p_matrix(true);
    if(laplacebeltrami_!=INPAR::TWOPHASE::matrix_smoothed)
    {
      for(unsigned int i =0; i<nsd; i++)
      {
        for(unsigned int j =0; j<nsd; j++)
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

      for(unsigned int i =0; i<nsd; i++)
      {
        for(unsigned int j =0; j<nsd; j++)
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
      double&   surf_coeff,             ///< surface tension coefficient gamma
      bool&     is_traction_jump        ///< is it a normal traction jump or done through laplace-beltrami
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

      is_traction_jump=false;

      // P = (\delta_{ij} - n_{i}n{j})
      // Calculations for evaluation with Laplace-Beltrami operator
      EvalProjectionMatrix<DISTYPE>(itraction_jump_matrix,eid, funct, derxy, normal);
      itraction_jump_matrix.Scale(surf_coeff);
    }

    return;
  }

  void ReadRestart(
      const int step
  );

private:
  INPAR::TWOPHASE::SurfaceTensionApprox surftensapprox_;
  INPAR::TWOPHASE::LaplaceBeltramiCalc  laplacebeltrami_;
  bool surfacetension_init_;

};

} // namespace XFEM

#endif /*XFEM_COUPLING_LEVELSET_H_*/

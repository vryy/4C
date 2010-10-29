/*!----------------------------------------------------------------------
\file strugenalpha_cmt.H

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STRUGENALPHA_CMT_H
#define STRUGENALPHA_CMT_H

#include "../drt_structure/strugenalpha.H"

// forward declarations
namespace MORTAR
{
  class ManagerBase;
}

namespace CONTACT
{

// forward declarations

/*!
\brief Generalized Alpha time integration for structural problems with meshtying or contact

\author popp (popp@lnm.mw.tum.de)
*/
class CmtStruGenAlpha : public virtual StruGenAlpha
{
public:

  /*!
  \brief Standard Constructor

  */
  CmtStruGenAlpha(ParameterList& params,
                  DRT::Discretization& dis,
                  LINALG::Solver& solver,
                  IO::DiscretizationWriter& output);

  /*!
  \brief Destructor

  */
  virtual ~CmtStruGenAlpha() {};

  /*!
  \brief Do time integration

  */
  void Integrate();

  /*!
  \brief Prepare contact for new time step

  */
  void PrepareStepContact();

  /*!
  \brief do constant predictor step

  */
  void ConstantPredictor();

  /*!
  \brief do consistent constant predictor step

  */
  void ConsistentPredictor();


  /// setup equilibrium with additional external forces
  void ApplyExternalForce(const STR::UTILS::MapExtractor& extractor,
                          Teuchos::RCP<Epetra_Vector> iforce);

  /// compute F (FSI)
  void computeF(const Epetra_Vector& x,
                Epetra_Vector& F);


  /// compute J at x (FSI)
  void computeJacobian(const Epetra_Vector& x);


  /// build linear system matrix and rhs (FSI)
  void Evaluate(Teuchos::RCP<const Epetra_Vector> disp);


  /// linear relaxation solve (FSI)
  Teuchos::RCP<Epetra_Vector> LinearRelaxationSolve(Teuchos::RCP<Epetra_Vector> relax);

  /*!
  \brief nonlinear solver for contact and meshtying
  
  This is a control routine which chooses the nonlinear solution method
  specified in the input file for the current contact or meshtying problem.
  It will either be called directly from Integrate() or from the Solve()
  method in the structural ADAPTER (adapter_structure_strugenalpha.cpp) in
  multiphysics problems (e.g. FSI simulations).

  */
  void CmtNonlinearSolve();
  
  /*!
  \brief do full (semi-smooth) Newton iteration

  This is a full Newton method for:
  
  (1) Penalty Meshtying
  (2) Lagrange Meshtying
  (3) Penalty Contact
  (4) Lagrange Contact (only for fixed-point PDASS, see below)
  
  Note for case (4): Here, it is only a standard full Newton if the PDASS
  (primal-dual active set strategy) is interpreted as a fixed-point
  method: The search for the correct active set (=contact nonlinearity)
  is represented by a fixed-point approach, whereas the large deformation
  linearization (=geometrical nonlinearity) is treated by a standard
  Newton scheme. This yields TWO nested iteration loops! If the PDASS
  is interpreted as a semi-smooth Newton method, this routine automatically
  yields another nonlinear solution method, namely semi-smooth Newton.
  
  This is a semi-smooth Newton method for:
  
  (1) Lagrange Contact (only for semi-smooth PDASS, see above)
  
  Note for case (1): Here, it is only a semi-smooth Newton if the PDASS
  (primal-dual active set strategy) is interpreted as a semi-smooth Newton
  method: The search for the correct active set (=contact nonlinearity)
  and the large deformation linearization (=geometrical nonlinearity) are
  merged into one semi-smooth Newton method and can therefore besolved
  within ONE iteration loop! If the PDASS is interpreted as a fixed-point
  method, this routine automatically yields another nonlinear solution
  method, namely standard full Newton.

  */
  void FullNewton();
    
  /*!
  \brief Line search version of full (semi-smooth) Newton iteration

  */
  void FullNewtonLineSearch();
  
  /*!
  \brief do pseudo transient continuation nonlinear iteration

  Pseudo transient continuation is a variant of full newton that has a
  larger convergence radius than newton and is therefore more stable and/or
  can do larger time steps. In this context, it is applicable wherever
  FullNewton() is applicable, thus in all cases except when treating
  Lagrange Contact with a semi-smooth Newton based PDASS.

  */
  void PTC();

  /*!
  \brief Check convergence of nonlinear iteration

  */
  bool Converged(const string type, const double disinorm,
                 const double resnorm, const double toldisp,
                 const double tolres);
  
  /*!
  \brief linear solve for meshtying and contact
  
  In most cases, the final linear system can either be condensed to displacement
  degrees of freedom (dual Lagrange strategy) only or contains only displacement
  degrees of freedom anyway (Penalty / Augmented Lagrange strategy). Yet, if a
  standard Lagrange strategy is applied, we obtain a saddle point system to be
  solved for BOTH displacement and Lagrange multiplier degrees of freedom. Thus,
  in that case, this function first builds the saddle point system from the
  respective matrix and vector blocks, then solves the saddle point system and
  finally extracts the displacement increments and the Lagrange multiplier solution.  
  In all other cases this function only performs the standard solver call.
  
  \param numiter(in): current nonlinear iteration step
  \param wanted(in): residual tolerance (for adapting solver tolerance)
  \param worst(in): current fresmnorm (for adapting solver tolerance)
  
  */
  virtual void LinearSolve(int numiter, double wanted, double worst);
  
  /*!
  \brief update configuration and output to file/screen

  */
  virtual void UpdateandOutput();

  /// update configuration
  virtual void Update();
  
  /// update elements
  virtual void UpdateElement();

  /// output to file/screen
  virtual void Output();

  /*!
  \brief output of linear+angular momentum and energies

  */
  void OutputEnergyMomentum();

  /*!
  \brief output of solution error norms

  */
  void OutputErrorNorms();
  
  /*!
  \brief restart calculation from a given step

  */
  void ReadRestart(int step);

  /*!
  \brief get manager

  */
  RCP<MORTAR::ManagerBase> GetManager() {return cmtmanager_; }
  
  /*!
  \brief return bool indicating if contact or meshtying are defined

  */
  bool HaveContactMeshtying() { return true; }
    
protected:

  // don't want = operator and cctor
  CmtStruGenAlpha operator = (const CmtStruGenAlpha& old);
  CmtStruGenAlpha(const CmtStruGenAlpha& old);

  RCP<MORTAR::ManagerBase>    cmtmanager_;

}; // class CmtStruGenAlpha
} // namespace CONTACT

#endif  // #ifndef STRUGENALPHA_CMT_H
#endif  // #ifdef CCADISCRET

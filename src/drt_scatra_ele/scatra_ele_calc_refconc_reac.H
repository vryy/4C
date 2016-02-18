/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_refconc_reac.H

\brief main file containing routines for calculation of scatra element formulated in reference concentrations
  and with advanced reaction terms

 <pre>
   Maintainer: Moritz Thon
               thon@mhpc.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-10364
 </pre>
 *----------------------------------------------------------------------*/


#ifndef SCATRA_ELE_CALC_REFCONC_REAC_H_
#define SCATRA_ELE_CALC_REFCONC_REAC_H_

#include "scatra_ele_calc_advanced_reaction.H"


namespace DRT
{
namespace ELEMENTS
{

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcRefConcReac : public ScaTraEleCalcAdvReac<distype>
{
private:
  /// private constructor, since we are a Singleton.
  ScaTraEleCalcRefConcReac(const int numdofpernode,const int numscal,const std::string& disname);

  typedef ScaTraEleCalc<distype> my;
  typedef ScaTraEleCalcAdvReac<distype> advreac;

public:
  virtual ~ScaTraEleCalcRefConcReac() {}

  /// Singleton access method
  static ScaTraEleCalcRefConcReac<distype> * Instance(
    const int numdofpernode,
    const int numscal,
    const std::string& disname,
    const ScaTraEleCalcRefConcReac* delete_me = NULL );

  /// called upon destruction
  virtual void Done();

protected:

  virtual double CalcReaCoeffFac(
      const std::vector<int>                    stoich,                  //!<stoichometrie of current condition
      const MAT::PAR::reaction_coupling         couplingtype,            //!<type of coupling the stoichiometry coefficients
      const std::vector<double>                 couprole,                //!<type of coupling role
      const double                              reacstart,               //!<reaction start coefficient
      const int                                 k,                       //!< id of current scalar
      const double                              scale = 1.0              //!< scale factor
  );

  virtual double CalcReaCoeffDerivFac(
      const std::vector<int>                  stoich,                  //!<stoichometrie of current condition
      const MAT::PAR::reaction_coupling       couplingtype,            //!<type of coupling the stoichiometry coefficients
      const std::vector<double>               couprole,                //!<type of coupling role
      const double                            reacstart,               //!<reaction start coefficient
      const int                               toderive,                //!<concentration to be derived to
      const int                               k,                       //!< id of current scalar
      const double                            scale = 1.0              //!< scale factor
  );

  virtual double CalcReaBodyForceTermFac(
      const std::vector<int>                      stoich,                 //!<stoichometrie of current condition
      const MAT::PAR::reaction_coupling           couplingtype,           //!<type of coupling the stoichiometry coefficients
      const std::vector<double>                   couprole,               //!<type of coupling role
      const double                                reacstart,              //!<reaction start coefficient
      const int                                   k,                      //!< id of current scalar
      const double                                scale = 1.0             //!< scale factor
  );

  virtual double CalcReaBodyForceDerivFac(
      const std::vector<int>                    stoich,                  //!<stoichometrie of current condition
      const MAT::PAR::reaction_coupling         couplingtype,            //!<type of coupling the stoichiometry coefficients
      const std::vector<double>                 couprole,                //!<type of coupling role
      const double                              reacstart,               //!<reaction start coefficient
      const int                                 toderive,                //!<concentration to be derived to
      const int                                 k,                       //!< id of current scalar
      const double                              scale = 1.0              //!< scale factor
  );

  //! calculation of convective element matrix: add conservative contributions
  virtual void CalcMatConvAddCons(
    Epetra_SerialDenseMatrix&     emat,        //!< element matrix to be filled
    const int                     k,           //!< index of current scalar
    const double                  timefacfac,  //!< domain-integration factor times time-integration factor
    const double                  vdiv,        //!< velocity divergence
    const double                  densnp       //!< density at time_(n+1)
  );

  //! set internal variables
  virtual void SetInternalVariablesForMatAndRHS();

  //! calculation of diffusive element matrix
  virtual void CalcMatDiff(
    Epetra_SerialDenseMatrix&     emat,        //!< element matrix to be filled
    const int                     k,           //!< index of current scalar
    const double                  timefacfac   //!< domain-integration factor times time-integration factor
  );

  //! calculate the Laplacian (weak form)
  void GetLaplacianWeakForm(
    double&                                  val,        //!< ?
    const LINALG::Matrix<my::nsd_,my::nsd_>& difftensor, //!< ?
    const int                                vi,         //!< ?
    const int                                ui          //!< ?
  )
  {
    val = 0.0;
    for (int j = 0; j<my::nsd_; j++)
    {
      for (int i = 0; i<my::nsd_; i++)
      {
        val += my::derxy_(j, vi)*difftensor(j,i)*my::derxy_(i, ui);
      }
    }
    return;
  };

  //! standard Galerkin diffusive term on right hand side
  virtual void CalcRHSDiff(
    Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
    const int                     k,           //!< index of current scalar
    const double                  rhsfac       //!< time-integration factor for rhs times domain-integration factor
  );

  //! calculate the Laplacian (weak form)
  void GetLaplacianWeakFormRHS(
    double&                                  val,        //!< ?
    const LINALG::Matrix<my::nsd_,my::nsd_>& difftensor, //!< ?
    const LINALG::Matrix<my::nsd_,1>&        gradphi,    //!< ?
    const int                                vi         //!< ?
  )
  {
    val = 0.0;
    for (int j = 0; j<my::nsd_; j++)
    {
      for (int i = 0; i<my::nsd_; i++)
      {
        val += my::derxy_(j,vi)*difftensor(j,i)*gradphi(i);
      }
    }
    return;
  };

  // add nodal displacements to point coordinates
  virtual void UpdateNodeCoordinates(){/*nothing to to since we want reference coordinates*/return;};

private:

  ///determinante of deformation gradient
  double J_;

  ///inverse of cauchy-green deformation gradient
  LINALG::Matrix<my::nsd_,my::nsd_> C_inv_;

  ///derivative dJ/dX by finite differences
  LINALG::Matrix<my::nsd_,1> dJdX_;

};  //end ScaTraEleCalcRefConcReac


}

}


#endif /* SCATRA_ELE_CALC_REFCONC_REAC_H_ */

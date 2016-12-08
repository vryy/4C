/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_multiporo_reac.H

 \brief evaluation class containing routines for calculation of scalar transport
        within multiphase porous medium

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_SCATRA_ELE_SCATRA_ELE_CALC_MULTIPORO_REAC_H_
#define SRC_DRT_SCATRA_ELE_SCATRA_ELE_CALC_MULTIPORO_REAC_H_


#include "scatra_ele_calc_poro_reac.H"

//forward declarations
namespace MAT
{
class ScatraMat;
}

namespace DRT
{
namespace ELEMENTS
{

template<int NSD,int NEN> class ScaTraEleInternalVariableManagerMultiPoro;

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcMultiPoroReac : public ScaTraEleCalcPoroReac<distype>
{
private:
  /// private constructor, since we are a Singleton.
  ScaTraEleCalcMultiPoroReac(const int numdofpernode,const int numscal,const std::string& disname);

  typedef ScaTraEleCalc<distype> my;
  typedef ScaTraEleCalcPoroReac<distype> pororeac;
  typedef ScaTraEleCalcPoro<distype> poro;
  typedef ScaTraEleCalcAdvReac<distype> advreac;

public:
  virtual ~ScaTraEleCalcMultiPoroReac() {}

  /// Singleton access method
  static ScaTraEleCalcMultiPoroReac<distype> * Instance(
    const int numdofpernode,
    const int numscal,
    const std::string& disname,
    const ScaTraEleCalcMultiPoroReac* delete_me = NULL );

  /// called upon destruction
  virtual void Done();

  /// Setup element evaluation
  virtual int SetupCalc(
      DRT::Element*               ele,
      DRT::Discretization&        discretization
      );

  protected:

  //! extract element based or nodal values
  //  return extracted values of phinp
  virtual void ExtractElementAndNodeValues(
      DRT::Element*                 ele,
      Teuchos::ParameterList&       params,
      DRT::Discretization&          discretization,
      DRT::Element::LocationArray&  la
  );

  //! set internal variables
  virtual void SetInternalVariablesForMatAndRHS();

  //! evaluate material
  virtual void Materials(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

  //! material MatMultiPoro
  virtual void MatMultiPoro(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

  //! Set advanced reaction terms and derivatives
  virtual void SetAdvancedReactionTerms(
      const int                                 k,           //!< index of current scalar
      const Teuchos::RCP<MAT::MatListReactions> matreaclist, //!< index of current scalar
      const double* gpcoord                                  //!< current Gauss-point coordinates
  );

  //! compute pore pressure
  virtual double ComputePorePressure(
    );

  //! get internal variable manager for multiporo formulation
  Teuchos::RCP<ScaTraEleInternalVariableManagerMultiPoro<my::nsd_,my::nen_> > VarManager()
      {return Teuchos::rcp_static_cast<ScaTraEleInternalVariableManagerMultiPoro<my::nsd_,my::nen_> >(my::scatravarmanager_);};

  //! calculation of convective element matrix in convective form
  //! the only difference to the base class version is, that there is no scaling with the density
  virtual void CalcMatConv(
    Epetra_SerialDenseMatrix&         emat,        //!< element matrix to be filled
    const int                         k,           //!< index of current scalar
    const double                      timefacfac,  //!< domain-integration factor times time-integration factor
    const double                      densnp,      //!< density at time_(n+1)
    const LINALG::Matrix<my::nen_,1>& sgconv       //!< subgrid-scale convective operator
  );

  //! adaption of convective term for rhs
  //! the only difference to the base class version is, that there is no scaling with the density
  virtual void RecomputeConvPhiForRhs(
    const int                          k,          //!< index of current scalar
    const LINALG::Matrix<my::nsd_,1>&  sgvelint,  //!< subgrid-scale velocity at Gauss point
    const double                       densnp,     //!< density at time_(n+1)
    const double                       densn,      //!< density at time_(n)
    const double                       vdiv       //!< velocity divergence
  );

  //! calculation of convective element matrix: add conservative contributions
  //! the only difference to the base class version is, that there is no scaling with the density
  virtual void CalcMatConvAddCons(
    Epetra_SerialDenseMatrix&     emat,        //!< element matrix to be filled
    const int                     k,           //!< index of current scalar
    const double                  timefacfac,  //!< domain-integration factor times time-integration factor
    const double                  vdiv,        //!< velocity divergence
    const double                  densnp       //!< density at time_(n+1)
  );

  void FillCouplingVector();

  protected:
  //! nodal flux values at t_(n+1)
  std::vector<LINALG::Matrix<my::nsd_,my::nen_> > efluxnp_;
  //! nodal phase pressures at t_(n+1) or t_(n+alpha_F)
  std::vector<LINALG::Matrix<my::nen_,1> > epresnp_;
  //! nodal phase saturation at t_(n+1) or t_(n+alpha_F)
  std::vector<LINALG::Matrix<my::nen_,1> > esatnp_;
  //! nodal solid pressure at t_(n+1) or t_(n+alpha_F)
  LINALG::Matrix<my::nen_,1> esolidpresnp_;

  //! a vector containing all quantities, the equation is coupled with
  //! (i.e. pressures, saturations and porosity)
  std::vector<std::pair<std::string,double> > couplingvalues_;

};


template<int NSD, int NEN>
class ScaTraEleInternalVariableManagerMultiPoro : public ScaTraEleInternalVariableManager<NSD, NEN>
{
  typedef ScaTraEleInternalVariableManager<NSD, NEN> my;

public:

  ScaTraEleInternalVariableManagerMultiPoro(int numscal):
    ScaTraEleInternalVariableManager<NSD, NEN>(numscal),
    multiconvelint_(0),
    multiconv_(0),
    pressure_(0),
    saturation_(0),
    solidpressure_(0.0),
    scalartophaseID_(numscal,-1)
  {return;}

  virtual ~ScaTraEleInternalVariableManagerMultiPoro(){};

  // compute and set internal variables
   void SetInternalVariablesMultiPoro(
       const LINALG::Matrix<NEN,1>&                 funct,       //! array for shape functions
       const LINALG::Matrix<NSD,NEN>&               derxy,       //! global derivatives of shape functions w.r.t x,y,z
       const std::vector<LINALG::Matrix<NEN,1> >&   ephinp,      //! scalar at t_(n+1) or t_(n+alpha_F)
       const std::vector<LINALG::Matrix<NEN,1> >&   ephin,       //! scalar at t_(n)
       const std::vector<LINALG::Matrix<NSD,NEN> >& econvelnp,   //! nodal convective velocity values at t_(n+1) or t_(n+alpha_F)
       const std::vector<LINALG::Matrix<NEN,1> >&   epresnp,
       const std::vector<LINALG::Matrix<NEN,1> >&   esatnp,
       const LINALG::Matrix<NEN,1>&                 esolidpresnp,
       const std::vector<LINALG::Matrix<NEN,1> >&   ehist        //! history vector of transported scalars
       )
   {
     // call base class (scatra)
     my::SetInternalVariables(funct,derxy,ephinp,ephin,econvelnp[0],ehist);

     const int numphases = econvelnp.size();

     // resize all phase related vectors
     multiconvelint_.resize(numphases);
     multiconv_.resize(numphases);
     pressure_.resize(numphases);
     saturation_.resize(numphases);

     for (int k = 0; k < numphases; ++k)
     {
       // convective velocity
       multiconvelint_[k].Multiply(econvelnp[k],funct);
       // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
       multiconv_[k].MultiplyTN(derxy,multiconvelint_[k]);

       // phase pressure
       pressure_[k] = epresnp[k].Dot(funct);
       // phase saturation
       saturation_[k] = esatnp[k].Dot(funct);
     }

     // solid pressure
     solidpressure_ = esolidpresnp.Dot(funct);

     // check if the IDs make sense
     for (std::vector<int>::iterator it = scalartophaseID_.begin() ; it != scalartophaseID_.end(); ++it)
     {
        if( *it<0 or *it>=numphases )
          dserror("Invalid phase ID %i",*it);
     }

     // overwrite convective term
     for (int k = 0; k < my::numscal_; ++k)
     {
       my::conv_phi_[k] = multiconvelint_[scalartophaseID_[k]].Dot(my::gradphi_[k]);
     }

     return;
   };

  /*========================================================================*/
  //! @name return methods for internal variables
  /*========================================================================*/

  //! return convective velocity
  virtual const LINALG::Matrix<NSD,1>& ConVel(const int k) const
  {
    return multiconvelint_[scalartophaseID_[k]];
  };

  //! return convective part in convective form
  virtual const LINALG::Matrix<NEN,1>& Conv(const int k) const
  {
    return multiconv_[scalartophaseID_[k]];
  };

  //! return pressure associated with scalar k
  double Pressure(const int k) const
  {
    return pressure_[scalartophaseID_[k]];
  };

  //! return pressures
  const std::vector<double>& Pressure() const
  {
    return pressure_;
  };

  //! return saturation associated with scalar k
  double Saturation(const int k) const
  {
    return saturation_[scalartophaseID_[k]];
  };

  //! return saturation associated with scalar k
  const std::vector<double>& Saturation() const
  {
    return saturation_;
  };

  //! return solid pressure
  double SolidPressure() const
  {
    return solidpressure_;
  };

  //! set scalar ID to phase ID mapping
  void SetPhaseID(const int scalarID, const int phaseID)
  {
    scalartophaseID_[scalarID] = phaseID;
  };

private:

  //! convective velocity
  std::vector<LINALG::Matrix<NSD,1> > multiconvelint_;
  //! convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
  std::vector<LINALG::Matrix<NEN,1> > multiconv_;

  //! phase pressure
  std::vector<double> pressure_;
  //! phase saturation
  std::vector<double> saturation_;
  //! solid pressure
  double solidpressure_;

  //! mapping from scalar ID to phase ID
  std::vector<int > scalartophaseID_;

};

}
}

#endif /* SRC_DRT_SCATRA_ELE_SCATRA_ELE_CALC_MULTIPORO_REAC_H_ */

/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_elch.H

\brief evalution of ScaTra elements for ion-transport equation

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef SCATRA_ELE_CALC_ELCH_H
#define SCATRA_ELE_CALC_ELCH_H

#include "scatra_ele_calc.H"

#include "../drt_inpar/inpar_elch.H"

namespace DRT
{
namespace ELEMENTS
{
  // forward declarations
  class ScaTraEleParameterElch;
  class ScaTraEleDiffManagerElch;

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcElch : public ScaTraEleCalc<distype>
{
  /// private constructor, since we are a Singleton.
  ScaTraEleCalcElch(const int numdofpernode,const int numscal);

  typedef ScaTraEleCalc<distype> my;

public:
  virtual ~ScaTraEleCalcElch() {}

  /// Singleton access method
  static ScaTraEleCalcElch<distype> * Instance(
    const int numdofpernode,
    const int numscal,
    bool create=true );

  /// called upon destruction
  virtual void Done();

  //! evaluate the element
  int Evaluate(DRT::ELEMENTS::Transport*     ele,
               Teuchos::ParameterList&       params,
               DRT::Discretization &         discretization,
               const std::vector<int> &      lm,
               Epetra_SerialDenseMatrix&     elemat1_epetra,
               Epetra_SerialDenseMatrix&     elemat2_epetra,
               Epetra_SerialDenseVector&     elevec1_epetra,
               Epetra_SerialDenseVector&     elevec2_epetra,
               Epetra_SerialDenseVector&     elevec3_epetra);

  int EvaluateService(DRT::ELEMENTS::Transport*     ele,
                      Teuchos::ParameterList&       params,
                      DRT::Discretization &         discretization,
                      const std::vector<int> &      lm,
                      Epetra_SerialDenseMatrix&     elemat1_epetra,
                      Epetra_SerialDenseMatrix&     elemat2_epetra,
                      Epetra_SerialDenseVector&     elevec1_epetra,
                      Epetra_SerialDenseVector&     elevec2_epetra,
                      Epetra_SerialDenseVector&     elevec3_epetra);

private:

  /*========================================================================*/
  //! @name general framework
  /*========================================================================*/

  //! calculate matrix and rhs. Here the whole thing is hidden.
  void Sysmat(
    DRT::Element*                         ele,       //!< the element we are dealing with
    Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs,      //!< element rhs to calculate
    Epetra_SerialDenseVector&             subgrdiff  //!< subgrid-diff.-scaling vector
  );

  /*========================================================================*/
  //! @name material and related and related functions
  /*========================================================================*/

  //! get the material parameters
  virtual void GetMaterialParams(
    const DRT::Element* ele,       //!< the element we are dealing with
    double&             densn,     //!< density at t_(n)
    double&             densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
    double&             densam,    //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManager> diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>  reamanager,   //!< reaction manager
    double&             visc,       //!< fluid viscosity
    const int           iquad=-1    //!< id of current gauss point (default = -1)
  );


  //! evaluate material
  virtual void Materials(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

  //! material ScaTra
  virtual void MatNewman(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

  //! material ScaTra
  virtual void MatDiffCond(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

  //! material ScaTra
  virtual void MatPhase(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               iphase,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

  void MatIon(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
    double&                                 visc,     //!< fluid viscosity
    const int                               iquad     //!< id of current gauss point
  );

  /*========================================================================*/
  //! @name methods for evaluation of individual terms
  /*========================================================================*/


  /*========================================================================*/
  //! @name additional service routines
  /*========================================================================*/

  void CalculateFlux(
  LINALG::Matrix<3,my::nen_>&         flux,
  const DRT::Element*             ele,
  INPAR::ELCH::ElchType           elchtype,
  const INPAR::SCATRA::FluxType   fluxtype,
  const int                       k
  );

  void CalErrorComparedToAnalytSolution(
    const DRT::Element*                   ele,
    const enum INPAR::ELCH::ElchType      elchtype,
    Teuchos::ParameterList&               params,
    Epetra_SerialDenseVector&             errors
  );

  void CalculateConductivity(
    const DRT::Element*               ele,
    const enum INPAR::ELCH::ElchType  elchtype,
    Epetra_SerialDenseVector&         sigma
  );

  void CalculateElectricPotentialField(
    const DRT::Element*               ele,
    const enum INPAR::ELCH::ElchType  elchtype,
    Epetra_SerialDenseMatrix&         emat,
    Epetra_SerialDenseVector&         erhs
  );

  /*========================================================================*/
  //! @name scalar degrees of freedom and related
  /*========================================================================*/

  //! electric potential
  LINALG::Matrix<my::nen_,1> epotnp_;

  //! nodal current density values at t_(n+1) or t_(n+alpha_F)
  LINALG::Matrix<my::nsd_,my::nen_> ecurnp_;

  DRT::ELEMENTS::ScaTraEleParameterElch* elchpara_;

  Teuchos::RCP<ScaTraEleDiffManagerElch> dme_;

  INPAR::ELCH::DiffCondMat diffcondmat_;

  bool cursolvar_;

  enum INPAR::ELCH::EquPot equpot_;

  //! flag for migration velocity in the stabilization operator (ELCH)
  bool migrationstab_;
  //! flag for migration velocity in formula for stabilization parameter (ELCH)
  bool migrationintau_;
  //! flag for migration velocity in residual of strong form (ELCH)
  bool migrationinresidual_;

  //TODO: BACI_ELE_CLEANING

  //! transference numbers
  double transelim_;

  std::vector<double> diffuselimderiv_;
  double diffuselim_;

  std::vector<double> eps_;
  std::vector<double> tort_;
  std::vector<double> epstort_;
};


/// ScaTraEleDiffManagerElch implementation
/*!
  This class keeps all elch-specific transport parameter needed for the evaluation of an element. The
  ScaTraEleDiffManagerElch is derived from the standard ScaTraEleDiffManager.
  for
*/
class ScaTraEleDiffManagerElch : public ScaTraEleDiffManager
{
  public:

    ScaTraEleDiffManagerElch(int numscal):
    ScaTraEleDiffManager(numscal),
    valence_(numscal,0.0),
    derivdiffcoef_(numscal,std::vector<double>(numscal,0.0)),
    transnum_(numscal,0.0),
    derivtransnum_(numscal,std::vector<double>(numscal,0.0)),
    cond_(0.0),
    derivcond_(numscal,0.0),
    thermfac_(1.0),
    derivthermfac_(numscal,0.0)
    {return;}

    ~ScaTraEleDiffManagerElch(){};

    //! Set valence of the single ionic species
    void SetValence(const double valence,const int k)
    {
      valence_[k]=valence;
      return;
    };

    //! Access routine for valence of the single ionic species
    const double GetValence(const int k){return valence_[k];};

    //! Set derivative of diffusion coefficients with respect to concentrations
    void SetDerivIsoDiffCoef(const double derivdiffcoef,const int k,const int iscal)
    {
      (derivdiffcoef_[k])[iscal] = derivdiffcoef;
      return;
    };

    //! Access routine for derivative of diffusion coefficients with respect to concentrations
    const double GetDerivIsoDiffCoef(const int k,const int iscal){return (derivdiffcoef_[k])[iscal];};

    //! Set transference numbers with respect to single ionic species
    void SetTransNum(const double transnum,const int k)
    {
      transnum_[k] = transnum;
      return;
    };

    //! Access routine for transference numbers with respect to single ionic species
    const double GetTransNum(const int k){return transnum_[k];};

    //! Set derivative of transference numbers with respect to concentrations
    void SetDerivTransNum(const double derivtransnum,const int k,const int iscal)
    {
      (derivtransnum_[k])[iscal] = derivtransnum;
      return;
    };

    //! Access routine for derivative of transference numbers with respect to concentrations
    const double GetDerivTransNum(const int k,const int iscal){return (derivtransnum_[k])[iscal];};

    //! Set conductivity of the electrolyte solution
    void SetCond(const double cond)
    {
      cond_ = cond;
      return;
    };

    //! Access routine for conductivity of the electrolyte solution
    const double GetCond(){return cond_;};

    //! Set derivative of the conductivity with respect to concentrations
    void SetDerivCond(const double derivcond,const int k)
    {
      derivcond_[k] = derivcond;
      return;
    };

    //! Access routine for derivative of the conductivity with respect to concentrations
    const double GetDerivCond(const int k){return derivcond_[k];};

    //! Set thermodynamic factor for a specific electrolyte solution
    void SetThermFac(const double thermfac)
    {
      thermfac_ = thermfac;
      return;
    };

    //! Access routine for the thermodynamic factor for a specific electrolyte solution
    const double GetThermFac(){return thermfac_;};

    //! Set derivative of thermodynamic factor with respect to concentrations
    void SetDerivThermFac(const double derivthermfac,const int k)
    {
      derivthermfac_[k] = derivthermfac;
      return;
    };

    //! Access routine for derivative of thermodynamic factor with respect to concentrations
    const double GetDerivThermFac(const int k){return derivthermfac_[k];};

    //! Output of transport parameter (to screen)
    void OutputTransportParams(const int numscal)
    {
      for(int k=0;k<numscal;++k)
        std::cout << "valence " << k << ":   " << valence_[k] << std::endl;

      for(int k=0;k<numscal;++k)
        std::cout << "diffusion coefficient " << k << ":   " << diff_[k] << std::endl;

      for(int k=0;k<numscal;++k)
      {
        for(int iscal=0;iscal<numscal;++iscal)
        {
          std::cout << "derivation diffusion coefficient (" << k << "," << iscal << "):  " << (derivdiffcoef_[k])[iscal] << std::endl;
        }
      }
      std::cout << std::endl;

      for(int k=0;k<numscal;++k)
        std::cout << "transference number " << k << ":   " << transnum_[k] << std::endl;

      for(int k=0;k<numscal;++k)
      {
        for(int iscal=0;iscal<numscal;++iscal)
        {
          std::cout << "derivation transference number (" << k << "," << iscal << "):  " << (derivtransnum_[k])[iscal] << std::endl;
        }
      }
      std::cout << std::endl;

      std::cout << "conductivity:   " << cond_ << std::endl;

      for(int k=0;k<numscal;++k)
        std::cout << "derivation of conductivity  " << k << ":   " << derivcond_[k] << std::endl;
      std::cout << std::endl;

      std::cout << "thermodynamic factor:   " << thermfac_ << std::endl;

      for(int k=0;k<numscal;++k)
        std::cout << "derivation of thermodynamic factor " << k << ":   " << derivthermfac_[k] << std::endl;
      std::cout << std::endl;

      return;
    };

    // return
    unsigned int LengthVector(){return diff_.size();};

    void IncreaseLengthVector(const int k, const int numscal)
    {
      if(diff_.size()== (unsigned) numscal)
      {
        valence_.push_back(0.0);
        diff_.push_back(0.0);
      }

      return;
    };


  private:

    //! valence of the single ionic species
    std::vector<double> valence_;

    //! derivative of diffusion coefficients with respect to concentrations
    std::vector< std::vector<double> > derivdiffcoef_;

    //! transference numbers for single ionic species
    std::vector<double> transnum_;

    //! derivative of transference numbers with respect to concentrations
    std::vector< std::vector<double> > derivtransnum_;

    //! conductivity of the electrolyte solution
    double cond_;

    //! derivative of the conductivity with respect to concentrations
    std::vector<double> derivcond_;

    //! thermodynamic factor for a specific electrolyte solution
    //! transport parameter used only in the Newman model
    double thermfac_;

    //! derivative of thermodynamic factor with respect to concentrations
    std::vector<double> derivthermfac_;

};

}
}

#endif


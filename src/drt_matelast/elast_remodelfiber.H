/*----------------------------------------------------------------------*/
/*!
\file elast_remodelfiber.H
\brief


the input line should read
  MAT 1 ELAST_RemodelFiber NUMMAT 1 MATIDS 100 TDECAY 1.0 SIGMAPRE 1.0 COLMASSFRAC 0.5

<pre>
Maintainer: Fabian Bräu
            braeu@lnm.mw.tum.de
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_REMODELFIBER_H_
#define ELAST_REMODELFIBER_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
namespace ELASTIC
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_RemodelFiber NUMMAT 1 MATIDS 100 TDECAY 1.0 SIGMAPRE 1.0
    class RemodelFiber : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      RemodelFiber(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~RemodelFiber() { ; }

      /// length of material list
      const int nummat_;

      /// the list of material IDs
      const std::vector<int>* matids_;

      /// @name material parameters
      //@{

      /// decay time of Poisson (degradation) process
      const double tdecay_;

      /// Cauchy prestress of new mass which is deposited during G&R
      const double sigmapre_;

      /// time constant for collagen growth
      const double k_growth_;

      /// Mass fraction of each fiber family in constraint mixture
      std::vector<double>* cur_w_collagen_;

      //@}

      /// create material instance of matching type with my parameters
      ///
      /// Remark:                                                  AN 03/2014
      /// Function has to be implemented, since it is a pure virtual
      /// function in the MAT::PAR::Material class.
      /// Function stays empty, since the call within the material factory is done only
      /// respective the elasthyper class
      ///
      virtual Teuchos::RCP<MAT::Material> CreateMaterial(){ return Teuchos::null; };

    };  // class RemodelFiber

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// General interface for fibers which remodel

  class RemodelFiber : public Summand
  {

  public:

    /// constructor with given material parameters
    RemodelFiber(MAT::ELASTIC::PAR::RemodelFiber* params);

    ///@name Packing and Unpacking
    //@{

    virtual void PackSummand(DRT::PackBuffer& data) const;

    virtual void UnpackSummand(const std::vector<char>& data, std::vector<char>::size_type& position);

    //@}

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_remodelfiber; };

    /// Return number of fiber families
    virtual unsigned GetNumFibers() { return potsumfiber_.size(); };

    //@}

    /// Update of summand
    virtual void Update();

    /// Setup of summand
    virtual void Setup(int numgp,DRT::INPUT::LineDefinition* linedef);

    /// Retrieve stress and cmat
    virtual void AddStressCmatRemodel(
      const LINALG::Matrix<3,3>* defgrd,   ///< deformation gradient
      Teuchos::ParameterList& params,      ///< Container for additional information
      LINALG::Matrix<6,6>& cmat,           ///< material stiffness matrix
      LINALG::Matrix<6,1>& stress,         ///< 2nd PK-stress
      const int eleGID                     ///< element ID
      );

    /// Return names of visualization data
    virtual void VisNames(std::map<std::string,int>& names);

    /// Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleId);

    /// Indicator for the chosen formulations
    virtual void SpecifyFormulation(
        bool& isoprinc, ///< global indicator for isotropic principal formulation
        bool& isomod, ///< global indicator for isotropic splitted formulation
        bool& anisoprinc, ///< global indicator for anisotropic principal formulation
        bool& anisomod, ///< global indicator for anisotropic splitted formulation
        bool& viscogeneral       ///< global indicator, if one viscoelastic formulation is used
    )
    {
      dserror("RemodelFiber cannot be used with ElastHyper. Use GrowthRemodel_ElastHyper!");
      return;
    }

  private:

    /// my material parameters
    MAT::ELASTIC::PAR::RemodelFiber* params_;

    /// last inverse inelastic stretch in fiberdirection for each Gauss point at converged state
    std::vector<std::vector<double> > last_ilambda_r_;

    /// current inverse inelastic stretch in fiberdirection for each Gauss point at converged state
    std::vector<std::vector<double> > current_ilambda_r_;

    /// current mass fractions of each fiber family
    std::vector<std::vector<double> > current_w_collagen_;

    /// current fiber Cauchy stress
    std::vector<std::vector<double> > stress_;

    /// map to materials/potential summands (fibers) which remodel
    std::vector<Teuchos::RCP<MAT::ELASTIC::Summand> > potsumfiber_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // ELAST_REMODELFIBER_H_

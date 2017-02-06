/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_aniso.H

\brief scatra_ele_calc_aniso.H

\level 3

<pre>
  \maintainer Lasse Jagschies
              jagschies@mhpc.mw.tum.de
              http://www.lnm.mw.tum.de
              089 - 289-10365
</pre>
 *----------------------------------------------------------------------*/


#ifndef SCATRA_ELE_CALC_ANISO_H_
#define SCATRA_ELE_CALC_ANISO_H_

#include "scatra_ele_calc.H"

namespace DRT
{
  namespace ELEMENTS
  {
    template <int NSD>
    class ScaTraEleDiffManagerAniso;

    template<DRT::Element::DiscretizationType distype,int probdim=DRT::UTILS::DisTypeToDim<distype>::dim>
    class ScaTraEleCalcAniso : public virtual ScaTraEleCalc<distype,probdim>
    {
    protected:

      /// (private) protected constructor, since we are a Singleton.
      /// this constructor is called from a derived class
      /// -> therefore, it has to be protected instead of private
      ScaTraEleCalcAniso(const int numdofpernode,const int numscal,const std::string& disname);

    private:

      typedef ScaTraEleCalc<distype,probdim> my;
      typedef ScaTraEleInternalVariableManager<my::nsd_,my::nen_> varmanager;

    public:

      virtual ~ScaTraEleCalcAniso() {}


      /// Singleton access method
      static ScaTraEleCalcAniso<distype,probdim> * Instance(
        const int numdofpernode,
        const int numscal,
        const std::string& disname,
        const ScaTraEleCalcAniso* delete_me = NULL );

      /// called upon destruction
      virtual void Done();

    protected:

      /*========================================================================*/
      //! @name material and related and related functions
      /*========================================================================*/

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

      //! material ScaTra
      void MatScaTraAniso(
        const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
        const int                               k,        //!< id of current scalar
        double&                                 densn,    //!< density at t_(n)
        double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
        double&                                 densam,   //!< density at t_(n+alpha_M)
        double&                                 visc,      //!< fluid viscosity
        const int                               iquad=-1   //!< id of current gauss point (default = -1)
      );

      //! get diffusion manager for anisotropic diffusivity
      Teuchos::RCP<ScaTraEleDiffManagerAniso<my::nsd_> > DiffManager(){return Teuchos::rcp_static_cast<ScaTraEleDiffManagerAniso<my::nsd_> >(my::diffmanager_);};

      /*========================================================================*/
      //! @name methods for evaluation of individual terms
      /*========================================================================*/

      //! calculate the Laplacian for all shape functions(strong form)
      void GetLaplacianStrongForm(
        LINALG::Matrix<my::nen_,1>&        diff  //!< laplace term to be computed
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
      for (unsigned j = 0; j<my::nsd_; j++)
      {
        for (unsigned i = 0; i<my::nsd_; i++)
        {
          val += my::derxy_(j, vi)*difftensor(j,i)*my::derxy_(i, ui);
        }
      }
      return;
    };

    //! calculate the Laplacian (weak form)
    void GetLaplacianWeakFormRHS(
      double&                                  val,        //!< ?
      const LINALG::Matrix<my::nsd_,my::nsd_>& difftensor, //!< ?
      const LINALG::Matrix<my::nsd_,1>&        gradphi,    //!< ?
      const int                                vi         //!< ?
    )
    {
      val = 0.0;
      for (unsigned j = 0; j<my::nsd_; j++)
      {
        for (unsigned i = 0; i<my::nsd_; i++)
        {
          val += my::derxy_(j,vi)*difftensor(j,i)*gradphi(i);
        }
      }
      return;
    };

    //! calculation of diffusive element matrix
    virtual void CalcMatDiff(
      Epetra_SerialDenseMatrix&             emat,        //!< element matrix to be filled
      const int                             k,           //!< index of current scalar
      const double                          timefacfac   //!< domain-integration factor times time-integration factor
    );

    //! standard Galerkin diffusive term on right hand side
    virtual void CalcRHSDiff(
        Epetra_SerialDenseVector&           erhs,        //!< element vector to be filled
        const int                           k,           //!< index of current scalar
        const double                        rhsfac       //!< time-integration factor for rhs times domain-integration factor
    );

    };


    /// Scatra anisotropic diffusion manager
    /*!
    This is a class to handle anisotropic diffusion. It is derived from the
    basic isotropic diffusion manager class in scatra_ele_calc.H. It doesn't
    contain subgrid diffusion.
    */
    template<int NSD>
    class ScaTraEleDiffManagerAniso : public ScaTraEleDiffManager
    {
    public:

      ScaTraEleDiffManagerAniso(int numscal):
        ScaTraEleDiffManager(numscal),
        difftensor_(numscal)
      {return;}

      virtual ~ScaTraEleDiffManagerAniso(){};

      //! Set the anisotropic diffusion coefficient
      virtual void SetAnisotropicDiff(const LINALG::Matrix<NSD,NSD> difftensor, const int k)
      {
        difftensor_[k] = difftensor;
        return;
      }

      //! Return the stored anisotropic diffusion coefficient
      virtual LINALG::Matrix<NSD,NSD> GetAnisotropicDiff(const int k){return difftensor_[k];}

    protected:

      //! tensor valued diffusion coefficient
      std::vector<LINALG::Matrix<NSD,NSD> > difftensor_;

    };

  }
}


#endif /* SCATRA_ELE_CALC_ANISO_H_ */

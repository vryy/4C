/*----------------------------------------------------------------------*/
/*!

 \brief evaluation class containing routines for calculation of scalar transport
        within porous medium

\level 2

\maintainer Anh-Tu Vuong
 *----------------------------------------------------------------------*/


#ifndef SCATRA_ELE_CALC_PORO_H_
#define SCATRA_ELE_CALC_PORO_H_

#include "scatra_ele_calc.H"

// forward declarations
namespace MAT
{
  class ScatraMat;
}

namespace DRT
{
  namespace ELEMENTS
  {
    class ScaTraEleDiffManagerPoro;


    template <DRT::Element::DiscretizationType distype>
    class ScaTraEleCalcPoro : public virtual ScaTraEleCalc<distype>
    {
     protected:
      /// (private) protected constructor, since we are a Singleton.
      ScaTraEleCalcPoro(const int numdofpernode, const int numscal, const std::string& disname);

     private:
      typedef ScaTraEleCalc<distype> my;

     public:
      virtual ~ScaTraEleCalcPoro() {}

      /// Singleton access method
      static ScaTraEleCalcPoro<distype>* Instance(const int numdofpernode, const int numscal,
          const std::string& disname, const ScaTraEleCalcPoro* delete_me = NULL);

      /// called upon destruction
      virtual void Done();

      /// Evaluate the element
      /*!
        Generic virtual interface function. Called via base pointer.
       */
      //   virtual int Evaluate(DRT::Element*                 ele,
      //                        Teuchos::ParameterList&       params,
      //                        DRT::Discretization &         discretization,
      //                        const std::vector<int> &      lm,
      //                        Epetra_SerialDenseMatrix&     elemat1_epetra,
      //                        Epetra_SerialDenseMatrix&     elemat2_epetra,
      //                        Epetra_SerialDenseVector&     elevec1_epetra,
      //                        Epetra_SerialDenseVector&     elevec2_epetra,
      //                        Epetra_SerialDenseVector&     elevec3_epetra);

     protected:
      //! evaluate action
      virtual int EvaluateAction(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const SCATRA::Action& action,
          DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra);

      //   int EvaluateODMesh(  DRT::Element*                 ele,
      //                        Teuchos::ParameterList&       params,
      //                        DRT::Discretization &         discretization,
      //                        const std::vector<int> &      lm,
      //                        Epetra_SerialDenseMatrix&     elemat1_epetra,
      //                        Epetra_SerialDenseMatrix&     elemat2_epetra,
      //                        Epetra_SerialDenseVector&     elevec1_epetra,
      //                        Epetra_SerialDenseVector&     elevec2_epetra,
      //                        Epetra_SerialDenseVector&     elevec3_epetra);
      //
      //   int EvaluateODFluid(  DRT::Element*                 ele,
      //                        Teuchos::ParameterList&       params,
      //                        DRT::Discretization &         discretization,
      //                        const std::vector<int> &      lm,
      //                        Epetra_SerialDenseMatrix&     elemat1_epetra,
      //                        Epetra_SerialDenseMatrix&     elemat2_epetra,
      //                        Epetra_SerialDenseVector&     elevec1_epetra,
      //                        Epetra_SerialDenseVector&     elevec2_epetra,
      //                        Epetra_SerialDenseVector&     elevec3_epetra);
      //
      //   //! calculate matrix and rhs. Here the whole thing is hidden.
      //   virtual void SysmatODMesh(
      //     DRT::Element*                         ele,       //!< the element we are dealing with
      //     Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
      //     const int                     numdofpernode
      //   );
      //
      //   //! calculate matrix and rhs. Here the whole thing is hidden.
      //   virtual void SysmatODFluid(
      //     DRT::Element*                         ele,       //!< the element we are dealing with
      //     Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
      //     const int                     numdofpernode
      //   );

      //! read element coordinates
      virtual void ReadElementCoordinates(const DRT::Element* ele);

      //! extract element based or nodal values
      //  return extracted values of phinp
      virtual void ExtractElementAndNodeValues(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la);

      //! extract element based or nodal values
      //  return extracted values of phinp
      virtual void ExtractElementAndNodeValuesPoro(DRT::Element* ele,
          Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Element::LocationArray& la);

      //! get the material parameters
      virtual void GetMaterialParams(const DRT::Element* ele,  //!< the element we are dealing with
          std::vector<double>& densn,                          //!< density at t_(n)
          std::vector<double>& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
          std::vector<double>& densam,  //!< density at t_(n+alpha_M)
          double& visc,                 //!< fluid viscosity
          const int iquad = -1          //!< id of current gauss point (default = -1)
      );

      //! compute porosity based on solid, fluid and (potentially) scatra solution
      virtual void ComputePorosity(const DRT::Element* ele  //!< the element we are dealing with
      );

      //! compute pore pressure
      virtual double ComputePorePressure();

      //! material ScaTra
      virtual void MatScaTra(
          const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
          const int k,                                       //!< id of current scalar
          double& densn,                                     //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
      );


      //! set diffusivity for poro scatra problem (i.e. scale by porosity)
      virtual void SetDiffusivity(
          const Teuchos::RCP<const MAT::ScatraMat>& material, const int k, const double scale);

      //! set densisties for poro scatra problem (i.e. scale by porosity)
      virtual void SetDensities(double porosity,
          double& densn,   //!< density at t_(n)
          double& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam   //!< density at t_(n+alpha_M));
      );

      //! calculate scalar(s) and domain integral
      virtual void CalculateScalars(const DRT::Element* ele,  //!< the element we are dealing with
          Epetra_SerialDenseVector& scalars,                  //!< scalar to be computed
          const bool inverting                                //!< bool indicating inversion
      );

      //! get poro diffusion manager
      Teuchos::RCP<ScaTraEleDiffManagerPoro> DiffManager()
      {
        return Teuchos::rcp_static_cast<ScaTraEleDiffManagerPoro>(my::diffmanager_);
      };

      /*========================================================================*/
      //! @name Galerkin approximation and related
      /*========================================================================*/

      //! initial node coordinates
      LINALG::Matrix<my::nsd_, my::nen_> xyze0_;

      //! nodal porosity values at t_(n+1)
      LINALG::Matrix<my::nen_, 1> eporosity_;

      //! flag indacting a node based porosity
      bool isnodalporosity_;
    };

    /// ScaTraEleDiffManagerPoro implementation
    /*!
      This class keeps all poro-specific transport parameter needed for the evaluation of an
      element. The ScaTraEleDiffManagerPoro is derived from the standard ScaTraEleDiffManager.
    */

    ////TODO: HACK
    // const int NO_CONVECTION_NR = 6;

    class ScaTraEleDiffManagerPoro : public ScaTraEleDiffManager
    {
     public:
      ScaTraEleDiffManagerPoro(int numscal) : ScaTraEleDiffManager(numscal), porosity_(0.0)
      {
        return;
      }

      ~ScaTraEleDiffManagerPoro(){};

      void SetPorosity(double porosity) { porosity_ = porosity; }

      double GetPorosity(const int k) const
      {
        //      if(k<NO_CONVECTION_NR)
        return porosity_;
        //      else
        //        return 1.0;
      }

     protected:
      //! porosity at gauss point
      double porosity_;
    };

    template <int NSD, int NEN>
    class ScaTraEleInternalVariableManagerPoro : public ScaTraEleInternalVariableManager<NSD, NEN>
    {
      typedef ScaTraEleInternalVariableManager<NSD, NEN> my;

     public:
      ScaTraEleInternalVariableManagerPoro(int numscal)
          : ScaTraEleInternalVariableManager<NSD, NEN>(numscal),
            zeroconvel_(true),
            zeroconv_(true),
            zero_(0.0)
      {
        return;
      }

      virtual ~ScaTraEleInternalVariableManagerPoro(){};

      /*========================================================================*/
      //! @name return methods for internal variables
      /*========================================================================*/

      //! return convective velocity
      virtual const LINALG::Matrix<NSD, 1>& ConVel(const int k) const {//    if(k<NO_CONVECTION_NR)
          {return my::convelint_;
    }
    //    else
    //      return zeroconvel_;
  };  // namespace ELEMENTS
  //! return convective part in convective form
  virtual const LINALG::Matrix<NEN, 1>& Conv(const int k) const {//    if(k<NO_CONVECTION_NR)
      {return my::conv_;
}  // namespace DRT
   //    else
   //      return zeroconv_;
}
;

//! return convective term of current scalar value
virtual const double& ConvPhi(const int k) const {//    if(k<NO_CONVECTION_NR)
    {return my::conv_phi_[k];
}
//    else
//      return zero_;
}
;

private:
LINALG::Matrix<NSD, 1> zeroconvel_;
LINALG::Matrix<NEN, 1> zeroconv_;
double zero_;
}
;
}
}


#endif /* SCATRA_ELE_CALC_PORO_H_ */

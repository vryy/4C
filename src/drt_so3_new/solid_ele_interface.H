/*----------------------------------------------------------------------*/
/*! \file

\brief Interface of solid elements

\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef SOLID_ELE_INTERFACE_H
#define SOLID_ELE_INTERFACE_H

#include "Teuchos_ParameterList.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "../drt_lib/drt_dserror.H"

namespace STR::ELEMENTS
{
  class ParamsInterface;
}  // namespace STR::ELEMENTS

namespace DRT
{
  // forward declaration
  class Discretization;
  class Condition;
  namespace UTILS
  {
    class GaussIntegration;
  }
  namespace INPUT
  {
    class LineDefinition;
  }

  namespace ELEMENTS
  {
    // forward declaration
    class Solid;

    class SolidEleInterface
    {
     public:
      /*!
      \brief Evaluate an element

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param elemat1 (out)  : matrix to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elemat2 (out)  : matrix to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elevec1 (out)  : vector to be filled by element. If NULL on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If NULL on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If NULL on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      virtual int Evaluate(DRT::ELEMENTS::Solid* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
          Epetra_SerialDenseVector& elevec3) = 0;

      /// Evaluate the element
      /*!
        This class does not provide a definition for this function; it
        must be defined in SolidEleCalc.
       */
      virtual int nln_force_stiff_mass(DRT::ELEMENTS::Solid* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params, Epetra_SerialDenseMatrix* elemat1_epetra,
          Epetra_SerialDenseMatrix* elemat2_epetra, Epetra_SerialDenseVector* elevec1_epetra,
          Epetra_SerialDenseVector* elevec2_epetra, Epetra_SerialDenseVector* elevec3_epetra,
          DRT::UTILS::GaussIntegration* intpoints = nullptr) = 0;

      /// Evaluate the element for GEMM time integration
      /*!
        This class does not provide a definition for this function; it
        must be defined in SolidEleCalc.
       */
      virtual int nln_force_stiff_mass_gemm(DRT::ELEMENTS::Solid* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params, Epetra_SerialDenseMatrix* elemat1_epetra,
          Epetra_SerialDenseMatrix* elemat2_epetra, Epetra_SerialDenseVector* elevec1_epetra,
          Epetra_SerialDenseVector* elevec2_epetra, Epetra_SerialDenseVector* elevec3_epetra,
          DRT::UTILS::GaussIntegration* intpoints = nullptr)
      {
        dserror("Not implemented for this element evaluation");
        return 0;
      }

      /// Evaluate Neumann condition on this element
      virtual int EvaluateNeumann(DRT::ELEMENTS::Solid* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          const DRT::UTILS::GaussIntegration* intpoints, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseMatrix* elemat1 = nullptr) = 0;

      /// calculate stress for output
      virtual int CalcStress(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra) = 0;

      /// Update the element at converged time step
      virtual int UpdateElement(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra) = 0;

      virtual int ResetStep(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra) = 0;

      virtual int StoreStep(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra) = 0;

      virtual int RecoverStep(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra) = 0;

      virtual int ResetIstep(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra) = 0;

      virtual int CalcEnergy(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra)
      {
        dserror("Not implemented for this element evaluation");
        return 0;
      }

      virtual int CalcErrorNorms(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra)
      {
        dserror("Not implemented for this element evaluation");
        return 0;
      }

      virtual int SetMultiEAS(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra)
      {
        dserror("Not implemented for this element evaluation");
        return 0;
      }

      virtual int ReadMultiRestart(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra)
      {
        dserror("Not implemented for this element evaluation");
        return 0;
      }

      virtual int CalcMultiDens(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra)
      {
        dserror("Not implemented for this element evaluation");
        return 0;
      }

      virtual int PostProcessStress(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra)
      {
        dserror("Not implemented for this element evaluation");
        return 0;
      }

      virtual int UpdateInverseDesign(DRT::ELEMENTS::Solid* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params, Epetra_SerialDenseMatrix* elemat1_epetra,
          Epetra_SerialDenseMatrix* elemat2_epetra, Epetra_SerialDenseVector* elevec1_epetra,
          Epetra_SerialDenseVector* elevec2_epetra, Epetra_SerialDenseVector* elevec3_epetra)
      {
        dserror("Not implemented for this element evaluation");
        return 0;
      }

      virtual int UpdatePrestress(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra)
      {
        dserror("Not implemented for this element evaluation");
        return 0;
      }

      virtual int SwitchInverseDesign(DRT::ELEMENTS::Solid* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params, Epetra_SerialDenseMatrix* elemat1_epetra,
          Epetra_SerialDenseMatrix* elemat2_epetra, Epetra_SerialDenseVector* elevec1_epetra,
          Epetra_SerialDenseVector* elevec2_epetra, Epetra_SerialDenseVector* elevec3_epetra)
      {
        dserror("Not implemented for this element evaluation");
        return 0;
      }

      virtual int GlobalGPstressMap(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra) = 0;

      virtual int InterpolatePoroVelocityToPoint(DRT::ELEMENTS::Solid* ele,
          DRT::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params, Epetra_SerialDenseMatrix* elemat1_epetra,
          Epetra_SerialDenseMatrix* elemat2_epetra, Epetra_SerialDenseVector* elevec1_epetra,
          Epetra_SerialDenseVector* elevec2_epetra, Epetra_SerialDenseVector* elevec3_epetra)
      {
        dserror("Not implemented for this element evaluation");
        return 0;
      }

      virtual int CalcMassVolume(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra)
      {
        dserror("Not implemented for this element evaluation");
        return 0;
      }

      /// recover variables condensed at element level
      virtual int RecoverCondensed(DRT::ELEMENTS::Solid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Epetra_SerialDenseMatrix* elemat1_epetra, Epetra_SerialDenseMatrix* elemat2_epetra,
          Epetra_SerialDenseVector* elevec1_epetra, Epetra_SerialDenseVector* elevec2_epetra,
          Epetra_SerialDenseVector* elevec3_epetra) = 0;

      /// lump mass matrix
      virtual void LumpMass(Epetra_SerialDenseMatrix* mass) = 0;

      /// setup element data
      virtual void Setup(DRT::ELEMENTS::Solid* ele, DRT::INPUT::LineDefinition* linedef) = 0;

      virtual void MaterialPostSetup(const DRT::ELEMENTS::Solid& ele) = 0;

      /// get stuff from the STR::ELEMENTS::ParamsInterface and write it in the parameter list
      /// todo: for now, we need this, since the solid materials still rely on it.
      ///       in a next step they should start using the paramsInterface as well
      virtual void ParamsInterfaceToList(
          STR::ELEMENTS::ParamsInterface& pi, Teuchos::ParameterList& pl);

     protected:
      //! ctor
      SolidEleInterface() = default;

    };  // class SolidEleInterface
  }     // namespace ELEMENTS

}  // namespace DRT

#endif  // SOLID_ELE_INTERFACE_H

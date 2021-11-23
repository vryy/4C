/*----------------------------------------------------------------------------*/
/*! \file
\brief three dimensional total Lagrange truss element used for scalar transport coupling

\level 3

*/
/*---------------------------------------------------------------------------*/
#ifndef TRUSS3_SCATRA_H
#define TRUSS3_SCATRA_H

#include "truss3.H"

#include "../drt_inpar/inpar_scatra.H"

namespace DRT
{
  namespace ELEMENTS
  {
    class Truss3ScatraType : public Truss3Type
    {
     public:
      DRT::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<DRT::Element> Create(const std::string eletype, const std::string eledistype,
          const int id, const int owner) override;

      Teuchos::RCP<DRT::Element> Create(const int id, const int owner) override;

      static Truss3ScatraType& Instance();

      std::string Name() const override { return "Truss3ScatraType"; }

      void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
          override;

     private:
      static Truss3ScatraType instance_;
    };

    /*!
     \brief three dimensional total Lagrange truss element for scalar transport coupling

     */
    class Truss3Scatra : public Truss3
    {
     public:
      /*!
       \brief Standard Constructor

       \param id    (in): A globally unique element id
       \param owner (in): owner processor of the element
       */
      Truss3Scatra(int id, int owner);

      Truss3Scatra(const Truss3Scatra& old);

      DRT::Element* Clone() const override;

      DRT::ElementType& ElementType() const override { return Truss3ScatraType::Instance(); }

      /// return SCATRA::ImplType
      const INPAR::SCATRA::ImplType& ImplType() const { return impltype_; };

      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef) override;

      int UniqueParObjectId() const override
      {
        return Truss3ScatraType::Instance().UniqueParObjectId();
      }

      void CalcInternalForceStiffTotLag(
          Epetra_SerialDenseVector& forcevec, Epetra_SerialDenseMatrix& stiffmat) override;

      void Pack(DRT::PackBuffer& data) const override;
      void Unpack(const std::vector<char>& data) override;

     protected:
      void ExtractElementalVariables(LocationArray& la, DRT::Discretization& discretization,
          Teuchos::ParameterList& params) override;

      void t3_energy(Teuchos::ParameterList& params, Epetra_SerialDenseVector& intenergy) override;

     private:
      //! scalar transport implementation type (physics)
      INPAR::SCATRA::ImplType impltype_;

      //! evaluate elemental specific values
      //!
      //! \param[out] truss_disp       nodal displacement
      //! \param[out] dtruss_disp_du   derivative of truss displacement w.r.t. global displacement
      //! \param[out] dN_dx            derivative of shape functions
      //! \param[out] nodal_concentration  nodal concentrations
      void PrepCalcInternalForceStiffTotLagScaTra(LINALG::Matrix<6, 1>& truss_disp,
          LINALG::Matrix<6, 6>& dtruss_disp_du, LINALG::Matrix<6, 1>& dN_dx,
          LINALG::Matrix<2, 1>& nodal_concentration);

      //! calculation of concentration at Gauß Points, given concentration at nodes
      double ProjectScalarToGaussPoint(double xi, const LINALG::Matrix<2, 1>& c) const;

      //! compute polynomial growth factor
      double GetGrowthFactorConcProp(
          int poly_num, double c_GP, double c_0, const std::vector<double>& poly_params) const;

      //! compute polynomial growth factor based on amount of substance
      double GetGrowthFactorAoSProp(
          int poly_num, double c_GP, double c_0, const std::vector<double>& poly_params) const;

      //! nodal concentration for this element
      std::vector<double> phi_ele_;

      // don't want = operator
      Truss3Scatra& operator=(const Truss3Scatra& old);
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif  // #ifndef TRUSS3_SCATRA_H

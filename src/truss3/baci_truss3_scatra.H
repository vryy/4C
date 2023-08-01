/*----------------------------------------------------------------------------*/
/*! \file
\brief three dimensional total Lagrange truss element used for scalar transport coupling

\level 3

*/
/*---------------------------------------------------------------------------*/
#ifndef BACI_TRUSS3_SCATRA_H
#define BACI_TRUSS3_SCATRA_H

#include "baci_inpar_scatra.H"
#include "baci_truss3.H"

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

      void CalcInternalForceStiffTotLag(const std::map<std::string, std::vector<double>>& ele_state,
          CORE::LINALG::SerialDenseVector& forcevec,
          CORE::LINALG::SerialDenseMatrix& stiffmat) override;

      void CalcGPStresses(Teuchos::ParameterList& params,
          const std::map<std::string, std::vector<double>>& ele_state) override;

      void Pack(DRT::PackBuffer& data) const override;
      void Unpack(const std::vector<char>& data) override;

     protected:
      void ExtractElementalVariables(LocationArray& la, const DRT::Discretization& discretization,
          const Teuchos::ParameterList& params,
          std::map<std::string, std::vector<double>>& ele_state) override;

      void Energy(const std::map<std::string, std::vector<double>>& ele_state,
          Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector& intenergy) override;

     private:
      //! scalar transport implementation type (physics)
      INPAR::SCATRA::ImplType impltype_;

      //! evaluate elemental specific values
      //!
      //! @param[in] ele_state              elemental states (depending on the instantiated element)
      //! @param[out] curr_nodal_coords     nodal displacement
      //! @param[out] dcurr_nodal_coords_du derivative of truss displacement w.r.t. global
      //! displacement
      //! @param[out] dN_dx               derivative of shape functions
      //! @param[out] nodal_concentration nodal concentrations
      void PrepCalcInternalForceStiffTotLagScaTra(
          const std::map<std::string, std::vector<double>>& ele_state,
          CORE::LINALG::Matrix<6, 1>& curr_nodal_coords,
          CORE::LINALG::Matrix<6, 6>& dcurr_nodal_coords_du, CORE::LINALG::Matrix<6, 1>& dN_dx,
          CORE::LINALG::Matrix<2, 1>& nodal_concentration);

      //! calculation of concentration at Gauss Points, given concentration at nodes
      double ProjectScalarToGaussPoint(double xi, const CORE::LINALG::Matrix<2, 1>& c) const;

      // don't want = operator
      Truss3Scatra& operator=(const Truss3Scatra& old);
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif  // TRUSS3_SCATRA_H

/*----------------------------------------------------------------------*/
/*! \file

\brief a common base class for all solid elements

\level 2


 *----------------------------------------------------------------------*/


#ifndef FOUR_C_SO3_BASE_HPP
#define FOUR_C_SO3_BASE_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_inpar_structure.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration ...
namespace STR
{
  namespace ELEMENTS
  {
    class ParamsInterface;
    enum EvalErrorFlag : int;
  }  // namespace ELEMENTS
}  // namespace STR
namespace Mat
{
  class So3Material;
}  // namespace Mat
namespace Discret
{
  namespace ELEMENTS
  {
    //! A wrapper for structural elements
    class SoBase : public Core::Elements::Element
    {
     public:
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id    (in): A globally unique element id
      \param owner (in): owner processor of the element
      */
      SoBase(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      SoBase(const SoBase& old);

      /*!
      \brief Default Constructor must not be called

      */
      SoBase() = delete;
      //@}

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Pack(Core::Communication::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element

      */
      void Unpack(const std::vector<char>& data) override;

      // get the kinematic type from the element
      Inpar::STR::KinemType KinematicType() const { return kintype_; }

      // get the kinematic type from the element
      void SetKinematicType(Inpar::STR::KinemType kintype) { kintype_ = kintype; }

      /*!
      \brief Does this element use EAS?

      ToDo: This function can be declared as pure virtual and each concrete derived
            class has to implement this function. This can be done during the up-coming
            cleaning procedure.                                      hiermeier 09/15
      */
      virtual bool HaveEAS() const { return false; };

      /*!
      \brief Return the material of this element

      Note: The input parameter nummat is not the material number from input file
            as in SetMaterial(int matnum), but the number of the material within
            the vector of materials the element holds

      \param nummat (in): number of requested material
      */
      virtual Teuchos::RCP<Mat::So3Material> SolidMaterial(int nummat = 0) const;

      /*!
       * @brief Evaluate Cauchy stress contracted with the normal vector and another direction
       * vector at given point in parameter space and calculate linearizations
       *
       * @param[in] xi            position in parameter space xi
       * @param[in] disp          vector of displacements
       * @param[in] n             normal vector (\f[\mathbf{n}\f])
       * @param[in] dir           direction vector (\f[\mathbf{v}\f]), can be either normal or
       *                          tangential vector
       * @param[out] cauchy_n_dir  cauchy stress tensor contracted using the vectors n and dir
       *                           (\f[ \mathbf{\sigma} \cdot \mathbf{n} \cdot \mathbf{v} \f])
       * @param[out] d_cauchyndir_dd    derivative of cauchy_n_dir w.r.t. displacements
       *                         (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}} { \mathrm{d} \mathbf{d}} \f])
       * @param[out] d2_cauchyndir_dd2  second derivative of cauchy_n_dir w.r.t. displacements
       *                      (\f[ \frac{ \mathrm{d}^2 \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}} { \mathrm{d} \mathbf{d}^2 } \f])
       * @param[out] d2_cauchyndir_dd_dn  second derivative of cauchy_n_dir w.r.t. displacements and
       *                                  vector n
       *                       (\f[ \frac{\mathrm{d}^2 \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}}
       *                                 {\mathrm{d} \mathbf{d} \mathrm{d} \mathbf{n} } \f])
       * @param[out] d2_cauchyndir_dd_ddir  second derivative of cauchy_n_dir w.r.t. displacements
       *                                    and direction vector v
       *                       (\f[ \frac{\mathrm{d}^2 \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}}
       *                                 {\mathrm{d} \mathbf{d} \mathrm{d} \mathbf{v} } \f])
       * @param[out] d2_cauchyndir_dd_dxi  second derivative of cauchy_n_dir w.r.t. displacements
       *                                   and local parameter coordinate xi
       *                       (\f[ \frac{\mathrm{d}^2 \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}}
       *                                 {\mathrm{d} \mathbf{d} \mathrm{d} \mathbf{\xi} } \f])
       * @param[out] d_cauchyndir_dn   derivative of cauchy_n_dir w.r.t. vector n
       *                        (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}} { \mathrm{d} \mathbf{n}} \f])
       * @param[out] d_cauchyndir_ddir  derivative of cauchy_n_dir w.r.t. direction vector v
       *                        (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}} { \mathrm{d} \mathbf{v}} \f])
       * @param[out] d_cauchyndir_dxi  derivative of cauchy_n_dir w.r.t. local parameter coord. xi
       *                        (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}} { \mathrm{d} \mathbf{\xi}} \f])
       * @param[in] temp                temperature
       * @param[out] d_cauchyndir_dT    derivative of cauchy_n_dir w.r.t. temperature
       *                        (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}} { \mathrm{d} T} \f])
       * @param[out] d2_cauchyndir_dd_dT  second derivative of cauchy_n_dir w.r.t. displacements
       *                                   and temperature
       *                       (\f[ \frac{\mathrm{d}^2 \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}}
       *                                 {\mathrm{d} \mathbf{d} \mathrm{d} T } \f])
       * @param[in] concentration     concentration
       * @param[out] d_cauchyndir_dc  derivative of cauchy_n_dir w.r.t. concentration
       *                        (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}} { \mathrm{d} c} \f])
       *
       * @note At the moment this method is only used for the nitsche contact formulation
       */
      virtual void get_cauchy_n_dir_and_derivatives_at_xi(const Core::LinAlg::Matrix<3, 1>& xi,
          const std::vector<double>& disp, const Core::LinAlg::Matrix<3, 1>& n,
          const Core::LinAlg::Matrix<3, 1>& dir, double& cauchy_n_dir,
          Core::LinAlg::SerialDenseMatrix* d_cauchyndir_dd,
          Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd2,
          Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd_dn,
          Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd_ddir,
          Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd_dxi,
          Core::LinAlg::Matrix<3, 1>* d_cauchyndir_dn,
          Core::LinAlg::Matrix<3, 1>* d_cauchyndir_ddir,
          Core::LinAlg::Matrix<3, 1>* d_cauchyndir_dxi, const std::vector<double>* temp,
          Core::LinAlg::SerialDenseMatrix* d_cauchyndir_dT,
          Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd_dT, const double* concentration,
          double* d_cauchyndir_dc)
      {
        FOUR_C_THROW("not implemented for chosen solid element");
      }

      /*!
       * @brief Evaluate Cauchy stress contracted with the normal vector and another direction
       * vector at given point in parameter space and calculate linearizations
       *
       * @param[in] xi            position in parameter space xi
       * @param[in] disp          vector of displacements
       * @param[in] n             normal vector (\f[\mathbf{n}\f])
       * @param[in] dir           direction vector (\f[\mathbf{v}\f]), can be either normal or
       *                          tangential vector
       * @param[out] cauchy_n_dir  cauchy stress tensor contracted using the vectors n and dir
       *                           (\f[ \mathbf{\sigma} \cdot \mathbf{n} \cdot \mathbf{v} \f])
       * @param[out] d_cauchyndir_dd    derivative of cauchy_n_dir w.r.t. displacements
       *                         (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}} { \mathrm{d} \mathbf{d}} \f])
       * @param[out] d2_cauchyndir_dd2  second derivative of cauchy_n_dir w.r.t. displacements
       *                      (\f[ \frac{ \mathrm{d}^2 \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}} { \mathrm{d} \mathbf{d}^2 } \f])
       * @param[out] d2_cauchyndir_dd_dn  second derivative of cauchy_n_dir w.r.t. displacements and
       *                                  vector n
       *                       (\f[ \frac{\mathrm{d}^2 \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}}
       *                                 {\mathrm{d} \mathbf{d} \mathrm{d} \mathbf{n} } \f])
       * @param[out] d2_cauchyndir_dd_ddir  second derivative of cauchy_n_dir w.r.t. displacements
       *                                    and direction vector v
       *                       (\f[ \frac{\mathrm{d}^2 \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}}
       *                                 {\mathrm{d} \mathbf{d} \mathrm{d} \mathbf{v} } \f])
       * @param[out] d2_cauchyndir_dd_dxi  second derivative of cauchy_n_dir w.r.t. displacements
       *                                   and local parameter coordinate xi
       *                       (\f[ \frac{\mathrm{d}^2 \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}}
       *                                 {\mathrm{d} \mathbf{d} \mathrm{d} \mathbf{\xi} } \f])
       * @param[out] d_cauchyndir_dn   derivative of cauchy_n_dir w.r.t. vector n
       *                        (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}} { \mathrm{d} \mathbf{n}} \f])
       * @param[out] d_cauchyndir_ddir  derivative of cauchy_n_dir w.r.t. direction vector v
       *                        (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}} { \mathrm{d} \mathbf{v}} \f])
       * @param[out] d_cauchyndir_dxi  derivative of cauchy_n_dir w.r.t. local parameter coord. xi
       *                        (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}} { \mathrm{d} \mathbf{\xi}} \f])
       * @param[in] temp                temperature
       * @param[out] d_cauchyndir_dT    derivative of cauchy_n_dir w.r.t. temperature
       *                        (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}} { \mathrm{d} T} \f])
       * @param[out] d2_cauchyndir_dd_dT  second derivative of cauchy_n_dir w.r.t. displacements
       *                                   and temperature
       *                       (\f[ \frac{\mathrm{d}^2 \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}}
       *                                 {\mathrm{d} \mathbf{d} \mathrm{d} T } \f])
       * @param[in] concentration     concentration
       * @param[out] d_cauchyndir_dc  derivative of cauchy_n_dir w.r.t. concentration
       *                        (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
       * \mathbf{v}} { \mathrm{d} c} \f])
       *
       * @note At the moment this method is only used for the nitsche contact formulation
       */
      virtual void get_cauchy_n_dir_and_derivatives_at_xi(const Core::LinAlg::Matrix<2, 1>& xi,
          const std::vector<double>& disp, const Core::LinAlg::Matrix<2, 1>& n,
          const Core::LinAlg::Matrix<2, 1>& dir, double& cauchy_n_dir,
          Core::LinAlg::SerialDenseMatrix* d_cauchyndir_dd,
          Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd2,
          Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd_dn,
          Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd_ddir,
          Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd_dxi,
          Core::LinAlg::Matrix<2, 1>* d_cauchyndir_dn,
          Core::LinAlg::Matrix<2, 1>* d_cauchyndir_ddir,
          Core::LinAlg::Matrix<2, 1>* d_cauchyndir_dxi, const std::vector<double>* temp,
          Core::LinAlg::SerialDenseMatrix* d_cauchyndir_dT,
          Core::LinAlg::SerialDenseMatrix* d2_cauchyndir_dd_dT, const double* concentration,
          double* d_cauchyndir_dc)
      {
        FOUR_C_THROW("not implemented for chosen solid element");
      }

      /** \brief set the parameter interface ptr for the solid elements
       *
       *  \param p (in): Parameter list coming from the time integrator.
       *
       *  \author hiermeier
       *  \date 04/16 */
      void set_params_interface_ptr(const Teuchos::ParameterList& p) override;

      /** \brief returns true if the parameter interface is defined and initialized, otherwise false
       *
       *  \author hiermeier
       *  \date 04/16 */
      inline bool IsParamsInterface() const override { return (not interface_ptr_.is_null()); }

      /** \brief get access to the parameter interface pointer
       *
       *  \author hiermeier
       *  \date 04/16 */
      Teuchos::RCP<Core::Elements::ParamsInterface> ParamsInterfacePtr() override;

     protected:
      /** \brief get access to the interface
       *
       *  \author hiermeier
       *  \date 04/16 */
      inline Core::Elements::ParamsInterface& params_interface()
      {
        if (not IsParamsInterface()) FOUR_C_THROW("The interface ptr is not set!");
        return *interface_ptr_;
      }

      /** \brief get access to the structure interface
       *
       *  \author vuong
       *  \date 11/16 */
      STR::ELEMENTS::ParamsInterface& str_params_interface();

      /** \brief error handling for structural elements
       *
       *  \author hiermeier \date 09/18 */
      void error_handling(const double& det_curr, Teuchos::ParameterList& params, int line_id,
          STR::ELEMENTS::EvalErrorFlag flag);

     protected:
      /*!
       * \brief This method executes the material_post_setup if not already executed.
       *
       * This method should be placed in the Evaluate call. It will internally check, whether the
       * material post_setup() routine was already called in if not, it invokes this call directly.
       *
       * @param params Container for additional information
       */
      void ensure_material_post_setup(Teuchos::ParameterList& params);

      /*!
       * \brief This method calls the post_setup routine of all materials.
       *
       * It can be used to pass information from the element to the materials after everything
       * is set up. For a simple element, the ParameterList is passed unchanged to the materials.
       */
      virtual void material_post_setup(Teuchos::ParameterList& params);

      //! kinematic type
      Inpar::STR::KinemType kintype_;

     private:
      /** \brief interface ptr
       *
       *  data exchange between the element and the time integrator. */
      Teuchos::RCP<Core::Elements::ParamsInterface> interface_ptr_;

      //! Flag of the status of the material post setup routine
      bool material_post_setup_;

    };  // class So_base

  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif

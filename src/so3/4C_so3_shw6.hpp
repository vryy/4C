/*----------------------------------------------------------------------*/
/*! \file
\brief
\level 1


*----------------------------------------------------------------------*/
#ifndef FOUR_C_SO3_SHW6_HPP
#define FOUR_C_SO3_SHW6_HPP

#include "4C_config.hpp"

#include "4C_linalg_serialdensematrix.hpp"
#include "4C_so3_weg6.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace ELEMENTS
  {
    class SoShw6Type : public Core::Elements::ElementType
    {
     public:
      std::string Name() const override { return "So_shw6Type"; }

      static SoShw6Type& Instance();

      Core::Communication::ParObject* Create(const std::vector<char>& data) override;

      Teuchos::RCP<Core::Elements::Element> Create(const std::string eletype,
          const std::string eledistype, const int id, const int owner) override;

      Teuchos::RCP<Core::Elements::Element> Create(const int id, const int owner) override;

      int initialize(Core::FE::Discretization& dis) override;

      void nodal_block_information(
          Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np) override;

      Core::LinAlg::SerialDenseMatrix ComputeNullSpace(
          Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp) override;

      void setup_element_definition(
          std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
          override;

     private:
      static SoShw6Type instance_;

      std::string get_element_type_string() const { return "SOLIDSHW6"; }
    };

    /*!
    \brief A C++ 6-node wedge Solid-Shell element inherited from so_weg6

    This element adopts the ANS approach to improve shear and curvature-thickness
    (trapezoidal) locking. However, being triangular inplane, the ANS evaluation
    points of the shear strains are choosen just on two of three surfaces (See
    also Discret::ELEMENTS::So_shw6::soshw6_anssetup).
    Thus, the element becomes mesh dependent (where in material space lies my
    parameter space), but this is in favour of less pronounced shear locking
    compared to a invariant formulation (locking free just for certain mesh configurations).

    See Diss. Frank Koschnick for details.

    */
    class SoShw6 : public SoWeg6
    {
     public:
      //! @name Friends
      friend class SoShw6Type;
      friend class Soweg6Surface;
      friend class Soweg6Line;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner : elements owning processor
      */
      SoShw6(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      SoShw6(const SoShw6& old);

      /*!
      \brief Deep copy this instance of Solid3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      Core::Elements::Element* Clone() const override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      int UniqueParObjectId() const override { return SoShw6Type::Instance().UniqueParObjectId(); }

      /*!
      \brief Pack this class so it can be communicated

      \ref pack and \ref unpack are used to communicate this element

      */
      void pack(Core::Communication::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref pack and \ref unpack are used to communicate this element

      */
      void unpack(const std::vector<char>& data) override;


      //@}

      //! @name Acess methods

      /*!
       \brief Does this element use EAS?
       */
      bool HaveEAS() const override { return (eastype_ != soshw6_easnone); };


      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const override;

      Core::Elements::ElementType& ElementType() const override { return SoShw6Type::Instance(); }

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          Input::LineDefinition* linedef) override;

      //@}

      //! @name Evaluation



     private:
      /*! \brief EAS technology enhancement types of so_shw6
       * Solid-Shell Wedge6 has EAS enhancement of GL-strains to avoid locking
       * \param soshw6_poissonthick: 1 parameter for correct kinematic of solid-shell,
       *                             allows linear strain in thick-dir
       * \param soshw6_none:         no EAS
       */
      enum EASType
      {
        soshw6_easpoisthick = 1,
        soshw6_easnone = 0
      };

      //! type of EAS technology
      EASType eastype_;

      struct EASData
      {
        Core::LinAlg::SerialDenseMatrix alpha;
        Core::LinAlg::SerialDenseMatrix alphao;
        Core::LinAlg::SerialDenseMatrix feas;
        Core::LinAlg::SerialDenseMatrix invKaa;
        Core::LinAlg::SerialDenseMatrix Kda;
        Core::LinAlg::SerialDenseMatrix eas_inc;
      };
      //! EAS data
      EASData easdata_;

      //! number of EAS parameters (alphas), defined by 'EASType"
      int neas_{};
      // line search parameter (old step length)
      double old_step_length_{};

      //! reorder nodes to optimally align material space with parameter space
      bool optimal_parameterspace_map_{};

      //! in case of changed nodal order this is true
      bool nodes_rearranged_{};

      //! don't want = operator
      SoShw6& operator=(const SoShw6& old);

      /*!
      \brief Evaluate an element

      Evaluate so_hex8 element stiffness, mass, internal forces, etc.

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param discretization : pointer to discretization for de-assembly
      \param lm (in)        : location matrix for de-assembly
      \param elemat1 (out)  : (stiffness-)matrix to be filled by element. If nullptr on input,
                              the controling method does not expect the element to fill
                              this matrix.
      \param elemat2 (out)  : (mass-)matrix to be filled by element. If nullptr on input,
                              the controling method does not expect the element to fill
                              this matrix.
      \param elevec1 (out)  : (internal force-)vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If nullptr on input,
                              the controlling method does not expect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      int evaluate(Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          std::vector<int>& lm, Core::LinAlg::SerialDenseMatrix& elemat1,
          Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
          Core::LinAlg::SerialDenseVector& elevec2,
          Core::LinAlg::SerialDenseVector& elevec3) override;

      //! Compute stiffness and mass matrix
      virtual void soshw6_nlnstiffmass(std::vector<int>& lm,            // location matrix
          std::vector<double>& disp,                                    // current displacements
          std::vector<double>& residual,                                // current residual displ
          Core::LinAlg::Matrix<NUMDOF_WEG6, NUMDOF_WEG6>* stiffmatrix,  // element stiffness matrix
          Core::LinAlg::Matrix<NUMDOF_WEG6, NUMDOF_WEG6>* massmatrix,   // element mass matrix
          Core::LinAlg::Matrix<NUMDOF_WEG6, 1>* force,      // element internal force vector
          Core::LinAlg::Matrix<NUMDOF_WEG6, 1>* force_str,  // structure force
          Core::LinAlg::Matrix<NUMGPT_WEG6, Mat::NUM_STRESS_3D>* elestress,  // stresses at GP
          Core::LinAlg::Matrix<NUMGPT_WEG6, Mat::NUM_STRESS_3D>* elestrain,  // strains at GP
          Teuchos::ParameterList& params,          // algorithmic parameters e.g. time
          const Inpar::STR::StressType iostress,   // stress output option
          const Inpar::STR::StrainType iostrain);  // strain output option

      static constexpr int num_sp = 5;   ///< number of ANS sampling points
      static constexpr int num_ans = 3;  ///< number of modified ANS strains (E_rt,E_st,E_tt)

      //! Setup ANS interpolation (shear and trapezoidal)
      void soshw6_anssetup(
          const Core::LinAlg::Matrix<NUMNOD_WEG6, NUMDIM_WEG6>& xrefe,  ///< material element coords
          const Core::LinAlg::Matrix<NUMNOD_WEG6, NUMDIM_WEG6>& xcurr,  ///< current element coords
          std::vector<Core::LinAlg::Matrix<NUMDIM_WEG6, NUMNOD_WEG6>>**
              deriv_sp,  ///< derivs eval. at all sampling points
          std::vector<Core::LinAlg::Matrix<NUMDIM_WEG6, NUMDIM_WEG6>>&
              jac_sps,  ///< jac at all sampling points
          std::vector<Core::LinAlg::Matrix<NUMDIM_WEG6, NUMDIM_WEG6>>&
              jac_cur_sps,  ///< current jac at all sampling points
          Core::LinAlg::Matrix<num_ans * num_sp, NUMDOF_WEG6>& B_ans_loc);  ///< modified B

      //! Transformation matrix parameter->material space
      void soshw6_evaluate_t(const Core::LinAlg::Matrix<NUMDIM_WEG6, NUMDIM_WEG6>& jac,
          Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, Mat::NUM_STRESS_3D>& TinvT);

      //! EAS technology, init
      void soshw6_easinit();

      //! EAS technology, setup necessary data
      void soshw6_eassetup(
          std::vector<Core::LinAlg::SerialDenseMatrix>** M_GP,  // M-matrix evaluated at GPs
          double& detJ0,                                        // det of Jacobian at origin
          Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, Mat::NUM_STRESS_3D>&
              T0invT,  // maps M(origin) local to global
          const Core::LinAlg::Matrix<NUMNOD_WEG6, NUMDIM_WEG6>& xrefe);  // material element coords

      //! Evaluate Cauchy stress (mapping with disp_based F per default)
      void soshw6_cauchy(Core::LinAlg::Matrix<NUMGPT_WEG6, Mat::NUM_STRESS_3D>* elestress,
          const int gp, const Core::LinAlg::Matrix<NUMDIM_WEG6, NUMDIM_WEG6>& defgrd,
          const Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1>& glstrain,
          const Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1>& stress);

      //! Reorient wedge to optimally align material with parameter space
      int soshw6_findoptparmap();

      //! recover condensed EAS variables
      void soshw6_recover(const std::vector<double>& residual);

      void pack_eas_data(Core::Communication::PackBuffer& data) const
      {
        add_to_pack(data, easdata_.alpha);
        add_to_pack(data, easdata_.alphao);
        add_to_pack(data, easdata_.feas);
        add_to_pack(data, easdata_.invKaa);
        add_to_pack(data, easdata_.Kda);
        add_to_pack(data, easdata_.eas_inc);
      };

      void unpack_eas_data(std::vector<char>::size_type& position, const std::vector<char>& data)
      {
        extract_from_pack(position, data, easdata_.alpha);
        extract_from_pack(position, data, easdata_.alphao);
        extract_from_pack(position, data, easdata_.feas);
        extract_from_pack(position, data, easdata_.invKaa);
        extract_from_pack(position, data, easdata_.Kda);
        extract_from_pack(position, data, easdata_.eas_inc);
      };

     private:
      std::string get_element_type_string() const { return "SOLIDSHW6"; }
    };  // class So_shw6


  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif

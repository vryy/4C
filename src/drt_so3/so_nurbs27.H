/*!----------------------------------------------------------------------
\file so_nurbs27.H

\brief implementation of the quadratic NURBS 27 element

\level 2

\maintainer Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271

*----------------------------------------------------------------------*/
#ifndef SONURBS27_H
#define SONURBS27_H


#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_node.H"

#include "so_base.H"

//@=
#include "../drt_inpar/inpar_structure.H"

#include "../drt_mat/material.H"

#include "../linalg/linalg_fixedsizematrix.H"

/// Several parameters which are fixed for Solid nurbs27
const int NUMNOD_SONURBS27 = 27;  ///< number of nodes
const int NODDOF_SONURBS27 = 3;   ///< number of dofs per node
const int NUMDOF_SONURBS27 = 81;  ///< total dofs per element
const int NUMGPT_SONURBS27 = 27;  ///< total gauss points per element
const int NUMDIM_SONURBS27 = 3;   ///< number of dimensions


namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    /*!
    \brief NURBS: namespace of all isogeometric analysis tools

    */
    namespace NURBS
    {
      class So_nurbs27Type : public DRT::ElementType
      {
       public:
        std::string Name() const { return "So_nurbs27Type"; }

        static So_nurbs27Type& Instance();

        virtual DRT::ParObject* Create(const std::vector<char>& data);

        virtual Teuchos::RCP<DRT::Element> Create(
            const std::string eletype, const std::string eledistype, const int id, const int owner);

        virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

        virtual int Initialize(DRT::Discretization& dis);

        virtual void NodalBlockInformation(
            DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

        virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
            const double* x0, int numdf, int dimns);

        virtual void SetupElementDefinition(
            std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

       private:
        static So_nurbs27Type instance_;
      };


      /*!
      \brief A C++ version of a 27-control point nurbs solid element

      A structural 27-control points nurbs solid displacement element for large deformations.
      As its discretization is fixed many data structures are evaluated just once and kept
      for performance.

      \author gammi
      */
      class So_nurbs27 : public So_base
      {
       public:
        //! @name Friends
        friend class So_nurbs27Type;

        //@}
        //! @name Constructors and destructors and related methods

        /*!
        \brief Standard Constructor

        \param id : A unique global id
        \param owner : elements owner
        */
        So_nurbs27(int id, int owner);

        /*!
        \brief Copy Constructor

        Makes a deep copy of a Element

        */
        So_nurbs27(const So_nurbs27& old);

        /*!
        \brief Deep copy this instance of Solid3 and return pointer to the copy

        The Clone() method is used from the virtual base class Element in cases
        where the type of the derived class is unknown and a copy-ctor is needed

        */
        DRT::Element* Clone() const;

        /*!
        \brief Get shape type of element
        */
        virtual DiscretizationType Shape() const;

        /*!
        \brief Return number of volumes of this element
        */
        virtual int NumVolume() const { return 1; }

        /*!
        \brief Return number of surfaces of this element
        */
        virtual int NumSurface() const { return 6; }

        /*!
        \brief Return number of lines of this element
        */
        virtual int NumLine() const { return 12; }

        /*!
        \brief Get vector of Teuchos::RCPs to the lines of this element

        */
        virtual std::vector<Teuchos::RCP<DRT::Element>> Lines();

        /*!
        \brief Get vector of Teuchos::RCPs to the surfaces of this element

        */
        virtual std::vector<Teuchos::RCP<DRT::Element>> Surfaces();

        /*!
        \brief Get vector of Teuchos::RCPs to the volumes of this element

        */
        virtual std::vector<Teuchos::RCP<DRT::Element>> Volumes();

        /*!
        \brief Return unique ParObject id

        every class implementing ParObject needs a unique id defined at the
        top of this file.
        */
        virtual int UniqueParObjectId() const
        {
          return So_nurbs27Type::Instance().UniqueParObjectId();
        }

        /*!
        \brief Pack this class so it can be communicated

        \ref Pack and \ref Unpack are used to communicate this element

        */
        virtual void Pack(DRT::PackBuffer& data) const;

        /*!
        \brief Unpack data from a char vector into this class

        \ref Pack and \ref Unpack are used to communicate this element

        */
        virtual void Unpack(const std::vector<char>& data);

        /*!
        \brief Destructor

        */
        virtual ~So_nurbs27();

        //@}

        //! @name Access methods


        /*!
        \brief Get number of degrees of freedom of a certain node
               (implements pure virtual DRT::Element)

        The element decides how many degrees of freedom its nodes must have.
        As this may vary along a simulation, the element can redecide the
        number of degrees of freedom per node along the way for each of it's nodes
        separately.
        */
        virtual int NumDofPerNode(const DRT::Node& node) const { return 3; }

        /*!
        \brief Get number of degrees of freedom per element
               (implements pure virtual DRT::Element)

        The element decides how many element degrees of freedom it has.
        It can redecide along the way of a simulation.

        \note Element degrees of freedom mentioned here are dofs that are visible
              at the level of the total system of equations. Purely internal
              element dofs that are condensed internally should NOT be considered.
        */
        virtual int NumDofPerElement() const { return 0; }

        /*!
        \brief Print this element
        */
        void Print(std::ostream& os) const;

        virtual DRT::ElementType& ElementType() const { return So_nurbs27Type::Instance(); }

        //! @name Input and Creation

        /*!
        \brief Read input for this element
        */
        virtual bool ReadElement(const std::string& eletype, const std::string& distype,
            DRT::INPUT::LineDefinition* linedef);

        //@}

        //! @name Evaluation

        /*!
        \brief Evaluate an element

        Evaluate So_nurbs27 element stiffness, mass, internal forces, etc.

        If NULL on input, the controling method does not expect the element
        to fill these matrices or vectors.

        \return 0 if successful, negative otherwise
        */
        int Evaluate(
            Teuchos::ParameterList&
                params,  ///< ParameterList for communication between control routine and elements
            DRT::Discretization& discretization,  ///< pointer to discretization for de-assembly
            std::vector<int>& lm,                 ///< location matrix for de-assembly
            Epetra_SerialDenseMatrix& elemat1,    ///< (stiffness-)matrix to be filled by element.
            Epetra_SerialDenseMatrix& elemat2,    ///< (mass-)matrix to be filled by element.
            Epetra_SerialDenseVector& elevec1,  ///< (internal force-)vector to be filled by element
            Epetra_SerialDenseVector& elevec2,  ///< vector to be filled by element
            Epetra_SerialDenseVector& elevec3   ///< vector to be filled by element
        );


        /*!
        \brief Evaluate a Neumann boundary condition

        this method evaluates a surface Neumann condition on the solid3 element

        \param params (in/out)    : ParameterList for communication between control routine
                                    and elements
        \param discretization (in): A reference to the underlying discretization
        \param condition (in)     : The condition to be evaluated
        \param lm (in)            : location vector of this element
        \param elevec1 (out)      : vector to be filled by element. If NULL on input,

        \return 0 if successful, negative otherwise
        */
        virtual int EvaluateNeumann(Teuchos::ParameterList& params,
            DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
            Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1 = NULL);

        //@}

        /*!
         \brief calculate a. scaled director matrix for thin shell-like structures
           or its inverse

         \return void
         */
        void CalcSTCMatrix(LINALG::Matrix<81, 81>& elemat1, const INPAR::STR::STC_Scale stc_scaling,
            const int stc_layer, std::vector<int>& lm, DRT::Discretization& discretization,
            bool do_inverse);

       protected:
        //! action parameters recognized by So_nurbs27
        enum ActionType
        {
          none,
          calc_struct_linstiff,
          calc_struct_nlnstiff,
          calc_struct_internalforce,
          calc_struct_linstiffmass,
          calc_struct_nlnstiffmass,
          calc_struct_eleload,
          calc_struct_fsiload,
          calc_struct_update_istep,
          calc_stc_matrix,  //! calculate scaled director matrix for thin shell-like structures
          calc_stc_matrix_inverse,  //! calculate inverse scaled director matrix for thin shell-like
                                    //! structures
          calc_struct_reset_istep,
          calc_struct_energy,
          calc_struct_stifftemp,  //!< TSI specific: mechanical-thermal stiffness
          calc_struct_nlnstifflmass
        };

        //! container containing all kind of solid3 thisandthat
        DRT::Container data_;
        //! vector of inverses of the jacobian in material frame
        std::vector<LINALG::Matrix<3, 3>> invJ_;
        //! determinant of Jacobian in material frame
        std::vector<double> detJ_;


        // internal calculation methods

        //! don't want = operator
        So_nurbs27& operator=(const So_nurbs27& old);

        //! init the inverse of the jacobian and its determinant in
        //! the material configuration
        virtual void InitJacobianMapping(DRT::Discretization& dis);

        //! Calculate nonlinear stiffness and mass matrix
        virtual void sonurbs27_nlnstiffmass(std::vector<int>& lm,  ///< location matrix
            DRT::Discretization& discretization,  ///< discretisation to extract knot vector
            std::vector<double>& disp,            ///< current displacements
            std::vector<double>& residual,        ///< current residual displ
            LINALG::Matrix<81, 81>* stiffmatrix,  ///< element stiffness matrix
            LINALG::Matrix<81, 81>* massmatrix,   ///< element mass matrix
            LINALG::Matrix<81, 1>* force,         ///< element internal force vector
            Teuchos::ParameterList& params);

        //! Evaluate Nurbs27 Shapefcts to keep them static
        const std::vector<LINALG::Matrix<27, 1>> sonurbs27_shapefcts(
            const std::vector<Epetra_SerialDenseVector>& myknots,
            const LINALG::Matrix<27, 1>& weights);
        //! Evaluate Nurbs27 Derivs to keep them static
        const std::vector<LINALG::Matrix<3, 27>> sonurbs27_derivs(
            const std::vector<Epetra_SerialDenseVector>& myknots,
            const LINALG::Matrix<27, 1>& weights);
        //! Evaluate Nurbs27 Weights to keep them static
        const std::vector<double> sonurbs27_gpweights();

        //! calculate internal energy
        double CalcIntEnergy(DRT::Discretization& discretization, std::vector<double>& disp,
            Teuchos::ParameterList& params);

        //! Lump mass matrix (bborn 07/08)
        void lumpmass(LINALG::Matrix<NUMDOF_SONURBS27, NUMDOF_SONURBS27>* emass);

      };  // class So_nurbs27



    }  // namespace NURBS
  }    // namespace ELEMENTS
}  // namespace DRT

#endif  // #ifndef SONURBS27_H

/*----------------------------------------------------------------------*/
/*!
\file so3_thermo.H
\brief 3d TSI solid element
\level 1
\maintainer Christoph Meier
*/


/*----------------------------------------------------------------------*
 | definitions                                               dano 11/12 |
 *----------------------------------------------------------------------*/
#ifndef SO3_THERMO_H
#define SO3_THERMO_H

/*----------------------------------------------------------------------*
 | headers                                                   dano 11/12 |
 *----------------------------------------------------------------------*/
#include "so3_thermo_eletypes.H"

// include thermal header because number of Gauss points are determined
// dependent on distype
#include "../drt_thermo/thermo_ele_impl_utils.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"

/*----------------------------------------------------------------------*
 |                                                           dano 11/12 |
 *----------------------------------------------------------------------*/
namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    //! A C++ version of a 3 dimensional solid element with modifications for
    //! thermomechanics
    //!
    //! A structural 3 dimensional solid displacement element for large deformations
    //! and with small and large strains
    template <class so3_ele, DRT::Element::DiscretizationType distype>
    class So3_Thermo : public so3_ele
    {
      //! @name Friends
      friend class So_hex8ThermoType;
      friend class So_hex8fbarThermoType;
      friend class So_tet4ThermoType;
      friend class So_tet10ThermoType;
      friend class So_hex27ThermoType;
      friend class So_hex20ThermoType;
      friend class So_nurbs27ThermoType;

     public:
      //@}
      //! @name Constructors and destructors and related methods


      //! Standard Constructor
      So3_Thermo(int id,  //!< (i) this element's global id
          int owner       //!< elements owner
      );

      //! Copy Constructor
      //! Makes a deep copy of a Element
      So3_Thermo(const So3_Thermo& old);

      //! Destructor
      virtual ~So3_Thermo(){};

      //@}

      //! number of element nodes
      static const int nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
      //! number of space dimensions
      static const int nsd_ = 3;
      //! number of dofs per node
      static const int numdofpernode_ = 3;
      //! total dofs per element
      static const int numdofperelement_ = numdofpernode_ * nen_;
      //! number of strains/stresses
      static const int numstr_ = 6;
      //! number of Gauss points per element (value is added in so3_thermo.cpp)
      int numgpt_;
      //! static const is required for fixedsizematrices
      //! TODO maybe more beauty is possible
      static const int numgpt_post = THR::DisTypeToSTRNumGaussPoints<distype>::nquad;

      //! @name Acess methods

      //! Deep copy this instance of Solid3 and return pointer to the copy
      //!
      //! The Clone() method is used from the virtual base class Element in cases
      //! where the type of the derived class is unknown and a copy-ctor is needed
      DRT::Element* Clone() const;

      //! Return unique ParObject id
      //!
      //! every class implementing ParObject needs a unique id defined at the top of
      //! this file.
      virtual int UniqueParObjectId() const;

      //! Pack this class so it can be communicated
      //! Pack and \ref Unpack are used to communicate this element
      virtual void Pack(DRT::PackBuffer& data) const;

      //! Unpack data from a char vector into this class
      //! Pack and \ref Unpack are used to communicate this element
      virtual void Unpack(const std::vector<char>& data);

      //@}

      //! @name Access methods

      //! Print this element
      void Print(std::ostream& os) const;

      //! return elementtype thermo element
      virtual DRT::ElementType& ElementType() const;

      //@}

      //! @name Input and Creation

      //! Query names of element data to be visualized using BINIO
      //!
      //! The element fills the provided map with key names of
      //! visualization data the element wants to visualize AT THE CENTER
      //! of the element geometry. The values is supposed to be dimension of the
      //! data to be visualized. It can either be 1 (scalar), 3 (vector), 6 (sym. tensor)
      //! or 9 (nonsym. tensor)
      //!
      //! Example:
      //! \code
      //!   // Name of data is 'Owner', dimension is 1 (scalar value)
      //!   names.insert(std::pair<std::string,int>("Owner",1));
      //!   // Name of data is 'StressesXYZ', dimension is 6 (sym. tensor value)
      //!   names.insert(std::pair<std::string,int>("StressesXYZ",6));
      //! \endcode
      //!
      //!  names (out): On return, the derived class has filled names with key
      //!               names of data it wants to visualize and with int dimensions
      //!               of that data.
      virtual void VisNames(std::map<std::string, int>& names);

      //!  Query data to be visualized using BINIO of a given name
      //!
      //! The method is supposed to call this base method to visualize the owner of
      //! the element.
      //! If the derived method recognizes a supported data name, it shall fill it
      //! with corresponding data.
      //! If it does NOT recognizes the name, it shall do nothing.
      //!
      //! warning: the method must not change size of data
      //!
      //!  name (in):   Name of data that is currently processed for visualization
      //! \param data (out):  data to be filled by element if element recognizes the name
      virtual bool VisData(const std::string& name, std::vector<double>& data);

      //! read input for this element
      bool ReadElement(const std::string& eletype,  //!< so3thermo
          const std::string& eledistype,            //!< hex8,tet4,...
          DRT::INPUT::LineDefinition* linedef       //!< what parameters have to be read
      );

      //@}

      //! @name Evaluation

      //! evaluate an element
      //! evaluate element stiffness, mass, internal forces, etc.
      //!
      //! if NULL on input, the controlling method does not expect the element
      //!  to fill these matrices or vectors.
      //!
      //!  \return 0 if successful, negative otherwise
      int Evaluate(
          Teuchos::ParameterList&
              params,  //!< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  //!< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la,      //!< location array for de-assembly
          Epetra_SerialDenseMatrix& elemat1,    //!< (stiffness-)matrix to be filled by element.
          Epetra_SerialDenseMatrix& elemat2,    //!< (mass-)matrix to be filled by element.
          Epetra_SerialDenseVector& elevec1,    //!< (internal force-)vector to be filled by element
          Epetra_SerialDenseVector& elevec2,    //!< vector to be filled by element
          Epetra_SerialDenseVector& elevec3     //!< vector to be filled by element
      );

      // PreEvaluate undertakes the task to calculate coupling term of matrix
      //! evaluation
      void PreEvaluate(
          Teuchos::ParameterList&
              params,  //!< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  //!< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la       //!< location array for de-assembly
      );

      //@}

      //! init the inverse of the jacobian and its determinant in the material
      //! configuration
      virtual void InitJacobianMapping(DRT::Discretization& dis);

      //@}

     protected:
      //! don't want = operator
      So3_Thermo& operator=(const So3_Thermo& old);

      //! action parameters recognized by so_hex8
      enum ActionType
      {
        none,
        calc_struct_linstiff,       //!< geometrical linear analysis: stiffness matrix
        calc_struct_nlnstiff,       //!< stiffness matrix
        calc_struct_internalforce,  //!< internal force
        calc_struct_linstiffmass,   //!< geometrical linear analysis: internal force,
                                    //!< its stiffness and mass matrix
        calc_struct_nlnstiffmass,   //!< internal force, its stiffness and mass matrix
        calc_struct_nlnstifflmass,  //!< internal force, its stiffness and lumped mass matrix
        calc_struct_stress,         //!< calculate stresses
        calc_struct_stifftemp,      //!< calculate coupling term k_dT for monolithic TSI
        calc_struct_update_istep,   //!< update all at element level
        calc_struct_reset_istep,    //!< reset elementwise internal variables
                                    //!< during iteration to last converged state
                                    //!< needed for predictor TangDis
        postprocess_stress,         //!< postprocess stresses
        calc_struct_energy          //!< compute internal energy
      };

      //! vector of inverses of the jacobian in material frame
      std::vector<LINALG::Matrix<nsd_, nsd_>> invJ_;
      //! determinant of Jacobian in material frame
      std::vector<double> detJ_;
      //! container containing all kind of solid3 thisandthat
      DRT::Container data_;
      //! vector of coordinates of current integration point in reference coordinates
      std::vector<LINALG::Matrix<nsd_, 1>> xsi_;

      DRT::UTILS::GaussIntegration intpoints_;

      //! flag indicating higher order element
      bool ishigherorder_;

      //! @name TSI related stuff

      //! evaluate an element
      //!
      //! evaluate So3_Thermo element stiffness, mass, internal forces, etc.
      //! templated evaluate routine of element matrixes
      //!
      //! if NULL on input, the controlling method does not expect the element
      //! to fill these matrices or vectors.
      //!
      //! \return 0 if successful, negative otherwise
      int EvaluateCouplWithThr(
          Teuchos::ParameterList&
              params,  //!< ParameterList for communication between control routine and elements
          DRT::Discretization& discretization,  //!< pointer to discretization for de-assembly
          DRT::Element::LocationArray& la,      //!< location array for de-assembly
          Epetra_SerialDenseMatrix& elemat1,    //!< (stiffness-)matrix to be filled by element.
          Epetra_SerialDenseMatrix& elemat2,    //!< (mass-)matrix to be filled by element.
          Epetra_SerialDenseVector& elevec1,    //!< (internal force-)vector to be filled by element
          Epetra_SerialDenseVector& elevec2,    //!< vector to be filled by element
          Epetra_SerialDenseVector& elevec3     //!< vector to be filled by element
      );

      //! Calculate temperature coupling term for the internal force (geometric linear)
      virtual void lin_fint_tsi(DRT::Element::LocationArray& la,  //!< location array
          std::vector<double>& disp,                              //!< current displacements
          std::vector<double>& temp,                              //!< current temperature
          LINALG::Matrix<numdofperelement_, 1>* force,            //!< element internal force vector
          LINALG::Matrix<numgpt_post, numstr_>* elestress,        //!< stresses at GP
          Teuchos::ParameterList& params,        //!< algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress  //!< stress output option
      );

      //! Calculate mechanical thermal stiffness term needed for monolithic TSI K_dT
      virtual void lin_kdT_tsi(DRT::Element::LocationArray& la,
          std::vector<double>& disp,  //!< (i): current displacement
          std::vector<double>& temp,  // current temperatures
          LINALG::Matrix<numdofperelement_, nen_>*
              stiffmatrix_kdT,  //!< (o): mechanical thermal stiffness term at current gp
          Teuchos::ParameterList& params);

      //! Calculate nonlinear stiffness and mass matrix with temperature fraction
      virtual void nln_stifffint_tsi(DRT::Element::LocationArray& la,  //!< location array
          DRT::Discretization& discretization,  ///< discretisation to extract knot vector
          std::vector<double>& disp,            //!< current displacements
          std::vector<double>& temp,            //!< current temperature
          LINALG::Matrix<numdofperelement_, numdofperelement_>*
              stiffmatrix,                                  // element stiffness matrix
          LINALG::Matrix<numdofperelement_, 1>* force,      //!< element internal force vector
          LINALG::Matrix<numgpt_post, numstr_>* elestress,  //!< stresses at GP
          Teuchos::ParameterList& params,                   //!< algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress             //!< stress output option
      );

      //! Calculate mechanical thermal stiffness term needed for monolithic TSI K_dT
      virtual void nln_kdT_tsi(DRT::Element::LocationArray& la,
          DRT::Discretization& discretization,  ///< discretisation to extract knot vector
          std::vector<double>& disp,            //!< (i): current displacement
          std::vector<double>& temp,            //!< current temperature
          LINALG::Matrix<numdofperelement_, nen_>*
              stiffmatrix_kdT,  //!< (o): mechanical thermal stiffness term at current gp
          Teuchos::ParameterList& params);

      //! \brief Interface to the temperature dependent material law working with soh8.
      //!
      //! Here the interface to the first solid material with temperature takes place.
      //! Stress and material tangent must be retrieved whereas all necessary data
      //! is passed.
      //! Add whatever your material needs, but make sure that exchange is not
      //! overdone performance-wise.
      void Materialize(
          LINALG::Matrix<numstr_, 1>* stresstemp,  //!< (o): Voigt-Vector of stresses at current gp
          LINALG::Matrix<numstr_, 1>*
              ctemp,  //!< (o): temperature dependent tangent matrix at current gp
          LINALG::Matrix<1, 1>*
              Ntemp,  // element temperature: (shapefcts . element temperature) at current gp
          LINALG::Matrix<numstr_, numstr_>* cmat,  //!< (o): material tangent matrix at current gp
          LINALG::Matrix<numstr_, 1>* glstrain,    //!< (o): total strain
          Teuchos::ParameterList& params  //!< parameter list to access time, etc. in materials
      );

      //! calculate the constant temperature tangent for stresstemp
      void Ctemp(LINALG::Matrix<numstr_, 1>* ctemp,  //!< temperature dependent material tangent
          Teuchos::ParameterList& params  //!< parameter list to access time, etc. in materials
      );
      // TODO this should really not be necessary if we use one consistent GP definition
      // throughout baci
      //! map the Intrepid gp numbering to the so_hex8 numbering, do nothing if not hex8
      int MapMyGpToSoHex8(int myGp);
      //@}

      //! calculate nonlinear B-operator (6x24)
      void CalculateBop(LINALG::Matrix<numstr_, numdofperelement_>* bop,
          LINALG::Matrix<nsd_, nsd_>* defgrd, LINALG::Matrix<nsd_, nen_>* N_XYZ);

      //! calculates nonlinear B-operator in vector notation (1x24)
      void CalculateBopVec(LINALG::Matrix<1, numdofperelement_>& bopvec,
          LINALG::Matrix<nsd_, nsd_>& defgrd, LINALG::Matrix<nsd_, nen_>& N_XYZ);

      //! calculate linear B-operator
      void CalculateBoplin(
          LINALG::Matrix<numstr_, numdofperelement_>* boplin, LINALG::Matrix<nsd_, nen_>* N_XYZ);

      //! push forward of material stresses to the current, spatial configuration
      void PK2toCauchy(LINALG::Matrix<numstr_, 1>* stress, LINALG::Matrix<nsd_, nsd_>* defgrd,
          LINALG::Matrix<nsd_, nsd_>* cauchystress);

      //! push forward of Green-Lagrange strain to Euler-Almansi strains
      void GLtoEA(LINALG::Matrix<numstr_, 1>* glstrain, LINALG::Matrix<nsd_, nsd_>* defgrd,
          LINALG::Matrix<nsd_, nsd_>* euler_almansi);

      //! @name TSI and thermoplasticity related stuff

      //! Calculate nonlinear stiffness and mass matrix with temperature fraction
      //! implementation for hex8fbar elements differs from standard implementation
      virtual void nln_stifffint_tsi_fbar(DRT::Element::LocationArray& la,  //!< location array
          std::vector<double>& disp,  //!< current displacements
          std::vector<double>& temp,  //!< current temperature
          LINALG::Matrix<numdofperelement_, numdofperelement_>*
              stiffmatrix,                                  // element stiffness matrix
          LINALG::Matrix<numdofperelement_, 1>* force,      //!< element internal force vector
          LINALG::Matrix<numgpt_post, numstr_>* elestress,  //!< stresses at GP
          Teuchos::ParameterList& params,                   //!< algorithmic parameters e.g. time
          const INPAR::STR::StressType iostress             //!< stress output option
      );

      //! Calculate mechanical thermal stiffness term needed for monolithic TSI K_dT
      virtual void nln_kdT_tsi_fbar(DRT::Element::LocationArray& la,
          std::vector<double>& disp,  //!< (i): current displacement
          std::vector<double>& temp,  //!< current temperature
          LINALG::Matrix<numdofperelement_, nen_>*
              stiffmatrix_kdT,  //!< (o): mechanical thermal stiffness term at current gp
          Teuchos::ParameterList& params);

      //@}

     private:
      DRT::Node** Nodes();

      Teuchos::RCP<MAT::Material> Material() const;

      int Id() const;

    };  // class So3_Thermo

  }  // namespace ELEMENTS
}  // namespace DRT


/*----------------------------------------------------------------------*/
#endif  // SO3_THERMO_H

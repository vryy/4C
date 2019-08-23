/*-----------------------------------------------------------*/
/*! \file
\brief discsh3 element

\level 3

\maintainer Christoph Meier
*/
/*-----------------------------------------------------------*/
// header file only included if not yet included!
#ifndef DISCSH3_H
#define DISCSH3_H


#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "../linalg/linalg_serialdensevector.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_discret_faces.H"

#include "Sacado.hpp"
typedef Sacado::Fad::DFad<double> FAD;


// forward declarations
struct _MATERIAL;
struct _DISCSH3_DATA;
// Several parameters which are fixed for Solid Hex8
const int NUMNOD_DISCSH3 = 3;  ///< number of nodes
const int NODDOF_DISCSH3 = 3;  ///< number of dofs per node
const int NUMDOF_DISCSH3 = 9;  ///< total dofs per element
const int NUMGPT_DISCSH3 = 3;  ///< total gpts per element
const int NUMDIM_DISCSH3 = 3;  ///< number of dimensions

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    // forward declarations
    class DiscSh3Line;

    class DiscSh3Type : public DRT::ElementType
    {
     public:
      std::string Name() const { return "DiscSh3Type"; }

      static DiscSh3Type& Instance();

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
      static DiscSh3Type instance_;

      //  std::map<int,LINALG::Matrix<3,1> > MapToNeighbours_;
    };

    /*!
    \brief A C++ wrapper for the DiscSh3 element

    \author mukherjee (mukherjee@lnm.mw.tum.de)
    */
    class DiscSh3 : public DRT::Element
    {
     public:
      //! @name Friends
      friend class DiscSh3Type;
      friend class DiscSh3Line;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      */
      DiscSh3(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      DiscSh3(const DiscSh3& old);

      /*!
      \brief Deep copy this instance of DiscSh3 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const;

      /*!
      \brief Get shape type of element
      */
      virtual DiscretizationType Shape() const;

      /*!
      \brief Return number of lines (Edges) of this element
      */
      virtual int NumLine() const
      {
        return 3;
      }  // Currenty implemented only for triangular elements

      /*!
      \brief Return number of surfaces of this element
      */
      virtual int NumSurface() const { return 1; }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element

      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Lines();

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element

      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Surfaces();

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      virtual int UniqueParObjectId() const { return DiscSh3Type::Instance().UniqueParObjectId(); }

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
      virtual ~DiscSh3();

      //@}

      //! @name Acess methods


      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual DRT::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      virtual int NumDofPerNode(const DRT::Node& node) const { return NODDOF_DISCSH3; }

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

      virtual DRT::ElementType& ElementType() const { return DiscSh3Type::Instance(); }

      /*!
      \brief Query names of element data to be visualized using BINIO

      The element fills the provided map with key names of
      visualization data the element wants to visualize AT THE CENTER
      of the element geometry. The values is supposed to be dimension of the
      data to be visualized. It can either be 1 (scalar), 3 (vector), 6 (sym. tensor)
      or 9 (nonsym. tensor)

      Example:
      \code
        // Name of data is 'Owner', dimension is 1 (scalar value)
        names.insert(std::pair<string,int>("Owner",1));
        // Name of data is 'StressesXYZ', dimension is 6 (sym. tensor value)
        names.insert(std::pair<string,int>("StressesXYZ",6));
      \endcode

      \param names (out): On return, the derived class has filled names with
                          key names of data it wants to visualize and with int dimensions
                          of that data.
      */
      //  virtual void VisNames(std::map<std::string,int>& names);

      /*!
      \brief Query data to be visualized using BINIO of a given name

      The method is supposed to call this base method to visualize the owner of
      the element.
      If the derived method recognizes a supported data name, it shall fill it
      with corresponding data.
      If it does NOT recognizes the name, it shall do nothing.

      \warning The method must not change size of data

      \param name (in):   Name of data that is currently processed for visualization
      \param data (out):  data to be filled by element if element recognizes the name
      */
      //  virtual bool VisData(const std::string& name, std::vector<double>& data);

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef);

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate shell8 element stiffness, mass, internal forces etc

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
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3);


      int EvaluateCurvatureEnergy(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
          Epetra_SerialDenseVector& elevec3);

      virtual int HowManyRandomNumbersINeed();

      virtual int HowManyRandomNumbersPerNode();

      void EvaluatePTC(Teuchos::ParameterList& params, Epetra_SerialDenseMatrix& elemat1);
      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a surfaces Neumann condition on the shell element

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


      /*!
      \brief Get Teuchos::RCP to the internal face adjacent to this element as master element and
      the parent_slave element
      */
      virtual Teuchos::RCP<DRT::Element> CreateFaceElement(
          DRT::Element* parent_slave,  //!< parent slave fluid3 element
          int nnode,                   //!< number of surface nodes
          const int* nodeids,          //!< node ids of surface element
          DRT::Node** nodes,           //!< nodes of surface element
          const int lsurface_master,   //!< local surface number w.r.t master parent element
          const int lsurface_slave,    //!< local surface number w.r.t slave parent element
          const std::vector<int>& localtrafomap  //! local trafo map
      );


      // calculate surface area of element (private)
      double CalcSurfArea(
          DRT::Discretization& discretization, const Epetra_Vector& discol, bool refconfig);

      // calculate volume enclosed of element (private)
      double CalcVolume(
          DRT::Discretization& discretization, const Epetra_Vector& discol, bool refconfig);

      //@}

      //! @name Other

      /*!
      \brief Get ptr to nodal thicknesses
      */
      const std::vector<double>* GetThickness() { return data_.Get<std::vector<double>>("thick"); }

      /*!
      \brief Get reference to directors
      */
      const Epetra_SerialDenseMatrix* GetDirectors()
      {
        return data_.Get<Epetra_SerialDenseMatrix>("a3ref");
      }

      //@}


     private:
      //! action parameters recognized by shell8
      enum ActionType
      {
        none,
        calc_struct_linstiff,
        calc_struct_nlnstiff,
        calc_struct_internalforce,
        calc_struct_linstiffmass,
        calc_struct_nlnstiffmass,
        calc_struct_nlnstifflmass,
        calc_struct_ptcstiff,
        calc_struct_stress,
        calc_struct_eleload,
        calc_struct_fsiload,
        calc_struct_refvol,
        calc_struct_currvol,
        calc_struct_refCG,
        calc_struct_currCG,
        calc_struct_refarea,
        calc_struct_currarea,
        calc_struct_update_istep,
        calc_struct_reset_istep  //!< reset elementwise internal variables
                                 //!< during iteration to last converged state
      };

      //! type of forces to be calculated in postprocessing
      enum ForceType
      {
        discsh3_none,
        discsh3_xyz,
        discsh3_rst,
        discsh3_rst_ortho
      };

      //! type of force calculation
      ForceType forcetype_;
      //! shell thickness
      double thickness_;
      //! shell thickness
      LINALG::Matrix<12, 12, FAD> stiff_dummy3nodes_;
      //! no. gaussian points in each spatial dimension
      int ngp_[3];
      //! no. gaussian points for triangle elements
      int ngptri_;
      //! no. of internal eas dofs
      int nhyb_;
      //! eas parameters
      int eas_[5];
      //! type of ANS: 0=none 1=Q 2=T 3=QT
      int ans_;
      //! scaled director conditioning
      double sdc_;
      //! number of the material law
      int material_;

      // Nodal positions at current time step
      LINALG::Matrix<1, NUMDOF_DISCSH3> x_n_;

      // Nodal positions at previous time step
      LINALG::Matrix<1, NUMDOF_DISCSH3> x_n_1_;

      DRT::UTILS::GaussRule2D gaussrule_;  ///< Gaussrule

      //! container containing all kind of discsh3 thisandthat
      DRT::Container data_;

      //! contact discretization (basically a copy)
      Teuchos::RCP<DRT::Discretization> discret_;

      bool CalcDampingForces_;

      // internal calculation methods

      // don't want = operator
      DiscSh3& operator=(const DiscSh3& old);

      //  evaluate the element (private)
      void sh3_nlnstiffmass(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, const std::vector<double>& vel, const std::vector<double>& disp,
          Epetra_SerialDenseMatrix* stiffmatrix, Epetra_SerialDenseMatrix* massmatrix,
          Epetra_SerialDenseVector* force);

      inline void CalcBrownian(Teuchos::ParameterList& params,
          const std::vector<double>& vel,         //!< element velocity vector
          const std::vector<double>& disp,        //!< element displacement vector
          Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
          Epetra_SerialDenseVector* force);

      //  evaluate the stiffness contribution from Area constraint (private)
      void CGConstrtStiffmass(Teuchos::ParameterList& params, const std::vector<double>& disp,
          const std::vector<double>& vel, Epetra_SerialDenseMatrix* stiffmatrix,
          Epetra_SerialDenseVector* force, const int NumGlobalNodes);

      //  evaluate the stiffness contribution from Area constraint (private)
      void CGGlobalConstrtStiffmass(Teuchos::ParameterList& params, const std::vector<double>& disp,
          const std::vector<double>& vel, Epetra_SerialDenseMatrix* stiffmatrix,
          Epetra_SerialDenseVector* force, const int NumGlobalNodes);

      //  evaluate the stiffness contribution from Area constraint (private)
      void BaryConstrtStiffmass(Teuchos::ParameterList& params, const std::vector<double>& disp,
          const std::vector<double>& vel, Epetra_SerialDenseMatrix* stiffmatrix,
          Epetra_SerialDenseVector* force, const int NumGlobalNodes);

      //  evaluate the stiffness contribution from Area constraint (private)
      void AreaConstrtStiffmass(Teuchos::ParameterList& params, const std::vector<double>& disp,
          const std::vector<double>& vel, Epetra_SerialDenseMatrix* stiffmatrix,
          Epetra_SerialDenseVector* force);

      //  evaluate the stiffness contribution from Area constraint (private)
      void AreaConstrtQuadStiffmass(Teuchos::ParameterList& params, const std::vector<double>& disp,
          const std::vector<double>& vel, Epetra_SerialDenseMatrix* stiffmatrix,
          Epetra_SerialDenseVector* force);

      //  evaluate the stiffness contribution from vol constraint (private)
      void VolConstrtStiffmass(Teuchos::ParameterList& params, const std::vector<double>& disp,
          Epetra_SerialDenseMatrix* stiffmatrix, Epetra_SerialDenseVector* force);

      //  evaluate the stiffness contribution from vol constraint (private)
      void VolConstrtGlobalStiff(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const std::vector<double>& disp,
          Epetra_SerialDenseMatrix* stiffmatrix, Epetra_SerialDenseVector* force);

      //  evaluate the stiffness contribution from vol constraint (private)
      void AreaConstrtGlobalStiff(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const std::vector<double>& disp,
          Epetra_SerialDenseMatrix* stiffmatrix, Epetra_SerialDenseVector* force);

      /*!
      \brief Submethod to compute interfacial area and its first and second
             derivatives w.r.t. the displacements (this is needed for surface energy problems)
      */
      void ComputeAreaDeriv(const LINALG::SerialDenseMatrix& x,  ///< spatial configuration
          const int numnode,                                     ///< number of nodes
          const int ndof,                                        ///< number of degrees of freedom
          double& A,                                             ///< area
          Teuchos::RCP<Epetra_SerialDenseVector> Adiff,          ///< first derivative
          Teuchos::RCP<Epetra_SerialDenseMatrix> Adiff2          ///< second derivative
      );

      void ComputeAreaRef(
          const LINALG::SerialDenseMatrix& x0, const int numnode, const int ndof, double& A);

      /*!

      this method evaluates normal and detA at gaussian point

      \param detA   (out) : are at gaussian point
      \param normal (out) : The normal at gaussian point, length is detA!
      \param x      (in)  : nodal coords in either material or spatial frame
      \param deriv  (in)  : derivatives of shape functions
      */
      void SurfaceIntegration(double& detA, std::vector<double>& normal,
          const Epetra_SerialDenseMatrix& x, const Epetra_SerialDenseMatrix& deriv);


      // calculate behavior function arising from area constraint of element (private)
      void CalcBehaviorFunctArea(
          Teuchos::ParameterList& params, const std::vector<double>& disp, FAD& EnergyAreaConstrt);


      // computes translational damping forces and stiffness (private)
      inline void MyDampingForces(Teuchos::ParameterList& params, const std::vector<double>& vel,
          const std::vector<double>& disp, Epetra_SerialDenseMatrix* stiffmatrix,
          Epetra_SerialDenseVector* force);


      // computes translational damping forces and stiffness (private)
      inline void MyTranslationalDamping(Teuchos::ParameterList& params,
          const std::vector<double>& vel, const std::vector<double>& disp,
          Epetra_SerialDenseMatrix* stiffmatrix, Epetra_SerialDenseVector* force);


      // computes translational damping forces  (private)
      void MyDampingForcesArea(Teuchos::ParameterList& params,  //!< parameter list
          const std::vector<double>& vel,   //!< vector containing first order time derivative of
                                            //!< nodal positions and nodal tangents of an element
          const std::vector<double>& disp,  //!< vector containing change in nodal positions and
                                            //!< nodal tangents of an element w.r.t. inital config.
          Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
          Epetra_SerialDenseVector* force);       //!< element internal force vector


      void MyDampingForcesVol(Teuchos::ParameterList& params,  //!< parameter list
          const std::vector<double>& vel,   //!< vector containing first order time derivative of
                                            //!< nodal positions and nodal tangents of an element
          const std::vector<double>& disp,  //!< vector containing change in nodal positions and
                                            //!< nodal tangents of an element w.r.t. inital config.
          Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
          Epetra_SerialDenseVector* force);       //!< element internal force vector

      inline void MyStochasticForces(Teuchos::ParameterList& params,  //!< parameter list
          const std::vector<double>& vel,                             //!< element velocity vector
          const std::vector<double>& disp,                            //!< element disp vector
          Epetra_SerialDenseMatrix* stiffmatrix,                      //!< element stiffness matrix
          Epetra_SerialDenseVector* force);  //!< element internal force vector


      // calculate behavior function arising from area constraint of element (private)
      void CalcBehaviorFunctVol(
          Teuchos::ParameterList& params, const std::vector<double>& disp, FAD& BehaviorFunctVol);

      // calculate surface area of element (private)
      FAD CalcSurfaceArea(const std::vector<double>& disp, bool refconfig);

      // calculate surface area of element (private)
      FAD CalcSurfaceAreaFAD(const std::vector<double>& disp, bool refconfig);

      // calculate surface area of element (private)
      FAD CalcSurfaceArea(DRT::Discretization& discretization, bool refconfig);


      // calculate surface area of triangle (private)
      FAD CalcSurfaceArea(LINALG::Matrix<1, 3, FAD>& vertex1, LINALG::Matrix<1, 3, FAD>& vertex2,
          LINALG::Matrix<1, 3, FAD>& vertex3);

      double CalcSurfaceArea(LINALG::Matrix<1, 3>& vertex1, LINALG::Matrix<1, 3>& vertex2,
          LINALG::Matrix<1, 3>& vertex3);

      // calculate volume enclosed of element (private)
      void CalcVolume(FAD& volume, const std::vector<double>& disp, bool refconfig);

      // calculate surface area of element (private)
      LINALG::Matrix<1, 3, FAD> CalcSurfaceNormal(std::vector<FAD>& x_FAD);

      void CheckIfOutwardsNormal(Teuchos::ParameterList& params, const int NumGElements);

      // Calculate cross product of 2 FAD vectors
      inline LINALG::Matrix<1, 3, FAD> CalcCrossProduct(
          LINALG::Matrix<1, 3, FAD>& vector1, LINALG::Matrix<1, 3, FAD>& vector2) const
      {
        LINALG::Matrix<1, 3, FAD> crossprod(true);
        // Cross Product
        crossprod(0) = vector1(1) * vector2(2) - vector1(2) * vector2(1);
        crossprod(1) = vector1(2) * vector2(0) - vector1(0) * vector2(2);
        crossprod(2) = vector1(0) * vector2(1) - vector1(1) * vector2(0);

        return crossprod;
      }

      /*!
       \brief Create matrix with material configuration

       \param x  (out)  : nodal coords in material frame
       */
      inline LINALG::Matrix<1, 9> MaterialConfiguration() const
      {
        LINALG::Matrix<1, 9> x(true);
        const int numnode = NumNode();
        //     std::cout<<"numnode="<<numnode<<std::endl;
        for (int i = 0; i < numnode; ++i)
        {
          //       std::cout<<"Material"<<std::endl;
          //       ((Nodes())[i])->Print(std::cout); std::cout<<std::endl;
          x(3 * i) = Nodes()[i]->X()[0];
          x(3 * i + 1) = Nodes()[i]->X()[1];
          x(3 * i + 2) = Nodes()[i]->X()[2];
        }
        return x;
      }

      /*!
       \brief Create matrix with material configuration

       \param x  (out)  : nodal coords in material frame
       */
      inline void MaterialConfiguration(LINALG::SerialDenseMatrix& x) const
      {
        const int numnode = NumNode();
        for (int i = 0; i < numnode; ++i)
        {
          x(i, 0) = Nodes()[i]->X()[0];
          x(i, 1) = Nodes()[i]->X()[1];
          x(i, 2) = Nodes()[i]->X()[2];
        }
        return;
      }

      /*! \brief Check if viscous forces need to be calculated
       *  Relevant for edge based integration
       */
      bool IfHaveDamping() { return CalcDampingForces_; }

      /*!
      \brief Create matrix with spatial configuration

      \param x     (out)  : nodal coords in spatial frame
      \param disp  (int)  : displacements
      */
      LINALG::Matrix<1, 9> SpatialConfiguration(const std::vector<double>& disp) const;

      /*!
      \brief Create matrix with spatial configuration

      \param x     (out)  : nodal coords in spatial frame
      \param disp  (int)  : displacements
      */
      inline void SpatialConfiguration(
          LINALG::SerialDenseMatrix& x, const std::vector<double> disp) const
      {
        const int numnode = NumNode();
        for (int i = 0; i < numnode; ++i)
        {
          x(i, 0) = Nodes()[i]->X()[0] + disp[i * 3 + 0];
          x(i, 1) = Nodes()[i]->X()[1] + disp[i * 3 + 1];
          x(i, 2) = Nodes()[i]->X()[2] + disp[i * 3 + 2];
        }
        return;
      }

      LINALG::Matrix<1, 9> SpatialConfiguration(DRT::Discretization& dis) const;

      LINALG::Matrix<1, 9> SpatialConfiguration(
          DRT::Discretization& dis, const Epetra_Vector& discol) const;

      LINALG::Matrix<1, 9> GetVel(DRT::Discretization& dis) const;

      void AddPrimaryDOFsMaster(DRT::Element& master, DRT::Element& neighbour,
          std::vector<int>& connectivity, std::vector<FAD>& x_FAD, DRT::Discretization& dis,
          bool refconfig);

      void AddPrimaryDOFsNeighbour(DRT::Element& master, DRT::Element& neighbour,
          std::vector<FAD>& x_FAD, DRT::Discretization& dis, bool refconfig);


      void SortPrimaryDOFs(DRT::Element& master, DRT::Element& neighbour, std::vector<FAD>& x_FAD,
          std::vector<FAD>& x_FAD_master, std::vector<FAD>& x_FAD_neighbour);

      inline void MyDampingConstants(Teuchos::ParameterList& params, LINALG::Matrix<3, 1>& gamma);


      /*!
      \brief Create matrix with spatial configuration

      \param x     (out)  : nodal coords in spatial frame
      \param disp  (int)  : displacements
      */


      //! lump mass matrix (bborn 07/08)
      void sh3_lumpmass(const std::vector<double>& disp,
          Epetra_SerialDenseMatrix* massmatrix);  //!< element mass matrix

      FAD CalcTheta(LINALG::Matrix<1, 3, FAD>& vector1, LINALG::Matrix<1, 3, FAD>& vector2);


    };  // class DiscSh3


    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================


    class DiscSh3LineType : public DRT::ElementType
    {
     public:
      std::string Name() const { return "DiscSh3LineType"; }

      static DiscSh3LineType& Instance();

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
      {
      }

      virtual void ComputeNullSpace(
          DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
      {
      }

     private:
      static DiscSh3LineType instance_;
    };

    /*!
    \brief An element representing a line edge of a shell8 element

    \note This is a pure Neumann boundary condition element. It's only
          purpose is to evaluate line Neumann boundary conditions that might be
          adjacent to a parent shell8 element. It therefore does not implement
          the DRT::Element::Evaluate method and does not have its own ElementRegister class.

    \author gee (gee@lnm.mw.tum.de)
    */
    class DiscSh3Line : public DRT::FaceElement
    {
     public:
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner: Processor owning this line
      \param nnode: Number of nodes attached to this element
      \param nodeids: global ids of nodes attached to this element
      \param nodes: the discretizations map of nodes to build ptrs to nodes from
      \param parent: The parent shell element of this line
      \param lline: the local line number of this line w.r.t. the parent element
      */
      DiscSh3Line(int id, int owner, int nnode, const int* nodeids, DRT::Node** nodes,
          DRT::ELEMENTS::DiscSh3* parent, const int lline);

      // overloading of Constructor
      DiscSh3Line(int id,      ///< element id
          int owner,           ///< owner (= owner of parent element with smallest gid)
          int nnode,           ///< number of nodes
          const int* nodeids,  ///< node ids
          DRT::Node** nodes,   ///< nodes of surface
          DRT::ELEMENTS::DiscSh3* parent_master,  ///< master parent element
          DRT::ELEMENTS::DiscSh3* parent_slave,   ///< slave parent element
          const int lsurface_master,  ///< local surface index with respect to master parent element
          const int lsurface_slave,   ///< local surface index with respect to slave parent element
          const std::vector<int>
              localtrafomap  ///< get the transformation map between the local coordinate systems of
                             ///< the face w.r.t the master parent element's face's coordinate
                             ///< system and the slave element's face's coordinate system
      );

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      DiscSh3Line(const DiscSh3Line& old);

      /*!
      \brief Deep copy this instance of an element and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of the drt_parobject.H file.
      */
      virtual int UniqueParObjectId() const
      {
        return DiscSh3LineType::Instance().UniqueParObjectId();
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
      virtual ~DiscSh3Line();

      //@}

      //! @name Acess methods

      /*!
      \brief Get shape type of element
      */
      virtual DiscretizationType Shape() const;

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
      virtual int NumDofPerElement() const { return 9; }

      /*!
       * \brief Return pointer to the parent element
       */
      virtual DRT::ELEMENTS::DiscSh3* ParentElement() const
      {
        DRT::Element* parent = this->DRT::FaceElement::ParentElement();
        // make sure the static cast below is really valid
        dsassert(dynamic_cast<DRT::ELEMENTS::DiscSh3*>(parent) != NULL,
            "Parent element is no shell element");
        return static_cast<DRT::ELEMENTS::DiscSh3*>(parent);
      }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      // Calculate cross product of 2 FAD vectors
      inline LINALG::Matrix<1, 3, FAD> CalcCrossProduct(
          LINALG::Matrix<1, 3, FAD>& vector1, LINALG::Matrix<1, 3, FAD>& vector2) const
      {
        LINALG::Matrix<1, 3, FAD> crossprod(true);
        // Cross Product
        crossprod(0) = vector1(1) * vector2(2) - vector1(2) * vector2(1);
        crossprod(1) = vector1(2) * vector2(0) - vector1(0) * vector2(2);
        crossprod(2) = vector1(0) * vector2(1) - vector1(1) * vector2(0);

        return crossprod;
      }

      inline LINALG::Matrix<3, 3, FAD> SkewSymmetricTrans(LINALG::Matrix<1, 3, FAD>& v) const
      {
        LINALG::Matrix<3, 3, FAD> w(true);
        w(0, 0) = 0;
        w(0, 1) = -v(2);
        w(0, 2) = v(1);
        w(1, 0) = v(2);
        w(1, 1) = 0;
        w(1, 2) = -v(0);
        w(2, 0) = -v(1);
        w(2, 1) = v(0);
        w(2, 2) = 0;

        return w;
      }

      /*!
       \brief Create matrix with material configuration

       \param x  (out)  : nodal coords in material frame
       */
      inline LINALG::Matrix<1, 9> MaterialConfiguration() const
      {
        LINALG::Matrix<1, 9> x(true);
        const int numnode = NumNode();
        for (int i = 0; i < numnode; ++i)
        {
          x(3 * i) = Nodes()[i]->X()[0];
          x(3 * i + 1) = Nodes()[i]->X()[1];
          x(3 * i + 2) = Nodes()[i]->X()[2];
        }
        return x;
      }


      void AddPrimaryDOFsMaster(DRT::ELEMENTS::DiscSh3& master, DRT::ELEMENTS::DiscSh3& neighbour,
          std::vector<int>& connectivity, std::vector<FAD>& x_FAD, DRT::Discretization& dis,
          bool refconfig);

      void AddPrimaryDOFsSlave(DRT::ELEMENTS::DiscSh3& master, DRT::ELEMENTS::DiscSh3& slave,
          std::vector<FAD>& x_FAD, DRT::Discretization& dis, bool refconfig);

      void AddCurrVel(DRT::ELEMENTS::DiscSh3& master, DRT::ELEMENTS::DiscSh3& slave,
          std::vector<int>& connectivity, std::vector<FAD>& v_FAD, DRT::Discretization& dis);


      void SortPrimaryDOFs(DRT::Element& master, DRT::Element& neighbour, std::vector<FAD>& x_FAD,
          std::vector<FAD>& x_FAD_master, std::vector<FAD>& x_FAD_neighbour);

      void ReassembleMATBlock(const int row_block,  ///< row block
          const int col_block,                      ///< column block
          Epetra_SerialDenseMatrix& mat_block,      ///< matrix block
          LINALG::Matrix<NUMDOF_DISCSH3, NUMDOF_DISCSH3>&
              elematrix_mm,  ///< element matrix master-master block
          LINALG::Matrix<NUMDOF_DISCSH3, NUMDOF_DISCSH3>&
              elematrix_ms,  ///< element matrix master-slave block
          LINALG::Matrix<NUMDOF_DISCSH3, NUMDOF_DISCSH3>&
              elematrix_sm,  ///< element matrix slave-master block
          LINALG::Matrix<NUMDOF_DISCSH3, NUMDOF_DISCSH3>&
              elematrix_ss,  ///< element matrix slave-slave block
          std::vector<int>&
              lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
          std::vector<int>&
              lm_slaveNodeToPatch  ///< local map between slave nodes and nodes in patch
      );

      void ReassembleRHSBlock(const int row_block,         ///< row block
          Epetra_SerialDenseVector& rhs_block,             ///< rhs block
          LINALG::Matrix<NUMDOF_DISCSH3, 1>& elevector_m,  ///< element vector master block
          LINALG::Matrix<NUMDOF_DISCSH3, 1>& elevector_s,  ///< element vector slave block
          std::vector<int>&
              lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
          std::vector<int>&
              lm_slaveNodeToPatch  ///< local map between slave nodes and nodes in patch
      );


      // Calculate gradient of normal of a surface
      void CalcGradienNormal(std::vector<FAD>& x_FAD, LINALG::Matrix<3, 3, FAD>& DnDx1,
          LINALG::Matrix<3, 3, FAD>& DnDx2, LINALG::Matrix<3, 3, FAD>& DnDx3);

      void LengthConstrtStiffmass(Teuchos::ParameterList& params,
          DRT::ELEMENTS::DiscSh3Line* edge,  ///< internal face element
          std::vector<FAD>& x_FAD_curr,
          DRT::Discretization& discretization,  ///< discretization
          LINALG::Matrix<1, 12, FAD>& ana_force_aux, LINALG::Matrix<12, 12, FAD>& stiff_dummy);

      // Spatial Configuration
      LINALG::Matrix<1, 6> SpatialConfiguration(DRT::Discretization& dis) const;

      LINALG::Matrix<1, 6> SpatialConfiguration(
          DRT::Discretization& dis, const Epetra_Vector& discol) const;


      FAD CalcTheta(LINALG::Matrix<1, 3, FAD>& vector1, LINALG::Matrix<1, 3, FAD>& vector2);

      FAD GetRefEdgeLength();

      FAD GetEdgeLengthPrevTimeStep();

      FAD GetCurrEdgeLength(DRT::Discretization& dis);

      double GetCurrEdgeLength(DRT::Discretization& dis, const Epetra_Vector& discol);

      FAD GetCurrEdgeLength(std::vector<FAD>& x_FAD_curr);

      // calculate surface area of element (private)
      LINALG::Matrix<1, 3, FAD> CalcSurfaceNormalMaster(std::vector<FAD>& x_FAD);

      LINALG::Matrix<1, 3, FAD> CalcSurfaceNormalSlave(std::vector<FAD>& x_FAD);


      virtual DRT::ElementType& ElementType() const { return DiscSh3LineType::Instance(); }

      //@}

      //! @name Evaluate methods

      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a line Neumann condition on the shell element

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

      virtual int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseMatrix* elemat1 = NULL);

      //    evaluate implementation for internal surface stabilization
      virtual int EvaluateEdges(const Teuchos::ParameterList& params,
          DRT::ELEMENTS::DiscSh3Line* edge,     ///< internal face element
          DRT::Discretization& discretization,  ///< discretization
          std::vector<int>& patchlm,            ///< patch local map
          std::vector<int>& lm_masterToPatch,   ///< local map between master dofs and patchlm
          std::vector<int>& lm_slaveToPatch,    ///< local map between slave dofs and patchlm
          std::vector<int>& lm_faceToPatch,     ///< local map between face dofs and patchlm
          std::vector<int>&
              lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
          std::vector<int>&
              lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
          std::vector<Epetra_SerialDenseMatrix>& elemat_blocks,  ///< element matrix blocks
          std::vector<Epetra_SerialDenseVector>& elevec_blocks   ///< element vector blocks
      );


      void AssembleInternalFacesUsingNeighborData(const Teuchos::ParameterList& params,
          DRT::ELEMENTS::DiscSh3Line* intface,              ///< internal face element
          Teuchos::RCP<MAT::Material>& material,            ///< material for face stabilization
          std::vector<int>& nds_master,                     ///< nodal dofset w.r.t. master element
          std::vector<int>& nds_slave,                      ///< nodal dofset w.r.t. slave element
          DRT::DiscretizationFaces& discretization,         ///< faces discretization
          Teuchos::RCP<LINALG::SparseMatrix> systemmatrix,  ///< systemmatrix
          Teuchos::RCP<Epetra_Vector> systemvector          ///< systemvector
      );


      void PatchLocationVector(DRT::Discretization& discretization,  ///< discretization
          std::vector<int>& nds_master,        ///< nodal dofset w.r.t master parent element
          std::vector<int>& nds_slave,         ///< nodal dofset w.r.t slave parent element
          std::vector<int>& patchlm,           ///< local map for gdof ids for patch of elements
          std::vector<int>& lm_masterToPatch,  ///< local map between lm_master and lm_patch
          std::vector<int>& lm_slaveToPatch,   ///< local map between lm_slave and lm_patch
          std::vector<int>& lm_faceToPatch,    ///< local map between lm_face and lm_patch
          std::vector<int>&
              lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
          std::vector<int>&
              lm_slaveNodeToPatch  ///< local map between slave nodes and nodes in patch
      );

      void RemapForceVectoMasterSlave(
          LINALG::Matrix<1, 12, FAD>& force,  ///< nodal dofset w.r.t master parent element
          LINALG::Matrix<NUMDOF_DISCSH3, 1>& elevector_m,
          LINALG::Matrix<NUMDOF_DISCSH3, 1>& elevector_s, std::vector<int>& connectivity,
          std::vector<int>& lm_slaveNodeToPatch);

      void RemapStiffmattoMasterSlave(
          LINALG::Matrix<12, 12, FAD>& stiffmatrix,  ///< nodal dofset w.r.t master parent element
          LINALG::Matrix<NUMDOF_DISCSH3, NUMDOF_DISCSH3>& elematrix_mm,
          LINALG::Matrix<NUMDOF_DISCSH3, NUMDOF_DISCSH3>& elematrix_ms,
          LINALG::Matrix<NUMDOF_DISCSH3, NUMDOF_DISCSH3>& elematrix_sm,
          LINALG::Matrix<NUMDOF_DISCSH3, NUMDOF_DISCSH3>& elematrix_ss,
          std::vector<int>&
              lm_masterNodeToPatch,  ///< local map between master nodes and nodes in patch
          std::vector<int>&
              lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
          std::vector<int>& connectivity);
      //@}

     private:
      // don't want = operator
      DiscSh3Line& operator=(const DiscSh3Line& old);

    };  // class DiscSh3Line



  }  // namespace ELEMENTS
}  // namespace DRT


#endif  // #ifndef DISCSH3_H

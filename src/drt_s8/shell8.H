/*----------------------------------------------------------------------------*/
/*!
\brief shell8

\level 1

\maintainer Christoph Meier

*/
/*---------------------------------------------------------------------------*/
#ifndef SHELL8_H
#define SHELL8_H


#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_node.H"

#ifndef MAXNOD_SHELL8
#define MAXNOD_SHELL8 (9) /* maximum of element nodes : max quad9 */
#endif

// forward declarations
struct _MATERIAL;
struct _S8_DATA;

namespace STR
{
  namespace ELEMENTS
  {
    class ParamsInterface;
  }  // namespace ELEMENTS
}  // namespace STR
namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    // forward declarations
    class Shell8Line;

    class Shell8Type : public DRT::ElementType
    {
     public:
      std::string Name() const { return "Shell8Type"; }

      static Shell8Type& Instance();

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
      static Shell8Type instance_;
    };

    /*!
    \brief A C++ wrapper for the shell8 element

    \author gee (gee@lnm.mw.tum.de)
    */
    class Shell8 : public DRT::Element
    {
     public:
      //! @name Friends
      friend class Shell8Type;
      friend class Shell8Line;

      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      */
      Shell8(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      Shell8(const Shell8& old);

      /*!
      \brief Deep copy this instance of Shell8 and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const;

      /*!
      \brief Get shape type of element
      */
      virtual DiscretizationType Shape() const;

      /*!
      \brief Return number of lines of this element
      */
      virtual int NumLine() const
      {
        if (NumNode() == 4 || NumNode() == 9)
          return 4;
        else
          return 3;
      }

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
      virtual int UniqueParObjectId() const { return Shell8Type::Instance().UniqueParObjectId(); }

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
      virtual ~Shell8();

      //@}

      //! @name Access methods

      /*!
      \brief Does this element use EAS?

      ToDo: This function can be declared as pure virtual and each concrete derived
            class has to implement this function. This can be done during the up-coming
            cleaning procedure.                                      hiermeier 09/15
      */
      virtual bool HaveEAS() const { return (nhyb_ > 0); };

      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual DRT::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      virtual int NumDofPerNode(const DRT::Node& node) const { return 6; }

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

      virtual DRT::ElementType& ElementType() const { return Shell8Type::Instance(); }

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
      virtual void VisNames(std::map<std::string, int>& names);

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
      virtual bool VisData(const std::string& name, std::vector<double>& data);

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

      /** \brief set the parameter interface ptr for the solid elements
       *
       *  \param p (in): Parameter list coming from the time integrator.
       *
       *  \author hiermeier
       *  \date 04/16 */
      void SetParamsInterfacePtr(const Teuchos::ParameterList& p);

      /** \brief returns true if the parameter interface is defined and initialized, otherwise false
       *
       *  \author hiermeier
       *  \date 04/16 */
      inline bool IsParamsInterface() const { return (not interface_ptr_.is_null()); }

      /** \brief get access to the parameter interface
       *
       *  \author hiermeier
       *  \date 04/16 */
      inline STR::ELEMENTS::ParamsInterface& ParamsInterface()
      {
        if (not IsParamsInterface()) dserror("The interface ptr is not set!");
        return *interface_ptr_;
      }

     private:
      //! type of forces to be calculated in postprocessing
      enum ForceType
      {
        s8_none,
        s8_xyz,
        s8_rst,
        s8_rst_ortho
      };
      //! type of force calculation
      ForceType forcetype_;
      //! shell thickness
      double thickness_;
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
      //! container containing all kind of shell8 thisandthat
      DRT::Container data_;
      //! old step length for elementwise internal variables (e.g. eas)
      double old_step_length_;
      //! pointer to the ELEMENTS::ParamsInterface object
      Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_;

      // internal calculation methods

      // don't want = operator
      Shell8& operator=(const Shell8& old);

      void s8_recover(
          std::vector<int>& lm, std::vector<double>& disp, std::vector<double>& residual);

      void s8_nlnstiffmass(std::vector<int>& lm, std::vector<double>& disp,
          std::vector<double>& residual, Epetra_SerialDenseMatrix* stiffmatrix,
          Epetra_SerialDenseMatrix* massmatrix, Epetra_SerialDenseVector* force,
          struct _MATERIAL* material);

      //! lump mass matrix (bborn 07/08)
      void s8_lumpmass(Epetra_SerialDenseMatrix* emass);  //!< element mass matrix

      void s8stress(struct _MATERIAL* material, std::vector<double>& mydisp);

      void s8_surfaceload(Teuchos::ParameterList& params, Epetra_SerialDenseVector& elevec1,
          const std::vector<double>& disp);

      void s8_integration_points(struct _S8_DATA& data);

      void s8BtDB(Epetra_SerialDenseMatrix& estif, const Epetra_SerialDenseMatrix& bop,
          const double D[][12], const int iel, const int numdf, const double weight);

      void s8tvkg(Epetra_SerialDenseMatrix& estif, double stress_r[],
          const std::vector<double>& funct, const Epetra_SerialDenseMatrix& deriv, const int numdf,
          const int iel, const double weight, const double e1, const double e2);

      void s8anstvkg(Epetra_SerialDenseMatrix& estif, double stress_r[],
          const std::vector<double>& funct, const Epetra_SerialDenseMatrix& deriv, const int numdf,
          const int iel, const double weight, const double e1, const double e2, const double frq[],
          const double fsq[], const std::vector<double> funct1q[],
          const std::vector<double> funct2q[], const Epetra_SerialDenseMatrix deriv1q[],
          const Epetra_SerialDenseMatrix deriv2q[], const int ansq, const int nsansq);

      void s8intforce(Epetra_SerialDenseVector& intforce, const double stress_r[],
          const Epetra_SerialDenseMatrix& bop, const int iel, const int numdf, const int nstress_r,
          const double weight);

      void s8tmas(const std::vector<double>& funct, const std::vector<double>& thick,
          Epetra_SerialDenseMatrix& emass, const int iel, const int numdf, const double facv,
          const double facw, const double facvw);

      void s8tmat(struct _MATERIAL* material, double stress[], double strain[], double** C,
          double gmkovc[][3], double gmkonc[][3], double gmkovr[][3], double gmkonr[][3],
          double gkovc[][3], double gkonc[][3], double gkovr[][3], double gkonr[][3],
          const double detc, const double detr, const double e3, const int option,
          const int ngauss);

      void s8tvma(double D[][12], double** C, double stress[], double stress_r[], const double e3,
          const double fact, const double condfac);

      void s8_YpluseqAx(std::vector<double>& y, const Epetra_SerialDenseMatrix& A,
          const std::vector<double>& x, const double factor, const bool init);
      void s8_YpluseqAx(Epetra_SerialDenseVector& y, const Epetra_SerialDenseMatrix& A,
          const std::vector<double>& x, const double factor, const bool init);

      void s8inv3(double a[][3], double* det);

      void s8trans3(double a[][3]);

      void s8unvc(double* enorm, double vec[], const int n);

      void s8matmatdense(Epetra_SerialDenseMatrix& R, const Epetra_SerialDenseMatrix& A,
          const Epetra_SerialDenseMatrix& B, const int ni, const int nk, const int nj,
          const int init, const double factor);
      void s8matmatdense(Epetra_SerialDenseMatrix& R, const double A[][12],
          const Epetra_SerialDenseMatrix& B, const int ni, const int nk, const int nj,
          const int init, const double factor);
      void s8mattrnmatdense(Epetra_SerialDenseMatrix& R, const Epetra_SerialDenseMatrix& A,
          const Epetra_SerialDenseMatrix& B, const int ni, const int nk, const int nj,
          const int init, const double factor);
      void s8mattrnvecdense(std::vector<double>& r, const Epetra_SerialDenseMatrix& A,
          const double b[], const int ni, const int nk, const int init, const double factor);

      double s8_localcoordsofnode(const int node, const int flag, const int numnode) const;

      void s8_shapefunctions(std::vector<double>& funct, Epetra_SerialDenseMatrix& deriv,
          const double r, const double s, const int numnode, const int doderiv) const;

      void s8tvbo(const double e1, const double e2, Epetra_SerialDenseMatrix& bop,
          const std::vector<double>& funct, const Epetra_SerialDenseMatrix& deriv, const int iel,
          const int numdf, const double akov[][3], const double a3kvp[][2], const int nsansq);

      void s8tvmr(const double x[][MAXNOD_SHELL8], const double a3[][MAXNOD_SHELL8],
          double akov[][3], double akon[][3], double amkov[][3], double amkon[][3], double* det,
          const std::vector<double>& funct, const Epetra_SerialDenseMatrix& deriv, const int iel,
          double a3kvp[][2], const int flag);

      void s8tmtr(const double x[][MAXNOD_SHELL8], const double a3[][MAXNOD_SHELL8],
          const double e3, double gkov[][3], double gkon[][3], double gmkov[][3], double gmkon[][3],
          double* det, const std::vector<double>& funct, const Epetra_SerialDenseMatrix& deriv,
          const int iel, const double condfac, const int flag);

      void s8tvhe(double gmkovr[][3], double gmkovc[][3], double gmkonr[][3], double gmkonc[][3],
          double gkovr[][3], double gkovc[][3], double* detr, double* detc, double amkovc[][3],
          double amkovr[][3], double akovc[][3], double akovr[][3], double a3kvpc[][2],
          double a3kvpr[][2], const double e3, const double condfac);

      void s8_jaco(const std::vector<double>& funct, const Epetra_SerialDenseMatrix& deriv,
          const double x[][MAXNOD_SHELL8], double xjm[][3], const std::vector<double>& hte,
          const double a3ref[][MAXNOD_SHELL8], const double e3, const int iel, double* det,
          double* deta);

      void s8_ans_colloquationpoints(const int nsansq, const int iel, const int ans, double xr1[],
          double xs1[], double xr2[], double xs2[], std::vector<double> funct1q[],
          Epetra_SerialDenseMatrix deriv1q[], std::vector<double> funct2q[],
          Epetra_SerialDenseMatrix deriv2q[], const double xrefe[][MAXNOD_SHELL8],
          const double a3r[][MAXNOD_SHELL8], const double xcure[][MAXNOD_SHELL8],
          const double a3c[][MAXNOD_SHELL8], double akovr1q[][3][3], double akonr1q[][3][3],
          double amkovr1q[][3][3], double amkonr1q[][3][3], double a3kvpr1q[][3][2],
          double akovc1q[][3][3], double akonc1q[][3][3], double amkovc1q[][3][3],
          double amkonc1q[][3][3], double a3kvpc1q[][3][2], double akovr2q[][3][3],
          double akonr2q[][3][3], double amkovr2q[][3][3], double amkonr2q[][3][3],
          double a3kvpr2q[][3][2], double akovc2q[][3][3], double akonc2q[][3][3],
          double amkovc2q[][3][3], double amkonc2q[][3][3], double a3kvpc2q[][3][2], double* detr,
          double* detc);

      void s8_ans_colloquationcoords(
          double xqr1[], double xqs1[], double xqr2[], double xqs2[], const int iel, const int ans);

      void s8_ansqshapefunctions(double frq[], double fsq[], const double r, const double s,
          const int iel, const int nsansq);

      void s8ansbbarq(Epetra_SerialDenseMatrix& bop, const double frq[], const double fsq[],
          const std::vector<double> funct1q[], const std::vector<double> funct2q[],
          const Epetra_SerialDenseMatrix deriv1q[], const Epetra_SerialDenseMatrix deriv2q[],
          const double akovc1q[][3][3], const double akovc2q[][3][3], const double a3kvpc1q[][3][2],
          const double a3kvpc2q[][3][2], const int& iel, const int& numdf, const int& nsansq);

      void s8anstvheq(double gmkovr[][3], double gmkovc[][3], double gmkonr[][3],
          double gmkonc[][3], double gkovr[][3], double gkovc[][3], double amkovc[][3],
          double amkovr[][3], double akovc[][3], double akovr[][3], double a3kvpc[][2],
          double a3kvpr[][2], double* detr, double* detc, double amkovr1q[][3][3],
          double amkovc1q[][3][3], double akovr1q[][3][3], double akovc1q[][3][3],
          double a3kvpr1q[][3][2], double a3kvpc1q[][3][2], double amkovr2q[][3][3],
          double amkovc2q[][3][3], double akovr2q[][3][3], double akovc2q[][3][3],
          double a3kvpr2q[][3][2], double a3kvpc2q[][3][2], double frq[], double fsq[],
          const double e3, const int nansq, const int iel, const double condfac);

      void s8eas(const int nhyb, const double e1, const double e2, const int iel, const int* eas,
          Epetra_SerialDenseMatrix& P);

      void s8transeas(Epetra_SerialDenseMatrix& P, Epetra_SerialDenseMatrix& transP,
          Epetra_SerialDenseMatrix& T, const double akovr[][3], const double akonr0[][3],
          const double detr, const double detr0, const int nhyb);

      void s8vthv(double gmkovc[][3], double gmkonc[][3], const std::vector<double>& epsh,
          double* detc, const double e3, const double condfac);

    };  // class Shell8


    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================


    class Shell8LineType : public DRT::ElementType
    {
     public:
      std::string Name() const { return "Shell8LineType"; }

      static Shell8LineType& Instance();

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
      static Shell8LineType instance_;
    };


    /*!
    \brief An element representing a line edge of a shell8 element

    \note This is a pure Neumann boundary condition element. It's only
          purpose is to evaluate line Neumann boundary conditions that might be
          adjacent to a parent shell8 element. It therefore does not implement
          the DRT::Element::Evaluate method and does not have its own ElementRegister class.

    \author gee (gee@lnm.mw.tum.de)
    */
    class Shell8Line : public DRT::FaceElement
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
      Shell8Line(int id, int owner, int nnode, const int* nodeids, DRT::Node** nodes,
          DRT::ELEMENTS::Shell8* parent, const int lline);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      Shell8Line(const Shell8Line& old);

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
        return Shell8LineType::Instance().UniqueParObjectId();
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
      virtual ~Shell8Line();

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
      virtual int NumDofPerNode(const DRT::Node& node) const { return 6; }

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
       * \brief Return pointer to the parent element
       */
      virtual DRT::ELEMENTS::Shell8* ParentElement() const
      {
        DRT::Element* parent = this->DRT::FaceElement::ParentElement();
        // make sure the static cast below is really valid
        dsassert(dynamic_cast<DRT::ELEMENTS::Shell8*>(parent) != NULL,
            "Parent element is no shell element");
        return static_cast<DRT::ELEMENTS::Shell8*>(parent);
      }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      virtual DRT::ElementType& ElementType() const { return Shell8LineType::Instance(); }

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

      //@}

     private:
      // don't want = operator
      Shell8Line& operator=(const Shell8Line& old);

    };  // class Shell8Line



  }  // namespace ELEMENTS
}  // namespace DRT


#endif  // #ifndef SHELL8_H

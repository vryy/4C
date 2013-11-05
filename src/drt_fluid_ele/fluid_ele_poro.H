/*----------------------------------------------------------------------*/
/*!
 \file fluid_ele_poro.H

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/


#ifndef FLUID_ELE_PORO_H_
#define FLUID_ELE_PORO_H_

#include "fluid_ele.H"

namespace DRT
{

class Discretization;

namespace ELEMENTS
{

class FluidPoroEleType : public FluidType
{
public:

  std::string Name() const { return "FluidPoroEleType"; }

  static FluidPoroEleType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

  /// pre-evaluation
  virtual void PreEvaluate(DRT::Discretization&                dis,
                          Teuchos::ParameterList&              p,
                          Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                          Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                          Teuchos::RCP<Epetra_Vector>          systemvector1,
                          Teuchos::RCP<Epetra_Vector>          systemvector2,
                          Teuchos::RCP<Epetra_Vector>          systemvector3);

private:

  static FluidPoroEleType instance_;
};

/*!
\brief A C++ wrapper for the fluid element
*/
class FluidPoro : public Fluid
{
public:

  //@}
  //! @name constructors and destructors and related methods

  /*!
  \brief standard constructor
  */
  FluidPoro(
    int id,    ///< A unique global id
    int owner  ///< number of owning processor
    );

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  FluidPoro(const FluidPoro& old);

  /*!
  \brief Destructor

  */
  virtual ~FluidPoro(){};

  /*!
  \brief Deep copy this instance of fluid and return pointer to the copy

  The Clone() method is used from the virtual base class Element in cases
  where the type of the derived class is unknown and a copy-ctor is needed

  */
  DRT::Element* Clone() const;

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of this file.
  */
  virtual int UniqueParObjectId() const { return FluidPoroEleType::Instance().UniqueParObjectId(); }

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
  \brief Get vector of RCPs to the lines of this element
  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();

  /*!
  \brief Get vector of RCPs to the surfaces of this element
  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Surfaces();

  /*!
  \brief Get vector of RCPs to the volume of this element
  */
  virtual std::vector<Teuchos::RCP<DRT::Element> > Volumes();

  //@}

  //! @name Acess methods

  /*!
  \brief Print this element
  */
  void Print(std::ostream& os) const;

  virtual DRT::ElementType & ElementType() const
  { return FluidPoroEleType::Instance(); }

  /*!
  \brief Get number of degrees of freedom of a certain node in case of multiple
         dofsets, i.e. monolithic solution approaches! (implements pure virtual DRT::Element)

  The element decides how many degrees of freedom its nodes must have.
  As this may vary along a simulation, the element can redecide the
  number of degrees of freedom per node along the way for each of it's nodes
  separately.

  */
  virtual int NumDofPerNode(const unsigned nds, const DRT::Node& node, const std::string disname) const;

  //@}

  //! @name Evaluation

  /*!
  \brief Evaluate an element, that is, call the element routines to evaluate fluid
  element matrices and vectors or evaluate errors, statistics or updates etc. directly.

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
  int Evaluate(Teuchos::ParameterList&   params,
               DRT::Discretization&      discretization,
               std::vector<int>&         lm,
               Epetra_SerialDenseMatrix& elemat1,
               Epetra_SerialDenseMatrix& elemat2,
               Epetra_SerialDenseVector& elevec1,
               Epetra_SerialDenseVector& elevec2,
               Epetra_SerialDenseVector& elevec3);

  //@}

private:

  // don't want = operator
  FluidPoro& operator = (const FluidPoro& old);
};

} // namespace ELEMENTS
} // namespace DRT

#endif /* FLUID_ELE_PORO_H_ */

/*----------------------------------------------------------------------*/
/*!
 \file wall1_scatra.H

 \brief a 2D solid-wall element with ScaTra coupling

 \level 2

 \maintainer Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef WALL1_SCATRA_H_
#define WALL1_SCATRA_H_

#include "wall1.H"

namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{

class Wall1ScatraType : public Wall1Type
{
public:

  std::string Name() const { return "Wall1ScatraType"; }

  static Wall1ScatraType & Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

private:

  static Wall1ScatraType instance_;
};

class Wall1_Scatra : public Wall1
{
public:

  /// @name Constructors and destructors and related methods

  /// Standard Constructor
  Wall1_Scatra(
    int id,  ///< A unique global id
    int owner
  );

  /// Destructor
  virtual ~Wall1_Scatra(){};

  /// Copy Constructor
  ///
  /// Makes a deep copy of a Element
  Wall1_Scatra(const Wall1_Scatra& old);

  /// Deep copy this instance of Wall1 and return pointer to the copy
  ///
  /// The Clone() method is used from the virtual base class Element in cases
  /// where the type of the derived class is unknown and a copy-ctor is needed
  DRT::Element* Clone() const;

  /// Return unique ParObject id
  ///
  /// every class implementing ParObject needs a unique id defined at the
  /// top of this file.
  virtual int UniqueParObjectId() const { return Wall1ScatraType::Instance().UniqueParObjectId(); }

  /// Pack this class so it can be communicated
  ///
  /// \ref Pack and \ref Unpack are used to communicate this element
  virtual void Pack(DRT::PackBuffer& data) const;

  /// Unpack data from a char vector into this class
  ///
  /// \ref Pack and \ref Unpack are used to communicate this element
  virtual void Unpack(const std::vector<char>& data);

  //@}

  virtual int NumDofPerNode(const unsigned nds, const DRT::Node& node, const std::string disname) const
  {
    if (nds!=0)
      return 1;
    return Wall1::NumDofPerNode(node);
  };

  /// Print this element
  void Print(std::ostream& os) const;

  virtual DRT::ElementType & ElementType() const
  { return Wall1ScatraType::Instance(); }

  //@}

  /// @name Input and Creation
  //@{

  /// Read input for this element
//  virtual bool ReadElement(const std::string& eletype,
//                   const std::string& distype,
//                   DRT::INPUT::LineDefinition* linedef);

  //@}

  /// @name Evaluation
  //@{

  void PreEvaluate
    (
      Teuchos::ParameterList& params,  ///< ParameterList for communication between control routine and elements
      DRT::Discretization&      discretization,  ///< pointer to discretization for de-assembly
      DRT::Element::LocationArray& la ///< location array for de-assembly
    );

  /// Evaluate an element
  ///
  /// Evaluate Wall1 element stiffness, mass, internal forces etc
  ///
  /// \return 0 if successful, negative otherwise
  virtual int Evaluate
  (
    Teuchos::ParameterList& params,  ///< ParameterList for communication between control routine and elements
    DRT::Discretization&      discretization,  ///< pointer to discretization for de-assembly
    DRT::Element::LocationArray& la, ///< location array for de-assembly
    Epetra_SerialDenseMatrix& elemat1,  ///< (stiffness-)matrix to be filled by element.
    Epetra_SerialDenseMatrix& elemat2,  ///< (mass-)matrix to be filled by element.
    Epetra_SerialDenseVector& elevec1,  ///< (internal force-)vector to be filled by element
    Epetra_SerialDenseVector& elevec2,  ///< vector to be filled by element
    Epetra_SerialDenseVector& elevec3  ///< vector to be filled by element
  );

protected:

  //! don't want = operator
  Wall1_Scatra& operator = (const Wall1_Scatra& old);

  int MyEvaluate
  (
    Teuchos::ParameterList& params,  ///< ParameterList for communication between control routine and elements
    DRT::Discretization&      discretization,  ///< pointer to discretization for de-assembly
    DRT::Element::LocationArray& la, ///< location array for de-assembly
    Epetra_SerialDenseMatrix& elemat1,  ///< (stiffness-)matrix to be filled by element.
    Epetra_SerialDenseMatrix& elemat2,  ///< (mass-)matrix to be filled by element.
    Epetra_SerialDenseVector& elevec1,  ///< (internal force-)vector to be filled by element
    Epetra_SerialDenseVector& elevec2,  ///< vector to be filled by element
    Epetra_SerialDenseVector& elevec3  ///< vector to be filled by element
  );
};


}
}


#endif /* WALL1_SCATRA_H_ */

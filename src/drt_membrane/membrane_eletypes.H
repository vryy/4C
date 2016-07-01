/*!----------------------------------------------------------------------
\file membrane_eletypes.H

\level 3

<pre>
\maintainer Fabian Bräu
            braeu@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

\brief Nonlinear Membrane Finite Element Type

*----------------------------------------------------------------------*/
#ifndef MEMBRANE_ELETYPES_H
#define MEMBRANE_ELETYPES_H

#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_linedefinition.H"


namespace DRT
{
// forward declarations
class Discretization;

namespace ELEMENTS
{

/*----------------------------------------------------------------------*
 |  TRI 3 Element                                          fbraeu 06/16 |
 *----------------------------------------------------------------------*/
class Membrane_tri3Type : public DRT::ElementType
{
public:

  std::string Name() const { return "Membrane_tri3Type"; }

  static Membrane_tri3Type & Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual int Initialize(DRT::Discretization& dis) {return 0;};

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

private:

  static Membrane_tri3Type instance_;
};

/*----------------------------------------------------------------------*
 |  TRI 6 Element                                          fbraeu 06/16 |
 *----------------------------------------------------------------------*/
class Membrane_tri6Type : public DRT::ElementType
{
public:

  std::string Name() const { return "Membrane_tri6Type"; }

  static Membrane_tri6Type & Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual int Initialize(DRT::Discretization& dis) {return 0;};

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

private:

  static Membrane_tri6Type instance_;
};

/*----------------------------------------------------------------------*
 |  QUAD 4 Element                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
class Membrane_quad4Type : public DRT::ElementType
{
public:

  std::string Name() const { return "Membrane_quad4Type"; }

  static Membrane_quad4Type & Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual int Initialize(DRT::Discretization& dis) {return 0;};

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

private:

  static Membrane_quad4Type instance_;
};

/*----------------------------------------------------------------------*
 |  QUAD 9 Element                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
class Membrane_quad9Type : public DRT::ElementType
{
public:

  std::string Name() const { return "Membrane_quad9Type"; }

  static Membrane_quad9Type & Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual int Initialize(DRT::Discretization& dis) {return 0;};

  virtual void NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

private:

  static Membrane_quad9Type instance_;
};

}
}

#endif  // #ifndef MEMBRANE_ELETYPES_H

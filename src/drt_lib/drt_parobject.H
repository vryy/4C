/*!----------------------------------------------------------------------
\file drt_parobject.H

drt_parser.cpp<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef PAROBJECT_H
#define PAROBJECT_H

#include <vector>
#include <string>
#include "linalg_serialdensematrix.H"
#include <Epetra_SerialDenseVector.h>
#include "linalg_fixedsizematrix.H"

using namespace std;


// These defines are used to have a unique identifier for every class
// that implements the pure virtual ParObject class.
// The unique identifier is stored in a communcation array as an integer
// so the receiving processor knows what type of class has been sent.
// the method DRT::Factory(data) in utils.cpp can then be used to
// create an Unpack the correct instance of ParObject
#define ParObject_Container             (int)0
#define ParObject_Condition             (int)1
#define ParObject_DofSet                (int)2
#define ParObject_Node                  (int)3
#define ParObject_DesignNode            (int)4
#define ParObject_Element               (int)5
#define ParObject_ElementRegister       (int)6
#define ParObject_Shell8                (int)7
#define ParObject_Shell8Register        (int)8
#define ParObject_Shell8Line            (int)9
#define ParObject_Fluid3                (int)10
#define ParObject_Fluid3Register        (int)11
#define ParObject_Fluid3Surface         (int)12
#define ParObject_Wall1                 (int)13
#define ParObject_Wall1Register         (int)14
#define ParObject_Wall1Line             (int)15
#define ParObject_Fluid2                (int)16
#define ParObject_Fluid2Register        (int)17
#define ParObject_Fluid2Line            (int)18
#define ParObject_So_hex8               (int)19
#define ParObject_Soh8Register          (int)20
#define ParObject_Ale3                  (int)24
#define ParObject_Ale3Register          (int)25
#define ParObject_Ale3Surface           (int)26
// missing ints result from removed Condif2 elements
#define ParObject_Sosh8Register         (int)34
#define ParObject_Ale2                  (int)35
#define ParObject_Ale2Register          (int)36
#define ParObject_Ale2Line              (int)37
#define ParObject_So_tet10              (int)38
#define ParObject_Sotet10Register       (int)39
#define ParObject_Ptet                  (int)41
#define ParObject_PtetRegister          (int)42
#define ParObject_NewtonianFluid        (int)43
#define ParObject_MicroMaterial         (int)44
#define ParObject_StVenantKirchhoff     (int)45
#define ParObject_So_sh8                (int)46
#define ParObject_Fluid3Line            (int)47
#define ParObject_LungPenalty           (int)48
#define ParObject_NeoHooke              (int)49
#define ParObject_ConvecDiffus          (int)50
#define ParObject_So_weg6               (int)51
#define ParObject_Sow6Register          (int)52
#define ParObject_SoDisp                (int)56
#define ParObject_SoDispRegister        (int)57
#define ParObject_So_tet4               (int)61
#define ParObject_Sotet4Register        (int)62
#define ParObject_So_hex27              (int)63
#define ParObject_Soh27Register         (int)64
#define ParObject_So_hex20              (int)65
#define ParObject_Soh20Register         (int)66
#define ParObject_So_sh8p8              (int)67
#define ParObject_Sosh8p8Register       (int)68
// missing ints result from removed So_ctet element
#define ParObject_CNode                 (int)71
#define ParObject_CElement              (int)72
#define ParObject_XFluid3               (int)73
#define ParObject_XFluid3Register       (int)74
#define ParObject_XFluid3Surface        (int)75
#define ParObject_XFluid3Line           (int)76
#define ParObject_So_shw6               (int)77
#define ParObject_Soshw6Register        (int)78
#define ParObject_AnisotropicBalzani    (int)79
#define ParObject_Beam2	                (int)80
#define ParObject_DGFluid2              (int)81
#define ParObject_DGFluid2Register      (int)82
#define ParObject_DGFluid2Interface     (int)83
#define ParObject_DGFluid2InterfaceRegister  (int)84
#define ParObject_DGFluid3              (int)85
#define ParObject_DGFluid3Register      (int)86
#define ParObject_DGFluid3Interface     (int)87
#define ParObject_DGFluid3InterfaceRegister  (int)88
#define ParObject_Beam2Register	        (int)89
#define ParObject_Beam3	                (int)90
#define ParObject_Beam3Register	        (int)91
#define ParObject_AAAneohooke           (int)92
#define ParObject_CarreauYasuda         (int)93
#define ParObject_Bele3                 (int)94
#define ParObject_Bele3Register         (int)95
#define ParObject_Bele3Line             (int)96
#define ParObject_ModPowerLaw           (int)97
#define ParObject_ConstraintElement2    (int)98
#define ParObject_ConstraintElement2Register (int)99
#define ParObject_StructuralSurface     (int)100
#define ParObject_StructuralLine        (int)101
#define ParObject_MooneyRivlin          (int)102
#define ParObject_LungOgden             (int)103
#define ParObject_ViscoNeoHooke         (int)104
#define ParObject_Yeoh                  (int)105
#define ParObject_AAAraghavanvorp_damage (int)106
// missing ints result from removed Condif3 elements
#define ParObject_Knotvector            (int)108
#define ParObject_ControlPoint          (int)109
#define ParObject_Fluid2Nurbs           (int)110
#define ParObject_ContChainNetw         (int)111
#define ParObject_Fluid3Nurbs           (int)112
#define ParObject_ArtWallRemod          (int)113
#define ParObject_PreStress             (int)114
#define ParObject_MatList               (int)115
#define ParObject_Truss3                (int)116
#define ParObject_Truss3Register        (int)117
#define ParObject_BioCell               (int)118
#define ParObject_InvDesign             (int)119
#define ParObject_Combust3              (int)120
#define ParObject_Combust3Register      (int)121
#define ParObject_Combust3Surface       (int)122
#define ParObject_Combust3Line          (int)123
#define ParObject_Ion                   (int)124
#define ParObject_ViscoAnisotropic      (int)125
#define ParObject_ConstraintElement3    (int)126
#define ParObject_ConstraintElement3Register (int)127
#define ParObject_Transport             (int)128
#define ParObject_TransportBoundary     (int)129
#define ParObject_SutherlandFluid       (int)130
#define ParObject_SutherlandCondif      (int)131
#define ParObject_CompOgden             (int)132
#define ParObject_ParMaterial           (int)133
#define ParObject_CHARMM                (int)134
#define ParObject_Fluid3NurbsSurface    (int)135
#define ParObject_Wall1Nurbs            (int)136
#define ParObject_DofKey                (int)137
#define ParObject_Ale2_Nurbs            (int)138
#define ParObject_Ale3_Nurbs            (int)139
#define ParObject_So_Nurbs27            (int)140
#define ParObject_Itskov                (int)141
#define ParObject_LogNeoHooke           (int)142
#define ParObject_So_Nurbs27_Surface    (int)143
#define ParObject_So_Nurbs27_Register   (int)144
#define ParObject_Beam2r                (int)145
#define ParObject_Beam2rRegister        (int)146
#define ParObject_Protein               (int)147
#define ParObject_Truss2                (int)148
#define ParObject_Truss2Register        (int)149
#define ParObject_Torsion2              (int)150
#define ParObject_Torsion2Register      (int)151
#define ParObject_Torsion3              (int)152
#define ParObject_Torsion3Register      (int)153
#define ParObject_ElastHyper            (int)154
#define ParObject_Vele3                 (int)155
#define ParObject_Vele3Register         (int)156
#define ParObject_Vele3Surface          (int)157
#define ParObject_Vele3Line             (int)158
#define ParObject_ArteryRegister        (int)159
#define ParObject_Artery                (int)160
#define ParObject_Cnst_1d_art           (int)161
#define ParObject_ArrheniusSpec         (int)162
#define ParObject_ArrheniusTemp         (int)163
#define ParObject_So_Hex8P1J1           (int)164
#define ParObject_SoHex8P1J1Register    (int)165
#define ParObject_XDiff3                (int)166
#define ParObject_XDiff3Register        (int)167
#define ParObject_XDiff3Surface         (int)168
#define ParObject_XDiff3Line            (int)169
#define ParObject_ArrheniusPV           (int)170

namespace DRT
{

/*!
\brief A virtual class with functionality to pack, unpack and communicate
       classes in parallel

This class is used to pack information usually stored in a class in a vector<char>.
This vector<char> can then be used to communicate the contents of a class and to
read/write binary io. Every class derived from ParObject must basically implement
the Pack and the Unpack method. There are several methods (mpst of them template specializations)
to ease the work of packing/unpacking<br><br>
Here is an example:
\code
// stuff in a class 'Fool' that needs to be packed
int                      i;
double                   b;
double*                  vec = new double[50];
vector<char>             bla;
Epetra_SerialDenseMatrix matrix;
\endcode
This is how it is packed into a vector<char>& data:<br>
\code
Fool::Pack (vector< char > &data) const
{
  data.resize(0);                       // resize data to zero
  int tmp = UniqueParObjectId()         // get the unique parobject id
  AddtoPack(data,tmp);                  // pack this id
  AddtoPack(data,i);                    // pack i
  AddtoPack(data,b);                    // pack b
  AddtoPack(data,vec,50*sizeof(double); // pack vec
  AddtoPack(data,bla);                  // pack bla
  AddtoPack(data,matrix);               // pack matrix
  return;
}
\endcode
Here is how this data can be unpacked again:<br>
\code
Fool::Unpack(const vector< char > &data)
{
  int position = 0;                      // used to mark current reading position in data
  int tmp;
  ExtractfromPack(position,data,tmp);    // unpack the unique id
  if (tmp != UniqueParObjectId()) dserror("data does not belong to this class");
  ExtractfromPack(position,data,i);      // extract i
  ExtractfromPack(position,data,b);      // extract b
  ExtractfromPack(position,data,bla);    // extract bla
  ExtractfromPack(position,data,matrix); // extract matrix
  if (position != (int)data.size()) dserror("Mismatch in size of data");
  return;
}
\endcode
<br>
Some remarks:

- Data has to be unpacked the order it was packed

- The first object in every packed data set has to be the unique parobject id, see head of file drt_parobject.H

- The size of data ( data.size() ) must 'fit' exactly

- A class should pack everything it needs to be exatcly recreated on a different processor.
  this specifically holds for classes used in a Discretization where data might be shifted around processors.

- Every object that carefully implements ParObject can very easily be communicated using the Exporter.

- Every class that carefully implements ParObject can pretty easily be written/read to/from binary io

- A class derived from or more base classes is responsible to also pack and unpack the base classes' data
  by calling the base class' implementation of Pack/Unpack

- The intention of this class is to pack and communicate rather 'small' units of data. Though possible,
  it is not meant to be used at the system level to communicate huge data sets such as sparse matrices or
  vectors of system length. It does therefore not support any Epetra_Vector
  or Epetra_CrsMatrix objects and is not supposed to in the future either.

<br>
Here is a list of data types that are currently supported by the existing AddtoPack and ExtractfromPack methods:
\code
bool, bool*
char, char*
enum, enum*
int, int*
double, double*
float, float*
string
vector<char>  // especially useful to pack other packs into a pack, e.g. a class packing its own base class
vector<int>
vector<double>
Epetra_SerialDenseMatrix
\endcode

Note that trying to pack an unsupported type of data (such as e.g. map<key,T> ) might compile and link but will
result in the most (or least) funny behavior. Also, this type of bug might be extremely hard to find....
<br><br>
Of course, you are welcome to add more specializations to the existing AddtoPack and ExtractfromPack templates.
If you do so, please update this documentation.

\author gee (gee@lnm.mw.tum.de)
*/
class ParObject
{
public:



  /*!
  \brief Standard Constructor

  */
  ParObject();

  /*!
  \brief Copy Constructor

  Makes a deep copy of a ParObject

  */
  ParObject(const DRT::ParObject& old);

  /*!
  \brief Destructor

  */
  virtual ~ParObject();


  //! @name pure virtual Packing and Unpacking

  /*!
  \brief Return unique ParObject id

  every class implementing ParObject needs a unique id defined at the
  top of drt_parobject.H (this file) and should return it in this method.
  */
  virtual int UniqueParObjectId() const = 0;

  /*!
  \brief Pack this class so it can be communicated

  Resizes the vector data and stores all information of a class in it.
  The first information to be stored in data has to be the
  unique parobject id delivered by UniqueParObjectId() which will then
  identify the exact class on the receiving processor.

  \param data (in/out): char vector to store class information

  */
  virtual void Pack(vector<char>& data) const = 0;

  /*!
  \brief Unpack data from a char vector into this class

  The vector data contains all information to rebuild the
  exact copy of an instance of a class on a different processor.
  The first entry in data has to be an integer which is the unique
  parobject id defined at the top of this file and delivered by
  UniqueParObjectId().

  \param data (in) : vector storing all data to be unpacked into this
                     instance.
  */
  virtual void Unpack(const vector<char>& data) = 0;


  //@}

public:

  //! @name Routines to help pack stuff into a char vector

  /*!
  \brief Add stuff to the end of a char vector data

  This method is templated for all basic types int char double enum bool etc.
  \param data (in/out) : char string stuff shall be added to
  \param stuff (in)    : basic data type (float int double char etc) that get's added to stuff

  \note To be more precise, you can use this template for types of data that sizeof(kind)
        works for. Do not use for classes or structs or stl containers!

  */
  template<typename kind> static void AddtoPack(vector<char>& data, const kind& stuff)
  {
    int size = sizeof(kind);
    int oldsize = data.size();
    data.resize(oldsize+size);
    memcpy(&data[oldsize],&stuff,size);
    return;
  }

  /*!
  \brief Add stuff to the end of a char vector data

  This method is an overload of the above template for vector<int>
  \param data (in/out) : char string stuff shall be added to
  \param stuff (in)    : vector<int> that get's added to stuff

  */
  static void AddtoPack(vector<char>& data, const vector<int>& stuff);

  /*!
  \brief Add stuff to the end of a char vector data

  This method is an overload of the above template for vector<double>
  \param data (in/out) : char string stuff shall be added to
  \param stuff (in)    : vector<double> that get's added to stuff

  */
  static void AddtoPack(vector<char>& data, const vector<double>& stuff);

  /*!
  \brief Add stuff to the end of a char vector data

  This method is an overload of the above template for vector<char>
  \param data (in/out) : char string stuff shall be added to
  \param stuff (in)    : vector<char> that get's added to stuff

  */
  static void AddtoPack(vector<char>& data, const vector<char>& stuff);

  /*!
  \brief Add stuff to the end of a char vector data

  This method is an overload of the above template for Epetra_SerialDenseMatrix
  \param data (in/out) : char string stuff shall be added to
  \param stuff (in)    : Epetra_SerialDenseMatrix that get's added to stuff

  */
  static void AddtoPack(vector<char>& data, const Epetra_SerialDenseMatrix& stuff);

  /*!
  \brief Add stuff to the end of a char vector data

  This method is an overload of the above template for Epetra_SerialDenseMatrix
  \param data (in/out) : char string stuff shall be added to
  \param stuff (in)    : Epetra_SerialDenseVector that get's added to stuff

  */
  static void AddtoPack(vector<char>& data, const Epetra_SerialDenseVector& stuff);

  /*!
  \brief Add stuff to the end of a char vector data

  This method is an overload of the above template for Matrix
  \param data (in/out) : char string stuff shall be added to
  \param stuff (in)    : Matrix that get's added to data

  */
  template<unsigned i, unsigned j>
  static void AddtoPack(vector<char>& data, const LINALG::Matrix<i,j>& stuff)
  {
    AddtoPack(data, i);
    AddtoPack(data, j);
    AddtoPack(data, stuff.A(), stuff.M()*stuff.N()*sizeof(double));
  }

  /*!
  \brief Add stuff to the end of a char vector data

  This method is an overload of the above template for string
  \param data (in/out) : char string stuff shall be added to
  \param stuff (in)    : string that get's added to stuff

  */
  static void AddtoPack(vector<char>& data, const string& stuff);

  /*!
  \brief Add stuff to a char vector data

  This method adds an array to data
  \param data (in/out)  : char vector stuff shall be added to
  \param stuff (in)     : ptr to stuff that has length stuffsize (in byte)
  \param stuffsize (in) : length of stuff in byte
  */
  template<typename kind> static void AddtoPack(vector<char>& data, const kind* stuff, const int stuffsize)
  {
    int oldsize = data.size();
    data.resize(oldsize+stuffsize);
    memcpy(&data[oldsize],stuff,stuffsize);
    return;
  }

  /*!
  \brief Extract stuff from a char vector data and increment position

  This method is templated for all basic type like int char double enum.
  To be precise, it will work for all objects where sizeof(kind) is well defined.

  \param position (in/out) : place in data where to extract stuff. Position will be incremented by sizeof(kind)
  \param data (in)         : char vector where stuff is extracted from
  \param stuff (out)       : basic data type (float int double char ...) to extract from data

  */
  template<typename kind> static void ExtractfromPack(int& position, const vector<char>& data, kind& stuff)
  {
    int size = sizeof(kind);
    memcpy(&stuff,&data[position],size);
    position += size;
    return;
  }

  /*!
  \brief Extract stuff from a char vector data and increment position

  This method is an overload of the above template for stuff of type vector<int>

  \param position (in/out) : place in data where to extract stuff. Position will be incremented by this method
  \param data (in)         : char vector where stuff is extracted from
  \param stuff (out)       : vector<int> to extract from data

  */
  static void ExtractfromPack(int& position, const vector<char>& data, vector<int>& stuff);

  /*!
  \brief Extract stuff from a char vector data and increment position

  This method is an overload of the above template for stuff of type vector<double>

  \param position (in/out) : place in data where to extract stuff. Position will be incremented by this method
  \param data (in)         : char vector where stuff is extracted from
  \param stuff (out)       : vector<double> to extract from data

  */
  static void ExtractfromPack(int& position, const vector<char>& data, vector<double>& stuff);

  /*!
  \brief Extract stuff from a char vector data and increment position

  This method is an overload of the above template for stuff of type vector<char>

  \param position (in/out) : place in data where to extract stuff. Position will be incremented by this method
  \param data (in)         : char vector where stuff is extracted from
  \param stuff (out)       : vector<char> to extract from data

  */
  static void ExtractfromPack(int& position, const vector<char>& data, vector<char>& stuff);

  /*!
  \brief Extract stuff from a char vector data and increment position

  This method is an overload of the above template for stuff of type Epetra_SerialDenseMatrix

  \param position (in/out) : place in data where to extract stuff. Position will be incremented by this method
  \param data (in)         : char vector where stuff is extracted from
  \param stuff (out)       : Epetra_SerialDenseMatrix to extract from data

  */
  static void ExtractfromPack(int& position, const vector<char>& data, Epetra_SerialDenseMatrix& stuff);


  /*!
  \brief Extract stuff from a char vector data and increment position

  This method is an overload of the above template for stuff of type Epetra_SerialDenseVector

  \param position (in/out) : place in data where to extract stuff. Position will be incremented by this method
  \param data (in)         : char vector where stuff is extracted from
  \param stuff (out)       : Epetra_SerialDenseVector to extract from data

  */
  static void ExtractfromPack(int& position, const vector<char>& data, Epetra_SerialDenseVector& stuff);

  /*!
  \brief Extract stuff from a char vector data and increment position

  This method is an overload of the above template for stuff of type Matrix

  \param position (in/out) : place in data where to extract stuff. Position will be incremented by this method
  \param data (in)         : char vector where stuff is extracted from
  \param stuff (out)       : Matrix to extract from data

  */
  template <unsigned int i, unsigned int j>
  static void ExtractfromPack(int& position, const vector<char>& data, LINALG::Matrix<i,j>& stuff)
  {
    int m = 0;
    ExtractfromPack(position,data,m);
    if (m!=i)
      dserror("first dimension mismatch");
    int n = 0;
    ExtractfromPack(position,data,n);
    if (n!=j)
      dserror("second dimension mismatch");
    ExtractfromPack(position, data, stuff.A(), stuff.M()*stuff.N()*sizeof(double));
  }

  /*!
  \brief Extract stuff from a char vector data and increment position

  This method is an overload of the above template for stuff of type string

  \param position (in/out) : place in data where to extract stuff. Position will be incremented by this method
  \param data (in)         : char vector where stuff is extracted from
  \param stuff (out)       : string to extract from data

  */
  static void ExtractfromPack(int& position, const vector<char>& data, string& stuff);

  /*!
  \brief Extract stuff from a char vector data and increment position

  This method extracts an array from data
  \param position (in/out) : place in data where to extract stuff. Position will be incremented by stuffsize
  \param data (in)         : char string where stuff is extracted from
  \param stuff (out)       : array of total length stuffsize (in byte)
  \param stuffsize (in)    : length of stuff in byte
  */
  template<typename kind> static void ExtractfromPack(int& position, const vector<char>& data, kind* stuff, const int stuffsize)
  {
    memcpy(stuff,&data[position],stuffsize);
    position += stuffsize;
    return;
  }

  //@}



}; // class ParObject
} // namespace DRT











#endif  // #ifndef PAROBJECT_H
#endif  // #ifdef CCADISCRET

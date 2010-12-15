
#include <sstream>
#include <string>
#include "cut_test_generator.H"

#if 0

void OutputGenerator::Generate( GEO::CUT::Element* element, const tetgenio & out )
{
  std::stringstream str;
  str << "tet." << counter_;
  counter_ += 1;
  std::string name = str.str();

  const_cast<tetgenio &>( out ).save_nodes( const_cast<char*>( name.c_str() ) );
  //const_cast<tetgenio &>( out ).save_poly( const_cast<char*>( name.c_str() ) );
  const_cast<tetgenio &>( out ).save_elements( const_cast<char*>( name.c_str() ) );
  const_cast<tetgenio &>( out ).save_faces( const_cast<char*>( name.c_str() ) );
}

#endif

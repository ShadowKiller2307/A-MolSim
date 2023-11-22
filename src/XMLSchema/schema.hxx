// Copyright (c) 2005-2014 Code Synthesis Tools CC
//
// This program was generated by CodeSynthesis XSD, an XML Schema to
// C++ data binding compiler.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License version 2 as
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
//
// In addition, as a special exception, Code Synthesis Tools CC gives
// permission to link this program with the Xerces-C++ library (or with
// modified versions of Xerces-C++ that use the same license as Xerces-C++),
// and distribute linked combinations including the two. You must obey
// the GNU General Public License version 2 in all respects for all of
// the code used other than Xerces-C++. If you modify this copy of the
// program, you may extend this exception to your version of the program,
// but you are not obligated to do so. If you do not wish to do so, delete
// this exception statement from your version.
//
// Furthermore, Code Synthesis Tools CC makes a special exception for
// the Free/Libre and Open Source Software (FLOSS) which is described
// in the accompanying FLOSSE file.
//

#ifndef SCHEMA_HXX
#define SCHEMA_HXX

#ifndef XSD_USE_CHAR
#define XSD_USE_CHAR
#endif

#ifndef XSD_CXX_TREE_USE_CHAR
#define XSD_CXX_TREE_USE_CHAR
#endif

// Begin prologue.
//
//
// End prologue.

#include <xsd/cxx/config.hxx>

#if (XSD_INT_VERSION != 4000000L)
#error XSD runtime version mismatch
#endif

#include <xsd/cxx/pre.hxx>

#include <xsd/cxx/xml/char-utf8.hxx>

#include <xsd/cxx/tree/exceptions.hxx>
#include <xsd/cxx/tree/elements.hxx>
#include <xsd/cxx/tree/types.hxx>

#include <xsd/cxx/xml/error-handler.hxx>

#include <xsd/cxx/xml/dom/auto-ptr.hxx>

#include <xsd/cxx/tree/parsing.hxx>
#include <xsd/cxx/tree/parsing/byte.hxx>
#include <xsd/cxx/tree/parsing/unsigned-byte.hxx>
#include <xsd/cxx/tree/parsing/short.hxx>
#include <xsd/cxx/tree/parsing/unsigned-short.hxx>
#include <xsd/cxx/tree/parsing/int.hxx>
#include <xsd/cxx/tree/parsing/unsigned-int.hxx>
#include <xsd/cxx/tree/parsing/long.hxx>
#include <xsd/cxx/tree/parsing/unsigned-long.hxx>
#include <xsd/cxx/tree/parsing/boolean.hxx>
#include <xsd/cxx/tree/parsing/float.hxx>
#include <xsd/cxx/tree/parsing/double.hxx>
#include <xsd/cxx/tree/parsing/decimal.hxx>

namespace xml_schema
{
  // anyType and anySimpleType.
  //
  typedef ::xsd::cxx::tree::type type;
  typedef ::xsd::cxx::tree::simple_type< char, type > simple_type;
  typedef ::xsd::cxx::tree::type container;

  // 8-bit
  //
  typedef signed char byte;
  typedef unsigned char unsigned_byte;

  // 16-bit
  //
  typedef short short_;
  typedef unsigned short unsigned_short;

  // 32-bit
  //
  typedef int int_;
  typedef unsigned int unsigned_int;

  // 64-bit
  //
  typedef long long long_;
  typedef unsigned long long unsigned_long;

  // Supposed to be arbitrary-length integral types.
  //
  typedef long long integer;
  typedef long long non_positive_integer;
  typedef unsigned long long non_negative_integer;
  typedef unsigned long long positive_integer;
  typedef long long negative_integer;

  // Boolean.
  //
  typedef bool boolean;

  // Floating-point types.
  //
  typedef float float_;
  typedef double double_;
  typedef double decimal;

  // String types.
  //
  typedef ::xsd::cxx::tree::string< char, simple_type > string;
  typedef ::xsd::cxx::tree::normalized_string< char, string > normalized_string;
  typedef ::xsd::cxx::tree::token< char, normalized_string > token;
  typedef ::xsd::cxx::tree::name< char, token > name;
  typedef ::xsd::cxx::tree::nmtoken< char, token > nmtoken;
  typedef ::xsd::cxx::tree::nmtokens< char, simple_type, nmtoken > nmtokens;
  typedef ::xsd::cxx::tree::ncname< char, name > ncname;
  typedef ::xsd::cxx::tree::language< char, token > language;

  // ID/IDREF.
  //
  typedef ::xsd::cxx::tree::id< char, ncname > id;
  typedef ::xsd::cxx::tree::idref< char, ncname, type > idref;
  typedef ::xsd::cxx::tree::idrefs< char, simple_type, idref > idrefs;

  // URI.
  //
  typedef ::xsd::cxx::tree::uri< char, simple_type > uri;

  // Qualified name.
  //
  typedef ::xsd::cxx::tree::qname< char, simple_type, uri, ncname > qname;

  // Binary.
  //
  typedef ::xsd::cxx::tree::buffer< char > buffer;
  typedef ::xsd::cxx::tree::base64_binary< char, simple_type > base64_binary;
  typedef ::xsd::cxx::tree::hex_binary< char, simple_type > hex_binary;

  // Date/time.
  //
  typedef ::xsd::cxx::tree::time_zone time_zone;
  typedef ::xsd::cxx::tree::date< char, simple_type > date;
  typedef ::xsd::cxx::tree::date_time< char, simple_type > date_time;
  typedef ::xsd::cxx::tree::duration< char, simple_type > duration;
  typedef ::xsd::cxx::tree::gday< char, simple_type > gday;
  typedef ::xsd::cxx::tree::gmonth< char, simple_type > gmonth;
  typedef ::xsd::cxx::tree::gmonth_day< char, simple_type > gmonth_day;
  typedef ::xsd::cxx::tree::gyear< char, simple_type > gyear;
  typedef ::xsd::cxx::tree::gyear_month< char, simple_type > gyear_month;
  typedef ::xsd::cxx::tree::time< char, simple_type > time;

  // Entity.
  //
  typedef ::xsd::cxx::tree::entity< char, ncname > entity;
  typedef ::xsd::cxx::tree::entities< char, simple_type, entity > entities;

  typedef ::xsd::cxx::tree::content_order content_order;
  // Flags and properties.
  //
  typedef ::xsd::cxx::tree::flags flags;
  typedef ::xsd::cxx::tree::properties< char > properties;

  // Parsing/serialization diagnostics.
  //
  typedef ::xsd::cxx::tree::severity severity;
  typedef ::xsd::cxx::tree::error< char > error;
  typedef ::xsd::cxx::tree::diagnostics< char > diagnostics;

  // Exceptions.
  //
  typedef ::xsd::cxx::tree::exception< char > exception;
  typedef ::xsd::cxx::tree::bounds< char > bounds;
  typedef ::xsd::cxx::tree::duplicate_id< char > duplicate_id;
  typedef ::xsd::cxx::tree::parsing< char > parsing;
  typedef ::xsd::cxx::tree::expected_element< char > expected_element;
  typedef ::xsd::cxx::tree::unexpected_element< char > unexpected_element;
  typedef ::xsd::cxx::tree::expected_attribute< char > expected_attribute;
  typedef ::xsd::cxx::tree::unexpected_enumerator< char > unexpected_enumerator;
  typedef ::xsd::cxx::tree::expected_text_content< char > expected_text_content;
  typedef ::xsd::cxx::tree::no_prefix_mapping< char > no_prefix_mapping;

  // Error handler callback interface.
  //
  typedef ::xsd::cxx::xml::error_handler< char > error_handler;

  // DOM interaction.
  //
  namespace dom
  {
    // Automatic pointer for DOMDocument.
    //
    using ::xsd::cxx::xml::dom::auto_ptr;

#ifndef XSD_CXX_TREE_TREE_NODE_KEY__XML_SCHEMA
#define XSD_CXX_TREE_TREE_NODE_KEY__XML_SCHEMA
    // DOM user data key for back pointers to tree nodes.
    //
    const XMLCh* const tree_node_key = ::xsd::cxx::tree::user_data_keys::node;
#endif
  }
}

// Forward declarations.
//
class doubleList;
class unsignedIntList;
class arrayOfThreeDoubles;
class arrayOfThreeUnsignedInts;
class Cuboid;
class simulationConfig;

#include <memory>    // ::std::auto_ptr
#include <limits>    // std::numeric_limits
#include <algorithm> // std::binary_search

#include <xsd/cxx/xml/char-utf8.hxx>

#include <xsd/cxx/tree/exceptions.hxx>
#include <xsd/cxx/tree/elements.hxx>
#include <xsd/cxx/tree/containers.hxx>
#include <xsd/cxx/tree/list.hxx>

#include <xsd/cxx/xml/dom/parsing-header.hxx>

class doubleList: public ::xml_schema::simple_type,
  public ::xsd::cxx::tree::list< ::xml_schema::double_, char, ::xsd::cxx::tree::schema_type::double_ >
{
  public:
  doubleList ();

  doubleList (size_type n, const ::xml_schema::double_& x);

  template < typename I >
  doubleList (const I& begin, const I& end)
  : ::xsd::cxx::tree::list< ::xml_schema::double_, char, ::xsd::cxx::tree::schema_type::double_ > (begin, end, this)
  {
  }

  doubleList (const ::xercesc::DOMElement& e,
              ::xml_schema::flags f = 0,
              ::xml_schema::container* c = 0);

  doubleList (const ::xercesc::DOMAttr& a,
              ::xml_schema::flags f = 0,
              ::xml_schema::container* c = 0);

  doubleList (const ::std::string& s,
              const ::xercesc::DOMElement* e,
              ::xml_schema::flags f = 0,
              ::xml_schema::container* c = 0);

  doubleList (const doubleList& x,
              ::xml_schema::flags f = 0,
              ::xml_schema::container* c = 0);

  virtual doubleList*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  virtual 
  ~doubleList ();
};

class unsignedIntList: public ::xml_schema::simple_type,
  public ::xsd::cxx::tree::list< ::xml_schema::unsigned_int, char >
{
  public:
  unsignedIntList ();

  unsignedIntList (size_type n, const ::xml_schema::unsigned_int& x);

  template < typename I >
  unsignedIntList (const I& begin, const I& end)
  : ::xsd::cxx::tree::list< ::xml_schema::unsigned_int, char > (begin, end, this)
  {
  }

  unsignedIntList (const ::xercesc::DOMElement& e,
                   ::xml_schema::flags f = 0,
                   ::xml_schema::container* c = 0);

  unsignedIntList (const ::xercesc::DOMAttr& a,
                   ::xml_schema::flags f = 0,
                   ::xml_schema::container* c = 0);

  unsignedIntList (const ::std::string& s,
                   const ::xercesc::DOMElement* e,
                   ::xml_schema::flags f = 0,
                   ::xml_schema::container* c = 0);

  unsignedIntList (const unsignedIntList& x,
                   ::xml_schema::flags f = 0,
                   ::xml_schema::container* c = 0);

  virtual unsignedIntList*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  virtual 
  ~unsignedIntList ();
};

class arrayOfThreeDoubles: public ::doubleList
{
  public:
  // Constructors.
  //
  arrayOfThreeDoubles ();

  arrayOfThreeDoubles (const ::doubleList&);

  arrayOfThreeDoubles (const ::xercesc::DOMElement& e,
                       ::xml_schema::flags f = 0,
                       ::xml_schema::container* c = 0);

  arrayOfThreeDoubles (const ::xercesc::DOMAttr& a,
                       ::xml_schema::flags f = 0,
                       ::xml_schema::container* c = 0);

  arrayOfThreeDoubles (const ::std::string& s,
                       const ::xercesc::DOMElement* e,
                       ::xml_schema::flags f = 0,
                       ::xml_schema::container* c = 0);

  arrayOfThreeDoubles (const arrayOfThreeDoubles& x,
                       ::xml_schema::flags f = 0,
                       ::xml_schema::container* c = 0);

  virtual arrayOfThreeDoubles*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  virtual 
  ~arrayOfThreeDoubles ();
};

class arrayOfThreeUnsignedInts: public ::unsignedIntList
{
  public:
  // Constructors.
  //
  arrayOfThreeUnsignedInts ();

  arrayOfThreeUnsignedInts (const ::unsignedIntList&);

  arrayOfThreeUnsignedInts (const ::xercesc::DOMElement& e,
                            ::xml_schema::flags f = 0,
                            ::xml_schema::container* c = 0);

  arrayOfThreeUnsignedInts (const ::xercesc::DOMAttr& a,
                            ::xml_schema::flags f = 0,
                            ::xml_schema::container* c = 0);

  arrayOfThreeUnsignedInts (const ::std::string& s,
                            const ::xercesc::DOMElement* e,
                            ::xml_schema::flags f = 0,
                            ::xml_schema::container* c = 0);

  arrayOfThreeUnsignedInts (const arrayOfThreeUnsignedInts& x,
                            ::xml_schema::flags f = 0,
                            ::xml_schema::container* c = 0);

  virtual arrayOfThreeUnsignedInts*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  virtual 
  ~arrayOfThreeUnsignedInts ();
};

class Cuboid: public ::xml_schema::type
{
  public:
  // llfc
  //
  typedef ::arrayOfThreeDoubles llfc_type;
  typedef ::xsd::cxx::tree::traits< llfc_type, char > llfc_traits;

  const llfc_type&
  llfc () const;

  llfc_type&
  llfc ();

  void
  llfc (const llfc_type& x);

  void
  llfc (::std::auto_ptr< llfc_type > p);

  // particlePerDimension
  //
  typedef ::arrayOfThreeUnsignedInts particlePerDimension_type;
  typedef ::xsd::cxx::tree::traits< particlePerDimension_type, char > particlePerDimension_traits;

  const particlePerDimension_type&
  particlePerDimension () const;

  particlePerDimension_type&
  particlePerDimension ();

  void
  particlePerDimension (const particlePerDimension_type& x);

  void
  particlePerDimension (::std::auto_ptr< particlePerDimension_type > p);

  // particleVelocity
  //
  typedef ::arrayOfThreeDoubles particleVelocity_type;
  typedef ::xsd::cxx::tree::traits< particleVelocity_type, char > particleVelocity_traits;

  const particleVelocity_type&
  particleVelocity () const;

  particleVelocity_type&
  particleVelocity ();

  void
  particleVelocity (const particleVelocity_type& x);

  void
  particleVelocity (::std::auto_ptr< particleVelocity_type > p);

  // h
  //
  typedef ::xml_schema::double_ h_type;
  typedef ::xsd::cxx::tree::traits< h_type, char, ::xsd::cxx::tree::schema_type::double_ > h_traits;

  const h_type&
  h () const;

  h_type&
  h ();

  void
  h (const h_type& x);

  // mass
  //
  typedef ::xml_schema::double_ mass_type;
  typedef ::xsd::cxx::tree::traits< mass_type, char, ::xsd::cxx::tree::schema_type::double_ > mass_traits;

  const mass_type&
  mass () const;

  mass_type&
  mass ();

  void
  mass (const mass_type& x);

  // generateNumber
  //
  typedef ::xml_schema::int_ generateNumber_type;
  typedef ::xsd::cxx::tree::traits< generateNumber_type, char > generateNumber_traits;

  const generateNumber_type&
  generateNumber () const;

  generateNumber_type&
  generateNumber ();

  void
  generateNumber (const generateNumber_type& x);

  // Constructors.
  //
  Cuboid (const llfc_type&,
          const particlePerDimension_type&,
          const particleVelocity_type&,
          const h_type&,
          const mass_type&,
          const generateNumber_type&);

  Cuboid (const ::xercesc::DOMElement& e,
          ::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0);

  Cuboid (const Cuboid& x,
          ::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0);

  virtual Cuboid*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  Cuboid&
  operator= (const Cuboid& x);

  virtual 
  ~Cuboid ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< llfc_type > llfc_;
  ::xsd::cxx::tree::one< particlePerDimension_type > particlePerDimension_;
  ::xsd::cxx::tree::one< particleVelocity_type > particleVelocity_;
  ::xsd::cxx::tree::one< h_type > h_;
  ::xsd::cxx::tree::one< mass_type > mass_;
  ::xsd::cxx::tree::one< generateNumber_type > generateNumber_;
};

class simulationConfig: public ::xml_schema::type
{
  public:
  // baseName
  //
  typedef ::xml_schema::string baseName_type;
  typedef ::xsd::cxx::tree::traits< baseName_type, char > baseName_traits;

  const baseName_type&
  baseName () const;

  baseName_type&
  baseName ();

  void
  baseName (const baseName_type& x);

  void
  baseName (::std::auto_ptr< baseName_type > p);

  // writeFrequency
  //
  typedef ::xml_schema::int_ writeFrequency_type;
  typedef ::xsd::cxx::tree::traits< writeFrequency_type, char > writeFrequency_traits;

  const writeFrequency_type&
  writeFrequency () const;

  writeFrequency_type&
  writeFrequency ();

  void
  writeFrequency (const writeFrequency_type& x);

  // t_end
  //
  typedef ::xml_schema::double_ t_end_type;
  typedef ::xsd::cxx::tree::traits< t_end_type, char, ::xsd::cxx::tree::schema_type::double_ > t_end_traits;

  const t_end_type&
  t_end () const;

  t_end_type&
  t_end ();

  void
  t_end (const t_end_type& x);

  // delta_t
  //
  typedef ::xml_schema::double_ delta_t_type;
  typedef ::xsd::cxx::tree::traits< delta_t_type, char, ::xsd::cxx::tree::schema_type::double_ > delta_t_traits;

  const delta_t_type&
  delta_t () const;

  delta_t_type&
  delta_t ();

  void
  delta_t (const delta_t_type& x);

  // logLevel
  //
  typedef ::xml_schema::int_ logLevel_type;
  typedef ::xsd::cxx::tree::traits< logLevel_type, char > logLevel_traits;

  const logLevel_type&
  logLevel () const;

  logLevel_type&
  logLevel ();

  void
  logLevel (const logLevel_type& x);

  // inputGenerator
  //
  typedef ::xml_schema::boolean inputGenerator_type;
  typedef ::xsd::cxx::tree::traits< inputGenerator_type, char > inputGenerator_traits;

  const inputGenerator_type&
  inputGenerator () const;

  inputGenerator_type&
  inputGenerator ();

  void
  inputGenerator (const inputGenerator_type& x);

  // inputText
  //
  typedef ::xml_schema::boolean inputText_type;
  typedef ::xsd::cxx::tree::traits< inputText_type, char > inputText_traits;

  const inputText_type&
  inputText () const;

  inputText_type&
  inputText ();

  void
  inputText (const inputText_type& x);

  // inputPicture
  //
  typedef ::xml_schema::boolean inputPicture_type;
  typedef ::xsd::cxx::tree::traits< inputPicture_type, char > inputPicture_traits;

  const inputPicture_type&
  inputPicture () const;

  inputPicture_type&
  inputPicture ();

  void
  inputPicture (const inputPicture_type& x);

  // inputXML
  //
  typedef ::xml_schema::boolean inputXML_type;
  typedef ::xsd::cxx::tree::traits< inputXML_type, char > inputXML_traits;

  const inputXML_type&
  inputXML () const;

  inputXML_type&
  inputXML ();

  void
  inputXML (const inputXML_type& x);

  // Cuboid
  //
  typedef ::Cuboid Cuboid_type;
  typedef ::xsd::cxx::tree::sequence< Cuboid_type > Cuboid_sequence;
  typedef Cuboid_sequence::iterator Cuboid_iterator;
  typedef Cuboid_sequence::const_iterator Cuboid_const_iterator;
  typedef ::xsd::cxx::tree::traits< Cuboid_type, char > Cuboid_traits;

  const Cuboid_sequence&
  Cuboid () const;

  Cuboid_sequence&
  Cuboid ();

  void
  Cuboid (const Cuboid_sequence& s);

  // Constructors.
  //
  simulationConfig (const baseName_type&,
                    const writeFrequency_type&,
                    const t_end_type&,
                    const delta_t_type&,
                    const logLevel_type&,
                    const inputGenerator_type&,
                    const inputText_type&,
                    const inputPicture_type&,
                    const inputXML_type&);

  simulationConfig (const ::xercesc::DOMElement& e,
                    ::xml_schema::flags f = 0,
                    ::xml_schema::container* c = 0);

  simulationConfig (const simulationConfig& x,
                    ::xml_schema::flags f = 0,
                    ::xml_schema::container* c = 0);

  virtual simulationConfig*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  simulationConfig&
  operator= (const simulationConfig& x);

  virtual 
  ~simulationConfig ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< baseName_type > baseName_;
  ::xsd::cxx::tree::one< writeFrequency_type > writeFrequency_;
  ::xsd::cxx::tree::one< t_end_type > t_end_;
  ::xsd::cxx::tree::one< delta_t_type > delta_t_;
  ::xsd::cxx::tree::one< logLevel_type > logLevel_;
  ::xsd::cxx::tree::one< inputGenerator_type > inputGenerator_;
  ::xsd::cxx::tree::one< inputText_type > inputText_;
  ::xsd::cxx::tree::one< inputPicture_type > inputPicture_;
  ::xsd::cxx::tree::one< inputXML_type > inputXML_;
  Cuboid_sequence Cuboid_;
};

#include <iosfwd>

#include <xercesc/sax/InputSource.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMErrorHandler.hpp>

// Parse a URI or a local file.
//

::std::auto_ptr< ::simulationConfig >
Configuration (const ::std::string& uri,
               ::xml_schema::flags f = 0,
               const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::simulationConfig >
Configuration (const ::std::string& uri,
               ::xml_schema::error_handler& eh,
               ::xml_schema::flags f = 0,
               const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::simulationConfig >
Configuration (const ::std::string& uri,
               ::xercesc::DOMErrorHandler& eh,
               ::xml_schema::flags f = 0,
               const ::xml_schema::properties& p = ::xml_schema::properties ());

// Parse std::istream.
//

::std::auto_ptr< ::simulationConfig >
Configuration (::std::istream& is,
               ::xml_schema::flags f = 0,
               const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::simulationConfig >
Configuration (::std::istream& is,
               ::xml_schema::error_handler& eh,
               ::xml_schema::flags f = 0,
               const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::simulationConfig >
Configuration (::std::istream& is,
               ::xercesc::DOMErrorHandler& eh,
               ::xml_schema::flags f = 0,
               const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::simulationConfig >
Configuration (::std::istream& is,
               const ::std::string& id,
               ::xml_schema::flags f = 0,
               const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::simulationConfig >
Configuration (::std::istream& is,
               const ::std::string& id,
               ::xml_schema::error_handler& eh,
               ::xml_schema::flags f = 0,
               const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::simulationConfig >
Configuration (::std::istream& is,
               const ::std::string& id,
               ::xercesc::DOMErrorHandler& eh,
               ::xml_schema::flags f = 0,
               const ::xml_schema::properties& p = ::xml_schema::properties ());

// Parse xercesc::InputSource.
//

::std::auto_ptr< ::simulationConfig >
Configuration (::xercesc::InputSource& is,
               ::xml_schema::flags f = 0,
               const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::simulationConfig >
Configuration (::xercesc::InputSource& is,
               ::xml_schema::error_handler& eh,
               ::xml_schema::flags f = 0,
               const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::simulationConfig >
Configuration (::xercesc::InputSource& is,
               ::xercesc::DOMErrorHandler& eh,
               ::xml_schema::flags f = 0,
               const ::xml_schema::properties& p = ::xml_schema::properties ());

// Parse xercesc::DOMDocument.
//

::std::auto_ptr< ::simulationConfig >
Configuration (const ::xercesc::DOMDocument& d,
               ::xml_schema::flags f = 0,
               const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::simulationConfig >
Configuration (::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d,
               ::xml_schema::flags f = 0,
               const ::xml_schema::properties& p = ::xml_schema::properties ());

#include <xsd/cxx/post.hxx>

// Begin epilogue.
//
//
// End epilogue.

#endif // SCHEMA_HXX

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

// Begin prologue.
//
//
// End prologue.

#include <xsd/cxx/pre.hxx>

#include "schema.h"

// doubleList
//

doubleList::
doubleList ()
: ::xsd::cxx::tree::list< ::xml_schema::double_, char, ::xsd::cxx::tree::schema_type::double_ > (this)
{
}

doubleList::
doubleList (size_type n, const ::xml_schema::double_& x)
: ::xsd::cxx::tree::list< ::xml_schema::double_, char, ::xsd::cxx::tree::schema_type::double_ > (n, x, this)
{
}

doubleList::
doubleList (const doubleList& o,
            ::xml_schema::flags f,
            ::xml_schema::container* c)
: ::xml_schema::simple_type (o, f, c),
  ::xsd::cxx::tree::list< ::xml_schema::double_, char, ::xsd::cxx::tree::schema_type::double_ > (o, f, this)
{
}

// unsignedIntList
//

unsignedIntList::
unsignedIntList ()
: ::xsd::cxx::tree::list< ::xml_schema::unsigned_int, char > (this)
{
}

unsignedIntList::
unsignedIntList (size_type n, const ::xml_schema::unsigned_int& x)
: ::xsd::cxx::tree::list< ::xml_schema::unsigned_int, char > (n, x, this)
{
}

unsignedIntList::
unsignedIntList (const unsignedIntList& o,
                 ::xml_schema::flags f,
                 ::xml_schema::container* c)
: ::xml_schema::simple_type (o, f, c),
  ::xsd::cxx::tree::list< ::xml_schema::unsigned_int, char > (o, f, this)
{
}

// arrayOfThreeDoubles
// 


// arrayOfThreeUnsignedInts
// 


// Cuboid
// 

const Cuboid::llfc_type& Cuboid::
llfc () const
{
  return this->llfc_.get ();
}

Cuboid::llfc_type& Cuboid::
llfc ()
{
  return this->llfc_.get ();
}

void Cuboid::
llfc (const llfc_type& x)
{
  this->llfc_.set (x);
}

void Cuboid::
llfc (::std::unique_ptr< llfc_type > x)
{
  this->llfc_.set (std::move (x));
}

const Cuboid::particlePerDimension_type& Cuboid::
particlePerDimension () const
{
  return this->particlePerDimension_.get ();
}

Cuboid::particlePerDimension_type& Cuboid::
particlePerDimension ()
{
  return this->particlePerDimension_.get ();
}

void Cuboid::
particlePerDimension (const particlePerDimension_type& x)
{
  this->particlePerDimension_.set (x);
}

void Cuboid::
particlePerDimension (::std::unique_ptr< particlePerDimension_type > x)
{
  this->particlePerDimension_.set (std::move (x));
}

const Cuboid::particleVelocity_type& Cuboid::
particleVelocity () const
{
  return this->particleVelocity_.get ();
}

Cuboid::particleVelocity_type& Cuboid::
particleVelocity ()
{
  return this->particleVelocity_.get ();
}

void Cuboid::
particleVelocity (const particleVelocity_type& x)
{
  this->particleVelocity_.set (x);
}

void Cuboid::
particleVelocity (::std::unique_ptr< particleVelocity_type > x)
{
  this->particleVelocity_.set (std::move (x));
}

const Cuboid::h_type& Cuboid::
h () const
{
  return this->h_.get ();
}

Cuboid::h_type& Cuboid::
h ()
{
  return this->h_.get ();
}

void Cuboid::
h (const h_type& x)
{
  this->h_.set (x);
}

const Cuboid::mass_type& Cuboid::
mass () const
{
  return this->mass_.get ();
}

Cuboid::mass_type& Cuboid::
mass ()
{
  return this->mass_.get ();
}

void Cuboid::
mass (const mass_type& x)
{
  this->mass_.set (x);
}

const Cuboid::generateNumber_type& Cuboid::
generateNumber () const
{
  return this->generateNumber_.get ();
}

Cuboid::generateNumber_type& Cuboid::
generateNumber ()
{
  return this->generateNumber_.get ();
}

void Cuboid::
generateNumber (const generateNumber_type& x)
{
  this->generateNumber_.set (x);
}


// Sphere
// 

const Sphere::centerCoordinates_type& Sphere::
centerCoordinates () const
{
  return this->centerCoordinates_.get ();
}

Sphere::centerCoordinates_type& Sphere::
centerCoordinates ()
{
  return this->centerCoordinates_.get ();
}

void Sphere::
centerCoordinates (const centerCoordinates_type& x)
{
  this->centerCoordinates_.set (x);
}

void Sphere::
centerCoordinates (::std::unique_ptr< centerCoordinates_type > x)
{
  this->centerCoordinates_.set (std::move (x));
}

const Sphere::initialVelocity_type& Sphere::
initialVelocity () const
{
  return this->initialVelocity_.get ();
}

Sphere::initialVelocity_type& Sphere::
initialVelocity ()
{
  return this->initialVelocity_.get ();
}

void Sphere::
initialVelocity (const initialVelocity_type& x)
{
  this->initialVelocity_.set (x);
}

void Sphere::
initialVelocity (::std::unique_ptr< initialVelocity_type > x)
{
  this->initialVelocity_.set (std::move (x));
}

const Sphere::radius_type& Sphere::
radius () const
{
  return this->radius_.get ();
}

Sphere::radius_type& Sphere::
radius ()
{
  return this->radius_.get ();
}

void Sphere::
radius (const radius_type& x)
{
  this->radius_.set (x);
}

const Sphere::distance_type& Sphere::
distance () const
{
  return this->distance_.get ();
}

Sphere::distance_type& Sphere::
distance ()
{
  return this->distance_.get ();
}

void Sphere::
distance (const distance_type& x)
{
  this->distance_.set (x);
}

const Sphere::mass_type& Sphere::
mass () const
{
  return this->mass_.get ();
}

Sphere::mass_type& Sphere::
mass ()
{
  return this->mass_.get ();
}

void Sphere::
mass (const mass_type& x)
{
  this->mass_.set (x);
}


// simulationConfig
// 

const simulationConfig::baseName_type& simulationConfig::
baseName () const
{
  return this->baseName_.get ();
}

simulationConfig::baseName_type& simulationConfig::
baseName ()
{
  return this->baseName_.get ();
}

void simulationConfig::
baseName (const baseName_type& x)
{
  this->baseName_.set (x);
}

void simulationConfig::
baseName (::std::unique_ptr< baseName_type > x)
{
  this->baseName_.set (std::move (x));
}

const simulationConfig::writeFrequency_type& simulationConfig::
writeFrequency () const
{
  return this->writeFrequency_.get ();
}

simulationConfig::writeFrequency_type& simulationConfig::
writeFrequency ()
{
  return this->writeFrequency_.get ();
}

void simulationConfig::
writeFrequency (const writeFrequency_type& x)
{
  this->writeFrequency_.set (x);
}

const simulationConfig::t_end_type& simulationConfig::
t_end () const
{
  return this->t_end_.get ();
}

simulationConfig::t_end_type& simulationConfig::
t_end ()
{
  return this->t_end_.get ();
}

void simulationConfig::
t_end (const t_end_type& x)
{
  this->t_end_.set (x);
}

const simulationConfig::delta_t_type& simulationConfig::
delta_t () const
{
  return this->delta_t_.get ();
}

simulationConfig::delta_t_type& simulationConfig::
delta_t ()
{
  return this->delta_t_.get ();
}

void simulationConfig::
delta_t (const delta_t_type& x)
{
  this->delta_t_.set (x);
}

const simulationConfig::logLevel_type& simulationConfig::
logLevel () const
{
  return this->logLevel_.get ();
}

simulationConfig::logLevel_type& simulationConfig::
logLevel ()
{
  return this->logLevel_.get ();
}

void simulationConfig::
logLevel (const logLevel_type& x)
{
  this->logLevel_.set (x);
}

const simulationConfig::domainSize_type& simulationConfig::
domainSize () const
{
  return this->domainSize_.get ();
}

simulationConfig::domainSize_type& simulationConfig::
domainSize ()
{
  return this->domainSize_.get ();
}

void simulationConfig::
domainSize (const domainSize_type& x)
{
  this->domainSize_.set (x);
}

void simulationConfig::
domainSize (::std::unique_ptr< domainSize_type > x)
{
  this->domainSize_.set (std::move (x));
}

const simulationConfig::containerType_type& simulationConfig::
containerType () const
{
  return this->containerType_.get ();
}

simulationConfig::containerType_type& simulationConfig::
containerType ()
{
  return this->containerType_.get ();
}

void simulationConfig::
containerType (const containerType_type& x)
{
  this->containerType_.set (x);
}

void simulationConfig::
containerType (::std::unique_ptr< containerType_type > x)
{
  this->containerType_.set (std::move (x));
}

const simulationConfig::boundaries_type& simulationConfig::
boundaries () const
{
  return this->boundaries_.get ();
}

simulationConfig::boundaries_type& simulationConfig::
boundaries ()
{
  return this->boundaries_.get ();
}

void simulationConfig::
boundaries (const boundaries_type& x)
{
  this->boundaries_.set (x);
}

void simulationConfig::
boundaries (::std::unique_ptr< boundaries_type > x)
{
  this->boundaries_.set (std::move (x));
}

const simulationConfig::cutOffRadius_type& simulationConfig::
cutOffRadius () const
{
  return this->cutOffRadius_.get ();
}

simulationConfig::cutOffRadius_type& simulationConfig::
cutOffRadius ()
{
  return this->cutOffRadius_.get ();
}

void simulationConfig::
cutOffRadius (const cutOffRadius_type& x)
{
  this->cutOffRadius_.set (x);
}

const simulationConfig::Cuboid_sequence& simulationConfig::
Cuboid () const
{
  return this->Cuboid_;
}

simulationConfig::Cuboid_sequence& simulationConfig::
Cuboid ()
{
  return this->Cuboid_;
}

void simulationConfig::
Cuboid (const Cuboid_sequence& s)
{
  this->Cuboid_ = s;
}

const simulationConfig::Sphere_sequence& simulationConfig::
Sphere () const
{
  return this->Sphere_;
}

simulationConfig::Sphere_sequence& simulationConfig::
Sphere ()
{
  return this->Sphere_;
}

void simulationConfig::
Sphere (const Sphere_sequence& s)
{
  this->Sphere_ = s;
}


#include <xsd/cxx/xml/dom/parsing-source.hxx>

// doubleList
//

doubleList::
doubleList (const ::xercesc::DOMElement& e,
            ::xml_schema::flags f,
            ::xml_schema::container* c)
: ::xml_schema::simple_type (e, f, c),
  ::xsd::cxx::tree::list< ::xml_schema::double_, char, ::xsd::cxx::tree::schema_type::double_ > (e, f, this)
{
}

doubleList::
doubleList (const ::xercesc::DOMAttr& a,
            ::xml_schema::flags f,
            ::xml_schema::container* c)
: ::xml_schema::simple_type (a, f, c),
  ::xsd::cxx::tree::list< ::xml_schema::double_, char, ::xsd::cxx::tree::schema_type::double_ > (a, f, this)
{
}

doubleList::
doubleList (const ::std::string& s,
            const ::xercesc::DOMElement* e,
            ::xml_schema::flags f,
            ::xml_schema::container* c)
: ::xml_schema::simple_type (s, e, f, c),
  ::xsd::cxx::tree::list< ::xml_schema::double_, char, ::xsd::cxx::tree::schema_type::double_ > (s, e, f, this)
{
}

doubleList* doubleList::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class doubleList (*this, f, c);
}

doubleList::
~doubleList ()
{
}

// unsignedIntList
//

unsignedIntList::
unsignedIntList (const ::xercesc::DOMElement& e,
                 ::xml_schema::flags f,
                 ::xml_schema::container* c)
: ::xml_schema::simple_type (e, f, c),
  ::xsd::cxx::tree::list< ::xml_schema::unsigned_int, char > (e, f, this)
{
}

unsignedIntList::
unsignedIntList (const ::xercesc::DOMAttr& a,
                 ::xml_schema::flags f,
                 ::xml_schema::container* c)
: ::xml_schema::simple_type (a, f, c),
  ::xsd::cxx::tree::list< ::xml_schema::unsigned_int, char > (a, f, this)
{
}

unsignedIntList::
unsignedIntList (const ::std::string& s,
                 const ::xercesc::DOMElement* e,
                 ::xml_schema::flags f,
                 ::xml_schema::container* c)
: ::xml_schema::simple_type (s, e, f, c),
  ::xsd::cxx::tree::list< ::xml_schema::unsigned_int, char > (s, e, f, this)
{
}

unsignedIntList* unsignedIntList::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class unsignedIntList (*this, f, c);
}

unsignedIntList::
~unsignedIntList ()
{
}

// arrayOfThreeDoubles
//

arrayOfThreeDoubles::
arrayOfThreeDoubles ()
: ::doubleList ()
{
}

arrayOfThreeDoubles::
arrayOfThreeDoubles (const ::doubleList& _xsd_doubleList_base)
: ::doubleList (_xsd_doubleList_base)
{
}

arrayOfThreeDoubles::
arrayOfThreeDoubles (const arrayOfThreeDoubles& x,
                     ::xml_schema::flags f,
                     ::xml_schema::container* c)
: ::doubleList (x, f, c)
{
}

arrayOfThreeDoubles::
arrayOfThreeDoubles (const ::xercesc::DOMElement& e,
                     ::xml_schema::flags f,
                     ::xml_schema::container* c)
: ::doubleList (e, f, c)
{
}

arrayOfThreeDoubles::
arrayOfThreeDoubles (const ::xercesc::DOMAttr& a,
                     ::xml_schema::flags f,
                     ::xml_schema::container* c)
: ::doubleList (a, f, c)
{
}

arrayOfThreeDoubles::
arrayOfThreeDoubles (const ::std::string& s,
                     const ::xercesc::DOMElement* e,
                     ::xml_schema::flags f,
                     ::xml_schema::container* c)
: ::doubleList (s, e, f, c)
{
}

arrayOfThreeDoubles* arrayOfThreeDoubles::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class arrayOfThreeDoubles (*this, f, c);
}

arrayOfThreeDoubles::
~arrayOfThreeDoubles ()
{
}

// arrayOfThreeUnsignedInts
//

arrayOfThreeUnsignedInts::
arrayOfThreeUnsignedInts ()
: ::unsignedIntList ()
{
}

arrayOfThreeUnsignedInts::
arrayOfThreeUnsignedInts (const ::unsignedIntList& _xsd_unsignedIntList_base)
: ::unsignedIntList (_xsd_unsignedIntList_base)
{
}

arrayOfThreeUnsignedInts::
arrayOfThreeUnsignedInts (const arrayOfThreeUnsignedInts& x,
                          ::xml_schema::flags f,
                          ::xml_schema::container* c)
: ::unsignedIntList (x, f, c)
{
}

arrayOfThreeUnsignedInts::
arrayOfThreeUnsignedInts (const ::xercesc::DOMElement& e,
                          ::xml_schema::flags f,
                          ::xml_schema::container* c)
: ::unsignedIntList (e, f, c)
{
}

arrayOfThreeUnsignedInts::
arrayOfThreeUnsignedInts (const ::xercesc::DOMAttr& a,
                          ::xml_schema::flags f,
                          ::xml_schema::container* c)
: ::unsignedIntList (a, f, c)
{
}

arrayOfThreeUnsignedInts::
arrayOfThreeUnsignedInts (const ::std::string& s,
                          const ::xercesc::DOMElement* e,
                          ::xml_schema::flags f,
                          ::xml_schema::container* c)
: ::unsignedIntList (s, e, f, c)
{
}

arrayOfThreeUnsignedInts* arrayOfThreeUnsignedInts::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class arrayOfThreeUnsignedInts (*this, f, c);
}

arrayOfThreeUnsignedInts::
~arrayOfThreeUnsignedInts ()
{
}

// Cuboid
//

Cuboid::
Cuboid (const llfc_type& llfc,
        const particlePerDimension_type& particlePerDimension,
        const particleVelocity_type& particleVelocity,
        const h_type& h,
        const mass_type& mass,
        const generateNumber_type& generateNumber)
: ::xml_schema::type (),
  llfc_ (llfc, this),
  particlePerDimension_ (particlePerDimension, this),
  particleVelocity_ (particleVelocity, this),
  h_ (h, this),
  mass_ (mass, this),
  generateNumber_ (generateNumber, this)
{
}

Cuboid::
Cuboid (const Cuboid& x,
        ::xml_schema::flags f,
        ::xml_schema::container* c)
: ::xml_schema::type (x, f, c),
  llfc_ (x.llfc_, f, this),
  particlePerDimension_ (x.particlePerDimension_, f, this),
  particleVelocity_ (x.particleVelocity_, f, this),
  h_ (x.h_, f, this),
  mass_ (x.mass_, f, this),
  generateNumber_ (x.generateNumber_, f, this)
{
}

Cuboid::
Cuboid (const ::xercesc::DOMElement& e,
        ::xml_schema::flags f,
        ::xml_schema::container* c)
: ::xml_schema::type (e, f | ::xml_schema::flags::base, c),
  llfc_ (this),
  particlePerDimension_ (this),
  particleVelocity_ (this),
  h_ (this),
  mass_ (this),
  generateNumber_ (this)
{
  if ((f & ::xml_schema::flags::base) == 0)
  {
    ::xsd::cxx::xml::dom::parser< char > p (e, true, false, false);
    this->parse (p, f);
  }
}

void Cuboid::
parse (::xsd::cxx::xml::dom::parser< char >& p,
       ::xml_schema::flags f)
{
  for (; p.more_content (); p.next_content (false))
  {
    const ::xercesc::DOMElement& i (p.cur_element ());
    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (i));

    // llfc
    //
    if (n.name () == "llfc" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< llfc_type > r (
        llfc_traits::create (i, f, this));

      if (!llfc_.present ())
      {
        this->llfc_.set (::std::move (r));
        continue;
      }
    }

    // particlePerDimension
    //
    if (n.name () == "particlePerDimension" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< particlePerDimension_type > r (
        particlePerDimension_traits::create (i, f, this));

      if (!particlePerDimension_.present ())
      {
        this->particlePerDimension_.set (::std::move (r));
        continue;
      }
    }

    // particleVelocity
    //
    if (n.name () == "particleVelocity" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< particleVelocity_type > r (
        particleVelocity_traits::create (i, f, this));

      if (!particleVelocity_.present ())
      {
        this->particleVelocity_.set (::std::move (r));
        continue;
      }
    }

    // h
    //
    if (n.name () == "h" && n.namespace_ ().empty ())
    {
      if (!h_.present ())
      {
        this->h_.set (h_traits::create (i, f, this));
        continue;
      }
    }

    // mass
    //
    if (n.name () == "mass" && n.namespace_ ().empty ())
    {
      if (!mass_.present ())
      {
        this->mass_.set (mass_traits::create (i, f, this));
        continue;
      }
    }

    // generateNumber
    //
    if (n.name () == "generateNumber" && n.namespace_ ().empty ())
    {
      if (!generateNumber_.present ())
      {
        this->generateNumber_.set (generateNumber_traits::create (i, f, this));
        continue;
      }
    }

    break;
  }

  if (!llfc_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "llfc",
      "");
  }

  if (!particlePerDimension_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "particlePerDimension",
      "");
  }

  if (!particleVelocity_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "particleVelocity",
      "");
  }

  if (!h_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "h",
      "");
  }

  if (!mass_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "mass",
      "");
  }

  if (!generateNumber_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "generateNumber",
      "");
  }
}

Cuboid* Cuboid::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class Cuboid (*this, f, c);
}

Cuboid& Cuboid::
operator= (const Cuboid& x)
{
  if (this != &x)
  {
    static_cast< ::xml_schema::type& > (*this) = x;
    this->llfc_ = x.llfc_;
    this->particlePerDimension_ = x.particlePerDimension_;
    this->particleVelocity_ = x.particleVelocity_;
    this->h_ = x.h_;
    this->mass_ = x.mass_;
    this->generateNumber_ = x.generateNumber_;
  }

  return *this;
}

Cuboid::
~Cuboid ()
{
}

// Sphere
//

Sphere::
Sphere (const centerCoordinates_type& centerCoordinates,
        const initialVelocity_type& initialVelocity,
        const radius_type& radius,
        const distance_type& distance,
        const mass_type& mass)
: ::xml_schema::type (),
  centerCoordinates_ (centerCoordinates, this),
  initialVelocity_ (initialVelocity, this),
  radius_ (radius, this),
  distance_ (distance, this),
  mass_ (mass, this)
{
}

Sphere::
Sphere (const Sphere& x,
        ::xml_schema::flags f,
        ::xml_schema::container* c)
: ::xml_schema::type (x, f, c),
  centerCoordinates_ (x.centerCoordinates_, f, this),
  initialVelocity_ (x.initialVelocity_, f, this),
  radius_ (x.radius_, f, this),
  distance_ (x.distance_, f, this),
  mass_ (x.mass_, f, this)
{
}

Sphere::
Sphere (const ::xercesc::DOMElement& e,
        ::xml_schema::flags f,
        ::xml_schema::container* c)
: ::xml_schema::type (e, f | ::xml_schema::flags::base, c),
  centerCoordinates_ (this),
  initialVelocity_ (this),
  radius_ (this),
  distance_ (this),
  mass_ (this)
{
  if ((f & ::xml_schema::flags::base) == 0)
  {
    ::xsd::cxx::xml::dom::parser< char > p (e, true, false, false);
    this->parse (p, f);
  }
}

void Sphere::
parse (::xsd::cxx::xml::dom::parser< char >& p,
       ::xml_schema::flags f)
{
  for (; p.more_content (); p.next_content (false))
  {
    const ::xercesc::DOMElement& i (p.cur_element ());
    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (i));

    // centerCoordinates
    //
    if (n.name () == "centerCoordinates" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< centerCoordinates_type > r (
        centerCoordinates_traits::create (i, f, this));

      if (!centerCoordinates_.present ())
      {
        this->centerCoordinates_.set (::std::move (r));
        continue;
      }
    }

    // initialVelocity
    //
    if (n.name () == "initialVelocity" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< initialVelocity_type > r (
        initialVelocity_traits::create (i, f, this));

      if (!initialVelocity_.present ())
      {
        this->initialVelocity_.set (::std::move (r));
        continue;
      }
    }

    // radius
    //
    if (n.name () == "radius" && n.namespace_ ().empty ())
    {
      if (!radius_.present ())
      {
        this->radius_.set (radius_traits::create (i, f, this));
        continue;
      }
    }

    // distance
    //
    if (n.name () == "distance" && n.namespace_ ().empty ())
    {
      if (!distance_.present ())
      {
        this->distance_.set (distance_traits::create (i, f, this));
        continue;
      }
    }

    // mass
    //
    if (n.name () == "mass" && n.namespace_ ().empty ())
    {
      if (!mass_.present ())
      {
        this->mass_.set (mass_traits::create (i, f, this));
        continue;
      }
    }

    break;
  }

  if (!centerCoordinates_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "centerCoordinates",
      "");
  }

  if (!initialVelocity_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "initialVelocity",
      "");
  }

  if (!radius_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "radius",
      "");
  }

  if (!distance_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "distance",
      "");
  }

  if (!mass_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "mass",
      "");
  }
}

Sphere* Sphere::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class Sphere (*this, f, c);
}

Sphere& Sphere::
operator= (const Sphere& x)
{
  if (this != &x)
  {
    static_cast< ::xml_schema::type& > (*this) = x;
    this->centerCoordinates_ = x.centerCoordinates_;
    this->initialVelocity_ = x.initialVelocity_;
    this->radius_ = x.radius_;
    this->distance_ = x.distance_;
    this->mass_ = x.mass_;
  }

  return *this;
}

Sphere::
~Sphere ()
{
}

// simulationConfig
//

simulationConfig::
simulationConfig (const baseName_type& baseName,
                  const writeFrequency_type& writeFrequency,
                  const t_end_type& t_end,
                  const delta_t_type& delta_t,
                  const logLevel_type& logLevel,
                  const domainSize_type& domainSize,
                  const containerType_type& containerType,
                  const boundaries_type& boundaries,
                  const cutOffRadius_type& cutOffRadius)
: ::xml_schema::type (),
  baseName_ (baseName, this),
  writeFrequency_ (writeFrequency, this),
  t_end_ (t_end, this),
  delta_t_ (delta_t, this),
  logLevel_ (logLevel, this),
  domainSize_ (domainSize, this),
  containerType_ (containerType, this),
  boundaries_ (boundaries, this),
  cutOffRadius_ (cutOffRadius, this),
  Cuboid_ (this),
  Sphere_ (this)
{
}

simulationConfig::
simulationConfig (const simulationConfig& x,
                  ::xml_schema::flags f,
                  ::xml_schema::container* c)
: ::xml_schema::type (x, f, c),
  baseName_ (x.baseName_, f, this),
  writeFrequency_ (x.writeFrequency_, f, this),
  t_end_ (x.t_end_, f, this),
  delta_t_ (x.delta_t_, f, this),
  logLevel_ (x.logLevel_, f, this),
  domainSize_ (x.domainSize_, f, this),
  containerType_ (x.containerType_, f, this),
  boundaries_ (x.boundaries_, f, this),
  cutOffRadius_ (x.cutOffRadius_, f, this),
  Cuboid_ (x.Cuboid_, f, this),
  Sphere_ (x.Sphere_, f, this)
{
}

simulationConfig::
simulationConfig (const ::xercesc::DOMElement& e,
                  ::xml_schema::flags f,
                  ::xml_schema::container* c)
: ::xml_schema::type (e, f | ::xml_schema::flags::base, c),
  baseName_ (this),
  writeFrequency_ (this),
  t_end_ (this),
  delta_t_ (this),
  logLevel_ (this),
  domainSize_ (this),
  containerType_ (this),
  boundaries_ (this),
  cutOffRadius_ (this),
  Cuboid_ (this),
  Sphere_ (this)
{
  if ((f & ::xml_schema::flags::base) == 0)
  {
    ::xsd::cxx::xml::dom::parser< char > p (e, true, false, false);
    this->parse (p, f);
  }
}

void simulationConfig::
parse (::xsd::cxx::xml::dom::parser< char >& p,
       ::xml_schema::flags f)
{
  for (; p.more_content (); p.next_content (false))
  {
    const ::xercesc::DOMElement& i (p.cur_element ());
    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (i));

    // baseName
    //
    if (n.name () == "baseName" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< baseName_type > r (
        baseName_traits::create (i, f, this));

      if (!baseName_.present ())
      {
        this->baseName_.set (::std::move (r));
        continue;
      }
    }

    // writeFrequency
    //
    if (n.name () == "writeFrequency" && n.namespace_ ().empty ())
    {
      if (!writeFrequency_.present ())
      {
        this->writeFrequency_.set (writeFrequency_traits::create (i, f, this));
        continue;
      }
    }

    // t_end
    //
    if (n.name () == "t_end" && n.namespace_ ().empty ())
    {
      if (!t_end_.present ())
      {
        this->t_end_.set (t_end_traits::create (i, f, this));
        continue;
      }
    }

    // delta_t
    //
    if (n.name () == "delta_t" && n.namespace_ ().empty ())
    {
      if (!delta_t_.present ())
      {
        this->delta_t_.set (delta_t_traits::create (i, f, this));
        continue;
      }
    }

    // logLevel
    //
    if (n.name () == "logLevel" && n.namespace_ ().empty ())
    {
      if (!logLevel_.present ())
      {
        this->logLevel_.set (logLevel_traits::create (i, f, this));
        continue;
      }
    }

    // domainSize
    //
    if (n.name () == "domainSize" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< domainSize_type > r (
        domainSize_traits::create (i, f, this));

      if (!domainSize_.present ())
      {
        this->domainSize_.set (::std::move (r));
        continue;
      }
    }

    // containerType
    //
    if (n.name () == "containerType" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< containerType_type > r (
        containerType_traits::create (i, f, this));

      if (!containerType_.present ())
      {
        this->containerType_.set (::std::move (r));
        continue;
      }
    }

    // boundaries
    //
    if (n.name () == "boundaries" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< boundaries_type > r (
        boundaries_traits::create (i, f, this));

      if (!boundaries_.present ())
      {
        this->boundaries_.set (::std::move (r));
        continue;
      }
    }

    // cutOffRadius
    //
    if (n.name () == "cutOffRadius" && n.namespace_ ().empty ())
    {
      if (!cutOffRadius_.present ())
      {
        this->cutOffRadius_.set (cutOffRadius_traits::create (i, f, this));
        continue;
      }
    }

    // Cuboid
    //
    if (n.name () == "Cuboid" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< Cuboid_type > r (
        Cuboid_traits::create (i, f, this));

      this->Cuboid_.push_back (::std::move (r));
      continue;
    }

    // Sphere
    //
    if (n.name () == "Sphere" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< Sphere_type > r (
        Sphere_traits::create (i, f, this));

      this->Sphere_.push_back (::std::move (r));
      continue;
    }

    break;
  }

  if (!baseName_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "baseName",
      "");
  }

  if (!writeFrequency_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "writeFrequency",
      "");
  }

  if (!t_end_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "t_end",
      "");
  }

  if (!delta_t_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "delta_t",
      "");
  }

  if (!logLevel_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "logLevel",
      "");
  }

  if (!domainSize_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "domainSize",
      "");
  }

  if (!containerType_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "containerType",
      "");
  }

  if (!boundaries_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "boundaries",
      "");
  }

  if (!cutOffRadius_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "cutOffRadius",
      "");
  }
}

simulationConfig* simulationConfig::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class simulationConfig (*this, f, c);
}

simulationConfig& simulationConfig::
operator= (const simulationConfig& x)
{
  if (this != &x)
  {
    static_cast< ::xml_schema::type& > (*this) = x;
    this->baseName_ = x.baseName_;
    this->writeFrequency_ = x.writeFrequency_;
    this->t_end_ = x.t_end_;
    this->delta_t_ = x.delta_t_;
    this->logLevel_ = x.logLevel_;
    this->domainSize_ = x.domainSize_;
    this->containerType_ = x.containerType_;
    this->boundaries_ = x.boundaries_;
    this->cutOffRadius_ = x.cutOffRadius_;
    this->Cuboid_ = x.Cuboid_;
    this->Sphere_ = x.Sphere_;
  }

  return *this;
}

simulationConfig::
~simulationConfig ()
{
}

#include <istream>
#include <xsd/cxx/xml/sax/std-input-source.hxx>
#include <xsd/cxx/tree/error-handler.hxx>

::std::unique_ptr< ::simulationConfig >
Configuration (const ::std::string& u,
               ::xml_schema::flags f,
               const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xsd::cxx::tree::error_handler< char > h;

  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      u, h, p, f));

  h.throw_if_failed< ::xsd::cxx::tree::parsing< char > > ();

  return ::std::unique_ptr< ::simulationConfig > (
    ::Configuration (
      std::move (d), f | ::xml_schema::flags::own_dom, p));
}

::std::unique_ptr< ::simulationConfig >
Configuration (const ::std::string& u,
               ::xml_schema::error_handler& h,
               ::xml_schema::flags f,
               const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      u, h, p, f));

  if (!d.get ())
    throw ::xsd::cxx::tree::parsing< char > ();

  return ::std::unique_ptr< ::simulationConfig > (
    ::Configuration (
      std::move (d), f | ::xml_schema::flags::own_dom, p));
}

::std::unique_ptr< ::simulationConfig >
Configuration (const ::std::string& u,
               ::xercesc::DOMErrorHandler& h,
               ::xml_schema::flags f,
               const ::xml_schema::properties& p)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      u, h, p, f));

  if (!d.get ())
    throw ::xsd::cxx::tree::parsing< char > ();

  return ::std::unique_ptr< ::simulationConfig > (
    ::Configuration (
      std::move (d), f | ::xml_schema::flags::own_dom, p));
}

::std::unique_ptr< ::simulationConfig >
Configuration (::std::istream& is,
               ::xml_schema::flags f,
               const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xsd::cxx::xml::sax::std_input_source isrc (is);
  return ::Configuration (isrc, f, p);
}

::std::unique_ptr< ::simulationConfig >
Configuration (::std::istream& is,
               ::xml_schema::error_handler& h,
               ::xml_schema::flags f,
               const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xsd::cxx::xml::sax::std_input_source isrc (is);
  return ::Configuration (isrc, h, f, p);
}

::std::unique_ptr< ::simulationConfig >
Configuration (::std::istream& is,
               ::xercesc::DOMErrorHandler& h,
               ::xml_schema::flags f,
               const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::sax::std_input_source isrc (is);
  return ::Configuration (isrc, h, f, p);
}

::std::unique_ptr< ::simulationConfig >
Configuration (::std::istream& is,
               const ::std::string& sid,
               ::xml_schema::flags f,
               const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xsd::cxx::xml::sax::std_input_source isrc (is, sid);
  return ::Configuration (isrc, f, p);
}

::std::unique_ptr< ::simulationConfig >
Configuration (::std::istream& is,
               const ::std::string& sid,
               ::xml_schema::error_handler& h,
               ::xml_schema::flags f,
               const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xsd::cxx::xml::sax::std_input_source isrc (is, sid);
  return ::Configuration (isrc, h, f, p);
}

::std::unique_ptr< ::simulationConfig >
Configuration (::std::istream& is,
               const ::std::string& sid,
               ::xercesc::DOMErrorHandler& h,
               ::xml_schema::flags f,
               const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::sax::std_input_source isrc (is, sid);
  return ::Configuration (isrc, h, f, p);
}

::std::unique_ptr< ::simulationConfig >
Configuration (::xercesc::InputSource& i,
               ::xml_schema::flags f,
               const ::xml_schema::properties& p)
{
  ::xsd::cxx::tree::error_handler< char > h;

  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      i, h, p, f));

  h.throw_if_failed< ::xsd::cxx::tree::parsing< char > > ();

  return ::std::unique_ptr< ::simulationConfig > (
    ::Configuration (
      std::move (d), f | ::xml_schema::flags::own_dom, p));
}

::std::unique_ptr< ::simulationConfig >
Configuration (::xercesc::InputSource& i,
               ::xml_schema::error_handler& h,
               ::xml_schema::flags f,
               const ::xml_schema::properties& p)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      i, h, p, f));

  if (!d.get ())
    throw ::xsd::cxx::tree::parsing< char > ();

  return ::std::unique_ptr< ::simulationConfig > (
    ::Configuration (
      std::move (d), f | ::xml_schema::flags::own_dom, p));
}

::std::unique_ptr< ::simulationConfig >
Configuration (::xercesc::InputSource& i,
               ::xercesc::DOMErrorHandler& h,
               ::xml_schema::flags f,
               const ::xml_schema::properties& p)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      i, h, p, f));

  if (!d.get ())
    throw ::xsd::cxx::tree::parsing< char > ();

  return ::std::unique_ptr< ::simulationConfig > (
    ::Configuration (
      std::move (d), f | ::xml_schema::flags::own_dom, p));
}

::std::unique_ptr< ::simulationConfig >
Configuration (const ::xercesc::DOMDocument& doc,
               ::xml_schema::flags f,
               const ::xml_schema::properties& p)
{
  if (f & ::xml_schema::flags::keep_dom)
  {
    ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
      static_cast< ::xercesc::DOMDocument* > (doc.cloneNode (true)));

    return ::std::unique_ptr< ::simulationConfig > (
      ::Configuration (
        std::move (d), f | ::xml_schema::flags::own_dom, p));
  }

  const ::xercesc::DOMElement& e (*doc.getDocumentElement ());
  const ::xsd::cxx::xml::qualified_name< char > n (
    ::xsd::cxx::xml::dom::name< char > (e));

  if (n.name () == "Configuration" &&
      n.namespace_ () == "")
  {
    ::std::unique_ptr< ::simulationConfig > r (
      ::xsd::cxx::tree::traits< ::simulationConfig, char >::create (
        e, f, 0));
    return r;
  }

  throw ::xsd::cxx::tree::unexpected_element < char > (
    n.name (),
    n.namespace_ (),
    "Configuration",
    "");
}

::std::unique_ptr< ::simulationConfig >
Configuration (::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d,
               ::xml_schema::flags f,
               const ::xml_schema::properties&)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > c (
    ((f & ::xml_schema::flags::keep_dom) &&
     !(f & ::xml_schema::flags::own_dom))
    ? static_cast< ::xercesc::DOMDocument* > (d->cloneNode (true))
    : 0);

  ::xercesc::DOMDocument& doc (c.get () ? *c : *d);
  const ::xercesc::DOMElement& e (*doc.getDocumentElement ());

  const ::xsd::cxx::xml::qualified_name< char > n (
    ::xsd::cxx::xml::dom::name< char > (e));

  if (f & ::xml_schema::flags::keep_dom)
    doc.setUserData (::xml_schema::dom::tree_node_key,
                     (c.get () ? &c : &d),
                     0);

  if (n.name () == "Configuration" &&
      n.namespace_ () == "")
  {
    ::std::unique_ptr< ::simulationConfig > r (
      ::xsd::cxx::tree::traits< ::simulationConfig, char >::create (
        e, f, 0));
    return r;
  }

  throw ::xsd::cxx::tree::unexpected_element < char > (
    n.name (),
    n.namespace_ (),
    "Configuration",
    "");
}

#include <xsd/cxx/post.hxx>

// Begin epilogue.
//
//
// End epilogue.


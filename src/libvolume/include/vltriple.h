/***************************************************************************
                          vltriple.h  -  description
                             -------------------
    begin                : Thu Apr 1 2004
    copyright            : (C) 2004 by Jose Ignacio Garzon
    email                :
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/


   #ifndef _vlTriple_h
#define _vlTriple_h

#include <stdio.h>

#include <ostream>
#include <math.h>

#include "vlconstants.h"

/**
 * This class repesents a datatype which stores 3 entries of the same type. The
 * entries can be accessed using x(), y() and z() functions.
 *
 * @author Sarang Lakare <sarang@users.sf.net>
 */
template <class T>
class vlTriple
{
public:
  /// default constructor
  vlTriple<T>(T const & x = 0, T const & y = 0, T const & z = 0);

  /// default destructor - no virtual.. it adds 4 bytes everytime
  ~vlTriple<T>();

  /// get value of x
  T x() const { return m_x; }

  /// get value of y
  T y() const { return m_y; }

  /// get value of z
  T z() const { return m_z; }

  /// Returns a reference to the vriable x - prefer x() to this
  T const & xRef() const { return m_x; }

  /// Returns a reference to the vriable y - prefer y() to this
  T const & yRef() const { return m_y; }

  /// Returns a reference to the vriable z - prefer z() to this
  T const & zRef() const { return m_z; }

  /// override for = operator
  template <class S>
  vlTriple<T> & operator = (vlTriple<S> const & other);

  /// set value of x
  template <class S>
  void x(const S & value) { m_x = static_cast<S>(value); }

  /// set value of y
  template <class S>
  void y(const S & value) { m_y = static_cast<S>(value); }

  /// set value of z
  template <class S>
  void z(const S & value) { m_z = static_cast<S>(value); }

  /// checks if the two triple are equal
  bool operator == (vlTriple<T> const & other) const;

  /// checks if the two triple are unequal
  bool operator != (vlTriple<T> const & other) const;

  /// multiples triple elements with +1
  vlTriple<T> operator +  () const;

  /// multiples triple elements with -1
  vlTriple<T> operator -  () const;

  /// adds another triple to self and returns the sum
  template <class S>
  vlTriple<T> operator + (vlTriple<S> const & other) const
  {
    return vlTriple<T>(m_x+other.x(), m_y+other.y(), m_z+other.z());
  };

  /// adds a const value to self and returns the sum
  vlTriple<T> operator + (T const & value) const
  {
    return vlTriple<T>(m_x+value, m_y+value, m_z+value);
  };

  /// subtracts another triple from self and returns the result
  template <class S>
  vlTriple<T> operator - (vlTriple<S> const & other) const
  {
    return vlTriple<T>(m_x-other.x(), m_y-other.y(), m_z-other.z());
  };

  /// subtracts a const value from self and returns the result
  template <class S>
  vlTriple<T> operator - (S const & value) const
  {
    return vlTriple<T>(m_x-value, m_y-value, m_z-value);
  };

  /// adds another triple to self
  template <class S>
  vlTriple<T> & operator += (vlTriple<S> const & other)
  {
    m_x+=other.x(); m_y+=other.y(); m_z+=other.z();
    return *this;
  };

  /// subtracts another triple from self
  template <class S>
  vlTriple<T> & operator -= (vlTriple<S> const & other)
  {
    m_x-=other.x(); m_y-=other.y(); m_z-=other.z();
    return *this;
  };

  /// multiples self with other and returns the result
  template <class S>
  vlTriple<T> operator * (S const & value) const
  {
    return vlTriple(m_x*value, m_y*value, m_z*value);
  };

  /// divides other from self and returns the result
  template <class S>
  vlTriple<T> operator / (S const & value) const
  {
    S factor =
    (value == 0.0f) ? 1.0 : (1.0/value);

    return vlTriple<T>(m_x*factor, m_y*factor, m_z*factor);
  };

  /// multiples self with other (self changes)
  template <class S>
  void operator *= (S const & value)
  {
    m_x*=value;m_y*=value;m_z*=value;
  };

  /// divides other from self (self changes)
  template <class S>
  void operator /= (S const & value)
  {
    float factor = (value == 0.0f) ? 1.0 : (1.0f/value);
    m_x = m_x*factor; m_y = m_y*factor; m_z = m_z*factor;
  };


  /// overriding << to enable writing triple to a stream
  friend std::ostream & operator << (std::ostream & os, vlTriple<T> const * const v) {
    return (os << *v);
  }

  /// overriding << to enable writing triple to a stream
  friend std::ostream & operator << (std::ostream & os, vlTriple<T> const & v) {
    return os<<"("<<v.m_x<<","<<v.m_y<<","<<v.m_z<<")";
  }

protected:
  /// member variables storing the 3 entries in the triple
  T m_x, m_y, m_z;
};

// include the inline definitions

/**
 * Default constructor.
 *
 * @param T      data type.
 * @param x      the x value.
 * @param y      the y value.
 * @param z      the z value
 */
template <class T>
vlTriple<T>::vlTriple(T const & x, T const & y, T const & z)
  : m_x(x),
    m_y(y),
    m_z(z)
{

}


/**
 * Default destructor.
 * Note that this is not a virtual function. The reason is that a virtual function
 * adds 4 bytes to the size of the class. Thus, a sizeof(vlTriple<char>) would be 7 bytes
 * if it were virtual.
 *
 * @param T      data type.
 */
template <class T>
vlTriple<T>::~vlTriple()
{

}


/**
 * copy operator for triple.
 *
 * @param T      data type.
 * @param S     Data type of the other triple to be copied
 * @param other  the other triple to copy.
 * @return pointer to self.
 */
template <class T> template <class S>
inline vlTriple<T> & vlTriple<T>::operator =(vlTriple<S> const & other) {
  m_x = other.x();
  m_y = other.y();
  m_z = other.z();
  return (*this);
}


/**
 * checks if the two triples are equal.
 *
 * @param T      data type.
 * @param other  the other triple to compare with.
 * @return true if the two triples are equal.
 */
template <class T>
inline bool vlTriple<T>::operator == (vlTriple<T> const & other) const
{
  return (m_x==other.m_x && m_y==other.m_y && m_z==other.m_z);
}


/**
 * compares two float triples if they are equal.
 *
 * @param other  the other triple to compare with.
 * @return true if the two triples are equal.
 */
template <>
inline bool vlTriple<float>::operator == (vlTriple<float> const & other) const
{
  const float diffDelta(static_cast<float>(0.000001));
  return ((m_x-other.m_x < diffDelta) &&
    (m_y-other.m_y < diffDelta) &&
    (m_z-other.m_z < diffDelta));
}


/**
 * compares two double triples if they are equal.
 *
 * @param other  the other triple to compare with.
 * @return true if the two triples are equal.
 */
template <>
inline bool vlTriple<double>::operator == (vlTriple<double> const & other) const
{
  const double diffDelta(static_cast<double>(0.0000000001));
  return ((m_x-other.m_x < diffDelta) &&
    (m_y-other.m_y < diffDelta) &&
    (m_z-other.m_z < diffDelta));
}


/**
 * checks if the two triples are unequal.
 *
 * @param T      data type.
 * @param other  the other triple to compare with.
 * @return true if the two triples are unequal.
 */
template <class T>
inline bool vlTriple<T>::operator != (vlTriple<T> const & other) const
{

  return (m_x!=other.m_x || m_y!=other.m_y || m_z!=other.m_z);
}


/**
 * compares two float triples if they are unequal.
 *
 * @param other  the other triple to compare with.
 * @return true if the two triples are unequal.
 */
template <>
inline bool vlTriple<float>::operator != (vlTriple<float> const & other) const
{
  const float diffDelta(static_cast<float>(0.000001));
  return ( (m_x-other.m_x > diffDelta) ||
    (m_y-other.m_y > diffDelta) ||
    (m_z-other.m_z > diffDelta) );
}


/**
 * compares two double triples if they are unequal.
 *
 * @param other  the other triple to compare with.
 * @return true if the two triples are unequal.
 */
template <>
inline bool vlTriple<double>::operator != (vlTriple<double> const & other) const
{
  const double diffDelta(static_cast<double>(0.0000000001));
  return ((m_x-other.m_x > diffDelta) ||
    (m_y-other.m_y > diffDelta) ||
    (m_z-other.m_z > diffDelta));
}


/**
 * Multiplies elements with +1.
 *
 * @param T      data type.
 * @return another triple with the result.
 */
template <class T>
inline vlTriple<T> vlTriple<T>::operator +  () const {
  return vlTriple<T>(*this);
}


/**
 * Multiplies elements with -1.
 *
 * @param T      data type
 * @return another triple with the result.
 */
template <class T>
inline vlTriple<T> vlTriple<T>::operator - () const {
  return vlTriple<T>(-m_x, -m_y, -m_z);
}

#endif // _vlTriple_h

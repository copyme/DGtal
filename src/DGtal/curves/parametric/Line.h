/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#pragma once

/**
 * @file Line.h
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, A3SI, France
 *
 * @date 2014/11/05
 *
 * Header file for module Line.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(Line_RECURSES)
#error Recursive header files inclusion detected in Line.h
#else // defined(Line_RECURSES)
/** Prevents recursive inclusion of headers. */
#define Line_RECURSES

#if !defined Line_h
/** Prevents repeated inclusion of headers. */
#define Line_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

/////////////////////////////////////////////////////////////////////////////
// class Line
/**
 * Description of class 'Line' <p>
 * \brief Aim:
 */
template <typename TSpace>
class Line
{
    // ----------------------- Standard services ------------------------------
public:

    typedef TSpace Space;
    typedef typename TSpace::RealPoint RealPoint;
    typedef typename TSpace::Point Point;

    Line ( bool, bool, bool );

    /**
     * Destructor.
     */
    ~Line() { }

    // ----------------------- Interface --------------------------------------
public:

    /**
     * @param t any angle between 0 and k*Pi.
     *
     * @return the vector (x(t),y(t), z(t))
     */
    RealPoint x ( double t ) const;

    /**
     * @param t any angle between 0 and k*Pi.
     *
     * @return the vector (x(t)',y(t)', z(t)')
     */
    RealPoint xp ( double ) const;

    /**
     * @brief inverse function of x
     * @param p = x(t)
     * @return t
     */
    double f ( const RealPoint & p ) const;

    /**
     * @brief inverse function of y
     * @param p = x(t)[
     * @return t
     */
    double g ( const RealPoint & p ) const;

    /**
     * @brief inverse function of z
     * @param p = x(t)
     * @return t
     */
    double h ( const RealPoint & p ) const;

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    // ------------------------- Protected Datas ------------------------------
private:
    // ------------------------- Private Datas --------------------------------
private:
    bool mask[3];
    // ------------------------- Hidden services ------------------------------
protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    Line();

private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    Line ( const Line & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    Line & operator= ( const Line & other );

    // ------------------------- Internals ------------------------------------
private:

}; // end of class Line


/**
 * Overloads 'operator<<' for displaying objects of class 'Line'.
 * @param out the output stream where the object is written.
 * @param object the object of class 'Line' to write.
 * @return the output stream after the writing.
 */
template <typename TSpace>
std::ostream&
operator<< ( std::ostream & out, const Line < TSpace > & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#if !defined(BUILD_INLINE)
#include "DGtal/curves/parametric/Line.ih"
#endif


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined Line_h

#undef Line_RECURSES
#endif // else defined(Line_RECURSES)

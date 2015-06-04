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
 * @file EllipticHelix.h
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, A3SI, France
 *
 * @date 2014/10/01
 *
 * Header file for module EllipticHelix.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(EllipticHelix_RECURSES)
#error Recursive header files inclusion detected in EllipticHelix.h
#else // defined(EllipticHelix_RECURSES)
/** Prevents recursive inclusion of headers. */
#define EllipticHelix_RECURSES

#if !defined EllipticHelix_h
/** Prevents repeated inclusion of headers. */
#define EllipticHelix_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <cmath>
#include "DGtal/base/Common.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

/////////////////////////////////////////////////////////////////////////////
// class EllipticHelix
/**
 * Description of class 'EllipticHelix' <p>
 * \brief Aim:
 */
template <typename TSpace>
class EllipticHelix
{
    // ----------------------- Standard services ------------------------------
public:
    typedef TSpace Space;
    typedef typename TSpace::RealPoint RealPoint;
    typedef typename TSpace::Point Point;

    /**
     * Destructor.
     */
    ~EllipticHelix() {}
    /**
     * Constructor.
     * @param rr small radius of the helix
     * @param rl big radius of the helix
     * @param bb distance between each turn.
     */
    EllipticHelix ( double rr, double rl,  double bb );

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
    RealPoint xp ( double t ) const;

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
    double r1, r2, b;
    // ------------------------- Hidden services ------------------------------
protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    EllipticHelix();

private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    EllipticHelix ( const EllipticHelix & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    EllipticHelix & operator= ( const EllipticHelix & other );

    // ------------------------- Internals ------------------------------------
private:

}; // end of class EllipticHelix


/**
 * Overloads 'operator<<' for displaying objects of class 'EllipticHelix'.
 * @param out the output stream where the object is written.
 * @param object the object of class 'EllipticHelix' to write.
 * @return the output stream after the writing.
 */
template <typename T>
std::ostream&
operator<< ( std::ostream & out, const EllipticHelix<T> & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#if !defined(BUILD_INLINE)
#include "DGtal/curves/parametric/EllipticHelix.ih"
#endif


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined EllipticHelix_h

#undef EllipticHelix_RECURSES
#endif // else defined(EllipticHelix_RECURSES)

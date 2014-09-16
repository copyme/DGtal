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
 * @file QuasiAffineTransformation2D.h
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, France
 *
 * @date 2014/09/16
 *
 * Header file for module QuasiAffineTransformation2D.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(QuasiAffineTransformation2D_RECURSES)
#error Recursive header files inclusion detected in QuasiAffineTransformation2D.h
#else // defined(QuasiAffineTransformation2D_RECURSES)
/** Prevents recursive inclusion of headers. */
#define QuasiAffineTransformation2D_RECURSES

#if !defined QuasiAffineTransformation2D_h
/** Prevents repeated inclusion of headers. */
#define QuasiAffineTransformation2D_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/math/linalg/SimpleMatrix.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

/////////////////////////////////////////////////////////////////////////////
// class QuasiAffineTransformation2D
/**
 * Description of class 'QuasiAffineTransformation2D' <p>
 * \brief Aim:
 */
class QuasiAffineTransformation2D
{
    // ----------------------- Standard services ------------------------------
public:

    /**
     * Destructor.
     */
    ~QuasiAffineTransformation2D();

    // ----------------------- Interface --------------------------------------
public:

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
    typedef SimpleMatrix < int, 2, 2 > Matrix;
    // ------------------------- Private Datas --------------------------------
private:
    Matrix m_matrix;
    int m_omega;
    Z2i::Vector m_vector;

    int a0, b0, c0, d0, e0, f0;
    int a1, b1, c1, d1;
    int b10, b11;
    Matrix H, H0;
    int alpha0, alpha1, beta1;
    Z2i::Vector U0, U1;
    std::vector<Z2i::Point> pavings;
    std::vector<Z2i::Point> pavingsRemainder; // remainders of the canonical paving points

    // ------------------------- Hidden services ------------------------------
protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    QuasiAffineTransformation2D();

private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    QuasiAffineTransformation2D ( const QuasiAffineTransformation2D & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    QuasiAffineTransformation2D & operator= ( const QuasiAffineTransformation2D & other );

    // ------------------------- Internals ------------------------------------
private:

}; // end of class QuasiAffineTransformation2D


/**
 * Overloads 'operator<<' for displaying objects of class 'QuasiAffineTransformation2D'.
 * @param out the output stream where the object is written.
 * @param object the object of class 'QuasiAffineTransformation2D' to write.
 * @return the output stream after the writing.
 */
std::ostream&
operator<< ( std::ostream & out, const QuasiAffineTransformation2D & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#if !defined(BUILD_INLINE)
#include "DGtal/images/QuasiAffineTransformation2D.ih"
#endif


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined QuasiAffineTransformation2D_h

#undef QuasiAffineTransformation2D_RECURSES
#endif // else defined(QuasiAffineTransformation2D_RECURSES)

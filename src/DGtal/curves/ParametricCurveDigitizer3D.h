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
 * @file CurveDigitizer.h
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, A3SI, France
 *
 * @date 2014/09/26
 *
 * Header file for module CurveDigitizer.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(CurveDigitizer_RECURSES)
#error Recursive header files inclusion detected in CurveDigitizer.h
#else // defined(CurveDigitizer_RECURSES)
/** Prevents recursive inclusion of headers. */
#define CurveDigitizer_RECURSES

#if !defined CurveDigitizer_h
/** Prevents repeated inclusion of headers. */
#define CurveDigitizer_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
#include <vector>
#include "DGtal/curves/parametric/C3DParametricCurve.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

/////////////////////////////////////////////////////////////////////////////
// class CurveDigitizer
/**
 * Description of class 'CurveDigitizer' <p>
 * \brief Aim:
 */
template <typename TParametricCurve>
class ParametricCurveDigitizer3D
{
    BOOST_CONCEPT_ASSERT(( concepts::C3DParametricCurve < TParametricCurve > ));
    // ----------------------- Standard services ------------------------------
public:
    typedef std::vector<typename TParametricCurve::Space::Point> DigitalCurve;
    typedef typename TParametricCurve::Space::Point Point;
    typedef typename TParametricCurve::Space::RealPoint RealPoint;
    typedef typename TParametricCurve::Space::RealVector RealVector;

    ParametricCurveDigitizer3D();
    /**
     * Destructor.
     */
    ~ParametricCurveDigitizer3D() {}

    // ----------------------- Interface --------------------------------------
public:
    void attach ( const TParametricCurve & t_curve );
    void init ( DigitalCurve & digiCurve );

    void digitize ( const double & tmin, double & tmax, double step );

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

    // ------------------------- Private Datas --------------------------------
private:
    // ------------------------- Protected Datas ------------------------------
protected:

    DigitalCurve * digitalCurve;
    const TParametricCurve * curve;
    bool adaptive;
    unsigned int attemptNumber;
    bool is26Connected ( const Point &x, const Point &y );
    // ------------------------- Hidden services ------------------------------
protected:
  virtual unsigned char findMainAxis ( const char & blockAxis, const double & i );
  virtual Point CurvePoint ( const char & mainAxis, const double & t, double & t_time );
private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    ParametricCurveDigitizer3D ( const ParametricCurveDigitizer3D & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    ParametricCurveDigitizer3D & operator= ( const ParametricCurveDigitizer3D & other );

    // ------------------------- Internals ------------------------------------

}; // end of class CurveDigitizer


/**
 * Overloads 'operator<<' for displaying objects of class 'CurveDigitizer'.
 * @param out the output stream where the object is written.
 * @param object the object of class 'CurveDigitizer' to write.
 * @return the output stream after the writing.
 */
template <typename T>
std::ostream&
operator<< ( std::ostream & out, const ParametricCurveDigitizer3D<T> & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#if !defined(BUILD_INLINE)
#include "DGtal/curves/ParametricCurveDigitizer3D.ih"
#endif


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined CurveDigitizer_h

#undef CurveDigitizer_RECURSES
#endif // else defined(CurveDigitizer_RECURSES)

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

/**
 * @file testLambdaMST3D.cpp
 * @ingroup Tests
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, A3SI, France
 *
 * @date 2014/10/03
 *
 * Functions for testing class LambdaMST3D.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <vector>
#include "DGtal/base/Common.h"
#include "ConfigTest.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/geometry/curves/Naive3DDSSComputer.h"
#include "DGtal/geometry/curves/LambdaMST3D.h"
#include "DGtal/curves/ParametricCurveDigitizer3D.h"
#include "DGtal/curves/parametric/EllipticHelix.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"

#ifdef __GNUC__
   #ifndef NDEBUG
      #include <fenv.h>
   #endif
#endif
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class LambdaMST3D.
///////////////////////////////////////////////////////////////////////////////

class testLambdaMST3D
{
    typedef EllipticHelix < Z3i::Space > MyHelix;
    typedef DGtal::Z3i::Point Point;
    typedef std::vector < Point > Range;
    typedef Range::const_iterator ConstIterator;
    typedef ParametricCurveDigitizer3D < MyHelix >  Digitizer;
    typedef ParametricCurveDigitizer3D < MyHelix >::DigitalCurve MyDigitalCurve;
    typedef Naive3DDSSComputer < ConstIterator, int, 8 > SegmentComputer;
    typedef SaturatedSegmentation<SegmentComputer> Segmentation;
private:
    Digitizer digitizer;
    MyHelix helix;
    MyDigitalCurve digitalCurve;
public:
    testLambdaMST3D () : helix ( 10, 10, 20 ) {
        digitizer.attach ( helix );
        digitizer.init ( digitalCurve );
	double maxt = 2 * M_PI;
        digitizer.digitize ( 0, maxt, 0.001 );
    }
    bool lambda64ByPoint ()
    {
        Segmentation segmenter ( digitalCurve.cbegin(), digitalCurve.cend(), SegmentComputer() );
        LambdaMST3D < Segmentation > lmst64;
        lmst64.attach ( segmenter );
        lmst64.init ( digitalCurve.begin(), digitalCurve.end() );
        for ( unsigned int i = 0; i < digitalCurve.size(); i++ )
        {
            lmst64.eval ( digitalCurve[i] );
        }
        return true;
    }
    bool lambda64()
    {
        Segmentation segmenter ( digitalCurve.cbegin(), digitalCurve.cend(), SegmentComputer() );
        LambdaMST3D < Segmentation > lmst64;
        lmst64.attach ( segmenter );
        lmst64.init ( digitalCurve.begin(), digitalCurve.end() );
        std::vector < LambdaMST3D < Segmentation >::RealVector > tangent;
        lmst64.eval ( tangent );
        return true;
    }
    bool lambda64Both()
    {
        Segmentation segmenter ( digitalCurve.cbegin(), digitalCurve.cend(), SegmentComputer() );
        LambdaMST3D < Segmentation > lmst64;
        lmst64.attach ( segmenter );
        lmst64.init ( digitalCurve.begin(), digitalCurve.end() );
        std::vector < LambdaMST3D < Segmentation >::RealVector > tangent;
        lmst64.eval ( tangent );
        for ( unsigned int i = 0; i < digitalCurve.size(); i++ )
        {
            if ( lmst64.eval ( digitalCurve[i] ) != tangent[i] )
                return false;
        }
        return true;
    }
};


///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

int main( int , char**  )
{
#ifdef __GNUC__
   #ifndef NDEBUG
    fetestexcept ( FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW );
   #endif
#endif
    bool res = true;
    testLambdaMST3D testLMST;
    trace.beginBlock ( "Testing LambdaMST3D" );
        trace.beginBlock ( "Testing point only calculation" );
          res &= testLMST.lambda64ByPoint();
        trace.endBlock();
        trace.beginBlock ( "Testing calculation for whole curve" );
           res &= testLMST.lambda64();
        trace.endBlock();
        trace.beginBlock ( "Testing values obtined from the both methods." );
           res &= testLMST.lambda64Both();
        trace.endBlock();
    trace.emphase() << ( res ? "Passed." : "Error." ) << endl;
    trace.endBlock();
    return res ? 0 : 1;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

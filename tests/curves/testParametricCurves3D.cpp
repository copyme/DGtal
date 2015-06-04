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
 * @file testParametricCurves3D.cpp
 * @ingroup Tests
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, A3SI, France
 *
 * @date 2014/09/26
 *
 * Functions for testing class ParametricCurveDigitizer3D.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"
#include "ConfigTest.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/curves/parametric/EllipticHelix.h"
#include "DGtal/curves/parametric/Line.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/curves/ParametricCurveDigitizer3D.h"
#include "DGtal/curves/parametric/DecoratorCurveTransformation.h"
#include "DGtal/images/RigidTransformation3D.h"
#include "DGtal/io/writers/PGMWriter.h"
#include <DGtal/io/boards/Board3D.h>
#include <functional>
#ifdef __GNUC__
   #ifndef NDEBUG
      #include <fenv.h>
   #endif
#endif
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace Z3i;
using namespace functors;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class ParametricCurveDigitizer3D.
///////////////////////////////////////////////////////////////////////////////
/**
 * Example of a test. To be completed.
 *
 */

class testParametricCurves3D
{
    typedef EllipticHelix < Z3i::Space > MyEllipticHelix;
    typedef EllipticHelix < Z3i::Space > MyHelix;
    typedef EllipticHelix < Z3i::Space > MyDisc;
    typedef Line < Z3i::Space > MyLine;
    typedef Board3D<> MyBoard;
private:
    MyEllipticHelix ellipticHelix;
    MyHelix helix;
    MyDisc disc;
    MyBoard board;
    MyLine line;
    template <typename TDCurve >
    void project ( const TDCurve & curve, const char * filename, const TDCurve * errors = 0 )
    {

        board.clear();
        for ( unsigned int i = 0; i < curve.size(); i++ )
        {
            board.setFillColor ( Color ( 0, 0, 255, 255 ) );
            board << SetMode3D(curve[i].className(), "");
            board << curve[i];

            board << SetMode3D(curve[i].className(), "PavingWired");
            board << curve[i];
//             board.setFillColor ( Color ( 255, 0, 0, 255 ) );
//             board << SetMode3D(curve[i].className(), "");
//             board << Z3i::Point ( curve.at(i)[0], curve.at(i)[1], -50 );
// 
//             board << SetMode3D ( curve[i].className(), "PavingWired");
//             board << Z3i::Point ( curve.at(i)[0], curve.at(i)[1], -50);
// 
//             board.setFillColor ( Color (  0, 255, 0 ) );
//             board << SetMode3D(curve[i].className(), "");
//             board << Z3i::Point ( curve.at(i)[0], -50, curve.at(i)[2] );
// 
//             board << SetMode3D ( curve[i].className(), "PavingWired");
//             board << Z3i::Point ( curve.at(i)[0], -50, curve.at(i)[2] );
// 
//             board.setFillColor ( Color (  255, 255, 0 ) );
//             board << SetMode3D(curve[i].className(), "");
//             board << Z3i::Point ( -50, curve.at(i)[1], curve.at(i)[2] );
// 
//             board << SetMode3D ( curve[i].className(), "PavingWired");
//             board << Z3i::Point ( -50, curve.at(i)[1], curve.at(i)[2] );
	}
	if (errors)
	{
	  for ( unsigned int i = 0; i < errors->size(); i++ )
	  {
	    board.setFillColor ( Color ( 255, 0, 0, 255 ) );
	    board << SetMode3D(errors->at(i).className(), "");
	    board << errors->at(i);
	    
	    board << SetMode3D(errors->at(i).className(), "PavingWired");
	    board << errors->at(i);
	  }
	}
	board.saveOBJ ( filename );
    }

public:
    testParametricCurves3D() : ellipticHelix ( 2, 36, 2 ), helix ( 6, 16, 10 ), disc ( 3 * 10 , 8 * 10, 0 ), line (true, true, true) { }

    bool testDigitizeEllipticHelix ()
    {
      typedef ParametricCurveDigitizer3D < MyEllipticHelix >  Digitizer;
      typedef ParametricCurveDigitizer3D < MyEllipticHelix >::DigitalCurve MyDigitalCurve;
      Digitizer digitizer;
      MyDigitalCurve digitalCurve;
      digitizer.attach ( ellipticHelix );
      digitizer.init ( digitalCurve );
      double maxt = 2. * M_PI;
      digitizer.digitize ( 0, maxt, 0.1 );
      project < MyDigitalCurve > ( digitalCurve, "EllipticHelix" );
      return true;
    }
    
    bool testDigitizeHelix ()
    {
      typedef ParametricCurveDigitizer3D < MyHelix >  Digitizer;
      typedef ParametricCurveDigitizer3D < MyHelix >::DigitalCurve MyDigitalCurve;
      Digitizer digitizer;
      MyDigitalCurve digitalCurve;
      digitizer.attach ( helix );
      digitizer.init ( digitalCurve );
      double maxt = 2. * M_PI + 0.01;
      digitizer.digitize ( 0.0, maxt, 0.001 );
      project < MyDigitalCurve > ( digitalCurve, "helix" );
      return true;
    }
    
    bool testDigitizeLine ()
    {
      typedef ParametricCurveDigitizer3D < MyLine >  Digitizer;
      typedef ParametricCurveDigitizer3D < MyLine >::DigitalCurve MyDigitalCurve;
      Digitizer digitizer;
      MyDigitalCurve digitalCurve;
      digitizer.attach ( line );
      digitizer.init ( digitalCurve );
      double maxt = 100;
      digitizer.digitize ( 0, maxt, 0.5 );
      project < MyDigitalCurve > ( digitalCurve, "line" );
      return true;
    }
    
    bool testDigitizeDisc ()
    {
      typedef ForwardRigidTransformation3D < Space, Identity, Space::RealPoint, Space::RealPoint > ForwardTrans;
      typedef BackwardRigidTransformation3D < Space, Identity, Space::RealPoint, Space::RealPoint > BackwardTrans;
      typedef DecoratorCurveTransformation < MyDisc, ForwardTrans, BackwardTrans > MyRotatedDisc;
      typedef ParametricCurveDigitizer3DDecorator < MyRotatedDisc >  Digitizer;
      typedef ParametricCurveDigitizer3D < MyRotatedDisc >::DigitalCurve MyDigitalCurve;
      

      ForwardTrans trans ( RealPoint ( 0, 0, 0 ), RealVector ( 0, 1, 0 ), M_PI_2, RealVector( 0, 0, 0) );
      BackwardTrans inverse ( RealPoint ( 0, 0, 0 ), RealVector ( 0, 1, 0 ), M_PI_2, RealVector( 0, 0, 0 ) );
      MyRotatedDisc rotDisc ( disc, trans, inverse );
      Digitizer digitizer;
      MyDigitalCurve digitalCurve;
      MyDigitalCurve digitalCurveErrors;
      digitizer.attach ( rotDisc );
//       digitizer.fastAdaptiveMethod ( true );
      digitizer.init ( digitalCurve );
      try
      {
	double maxt = 2 * M_PI - 0.04;
	digitizer.digitize ( 0.04, maxt, 0.00001 );
      } catch ( std::exception & e) {
// 	if (  std::abs ( digitalCurve.back()[0] - digitalCurve.front()[0] ) >= 2
// 	  || std::abs ( digitalCurve.back()[1] - digitalCurve.front()[1] ) >= 2
// 	  || std::abs ( digitalCurve.back()[2] - digitalCurve.front()[2] ) >= 2 )
// 	{
	;//return false;
// 	}
      }
      project < MyDigitalCurve > ( digitalCurve, "disc", &digitalCurveErrors );
      return true;
    }
};

///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

int main( int , char** )
{

// #ifdef __GNUC__
//    #ifndef NDEBUG
//     feenableexcept ( FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW );
//    #endif
// #endif

    bool res = true;
    testParametricCurves3D parametricCurveTest;
    trace.beginBlock ( "Testing Digitization of parametric 3D curves" );
    std::cout << std::asin(1) << std::endl;
//     res &= parametricCurveTest.testDigitizeEllipticHelix();
//     res &= parametricCurveTest.testDigitizeHelix();
//     res &= parametricCurveTest.testDigitizeLine();
    res &= parametricCurveTest.testDigitizeDisc();
    trace.emphase() << ( res ? "Passed." : "Error." ) << endl;
    trace.endBlock();
    return res ? 0 : 1;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

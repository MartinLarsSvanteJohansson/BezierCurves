/*************************/
/*
 * File: main.cpp
 * Author: Martin Johansson
 *
 * 
 *
 */
/**************************/


#include "CubicCurve.h"
#include "DirectionRad.h"

#include <iostream>


int main(int argc, char * argv[]){

	CCubicCurve cc(CMetricMapPoint(0,0), CMetricMapPoint(0,3), CMetricMapPoint(3,3), CMetricMapPoint(4,0));
	
	/* Elevate Degree of the curve */
	std::cout << "EXAMPLE #1. Elevation\nInput Curve: ";
	for (auto i : cc.getCurve())
		std::cout << "{" << i.getX() << "," << i.getY() << "}, ";

	cc.elevate();
	std::cout << "\nResulting elevated curve: ";
	for (auto i : cc.getCurve())
		std::cout << "{" << i.getX() << "," << i.getY() << "}, ";

	/* Split bezier curve into two */
	std::cout << "\n\nEXAMPLE #2. Splitting\n";
	CCubicCurve firstCurve;
	CCubicCurve secondCurve;

	// Split at t = 0.3
	cc.splitCurve(firstCurve, secondCurve, 0.3);


	/* Generate a cubic curve from a arc segment */
	std::cout << "\n\nEXAMPLE #3. Circular arc approximation\n";
	double radius = 1;
	double radStart = 0.0;
	double radSpan = CDirectionRad::getPi() / 2;
	CMetricMapPoint origin = CMetricMapPoint(0, 0);
	cc.arc2curve(radius, radStart, radSpan, origin);

    /* Generate a circle*/
    std::cout << "\nEXAMPLE #4. Circle approximation\n";
    cc.cubicCircle();
    std::cout << "Resulting control points:  ";
    for (auto i : cc.getCurve())
        std::cout << "{" << i.getX() << "," << i.getY() << "}, ";

	return 0;
}

/*
 * CubicCurve.cpp
 *
 *
 * Martin Johansson
 * martin.lars.svante.johansson@gmail.com
 *
 */



#include "CubicCurve.h"
#include <stdlib.h>
#include <cmath>

CCubicCurve::CCubicCurve(void)
{
	_length = -1;
}

CCubicCurve::CCubicCurve( const CMetricMapPoint start, const CMetricMapPoint ctrl1, const CMetricMapPoint ctrl2, const CMetricMapPoint end )
{
	_start = start;
	_ctrl1 = ctrl1;
	_ctrl2 = ctrl2;
	_end   = end;
	calculateLength();
	convert();
}

/**
* @brief Generates a cubic curve from two points (a line)
*/
CCubicCurve::CCubicCurve( CMetricMapPoint start, CMetricMapPoint end )
{
	_start = start;
	_end   = end;

	double dx = end.getX()-start.getX();
	double dy = end.getY()-start.getY();
	_length = sqrt(dx*dx + dy*dy);
	double theta = atan(dy/dx);
	
	// Assign the two control points to be on the line
	_ctrl1 = CMetricMapPoint(_start.getX() + _length*cos(theta) * 1 / 3, _start.getY() + _length*sin(theta) * 1 / 3);
	_ctrl2 = CMetricMapPoint(_start.getX() + _length*cos(theta) * 2 / 3, _start.getY() + _length*sin(theta) * 2 / 3);
	convert();

}

CCubicCurve::CCubicCurve( const CCubicCurve& cubicCurve )
{
	_start  = cubicCurve._start;
	_ctrl1  = cubicCurve._ctrl1;
	_ctrl2  = cubicCurve._ctrl2;
	_end    = cubicCurve._end;
	_length = cubicCurve._length;
}

CCubicCurve::~CCubicCurve(void)
{
}

CCubicCurve& CCubicCurve::operator=( const CCubicCurve& cubicCurve )
{
	_start  = cubicCurve._start;
	_ctrl1  = cubicCurve._ctrl1;
	_ctrl2  = cubicCurve._ctrl2;
	_end    = cubicCurve._end;
	_length = cubicCurve._length;
	return *this;
}

/**
* @brief Returns the \p X,Y point at the l position on the cubic curve
*/
CMetricMapPoint CCubicCurve::getPointOnCurve( double l )
{
	if (l < 0)
		l = 0;
	if (l > 1)
		l= 1;

	return CMetricMapPoint(calculatePos( l, _start.getX(), _ctrl1.getX(), _ctrl2.getX(), _end.getX() ), calculatePos( l, _start.getY(), _ctrl1.getY(), _ctrl2.getY(), _end.getY() ));
}


/**
* @returns The cubic curve
*/
std::vector<CMetricMapPoint> CCubicCurve::getCurve() const
{
	return _curve;
}

double CCubicCurve::getCurvatureOnCurve(double l)
{
	return 0;
}

/**
*  \todo Implement
*/
CDirectionRad CCubicCurve::getDirectionOnCurve(double l)
{
	return CDirectionRad(0);
}
 
/**
* @brief What is the heading at a certain point on the segment
*/
double CCubicCurve::getHeading(double l)
{
	CMetricMapPoint p = calculateDerivative(l);
	return std::atan2(p.getY(), p.getX());
}

/**
* @brief Retrieves the lenght of the curve.
*/
double CCubicCurve::getLenght()
{
	if (_length == -1)
		calculateLength();

	return _length;
}


CMetricMapPoint CCubicCurve::getMin()
{
	CMetricMapPoint min = _start;

	if (min.getX() > _ctrl1.getX()) min.setX(_ctrl1.getX());
	if (min.getY() > _ctrl1.getY()) min.setY(_ctrl1.getY());

	if (min.getX() > _ctrl2.getX()) min.setX(_ctrl2.getX());
	if (min.getY() > _ctrl2.getY()) min.setY(_ctrl2.getY());

	if (min.getX() > _end.getX()) min.setX(_end.getX());
	if (min.getY() > _end.getY()) min.setY(_end.getY());

	return min;
}

/**
* @brief Return maximum value in X and Y coordinates of the four points.
*/
CMetricMapPoint CCubicCurve::getMax()
{
	CMetricMapPoint max = _start;

	if (max.getX() < _ctrl1.getX()) max.setX(_ctrl1.getX());
	if (max.getY() < _ctrl1.getY()) max.setY(_ctrl1.getY());

	if (max.getX() < _ctrl2.getX()) max.setX(_ctrl2.getX());
	if (max.getY() < _ctrl2.getY()) max.setY(_ctrl2.getY());

	if (max.getX() < _end.getX()) max.setX(_end.getX());
	if (max.getY() < _end.getY()) max.setY(_end.getY());

	return max;
}

/**
* @brief The Bernstein polynomial equation
*/
double CCubicCurve::calculatePos(double l, double P0, double P1, double P2, double P3)
{
	return (1 - l)*(1 - l)*(1 - l)*P0 + 3 * (1 - l)*(1 - l)*l*P1 + 3 * (1 - l)*l*l*P2 + l*l*l*P3;
}

CMetricMapPoint CCubicCurve::calculateDerivative(double l)
{
	// Q(t) = 3(1-t)2(B-A) + 3t2(D-C) + 6t(1-t)(C-B), t = [0-1]
	CMetricMapPoint q = (_ctrl1 - _start) * 3 * (1 - l*l) + (getCtrl2() - getCtrl1()) * 6 * l*(1 - l) + (getEnd() - getCtrl2()) * 3 * l*l*l;
	return q;
}

double CCubicCurve::calculateXPos(double l)
{
	return calculatePos(l, _start.getX(), _ctrl1.getX(), _ctrl2.getX(), _end.getX());
}

double CCubicCurve::calculateYPos(double l)
{
	return calculatePos(l, _start.getY(), _ctrl1.getY(), _ctrl2.getY(), _end.getY());
}

CCubicCurve CCubicCurve::createRandomCubicCurve( CMetricMapPoint startPoint, double Lenght )
{
	CMetricMapPoint ctrl1Point	= createRandomCubicCurvePoint(startPoint, Lenght/4);
	CMetricMapPoint ctrl2Point	= createRandomCubicCurvePoint(ctrl1Point, Lenght/4);
	CMetricMapPoint endPoint	= createRandomCubicCurvePoint(ctrl2Point, Lenght/4);
	return CCubicCurve(startPoint, ctrl1Point, ctrl2Point, endPoint);
}

CCubicCurve CCubicCurve::createRandomCubicCurve( CMetricMapPoint startPoint, CDirectionRad direction, double lenght )
{
	CMetricMapPoint ctrl1Point = createRandomCubicCurvePoint(startPoint, direction, lenght/4);
	CMetricMapPoint ctrl2Point = createRandomCubicCurvePoint(ctrl1Point, direction, lenght/4);
	CMetricMapPoint endPoint = createRandomCubicCurvePoint(ctrl2Point,   direction, lenght/4);
	return CCubicCurve(startPoint, ctrl1Point, ctrl2Point, endPoint);
}

CMetricMapPoint CCubicCurve::createRandomCubicCurvePoint( CMetricMapPoint startPoint, double Lenght )
{
	CMetricMapPoint tmpStartPoint(startPoint.getX(),startPoint.getY());

	// Rand direction between 0 and 360 degrees
	CDirectionRad direction = (double)(((double)rand()/RAND_MAX) * 6.28);
	CMetricMapPoint tmpEndPoint(tmpStartPoint);

	// Set the direction to 90 degrees
	direction = 3.14/2;
	tmpEndPoint.set(tmpStartPoint,direction,Lenght);

	return CMetricMapPoint(tmpEndPoint.getX(), tmpEndPoint.getY());
}

CMetricMapPoint CCubicCurve::createRandomCubicCurvePoint( CMetricMapPoint startPoint, CDirectionRad startDirection, double lenght )
{
	CMetricMapPoint tmpStartPoint(startPoint.getX(),startPoint.getY());
	CDirectionRad direction = startDirection + (double)(((double)rand()/RAND_MAX) * 6.28/10);
	direction=3.14/2;
	CMetricMapPoint tmpEndPoint(tmpStartPoint);
	tmpEndPoint.set(tmpStartPoint,direction,lenght);

	return CMetricMapPoint(tmpEndPoint.getX(), tmpEndPoint.getY());
}

void CCubicCurve::convert()
{
	_curve.clear();
	_curve.push_back(_start);
	_curve.push_back(_ctrl1);
	_curve.push_back(_ctrl2);
	_curve.push_back(_end);
}

/**
* @brief Create one new cubic curve of two and return maximum error
*/
double CCubicCurve::join( CCubicCurve cubicCurveSourceFirst, CCubicCurve cubicCurveSourceLast, CCubicCurve *cubicCurveTarget )
{
	cubicCurveTarget->_start = cubicCurveSourceFirst._start;
	cubicCurveTarget->_ctrl1 = cubicCurveSourceFirst._ctrl1;
	cubicCurveTarget->_ctrl2 = cubicCurveSourceLast._ctrl2;
	cubicCurveTarget->_end   = cubicCurveSourceLast._end;

	/* Move Ctrl1 and Ctrl2 along their trajectories from start and end to adjust to the source curves to minimize the error */

	/* Calculate the maximum error */
	double    error = 0;
	double maxError = 0;
	CMetricMapPoint tmpPoint;
	CMetricMapPoint targetPoint;
	for (int i = 0; i < 2; i++)	{
		for (double l = 0; l<1; l+=0.1)	{
			if (i == 0)
			{
				tmpPoint    = cubicCurveSourceFirst.getPointOnCurve(l);
				targetPoint = cubicCurveTarget->getPointOnCurve(l/2);
				error       = targetPoint.distance(tmpPoint);
			}
			if (i == 1)
			{
				tmpPoint    = cubicCurveSourceLast.getPointOnCurve(l);
				targetPoint = cubicCurveTarget->getPointOnCurve(0.5 + l/2);
				error       = abs(targetPoint.distance(tmpPoint));
			}
			if(error > maxError) maxError = error;
		}
	}

	return maxError;
}

double CCubicCurve::minimalDistanceBetween(std::vector<CCubicCurve> cubicCurves, CMetricMapPoint point)
{
	if(cubicCurves.size() == 0) return -1;
	int selectedCuevIndex = -1;
	if (cubicCurves.size() == 1) selectedCuevIndex = 0;
	else
	{
		/* Check witch curve segment that is the nearest one! */
		double distanceToStartEnd = abs(point.distance(cubicCurves[0]._start) + point.distance(cubicCurves[0]._end)); 
		selectedCuevIndex = 0;
		for (unsigned int curvIndex = 0; curvIndex < cubicCurves.size(); curvIndex++)
		{
			if (abs(point.distance(cubicCurves[curvIndex]._start) + point.distance(cubicCurves[curvIndex]._end)) < distanceToStartEnd)
			{
				distanceToStartEnd = abs(point.distance(cubicCurves[curvIndex]._start) + point.distance(cubicCurves[curvIndex]._end)); 
				selectedCuevIndex = curvIndex;
			}
		}
	}

	if (selectedCuevIndex != -1)
	{
		return cubicCurves[selectedCuevIndex].minimalDistanceTo(point);
	}
	return -1;
}


/**
* @brief  Returns the minimum distance between the cubic curve and the input point
*/
double CCubicCurve::minimalDistanceTo(CMetricMapPoint point)
{
	if (point == _start) return 0;
	if (point == _end) return   0;

	//TODO: Make defines for CurvIndexStep and ...
	double curvIndex     = 0.25;
	double curvIndexStep = 0.1;
	double distance = point.distance(_start);
	while (abs(curvIndexStep) > 0.001)
	{
		while(point.distance(getPointOnCurve(curvIndex)) < distance)
		{
			distance = point.distance(getPointOnCurve(curvIndex));
			curvIndex += curvIndexStep;
		}
		curvIndexStep = curvIndexStep * -0.1;
	}
	return distance;
}


/**
* @brief  Same as \p minimalDistanceTo(CMetricMapPoint point) but also returns at what parametric value the minimum distance is found.
* @see minimalDistanceTo(CMetricMapPoint point)
*/
double CCubicCurve::minimalDistanceTo(CMetricMapPoint point, double &tMin)
{
	if (point == _start) return 0;
	if (point == _end) return   0;

	double stepSize = 0.01;
	double minDist = 9999999;

	for(double curvIndex = 0; curvIndex < (1+stepSize); curvIndex = curvIndex + stepSize)
	{
		CMetricMapPoint p = getPointOnCurve(curvIndex);
		double dist = point.distance(p);
		
		if(dist < minDist) 
		{
			minDist = dist;
			tMin = curvIndex;
		}
	}

	return minDist;
}


/**
* @brief The direction from start to ctrl1 and end to ctrl2 is set. Optimize the distance to minimize curve distance to point
*/
double CCubicCurve::optimizeCtrlPointsDistance(CMetricMapPoint point)
{
	CMetricMapPoint bestCtrl1      = _ctrl1;
	CMetricMapPoint bestCtrl2      = _ctrl2;
	CDirectionRad startDirection   = _start.direction(_ctrl1);
	CDirectionRad endDirection     = _end.direction(_ctrl2);
	double        startDistance    = _start.distance(_ctrl1);
	double        endDistance      = _end.distance(_ctrl2);
	double        pointDistance    = minimalDistanceTo(point);
	double        pointDistanceOld = pointDistance*2;

	double  stepSizeStartDistance         = _start.distance(_end)/12;
	double  stepSizeEndDistance           = _start.distance(_end)/12;
	double  stepSizeStartDistanceCurrent  = stepSizeStartDistance;
	double  stepSizeEndDistanceCurrent    = stepSizeEndDistance;

	/*
	Step direction
	---------------
	0: inc distance _start -> ctrl1
	1: inc distance _start -> ctrl1 & inc distance _end   -> ctrl2
	2: inc distance _start -> ctrl1 & dec distance _end   -> ctrl2
	3: inc distance _end   -> ctrl2
	4: dec distance _end   -> ctrl2
	5: dec distance _start -> ctrl1
	6: dec distance _start -> ctrl1 & inc distance _end   -> ctrl2
	7: dec distance _start -> ctrl1 & dec distance _end   -> ctrl2
	*/
	int     stepDirection         = 0; 

	/* Optimize CTRL1 & CTRL2 */
	while((abs(stepSizeStartDistance) > 0.01) || (abs(stepSizeEndDistance) > 0.01)) 
	{
		/* Check if need to change step direction */
		if (pointDistance < pointDistanceOld) // the distance decreases!
		{
			pointDistanceOld = pointDistance;
			bestCtrl1      = _ctrl1;
			bestCtrl2      = _ctrl2;
		}
		else // Change direction
		{
			stepDirection++;
			if (stepDirection == 8) stepDirection = 0;
			switch (stepDirection)
			{
				case 0:
					stepSizeStartDistance = stepSizeStartDistance * 0.1;
					stepSizeEndDistance   = stepSizeEndDistance * 0.1;
					
					stepSizeStartDistanceCurrent = stepSizeStartDistance;
					stepSizeEndDistanceCurrent   = 0;
					break;
				case 1:
					stepSizeStartDistanceCurrent = stepSizeStartDistance;
					stepSizeEndDistanceCurrent   = stepSizeEndDistance;
					break;
				case 2:
					stepSizeStartDistanceCurrent = stepSizeStartDistance;
					stepSizeEndDistanceCurrent   = -stepSizeEndDistance;
					break;
				case 3:
					stepSizeStartDistanceCurrent = 0;
					stepSizeEndDistanceCurrent   = stepSizeEndDistance;
					break;
				case 4:
					stepSizeStartDistanceCurrent = 0;
					stepSizeEndDistanceCurrent   = -stepSizeEndDistance;
					break;
				case 5:
					stepSizeStartDistanceCurrent = -stepSizeStartDistance;;
					stepSizeEndDistanceCurrent   = 0;
					break;
				case 6:
					stepSizeStartDistanceCurrent = -stepSizeStartDistance;;
					stepSizeEndDistanceCurrent   = stepSizeEndDistance;
					break;
				case 7:
					stepSizeStartDistanceCurrent = -stepSizeStartDistance;;
					stepSizeEndDistanceCurrent   = -stepSizeEndDistance;
					break;
				default:
					break;
			}
		}

		/* Test new ctrl points */
		startDistance += stepSizeStartDistanceCurrent;
		endDistance   += stepSizeEndDistanceCurrent;
		
		_ctrl1.set(_start, startDirection, startDistance);
		_ctrl2.set(_end,   endDirection,   endDistance);
		pointDistance = minimalDistanceTo(point);
	}
	return minimalDistanceTo(point);
}

/**
* @brief Create one cubic (this) curve and return maximum error
*/
double CCubicCurve::cubuicCurveFromMetricPoints(CMetricMapPoint before, std::vector<CMetricMapPoint> points, CMetricMapPoint after)
{
	if (points.size()<2) return -1; 
	_start = points.front();
	_end = points.back();
	CDirectionRad startDirection = before.direction(_start);
	CDirectionRad endDirection = after.direction(_end);
	double lenght = _start.distance(_end);
	_ctrl1.set(_start, startDirection, lenght/3);
	_ctrl2.set(_end, endDirection, lenght/3);

	double minimalDistance    = 0;
	double maximalDistance    = 0;
	double minimalDistanceOld = 0;
	double maximalDistanceOld = 0;
	double deltaError         = lenght;
	double stepSize           = lenght/12;
	double ctrlPoint1Lenght   = lenght/3;
	double ctrlPoint2Lenght   = lenght/3;

	/* Optimera ctrl1 */
	while(abs(stepSize) > 0.01) 
	{
		_ctrl1.set(_start, startDirection, ctrlPoint1Lenght);
		for (unsigned int pointIndex = 1; pointIndex < points.size()-1; pointIndex++)
		{
			minimalDistance = minimalDistanceTo(points[pointIndex]);
			if (minimalDistance > maximalDistance) maximalDistance = minimalDistance;
		}
		deltaError = maximalDistanceOld - maximalDistance;
		if (deltaError > 0) stepSize = stepSize * -0.01;
		else                ctrlPoint1Lenght = ctrlPoint1Lenght + stepSize;
	}

	/* Optimera ctrl2 */
	minimalDistance    = 0;
	maximalDistance    = 0;
	minimalDistanceOld = 0;
	maximalDistanceOld = 0;
	deltaError         = lenght;
	stepSize = lenght/12;
	while(abs(stepSize) > 0.01) 
	{
		_ctrl2.set(_end,   endDirection,   ctrlPoint2Lenght);
		for (unsigned int pointIndex = 1; pointIndex < points.size()-1; pointIndex++)
		{
			minimalDistance = minimalDistanceTo(points[pointIndex]);
			if (minimalDistance > maximalDistance) maximalDistance = minimalDistance;
		}
		deltaError = maximalDistanceOld - maximalDistance;
		if (deltaError > 0) stepSize = stepSize * -0.01;
		else                ctrlPoint2Lenght = ctrlPoint2Lenght + stepSize;
	}

	/* check error */
	minimalDistance    = 0;
	maximalDistance    = 0;
	for (unsigned int pointIndex = 1; pointIndex < points.size()-1; pointIndex++)
	{
		minimalDistance = minimalDistanceTo(points[pointIndex]);
		if (minimalDistance > maximalDistance) maximalDistance = minimalDistance;
	}

	return maximalDistance;
}

void CCubicCurve::calculateLength()
{
	CMetricMapPoint p1,p2;
	double length = 0;
	p1 = _start;
	double stepSize = 0.01;

	for(double l = stepSize; l < (1 - stepSize); l+=stepSize)
	{
		p2 = getPointOnCurve(l);
		length += p1.distance(p2);
		p1 = p2;
	}
	length += p1.distance(_end);
	_length = length;
}

double CCubicCurve::calculateLength( std::vector<CCubicCurve> track )
{
	double lenght = 0;
	for (int i = 0; i < track.size(); i++) 
		lenght += track[i].getLenght();
	return lenght;
}

/**
* @brief Given a point known to be on the curve, get the metric distance path to it
*/
double CCubicCurve::calculateLengthToPointOnCurve(CMetricMapPoint P)
{

	double t;
	minimalDistanceTo(P, t);

	/* Given a l value, find the distance */
	CMetricMapPoint p1,p2;
	p1 = _start;
	double length = 0;
	double stepSize = 0.01;

	double l;
	for(l = stepSize; l < (t - stepSize); l+=stepSize)
	{
		p2 = getPointOnCurve(l);
		length += p1.distance(p2);
		p1 = p2;
	}
	p2 = getPointOnCurve(l+stepSize);
	length += p1.distance(p2);
	return length;

}

void CCubicCurve::removeShortSegments(std::vector<CCubicCurve> &track, double maxLenght )
{
	double segmentLenght = (track)[0].getLenght();
	double longestSegmentLenght = (track)[0].getLenght();
	int shortestSegmentindex = 0;
	std::vector<CCubicCurve> tmpTrack;

	for (int segmentIndex = 0; segmentIndex < track.size(); segmentIndex++)
	{
		if ((track)[segmentIndex].getLenght() > longestSegmentLenght) longestSegmentLenght = (track)[segmentIndex].getLenght();
		if ((track)[segmentIndex].getLenght() < segmentLenght) 
		{
			segmentLenght = (track)[segmentIndex].getLenght();
			shortestSegmentindex = segmentIndex;
		}
	}

	if (longestSegmentLenght < maxLenght) {
		//throw std::logic_error(std::string(__FILE__) + _NL_ + std::string(__FUNCTION__) + _NL_ + std::string("No segment longer then max lenght!"));
	}

	/* clear short segments */
	while(segmentLenght < maxLenght)
	{
		segmentLenght = (track)[0].getLenght();
		shortestSegmentindex = 0;
		for (int segmentIndex = 0; segmentIndex < (track).size(); segmentIndex++)
		{
			if ((track)[segmentIndex].getLenght() < segmentLenght) 
			{
				segmentLenght = (track)[segmentIndex].getLenght();
				shortestSegmentindex = segmentIndex;
			}
		}

		for (int segmentIndex = 0; segmentIndex < shortestSegmentindex; segmentIndex++) tmpTrack.push_back((track)[segmentIndex]);
		for (int segmentIndex = shortestSegmentindex + 1; segmentIndex < (track).size(); segmentIndex++) tmpTrack.push_back((track)[segmentIndex]);
		// If the shortest segment is the first segment!
		if (shortestSegmentindex == 0)                      tmpTrack[0]._start                      = (track)[0]._start;
		// If the shortest segment is the last segment
		else if (shortestSegmentindex == (track).size()-1) tmpTrack[shortestSegmentindex-1]._end   = (track)[shortestSegmentindex]._end;
		// All other cases
		else	                                            tmpTrack[shortestSegmentindex - 1]._end = tmpTrack[shortestSegmentindex]._start;
		(track) = tmpTrack;
		tmpTrack.clear();
	}

	segmentLenght = (track)[0].getLenght();
	shortestSegmentindex = 0;
	for (int segmentIndex = 0; segmentIndex < (track).size(); segmentIndex++)
	{
		if ((track)[segmentIndex].getLenght() < segmentLenght) 
		{
			segmentLenght = (track)[segmentIndex].getLenght();
			shortestSegmentindex = segmentIndex;
		}
	}
}

double CCubicCurve::findClosesPointOnPath( std::vector<CCubicCurve> track, CMetricMapPoint point, unsigned int &index, double &parameter )
{
	
	double shortestDistanceToPoint, testDistanceToPoint;
	//*index     = 0;
	if(index > 1) index = index - 1;
	parameter = 0;
	shortestDistanceToPoint = testDistanceToPoint = point.distance(track[index].getPointOnCurve(0));
	if (shortestDistanceToPoint == 0) return shortestDistanceToPoint;
	for (unsigned int cubicCurvIndex = index; cubicCurvIndex < track.size(); cubicCurvIndex++)
	{
		for (double posOnCubicCurv = 0; posOnCubicCurv <= 1; posOnCubicCurv += 0.001)
		{
			testDistanceToPoint = point.distance(track[cubicCurvIndex].getPointOnCurve(posOnCubicCurv));
			if(testDistanceToPoint < shortestDistanceToPoint)
			{
				shortestDistanceToPoint = testDistanceToPoint;
				index     = cubicCurvIndex;
				parameter = posOnCubicCurv;
			}
		}
	}
	return shortestDistanceToPoint;
}


void CCubicCurve::addCurvePoint(double x, double y)
{
	_curve.push_back(CMetricMapPoint(x, y));
}

/**
* @brief Creates a smooth bezier spline through prescribed knot points
* \todo overload the + operator to make this more efficient
* \bug Thomas algorithm not solving correctly
*/
void CCubicCurve::createSplineFromPoints(std::vector<CMetricMapPoint> &points)
{
	double n = points.size()-1;

	std::vector<int> a;
	std::vector<int> b;
	std::vector<int> c;
	std::vector<int> r;
	std::vector<int> rY;
	std::vector<double> p1;
	std::vector<double> p2;
	std::vector<double> p1Y;
	std::vector<double> p2Y;
	std::vector<CMetricMapPoint> _p1;

	for(int i=0; i <n ; i++) {
		p1.push_back(-1);
		p2.push_back(-1);
		p1Y.push_back(-1);
		p2Y.push_back(-1);
	}
	
	// Left most segment
	a.push_back(0);
	b.push_back(2);
	c.push_back(1);
	r.push_back(static_cast<int>(points[0].getX()) + static_cast<int>(2 * points[1].getX()));
	rY.push_back( points[0].getY() + 2*points[1].getY() );

	std::vector<CMetricMapPoint>::iterator it = points.begin();
	// Skip the first point
	it++;
	
	// Internal segments
	for(int i = 1; i < (n-1) ; i++)
	{
		a.push_back(1);
		b.push_back(4);
		c.push_back(1);
		r.push_back( 4*points[i].getX() + 2*points[i+1].getX() );
		rY.push_back( 4*points[i].getY() + 2*points[i+1].getY() );
	}

	// Right segment
	a.push_back(2);
	b.push_back(7);
	c.push_back(0);
	r.push_back( 8*points[n-1].getX() + points[n].getX() );
	rY.push_back( 8*points[n-1].getY() + points[n].getY() );

	// Solve Ax=b with the Thomas algorithm 
	for(int i = 1; i < n; i++)
	{
		double m = a[i]/b[i-1];
		b[i] = b[i] - m*c[i-1];
		r[i] = r[i] - m*r[i-1];
		rY[i] = rY[i] - m*rY[i-1];
	}

	// This gives us p1
	p1[n-1] = r[n-1] / b[n-1];
	p1Y[n-1] = rY[n-1] / b[n-1];

	for(int i = n-2 ; i>= 0; --i)
	{
		p1[i] = (r[i] - c[i] * p1[i+1] ) / b[i];
		p1Y[i] = (rY[i] - c[i] * p1Y[i+1] )/ b[i];
	}

	// We now know p1 and can now compute p2
	for(int i=0; i < n-1; i++)
	{
		p2[i] = 2*points[i+1].getX() - p1[i+1];
		p2Y[i] = 2*points[i+1].getY() - p1Y[i+1];
	}

	p2[n-1] = 0.5*(points[n].getX() + p1[n-1]);
	p2Y[n-1] = 0.5*(points[n].getY() + p1Y[n-1]);

	for(int i = 0; i < n; i++)
	{
		_curve.push_back(CMetricMapPoint( points[i].getX(), points[i].getY() ));		
		_curve.push_back(CMetricMapPoint(p1[i], p1Y[i]));
		_curve.push_back(CMetricMapPoint(p2[i], p2Y[i]));
		_curve.push_back(CMetricMapPoint( points[i+1].getX(), points[i+1].getY() ));

	}
	
}

/**
 * @brief This function splits the curve into two parts at a given \p t value, 
 * @param[out] firstPart First cubic curve.
 * @param[out] lastPart Second cubic curve.
 * @param[in] t At what parametric distance to split 
 */
void CCubicCurve::splitCurve(CCubicCurve &firstPart, CCubicCurve &lastPart, double t)
{
	CCubicCurve orig; // original curve
	orig = *this;
	if(this->_curve.size() == 0)
	{
		this->_curve.push_back(this->_start);
		this->_curve.push_back(this->_ctrl1);
		this->_curve.push_back(this->_ctrl2);
		this->_curve.push_back(this->_end);
	}
	this->splitCurve(t);

	firstPart._start = this->_curve[0];
	firstPart._ctrl1 = this->_curve[1];
	firstPart._ctrl2 = this->_curve[2];
	firstPart._end   = this->_curve[3];

	lastPart._start = this->_curve[4];
	lastPart._ctrl1 = this->_curve[5];
	lastPart._ctrl2 = this->_curve[6];
	lastPart._end   = this->_curve[7];

	*this = orig;
}

void CCubicCurve::splitCurveKeepFirst( double t)
{
	if(this->_curve.size() == 0)
	{
		this->_curve.push_back(this->_start);
		this->_curve.push_back(this->_ctrl1);
		this->_curve.push_back(this->_ctrl2);
		this->_curve.push_back(this->_end);
	}
	this->splitCurve(t);

	this->_start = this->_curve[0];
	this->_ctrl1 = this->_curve[1];
	this->_ctrl2 = this->_curve[2];
	this->_end   = this->_curve[3];
}

void CCubicCurve::splitCurveKeepLast( double t)
{
	if(this->_curve.size() == 0)
	{
		this->_curve.push_back(this->_start);
		this->_curve.push_back(this->_ctrl1);
		this->_curve.push_back(this->_ctrl2);
		this->_curve.push_back(this->_end);
	}
	this->splitCurve(t);

	this->_start = this->_curve[4];
	this->_ctrl1 = this->_curve[5];
	this->_ctrl2 = this->_curve[6];
	this->_end   = this->_curve[7];
}

/**
* @brief This function splits the curve into two parts at a given t value and stores it in _curve
*/
void CCubicCurve::splitCurve(double t)
{
	std::vector<std::vector<double> > M;
	
	CMetricMapPoint P0, P1, P2, P3;
	CMetricMapPoint P4, P5, P6, P7;
	
	// Calculates the first cubic curve segment
	P0 = _curve[0];
	P1 = _curve[1]*t - _curve[0]*(t-1);
	P2 = _curve[2]*t*t - _curve[1]*2*t*(t-1) + _curve[0]*(t-1)*(t-1);
	P3 = _curve[3]*t*t*t - _curve[2]*3*t*t*(t-1) + _curve[1]*3*t*(t-1)*(t-1) - _curve[0]*(t-1)*(t-1)*(t-1);
	
	// Calculates the second cubic curve segment
	P4 = _curve[3]*t*t*t - _curve[2]*3*t*t*(t-1) + _curve[1]*3*t*(t-1)*(t-1) - _curve[0]*(t-1)*(t-1)*(t-1);
	P5 = _curve[3]*t*t - _curve[2]*2*t*(t-1) + _curve[1]*(t-1)*(t-1);
	P6 = _curve[3]*t - _curve[2]*(t-1);
	P7 = _curve[3];

	// Replace the four first points, and then add the last four points
	_curve[0] = P0;
	_curve[1] = P1;
	_curve[2] = P2;
	_curve[3] = P3;
	_curve.push_back(P4);
	_curve.push_back(P5);
	_curve.push_back(P6);
	_curve.push_back(P7);
}

/**
* @brief Finds the normalized tangent at a given parametric position t. Calculated from the first derivative of the cubic bezier
*/
CMetricMapPoint CCubicCurve::findTangent(double t)
{
	//\mathbf{B}'(t) = 3(1-t)^2(\mathbf{P}_1 - \mathbf{P}_0) + 6(1-t)t(\mathbf{P}_2 - \mathbf{P}_1) + 3t^2(\mathbf{P}_3 - \mathbf{P}_2) \,.
	CMetricMapPoint m = ( (_curve[1] - _curve[0])*3*(1-t)*(1-t) + (_curve[2] - _curve[1])*6*t*(1-t) + (_curve[3] - _curve[2])*3*t*t );
	double d = sqrt(m.getX()*m.getX() + m.getY()*m.getY());
	m.setX(m.getX()/d);
	m.setY(m.getY()/d);
	return m;
}

/**
* @brief Elevate the order of the Bezier curve by one
* @returns Updates \p _curve with the new control points
*/
void CCubicCurve::elevate()
{
	if (_length == -1)
		return;

	int L = _curve.size();
	if (L == 0) {
		convert();
		L = _curve.size();
	}

	std::vector<CMetricMapPoint> newCurve;
	newCurve.resize(L + 1);
	newCurve[0] = (_curve[0]);
	
	for (double i = 1; i < L; i++)
	{
		double x = (i / L) * _curve[i - 1].getX() + (L - i) / L *_curve[i].getX();
		double tmp1 = (i / L) *_curve[i - 1].getY();
		double y = (i / L) * _curve[i - 1].getY() + ( (L - i) / L )*_curve[i].getY();
		newCurve[i] =  CMetricMapPoint(x, y) ;
	}
	newCurve[L] = _curve[L-1];
	_curve = newCurve;


}

/**
* @brief  Move the cubic curve x and y steps. Used mainly to debug stuff
*/
void CCubicCurve::moveCurve(double x, double y)
{
	_start.set(	_start.getX()+ x, _start.getY() + y);
	_ctrl1.set(	_ctrl1.getX()+ x, _ctrl1.getY() + y);
	_ctrl2.set(	_ctrl2.getX()+ x, _ctrl2.getY() + y);
	_end.set( _end.getX() + x, _end.getY() + y);

	_curve.clear();
	_curve.push_back(_start);
	_curve.push_back(_ctrl1);
	_curve.push_back(_ctrl2);
	_curve.push_back(_end);

}

/**
* @brief Debugging function currently only for straight line curves
*/
bool CCubicCurve::isPointOnCurve(CMetricMapPoint P)
{
	
	double eps = 5;

	float a = (_end.getY()- _start.getY()) / (_end.getX() - _start.getX());
    float b = _start.getY() - a * _start.getX();
   if ( fabs(P.getY() - (a*P.getX()+b)) < eps)
	   return true;   
   else
	return false;

}

/**
* @brief Generates an approximation to a full circle
*/
void CCubicCurve::cubicCircle()
{
    double const c = 0.551915024494; // This gives a maximal radial drift of 0.019608%
    _curve.clear();
    _start = CMetricMapPoint(0, 0);
    _ctrl1 = CMetricMapPoint(c, 1);
    _ctrl2 = CMetricMapPoint(1, c);
    _end = CMetricMapPoint(1, 0);

    _curve.push_back(_start);
    _curve.push_back(_ctrl1);
    _curve.push_back(_ctrl2);
    _curve.push_back(_end);

    _curve.push_back(CMetricMapPoint(1, 0));
    _curve.push_back(CMetricMapPoint(1, -c));
    _curve.push_back(CMetricMapPoint(c, -1));
    _curve.push_back(CMetricMapPoint(-1, 0));

    _curve.push_back(CMetricMapPoint(0, -1));
    _curve.push_back(CMetricMapPoint(-c, -1));
    _curve.push_back(CMetricMapPoint(-1, -c));
    _curve.push_back(CMetricMapPoint(-1, 0));

    _curve.push_back(CMetricMapPoint(-1, 0));
    _curve.push_back(CMetricMapPoint(-1, c));
    _curve.push_back(CMetricMapPoint(-c, 1));
    _curve.push_back(CMetricMapPoint(0, 1));

}

void CCubicCurve::arc2curve(double R, double radStart, double radTotal, CMetricMapPoint origin)
{

	// If the segment is short, skip it
	if (radTotal < 0.01)
		return;

	// Rotational angel
	double rotAngle = radTotal / 2 + radStart;

	// Find the circle arc control points
	double t1 = radTotal / 2;
	double x0 = R*std::cos(t1);
	double y0 = -R * std::sin(t1);
	double x3 = R*std::cos(t1);
	double y3 = R*std::sin(t1);
	double x1 = (4 * R - R*cos(t1)) / 3;
	double y1 = ((R - R*cos(t1))*(R*cos(t1) - 3 * R)) / (3 * R*sin(t1));
	double x2 = x1;
	double y2 = -y1;

	_curve.clear();
	_curve.push_back(CMetricMapPoint(x0, y0));
	_curve.push_back(CMetricMapPoint(x1, y1));
	_curve.push_back(CMetricMapPoint(x2, y2));
	_curve.push_back(CMetricMapPoint(x3, y3));

	// Rotate and translate all the points
	for (int i = 0; i<4; i++) {

		_curve[i].setX(_curve[i].getX() + origin.getX());
		_curve[i].setY(_curve[i].getY() + origin.getY());

		// Translate back to origin
		_curve[i].setX(_curve[i].getX() - origin.getX());
		_curve[i].setY(_curve[i].getY() - origin.getY());

		// Rotate points
		double px_temp = _curve[i].getX()*cos(rotAngle) - _curve[i].getY()*sin(rotAngle);
		double py_temp = _curve[i].getX()*sin(rotAngle) + _curve[i].getY()*cos(rotAngle);

		// Translate back
		_curve[i].setX(px_temp + origin.getX());
		_curve[i].setY(py_temp + origin.getY());
	}

	_start = _curve[0];
	_ctrl1 = _curve[1];
	_ctrl2 = _curve[2];
	_end = _curve[3];
}


/* -- Currently being merged from Matlab/Python below
  
// Least square cubic curve fitting below -- 
void CCubicCurve::startFit(std::vector<CMetricMapPoint>& p)
{
	// Max allowed square distance
	double maxAllowedSqD = 0.001;
	// First and last point (-- should be 0 to end-1?)
	int ibi[] = {1, p.size()}	;

	// The resulting fitted cubic curves
	std::vector<CCubicCurve> K;
	// Final break point indices (-- optional?)
	std::vector<int> fbi; 
	std::vector<CCubicCurve> M = _bzapproxu(p, maxAllowedSqD, ibi), fbi; 


}


std::vector<CCubicCurve> CCubicCurve::_bzapproxu(std::vector<CMetricMapPoint>& p, double maxAllowedSqD, int ibi[], std::vector<int>& fbi )
{
	if(p.size() < 4) {
		// Print error - needs at least four points required
	}

	if(maxAllowedSqD < 0) {
		// Print error - Must be greater than zero
	}

	/// \todo - make sure first and last are included, and sort
	// ibi=[ibi; 1; size(Mat,1)];

	std::vector<CCubicCurve> K;
	std::vector<int> ti;
	bool uniform = true;

	FindBzCP4AllSeg(p, ibi, uniform, K, ti);
}
*/
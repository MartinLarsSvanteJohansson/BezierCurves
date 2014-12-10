#pragma once
#include "MetricMapPoint.h"
#include <vector>


/*! \class CCubicCurve 
 *  \brief     CCubicCurve is a 2-dimensional path implemented with 4 CMetricMapPoints
 *  \details   ...
 *
 *  \author    Martin Johansson
 *  \version   0.1
 *  \date      2014-05-04
 *
 * \par Status
 * - Declared
 * - \b Started
 * - Ready
 * - Tested
 * - Approved
 * - Closed
 * - Locked
 *  \bug       --
 *  \warning   No Error handling implemented
 *  \todo Upgrade the documenetation
 */
class CCubicCurve
{
public:
	/* Default Constructor */
	CCubicCurve(void);
	CCubicCurve(const CMetricMapPoint start, const CMetricMapPoint ctrl1,	const CMetricMapPoint ctrl2, const CMetricMapPoint end);
	CCubicCurve( CMetricMapPoint start, CMetricMapPoint end );
	/* Copy Constructor */
	CCubicCurve(const CCubicCurve& cubicCurve);
	/* Default Destructor */
	~CCubicCurve(void);

	/* Assignment operator, and makes them the same size */
	CCubicCurve& operator=(const CCubicCurve& cubicCurve); 
	
	

	/* Access curve data */
	CMetricMapPoint                 getPointOnCurve(double l);
	CDirectionRad                   getDirectionOnCurve(double l);
	double                          getCurvatureOnCurve(double l);
	double                          getHeading(double l); 
	double                          getLenght(); 
	std::vector<CMetricMapPoint>    getCurve() const;
	CMetricMapPoint                 getMax(); 
	CMetricMapPoint                 getMin(); 


	inline CMetricMapPoint getStart() { return _start; }
	inline CMetricMapPoint getCtrl1() { return _ctrl1; } 
	inline CMetricMapPoint getCtrl2() { return _ctrl2; }
	inline CMetricMapPoint getEnd()   { return _end;   }
	void setStart();
	void setCtrl1();
	void setCtrl2();
	void setEnd();

	
	double          calculatePos(double l, double P0, double P1, double P2, double P3); /// \todo Make private
	double          calculateXPos(double l);
	double          calculateYPos(double l);
	CMetricMapPoint calculateDerivative(double l);


	/* Cubic Curve transformations and manipulations */
	double  join(const CCubicCurve cubicCurveSourceFirst, const CCubicCurve cubicCurveSourceLast, CCubicCurve *cubicCurveTarget); 
	double  cubuicCurveFromMetricPoints(CMetricMapPoint before, std::vector<CMetricMapPoint> points, CMetricMapPoint after);
	void    addCurvePoint(double x, double y);
	void    convert();

	/* Path functions */
	static void     removeShortSegments(std::vector<CCubicCurve> &track, double maxLenght );
	static double   findClosesPointOnPath(std::vector<CCubicCurve> track, CMetricMapPoint point, unsigned int &index, double &parameter); 
	static double   calculateLength(std::vector<CCubicCurve> track);
	double          calculateLengthToPointOnCurve(CMetricMapPoint P);
	
	/* Debug functions */
	CCubicCurve     createRandomCubicCurve( CMetricMapPoint startPoint, double Lenght );
	CMetricMapPoint createRandomCubicCurvePoint( CMetricMapPoint startPoint, double Lenght );
	CCubicCurve     createRandomCubicCurve( CMetricMapPoint startPoint, CDirectionRad direction, double lenght );
	CMetricMapPoint createRandomCubicCurvePoint( CMetricMapPoint startPoint, CDirectionRad direction, double lenght );
	
	double minimalDistanceTo(CMetricMapPoint point);
	double minimalDistanceTo(CMetricMapPoint point, double &t);
	double minimalDistanceBetween(std::vector<CCubicCurve> cubicCurves, CMetricMapPoint point);
	double optimizeCtrlPointsDistance(CMetricMapPoint point);

	/* Various Algorithms*/
	void                createSplineFromPoints(std::vector<CMetricMapPoint> &points);
	void                elevate();
	void                splitCurve(double t);
	CMetricMapPoint     findTangent(double t);
	void                findNormal(double t);
	void                splitCurve(CCubicCurve &firstPart, CCubicCurve &lastPart, double t);
	void                splitCurveKeepFirst( double t);
	void                splitCurveKeepLast( double t);
	void                moveCurve(double x, double y);
	bool                isPointOnCurve(CMetricMapPoint P);
    void                cubicCircle();
	void                arc2curve(double R, double radStart, double radTotal, CMetricMapPoint origin);
	

	/* To be merged in*/
	/*
	double getDerivativeValue(double t) { return comp.getDerivative(1, t, x_values);  };
	CMetricMapPoint getNormal(double t)
	std::vector<CMetricMapPoint> generateBoundingBox()
	double getBoundingBoxArea();
	CCubicCurve tightBox();
	boolean hasBoxOverlap(CCubicCurve);
	CCubicCurve[] split() // Split in half
	CCubicCurve normalize()
	void startFit(std::vector<CMetricMapPoint>& p);

	*/
	
	

private:
	
	CMetricMapPoint                 _start;
	CMetricMapPoint                 _ctrl1;
	CMetricMapPoint                 _ctrl2;
	CMetricMapPoint	                _end;
	std::vector<CMetricMapPoint>    _curve; /*!< A extra datastructure to handle higher order bezier curves */
	double                          _length;         
	
	void    calculateLength();
	
	
	
	
	




};


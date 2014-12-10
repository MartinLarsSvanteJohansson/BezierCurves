#pragma once
#include "DirectionRad.h" 
/*! \class CMetricMapPoint 
 *  \brief     CMetricMapPoint is a 2-dimensional point implemented with 2 doubles X and Y (and a experimental direction dir)
 *  \details   
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
 *  \todo Upgrade the documentation
 */
class CMetricMapPoint
{
public:
	CMetricMapPoint(double X, double Y);
	CMetricMapPoint(double X, double Y, double D);
	CMetricMapPoint(void);
	~CMetricMapPoint(void);

	/* - Getters*/
	double getX() const;
	double getY() const;
	double getDir() const;

	/* - Setters */
	void setX(double X);
	void setY(double Y);
	void setDir(double D);
	void set(double X, double Y);
	void set(double X, double Y, double D);
	void set(CMetricMapPoint refPoint, CDirectionRad direction, double lenght);
	
	/* - Operator overloads*/
	CMetricMapPoint& operator=(const CMetricMapPoint& MapPoint);
	CMetricMapPoint	operator*(double z);
	CMetricMapPoint operator-(double z);
	CMetricMapPoint operator-(const CMetricMapPoint& MapPoint);
	CMetricMapPoint operator+(const CMetricMapPoint& MapPoint);
	bool			operator==(const CMetricMapPoint& MapPoint);
	bool			operator!=(const CMetricMapPoint& MapPoint);
	double			distance(CMetricMapPoint mapPoint);
	CDirectionRad	direction(CMetricMapPoint mapPoint);
	CMetricMapPoint move(double xDistande, double yDistande);
	CMetricMapPoint rotate(CDirectionRad r);


private:
	double _X;
	double _Y;
	double _dir;
};


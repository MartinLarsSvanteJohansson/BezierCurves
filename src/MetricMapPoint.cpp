/*
* MetricMapPoint.cpp
*
*
* Martin Johansson
* martin.lars.svante.johansson@gmail.com
*
*/

#include "MetricMapPoint.h"
#include <cmath>

CMetricMapPoint::CMetricMapPoint(double X, double Y) {_X = X; _Y = Y; _dir = -1;}
CMetricMapPoint::CMetricMapPoint(double X, double Y, double D) {_X = X; _Y = Y; _dir = D;}
CMetricMapPoint::CMetricMapPoint(void){}
CMetricMapPoint::~CMetricMapPoint(void){}

double CMetricMapPoint::getX() const {return _X;}
double CMetricMapPoint::getY() const {return _Y;}
double CMetricMapPoint::getDir() const {return _dir;}
void CMetricMapPoint::setX(double X){_X = X;}
void CMetricMapPoint::setY(double Y){_Y = Y;}
void CMetricMapPoint::setDir(double D){_dir = D;}
void CMetricMapPoint::set(double X, double Y){this->setX(X);this->setY(Y);}
void CMetricMapPoint::set(double X, double Y, double D){this->setX(X);this->setY(Y); this->setDir(D);}

void CMetricMapPoint::set( CMetricMapPoint refPoint, CDirectionRad direction, double lenght )
{
	this->_X = refPoint._X + lenght * cos(direction.getDouble());
	this->_Y = refPoint._Y + lenght * sin(direction.getDouble());
}

CMetricMapPoint& CMetricMapPoint::operator=( const CMetricMapPoint& MapPoint )
{
	this->setX(MapPoint._X);
	this->setY(MapPoint._Y); 
	this->setDir(MapPoint._dir);
	return *this;
}

CMetricMapPoint CMetricMapPoint::operator*(double z)
{
	return CMetricMapPoint(this->getX()*z, this->getY()*z, this->getDir()*z);
}
CMetricMapPoint CMetricMapPoint::operator-(double z)
{
	return CMetricMapPoint(this->getX() - z, this->getY() - z);
}
CMetricMapPoint CMetricMapPoint::operator-(const CMetricMapPoint& MapPoint)
{
	return CMetricMapPoint(this->_X - MapPoint._X, this->_Y - MapPoint._Y);
}
CMetricMapPoint CMetricMapPoint::operator+(const CMetricMapPoint& MapPoint)
{
	return CMetricMapPoint(this->_X + MapPoint._X, this->_Y + MapPoint._Y);
}

bool CMetricMapPoint::operator==( const CMetricMapPoint& MapPoint )
{
	//TODO: make a define for 0.1
	if((std::abs(_X - MapPoint._X) < 0.1) && (std::abs(_Y - MapPoint._Y) < 0.1))
		return true;
	else 
		return false;
}

double CMetricMapPoint::distance( CMetricMapPoint mapPoint )
{
	double dX = _X - mapPoint.getX();
	double dY = _Y - mapPoint.getY();
	return sqrt((dX*dX) + (dY*dY));
}

CDirectionRad CMetricMapPoint::direction( CMetricMapPoint mapPoint )
{
	double dX = mapPoint.getX() - _X;
	double dY = mapPoint.getY() - _Y;
	CDirectionRad angle = atan2(dY,dX);
	return angle;
}

bool CMetricMapPoint::operator!=( const CMetricMapPoint& MapPoint )
{
	//TODO: make a define for 0.1
	if((std::abs(_X - MapPoint._X) < 0.1) && (std::abs(_Y - MapPoint._Y) < 0.1)) return false;
	else return true;
}

CMetricMapPoint CMetricMapPoint::move( double xDistande, double yDistande )
{
	_X += xDistande;
	_Y += yDistande;
	return *this;
}

/* Rotates this point around origo! */

CMetricMapPoint CMetricMapPoint::rotate( CDirectionRad r )
{
	CDirectionRad newDirection = CMetricMapPoint(0,0).direction(*this) + r.getDouble();
	this->set(CMetricMapPoint(0,0),newDirection,CMetricMapPoint(0,0).distance(*this));
	return *this;
}

//CMetricMapPoint CMetricMapPoint::cubicCurvePoint()
//{
//	return CMetricMapPoint(_X,_Y);
//}



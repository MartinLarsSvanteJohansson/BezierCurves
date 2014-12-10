/*
* DirectionRad.cpp
*
*
* Martin Johansson
* martin.lars.svante.johansson@gmail.com
*
*/

#include "DirectionRad.h"

double const CDirectionRad::_pi = 3.141592653589793238462643383279502884;

CDirectionRad::CDirectionRad(void)  
{
	_direction = 0;
	_2pi =		2*_pi;
}

CDirectionRad::CDirectionRad(const CDirectionRad& directionRad)
{
	_direction = directionRad._direction;
	_2pi = 2*3.141592653589793238462643383279502884;
}

CDirectionRad::CDirectionRad( const double& directionRad )
{
	_2pi = 2*3.141592653589793238462643383279502884;
	_direction = checkOverflow(directionRad);
}

CDirectionRad::~CDirectionRad(void)
{
}

CDirectionRad& CDirectionRad::operator=(const CDirectionRad& directionRad)
{
	_direction = checkOverflow(directionRad._direction);
	return *this;
}
bool CDirectionRad::operator!=(const CDirectionRad& directionRad)
{
	return (_direction == directionRad._direction) ? false : true;
}
bool CDirectionRad::operator==(const CDirectionRad& directionRad)
{
	return (_direction == directionRad._direction) ? true : false;
}

CDirectionRad& CDirectionRad::operator+(const CDirectionRad& directionRad)
{
	_direction = checkOverflow(_direction + directionRad._direction);
	return *this;
}
CDirectionRad& CDirectionRad::operator-(const CDirectionRad& directionRad)
{
	_direction = checkOverflow(_direction - directionRad._direction);
	return *this;
}
CDirectionRad& CDirectionRad::operator*(const CDirectionRad& directionRad)
{
	_direction = checkOverflow(_direction * directionRad._direction);
	return *this;
}

CDirectionRad& CDirectionRad::operator/(const CDirectionRad& directionRad)
{
	_direction = checkOverflow(_direction / directionRad._direction);
	return *this;
}

CDirectionRad& CDirectionRad::operator+(const double& d){ _direction = checkOverflow(_direction + d); return *this; }
CDirectionRad& CDirectionRad::operator-(const double& d){ _direction = checkOverflow(_direction - d); return *this; }
CDirectionRad& CDirectionRad::operator*(const double& d){ _direction = checkOverflow(_direction * d); return *this; }
CDirectionRad& CDirectionRad::operator/(const double& d){ _direction = checkOverflow(_direction / d); return *this; }

double CDirectionRad::checkOverflow( double directionRad )
{
	while(directionRad >= _2pi) directionRad -= _2pi;
	while(directionRad < 0)     directionRad += _2pi;
	return directionRad;
}

double CDirectionRad::getDegrees() const
{
	return _direction * 180/_pi;
}

CDirectionRad& CDirectionRad::setDegrees(double degrees)
{
	_direction = checkOverflow(degrees * _pi/180);
	return *this;
}

double CDirectionRad::getDouble() const
{
	return _direction;
}

double CDirectionRad::getPi() 
{
	return _pi;
}

double CDirectionRad::get2Pi() const
{
	return _2pi;
}

double CDirectionRad::capRadian(double r)
{
		while (r >= 2.0*_pi)
			r -= 2.0*_pi;
		while (r < 0.0)
			r += 2.0*_pi;
		if (r < 0.00000001  ||  r > _2pi - 0.00000001)
			return 0.0;
		return r;
}
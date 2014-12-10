#pragma once
/*! \class CDirectionRad 
 *  \brief     CDirectionRad is a 
 *  \details   ....
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
class CDirectionRad
{
public:
	CDirectionRad(void);
	CDirectionRad(const double& directionRad);
	CDirectionRad(const CDirectionRad& directionRad);
	~CDirectionRad(void);

	CDirectionRad&  operator=(const CDirectionRad& directionRad);
	bool            operator!=(const CDirectionRad& directionRad);
	bool            operator==(const CDirectionRad& directionRad);
	CDirectionRad&  operator+(const CDirectionRad& directionRad);
	CDirectionRad&  operator-(const CDirectionRad& directionRad);
	CDirectionRad&  operator*(const CDirectionRad& directionRad);
	CDirectionRad&  operator/(const CDirectionRad& directionRad);
	CDirectionRad&  operator+(const double& d);
	CDirectionRad&  operator-(const double& d);
	CDirectionRad&  operator*(const double& d);
	CDirectionRad&  operator/(const double& d);

	double          getDegrees() const;
	static double   getPi();
	double          get2Pi() const;
	double          getDouble() const;
	CDirectionRad&  setDegrees(double degrees);
	double          capRadian(double r); 

private:
	double              _direction;
	static const double _pi;
	double              _2pi;

    double checkOverflow(double directionRad);

};


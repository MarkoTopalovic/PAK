#pragma once

#include <ctime>

class TSValue
{
public:
	static const double NO_DATA_VALUE;

	TSValue();
	TSValue(double dt, double value, double accuracy);
	TSValue(double dt);

	bool hasValue();
	double getValue();
	double getAccuracy();
	double getTime();

	bool HasValue;
	double value;
	double accuracy;
	double time;  
};
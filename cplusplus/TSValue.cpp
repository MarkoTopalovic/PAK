#include "TSValue.h"

const double TSValue::NO_DATA_VALUE=-9999.0;

TSValue::TSValue()
{
}

TSValue::TSValue(double dt, double value, double accuracy)
{
	this->time = dt;
	this->value = value;
	this->accuracy = accuracy;
	this->HasValue = true;
}

TSValue::TSValue(double dt)
{
	this->HasValue = false;
	this->time = dt;
}

bool TSValue::hasValue()
{
	return HasValue;
}

double TSValue::getValue()
{
	return value;
}

double TSValue::getAccuracy()
{
	return accuracy;
}

double TSValue::getTime()
{
	return time;
}
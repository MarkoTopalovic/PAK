#pragma once

#include <iostream>
#include <vector>
#include "TSValue.h"

#include "math.h"
#include <string>
#include <stdio.h>

using namespace std;

enum TimeUnit
{
	HOUR,
	DAY,
	MONTH
};

class TimeSeries
{
	public:
		TimeSeries(void);

		void addValue(double dt);
		void addValue(struct tm dt);

		void addValue(double dt, double value, double accuracy);
		void addValue(struct tm dt, double value, double accuracy);
		
		TSValue& getValue(int index);
		int getCount();
		std::vector<TSValue>& getValues();
		
		bool getY(double time, double *y);

		void toFEQ(string path);
		std::string toFEQ();

		TimeSeries fromFEQ(string path);
		TimeSeries fromFEQ(FILE *f);

		std::vector<TSValue> m_values;
		
		static struct tm* no_dts_localtime(time_t *rawtime);
};


#include "TimeSeries.h"

TimeSeries::TimeSeries(void){}

void TimeSeries::addValue(double dt)
{
	m_values.push_back(TSValue(dt));
}

void TimeSeries::addValue(struct tm dt)
{
	time_t t_t = mktime(&dt);
	m_values.push_back(TSValue(t_t));
}

void TimeSeries::addValue(double dt, double value, double accuracy)
{
	m_values.push_back(TSValue(dt, value, accuracy));
}

void TimeSeries::addValue(struct tm dt, double value, double accuracy)
{
	time_t t_t = mktime(&dt);

	m_values.push_back(TSValue(t_t, value, accuracy));
}

TSValue& TimeSeries::getValue(int index)
{
	return m_values[index];
}

int TimeSeries::getCount()
{
	return m_values.size();
}

std::vector<TSValue>& TimeSeries::getValues()
{
	return m_values;
}
		
bool TimeSeries::getY(double time, double *y)
{
	int i;
    int nPointCount;
    TSValue p1, p2;

    *y = 0.0;

    nPointCount = getCount();

	if(nPointCount >= 1)
		if(time < getValue(0).getTime())
			return false;

    if (nPointCount > 1)
    {	
        if (time >= getValue(nPointCount - 1).getTime())
        {
            if (nPointCount > 1)
            {
                p2 = getValue(nPointCount - 1);
				*y = p2.getValue();

                return true;
            }
            else
            {
                return false;
            }
        }
		
		int first = 1;
		int last = nPointCount - 1;
		int middle = (first+last)/2;

		double dwa = this->m_values[0].time;

		if(getValue(0).getTime() < getValue(1).getTime())
		{
			while (first <= last) 
			{
				if (getValue(middle).getTime() >= time && getValue(middle-1).getTime() <= time)
					break;    
				else if (getValue(middle).getTime() > time)
					last = middle - 1;
				else
					first = middle + 1;
 
				middle = (first + last)/2;
			}
			
			if(middle == 0)
				p1=getValue(middle);
			else
				p1=getValue(middle-1);

			p2=getValue(middle);
		}
		else
		{
			while (first < last) {
				if (getValue(middle).getTime() <= time && getValue(middle-1).getTime() >= time)
					break;    
				else if (getValue(middle).getTime() < time) {
					last = middle - 1;
				}
				else
					first = middle + 1;
 
				middle = (first + last)/2;
			}

			p2=getValue(middle-1);
			p1=getValue(middle);
		}

		if(p1.time != p2.time)
			*y = double(time - p1.getTime()) * (p2.getValue() - p1.getValue()) / (p2.getTime() - p1.getTime()) + p1.getValue();
		else
			*y = p2.getValue();

        return true;
    }

    *y = getValue(0).getValue();
    return false;
}

void TimeSeries::toFEQ(string path)
{
	FILE *f = fopen(path.c_str(), "w");

	string feqCont = toFEQ();

	fprintf(f, "%s", feqCont.c_str());

	fflush(f);

	fclose(f);
}

std::string TimeSeries::toFEQ()
{
	string header;

	header = "TOTAL FLOW\nREFL=0.0        FAC=1.0\nYEAR MN DY      HOUR DISCHARGE  INFLOW HYDROGRAPH AT UPSTREAM END\n";
	
	double hour = 0.0;
	double discharge = 0.0;

	TSValue val;

	for(int j = 0; j < getCount(); j++)
	{
		double hour = 0.0;
		double discharge = 0.0;

		val = getValue(j);

		time_t rawtime;
		struct tm *ptm;
		rawtime = val.getTime();

		ptm = no_dts_localtime(&rawtime);
		 
		char content[500];

		double minutes;
		if(ptm->tm_min != 0)
			minutes = (ptm->tm_min / 60.0);
		else
			minutes = 0;

		double seconds;
		if(ptm->tm_sec != 0)
			seconds = (ptm->tm_sec / 3600.0);
		else
			seconds = 0;

		hour = ptm->tm_hour;

		sprintf(content,  "%4d%3d%3d %9.5lf %9.2lf\n", ptm->tm_year + 1900, ptm->tm_mon + 1, ptm->tm_mday, hour + minutes + seconds, val.value);

		header += content;

		if(j == getCount() - 1)
		{
 			sprintf(content,  "%4d%3d%3d%10.5f\n", ptm->tm_year + 1900, ptm->tm_mon + 1, ptm->tm_mday, hour + minutes);
			header += content;
		}
	}

	return header;
}

TimeSeries TimeSeries::fromFEQ(string path)
{
	FILE *f = fopen(path.c_str(), "r");

	TimeSeries ts = ts.fromFEQ(f);
	fclose(f);

	return ts;
}

TimeSeries TimeSeries::fromFEQ(FILE *f)
{
	TimeSeries ts;

	char s[255];

	for(int i = 0; i < 3; i++)
		fgets(s, 255, f);


	double hour = -1, last_hour = -2;
	double val = 0;

	time_t rawtime = time(0);
	struct tm *ptm;
	ptm = no_dts_localtime(&rawtime);

	ptm->tm_min = 0;
	ptm->tm_sec = 0;

	int y, d = -2, m;

	while(!(d == ptm->tm_mday && hour <= last_hour))
	{
		last_hour = hour;
		d = ptm->tm_mday;

		int testScan = fscanf(f, "%d%d%d%lf%lf", &(ptm->tm_year), &(ptm->tm_mon), &(ptm->tm_mday), &hour, &val);
		if(testScan == 0)
		{
			fgets(s, 255, f);
			last_hour = -2;
			continue;
		}

		ptm->tm_mon -= 1;

		ptm->tm_year -= 1900;
		ptm->tm_hour = (int)floor(hour);

		ptm->tm_min = (int) floor(hour - ptm->tm_hour);

		double seconds = hour - (ptm->tm_hour + floor((ptm->tm_min / 60.0)));

		ptm->tm_sec = (int) floor((seconds * 3600) + 0.5);

		ptm->tm_isdst = 0;
		
		time_t t_t = mktime(ptm);

		if(!(d == ptm->tm_mday && hour <= last_hour))
			ts.addValue(*ptm, val, 0.0);
	}

	return ts;
}

struct tm* TimeSeries::no_dts_localtime(time_t *rawtime)
{
	struct tm *dt = localtime(rawtime);
	dt->tm_hour -= dt->tm_isdst;
	mktime(dt);

	dt->tm_isdst = 0;

	return dt;
}

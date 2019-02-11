#include "TimeSeries.h"
#include <stdint.h>
#include <cstring>
#include <iostream>
#include <vector>
#include "math.h"
#include <stdio.h>
#include "NamingStrategy.h"
//#include "TSValue.h"


extern "C" void save_series_(char *ime, double* vreme, double* serie, int64_t &serie_steps, int64_t &serie_pts, int &current, int &n, int &ind);

void save_series_(char *ime, double* vreme, double* serie, int64_t &serie_steps, int64_t &serie_pts, int &current, int &n, int &ind)
{//CALL save_series(MPOINT_ID(IMP),MP_VREME, MP_RESULTS_NIZ,     BRKORAKA,            MAX_MPOINTS,       IMP,      INDMPI,   isave)
    int64_t i,start;
    TimeSeries ts;

    if(ind == 1)
    {
        start=0;
    }else{
        start=1;
    }
    for (i=start ; i<serie_steps ; i++) 
    {
        ts.addValue(vreme[i]*86400,serie[i],0.0);
    }
	char imeC[100000]; //maksimalna duzina imena
	for (i = 15 - n; i<15; i++)
	{
		imeC[i - 15 + n] = ime[i];
	}

	std::string a(imeC);
	std::string st;
	st = a.substr(0, a.size() - (a.size() - n));
	std::string b;
	b = "_pak.tab";
	std::string imeCpak = st + b;
	ts.toFEQ(imeCpak);
}


#include <iostream>
#include <string>
#include <sstream>
#include <unistd.h>
#include <sys/times.h>
#include "BGAutils.hpp"

using namespace std;

ofstream& mktempofstream(ofstream& out, char *writefile)
{
    const int max_attempts = 10;
    char *pX = strstr(writefile, "XXXXXX");
    if (!pX)
    {
	throw IOError("mktempofstream(): XXXXXX missing in filename");
    }
    string alnum = "0123456789"
	"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    for (int i = 0; ; ++i)
    {
	// create random name
	for (int j = 0; j < 6; ++j)
	    pX[j] = alnum[ rand() % alnum.size() ];
	ifstream fid(writefile);
	if (!fid)
	    break;
	else
	    fid.close();
	if ( !(i < max_attempts) )
	    throw IOError("mktempofstream(): max_attempts exceeded");
    }
    // we should have a good name here
    out.open(writefile);
    if (!out)
    {
	ostringstream oss;
	oss << "mktempofstream(): unable to open " <<
				  writefile << "for writing";
	throw IOError(oss.str());
    }
    return out;
}

double BGA::CPUTime()
{
    tms tbuf;
    times(&tbuf);
    return 1.0*tbuf.tms_utime/sysconf(_SC_CLK_TCK);
}

void Counters_t::PrintRunStats()
{
    char hostname[255];
    gethostname(hostname, 255);
    cout << "Run statistics:" << endl;
    cout << "penalty_calls = " << penalty_calls << endl;
    cout << "distance_calls = " << distance_calls << endl;
    cout << "UserCPUtime = " << BGA::CPUTime() << 's' << endl;
    cout << "Host = " << hostname << endl;
}

// define global instance of Counters_t
Counters_t BGA::cnt;

double vdnorm(const valarray<double>& v)
{
    double s2 = 0.0;
    for (const double *pv = &v[0]; pv != &v[v.size()]; ++pv)
	s2 += (*pv)*(*pv);
    return sqrt(s2);
}

double vddot(const valarray<double>& v1, const valarray<double>& v2)
{
    // assuming v1.size() == v2.size()
    double dotprod = 0.0;
    const double *pv1 = &v1[0], *pv2 = &v2[0];
    for (; pv1 != &v1[v1.size()]; ++pv1, ++pv2)
	dotprod += (*pv1)*(*pv2);
    return dotprod;
}

valarray<double> vdcross(const valarray<double>& v1, const valarray<double>& v2)
{
    if (v1.size() != 3 || v2.size() != 3)
    {
	throw runtime_error("vdcross(): invalid valarray size");
    }
    valarray<double> cross(3);
    cross[0] = v1[1]*v2[2] - v1[2]*v2[1];
    cross[1] = v1[2]*v2[0] - v1[0]*v2[2];
    cross[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return cross;
}

valarray<double> vdrecipw0(const valarray<double>& v)
{
    // calculate reciprocal value, with checking for zeros
    const double zero_recip_gain = 10;
    double min_positive = 0.0;
    for (const double *p = &v[0]; p != &v[v.size()]; ++p)
    {
	if (*p > 0 && (!min_positive || *p < min_positive))
	    min_positive = *p;
    }
    if (!min_positive)
	min_positive = 1.0;
    valarray<double> recip(v.size());
    double *r = &recip[0];
    for (const double *p = &v[0]; p != &v[v.size()]; ++p, ++r)
	*r = (*p != 0) ? 1.0/(*p) : zero_recip_gain*1.0/min_positive;
    return recip;
}

// End of file

/***********************************************************************
* Short Title: class DistanceTable - definitions
*
* Comments:
*
* $Id$
*
* <license text>
***********************************************************************/

#include <sstream>
#include <cmath>
#include "Exceptions.hpp"
#include "DistanceTable.hpp"

RegisterSVNId DistanceTable_cpp_id("$Id$");

using namespace std;

////////////////////////////////////////////////////////////////////////
// class DistanceTable
////////////////////////////////////////////////////////////////////////

// constructors

DistanceTable::DistanceTable() : vector<double>()
{
    init();
}

DistanceTable::DistanceTable(const double* v, size_t sz) :
    vector<double>(v, v + sz)
{
    init();
}

DistanceTable::DistanceTable(const vector<double>& v) : vector<double>(v)
{
    init();
}

DistanceTable::DistanceTable(const DistanceTable& d0)
{
    *this = d0;
}

DistanceTable& DistanceTable::operator= (const vector<double>& v)
{
    vector<double>& this_vector = *this;
    this_vector = v;
    init();
    return *this;
}

DistanceTable& DistanceTable::operator= (const DistanceTable& d0)
{
    if (this == &d0)    return *this;
    assign(d0.begin(), d0.end());
    est_num_atoms = d0.est_num_atoms;
    count_unique = d0.count_unique;
    resolution = d0.resolution;
    return *this;
}

DistanceTable::const_iterator
DistanceTable::find_nearest(const double& dfind) const
{
    const_iterator ii = lower_bound(begin(), end(), dfind);
    if (    ( ii == end() && size() != 0 ) ||
	    ( ii != begin() && (dfind - *(ii-1)) < (*ii - dfind) )
       )
	--ii;
    return ii;
}

DistanceTable::iterator DistanceTable::return_back(const double& dback)
{
    iterator ii = lower_bound(begin(), end(), dback);
    return insert(ii, dback);
}

int DistanceTable::estNumAtoms() const
{
    return est_num_atoms;
}

int DistanceTable::countUnique() const
{
    if (count_unique >= 0)  return count_unique;
    count_unique = empty() ? 0 : 1;
    const_iterator ilo = begin();
    const_iterator ihi = (ilo == end()) ? ilo : ilo + 1;
    for (; ihi != end(); ++ilo, ++ihi)
    {
	if ( (*ihi - *ilo) > resolution )   ++count_unique;
    }
    return count_unique;
}

vector<double> DistanceTable::unique() const
{
    vector<double> dtu(countUnique());
    vector<double>::iterator dui = dtu.begin();
    double d0 = -1.0;
    for (const_iterator di = begin(); di != end() && dui != dtu.end(); ++di)
    {
	if ( (*di - d0) > resolution )
	{
	    *(dui++) = *di;
	    d0 = *di;
	}
    }
    return dtu;
}

double DistanceTable::getResolution() const
{
    return resolution;
}

void DistanceTable::setResolution(double res)
{
    resolution = res;
    count_unique = -1;
}

// private methods

void DistanceTable::init()
{
    double xn = ( 1.0 + sqrt(1.0 + 8.0*size()) )/2.0;
    est_num_atoms = int(xn);
    count_unique = -1;
    resolution = 1e-5;  // same as sqrt(BGA::eps_badness)
    if (empty())    return;
    // check for any negative element
    size_t minidx = min_element(begin(), end()) - begin();
    if (at(minidx) < 0.0)
    {
        ostringstream emsg;
        emsg << "DistanceTable::init() negative distance at " << minidx;
        throw InvalidDistanceTable(emsg.str());
    }
    // sort and count unique distance
    sort(begin(), end());
}

// non-member operators

istream& operator>>(istream& fid, DistanceTable& dstbl)
{
    vector<double> data_buffer;
    // skip header
    string word;
    while (!fid.eof() && fid >> word && !word.empty() && word[0] == '#')
    {}
    istringstream wstrm(word);
    double x;
    if (wstrm >> x)     data_buffer.push_back(x);
    while (fid >> x)    data_buffer.push_back(x);
    dstbl.swap(data_buffer);
    dstbl.init();
    return fid;
}

// End of file

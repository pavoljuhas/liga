/***********************************************************************
* Short Title: object definitions for Biosphere Genetic Algorithm
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/


//#include <iostream>
//#include <string>
//#include <sstream>
//#include <fstream>
//#include <valarray>
//#include <vector>
//#include <list>

#include <stdexcept>
#include "BGAlib.hpp"

// exceptions
struct InvalidDistanceTable { };
struct InvalidMolecule { };
struct IOError { };

// read lines that do not start with number
bool read_header(istream& fid, string& header)
{
    double x;
    string line;
    istringstream istrs;
    header.clear();
    for (   istream::pos_type p = fid.tellg();
	    !fid.eof() && getline(fid, line);
	    p = fid.tellg()
	)
    {
	istrs.clear();
	istrs.str(line);
	if (istrs >> x)
	{
	    fid.seekg(p);
	    break;
	}
	else
	{
	    header.append(line + '\n');
	}
    }
    return !(fid.rdstate() & ios::badbit);
}

inline bool read_header(istream& fid)
{
    string dummy;
    return read_header(fid, dummy);
}

// read as many numbers as possible
template<class T> bool read_data(istream& fid, vector<T>& v)
{
    // prepare v
    T x;
    while (fid >> x)
    {
	v.push_back(x);
    }
    return !(fid.rdstate() & ios::badbit);
}

////////////////////////////////////////////////////////////////////////
// SandSphere definitions
////////////////////////////////////////////////////////////////////////
SandSphere::SandSphere(int GridMax, const vector<double>& vt) :
    gridmax(GridMax)
{
    init(vt);
}

SandSphere::SandSphere(int GridMax, int s, const double *pt) :
    gridmax(GridMax)
{
    vector<double> vt(s);
    for (int i = 0; i < s; ++i) { vt[i] = pt[i]; }
    init(vt);
}

SandSphere::SandSphere(int GridMax, const char *file) :
    gridmax(GridMax)
{
    // open file for reading
    ifstream fid(file);
    if (!fid)
    {
	cerr << "E: unable to read '" << file << "'" << endl;
	throw IOError();
    }
    // read values to vt
    vector<double> vt;
    bool result = read_header(fid) && read_data(fid, vt);
    // check if everything was read
    if ( !result || !fid.eof() )
    {
	cerr << "E: " << file << ':' << fid.tellg() <<
	    ": error reading SandSphere" << endl;
	throw IOError();
    }
    fid.close();
    init(vt);
}

void SandSphere::init(const vector<double>& t)
{
    NDist = t.size();
    if (NDist == 0)
    {
	cerr << "E: target distance table is empty" << endl;
	throw InvalidDistanceTable();
    }
    // calculate and check NAtoms
    double xNAtoms = 0.5 + sqrt(1 + 8.0*NDist)/2.0;
    NAtoms = int(xNAtoms);
    if (double(NAtoms) != xNAtoms)
    {
	cerr << "E: incorrect length of target distance table, NAtoms=" <<
		xNAtoms << '\n';
	throw InvalidDistanceTable();
    }
    // fill in and check distance valarray d
    d.resize(NDist);
    for (int i = 0; i < NDist; ++i)
    {
	d[i] = t[i];
    }
    sort(&d[0], &d[d.size()]);
    if (d[0] <= 0)
    {
	cerr << "E: non-positive entry in target distance table, " <<
	    "d[0]=" << d[0] << '\n';
	throw InvalidDistanceTable();
    }
    // calculate grid parameters
    dmax = d[d.size()-1];
    delta = dmax/gridmax;
    // calculate d2, d2lo, d2hi
    d2.resize(NDist);
    d2lo.resize(NDist);
    d2hi.resize(NDist);
    d /= delta;
    d2 = d*d;
    // take care of grid tolerance:
    SetGridTol(defGridTol);
}

void SandSphere::SetGridTol(double t)
{
    vGridTol = t;
    for (int i = 0; i < NDist; ++i)
    {
	d2lo[i] = (d[i] < vGridTol) ? 0.0 : pow(d[i] - vGridTol, 2);
	d2hi[i] = pow(d[i] + vGridTol, 2);
    }
    for ( list<Molecule *>::const_iterator i = molecules.begin();
	    i != molecules.end(); ++i )
    {
	(*i)->UnCache();
    }
}

double SandSphere::GridTol()
{
    return vGridTol;
}

////////////////////////////////////////////////////////////////////////
// Molecule definitions
////////////////////////////////////////////////////////////////////////
Molecule::Molecule(SandSphere *SS) : ss(SS)
{
    init();
}

Molecule::Molecule(SandSphere *SS,
	int s, int *ph, int *pk) : ss(SS)
{
    h.resize(s);
    k.resize(s);
    for (int i = 0; i < s; ++i)
    {
	h[i] = ph[i];
	k[i] = pk[i];
    }
    init();
}

Molecule::Molecule(SandSphere *SS,
	const vector<int>& vh, const vector<int>& vk) : ss(SS)
{
    h = vh;
    k = vk;
    init();
}

Molecule::Molecule(SandSphere *SS,
	int s, double *px, double *py) : ss(SS)
{
    h.resize(s);
    k.resize(s);
    for (int i = 0; i < s; ++i)
    {
	h[i] = (int) round(px[i] / ss->delta);
	k[i] = (int) round(py[i] / ss->delta);
    }
    init();
}

Molecule::Molecule(SandSphere *SS,
	const vector<double>& vx, const vector<double>& vy) : ss(SS)
{
    h.resize(vx.size());
    k.resize(vx.size());
    for (size_t i = 0; i < h.size(); ++i)
    {
	h[i] = (int) round(vx[i] / ss->delta);
	k[i] = (int) round(vy[i] / ss->delta);
    }
    init();
}

Molecule::Molecule(const Molecule& M) :
    ss(M.ss)
{
    init();
    *this  = M;
}

Molecule& Molecule::operator=(const Molecule& M)
{
    if (this == &M) return *this;
    // data storage
    if (ss != M.ss)
    {
	ss->molecules.remove(this);
	ss = M.ss;
	ss->molecules.push_back(this);
    }
    // h, k assignment must preceed fix_size()
    h = M.h;			// x-coordinates
    k = M.k;			// y-coordinates
    // parameters
    if (NAtoms != M.NAtoms)
    {
	// this sets NAtoms, NDist, resizes all valarray's
	fix_size();
    }
    // badness evaluation
    cached = M.cached;
    if (M.cached)
    {
	abad = M.abad;		// individual atom badnesses
	abadMax = M.abadMax;	// maximum atom badness
	mbad = M.mbad;		// molecular badness
	d2 = M.d2;		// sorted table of squared distances 
	ssdIdxFree = M.ssdIdxFree;	// available elements in ss.dist
    }
    // IO helpers
    output_format = M.output_format;
    opened_file = M.opened_file;
    return *this;
}

void Molecule::init()
{
    ss->molecules.push_back(this);
    // check coordinate sizes
    OutFmtGrid();
    UnCache();
    fix_size();
}

void Molecule::fix_size()
{
    UnCache();
    if (h.size() != k.size())
    {
	cerr << "E: invalide coordinate vectors" << endl;
	throw InvalidMolecule();
    }
    NAtoms = h.size();
    NDist  = NAtoms*(NAtoms-1)/2;
    abad.resize(NAtoms, 0.0);
    abadMax = 0.0;
    d2.resize(NDist, 0);
//    ssdIdxUsed.clear();
    ssdIdxFree.clear();
}

Molecule::~Molecule()
{
    // debug: cout << "ss->molecules.size() = " << ss->molecules.size() << endl;
    ss->molecules.remove(this);
}

double Molecule::dist(const int& i, const int& j)
{
    return sqrt(1.0*dist2(i, j));
}

double Molecule::out_penalty(int i)
{
    double Rout = sqrt(h[i]*h[i] + k[i]*k[i] + 0.0) - ss->gridmax;
    return (Rout > 0.0) ? Rout : 0.0;
}

typedef struct {
    int d2;
    int i;
    int j;
} d2idx_type;

bool d2idx_cmp(const d2idx_type& p, const d2idx_type& q)
{
    return (p.d2 < q.d2);
}

void Molecule::calc_df()
{
    cached = true;
//    ssdIdxUsed.clear();
    ssdIdxFree.clear();
    d2idx_type d2idx[NDist];
    // check if molecule is not too large
    if (NDist > ss->NDist)
    {
	cerr << "E: molecule too large" << endl;
	throw InvalidMolecule();
    }
    // calculate and store distances
    int ij = 0;
    for (int i = 0; i < NAtoms; ++i)
    {
	for (int j = i + 1; j < NAtoms; ++j, ++ij)
	{
	    d2idx[ij].d2 = dist2(i, j);
	    d2idx[ij].i = i;
	    d2idx[ij].j = j;
	}
    }
    sort(d2idx, d2idx+NDist, d2idx_cmp);
    for (ij = 0; ij < NDist; ++ij)
    {
	d2[ij] = d2idx[ij].d2;
    }
    // evaluate abad[i]
    abad = 0.0;
    int ssdIdx = 0;
    for (ij = 0; ij < NDist; ++ij)
    {
	for(; ss->d2hi[ssdIdx] < d2idx[ij].d2 && ssdIdx < ss->NDist; ++ssdIdx)
	{
	    ssdIdxFree.push_back(ssdIdx);
	}
	// did we find a baddie?
	if (ssdIdx < ss->NDist  &&  ss->d2lo[ssdIdx] > d2idx[ij].d2)
	{
	    abad[d2idx[ij].i]++;
	    abad[d2idx[ij].j]++;
	}
	// otherwise it is a matching distance
	else
	{
	    ssdIdx++;
	}
    }
    // now add penalty for being outside the SandSphere
    for (int i = 0; i < NAtoms; ++i)
    {
	abad[i] += out_penalty(i);
    }
    mbad = abad.sum();
    abadMax = (NAtoms > 0) ? max(abad.max(), (double) NAtoms) : 0.0;
}

double Molecule::ABadness(int i)
{
    if (!cached)
    {
	calc_df();
    }
    return abad[i];
}

double Molecule::AFitness(int i)
{
    // this will update abadMax if necessary
    double badness_i = ABadness(i);
    return abadMax - badness_i;
}

double Molecule::MBadness()
{
    if (!cached)
    {
	calc_df();
    }
    return mbad;
}

double Molecule::MFitness()
{
    // this will update abadMax if necessary
    double mbadness = MBadness();
    return NAtoms*abadMax - mbadness;
}

////////////////////////////////////////////////////////////////////////
// Molecule operators
////////////////////////////////////////////////////////////////////////
Molecule& Molecule::Shift(int dh, int dk)
{
    for (int i = 0; i < NAtoms; ++i)
    {
	if (cached)
	{
	    abad[i] -= out_penalty(i);
	}
	h[i] += dh;
	k[i] += dk;
	if (cached)
	{
	    abad[i] += out_penalty(i);
	}
    }
    if (cached)
    {
	abadMax = (NAtoms > 0) ? max(abad.max(), (double) NAtoms) : 0.0;
	mbad = abad.sum();
    }
    return *this;
}

Molecule& Molecule::Center()
{
    double mean_h = 0.0;
    double mean_k = 0.0;
    for (int i = 0; i < NAtoms; ++i)
    {
	mean_h += h[i];
	mean_k += k[i];
    }
    mean_h /= NAtoms;
    mean_k /= NAtoms;
    Shift( (int) round(-mean_h), (int) round(-mean_k) );
    return *this;
}

Molecule& Molecule::Part(const Molecule& M, const int cidx)
{
    h.resize(1, M.h[cidx]);
    k.resize(1, M.k[cidx]);
    fix_size();
    return *this;
}

Molecule& Molecule::Part(const Molecule& M, const list<int>& cidx)
{
    vector<int> *h_new, *k_new;
    if (&M == this)
    {
	h_new = new vector<int>;
	k_new = new vector<int>;
    }
    else
    {
	h_new = &h;
	k_new = &k;
    }
    h_new->resize(cidx.size());
    k_new->resize(cidx.size());
    int idx = 0;
    for ( list<int>::const_iterator li = cidx.begin();
	    li != cidx.end(); ++li )
    {
	if (*li < 0 || *li >= M.NAtoms)
	{
	    throw range_error("in Molecule::Part()");
	}
	(*h_new)[idx] = M.h[*li];
	(*k_new)[idx] = M.k[*li];
	++idx;
    }
    if (&M == this)
    {
	h = *h_new;
	k = *k_new;
	delete h_new, k_new;
    }
    fix_size();
    return *this;
}

Molecule& Molecule::Pop(const Molecule& M, const int cidx)
{
    if (cidx < 0 || cidx >= M.NAtoms)
    {
	throw range_error("in Molecule::Pop(list<int>)");
    }
    if (this != &M)  *this = M;
    h.erase(h.begin() + cidx);
    k.erase(k.begin() + cidx);
    fix_size();
    return *this;
}

Molecule& Molecule::Pop(const Molecule& M, const list<int>& cidx)
{
    if (cidx.size() == 0)
    {
	if (this != &M)  *this = M;
	return *this;
    }
    list<int> sidx(cidx);
    sidx.sort();
    if (sidx.front() < 0 || sidx.back() >= M.NAtoms)
    {
	throw range_error("in Molecule::Pop(list<int>)");
    }
    sidx.push_back(M.NAtoms);
    vector<int> h_new, k_new;
    int j = 0;
    for ( list<int>::iterator li = sidx.begin();
	    li != sidx.end(); ++li )
    {
	for (; j < *li; ++j)
	{
	    h_new.push_back(M.h[j]);
	    k_new.push_back(M.k[j]);
	}
	j = *li + 1;
    }
    h = h_new;
    k = k_new;
    fix_size();
    return *this;
}

Molecule& Molecule::Clear()
{
    h.clear();
    k.clear();
    fix_size();
    return *this;
}

Molecule& Molecule::Add(Molecule& m)
{
    for (int i = 0; i < m.NAtoms; ++i)
    {
	h.push_back(m.h[i]);
	k.push_back(m.k[i]);
    }
    fix_size();
    return *this;
}

Molecule& Molecule::Add(int nh, int nk)
{
    h.push_back(nh);
    k.push_back(nk);
    fix_size();
    return *this;
}

Molecule& Molecule::MoveAtomTo(int idx, int nh, int nk)
{
    if (idx >= NAtoms)
    {
	throw range_error("in Molecule::MoveAtomTo()");
    }
    UnCache();
    h[idx] = nh;
    k[idx] = nk;
    return *this;
}

////////////////////////////////////////////////////////////////////////
// Molecule IO functions
////////////////////////////////////////////////////////////////////////

Molecule::ParseHeader::ParseHeader(const string& s) : header(s)
{
    // parse format
    string fmt;
    // initialize members:
    state = 
	read_token("BGA molecule format", fmt) &&
	read_token("NAtoms", NAtoms) &&
	read_token("delta", delta);
    if (!state)
    {
	return;
    }
    if (fmt == "grid")
	format = GRID;
    else if (fmt == "xy")
	format = XY;
    else if (fmt == "atomeye")
	format = ATOMEYE;
    else
    {
	state = false;
	return;
    }
}

template<class T> bool Molecule::ParseHeader::read_token(
	const char *token, T& value
	)
{
    const char *fieldsep = ":= ";
    int ltoken = strlen(token);
    string::size_type sp;
    const string::size_type npos = string::npos;
    if ( 
	    npos == (sp = header.find(token)) || 
	    npos == (sp = header.find_first_not_of(fieldsep, sp+ltoken))
       )
    {
	return false;
    }
    istringstream istrs(header.substr(sp));
    bool result = (istrs >> value);
    return result;
}

istream& Molecule::ReadGrid(istream& fid)
{
    // read values to integer vector vhk
    string header;
    vector<int> vhk;
    bool result = read_header(fid, header) && read_data(fid, vhk);
    if (!result) return fid;
    // parse header
    double vhk_scale = 1.0;
    int vhk_NAtoms = vhk.size()/2;
    ParseHeader ph(header);
    if (ph)
    {
	vhk_scale = ph.delta / ss->delta;
	if ( vhk_NAtoms != ph.NAtoms )
	{
	    cerr << "E: " << opened_file << ": expected " << ph.NAtoms <<
		" atoms, read " << vhk_NAtoms << endl;
	    throw IOError();
	}
    }
    // check if all coordinates have been read
    if ( vhk.size() % 2 )
    {
	cerr << "E: " << opened_file << ": incomplete data" << endl;
	throw IOError();
    }
    Clear();
    h.resize(vhk_NAtoms);
    k.resize(vhk_NAtoms);
    for (int i = 0, iv = 0; i < vhk_NAtoms; ++i)
    {
	h[i] = vhk[iv++];
	k[i] = vhk[iv++];
    }
    fix_size();
    return fid;
}

bool Molecule::ReadGrid(const char* file)
{
    // open file for reading
    ifstream fid(file);
    if (!fid)
    {
	cerr << "E: unable to read '" << file << "'" << endl;
	throw IOError();
    }
    opened_file = file;
    bool result = ReadGrid(fid);
    opened_file.clear();
    fid.close();
    return result;
}

istream& Molecule::ReadXY(istream& fid)
{
    // read values to integer vector vxy
    string header;
    vector<double> vxy;
    bool result = read_header(fid, header) && read_data(fid, vxy);
    if (!result) return fid;
    int vxy_NAtoms = vxy.size()/2;
    // check how many numbers were read
    ParseHeader ph(header);
    if (ph && vxy_NAtoms != ph.NAtoms)
    {
	cerr << "E: " << opened_file << ": expected " << ph.NAtoms <<
	    " atoms, read " << vxy_NAtoms << endl;
	throw IOError();
    }
    if ( vxy.size() % 2 )
    {
	cerr << "E: " << opened_file << ": incomplete data" << endl;
	throw IOError();
    }
    Clear();
    h.resize(vxy.size()/2);
    k.resize(vxy.size()/2);
    for (size_t i = 0, iv = 0; i < vxy.size()/2; ++i)
    {
	h[i] = (int) round(vxy[iv++] / ss->delta);
	k[i] = (int) round(vxy[iv++] / ss->delta);
    }
    fix_size();
    return fid;
}

bool Molecule::ReadXY(const char* file)
{
    // open file for reading
    ifstream fid(file);
    if (!fid)
    {
	cerr << "E: unable to read '" << file << "'" << endl;
	throw IOError();
    }
    opened_file = file;
    bool result = ReadXY(fid);
    opened_file.clear();
    fid.close();
    return result;
}

bool write_file(const char* file, Molecule& m)
{
    // open file for writing
    ofstream fid(file);
    if (!fid)
    {
	cerr << "E: unable to write to '" << file << "'" << endl;
	throw IOError();
    }
    bool result = (fid << m);
    fid.close();
    return result;
}

bool Molecule::WriteGrid(const char* file)
{
    file_fmt_type org_ofmt = output_format;
    OutFmtGrid();
    bool result = write_file(file, *this);
    output_format = org_ofmt;
    return result;
}

bool Molecule::WriteXY(const char* file)
{
    file_fmt_type org_ofmt = output_format;
    OutFmtXY();
    bool result = write_file(file, *this);
    output_format = org_ofmt;
    return result;
}

bool Molecule::WriteAtomEye(const char* file)
{
    file_fmt_type org_ofmt = output_format;
    OutFmtAtomEye();
    bool result = write_file(file, *this);
    output_format = org_ofmt;
    return result;
}

Molecule& Molecule::OutFmtGrid()
{
    output_format = GRID;
    return *this;
}

Molecule& Molecule::OutFmtXY()
{
    output_format = XY;
    return *this;
}

Molecule& Molecule::OutFmtAtomEye()
{
    output_format = ATOMEYE;
    return *this;
}

istream& operator>>(istream& fid, Molecule& m)
{
    string header;
    istream::pos_type p = fid.tellg();
    if( !read_header(fid, header) )
    {
	fid.setstate(ios_base::failbit);
	return fid;
    }
    fid.seekg(p);
    Molecule::ParseHeader ph(header);
    if (!ph)
    {
	fid.setstate(ios_base::failbit);
	return fid;
    }
    bool result;
    switch (ph.format)
    {
	case m.GRID:
	    result = m.ReadGrid(fid);
	    break;
	case m.XY:
	    result = m.ReadXY(fid);
	    break;
	case m.ATOMEYE:
	    throw runtime_error("reading of atomeye files not implemented");
	    break;
    }
    if (!result)
    {
	fid.setstate(ios_base::failbit);
    }
    return fid;
}

ostream& operator<<(ostream& fid, Molecule& m)
{
    switch (m.output_format)
    {
	case m.GRID:
	    fid << "# BGA molecule format = grid" << endl;
	    fid << "# NAtoms = " << m.NAtoms << endl;
	    fid << "# delta = " << m.ss->delta << endl;
	    for (int i = 0; i < m.NAtoms; ++i)
	    {
		fid << m.h[i] << '\t' << m.k[i] << endl;
	    }
	    break;
	case m.XY:
	    fid << "# BGA molecule format = xy" << endl;
	    fid << "# NAtoms = " << m.NAtoms << endl;
	    fid << "# delta = " << m.ss->delta << endl;
	    for (int i = 0; i < m.NAtoms; ++i)
	    {
		fid <<
		    m.ss->delta * m.h[i] << '\t' <<
		    m.ss->delta * m.k[i] << endl;
	    }
	    break;
	case m.ATOMEYE:
	    fid << "# BGA molecule format = atomeye" << endl;
	    fid << "# NAtoms = " << m.NAtoms << endl;
	    fid << "# delta = " << m.ss->delta << endl;
	    fid << "Number of particles = " << m.NAtoms << endl;
	    fid << "A = 1.0 Angstrom (basic length-scale)" << endl;
	    fid << "H0(1,1) = " << 2.0 * m.ss->dmax << " A" << endl;
	    fid << "H0(1,2) = 0 A" << endl;
	    fid << "H0(1,3) = 0 A" << endl;
	    fid << "H0(2,1) = 0 A" << endl;
	    fid << "H0(2,2) = " << 2.0 * m.ss->dmax << " A" << endl;
	    fid << "H0(2,3) = 0 A" << endl;
	    fid << "H0(3,1) = 0 A" << endl;
	    fid << "H0(3,2) = 0 A" << endl;
	    fid << "H0(3,3) = " << 2.0 * m.ss->dmax << " A" << endl;
	    fid << ".NO_VELOCITY." << endl;
	    // 4 entries: x, y, z, Uiso
	    fid << "entry_count = 4" << endl;
	    fid << "auxiliary[0] = abad [au]" << endl;
	    fid << endl;
	    // pj: now it only works for a single Carbon atom in the molecule
	    fid << "12.0111" << endl;
	    fid << "C" << endl;
	    for (int i = 0; i < m.NAtoms; i++)
	    {
		fid <<
		    (m.h[i]*m.ss->delta + m.ss->dmax)/(2*m.ss->dmax) << " " <<
		    (m.k[i]*m.ss->delta + m.ss->dmax)/(2*m.ss->dmax) << " " <<
		    0.5 << " " <<
		    m.abad[i] << endl;
	    }
	    break;
    }
    return fid;
}

void Molecule::PrintBadness()
{
    // call to MBadness() will update abadMax if necessary
    cout << "MBadness() = " << MBadness() << endl;
    cout << "ABadness() =";
    double mx = (NAtoms > 0) ? abad.max() : 0.0;
    for (int i = 0; i < NAtoms; ++i)
    {
	cout << ' ' << ABadness(i);
	if (ABadness(i) == mx)
	{
	    cout << '*';
	    mx += 1.0;
	}
    }
    cout << endl;
}

void Molecule::PrintFitness()
{
    cout << "MFitness() = " << MFitness() << endl;
    cout << "AFitness() =";
    double mx = (NAtoms > 0) ? abad.max() : 0.0;
    for (int i = 0; i < NAtoms; ++i)
    {
	cout << ' ' << AFitness(i);
	if (ABadness(i) == mx)
	{
	    cout << '*';
	    mx += 1.0;
	}
    }
    cout << endl;
}

/***********************************************************************
* Short Title: declaration of Lattice class
*
* Comments: class for general coordinate system
*
* $Id$
*
* <license text>
***********************************************************************/

#ifndef LATTICE_HPP_INCLUDED
#define LATTICE_HPP_INCLUDED

#include "R3linalg.hpp"

class Lattice
{
    public:

        // constructors

        Lattice();
        Lattice(double a0, double b0, double c0,
                double alpha0, double beta0, double gamma0);
        // create from base vectors
        template <class V>
            Lattice(const V& va0, const V& vb0, const V& vc0);

        // methods

        // set lattice parameters
        void setLatPar(
                double a0, double b0, double c0,
                double alpha0, double beta0, double gamma0 );
        // set lattice base vectors
        void setLatBase(
                const R3::Vector& va0,
                const R3::Vector& vb0,
                const R3::Vector& vc0 );
        template <class V>
            void setLatBase(const V& va0, const V& vb0, const V& vc0);
        // direct lattice
        inline double a() const         { return _a; }
        inline double b() const         { return _b; }
        inline double c() const         { return _c; }
        inline double alpha() const     { return _alpha; }
        inline double beta() const      { return _beta; }
        inline double gamma() const     { return _gamma; }
        inline double cosalpha() const  { return cosa; }
        inline double cosbeta() const   { return cosb; }
        inline double cosgamma() const  { return cosg; }
        inline double sinalpha() const  { return sina; }
        inline double sinbeta() const   { return sinb; }
        inline double singamma() const  { return sing; }
        // reciprocal lattice
        inline double ar() const        { return _ar; }
        inline double br() const        { return _br; }
        inline double cr() const        { return _cr; }
        inline double alphar() const    { return _alphar; }
        inline double betar() const     { return _betar; }
        inline double gammar() const    { return _gammar; }
        inline double cosalphar() const { return cosar; }
        inline double cosbetar() const  { return cosbr; }
        inline double cosgammar() const { return cosgr; }
        inline double sinalphar() const { return sinar; }
        inline double sinbetar() const  { return sinbr; }
        inline double singammar() const { return singr; }
        // vector operations using lattice coordinates
        template <class V>
            inline double dot(const V& u, const V& v) const;
        template <class V>
            inline double norm(const V& u) const;
        template <class V>
            inline double dist( const V& u, const V& v) const;
        // angle in degrees
        template <class V>
            inline double angledeg(const V& u, const V& v) const;
        // angle in radians
        template <class V>
            inline double anglerad(const V& u, const V& v) const;
        // conversion of coordinates and tensors
        const R3::Vector& cartesian(const R3::Vector& vl) const;
        const R3::Vector& fractional(const R3::Vector& vc) const;
        const R3::Matrix& cartesianMatrix(const R3::Matrix& Ml) const;
        const R3::Matrix& fractionalMatrix(const R3::Matrix& Mc) const;
        // lattice related tensors
        // metrics tensor
        const R3::Matrix& metrics() const   { return _metrics; }
        // matrix of base vectors, base() = stdbase() * baserot()
        const R3::Matrix& base() const      { return _base; }
        // standard base vectors
        const R3::Matrix& stdbase() const   { return _stdbase; }
        // base rotation matrix
        const R3::Matrix& baseRot() const   { return _baserot; }
        // inverse of base matrix
        const R3::Matrix& recbase() const   { return _recbase; }

    private:

        // data - direct lattice parameters
        double _a, _b, _c;
        double _alpha, _beta, _gamma;
        double cosa, cosb, cosg;
        double sina, sinb, sing;

        // data - reciprocal lattice parameters
        double _ar, _br, _cr;
        double _alphar, _betar, _gammar;
        double cosar, cosbr, cosgr;
        double sinar, sinbr, singr;

        // data - tensors
        // _base = _stdbase * _baserot
        R3::Matrix _metrics;        // metrics tensor
        R3::Matrix _base;           // lattice base
        R3::Matrix _stdbase;        // standard unit cell base
        R3::Matrix _baserot;        // base rotation matrix
        R3::Matrix _recbase;        // inverse of base matrix
        // base multiplied by magnitudes of reciprocal vectors
        R3::Matrix _normbase;
        R3::Matrix _recnormbase;    // inverse of _normbase


};

////////////////////////////////////////////////////////////////////////
// inline and template definitions
////////////////////////////////////////////////////////////////////////

template <class V>
Lattice::Lattice(const V& va0, const V& vb0, const V& vc0)
{
    setLatBase(va0, vb0, vc0);
}

template <class V>
void Lattice::setLatBase(const V& va0, const V& vb0, const V& vc0)
{
    R3::Vector va(va0[0], va0[1], va0[2]);
    R3::Vector vb(vb0[0], vb0[1], vb0[2]);
    R3::Vector vc(vc0[0], vc0[1], vc0[2]);
    setLatBase(va, vb, vc);
}

template <class V>
inline double Lattice::dot(const V& u, const V& v) const
{
    double dp = u[0]*v[0]*_metrics(0,0) +
        u[1]*v[1]*_metrics(1,1) +
        u[2]*v[2]*_metrics(2,2) +
        (u[0]*v[1] + u[1]*v[0])*_metrics(0,1) +
        (u[0]*v[2] + u[2]*v[0])*_metrics(0,2) +
        (u[2]*v[1] + u[1]*v[2])*_metrics(1,2);
    return dp;
}

template <class V>
inline double Lattice::norm(const V& u) const
{
    return sqrt(dot(u, u));
}

template <class V>
inline double Lattice::dist(const V& u, const V& v) const
{
    static R3::Vector duv;
    duv[0] = u[0] - v[0];
    duv[1] = u[1] - v[1];
    duv[2] = u[2] - v[2];
    return norm(duv);
}

template <class V>
inline double Lattice::angledeg(const V& u, const V& v) const
{
    return 180.0/M_PI * anglerad(u, v);
}

template <class V>
inline double Lattice::anglerad(const V& u, const V& v) const
{
    double ca = dot(u, v)/(norm(u) * norm(v));
    return acos(ca);
}

#endif  // LATTICE_HPP_INCLUDED

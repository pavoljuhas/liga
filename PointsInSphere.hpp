#include "dbprint.h"
/***********************************************************************
*
* pdffit2           by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2007 trustees of the Michigan State University
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
************************************************************************
*
* classes PointsInSphere, ReflectionsInQminQmax, ReflectionsInDmaxDmin
*
* Constructors:
*
*     PointsInSphere(Rmin, Rmax, a, b, c, alpha, beta, gamma)
*     ReflectionsInQminQmax(Qmin, Qmax, a, b, c, alpha, beta, gamma)
*     ReflectionsInDmaxDmin(Dmax, Dmin, a, b, c, alpha, beta, gamma)
*
*     template<class Lattice> PointsInSphere(Rmin, Rmax, const Lattice&)
*
*     where class Lattice must provide methods a(), b(), c(),
*     alpha(), beta(), gamma()
*
* Examples:
*
*     PointsInSphere sph(Rmin, Rmax, a, b, c, alpha, beta, gamma)
*     for (sph.rewind(); !sph.finished(); sph.next())
*     { 
*         // lattice indices are in sph.m, sph.n, sph.o or sph.mno[3]
*         // sph.r() is distance from origin,
*         // where Rmin < sph.r() < Rmax
*     }
*
*     ReflectionsInQminQmax refl(Qmin, Qmax, a, b, c, alpha, beta, gamma)
*     for (ReflectionsInQminQmax ref(Qmin, Qmax, a, b, c, alpha, beta, gamma);
*	   !ref.finished(); ref.next() )
*     { 
*         // Miller indices are in ref.h, ref.k, ref.l or ref.hkl[3]
*         // ref.Q() is magnitude of Q vector
*         // ref.d() is lattice plane spacing
*     }
*
* Tip: add epsilon to Rmax to avoid roundoff issues
*
* $Id$
*
***********************************************************************/

#ifndef POINTSINSPHERE_HPP_INCLUDED
#define POINTSINSPHERE_HPP_INCLUDED

// ensure math constants get defined for MSVC 
#define _USE_MATH_DEFINES
#include <cmath>


namespace NS_POINTSINSPHERE {

class LatticeParameters
{
    public:

        // data
        // input arguments
        double a, b, c, alpha, beta, gamma;
        // cosines and sines of direct lattice angles
        double ca, cb, cg, sa, sb, sg;
        // reciprocal lattice and its cosines and sines
        double ar, br, cr, alphar, betar, gammar;
        double car, cbr, cgr, sar, sbr, sgr;

        // constructor
        LatticeParameters(double _a, double _b, double _c,
                double _alpha, double _beta, double _gamma);

        // methods
        // calculate all properties from current lattice parameters
        void update();
        // return a reciprocal of this lattice
        LatticeParameters reciprocal() const;
        
    private:

        // methods
        static double _cosd(double x);
        static double _sind(double x);
};

}	// namespace NS_POINTSINSPHERE


class PointsInSphere
{
    public:

        // data
        // input arguments
        const double Rmin, Rmax;
        // results
        // mno array and m, n, o aliases are supposed to be read only
        int mno[3];
        int &m, &n, &o;

        // constructors
        PointsInSphere(double _Rmin, double _Rmax,
                const NS_POINTSINSPHERE::LatticeParameters& _latpar);
        PointsInSphere(double _Rmin, double _Rmax,
                double _a, double _b, double _c,
                double _alpha, double _beta, double _gamma);
        template <class L>
        PointsInSphere(double _Rmin, double _Rmax, const L&);

        // methods
        void rewind();
        inline void next()              { next_o(); }
        inline bool finished() const    { return !(m < hi_m); }
        double r() const;

    private:

        // data
        const NS_POINTSINSPHERE::LatticeParameters latpar;
        // calculated constants set by init()
        double RminSquare, RmaxSquare;
        // 2D reciprocal parameters and cosine in bc plane
        double b2r, c2r, ca2r;
        // reciprocal c
        double c1r;
        // offset of the nearest point to [0,0,0]
        double dn0dm, do0dm, do0dn;
        // loop variables
        double n0plane, o0plane, o0line;
        double mHalfSpan, nHalfSpan, oHalfSpan;
        // o indices excluded due to Rmin
        double oExclHalfSpan;
        int hi_m, hi_n, hi_o, outside_o;
        double RplaneSquare;

        // methods
        // loop advance
        void next_m();
        void next_n();
        void next_o();
        void init();
};

// template constructor

template <class L>
PointsInSphere::PointsInSphere(double _Rmin, double _Rmax, const L& lat) :
    Rmin(_Rmin), Rmax(_Rmax),
    m(mno[0]), n(mno[1]), o(mno[2]),
    latpar(lat.a(), lat.b(), lat.c(), lat.alpha(), lat.beta(), lat.gamma())
{
    init();
    rewind();
}


class ReflectionsInQminQmax
{
    private:

        // data - sph must be initialized before hkl and h, k, l
        const NS_POINTSINSPHERE::LatticeParameters latpar;
        PointsInSphere sph;

    public:

        // data
        // input arguments
        const double Qmin, Qmax;
        // results - hkl array and h, k, l aliases are read only
        int *hkl;
        int &h, &k, &l;

        // constructors
        ReflectionsInQminQmax(double _Qmin, double _Qmax,
                const NS_POINTSINSPHERE::LatticeParameters& _latpar);
        ReflectionsInQminQmax(double _Qmin, double _Qmax,
                double _a, double _b, double _c,
                double _alpha, double _beta, double _gamma);
        template <class L>
        ReflectionsInQminQmax(double _Qmin, double _Qmax, const L&);

        // methods
        inline void rewind()            { sph.rewind(); }
        inline void next()              { sph.next(); }
        inline bool finished() const    { return sph.finished(); }
        inline double Q() const         { return 2.0*M_PI*sph.r(); }
        inline double d() const         { return 1.0/sph.r(); }
};

// template constructor
template <class L>
ReflectionsInQminQmax::ReflectionsInQminQmax(
        double _Qmin, double _Qmax, const L& lat) :
	    latpar(lat.a(), lat.b(), lat.c(),
                    lat.alpha(), lat.beta(), lat.gamma()),
	    sph(Qmin*M_1_PI/2.0, Qmax*M_1_PI/2.0, latpar.reciprocal()),
	    Qmin(_Qmin), Qmax(_Qmax),
	    hkl(sph.mno), h(hkl[0]), k(hkl[1]), l(hkl[2])
{ }


class ReflectionsInDmaxDmin : public ReflectionsInQminQmax
{
    public:

        // data
        // input arguments
        const double Dmax, Dmin;

        // constructors
        ReflectionsInDmaxDmin(double _Dmax, double _Dmin,
                const NS_POINTSINSPHERE::LatticeParameters& _latpar);
        ReflectionsInDmaxDmin(double _Dmax, double _Dmin,
                double _a, double _b, double _c,
                double _alpha, double _beta, double _gamma);
        template <class L>
        ReflectionsInDmaxDmin(double _Dmax, double _Dmin, const L&);
};

// template constructor
template <class L>
ReflectionsInDmaxDmin::ReflectionsInDmaxDmin(
        double _Dmax, double _Dmin, const L& lat) :
	    ReflectionsInQminQmax(2.0*M_PI/_Dmax, 2.0*M_PI/_Dmin, lat),
	    Dmax(_Dmax), Dmin(_Dmin)
{ }

#endif	// POINTSINSPHERE_HPP_INCLUDED

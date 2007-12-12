/***********************************************************************
* Short Title: declarations of AtomFilter_t and derived classes
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef ATOMFILTER_T_HPP_INCLUDED
#define ATOMFILTER_T_HPP_INCLUDED

class Atom_t;
class Molecule;

////////////////////////////////////////////////////////////////////////
// class AtomFilter_t
////////////////////////////////////////////////////////////////////////

class AtomFilter_t
{
    public:

        // destructor
        virtual ~AtomFilter_t() { }

        // methods
        virtual bool Check(Atom_t*, Molecule* pm=NULL)  { return true; }

};

class BondAngleFilter_t : public AtomFilter_t
{
    public:

        // constructor and destructor
        BondAngleFilter_t(double _max_blen);
        virtual ~BondAngleFilter_t() { }

        // methods
        bool Check(Atom_t*, Molecule* pm);
        void setBondAngleRange(double lo, double hi);

    private:

        // data
        double max_blen;
        double lo_bangle;
        double hi_bangle;
};

class LoneAtomFilter_t : public AtomFilter_t
{
    public:

        // constructor and destructor
        LoneAtomFilter_t(double _max_dist);
        virtual ~LoneAtomFilter_t() { }

        // methods
        bool Check(Atom_t*, Molecule* pm);

    private:

        // data
        double max_dist;
};

#endif	// ATOMFILTER_T_HPP_INCLUDED

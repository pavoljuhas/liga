/***********************************************************************
* Short Title: useful functions for string manipulation
*
* Comments:
*
* $Id$
* 
* <license text>
***********************************************************************/

#ifndef STRINGUTILS_HPP_INCLUDED
#define STRINGUTILS_HPP_INCLUDED

#include <string>

template<typename StrType, typename Sequence>
std::string join(StrType sep, Sequence words)
{
    std::string joined;
    typename Sequence::iterator w = words.begin();
    typename Sequence::iterator last = words.end();
    if (w == last)	return joined;
    for (joined = *(w++); w != last; ++w)
    {
	joined += sep;
	joined += *w;
    }
    return joined;
}

#endif	// STRINGUTILS_HPP_INCLUDED

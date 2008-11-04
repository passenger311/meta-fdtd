#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <stdio.h>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>

using namespace std;

struct strCmp {
	bool operator()( const char* s1, const char* s2 ) const {
		return strcmp( s1, s2 ) < 0;
	}
};

struct stringCmp {
	bool operator()( string s1, string s2 ) const {
		return s1 < s2;
	}
};

typedef map<const char*, const char*, strCmp> dictmap;

#endif /*GLOBAL_H_*/

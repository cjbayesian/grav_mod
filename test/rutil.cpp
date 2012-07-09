#define MATHLIB_STANDALONE 1
#include <Rmath.h>

#include<iostream>
#include<cstdlib>
//typedef enum {
//    BUGGY_KINDERMAN_RAMAGE,
//    AHRENS_DIETER,
//    BOX_MULLER,
//    USER_NORM,
//    INVERSION,
//    KINDERMAN_RAMAGE
//} N01type;


using namespace std;

int main()
{
    qnorm(0.7, 0.0, 1.0, 0, 0);	
	for(int i=1;i<=10;i++)
		cout << rnorm(0,1) <<"\n";
	return 0;
}


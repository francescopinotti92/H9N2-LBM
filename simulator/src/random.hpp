//
//  random.hpp
//  SimulateSingleMarket
//
//  Created by user on 28/08/2021.
//

#ifndef random_hpp
#define random_hpp

#include <cstdint>
#include <random>
#include <math.h>
#include <stdio.h>

extern std::mt19937_64 m_mt;
extern std::uniform_real_distribution<double> rand_uniform;

bool getBool(const double& p);
double getUni();
double getUniPos();
double getGamma(const double& a, const double& b); // a is shape, b is scale (not rate)
double getBeta(const double& a, const double& b);
int getUniInt(const int& max);
int getBinom(double p, int n);
int getGeom1( const double& p );



double gamma_int (const int& a);
double gamma_large (const double& a);
double gamma_frac (const double& a);



#endif /* random_hpp */

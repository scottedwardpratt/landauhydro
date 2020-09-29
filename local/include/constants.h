#ifndef __CORAL_CONSTANTS_H__
#define __CORAL_CONSTANTS_H__

#include <cmath>
#include <complex>

const double ZERO         = 0.0;
const double HBARC        = 197.3269718;       // hbar times c
const double ALPHA        = 1.0/137.03599976;   // fine structure constant
const double PI           = 3.1415926535897932384626433832795028841972;
const double SQRTPI       = 1.7724538509055160272981674833411451827975;
const double SQRTFOURPI   = 3.5449077018110320545963349666822903655951;
const double DEGRAD       = 57.2957795130823208767981548141051703324055;
const double AMU          = 931.494;          // atomic mass unit
const double ProtonMass   = 938.272;
const double KaonMass     = 493.677;
const double PionMass     = 139.57018;
const double Pion0Mass    = 134.9766;
const double LambdaMass   = 1115.7;
const double NeutronMass  = 939.565;
const double RhoMass      = 771.1;
const double XiMass       = 1321.3;
const double XiStarMass   = 1530.0;
const std::complex< double > ci = std::complex< double >(0.0,1.0);

#endif

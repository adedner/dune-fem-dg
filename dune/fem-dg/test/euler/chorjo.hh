#ifndef EULERCHORIN_HH
#define EULERCHORIN_HH
/**************
 *   euler.c   *
 **************/


/*********************************************************************                      *
*     Bestimmung der 1-D Riemannloesung der Eulergleichungen         *
*     nach dem Verfahren von Chorin (liefert die "exakte" Loesung).  *
*     Literatur : Random Choice Solution of Hyperbolic Systems,      *
*                 A.J. Chorin, J. of Comp. Phys. 22, 517-533 (1976)  *
*********************************************************************/

#include<cstdio>
#include<cmath>

namespace EULERCHORIN {

double hilf(double,double,double);
double maximum(double,double);
void iteration(double,double,double,double,double,double,double *p,
               double *u,double);
void lsg(double,double,double *q_erg,double *u_erg,double *p_erg,
         double,double,double,double,double,double,double);
}
#endif

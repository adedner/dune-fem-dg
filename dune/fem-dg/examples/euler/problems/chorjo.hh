#ifndef DUNE_FEM_DG_EULERCHORIM_HH
#define DUNE_FEM_DG_EULERCHORIM_HH

#include <cmath>
#include <iostream>

namespace EULERCHORIN {

/*********************************************************************
*    Bestimmung von phi (siehe Artikel)                              *
*********************************************************************/
inline double hilf(double h1,double h2,double ny)
{
   double erg,h,w,z;

   w=h1/h2;
   if (w < 1)
     { z=(ny-1.0)/(2.0*ny);
       h=pow(w,z);
       erg=(1.0-w)/(1.0-h)*(ny-1.0)/(2.0*std::sqrt(ny));
     }
   else
       erg=std::sqrt((ny+1.0)/2.0*w+(ny-1.0)/2.0);
   return(erg);
}



/*********************************************************************
*   liefert Maximum zweier reeller Zahlen                            *
*********************************************************************/
inline double maximum(double z1,double z2)
{
   double max;

   if (z1 > z2)
     max=z1;
   else
     max=z2;
   return(max);
}




/*********************************************************************
*    iterative Bestimmung von p*, Bestimmung von u*,                 *
*    siehe Artikel                                                   *
*********************************************************************/
inline void iteration(double ql,double qr,double ul,double ur,double pl,\
               double pr,double *p,double *u,double ny)
{
   double  ml,mr,ml_v,mr_v;
   double  erg,p1,p2,eps,a;
   int     z,l,abbruch1,abbruch2;

   abbruch1=20;
   abbruch2=10;
   l=0;
   z=0;
   a=0.5;
   eps=10E-06;
   ml_v=10E20;
   mr_v=10E20;

   *p=(pl+pr)/2.0;
marke:
   erg=hilf(*p,pl,ny);
   ml=-std::sqrt(pl*ql)*erg;
   erg=hilf(*p,pr,ny);
   mr=std::sqrt(pr*qr)*erg;
   p1=(ul-ur+pr/mr-pl/ml)/(1.0/mr-1.0/ml);
   p2=maximum(p1,eps);
   *p=a*p2+(1.0-a)*(*p);
   erg=maximum(std::abs(mr-mr_v),std::abs(ml-ml_v));
   if (erg >= eps) {
       ml_v=ml;
       mr_v=mr;
       l=l+1;
       if (l <= abbruch1)
          goto marke;
       else {
         z=z+1;
         if (z <= abbruch2) {
             a=a/2.0;
             l=0;
             goto marke;
             }
         else ;
              /*printf(" Ueberpruefen der Konvergenz !! "); */
         }
       }
   *u=(*p-pl)/ml+ul;
}





/*********************************************************************
*    Dieses Unterprogramm bestimmt im Punkte(x,t) die Riemann-       *
*    loesung der Eulergleichungen.                                   *
*    Benoetigte Hilfsproceduren :                                    *
*          iteration,(maximum),(hilf),(hoch)                         *
*********************************************************************/
inline void lsg(double x,double t,double *q_erg,double *u_erg,double *p_erg,\
         double ql,double qr,double ul,double ur,double pl,double pr,\
         double ny)
{
   double  q,s,h1,h2,h3,p,u;
   double  links,rechts,cl,cr,c;

/*  Bestimmung von p* und u*  */
   iteration(ql,qr,ul,ur,pl,pr,&p,&u,ny);

/*  Punkt liegt links der Kontaktunstetigkeit  */
   if(x <= u*t) {
/*  linke Welle ist ein Schock  */
       if (pl < p) {
           h1=p/pl;
           h2=(ny+1.0)/(ny-1.0);
           q=(1.0+h1*h2)/(h1+h2)*ql;
           s=(ql*ul-q*u)/(ql-q);
/*  Punkt liegt links vom Schock  */
           if (x < s*t) {
               *q_erg=ql;
               *u_erg=ul;
               *p_erg=pl;
               }
           else {
/*  Punkt liegt rechts vom Schock  */
                *q_erg=q;
                *u_erg=u;
                *p_erg=p;
                }
/*  linke Welle ist eine Verduennungswelle  */
           }
       else {
           cl=std::sqrt(ny*pl/ql);
           c=cl+0.5*(ny-1.0)*(ul-u);
           q=ny*p/(c*c);
           rechts=u-c;
           links=ul-cl;
/*  Punkt liegt links der Verduennungswelle  */
           if (x < links*t) {
                  *q_erg=ql;
                  *u_erg=ul;
                  *p_erg=pl;
                  }
/*  Punkt liegt rechts der Verduennungswelle  */
           else {
                if (x > rechts*t) {
                     *q_erg=q;
                     *u_erg=u;
                     *p_erg=p;
                     }
/*  Punkt liegt in der Verduennungswelle  */
              else {
                  cl=std::sqrt(ny*pl/ql);
                  c=std::sqrt(ny*p/q);
                  *u_erg=(x/t+c+0.5*(ny-1.0)*u)/(1.0+0.5*(ny-1.0));
                  c=(*u_erg)-x/t;
                  h1=pow(q,ny);
                  h2=c*c/ny*h1/p;
                  h3=1.0/(ny-1.0);
                  *q_erg=pow(h2,h3);
                  *p_erg=c*c*(*q_erg)/ny;
                  }
           }
       }
    }
/*  Punkt liegt rechts der Kontaktunstetigkeit  */
   else {
/* rechte Welle ist ein Schock  */
        if (p > pr) {
             h1=pr/p;
             h2=(ny+1.0)/(ny-1.0);
             q=(h1+h2)/(1.0+h1*h2)*qr;
             s=(qr*ur-q*u)/(qr-q);
/*  Punkt liegt links vom Schock  */
             if (x <= s*t) {
                  *q_erg=q;
                  *u_erg=u;
                  *p_erg=p;
                  }
/*  Punkt liegt rechts vom Schock  */
             else {
                 *q_erg=qr;
                 *u_erg=ur;
                 *p_erg=pr;
                 }
          }
/*  rechte Welle ist eine Verduennungswelle  */
        else {
            cr=std::sqrt(ny*pr/qr);
            c=cr+0.5*(ny-1.0)*(u-ur);
            q=ny*p/(c*c);
            links=u+c;
            rechts=ur+cr;
/*  Punkt liegt links der Verduennungswelle  */
            if (x <= links*t) {
                 *q_erg=q;
                 *u_erg=u;
                 *p_erg=p;
                 }
            else {
/*  Punkt liegt rechts der Verduennungswelle  */
                 if (x >= rechts*t) {
                     *q_erg=qr;
                     *u_erg=ur;
                     *p_erg=pr;
                     }
/*  Punkt liegt in der Verduennungswelle  */
                else {
                   cr=std::sqrt(ny*pr/qr);
                   *u_erg=(x/t-cr+ur/2.0*(ny-1.0))/(1.0+(ny-1.0)/2.0);
                   c=x/t-(*u_erg);
                   h1=pow(qr,ny);
                   h2=c*c/ny*h1/pr;
                   h3=1.0/(ny-1.0);
                   *q_erg=pow(h2,h3);
                   *p_erg=c*c*(*q_erg)/ny;
                   }
              }
          }
     }
}

} // end namespace EULERCHORIN

#endif

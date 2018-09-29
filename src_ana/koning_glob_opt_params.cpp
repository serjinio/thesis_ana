#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv){

  double z, a, e;
  double n;
  int zint, aint, eint;

  double v1, v2, v3, v4, w1, w2, d1, d2, d3;
  double vso1, vso2, wso1, wso2, efp, vbarc;

  double rc;               // Coulomb
  double vr, rr, ar;       // Volume real
  double wv, rwv, awv;     // Volume imaginary
  double ws, rws, aws;     // Surface imarginary
  double vls, rls, als;    // LS real
  double wls, rwls, awls;  // LS imaginary

  char filename[100];
  FILE *fout;

  z = atof(argv[1]);
  a = atof(argv[2]);
  e = atof(argv[3]);
  n = a-z;

  zint = (int)z;
  aint = (int)a;
  eint = (int)e;


  sprintf(filename,"koning_z%da%de%d.pot",zint,aint,eint);

   if((fout=fopen(filename, "w")) == NULL)
   {
     printf("Cannot open the output file. \n");
     exit(1);
   }


  rc = 1.198+0.697*pow(a,-0.66667)+12.994*pow(a,-1.66667);


  v1 = 59.3+21.0*(n-z)/a-0.024*a;
  v2 = 0.007067+4.23*pow(10,-6)*a;
  v3 = 1.729*pow(10,-5)+1.136*pow(10,-8)*a;
  v4 = 0.000000007;
  w1 = 14.667+0.009629*a;
  w2 = 73.55+0.0795*a;
  d1 = 16.0+16.0*(n-z)/a;
  d2 = 0.018+0.003802/(1+exp((a-156.0)/8.0));
  d3 = 11.5;
  vso1 = 5.922+0.003*a;
  vso2 = 0.004;
  wso1 = -3.1;
  wso2 = 160.0;
  efp = -8.4075+0.01378*a;
  vbarc = 1.73/rc*z*pow(a,-0.33333);



  vr = v1*(1-v2*(e-efp)+v3*pow((e-efp),2)-v4*pow((e-efp),3))\
    +vbarc*v1*(v2-2*v3*(e-efp)+3*v4*pow((e-efp),2));

  rr = 1.3039-0.4054*pow(a,-0.33333);
  ar = 0.6778-1.487*pow(10,-4)*a;
  wv = w1*pow((e-efp),2)/(pow((e-efp),2)+pow(w2,2));
  rwv = rr;
  awv = ar;
  ws = d1*pow((e-efp),2)/(pow((e-efp),2)+d3*d3)*exp(-d2*(e-efp));
  rws = 1.3424-0.01585*pow(a,0.33333);
  aws = 0.5187+5.205*pow(10,-4)*a;
  vls = vso1*exp(-vso2*(e-efp));
  rls = 1.1854-0.647*pow(a,-0.33333);
  als = 0.59;
  wls = wso1*pow((e-efp),2)/(pow((e-efp),2)+pow(wso2,2));
  rwls = 1.1854-0.647*pow(a,-0.33333);
  awls = 0.59;


  printf("Z= %3.0f     N= %3.0f     A= %3.0f     E= %5.1f    r_c= %8.5f\n", z, n, a, e, rc);



  printf("%10.5f%10.5f%10.5f\n", vr, rr, ar);
  printf("%10.5f%10.5f%10.5f\n", wv, rwv, awv);
  printf("%10.5f%10.5f%10.5f\n", ws, rws, aws);
  printf("%10.5f%10.5f%10.5f\n", vls, rls, als);
  printf("%10.5f%10.5f%10.5f\n\n", wls, rwls, awls);

  fprintf(fout,"%10.5f%10.5f%10.5f\n", z, a, e);
  fprintf(fout,"%10.5f%10.5f%10.5f\n", vr, rr, ar);
  fprintf(fout,"%10.5f%10.5f%10.5f\n", wv, rwv, awv);
  fprintf(fout,"%10.5f%10.5f%10.5f\n", ws, rws, aws);
  fprintf(fout,"%10.5f%10.5f%10.5f\n", vls, rls, als);
  fprintf(fout,"%10.5f%10.5f%10.5f\n", wls, rwls, awls);


  exit(0);
}

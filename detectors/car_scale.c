/************************************************************
*
* Version 08/03/2001
************************************************************/
#include <stdio.h>
#include <math.h>
#define PI 3.14159

main(int argc, char **argv)
{
int day, month, year;
double xx1, yy1, xx2, yy2;
double rho1, theta1, rho2, theta2, rho45, rho10, rho50, rho75;
  printf("car_scale/Version 08-03-2001, from polar OXY to sky coordinates\n");
if(argc == 7 && (*argv[4] && *argv[4] != ' ')) argc = 4;
if(argc == 7 && (*argv[3] && *argv[3] != ' ')) argc = 3;
if(argc == 7 && (*argv[2] && *argv[2] != ' ')) argc = 2;
if(argc != 3)
  {
  printf("Error, syntax is: car_scale rho theta date (coordinates in pixels)\n");
  printf(" Example: car_scale 0.024 134.8 30-01-97 \n");
  exit(-1);
  }
sscanf(argv[1],"%lf",&rho1);
sscanf(argv[2],"%lf",&theta1);
sscanf(argv[3],"%d-%d-%d",&day,&month,&year);
printf(" OK: rho1=%.3f pixels, theta1=%.2f degrees ", rho1, theta1);
printf(" Date is %.2d-%.2d-%.2d\n", day, month, year);

/* To Cartesian coordinates: */
theta1 *= (double)(PI / 180.);
xx1 = rho1 * cos(theta1);
yy1 = rho1 * sin(theta1);

/* Distorsion correction: */
 xx2 = xx1 - tan(2.5*PI/180.) * yy1;
 yy2 = 1.01 * yy1 / cos(2.5*PI/180.);
printf("yy1=%f yy2=%f\n",(float)yy1, (float)yy2);

/* To Polar coordinates: */
theta2 = atan2(yy2,xx2);
theta2 *= (double)(180./PI);
rho2 = xx2 * xx2 + yy2 * yy2;
rho2 = sqrt(rho2);

printf("After distorsion correction: rho = %.3f pixels, theta = %.2f degrees\n",
        rho2,theta2);

switch(year)
{
/* September 1994 ************************************************/
  case 94:
    switch(month)
       {
       default:
            printf("Calibration for September 1994:");
            printf(" .00750\"/pix 4.5mm; (87.5+theta deg)\n");
            rho45 = rho2 * 0.00750;
            theta2 = 87.5 + theta2;
            printf("=> rho(4.5mm)=%.3f arcsec (ireduc=2x)", rho45);
            printf(" or %.3f arcsec (ireduc=4x)\n", (float)(2.*rho45));
       break;
       }
    break;
/* July 1995 ************************************************/
  case 95:
    switch(month)
       {
       default:
       switch(day)
          {
          case 20:
            printf("Calibration for 20 Jul 1995:");
            printf(" 0.00734\"/pix (93.7+theta deg)\n");
/* 0.00765 arcsec/pixel with 4.5 mm eyepiece:  (F=50.4m: 0.4093" per 0.1mm)*/
            rho45 = rho2 * 0.00734;
            theta2 = 93.7 + theta2;
            printf("=> rho(4.5mm)=%.3f arcsec (ireduc=2x)", rho45);
            printf(" or %.3f arcsec (ireduc=4x)\n", (float)(2.*rho45));
            break;
          case 21:
            printf("Calibration for 21 Jul 1995:");
            printf(" 0.00724\"/pix (93.7+theta deg)\n");
            rho45 = rho2 * 0.007245;
            theta2 = 93.7 + theta2;
            printf("=> rho(4.5mm)=%.3f arcsec (ireduc=2x)", rho45);
            printf(" or %.3f arcsec (ireduc=4x)\n", (float)(2.*rho45));
            break;
          case 22:
            printf("Calibration for 22 Jul 1995:");
            printf(" 0.00737\"/pix 4.5mm; 0.0168\"/pix 10mm (93.7+theta deg)\n");
            rho45 = rho2 * 0.00737;
            theta2 = 93.7 + theta2;
            printf("=> rho(4.5mm)=%.3f arcsec (ireduc=2x)", rho45);
            printf(" or %.3f arcsec (ireduc=4x)\n", (float)(2.*rho45));
            rho10 = rho2 * 2. * 0.0168;
            printf("=> rho(10mm)=%.3f arcsec (ireduc=4x)\n", rho10);
            break;
          case 23:
            printf("Calibration for 23 Jul 1995:");
            printf(" .00730\"/pix 4.5mm; .0171 oc10; .251 oc50 (93.7+theta deg)\n");
/* No calibration for 4.5, so I take the average between 21/07/95 and 22/07/95
*/
            rho45 = rho2 * 0.00730;
            theta2 = 93.7 + theta2;
            printf("=> rho(4.5mm)=%.3f arcsec (ireduc=2x)", rho45);
            printf(" or %.3f arcsec (ireduc=4x)\n", (float)(2.*rho45));
            rho10 = rho2 * 2. * 0.0171;
            printf("=> rho(10mm)=%.3f arcsec (ireduc=4x)\n", rho10);
            rho50 = rho2 * 2. *0.251;
            printf("=> rho(50mm)=%.3f arcsec (ireduc=4x)\n", rho50);
          case 24:
            printf("Calibration for 24 Jul 1995:");
            printf(" .00728\"/pix 4.5mm; .0170 oc10; (93.7+theta deg)\n");
            rho45 = rho2 * 0.00728;
            theta2 = 93.7 + theta2;
            printf("=> rho(4.5mm)=%.3f arcsec (ireduc=2x)", rho45);
            printf(" or %.3f arcsec (ireduc=4x)\n", (float)(2.*rho45));
            rho10 = rho2 * 2. * 0.0170;
            printf("=> rho(10mm)=%.3f arcsec (ireduc=4x)\n", rho10);
            break;
          default:
          break;
          }
       break;
       }
    break;
/* January 1997 ************************************************/
  case 97:
    printf("Calibration for 29-30 Jan 1997:");
    printf(" 0.0341\"/pix and 0.0247\"/pix (83.3deg-theta)\n");
/* 0.0341 arcsec/pixel with 10 mm eyepiece:  (F=50.4m)*/
    rho10 = rho2 * 0.0341;
/* 0.0247 arcsec/pixel with 7.5 mm eyepiece:  (F=50.4m)*/
    rho75 = rho2 * 0.0247;
    printf("=> rho(10mm)=%.4f arcsec, rho(7.5mm)=%.4f arcsec (ireduc=4x)\n",
           rho10, rho75);
/*    theta2 = 83.3 - theta2;
* Error detected on July 23rd 2001 !!!!
*/
    theta2 = 96.7 + theta2;
    break;
/* Default: ************************************************/
  default:
    rho10 = 0.;
    rho75 = 0.;
    theta2 = 0.;
    break;
}

if(theta2 < 0.) theta2 += 180.;

printf("=> theta=%.2f degrees, (NE orientation)\n", theta2);
exit(0);
}

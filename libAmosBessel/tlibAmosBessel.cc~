/* 
 * Copyright (C) 2007, 2008, 2009, 2010, 2011 M. T. Homer Reid
 * This file is part of libAmosBessel. 
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * tlibAmosBessel.cc -- test/demonstration program for libAmosBessel
 *
 * homer reid        -- 2/2011
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include <libhrutil.h>

#include "libAmosBessel.h"

#define CSIGN(z) (__imag__ (z) > 0.0 ? '+' : '-') 

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  double zr, zi, MinOrder;
  cdouble z;
  int NumOrders, Scale;
  char *WFString, WhichFunction;

  /*--------------------------------------------------------------*/
  /*- process command-line arguments -----------------------------*/
  /*--------------------------------------------------------------*/
  zr=zi=MinOrder=INFINITY;
  ArgStruct ASArray[]=
   { {"zr",            PA_DOUBLE, (void *)&zr,            0, "real part of argument"},
     {"zi",            PA_DOUBLE, (void *)&zi,            0, "imag part of argument"},
     {"MinOrder",      PA_DOUBLE, (void *)&MinOrder,  "0.0", "minimum order"},
     {"NumOrders",     PA_INT,    (void *)&NumOrders,   "3", "number of orders to compute"},
     {"WhichFunction", PA_STRING, (void *)&WFString,    "J", "which function to compute"},
     {"Scale",         PA_BOOL,   (void *)&Scale,       "0", "compute scaled functions"},
     {0,0,0,0,0}
   };
  ProcessArguments(argc, argv, ASArray);
  if ( isinf(zr) )
   ASUsage(argv[0],ASArray,"--zr option is mandatory");
  if ( isinf(zi) )
   ASUsage(argv[0],ASArray,"--zi option is mandatory");
  if ( !WFString )
   ASUsage(argv[0],ASArray,"--WhichFunction option is mandatory");
  WhichFunction=WFString[0];

  printf("at argument Z=%15.8e %c %15.8eI: \n",zr,CSIGN(1.0fi*zi),zi);
  z = zr + 1.0fi*zi;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (WhichFunction=='A' || WhichFunction=='a' || WhichFunction=='B' || WhichFunction=='b')
   { 
     cdouble f;
     AmosAiry(WhichFunction, z, Scale, &f);
     printf("%c(Z) = %+21.15e %c %20.14eI\n",WhichFunction, 
             __real__ f, CSIGN(f), fabs(__imag__ f));
   }
  else
   { cdouble f[NumOrders];
     double Order; 
     int no;

     AmosBessel(WhichFunction, z, MinOrder, NumOrders, Scale, f);

     if (Scale)
      { 
        _Complex long double ScaleFac;
        _Complex long double ff[NumOrders];

        switch(WhichFunction)
         { 
           case 'J': 
           case 'j': 
             ScaleFac = exp(fabs(zi)) + 0.0fi;
             break;

           case 'Y': 
           case 'y': 
             ScaleFac = exp(fabs(zi)) + 0.0fi;
             break;

           case 'I': 
           case 'i': 
             ScaleFac = exp(fabs(zr)) + 0.0fi;
             break;

           case 'K': 
           case 'k': 
             ScaleFac = cexp(-1.0*z);
             break;

           case 'O':
           case 'o':
             ScaleFac = cexp(1.0fi*z);
             break;

           case 'T': 
           case 't': 
             ScaleFac = cexp(-1.0*1.0fi*z);
             break;
         };

        for(no=0, Order=MinOrder; no<NumOrders; no++)
	 { ff[no]=ScaleFac * f[no];
           //f[no]=(cdouble)ff[no];
         };

      };

     for(no=0, Order=MinOrder; no<NumOrders; no++, Order+=1.0)
      printf("%c_{%g}(Z) = %+21.14e %c %20.14eI\n",WhichFunction,Order, 
           __real__ f[no], CSIGN(f[no]), fabs(__imag__ f[no]));
   };

  printf("\n");
  printf("Thank you for your support.\n");
}

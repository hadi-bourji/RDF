#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;
void function ( double * rx, double * ry, double * rz,  
                int N, double * box2, double rc2, double dr, int nB, int * H ) 
{

   int i,j;
   double dx,dy,dz,r2;
   int bin; 
   for (i=0;i<(N-1);i++) {
     for (j=i+1;j<N;j++) {
         dx = rx[i]-rx[j];
         dy = ry[i]-ry[j];
         dz = rz[i]-rz[j];
         dx = dx-(box2[0]*rint(dx/box2[0]));
         dy = dy-(box2[1]*rint(dy/box2[1]));
         dz = dz-(box2[2]*rint(dz/box2[2]));
         r2 = dx*dx + dy*dy + dz*dz;
         bin=(int)(((sqrt(r2)-(dr/2.))/dr)+1);
         H[bin] += 1;    
         
       }
     }
    }
 
extern "C" {
    void My_Function(double * rx, double * ry, double * rz,  
                     int N, double * box2, double rc2, double dr, int nB , int * H )
    {
        return  function(rx,ry,rz,N,box2,rc2,dr,nB,H);
    }
}

// How to use:
//g++ -fPIC -shared -o libTest.so rdf.cpp
//python calculations-CTypes.py

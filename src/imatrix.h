// Simple interval library from Luc JAULIN, with minor modifications from Fabrice LE BARS.

#ifndef __IMATRIX__
#define __IMATRIX__

#include <iostream>
#include <vector>
#include <math.h>
#include "interval.h"
#include "box.h"
#include "rmatrix.h"


#include "tnt.h"
using namespace TNT;

#include "jama_lu.h"
using namespace JAMA;


class imatrix
{
public:
	Array2D <interval> data;
//---------------------CONSTRUCTEURS--------------------------------------------
imatrix ();
imatrix (int,int);
imatrix (const rmatrix&);
imatrix (const box&);
imatrix (const imatrix&);
~imatrix ();
imatrix& operator=(const imatrix&);
//interval& imatrix::operator() (int i,int j)  const;
interval operator() (int i,int j) const;
interval GetVal (int,int) const;
void SetVal (int,int,interval);

int dim1 (void) const;
int dim2 (void) const;

//----------------------OPERATEURS----------------------------------------------
friend imatrix	operator+ (const imatrix&,const imatrix&);
friend imatrix	operator- (const imatrix&);
friend imatrix	operator- (const imatrix&,const imatrix&);
friend imatrix	operator* (const imatrix&,const imatrix&);
friend imatrix  operator* (const rmatrix& X,const imatrix& Y);
friend imatrix  operator*(const double a,const imatrix& X);
friend box operator*(const imatrix& A,const box& x);

friend ostream& operator<<(ostream&, const imatrix&);
//-------------------FONCTIONS NON MEMBRES---------------------------------------


friend rmatrix Center	(const imatrix&);
friend box ToBox(const imatrix&);
friend box Row   (const imatrix& B, int i);
friend box Column(const imatrix& B, int j);
friend imatrix Transpose(const imatrix& X);
friend imatrix iZeros (int n, int m);
friend imatrix iEye (int n);
friend imatrix RotationPhiThetaPsi(interval& phi, interval& theta, interval& psi);
friend void Cmult(imatrix& C, imatrix& A, imatrix& B);
friend void Cmult(box& c, imatrix& A, box& b);
friend void Crot(imatrix&);
friend void Cantisym(imatrix& A);

/*friend void         Update            (imatrix&);
friend imatrix          Rand              (const imatrix& X);
friend interval     Determinant       (imatrix&, imatrix&);
friend bool         Emptyimatrix          (const imatrix&);
friend bool	    Disjoint	      (const imatrix&,const imatrix&);
friend iboolean	    In		      (imatrix,imatrix);
friend imatrix	    Inf 	      (imatrix);
friend imatrix	    Inter 	      (const imatrix&,const imatrix&);
friend imatrix          Concat            (const imatrix&, const imatrix&);
friend interval     Norm              (imatrix);
friend interval     NormEuclid        (imatrix, imatrix);
friend interval     NormInf           (imatrix, imatrix);
friend bool	    Subset	      (imatrix&,imatrix&);
friend bool	    Subset	      (imatrix&,imatrix&,double);
friend imatrix	    Sup 	      (imatrix);
friend imatrix	    Union	      (imatrix&,imatrix&);
friend double	    Width	      (imatrix&);
friend double	    Width	      (imatrix&,vector<int>&);
friend double	    Width	      (imatrix&,imatrix&);
friend imatrix          Zeros             (int,int);
*/

};
#endif


// Simple interval library from Luc JAULIN, with minor modifications from Fabrice LE BARS.

#ifndef __RMATRIX__
#define __RMATRIX__

#include "interval.h"
#include "box.h"


#include "tnt.h"
using namespace TNT;

#include "jama_lu.h"
using namespace JAMA;


class rmatrix
{
public:
	Array2D <double> data;
//---------------------CONSTRUCTEURS--------------------------------------------
rmatrix ();
rmatrix (int,int);
rmatrix (const rmatrix&);
rmatrix (const box &);
~rmatrix ();
rmatrix& operator=(const rmatrix&);

double GetVal (int,int) const;
void SetVal (int,int,double);

int dim1 (void) const;
int dim2 (void) const;

//----------------------OPERATEURS----------------------------------------------
friend rmatrix	operator+(const rmatrix&,const rmatrix&);
friend rmatrix	operator-(const rmatrix&);
friend rmatrix	operator-(const rmatrix&,const rmatrix&);
friend rmatrix	operator*(const rmatrix&,const rmatrix&);
friend rmatrix	operator*(const double,const rmatrix&);
//friend rmatrix	operator*(const rmatrix&,const interval&);
//friend ostream& operator<<(ostream&, const rmatrix&);
//-------------------FONCTIONS MEMBRES---------------------------------------

friend rmatrix	Zeros (int, int);
friend rmatrix	Eye (int);
friend double	Det (rmatrix&);
friend rmatrix	Inv (rmatrix&);
friend rmatrix RotationPhiThetaPsi(double phi,double theta,double psi);



/*friend double       Angle             (rmatrix& ,rmatrix&); // Il faut des vecteurs de dim 2
friend int          Size              (const rmatrix&);
friend int	    AxePrincipal      (rmatrix&);
friend int	    AxePrincipal      (rmatrix&, vector<int>&);
friend int	    AxePrincipal      (rmatrix&, rmatrix&);
friend void         Update            (rmatrix&);
friend rmatrix          Rand              (const rmatrix& X);
friend rmatrix	    Center	      (const rmatrix&);
friend rmatrix          Center            (const rmatrix&, vector<int>&);
//friend void	    CheckRange	      (rmatrix&,rmatrix&);
friend interval     Determinant       (rmatrix&, rmatrix&);
friend bool         Emptyrmatrix          (const rmatrix&);
friend bool	    Disjoint	      (const rmatrix&,const rmatrix&);
friend double         decrease          (const rmatrix&, const rmatrix&);
friend double         decrease          (const rmatrix&, const rmatrix&, vector<int>);
friend double         Eloignement       (rmatrix&,rmatrix&);
friend double         Eloignement2      (rmatrix&,rmatrix&);
friend double         EloignementRelatif2(rmatrix&,rmatrix&);
friend double         Marge             (rmatrix,rmatrix);
friend iboolean	    In		      (rmatrix,rmatrix);
friend rmatrix	    Inf 	      (rmatrix);
friend rmatrix	    Inflate 	      (rmatrix&,double);
friend rmatrix	    Inter 	      (const rmatrix&,const rmatrix&);
friend rmatrix          Concat            (const rmatrix&, const rmatrix&);
friend rmatrix          Proj              (const rmatrix&, int, int);
//friend void       Inter1            (rmatrix&,rmatrix&,const rmatrix&,const rmatrix&,const rmatrix&);
friend interval     Norm              (rmatrix);
friend interval     NormEuclid        (rmatrix, rmatrix);
friend interval     NormInf           (rmatrix, rmatrix);
friend void         Phi0              (rmatrix&,rmatrix*);             //0-intersection (voir qminimax)
friend void         Phi1              (rmatrix&,rmatrix&,rmatrix*);       //1-intersection (voir qminimax)
friend void         Phi2              (rmatrix&,rmatrix&,rmatrix&,rmatrix*);  //2-intersection
friend interval     ProduitScalaire   (rmatrix&,rmatrix&);
friend bool         Prop              (rmatrix&,rmatrix&);
friend bool	    Subset	      (rmatrix&,rmatrix&);
friend bool	    Subset	      (rmatrix&,rmatrix&,double);
friend bool         Isrmatrix             (rmatrix);
friend void         Sucre             (rmatrix&,rmatrix&);
friend rmatrix	    Sup 	      (rmatrix);
friend rmatrix	    Union	      (rmatrix&,rmatrix&);
friend double         Volume            (rmatrix&);
friend double	    Width	      (rmatrix&);
friend double	    Width	      (rmatrix&,vector<int>&);
friend double	    Width	      (rmatrix&,rmatrix&);
//friend rmatrix          Empty             (int);
//friend rmatrix          Infinity          (int);   */
};
#endif

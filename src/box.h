// Simple interval library from Luc JAULIN, with minor modifications from Fabrice LE BARS.

#ifndef __BOX__
#define __BOX__

#include "interval.h"

using namespace std;

class box
{
public:
	interval* data;
	int dim;
	//---------------------CONSTRUCTEURS--------------------------------------------
	box();
	box(int);
	box(interval, int);
	box(interval x);
	box(interval x, interval y);
	box(interval x, interval y, interval z);
	box(const box&);
	~box();
	interval& operator[] (int) const;
	box& operator=(const box&);
	void Resize(int);
	box& Intersect(const box& Y);
	bool IsEmpty(void) const;
	double Width(void);
	double SumWidth(void);
	//----------------------OPERATEURS----------------------------------------------
	friend box operator&(const box&, const box&);
	friend box operator|(const box&, const box&);
	friend box operator+(const box&, const box&);
	friend box operator-(const box&);
	friend box operator-(const box&, const box&);
	friend box operator*(const interval&, const box&);
	friend box operator*(const double, const box&);
	friend box operator*(const box&, const interval&);
	friend bool operator==(const box&, const box&);
	friend ostream& operator<<(ostream&, const box&);
#ifdef QT_VERSION 
	friend QDebug operator<<(QDebug, const box&);
#endif // QT_VERSION 
	//-------------------FONCTIONS NON MEMBRES---------------------------------------
	friend double       Angle(box&, box&); // Il faut des vecteurs de dim 2
	friend int          Size(const box&);
	friend int	    AxePrincipal(box&);
	friend int	    AxePrincipal(box&, vector<int>&);
	friend int	    AxePrincipal(box&, box&);
	friend void         Update(box&);
	friend void	    Bisect(box&, box&, box&);
	friend void	    Bisect(box&, box&, box&, vector<int>&);
	friend void	    Decoup(box&, box&, box&);
	friend void         Trisect(box&, box&, box&, box&);
	friend void	    Bisect(box&, box&, box&, box&);
	friend void	    BisectAlong(box&, box&, box&, int);
	friend void	    DecoupAlong(box&, box&, box&, int);
	friend void         TrisectAlong(box&, box&, box&, box&, int);
	friend void         BisectHere(box&, box&, box&, int, double);
	friend box          Rand(const box& X);
	friend box	    Center(const box&);
	friend box          Center(const box&, vector<int>&);
	//friend void CheckRange(box&,box&);
	friend interval     Determinant(box&, box&);
	friend bool	    Disjoint(const box&, const box&);
	friend double         decrease(const box&, const box&);
	friend double         decrease(const box&, const box&, vector<int>);
	friend box EmptyBox(const box&);
	friend box EmptyBox(int);
	friend double         Eloignement(box&, box&);
	friend double         Eloignement2(box&, box&);
	friend double         EloignementRelatif2(box&, box&);
	friend double         Marge(box, box);
	friend iboolean	    In(box, box);
	friend box	    Inf(box);
	friend box	    Inflate(box&, double);
	friend box	    Inter(const box&, const box&);
	friend box          Inter(vector<box>&);
	friend box          Union(vector<box>&);
	friend box          Concat(const box&, const box&);
	friend box          Proj(const box&, int, int);
	//friend void       Inter1(box&,box&,const box&,const box&,const box&);
	friend interval     Norm(box);
	friend interval     NormEuclid(box, box);
	friend interval     NormInf(box, box);
	friend interval     ProduitScalaire(box&, box&);
	friend bool         Prop(box&, box&);
	friend bool	    Subset(box&, box&);
	friend bool         SubsetStrict(box&, box&);
	friend bool	    Subset(box&, box&, double);
	friend bool         IsBox(box);
	friend void         Sucre(box&, box&);
	friend box	    Sup(box);
	friend box	    Union(const box&, const box&);
	friend double         Volume(box&);
	friend double	    Width(box&);
	friend double	    Width(box&, vector<int>&);
	friend double	    Width(box&, box&);
	friend box Zeros(int);
	friend box Empty(int);
	friend box Infinity(int);

	friend void Cplus(box&, box&, box&, int sens);
	friend void Cmoins(box&, box&, box&, int sens);
	friend void C_q_in(box&, int, vector<box>&);
	friend void Cnorm(interval&R, box& X);
	friend void Cdistance(interval& R, box& X, box& Y);
	friend void CProdScalaire(interval& R, box& X, box& Y);
	friend void COrtho(box& X, box& Y);
	friend void CProd(box& Y, interval& a, box& X, int sens);

	// Operation sur les boites
	friend vector<box>* diff(box x, box y);
};

//box Empty(int);

#endif // __BOX__

/* Copyright (c) 2012, Christoph Reimann */  
 
/* change the following macros if orbitals are to be arranged in a different order */
#define L_QN(l) (l)
/* for matrices use indices ranging from m_MIN to m_MAX */
#define M_INDEX_MIN(l) (-(l))
#define M_INDEX_MAX(l) (+(l))
/* map mIndex to the quantum number m */
#define M_QN(l,mIndex) M_INDEX_MIN(l)+(mIndex)
/* determine the index from the quantum number */
#define M_INDEX(l,m) (m)-M_INDEX_MIN(l)

/* angular momentum quantum numbers l,m in arrays are combined
   into one dimension by this macro
   l:        0  1  2  3  4  ...
   L_DIM:    1  4  9 16 25  ...
   M_DIM:    1  3  5  7  9  ... */   
#define L_DIM(l) ((l)+1)*((l)+1)
#define M_DIM(l) (2*(l)+1)
#define LM_INDEX(l,m) ((l)*(l)+(m))

/* gaussians: cartesian representation 
   l:        0  1  2  3  4  ...
   IJK_DIM:  1  3  6 10 15  ...
   C_DIM:    1  4 10 20 35  ... */
#define C_DIM(l) ((l)+1)*((l)+2)*((l)+3)/6
#define C_INDEX(l,c) (C_DIM((l)-1)+(c))
#define IJK_DIM(l) ((l)+1)*((l)+2)/2
#define CIJK_INDEX(l,c) (C_DIM((l)-1)+(c))*3

/* exponents nx,ny,nz for cartesian gaussian shells with
   angular momentum up to am
   default: standard LIBINT ordering, i.e.
            p : px , py , pz
	    d : dxx , dxy , dxz , dyy , dyz , dzz
	    f : fxxx , fxxy , fxxz , fxyy , fxyz , fxzz , fyyy , fyyz , fyzz , fzzz
	    etc. 
   parameter: max. angular momentum */
int * cartesianShellOrder(const int am);
int * cartesianShellOrderIndex(const int am, int *ijk);

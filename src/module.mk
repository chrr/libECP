
libecp_SRC =	ecp.c\
		type2.c\
		type1.c\
		dimensions.c\
		spherical_harmonics.c\
		transformations.c\
		angular_integrals.c\
		util.c\
		bessel.c\
		gc_integrators.c\
		ecp_array.c\
		libecp.c\
		getIntegrals.c

SRC += $(patsubst %, src/%, $(libecp_SRC)) 

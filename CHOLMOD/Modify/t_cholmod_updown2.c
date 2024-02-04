//------------------------------------------------------------------------------
// CHOLMOD/Modify/t_cholmod_updown2: template for cholmod_updown2
//------------------------------------------------------------------------------

// Author: Roger Hsiao

//------------------------------------------------------------------------------

/* Updates/downdates the LDL' factorization, by computing a new factorization of
 *
 *   	Lnew * Dnew * Lnew' = Lold * Dold * Lold' +/- C*C'
 *
 * This file is not compiled separately.  It is included into
 * cholmod_updown2.c.  There are no user-callable routines in this file.
 *
 * The next include statements, below, create the numerical update/downdate
 * kernels from t_cholmod_updown2_numkr.c.  There are 4 compiled versions of this
 * file, one for each value of WDIM in the set 1, 2, 4, and 8.  Each calls
 * multiple versions of t_cholmod_updown2_numkr; the number of versions of each
 * is equal to WDIM.  Each t_cholmod_updown2_numkr version is included as a
 * static function within its t_cholmod_updown2.c caller routine.  Thus:
 *
 *	t*_updown2.c	creates these versions of t_cholmod_updown2_numkr.c:
 *	---------	---------------------------------------------------
 *
 *	updown2_1_r	updown2_1_1
 *
 *	updown2_2_r	updown2_2_1     updown2_2_2
 *
 *	updown2_4_r	updown2_4_1     updown2_4_2     updown2_4_3     updown2_4_4
 *
 *	updown2_8_r	updown2_8_1     updown2_8_2     updown2_8_3     updown2_8_4
 *			updown2_8_5     updown2_8_6     updown2_8_7     updown2_8_8
 *
 * workspace: Xwork (nrow*wdim)
 */

/* ========================================================================== */
/* === routines for numeric update/downdate along one path ================== */
/* ========================================================================== */

// #include <stdint.h>
// #include "../Include/cholmod.h"
// #include "../Include/cholmod_internal.h"
// #include "../Include/cholmod_template.h"
// #include "../Include/cholmod_types.h"
// #ifndef WDIM
// #define WDIM 8
// #endif


#undef FORM_NAME
#undef NUMERIC

#define FORM_NAME(k,rank) updown2_ ## k ## _ ## rank
#define NUMERIC(k,rank) FORM_NAME(k,rank)

#define RANK 1
#include "t_cholmod_updown2_numkr.c"

#if WDIM >= 2
#define RANK 2
#include "t_cholmod_updown2_numkr.c"
#endif

#if WDIM >= 4
#define RANK 3
#include "t_cholmod_updown2_numkr.c"
#define RANK 4
#include "t_cholmod_updown2_numkr.c"
#endif

#if WDIM == 8
#define RANK 5
#include "t_cholmod_updown2_numkr.c"
#define RANK 6
#include "t_cholmod_updown2_numkr.c"
#define RANK 7
#include "t_cholmod_updown2_numkr.c"
#define RANK 8
#include "t_cholmod_updown2_numkr.c"
#endif

/* ========================================================================== */
/* === numeric update/downdate for all paths ================================ */
/* ========================================================================== */

static void NUMERIC (WDIM, r)
(
    cholmod_sparse *C,	/* in packed or unpacked, and sorted form */
			/* no empty columns */
    cholmod_sparse *D,	/* in packed or unpacked, and sorted form */
			/* no empty columns */
    Int rank,		/* rank of the update/downdate */
    cholmod_factor *L,	/* with unit diagonal (diagonal not stored) */
			/* temporary workspaces: */
    double WC [ ],	/* n-by-WDIM dense matrix, initially zero */
    double WD [ ],	/* n-by-WDIM dense matrix, initially zero */
    Path_type Path [ ],
    Int npaths,
    Int mask [ ],	/* size n */
    Int maskmark,
    cholmod_common *Common
)
{
    // Roger notes: 
    // 1. The workspace W stores the columns of C that we are working on (in row major order). 
    //    We allocated 2x the workspace to store the columns of D. For ease of implementation
    //    right now, we will split W into WC and WD, where WD is appended to WC
    // 2. There are up to WDIM paths in W. 
    //    WCpath/WDpath denotes the start of the column of path
    // 3. WC1 and WD1 point to the wfirst column of C and D
    // 4. aC and aD point to the current AlphaC and alphaD, respectively

    double AlphaC [8] ;
    double AlphaD [8] ;
    double *Cx, *Dx, *WCpath, *WDpath, *WC1, *WD1, *aC, *aD ;
    Int i, j, p, ccol, pend, wfirst, e, path, packed ;
    Int *Ci, *Cp, *Cnz, *Di, *Dp, *Dnz ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    Ci = C->i ;
    Cx = C->x ;
    Cp = C->p ;
    Cnz = C->nz ;
    packed = C->packed ;
    ASSERT (IMPLIES (!packed, Cnz != NULL)) ;
    ASSERT (L->n == C->nrow) ;
    DEBUG (CHOLMOD(dump_real) ("num_d: in W:", W, WDIM, L->n, FALSE, 1,Common));

    Di = D->i ;
    Dx = D->x ;
    Dp = D->p ;
    Dnz = D->nz ;

    /* ---------------------------------------------------------------------- */
    /* scatter C into W */
    /* ---------------------------------------------------------------------- */

    __itt_frame_begin_v3(CHOLMOD(scatter_domain), NULL);

    for (path = 0 ; path < rank ; path++)
    {
	/* W (:, path) = C (:, Path [path].col) */
	ccol = Path [path].ccol ;
	WCpath = WC + path ;
	WDpath = WD + path ;
	PRINT1 (("Ordered Columns [path = "ID"] = "ID"\n", path, ccol)) ;
	p = Cp [ccol] ;
	pend = (packed) ? (Cp [ccol+1]) : (p + Cnz [ccol]) ;
	/* column C can be empty */
	for ( ; p < pend ; p++)
	{
	    i = Ci [p] ;
	    ASSERT (i >= 0 && i < (Int) (C->nrow)) ;
	    if (mask == NULL || mask [i] < maskmark)
	    {
		WCpath [WDIM * i] = Cx [p] ;
		WDpath [WDIM * i] = Dx [p] ;
	    }
	    PRINT1 (("    row "ID" : %g mask "ID"\n", i, Cx [p],
		    (mask) ? mask [i] : 0)) ;
	    PRINT1 (("    row "ID" : %g mask "ID"\n", i, Dx [p],
		    (mask) ? mask [i] : 0)) ;
	}
	AlphaC [path] = 1.0 ;
	AlphaD [path] = 1.0 ;
    }
    DEBUG (CHOLMOD(dump_real) ("num_d: WC:", WC, WDIM, L->n, FALSE, 1,Common)) ;
    DEBUG (CHOLMOD(dump_real) ("num_d: WD:", WD, WDIM, L->n, FALSE, 1,Common)) ;

    __itt_frame_end_v3(CHOLMOD(scatter_domain), NULL);

    /* ---------------------------------------------------------------------- */
    /* numeric update/downdate of the paths */
    /* ---------------------------------------------------------------------- */

    /* for each disjoint subpath in Tbar in DFS order do */
    for (path = rank ; path < npaths ; path++)
    {

	/* determine which columns of W to use */
	wfirst = Path [path].wfirst ;
	e = Path [path].end ;
	j = Path [path].start ;
	ASSERT (e >= 0 && e < (Int) (L->n)) ;
	ASSERT (j >= 0 && j < (Int) (L->n)) ;

	WC1 = WC + wfirst ;	/* pointer to row 0, column wfirst of W */
	WD1 = WD + wfirst ;	/* pointer to row 0, column wfirst of W */
	aC = AlphaC + wfirst ;	/* pointer to AlphaC [wfirst] */
	aD = AlphaD + wfirst ;	/* pointer to AlphaD [wfirst] */

	PRINT1 (("Numerical update/downdate of path "ID"\n", path)) ;
	PRINT1 (("start "ID" end "ID" wfirst "ID" rank "ID" ccol "ID"\n", j, e,
		wfirst, Path [path].rank, Path [path].ccol)) ;

#if WDIM == 1
	NUMERIC (WDIM,1) (j, e, aC, aD, WC1, WD1, L, Common) ;
#else

	switch (Path [path].rank)
	{
	    case 1:
		NUMERIC (WDIM,1) (j, e, aC, aD, WC1, WD1, L, Common) ;
		break ;

#if WDIM >= 2
	    case 2:
		NUMERIC (WDIM,2) (j, e, aC, aD, WC1, WD1, L, Common) ;
		break ;
#endif

#if WDIM >= 4
	    case 3:
		NUMERIC (WDIM,3) (j, e, aC, aD, WC1, WD1, L, Common) ;
		break ;
	    case 4:
		NUMERIC (WDIM,4) (j, e, aC, aD, WC1, WD1, L, Common) ;
		break ;
#endif

#if WDIM == 8
	    case 5:
		NUMERIC (WDIM,5) (j, e, aC, aD, WC1, WD1, L, Common) ;
		break ;
	    case 6:
		NUMERIC (WDIM,6) (j, e, aC, aD, WC1, WD1, L, Common) ;
		break ;
	    case 7:
		NUMERIC (WDIM,7) (j, e, aC, aD, WC1, WD1, L, Common) ;
		break ;
	    case 8:
        NUMERIC (WDIM,8) (j, e, aC, aD, WC1, WD1, L, Common) ;
        break ;
#endif

	}
#endif

    }
}

/* prepare for the next inclusion of this file in cholmod_updown2.c */
#undef WDIM

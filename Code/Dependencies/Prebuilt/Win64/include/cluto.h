/*
 * Copyright 2003 Regents of the University of Minnesota. All Rights Reserved
 *
 * cluto.h 
 *
 * This file contains function prototypes and constant definitions 
 * for CLUTO
 *
 * Started 11/23/01
 * George
 *
 */

#ifndef __cluto_h__
#define __cluto_h__


#if !defined __cdecl
#define __cdecl
#endif

#if !defined _declspec
#define _declspec(x)
#endif

/*------------------------------------------------------------------------
* Function prototypes 
*-------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

_declspec(dllexport) void __cdecl 
     __CLUTO_VP_ClusterDirect(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int crfun, int rowmodel, int colmodel, 
              float colprune, int ntrials, int niter, int seed, int dbglvl, 
              int nparts, int *part, int *xcode);
_declspec(dllexport) void __cdecl 
     CLUTO_VP_ClusterDirect(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int crfun, int rowmodel, int colmodel, 
              float colprune, int ntrials, int niter, int seed, int dbglvl, 
              int nparts, int *part);


_declspec(dllexport) void __cdecl 
     __CLUTO_VP_ClusterRB(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int crfun, int rowmodel, int colmodel, 
              float colprune, int ntrials, int niter, int seed, int rbtype, 
              int kwayrefine, int dbglvl, int nparts, int *part, int *xcode);
_declspec(dllexport) void __cdecl 
     CLUTO_VP_ClusterRB(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int crfun, int rowmodel, int colmodel, 
              float colprune, int ntrials, int niter, int seed, int rbtype, 
              int kwayrefine, int dbglvl, int nparts, int *part);


_declspec(dllexport) void __cdecl 
     __CLUTO_VP_ClusterRBTree(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int crfun, int rowmodel, int colmodel, 
              float colprune, int ntrials, int niter, int seed, int dbglvl, 
              int nparts, int *part, int *ptree, int *xcode);
_declspec(dllexport) void __cdecl 
     CLUTO_VP_ClusterRBTree(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int crfun, int rowmodel, int colmodel, 
              float colprune, int ntrials, int niter, int seed, int dbglvl, 
              int nparts, int *part, int *ptree);


_declspec(dllexport) void __cdecl 
     __CLUTO_VA_Cluster(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int crfun, int rowmodel, int colmodel, 
              float colprune, int dbglvl, int nparts, int *part, int *ptree, 
              float *tsims, float *gains, int *xcode);
_declspec(dllexport) void __cdecl 
     CLUTO_VA_Cluster(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int crfun, int rowmodel, int colmodel, 
              float colprune, int dbglvl, int nparts, int *part, int *ptree, 
              float *tsims, float *gains);


_declspec(dllexport) void __cdecl 
     __CLUTO_VA_ClusterBiased(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int crfun, int rowmodel, int colmodel, 
              float colprune, int dbglvl, int pnparts, int nparts, int *part, 
              int *ptree, float *tsims, float *gains, int *xcode);
_declspec(dllexport) void __cdecl 
     CLUTO_VA_ClusterBiased(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int crfun, int rowmodel, int colmodel, 
              float colprune, int dbglvl, int pnparts, int nparts, int *part, 
              int *ptree, float *tsims, float *gains);


_declspec(dllexport) void __cdecl 
     __CLUTO_SP_ClusterRB(int nrows, int *xadj, int *adjncy, float *adjwgt, 
              int crfun, int ntrials, int niter, int seed, int cstype, int kwayrefine, 
              int dbglvl, int nparts, int *part, int *xcode);
_declspec(dllexport) void __cdecl 
     CLUTO_SP_ClusterRB(int nrows, int *xadj, int *adjncy, float *adjwgt, 
              int crfun, int ntrials, int niter, int seed, int cstype, int kwayrefine, 
              int dbglvl, int nparts, int *part);


_declspec(dllexport) void __cdecl 
    __CLUTO_SP_ClusterDirect(int nrows, int *xadj, int *adjncy, float *adjwgt, 
              int crfun, int ntrials, int niter, int seed, int dbglvl, int nparts, 
              int *part, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_SP_ClusterDirect(int nrows, int *xadj, int *adjncy, float *adjwgt, 
              int crfun, int ntrials, int niter, int seed, int dbglvl, int nparts, 
              int *part);


_declspec(dllexport) void __cdecl 
    __CLUTO_SA_Cluster(int nvtxs, int *xadj, int *adjncy, float *adjwgt,
              int crfun, int dbglvl, int nparts, int *part, int *ptree, 
              float *tsims, float *gains, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_SA_Cluster(int nvtxs, int *xadj, int *adjncy, float *adjwgt,
              int crfun, int dbglvl, int nparts, int *part, int *ptree, 
              float *tsims, float *gains);


_declspec(dllexport) int  __cdecl 
    __CLUTO_VP_GraphClusterRB(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int rowmodel, int colmodel, float colprune, 
              int grmodel, int nnbrs, float edgeprune, float vtxprune, int mincmp, 
              int ntrials, int seed, int rbtype, int dbglvl, int nparts, int *part, 
              float *r_crvalue, int *xcode);
_declspec(dllexport) int  __cdecl 
    CLUTO_VP_GraphClusterRB(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int rowmodel, int colmodel, float colprune, 
              int grmodel, int nnbrs, float edgeprune, float vtxprune, int mincmp, 
              int ntrials, int seed, int rbtype, int dbglvl, int nparts, int *part, 
              float *r_crvalue);


_declspec(dllexport) int __cdecl  
    __CLUTO_SP_GraphClusterRB(int nvtxs, int *xadj, int *adjncy, float *adjwgt,
              int nnbrs, float edgeprune, float vtxprune, int mincmp, int ntrials, 
              int seed, int cstype, int dbglvl, int nparts, int *part, 
              float *r_crvalue, int *xcode);
_declspec(dllexport) int __cdecl  
    CLUTO_SP_GraphClusterRB(int nvtxs, int *xadj, int *adjncy, float *adjwgt,
              int nnbrs, float edgeprune, float vtxprune, int mincmp, int ntrials, 
              int seed, int cstype, int dbglvl, int nparts, int *part, float *r_crvalue);


_declspec(dllexport) void __cdecl 
    __CLUTO_V_GetGraph(int nrows, int ncols, int *rowptr, int *rowind, float *rowval, 
              int simfun, int rowmodel, int colmodel, float colprune, int grmodel, int nnbrs, 
              int dbglvl, int **r_xadj, int **r_adjncy, float **r_adjwgt, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_V_GetGraph(int nrows, int ncols, int *rowptr, int *rowind, float *rowval, 
              int simfun, int rowmodel, int colmodel, float colprune, int grmodel, int nnbrs, 
              int dbglvl, int **r_xadj, int **r_adjncy, float **r_adjwgt);


_declspec(dllexport) void __cdecl 
    __CLUTO_S_GetGraph(int nvtxs, int *xadj, int *adjncy, float *adjwgt, int grmodel, 
               int nnbrs, int dbglvl, int **r_xadj, int **r_adjncy, float **r_adjwgt,
               int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_S_GetGraph(int nvtxs, int *xadj, int *adjncy, float *adjwgt, int grmodel, 
               int nnbrs, int dbglvl, int **r_xadj, int **r_adjncy, float **r_adjwgt);


_declspec(dllexport) float __cdecl 
    __CLUTO_V_GetSolutionQuality(int nrows, int ncols, int *rowptr, int *rowind,  
               float *rowval, int simfun, int crfun, int rowmodel, int colmodel, 
               float colprune, int nparts, int *part, int *xcode);
_declspec(dllexport) float __cdecl 
    CLUTO_V_GetSolutionQuality(int nrows, int ncols, int *rowptr, int *rowind,  
               float *rowval, int simfun, int crfun, int rowmodel, int colmodel, 
               float colprune, int nparts, int *part);


_declspec(dllexport) float __cdecl 
    __CLUTO_S_GetSolutionQuality(int nvtxs, int *xadj, int *adjncy,  
               float *adjwgt, int crfun, int nparts, int *part, int *xcode);
_declspec(dllexport) float __cdecl 
    CLUTO_S_GetSolutionQuality(int nvtxs, int *xadj, int *adjncy,  
               float *adjwgt, int crfun, int nparts, int *part);


_declspec(dllexport) void __cdecl 
    __CLUTO_V_GetClusterStats(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int rowmodel, int colmodel, float colprune, 
              int nparts, int *part, int *pwgts, 
	      float *cintsim, float *cintsdev, float *izscores, 
	      float *cextsim, float *cextsdev, float *ezscores, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_V_GetClusterStats(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int rowmodel, int colmodel, float colprune, 
              int nparts, int *part, int *pwgts, 
	      float *cintsim, float *cintsdev, float *izscores, 
	      float *cextsim, float *cextsdev, float *ezscores);


_declspec(dllexport) void __cdecl 
    __CLUTO_S_GetClusterStats(int nvtxs, int *xadj, int *adjncy, float *adjwgt, 
              int nparts, int *part, int *pwgts, 
	      float *cintsim, float *cintsdev, float *izscores, 
	      float *cextsim, float *cextsdev, float *ezscores, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_S_GetClusterStats(int nvtxs, int *xadj, int *adjncy, float *adjwgt, 
              int nparts, int *part, int *pwgts, 
	      float *cintsim, float *cintsdev, float *izscores, 
	      float *cextsim, float *cextsdev, float *ezscores);


_declspec(dllexport) void __cdecl 
    __CLUTO_V_GetClusterFeatures(int nrows, int ncols, int *rowptr,
              int *rowind, float *rowval, int simfun, int rowmodel, int colmodel, 
              float colprune, int nparts, int *part, int nfeatures,
              int *internalids, float *internalwgts, int *externalids,
              float *externalwgts, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_V_GetClusterFeatures(int nrows, int ncols, int *rowptr,
              int *rowind, float *rowval, int simfun, int rowmodel, int colmodel, 
              float colprune, int nparts, int *part, int nfeatures,
              int *internalids, float *internalwgts, int *externalids,
              float *externalwgts);


_declspec(dllexport) void __cdecl 
    __CLUTO_V_BuildTree(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int crfun, int rowmodel, int colmodel, 
              float colprune, int treetype, int dbglvl, int nparts, int *part, 
              int *ptree, float *tsims, float *gains, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_V_BuildTree(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int crfun, int rowmodel, int colmodel, 
              float colprune, int treetype, int dbglvl, int nparts, int *part, 
              int *ptree, float *tsims, float *gains);


_declspec(dllexport) void __cdecl 
    __CLUTO_S_BuildTree(int nvtxs, int *xadj, int *adjncy, float *adjwgt, 
              int crfun, int treetype, int dbglvl, int nparts, int *part, int *ptree, 
              float *tsims, float *gains, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_S_BuildTree(int nvtxs, int *xadj, int *adjncy, float *adjwgt, 
              int crfun, int treetype, int dbglvl, int nparts, int *part, int *ptree, 
              float *tsims, float *gains);


_declspec(dllexport) void __cdecl 
    __CLUTO_V_GetTreeStats(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int rowmodel, int colmodel, float colprune, 
              int nparts, int *part, int *ptree, int *pwgts, float *cintsim, 
              float *cextsim, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_V_GetTreeStats(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int rowmodel, int colmodel, float colprune, 
              int nparts, int *part, int *ptree, int *pwgts, float *cintsim, 
              float *cextsim);


_declspec(dllexport) void __cdecl 
    __CLUTO_V_GetTreeFeatures(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int rowmodel, int colmodel, float colprune, 
              int nparts, int *part, int *ptree, int nfeatures, int *internalids,
              float *internalwgts, int *externalids, float *externalwgts, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_V_GetTreeFeatures(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int rowmodel, int colmodel, float colprune, 
              int nparts, int *part, int *ptree, int nfeatures, int *internalids,
              float *internalwgts, int *externalids, float *externalwgts);


_declspec(dllexport) void __cdecl 
    __CLUTO_InternalizeMatrix(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int rowmodel, int colmodel, float colprune, 
              int *part, int *r_nrows, int *r_ncols, int **r_rowptr, int **r_rowind, 
              float **r_rowval, int **r_rimap, int **r_cimap, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_InternalizeMatrix(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int rowmodel, int colmodel, float colprune, 
              int *part, int *r_nrows, int *r_ncols, int **r_rowptr, int **r_rowind, 
              float **r_rowval, int **r_rimap, int **r_cimap);


_declspec(dllexport) void __cdecl 
    __CLUTO_S_TreeReorderInternal(int nrows, int *rwgts, float *smat, int memflag,
              int dbglvl, int *ptree, int **ftree, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_S_TreeReorderInternal(int nrows, int *rwgts, float *smat, int memflag,
              int dbglvl, int *ptree, int **ftree);


_declspec(dllexport) void __cdecl 
    __CLUTO_V_TreeReorder(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int rowmodel, int colmodel, float colprune, 
              int dbglvl, int *ptree, int **ftree, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_V_TreeReorder(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int rowmodel, int colmodel, float colprune, 
              int dbglvl, int *ptree, int **ftree);


_declspec(dllexport) void __cdecl 
    __CLUTO_S_TreeReorder(int nvtxs, int *xadj, int *adjncy, float *adjwgt, 
              int dbglvl, int *ptree, int **ftree, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_S_TreeReorder(int nvtxs, int *xadj, int *adjncy, float *adjwgt, 
              int dbglvl, int *ptree, int **ftree);


_declspec(dllexport) void __cdecl 
    __CLUTO_V_ClusterTreeReorder(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int rowmodel, int colmodel, float colprune, 
              int dbglvl, int nparts, int *part, int *ptree, int **ftree, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_V_ClusterTreeReorder(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int rowmodel, int colmodel, float colprune, 
              int dbglvl, int nparts, int *part, int *ptree, int **ftree);


_declspec(dllexport) void __cdecl 
    __CLUTO_S_ClusterTreeReorder(int nvtxs, int *xadj, int *adjncy, float *adjwgt, 
              int dbglvl, int nparts, int *part, int *ptree, int **ftree, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_S_ClusterTreeReorder(int nvtxs, int *xadj, int *adjncy, float *adjwgt, 
              int dbglvl, int nparts, int *part, int *ptree, int **ftree);


_declspec(dllexport) void __cdecl 
    __CLUTO_V_ReorderPartitions(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int rowmodel, int colmodel, float colprune,
              int nparts, int *part, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_V_ReorderPartitions(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int rowmodel, int colmodel, float colprune,
              int nparts, int *part);


_declspec(dllexport) void __cdecl 
    __CLUTO_S_ReorderPartitions(int nvtxs, int *xadj, int *adjncy, float *adjwgt,
              int nparts, int *part, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_S_ReorderPartitions(int nvtxs, int *xadj, int *adjncy, float *adjwgt,
              int nparts, int *part);


_declspec(dllexport) void __cdecl 
    __CLUTO_V_GetClusterDistanceMatrix(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int rowmodel, int colmodel, float colprune,
              int dbglvl, int nparts, int *part, float *distmat, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_V_GetClusterDistanceMatrix(int nrows, int ncols, int *rowptr, int *rowind,
              float *rowval, int simfun, int rowmodel, int colmodel, float colprune,
              int dbglvl, int nparts, int *part, float *distmat);


_declspec(dllexport) void __cdecl 
    __CLUTO_S_GetClusterDistanceMatrix(int nvtxs, int *xadj, int *adjncy, float *adjwgt, 
              int dbglvl, int nparts, int *part, float *distmat, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_S_GetClusterDistanceMatrix(int nvtxs, int *xadj, int *adjncy, float *adjwgt, 
              int dbglvl, int nparts, int *part, float *distmat);


_declspec(dllexport) void __cdecl 
    __CLUTO_V_GetClusterSummaries(int nrows, int ncols, int *rowptr, int *rowind, 
              float *rowval, int simfun, int rowmodel, int colmodel, float colprune, int nparts, 
              int *part, int sumtype, int nfeatures, int *r_nsum, int **r_spid, float **r_swgt, 
              int **r_sumptr, int **r_sumind, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_V_GetClusterSummaries(int nrows, int ncols, int *rowptr, int *rowind, 
              float *rowval, int simfun, int rowmodel, int colmodel, float colprune, int nparts, 
              int *part, int sumtype, int nfeatures, int *r_nsum, int **r_spid, float **r_swgt, 
              int **r_sumptr, int **r_sumind);


_declspec(dllexport) void __cdecl 
    __CLUTO_V_GetTopNSimilarities(int nrows, int ncols, int *rowptr, int *rowind, float *rowval, 
              int simfun, int rowmodel, int colmodel, float colprune, int nparts, int *part,
              int topn, int *clids, float *clsims, int *xcode);
_declspec(dllexport) void __cdecl 
    CLUTO_V_GetTopNSimilarities(int nrows, int ncols, int *rowptr, int *rowind, float *rowval, 
              int simfun, int rowmodel, int colmodel, float colprune, int nparts, int *part,
              int topn, int *clids, float *clsims);
 

#ifdef __cplusplus
}
#endif




/*------------------------------------------------------------------------
* Constant definitions 
*-------------------------------------------------------------------------*/

/* Different choices for RowModel */
#define CLUTO_ROWMODEL_NONE           1       /* Use the rows as is */
#define CLUTO_ROWMODEL_MAXTF          2       /* TF = .5 + .5(TF/max(|TF|)) */
#define CLUTO_ROWMODEL_SQRT           3       /* TF = 1+sign(TF, sqrt(|TF|)) */
#define CLUTO_ROWMODEL_LOG            4       /* TF = 1+sign(TF, log2(|TF|)) */


/* Different choices for ColModel */
#define CLUTO_COLMODEL_NONE           1       /* Use the columns as is */
#define CLUTO_COLMODEL_IDF            2       /* Scale-the columns following a binary IDF model */


/* Different cluster criterion functions */
#define CLUTO_CLFUN_I1       1       /* The I1 from the paper */
#define CLUTO_CLFUN_I2       2       /* The I2 from the paper */
#define CLUTO_CLFUN_E1       3       /* The E1 from the paper */
#define CLUTO_CLFUN_G1       4       /* The G1 from the paper */
#define CLUTO_CLFUN_G1P      5       /* The G1' from the paper */
#define CLUTO_CLFUN_H1       6       /* The H1 from the paper */
#define CLUTO_CLFUN_H2       7       /* The H2 from the paper */
#define CLUTO_CLFUN_SLINK    8       /* The traditional single-link / MST method */
#define CLUTO_CLFUN_SLINK_W  9       /* The traditional single-link / MST method, cluster size weighted */
#define CLUTO_CLFUN_CLINK    10      /* The traditional complete-link method */
#define CLUTO_CLFUN_CLINK_W  11      /* The traditional complete-link method, cluster size weighted */
#define CLUTO_CLFUN_UPGMA    12      /* The traditional UPGMA method */
#define CLUTO_CLFUN_UPGMA_W  13      /* The traditional weighted UPGMA method */

/* The following are criterion functions for graph-based clustering */
#define CLUTO_CLFUN_CUT     15      /* Edge-cut based */
#define CLUTO_CLFUN_RCUT    16      /* Ratio cut */
#define CLUTO_CLFUN_NCUT    17      /* Normalized cut */
#define CLUTO_CLFUN_MMCUT   18      /* Min-Max cut */


/* Different cluster selection schemes for RB */
#define CLUTO_CSTYPE_LARGEFIRST         1       /* Select the largest cluster to bisect */
#define CLUTO_CSTYPE_BESTFIRST          2       /* Select the cluster that leads the best value of the criterion function */
#define CLUTO_CSTYPE_LARGESUBSPACEFIRST 3       /* Selects the cluster that leads to the largest subspace reduction */

/* Different dbglvl options */
#define CLUTO_DBG_PROGRESS         1       /* Show simple progress statistics */
#define CLUTO_DBG_RPROGRESS        2       /* Show simple progress statistics during refinement */
#define CLUTO_DBG_APROGRESS        4       /* Show progress during the agglomeration */
#define CLUTO_DBG_CPROGRESS        8       /* Show progress statistics during coarsening */
#define CLUTO_DBG_MPROGRESS        16      /* Show vertex movement information during refinement */
#define CLUTO_DBG_CCMPSTAT         32      /* Show stats during cc elimination */


/* Different option for memory re-use for the SA-routines */
#define CLUTO_MEM_NOREUSE       1       /* Preserves the supplied information */
#define CLUTO_MEM_REUSE         2       /* Does not preserve the supplied information */

/* Different types of trees that CLUTO can build on top of the clustering solution */
#define CLUTO_TREE_TOP       1       /* Builds a tree on top of the supplied clustering */
#define CLUTO_TREE_FULL      2       /* Builds the entire tree that preserves the clustering */

/* Different similarity functions that CLUTO supports */
#define CLUTO_SIM_COSINE        1       /* Similarity is measured using the cosine function */
#define CLUTO_SIM_CORRCOEF      2       /* Similarity is measured using Pearson's correlatio coefficient */
#define CLUTO_SIM_EDISTANCE     3       /* Similarity is measured using the negative Euclidean distance */
#define CLUTO_SIM_EJACCARD      4       /* Similarity is measured using the extended Jaccard */

/* Different types of optimizers implemeted by CLUTO */
#define CLUTO_OPTIMIZER_SINGLELEVEL     1       /* Traditional single-level optimizer */
#define CLUTO_OPTIMIZER_MULTILEVEL      2       /* Better multi-level optimizer */

/* Different ways for performing the graph coarsening */
#define CLUTO_MTYPE_HEDGE       1       /* Heavy-edge matching */
#define CLUTO_MTYPE_HSTAR       2       /* Heavy-star matching */
#define CLUTO_MTYPE_HSTAR2      3       /* Heavy-star matching */


/* Different type of neighborhood graph models */
#define CLUTO_GRMODEL_EXACT_ASYMETRIC_DIRECT    1       /* Computes similarity exactly, and 
                                                           includes edges for all of them */
#define CLUTO_GRMODEL_EXACT_SYMETRIC_DIRECT     2       /* Computes similarity exactly, and 
                                                           includes edges only if they are shared */
#define CLUTO_GRMODEL_INEXACT_ASYMETRIC_DIRECT  3       /* Computes most similar vertices inexactly, 
                                                           and includes edges for all of them */
#define CLUTO_GRMODEL_INEXACT_SYMETRIC_DIRECT   4       /* Computes most similar vertices inexactly, 
                                                           includes edges only if they are shared */
#define CLUTO_GRMODEL_EXACT_ASYMETRIC_LINKS     5       /* Computes similarity exactly, and 
                                                           includes edges for all of them */
#define CLUTO_GRMODEL_EXACT_SYMETRIC_LINKS      6       /* Computes similarity exactly, and 
                                                           includes edges only if they are shared */
#define CLUTO_GRMODEL_INEXACT_ASYMETRIC_LINKS   7       /* Computes most similar vertices inexactly, 
                                                           and includes edges for all of them */
#define CLUTO_GRMODEL_INEXACT_SYMETRIC_LINKS    8       /* Computes most similar vertices inexactly, 
                                                           includes edges only if they are shared */

#define CLUTO_GRMODEL_ASYMETRIC_DIRECT          CLUTO_GRMODEL_EXACT_ASYMETRIC_DIRECT
#define CLUTO_GRMODEL_ASYMETRIC_LINKS           CLUTO_GRMODEL_EXACT_ASYMETRIC_LINKS
#define CLUTO_GRMODEL_SYMETRIC_DIRECT           CLUTO_GRMODEL_EXACT_SYMETRIC_DIRECT
#define CLUTO_GRMODEL_SYMETRIC_LINKS            CLUTO_GRMODEL_EXACT_SYMETRIC_LINKS
#define CLUTO_GRMODEL_NONE                      9


/* Summary Types */
#define CLUTO_SUMMTYPE_MAXCLIQUES               1       /* Find maximal feature cliques */
#define CLUTO_SUMMTYPE_MAXITEMSETS              2       /* Find maximal itemset cliques */

/* Cluto's exit codes */
#define CLUTO_EXIT_NORMAL               0
#define CLUTO_EXIT_NOTENOUGHMEMORY      1
#define CLUTO_EXIT_ERROR                2

/* Cluto's version number */
#define CLUTO_VER_MAJOR         2
#define CLUTO_VER_MINOR         1
#define CLUTO_VER_SUBMINOR      2



#endif

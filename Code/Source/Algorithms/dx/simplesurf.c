/***********************************************************************/
/* Open Visualization Data Explorer                                    */
/* (C) Copyright IBM Corp. 1989,1999                                   */
/* ALL RIGHTS RESERVED                                                 */
/* This code licensed under the                                        */
/*    "IBM PUBLIC LICENSE - Open Visualization Data Explorer"          */
/***********************************************************************/

/* SimplifySurface

   is a module that approximates a triangulated surface and resamples data
   
   */

/*********************************************************************/
/*

  Author: Andre GUEZIEC, March/April 1997

 */
/*********************************************************************/


#include "dx.h"
#include "simplesurf.h"
#include <stdlib.h>
#include <string.h>

#define dx_bzero(b, len) (memset((b), '\0', (len)), (void)0)

#if !defined(MIN)
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#if !defined(MAX)
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

/* 
   this file contains the "guts" of the surface simplification procedure
   
   As much as possible I am using simple types 

   For best memory management, I have converted my routines so that they
   use the data explorer DXAllocate and DXFree calls

   */

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfSimplifySurfaceDriver                                              |
  |                                                                          |
  +--------------------------------------------------------------------------+*/


int _dxfSimplifySurfaceDriver(int nV, float *v, 
			      int nT, int *t,
			      int dim_data,
			      float *vertex_data,
			      float tolerance, 
			      float data_tolerance,
			      int *old_nV_after_manifold_conversion,
			      int *old_nT_after_manifold_conversion,
			      int *new_nV, float **new_v,
			      float **new_vertex_data,
			      int *new_nT, int   **new_t,
			      int   **vertex_parents,
			      int   **face_parents,
			      int   **vertex_lut, /* for fixing non-manifolds */
			      int   **face_lut,   /* same thing: for eliminating degenerate faces */
			      float *old_positional_error,
			      float **new_positional_error, 
			      float **face_normals,
			      float **old_face_areas,
			      float **new_face_areas,
			      int preserve_volume,
			      int simplify_boundary,
			      int preserve_boundary_length)

     /* all of the arrays that are ** must be allocated by the program
	if any of such arrays is not allocated properly, than the
	routine must return ERROR (zero). However, in case of an
	error, there is no need to DXFree() such arrays, as they are
	freed by the caller routine */

{
  int success = 0;

  int nE;

  EdgeS *e = NULL;

  Htable edge_table[1];

  /* create a simplification data structure where all the relevant arrays
     are copied */

  SimpData simp_data[1];
    
  dx_bzero(simp_data, sizeof (SimpData));
  
  dx_bzero(edge_table, sizeof (Htable));

  /* the first step is to convert the surface to a manifold, in case 
     Isosurface would have created non-manifold vertices (I know that happens sometimes)
     
     It is also possible that a surface provided directly by the user would not
     be a manifold */

  if (_dxfToOrientedManifold(nV, v, nT, t, new_nV, new_v, new_nT, new_t, vertex_lut, face_lut, 
			     &nE, edge_table, &e))
    {
      *old_nV_after_manifold_conversion = *new_nV;
      *old_nT_after_manifold_conversion = *new_nT;

      simp_data->nT       = *new_nT;       /* the number of vertices, edges and triangles  */
      simp_data->nE       = nE;            /* after conversion to manifold and before simplification */
      simp_data->nV       = *new_nV;
      simp_data->data_dim = (vertex_data) ? dim_data: 0;

      simp_data->preserve_volume          = preserve_volume;
      simp_data->simplify_boundary        = simplify_boundary;
      simp_data->preserve_boundary_length = preserve_boundary_length;

      simp_data->tolerance                = tolerance;
      simp_data->data_tolerance           = data_tolerance;

      *old_face_areas = (float *) DXAllocateZero( simp_data->nT * sizeof (float));

      *vertex_parents = (int   *) DXAllocateZero( simp_data->nV * sizeof (int));

      *face_parents   = (int   *) DXAllocateZero( simp_data->nT * sizeof (int));

      /* part of the simplification data structure that references existing arrays */

      simp_data->edge_array      = e;
      simp_data->e_table         = edge_table;
      simp_data->vert            = (Vertex *) *new_v;
      simp_data->triangles       = *new_t;
      simp_data->vertex_father   = *vertex_parents;
      simp_data->triangle_father = *face_parents;
      simp_data->old_face_areas  = *old_face_areas;
      simp_data->vertex_lut      = *vertex_lut;

      if (*old_face_areas && *vertex_parents && *face_parents &&
	  
	  (_dxfCreateSimplificationDataStructure(simp_data, vertex_data, old_positional_error)) &&

	  (_dxfSimplifyManifoldSurface(simp_data)))

	{
	  /* free only the simplification data structures that are not exported by the driver routine,
	     to avoid running out of memory in _dxfBuildSimplifiedSurface(): */

	  _dxfPartialFreeSimplificationDataStructure(simp_data);

	  /* build a simplified surface using the vertex_father and triangle_father arrays */

	  if (_dxfBuildSimplifiedSurface(simp_data, new_nV, new_v, new_nT, new_t,
					 new_face_areas, face_normals,
					 new_positional_error, new_vertex_data))
	      
	    success = 1;

	  _dxfFreeSimplificationDataStructure(simp_data, vertex_data, old_positional_error);
	
	}

        
      _dxfFreeHashTable(edge_table);
      if (e) DXFree((Pointer)e);

     
    }
  /* if the conversion to oriented manifold fails, vertex_lut, new_v, new_t ... etc.
     are freed in the caller routine */
  
  return success;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfToOrientedManifold                                                 |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfToOrientedManifold(int nV,   float *v,   int nT,   int *t, 
			   int *new_nV, float **new_v, int *new_nT, int **new_t, int **vlut, int **flut,
			   int *n_edges, Htable *e_table, EdgeS **edges)
	    
{
  int success = 0;

  /* we first eliminate the degenerate triangles */

  if (_dxfEliminateDegenerateTriangles(nV, nT, t, new_nT, new_t, flut))
    {
      /* after that the number of triangles does not change any more */

      nT = *new_nT;

      /* we then eliminate the standalone vertices */
      
      if (_dxfEliminateStandaloneVertices(nV,v,new_nV,new_v,nT,*new_t,vlut))
	{
	  /* we create a fathers array that is large enough to contain all
	     the triangle corners */

	  int *fathers = (int *) DXAllocateZero (3 * nT * sizeof (int));
	 
	  if (fathers)
	    {
	      if (_dxfJoinTriangleCorners(nT, *new_t, e_table, n_edges, edges, fathers))
		{  
		 
		  if (_dxfBuildOrientedManifold(new_nV, new_v, nT, *new_t, vlut, n_edges,
						edges, e_table, fathers))

		    { success = 1;}
		}
	    
	      DXFree((Pointer)fathers);
	    }
	}
     
    }

  return success;

}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfJoinTriangleCorners                                                |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfJoinTriangleCorners(int nT, int *t, Htable *e_table, int *n_edges, EdgeS **edges, int *fathers)

     /* we join two triangle corners only if they share a manifold edge 

	we define a manifold edge as an edge that has at most two incident triangles
	and such that the two triangles are consistently oriented. Also, the two triangles
	should not have the same three vertices */

{
  int success = 0;

  /* we first determine the edges that are manifold 

     in a second pass, we take all manifold edges and we join the corresponding triangle corners */

    
  register int i;

  register u_char j,k;

  int nT3 = 3 * nT;
      
  /* allocate an array big enough to contain all possible edges */
   
  EdgeS *edge_array = (EdgeS *) DXAllocateZero (nT3 * sizeof (EdgeS));

  if (edge_array) 
    {
       
      EdgeS 
	the_edge[1], 
	*cur_edge,
	*edge_in_table;
	
      int 
	*incident_tri  = SurfEdgeGetTriangles(the_edge),
	*endpoints     = SurfEdgeGetVertices (the_edge),
	nE             = 0, 
	total_nE       = nT3,

	/* initialize a pointer to the current triangle */

	*the_triangle = t,
	first_tri,  
	second_tri,
	corner_v0first, 
	corner_v1first, 
	corner_v0second,
	corner_v1second,
	edge_paired_as_manifold = 0,
	v0,                           /* surface edge end-points */
	v1,
	*triangles_in_table;
			
      /* initialize the fathers array */

      _dxfInitFather(nT3, fathers);
           
      dx_bzero((char *)e_table, sizeof(Htable)); /* reset the hash table */
     
      /* the size of the table should be chosen such that 2**tbl_size ~ 2 * number of  edges 
	 for best performances. Here taking nT is plenty enough*/
      
      if (!_dxfInitHashTable(e_table, 
			     _dxfFlog2(nT), 
			     (int (*)(KEYTYPE))_dxfEdgeSKey, 
			     (int (*)(RECTYPE, RECTYPE))_dxfEdgeSFilter, 
			     NULL))
	goto error;
		

      /* initialize the number of faces per edge to 1 */
       
      SurfEdgeLabel(the_edge) = 1; 

      /* for each triangle, add the three edges to the hash table if they are not already
	 indexed in the hash table */

      for (the_triangle = t, i=0; i<nT; i++, the_triangle += 3)
	{
	  
	  for (k=0;k<3;k++)
		    
	    {
	      edge_paired_as_manifold = 0;
	      v0 = the_triangle[k];       /* surface edge end-points */
	      v1 = the_triangle[(k+1)%3];

	      if (v0>=v1) 
		{ 
		  endpoints[0]     = v0;       endpoints[1]    = v1; 
		  incident_tri[0]  =  i;       incident_tri[1] = -1;
		}
	      else
		{ 
		  endpoints[0]     = v1;       endpoints[1]    = v0; 
		  incident_tri[0]  = -1;       incident_tri[1] = i;
		}  

	      edge_in_table = (EdgeS *)_dxfFindRecord(e_table,(KEYTYPE)the_edge);

	      /* if the edge is not already in the table we add it */

	      if (!edge_in_table)
		{
		  /* copy this edge to the edge array */
		  
		  memcpy(edge_array+nE, the_edge, sizeof(EdgeS));

		  /* add the edge to the hash table */

		  if (!_dxfAddRecord(e_table, (RECTYPE)(edge_array+nE), (KEYTYPE)(edge_array+nE)))
		    goto error;

		  /* increment the number of edges */

		  nE++;

		}

	      else
		{ /* there is already an edge (v0,v1) in the table */

		  /* either one of the following two cases apply:

		     CASE ONE: we have already determined that the edge is not a manifold edge,

		     by giving it a label of zero. In that case, 
		     we do not even add the edge to the list of edges, 
		     but we put a label of zero in the edge slot.

		     and we increment the edge counter by one 

		     */

		  if (SurfEdgeLabel(edge_in_table) == 0)
		    {
		      
		      dx_bzero(edge_array+nE, sizeof(EdgeS));

		      nE++;

		    }

		  /* CASE TWO: the edge is supposed to be a manifold edge with a label == 2
		     
		   *but now* there is a third triangle in the picture.
		     
		   so the edge cannot be manifold. 
		   so we assign it a label of zero.

		   the problem is, we thought that we had a single edge, and we have
		   two, plus the new one is three. So we increment the edge counter by two.
		   */

		  else if (SurfEdgeLabel(edge_in_table) == 2)
		    {
		      /* again none of the three edges is a manifold edge */

		      SurfEdgeLabel(edge_in_table) = 0;

		      dx_bzero(edge_array+nE, sizeof(EdgeS)); nE++;

		      dx_bzero(edge_array+nE, sizeof(EdgeS)); nE++;

		      total_nE++; /*reinstate the edge that we thought was paired */
		    }
		  
		  /* CASE THREE: we did not establish whether  the edge is a manifold edge :

		     label == 1

		     but we have two triangles and we can establish it */

		  else if (SurfEdgeLabel(edge_in_table) == 1)
		    {
		      first_tri = i;

		      triangles_in_table = SurfEdgeGetTriangles(edge_in_table);
			    
		      /* if the orientations of the faces do match,
			 we update the edge record to recognize that it is a valid edge,
			 and we will *later*  join the corresponding vertices using the father array,
			 if it is established in the edge that the edge (v0,v1) *is indeed* manifold */
			  
		      if ((v0>=v1) && triangles_in_table[0] == -1)
			{

			  /* do the triangle orientations match ? */

			  /* so v0 and v1 are visited in that order in the_triangle.
			     they have to be visited in the opposite order in triangle_in_table */

			  second_tri = triangles_in_table[1];
			      
			  /* we know that v1 appears first in second_tri (actually I am not sure what
			     this means because it is only true modulo circular permutation. 

			     We can locate v1 in t1: */
			      
			  j = 0; while (j<3 && t[second_tri*3+j] != v1) j++;

			  /* for the orientation to match, v0 must appear next, right after v1
			      
			     also, we do not want the two triangles to have exactly the same three
			     vertices */
			  
			  if (j < 3 && t[second_tri*3+(j+1)%3] == v0 && 
			      t[second_tri*3+(j+2)%3] != t[first_tri*3+(k+2)%3])
			    {
			      edge_paired_as_manifold = 1;
			      triangles_in_table[0]   = first_tri;
			    }

			}

		      else if ((v0<v1) && triangles_in_table[1] == -1)
			{
			  second_tri = triangles_in_table[0];

			  j = 0; while (j<3 && t[second_tri*3+j] != v1) j++;

			  if (j < 3 && t[second_tri*3+(j+1)%3] == v0 && 
			      t[second_tri*3+(j+2)%3] != t[first_tri*3+(k+2)%3]) 
			    {
			      edge_paired_as_manifold = 1;
			      triangles_in_table[1]   = first_tri;
			    }

			}
		      
		      if (edge_paired_as_manifold)
			{
			  /* so it worked, and we indeed have a manifold edge (for the moment)

			     we set the edge label to 2 (manifold edge) in the table */

			  SurfEdgeLabel(edge_in_table) = 2;

			  total_nE--; /* substract the edge that is paired */
			}
		      else
			{
			  /* so there are (at least) two triangles sharing that edge, but they
			     are not matchable. so we know that the edge is not manifold

			     we set the edge label accordingly
			     */

			  SurfEdgeLabel(edge_in_table) = 0;
			  
			  /* introduce a dummy edge in the edge array */

			  dx_bzero(edge_array+nE, sizeof(EdgeS)); 

			  /* increment the number of edges */

			  nE ++;
			}
		    } /* end of CASE 3 */
		}
	    }
	}

      /* Now the manifold edges have been identified.
	 For each manifold edge, I join the triangle corners associated with this edge */

      /* so I loop on the manifold edges */

      for (cur_edge = edge_array, i=0 ; i<nE; i++, cur_edge++)
	{
	  /* if the edge was established to be a manifold edge */

	  if (SurfEdgeLabel(cur_edge) == 2)
	    {
	      
	      /* retrieve v0, v1, t0, t1 */

	      incident_tri = SurfEdgeGetTriangles(cur_edge);
	      endpoints    = SurfEdgeGetVertices(cur_edge);


	      v0 = endpoints[0];

	      v1 = endpoints[1];

	      first_tri  = incident_tri[0];
	      second_tri = incident_tri[1];
	      
	      /* what are the corners associated with the first triangle incident_tri[0] ? */

	      /* ANSWER we loop on k until v0 is found */

	      k = 0 ; while (k<3 && t[first_tri*3+k] != v0) k++;

	      if (k < 3)
		{
		  corner_v0first = first_tri * 3 + k;
		  corner_v1first = first_tri * 3 + (k+1)%3;
			      
		  /* corners associated with the second tri ? */

		  /* ANSWER we loop on j until v1 is found */

		  
		  j = 0; while (j<3 && t[second_tri*3+j] != v1) j++;

		  if (j < 3)
		    {
		      /* note that we have already verified that first_tri and second_tri 
			 have different vertices when we established that the edge cur_edge was
			 manifold, so there is no need to check for this again */


		      corner_v1second = second_tri * 3 + j;
		      corner_v0second = second_tri * 3 + (j+1)%3;
			    
		      /* join the two vertex pairs */
		      _dxfJoin(corner_v1second, corner_v1first, fathers);
		      _dxfJoin(corner_v0second, corner_v0first, fathers);
		    }
		}
	    }
	}
      
      /* the hash table is not needed once the corners have been joined 

	 because we are going to allocate many more arrays,
	 it is better to free the space now */

      _dxfFreeHashTable(e_table);
      dx_bzero((char *)e_table, sizeof (Htable));

      /* return the array of edges. 
	 Edges need to be recomputed later in _dxfIndexEdges
	 we reallocate the array of edges at that time 
	 */
      
      *edges   = edge_array;

      /* return the number of edges */
      
      *n_edges = nE;

      success = 1;
    }

  return success;

error:

  if (edge_array) DXFree((Pointer)edge_array);
  _dxfFreeHashTable(e_table);
		 
  return 0;

}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfBuildOrientedManifold                                              |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfBuildOrientedManifold(int *nV, float **v, int nT, int *t, int **vlut, 
			      int *nE, EdgeS **edges, Htable *e_table, int *fathers)
{
 
  int
    success = 0;

  int
    nT3    = 3 * nT,
    new_nV = 0;

  float 
    *new_v = NULL;

  int 
    *new_t     = NULL,
    *old_vlut  = *vlut,
    *new_vlut  = NULL;

  register int i;

  /* create a conversion table from old vertices to new vertices (triangle corners)
     this is exactly the new set of triangles */


  new_t = (int *)   DXAllocateZero (nT3 * sizeof (int));

  if (!new_t) goto error;

  for (new_nV = i = 0; i < nT3; i++)
    {
      int corner_father;

      /* perform path compression of the fathers at the same time */
		  
      if ((corner_father = _dxfFather(i,fathers)) == i) 
	{
	  new_t[i] = new_nV;
	  new_nV ++;
	} 

      /* this works because the father always has lower index than the child */
      
      else new_t[i] = new_t[corner_father];

      /* corner_father should always be <= i */
    }

  /* new_nV determines the size of the new vertex array and new vlut array */

  new_vlut = (int *)   DXAllocateZero (new_nV * sizeof (int));

  if (!new_vlut) goto error;

  /* compute the new vertex look-up table */

  /* loop through new_t. For each different vertex, look-up the old
     vertex, and the look-up table of the old vertex */

  for (i = 0; i < nT3; i++) if (i == fathers[i])

    /* for each corner that is a corner group representative,
       we have already assigned a new consecutive number.
       We now update the entry in the new vertex look-up table by
       retrieving the old vertex lookup entry of the vertex corresponding 
       to the corner */

    new_vlut[new_t[i]] = old_vlut[t[i]];

  /* we discard the old look-up table and replace the old look-up table
     with the new lookup table */

  if (*vlut) 
    {DXFree((Pointer)(*vlut));}
 
  *vlut = new_vlut; 

  /* be careful, here, in case there is an error, new_vlut will be freed */


  /* reshuffle the new vertex array */

  new_v    = (float *) DXAllocateZero (new_nV * sizeof (Vertex));

  if (!new_v) goto error;

  /* copy the vertices only if the corner is a father: only for each actual new vertex */

  for (i = 0; i< nT3; i++)
    
    if (i == fathers[i])
      
      memcpy(new_v + 3 * new_t[i], *v + 3 * t[i], sizeof (Vertex));
    

  /* discard the old vertex array and replace with the new vertex array */

  if (*v) {DXFree((Pointer)(*v));}
   

  *v  = new_v;

  /* be careful, here, in case there is an error, new_v will be freed 
     except that there cannot be an error after this point */

  *nV = new_nV;


  /* we overwrite the old triangle with the new triangles
     and free the new triangles */

  memcpy(t, new_t, nT3 * sizeof (int));

  DXFree((Pointer) new_t);

  new_t = NULL;
  
  /* refill the edge hash table with the new vertex indices */

  /* to start, reinitialize the hash table */
  dx_bzero((char *)e_table, sizeof (Htable));

  if (_dxfIndexEdges(nE, edges, e_table, nT, t))
    /* if this fails, we do not go to error, because there is nothing to
       free (the hash table is freed in the error: of _dxfIndexEdges and so are the edges) 
       *v, *vlut ... will be freed later since they have been properly allocated 
       */
    success = 1;

  return success;

error:
  if (new_t)    {DXFree((Pointer)new_t);}
  
  if (new_vlut) 
    {
      DXFree((Pointer)new_vlut); 
      /* new_vlut was computed and vlut freed and replaced with new_vlut,
	 so we must set *vlut = NULL */
      *vlut = NULL;}
 
  if (new_v)    
    {
      DXFree((Pointer)new_v); 
      /* new_v was computed and *v freed and replaced with new_v,
	 so we must set *v= NULL */
      *v = NULL;
    }
  
  return 0;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfIndexEdges                                                         |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfIndexEdges(int *n_edges, EdgeS **new_edges, Htable *e_table, int nT, int *t)
{
 /*
    It is very important to count the edges again!
     
     particularly, we may have implicitely joined edges that we thought were not joined */

  EdgeS 
    *edges = NULL,
    the_edge[1], 
    *edge_in_table;
	     
  int
    nE                   = *n_edges,
    *the_triangle        = t,
    *incident_tri        = SurfEdgeGetTriangles(the_edge),
    *endpoints           = SurfEdgeGetVertices (the_edge),
    second_tri,
    v0,                    /* surface edge end-points */
    v1,
    *triangles_in_table;

   
  register int i;

  register u_char k;

  /* this routine works for both a manifold and a non manifold surface
     but the number of edges must be known in advance */

  dx_bzero((char *)e_table, sizeof(Htable)); /* reset the hash table */


  /* initialize the number of faces per edge to 1 */
       
  SurfEdgeLabel(the_edge) = 1; 

  /* reallocate the array of edges */

  /* there may be fewer edges than nE, especially if some edges were
     implicitely merged. However, there cannot be more than nE edges
     So we are safe. However, it is *critical* to *not reallocate* after
     the final hash_table is build in this procedure. Otherwise,
     the edges may be relocated, and the hash-table entries will point
     to NIL. */

  *new_edges =  (EdgeS *) DXReAllocate(*new_edges, nE * sizeof (EdgeS));

  if (!(*new_edges)) goto error;

  edges = *new_edges;

  /* initialize the hash table with roughly as much entries as
     the number of edges */

   if (!_dxfInitHashTable(e_table, 
			     _dxfFlog2(nE), 
			     (int (*)(KEYTYPE))_dxfEdgeSKey, 
			     (int (*)(RECTYPE, RECTYPE))_dxfEdgeSFilter, 
			     NULL)) goto error;
  
    for (the_triangle = t, nE = i = 0; i < nT; i++, the_triangle += 3)
	{
	  
	  for (k=0;k<3;k++)
		    
	    {
	      v0 = the_triangle[k];        /* surface edge end-points */
	      v1 = the_triangle[(k+1)%3];

   
	      if (v0>=v1) 
		{ 
		  endpoints[0]     = v0;           endpoints[1] = v1; 
		  incident_tri[0]  =  i;        incident_tri[1] = -1;
		  second_tri       =  0;
		}
	      else
		{ 
		  endpoints[0]     = v1;           endpoints[1] = v0; 
		  incident_tri[0]  = -1;        incident_tri[1] = i;
		  second_tri       =  1;
		}  

	      edge_in_table = (EdgeS *)_dxfFindRecord(e_table,(KEYTYPE)the_edge);

	      if (edge_in_table)
		{
		  
		  triangles_in_table = SurfEdgeGetTriangles(edge_in_table);
			
		  SurfEdgeLabel(edge_in_table) += 1;

		  triangles_in_table[second_tri] = i;

		}
	      else
		{
		  
		  memcpy(edges+nE, the_edge, sizeof(EdgeS));


		  /* add the edge to the hash table */

		  if (!_dxfAddRecord(e_table, (RECTYPE)(edges+nE), (KEYTYPE)(edges+nE)))
		    goto error;

		  /* increment the number of edges */
		  nE++;

		}
	    }
	}

    /* now the final number of edges is known. 
	I return this number to the caller routine */

    *n_edges     = nE;


  return 1;

error:
  
  _dxfFreeHashTable(e_table);

  if (*new_edges) DXFree((Pointer)(*new_edges)); *new_edges = NULL;
	
  return 0;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfEliminateDegenerateTriangles                                       |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfEliminateDegenerateTriangles(int nV, int nT, int *t, int *new_nT, int **new_t, int **flut)
{
  int success = 0; /* variable converted to one and returned if the procedure executes
                      successfully */
  
  int 
    trisize        = 3 * sizeof (int),
    *t1            = t,
    *new_tris      = NULL,
    *face_lut      = NULL;

  register int i;

  /* initially, set the new number of triangles to zero */

  *new_nT = 0;

  /* 4/20/97: with respect to the original implementation, we now have the difference
     that there may be a pre-existing "flut" face look-up table, which was created
     for instance when quads were converted to triangles */


  if (*flut == NULL)
    {
      face_lut = *flut = (int *) DXAllocateZero(nT * sizeof (int));
      
      if (!(*flut)) goto error;

      /* initialize the face_lut to point to the identity */
      
      for (i=0;i<nT;i++) face_lut[i] = i;
    }
  else face_lut  = *flut;


  if (!(*flut)) goto error;

  /* allocate more than enough new triangles and space for the new "flut": 
     we can't add any triangle above nT but we may eliminate some */

  new_tris   = *new_t = (int *) DXAllocateZero(nT * trisize);

  if (new_tris && (*new_t))
    {
      
      /* take each old triangle in turn */

      for (i=0 ; i<nT; i++, t1+=3)

	/* if all three triangle indices are different, add the triangle to the new triangle list */
	
	/* of course, all the indices should be  0 <=   < nV

	   (we never know what sort of strange input people could create) */

	if (t1[0] >= 0 && t1[0] < nV && t1[1] >= 0 && t1[1] < nV && t1[2] >= 0 && t1[2] < nV)
	
	  if (t1[0] != t1[1] && t1[1] != t1[2] && t1[0] != t1[2])
	    {
	      if (*new_nT < i)
		face_lut[*new_nT] = face_lut[i]; /* i is always larger (>=) than *new_nT
						    so we will never overwrite a face_lut
						    entry that we will later need */

	      /* copy this triangle to the new list */

	      memcpy(new_tris + 3 * (*new_nT), t1, trisize);
	    
	      *new_nT +=1;
	    }

      /* now reallocate the new triangle array using the new number of triangles 

       reallocate also the face look up table */

      
      if (nT != *new_nT)
	{
	  
	  *new_t = (int *) DXReAllocate(*new_t,  *new_nT * trisize);
	  
	  if (!(*new_t)) goto error;

	  *flut  = (int *) DXReAllocate(*flut,    *new_nT * sizeof (int));

	  if (!(*flut)) goto error;

	}

      /* don't forget to return "success" if there was no problem */

      success = 1;
    }

  else goto error;
  

  return success;

error:

  if (*new_t)   {DXFree((Pointer) (*new_t)); *new_t = NULL;}
  if (*flut)    {DXFree((Pointer) (*flut));  *flut  = NULL;}

  return 0;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfEliminateStandaloneVertices                                        |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfEliminateStandaloneVertices(int nV, float *v, int *new_nV, float **new_v_exported, 
				    int nT, int *t, int **vlut)
{
  /* the *new_v_exported notation stands for "Manifold set of Vertices" */

  register int i, i1;
  register u_char k;

  float *new_v = NULL;

  int   
    *new_tris                = t,
    *vertex_lut              = NULL,
    *old_to_new_vertex_index = NULL;

  /* allocate the new arrays using the old counts since we potentially
     eliminate standalone vertices, the arrays cannot grow larger
     during this procedure (the vertex array may grow larger later on */

  int *valence = (int *) DXAllocateZero(nV * sizeof (int));

  if (!valence) goto error;

  old_to_new_vertex_index = (int *) DXAllocateZero(nV * sizeof (int));

  if (!old_to_new_vertex_index) goto error;

  new_v = *new_v_exported = (float *) DXAllocateZero(3 * nV * sizeof (float));

  if (!(*new_v_exported)) goto error;

  vertex_lut = *vlut = (int *) DXAllocateZero(nV *sizeof (int));

  if (!(*vlut)) goto error;

  /* using the set of triangles, compute the valence of each vertex */

  _dxfComputeVertexValence(nT, t, valence);
  
  *new_nV = 0; /* reset the number of new vertices */
  
  /* build the look-up table of vertices with non-zero valence: the
     remaining vertices */

  
  for (i=i1=0; i<nV; i++, i1+=3)

    if (valence[i]>0)
      {
	vertex_lut[*new_nV] = i;
	old_to_new_vertex_index[i]     = *new_nV;

	/* copy the vertices that are not standalone in the new array */
            
	memcpy(new_v + (*new_nV)*3, v + i1, 3 * sizeof (float));

	*new_nV += 1; /* increment the count of new vertices */
      }

  /* change the vertex indices in the triangle array */
  
  for (i=i1=0; i<nT; i++, i1+=3)
                  
    /* for each triangle, change the index to its three vertices */
               
    for (k=0; k<3; k++) new_tris[i1+k] = old_to_new_vertex_index[t[i1+k]];
             

  /* if there are fewer vertices than before, we reallocate the vertex array and the vertex lut */

  if (*new_nV != nV)
    {
      *new_v_exported   = (float *)  DXReAllocate (*new_v_exported, *new_nV * sizeof (Vertex));

      if (!(*new_v_exported))   goto error;

      *vlut = (int *)    DXReAllocate(*vlut, *new_nV * sizeof (int));

      if (!(*vlut)) goto error;
    }

  DXFree((Pointer) old_to_new_vertex_index);
  DXFree((Pointer) valence);
  
  return 1;

error:

  if (*new_v_exported)         DXFree((Pointer) (*new_v_exported));    *new_v_exported    = NULL;
  if (*vlut)                   DXFree((Pointer) (*vlut));              *vlut              = NULL;
  if (valence)                 DXFree((Pointer) valence);
  
  if (old_to_new_vertex_index) DXFree((Pointer) old_to_new_vertex_index);

  return 0;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   dxfInitHashTable                                                       |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfInitHashTable(Htable *table, int numbits, int (*fcode)(KEYTYPE), 
		      int (*filter)( RECTYPE, RECTYPE), KEYTYPE (*compute_key)(RECTYPE))
{
  register u_int i;

  /* relation between size and number of bits */
  table->size = 1 << numbits;

  /* the codes will be used as hash entries modulo table->mask */
  table->mask = table->size-1;

  /* alloc the table of entries */

  table->bucket = (Nodetype **) DXAllocateZero (sizeof(Nodetype *) * table->size);

  if (!table->bucket) goto error;

  for(i=0;i<table->size;i++)
    table->bucket[i] = (Nodetype *) NULL;

  table->fcode  = fcode;  /* function that computes the code from the key */
  table->filter = filter;/* filtering function:
	returns 1 if the filter accepts the record
	returns 0 if the record is rejected, in which case the next
	      record must be found and processed */

  
  if (compute_key == NULL) /* function that computes the key from the record */
    table->compute_key = _dxfIdentifyRecord2Key;
  else 
    table->compute_key = compute_key;

  /* initialize the first block of entries */

  table->nodes_per_block = 1 << 10 ;

  table->block    = (Block *) NULL;

  /* initialize the free node */

  table->freenode = (Nodetype *) NULL;

  return 1;

error:

  dx_bzero(table, sizeof (Htable));
  return 0;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfNewBlock                                                           |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfNewBlock(Htable *table)
{
  Block    *newblock = NULL;

  newblock = (Block *) DXAllocateZero (sizeof (Block));

  if (!newblock) goto error;

  newblock->next = table->block;
  newblock->rec  = NULL;

  table->block = newblock;

  return 1;

error:

  return 0;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   *_dxfgetnode                                                           |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

Nodetype *_dxfgetnode(Htable *table)
{
  Nodetype *retval;

  register u_int i;

 if (table->freenode)
   {
     retval          = table->freenode;
     table->freenode = table->freenode->next;
   }
  else
    {
      /* create a new block */

      /* this can potentially fail if there is not enough memory */

      if (!_dxfNewBlock(table)) goto error;

      /* allocate a new array of nodes */

      table->block->rec= (RECTYPE) DXAllocateZero (sizeof(Nodetype)
						   * table->nodes_per_block);

      if (!table->block->rec) goto error;

      /* link the nodes newly created  */

      for(retval= (Nodetype *) table->block->rec,
	  i = 0; i < table->nodes_per_block-1; i++, retval++)
	{
	  retval->next = retval+1 ;
	}

      /* the last node has to point to NULL */

      retval->next    = (Nodetype *) NULL;

      retval = (Nodetype *) table->block->rec ;
      table->freenode = retval+1;
    }

  return(retval);

error:

  return (NULL);
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   *_dxfnewnode                                                           |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

Nodetype *_dxfnewnode(Htable *table, RECTYPE rec, int code)
{
  
  Nodetype *retval = _dxfgetnode(table);

  /* if _dxfgetnode(table) returns a NULL pointer, we should not try to access it */

  if (retval)
    {
      retval->code     = code;
      retval->rec      = rec;
    }
  return(retval);
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfHashTableReset                                                     |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

void _dxfHashTableReset(Htable *table)
{
  /* resets the hash table by cleaning all buckets and
     deleting all records */

  Block 
    *block=table->block,
    *next;

  while(block)
    {
      next=block->next;
      if (block->rec) DXFree((Pointer) block->rec);
      DXFree((Pointer) block);
      block=next;
    }
  table->block = NULL;

  /* clear  buckets */

  if (table->size > 0 && table->bucket) 
    dx_bzero((char *)table->bucket, table->size * sizeof(Nodetype *)) ;

  return;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfAddRecord                                                          |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfAddRecord(Htable *table, RECTYPE rec, KEYTYPE key)
{
  Nodetype  *node = NULL;

  register int i;

  int code= table->fcode(key);

  node=_dxfnewnode(table,rec,code);

  if (!node) goto error;

  /* add the node containing the record to the appropriate bucket */

  i=node->code & table->mask;

  node->next = table->bucket[i];
  table->bucket[i]=node;

  return 1;

error:

  return 0;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   *_dxfFindRecord                                                        |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

char *_dxfFindRecord(Htable *table, KEYTYPE key)
{
  register int i;
  
  table->search_code = table->fcode (key);

  i = table->search_code & table->mask;

  table->search_key = key;

  table->search_node= table->bucket[i];

  return (_dxfNextRecord(table));
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   *_dxfNextRecord                                                        |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

char *_dxfNextRecord(Htable *table)
{
  u_char found=0;

  Nodetype  *node=table->search_node;

  /* if a filter function was provided, the first node with the same code AND such
     that the filter function comparing the node's record and the search key "agrees"
     is selected */

  if (table->filter)
    {
      while (node && !found)
      {
	if (node->code == table->search_code &&
	    table->filter(node->rec,table->compute_key(table->search_key)))
	  found=1;
	else node = node->next;
      }
    }

  /* if no filtering function was provided, the next node with the same code is
     retrieved */


  else
    while (node && !found)
      {
	if (node->code == table->search_code)
	  found=1;
	else node = node->next;
      }

  if (found)
    {
      table->search_node = node->next;
      return ((char *)node->rec);
    }

  else return ((char *) NULL);

}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfIdentifyRecord2Key                                                 |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

/* when a filtering function IS provided, we need to compare a current record with
   the search key. In the general case, a "compute_key" function must be provided
   to convert the key to a record or something that can be directly compared with
   a record
   in many cases the search key is the same as a record, or can be used without conversion.
   This is when the "dummy" function _dxfIdentifyRecord2Key is used.
   */

KEYTYPE _dxfIdentifyRecord2Key(RECTYPE rec)
{
  return (KEYTYPE) rec;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfFlog2                                                              |
  |                                                                          |
  +--------------------------------------------------------------------------+*/
/* 
   Log base 2 (integer form)
   
   The following relation holds exactly

   2 ^ Log2(N) <= N < 2 * 2 ^ Log2(N)

   if N==0 then the result is by definition -1, but the true result
   is of course -infty
*/

int _dxfFlog2(int n)
{
   int d=0;
   if (n==0) return (-1);
   while(n>1)
        {d++;
         n/=2;
        }
   return(d);
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfEdgeSKey                                                           |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

/* compute the Hash-Table code from the key in the case of a surface edge */

int _dxfEdgeSKey(EdgeS *e)
{
  return(e->v[1]^(e->v[0] << 8));
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfEdgeSFilter                                                        |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfEdgeSFilter(EdgeS *e1, EdgeS *e2)
{
  if ((e1->v[0] == e2->v[0]) &&
      (e1->v[1] == e2->v[1])) return(1);
  else
    return(0);
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfEdgeSFilter                                                        |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

EdgeS *_dxfFindEdge(Htable *e, int vA, int vB)
{
  EdgeS edge;
  
  if (vA > vB)
    { 
      edge.v[0]=vA;edge.v[1]=vB;
    }
  else
    { 
      edge.v[0]=vB;edge.v[1]=vA;
    }

  return (EdgeS *) _dxfFindRecord(e,(KEYTYPE) &edge);
 
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfInitArray                                                          |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

void _dxfInitArray(char *array, char *data, int n, int size)
{
  
  register int i;

  for(i=0;i<n;i++, array += size) 
    memcpy(array, data, size);

  return;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfInitFather                                                         |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

void _dxfInitFather(int n, int *father)
{
  register int i;

  for(i=0;i<n;i++) father[i] = i; 
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfFather                                                             |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfFather(int i, int *father)
{
  static int Ci, Fi;

  /*--- path traversal ---*/
  for(Ci=i;(Fi=father[Ci])!=Ci;Ci=Fi);

  /*--- path compression ---*/
  for(;(Fi=father[i])!=Ci;i=Fi) father[i] = Ci;

  return(Ci);
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfJoin                                                               |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfJoin(int i, int j, int *father)
{
  int rep_i, rep_j;

  if((rep_i = _dxfFather(i,father)) != (rep_j = _dxfFather(j,father)))
    {
      /* the smallest representative becomes the father of
	 the other representative */

      if (rep_i < rep_j)

	return (father[rep_j] = rep_i);

      else if (rep_j < rep_i)

	return (father[rep_i] = rep_j);

    }

  return -1;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfCreateSimplificationDataStructure                                  |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfCreateSimplificationDataStructure(SimpData *simp_data, float *vertex_data, float *old_positional_error)
{
  /* initialize all pointers except vx_data and err_volume to NULL */

  simp_data->vx_data                  = vertex_data;
  simp_data->err_volume               = old_positional_error;

  simp_data->edge2index               = NULL;
  simp_data->index2edge               = NULL;
  simp_data->edge_father              = NULL;
  simp_data->normal                   = NULL;
  simp_data->area                     = NULL;
  simp_data->compactness              = NULL;
  simp_data->valence                  = NULL;
  simp_data->boundary_vert            = NULL;
  simp_data->tol_volume               = NULL;
  simp_data->vx_data_error            = NULL;
  simp_data->vx_data_potential_values = NULL;
  simp_data->edge_heap                = NULL;


  simp_data->edge2index = (int *) DXAllocateZero(simp_data->nE * sizeof (int));

  if (!simp_data->edge2index) goto error;

  simp_data->index2edge = (int *) DXAllocateZero(simp_data->nE * sizeof (int));
  
  if (!simp_data->index2edge) goto error;

  simp_data->edge_father = (int *) DXAllocateZero(simp_data->nE * sizeof (int));
  
  if (!simp_data->edge_father) goto error;

  /* initialize the father arrays */

  _dxfInitFather(simp_data->nV, simp_data->vertex_father);

  _dxfInitFather(simp_data->nE, simp_data->edge_father);

  _dxfInitFather(simp_data->nT, simp_data->triangle_father);


  simp_data->normal = (Vertex *) DXAllocateZero(simp_data->nT * sizeof (Vertex));

  if (!simp_data->normal) goto error;

  simp_data->area   = (float  *) DXAllocateZero(simp_data->nT * sizeof (float));

  if (!simp_data->area) goto error;

  simp_data->compactness   = (float  *) DXAllocateZero(simp_data->nT * sizeof (float));

  if (!simp_data->compactness) goto error;

  /* compute the triangle normals, areas and compactness */

  _dxfTrianglesNormalAreaCompactness(simp_data->nT, (Face *)simp_data->triangles, simp_data->vert,
				     simp_data->normal, simp_data->area, simp_data->compactness);

  /* save the areas in the "old_face_areas" array */
					 
  memcpy(simp_data->old_face_areas, simp_data->area, simp_data->nT * sizeof (float));


  /* use default values for the maximum angle between face normals 
     and maximum ratio of compactness values */

  simp_data->min_scalprod_nor  = (float) cos ((double)SIMPLIFY_MAX_ANGLE_NOR);

  simp_data->compactness_ratio = SIMPLIFY_COMPACTNESS_RATIO;



  /* allocate and compute the vertex valence */

  simp_data->valence = (u_short *) DXAllocateZero(simp_data -> nV * sizeof (u_short));

  if (!simp_data->valence) goto error;

  _dxfComputeVertexValence(simp_data->nT, (int *)simp_data->triangles, simp_data->valence);
  

  /* allocate and compute the "boundary_vert" array that flags boundary vertices with a flag value of 1 */

  simp_data->boundary_vert = (u_char *)  DXAllocateZero(simp_data -> nV * sizeof (u_char));

  if (!simp_data->boundary_vert) goto error;

  _dxfFlagBoundaryVertices(simp_data->nE, simp_data->edge_array, simp_data->boundary_vert);

  
  /* allocate and initialize the positional error */

  /* there may or may not be positional error defined */

  if (simp_data->err_volume)
    {
      /* there is positional error already defined */

      float *old_error_volume = simp_data->err_volume;

      simp_data->err_volume = (float *) DXAllocateZero(simp_data->nV * sizeof (float));

      if (!simp_data->err_volume) goto error;

      memcpy(simp_data->err_volume, old_error_volume, simp_data->nV * sizeof (float));

    }
  else
    {
      simp_data->err_volume = (float *) DXAllocateZero(simp_data->nV * sizeof (float));
 
      if (!simp_data->err_volume) goto error;
   }

  simp_data->tol_volume = (float *) DXAllocateZero(simp_data->nV * sizeof (float));
  
  if (!simp_data->tol_volume) goto error;

  /* otherwise initialize the tolerance array */

   _dxfInitArray((char *)simp_data->tol_volume,(char *) &simp_data->tolerance, simp_data->nV, sizeof (float));


  /* treat the vertex-related data */

  /* there may or may not be vertex related data defined */
  
  if (simp_data->vx_data && simp_data->data_dim > 0)
    {
      register int i;

      float *old_vertex_data = simp_data->vx_data;

      simp_data->vx_data = (float *) DXAllocateZero(simp_data->nV * simp_data->data_dim * sizeof (float));

      if (!simp_data->vx_data) goto error;

      simp_data->vx_data_error = (float *) DXAllocateZero(simp_data->nV * sizeof (float));

      simp_data->vx_data_potential_values = (float *) DXAllocateZero(3 * simp_data->data_dim * sizeof (float));

      if (!simp_data->vx_data_potential_values) goto error;

      simp_data->vx_old_data = simp_data->vx_data_potential_values + simp_data->data_dim;
      simp_data->vx_new_data = simp_data->vx_old_data              + simp_data->data_dim;

      /* copy the vertex_data using the vertex_lut */

      for (i=0;i<simp_data->nV;i++)

	memcpy(simp_data->vx_data + i * simp_data->data_dim, 
	       old_vertex_data + simp_data->vertex_lut[i] * simp_data->data_dim, 
	       simp_data->data_dim* sizeof (float));

    }
  else simp_data->vx_data = NULL;

  /* allocate and initialize the edge heap */

  simp_data->edge_heap = _dxfNewHeap(simp_data->nE);

  if (!simp_data->edge_heap) goto error;

  simp_data->get_lowest_weight_edge = _dxfRemoveLowestWeightEdgeFromHeap;
  simp_data->add_edge               = _dxfAddEdge2Heap;
  simp_data->remove_edge            = _dxfRemoveEdgeFromHeap;

  simp_data->heap_initial_index     = HeapInitialIndex(simp_data->edge_heap);
  simp_data->heap_outside_index     = HeapOutsideIndex(simp_data->edge_heap);
      
  /* initialize the heap indices */
  
  _dxfInitArray((char *)simp_data->edge2index,(char *) &simp_data->heap_initial_index, 
		simp_data->nE, sizeof (int));
 
  /* initialize various counters: */
  simp_data->num_edg_collapsed           = 0;
  simp_data->num_edg_weights             = 0;
  simp_data->edg_rejected_4_topology     = 0;
  simp_data->edg_rejected_4_geometry     = 0;
  simp_data->edg_rejected_4_tolerance    = 0;

  /* initialize the last edge added to the heap: */
  simp_data->last_edge_added2heap        = 0;

  /* initialize the number of edges left to be tested: */
  simp_data->num_edg_remaining_4_testing = simp_data->nE;

  simp_data->valence_max = SIMPLIFY_VALENCE_MAX;
  
  return 1;

error:

  if (simp_data->edge2index)               
    {
      DXFree((Pointer) simp_data->edge2index);
      simp_data->edge2index = NULL;
    }
 
  if (simp_data->index2edge)
    {
      DXFree((Pointer) simp_data->index2edge);
      simp_data->index2edge = NULL;
    }

  if (simp_data->edge_father) 
    {
      DXFree((Pointer) simp_data->edge_father);
      simp_data->edge_father = NULL;
    }

  if (simp_data->normal) 
    {
      DXFree((Pointer) simp_data->normal);
      simp_data->normal = NULL;
    }

  if (simp_data->area)
    {
      DXFree((Pointer) simp_data->area);
      simp_data->area = NULL;
    }

  if (simp_data->compactness) 
    {
      DXFree((Pointer) simp_data->compactness);
      simp_data->compactness = NULL;
    }
 
  if (simp_data->valence)   
    {
      DXFree((Pointer) simp_data->valence);
      simp_data->valence = NULL;
    }

  if (simp_data->boundary_vert)  
    {
      DXFree((Pointer) simp_data->boundary_vert);
      simp_data->boundary_vert = NULL;
    }

  if (simp_data->err_volume && simp_data->err_volume != old_positional_error)
    {
      DXFree((Pointer) simp_data->err_volume);
      simp_data->err_volume = old_positional_error;
    }

  if (simp_data->tol_volume) 
    {
      DXFree((Pointer) simp_data->tol_volume);
      simp_data->tol_volume = NULL;
    }

  if (simp_data->vx_data && simp_data->vx_data != vertex_data)                  
    {
      DXFree((Pointer) simp_data->vx_data);
      simp_data->vx_data = vertex_data;
    }
  
  if (simp_data->vx_data_error) 
    {
      DXFree((Pointer) simp_data->vx_data_error);
      simp_data->vx_data_error = NULL;
    }
  
  if (simp_data->vx_data_potential_values)
    { 
      DXFree((Pointer)simp_data->vx_data_potential_values); 
      simp_data->vx_data_potential_values = NULL;
    }
  
  if (simp_data->edge_heap)   
    {
      DXFree((Pointer) simp_data->edge_heap);
      simp_data->edge_heap = NULL;
    }

  return 0;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfFreeSimplificationDataStructure                                    |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfFreeSimplificationDataStructure(SimpData *simp_data, float *vertex_data, float *old_positional_error)
{

  if (simp_data->edge2index)               
    {
      DXFree((Pointer) simp_data->edge2index);
      simp_data->edge2index = NULL;
    }
 
  if (simp_data->index2edge)
    {
      DXFree((Pointer) simp_data->index2edge);
      simp_data->index2edge = NULL;
    }

  if (simp_data->edge_father) 
    {
      DXFree((Pointer) simp_data->edge_father);
      simp_data->edge_father = NULL;
    }

  if (simp_data->normal) 
    {
      DXFree((Pointer) simp_data->normal);
      simp_data->normal = NULL;
    }

  if (simp_data->area)
    {
      DXFree((Pointer) simp_data->area);
      simp_data->area = NULL;
    }

  if (simp_data->compactness) 
    {
      DXFree((Pointer) simp_data->compactness);
      simp_data->compactness = NULL;
    }
 
  if (simp_data->valence)   
    {
      DXFree((Pointer) simp_data->valence);
      simp_data->valence = NULL;
    }

  if (simp_data->boundary_vert)  
    {
      DXFree((Pointer) simp_data->boundary_vert);
      simp_data->boundary_vert = NULL;
    }

  if (simp_data->err_volume && simp_data->err_volume != old_positional_error)
    {
      DXFree((Pointer) simp_data->err_volume);
      simp_data->err_volume = old_positional_error;
    }

  if (simp_data->tol_volume) 
    {
      DXFree((Pointer) simp_data->tol_volume);
      simp_data->tol_volume = NULL;
    }

  if (simp_data->vx_data && simp_data->vx_data != vertex_data)                  
    {
      DXFree((Pointer) simp_data->vx_data);
      simp_data->vx_data = vertex_data;
    }
  
  if (simp_data->vx_data_error) 
    {
      DXFree((Pointer) simp_data->vx_data_error);
      simp_data->vx_data_error = NULL;
    }
  
  if (simp_data->vx_data_potential_values)
    { 
      DXFree((Pointer)simp_data->vx_data_potential_values); 
      simp_data->vx_data_potential_values = NULL;
    }
  
  if (simp_data->edge_heap)   
    {
      DXFree((Pointer) simp_data->edge_heap);
      simp_data->edge_heap = NULL;
    }
  
  return 1;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfPartialFreeSimplificationDataStructure                             |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfPartialFreeSimplificationDataStructure(SimpData *simp_data)
{

  if (simp_data->edge2index)               
    {
      DXFree((Pointer) simp_data->edge2index);
      simp_data->edge2index = NULL;
    }
 
  if (simp_data->index2edge)
    {
      DXFree((Pointer) simp_data->index2edge);
      simp_data->index2edge = NULL;
    }

  if (simp_data->edge_father) 
    {
      DXFree((Pointer) simp_data->edge_father);
      simp_data->edge_father = NULL;
    }

  if (simp_data->compactness) 
    {
      DXFree((Pointer) simp_data->compactness);
      simp_data->compactness = NULL;
    }
 
  if (simp_data->valence)   
    {
      DXFree((Pointer) simp_data->valence);
      simp_data->valence = NULL;
    }

  if (simp_data->boundary_vert)  
    {
      DXFree((Pointer) simp_data->boundary_vert);
      simp_data->boundary_vert = NULL;
    }

  if (simp_data->tol_volume) 
    {
      DXFree((Pointer) simp_data->tol_volume);
      simp_data->tol_volume = NULL;
    }

  if (simp_data->vx_data_error) 
    {
      DXFree((Pointer) simp_data->vx_data_error);
      simp_data->vx_data_error = NULL;
    }
  
  if (simp_data->vx_data_potential_values)
    { 
      DXFree((Pointer)simp_data->vx_data_potential_values); 
      simp_data->vx_data_potential_values = NULL;
    }
  
  if (simp_data->edge_heap)   
    {
      DXFree((Pointer) simp_data->edge_heap);
      simp_data->edge_heap = NULL;
    }
  
  return 1;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfTrianglesNormalAreaCompactness                                     |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfTrianglesNormalAreaCompactness( int nT, Face *t, Vertex *v, Vertex *t_normal,
					float *t_area, float *t_compactness)
{
  int success = 1;

  static VertexD triVert[3];
  static double n[4];

  register int i;
  register u_char j,k;

  /* compute the face normals */
   
  for(i=0;i<nT;i++,t++)
    {
      double
	sum_lengths_sq = 0.;

      for (j=0;j<3;j++)
	  for (k=0;k<3;k++)
	    triVert[j][k] = v[t[0][j]][k];
      
       _dxfTriangleNormalQR2D(triVert, n);

       for (k=0;k<3;k++)
	 t_normal[i][k] = n[k];

       t_area[i] = n[3];

       /* compute the squared lenghts of the three edges: */
  
       for(j=0;j<3;j++)
	 {
	   double len_edg_sq = 0, tmp;

	   for(k=0;k<3;k++)
	     {
	       tmp         = triVert[(j+1)%3][k] - triVert[j][k];
	       len_edg_sq += tmp * tmp;
	     }

	   sum_lengths_sq += len_edg_sq;

	 }

       t_compactness[i] = FOUR_SQRT_3 * t_area[i] / sum_lengths_sq;

    }

  return success;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _DxfTriangleNormalQR2D                                                 |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

/* computation of the triangle normal in double precision

   n is an array 

   double n[4]

   in n[3] the triangle area is returned 

   
   this method is much more numerically stable than the cross product of the vectors.

   it uses Householder orthogonalization of a set of two basis vectors chosen for the
   triangle. Experimental studies show that it is better to take both the shortest
   and longest edges of the triangle to form the basis

   */

/*  Turn off set but never used warnings for ("origin" below)  */
#ifdef sgi
#  pragma set woff 1552
#endif

int _dxfTriangleNormalQR2D(VertexD *tri, double *n)
{

  static VertexD x1,x2;

  int origin = 0;
   
  _dxfTriangleBasisVectors(tri, x1, x2, origin);

  return _dxfVectorProductQRD(x1,x2,n);

}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfTriangleNormalQR2                                                  |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

/* computation of the triangle normal in single precision

   the triangle area is *not* computed */

int _dxfTriangleNormalQR2(Vertex *tri, Vertex n)
{

  static float eps = 1e-8;

  static Vertex x1,x2;
 
  int origin = 0;

  dx_bzero((char *)n, 3 * sizeof (float));

  _dxfTriangleBasisVectors(tri, x1, x2, origin);

  {
  
    float 
      sign1 = (x1[0] >= 0.) ? 1. : -1.,
      sign2,
      norm1,
      norm2;

    float beta1=0, beta;

    static Vertex v1,v2;

    register u_char k;

    norm1 = NORM3(x1);

    memcpy(v1,x1,sizeof(Vertex));

    v1[0] += sign1 *  norm1;

    /* we have to record the sign of H1 x1: -sign1 */

    if (norm1 >= eps)
      
      /* we apply H1 to the vector x2 and replace x2 by H1 x2 
	 if the column is already zero, it is not necessary to use h1 */
      {

	beta1 = - 2. / SCALPROD3(v1,v1);
	
	beta  = beta1 * SCALPROD3(v1,x2);
     
	/* x2 = x2 + beta v1 */
	
	for(k=0;k<3;k++)
	  x2[k] += beta * v1[k];
      }

    /* we create the v2 householder vector;
       it is by the way only a 2-vector */
    
    {double norm2_d = x2[1]*x2[1] + x2[2] * x2[2];
      
     if (norm2_d > 0.0) norm2 = (float) sqrt(norm2_d);
  
     else norm2 = 0.0;}

    /* we record the sign of H2 x2: -sign2 */

    sign2 = (x2[1] >= 0.) ? 1. : -1.;

    v2[0] = 0; 
    v2[1] = x2[1] + sign2 * norm2;
    v2[2] = x2[2];
      
  
    if (norm2 >=eps)

      {
	float beta2 = -2. *v2[2] / (v2[1]*v2[1] + v2[2]*v2[2]);

	n[1] =      beta2 * v2[1];
	n[2] = 1. + beta2 * v2[2];
      }


    if (norm1 >= eps)
      {

	/* we apply H1 to the H2 e3 vector, and that's the normal
	   up to a -1 factor:*/
    
	beta = beta1 * SCALPROD3(v1, n);

	/* n = n + beta v1 */
	
	for(k=0;k<3;k++)
	  n[k] += beta * v1[k];
      }

    
    if (sign1 * sign2 <0.) 

      for(k=0;k<3;k++)

	n[k] = - n[k];
    
  }


  return 1;
}

#ifdef sgi
#  pragma reset woff 1552
#endif

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfVectorProductQRD                                                   |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

/* compute the vector product using QR decomposition (which, to me, is the
   same as Householder orthogonalization. I guess the correct name is QR, since
   a Householder transformation represents a particular symmetry, and
   several Householder transformation must be combined to produce the final
   rotation + symmetry.  

   In the case of a 3 * 2 system, either one or two Householder transformations are required.

   in n[3] store the area of the triangle defined by the origin
   and the two vectors x1 and x2 */


int _dxfVectorProductQRD(VertexD x1, VertexD x2, double *n)
{  

  dx_bzero((char *)n, 4 * sizeof (double));

  {
    static double eps = 1e-15;

    /* do a QR decomposition of the matrix (x1,x2) */

    /* the first householder vector x1 is equal to : 
       x1 + sign(x11) |x1| (1 0 0)^t
       we will call H1 the corresponding Householder reflexion
       */
    

    double 
      sign1 = (x1[0] >= 0.) ? 1. : -1.,
      sign2,
      norm1,
      norm2;

    double beta1=0, beta;

    static VertexD v1,v2;

    register u_char k;

    norm1 = NORM3(x1);

    memcpy(v1,x1,sizeof(VertexD));

    v1[0] += sign1 *  norm1;

    /* we have to record the sign of H1 x1: -sign1 */

    if (norm1 >= eps)
      
      /* we apply H1 to the vector x2 and replace x2 by H1 x2 
	 if the column is already zero, it is not necessary to use h1 */
      {
	/* beta = -2/ (v1^t v1)* (v1^x2) */


	beta1 = - 2. / SCALPROD3(v1,v1);
	
	beta  = beta1 * SCALPROD3(v1,x2);
     
	/* x2 = x2 + beta v1 */
	
	for(k=0;k<3;k++)
	  x2[k] += beta * v1[k];
      }

    /* we create the v2 householder vector;
       it is by the way only a 2-vector */
    
    norm2 = x2[1]*x2[1] + x2[2] * x2[2];
      
    if (norm2 > 0.0) norm2 = sqrt(norm2);

    /* we record the sign of H2 x2: -sign2 */

    sign2 = (x2[1] >= 0.) ? 1. : -1.;

    v2[0] = 0; 
    v2[1] = x2[1] + sign2 * norm2;
    v2[2] = x2[2];
      
    /* compute the area of the triangle */

    n[3] = norm1 * norm2 /2.;
    
  
    if (norm2 >=eps)

      /* we apply H2 to the e3 vector and put the result directly in n */
      {
	/* beta2 = -2 /v2^t v2 * v2^t e3 */
	double beta2 = -2. *v2[2] / (v2[1]*v2[1] + v2[2]*v2[2]);

	n[1] =      beta2 * v2[1];
	n[2] = 1. + beta2 * v2[2];
      }


    if (norm1 >= eps)
      {

	/* we apply H1 to the H2 e3 vector, and that's the normal
	   up to a -1 factor:*/
    
	beta = beta1 * SCALPROD3(v1, n);

	/* n = n + beta v1 */
	
	for(k=0;k<3;k++)
	  n[k] += beta * v1[k];
      }

    
    /* There is a unique QR decomposition such that the diagonal of R
       is made of positive elements and Q is a rotation.
       The last vector of Q is going to be n.
       if ( R(1,1)*R(2,2) < 0 ) we have to invert n */

    if (sign1 * sign2 <0.) 

      for(k=0;k<3;k++)

	n[k] = - n[k];
    
  }

  return 1;


}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfFlagBoundaryVertices                                               |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

/* as the name indicates, mark the boundary vertices with a 
value of 1 in the array "boundary_vert"
the interior vertices are marked with 0 */

int _dxfFlagBoundaryVertices(int nE, EdgeS *edges, u_char *boundary_vert)
{
  /* loop through edges. 

     if an edge is a boundary edge, then mark both vertices as boundary vertices */

  register int i;

  for (i=0; i<nE; i++, edges++)
    
    
    if (SurfEdgeLabel(edges) == 1)
      {
	int *endpoints = SurfEdgeGetVertices(edges);

	boundary_vert[endpoints[0]] = boundary_vert[endpoints[1]] = 1;

      }
    
  return 1;
}

/*

  static functions for Heap processing 

 */

static void HeapSwitch(int i, int j, Heap *h);
static void HeapDown(int father, Heap *h);
static void HeapUp(Heap *h);
static int  HeapDel(int i, Heap *h);

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   HeapSwitch                                                             |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

static void HeapSwitch(int i, int j, Heap *h)
{
	/* switch two elements of the Heap, and update the permutation
	and inverse permutation arrays */

  if(i!=j && i<h->n && j<h->n)
    {
      float key;
      int   pos;
      
      key = h->key[i];  h->key[i]  = h->key[j];  h->key[j] = key;
      pos = h->invp[i]; h->invp[i] = h->invp[j]; h->invp[j] = pos;

      h->perm[h->invp[i]] = i;
      h->perm[h->invp[j]] = j;
    }
}


/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   HeapDown                                                               |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

static void HeapDown(int father, Heap *h)
{

	/* trickle an element down the Heap until 
	the Heap condition is verified 

	used when deleting an element from the Heap
	*/


  if(0<=father && father<h->last-1)
    {
      int left,right,child;
      
      while((left=LeftSon(father))<h->last)
	{
	  child = father;
	  if(h->key[left]<h->key[child]) child = left;
	  if((right=RightSon(father))<h->last &&
	     h->key[right]<h->key[child]) child = right;
	  if(child==father)
	    break;
	  HeapSwitch(father,child,h);
	  father = child;
	}
    }
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   HeapUp                                                                 |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

static void HeapUp(Heap *h)
{
	
/* trickle the last heap element up the Heap until the heap
	condition is verified 
	
	used when adding an element to the heap */

  int child,father;

  child = h->last-1;

  while(child>0 && h->key[(father=(child-1)/2)]>=h->key[child])
    {
      HeapSwitch(father,child,h);
      child = father;
    }
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   HeapDel                                                                |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

static int HeapDel(int i, Heap *h)
{

/* delete an element inside a Heap using its actual Heap index
   this procedure is not meant to be called outside this file */
	
  if (0<=i && i<h->last)
    {
      h->last--;
      HeapSwitch(i,h->last,h);
      HeapDown(i,h);
      return 1;
    }
  else return 0;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfHeapDelete                                                         |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int  _dxfHeapDelete(int i, Heap *h)
{
  /* delete a Heap element using the original index it was assigned
     when first added to the Heap.
	 This is the preferred procedure for deleting Heap elements
	 (together with HeapDelMin())*/

  if(0<=i && i<h->n)
    return HeapDel(h->perm[i],h);

  else return 0;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfHeapDelMin                                                         |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfHeapDelMin(Heap *h)
{
	
  /* retrieve the Heap element with the smallest key
   (first element in the Heap), take it off the Heap,
   and reorganize the Heap */

  if(h->last>0)
    {
      HeapDel(0,h);
      return(h->invp[h->last]);
    }
  else
    return(-1);			/* heap is empty */
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   *_dxfNewHeap                                                           |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

Heap *_dxfNewHeap(int n)
{

	/* create a new Heap: allocate and initialize */

  Heap *heap = (Heap*)0;

  if(n>0)
    {
      unsigned long bytes = _dxfHeapSize(n);
      char         *block = (char *) DXAllocateZero(bytes);

      if (!block) goto error;
     
      heap       = (Heap*)block;
      heap->n    = n;
      heap->key  = (float*)(block+sizeof(Heap));
      heap->perm = (int*)(block+sizeof(Heap)+n*sizeof(float));
      heap->invp = (int*)(block+sizeof(Heap)+n*sizeof(float)+n*sizeof(int));

      _dxfResetHeap(heap);
    }

  return(heap);

error:

  return NULL;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfHeapSize                                                           |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

unsigned long _dxfHeapSize(int n)
{
	/* used by NewHeap(), compute the size that is needed
	   given the maximum number of elements that can be contained
	   by the Heap */


  unsigned long bytes = sizeof(Heap)+n*sizeof(float)+2*n*sizeof(int);
      
  return(bytes);
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfResetHeap                                                          |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

void _dxfResetHeap(Heap *heap)
{

	/* without freeing the memory, re-initialize a Heap for
	   subsequent use */


  heap->last = 0; 
      
  if(heap->n>0)
    {
      register int i;
      
      for(i=0;i<heap->n;i++)
	{
	  heap->key[i]=0.0;
	  heap->perm[i] = heap->invp[i] = i;
	}
    }

  return;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfHeapAdd                                                            |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfHeapAdd(float key, Heap *h)
{
  int indx = -1;

  if(h->last<h->n) /* don't want to enlarge the heap if h->last==h->n ? */
    {
      int last = h->last++;

      indx = h->invp[last];
      h->key[last] = key;
      HeapUp(h);
    }

  return(indx);
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfRemoveEdgeFromHeap                                                 |
  |                                                                          |
  +--------------------------------------------------------------------------+*/


/* remove a surface edge from the Heap using the edge's original
   index */


int _dxfRemoveEdgeFromHeap(SimpData *simp_data, int index, int the_edge)
{
  return _dxfHeapDelete(index, simp_data->edge_heap);
} 

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfAddEdge2Heap                                                       |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

/* add a surface edge to the heap and record the index associated
   with the edge */

int _dxfAddEdge2Heap(SimpData *simp_data, int the_edge)
{
  int 
    added      = 0;


  /* verify that the Heap is not already full and that 
     the edge does not have an index that forbids adding to the Heap */


  if (!HeapFull(simp_data->edge_heap) && simp_data->edge2index[the_edge] >=0)
    {
      /* compute the representatives of the edge vertices */
      int 
	*v = SurfEdgeGetVertices(simp_data->edge_array+the_edge),
	v0 = FATHER(v[0], simp_data->vertex_father),
	v1 = FATHER(v[1], simp_data->vertex_father);

      /* 1) to compute the key to the heap, we choose the edge length to start with */
      float 
	weight = DIST3(simp_data->vert[v0], simp_data->vert[v1]);

      /* 1.1) we add the sum of the errors at the two vertices */

      weight += simp_data->err_volume[v0] + simp_data->err_volume[v1];

      /* 2) add the edge to the heap */
      
	
      simp_data->index2edge[
			    simp_data->edge2index[the_edge]=
			    _dxfHeapAdd(
					weight,
					simp_data->edge_heap)
                           ] = the_edge;
      added = 1;
	
    }

  return added;

}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfRemoveLowestWeightEdgeFromHeap                                     |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfRemoveLowestWeightEdgeFromHeap(SimpData *simp_data)
{
  int 
    idx  = _dxfHeapDelMin(simp_data->edge_heap);
	
  if (idx != -1)
  
    return simp_data->index2edge[idx];

  else
    return -1;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfBuildSimplifiedSurface                                             |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfBuildSimplifiedSurface(SimpData *simp_data, int *new_nV, float **new_v, 
			       int *new_nT, int **new_t,
			       float **new_face_areas, float **face_normals,
			       float **new_positional_error, float **new_vertex_data)

/*...........................................................................*/
/*
   
    Use the information contained in the vertex and triangle fathers
    to construct a simplified surface:
    
    only the vertices and triangles that are their own "father" are carried to
    the simplified surface. The vertex_parent and face_parent array are reshuffled
    to be used later for the purpose of interpolating data (positions-dependent
    or connections-dependent) from the original surface to the simplified surface)

    I use indifferently the word "parent" or "father" in this
    implementation. I think "parent" is better than "father", which
    refers to the most "direct" parent. But the father of the father
    is a "parent", and it sounds better to me than a
    "father". However, there was much existing code using the word
    "father" so I kept it.
    

    */
     /*...........................................................................*/

{
 
  /* count the new number of vertices */

  int 
    *vertex_parents = simp_data->vertex_father,
    *face_parents   = simp_data->triangle_father;

  register int i,j;

  float 
    *v            = NULL;

  int   
    *t            = NULL,
    *new_v_number = NULL,
    *new_t_number = NULL;
  


  /* I) build the new vertices and copy the vertex errors and data */




  new_v_number = (int *) DXAllocateZero(simp_data->nV * sizeof (int));

  if (!new_v_number) goto error;

  *new_nV = 0;

  for (i=0; i< simp_data->nV; i++) if (i == FATHER(i,vertex_parents)) 
    
    {
      new_v_number[i] = *new_nV;
      *new_nV += 1;
    }


  /* allocate the new vertex array */

  v               = (float *) DXAllocateZero (*new_nV * sizeof (Vertex));

  if (!v) goto error;

  *new_positional_error = (float *)  DXAllocateZero (*new_nV * sizeof (float));

  if (!(*new_positional_error)) goto error;
  

  if (simp_data->vx_data)
    {
      *new_vertex_data = (float *) DXAllocateZero (*new_nV * simp_data->data_dim * sizeof (float));

      if (!(*new_vertex_data)) goto error;
    }


  /* copy the vertices and other vertex attributes and update the vertex_parents array */

  
  for(j=i=0;i<simp_data->nV;i++) if (i == vertex_parents[i]) 
    {
      memcpy(v + 3 * j,simp_data->vert[i],sizeof(Vertex));

      (*new_positional_error)[j] = simp_data->err_volume[i];

      if (simp_data->vx_data) 
	memcpy(*new_vertex_data   + j * simp_data->data_dim,
	       simp_data->vx_data + i * simp_data->data_dim, simp_data->data_dim * sizeof (float));
      
      j++;
    }
  
  for(i=0;i<simp_data->nV;i++) vertex_parents[i] = new_v_number[vertex_parents[i]];
 
  DXFree((Pointer) new_v_number); new_v_number = NULL;

  
  /* copy the new vertices into *new_v and reallocate *new_v */

  
  DXFree((Pointer) (*new_v));

  *new_v = (float *) DXAllocateZero(*new_nV * sizeof (Vertex));

  if (!(*new_v)) goto error;

  memcpy(*new_v, v, *new_nV * sizeof (Vertex));

  DXFree((Pointer) v); v = NULL;




  /* II) build the new triangles and copy triangle normals and areas */




  new_t_number = (int *) DXAllocateZero(simp_data->nT * sizeof (int));

  if (!new_t_number) goto error;

  *new_nT = 0;

  for (i=0; i< simp_data->nT; i++) if (i == FATHER(i,face_parents)) 
    
    {
      new_t_number[i] = *new_nT;
      *new_nT += 1;
    }

  t               = (int *)   DXAllocateZero (*new_nT * sizeof (Face));

  if (!t) goto error;

  *new_face_areas = (float *) DXAllocateZero (*new_nT * sizeof (float));

  if (!(*new_face_areas)) goto error;

  *face_normals   = (float *) DXAllocateZero (*new_nT * sizeof (Vertex));

  if (!(*face_normals)) goto error;
  
  /* copy the triangles and other triangle attributes and update the face_parents array */
  
  for(j=i=0;i<simp_data->nT;i++) if (i == face_parents[i]) 
    {
      register u_char k;

      for (k=0;k<3;k++) t[j*3+k] = vertex_parents[ (*new_t)[i * 3 + k]];

      (*new_face_areas)[j] = simp_data->area[i];

      memcpy(*face_normals + 3 * j, simp_data->normal+i, sizeof(Vertex));

      j++;
    }
  
  
  for(i=0;i<simp_data->nT;i++) face_parents[i] = new_t_number[face_parents[i]];
 

  DXFree((Pointer) new_t_number); new_t_number = NULL;


  /* copy the new triangles into *new_t and reallocate *new_t */
 
  DXFree((Pointer) (*new_t));

  *new_t = (int *) DXAllocateZero(*new_nT * sizeof (Face));

  if (!(*new_t)) goto error;

  memcpy(*new_t, t, *new_nT * sizeof (Face));

  DXFree((Pointer) t); t = NULL;


  return 1;


error:

  if (new_t_number) DXFree((Pointer) new_t_number);
  if (new_v_number) DXFree((Pointer) new_v_number);

  if (t)            DXFree((Pointer) t);
  if (v)            DXFree((Pointer) v);

  
  return 0;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfMarkEdgesAdjacent2Boundary                                         |
  |                                                                          |
  +--------------------------------------------------------------------------+*/


/*
  As the name indicates, I am just using two different flags for edges adjacent to
  the surface boundary. 

  An edge that is a boundary edge, which is the same as an edge incident to two 
  boundary vertex, is said to be a boundary edge of the second kind,
  with a flag value of -2

  Otherwise, an edge incident to a single boundary vertex is said to
  be a boundary edge of the first kind, with a flag value of -1.

  This routine is executed *only* if the user does not wish to simplify the boundary.
  
  Once an edge is flagged with a negative flag value, it is impossible to add it
  later to the edge_heap.

 */

int _dxfMarkEdgesAdjacent2Boundary(SimpData *simp_data)
{

  int 
    e  = 0,
    nE = simp_data->nE,
    edges_adjacent_2_boundary;

  u_char 
    *boundary_vertex     = simp_data->boundary_vert,
    boundary1,
    boundary2;

  EdgeS *current_edge    = simp_data->edge_array;

  /* for each edge of the surface, determine whether it is a boundary
     edge of the first or second type. If yes, label it accordingly
     using the edge2index array */

  for (edges_adjacent_2_boundary = nE,
       e=0; e<nE; e++, current_edge++)
    {
      
      /* if the edge label is 1, then it is a boundary edge of the second kind */

      if (SurfEdgeLabel(current_edge) == 1) simp_data->edge2index[e] = -2;

      else
	{

	  int 
	    *v = SurfEdgeGetVertices(current_edge);

	  boundary1        = boundary_vertex[v[0]];
	  boundary2        = boundary_vertex[v[1]];

	  if (boundary1 || boundary2) simp_data->edge2index[e] = -1;
	  
	  else  edges_adjacent_2_boundary--;
	}
      
    }


  return edges_adjacent_2_boundary;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfBuildEdgeHeap                                                      |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

/* add to the Heap edges all the edges that can be added */

int _dxfBuildEdgeHeap(SimpData *simp_data)
{
  
  int 
    can_add_edge  = 1,
    n_edges_added = 0,
    nE            = simp_data->nE;

  while (can_add_edge && simp_data->last_edge_added2heap < nE)
    { 
      if ( simp_data->edge2index[simp_data->last_edge_added2heap]  == 
	   simp_data->heap_initial_index )

	   /* *new edges* are entered at this stage.  if an edge has
	      been visited before, it can only be reentered if a
	      neighbor was simplified.  _dxfReinstateEdgesInHeap does this
	      operation */
	
	n_edges_added += 
	  (can_add_edge = 
	   simp_data->add_edge(simp_data, simp_data->last_edge_added2heap));


      simp_data->last_edge_added2heap ++;

    } /* end while */

  return (n_edges_added);
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfRotateParentTriangle                                               |
  |                                                                          |
  +--------------------------------------------------------------------------+*/


/* performs a circular permutation (rotation) of the 
   vertex indices of a triangle such that the parent of the 
   first index is the same as "parent_vertex" */

int _dxfRotateParentTriangle(int  parent_vertex, int *vertex_fathers, 
			     int *tri_v, int *tri_v_parents, int *tri_v_rotated)
{

  /* the triangle must be a parent triangle itself (otherwise
	 all three triangle vertices are not necessarily different */

  tri_v_parents[0] = FATHER(tri_v[0], vertex_fathers);
  tri_v_parents[1] = FATHER(tri_v[1], vertex_fathers);
  tri_v_parents[2] = FATHER(tri_v[2], vertex_fathers);

  if (tri_v_parents[0] == parent_vertex)
    {
      tri_v_rotated[0] = tri_v[0];
      tri_v_rotated[1] = tri_v[1];
      tri_v_rotated[2] = tri_v[2];
    }
  else if (tri_v_parents[1] == parent_vertex)
    {
      /* rotate the vertex parents */
     
	  tri_v_rotated[0] = tri_v_parents[0]; /* use tri_v_rotated[0] as a temporary buffer */
      
      tri_v_parents[0] = tri_v_parents[1];
      tri_v_parents[1] = tri_v_parents[2];
      tri_v_parents[2] = tri_v_rotated[0];  

      /* rotate the vertices */
      tri_v_rotated[0] = tri_v[1];
      tri_v_rotated[1] = tri_v[2];
      tri_v_rotated[2] = tri_v[0];
    }
  else 
    {
      /* rotate the vertex parents */
      tri_v_rotated[0] = tri_v_parents[0]; /* use tri_v_rotated[0] as a temporary buffer */

      tri_v_parents[0] = tri_v_parents[2];
      tri_v_parents[2] = tri_v_parents[1];  
      tri_v_parents[1] = tri_v_rotated[0];
	      
      /* rotate the vertices */
      tri_v_rotated[0] = tri_v[2];
      tri_v_rotated[1] = tri_v[0];
      tri_v_rotated[2] = tri_v[1];
    }

  return 1;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfRotateParentEdge                                                   |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfRotateParentEdge(int parent_vertex, int *vertex_fathers,
			 int *edge_v, int *edge_t, int *edge_v_parents, 
			 int *edge_v_rotated, int *edge_t_rotated)
{
  

  /* rotate a parent edge with respect to a parent vertex, so that
     the child of the parent vertex will be first (index [0])
     among the two edge vertices 
     ( currently, by default the first vertex is the vertex with the
       highest number )
  */

  edge_v_parents[0] = FATHER(edge_v[0], vertex_fathers);
  edge_v_parents[1] = FATHER(edge_v[1], vertex_fathers);

  /* if the vertices are in the right order, just copy the vertices */

  if (edge_v_parents[0] == parent_vertex)
    {
      edge_v_rotated[0] = edge_v[0];
      edge_v_rotated[1] = edge_v[1];
      edge_t_rotated[0] = edge_t[0];
      edge_t_rotated[1] = edge_t[1];
    }
  else
    {
      edge_v_rotated[0] = edge_v_parents[0];  /* use edge_v_rotated[0] as a temporary buffer */

      edge_v_parents[0] = edge_v_parents[1];
      edge_v_parents[1] = edge_v_rotated[0];
      
      edge_v_rotated[0] = edge_v[1];
      edge_v_rotated[1] = edge_v[0];
      edge_t_rotated[0] = edge_t[1];
      edge_t_rotated[1] = edge_t[0];

    }

  return 1;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfDirectedParentVertexStar                                           |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

/* 
   compute a vertex star and store the triangles in either
   CLOCKWISE or COUNTER-CLOCKWISE direction
 */

int _dxfDirectedParentVertexStar(int parent_vertex, int first_triangle,
				 int valence, int *star, SimpData *simp_data, int direction)
{   


  /* labeling system for directed triangle stars 

	          v*1\      /v*0  
	           e*1\  *1/e*0   
	       *val0-1 \  / *0    
	     ___________\/___________ we omit on purpose the first vertex, because in the case
	     v*val0-1   v0            of a boundary edge, we know it already


	     */

  int 
    *vstar = star   + valence,
    *estar = vstar  + valence,
    /* rotated triangle vertices: */
    tri_v_rotated[3],  tri_v_parents[3], 
    /* rotated edge vertices and adjacent triangles: */
    edge_v_rotated[2], edge_v_parents[2], edge_t_rotated[2],
    star_triangle = first_triangle,
    i;

     
  for (i = 0; i < valence; i++)
    {
      star[i] = star_triangle;

      _dxfRotateParentTriangle(parent_vertex, simp_data->vertex_father,
			       simp_data->triangles + 3 * star_triangle,
			       tri_v_parents, tri_v_rotated);

      /* find previous or next edges and get their parents */

      estar[i] = _dxfNextParentEdge(tri_v_rotated, direction);

      _dxfRotateParentEdge(parent_vertex, simp_data->vertex_father,
			   (int *) SurfEdgeGetVertices (simp_data->edge_array + estar[i]),
			   (int *) SurfEdgeGetTriangles(simp_data->edge_array + estar[i]),
			   edge_v_parents, edge_v_rotated, edge_t_rotated);

      vstar[i] = edge_v_parents[1];
      /* this vertex number should be the same as tri_v_parents[direction] */

      if (i < valence-1) 
	/* otherwise, the next triangle would be -1 and finding its father 
	   would produce a segmentation fault or another bad thing */
	star_triangle = _dxfNextParentTriangle(edge_t_rotated, direction);
    } 


  return 1;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfManifoldLinkTest                                                   |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

/*
   detect duplicate vertices in the link of an edge.
   In case duplicate vertices are detected, collapsing an edge
   would generate a singular vertex.
  */

int _dxfManifoldLinkTest(int *link0, int *link1, int val0, int val1)
{
  /* among the vertices of the link of v0, link0 and of the link of v1, link1,
     is there any duplicate vertex ? */

  int 
    link_is_manifold = 1,
    val              = val0 + val1;

  /* since 
	val0 + val1 <= SIMPLIFY_VALENCE_MAX4,
	I can use a static array here. */

  static int
    link[SIMPLIFY_VALENCE_MAX4];

  {
    register int i;
      
    /* gather vertices of link0 and link1 in a single array, sort
       the array, and determine whether there are any duplicate
       vertices */

    memcpy(link,        link0, val0 * sizeof (int));
    memcpy(link + val0, link1, val1 * sizeof (int));

    qsort(link, val, sizeof(int), (int (*)(const void*, const void*)) _dxfCmpIntSmallFirst);

    i = 0;

    while ((i < val -1 ) && (link_is_manifold))
      {
	if (link[i+1] == link[i]) /* if there are duplicate vertices,
				     they must be contiguous in the
				     sorted array */
				       
	  link_is_manifold = 0;
	i++;
      }
     
  }

  return link_is_manifold; 

}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfCmpIntSmallFirst                                                   |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

/*
    compare integers. Using this function in qsort will
	result in listing the smaller integers first
  */

int _dxfCmpIntSmallFirst(char *s1, char *s2)
{
  int 
    *f1 = (int *)s1,
    *f2 = (int *)s2;

 if (*f2 > *f1)
    return(-1);

  else if (*f2 < *f1)
    return(1);

  else
    return(0);
  
}
/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfMakeEdgeBarycenter                                                 |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

/* 
    compute the barycenter of two vertices 
  */

void _dxfMakeEdgeBarycenter(Vertex v0, Vertex v1, float alpha0, Vertex v)
{  

  register u_char j;

  float alpha1 = 1. - alpha0;

  for(j=0; j<3; j++)

    v[j] = alpha0 * v0[j] + alpha1 * v1[j];
   
  return;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfClosestPointOnEdge                                                 |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

/* 
    compute the closest point XP to a point X on an edge (A,B) 

    return the barycentric cordinnates of X with respect to A and B
	
  */


float _dxfClosestPointOnEdge(Vertex X, Vertex XP, Vertex A, Vertex B, 
			     float *t)/* barycentric coordinates: t, 1-t */
{
  float d;
  Vertex BA,BX;
  float BA_BX,BA_2;

  MakeVec(B,A,BA);
  MakeVec(B,X,BX);
  EdgeBarycentricCoordinates(BA,BX,BA_BX,BA_2);

  if (BA_BX <= 0.)
    {
      *t = 0.;
      d = NORM3(BX);
      VecCopy(B,XP);
    }
  else if (BA_2 <= BA_BX)
    {
      *t = 1.;
      d = DIST3(A,X);
      VecCopy(A,XP);
    }
  else
    {
      /* the following is not robust numerically */
      /*
	d= SCALPROD3(BX,BX) - BA_BX*BA_BX / BA_2;
	if (d<=0.0) d = 0.0;
	else d =(float) sqrt ((double) d);
	*/

      /* it is more robust to compute d= NORM3(Vecprod(BA,BC))/BA 
	 or to just use the distance from X to XP */

      *t = BA_BX/BA_2;
	
      _dxfMakeEdgeBarycenter(A,B,*t,XP);
      
      d = DIST3(X,XP);

    }


  return(d);
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfNormalize3Vector                                                   |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

float _dxfNormalize3Vector(Vertex v)
{
  register u_char j;
  register double norm;

  for (norm=0.0,j=0;j<3;j++) norm+=v[j]*v[j];


  if (norm>0.0)
   {
    norm=sqrt(norm);
    for (j=0;j<3;j++) v[j] /= norm;
   } 

  return (float)norm ;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfSolveQuadraticEquation                                             |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfSolveQuadraticEquation(double A, double B, double C, 
			       double *sol1, double *sol2, double eps)
{
  /*
  
    returns the number of real solutions of a second degree equation as well as the solutions

   everything below eps is supposed to be 0
   
   */

  int num_sol;

  A *= 2.;

  if ( fabs( A ) < eps)
    {
      /* B x + C = 0 */
      
      if ( fabs( B ) < eps)
	num_sol = 0;

      else 
	{
	  num_sol = 1;
	  *sol1   = - C / B;
	}
    }

  else
    {
      double delta  =  B * B - 2. * A * C;

      if ( delta < - eps )
	num_sol = 0;

     else
       {
	 /* the discriminant is zero or positive */
	
	 if (delta < eps)
	   {
	     num_sol = 1;
	      
	     *sol1 = -B / A;
	   }
	 else
	   { 
	     delta = sqrt (delta);

	     num_sol = 2;

	      *sol1 = (-B - delta ) /A;

	      *sol2 = (-B + delta ) /A;
	    }
       }
    }

  return num_sol;

}
/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfSimplifiedBoundaryVertexLocation                                   |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

/* 
   determine the location of a potential simplified boundary vertex. 
   applies only when a boundary edge is collapsed, which implies that
   the two endpoints of the edge are boundary vertices
 */


int _dxfSimplifiedBoundaryVertexLocation(SimpData *simp_data, int v0, int v1,
					 int *vstar0, int *vstar1, 
					 int val0, int val1, int method)
{
  int success = 1;

  if (method == 0)

    memcpy(simp_data->simplified_vertex, simp_data->vert[v0], sizeof (Vertex));

  else if (method == 1)
    {

      /* place the simplified vertex so that the length of the boundary
	 will be preserved.  The vertex must be situated on an ellipse 

      ___________\/___________\/__________
      v*val0-1   v0           v1          v*val1-1
      == vl                                  == vr

      Firstly, we retrieve the vertices vl, v0, v1, vr from the stars.

      we compute the length of the polygonal arc vl--v0--v1--vr, and
      call 2a this length
     
      The focal points of the ellipsoid are vr and vl, so we note 2c the distance
      vr-vl.

      */

      float 
	a, a2, c, eps1 = (float)1e-6,
	*ve = (float *) simp_data->simplified_vertex;

      Vertex 
	v_zero, v_one, vr, vl, vm;
  
      memcpy(v_zero, simp_data->vert[v0], sizeof (Vertex));
      memcpy(v_one,  simp_data->vert[v1], sizeof (Vertex));
  
      memcpy(vl, simp_data->vert[vstar0[val0-1]], sizeof (Vertex));
      memcpy(vr, simp_data->vert[vstar1[val1-1]], sizeof (Vertex));
      
      /* compute the parameters of the ellipsoid */

      a  = 
	(
	 DIST3(simp_data->vert[v0], simp_data->vert[v1]) +
	 DIST3(simp_data->vert[vstar0[val0-1]], simp_data->vert[v0]) +
	 DIST3(simp_data->vert[vstar1[val1-1]], simp_data->vert[v1])
	 ) 
	/ 2.;

      a2 = a * a;
  
      c  = DIST3(vl,vr)/2.;
  
      /* Compute the middle of v_zero and v_one */

      vm[0] = (v_zero[0] + v_one[0]) / 2.;
      vm[1] = (v_zero[1] + v_one[1]) / 2.;
      vm[2] = (v_zero[2] + v_one[2]) / 2.;
 
      /* Then, our first next test is to determine whether a and c are
	 sufficiently different, otherwise, taking the middle of v_zero and
	 v_one will certainly work since we will have the same cord length
	 along the boundary */

      if (a == 0.0 || /* if a == 0, the length to preserve is zero, then the segment middle is fine! */
	  1. - c/a < eps1) 

	/* limit of the floating point precision */

	{
	  memcpy(ve, vm, sizeof(Vertex));
	}

      else /* if (a > c ) */
	{
	  float 
	    t,  x,  y, 
	    ye;

	  Vertex 
	    vp;
  
	  /* then we consider the middle of v0 and v1, vm,
	     we project vm on the segment vr vl and obtain vp

	     finally, we project vp on the ellipsoid defined by vr, vl, a
	     and c along the direction vp--->vm


                                  ve


                                  vm
                                  |y
           vl________________0____|_____________vr
           -c                     x              c
                                   vp


			     by projecting vm onto the segment (vl,vc), we obtain
			     vp, and the coordinates of vm in the coordinate system formed by
			     0(the middle of vr and vl) vr and vm.

           x and y are the coordinates of vm in this coordinate system
           We keep x unchanged, and compute ye such that

           x^2/a^2 + ye^2/b^2 = 1, with b^2 = a^2-c^2

                                           x^2
           hence        ye^2    =  (  1 -  ---  ) ( a^2 - c^2 )
                                           a^2

           We take for ye the positive square root of ye^2 

	   and finally ve is obtained by

	   ve-vp = ye/y (vm-vp) 

	   ve    = vp + ye/y (vm-vp)


	   5/7/97 bug of Floating Point exception noticed on DEC
	   A special case is encountered if c == 0,
	   then we also have x = 0 and the formula is simply ye=a
	   this case is handled with the general case without particular
	   attention, except that we can't divide by c and
	   also we should not execute the part of the code relative to t == 0 or t == 1
	   */

 
	  /* compute x and y of vm with respect to vr and vl */
	  
	  y  = _dxfClosestPointOnEdge(vm, vp, vl, vr, &t);


	  /* what do we do if y == 0, in which case vm is accidentally
	     exactly on the line segment vr-vl?

	     In that case, we still have a non-degenerate ellipse since we
	     ruled out the fact that a and c are the same. Then either
	     v_zero and v_one cannot be on the line segment vr-vl.

	     So we project instead v0 on the segment vr vl 
	 
	          v0
	      ___|________________
            vl   vp0      vp      vr

                                ----->
		   and we use the direction vp0 v0 for projecting the point vp on the  ellipse

		   */


	  if (y <= eps1 * c) /* it can be that c == 0 and y == 0 */
	    {
	      /* limit of floating point precision */
	      /* In that case, either v_zero or v_one cannot be on the line vl, vr,
		 so we replace vm with v_zero */
	  
	      memcpy(vm, v_zero, sizeof(Vertex));
	  
	      y  = _dxfClosestPointOnEdge(vm, vp, vl, vr, &t);

	      if (y <= eps1 *c)
		{
		  memcpy(vm, v_one, sizeof(Vertex));
	  
		  y  = _dxfClosestPointOnEdge(vm, vp, vl, vr, &t);
		}
	    }

	  /* what would it mean for y to be equal to zero after these
	     various projections: that all points v_zero, v_one, vm are on vl_vr,
	     meaning that a = c. Just in case, we use an if (y!=0.0) in the end */

	  /* in case vm does not project orthogonally on the segment (vl,vr)
	     the routine ClosestPointOnEdge projects vm either on vl and assigns
	     1 to t or on vr and assigns 0 to t 
	     we need to determine the angle that vm,vp makes with (vl,vr)
	     */

	  if (c > 0.0 /* c == 0 is handled as the general case */
	      && (t == 0. || t == 1.))
	    {
	      Vertex u,v;
	      double 
		cos_uv, cos_uv_2, sin_uv_2,
		A, B, C, sol1 = 0, sol2 = 0, eps2 = 1e-12;

	      MakeVec(vl,vr,u); /* if the segment (vr,vl) is of zero length, which is
				   possible in nasty cases, u will be a vector of zero
				   lenght, and this invalidates the following computations
				   we need to find another solution */
	      _dxfNormalize3Vector(u);
	  
	      MakeVec(vp,vm,v);
	      _dxfNormalize3Vector(v);

	      cos_uv   = u[0] * v[0] + u[1]*v[1] + u[2]*v[2]; /* FScalProd(u,v); */

	      if (cos_uv < 0.)
		cos_uv = -cos_uv;

	      cos_uv_2 = cos_uv * cos_uv;

	      sin_uv_2 = 1. - cos_uv_2;

	      /* ye is the solution of the following equation of the second
		 degree:

		 x^2/a^2 + y^2/(a^2-c^2) = 1 <=>
		 
		 (ye cos_uv + c)^2/a^2 + (ye sin_uv)^2/(a^2-c^2) = 1 <=>

		 ye^2( cos_uv_2 + a^2/(a^2-c^2) sin_uv_2) +
		 ye  (2 c cos_uv) +
		 c^2-a^2                                  = 0


		 A ye^2 + cb ye + C = 0 */

	      A = cos_uv_2 + sin_uv_2 * a2 / (a+c)/(a-c);
	      B = 2. * cos_uv * c;
	      C = (c+a)*(c-a); /* c^2-a^2 */

	      _dxfSolveQuadraticEquation(A, B, C, &sol1, &sol2, eps2);

	      /* there is only one positive solution */

	      if (sol1 > 0.)
		ye = sol1;
	      else if (sol2 > 0.)
		ye = sol2;
	      else 
		{
		  ye = 0;

		  /* The quadratic equation has no positive root */
		}
	    } /*endif vm project on an extremity */

	  else /* vm projects normally inside the segment */
	    {
              double tmp;

	      /* determine x from t: */

	      /* vp = t * vl + (1-t) * vr 
		 t == 1 <==> x = -c
		 t == 0 <==> x =  c
	 
		 x = c - 2 c t = c ( 1-2t) */
	      
	      x = c * (1. - 2. * t);
  

	      /* ye = sqrt ((double) (a2 - c * c)*( 1. - x*x/a2)); */

	      /* since a and c are likely to be very close, a more robust
		 computation is to avoid computations such as a2 - c^2
		 and compute sqrt(a2) is a is readily available */

	      tmp = (double) (a -c) * (a+c) * (a-x) * (a+x);

	      /* this is to avoid floating point exceptions when computing the square root */

              if (tmp >0.0) ye = sqrt (tmp) / a;

	      else ye = 0.0;

	    }
	
	  if (y != 0.0)
	    ye /= y;

	  /* do the following vector composition without
	     a function call ve    = vp + ye/y (vm-vp) */
    
	  ve[0] = vp[0] + ye * (vm[0] - vp[0]);
	  ve[1] = vp[1] + ye * (vm[1] - vp[1]);
	  ve[2] = vp[2] + ye * (vm[2] - vp[2]);

	  /* test whether the boundary length is preserved */
	  {

	    float 
	      l_before = 2. * a,
	      l_after  = DIST3(vl,ve) + DIST3(ve,vr),

	      eps3 = simp_data->tolerance * eps1; /* this is a relative accuracy
						     tolerance * floating point accuracy
						     using eps1 here would not be general enough */

	    if (l_before -l_after > eps3 || l_before - l_after < -eps3)

	      {
		success = 0;
	      }

	  }

	} /* if (a>c) otherwise, vm, the middle is good for ve */

    } /* if (method = 1, length preservation) otherwise */
  
  return success;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfFastCompactness3DTriangle2                                         |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

/* compute the triangle compactness and area
   this is a fast method, in comparison with the method that 
   orthogonalizes the triangle and computes the area afterwards 
   
   this fast method uses a cross product 
 */

float _dxfFastCompactness3DTriangle2(Vertex *tri, float *triangle_area)
{
  
  /* compute the squared lenghts of the three edges: */
  
  double sum_lengths_sq = 0., tmp;
   
  {

    register u_char i,j;

    for(i=0;i<3;i++)
      {
	double len_edg_sq = 0;

	for(j=0;j<3;j++)
	  {
	    tmp         = tri[(i+1)%3][j] - tri[i][j];
	    len_edg_sq += tmp * tmp;
	  }

	sum_lengths_sq += len_edg_sq;

      }
  }
      

  {
    double 
      m0 =    
	       tri[1][1] - tri[0][1],
      m1 =
	       tri[2][2] - tri[1][2],
      m2 =
	       tri[1][2] - tri[0][2],
      m3 =
	       tri[2][1] - tri[1][1],
      m4 =
	       tri[2][0] - tri[1][0],
      m5 =
	       tri[1][0] - tri[0][0],
      n0  =
	       m0 * m1 -
	       m2 * m3,
		 
      n1  =
	       m2 * m4 -
	       m5 * m1,
			 
      n2  =
	       m5 * m3 -
	       m0 * m4;
    
    tmp            = n0 * n0 + n1 * n1 + n2 * n2;

    if (tmp > 0.0) *triangle_area = sqrt (tmp) / 2.;

    else *triangle_area = 0;
     
	       
  }       

  
  return (float) (FOUR_SQRT_3 * *triangle_area / sum_lengths_sq);

}


/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfBoundaryCollapse2KndGeometricallyOk                                |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfBoundaryCollapse2KndGeometricallyOk(SimpData *simp_data, int v0, int v1,
					
					    /* this procedure must be run twice, once for the
					       star of the vertex v0 and once for the star of the vertex v1
					       Here, v1 rather stands for the first vertex of the star
					       */

					int *star0, int *vstar0, int val0, Vertex *s_normal0, float *s_area0,
					int direction, float min_scalprod_nor, float compactness_ratio) 

{ 
  int collapsible = 1, i = 0;

  Vertex faceVert[3];
  
  float 
    *s_comp0 = s_area0 + val0,
    min_compactness_before,
    min_compactness_after;

  /* labeling system for directed triangle stars 

	          v*1\      /v*0  
	           e*1\  *1/e*0   
	       *val0-1 \  / *0    
	     ___________\/___________   we omit on purpose the first vertex, because in the case
	     v*val0-1   v0          v1  of a boundary edge, we know it already


	     */

  /* retrieve the compactness of the first triangle */

  min_compactness_before = simp_data->compactness[star0[0]];

  min_compactness_after  = 1; /* initialization of the compactness
				 after simplification */
  
  while (collapsible && i < val0-1)
    {
   
      /* retrieve the face compactness before simplification */
   
      min_compactness_before = MIN ( simp_data->compactness[star0[i+1]], min_compactness_before);

      /* recompute the compactness after changing the first vertex */

      memcpy(faceVert[0], simp_data->simplified_vertex, sizeof (Vertex));
    
      if (direction == DIRECTION_CLOCKWISE)
	{
	  /* reverse the order of triangle vertices, so that the normal
	     that we will compute will respect the surface orientation */
	  memcpy(faceVert[2], simp_data->vert[vstar0[i]], sizeof (Vertex));
	  memcpy(faceVert[1], simp_data->vert[vstar0[i+1]], sizeof (Vertex));
	}

      else
	{
	  memcpy(faceVert[1], simp_data->vert[vstar0[i]], sizeof (Vertex));
	  memcpy(faceVert[2], simp_data->vert[vstar0[i+1]], sizeof (Vertex));
	}

      s_comp0[i+1] = _dxfFastCompactness3DTriangle2(faceVert, s_area0+i+1);

      min_compactness_after = MIN ( s_comp0[i+1], min_compactness_after);

      /* recompute the face normal assuming that the vertex v0 has
         moved and save this normal for future use*/

      _dxfTriangleNormalQR2(faceVert, s_normal0[i+1]);


      /* check whether the triangle normal orientation is consistent */

      collapsible = ((SCALPROD3(s_normal0[i+1], simp_data->normal[star0[i+1]])) > min_scalprod_nor);

      i++;
    }

  if (collapsible)

    /* check whether the minimum triangle compactness is too much degraded */

    collapsible = min_compactness_after > (compactness_ratio * min_compactness_before);
  
  return collapsible;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfErrorTestDataOnSurfaceUpdate                                       |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

/* 
   given a new vertex of the mapping between simplified and original surfaces,
   update the potential data error value at the potential simplified vertex location 
   
   */
   


int _dxfErrorTestDataOnSurfaceUpdate(SimpData *simp_data, int new_v, int old_v, 
				      int new_v2, int new_v3,
				      float alpha_sv, float alpha_n_v2, float alpha_n_v3,
				      int old_v1, int old_v2, int old_v3, 
				      float alpha_o_v1, float alpha_o_v2, float alpha_o_v3)
{

  int do_collapse = 1;

  float 
    o_alpha[3],
    n_alpha[3],
    *o_data[3],
    *n_data[3],
    o_data_error[3],
    data_difference,
    new_error_link,
    tmp_data_diff=0,
    new_error_sv,
    old_error;
  	
  register int j;
  register u_char k;

  o_data[0]       = simp_data->vx_data + old_v1 * simp_data->data_dim;
  o_data[1]       = simp_data->vx_data + old_v2 * simp_data->data_dim;
  o_data[2]       = simp_data->vx_data + old_v3 * simp_data->data_dim;

  n_data[0]       = simp_data->vx_data_potential_values;
  n_data[1]       = simp_data->vx_data + new_v2 * simp_data->data_dim;
  n_data[2]       = simp_data->vx_data + new_v3 * simp_data->data_dim;
	
  o_alpha[0]      = alpha_o_v1;
  o_alpha[1]      = alpha_o_v2;
  o_alpha[2]      = alpha_o_v3;
 	
  n_alpha[0]      = alpha_sv;
  n_alpha[1]      = alpha_n_v2;
  n_alpha[2]      = alpha_n_v3;
 
  o_data_error[0] = simp_data->vx_data_error[old_v1];
  if (old_v > 1) 
    o_data_error[1] = simp_data->vx_data_error[old_v2];
  if (old_v > 2) 
    o_data_error[2] = simp_data->vx_data_error[old_v3];
  


  /* interpolate vertex data in the old and new surfaces 
     and decide if they are sufficiently similar */

  for (data_difference = 0., j=0; j<simp_data->data_dim; j++)
    {
      simp_data->vx_old_data[j] = 0.0;

      for (k=0 ; k < old_v; k++) simp_data->vx_old_data[j] += o_alpha[k] * o_data[k][j];

      simp_data->vx_new_data[j] = 0.0;
	  
      for (k=0 ; k < new_v; k++) simp_data->vx_new_data[j] += n_alpha[k] * n_data[k][j];

      tmp_data_diff    = simp_data->vx_new_data[j] - simp_data->vx_old_data[j];
	  
      data_difference += tmp_data_diff * tmp_data_diff;

    }


  /* there are two cases, either the data dimension is one or it is larger than one
     in the first case, no square root is required for computing the discrepancy
     between data values, and since I am very stingy about expensive
     floating point computations, I am avoiding taking this square root if I can.
     */

  if (simp_data->data_dim == 1) {data_difference = (tmp_data_diff > 0.0) ? tmp_data_diff : -tmp_data_diff;}

  else
     if (data_difference > 0.0) data_difference = sqrt (data_difference);
	  
  old_error = 0.0;

  for (k=0 ; k < old_v; k++) old_error += o_alpha[k] * o_data_error[k];

  /* can the sum of the data difference and the old error at the
     current mesh point be explained by a new error at the simplified
     vertex, and if it is the case, what is this new error ? */

  new_error_link = (new_v == 3) ?
    alpha_n_v2 * simp_data->vx_data_error[new_v2] + alpha_n_v3 * simp_data->vx_data_error[new_v3]:
    (new_v == 2) ? alpha_n_v2 * simp_data->vx_data_error[new_v2]:
    0.0;


  /* there is a special case if the simplified vertex is associated with an
     alpha that's zero: in this case, the error must be explained with the error at
     the link */

  if (alpha_sv == 0. )
    {
      if (new_error_link < old_error + data_difference)
	do_collapse = 0;

      /* but otherwise no potential error is computed */

    }
  else /* alpha_sv * error_sv + new_error_link > old_error + data_difference */
    {
      new_error_sv = (old_error + data_difference - new_error_link) / alpha_sv;

      if (new_error_sv > simp_data->vx_data_potential_err) simp_data->vx_data_potential_err = new_error_sv;

      do_collapse = (simp_data->vx_data_potential_err <= simp_data->data_tolerance);
    }


  return do_collapse;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfSegmentIntersection                                                |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

float _dxfSegmentIntersection(Vertex A, Vertex B, Vertex C, Vertex D,
			      float *lambda, float *mu, 
			      int changed_AorB, int changed_C, int changed_D)
{
  /* check whether for two segments in R3, the intersection is in the
    interior of each of them.

    lamda values of 1 and 0 qualify as being in the interior, except
    when one of two edges is of length zero, in which case we set both
    lambda and mu to be -1, to reflect the fact that there is no
    intersection */

  float dist = 0;

  static VertexD w[3];
  static double  work[2], max_z, min_z, sign = 0., scale;

  register u_char i,j;


 
  /* if both Vertices A and B did not change since the previous call
     of the routine, the scaling factor, householder vector and sign
     are still correct and do not need to be re-evaluated (this spares
     us the computation of one square root */

  /* if C or D changed, we only have to recompute the portion
     corresponding to C or D 
     in surface simplification, often the only variable that changes
     is D */

  if (changed_AorB)
    {

      for (j=0;j<3;j++)
	{
	  /* substract A from the rest of the vertices */

	  w[0][j] = B[j]-A[j];
	  w[1][j] = C[j]-A[j];
	  w[2][j] = D[j]-A[j];
	}

      /* record the norm of the first vector and scale the 
	 vertices accordingly */

      scale = w[0][0] * w[0][0] + w[0][1] * w[0][1] + w[0][2] * w[0][2];
      if (scale > 0.0) scale = sqrt (scale);

      /* particular case if the scale is zero */

      if (scale > 0.)
	{
	  /* apply a scaling factor so that w[0][2] will be negative */

	  sign = (w[0][2] > 0.) ? -1. : 1.;

	  sign /= scale;

	  for (i=0;i<3;i++) 
	    for (j=0;j<3;j++)
	      w[i][j] *= sign; /* apply the scaling */
    
      /* determine the Householder vector for transforming
	 w[0] to z */
	  
	  w[0][2] -= 1.;
     
	  _dxfHouseholderPreMultiplicationDbl((double *) w[1], 3,
					  (double *) w[0], 3, 2, work);
   

	  /* now the first segment is the same as (0,0,1) */
	}
    }

  else 
    {
      if (changed_C)
	{
	  for (j=0;j<3;j++)
	    w[1][j] = (C[j]-A[j]) * sign; /* apply scaling */

	  _dxfHouseholderPreMultiplicationDbl((double *) w[1], 3, /* apply Householder */
					  (double *) w[0], 3, 1, work);
	}

      if (changed_D)
	{
	  for (j=0;j<3;j++)
	    w[2][j] = (D[j]-A[j]) * sign; /* apply scaling */

	  _dxfHouseholderPreMultiplicationDbl((double *) w[2], 3, /* apply Householder */
					  (double *) w[0], 3, 1, work); 
	}
    }
	

#define xC_ w[1][0]
#define yC_ w[1][1]
#define zC_ w[1][2]
#define xD_ w[2][0]
#define yD_ w[2][1]
#define zD_ w[2][2]
     
  if (scale > 0.)
    {
     /* there are two cases in which there is no intersection 
	if zC and zD >1 or if zC and zD < 0 */

     if (zC_> zD_) {max_z = zC_; min_z = zD_;}
     else          {max_z = zD_; min_z = zC_;}
    
     if ((min_z > 1.) || (max_z  < 0.))
	
       {
	  *lambda = *mu = -1.;
       }
     
     else
       {
	 static double xDC, yDC, zDC, dlambda, dmu, a, half_b, c, tmp;
	 /* x = xD + mu xDC 
	    y = yD + mu yDC */

	 xDC = xC_ - xD_ ;
	 yDC = yC_ - yD_ ;
	 zDC = zC_ - zD_ ;
	
	 c      = xD_ * xD_ + yD_ * yD_;
	 a      = xDC * xDC + yDC * yDC;
	 
	 if (a > 0.)
	   {
	     half_b  = xD_  * xDC + yD_ * yDC;

	     dmu     = -half_b / a;
	     *mu     = dmu;                     /* mu = -b/(2a) */

             tmp     =  c + half_b * dmu;
             if (tmp >0.0) tmp = sqrt (tmp);
             else tmp= 0.0;

	     dist    =  scale * tmp; /* c - b^2/4a   */

	     /* lambda is given by the 1-z_value at that point */

	     *lambda = 1. - zD_ - dmu * zDC;
	   }

	 /* otherwise, either the two segment are parallel or the
	    second segment is of zero length */
	 
	 else if (min_z < max_z)
	   {

	     /* the second segment is parallel to the first *and*
		there is an intersection */

	     min_z   = MAX(min_z,0.);
	     max_z   = MIN(max_z,1.);

	     dlambda = (max_z + min_z)/2.; 
	   
	     /* this is definitely between [0,1] and corresponds to
		a point inside CD to */

	     *mu = (dlambda - zD_) / zDC;
	     
	     *lambda = 1. - dlambda;
	     dist    = scale * sqrt (c);
	     
	   }
	 else
	   {
	     /* the second segment is degenerate */
	     *lambda = *mu = -1.;
	   }

       }
    }/* endif (scale > 0) */

#undef xC_ 
#undef yC_ 
#undef zC_ 
#undef xD_ 
#undef yD_ 
#undef zD_ 

    
  else
    {
      /* the first edge is degenerate. There is no mutual intersection */

      *lambda = *mu = -1.;

    }

  return dist;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfHouseholderPreMultiplication                                       |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfHouseholderPreMultiplication(
			     float *A,  /* input matrix, given column
				           by column */
			     int mA,    /* leading dimension of matrix a */
			     float *v,  /* householder vector    */
			     int m,     /* dimension of v        */
			     int n,     /* second dimension of A */
			     float *w   /* workspace, dimension n */
				 ) 
{
  int retval = 0;
  
  register int i,j  ;
  
  float *w1 = w;

  /* 1) compute beta = -2 vT v: */

  float beta = 0; 

  for (i=0;i<m;i++) beta += v[i]*v[i];

  beta = -2./beta;

  /* 2) compute w = beta AT v */
  
  if (w1 == NULL)
    {
      w1 = (float *) DXAllocateZero (n * sizeof (float));
    }
  
  if (w1)
    { 
      float *A1 = A;

      for (i=0;i< n; i++, A1 += mA) 
      
	{
	  for (w1[i]=0., j=0; j<m;j++)
	
	    w1[i] += A1[j]*v[j];

	  w1[i] *= beta;
	}

      /* 3) compute A = A + w vT */

      A1 = A;

      for (i=0;i< n; i++, A1 += mA) 
	
	for (j=0; j<m;j++)
	  
	  A1[j] += w1[i] * v[j];
	
      retval = 1;

      if (w1 && w != w1)
	DXFree((Pointer) w1);
    
    }
 
  return retval;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfHouseholderPreMultiplicationDbl                                    |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfHouseholderPreMultiplicationDbl(
			     double *A,  /* input matrix, given column
				           by column */
			     int mA,    /* leading dimension of matrix a */
			     double *v,  /* householder vector    */
			     int m,     /* dimension of v        */
			     int n,     /* second dimension of A */
			     double *w   /* workspace, dimension n */
				 ) 
{
  int retval = 0;
  
  register int i,j  ;
  
  double *w1 = w;

  /* 1) compute beta = -2 vT v: */

  double beta = 0; 

  for (i=0;i<m;i++) beta += v[i]*v[i];

  beta = -2./beta;

  /* 2) compute w = beta AT v */
  
  if (w1 == NULL)
    {
      w1 = (double *) DXAllocateZero (n * sizeof (double));
    }
  
  if (w1)
    { 
      double *A1 = A;

      for (i=0;i< n; i++, A1 += mA) 
      
	{
	  for (w1[i]=0., j=0; j<m;j++)
	
	    w1[i] += A1[j]*v[j];

	  w1[i] *= beta;
	}

      /* 3) compute A = A + w vT */

      A1 = A;

      for (i=0;i< n; i++, A1 += mA) 
	
	for (j=0; j<m;j++)
	  
	  A1[j] += w1[i] * v[j];
	
      retval = 1;

      if (w1 && w != w1)
	DXFree((Pointer) w1);
    
    }
  
  return retval;
}


/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfErrorWithinToleranceVBoundary2ndKnd                                |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfErrorWithinToleranceVBoundary2ndKnd(SimpData *simp_data, int v0, int v1, 
					    int *vstar0, int *vstar1, int val0, int val1)
{
  int edge_collapsible = 1;

  /* we construct a tiling of (e*) and a tiling of (vs*)
     such that there is a one to one mapping between the two tilings */
  
 
  /*   e*:
	          v*1\      /v*0  \      /
	           e*1\  *1/e*0    \ *1 /
	       *val0-1 \  / *0   *0 \  / *val1-1
	     ___________\/___________\/__________
	     v*val0-1   v0    vs      v1
             vl                                  vr
	     */

  /* The difference between the interior vertex and boundary vertex
    approaches is that vs v0 and v1 are projected on edges exclusively.
    That limits the number of edge intersections that must be computed
   */

  /* 1) determine the closest point to the vertex vs on the edge v0,v1 */

  /* 2) determine the closest point to the vertex v0 on the edge 
     (vl,vs) */

  /* 3) determine the closest point to the vertex v1 on the edge (vs,vr) */

  /* 4) take in turn each edge (vs,v*1)    (vs,v*val0-1) 
     and intersect them with each edge (v0,v*0),  (v0, v*val0-1) */


  /* 5) do the same as 4) for the star of the vertex v1 
   */

  register int i;


  float 
    potential_err_pivot,
    largest_err_pivot = 0.,

    /* set the maximum tolerance for the star */
  
    tol_pivot = MIN(simp_data->tol_volume[v0], simp_data->tol_volume[v1]);

  for (i=0;i<val0;i++)
    if (simp_data->tol_volume[vstar0[i]] < tol_pivot)
      tol_pivot = simp_data->tol_volume[vstar0[i]];

  for (i=1;i<val1;i++) 

    /* the first vertex of the star of v1 is the same as the first
       vertex of the star of v0: vstar1[0] = vstar0[0] */

    if (simp_data->tol_volume[vstar1[i]] < tol_pivot)
      tol_pivot = simp_data->tol_volume[vstar1[i]];

  /* verify that the maximum tolerance is not already exceeded at
     the vertices */

  i = 0;

  while (edge_collapsible && i < val0)
  /* the equality (star_err == tol) is allowed */  
    if (simp_data->err_volume[vstar0[i++]] > tol_pivot)
      edge_collapsible = 0;

  i = 1;
  while (edge_collapsible && i < val1)
    if (simp_data->err_volume[vstar1[i++]] > tol_pivot)
      edge_collapsible = 0;
  
  if (simp_data->err_volume[v0] > tol_pivot || simp_data->err_volume[v1] > tol_pivot)
    edge_collapsible = 0;

  if (edge_collapsible)
    {
      float 
	t, dist;

      Vertex wp[2];

      dx_bzero ((char *)wp, 2 *sizeof (Vertex));

      _dxfErrorTestDataOnSurfaceInit(simp_data);

      /* 1 determine the closest point to the vertex vs on the edge v0,v1 */
     
      dist = _dxfClosestPointOnEdge(simp_data->simplified_vertex, wp[0],
				    simp_data->vert[v0], simp_data->vert[v1],&t);

      _dxfErrorTestDataOnSurfaceInterpolate2(simp_data, v0, v1, t);

      potential_err_pivot = dist + t * simp_data->err_volume[v0] + (1.-t) * simp_data->err_volume[v1];
     
      UPDATE_ERR_PIVOT;
     

      if (edge_collapsible)
	{
	  
	  /* 2) determine the closest point to the vertex v0 on the edge 
	     (vl,vs) */

	 
	  dist = _dxfClosestPointOnEdge(simp_data->vert[v0], wp[0],
					simp_data->vert[vstar0[val0-1]], simp_data->simplified_vertex, 
					&t);

	  /* t * err[vl] + 1-t * err[vs] >= dist + err[v0] */

	  if (t < 1.)
	    {
	      potential_err_pivot = 
		
		(dist + simp_data->err_volume[v0] -t * simp_data->err_volume[vstar0[val0-1]] )/ (1.-t);

	      UPDATE_ERR_PIVOT;
	    }

	  /*
	    v0
	    if (1-t == 0 <==> t = 1)                             |
	    that means that v0 projects orthogonally on vl      vl____________vs

	    we can't treat this constraint and the edge is not collapsible 
	    unless err[vl] >= dist + err[v0]
	    */
	  
	  else if (simp_data->err_volume[vstar0[val0-1]] < dist + simp_data->err_volume[v0])
	    
	    edge_collapsible = 0;

	  if (edge_collapsible)
	    {
	      if (simp_data->vx_data)

		edge_collapsible =
		  _dxfErrorTestDataOnSurfaceUpdate(simp_data, 
						   2, 1, vstar0[val0-1], 0, 1.-t, t, 0.,  /* vs, vl */
						   v0, 0, 0, 1., 0., 0.);    

	      /* 3) determine the closest point to the vertex v1 on the edge (vs,vr) */

	      dist = _dxfClosestPointOnEdge(simp_data->vert[v1], wp[0],
					    simp_data->vert[vstar1[val1-1]], simp_data->simplified_vertex, &t);

	      /* t * err[vr] + 1-t * err[vs] >= dist + err[v1] */
                                
	      if (t < 1.)
		{
		  potential_err_pivot = 
		    (dist + simp_data->err_volume[v1] -t * simp_data->err_volume[vstar1[val1-1]])
		    / (1.-t);

		  UPDATE_ERR_PIVOT;
		}
	  
	      else if (simp_data->err_volume[vstar1[val1-1]] < dist + simp_data->err_volume[v1])
		
		edge_collapsible = 0;

	
	      if (simp_data->vx_data && edge_collapsible)

		edge_collapsible = _dxfErrorTestDataOnSurfaceUpdate(simp_data, 
							   2, 1, vstar1[val1-1], 0, 1.-t, t, 0.,  /* vs, vl */
							   v1, 0, 0, 1., 0., 0.);    
      

	      if (edge_collapsible)
		{
		  float 
		    lambda1, lambda2;

		  int j;

		  /* 4) take in turn each edge (vs,v*1)    (vs,v*val0-1) 
		     and intersect them with each edge (v0,v*0),  (v0, v*val0-1) */

		  i = 1;
		  while (edge_collapsible && i<val0)
		    {
		 
		      j = 0;
		      while (edge_collapsible && j<val0)
			{
		
			  dist = _dxfSegmentIntersection(simp_data->simplified_vertex,
							 simp_data->vert[vstar0[i]], 
							 simp_data->vert[v0], 
							 simp_data->vert[vstar0[j]],
							 &lambda1, &lambda2, !j, 0, 1);
					       
			  /* If the closest point is located in the
			     inside of either edge, determine a crossing point. */

			  if (lambda2 >= 0. && lambda2 <= 1.) /* inside 2 */
			    {
			      dist += lambda2   * simp_data->err_volume[v0] + 
				(1. -lambda2)  * simp_data->err_volume[vstar0[j]];


			      if (simp_data->vx_data && lambda1 >= 0. && lambda1 <= 1.)
				
				edge_collapsible = 
				  _dxfErrorTestDataOnSurfaceUpdate(simp_data, 2, 2, 
                                            vstar0[i], 0, lambda1, 1. -lambda1, 0.,
				        v0, vstar0[j], 0, lambda2, 1. -lambda2, 0.);


			      if (lambda1 > 0. && lambda1 <= 1.) /* inside 1 */
		      
				{
				  potential_err_pivot = 
				    (dist - (1. -lambda1)  * simp_data->err_volume[vstar0[i]])
				    / lambda1;

	    
				  UPDATE_ERR_PIVOT;
			      
				} 
			  
			      /* if lambda1 = 0 we can't accomodate
				 that constraint with the error value at the simplified vertex
				 and the edge might be declared "non-collapsible"

				 vs_______v*1
                                          |
				       v0_|_v*0

				       dist == 0 can occur when the two vertices
				       vstar0[j] and vstar0[i] are the same
				       but we prohibit that */

			      else if (lambda1 == 0.) 
				{
				  /* the error at vs cannot accomodate
				     that constraint, but maybe it is
				     already satisfied */

				  if (simp_data->err_volume[vstar0[i]] < dist)
			    
				    edge_collapsible = 0;
				}
			    }

		    

			  j++;
			  if (i == j) j++; /* we prohibit (i == j) */
			}
		      i++;
		    }

		  if (edge_collapsible)
		    {
		      /* 5) do the same as 4) for the star of the vertex v1:

			 take in turn each edge (vs,v*1)    (vs,v*val1-1) 
			 and intersect them with each edge (v1,v*0),  (v1, v*val1-1) */

		      i = 1;
		      while (edge_collapsible && i<val1)
			{
		 
			  j = 0;
			  while (edge_collapsible && j<val1)
			    {
		
			      dist = _dxfSegmentIntersection(simp_data->simplified_vertex,
							     simp_data->vert[vstar1[i]], 
							     simp_data->vert[v1], 
							     simp_data->vert[vstar1[j]],
							     &lambda1, &lambda2, !j, 0, 1);
					       
			      			      
			      if (lambda2 >= 0. && lambda2 <= 1.) /* inside 2 */
				{

				  dist += 
				    lambda2        * simp_data->err_volume[v1] +
				    (1. -lambda2)  * simp_data->err_volume[vstar1[j]];


				  if (simp_data->vx_data && lambda1 >= 0. && lambda1 <= 1.)

				    edge_collapsible = 
				      _dxfErrorTestDataOnSurfaceUpdate(simp_data, 2, 2,
				            vstar1[i], 0, lambda1, 1. -lambda1, 0.,
				        v1, vstar1[j], 0, lambda2, 1. -lambda2, 0.);

				  if (lambda1 > 0. && lambda1 <= 1.) /* inside 1 */
		      
				    {
				      potential_err_pivot = 
					(dist - (1. -lambda1)  * simp_data->err_volume[vstar1[i]])
					/ lambda1;

				      
				      UPDATE_ERR_PIVOT;
			      
				    } 
			  
				  else if (lambda1 == 0.) 
				    {
				      /* the error value at vs cannot accomodate that constraint,
					 but maybe it is already satisfied */
				
				      if (simp_data->err_volume[vstar1[i]] < dist)
					edge_collapsible = 0;
				    }
				}

			      j++;
			      if (i == j) j++; /* we prohibit (i == j) */
			    }
			  i++;
			}
		    }
		}
	    }
	}
    }

  /* update the tolerance and error volumes in case the edge will be simplified */

  if (edge_collapsible) 
    {

      simp_data->err_volume[v0] = largest_err_pivot;

      /* update the tolerance volume */
 	  
      simp_data->tol_volume[v0] =  MIN(simp_data->tol_volume[v0], simp_data->tol_volume[v1]);
	
      _dxfErrorTestDataOnSurfaceEnd(simp_data, v0);

    }
  
  return edge_collapsible;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfRemoveEdgesFromHeap                                                |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfRemoveEdgesFromHeap(int *estar, u_short val, int *edge2index, 
			    SimpData *simp_data)
{
  int 
    num_removed = 0;

  register u_short i;
 
  /* remove all edges of the star */
 
  for(i=0;i<val;i++)
    {
      int the_edge = estar[i];
      int idx      = edge2index[the_edge];
      num_removed += simp_data->remove_edge(simp_data, idx, the_edge);
    }

  return num_removed;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfReinstateEdgesInHeap                                               |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int  _dxfReinstateEdgesInHeap(int *estar, u_short val, SimpData *simp_data, int *num_edg_added)

     /* return the number of edges that effectively have been reinstated
	(were in the queue before, were deleted and are now reintered */
{
  
  int 
    num_edg_reinstated = 0;

  register int i;

  *num_edg_added = 0;

  /* add all edges of the star to the heap  */
 
  for(i=0;i<val;i++)
    {
      int
	outside = (simp_data->edge2index[estar[i]] == simp_data->heap_outside_index),
	added    = simp_data->add_edge(simp_data, estar[i]);
		 
      if (added) 
	{
	  *num_edg_added += 1;
	  if (outside) num_edg_reinstated++;
	  
	}
    }

  return num_edg_reinstated;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfCollapseBoundaryEdge2ndKnd                                         |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfCollapseBoundaryEdge2ndKnd(SimpData *simp_data, 
				   int edg_num, int v0, int v1, int val0, int val1, 
				   int *star0, int *star1, int *vstar0, int *vstar1, int *estar0, int *estar1,
				   Vertex *s_normal0, Vertex *s_normal1, float *s_area0, float *s_area1)
{
  int 
    n_edges_added[2],
    n_edges_reinstated = 0;
 

    /*   e*:
	          v*1\      /v*0  \      /
	           e*1\  *1/e*0    \ *1 /
	       *val0-1 \  / *0   *0 \  / *val1-1
	     ___________\/___________\/__________
	     v*val0-1   v0    vs      v1
             vl                                  vr
	     */


  /* father relationships:

                    /\
              e*0[0]  \e*1[0]
                  /    \
                 /<-----\		  
                /     ---\--> *1[1]
		   *1[0]
            v0  <-------  v1


	    This is only true if val1 > 1
	    otherwise if val1 ==1 we must have val0 > 1
     */

  /* 1) update the vertex, triangle and edge fathers */

  simp_data->vertex_father[v1]         = v0;

  if (val1 > 1)
    {
      simp_data->edge_father   [estar1[0]] = simp_data->edge_father   [estar0[0]];
      
      simp_data->triangle_father[star1[0]] = simp_data->triangle_father[star1[1]];
    }
  else
    {
      simp_data->edge_father   [estar0[0]] = simp_data->edge_father   [estar1[0]];
      
      simp_data->triangle_father[star0[0]] = simp_data->triangle_father[star0[1]];
    }

  /* 2) change the valences of the vertices */

  simp_data->valence[v0]               = val0 + val1 - 2;
  simp_data->valence[vstar0[0]] --;

  /* 3) remove from the edge heap all edges that have changed   */
  
  _dxfRemoveEdgesFromHeap(estar0, val0, simp_data->edge2index, simp_data);

  _dxfRemoveEdgesFromHeap(estar1, val1, simp_data->edge2index, simp_data);

  /* mark the 2 discarded edges. 
     
     verify that estar1[0] was not already counted as outside the heap */

  if (simp_data->edge2index[estar1[0]] != simp_data->heap_outside_index)
    simp_data->num_edg_remaining_4_testing --;
  
  simp_data->edge2index[estar1[0]] = -4;
  simp_data->edge2index[edg_num]   = -7;
  

  /* 5) reinstate edges in the heap */

  n_edges_reinstated = 
    _dxfReinstateEdgesInHeap(estar0, val0, simp_data, n_edges_added);
  
  /* reinstate all but the first edge (that collapses) in the star of v1 */

  n_edges_reinstated += 
    _dxfReinstateEdgesInHeap(estar1+1, val1-1, simp_data, n_edges_added+1);
 
  /* 6) update the number of edges left to be tested */

  simp_data->num_edg_remaining_4_testing += n_edges_reinstated;

  /* 7) update the coordinates of the simplified vertex */

  memcpy(simp_data->vert[v0], simp_data->simplified_vertex, sizeof(Vertex));

  /* replace the triangle normals, areas, and compactness in the simplified stars */

  {
    register int i;

    float 
      *s_comp0 = s_area0 + val0,
      *s_comp1 = s_area1 + val1;

    /* replace for v0 */
    for(i=1;i<val0;i++)
      {
	simp_data->area[star0[i]]         = s_area0[i];
	simp_data->compactness[star0[i]]  = s_comp0[i];
	memcpy(simp_data->normal[star0[i]], s_normal0[i], sizeof(Vertex));
      }
    /* replace for v1 */

    for(i=1;i<val1;i++)
      {
	simp_data->area[star1[i]]         = s_area1[i];
	simp_data->compactness[star1[i]]  = s_comp1[i];
	memcpy(simp_data->normal[star1[i]], s_normal1[i], sizeof(Vertex));
      }
  }
  


  return n_edges_added[0] + n_edges_added[1];
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfCollapsibilityTestsBoundaryEdge2ndKind                             |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfCollapsibilityTestsBoundaryEdge2ndKind(SimpData *simp_data, 
					       int edge_num, int v0, int v1, int val0, int val1, 
					       int move_simplified_vertex)
{
  int 
    collapsible = 0,
    *edge_t     = (int *) SurfEdgeGetTriangles (simp_data->edge_array+edge_num),
    *star0      = NULL;

  Vertex *s_normal0 = NULL;

  float *s_area0    = NULL;

  /* perform collapsibility tests particular to a boundary edge of the
     second kind, such that both edge vertices are boundary vertices

     In addition of the two vertices touching a boundary, the edge must
     effectively belong to the boundary, which means that exactly one
     triangle must be missing, with an index of -1. If both triangles had
     a -1 index, though, we would not have a valid edge */

  if (simp_data->boundary_vert[v0] && simp_data->boundary_vert[v1] &&
      (edge_t[0] == -1 || edge_t[1] == -1))
    {
 
      /* we only accept the edges such that both end vertices are located on the
	 boundary */
      
      int 
	first_triangle = edge_t[0];
      
      /* perform the valence test. For this particular case, the valence
	 test is 1 <= v(v0) + v(v1) -2 <= valence_max */

      if ((val0 + val1 -2 >= 1) && (val0 + val1 -2 < simp_data->valence_max))
	{
	 
	  /* build the edge star:
	     contrary to the procedure applied to the surface interior vertices,
	     we compute the vertex star in a directed fashion.
	     Here is the labeling system:

	          v*1\      /v*0  \      /
	           e*1\  *1/e*0    \ *1 /
	       *val0-1 \  / *0   *0 \  / *val1-1
	     ___________\/___________\/__________
	     v*val0-1   v0           v1

	     */

	  /* allocate storage space for the vertex stars */

	  star0 = (int *) simp_data->vertex_star_buffer;

	  dx_bzero (star0, (val0 + val1) * 3 * sizeof (int));

	  /*
	    
	    we next rotate the edge so that edge_t[0] is a valid triangle,     t[0]
	    i.e., so that the representation above holds true              v[0]-----v[1]
                                                                                t[1]
										*/  
	  if (first_triangle == -1)
	    {
	      /* rotate the edge, i.e. echange v0, v1, val0 and val1 
		 using first_triangle as a temporary storage location
		 
		 assign edge_t[1] to first_triangle
		 */

	      first_triangle = v0;
	      v0             = v1;
	      v1             = first_triangle;
	      
	      first_triangle = val0;
	      val0           = val1;
	      val1           = first_triangle;

	      first_triangle = edge_t[1];
	    }

	  
	  {int 
	      *vstar0        = star0  + val0,
	      *estar0        = vstar0 + val0,
	      *star1         = estar0 + val0,
	      *vstar1        = star1 + val1,
	      *estar1        = vstar1 + val1;
	      
	    first_triangle = FATHER(first_triangle, simp_data->triangle_father);

	    _dxfDirectedParentVertexStar(v0, first_triangle, val0, star0, simp_data, 
					 DIRECTION_COUNTER_CLOCKWISE);

	    /* compute the star of v1 in reverse order */

	    _dxfDirectedParentVertexStar(v1, first_triangle, val1, star1, simp_data, DIRECTION_CLOCKWISE);

	   
	    /* 
	       up to now, we found a parent triangle 
	       and we located the next edge, "n_edgeCCW" or "n_edgeCW" 
	       depending on whether we are rotating counterclockwise
	       or clockwise around "parent vertex".
	       
	       The parent of this edge must be stored in the e_star. 

	                              _/                                                          
	                    n_edgeCCW/\
		   parent vertex    /  \  
	                        \_ /____\                                                         
	                           n_edge\_                                                      
	                             CW       

	       Similarly to before, we need to orient the edge so that 
	       the parent of the first edge vertex is precisely "parent vertex"
	       */  
	
	    /* the manifold test is very similar to the test we
	       perform for an interior edge case we compare the
	       vertices of the links and determine whether duplicate
	       vertices exist */
	      
	    if (_dxfManifoldLinkTest(vstar0+1, vstar1+1, val0-1, val1-1))
	      {
		  Vertex 
		    *s_normal1;
		  float 
		    *s_area1 ;

		  s_normal0 = (Vertex *) simp_data->edge_star_buffer_vx;

		  dx_bzero (s_normal0, ( val0 + val1 ) * sizeof (Vertex));
		  
		  s_normal1 = s_normal0 + val0;
		       
		  s_area0 = (float *) simp_data->edge_star_buffer_fl;

		  dx_bzero (s_area0, 2 * ( val0 + val1 ) * sizeof (float));

		  s_area1 = s_area0 + 2 * val0;

		  /* position the simplified vertex so that the length of the boundary
		     will be preserved.  

		     The vertex must be situated on an ellipse */

		  if ((_dxfSimplifiedBoundaryVertexLocation(simp_data, v0, v1, vstar0, vstar1, 
							    val0, val1, move_simplified_vertex) )
		      &&
			    
		      /* compute the difference in orientation for the
			 triangle normals before and after simplification */
		      
		      /* compute the minimum compactnesses before and after simplification */

		      (_dxfBoundaryCollapse2KndGeometricallyOk(simp_data, v0, v1, star0, vstar0, val0, 
							       s_normal0, s_area0,
							       DIRECTION_COUNTER_CLOCKWISE,
							       simp_data->min_scalprod_nor, 
							       simp_data->compactness_ratio) )
		      &&
		      (_dxfBoundaryCollapse2KndGeometricallyOk(simp_data, v1, v0, star1, vstar1, val1, 
							       s_normal1, s_area1, DIRECTION_CLOCKWISE,
							       simp_data->min_scalprod_nor, 
							       simp_data->compactness_ratio) )
		      ) /* geometry test */


		    {
		      /* compute the tolerance after simplification */
			    
		      if (_dxfErrorWithinToleranceVBoundary2ndKnd(simp_data, v0, v1, vstar0, vstar1,
								  val0, val1))
			{
				
			  simp_data->num_edg_collapsed +=2;

			  simp_data->num_edg_weights += 

			    _dxfCollapseBoundaryEdge2ndKnd(simp_data, edge_num, v0, v1, val0, val1,
							   star0,  star1,  vstar0, vstar1, 
							   estar0, estar1, s_normal0, s_normal1,
							   s_area0, s_area1); 

			  collapsible = 1;

			}
		      else /* reject for tolerance off limits */
			simp_data->edg_rejected_4_tolerance ++;
			   
		    }
		  else /* reject for geometry */
		    simp_data->edg_rejected_4_geometry ++;
	      }
	    else /* reject for topology */
	    
	      {
		simp_data->edge2index[edge_num] = -3; 
			
		/* this condition will not be modified even if the
		   surface is simplified around the edge */

		simp_data->edg_rejected_4_topology ++;
	      }
	  }
	}
      else   
	simp_data->edg_rejected_4_topology ++;

    } /* end if (edge_t[0] == -1 || edge_t[1] == -1) */

  return collapsible;
      
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfManifoldLinkTest1stKnd                                             |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfManifoldLinkTest1stKnd(int v0, int val0, int *star1, int *link1, int val1, SimpData *simp_data)
{

  int link_is_manifold = 1;

  /* 
     v0 is a  boundary vertex
     v1 is an interior vertex


     | /\vt
     |/  \<-----this triangle is *1[0]
     v0___\_v1
     |\   /
     | \ /<-----this triangle is *1[1]
     |  vb


     this test works by taking 
     all vertices in the link of v0, except v1, and verify that they are
     all different from the vertices of the link of v1, except vt, v0 and vb
     
     we take from vstar0 as many vertices as val0 (remember that vstar0 has one element
     more than star0, because v0 is a boundary vertex)

     we skip the vertices 0, 1, and val1-1 from the link of v1, i.e.,
     we take from vstar1 as many vertices as val1-3
     
   */
  
  int 
    val = val0 + val1 -3;

  /* since 
	val0 + val1 <= SIMPLIFY_VALENCE_MAX4,
	I can use a static array here. */

  static int
    link[SIMPLIFY_VALENCE_MAX4];

    {
    
      int 
	i  = 0,
	
	tri_v_rotated[3],  tri_v_parents[3], 
    
	/* rotated edge vertices and adjacent triangles: */
	edge_v_rotated[2], edge_v_parents[2], edge_t_rotated[2],
    
	triangle = star1[0],
	edge;

      do
	{
	  triangle = FATHER (triangle,  simp_data->triangle_father);

	  _dxfRotateParentTriangle(v0, simp_data->vertex_father,
				   simp_data->triangles + 3 * triangle,
				   tri_v_parents, tri_v_rotated);

	  edge = _dxfNextParentEdge(tri_v_rotated, DIRECTION_COUNTER_CLOCKWISE);

	  _dxfRotateParentEdge(v0, simp_data->vertex_father,
			       (int *) SurfEdgeGetVertices (simp_data->edge_array + edge),
			       (int *) SurfEdgeGetTriangles(simp_data->edge_array + edge),
			       edge_v_parents, edge_v_rotated, edge_t_rotated);

	  link[i++] = edge_v_parents[1];

	  triangle = edge_t_rotated[0];

	}
      while (triangle != -1);


      /* now repeat the same operation in the clockwise direction
	 starting at vb (for Vertex_at_Bottom) */

      triangle = star1[1];

      do
	{
	  triangle = FATHER (triangle, simp_data->triangle_father);

	  _dxfRotateParentTriangle(v0, simp_data->vertex_father,
				   simp_data->triangles + 3 * triangle,
				   tri_v_parents, tri_v_rotated);

	  edge = _dxfNextParentEdge(tri_v_rotated, DIRECTION_CLOCKWISE);

	  _dxfRotateParentEdge(v0, simp_data->vertex_father,
			       (int *) SurfEdgeGetVertices (simp_data->edge_array + edge),
			       (int *) SurfEdgeGetTriangles(simp_data->edge_array + edge),
			       edge_v_parents, edge_v_rotated, edge_t_rotated);

	  link[i++] = edge_v_parents[1];

	  triangle = edge_t_rotated[1];
	}
      while (triangle != -1);


      /* gather vertices of link0 and link1 in a single array, sort
	 the array, and determine whether there are any duplicate
	 vertices */

      memcpy(&link[i], link1+2, (val1 -3) * sizeof (int));

      qsort(link, val, sizeof(int), (int (*)(const void*, const void*)) _dxfCmpIntSmallFirst);

      i = 0;

      while ((i < val -1 ) && (link_is_manifold))
	{
	
  
	  if (link[i+1] == link[i]) /* if there are duplicate vertices,
				       they must be contiguous in the
				       sorted array */
				       
	    link_is_manifold = 0;
	  i++;
	}

    }

  return (link_is_manifold);
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfBoundaryCollapse1stKndGeometricallyOk                              |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfBoundaryCollapse1stKndGeometricallyOk(SimpData *simp_data, int v0, int *star0, int *vstar0, 
					      int val0, Vertex *s_normal0, float *s_area0,
					      float min_scalprod_nor, float compactness_ratio) 

{ 
  int collapsible = 1, i = 1;

  Vertex faceVert[3];
  
  float 
    *s_comp0 = s_area0 + val0,
    min_compactness_before,
    min_compactness_after;

  /*
    labeling system          *[0]\  /
                         v*[0]____\/____
			      *[1]/\
                                 /  \
				 *[2]

    both triangles *[0] and *[1] disappear after edge collapse

   */


  /* compute the compactness of the collapsed triangles */
  
  min_compactness_before = MIN ( simp_data->compactness[star0[1]], simp_data->compactness[star0[0]]);

  min_compactness_after  = 1; /* initialization of the compactness
				 after simplification */
  
  while (collapsible && i < val0-1)
    {
 
      /* retrieve the face compactness before simplification */
      min_compactness_before = MIN ( simp_data->compactness[star0[i+1]], min_compactness_before);

      /* recompute the compactness, and area, after changing the first vertex */
  
      memcpy(faceVert[0], simp_data->simplified_vertex, sizeof (Vertex));
      memcpy(faceVert[1], simp_data->vert[vstar0[i]], sizeof (Vertex));
      memcpy(faceVert[2], simp_data->vert[vstar0[i+1]], sizeof (Vertex));

      s_comp0[i+1] = _dxfFastCompactness3DTriangle2(faceVert, s_area0+i+1);

      min_compactness_after = MIN ( s_comp0[i+1], min_compactness_after);

      /* recompute the face normal assuming that the vertex v0 has
         moved and save this normal for future use*/

      _dxfTriangleNormalQR2(faceVert, s_normal0[i+1]);

      /* check whether the triangle normal orientation is consistent */

      collapsible = (SCALPROD3(s_normal0[i+1], simp_data->normal[star0[i+1]]) > min_scalprod_nor);

      i++;
    }

  if (collapsible)

    /* check whether the minimum triangle compactness is too much degraded */

    collapsible = min_compactness_after > (compactness_ratio * min_compactness_before);
  
  return collapsible;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfMakeBarycenter                                                     |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

void _dxfMakeBarycenter(int nV, Vertex *vv, Vertex bary, float *bary_coord)
{
  register int i;
  register u_char j;
 
  /* initialize barycenter */

  for(j=0;j<3;j++)
    bary[j]=0.0;

  for(i=0;i<nV;i++)
    {
      for(j=0;j<3;j++)
	bary[j]+= bary_coord[i]*vv[i][j];
    }

  return;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfSolveSymmetric2x2Eqn                                               |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfSolveSymmetric2x2Eqn(double a,  double b, double d, double e, double f, double *x, double *y)
{
  int retval = 0;

  static double eps = 1e-13;
  
  if (fabs(b) < eps)
    {
      *x = (fabs(a) < eps) ? 0. : e / a;
      *y = (fabs(d) < eps) ? 0. : f / d;

      retval = 2;
    }
  else
    {
      double 
	theta = (d - a) / 2. / b,
	sign  = (theta > 0.) ? 1. : -1,
	t     = sign / ( theta * sign + sqrt (theta * theta + 1.0)),
	/* we do not check for overflows of theta *theta */
	c     = 1 / (sqrt (1.0 + t*t)),
        s     = t * c,

        /*    |  a    b  |     |  c    s  | |  A       |  |  c   -s  | 
              |          |  =  |          | |          |  |          |
              |  b    d  |     | -s    c  | |       C  |  |  s    c  |
        */
	
	A    = a - t * b,
	C    = d + t * b,
	Rx   = (fabs(A) < eps) ? 0. : (c * e - s * f) / A,
	Ry   = (fabs(C) < eps) ? 0. : (s * e + c * f) / C;
      
      *x =  c * Rx + s * Ry ; /* apply the inverse rotation */
      *y = -s * Rx + c * Ry ; /*          "    "            */
	
      retval = 1;
    }

  return retval;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfTriangle3DBarycentricCoordinates2                                  |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfTriangle3DBarycentricCoordinates2(Vertex *tri, Vertex w, Vertex wp, Vertex bary, float *residual)

/* compute the barycentric coordinates of a vector with respect to a point
   the projection wp and the residual are only valid if the function
   returns 1 (vertex projects inside)
*/
{
  int inside;

  static Vertex x1, x2;

  static float eps_vertex_outside = EPS_VERTEX_OUTSIDE_TRIANGLE;

  double a = 0., b = 0., c = 0., d = 0., e = 0., x3, x, y;

  register u_char j;

  int origin;

  /*
     x1 is the vector representing the shortest edge of the triangle,
     which is best numerically

   Min 
            
   ||  |  x1[0]  x2[0] |               | w[0] - tri[origin][0] |  ||^2
   ||  |  x1[1]  x2[1] |  *  | b1 | -  | w[1] - tri[origin][1] |  ||  
   ||  |  x1[2]  x2[2] |     | b2 |    | w[2] - tri[origin][2] |  ||  


   using the normal equations 
   */

  _dxfTriangleBasisVectors(tri, x1, x2, origin);

  for (j=0; j<3;j++)
    {
      x3 = w[j] - tri[origin][j];

      a += x1[j] * x1[j];
      b += x1[j] * x2[j];
      c += x2[j] * x2[j];
      d += x1[j] * x3;
      e += x2[j] * x3;
    }
  
  _dxfSolveSymmetric2x2Eqn(a,b,c,d,e, &x, &y);

  bary [origin] = 1. - x - y; 
    
  bary [(origin + 1)%3] = x; 
     
  bary [(origin + 2)%3] = y; 
     
  _dxfMakeBarycenter(3, tri, wp, bary);
      
  /* determine whether w projects inside or outside tri */

  inside = !(bary[0] < eps_vertex_outside || 
	     bary[1] < eps_vertex_outside ||  
	     bary[2] < eps_vertex_outside);

  /* we actually do not need to create a variable eps_vertex_outside
     the following test would work all the same 

  inside = !(bary[0] < EPS_VERTEX_OUTSIDE_TRIANGLE || 
	     bary[1] < EPS_VERTEX_OUTSIDE_TRIANGLE ||  
	     bary[2] < EPS_VERTEX_OUTSIDE_TRIANGLE);

	     */

  *residual   = DIST3(w, wp);

  return inside;

}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfClosestPointOnTriangle                                             |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfClosestPointOnTriangle(Vertex *tri, Vertex w, Vertex wp, Vertex bary, float *dist)
{
  int inside = _dxfTriangle3DBarycentricCoordinates2(tri, w, wp, bary, dist);
	
  if (!inside)
    {
      static float eps= EPS_VERTEX_OUTSIDE_TRIANGLE;
  
      /* either 2 barycentric coordinates are negative */

      if (bary[0] < eps && bary[1] < eps)
	{
	  bary[0] = bary[1] = 0.;
	  bary[2] = 1.;
	  memcpy(wp, tri[2], sizeof(Vertex));
	  *dist = DIST3(w,wp);
	}
	
      else if (bary[0] < eps && bary[2] < eps)
	{
	  bary[0] = bary[2] = 0.;
	  bary[1] = 1.;
	  memcpy(wp, tri[1], sizeof(Vertex));
	  *dist = DIST3(w,wp);
	}

      else if (bary[1] < eps && bary[2] < eps)
	{
	  bary[1] = bary[2] = 0.;
	  bary[0] = 1.;
	  memcpy(wp, tri[0], sizeof(Vertex));
	  *dist = DIST3(w,wp);
	}

      /* or 1 barycentric coordinate is negative */

      else if (bary[0] < eps)
	{
	  *dist =  _dxfClosestPointOnEdge(w, wp, tri[1], tri[2], &bary[1]);
	  bary[2] = 1. - bary[1];
	  bary[0] = 0.;
	} 

      else if (bary[1] < eps)
	{
	  *dist =  _dxfClosestPointOnEdge(w, wp, tri[0], tri[2], &bary[0]);
	  bary[2] = 1. - bary[0];
	  bary[1] = 0.;
	}
      
      else if (bary[2] < eps)
	{
	  *dist =  _dxfClosestPointOnEdge(w, wp, tri[0], tri[1], &bary[0]);
	  bary[1] = 1. - bary[0];
	  bary[2] = 0.;
	}
    }

  return (inside);
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfErrorWithinToleranceVBoundary1stKnd                                |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfErrorWithinToleranceVBoundary1stKnd(SimpData *simp_data, int v0, int v1, int *vstar1, int val1)
{
  int edge_collapsible = 1;

  register int i;
    
  float
    tol_pivot, potential_err_pivot, largest_err_pivot = 0., 
    dist,
    lambda1, lambda2;

  /* set the maximum tolerance for the star */
  
  tol_pivot = simp_data->tol_volume[v1];
  
  for (i=0;i<val1;i++)
    if (simp_data->tol_volume[vstar1[i]] < tol_pivot)
      tol_pivot = simp_data->tol_volume[vstar1[i]];

  /* verify that the maximum tolerance is not already exceeded at
     the vertices */

  i = 0; 
  while (edge_collapsible && i < val1)
    if (simp_data->err_volume[vstar1[i++]] > tol_pivot)
      edge_collapsible = 0;
  
  if (simp_data->err_volume[v1] > tol_pivot)
    edge_collapsible = 0;
  
  if (edge_collapsible)
    {
      Vertex
	triangle[3],
	vs[1],/* local copy of the coordinates of the simplified vertex */
	wp[2],
	alpha[SIMPLIFY_VALENCE_MAX1];  /* val0 + val1 -4 <= valence max <=> val0, val1 <= valence_max +1 
					  because val0 >=3 and val1 >=3 */
      float
	dp;

      int 
	i_min         = -1;

      dx_bzero((char *)triangle, 3 * sizeof (Vertex));

      dx_bzero((char *)vs, sizeof (Vertex));

      dx_bzero((char *)wp, 2 * sizeof (Vertex));

      dx_bzero((char *)alpha, SIMPLIFY_VALENCE_MAX1 * sizeof (Vertex));

      _dxfErrorTestDataOnSurfaceInit(simp_data);

      /* the simplified is already positioned, since it must be v0 */
      memcpy(vs[0], simp_data->simplified_vertex, sizeof (Vertex));

      /* if there is data attached to the surface, we also copy the data of v0
	 to the new potential data */

      if (simp_data->vx_data) memcpy(simp_data->vx_data_potential_values, 
				     simp_data->vx_data + v0*simp_data->data_dim,
				     simp_data->data_dim * sizeof(float));

      /* 1) determine the closest point to the vertex v1 on the triangles of star1 after simplification
       */
      
      memcpy(triangle[0], vs[0], sizeof (Vertex));

      /*                  v*[val1-1]
    labeling system          *[0]\  /
                         v*[0]____\/____
			      *[1]/\
                            v*[1]/  \
				 *[2]

				 both triangles *[0] and *[1] disappear after edge collapse

				 */
      i    = 1;

      dist = 0.;

      while (i<val1-1)
	{
	  memcpy(triangle[1], simp_data->vert[vstar1[i]],   sizeof(Vertex));
	  memcpy(triangle[2], simp_data->vert[vstar1[i+1]], sizeof(Vertex));
	
	  _dxfClosestPointOnTriangle(triangle, simp_data->vert[v1], wp[0], alpha[i], &dp);
      
	  dp += simp_data->err_volume[v1];

	  if (alpha[i][0]>0.)
	    {
	      dp = (dp - 
		alpha[i][1] * simp_data->err_volume[vstar1[i]] -
		alpha[i][2] * simp_data->err_volume[vstar1[i+1]])
		/ alpha[i][0];

	      if (i_min == -1 || dp < dist)  {dist = dp; i_min = i;}
		
	    }
	  else
	    {
	      /* otherwise v1 projects onto the link edge of "triangle"
		 the epsilon value at the simplified vertex cannot compensate
		 but it is possible that the error values at the link are already sufficient */

	      if (alpha[i][1] * simp_data->err_volume[vstar1[i]]   + 
		  alpha[i][2] * simp_data->err_volume[vstar1[i+1]] >= dp)
		{
		  /* look no further for the closest triangle */
		  i_min = i;
		  i = val1-1; /* so that the loop will immediately stop*/
		  dist = 0;
		}
	    }
	 
	  i++;
	}
	  
      /* the closest triangle wins */
      
      if (i_min >=0)
	       
	{
	  potential_err_pivot = dist;

	  UPDATE_ERR_PIVOT;	
   
	  if (edge_collapsible && simp_data->vx_data)

	    edge_collapsible = 
	      _dxfErrorTestDataOnSurfaceUpdate(simp_data, 3, 1,
		       vstar1[i_min], vstar1[i_min+1], 
		       alpha[i_min][0], alpha[i_min][1], alpha[i_min][2],
		       v1, 0, 0, 1., 0., 0.);

	}

      else 
	edge_collapsible = 0; /* v1 did not project inside of the new star */


      if (edge_collapsible)
	{
	  /* 2)
	     
                /                 v*[val1-1]     v*[val1-2]
               /    __/                    \    /
	      /  __/                        \  /
	   v0/__/________                ____\/____
           ==\  \__                     v*0  /\
           vs \    \__                      /  \
               \      \                    /    \
                \                         v*1    v*2

		take in turn each edge (vs,v*2) ..., (vs, v*[val1-2])
		
		and "intersect" them with each edge (v1, v*[1]) ... (v1, v*[val1-1])
		*/
	  int j;

	  i = 2;

	  while (edge_collapsible && i < val1-1)
	    {
	      j = 1;
	      
	      while (edge_collapsible && j < val1)
		{         
		  dist =  _dxfSegmentIntersection(simp_data->simplified_vertex,
						  simp_data->vert[vstar1[i]], 
						  simp_data->vert[v1], 
						  simp_data->vert[vstar1[j]],
						  &lambda1, &lambda2, (j==1), 0, 1);
		     
		  /* If the closest point is located in the
		     inside of either edge, determine a crossing point. */

		  if (lambda2 >= 0. && lambda2 <= 1.) /* inside 2 */
		    {
		      dist += 
			lambda2        * simp_data->err_volume[v1] +
		        (1. -lambda2)  * simp_data->err_volume[vstar1[j]];


		      if (simp_data->vx_data && lambda1 >= 0. && lambda1 <= 1.)
			
			edge_collapsible = 
			  _dxfErrorTestDataOnSurfaceUpdate(simp_data, 2, 2,
				            vstar1[i], 0, lambda1, 1. -lambda1, 0.,
				        v1, vstar1[j], 0, lambda2, 1. -lambda2, 0.);


		      if (lambda1 > 0. && lambda1 <= 1.) /* inside 1 */
		      
			{
			  potential_err_pivot = 
			    (dist + (lambda1 -1.) * simp_data->err_volume[vstar1[i]])
			    / lambda1;
	    
			  UPDATE_ERR_PIVOT;
			}
		      else if (lambda1 == 0.)
			{
			  /* the error at vs cannot accomodate that constraint,
			     but maybe it is already satisfied */

			  if (simp_data->err_volume[vstar1[i]] < dist)
			    edge_collapsible = 0;
			}
		    }
		  j++;

		  if (i == j) j++; /* we prohibit (i == j) */
		
		}
	      i++;
	    }
	}
	
    }

  if (edge_collapsible)
    {
      /* update the tolerance and error values */

      /* only the error at the pivot will change */
      
      simp_data->err_volume[v0] = largest_err_pivot;

      /* update the tolerance volume */
 	  
      simp_data->tol_volume[v0] =  MIN(simp_data->tol_volume[v0], simp_data->tol_volume[v1]);	

      _dxfErrorTestDataOnSurfaceEnd(simp_data, v0);

    }

  return edge_collapsible;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfCollapseBoundaryEdge1stKnd                                         |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfCollapseBoundaryEdge1stKnd(SimpData *simp_data, int edg_num, int v0, int v1, int val0, int val1, 
				   int *star1, int *vstar1, int *estar1,
				   Vertex *s_normal1, float *s_area1)
{
  int 
    n_edges_added      = 0,
    n_edges_reinstated = 0,
    et, eb,
    tri_v_rotated[3],  tri_v_parents[3];
 
   /*                     v*[val1-1]
       labeling system          *[0]\  /
                            v*[0]____\/____
	   		         *[1]/\
                               v*[1]/  \
				    *[2]

				 both triangles *[0] and *[1] disappear after edge collapse

   */


  /* father relationships:

                    /\
                et /  \e*1[val1-1]
                  /    \
                 /<-----\		  
                /     ---\--> *1[val1-1]
		   *1[0]
            v0  <-------  v1
               \   *1[1]/
            eb  \    --/------> *1[2]
                 \<---/e*1[1]
                  \  /
                   \/
     */

  /* 1) update the vertex, triangle and edge fathers */

  /* 1.1) First retrieve the edges et and eb on "top" and "bottom" of
     v0, since these edges were not pre-stored */

  /* take the triangle star1[0], orient it with respect to v0
     and use the oriented triangle vertices to determine the next edge counter clockwise starting
     from v0: et */

  _dxfRotateParentTriangle(v0, simp_data->vertex_father, simp_data->triangles + 3 * star1[0],
			   tri_v_parents, tri_v_rotated);

  et =  (int) (((EdgeS *) _dxfFindEdge(simp_data->e_table, tri_v_rotated[0], tri_v_rotated[2])) - 
	       simp_data->edge_array);

  et =  FATHER (et, simp_data->edge_father);

  /* do the same with the triangle star1[1] to determine the next edge clockwise starting
     from v0: eb */

  _dxfRotateParentTriangle(v0, simp_data->vertex_father, simp_data->triangles + 3 * star1[1],
			   tri_v_parents, tri_v_rotated);

  eb =  (int) (((EdgeS *) _dxfFindEdge(simp_data->e_table, tri_v_rotated[0], tri_v_rotated[1])) - 
	       simp_data->edge_array);

  eb =  FATHER (eb, simp_data->edge_father);

  /* 1.2) */

  simp_data->vertex_father[v1]   = v0;

  simp_data->edge_father   [estar1[0]]      = simp_data->edge_father   [estar1[1]] = eb;
  simp_data->edge_father   [estar1[val1-1]] = et;

  simp_data->triangle_father[star1[0]] = star1[val1-1];
  simp_data->triangle_father[star1[1]] = star1[2];
   
  /* 2) change the valences of the vertices */

  simp_data->valence[v0]               = val0 + val1 - 4;
  simp_data->valence[vstar1[1]]      --;
  simp_data->valence[vstar1[val1-1]] --;

  /* 3) remove from the edge heap all edges that have changed   */
  
  _dxfRemoveEdgesFromHeap(estar1, val1, simp_data->edge2index, simp_data);

  /* mark the 2 discarded edges. 
     
     verify that they were  not already counted as outside the heap */

  if (simp_data->edge2index[estar1[1]] != simp_data->heap_outside_index)
    simp_data->num_edg_remaining_4_testing --;
   
  if (simp_data->edge2index[estar1[val1-1]] != simp_data->heap_outside_index)
    simp_data->num_edg_remaining_4_testing --;
  
  simp_data->edge2index[estar1[1]] = simp_data->edge2index[estar1[val1-1]] = -4;
  simp_data->edge2index[edg_num]   = -7;
  

  /* 5) reinstate all but the two edges and last edge  in the star of v1 */

  n_edges_reinstated += 
    _dxfReinstateEdgesInHeap(estar1+2, val1-3, simp_data, &n_edges_added);
 
  /* 6) update the number of edges left to be tested */

  simp_data->num_edg_remaining_4_testing += n_edges_reinstated;

  /* 7) update the coordinates of the simplified vertex" no need since
   its coordinates do not change*/

  /* replace the triangle normals, areas and compactness in the simplified star of v1 */

  {
    register int i;

    float *s_comp1 = s_area1 + val1;
    
    /* replace for v1 */

    for(i=2;i<val1;i++)
      {
	simp_data->area[star1[i]]         = s_area1[i];
	simp_data->compactness[star1[i]]  = s_comp1[i];
	memcpy(simp_data->normal[star1[i]], s_normal1[i], sizeof(Vertex));
      }
  }

  return n_edges_added;
}


/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfCollapsibilityTestsBoundaryEdge1stKind                             |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfCollapsibilityTestsBoundaryEdge1stKind(SimpData *simp_data, 
					       int edge_num, int v0, int v1, int val0, int val1)
{
   /* perform collapsibility tests particular to a boundary edge of
     the first kind, such that exactly one of the edge vertices is a
     boundary vertex */

  int 
    collapsible = 0;

  int 
    *star1 = NULL;

  Vertex 
    *s_normal1 = NULL;

  float 
    *s_area1 = NULL;

  /* the valence of v0 becomes val0 + val1 - 4 */
  
  if ((val0 + val1 -4 >= 1) && (val0 + val1 -4 < simp_data->valence_max))
    {
      int
	*edge_t     = (int *) SurfEdgeGetTriangles (simp_data->edge_array+edge_num),
	first_tri   = edge_t[0],
	second_tri  = edge_t[1],
	/* this is the counter-clockwise order with respect to v1 */
	*vstar1, 
	*estar1;

      /* permute v0 and v1  so that v0 is the boundary vertex */

      if (simp_data->boundary_vert[v1])
	{
	  int 
	    tmp = v0;
	  v0 = v1;
	  v1 = tmp;

	  tmp = val0; val0 = val1; val1 = tmp;
	  first_tri = edge_t[1];      second_tri = edge_t[0];
	}
 
      /* initialize the storage pointer for the buffer star1 */

      star1 = (int *) simp_data->vertex_star_buffer;

      dx_bzero (star1, 3 * val1 * sizeof (int));

      /*
                   | / first_triangle
                   |/
                   v0----------------v1
                    \
                     \
                   boundary
                   */

    
      estar1 = star1  + val1;
      vstar1 = estar1 + val1;

      first_tri = FATHER(first_tri, simp_data->triangle_father);

      /* first star elements (index 0)*/
      star1[0]  = first_tri;
      estar1[0] = edge_num;
      vstar1[0] = v0;
	  
      /* complete the star from indices 1 to val1-1 */
      _dxfParentVertexStar(v1, second_tri, (u_short) val1, star1, simp_data); 

      if (_dxfManifoldLinkTest1stKnd(v0, val0, star1, vstar1, val1, simp_data))
	{

	  /* the simplified vertex is already positioned, since it must be v0 */
	  
	  memcpy(simp_data->simplified_vertex, simp_data->vert[v0], sizeof (Vertex));

	  /* initialize the storage pointer for the buffer s_normal1 */

	  s_normal1 = (Vertex *) simp_data->edge_star_buffer_vx;

	  dx_bzero(s_normal1, val1 * sizeof (Vertex));


	  /* initialize the storage pointer for the buffer s_area1 */

	  s_area1 = (float *) simp_data->edge_star_buffer_fl;

	  dx_bzero (s_area1, 2 * val1 * sizeof (float));

		  
	  if (_dxfBoundaryCollapse1stKndGeometricallyOk(simp_data, v1, star1, vstar1, val1, 
							s_normal1, s_area1, simp_data->min_scalprod_nor, 
							simp_data->compactness_ratio))

	    {
	      /* perform the error tolerance test */
		      
	      if (_dxfErrorWithinToleranceVBoundary1stKnd(simp_data, v0, v1, vstar1, val1))
		{
		  simp_data->num_edg_collapsed +=3;
		  
		  simp_data->num_edg_weights += 

		    _dxfCollapseBoundaryEdge1stKnd(simp_data, edge_num, v0, v1, val0, val1, 
						   star1, vstar1, estar1, s_normal1, s_area1);
			 
		  collapsible = 1;

		}
	      else /* reject for tolerance off limits */
		simp_data->edg_rejected_4_tolerance ++;
	    }

	  else /* reject for geometry */
	    simp_data->edg_rejected_4_geometry ++;

	}
      else
	{ 
	  /* this condition will not change if the surface is
	     further simplified */

	  simp_data->edge2index[edge_num] = -3; 
			
	  simp_data->edg_rejected_4_topology ++;
	}
    }
  else   
    simp_data->edg_rejected_4_topology ++;

  return collapsible;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfCollapseTopologicallyFeasible                                      |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int  _dxfCollapseTopologicallyFeasible(SimpData *simp_data, int *vstar0, int *vstar1, 
				       u_short val0, u_short val1)
{

  register u_short i, j;
  int feasible = 1;

  /* count each vertex only once 
                        *1[val1-1]  /\ *0[1]  counted twice
                           *1[0]   /  \ *0[0]
                                   \  /
			      *1[1] \/ *0[val0-1] counted twice
   */

  {
    /* verify that no two vertices of the link are the same */

    val0 -= 2;
    val1 -= 2;

    i = 2;
    while(( i <= val0) && (feasible))
      {

	j = 2;
	while ((j <= val1) && (feasible))
	  {

	    if (vstar0[i] == vstar1[j])
	      feasible = 0;
	    j++;
	  }
	i++;
      }
  }

  return(feasible);

}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfBuildParentEdgeStars                                               |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

void _dxfBuildParentEdgeStars(EdgeS *edg, int v0, int v1, u_short val0,
			      u_short val1, int *star0, int *star1, SimpData *simp_data)


/*...........................................................................*/
/*
   The following procedure computes the star of v0 and v1 belonging to
   an "active" surface data structure, i.e. the surface is undergoing a
   simplification process.

   the ordering of the stars is done such that:

         e2    e1                     e5    e4
          \ t2 /  t1                  \ t5 /
      t3   \  /                   t0   \  /  t4
  e3 _______\/_______e0       e0 _______\/_______e3
            /\v0     v1        v0       /\v1    
      t4   /  \  t0               t1   /  \  t3
          / t5 \                      / t2 \
         e4     e5                   e1     e2
*/
/*...........................................................................*/

{
  int *estar0 = star0  + val0, *vstar0 = estar0 + val0;
  int *estar1 = star1  + val1, *vstar1 = estar1 + val1,
      e_num   = (int) (edg - simp_data->edge_array);

  int *edg_triangles = SurfEdgeGetTriangles(edg);
  /* complete first the straighforward portion of the stars
         /\
        /  \
       /    \    father(edg_triangles[0]);
      /      \           
   v0/edg->num\ v1         
     \        /          
      \      /           
       \    / father(edg_triangles[1]);
        \  / 
	 \/ 
	 */

  star0 [0] = FATHER(edg_triangles[1], simp_data->triangle_father);
  vstar0[0] = v1;
  estar0[0] = e_num; /* the simplifiable edge is always a root edge */

  _dxfParentVertexStar(v0, edg_triangles[0], val0, star0, simp_data);
 
  star1 [0] = FATHER(edg_triangles[0], simp_data->triangle_father);
  vstar1[0] = v0;
  estar1[0] = e_num;
  
  _dxfParentVertexStar(v1, edg_triangles[1], val1, star1, simp_data);
 
  return;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfParentVertexStar                                                   |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

void _dxfParentVertexStar(int vf, int t0, u_short val, int *star, SimpData *simp_data)

/*...........................................................................*/
/*
   an "active" vertex is a vertex belonging to the inside of a surface
   ( not a boundary vertex ) undergoing a simplification process
*/
/*...........................................................................*/
{
  int *estar = star  + val;
  int *vstar = estar + val;
  register int i;
  int t = t0;
  
  /* complete the star from 1 to the valence,
     the 0-elements have been already completed */

  for(i=1;i<val;i++)
    {
      register u_char j=0;
      int vA,vAf,vB,vBf,ef, e_num, *endpoints, *facing_tris;
      EdgeS *edg, *edgf;

      t = FATHER(t, simp_data->triangle_father);

      /* complete the triangle part of the vertex star */
      star[i] = t;
      do 
	{
	  vA = simp_data->triangles[3 * t + j];
	  vB = simp_data->triangles[3 * t + (j+2)%3];
	  vAf  = FATHER(vA, simp_data->vertex_father);
	}
      while((j++ < 3) && (vAf != vf));
       
      /* now complete the edge portion of the vertex star */

      edg   = _dxfFindEdge(simp_data->e_table,vA,vB);

      e_num = (int) (edg - simp_data->edge_array); /* differential pointer is
					also the edge number */

      ef    = FATHER (e_num, simp_data->edge_father);

      estar[i] = ef;

      /* find out the parents of the vertices of ef */

      edgf = simp_data->edge_array+ef;

      endpoints   = SurfEdgeGetVertices(edgf);
      facing_tris = SurfEdgeGetTriangles(edgf);

      vAf  = FATHER(endpoints[0], simp_data->vertex_father);
      vBf  = FATHER(endpoints[1], simp_data->vertex_father);
      
      /* and the vertex portion of the vertex star */

      if (vAf == vf) 
	{
	  vstar[i] = vBf;
	  /*                          /\f0  ^ccw direction
             the next triangle is  v0/__\v1 |          v0 = vf, v0 > v1 
                                     \  /   |
                                      \/f1   */
 
	  t        = facing_tris[0];
	}
      else
	{
	  vstar[i] = vAf;
	  t        = facing_tris[1];
	 }
    
    }

  return;
}

static void  CentroidEdgeLink(Vertex centroid, Vertex *vert, int val0, int val02, int val014);

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   CentroidEdgeLink                                                       |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

static void  CentroidEdgeLink(Vertex centroid, Vertex *vert, int val0, int val02, int val014)
{
  register int i;
  register u_char j;

  dx_bzero((char *)centroid, sizeof(Vertex));

  for(j=0;j<3;j++)
    {
      for(i=0; i<=val02; i++)
   
	centroid[j] += vert[i][j];
      
      for(i=val0; i<=val014; i++)

	centroid[j] += vert[i][j];
      
      centroid[j] /= val014; 
    }

  return;
}

static void  BuildEdgeStarVal1Is3(Vertex *, Face *, int, int, int, int, int, int *, Vertex *);

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   BuildEdgeStarVal1Is3                                                   |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

static void  BuildEdgeStarVal1Is3(Vertex *vert, Face *face, 
				  int v0, int v1, int val0, int val01, 
				  int val02, int *vstar0, Vertex *surfvert)
{
  register int i,i1;
  
  /* example of numbering for the case when val0 = 6, val1 = 3 

          v1         _v0
          \        _//|
           \ t1  _/ / |  The interior vertices come last
        t2  \  _/t0/  |  The indices of triangles are shifted by one
     v2______\/___/ t6|  with respect to the vertex star
           v5/\_  \v6 |
        t3  /   \t5\  |
           /  t4  \_\ |
          /         \\|
         v3          v4      */

        
  /* 1) copy the locations of the vertices */

  /* 1.1) star of v0 */
  for (i=0, i1=1; i< val01; i++, i1++)
    memcpy(vert[i], surfvert[vstar0[i1]], sizeof(Vertex));
  memcpy(vert[val01], surfvert[v0], sizeof(Vertex));

  /* 1.2) star of v1 */
  memcpy(vert[val0], surfvert[v1], sizeof(Vertex));

  /* 2) create the faces or connections between vertices */

 face[0][0] = val01; face[0][1] = val0; face[0][2] = 0;
  
  for (i = 1, i1 = 0; i <= val02; i++, i1++)
    {
      face[i][0] =  val01; face[i][1] =  i1; face[i][2] =  i;
    }
  face[val01][0] = val01; 
  face[val01][1] = val02;
  face[val01][2] = val0;
  
  face[val0][0] = val0;
  face[val0][1] = val02;
  face[val0][2] = 0; 

  return;
}


/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfBuildEdgeStar                                                      |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

void _dxfBuildEdgeStar(Vertex *star_vert, Face *star_face, Plane *star_plane, float *star_areas, 
		       float *star_comp, int v0, int v1, u_short val0, u_short val1, 
		       int *star0, int *star1, SimpData *s)

/*...........................................................................*/
/*
   
   represent the edge star with origin at the centroid of the edge link

   copy the vertices and faces of both vertex stars in a single array:
   here is the new numbering of the vertices:

       1      0                    val0+val1-4 
        \    /               \    /   ..
         \  /                 \  /
   .. ____\/____           ____\/____val0+1    
    val0-1/\v0      val0+val1-3/\
         /  \                 /  \
        /    \               /    \val0
       ..    val0 -2             

    The number of triangles in star_val is val0 + val1 -2
    The number of vertices  is             val0 + val1 -4 + 2,
    since v0 and v1 belong to the edge star
    In a vertex star, there is one more vertex (the center) than 
    there are triangles.
*/
/*...........................................................................*/

{
  /* there are two cases: either val1 is exactly 3 or it is > 3 */
 
  int  star_val    = val0 + val1 - 2,
           val013  = star_val    -1,
           val014  = val013      -1,
           val01   = val0        -1,
           val02   = val01       -1,
          *vstar0  = star0       + 2 * val0, 
          *vstar1  = star1       + 2 * val1;

  Vertex translation;

  if (val1 == 3)
    {
      BuildEdgeStarVal1Is3(star_vert, star_face, v0, v1, val0, val01, val02, vstar0, s->vert);
    }
  else
    {
  
      register int i,i1;
 
      /* 0) copy the locations of the vertices */

      /* 0.1) star of v0 */
      
      for (i=0, i1=1; i< val01; i++, i1++)
	memcpy(star_vert[i], s->vert[vstar0[i1]], sizeof(Vertex));

      memcpy(star_vert[val01], s->vert[v0], sizeof(Vertex));

      /* 0.2) star of v1 */

      for (i=val0, i1=2; i<= val014; i++, i1++)
	memcpy(star_vert[i], s->vert[vstar1[i1]], sizeof(Vertex)); 
    
      memcpy(star_vert[val013], s->vert[v1], sizeof(Vertex));

      /* 1) copy the faces */
  
      /* numbering of the faces :
                         val0 +
    \ 1  /               \val1/
   2 \  / 0               \-3/
  ____\/____           ____\/____
      /\                   /\
  .. /  \ val0-1          /  \val0+1
    /    \               /val0\    */

      star_face[0][0] = val01; star_face[0][1] = val013; star_face[0][2] = 0;
  
      for (i = 1, i1 = 0; i <= val02; i++, i1++)
	{
	  star_face[i][0] =  val01; star_face[i][1] =  i1; star_face[i][2] =  i;
	}
      star_face[val01][0] = val01; 
      star_face[val01][1] = val02;
      star_face[val01][2] = val013;
  
      star_face[val0][0] = val013;
      star_face[val0][1] = val02;
      star_face[val0][2] = val0; 
 
  
      for (i = val0+1, i1 = val0; i <= val014; i++, i1++)
	{
	  star_face[i][0] = val013; star_face[i][1] = i1; star_face[i][2] = i;
	}
      
      star_face[val013][0] = val013;
      star_face[val013][1] = val014;
      star_face[val013][2] = 0;
    
    }
  /*  2) copy the normals of the faces and plug them into the 
      equations of the planes  */

  _dxfCopyFaceNormalsAreasCompactness(star_plane, star_areas, star_comp, 
				      star0, star1, val0, val013, s->normal, s->area, s->compactness);

  /* 3) compute the centroid of the edge link, and put it in the
     end of the star_vert array (the vertices v0 and v1 are excluded) */

  CentroidEdgeLink(star_vert[star_val], star_vert, val0, val02, val014);
 
  /* 4) change the coordinate system such that the origin will be the
     centroid of the edge link */
  
  
  _dxfOppositeVector(star_vert[star_val], translation);
  
  _dxfApplyTranslation(star_val, star_vert, translation);
 
  return;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfCopyFaceNormalsAreasCompactness                                    |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

void _dxfCopyFaceNormalsAreasCompactness(Plane *plane, float *s_area, float *s_comp, 
					 int *star0, int *star1, int val0, int val013, 
					 Vertex *t_normal, float *t_area, float *t_comp)
{
  register int i, i1;
  int val01 = val0-1;

  /*  copy the normals of the faces and plug them into the 
     equations of the planes  */
   
  for (i=0, i1=1; i< val01; i++, i1++)
    {
      memcpy(plane[i], t_normal[star0[i1]], sizeof(Vertex));
      s_area[i] = t_area[star0[i1]];
      s_comp[i] = t_comp[star0[i1]];
    }

  memcpy(plane[val01], t_normal[star0[0]], sizeof(Vertex));
  s_area[val01] = t_area[star0[0]];
  s_comp[val01] = t_comp[star0[0]];

  for (i=val0, i1=2; i<= val013; i++, i1++)
    {
      memcpy(plane[i], t_normal[star1[i1]], sizeof(Vertex));
      s_area[i] = t_area[star1[i1]];
      s_comp[i] = t_comp[star1[i1]];
    }

  return;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfApplyTranslation                                                   |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

void _dxfApplyTranslation(int nV, Vertex *vv, Vertex T)
{
  register int i;
  register u_char j;

  for(i=0;i<nV;i++)
  
    for(j=0;j<3;j++)
	
      vv[i][j] = vv[i][j] + T[j];

  return;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfOppositeVector                                                     |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

void _dxfOppositeVector(Vertex u, Vertex v)
{
  register u_char j;

  for(j=0;j<3;j++)

      v[j] = -u[j];

  return;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfMinFloatArray                                                      |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

float _dxfMinFloatArray(float *array, int n)
{
  float min = array[0];
 
  register int i;
 
  for(i=1;i<n;i++) 
    if (array[i] < min)
      min = array[i];
  
  return min;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfStarPlaneEquationsAndVolume                                        |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

double _dxfStarPlaneEquationsAndVolume(Vertex *star_vert,  Face *star_face, Plane  *star_plane, 
				       float *star_areas, int   star_val)
/*...........................................................................*/
/*
   compute the equations of the planes of the star, supposing that the
   normals are known and are stored as the 3 first float values of the
   equations of the planes.
   compute the volume of the star
*/
/*...........................................................................*/
{
  double star_volume_x_6 = 0.;

  register u_short j;

  for(j=0;j<star_val;j++)
    {
      /* determine the equation of the plane when the normal is known.

	 (n,v)  + d  = 0  is the equation, n known vector, d unknown value.
	 (n,v1) + d1 = 0;
	 (n,v2) + d2 = 0; Three data points to determine the optimum d.
	 (n,v3) + d3 = 0;
	 d           = -[ (n,v1) + (n,v2) + (n,v3) ]/3.;
	 */

      register u_char k;

      star_plane[j][3] = 0.;

      for(k=0;k<3;k++)
	
	star_plane[j][3] -= SCALPROD3(star_plane[j], star_vert[star_face[j][k]]);
	
      star_plane[j][3] /= 3.;

      
      /* for volume computation, retrieve the face area from the global array
	 the volume is equal to one third of the triangle area times the
	 height of the origin with respect to the triangle */

      star_volume_x_6  -= star_areas[j] * star_plane[j][3] / 3.;
    }

  return (star_volume_x_6 * 6.);
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfComposeVectors                                                     |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

void _dxfComposeVectors(Vertex u, Vertex v, float lambda, Vertex w)

/*...........................................................................*/
/*
  w = u + lambda v
*/
/*...........................................................................*/

{
  register u_char j;

  for(j=0;j<3;j++)

      w[j] = u[j] + lambda * v[j];

  return;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfPositionVertexLSConstraints2                                       |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int  _dxfPositionVertexLSConstraints2(Plane *star_plane, u_short star_val, Vertex simplified_vertex)

/*...........................................................................*/
/*

     Position the vertex such that 
     1) the volume of the star is preserved
     2) the sum of distances plane-vertex is minimal plus
     3) the distance vertex-origin (minimize the sum of square distances
     vertex - vertices_of_the_link)

     This version uses the Normal Least Squares Equations which have a
     closed form solution

     We have changed the coordinate system such that the normal
     of the plane that guarantees to preserve the volume
     is given either by (0,0,1) or by (0,0,-1)

*/
/*...........................................................................*/
{
  int retval = 1;

  double 
    A1[SIMPLIFY_VALENCE_MAX4],
    A2[SIMPLIFY_VALENCE_MAX4],
    B [SIMPLIFY_VALENCE_MAX4];

  register int i;
  
  /*1) each of the planes provide a LS equation, plus the minimum
    distance to the origin which provides 2 equations */

  double 
    
    /* height of the plane (+/-) z + delta = 0 */

    delta = 
    star_plane[star_val][3] *
    star_plane[star_val][2]; 

  /* 2) write the quadratic constraints associated with the planes */

  for (i = 0; i < star_val; i++)
    {
      A1[i] = star_plane[i][0];
      A2[i] = star_plane[i][1];
      B [i] = delta * star_plane[i][2] - star_plane[i][3] ;
    }

  /* 3) write the quadratic constraints associated with 
     the center of the link which is the current origin */

  A1[star_val    ] = 1; /* x = 0 is the quadratic constraint */
  A1[star_val + 1] = 0; 
  A2[star_val    ] = 0; 
  A2[star_val + 1] = 1; /* y = 0 is the quadratic constraint */
  B[star_val     ] = 
    B[star_val+1] = 0.;

  {
    /* write the normal equations */
    static double a,b,c,d,e,x,y;

    for (a=b=c=d=e=0.,i=0;i<star_val+2 ;i++)
      {
	a += A1[i] * A1[i];
	b += A1[i] * A2[i];
	c += A2[i] * A2[i];
	d += A1[i] * B [i];
	e += A2[i] * B [i];
      }

    /* solve the normal equations */

    _dxfSolveSymmetric2x2Eqn(a,b,c,d,e, &x, &y);
     
    simplified_vertex[0] = (float)  x;
    
    simplified_vertex[1] = (float)  y;
    
    simplified_vertex[2] = (float) -delta;
    
  }

  return retval;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfFPlaneForIdenticalVolume                                           |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

void     _dxfFPlaneForIdenticalVolume(u_short val0, u_short val1, Plane p, Vertex *star_vert, 
				      Face *star_face, float star_volume)
{
  static Vertex  n;
      
  u_short i, val02 = val0 - 2, val013 = val0 + val1 - 3;
  
  dx_bzero((char *)p, sizeof(Vertex));

  for(i=1; i<=val02; i++)
    {
      u_short v1  = star_face[i][1],
              v2  = star_face[i][2];
      VecProd3(star_vert[v1], star_vert[v2], n);
      AddVec(p, n);
    }

  for(i=val0; i<=val013; i++)
    {
      u_short v1  = star_face[i][1],
              v2  = star_face[i][2];
      VecProd3(star_vert[v1], star_vert[v2], n);
      AddVec(p, n);
    }
 
  p[3] = -star_volume;

  /* normalize the equation of the plane */
  {
    double norm = NORM3(p);

    if (norm > 0.0)
      for (i=0; i < 4; i++)
	p[i] /= norm;
  }
   
  return;
}


/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfSimplifiedVertexLocation                                           |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

void _dxfSimplifiedVertexLocation(
			      SimpData *s, 
			      u_short   val0, 
			      u_short   val1, 
			      u_short   star_val,
			      Vertex   *star_vert, 
			      Face     *star_face,
			      Plane    *star_plane, 
			      float    *star_areas,

		    /* workspace needed for the householder transform: */

			      float *work, 
			      u_char method)
/*
  various methods for computing the position of the simplified vertex
*/
{
  if (method == 0)
   
    /* v0 is the simplified vertex */
    memcpy(s->simplified_vertex, star_vert[val0-1], sizeof(Vertex));
  

  else if (method == 1)
    {
      /* take the vertex that minimizes the sum of distances to the
	 vertices of the edge star plus the planes defined by the
	 triangles of the edge star and that preserves the volume of
	 the edge star */

      float  
	star_volume,   /* volume of the edge star (x 6) */
	sign;
	
      int    
	star_val_1      = star_val + 1,
	size_planes     = (star_val_1) * sizeof (Plane);

      float  *star_plane_sav = (float *) (star_plane + star_val_1);

      static Vertex 
	v, z     = {0,0,1.};

      star_volume = _dxfStarPlaneEquationsAndVolume(star_vert, star_face, star_plane, star_areas, star_val);

      /* compute the equation of the plane p such that the star of any
	 vertex on that plane would have the same volume as star_volume */

      _dxfFPlaneForIdenticalVolume(val0, val1, star_plane[star_val], star_vert, star_face,  star_volume);

      /* now change again the coordinate system such that either
	 (0,0,1) or (0,0,-1) corresponds to the normal of P, 

	 change the equations of the planes. THE VERTICES ARE LEFT UNCHANGED */

      /* compute the householder vector */

      sign = (star_plane[star_val][2] > 0.) ? 1. : -1.;

      _dxfComposeVectors(star_plane[star_val], z, sign, v);

      /* Apply the Householder transform to the planes
	 while saving their previous positions */
      
      memcpy(star_plane_sav, star_plane, size_planes);
		     
      _dxfHouseholderPreMultiplication((float *)star_plane, 4, v, 3, star_val_1, work);
	
      /* this version of the positioning solves the Normal Equations
	 rather than the LS equations and has a closed form solution
	 (does not require usage of a numerical library) */

      _dxfPositionVertexLSConstraints2( star_plane, star_val, s->simplified_vertex);

      /* change back the position of the simplified vertex to
	 the original coordinate system */

      _dxfHouseholderPreMultiplication(s->simplified_vertex, 3, v, 3, 1, work);

      /* idem with the planes */
     
      memcpy(star_plane, star_plane_sav, size_planes);
		
    } 

  else if (method == 2)
    {
      /* take the middle of v0 and v1 for the position of the
         simplified vertex */
      u_short val01 = val0-1, val013 = val0+val1-3;

      register u_char j;

      for (j=0;j<3;j++)

	s->simplified_vertex[j] = (star_vert[val01 ][j]  + 
				   star_vert[val013][j] ) /2.;
    }

  return; 
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfCollapseGeometricallyFeasible                                      |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfCollapseGeometricallyFeasible(Vertex *svert, Face *sface, 
				      Plane *splane, float *old_comp, float new_comp_min,
				      Vertex *snor0, Vertex *snor1,   
				      u_short val0, u_short sval, 
				      float min_scalprod, float compactness_ratio)
{
  int feasible = 1;
  
  /* check whether any triangle normal might have switched its orientation */

  if (_dxfNoFlipOverCheck(1,   val0-1, snor0, splane, min_scalprod) == 1 &&
      _dxfNoFlipOverCheck(val0,  sval, snor1, splane, min_scalprod) == 1)
    {
      float old_comp_min = _dxfMinFloatArray(old_comp,sval);
     
      feasible = (new_comp_min  > (compactness_ratio * old_comp_min));
    }
  else   feasible = 0;
  

  return feasible;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfNoFlipOverCheck                                                    |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfNoFlipOverCheck(u_short first_face, u_short last_face, 
			Vertex *nor, Plane *plane, float min_scalprod_nor)
{
  
  register u_short i1=0, i=first_face;

  int orientation_consistent;

  do
    {
      orientation_consistent = 
	(SCALPROD3(nor[i1],plane[i]) > min_scalprod_nor);

      i++;
      i1++;
    }
  while ((orientation_consistent ) && (i < last_face));

  return orientation_consistent;

}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfErrorWithinToleranceV                                              |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfErrorWithinToleranceV(SimpData *simp_data, Vertex *star_vert,
			      Face *star_face, Plane *star_plane, 
			      int *star0, int *star1, int val0, int val1)

{
  /* The purpose of Method V is to construct in parallel a tiling of
     e* (edge undergoing collapsibility tests) denoted til(e*) and a
     tiling of vs* (potential simplified vertex) denoted til(vs*), such
     that there is a one to one mapping between
     til(e*) and til(vs*) 
     */

  /* refresher: numbering strategies for e*:
     
     numbering of the vertices :
     ~~~~~~~~~~~~~~~~~~~~~~~~~
      
       1      0                    val0+val1-4 
        \    /               \    /   ..
         \  /                 \  /
   .. ____\/____           ____\/____val0+1    
  val0-1  /\        val0+val1-3/\
         /  \                 /  \
        /    \               /    \
       ..     val0-2               val0
 
     numbering of the faces : t0,t1,...t_val0+val1-3 triangles of e*
     ~~~~~~~~~~~~~~~~~~~~~~
                             val0 +
        \ 1  /               \val1/
       2 \  / 0               \-3/
      ____\/____           ____\/____
          /\                   /\
      .. /  \ val0-1          /  \val0+1
        /    \               /val0\      

	q1,...,qval0-2,....,qval0,...,q_val0+val1-3 triangles of vs*
	*/

  /* 1) determine the closest point to the vertex vs on the triangles t0
     and t_val0-1 */

  /* 2) determine the closest point to the vertex v0 on the
     triangles q1,...,qval0-2. */

  /* 3) determine the closest point to the vertex v1 on the triangles
     qval0,...,q_val0+val1-3 */

  /* 4) Take in turn each edge (vs,v0),(vs,v1),..., (vs,v_val0-2)
     call the ith edge ei
     determine the closest point on each edge incident to v_val0-1:

     (v_val0-1, v0),..., (v_val0-1, v_val0-2), and (v_val0-1, v_val0+val1-3)
     If the closest point is located in the inside of either edge, 
     determine a crossing point.
     */


  /* 5) Take in turn each edge (vs,v_val0-2),...,(vs,v_val0 + val1-4) and
     (vs,v0),

     call the ith edge ei
     determine the closest point on each edge incident to v_val0+val1-3:

     */


  /* Firstly determine numbers of elements in old and new configurations */
   
  int 
    edge_collapsible = 1,

    /* vertex stars */
    *vstar0 = star0 + 2 * val0, 
    *vstar1 = star1 + 2 * val1,
    i,j, i_min = -1;

  static  Vertex 
    vs[1],             /* local copy of the simplified vertex */
    a_triangle[3],
    wp[2],
    alpha[SIMPLIFY_VALENCE_MAX1]; 
  /* val0 + val1 -4 <= valence max <=> val0, val1 <= valence_max +1 
     because val0 >=3 and val1 >=3 */

  static float
    dp[2];

  float 
    tol_pivot, potential_err_pivot, largest_err_pivot = 0., d_min,
    lambda1, lambda2;

  _dxfErrorTestDataOnSurfaceInit(simp_data);

  /* set the maximum tolerance for the star */

  tol_pivot = simp_data->tol_volume[vstar0[0]];

  for (i=1;i<val0;i++)
    if (simp_data->tol_volume[vstar0[i]] < tol_pivot)
      tol_pivot = simp_data->tol_volume[vstar0[i]];

  for (i=0;i<val1;i++)
    if (simp_data->tol_volume[vstar1[i]] < tol_pivot)
      tol_pivot = simp_data->tol_volume[vstar1[i]];

  
  /* verify that the maximum tolerance is not already exceeded at
     the vertices */

  i = 0;
  while (edge_collapsible && i < val0)
  /* the equality (star_err = tol) is allowed */  
    if (simp_data->err_volume[vstar0[i++]] > tol_pivot)
      edge_collapsible = 0;

  i = 0; 
  while (edge_collapsible && i < val1)
    if (simp_data->err_volume[vstar1[i++]] > tol_pivot)
      edge_collapsible = 0;


  if (edge_collapsible)
    {
 
      Vertex 
	*qi, *t0, *t1;
 
      /* local copy of the simplified vertex, translated by the due amount */
      memcpy(vs[0], simp_data->simplified_vertex, sizeof (Vertex));
      AddVec(vs[0], star_vert[val0+val1-2]);

      /* 1) */  

      /* then fabricate the triangle t0 = v1 v0 vval0-1 */
      memcpy(a_triangle[0], simp_data->vert[vstar0[0]], sizeof (Vertex));
      t0 =  ExtractTriangleFromVertexStar(vstar1,0);

      
      /* find the point closest to vs on t0 and record the projection as
	 well as the distance */
       
      _dxfClosestPointOnTriangle(t0, vs[0], wp[0], alpha[0], &dp[0]);
      
      dp[0] += 

	  alpha[0][0] * simp_data->err_volume[vstar0[0]] +
	  alpha[0][1] * simp_data->err_volume[vstar1[0]] +
	  alpha[0][2] * simp_data->err_volume[vstar1[1]];

      /*  find the point closest to vs on t1 and record the same info 
	  triangle t1 = (v0, v1, v2) */
 
      memcpy(a_triangle[0], simp_data->vert[vstar1[0]], sizeof (Vertex));
      t1 =  ExtractTriangleFromVertexStar(vstar0,0);

      _dxfClosestPointOnTriangle(t1, vs[0], wp[1], alpha[1], &dp[1]);
   
      dp[1] +=
	  
	  alpha[1][0] * simp_data->err_volume[vstar1[0]] +
	  alpha[1][1] * simp_data->err_volume[vstar0[0]] +
	  alpha[1][2] * simp_data->err_volume[vstar0[1]];

      /* the triangle with the smallest potential error wins: */

      if (dp[0] < dp[1])

	{
	  /* use the barycentric coordinates to interpolate the potential data values */
	  
	  _dxfErrorTestDataOnSurfaceInterpolate3(simp_data,   
						 vstar0[0],   vstar1[0], vstar1[1], 
						 alpha[0][0], alpha[0][1], alpha[0][2]);

	  potential_err_pivot = dp[0];

	  /* no need to test the data since we've just use this map point to
	     compute new data values, so they have to match ! */
	  /*
	  if (simp_data->vx_data)

	    edge_collapsible = *//* SimpData simp_data,        number_new_vertices_in_map, number_old_vertices,
	                                         new_v2,     new_v3, alpha_simp_vertex, alpha_n_v2, alpha_n_v3,
			  old_v1,                old_v2,     old_v3, alpha_o_v1,        alpha_o_v2, alpha_o_v3
			  *//*
	    _dxfErrorTestDataOnSurfaceUpdate(simp_data, 1, 3,
                                 0, 0, 1., 0., 0., 
		      vstar0[0], vstar1[0], vstar1[1], alpha[0][0], alpha[0][1], alpha[0][2]);
		      */
	}
      else
	{
	   /* use the barycentric coordinates to interpolate the potential data values */
	  
	  _dxfErrorTestDataOnSurfaceInterpolate3(simp_data,   
						 vstar1[0],   vstar0[0], vstar0[1], 
						 alpha[1][0], alpha[1][1], alpha[1][2]);

	  potential_err_pivot = dp[1];

	  /* no need to test the data since we've just use this map point to
	     compute new data values, so they have to match ! */
	 
	}

      UPDATE_ERR_PIVOT;

      /* 2)  determine the closest point to the vertex v0 on the
	 triangles q1,...,qval0-2. */

      if (edge_collapsible)
	{
	 

	  /* the triangles q1,...,val0-2 are constructed explicitely
	     with the triples of vertices:
	     (vs,v0,v1), vs(v_val03, v_val02) */

	  memcpy(a_triangle[0], vs[0], sizeof(Vertex));

	  i = 1;
	  d_min = 0;

	  while (i< val0-1)
	    {
	      qi = ExtractTriangleFromVertexStar(vstar0,i);

	      _dxfClosestPointOnTriangle(qi, simp_data->vert[vstar1[0]], wp[0], alpha[i], &dp[0]);

	      if (alpha[i][0]>0.)

		/* If v0 does not project on the link, then the projection is useful and can produce
		   a valid constraint */

		{
		  dp[0]    = 
		    ( dp[0]     + simp_data->err_volume[vstar1[0]]
		      - alpha[i][1] * simp_data->err_volume[vstar0[i]]
		      - alpha[i][2] * simp_data->err_volume[vstar0[i+1]])
		    / alpha[i][0] ;
		    
		  if (i_min  == -1 || dp[0] < d_min) {d_min = dp[0]; i_min = i;}
		}
	      else
		{
		  /* otherwise v0 projects onto the link edge of "triangle"
		     the epsilon value at the simplified vertex cannot compensate
		     but it is possible that the error values at the link are high enough
		     to treat that constraint */

		  if (alpha[i][1] * simp_data->err_volume[vstar0[i]]   + 
		      alpha[i][2] * simp_data->err_volume[vstar0[i+1]] >= dp[0])
		    {
		      /* look no further for the closest triangle */
		      i_min = i;
		      i = val0-1; /* next i will be = val0 and the loop will immediately stop */
		      d_min = 0.;
		    }
		}
	      i++;
	    }
	  
	  /* the triangle corresponding to the smallest error wins */

	  if (i_min >=0)
	    {
	      /* if vertex0 does not project on the rim,

		 i.e., 

		 if (alpha[i_min][0] > 0.)


		 a potential new error value is generated at the pivot
		 as follows: */
	  
	     
		  potential_err_pivot = d_min;
		  
		  UPDATE_ERR_PIVOT;

		  if (simp_data->vx_data && edge_collapsible)

		    edge_collapsible = 
		      _dxfErrorTestDataOnSurfaceUpdate(simp_data, 3, 1, 
						       vstar0[i_min], vstar0[i_min+1], 
						       alpha[i_min][0], alpha[i_min][1], alpha[i_min][2],
						       vstar1[0], 0, 0, 1., 0., 0.);

		  /* even though v0 == vstar1[0] is the same vertex, 
		     the data is allowed to be modified at that vertex, 
		     we still have to perform this test in case the simplified
		     vertex would have been allowed to move: we would then verify
		     that the old data is still compatible with the new directions */

	    }

	  else edge_collapsible = 0; /* no projection was useful */
	 
	  /* 4) */

	  if (edge_collapsible)
	    {
	      /* take each edge (vs,v1) .. (vs,v_val0-1)
		 whereby v1, v_val0-1 refer to indices to
		 the star the first vertex vstar0 */
	
	      i = 1;
	      while (edge_collapsible && i<val0)
		{
		 
		  j = 0;
		  while (edge_collapsible && j<val0)
		    {
		
		      d_min = _dxfSegmentIntersection(vs[0], 
						      simp_data->vert[vstar0[i]], 
						      simp_data->vert[vstar1[0]],
						      simp_data->vert[vstar0[j]], 
						      &lambda1, &lambda2, !j, 0, 1);

		      
		      /* If the closest point is located in the
			 inside of either edge, determine a crossing point. */

		      if (lambda2 >= 0. && lambda2 <= 1.) /* inside 2 */
			{


			  d_min += lambda2  * simp_data->err_volume[vstar1[0]] + 
			    (1. -lambda2)  * simp_data->err_volume[vstar0[j]];

			  if (simp_data->vx_data && lambda1 >= 0. && lambda1 <= 1.)

			    edge_collapsible = 

			      _dxfErrorTestDataOnSurfaceUpdate(simp_data, 2, 2,
							       vstar0[i], 0, lambda1, 1. -lambda1, 0.,
						    vstar1[0], vstar0[j], 0, lambda2, 1. -lambda2, 0.);

			  if (lambda1 > 0. && lambda1 <= 1.) /* inside 1 */
		      
			    {
			      potential_err_pivot = 
				(d_min + ( lambda1 - 1.)  * simp_data->err_volume[vstar0[i]])
				/ lambda1;

	    
			      UPDATE_ERR_PIVOT;
			      
			    } 
			  
			  /* LAMBDA1 = 0 IS AGAIN A SPECIAL CASE 
			     if lambda1 = 0
			     we can't treat that constraint and the edge is
			     not collapsible */

			  else if (lambda1 == 0.) 
			    {
			      /* the error at vs cannot accomodate that constraint,
				 but maybe it is already satisfied */

			      if (simp_data->err_volume[vstar0[i]] < d_min)
				edge_collapsible = 0;
			    }
			}

		    

		      j++; if (i==j) j++;
		    }
		  i++;
		}

	      /* 3) determine the closest point to the vertex v1 on
		   the triangles qval0,...,q_val0+val1-3 */


	      if (edge_collapsible)
		{
	
		  i_min = -1;

		  /* the triangles qval0,...,q_val0+val1-3 are
		     constructed explicitely using the vertices of
		     vstar1 and the triples of vertices with the
		     triples of vertices: (vs,v1,v2), (vs,v_val1-2,
		     v_val1-1) */

		  i = 1;

		  while (i<val1-1)
		    {
		      
		      qi = ExtractTriangleFromVertexStar(vstar1,i);

		      _dxfClosestPointOnTriangle(qi, simp_data->vert[vstar0[0]], wp[0], alpha[i], &dp[0]);

		      if  (alpha[i][0]>0.) 

			{  
			  dp[0] =  
			    ( dp[0] + simp_data->err_volume[vstar0[0]]
			    - alpha[i][1] * simp_data->err_volume[vstar1[i]]
			    - alpha[i][2] * simp_data->err_volume[vstar1[i+1]]) 
			    / alpha[i][0];
		   
			  if (i_min == -1 || dp[0] < d_min) {d_min = dp[0]; i_min = i;}
			}
		      else
			{
			  /* otherwise v1 projects onto the link edge of "triangle"
			     the epsilon value at the simplified vertex cannot compensate
			     but it is possible that the error values at the link are high enough
			     to treat that constraint */
			  if (alpha[i][1] * simp_data->err_volume[vstar1[i]]   + 
			      alpha[i][2] * simp_data->err_volume[vstar1[i+1]] >= dp[0])
			    {
			      /* look no further for the closest triangle */
			      i_min = i;
			      i = val1-1; /* and the loop will immediately stop  */
			      d_min = 0.;
			    }

			}
		      i++;
		    }
	  
		  /* the closest triangle wins */
		  
		  if (i_min >=0)
		    
		    {
		      potential_err_pivot = d_min;
		      
		      UPDATE_ERR_PIVOT;
		   

		      if (simp_data->vx_data && edge_collapsible)

			edge_collapsible = 
			  _dxfErrorTestDataOnSurfaceUpdate(simp_data, 3, 1, 
							   vstar1[i_min], vstar1[i_min+1], 
							   alpha[i_min][0], alpha[i_min][1], alpha[i_min][2],
							   vstar0[0], 0, 0, 1., 0., 0.);
		    }
		  else edge_collapsible = 0; /* no projection was useful */
		  	 
		  /* 5) */

		  if (edge_collapsible)
		    {
		      /* Take in turn each edge (vs,v2),(vs,v1),...,
			 (vs,v_val1-2) wherein v2... refer to star 1

			 original vertices of the star:
			 0,1,2,...,val1-1
			 we take
			     2,...,val1-2 */

		      i = 2;
		    
		      while (edge_collapsible && i<val1-1)
			{
		 
			  j = 0;
			  while (edge_collapsible && j<val1)
			    {
		  
			      /* Determine the closest point on each
				 edge incident to "v1" */

			      d_min = _dxfSegmentIntersection(vs[0], 
							      simp_data->vert[vstar1[i]],
							      simp_data->vert[vstar0[0]],
							      simp_data->vert[vstar1[j]],
							      &lambda1, &lambda2, !j, 0, 1);
		      
		      
			      if (lambda2 >= 0. && lambda2 <= 1.) /* inside 2 */
				{
				  d_min += lambda2 * simp_data->err_volume[vstar0[0]] +
				    (1.-lambda2)  * simp_data->err_volume[vstar1[j]];

				  if (simp_data->vx_data && lambda1 >= 0. && lambda1 <= 1.)

				    edge_collapsible = 
				      _dxfErrorTestDataOnSurfaceUpdate(simp_data, 2, 2,
								       vstar1[i], 0, lambda1, 1. -lambda1, 0.,
				                            vstar0[0], vstar1[j], 0, lambda2, 1. -lambda2, 0.);

				  if (lambda1 > 0. && lambda1 <= 1.) /* inside 1 */
		      
				    {
				      potential_err_pivot = 
					(d_min - (1.-lambda1)  * simp_data->err_volume[vstar1[i]])
					/ lambda1;

	    
				      UPDATE_ERR_PIVOT;
	      
				    }
				  else if (lambda1 == 0.)
				    {
				      /* the error at vs cannot accomodate that constraint,
					 but maybe it is already satisfied */

				      if (simp_data->err_volume[vstar1[i]] < d_min)
					edge_collapsible = 0;
				    }
				}
			      j++;if (i==j) j++;
			    }
			  i++;
			}
		    }
		}
	    }
	}
    }
    


  if (edge_collapsible)
    {
      /* update the tolerance and error values */

      /* only the error at the pivot will change */
      
      simp_data->err_volume[vstar1[0]] = largest_err_pivot;

      /* update the tolerance volume */
 	  
      simp_data->tol_volume[vstar1[0]] =  MIN(simp_data->tol_volume[vstar1[0]], 
					      simp_data->tol_volume[vstar0[0]]);

      _dxfErrorTestDataOnSurfaceEnd(simp_data, vstar1[0]);

    }

  return edge_collapsible;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfSimplifiedStarNormals                                              |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

void _dxfSimplifiedStarNormals(Face *sface, Vertex *svert, Vertex *snor, 
			       float *sarea, float *scomp, Vertex simpvert, 
			       u_short first_face, u_short last_face)

/*...........................................................................*/
/*
   compute the face normals in the simplified star as well as the compactness
   and area of the faces.
*/
/*...........................................................................*/
{
  register u_short i,i1;
  register u_char j,k;

  static VertexD triVert[3];
  static double n[4];

  /* introduce the simplified vertex */
  
  CopyType2Type(float, double, simpvert, triVert[0], 3);
      
  for(i=first_face, i1=0; i < last_face; i++, i1++)
    {
      double sum_lengths_sq = 0.0;

      /* get the second and third vertices of the face.. */
      u_short v1  = sface[i][1],
              v2  = sface[i][2];

      /* ..and copy them at the 2nd and 3rd positions of faceVert */

      CopyType2Type(float, double, svert[v1], triVert[1], 3);
      
      CopyType2Type(float, double, svert[v2], triVert[2], 3);

	      
      _dxfTriangleNormalQR2D(triVert, n);

   
      CopyType2Type(double, float, n, snor[i1], 3);

      sarea[i1] = n[3];

       /* compute the squared lenghts of the three edges: */
  
       for(j=0;j<3;j++)
	 {
	   double len_edg_sq = 0, tmp;

	   for(k=0;k<3;k++)
	     {
	       tmp         = triVert[(j+1)%3][k] - triVert[j][k];
	       len_edg_sq += tmp * tmp;
	     }

	   sum_lengths_sq += len_edg_sq;

	 }

       scomp[i1] = FOUR_SQRT_3 * sarea[i1] / sum_lengths_sq;

    }

  return;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfReplaceFaceNormals                                                 |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

void _dxfReplaceFaceNormals(int *star, Vertex *new_nor, float *new_area, float *new_comp, Vertex *nor, 
			    float *area, float *comp, u_short val)
{

/*...........................................................................*/
/*
  
   NOTE the count of the star goes from 0 to val-1
        the corresponding new normals go from 2 to val-1 
	(triangle 0 and 1 are omitted since they will be removed)
*/
/*...........................................................................*/
  
  register u_short i,i1;
 
  for(i=2, i1=0; i<val; i++, i1++)
    {
      memcpy(nor[star[i]], new_nor[i1], sizeof(Vertex));
      area[star[i]] = new_area[i1];
      comp[star[i]] = new_comp[i1];
    }

  return;
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfCollapseEdge                                                       |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfCollapseEdge(SimpData *simp_data, int v0, int v1, int *star0, int *star1, 
		     u_short val0, u_short val1, Vertex *nor0, Vertex *nor1, float *area0, float *area1,
		     float *comp0, float *comp1)

/*...........................................................................*/
/*
    
    Perform the simplification of an interior edge, and update the edge heap, by removing
    all edges of the star and adding back the new ones, keyed with the
    new distances.
 
 father(*0[1])  /\
             /\/  \
              /\   \
             / *0[1]\
            /________\       
         v0 \        /v1
             \ *O[0]/
              \   \/
               \  /\/
                \/  father(*0[1])
  

                           /                             \
       The father of both /                      is       \estar1[val1-1]
                         /estar0[1]                        \
                        /                                   \
                       /________ estar0[0]                   \



       The father  of     /     is  \
                         /           \
                        /estar1[1]    \estar0[val0-1]
                       /               \

    NOTE:
    In the *0 or star0 all triangles are parents of themselves 
    (root of the disjoint-set forest)

    in the surface, all vertices, triangles and edges are their original.
    We have to extract the parent each time.
*/
/*...........................................................................*/
{
  int 
    *estar0 = star0  + val0, *vstar0 = estar0 + val0,
    *estar1 = star1  + val1, *vstar1 = estar1 + val1,
    n_edges_added[2],
    n_edges_reinstated = 0;
 
  /* 1) perform the topological transformations */

  simp_data->triangle_father[star1[0]] = star1[val1-1];
  simp_data->triangle_father[star1[1]] = star1[2];

  simp_data->vertex_father[v1]         = v0; 

  /* the remaining edges must be valid, and the way we choose to pick
     the triangle fathers dictates the way we choose edge father */
    
  simp_data->edge_father[estar0[0]] = 
    simp_data->edge_father[ estar1[val1-1] ] = estar0[1];

  simp_data->edge_father[estar1[1]] = estar0[val0-1];

  /* 2) change the valences of the vertices */

  simp_data->valence[v0]            = val0 + val1 - 4;

  /* v0 is not the only vertex whose valence changes:
     the valence of vstar0[1] and vstar1[1] decreases by one 

       /\vstar0[1]            |
      /  \                    |
     /____\           --->    |
     \    /                   |
      \  /                    |
       \/vstar1[1]            | */

  simp_data->valence[vstar0[1]] --;
  simp_data->valence[vstar1[1]] --;
  

  /* 3 remove the old edge lengths from the heap and put back the new ones */

  /* 3.1 remove all edges in the star of v0 but v0 v1, because it was
   already removed from the heap (current edge under simplification) */

  _dxfRemoveEdgesFromHeap(estar0+1, val0-1, simp_data->edge2index, simp_data);

  /* 3.2 remove all edges but v0-v1 in the star of v1 */

  _dxfRemoveEdgesFromHeap(estar1+1, val1-1, simp_data->edge2index, simp_data);

  /* 4. mark the three edges that are to be discarded from the 
     structure as collapsed with a -4 index */

  /* update the number of edges left to be tested:
     the three collapsed edges wont be tested any more 
     one among them was removed already : estar0[1]
     for the last two, we want to verify that they were not
     already counted as outside the heap */

  if (simp_data->edge2index[estar1[val1-1]] != simp_data->heap_outside_index)
    simp_data->num_edg_remaining_4_testing --;
  if (simp_data->edge2index[estar1[1]] != simp_data->heap_outside_index)
    simp_data->num_edg_remaining_4_testing --;

  simp_data->edge2index[estar1[val1-1]] = simp_data->edge2index[estar1[1]] = -4;

  simp_data->edge2index[estar0[0]] = -7;

  /* 5.3 reinstate in the heap all the edges but 0  in estar0
     that are their own parent, that is, all the edges that are a root
     of the disjoint set forest) */

  n_edges_reinstated  = 
    _dxfReinstateEdgesInHeap(estar0+1, val0-1, simp_data, n_edges_added);
  
  /* 2.3 reinstate in the heap all edges but 0 and 1 and val1-1 in estar1 */

  n_edges_reinstated += 
    _dxfReinstateEdgesInHeap(estar1+2, val1-3, simp_data, n_edges_added+1);
 
  /* update the number of edges left to be tested */

  simp_data->num_edg_remaining_4_testing += n_edges_reinstated;

  /* 3. copy the simplified vertex in the right place */

  memcpy(simp_data->vert[v0], simp_data->simplified_vertex, sizeof(Vertex));

  /* 4. replace the face normals, areas and compactness for the simplified stars with the new ones */

  _dxfReplaceFaceNormals(star0, nor0, area0, comp0, 
			 simp_data->normal, simp_data->area, simp_data->compactness, val0);
  _dxfReplaceFaceNormals(star1, nor1,  area1, comp1, 
			 simp_data->normal, simp_data->area, simp_data->compactness, val1);
 
  return n_edges_added[0] + n_edges_added[1];
}

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfSimplifyManifoldSurface                                            |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

int _dxfSimplifyManifoldSurface(SimpData *simp_data)
{
  int 
    e,
    num_edg_tested                 = 0,                   /* number of edges extracted from the heap */
    num_edg_added_to_heap          = 0;

  /* allocate buffers for storing various information on the vertex and edge stars */

  /* since the valence of the simplified vertex can be as low as  v0 + v1 - 4 (for an interior vertex)
     and is less than SIMPLIFY_VALENCE_MAX, assuming that v0 and v1 >= 1 we
     can be sure that the valence of v0 and v1 that are accepted for processing is <= SIMPLIFY_VALENCE_MAX+3.
     conservatively, we use SIMPLIFY_VALENCE_MAX+4 as a bound */

  int 
    vertex_star_buffer_size = 6 * SIMPLIFY_VALENCE_MAX4 * sizeof (int),
    edge_star_buffer_flsize = 9 * SIMPLIFY_VALENCE_MAX4 * sizeof (float),
    edge_star_buffer_vxsize = 4 * SIMPLIFY_VALENCE_MAX4 * sizeof (Vertex),
    edge_star_buffer_plsize = 4 * SIMPLIFY_VALENCE_MAX4 * sizeof (Plane),
    edge_star_buffer_fasize = 2 * SIMPLIFY_VALENCE_MAX4 * sizeof (Face),

    *vertex_star_buffer = NULL;

  Pointer
    edge_star_buffer_fl = NULL,
    edge_star_buffer_vx = NULL,
    edge_star_buffer_pl = NULL,
    edge_star_buffer_fa = NULL;

  vertex_star_buffer = (int *)    DXAllocateZero(vertex_star_buffer_size);

  if (!vertex_star_buffer) goto error;

  edge_star_buffer_fl = (Pointer) DXAllocateZero(edge_star_buffer_flsize);

  if (!edge_star_buffer_fl) goto error;

  edge_star_buffer_vx = (Pointer) DXAllocateZero(edge_star_buffer_vxsize);

  if (!edge_star_buffer_vx) goto error;

  edge_star_buffer_pl = (Pointer) DXAllocateZero(edge_star_buffer_plsize);

  if (!edge_star_buffer_pl) goto error;

  edge_star_buffer_fa = (Pointer) DXAllocateZero(edge_star_buffer_fasize);

  if (!edge_star_buffer_fa) goto error;

  /* reference such storage arrays in the simp_data structure */

  simp_data->vertex_star_buffer  = vertex_star_buffer;
  simp_data->edge_star_buffer_fl = edge_star_buffer_fl;
  simp_data->edge_star_buffer_vx = edge_star_buffer_vx;


  /* 1) mark the edges touching the surface boundary as not simplifiable: */

  if (!simp_data -> simplify_boundary)
    /*num_edg_adjacent_2_boundary = */
      _dxfMarkEdgesAdjacent2Boundary(simp_data);

  /* 2) fill a Heap with edges, keyed with the weight of the edge */


  while ((num_edg_added_to_heap=_dxfBuildEdgeHeap(simp_data))){  

    simp_data->num_edg_weights += num_edg_added_to_heap;

    do
      {
	    
	e = simp_data->get_lowest_weight_edge(simp_data);

	if (e != -1) 
	  {	
		  
	    /* variable set to 1 if the edge is effectively collapsed */
		  
	    int 
	      collapsible = 0;
	    EdgeS 
	      *edg        = simp_data->edge_array+e;
	    int 
	      *endpoints  = SurfEdgeGetVertices(edg),
	      v0          = FATHER(endpoints[0], simp_data->vertex_father),
	      v1          = FATHER(endpoints[1], simp_data->vertex_father);
	    u_short 
	      val0        = simp_data->valence[v0],
	      val1        = simp_data->valence[v1];

	    num_edg_tested ++;
	    /* mark the edge as removed from the heap: give an index such that HeapDelete would fail + 1*/
		      
	    simp_data->edge2index[e] = simp_data->heap_outside_index;
		      
	    /* decrement the number of edges remaining for testing */
			
	    simp_data->num_edg_remaining_4_testing --;

	    /* treat a boundary vertex of the second kind */

	    if (simp_data->simplify_boundary && simp_data->boundary_vert[v0] && simp_data->boundary_vert[v1])
		    
	      collapsible =
		_dxfCollapsibilityTestsBoundaryEdge2ndKind(simp_data, e, v0, v1,
							   val0, val1, simp_data->preserve_boundary_length);

	    /* treat a boundary vertex of the fisrt kind */

	    else if (simp_data->simplify_boundary && (simp_data->boundary_vert[v0] ||
						      simp_data->boundary_vert[v1]))
		    
	      collapsible =
		_dxfCollapsibilityTestsBoundaryEdge1stKind(simp_data, e, v0, v1, val0, val1);

	    /* treat an interior edge */
	    else
	      {
		/* perform the edge valence test for an interior edge */
		      
		if ( (val0 + val1 - 4 >= 3) && (val0 + val1 - 4 < simp_data->valence_max))
		  {
		    
		    /* BUILD THE 2 VERTEX STARS */
		      
		    int
		      *star0  = (int *) vertex_star_buffer,
		      *estar0 = star0  + val0, *vstar0 = estar0 + val0,
		      *star1  = vstar0 + val0, *estar1 = star1  + val1,
		      *vstar1 = estar1 + val1;
		      
		    dx_bzero(star0, 3 * (val0 + val1) * sizeof (int));
			 
		      /* 3) construct the edge star */
		      
		    _dxfBuildParentEdgeStars(edg,v0,v1,val0,val1,star0,star1,simp_data);
		      
		    /* 3.1) test whether the simplification is topologically feasible */
		      
		    if (_dxfCollapseTopologicallyFeasible(simp_data, vstar0, vstar1, val0, val1))
		      {
			/* 3.2) compute a simplified vertex, and test whether the simplification is
			         geometrically possible, if not, go to END */
			  
			      /* the number of vertices is given by val0 + val1 - 2 and the number of
				 triangles is also val0 + val1 - 2 */
			  
			u_short 
			  sval  = val0 + val1 -2, 
			  sval1 = sval + 1,
			  val02 = val0-2, val12 = val1-2;
			
			/* BUILD THE EDGE STAR */
			  
			Vertex 
			  *svert       = (Vertex *) edge_star_buffer_vx,
			  *snor0       = svert + sval1,
			  *snor1       = snor0 + val02;
			Face 
			  *sface       = (Face  *) edge_star_buffer_fa;
			Plane 
			  *splane      = (Plane *) edge_star_buffer_pl;
			float 
			  *scomp0      = (float *) edge_star_buffer_fl,
			  *sarea0      = scomp0 + val02,
			  *scomp1      = sarea0 + val02,
			  *sarea1      = scomp1 + val12,
			  *scomp       = sarea1 + val12,
			  *sarea       = scomp  + sval,
			  *wspace      = sarea  + sval,
			  min_comp_after;

			dx_bzero((char *)svert,  edge_star_buffer_vxsize);
			dx_bzero((char *)sface,  edge_star_buffer_fasize);
			dx_bzero((char *)splane, edge_star_buffer_plsize);
			dx_bzero((char *)scomp0, edge_star_buffer_flsize);
			      
			/* the new origin is placed at the centroid of the edge link */
			      
			_dxfBuildEdgeStar(svert, sface, splane, sarea, scomp,
					  v0, v1, val0, val1, star0, star1, simp_data);
			  
			_dxfSimplifiedVertexLocation(simp_data, val0, val1, sval, svert, sface,
						     splane, sarea, wspace, simp_data->preserve_volume);
			  
			/* to determine whether the collapse is geometrically
			   possible, the new vertex normals are needed */
			  
			_dxfSimplifiedStarNormals(sface, svert, snor0, sarea0, scomp0, 
						  simp_data->simplified_vertex, 1, val0-1);
			_dxfSimplifiedStarNormals(sface, svert, snor1, sarea1, scomp1, 
						  simp_data->simplified_vertex, val0, sval);
			  
			/* minimum triangle compactness after simplification :*/
			min_comp_after = MIN( _dxfMinFloatArray( scomp0, val02 ), 
					      _dxfMinFloatArray( scomp1, val12 ) );
			  
			if (_dxfCollapseGeometricallyFeasible(svert, sface, splane, scomp, 
							      min_comp_after, snor0, snor1, val0, sval,
							      simp_data->min_scalprod_nor, 
							      simp_data->compactness_ratio))

			  {
			    /* 3.3) if collapsible, test whether the new errors are within the
			       tolerance limits  */
			
			    collapsible = _dxfErrorWithinToleranceV(simp_data, svert, sface, splane,
								    star0, star1, val0, val1);
			      
			    if (collapsible)
				
			      {
				  
				/* 3.4) simplify (collapse the edge ) and replace old edges in
				   the heap with new ones */
				  
				simp_data->num_edg_collapsed += 3;
				  
				/* translate back to the (0,0,0) origin */
				  
				AddVec( simp_data->simplified_vertex, svert[sval]);
				  
				simp_data->num_edg_weights += 
				  _dxfCollapseEdge(simp_data, v0, v1, star0, star1,
						   val0, val1, snor0, snor1, sarea0, sarea1,
						   scomp0, scomp1);
			      }
			    else
			      {
				/* rejected for tolerance off limits the edge can be
				   reinserted afterwards if its length changes*/
				      
				simp_data->edg_rejected_4_tolerance ++;
			      }
			      
			  }
			else
			  {/* geometrically infeasible */

			    simp_data->edg_rejected_4_geometry ++;
			  }
			  
		      }
		    else 
		      { /* topologically infeasible */

			simp_data->edge2index[e] = -3; 
			
			/* this condition will not be modified even if the surface is simplified around 
			   the edge */

			simp_data->edg_rejected_4_topology ++;  
		      } 
		      
			
		  }
		else 
		  {

		    /* If the valence is too high or too low, reject also for topology, However, this condition
		       might change */

		    simp_data->edg_rejected_4_topology ++; 
		  }
	      } /* end treat interior vertex */

	  } /* end if (e != -1) */
	
      }    while(e != -1);

	  
  } /* end while (_dxfBuildEdgeHeap) */
 
  DXFree((Pointer)vertex_star_buffer);

  DXFree(edge_star_buffer_fl);

  DXFree(edge_star_buffer_vx);

  DXFree(edge_star_buffer_pl);

  DXFree(edge_star_buffer_fa);


  return 1;

error:

  if (vertex_star_buffer)    DXFree((Pointer)vertex_star_buffer);

  if (edge_star_buffer_fl)   DXFree(edge_star_buffer_fl);

  if (edge_star_buffer_vx)   DXFree(edge_star_buffer_vx);

  if (edge_star_buffer_pl)   DXFree(edge_star_buffer_pl);

  if (edge_star_buffer_fa)   DXFree(edge_star_buffer_fa);

  return 0;
}

/*****************************************************************************/



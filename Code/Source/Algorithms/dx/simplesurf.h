/***********************************************************************/
/* Open Visualization Data Explorer                                    */
/* (C) Copyright IBM Corp. 1989,1999                                   */
/* ALL RIGHTS RESERVED                                                 */
/* This code licensed under the                                        */
/*    "IBM PUBLIC LICENSE - Open Visualization Data Explorer"          */
/***********************************************************************/

/*********************************************************************/
/*

  Author: Andre GUEZIEC, March/April 1997

  */
/*********************************************************************/


#ifndef _SIMPLESURF_H_
#define _SIMPLESURF_H_

/* SimplifySurface

   is a module that approximates a triangulated surface and resamples data

 */

/* partial bug log 
   ~~~~~~~~~~~~~~~

   this is a partial log of the bugs I discovered when developing this module.
   it is there mainly to help me, and other potential developers working 
   on this module spend less time hunting bugs, especially if we want to
   change the code and redo the same errors.


   3/29/97   DXCull, DXInvalidateConnections, DXInvalidateUnreferencedPositions

   I was trying to run something like:

   obj1 = DXInvalidateUnreferencedPositions(obj)
   obj2 = DXCull(obj1)

   That corresponded to what the User manual suggested
   but I was getting an error at DXGetObjectClass with the result of DXCull

   I decided to code things differently:

   new_in[0] = DXInvalidateConnections(in[0]);
   if (!DXInvalidateUnreferencedPositions(new_in[0])) goto error;
   if (!DXCull(new_in[0]))  goto error;
 
   
   ...........................................................................
   

   3/31/97   DXNewArrayV must be used in place of DXNewArray if a "shape" int * vector
   is passed to the routine rather than an int. Typically if rank > 1 but also if
   rank = 1 and shape is obtained from DXArrayInfo, which expects an int *.
   
   The problem with this bug is that the error message appeared in the end of the module
   something like Out of Memory Out of Memory, and I had no obvious way of tracking it down to the
   source.


   ...........................................................................


   3/31/97      _dxfVertexNormalsfromTriangleNormals(nor, new_nV, new_nT, new_connections,
		  				    face_normals, new_face_areas);
						   
   was executed *after* DXFree((Pointer)new_connections); which is a stupid thing to do!
   I was freeing after 
		      DXAddArrayData(con,  0, new_nT, new_connections);
   and when later on I added the vertex normals from triangle normals routine,
   I forgot that this array was freed. Part of the problem is that I wasnt used to the
   style of allocating and freeing with DXAllocate and DXFree whereby you do a goto error 
   if something bad happens. I was doing things differently:
 
   return code = 0
   buffer = (malloc)

   if (buffer)  { do the work;  free(buffer); set return code to 1}
   else         { write error message }
   
   return(return code)


   ...........................................................................

             

   3/31/97 DXConvertArray() sets an error code if the array is already of the same type that requires
   conversion. Same thing: it is not documented and it took me a while to realize what was going wrong
   since the error was showing up in the end of the module execution, and was not informative:

   Out of Memory, or Null Object, or something similar


   ...........................................................................



   4/2/97 MACRO DEFINITIONS: I had a hard time compiling (2 hours) because of syntax errors
   in the macro definitions that are not picked up by the compiler. They are detected much
   later and point so seemingly incomprehensible errors in the c-code.


   ...........................................................................


   4/2/97 SurfEdgeLabel: when building the "e_table" hash table of edges, I was forgetting 
   to initialize the label of the edge to 1 

   SurfEdgeLabel(edge) = 1;

   Accordingly, edges that were boundary edges had a label of 28800 or something like that
   (an arbitrary value) and were misclassified as interior edges.

   then, the routine that was looping around edges was crashing because it expected boundary
   edges to be interior edges.

   This bug was easy to locate in dbx: I started from the location where the program crashed,
   I printed the edge array and noticed right away that
   the labels were funny. I was leaded quickly to the routine that generated the edges



   ...........................................................................


   4/2/97 DXAllocate the valence instead of DXAllocateZero: stupid error. I had only two
   DXAllocate instead of DXAllocateZero in the entire program. 

   I noticed it in dbx because the module was simplifying surface so little.
   by stopping at the right line, I notices way too many topologically infeasible collapses.
   That made me want to check the valence

   simp_data->boundary_vert was also DXAllocated and that was probably another bug, fixed
   by hitting two birds with a stone.




   ...........................................................................



   4/3/97 rank and shape of a flag input. By default, the module builder was using the test

   if (type != TYPE_INT || category != CATEGORY_REAL || rank != 1 || shape != 2) 
      DXSetError(ERROR_DATA_INVALID, "input \"preserve_volume\"");
   
   the right test to use is 

   if (type != TYPE_INT || category != CATEGORY_REAL || !((rank == 0) || ((rank == 1)&&(shape == 1))))



   ...........................................................................



   4/3/97 DXGetComponentValue(field, "positional_error")

   the "positional_error" component was ignored by the simplification because
   of a faulty test for the "dep" attribute
   if was using 

   if (strcmp("positions", DXGetString((String)attr)))

   instead of 

   if (!strcmp("positions", DXGetString((String)attr)))

   which is the "c" translation of : are "positions" and DXGetString((String)attr) the same

   the (!) will presumably get me many more times in the future. This is the advantage of
   documenting the bugs. Hopefully I'll remember this one.



   ...........................................................................


   4/3/97 bug in the routine _dxfFloatDataBoundingBoxDiagonal. I noticed it when running
   an example that fed several isosurfaces to the module. The bounding box value was inconsistent.

   the code was as follows:

   for(j=0;j<data_dim;j++)
    {
      bbox[j] = bbox[j+data_dim] = data[j];
    }
    

  cur_data = data + data_dim;

  for (i=0;i<n_data;i++,cur_data+=data_dim)

    for(j=0;j<data_dim;j++)
      {
	if (cur_data[j] < bbox[j]) bbox[j] = cur_data[j];

	else if (cur_data[j] > bbox[j+data_dim]) bbox[j+data_dim] = cur_data[j];
      }

      The problem was that the first set of data values was used to initialize the bounding box
      and *then* after incrementing cur_data, I was examining n_data values.
      Hence I was peeking outside of the data array and picking any ridiculous value there or
      risking a core dump

      one solution is to do the loop as

      for (i=1;i<n_data;i++,cur_data+=data_dim)


   ...........................................................................


   4/3/97 GET_POSITIONAL_ERROR the existing positional error was not used. I noticed it because
   the module was apparently simplifying too much on a second simplification.

   what happened is that the portion of the code 

   starting with

   array = (Array)DXGetComponentValue(field, "positional_error");

   was just before a goto label GET_OTHER_INPUTS, and the code jumped to that goto in case
   no data "dep" positions was available, which concerns most of the cases. 

   so I changed the name of the goto label to  GET_POSITIONAL_ERROR and placed it before the
   relevant code.

   ...........................................................................


   4/3/97 this was an important omission rather than a bug.

   when running the debugged module on non-manifold input, I realize that if
   some faces in the input were degenerate, (if the three vertices of the face
   weren't all different), the conversion to manifold would not eliminate such faces.

   Curiously, this elimination of degenerate faces was done in the first version
   of the conversion to manifold, but when re-programming it I decided that it
   was not necessary. I think that it *is* necessary, because it perturbs 
   the joining of the corners. I am supposed to join corner c1 and c2 and c3 and c4
   and all the corners must be different. If a face is degenerate, then for
   instance c1 and c3 are the same and I don't know what will happen.

   This is a nasty omission, because a lot of code was depending on the number of
   faces not being changed in the conversion to manifold. Now there is a "face_lut"
   and a lot of code needed to be updated. This took approximately 3 hours.
   
   However, after incorporating this capability, problems still subsist,
   leading us to the next bug:
   

   ...........................................................................


   4/3/97 - 4/4/97 dealing with erroneous input
  
   in order to convert non manifold input surfaces to manifold surfaces,
   it is necessary to identify the manifold edges and to join triangles
   that share a manifold edge. If the edge is not manifold, no triangles can
   be joined there.

   ...........................................................................

   4/4/97 with new method for fixing non-manifold input, the edge count
   was not correct. I was counting too many edges and hence there were some
   additional spurious edges.

   In fact, in the routine _dxfLabelEdges() it is absolutely *necessary*
   to compute the number of edges again, as some triangles might be
   implicitely joined and there may be fewer edges than we previously thought.

   This bug was very difficult to debug because the -g version was running fine
   in the debugger, but I was convinced that it was a problem of memory being overwritten.
   In practice, no memory was overwritten but there were fewer valid elements in the
   edges array than the counter predicted leading in any catastrophic behavior when
   such edges were processed.

   ...........................................................................

   4/4/97 SIGH! as if it wasn't enough for today: I fixed the counting problem for
   the number of edges nE, and thus passed the pointer &nE to a bunch of routines

   I forgot that _dxfInitHashTable(e_table, _dxfFlog2(nE) + 1,...
   requires an estimate of the number of edges ! nE.
   It was receiving a NE of 0!!!, resulting in all the hash entries being added to the
   same list. 
   
   I noticed it because the program was soooo slooooow! 

   ...........................................................................

   4/6/97 definitions of uint, u_int, ushort, u_short ...

	the problem is that if a compiler knows about them on a particular
	platform, another compiler on a different platform might not know!

	below is a series of directives ifdef and include that seem to
	work

   ...........................................................................

   4/7/97 DXArrayConvert again
   
   I was getting an error ERROR: SimplifySurface: Invalid data: data ranks dont match
   that after one hour of debugging, I attributed to DXArrayConvert again
   so now I call instead of DXArrayConvert(array, TYPE_FLOAT, CATEGORY_REAL, 0)
   
   data_array = DXArrayConvert(array, TYPE_FLOAT, CATEGORY_REAL, rank);

   ...........................................................................

    4/7/97 

	line 6226 of the _simplifysurface.c file
    
    the line was

    s_comp[i] = s_comp[star1[i1]];

    instead of:

    s_comp[i] = t_comp[star1[i1]];

   resulting in reading the memory in some strange place
   it took me 3 hours to hunt this bug which was causing a segmentation
   fault on a dec machine (but not on other machines apparently)
   the debugger dbx was very useful for that, as I could notice that
   the array s_comp was not filled properly. It was surprisingly hard
   to finally find this typo, though, who easily survived compiling.

   ...........................................................................

    4/7/97 

     *new_edges =  (EdgeS *) DXReAllocate(*new_edges, nE * sizeof (EdgeS));
     line 760 of _simplifysurface.c

    the problem is that I was reallocating *after* filling up the h-table
    of edges. The edges were copied to another place and the edge entries
    in the h-table were pointing to nothing=or to a completely different
    piece of memory.

    This bug only appeared in the purified version of the dxexec because
    the reallocate was a realloc with a smaller size, which results in
    leaving the memory where it is with most implementations, except
    the implementation of Paula for the memorystubs.o which mallocs,
    copies and frees (places the information in a completely different memory
    location.

    The most annoying part of this bug is that I am sure to have encountered
    it before. I should have diagnosed it much sooner.


   ...........................................................................

   4/8/97

   SQRT: DOMAIN ERROR observed on hp and on sun. What I did was to
   change most of the sqrt calls, except those embedded into macros,
   and make sure that do not call sqrt on 0 or negative input.

   ...........................................................................
   
   4/9/97

   The "positional_error" component was probably resampled like
   the other components. Then there should have been two components
   called  "positional_error" dep on "positions": the newly created component
   + the resampled component. Apparently there was only one??

   What I did was to add a test to make sure that "positional_error"
   would not be resampled.

   ...........................................................................

   4/10/97

   function 
   Nodetype *_dxfnewnode()

   I wasnt checking that what was returned was different from NULL
   before doing

   Nodetype *retval = _dxfgetnode(table);

   retval->code     = code;
   retval->rec      = rec;

   now it reads 
   if (retval)
    {
      retval->code     = code;
      retval->rec      = rec;
    }
    return(retval);

    I hope that it will get rid of the "Memory Corrupt" signals 
  
   ...........................................................................

   4/11/97 improvement rather than bug fix

   I added a "data" switch in order to turn on or off the use of
   vertex related data to decide whether to collapse edges
   ...........................................................................

   4/12/97-4/13/97 improvement rather than bug fix

   I added a "stats" switch in order to provide information on
   how many vertices and triangles the 
   module starts with and how many vertices and triangles the module ends up
   with


   ...........................................................................


   4/15/97 Lloydt reported that after SimplifySurface, various attributes
   of dx fields were lost

   --> call DXCopyAttributes() after the components are created

   when printing the objects resulting from SimplifySurface, I noticed
   that sometimes the Bounding Box and Neighbors were created and sometimes
   they werent.

   --> call DXEndField() after the field is created (after the worker 
   routine has completed


   
   ...........................................................................


   4/30/97 Lloydt noted that SimplifySurface passes "quads" through
   without processing them

   --> implemented a _dxfQuads2Triangles internal conversion routine
   from quads to triangles. The same connection dependent
   data mapping can be used as before using the "face_lut" array.

   

   ...........................................................................

   4/30/97 treat the following classes in the doLeaf() handler

    case CLASS_PRIVATE:
    case CLASS_LIGHT:
    case CLASS_CAMERA:

    (before, I would have set an error "ERROR_BAD_CLASS", "encountered in object traversal");

   ...........................................................................
 
   4/30/97 series are handled better

   --> set is_series = (DXGetGroupClass((Group)in[0]) == CLASS_SERIES);


   --> use

   new_in[0] = DXGetSeriesMember((Series)in[0], i, &position);
 
   DXSetSeriesMember((Series)new_group, i, (double) position, new_out[0]);

   ...........................................................................

   note: most of the 4/29/96 changes were done during the vacation period
   4/15--4/28. They were checked on 4/30/96

   ...........................................................................

   5/1/97 

   I discovered that the "memory corrupt" signals are set by the DXFree routine
   when it frees a portion of memory that was already freed.

   in the routine _dxfBuildOrientedManifold I had the following piece of code:

   if (*vlut) 
    {DXFree((Pointer)(*vlut));}
 
    *vlut = new_vlut; 

    be careful, here, in case there is an error, new_vlut will be freed


   I had to add the following after the error: flag

   if (new_vlut) 
    {
      DXFree((Pointer)new_vlut); 
      new_vlut was computed and vlut freed and replaced with new_vlut,
      so we must set *vlut = NULL 
      *vlut = NULL;}


      also I redid the routines _dxfCreateSimplificationDataStructure
      _dxfFreeSimplificationDataStructure 
      _dxfPartialFreeSimplificationDataStructure to make sure that pointers were set
     to NULL after being freed.

     _dxfPartialFreeSimplificationDataStructure was created to make space and to avoid running
     out of memory in _dxfBuildSimplifiedSurface

     also. In
     _dxfBuildSimplifiedSurface I made sure pointers were set to NULL after being
     freed, to avoid freeing them twice in the error: part of the code

     the building of new surface vertices and new surface triangles was separated
     for clarity and for freeing un-used memory before allocating new memory

   ...........................................................................

      5/7/97 bug of Floating Point exception noticed on DEC

      in the function  _dxfSimplifiedBoundaryVertexLocation()

      around lines 3700-3800 of _simplifysurface.c

	   A special case is encountered if c == 0,
	   then we also have x = 0 and the formula is simply ye=a
	   this case is handled with the general case without particular
	   attention, except that we can't divide by c and
	   also we should not execute the  portion of code relative to 
	   t == 0 or t == 1
 
   ...........................................................................
   */ 


#include "dx.h"
#include <math.h>

#if defined(HAVE_VALUES_H)
#include <values.h>
#endif

#if defined(HAVE_SYS_BSD_TYPES_H)
#include  <sys/bsd_types.h>
#endif

#if defined(HAVE_SYS_TYPES_H)
#include <sys/types.h>
#endif

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

/* SC: defs for Visual C++ */
#ifdef _MSC_VER
typedef unsigned char u_char;
typedef unsigned short u_short;
typedef unsigned int u_int;
#endif

/* some vector operations: */

#undef  SQUARED_DISTANCE3
#define SQUARED_DISTANCE3(u,v) (((u)[0]-(v)[0])*((u)[0]-(v)[0]) + ((u)[1]-(v)[1])*((u)[1]-(v)[1]) +\
                                ((u)[2]-(v)[2])*((u)[2]-(v)[2]))


     /* At some point, I might wast to make sure that sqrt is called on non-zero input.
	I will have to replace the macro with a function call */
	
#undef  DIST3
#define DIST3(u,v)             ( sqrt ( SQUARED_DISTANCE3((u),(v))))

#undef  SCALPROD3
#define SCALPROD3(u,v)         ( (u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2] )

#undef  NORM3
#define NORM3(u)               ( sqrt ( SCALPROD3((u),(u)) ) )

#undef  MakeVec
#define MakeVec(u,v,w) {(w)[0] = (v)[0] - (u)[0]; (w)[1] = (v)[1] - (u)[1]; (w)[2] = (v)[2] - (u)[2];}

#undef  AddVec
#define AddVec(u,v)    {/* u = u + v */ (u)[0] += (v)[0]; (u)[1] += (v)[1]; (u)[2] += (v)[2];}

#undef  VecCopy
#define VecCopy(u,v)   {(v)[0] = (u)[0]; (v)[1] = (u)[1]; (v)[2] = (u)[2];}

#undef  VecProd3
#define VecProd3(u,v,w) { (w)[0] = (u)[1]*(v)[2]-(u)[2]*(v)[1]; (w)[1] = (u)[2]*(v)[0]-(u)[0]*(v)[2]; \
                          (w)[2] = (u)[0]*(v)[1]-(u)[1]*(v)[0];}

#undef  EdgeBarycentricCoordinates
#define EdgeBarycentricCoordinates(u,v,a,b) { a = SCALPROD3(u,v); b = SCALPROD3(u,u); }


/* prototype of the driver routine for simplification */

int _dxfSimplifySurfaceDriver(int nV, float *v, int nT, int *t,
			     int dim_data,
			     float *vertex_data,
			     float tolerance, 
			     float data_tolerance,
			     /* to be able to resample data defined on a non-manifold
				surface, we need look-up tables between vertices before
				and after manifold conversion and between triangles before 
				and after manifold conversion */
			     int *old_nV_after_manifold_conversion,
			     int *old_nT_after_manifold_conversion,
			     int *new_nV,             /* new number of vertices after simplification*/
			     float **new_v,           /* new vertex array (positions) after 
							 simplification */
			     float **new_vertex_data, /* new data values after simplification */
			     int *new_nT,             /* new number of triangles after simplification*/
			     int   **new_t,           /* new triangle array (connections) after
							 simplification */
			     int   **vertex_parents,  /* relationship between simplified vertices
							 and manifold converted vertices */
			     int   **face_parents,    /* relationship between simplified triangles
							 and manifold converted triangles */
			     int   **vertex_lut,      /* for conversion from non-manifold */
			     int   **face_lut,        /* to manifold surfaces */

			     float *old_positional_error, 

	     /* in case we resimplify a surface that was already simplified,
		we use the old positional error as a starting point. The old positional
		error may also have been provided by uncertainty measurements */
								     
			     float **new_positional_error, 
			     float **face_normals,   /* normals of the simplified faces */
			     float **old_face_areas, /* face areas after conversion to manifold */
			     float **new_face_areas, /* face areas after simplification */
			     int preserve_volume,    /* flag = 1 is the volume should be preserved */
			     int simplify_boundary,  /* flag = 1 if the boundary should be
							simplified */
			     /* flag = 1 if the boundary length should be preserved, in case
				the boundary is simplified (simplify_boundary = 1) */
			     int preserve_boundary_length);


				    
float _dxfFloatDataBoundingBoxDiagonal(float *data, int n_data, int data_dim);
/* compute the bounding box of any multidimensional array of data */
    
#ifdef DX_HAVE_ARRAY

int   _dxfVertexNormalsfromTriangleNormals(Array nor, int nV, int nT, int *triangles,
					  float *face_normals, float *face_areas);
/* convert the triangle normals to vertex normals */

/* resample components dep positions and dep connections after simplification 
   "resample" means: 1) create a new component with the same name for the simplified surface
                     2) compute new value at each new vertex or at each new triangle
		     by averaging the values at the corresponding old vertices and the corresponding
		     old triangles:
		     
		     the correspondences between old elements and new elements are obtained by
		     querying the following arrays:  *vertex_parents, *face_parents, *vertex_lut, 
		     *face_lut

		     the array "vertex_parents" gives for each old vertex after conversion to manifold
		     the new vertex to which it maps
		     the array "vertex_lut" gives for each vertex after conversion the corresponding
		     vertex before conversion.
 
		     the array "face_parents" gives for each old triangle after conversion to manifold
		     the new triangle to which it maps
		     the array "face_lut" gives for each triangle after conversion the corresponding
		     triangle before conversion.

		     When averaging vertex attributes, each vertex has weight 1
		     When averaging triangle attributes, each triangle has a weight equal to its area

		     the following components are not resampled
		     POSITIONS, CONNECTIONS, INVALID_POSITIONS INVALID_CONNECTIONS
		     NORMALS, NEIGHBORS, POSITIONAL_ERROR
		     and DATA if it is used by SimplifySurface.

   */

int   _dxfResampleComponentValuesAfterSimplification(Field original_surface, 
						    Field *simplified_surface,
						    int old_nV_after_conversion,
						    int old_nT_after_conversion,
						    int new_nV,
						    int new_nT,
						    int *vertex_parents, int *face_parents, 
						    int *vertex_lut, int *face_lut, float *face_areas);

/* routine specialized for resampling position dependent attributes */

Array _dxfResampleComponentDepPosAfterSimplification(Array array, 
						    int old_nV_after_conversion,
						    int new_nV, int *vertex_parents, 
						    int *vertex_lut);

/* routine specialized for resampling connection dependent attributes */

Array _dxfResampleComponentDepConAfterSimplification(Array array, 
						    int old_nT_after_conversion,
						    int new_nT,
						    int *face_parents, int *face_lut,
						    float *face_areas);

#endif  // DX_HAVE_ARRAY

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfComputeVertexValence                                               |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

/* using the set of triangles, compute the valence of each vertex */

#define _dxfComputeVertexValence(nT, t, valence)                               \
{                                                                              \
  register int i, i1;                                                          \
  register u_char k;                                                           \
                                                                               \
  /* for each triangle, increment the valence of its three vertices */         \
                                                                               \
  for (i=i1=0; i<(nT); i++, i1+=3) for (k=0; k<3; k++) (valence)[(t)[i1+k]]++; \
										 }


/* 
   following the suggestion of Greg Abram

   I am using macros to implement the resampling functions,
   so that they can be used on all the different scalar data types
   without having to copy the code every time.

   A few other macros are defined later on to reduce the number of 
   function calls on some computation-intensive routines
 */ 



/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfRESAMPLE_DEP_POS_AFTER_SIMPLIFICATION                              |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

#define _dxfRESAMPLE_DEP_POS_AFTER_SIMPLIFICATION(type)            \
{                                                                  \
  register int j;                                                  \
                                                                   \
  type                                                             \
    *old_data = (type *) DXGetArrayData(array),                    \
    *old_data_current,                                             \
    *new_data = (type *) new_array_data;                           \
                                                                   \
  float                                                            \
    *vertex_data_averages_current;                                 \
                                                                   \
  int                                                              \
    total_new_num_elements = new_nV * tensor_size;                 \
                                                                   \
  /*for each data element (including each tensor element)          \
    we compute the average using the vertex lut and vertex parent  \
                                                                   \
    we need a vertex_weights array allocated to zero beforehand    \
                                                                   \
    we need an array of floats to store temporarily the averages   \
    at each vertex                                                 \
    */                                                             \
                                                                   \
                                                                                       \
 for(i=0; i<old_nV_after_conversion; i++) /* this is the old number of vertices        \
					     after the conversion to manifolds */      \
    {                                                                                  \
      int                                                                              \
                                                                                       \
	old_vertex = vertex_lut[i],                                                    \
	new_vertex = vertex_parents[i];                                                \
                                                                                       \
      /* we add one to the vertex weights of the parents */                            \
                                                                                       \
      vertex_weights[new_vertex] += 1;                                                 \
                                                                                       \
      /* we add the tensor_size values to the average */                               \
                                                                                       \
      /* get the right offset in both old data and new data arrays */                  \
                                                                                       \
      old_data_current = old_data + old_vertex * tensor_size;                          \
                                                                                       \
      vertex_data_averages_current = vertex_data_averages + new_vertex * tensor_size;  \
                                                                                       \
      for (j=0;j<tensor_size;j++)                                                      \
                                                                                       \
	vertex_data_averages_current[j] += (float) old_data_current[j];                \
                                                                                       \
    }                                                                                  \
                                                                                       \
  /* then loop on all the averages, and divide them by the vertex_weights */           \
                                                                                       \
  vertex_data_averages_current = vertex_data_averages;                                 \
                                                                                       \
  for(i=0; i< new_nV; i++,vertex_data_averages_current+=tensor_size)                   \
                                                                                       \
    for (j=0;j<tensor_size;j++)                                                        \
                                                                                       \
      vertex_data_averages_current[j] /= vertex_weights[i];                            \
                                                                                       \
  /* finally convert the vertex_data_averages to the requested scalar type */          \
                                                                                       \
                                                                                       \
  for(i=0; i< total_new_num_elements; i++)                                             \
                                                                                       \
    new_data[i] = (type) vertex_data_averages[i];                                      \
											 }


/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfRESAMPLE_DEP_CON_AFTER_SIMPLIFICATION                              |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

#define _dxfRESAMPLE_DEP_CON_AFTER_SIMPLIFICATION(type)            \
{                                                                  \
  register int j;                                                  \
                                                                   \
  type                                                             \
    *old_data = (type *) DXGetArrayData(array),                    \
    *old_data_current = old_data,                                  \
    *new_data = (type *) new_array_data;                           \
                                                                   \
  float                                                            \
    *face_data_averages_current;                                   \
                                                                   \
  int                                                              \
    total_new_num_elements = new_nT * tensor_size;                 \
                                                                   \
  for(i=0; i< old_nT_after_conversion; i++, old_data_current += tensor_size)           \
      /* this is the old number of faces after conversion and before simplification */ \
    {                                                                                  \
      int                                                                              \
	old_face = face_lut[i],                                                        \
	new_face = face_parents[i];                                                    \
                                                                                       \
      /* we add the area of the face to the weight of the parent triangle */           \
                                                                                       \
      face_data_weights[new_face] += old_face_areas[old_face];                         \
                                                                                       \
      /* we next update the averages for all associated tensor_size values */          \
                                                                                       \
      /* get the right offset in  both old data and new data arrays */                 \
                                                                                       \
      old_data_current = old_data + old_face * tensor_size;                            \
                                                                                       \
      face_data_averages_current = face_data_averages + new_face * tensor_size;        \
                                                                                       \
      for (j=0;j<tensor_size;j++)                                                      \
                                                                                       \
	face_data_averages_current[j] += old_face_areas[old_face] * (float)            \
							     old_data_current[j];      \
                                                                                       \
    }                                                                                  \
                                                                                       \
  /* then loop on all the averages, and divide them by the face_weights */             \
                                                                                       \
  face_data_averages_current = face_data_averages;                                     \
                                                                                       \
  for(i=0; i< new_nT; i++,face_data_averages_current+=tensor_size)                     \
                                                                                       \
    for (j=0;j<tensor_size;j++)                                                        \
                                                                                       \
      face_data_averages_current[j] /= face_data_weights[i];                           \
                                                                                       \
  /* finally convert the vertex_data_averages to the requested scalar type */          \
                                                                                       \
                                                                                       \
  for(i=0; i< total_new_num_elements; i++)                                             \
                                                                                       \
    new_data[i] = (type) face_data_averages[i];                                        \
											 }


/* copy an array from one type to another type */

/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   CopyType2Type                                                          |
  |                                                                          |
  +--------------------------------------------------------------------------+*/

#define CopyType2Type(type1, type2, a, b, n)                       \
{                                                                  \
  register int i;                                                  \
                                                                   \
  type1 *a1 = (type1 *) (a);                                       \
                                                                   \
  type2 *b1 = (type2 *) (b);                                       \
                                                                   \
  for (i=0; i<(n) ; i++) b1[i] = (type2) a1[i];}


/* convert type "quad" connections to type "triangles" */

int _dxfQuads2Triangles(int *nT, int **connections, int **face_lut);


/* data structures for the surface edges and associated routines */


typedef struct _EdgeS_ {

  short            label;               /* used for a flag or for
					   classifying edges among boundary, regular and singular edges */
  int		    v[2];		/* end points */
  int		    t[2];		/* incident triangles */
} EdgeS;

/*
  surface edges have an orientation:
  v[0] >= v[1] and t[0] is before t[1] in clockwise order
  when we rotate around v0.

                      /|\v1
                     / | \
                    /  |  \              
                   / t0|t1 \             
                   \   |   /
                    \  |  /
		     \ | /
		      \|/v0    v0 >= v1
*/

/* routines for indexing edges in a hash table */

int         _dxfEdgeSKey(EdgeS *);
int         _dxfEdgeSFilter(EdgeS *, EdgeS *);


#define SurfEdgeLabel(edg)           (edg)->label
#define SurfEdgeGetVertices(edg)     (edg)->v
#define SurfEdgeGetTriangles(edg)    (edg)->t


/* my hash tables:

   they do not have the restriction that no more than
   16 elements should have the same code, that is documented in the DX hash tables: 
   
   if someone gives me a surface with 1,000,000 triangles, it should
   roughly have 1,500,000 edges and I don't know of a simple coding
   function that would never give 16 times the same code to edges (it
   may be possible to design such a function but I don't know how)

   */
   
/*...........................................................................*/

/* Following are includes needed for hash-tables using linear probing.
   In case of a code conflict, data items are chained */

/*...........................................................................*/

typedef char  *RECTYPE ;
typedef char  *KEYTYPE ;

typedef struct _nodetype_{
  int   code;
  RECTYPE rec;
  struct _nodetype_ *next;
} Nodetype;

typedef struct _Block_{
  RECTYPE rec;
  struct _Block_ *next;
} Block;

typedef struct _htable_{
  Nodetype **bucket;
  Block     *block;
  Nodetype  *freenode;

  int nodes_per_block ;
  int size, mask;
  int        search_code;
  Nodetype  *search_node;
  KEYTYPE     search_key;

  int     (*fcode)(KEYTYPE);
  int     (*filter)(RECTYPE, RECTYPE);
  KEYTYPE (*compute_key)(RECTYPE);

} Htable;

/*...........................................................................*/

/* prototypes for H-table procedures 
*/

/*...........................................................................*/

/* add a record to the hash-table */
int         _dxfAddRecord(Htable *, RECTYPE, KEYTYPE);
KEYTYPE     _dxfIdentifyRecord2Key(RECTYPE rec);


/* initialize the hash-table */                                                            

int         _dxfInitHashTable(Htable *, int, int (*fcode)(KEYTYPE), 
			     int (*filter)(RECTYPE,RECTYPE), 
			     KEYTYPE (*compute_key)(RECTYPE));

/* determine whether a particular record is in the hash-table and if yes, find that record 
   meaning, return a pointer to that record*/

char        *_dxfFindRecord(Htable *, KEYTYPE);
void        _dxfHashTableReset(Htable *);

int        _dxfNewBlock(Htable *);

/* once _dxfFindRecord() has been called, retrieve the next record with the same code and
   passing through the same filter() from the hash table */

char        *_dxfNextRecord(Htable *);
Nodetype    *_dxfgetnode(Htable *);
Nodetype    *_dxfnewnode(Htable *,RECTYPE, int);

int          _dxfFlog2(int n); /* log base two using integer arithmetic */

/* reset the hash table and free the buckets
   after that, either (1) inithashtable must be called again to refill the hash table with new records or
   completely different data
   
   or (2) free(table) to finish it off
   unless it was defined as Htable table[1] in which case nothing more has to be done */

#define       _dxfFreeHashTable(table) {_dxfHashTableReset(table);if ((table)->bucket) \
  {DXFree((Pointer)(table)->bucket);(table)->bucket = NULL;}}


/* find a particular surface edge in a Hash table indexing edges */

EdgeS *_dxfFindEdge(Htable *e, int vA, int vB);



/* convert a surface mesh to an oriented manifold surface mesh
   manifold means that each vertex has a single sheet of triangles touching it
   and each edge has at most two incident triangles */

int    _dxfToOrientedManifold(int nV,   float *v,   int nT,   int *t, 
				  int *nVm, float **vm, int *nTm, int **tm, int **vlut, int **flut,
				  int *n_edges, Htable *e_table, EdgeS **e);

/* this is the first operation of conversion to a manifold:
   triangles that are degenerate are eliminated
   a triangle is degenerate if either (1) one or more vertices are less then 0 or more than nV
   or (2) two or three vertices in the triangle are the same
   */

int    _dxfEliminateDegenerateTriangles(int nV, int nT, int *t, int *new_nT, int **new_t, int **flut);

/* second operation of conversion to manifold:
   the vertices that are not referenced after eliminating the degenerate triangles are
   eliminated from the list of vertices */

int    _dxfEliminateStandaloneVertices(int nV, float *v, int *nVm, float **vm, int nT, int *t, 
					      int **vlut);

/* then build a list of edges that are manifold edges, and join the corners of the triangles
   that share such an edge */

int   _dxfJoinTriangleCorners(int nT, int *t, Htable *e_table, int *n_edges, EdgeS **edges, 
				     int *fathers);

int    _dxfBuildOrientedManifold(int *nV, float **v, int nT, int *t, int **vlut, 
					int *nE, EdgeS **edges, Htable *e_table, int *fathers);

int    _dxfIndexEdges(int *nE, EdgeS **edges, Htable *e_table, int nT, int *t);

void   _dxfInitFather(int n, int *father);


/* perform path conpression on the representatives of a forest data structure
   the macro FATHER is more efficient than the _dxfFather(i,f) procedure
   it calls the procedure only if necessary */

#define FATHER(i, f) (((i) == (f)[i]) ? (i) : \
                      ((f)[i] == (f)[(f)[i]]) ? (f)[i] : _dxfFather(i,f))

int _dxfFather(int i, int *father);

int _dxfJoin(int i, int j, int *father);



/*...........................................................................*/
/* 
	default parameters for simplification
*/
/*...........................................................................*/


#define SIMPLIFY_MAX_ANGLE_NOR      1.5     /* maximum angle between normals of corresponding triangles
					       before and after an edge-collapse */
#define SIMPLIFY_COMPACTNESS_RATIO .2       /* maximum ratio between the minimum compactness
					       in the vertex star after simplification and the
					       minimum compactness in the edge star before simplification */
#define SIMPLIFY_VALENCE_MAX        30      /* maximum valence at a vertex after simplification */
#define SIMPLIFY_VALENCE_MAX1       31
#define SIMPLIFY_VALENCE_MAX4       34

#define FOUR_SQRT_3 6.928203230275508
#define EPS_VERTEX_OUTSIDE_TRIANGLE -1e5   /* threshold on the barycentric coordinate such that if 
					      a barycentric coordinate is less than that the corresponding
					      point is decided to be on the other side of the edge */

/* typedef for some simple geometric entities: */

typedef float  Vertex[3];
typedef float  Plane[4];
typedef double VertexD[3];
typedef int    Face[3];

#define DIRECTION_COUNTER_CLOCKWISE 2
#define DIRECTION_CLOCKWISE         1

/*...........................................................................*/

/*
  data structures for binary heaps
*/

/*...........................................................................*/


typedef struct _Heap_ {
  int		n,last;
  int		*perm;
  int		*invp;
  float		*key;
  
} Heap;

#define LeftSon(i)	2*(i)+1
#define RightSon(i)	2*(i)+2

Heap          *_dxfNewHeap(int n);
unsigned long  _dxfHeapSize(int n);
void           _dxfResetHeap(Heap *heap);
int            _dxfHeapAdd(float k, Heap *heap);
int            _dxfHeapDelMin(Heap *heap);
int            _dxfHeapDelete(int i, Heap *heap);

#define              HeapLength(heap)       (heap)->last
#define              HeapLengthMax(heap)    (heap)->n
#define              HeapFull(heap)         (HeapLength(heap) == HeapLengthMax(heap))
#define              HeapValidIndex(heap, index) \
                           ((index) < HeapLengthMax(heap) && (index) >= 0)
   
     /* you can't remove these elements from the heap  but you can add them in: */

#define              HeapOutsideIndex(heap) HeapLengthMax(heap) + 1
#define              HeapInitialIndex(heap) HeapLengthMax(heap)

  /* you can neither add nor remove these elements */

#define              HeapInvalidIndex(heap) -1


/* internal structure used for simplification: */

typedef struct _SimpData_ {
  
  /* data marked "to be allocated" is managed by the 

     _dxfCreateSimplification....

     _dxfFreeSimplification....

     routines. it is used only internally to the _dxfSimplifySurfaceDriver routine and is not exported

     */

  Heap   
        *edge_heap;                  /* to be allocated */

  int   
        nV, nT, nE, 
        data_dim, 
        *triangles,
        *edge2index,                 /* to be allocated */
        *index2edge,                 /* to be allocated */
        *vertex_father,
        *edge_father,                /* to be allocated */
        *triangle_father,
        *vertex_lut,
        preserve_volume,
	simplify_boundary,
	preserve_boundary_length,
        valence_max,

    /* various counters: */

        num_edg_remaining_4_testing,
        num_edg_collapsed,
        num_edg_weights,
        edg_rejected_4_topology,
        edg_rejected_4_geometry,
        edg_rejected_4_tolerance,
        last_edge_added2heap,     /* edge number for the last edge added to the heap */
        heap_initial_index,
        heap_outside_index;
 
  EdgeS 
        *edge_array;

  Htable 
        *e_table;

  Vertex
        *vert,
        *normal,                   /* to be allocated */
        simplified_vertex;
      
  u_char 
        *boundary_vert;            /* to be allocated */

  u_short 
        *valence;                  /* to be allocated */
  
  float 
        tolerance,
        data_tolerance,
        *err_volume,               /* to be allocated */
        *tol_volume,               /* to be allocated */
        *old_face_areas,
        *area,                     /* to be allocated */
        *compactness,              /* to be allocated */ 
        min_scalprod_nor,
        compactness_ratio,
        *vx_data,                  /* to be allocated */
        *vx_data_error,            /* to be allocated */
        *vx_data_potential_values, /* to be allocated */
                                   /* potential data values at the simplified vertex */
        *vx_old_data,
        *vx_new_data,              /* necessary to compute the discrepancy in data values
				      at the vertices (or corners) of the mapping between
				      original and simplified surfaces */

        vx_data_potential_err;     /* potential data error  at the simplified vertex */

  int   (*get_lowest_weight_edge)(struct _SimpData_ *); /* pointer to a function used
							   for extracting the edge
							   with lowest weight from
							   the heap */

  int   (*add_edge)              (struct _SimpData_ *, int);
  int   (*remove_edge)           (struct _SimpData_ *, int, int);

  /* other buffers used for storing information relative to vertex and edge stars */

  int 
      *vertex_star_buffer;

  Pointer 
      edge_star_buffer_fl,
      edge_star_buffer_vx;
 
} SimpData;



int  _dxfRemoveLowestWeightEdgeFromHeap(SimpData *simp_data);

int  _dxfAddEdge2Heap(SimpData *simp_data, int the_edge);

int  _dxfRemoveEdgeFromHeap(SimpData *simp_data, int index, int the_edge);

int  _dxfBuildSimplifiedSurface(SimpData *simp_data, int *new_nV, float **new_v, 
			       int *new_nT, int **new_t,
			       float **new_face_areas, float **face_normals,
			       float **new_positional_error, float **new_vertex_data);

int  _dxfSimplifyManifoldSurface(SimpData *simp_data);

int  _dxfCreateSimplificationDataStructure(SimpData *simp_data, float *data, float *old_pos_err);

/* free all the simplification data structures: */

int  _dxfFreeSimplificationDataStructure(SimpData *simp_data, float *data, float *old_pos_err);

/* free only the simplification data structures that are not exported by the driver routine,
   to avoid running out of memory in _dxfBuildSimplifiedSurface() */

int  _dxfPartialFreeSimplificationDataStructure(SimpData *simp_data);

int  _dxfTrianglesNormalAreaCompactness( int nT, Face *t, Vertex *v, Vertex *t_normal,
					float *t_area, float *t_compactness);

int  _dxfTriangleNormalQR2D(VertexD *tri, double *n);

int  _dxfTriangleNormalQR2(Vertex *tri, Vertex n);

int  _dxfVectorProductQRD(VertexD x1, VertexD x2, double *n);

int  _dxfFlagBoundaryVertices(int nE, EdgeS *edges, u_char *boundary_vert);
 
int  _dxfMarkEdgesAdjacent2Boundary(SimpData *simp_data);

int  _dxfBuildEdgeHeap(SimpData *simp_data);

void _dxfInitArray(char *array, char *data, int n, int size);

int  _dxfCollapsibilityTestsBoundaryEdge2ndKind(SimpData *simp_data, 
					       int edge_num, int v0, int v1, int val0, int val1, 
					       int move_simplified_vertex);

int  _dxfDirectedParentVertexStar(int parent_vertex, int first_triangle,
				 int valence, int *star, SimpData *simp_data, int direction);

int  _dxfRotateParentTriangle(int  parent_vertex, int *vertex_fathers, 
			     int *tri_v, int *tri_v_parents, int *tri_v_rotated);

int  _dxfRotateParentEdge(int parent_vertex, int *vertex_fathers,
			 int *edge_v, int *edge_t, int *edge_v_parents, 
			 int *edge_v_rotated, int *edge_t_rotated);

/* macros used for finding the next edge when rotating around a vertex in
   on a simplified surface */



#define _dxfNextParentEdge(oriented_tri, directionCW)                                                     \
	 _dxfFather(                                                                                      \
	  (int) (((EdgeS *)_dxfFindEdge(simp_data->e_table,(oriented_tri)[0],(oriented_tri)[directionCW]))\
		 -simp_data->edge_array), simp_data->edge_father)



#define _dxfNextParentTriangle(oriented_edge_t, directionCW)                                \
             ((directionCW) == 1)? FATHER((oriented_edge_t)[1], simp_data->triangle_father):\
				   FATHER((oriented_edge_t)[0], simp_data->triangle_father)

int   _dxfManifoldLinkTest(int *link0, int *link1, int val0, int val1);

int   _dxfCmpIntSmallFirst(char *s1, char *s2);

float _dxfClosestPointOnEdge(Vertex X, Vertex XP, Vertex A, Vertex B, 
			    float *t);/* barycentric coordinates: t, 1-t */

void  _dxfMakeEdgeBarycenter(Vertex v0, Vertex v1, float alpha0, Vertex v);

float _dxfNormalize3Vector(Vertex v);

int   _dxfSolveQuadraticEquation(double A, double B, double C, 
	 	 	        double *sol1, double *sol2, double eps);

int   _dxfBoundaryCollapse2KndGeometricallyOk(SimpData *simp_data, int v0, int v1,
			int *star0, int *vstar0, int val0, Vertex *s_normal0, float *s_area0,
			int direction, float min_scalprod_nor, float compactness_ratio);

float _dxfFastCompactness3DTriangle2(Vertex *tri, float *triangle_area);



/*+--------------------------------------------------------------------------+
  |                                                                          |
  |   _dxfTriangleBasisVectors                                               |
  |                                                                          |
  +--------------------------------------------------------------------------+*/



#define  _dxfTriangleBasisVectors(tri, x1, x2, origin)                         \
{                                                                              \
  /* compute two basis vectors of a triangle, such that the first basis vector \
     will correspond to the shortest edge,                                     \
     return the number of the vertex that is chosen as the triangle origin */  \
                                                                               \
                                                                               \
  int the_origin = 0;                                                          \
                                                                               \
  /* I- find the shortest edge of the triangle                                 \
                   2                                                           \
                  /\                                                           \
          edg1   /  \ edg0                                                     \
                /____\                                                         \
               0 edg2 1 */                                                     \
                                                                               \
  register u_char j;                                                           \
  double min_len = SQUARED_DISTANCE3((tri)[0],(tri)[2]);                       \
  u_char shortest_edge = 1;                                                    \
                                                                               \
  for (j=0;j<2;j++)                                                            \
    {                                                                          \
      double len = SQUARED_DISTANCE3((tri)[j],(tri)[j+1]);                     \
                                                                               \
      if (len < min_len)                                                       \
	{                                                                      \
	  min_len = len;                                                       \
	  shortest_edge = (j+2)%3;                                             \
	}                                                                      \
      }                                                                        \
                                                                               \
  the_origin = (shortest_edge + 1) %3;                                         \
                                                                               \
  (origin) = the_origin;                                                       \
                                                                               \
  /* II- assign the smallest edge to one of the vectors x1 and x2 */           \
  {                                                                            \
    u_char                                                                     \
      dest1 =  (shortest_edge+2)%3,                                            \
      dest2 =   shortest_edge;                                                 \
                                                                               \
    register u_char k;                                                         \
                                                                               \
    for(k=0;k<3;k++)                                                           \
      {                                                                        \
	(x1)[k] = (tri)[dest1][k] - (tri)[the_origin][k];                      \
	(x2)[k] = (tri)[dest2][k] - (tri)[the_origin][k];                      \
      }                                                                        \
  }                                                                            \
										 }

int _dxfErrorWithinToleranceVBoundary2ndKnd(SimpData *simp_data, int v0, int v1, 
					   int *vstar0, int *vstar1, int val0, int val1);


#define ExtractTriangleFromVertexStar( vstar, i) a_triangle; {                 \
  memcpy(a_triangle[1], simp_data->vert[(vstar)[(i)]],   sizeof(Vertex));      \
  memcpy(a_triangle[2], simp_data->vert[(vstar)[(i)+1]], sizeof(Vertex));}

/* macro for taking the minimum of the potential errors at the pivot vertex,
   also called the *simplified_vertex in my article */



#define UPDATE_ERR_PIVOT                                                       \
{	                                                                       \
              if (potential_err_pivot > largest_err_pivot)                     \
		largest_err_pivot = potential_err_pivot;                       \
	       edge_collapsible = (largest_err_pivot <= tol_pivot);}



#define    _dxfErrorTestDataOnSurfaceInit(data)                                                         \
{                                                                                                       \
 /* initialize the potential error value for the data at the potential simplified vertex location */    \
	       if ((data)->vx_data) {(data)->vx_data_potential_err = 0.0;}}



#define    _dxfErrorTestDataOnSurfaceEnd(data,v)                                                        \
{                                                                                                       \
                                                                                                        \
 /* replace the data error value at the simplified vertex with the potential data error value */        \
                                                                                                        \
    if ((data)->vx_data)                                                                                \
        {(data)->vx_data_error[(v)] = (data)->vx_data_potential_err;                                    \
                                                                                                        \
 /* replace the data value with the potential data value */                                             \
         memcpy((data)->vx_data + (v) * (data)->data_dim, (data)->vx_data_potential_values,             \
                (data)->data_dim * sizeof (float));                                                     \
													  }}


#define _dxfErrorTestDataOnSurfaceInterpolate2(data, v0, v1, t)                                             \
{                                                                                                           \
  if ((data)->vx_data)                                                                                      \
	{                                                                                                   \
	  /* interpolate the data value at the potential simplified vertex */                               \
                                                                                                            \
	  int coordinate = 0;                                                                               \
                                                                                                            \
	  float                                                                                             \
	    *data0 = (data)->vx_data + (v0) * (data)->data_dim,                                             \
	    *data1 = (data)->vx_data + (v1) * (data)->data_dim,                                             \
	    t1     = 1. -(t);                                                                               \
                                                                                                            \
	  for (coordinate = 0; coordinate < (data)->data_dim; coordinate++)                                 \
                                                                                                            \
	    (data)->vx_data_potential_values[coordinate] = data0[coordinate] * (t) + data1[coordinate] * t1;\
	}                                                                                                   \
													      }


#define _dxfErrorTestDataOnSurfaceInterpolate3(data, v0, v1, v2, t0, t1, t2)                              \
{                                                                                                         \
  if ((data)->vx_data)                                                                                    \
	{                                                                                                 \
	  /* interpolate the data value at the potential simplified vertex */                             \
                                                                                                          \
	  int coordinate = 0;                                                                             \
                                                                                                          \
	  float                                                                                           \
	    *data0 = (data)->vx_data + (v0) * (data)->data_dim,                                           \
	    *data1 = (data)->vx_data + (v1) * (data)->data_dim,                                           \
	    *data2 = (data)->vx_data + (v2) * (data)->data_dim;                                           \
                                                                                                          \
	  for (coordinate = 0; coordinate < (data)->data_dim; coordinate++)                               \
                                                                                                          \
	    (data)->vx_data_potential_values[coordinate] =                                                \
                  data0[coordinate] * (t0) + data1[coordinate] * (t1) + data2[coordinate] * (t2);         \
	}                                                                                                 \
													    }

int   _dxfErrorTestDataOnSurfaceUpdate(SimpData *simp_data, int new_v, int old_v, 
				      int new_v2, int new_v3,
				      float alpha_sv, float alpha_n_v2, float alpha_n_v3,
				      int old_v1, int old_v2, int old_v3, 
				      float alpha_o_v1, float alpha_o_v2, float alpha_o_v3);

float  _dxfSegmentIntersection(Vertex A, Vertex B, Vertex C, Vertex D,
			      float *lambda, float *mu, 
			      int changed_AorB, int changed_C, int changed_D);

int    _dxfHouseholderPreMultiplication(float *A, int mA, float *v, int m, int n, float *w );

int    _dxfHouseholderPreMultiplicationDbl(double *A, int mA, double *v, int m, int n, double *w );

int    _dxfCollapseBoundaryEdge2ndKnd(SimpData *simp_data, int edg_num, int v0, int v1, int val0, 
				     int val1, int *star0, int *star1, int *vstar0, int *vstar1, 
				     int *estar0, int *estar1, Vertex *s_normal0, Vertex *s_normal1, 
				     float *s_area0, float *s_area1);

int    _dxfRemoveEdgesFromHeap(int *estar, u_short val, int *edge2index, SimpData *simp_data);

int    _dxfReinstateEdgesInHeap(int *estar, u_short val, SimpData *simp_data, 
			       int *num_edg_added);

int    _dxfCollapsibilityTestsBoundaryEdge1stKind(SimpData *simp_data, 
						 int edge_num, int v0, int v1, int val0, int val1);

int    _dxfManifoldLinkTest1stKnd(int v0, int val0, int *star1, int *link1, int val1, 
				 SimpData *simp_data);

int    _dxfBoundaryCollapse1stKndGeometricallyOk(SimpData *simp_data, int v0, int *star0, int *vstar0, 
						int val0, Vertex *s_normal0, float *s_area0,
							float min_scalprod_nor, float compactness_ratio);

int    _dxfErrorWithinToleranceVBoundary1stKnd(SimpData *simp_data, int v0, int v1, int *vstar1, 
					      int val1);

int    _dxfClosestPointOnTriangle(Vertex *tri, Vertex w, Vertex wp, Vertex bary, float *dist);

int    _dxfTriangle3DBarycentricCoordinates2(Vertex *tri, Vertex w, Vertex wp, Vertex bary, 
					    float *residual);

void   _dxfMakeBarycenter(int nV, Vertex *vv, Vertex bary, float *bary_coord);

int    _dxfSolveSymmetric2x2Eqn(double a,  double b, double d, double e, double f,
			       double *x, double *y);


int    _dxfCollapseBoundaryEdge1stKnd(SimpData *simp_data, int edg_num, int v0, int v1, int val0, 
				     int val1, int *star1, int *vstar1, int *estar1,
				     Vertex *s_normal1, float *s_area1);

int    _dxfCollapseTopologicallyFeasible(SimpData *simp_data, int *vstar0, int *vstar1, 
					u_short val0, u_short val1);

void   _dxfBuildParentEdgeStars(EdgeS *edg, int v0, int v1, u_short val0,
			       u_short val1, int *star0, int *star1, SimpData *s);

void   _dxfParentVertexStar(int vf, int t0, u_short val, int *star, SimpData *simp_data);

void   _dxfApplyTranslation(int nV, Vertex *vv, Vertex T);

void   _dxfOppositeVector(Vertex u, Vertex v);

void   _dxfCopyFaceNormalsAreasCompactness(Plane *plane, float *s_area, float *s_comp, 
					  int *star0, int *star1, int val0, int val013, 
					  Vertex *t_normal, float *t_area, float *t_comp);

float  _dxfMinFloatArray(float *array, int n);

void   _dxfSimplifiedVertexLocation(SimpData *simpdata, u_short val0, u_short val1, u_short star_val,
				   Vertex *v, Face   *f, Plane  *p, float *area, float *w, 
				   u_char method);

double _dxfStarPlaneEquationsAndVolume(Vertex *star_vert,  Face *star_face, Plane  *star_plane, 
				      float *star_areas, int   star_val);

void   _dxfComposeVectors(Vertex u, Vertex v, float lambda, Vertex w);

int    _dxfPositionVertexLSConstraints2(Plane *star_plane, u_short star_val, Vertex simplified_vertex);

void   _dxfFPlaneForIdenticalVolume(u_short val0, u_short val1, Plane p, Vertex *star_vert, 
				   Face *star_face, float star_volume);

int    _dxfCollapseGeometricallyFeasible(Vertex *svert, Face *sface, 
					Plane *splane, float *old_comp, float new_comp_min,
					Vertex *snor0, Vertex *snor1,   
					u_short val0, u_short sval, 
					float min_scalprod, float compactness_ratio);

int    _dxfNoFlipOverCheck(u_short first_face, u_short last_face, 
			  Vertex *nor, Plane *plane, float min_scalprod_nor);

int    _dxfErrorWithinToleranceV(SimpData *simp_data, Vertex *star_vert, Face *star_face,
				Plane *star_plane, int *star0, int *star1, int val0, int val1);

void   _dxfSimplifiedStarNormals(Face *sface, Vertex *svert, Vertex *snor, float *sarea, float *scomp, 
				Vertex simpvert, u_short first_face, u_short last_face);

int    _dxfCollapseEdge(SimpData *simp_data, int v0, int v1, int *star0, int *star1, 
		       u_short val0, u_short val1, Vertex *nor0, Vertex *nor1, float *area0, 
		       float *area1, float *comp0, float *comp1);

void   _dxfReplaceFaceNormals(int *star, Vertex *new_nor, float *new_area, float *new_comp, 
			     Vertex *nor, float *area, float *comp, u_short val);

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#endif  /* _SIMPLESURF_H_ */

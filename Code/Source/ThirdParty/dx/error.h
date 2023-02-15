/***********************************************************************/
/* Open Visualization Data Explorer                                    */
/* (C) Copyright IBM Corp. 1989,1999                                   */
/* ALL RIGHTS RESERVED                                                 */
/* This code licensed under the                                        */
/*    "IBM PUBLIC LICENSE - Open Visualization Data Explorer"          */
/***********************************************************************/


/* TeX starts here. Do not remove this comment. */

#ifndef _DXI_ERROR_H_
#define _DXI_ERROR_H_

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

/*
\section{Error handling}
\label{errorsec}

In general, the routines in the Data Explorer library return either a
pointer, such as an object reference, or an integer error indication.
Success is indicated by returning a non-null pointer or by returning
the non-zero integer constant {\tt OK}.  Failure is indicated by
returning {\tt NULL}, or by returning {\tt ERROR} (which is defined to
be zero).

In case of failure, the library routine may have set the error code by
calling {\tt DXSetError()}, in which case the calling routine generally
should just return {\tt NULL} or {\tt ERROR} to propagate the error
up.  On the other hand, the library routine may not have set the error
code.  In this case it is up to the calling routine to decide whether
the null return indicates an error, in which case the calling routine
should set the error code by calling {\tt DXSetError()}, and then return
{\tt NULL} or {\tt ERROR}; or whether the null return was not an
error, in which case the calling routine should proceed.  This manual
documents for each routine whether it sets the error code when it
returns null.

The error codes are defined as follows: \index{errorcodes}
*/

typedef enum errorcode {
    ERROR_NONE,
    ERROR_INTERNAL,
    ERROR_UNEXPECTED,
    ERROR_ASSERTION,
    ERROR_NOT_IMPLEMENTED,
    ERROR_NO_MEMORY,
    ERROR_BAD_CLASS,
    ERROR_BAD_TYPE,
    ERROR_NO_CAMERA,
    ERROR_MISSING_DATA,
    ERROR_DATA_INVALID,
    ERROR_BAD_PARAMETER,
    ERROR_NO_HARDWARE_RENDERING,
    ERROR_MAX
} ErrorCode;

typedef int Error;
#ifndef ERROR
#define ERROR 0
#endif
#ifndef OK
#define OK 1
#endif

typedef void *Pointer;

#ifndef NULL
#define NULL 0
#endif

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#endif /* _DXI_ERROR_H_ */

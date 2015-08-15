/* Fortran to C interface to
 *     PO_metis_EdgeND,
 *     PO_metis_NodesND,
 *     PO_METIS_NodeWND 
 */

#include <metis.h>
#include "mph_part_metis_struct.h"
#include "mph_part_metis_proto.h"

/* PO_METIS_EdgeND */
void PO_METIS_EDGEND(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options,
		     idxtype *perm, idxtype *iperm, int *nbdom, SepInfoType *gbpartinfo)
{
  PO_METIS_EdgeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm, nbdom, gbpartinfo);
}
void po_metis_edgend(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options,
		     idxtype *perm, idxtype *iperm, int *nbdom, SepInfoType *gbpartinfo)
{
  PO_METIS_EdgeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm, nbdom, gbpartinfo);
}
void po_metis_edgend_(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options,
		      idxtype *perm, idxtype *iperm, int *nbdom, SepInfoType *gbpartinfo)
{
  PO_METIS_EdgeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm, nbdom, gbpartinfo);
}
void po_metis_edgend__(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options,
		       idxtype *perm, idxtype *iperm, int *nbdom, SepInfoType *gbpartinfo)
{
  PO_METIS_EdgeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm, nbdom, gbpartinfo);
}

/* PO_METIS_NodeND */
void PO_METIS_NODEND(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options,
		     idxtype *perm, idxtype *iperm, int *nbdom, SepInfoType *gbpartinfo)
{
  PO_METIS_NodeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm, nbdom, gbpartinfo);
}
void po_metis_nodend(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options,
		     idxtype *perm, idxtype *iperm, int *nbdom, SepInfoType *gbpartinfo)
{
  PO_METIS_NodeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm, nbdom, gbpartinfo);
}
void po_metis_nodend_(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options,
		      idxtype *perm, idxtype *iperm, int *nbdom, SepInfoType *gbpartinfo)
{
  PO_METIS_NodeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm, nbdom, gbpartinfo);
}
void po_metis_nodend__(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options,
		       idxtype *perm, idxtype *iperm, int *nbdom, SepInfoType *gbpartinfo)
{
  PO_METIS_NodeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm, nbdom, gbpartinfo);
}


/* PO_METIS_NodeWND  */
void PO_METIS_NODEWND(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, int *numflag, int *options,
		      idxtype *perm, idxtype *iperm, int *nbdom, SepInfoType *gbpartinfo)
{
  PO_METIS_NodeWND(nvtxs, xadj, adjncy, vwgt, numflag, options, perm, iperm, nbdom, gbpartinfo);
}
void po_metis_nodewnd(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, int *numflag, int *options,
		      idxtype *perm, idxtype *iperm, int *nbdom, SepInfoType *gbpartinfo)
{
  PO_METIS_NodeWND(nvtxs, xadj, adjncy, vwgt, numflag, options, perm, iperm, nbdom, gbpartinfo);
}
void po_metis_nodewnd_(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, int *numflag, int *options,
		       idxtype *perm, idxtype *iperm, int *nbdom, SepInfoType *gbpartinfo)
{
  PO_METIS_NodeWND(nvtxs, xadj, adjncy, vwgt, numflag, options, perm, iperm, nbdom, gbpartinfo);
}
void po_metis_nodewnd__(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, int *numflag, int *options,
			idxtype *perm, idxtype *iperm, int *nbdom, SepInfoType *gbpartinfo)
{
  PO_METIS_NodeWND(nvtxs, xadj, adjncy, vwgt, numflag, options, perm, iperm, nbdom, gbpartinfo);
}






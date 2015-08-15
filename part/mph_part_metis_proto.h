#ifndef MPH_PART_METIS_PROTO_H
#define MPH_PART_METIS_PROTO_H 1

/* ometis.c */
#define PO_MlevelNestedDissection		__PO_MlevelNestedDissection

void PO_METIS_EDGEND(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, int *,SepInfoType *); 
void po_metis_edgend(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, int *,SepInfoType *); 
void po_metis_edgend_(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, int *,SepInfoType *); 
void po_metis_edgend__(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, int *,SepInfoType *); 

void PO_METIS_NODEND(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, int *,SepInfoType *); 
void po_metis_nodend(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, int *,SepInfoType *); 
void po_metis_nodend_(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, int *,SepInfoType *); 
void po_metis_nodend__(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, int *,SepInfoType *); 

void PO_METIS_NODEWND(int *, idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, int *,SepInfoType *); 
void po_metis_nodewnd(int *, idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, int *,SepInfoType *); 
void po_metis_nodewnd_(int *, idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, int *,SepInfoType *); 
void po_metis_nodewnd__(int *, idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, int *,SepInfoType *); 

void PO_METIS_EdgeND(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, int *,SepInfoType *); 
void PO_METIS_NodeND(int *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, int *,SepInfoType *); 
void PO_METIS_NodeWND(int *, idxtype *, idxtype *, idxtype *, int *, int *, idxtype *, idxtype *, int *,SepInfoType *); 

void PO_MlevelNestedDissection(CtrlType *, GraphType *, idxtype *, float, int, int *, int, int, int *,SepInfoType *);

#endif






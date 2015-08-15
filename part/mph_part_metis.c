#include <metis.h>
#include "mph_part_metis_struct.h"
#include "mph_part_metis_proto.h"


/*************************************************************************
* This function is the entry point for OEMETIS
**************************************************************************/
void PO_METIS_EdgeND(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, 
                  idxtype *perm, idxtype *iperm, int *nbdom, SepInfoType *gbpartinfo ) 
{
  int i, j, allsep,keptsep,fatherid,sepcnter;
/*  SepInfoType *gbpartinfo;*/

  GraphType graph;
  CtrlType ctrl;

#if MAPHYS_DEBUG
   fprintf(stdout," \n");
   fprintf(stdout,"===============================================\n");
   fprintf(stdout,"               PO_METIS_EdgeND  \n");
   fprintf(stdout,"===============================================\n");
#endif

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  SetUpGraph(&graph, OP_OEMETIS, *nvtxs, 1, xadj, adjncy, NULL, NULL, 0);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = OEMETIS_CTYPE;
    ctrl.IType = OEMETIS_ITYPE;
    ctrl.RType = OEMETIS_RTYPE;
    ctrl.dbglvl = OEMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }
  ctrl.oflags  = 0;
  ctrl.pfactor = -1;
  ctrl.nseps   = 1;

  ctrl.optype = OP_OEMETIS;
  ctrl.CoarsenTo = 20;
  ctrl.maxvwgt = 1.5*(idxsum(*nvtxs, graph.vwgt)/ctrl.CoarsenTo);

  InitRandom(-1);

  AllocateWorkSpace(&ctrl, &graph, 2);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

    keptsep  = 0;
    keptsep  = *nbdom-1;
    /*fprintf(stdout,"keptsep  :%10d%10d\n",*nbdom,keptsep);*/
    if(keptsep<=0) return;
    allsep   = 0;
    fatherid = 0;
    sepcnter = 0;

  PO_MlevelNestedDissection(&ctrl, &graph, iperm, ORDER_UNBALANCE_FRACTION, *nvtxs, &allsep, keptsep,fatherid,&sepcnter,gbpartinfo);

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  for (i=0; i<*nvtxs; i++)
    perm[iperm[i]] = i;

  FreeWorkSpace(&ctrl, &graph);

  if (*numflag == 1)
    Change2FNumberingOrder(*nvtxs, xadj, adjncy, perm, iperm);
}





/*************************************************************************
* This function is the entry point for ONCMETIS
**************************************************************************/
void PO_METIS_NodeND(int *nvtxs, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, 
                  idxtype *perm, idxtype *iperm, int *nbdom, SepInfoType *gbpartinfo ) 
{
  int i, ii, j, l, wflag, nflag,allsep,keptsep,fatherid,sepcnter,debug;

/*  SepInfoType *gbpartinfo;*/
  idxtype *gucmpdispl;


  GraphType graph;
  CtrlType ctrl;
  idxtype *cptr, *cind, *piperm;

#if MAPHYS_DEBUG
   fprintf(stdout," \n");
   fprintf(stdout,"===============================================\n");
   fprintf(stdout,"               METIS_NodeND  \n");
   fprintf(stdout,"===============================================\n");
#endif

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType   = ONMETIS_CTYPE;
    ctrl.IType   = ONMETIS_ITYPE;
    ctrl.RType   = ONMETIS_RTYPE;
    ctrl.dbglvl  = ONMETIS_DBGLVL;
    ctrl.oflags  = ONMETIS_OFLAGS;
    ctrl.pfactor = ONMETIS_PFACTOR;
    ctrl.nseps   = ONMETIS_NSEPS;
  }
  else {
    ctrl.CType   = options[OPTION_CTYPE];
    ctrl.IType   = options[OPTION_ITYPE];
    ctrl.RType   = options[OPTION_RTYPE];
    ctrl.dbglvl  = options[OPTION_DBGLVL];
    ctrl.oflags  = options[OPTION_OFLAGS];
    ctrl.pfactor = options[OPTION_PFACTOR];
    ctrl.nseps   = options[OPTION_NSEPS];
  }
  if (ctrl.nseps < 1)
    ctrl.nseps = 1;

  ctrl.optype = OP_ONMETIS;
  ctrl.CoarsenTo = 100;

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  InitRandom(-1);

  /*=============================================================
  * AZZAM MODIF FOR THIS VERSION
  --=============================================================*/
   /*ctrl.oflags  = 0;*/
  if (ctrl.pfactor > 0) { 
   ctrl.pfactor = -1;
   fprintf(stdout," \n");
   fprintf(stdout,"===============================================\n");
   fprintf(stdout,"         pfactor forced to -1 in this version  \n");
   fprintf(stdout,"===============================================\n");
   }

 /*
 if (ctrl.oflags&OFLAG_COMPRESS){
      ctrl.oflags--; 
   fprintf(stdout," \n");
   fprintf(stdout,"===============================================\n");
   fprintf(stdout,"    test  ctrl.oflags forced to 0 in this version  \n");
   fprintf(stdout,"===============================================\n");
   }*/
  /*=============================================================
  * END AZZAM MODIF FOR THIS VERSION
  --=============================================================*/


  if (ctrl.pfactor > 0) { 
    /*============================================================
    * Prune the dense columns
    ==============================================================*/
    fprintf(stdout,"Prune the dense columns \n");
    piperm = idxmalloc(*nvtxs, "ONMETIS: piperm");

    PruneGraph(&ctrl, &graph, *nvtxs, xadj, adjncy, piperm, (float)(0.1*ctrl.pfactor));
  }  
  else if (ctrl.oflags&OFLAG_COMPRESS) {
    /*============================================================
    * Compress the graph 
    ==============================================================*/
    fprintf(stdout,"Compress the graph \n");

    cptr = idxmalloc(*nvtxs+1, "ONMETIS: cptr");
    cind = idxmalloc(*nvtxs, "ONMETIS: cind");

    CompressGraph(&ctrl, &graph, *nvtxs, xadj, adjncy, cptr, cind);

    if (graph.nvtxs >= COMPRESSION_FRACTION*(*nvtxs)) {
      ctrl.oflags--; /* We actually performed no compression */
      GKfree(&cptr, &cind, LTERM);
      fprintf(stdout,"We actually performed no compression \n");
    }else if (2*graph.nvtxs < *nvtxs && ctrl.nseps == 1){
      ctrl.nseps = 5;
      fprintf(stdout,"ctrl.nseps = %10d \n",ctrl.nseps );}
  }
  else {
    fprintf(stdout,"SetUpGraph \n");
    SetUpGraph(&graph, OP_ONMETIS, *nvtxs, 1, xadj, adjncy, NULL, NULL, 0);
  }


  /*=============================================================
  * Do the nested dissection ordering 
  --=============================================================*/
  ctrl.maxvwgt = 1.5*(idxsum(graph.nvtxs, graph.vwgt)/ctrl.CoarsenTo);
  AllocateWorkSpace(&ctrl, &graph, 2);

  if (ctrl.oflags&OFLAG_CCMP) {
    fprintf(stdout,"METIS_NodeND: coucou Azzam1 calling MlevelNestedDissectionCC \n");
    MlevelNestedDissectionCC(&ctrl, &graph, iperm, ORDER_UNBALANCE_FRACTION, graph.nvtxs);
  }
  else {
    /*fprintf(stdout,"METIS_NodeND: coucou Azzam2 calling PO_MlevelNestedDissection  \n");*/
    /*keptsep  = 31;*/
    keptsep  = 0;
    keptsep  = *nbdom-1;
    /*fprintf(stdout,"keptsep  :%10d%10d\n",*nbdom,keptsep);*/
    if(keptsep<=0) return;
    allsep   = 0;
    fatherid = 0;
    sepcnter = 0;


/*
    gbpartinfo = (SepInfoType *)calloc(keptsep,sizeof(SepInfoType) );
    if( ! gbpartinfo )
    {
    fprintf( stdout, "Erreur a l'allocation de gbpartinfo\n\n" );
    exit(1);
    }
*/
#if MAPHYS_DEBUG
      fprintf(stdout,"AZZAM : FIND %d SEPARATOR AT EACH STEP \n",ctrl.nseps);
#endif
    PO_MlevelNestedDissection(&ctrl, &graph, iperm, ORDER_UNBALANCE_FRACTION, \
                   graph.nvtxs, &allsep,keptsep,fatherid,&sepcnter,gbpartinfo);


  /*=============================================================
  * print out the structure constructed 
  --=============================================================*/
 /*
    fprintf(stdout,"=======================================================================================================\n"); 
    fprintf(stdout," number of subdomains is equal to :%10d\n",*nbdom);
    fprintf(stdout,"=======================================================================================================\n"); 
    fprintf(stdout,"%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s\n","Id","father","rson","lson","sepst","seped","rgst","rged","lgst","lged"); 
    fprintf(stdout,"=======================================================================================================\n"); 
    for (i=0; i<keptsep; i++){
    fprintf(stdout,"%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d\n",\
            i+1,gbpartinfo[i].father, gbpartinfo[i].rson, gbpartinfo[i].lson, gbpartinfo[i].sepst, gbpartinfo[i].seped, \
                            gbpartinfo[i].rgst, gbpartinfo[i].rged, gbpartinfo[i].lgst, gbpartinfo[i].lged);
    }
    fprintf(stdout,"======================================================================================================\n"); 
*/
  /*=============================================================
  * end of printing out the structure constructed 
  --=============================================================*/


/*############################## */
  } /* END IF */
/*############################## */


  /*fprintf(stdout,"LAST INTERIOR INDEX: %d\n",allsep);*/

  FreeWorkSpace(&ctrl, &graph);

  if (ctrl.pfactor > 0) { /* Order any prunned vertices */
    fprintf(stdout,"ctrl.pfactor > 0  Order any prunned vertices\n");
    if (graph.nvtxs < *nvtxs) { 
      idxcopy(graph.nvtxs, iperm, perm);  /* Use perm as an auxiliary array */
      for (i=0; i<graph.nvtxs; i++)
        iperm[piperm[i]] = perm[i];
      for (i=graph.nvtxs; i<*nvtxs; i++)
        iperm[piperm[i]] = i;
    }

    GKfree(&piperm, LTERM);
  }
  else if (ctrl.oflags&OFLAG_COMPRESS) { /* Uncompress the ordering */
    fprintf(stdout,"ctrl.pfactor < 0  \n");
    if (graph.nvtxs < COMPRESSION_FRACTION*(*nvtxs)) { 
       fprintf(stdout,"ctrl.pfactor < 0  Uncompress the ordering \n");

      /* allocate memory for a pointer that point the new displacement 
      in iperm of each index in the compressed graph. in other term
      each vertex displacement in my data structure has a new value 
      after uncompression so find this new value using this vector
       I USE THE C NUMBERING HERE starting from 0                */
      gucmpdispl = (idxtype *)calloc((graph.nvtxs+1),sizeof(idxtype) );
      if( ! gucmpdispl )
      {
      fprintf( stdout, "Erreur a l'allocation de graph uncompress displacement \n\n" );
      exit(1);
      }
      /* construct perm from iperm */
      for (i=0; i<graph.nvtxs; i++)
        perm[iperm[i]] = i; 

      for (l=ii=0; ii<graph.nvtxs; ii++) {
        gucmpdispl[ii]=l;

        i = perm[ii];
        for (j=cptr[i]; j<cptr[i+1]; j++){
          iperm[cind[j]] = l++;
         /* fprintf(stdout,"ii,cptr[i],cptr[i+1],cind[j],iperm[cind[j] gucmpdispl(ii)]%10d [ %10d%10d ] %10d%10d {%10d} \n",\
                 ii,cptr[i],cptr[i+1],cind[j],iperm[cind[j]],gucmpdispl[ii]);*/
         /* getchar();*/
        }
      } 
      gucmpdispl[ii]=l;
      /* endfor(ii) for creating the new displacement and new perm*/

    /* Now I should correct my data structure according 
    to the new permutation after uncompression. So I will 
    use the gucmpdispl pointer to allow me to find the 
    new displcement of my value in gbpartinfo   */
    debug=0;
    if(debug == 1){
    fprintf(stdout,"=======================================================================================================\n"); 
    fprintf(stdout,"                        NEW DATA STRUCTURE INFORMATION AZZAM :%10d\n",*nbdom);
    fprintf(stdout,"=======================================================================================================\n"); 
    fprintf(stdout,"%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s\n","Id","father","rson","lson","sepst","seped","rgst","rged","lgst","lged"); 
    fprintf(stdout,"=======================================================================================================\n"); 
    }
    for (i=0; i<keptsep; i++){
      /* find new displ for starting values st*/
      gbpartinfo[i].sepst = gucmpdispl[gbpartinfo[i].sepst-(1)] +(1); /* -(1) +(1) COZ FORTRAN NUMBERING)*/
      gbpartinfo[i].rgst  = gucmpdispl[gbpartinfo[i].rgst-(1) ] +(1);
      gbpartinfo[i].lgst  = gucmpdispl[gbpartinfo[i].lgst-(1) ] +(1);

      /* find new displ for ending values ed (that is displ of next index-1 */
      gbpartinfo[i].seped = gucmpdispl[(gbpartinfo[i].seped-(1)) +1]-1 +(1);
      gbpartinfo[i].rged  = gucmpdispl[(gbpartinfo[i].rged-(1) ) +1]-1 +(1);
      gbpartinfo[i].lged  = gucmpdispl[(gbpartinfo[i].lged-(1) ) +1]-1 +(1);
      if(debug == 1)fprintf(stdout,"%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d\n",\
            i+1,gbpartinfo[i].father, gbpartinfo[i].rson, gbpartinfo[i].lson, gbpartinfo[i].sepst, gbpartinfo[i].seped, \
                            gbpartinfo[i].rgst, gbpartinfo[i].rged, gbpartinfo[i].lgst, gbpartinfo[i].lged);

    }/* endfor for creating my new dastructure */

      if(debug == 1) fprintf(stdout,"=======================================================================================================\n"); 
    if(gbpartinfo[0].seped != *nvtxs) fprintf(stdout," ERROR INDICES PAY ATTENTION %10d%10d\n",gbpartinfo[0].seped,*nvtxs);
    else fprintf(stdout," GOOD FINISH PAY ATTENTION \n");
    free(gucmpdispl); 


    }/* endif (graph.nvtxs < COMPRESSION_FRACTION) that is graph is compressed*/

    GKfree(&cptr, &cind, LTERM);
  }

#if MAPHYS_DEBUG
  fprintf(stdout,"  FINISH  \n");
#endif

  for (i=0; i<*nvtxs; i++)
    perm[iperm[i]] = i;

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  if (*numflag == 1)
    Change2FNumberingOrder(*nvtxs, xadj, adjncy, perm, iperm);

}





/*************************************************************************
* This function is the entry point for ONWMETIS. It requires weights on the
* vertices. It is for the case that the matrix has been pre-compressed.
**************************************************************************/
void PO_METIS_NodeWND(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, int *numflag, 
                   int *options, idxtype *perm, idxtype *iperm, int *nbdom, SepInfoType *gbpartinfo ) 
{
  int i, j, tvwgt, allsep,keptsep,fatherid,sepcnter;
/*  SepInfoType *gbpartinfo;*/
  idxtype *gucmpdispl;

  GraphType graph;
  CtrlType ctrl;

   fprintf(stdout," \n");
   fprintf(stdout,"===============================================\n");
   fprintf(stdout,"               METIS_NodeWND  \n");
   fprintf(stdout,"===============================================\n");


  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  SetUpGraph(&graph, OP_ONMETIS, *nvtxs, 1, xadj, adjncy, vwgt, NULL, 2);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = ONMETIS_CTYPE;
    ctrl.IType = ONMETIS_ITYPE;
    ctrl.RType = ONMETIS_RTYPE;
    ctrl.dbglvl = ONMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }

  ctrl.oflags  = OFLAG_COMPRESS;
  ctrl.pfactor = 0;
  ctrl.nseps = 2;
  ctrl.optype = OP_ONMETIS;
  ctrl.CoarsenTo = 100;
  ctrl.maxvwgt = 1.5*(idxsum(*nvtxs, graph.vwgt)/ctrl.CoarsenTo);

  InitRandom(-1);

  AllocateWorkSpace(&ctrl, &graph, 2);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  /*=============================================================
  * Modif by Azzam 
  --=============================================================*/
    keptsep  = 0;
    keptsep  = *nbdom-1;
    if(keptsep<=0) return;
    allsep   = 0;
    fatherid = 0;
    sepcnter = 0;


    fprintf(stdout,"-------AZZAM-------- %10d\n",graph.nvtxs);

/*    MlevelNestedDissection(&ctrl, &graph, iperm, ORDER_UNBALANCE_FRACTION, \
                   graph.nvtxs);*/

    PO_MlevelNestedDissection(&ctrl, &graph, iperm, ORDER_UNBALANCE_FRACTION, \
                   graph.nvtxs, &allsep,keptsep,fatherid,&sepcnter,gbpartinfo);


  /*=============================================================
  * print out the structure constructed 
  --=============================================================*/
 /*
    fprintf(stdout,"=======================================================================================================\n"); 
    fprintf(stdout," number of subdomains is equal to :%10d\n",*nbdom);
    fprintf(stdout,"=======================================================================================================\n"); 
    fprintf(stdout,"%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s\n","Id","father","rson","lson","sepst","seped","rgst","rged","lgst","lged"); 
    fprintf(stdout,"=======================================================================================================\n"); 
    for (i=0; i<keptsep; i++){
    fprintf(stdout,"%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d\n",\
            i+1,gbpartinfo[i].father, gbpartinfo[i].rson, gbpartinfo[i].lson, gbpartinfo[i].sepst, gbpartinfo[i].seped, \
                            gbpartinfo[i].rgst, gbpartinfo[i].rged, gbpartinfo[i].lgst, gbpartinfo[i].lged);
    }
    fprintf(stdout,"======================================================================================================\n"); 
*/
  /*=============================================================
  * end of printing out the structure constructed 
  --=============================================================*/




  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  for (i=0; i<*nvtxs; i++)
    perm[iperm[i]] = i;

  FreeWorkSpace(&ctrl, &graph);

  if (*numflag == 1)
    Change2FNumberingOrder(*nvtxs, xadj, adjncy, perm, iperm);
}





/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
void PO_MlevelNestedDissection(CtrlType *ctrl, GraphType *graph, idxtype *order, float ubfactor, int lastvtx, \
                            int *allsep, int keptsep, int fatherid, int *sepcnter, SepInfoType *gbpartinfo)
{
  int i, j, nvtxs, nbnd, tvwgt, tpwgts2[2];
  idxtype *label, *bndind;
  GraphType lgraph, rgraph;
  int rkeptsep,lkeptsep,firstvtx,mysepid,debug;

  nvtxs = graph->nvtxs;

/*
   fprintf(stdout,"  \n");
   fprintf(stdout,"PO_MlevelNestedDissection  \n");
   fprintf(stdout,"voici nbnd number of bissector node before %d\n",graph->nbnd); 
*/

  /*=============================================================
  * AZZAM MODIF 
  --=============================================================*/
/*  if (graph->nvtxs <= MMDSWITCH)   {
    fprintf(stdout,"JUST REORDERING GRAPH %d\n",graph->nvtxs); 
    MMDOrder(ctrl, graph, order, lastvtx); 
   }
   return;
*/


  /* Determine the weights of the partitions */
  tvwgt = idxsum(nvtxs, graph->vwgt);
  tpwgts2[0] = tvwgt/2;
  tpwgts2[1] = tvwgt-tpwgts2[0];



  switch (ctrl->optype) {
    case OP_OEMETIS:
      MlevelEdgeBisection(ctrl, graph, tpwgts2, ubfactor);

      IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->SepTmr));
      ConstructMinCoverSeparator(ctrl, graph, ubfactor);
      IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->SepTmr));

      break;
    case OP_ONMETIS:
      MlevelNodeBisectionMultiple(ctrl, graph, tpwgts2, ubfactor);

     /*fprintf(stdout,"Nvtxs: %6d, [%6d %6d %6d]\n", graph->nvtxs, graph->pwgts[0], graph->pwgts[1], graph->pwgts[2]);*/

      IFSET(ctrl->dbglvl, DBG_SEPINFO, printf("Nvtxs: %6d, [%6d %6d %6d]\n", graph->nvtxs, graph->pwgts[0], graph->pwgts[1], graph->pwgts[2]));

      break;
  }

   /*fprintf(stdout,"voici nbnd number of bissector node after %d\n",graph->nbnd); */ 
  /* Order the nodes in the separator */
  nbnd = graph->nbnd;
  bndind = graph->bndind;
  label = graph->label;

  *allsep   = *allsep+nbnd;
  *sepcnter = *sepcnter+1;
  mysepid   = *sepcnter;
  rkeptsep  = (keptsep-1)/2;
  lkeptsep  = (keptsep-1)/2;

/*
  fprintf(stdout,"=======================================================================\n"); 
  fprintf(stdout,"%10s%10s%10s%10s%10s%10s%10s\n","allsep","keptsep","rkeptsep","lkeptsep","sepcnter","mysepid","fatherid"); 
  fprintf(stdout,"%10d%10d%10d%10d%10d%10d%10d\n",*allsep,keptsep,rkeptsep,lkeptsep,*sepcnter,mysepid,fatherid);
  fprintf(stdout,"=======================================================================\n"); 
*/




  for (i=0; i<nbnd; i++) 
   {
   /*fprintf(stdout,"sep indice[i]: %d\n",label[bndind[i]]+1); */
   /*fprintf(stdout,"%d\n",label[bndind[i]]+1);*/
   order[label[bndind[i]]] = --lastvtx;
   /*fprintf(stdout,"%10d%10d\n",nvtxs,order[label[bndind[i]]]);*/
   }
  SplitGraphOrder(ctrl, graph, &lgraph, &rgraph);


  /* Free the memory of the top level graph */
  GKfree(&graph->gdata, &graph->rdata, &graph->label, LTERM);

  /*=============================================================
  * AZZAM: Initialize my structure of father son index 
  --=============================================================*/
  /* find my father and put it into my structure then put myself into his structure */
  gbpartinfo[mysepid-1].father = fatherid;
  if (fatherid!= 0) { /* Je suis pas le premier separateur */
      if (gbpartinfo[fatherid-1].rson==0)  {/* right son not registred yet ==> register it */
         gbpartinfo[fatherid-1].rson = mysepid;
      }
      else {/* right son already enregistred ==> register left one */
         gbpartinfo[fatherid-1].lson = mysepid;
      }
  }
  /* if I have no son so put -1 on my [rl]son*/
  if(rkeptsep>0) {
     if (rgraph.nvtxs <= MMDSWITCH) gbpartinfo[mysepid-1].rson   = -1; /* Matrix too small to be partitioned */
  }else {gbpartinfo[mysepid-1].rson   = -1;}


  if(lkeptsep>0) {
     if (lgraph.nvtxs <= MMDSWITCH) gbpartinfo[mysepid-1].lson   = -1; /* Matrix too small to be partitioned */
  }else {gbpartinfo[mysepid-1].lson   = -1;}


/*fortran numbering */
  gbpartinfo[mysepid-1].sepst = lastvtx+1 ; 
  gbpartinfo[mysepid-1].seped = lastvtx+1+nbnd-1 ;
  gbpartinfo[mysepid-1].rgst  = lastvtx-rgraph.nvtxs+1;
  gbpartinfo[mysepid-1].rged  = lastvtx;
  gbpartinfo[mysepid-1].lgst  = lastvtx-rgraph.nvtxs-lgraph.nvtxs+1;
  gbpartinfo[mysepid-1].lged  = lastvtx-rgraph.nvtxs;

  debug=0;
  if(debug==1)fprintf(stdout,
		      "AZZAM Rgraph sepId fathId sepsize   : %20d%20d%10d%10d%10d\n",
		      rgraph.nvtxs,lgraph.nvtxs,mysepid,fatherid,nbnd);


  /*=============================================================
  * END Initialize structure 
  --=============================================================*/

  if(rkeptsep>0) {
     if (rgraph.nvtxs > MMDSWITCH) {
       /*fprintf(stdout,"RIGHT GRAPH CALLING ND %10d%10d%10d%10d\n",rgraph.nvtxs,lastvtx,mysepid,fatherid); */
       PO_MlevelNestedDissection(ctrl, &rgraph, order, ubfactor,lastvtx,allsep,rkeptsep,mysepid,sepcnter,gbpartinfo);}
     else {
       fprintf(stdout,"RIGHT Part small to be partitioned %d%10d%10d\n",rgraph.nvtxs,mysepid,fatherid); 
       fprintf(stdout,"\n############## Matrix too small to be partitioned ##############\n\n"); 
/*       exit(1);*/
/*       goto RIGHTORDER;*/
       firstvtx = lastvtx-rgraph.nvtxs;
       for (i=0; i<rgraph.nvtxs; i++) 
       order[rgraph.label[i]]= firstvtx++;
/*     fprintf(stdout,"%d\n",order[rgraph.label[i]] );*/
/*     order[rgraph.label[i]] = --lastvtx;*/
         }
  }
  else {
/*    RIGHTORDER:*/
    /*fprintf(stdout,"RIGHT GRAPH CALLING MMDORDER %d%10d%10d\n",rgraph.nvtxs,mysepid,fatherid); */
    MMDOrder(ctrl, &rgraph, order, lastvtx); 
    GKfree(&rgraph.gdata, &rgraph.rdata, &rgraph.label, LTERM);
  }




  if(lkeptsep>0) {
     if (lgraph.nvtxs > MMDSWITCH) {
      /* fprintf(stdout,"LEFT GRAPH CALLING ND %10d%10d%10d%10d\n",lgraph.nvtxs,lastvtx-rgraph.nvtxs,mysepid,fatherid); */
       PO_MlevelNestedDissection(ctrl, &lgraph, order, ubfactor, lastvtx-rgraph.nvtxs, allsep,lkeptsep,mysepid,sepcnter,gbpartinfo);}
     else {
       fprintf(stdout,"LEFT Part small to be partitioned %d%10d%10d\n",lgraph.nvtxs,mysepid,fatherid);
       fprintf(stdout,"\n############## Matrix too small to be partitioned ##############\n\n"); 
/*       exit(1);*/
/*       goto LEFTORDER;*/
       firstvtx = lastvtx-rgraph.nvtxs-lgraph.nvtxs;
       for (i=0; i<lgraph.nvtxs; i++) 
       order[lgraph.label[i]]= firstvtx++;
      }
  }
  else {
/*    LEFTORDER:*/
   /* fprintf(stdout,"LEFT GRAPH CALLING MMDORDER %d%10d%10d\n",lgraph.nvtxs,mysepid,fatherid); */
    MMDOrder(ctrl, &lgraph, order, lastvtx-rgraph.nvtxs); 
    GKfree(&lgraph.gdata, &lgraph.rdata, &lgraph.label, LTERM);
  }
}


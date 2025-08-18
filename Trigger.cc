int IsThreshold(int *tr_inp, double *base, double *vempeak, double threshold)
{
  int j,k;
  int threshx[3];

  for(k=0; k<3; k++)  threshx[k] = (int)(base[k]+threshold*vempeak[k]+.5);
  for(j=1; j<NSLOT; j++) {
    int nth[3];
    for(k=0; k<3; k++) nth[k] = (tr_inp[3*j+k]>threshx[k] || tr_inp[3*(j-1)+k]>threshx[k]);
    if(nth[0]+nth[1]+nth[2]>1) return j+1000*(nth[0]+2*nth[1]+3*nth[2]);
  }
  return 0;
}
/*-----------------------------------------------------------------------------------------------------*/
int IsToTd(int *tr_inp, double *base, double *vempeak, double threshold, int occupancy, int window)
{
  int j,k,n[3];
  int threshx[3],tdx[NSLOT][3];

  for(k=0; k<3; k++)  threshx[k] = (int)(base[k]+threshold*vempeak[k]+.5);

  // compute deconvoluted trace
  for(j=1; j<NSLOT; j++) {
    for(k=0; k<3; k++) {
      int Ai = tr_inp[3*j+k]*44;
      int Bi = tr_inp[3*j+k]*16-Ai/4;
      int Ci = Bi*52;
      int Di = Ci/16+8; 
      tdx[j][k] = Di/16;
    }
  }

  n[0] = n[1] = n[2] = 0;
  for(j=1; j<NSLOT; j++) {
    // update counts with the next slot
    for(k=0; k<3; k++)  n[k] += (tdx[j][k]>threshx[k]); 
    // trigger condition in this window: at least 2 PMTs with more than <occupancy> slots above <threshold>
    int p1 = n[0]>occupancy;
    int p2 = n[1]>occupancy; 
    int p3 = n[2]>occupancy; 
    if(p1+p2+p3>1) return j+1000*(p1+2*p2+3*p3);
    // remove from the count the slot before the window 
    if(j>window) {
      for(k=0; k<3; k++) n[k] -= (tdx[j-window][k]>threshx[k]);  // fixed point option 
    }
  }
  return 0;
}
/*-----------------------------------------------------------------------------------------------------*/
int IsMoPS(int *tr_inp, int low,int high,int window,int min,int nveto)
{
  int j,k,veto,lveto,jump,Jump,trace[3][NSLOT],Step[3][NSLOT],nStepWindow[3],p1,p2,p3;
  // position des "steps positifs"
  for(k=0; k<3; k++) {
    for(j=0; j<NSLOT; j++) {
      Step[k][j] = 0;
      trace[k][j] = tr_inp[3*j+k];
    }
    Jump = veto = 0;
    if(trace[k][1]>trace[k][0]) {jump = trace[k][1]-trace[k][0]; Jump = jump;}
    for(j=2; j<NSLOT; j++) {
      jump = trace[k][j]-trace[k][j-1];
      // fin de montee: on enregistre le step et on (re)definit nveto si necessaire
      if(jump<=0) {
        if(Jump>=low && veto<=0) Step[k][j-1] = Jump;

        if(nveto>=0 && Jump>=16) {
          lveto = 2;
          if(Jump>=32) lveto++;
          if(Jump>=64) lveto++;
          if(Jump>=128) lveto++;
          if(Jump>=256) lveto++;
          if(Jump>=512) lveto++;
          veto = max(veto,nveto+lveto);
        }
        Jump = 0;
      }
      else Jump += jump;
      if(veto>0) veto--;
    }
  }
  // condition de trigger sur fenetre glissante
  for(k=0; k<3; k++) nStepWindow[k] = 0;
  for(j=0; j<NSLOT; j++) {
    for(k=0; k<3; k++) {
     nStepWindow[k] += (Step[k][j]>=low && Step[k][j]<=high);
    }
    p1 = (nStepWindow[0]>=min);
    p2 = (nStepWindow[1]>=min);
    p3 = (nStepWindow[2]>=min);
    if(p1+p2+p3>1) return j+1000*(p1+2*p2+3*p3);
    if(j>=window)
      for(k=0; k<3; k++) nStepWindow[k] -= ((Step[k][j-window]>=low && Step[k][j-window]<=high));
  }
  return 0;
}
/*----------------------------------------------------------------------------------------*/

/*
 * Compile: Use Makefile
 * Version 5 changes pairing efficiency for most/all males. Changes in GetMate only
 */
#include "GD-FellerV5.h"

double ML(double pfcf){
  double x,val;
  x=log2(pfcf);
  if(x < -8) {val=1.;}
         else val=max(0.033,0.033-0.4835*(x+6));
  return val;
}

double GL(double pfcf){
  double x,val;
  x=log2(pfcf);
  if(x < -8) {val=0.;}
  else if(x < -6) {val=1.+ 0.5*(x+6.);}
	 else val=1.;
  return val;
}

double LE(double pfcf){
  double x;
  x=log2(pfcf);
  return max(0.025, (0.25269 + 0.031974 * x));
}

double FT(double pfcf){
  double x;
  x=log2(pfcf);
  return max(0.,0.127742*(x+8));
}

void GetPars(double T){
  
  ovi1=(0.03154*exp((T-4.7511)/10.580));                        /*oviposition*/
  ovic=2./(ovi1*FT(1.));                              /*egglaying coefficient */
  hatch=-0.0167105+0.03866*exp((T+7.26)/16.47906);              /*eclosion*/
  me= 0.01;                                                     /*egg mortality*/
  ma= 0.04;                                                     /* adult mortality */
  LoptT= 150;                              /* LoptT=BS*Lopt where BS=15 and Lopt=5 */
  pa=0.5787;                                                    /*emergence*/
  fem= 1./4300.;  /* female release fraction. carv15: 1F every 4300M */
  wipe=1;         /* No cleaning of breeding sites for treatment */
  u=1.;           /* Mean life of food*/

}

void GetRates(double dt){
  unsigned int i;
  unsigned int count;
  double nl,mor;
  double mortfac[SIZE*SIZE]; /*female mortality factors */

   mor=ma/(1-ML(pfcf));
   for(i=0;i<SIZE;i++){mortfac[i]=1.;};               /*wilds have no correction */
   for(i=SIZE;i<2*SIZE;i++){mortfac[i]=fm0LG;};
   for(i=2*SIZE;i<3*SIZE;i++){mortfac[i]=fm1LG;};
   for(i=3*SIZE;i<SIZE*SIZE;i++){mortfac[i]=fm2LG;};
   
   
   for(i=0;i<SIZE;i++){W[i]=me*((double)X[i]);};                      /* egg mortality */
   for(i=0;i<SIZE;i++){W[i+SIZE]=hatch*GL(pfcf)*((double)X[i]);};     /* egg hatching*/
   
   for(i=0;i<SIZE;i++){W[i+2*SIZE]=LE(pfcf)*ML(pfcf)*((double)X[i+SIZE]);};         /* larva mortality */
   for(i=0;i<SIZE;i++){W[i+3*SIZE]=LE(pfcf)*(1.-ML(pfcf))*((double)X[i+SIZE]);};    /* pupation */
   
   for(i=0;i<NOLET;i++){W[i+4*SIZE]=pa*ML(pfcf)*((double)X[i+2*SIZE]);};         /* pupa mortality*/
   for(i=NOLET;i<SIZE;i++){W[i+4*SIZE]=pa*(19./20.)*((double)X[i+2*SIZE]);};
   for(i=0;i<NOLET;i++){W[i+5*SIZE]=pa*(1.-ML(pfcf))*((double)X[i+2*SIZE]);};    /* pupa emergence */
   for(i=NOLET;i<SIZE;i++){W[i+5*SIZE]=(pa/20.)*((double)X[i+2*SIZE]);}; 

   for(i=0;i<NOLET+1;i++)W[i+6*SIZE]=mor*((double)X[i+3*SIZE]);                  /* male mortality  */   
   for(i=NOLET+1;i<SIZE;i++)W[i+6*SIZE]=ma*((double)X[i+3*SIZE]);                /* male mortality 2 lethal  */     
   for(i=0;i<SIZE*SIZE;i++)W[i+7*SIZE]=mor*mortfac[i]*((double)X[i+4*SIZE]);     /* female mortality */
   
   for(i=0;i<SIZE*SIZE;i++)W[i+SIZE*(7+SIZE)]=ovi1*X[i+4*SIZE];                  /* oviposition */

   return ;
 }

void UpdateFood(double dd){
  unsigned int i,larv;
  double LfCf;
   larv=0;                                                
   for(i=0;i<SIZE;i++){larv+=(X[i+SIZE]);};        
   LfCf=max(pfcf-1.,0.);   
   pfcf=LfCf*exp(-u*dd)+(LoptT+Xl)/((double)larv);
   Xl = Xl * exp(-u*dd);
   return ;
}

void GetMate(unsigned int female, unsigned int *class, double *Proph){
   unsigned int i, male, test;
   double countM,mc,maux[SIZE];
   double mm;

   test=0; mc=countM=0.;
/*
 * Put pairing-weighted males in an auxiliary array maux.
 * Only Wild males are not weighted by efi
 */
   maux[0]=(double)X[3*SIZE];
   maux[1]=efi0LG*(double)X[1+3*SIZE];
   maux[2]=efi1LG*(double)X[2+3*SIZE];
   maux[3]=efi2LG*(double)X[3+3*SIZE];
   
   for(i=0;i<SIZE;i++)countM+=maux[i];
/*
 * Random deviate
 */ 
   mm=countM*dsfmt_genrand_close_open(&dsfmt);
/*
 * Pick a male from class i
 */
   for(i=0; i<SIZE; i++){
       mc+=maux[i];
       if((mm <= mc)&&(test==0)){male=i; test=1;}  
       }     
   if(test==0){
       fprintf(stderr,"Error: No male chosen \n");
       exit(0);
     } /*found male */
   *class=SIZE*(4+female)+male; /* SIZE*4 is the female population offset */
   *Proph=(FProp[female]+Prop[male+3*SIZE])/2.;
   return ;
   }
 
void Oviposition(unsigned int ff){       /* ff index is 4*Female+Male */
             /* Egg indices: (0,Wild) (1,No lethal) (2,One lethal) (3,Two lethal) */
             /* Parental indices: Wild, 0, 1, 2 */
  unsigned int i,fi,legg,legg2,legg4;
  double countF;
  countF=0.; for(i=0;i<SIZE*SIZE;i++)countF+=(double)X[i+4*SIZE];
  eggs=(eggs*countF+ovic*FT(pfcf))/(countF+1.);
  legg=lround(eggs);
  legg2=lround(eggs/2.);
  legg4=lround(eggs/4.);  
  fi=4*SIZE+ff;               
  switch(ff){
    case 0 :                /* Female Wild Male Wild */
      X[0]+=legg; oo+=legg;
      break;
    case 1 :                /* Female Wild Male 0 */
      if(legg>0){
	Prop[1]=(legg*Prop[fi]+X[1]*Prop[1])/(double)(legg+X[1]);
        X[1]+=legg;
        oo+=legg;
      };
      break;
    case 2 :                /* Female Wild Male 1 */
      if(legg2>0){
	Prop[1]=(legg2*Prop[fi]+(X[1])*Prop[1])/(double)(legg2+X[1]);
        X[1]+=legg2;
        oo+=legg2;    /* Corrected from V1.0x */
        Prop[2]=(legg2*Prop[fi]+(X[2])*Prop[2])/(double)(legg2+X[2]);
        X[2]+=legg2;
        ofl+=legg2;
      };
      break;
    case 3 :                /* Female Wild Male 2 */
       if(legg>0){
	 Prop[2]=(legg*Prop[fi]+X[2]*Prop[2])/(double)(legg+X[2]);
         X[2]+=legg;
         ofl+=legg;
      };
      break;
    case 4 :                /* Female 0 Male Wild */
      if(legg>0){
	Prop[1]=(legg*Prop[fi]+X[1]*Prop[1])/(double)(legg+X[1]);
        X[1]+=legg;
        oo+=legg;
      };
      break;
    case 5 :                /* Female 0 Male 0 */
      if(legg>0){
	Prop[1]=(legg*Prop[fi]+X[1]*Prop[1])/(double)(legg+X[1]);
        X[1]+=legg;
        oo+=legg;
      };
      break;
    case 6 :                /* Female 0 Male 1 */
      if(legg2>0){
	Prop[1]=(legg2*Prop[fi]+(X[1])*Prop[1])/(double)(legg2+X[1]);
        X[1]+=legg2;
        oo+=legg2; 
        Prop[2]=(legg2*Prop[fi]+(X[2])*Prop[2])/(double)(legg2+X[2]);
        X[2]+=legg2;
        ofl+=legg2;
      };
    break;
    case 7 :                /* Female 0 Male 2 */
      if(legg>0){
	Prop[2]=(legg*Prop[fi]+X[2]*Prop[2])/(double)(legg+X[2]);
        X[2]+=legg;
        ofl+=legg;
      };
      break;
    case 8 :                /* Female 1 Male Wild */
      if(legg2>0){
	Prop[1]=(legg2*Prop[fi]+(X[1])*Prop[1])/(double)(legg2+X[1]);
        X[1]+=legg2;
        oo+=legg2;    /* Corrected from V1.0x */
        Prop[2]=(legg2*Prop[fi]+(X[2])*Prop[2])/(double)(legg2+X[2]);
        X[2]+=legg2;
        ofl+=legg2;
      };
      break;
    case 9 :                /* Female 1 Male 0 */
      if(legg2>0){
        Prop[1]=(legg2*Prop[fi]+(X[1])*Prop[1])/(double)(legg2+X[1]);
        X[1]+=legg2;
        oo+=legg2;
        Prop[2]=(legg2*Prop[fi]+(X[2])*Prop[2])/(double)(legg2+X[2]);
        X[2]+=legg2;
        ofl+=legg2;
      };
    break;
    case 10 :               /* Female 1 Male 1 */
      if(legg2>0){
	Prop[2]=(legg2*Prop[fi]+(X[2])*Prop[2])/(double)(legg2+X[2]);
        X[2]+=legg2;
        if(legg4>0){
	  Prop[1]=(legg4*Prop[fi]+(X[1])*Prop[1])/(double)(legg4+X[1]);
          X[1]+=legg4;
          oo+=legg4;
          Prop[3]=(legg4*Prop[fi]+(X[3])*Prop[3])/(double)(legg4+X[3]);
          X[3]+=legg4;
          ofl+=legg2+legg4;
	};
      };
    break;
    case 11 :               /* Female 1 Male 2 */
      if(legg2>0){
	Prop[2]=(legg2*Prop[fi]+(X[2])*Prop[2])/(double)(legg2+X[2]);
        X[2]+=legg2;
        Prop[3]=(legg2*Prop[fi]+(X[3])*Prop[3])/(double)(legg2+X[3]);
        X[3]+=legg2;
        ofl+=2*legg2;
      };
      break;
    case 12 :               /* Female 2 Male Wild */
      if(legg>0){
	Prop[2]=(legg*Prop[fi]+X[2]*Prop[2])/(double)(legg+X[2]);
        X[2]+=legg;
        ofl+=legg;
      };
    break;
    case 13 :               /* Female 2 Male 0 */
      if(legg>0){
	Prop[2]=(legg*Prop[fi]+X[2]*Prop[2])/(double)(legg+X[2]);
        X[2]+=legg;
        ofl+=legg;
      };
      break;
    case 14 :               /* Female 2 Male 1 */
      if(legg2>0){
	Prop[2]=(legg2*Prop[fi]+(X[2])*Prop[2])/(double)(legg2+X[2]);
        X[2]+=legg2;
        Prop[3]=(legg2*Prop[fi]+(X[3])*Prop[3])/(double)(legg2+X[3]);
        X[3]+=legg2;
        ofl+=2*legg2;
      };
    break;
    case 15 :               /* Female 2 Male 2 */
      if(legg>0){
	Prop[3]=(legg*Prop[fi]+X[3]*Prop[3])/(double)(legg+X[3]);
        X[3]+=legg;
        ofl+=legg;
      };
      break;
    default:
      fprintf(stderr, "Index error in oviposition\n"); exit(0);
  }
  return;
}

 void DoDelta(unsigned int j){
   unsigned int k,n,fs;
   unsigned int mr,fr;
   double tt,Proph=0.,dead;
   dead=FT(pfcf)*xl/FT(1.);
   if(j<SIZE){       /* egg mortality */
     X[j]--;
     }
   else if(j<2*SIZE){     /* egg hatching */
     X[j-SIZE]--;
     Prop[j]=(Prop[j-SIZE]+X[j]*Prop[j])/(double)(1+X[j]);
     X[j]++;
     }
   else if(j<3*SIZE){     /* larva mortality and recycling*/
     X[j-SIZE]--;
     Xl+=dead; 
      }
   else if(j<4*SIZE){     /* pupation */
     X[j-2*SIZE]--;
     Prop[j-SIZE]=(Prop[j-2*SIZE]+X[j-SIZE]*Prop[j-SIZE])/(double)(1+X[j-SIZE]);
     X[j-SIZE]++;
     }
   else if(j<5*SIZE){     /* pupa mortality and recycling*/
     X[j-2*SIZE]--;
     Xl+=dead; 
     }
   else if(j<6*SIZE){     /* pupa emergence and mating (if female emerges) */
     X[j-3*SIZE]--;
     tt=dsfmt_genrand_close_open(&dsfmt);
     if(tt<=0.5){   /*male emergence */
       Prop[j-2*SIZE]=(Prop[j-3*SIZE]+X[j-2*SIZE]*Prop[j-2*SIZE])/(double)(1+X[j-2*SIZE]);
       X[j-2*SIZE]++;
       }
     else{          /*female emergence */
       fs=0; for(k=0;k<SIZE;k++)fs+=X[k+SIZE*(4+j-5*SIZE)];
       FProp[j-5*SIZE]=(Prop[j-3*SIZE]+fs*FProp[j-5*SIZE])/(double)(1+fs);
       GetMate(j-5*SIZE,&n,&Proph); /* send female, get fecundated index and %ROCK */
       Prop[n]=(Proph+Prop[n]*X[n])/(double)(1+X[n]);
       X[n]++;
       };
     }
   else if(j<7*SIZE){     /* male mortality */
     X[j-3*SIZE]--;
     }
   else if(j<(7+SIZE)*SIZE){    /* female mortality */
     X[j-3*SIZE]--;
     }
   else if(j<(7+2*SIZE)*SIZE){ /*oviposition event */
     Oviposition(j-SIZE*(7+SIZE));
     }
   else {fprintf(stderr, "Wrong event index %d  \n",j); exit(0);};
  
   return ;
   }    

unsigned int BinS(double x){
  unsigned int i, first, last, middle;
  double ratesum[N], aux;
  aux=0;
  for (i=0;i<N;i++){
    aux+=W[i];
    ratesum[i]=aux;
  }
  first=0; last=N-1; 
  while (first < last) {
    middle=(first+last)/2;
    if (ratesum[middle] < x) first = middle+1;    
    else last=middle;
  };
  return first;   
}

int GetIndata(void )
{
  unsigned int j;
  int data, dum;
  FILE *infile;
  char mystring[1000];
  int offset,justread;
  data=0;
  justread=offset=0;
  
  if((infile=fopen("IC.dat","r")) != NULL){
  if(fgets(mystring,1000, infile)!=NULL){
      fclose(infile);
      sscanf(mystring+offset,"%d%n",&dum,&justread); offset+=justread;
      for(j=0;j<LIFE;j++){
	sscanf(mystring+offset,"%d%n",&X[j*SIZE],&justread); offset+=justread;
        fprintf(stderr,"%d ",X[j*SIZE]);
        }
      sscanf(mystring+offset,"%lf%lf%n",&Xl,&pfcf); 
      fprintf(stderr,"%lf %lf\n",Xl,pfcf);
      data=1;
      fprintf(stderr,"Read indata from IC.dat\n");
      }
    }  
  else fprintf(stderr, "No initial data IC\n");

  return data;
}

void SaveTransient(int dt, char fn[])
{
  unsigned int j;
  char name[20];
  FILE *ICfile, *OUTfile ;
  strcpy(name,"TR-");
  strcat(name,fn);
  if( ((ICfile=fopen("IC.dat","w")) != NULL)&& ((OUTfile=fopen(name,"w")) != NULL)){
     fprintf(ICfile, "%d ",(int)dt); 
     fprintf(OUTfile,"%d ",(int)dt); 
     for(j=0;j<LIFE;j++){
       fprintf(ICfile, "%d ",X[j*SIZE]);
       fprintf(OUTfile,"%d ",X[j*SIZE]);
       }
     fprintf(ICfile, "%lf %lf \n",Xl,pfcf); fflush(ICfile);
     fprintf(OUTfile,"%lf %lf \n",Xl,pfcf); fflush(OUTfile);
     fclose(ICfile);
     fclose(OUTfile);
     fprintf(stderr, "Transient written onto IC.dat and TR- file \n");
     }
 else fprintf(stderr, "Error writing transient\n");
}

double relfac(double xn, double xo){
  double fac,K,xp;
  K=exp(-ma*7.);
  xp=((xn-xo)*K+xn*(1.-xo))/((xn-xo)*K+(1.-xo));
  
if(((xp > Target) && (xn >= xo) ) || ((xp < Target) && (xn <= xo))){
  fac=(1.-xo)*((1.-xn)*(Target/(1.-Target))-K*xn)/(xn*(1.-xo)-K*xo*(1.-xn));
  }
else
  fac=1.;

  return min(max(fac,0.8),1.2);
 }

void DoRelease(double release){
  unsigned int mr, fr, fs,k,n;
  double Proph=0.;

  mr=(unsigned int)(release*(1.-fem)); /*males*/
  if(mr>0){
    Prop[4*SIZE-1]=(mr+Prop[4*SIZE-1]*X[4*SIZE-1])/(double)(mr+X[4*SIZE-1]);
    X[4*SIZE-1]+=mr;       
    relcount++;
    }
  else {
     fprintf(stderr,"Warning: Release due but release size is %d\n",mr); fflush(stderr);
     };
  fr=(unsigned int)(release*fem);
  if(fr>0){                         
    fs=0; for(k=0;k<SIZE;k++)fs+=X[k+7*SIZE]; /* sum over females 2 lethal */ 
    FProp[SIZE-1]=(fr+FProp[SIZE-1]*fs)/(double)(fr+fs);
    for(k=0;k<fr;k++){
	GetMate(3,&n,&Proph);          
        Prop[n]=(Proph+Prop[n]*X[n])/(double)(1+X[n]);
        X[n]++;
        }
       }
  };

void main(int argc, char *argv[])
{ 
  unsigned int n;
  int Weeks,Stat,ind;
  unsigned int i,j,test,testR, event;
  unsigned long eventcount;
  double dt, T, degs;
  double R,dd,rc,printT, Proph,aux,fac;
  char name[20],fn[20],fnr[20],number[4];
  FILE *file, *relf;

  /* Reading section */
  
if(argc<14){
 fprintf(stderr,"Usage: File1, Stat, Transient, Duration, After, xl, Target, seed, efi*3, fm*3\n");
 fprintf(stderr,"Note: File IC.dat with populations, Xl and Pf/Cf is also required. \n");
 exit(0);
}

sscanf(argv[1],"%s",&name);      fprintf(stderr,"File1: %s\n",name); fflush(stderr);
sscanf(argv[2],"%d",&Stat);      fprintf(stderr,"Stat=%i\n",Stat); fflush(stderr);
sscanf(argv[3],"%d",&Transient); fprintf(stderr,"Transient=%i\n",Transient); fflush(stderr);
sscanf(argv[4],"%d",&Duration);  fprintf(stderr,"Duration=%i\n",Duration); fflush(stderr);
sscanf(argv[5],"%d",&After);     fprintf(stderr,"After=%i\n",After); fflush(stderr);
sscanf(argv[6],"%lf",&xl);       fprintf(stderr,"xl=%lf\n",xl); fflush(stderr);
sscanf(argv[7],"%lf",&Target);   fprintf(stderr,"Target=%lf\n",Target); fflush(stderr);
sscanf(argv[8],"%d",&seed);      fprintf(stderr,"Seed for dsfmt =%i\n",seed); fflush(stderr);
sscanf(argv[9],"%lf",&efi0LG);   
sscanf(argv[10],"%lf",&efi1LG);  
sscanf(argv[11],"%lf",&efi2LG);  fprintf(stderr,"Male mating =%lf %lf %lf\n",efi0LG,efi1LG,efi2LG); fflush(stderr);
sscanf(argv[12],"%lf",&fm0LG);     
sscanf(argv[13],"%lf",&fm1LG);   
sscanf(argv[14],"%lf",&fm2LG);   fprintf(stderr,"Female death =%lf %lf %lf\n",fm0LG,fm1LG,fm2LG); fflush(stderr);

if(Stat>maxStat){
fprintf(stderr,"Stat: %d reset to %d \n",Stat,maxStat);
   Stat=maxStat;
 }
 
if(xl>MAX_XL || xl<0){
  xl=(double)(MAX_XL/2.);
fprintf(stderr,"Warning: xl reset to %lf \n",xl);
 }

Weeks=Transient+Duration+After;
if((Weeks>maxWeeks) || (Weeks<1)){
fprintf(stderr,"Improper no. of weeks: %d Maxweeks: %d Retry \n",Weeks,maxWeeks);
 exit(0);
 }
 fprintf(stderr,"Total nr of weeks: %d \n",Weeks);


/*Initialise temperature dependence*/
  degs=(double)Temp;
  GetPars(degs);

  /* initialize randon generator */
  dsfmt_init_gen_rand(&dsfmt,seed);
  
  /* initialize Target factor */
  if((Target>0.999999) || (Target<0.000001)){
    fprintf(stderr,"Target: %f reset to 1/2 \n",Target);
    Target=0.5;
    }  
  /* Repetitions loop */
for (i=0; i<Stat; i++){

    for(j=0;j<PoPS;j++)X[j]=0;
    for(j=0;j<PoPS;j++)Prop[j]=0.;
    for(j=0;j<SIZE;j++)FProp[j]=0.;
    Prop[4*SIZE-1]=1.;                         /*Released males are full Rockefeller */
    FProp[SIZE-1]=1.;                          /*Released females are full Rockefeller */
    ind=GetIndata();
    if(ind==0){
    fprintf(stderr,"Exit: No initial file(s) - Retry"); fflush(stderr);
    exit(0);
    }
  eggs=ovic*FT(pfcf);                /*initial eggs laid with input pfcf */

  for(j=0;j<N;j++){W[j]=0.;};

  dt=dd=printT=0.;
  testR=0;                                      /* For initial release of Oxitec */
  relcount=0;
  T=(double)(Weeks*7);

  eventcount=0;
  
  sprintf(number,"%d",i);
  strcpy(fn,name);
  strcat(fn,number);
  file=fopen(fn,"w");
  strcpy(fnr,"Rel");
  strcat(fnr,name);
  strcat(fnr,number);
  file=fopen(fn,"w");
  relf=fopen(fnr,"w");
  while (dt <= T){
    eventcount++;
    GetRates(dt);
    R=0.;
    for (j=0;j<N;j++){R+=(double)W[j];};
    dd=(double) (R*dsfmt_genrand_close_open(&dsfmt));
    
     event=BinS(dd); 

    dd= (-log((double)(dsfmt_genrand_close_open(&dsfmt))))/R;         /* Time for present event */
    UpdateFood(dd);
    DoDelta(event);
    dt+= dd;                                                          /* Elapsed time */

    printT+=dd;
    if(printT>7){
      printT-=7;
      fprintf(file,"%d ",(int)dt); 
      for(j=0;j<PoPS;j++){fprintf(file,"%d %f ",X[j],Prop[j]);};
      for(j=0;j<SIZE;j++)fprintf(file,"%f ",FProp[j]); 
      fprintf(file,"%f %f \n",Xl,pfcf); fflush(file);

      if((dt > (Transient*7)) && (dt< (Transient+Duration)*7)){
        if(testR==0){
          oo=ofl=0;
          testfl=testoldfl=0.;
          SaveTransient((int)dt,fn);
          X[0]=(unsigned int)(wipe*X[0]);
          X[SIZE]=(unsigned int)(wipe*X[SIZE]);
          X[2*SIZE]=(unsigned int)(wipe*X[2*SIZE]); /*Clean breeding sites before intervention*/
	  fac=1.;
	  aux=fac*(double)(11*X[3*SIZE]);               /*initial release 11 x adult_males */
          release=min(aux,(double)MAX_REL);
          testR=1;
          };
	if(oo+ofl>0){
          testfl=(double)ofl/((double)(oo+ofl));
	  fac=relfac(testfl,testoldfl);
          aux=fac*release;
	  release=min(aux,(double)MAX_REL);
	  testoldfl=testfl;
	  };
	
	DoRelease(release);
	fprintf(relf,"%d %d %f %f %d \n",(int)dt,(unsigned int)release,testfl/Target,fac,relcount); fflush(relf);
	oo=ofl=0;
	
        };
    };
  };
  fclose(file);	
  fclose(relf);	
  fprintf(stderr, "Events: %lu\n",eventcount);
 };
 exit(0);
}

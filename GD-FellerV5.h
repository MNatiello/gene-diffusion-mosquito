#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>

#include "../dSFMT-src-2.2.3/dSFMT.h"
uint32_t seed=17;   /*Should be unsigned int. No point in using -1 or any negative*/
dsfmt_t dsfmt;
static const int dsfmt_mexp = DSFMT_MEXP;

#define MAX_REL 100000  /*Maximal release*/
#define MAX_XL     1.0  /*Maximal xl*/
#define maxWeeks  1040  /*Max number of weeks */
#define maxStat	    50  /*Max number of repetitions */
#define LIFE         5  /*Wild life forms: Egg, Larva, Pupa, Males, Females */
#define LET          3  /*copies of lethal gene 0, 1 or 2*/
#define NOLET        2  /*No-lethal-gene-entries for each life form*/
#define SIZE         4  /*SIZE=1+LET*/
#define N           60  /*Total number of events N=7*SIZE+2*SIZE*SIZE */
#define PoPS        32  /*population size: SIZE*8*/
#define Temp        26  /*Celsius */
#define max(a,b)  ((a) > (b) ? (a) : (b))
#define min(a,b)  ((a) < (b) ? (a) : (b))

unsigned int Transient, Duration, After;     
double Target;                               /* Target */

/*Rate parameters*/
double hatch, me, pa, ovi1, ovic, ma, u, fem, xl;
unsigned int LoptT;
double eggs;
/*end rate parameters*/

/*Rate parameters discriminating genetics */
double efi0LG, efi1LG, efi2LG;  /*sexual mating efficiency non-wild males */
double fm0LG, fm1LG, fm2LG;     /*factor enhancing mortality female hybrids */
/*end */
      
unsigned int X[PoPS]; /* subpopulations X: LifeForm x (1+Lethal Gene)
                                             8      x (1+          3)
                           life-form: 
                            egg              [0-3]
                            larva            [4-7]
                            pupa             [8-11]
                            adult male       [12-15]
                            adult female     [16-31] adultfemales are already mated
                            SIZE is 1 wild lineage and 0,1,2 lethal Gene */
double Prop[PoPS];  /* Proportion of Rockefeller gene on each class */
double FProp[SIZE]; /* Total adult female proportion by class */
double W[N];        /* Event-rates */

double Xl,pfcf;                    /*global variables for food */
unsigned int relcount;             /*global variables for release */
double release,testfl,testoldfl,wipe;
unsigned int oo,ofl;

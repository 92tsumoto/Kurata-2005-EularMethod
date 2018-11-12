//#ifndef __SYSPARA_H_INCLUDE 
//#define __SYSPARA_H_INCLUDE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "mkl.h"
#include "./include/xhplot.h"

#define NN 16
#define BUF 100
#define NUM 60

//#define R 8314.472		// J/kmol/K
//#define F 96485.33771638995	// C/mol
//#define T 310.0		// K
#define R 8314400000000000		// J/kmol/K
#define F 96485000000000000	// C/mol
#define T 310.15		// K

#define VNMAX 400*5+1
#define dvm 5

struct varstruct {

    int datas;
    int line_wid[NUM];
	
	int n;

	// An invariant constant
	double FRT,RTF;

	// Cell tupe
	int celltype;

	// Cell Geometry
	double vcell;
	double vi,vup,vrel;
	double CapHV;

	// Ion Valences 
	double zna,zk,zca;

	// Reversal potential
	double Ena,Ek,Eks,Eca;
	double prnak;
			
	// Total Ion currents 
	double Itotal;

	// Net Ion Fluxes
	double jcanet,jnanet,jknet;

	// Troponin Ca Buffering (in Myoplasm)	
	double concTc,rftn,kftc,kbtc;

	// Calmodulin Ca Buffering (in Myoplasm) --- Rapid Buffering Approximation	
	double conccm,kdcm,bcm;

	// Calsequestrin Ca Buffering (in SR) --- Rapid Buffering Approximation
	double conccq,kdcq,bcq;

	// Extracellular ion concentrations
	double nao,ko,cao;

	// Intracellular ion concentrations
	double ki;

	// Base Currnt Stimulus
    double Istim,Istim_base;
	double omega;

	// test variable
	double dt;
	double gkr_rate,gks_rate;
	// Sttimulus parameters
	double BCL;  // Base cycle length = stimulus period
	int beat; // Number of stimulus

    int m;
    int l;

    double x0[NUM][NN];
    double tsign[NUM];
    double tend[NUM];

    int pflag;
    int write, graph;
    int write0;
    int half;

} var;

// Fast and Late sodium currnets
struct inastruct {

	double Gna,fast,pkna,Enaf;
	double am,bm,minf,taum,hinf,tauh;
	double *Tam,*Tbm,*Thinf,*Ttauh;

} ina;

// Transient Outward Current (Ito)
struct itostruct {

	double ik,Gto,Ek;
	double rinf,taur,qinf,tauq;
	double ar,br,aq1,bq1,aq2,bq2;
	double *Tar,*Tbr,*Taq1,*Tbq1,*Taq2,*Tbq2;
	double pnato,dVgq;

} ito;

// L-type Calcium channel current (IcaL)
struct icalstruct {

	double dinf,taud,finf,tauf;
	double *Tdinf,*Ttaud,*Tfinf,*Ttauf;
	double kmfca,gfCa_inf;
	double gca,Eca;
	double ica,icana,icak;
	double tmp;
} ical;

// Rapid activating potassium current (Ikr)
struct ikrstruct {

	double ik,Gkr;
	double pinf,ap,bp,taup,piinf;
	double *Tpinf,*Tap,*Tbp,*Tpiinf;

} ikr;

// Slowlactivating potassium current (Iks)
struct iksstruct {

	double ik,Gks,Ek;
	double ninf,taun;
	double *Tninf,*Ttaun;
		
} iks;

// Inward rectifier potassium current (Ik1)
struct ik1struct {

	double ik,Gk1;
	double k1inf,ak1,bk1;
	double *Tak1,*Tbk1;

} ik1;

struct ncxstruct {

	double kncx;
	double kmnaex,kmcaex,rncx,ksat;
	double jncx;
	double *Texp0,*Texp1;
	double exp0,exp1,c1,c2,c3;

} ncx;

// Sodium-Potassium Pump
struct inakstruct {

	double kmnap,kmkp,nna;
	double G,inak;
	double rhonak,*Trhonak;
	double *Texp0,exp0,*Texp1,exp1;

} inak;

// Sarcolemmal Ca Pump
struct ipcastruct {

	double G,km,ca;

} ipca;

// Na Background Current
struct inabstruct {

	double pnab,na;


} inab;

// K Background Current
struct ikbstruct {

	double G,k;
	double xkb,*Txkb;

} ikb;

// Ca Background Current
struct icabstruct {

	double pcab,gacai,gacao,ca;
	double exp,*Texp;

} icab;

// SR calcium release flux, via RyR (Jrel)
struct jrelstruct {

	double drss,tau_dr;
	double dfss,tau_df;
	double prel,nrel;
	double ca;

} jrel;

// Calcium uptake via SERCA pump
struct jupstruct {

	double pup,kup;
	double ca;

} jup;

// diffusion flux
struct jleakstruct {

	double pleak,ca;

} jleak;

// Translocation of Ca Ions from NSR to JSR
struct jtrstruct {

	double tau,ca;

} jtr;

void val_consts(FILE *);
void make_ExPTable();

//void eular(int n,double h,double x[],double t);
void runge(int n,double h,double x[],double t);
void function(double x[],double f[],double t);
void input_para(FILE *);
void initial_mem();
void closed_mem();

void eventloop(FILE *, int *mode, int *P, double m[]);
void orbit(int *mode, double m[], double x2);
void draw_p(int *mode, int P, double x[], double x2);
void mouse(int *mode, double x[], double x2);

void comp_reversal_potential(double x[]);
void comp_ina(double x[]);
void comp_ito(double x[]);
void comp_ical(double x[]);
void comp_ikr(double x[]);
void comp_iks(double x[]);
void comp_ik1(double x[]);
void comp_inaca(double x[]);
void comp_inak(double x[]);
void comp_ipca(double x[]);
void comp_ikb(double x[]);
void comp_icab(double x[]);
void comp_inab(double x[]);
void comp_cicr(double x[]);
void comp_jup(double x[]);
void comp_jleak(double x[]);
void comp_jtr (double x[]);
void comp_concentration (double x[]);

main(int argc,char **argv);

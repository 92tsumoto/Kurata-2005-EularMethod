#include "syspara.h"

void comp_ina(double x[])
{
	//MKL_INT iV=0;
	int iV=0;
	double V1,V2,d1,d2;
	
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;
	//printf("iV=%d,V1=%f,V2=%f,d1=%f,d2=%f\n",iV,V1,V2,d1,d2);
	
	ina.am = ina.Tam[iV]*d2 + ina.Tam[iV+1]*d1;
	ina.bm = ina.Tbm[iV]*d2 + ina.Tbm[iV+1]*d1;
	ina.minf = ina.am/(ina.am+ina.bm);
	ina.taum = 1.0/(ina.am+ina.bm);
	ina.hinf = ina.Thinf[iV]*d2 + ina.Thinf[iV+1]*d1;
	ina.tauh = ina.Ttauh[iV]*d2 + ina.Ttauh[iV+1]*d1;
	
	ina.fast = ina.Gna*(x[0]-ina.Enaf)*x[1]*x[1]*x[1]*x[2]*x[2];

}

// Ito Transient Outward Current
void comp_ito (double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ito.ar = ito.Tar[iV]*d2 + ito.Tar[iV+1]*d1;
	ito.br = ito.Tbr[iV]*d2 + ito.Tbr[iV+1]*d1;
	ito.rinf = ito.ar/(ito.ar+ito.br);
	ito.taur = 1.0/(ito.ar+ito.br);

	ito.aq1 = ito.Taq1[iV]*d2 + ito.Taq1[iV+1]*d1;
	ito.bq1 = ito.Tbq1[iV]*d2 + ito.Tbq1[iV+1]*d1;
	ito.aq2 = ito.Taq2[iV]*d2 + ito.Taq2[iV+1]*d1;
	ito.bq2 = ito.Tbq2[iV]*d2 + ito.Tbq2[iV+1]*d1;

	ito.qinf = ito.aq1/(ito.aq1+ito.bq1);
	ito.tauq = 1.0/(ito.aq2+ito.bq2);

	ito.ik = ito.Gto*(x[0]-ito.Ek)*x[5]*x[6];

}


// L-type calcium current
void comp_ical(double x[])
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

// VDA
	ical.dinf = ical.Tdinf[iV]*d2 + ical.Tdinf[iV+1]*d1;
	ical.taud = ical.Ttaud[iV]*d2 + ical.Ttaud[iV+1]*d1;
// VDI 
	ical.finf = ical.Tfinf[iV]*d2 + ical.Tfinf[iV+1]*d1;
	ical.tauf = ical.Ttauf[iV]*d2 + ical.Ttauf[iV+1]*d1;
// CDI 
	ical.gfCa_inf = ical.kmfca/(ical.kmfca+x[12]);

	ical.ica =ical.gca*(x[0]-ical.Eca)*x[3]*x[4]*ical.gfCa_inf;

}

// Rapidly Activating Potassium Current 
void comp_ikr (double x[])
{
	MKL_INT iV=0;	
	double V1,V2,d1,d2;
	
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ikr.pinf = ikr.Tpinf[iV]*d2 + ikr.Tpinf[iV+1]*d1;
	ikr.ap = ikr.Tap[iV]*d2 + ikr.Tap[iV+1]*d1;
	ikr.bp = ikr.Tbp[iV]*d2 + ikr.Tbp[iV+1]*d1;
	ikr.taup = 1.0/(ikr.ap+ikr.bp);

	ikr.piinf = ikr.Tpiinf[iV]*d2 + ikr.Tpiinf[iV+1]*d1;

	//printf("pinf=%f,ap=%f,bp=%f,taup=%f,piinf=%f\n",ikr.pinf,ikr.ap,ikr.bp,ikr.taup,ikr.piinf);
	ikr.ik = var.gkr_rate*ikr.Gkr*x[7]*ikr.piinf*(x[0]-var.Ek);

}

// Slowly Activating Potassium Current 
void comp_iks (double x[])
{
	
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	iks.ninf = iks.Tninf[iV]*d2 + iks.Tninf[iV+1]*d1;
	iks.taun = iks.Ttaun[iV]*d2 + iks.Ttaun[iV+1]*d1;
	
	iks.ik = var.gks_rate*iks.Gks*(x[0]-iks.Ek)*x[8]*x[8];

}

// Inward rectifier potassium current (Ik1)
void comp_ik1 (double x[])
{
        
	MKL_INT iV=0;   
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ik1.ak1 = ik1.Tak1[iV]*d2 + ik1.Tak1[iV+1]*d1;
	ik1.bk1 = ik1.Tbk1[iV]*d2 + ik1.Tbk1[iV+1]*d1;
	ik1.k1inf = ik1.ak1/(ik1.ak1+ik1.bk1);

	ik1.ik = ik1.Gk1*ik1.k1inf*(x[0]-var.Ek);

}

// Sodium-Calcium Exchanger V-S

void comp_inaca (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ncx.exp0=ncx.Texp0[iV]*d2 + ncx.Texp0[iV+1]*d1;
	ncx.exp1=ncx.Texp1[iV]*d2 + ncx.Texp1[iV+1]*d1;
	ncx.c2 = var.cao*x[15]*x[15]*x[15];
	ncx.c3 = var.nao*var.nao*var.nao*x[12];

	ncx.jncx=ncx.c1*(ncx.c2*ncx.exp0-ncx.c3*ncx.exp1)/(1.0+ncx.ksat*ncx.exp1);

}

// Sodium-Potassium Pump

void comp_inak (double x[])
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	inak.exp0 = inak.Texp0[iV]*d2 + inak.Texp0[iV+1]*d1;
	inak.exp1 = inak.Texp1[iV]*d2 + inak.Texp1[iV+1]*d1;

	inak.inak = inak.G*(var.ko/(var.ko+inak.kmkp))/(1.0+pow((inak.kmnap/x[15]),inak.nna))/(1.0+0.1245*inak.exp0+0.0365*inak.rhonak*inak.exp1);

}

// Sarcolemmal Ca Pump 

void comp_ipca (double x[])
{

	ipca.ca = ipca.G*x[12]/(ipca.km + x[12]);

}

// Ca Background Current 

void comp_icab (double x[])
{
	
	icab.ca = icab.pcab*(x[0]-var.Eca);

}

// Na Background Current 

void comp_inab (double x[])
{

	inab.na = inab.pnab*(x[0]-var.Ena);

}

void comp_cicr (double x[])
{

	double ica_total;
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (x[0]+200)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	jrel.drss = ical.Tdinf[iV]*d2 + ical.Tdinf[iV+1]*d1;
	jrel.dfss = ical.Tfinf[iV]*d2 + ical.Tfinf[iV+1]*d1;

	ica_total = ical.ica + icab.ca - 2.0*ncx.jncx + ipca.ca;
	jrel.ca = jrel.prel/(1.0+exp((ica_total+5.0)/0.9))*pow((x[9]*x[10]),jrel.nrel)*(x[13]-x[12]);

}

void comp_jup (double x[])
{

	jup.ca = jup.pup*x[12]*x[12]/(x[12]*x[12] + jup.kup*jup.kup);
}

void comp_jleak (double x[])
{
	jleak.ca = jleak.pleak*(x[14]-x[12]);

}
void comp_jtr (double x[])
{

	jtr.ca = (x[14]-x[13])/jtr.tau;

}

void comp_concentration (double x[])
{
	var.bcm=1.0/(1.0+var.conccm*var.kdcm/((var.kdcm+x[12])*(var.kdcm+x[12])));         // Dfcm=kfCM*Cai*(1-fcm)-kbCM*fcm;
	var.bcq = 1.0/(1.0+var.conccq*var.kdcq/((var.kdcq+x[13])*(var.kdcq+x[13])));	// Dfcq=kfCQ*Carel*(1-fcq)-kbCQ*fcq;

}


// Reversal potentials */

void comp_reversal_potential(double x[])
{
	var.Ena = (var.RTF/var.zna)*log(var.nao/x[15]);
	ina.Enaf = (var.RTF/var.zna)*log((var.nao+ina.pkna*var.ko)/(x[15]+ina.pkna*var.ki));
	iks.Ek = var.RTF*log((var.prnak*var.nao+var.ko)/(var.prnak*x[15]+var.ki));
	ito.Ek = var.RTF*log((ito.pnato*var.nao+var.ko)/(ito.pnato*x[15]+var.ki));
	var.Eca = (var.RTF/var.zca)*log(var.cao/x[12]);

	//printf("Ena=%lf, Ek=%lf, Eks=%lf\n",var.Ena,var.Ek,var.Eks);
}


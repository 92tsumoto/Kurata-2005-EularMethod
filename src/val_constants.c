
#include "syspara.h"

void val_consts(FILE *fp1)
{
	int i,w;
	double v_old,dvdt,dvdt_new;

	// frequency
		var.omega = 2.0*M_PI/var.BCL;

	// Cell Geometry */
		var.CapHV = 153.4;  // Membrane Capacitance (pF)
		var.vcell = 38.0*0.000000000001; // Cell Volume:3.801e-5 (uL)
		var.vi = var.vcell*0.68;      // Myoplasm volume (uL) = 68% for Cell volume
		var.vup = var.vcell*0.0552;    // NSR volume (uL)
		var.vrel = var.vcell*0.0048;    // JSR volume (uL)


	// Ion Valences
		var.zna = 1.0;  // Na valence
		var.zk = 1.0;   // K valence
		var.zca = 2.0;  // Ca valence

	// invariant constant
		var.FRT = F/(R*T);
		var.RTF = R*T/F;

	// Extracellular Concentrations
		var.nao = 140.0;    // (mM) 
		var.ko = 5.4;      // (mM)
		var.cao = 2.0;     // (mM)

	// Intracellular Concentrations
		//var.nai = 10.0;    // (mM) 
		var.ki = 140.0;      // (mM)
		//var.cai = 0.0001;     // [Ca] in myoplasm (mM)
		//var.carel = 1.0;     // [Ca] in SR Ca release site (junctional SR) (mM)
		//var.caup = 1.0;     // [Ca] in SR Ca uptake site (network SR) (mM)

		var.prnak = 0.01833;     // Na/K Permiability Ratio
		ina.pkna = 0.12;     // Na/K Permiability Ratio

	// Equilibrium Potentials (mV)
		var.Ek = (var.RTF/var.zk)*log(var.ko/var.ki);

	// Fast sodium current
		ina.Gna = 7.8;	// (mS/uF).
		
	// Transient outward current
		ito.pnato = 0.043;

		if(var.celltype==0){
			ito.dVgq= -12.0; //gating variables (q/r) % Endo (from Bernus 2002 model)
			ito.Gto = 0.13;	// Endo (from Bernus 2002 model).
		} else if(var.celltype==1){
			ito.dVgq= 0.0; //gating variables (q/r) % Epi (default for Kurata 2005 model)
			ito.Gto = 0.4;	// Epi (default for Kurata 2005 model).
		} else if(var.celltype==2){
			ito.dVgq= -4.0; //gating variables (q/r) % Mid (from Bernus 2002 model)
			ito.Gto = 0.35;	// Mid (from Bernus 2002 model).
		}

	// L-type calcium current
		ical.kmfca = 0.00035;
		ical.Eca = 52.8; // (mV)

		if(var.celltype == 0){ // Endo (based on O'Hara-Rudy 2011 model)
			ical.gca = 0.064*3.9*(1/1.2);
		} else if (var.celltype == 1){ // Epi (Default for Kurata 2005 model)
			ical.gca = 0.064*3.9;      	
		} else if (var.celltype == 2){ // Mid (based on O'Hara-Rudy 2011 model)
			ical.gca = 0.064*3.9*(2.5/1.2);
		}

	// Rapid delayed rectifier potassium current (Ikr)
		if(var.celltype == 0){
			ikr.Gkr = 0.012*(1.0/1.3);  //Endo (based on O'Hara-Rudy 2011 model)
		} else if(var.celltype == 1){
			ikr.Gkr = 0.012;  		// Epi (default for Kurata 2005 model)
		} else if(var.celltype == 2){
			ikr.Gkr = 0.012*(0.8/1.3);  // Mid (based on O'Hara-Rudy 2011 model)
		}

	// Slow delayed rectifier potassium current (Iks)
		if(var.celltype == 0){
			iks.Gks = 0.024*1.5*(1.0/1.4);	// Endo (based on O'Hara-Rudy 2011 model)
		} else if(var.celltype == 1){
			iks.Gks = 0.024*1.5;			// Epi (default for Kurata 2005 model)
		} else if(var.celltype == 2){
			iks.Gks = 0.024*1.5*(1.0/1.4);		// Mid (based on O'Hara-Rudy 2011 model)
		}
		
		//iks.Ek = var.RTF*log((var.prnak*var.nao+var.ko)/(var.prnak*var.nai+var.ki));
	
	// Inward rectifier K current: Ik1
		if(var.celltype == 0){
			ik1.Gk1 = 3.9*(1.0/1.2);		// Endo (based on O'Hara-Rudy 2011 model)
		} else if(var.celltype == 1){
			ik1.Gk1 = 3.9;					// Epi (default for Kurata 2005 model)
		} else if(var.celltype == 2){
			ik1.Gk1 = 3.9*(1.3/1.2);		// Mid (based on O'Hara-Rudy 2011 model)
		}

	// Sodium-Calcium Exchanger V-S
		ncx.kmnaex = 87.5;
		ncx.kmcaex = 1.38;
		ncx.rncx = 0.35;
		ncx.ksat = 0.1;

		if(var.celltype == 0){
			ncx.kncx = 1000.0*(1.0/1.1);	// Endo (based on O'Hara-Rudy 2011 model);
		} else {
			ncx.kncx = 1000.0;				// Epi (default for Kurata 2005 model) or Mid;
		}
		ncx.c1 = ncx.kncx/(pow(ncx.kmnaex,3)+pow(var.nao,3))/(ncx.kmcaex+var.cao);

	// Sodium-Potassium Pump
		inak.kmnap = 10.0;
		inak.kmkp = 1.5;
		inak.nna = 1.5;

		if(var.celltype == 0){
			inak.G = 1.3*0.68*(1.0/0.9);	// Endo (based on O'Hara-Rudy 2011 model)
		} else if(var.celltype == 1){
			inak.G = 1.3*0.68;				// Epi (default for Kurata 2005 model)
		} else if(var.celltype == 2){
			inak.G = 1.3*0.68;				// Mid (based on O'Hara-Rudy 2011 model)
		}
		inak.rhonak = (exp(var.nao/67.3)-1.0)/7.0;

	// Sarcolemmal Ca Pump
		ipca.G = 0.275*0.4;		// Max. Ca current through sarcolemmal Ca pump (mS/uF)
		ipca.km = 0.0005;		// Half-saturation concentration of sarcolemmal Ca pump (mM)

	// K Background Current 
		if(var.celltype != 1){
			ikb.G = 0.003;		// Max. conductance of K background (mS/uF)
		} else {
			ikb.G = 0.6*0.003;		// Max. conductance of K background (mS/uF)
		}

	// Ca Background Current 
		icab.pcab = 0.00085*0.6;	// (cm/s)

	// Na Background Current 
		inab.pnab = 0.001;

	// SR calcium release flux, via RyR (Jrel)
		jrel.nrel = 3.0;
		if(var.celltype!=2){
			jrel.prel = 50.0;	// Epi (default for Kurata 2005 model) or Endo
		} else {
			jrel.prel = 50.0*1.7;	// Mid (based on O'Hara-Rudy 2011 model)
		}
		jrel.tau_dr = 1.0*4.0;
		jrel.tau_df = 4.0*4.0;

	// calcium uptake via SERCA pump (Jup)
		jup.kup = 0.00025;
		if(var.celltype != 1){
			jup.pup = 0.000425*5.2*(1.0/1.3);	// Mid or Endo (based on O'Hara-Rudy 2011 model)
		} else {
			jup.pup = 0.000425*5.2;				// Epi (default for Kurata 2005 model)
		}

	// Translocation of Ca Ions from NSR to JSR
		jtr.tau = 180.0;      // Time constant of Ca transfer from NSR to JSR (ms)

	// Ca leakage from SR to Cytoplasm
		jleak.pleak = 0.00026;

	//Troponin Ca Buffering (in Myoplasm)
		var.concTc = 0.07;
		var.rftn = 10.0;
		var.kftc = 40.0*var.rftn;
		var.kbtc = 0.02*var.rftn;
		
	//Calmodulin Ca Buffering (in Myoplasm) --- Rapid Buffering Approximation
		if(var.celltype != 1){
			var.conccm = 0.05*(1.0/1.3);  //% Mid or Endo (based on O'Hara-Rudy 2011 model)
		} else {
			var.conccm = 0.05;   //% Epi (default for Kurata 2005 model)
		}
		var.kdcm = 0.00238;	// % Kd for Ca-binding to Calmoduline 

	//Calsequestrin Ca Buffering (in SR) --- Rapid Buffering Approximation	
		var.conccq = 10.0;     // Max. [Ca] buffered in CSQN (mM)
		var.kdcq = 0.8;     // % Kd for Ca-binding to Calsequestrin

}


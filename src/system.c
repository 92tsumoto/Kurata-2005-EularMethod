#include "syspara.h"

void function(double x[],double f[],double t)
{

	comp_reversal_potential(x);
	comp_ina(x);
	comp_ito(x);
	comp_ical(x);
	comp_ikr(x);
	comp_iks(x);
	comp_ik1(x);
	comp_inaca(x);
	comp_inak(x);
	comp_ipca(x);
	comp_icab(x);
	comp_inab(x);
	
	// Total Membrane Current (pA/pF)
	var.Itotal = ina.fast + ical.ica + ito.ik + ikr.ik + iks.ik + ik1.ik 
					+ inab.na + icab.ca + inak.inak + ncx.jncx + ipca.ca;

	// Net Ion Fluxes (umol/cell)
	var.jcanet = (ical.ica + icab.ca - 2.0*ncx.jncx + ipca.ca)*var.CapHV/(var.zca*F);
	var.jnanet = (ina.fast + inab.na + 3.0*inak.inak + 3.0*ncx.jncx)*var.CapHV/(var.zna*F);
	var.jknet = (ito.ik + ikr.ik + iks.ik + ik1.ik - 2.0*inak.inak)*var.CapHV/(var.zk*F);

	comp_cicr(x);
	comp_jup(x);
	comp_jleak(x);
	comp_jtr(x);
	comp_concentration(x);

	// Vm
	f[0] = (var.Istim-var.Itotal)/1.0;
	// Fast sodium current
	f[1] = ((ina.minf - x[1])/ina.taum)/1.0; // m
	f[2] = ((ina.hinf - x[2])/ina.tauh)/1.0; // h
	// L-Type calcium current
	f[3] = ((ical.dinf - x[3])/ical.taud)/1.0; // d
	f[4] = ((ical.finf - x[4])/ical.tauf)/1.0; // f 
	// Ito
	f[5] = ((ito.rinf - x[5])/ito.taur)/1.0; // r
	f[6] = ((ito.qinf - x[6])/ito.tauq)/1.0; // q
	// Ikr
	f[7] = ((ikr.pinf - x[7])/ikr.taup)/1.0; // p
	// IKs
	f[8] = ((iks.ninf - x[8])/iks.taun)/1.0; // n 
	// CICR 
	f[9] = ((jrel.drss - x[9])/jrel.tau_dr)/1.0; // dr
	f[10] = ((jrel.dfss - x[10])/jrel.tau_df)/1.0; // df
	// Troponin Ca Buffering (in Myoplasm)
	f[11] = (var.kftc*x[12]*(1.0-x[11])-var.kbtc*x[11])/1.0; // ftc
	// [Ca]i
	f[12] = (var.bcm*((-var.jcanet+jrel.ca*var.vrel - jup.ca*var.vup + jleak.ca*var.vup)/var.vi-var.concTc*(var.kftc*x[12]*(1.0-x[11])-var.kbtc*x[11])))/1.0;
	// [Ca]_rel
	f[13] = (var.bcq*(jtr.ca - jrel.ca))/1.0;
	// [Ca]_up
	f[14] = (jup.ca - jtr.ca*var.vrel/var.vup - jleak.ca)/1.0;
	// [Na]i
	f[15] = (-var.jnanet/var.vi)/1.0;
	// [K]i
	//f[16] = -var.jknet/var.vi;

}

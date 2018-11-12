#include "syspara.h"

void make_ExpTable()
{

	int vindex,kiindex;
	double v,ki;
    
	for(vindex=0;vindex<VNMAX;vindex++){

        v = (double)vindex/dvm-200.0;
		
        /** for ina **/
		ina.Tam[vindex] = 0.32*(v+47.13)/(1.0-exp(-(v+47.13)/10.0));
		ina.Tbm[vindex] = 0.08*exp(-v/11.0);

		ina.Thinf[vindex] = 0.5*(1.0-tanh(7.74+0.12*v));
		if(v<180){
			ina.Ttauh[vindex] = 0.25+2.24*(1.0-tanh(7.74+0.12*v))/(1.0-tanh(0.07*(v+92.4)));
		} else {
			ina.Ttauh[vindex] = 0.25;
		}

		// ito
		ito.Tar[vindex] = 0.5266*exp(-0.0166*(v-42.2912))/(1.0+exp(-0.0943*(v-42.2912)));
		ito.Tbr[vindex] = (0.5149*exp(-0.1344*(v-5.0027))+0.00005186*v)/(1.0+exp(-0.1348*(v-0.00005186)));

		ito.Taq1[vindex] = (0.0721*exp(-0.173*(v+34.2531-ito.dVgq))+0.00005612*(v-ito.dVgq))/(1.0+exp(-0.1732*(v+34.2531-ito.dVgq)));
		ito.Tbq1[vindex] = (0.0767*exp(-0.00000000166*(v+34.0235-ito.dVgq))+0.0001215*(v-ito.dVgq))/(1.0+exp(-0.1604*(v+34.0235-ito.dVgq)));

		ito.Taq2[vindex] = (0.0721*exp(-0.173*(v+34.2531))+0.00005612*v)/(1.0+exp(-0.1732*(v+34.2531)));
		ito.Tbq2[vindex] = (0.0767*exp(-0.00000000166*(v+34.0235))+0.0001215*v)/(1.0+exp(-0.1604*(v+34.0235)));

		// for ical
		ical.Tdinf[vindex] = 1.0/(1.0+exp(-(v+7.64)/6.32));
		ical.Ttaud[vindex] = (1.4/(1.0+exp(-(v+35.0)/13.0))+0.25)*1.4/(1.0+exp((v+5.0)/5.0))+1.0/(1.0+exp(-(v-50.0)/20.0));

		ical.Tfinf[vindex] = 1.0/(1.0+exp((v+24.6)/6.9));
		ical.Ttauf[vindex] = 23.9*0.75/(0.1389*exp(-(0.0358*(v-10.9))*(0.0358*(v-10.9)))+0.0519);

		// for ikr 
		ikr.Tpinf[vindex] = 1.0/(1.0+exp(-(v+14.0)/7.7));
		if(v==-14.0){
			ikr.Tap[vindex] = 0.0;
		} else {
			ikr.Tap[vindex] = 0.0003*(v+14.0)/(1.0-exp(-(v+14.0)/5.0));
		}
		ikr.Tbp[vindex] = 0.000073898*(v-3.4328)/(exp((v-3.4328)/5.1237)-1.0);
		ikr.Tpiinf[vindex] = 1.0/(1.0+exp((v+15.0)/22.4));

		// for iks 
		iks.Tninf[vindex] = 1.0/sqrt(1.0+exp(-(v-9.4)/11.8));
		iks.Ttaun[vindex] = 555.0/(1.0+exp(-(v+22.0)/11.3))+129.0;

		// ik1 
		ik1.Tak1[vindex] = 0.1/(1.0+exp(0.06*(v-var.Ek-200.0)));
		ik1.Tbk1[vindex] = (3.0*exp(0.0002*(v-var.Ek+100.0))+exp(0.1*(v-var.Ek-10.0)))/(1.0+exp(-0.5*(v-var.Ek)));

		// inaca
		ncx.Texp0[vindex] = exp(ncx.rncx*var.FRT*v);
		ncx.Texp1[vindex] = exp((ncx.rncx-1.0)*var.FRT*v);

		// inak 
		inak.Texp0[vindex] = exp(-0.1*var.FRT*v);
		inak.Texp1[vindex] = exp(-1.0*var.FRT*v);

	}
//exit(0);
}

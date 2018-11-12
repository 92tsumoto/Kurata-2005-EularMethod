/* produced by Tsumoto. K 2008.10.27 */

#include <string.h>
#include <stdlib.h>
#include "syspara.h"

FILE *fopen(), *fpin, *fp0, *fp1, *fp2, *fp3;
int mode = 1;
int P = 2;
int beats = 30;

typedef double Number;

main(argc,argv)
int argc;
char **argv;
{
	int i,w;
	int ii=0;
	double x[NN];
	double t = 0.0;
	double time=0.0;
	double h;
	double v_old,dvdt,dvdt_new;
	double t_stok;
	char *tmpname;
	char cmd[BUFSIZ];
	double tend;

/* Action Potential Duration and Max. Info */
//	double *vmax ; // Max. Voltage (mV)
//	double *dvdtmax ; // Max. dv/dt (mV/ms)
//	double *apd; // Action Potential Duration
//	double *toneapd; // Time of dv/dt Max.
//	double *ttwoapd; // Time of 90% Repolarization
//	double *rmbp; // Resting Membrane Potential
//	double *nair; // Intracellular Na At Rest
//	double *cair; // Intracellular Ca At Rest
//	double *kir ; // Intracellular K At Rest
//	double caimax [beats] ; // Peak Intracellular Ca

//	vmax=(Number *)calloc(beats,sizeof(Number));
//	dvdtmax=(Number *)calloc(beats,sizeof(Number));
//	apd=(Number *)calloc(beats,sizeof(Number));
//	toneapd=(Number *)calloc(beats,sizeof(Number));
//	ttwoapd=(Number *)calloc(beats,sizeof(Number));
//	rmbp=(Number *)calloc(beats,sizeof(Number));
//	nair=(Number *)calloc(beats,sizeof(Number));
//	cair=(Number *)calloc(beats,sizeof(Number));
//	kir=(Number *)calloc(beats,sizeof(Number));
//	if(vmax==NULL || dvdtmax==NULL || apd==NULL || toneapd==NULL || ttwoapd==NULL 
//		|| rmbp==NULL || nair==NULL || cair==NULL || kir==NULL
//		) exit(1);

//int i; // Stimulation Counter

	tmpname = "temp";

	sprintf(cmd, "/usr/bin/cpp -P %s > %s", argv[1],tmpname);
	if(system(cmd) == -1){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(1);
	}
	if((fpin=fopen(tmpname,"r"))==NULL){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(1);
	}
	if ((fp1 = fopen("para.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp2 = fopen("data.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp3 = fopen("nstate.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}

// parameter inputs
	input_para(fpin);

	if (var.write){
		if ((fp0 = fopen(argv[2],"w"))==NULL){
			fprintf(stderr, "%s cannot open.\n",argv[2]);
			exit(-1);
		}
	}

	xhplot(WINDOW, 700.0, 700.0, WHITE);
	xhplot(DIRECT, 0.0, 0.0, WHITE);

	for (ii = 0; ii < var.datas; ii++){
		long j;
		time = 0.0;
		tend = var.tend[ii];
		for (i = 0; i < NN; i++){ 
			x[i] = var.x0[ii][i];
		}
		if (var.half){
			h = 1.0*M_PI / var.m;
		}
		else {
			//h = 2.0*M_PI / var.m;
			h = 1.0 / var.m;
		}
		h *= var.tsign[ii];

		xddp.line_wid = var.line_wid[ii];
		xhplot(LINEATT,0,0,WHITE);

		
		// parameter values input.
		val_consts(fp1);
		printf("exit consts\n");
	
		// initial values input.
		initial_mem();
		printf("exit memory initialization\n");

		printf("Istim=%lf\n",var.Istim_base);

		// Tablize exp functions.	
		printf("start tablization\n");
		make_ExpTable();
		printf("finished tablization\n");

		// Initialization time
		//time -= h;
		//var.dt = h;
		var.beat = 0;

		for (var.beat=0;var.beat<beats;var.beat++){
			eventloop(fp1,&mode,&P,x);

			for (j = 0; j< (var.m * var.l * var.BCL); j++){
				t = h*j;
				time += h;
				if (fabs(time) > tend &&  tend != 0.0) break;

				if (time-(var.BCL*var.beat) >= 0.0 && time-(var.BCL*var.beat) < 1.0){
					var.Istim = var.Istim_base;
				} else {
					var.Istim = 0;
				}

				v_old = x[0];

				eular(NN,h,x,t);
				//runge(NN,h,x,t);
	
				/*printf("%lf ",time);
				for(w=0;w<NN;w++){
					if(w!=NN-1){
						printf("%lf ",x[w]);
					} else {
						printf("%lf\n",x[w]);
					}
				} */
				//printf("%lf %lf %lf %lf\n",time,ina.fast,ical.ica,ito.ik);

				dvdt_new = (x[0]-v_old)/h;

		/*		if(var.beat>=0){
					if (x[0] > vmax[var.beat] )
						vmax[var.beat] = x[0];
					if (x[12] > caimax[var.beat] )
						caimax[var.beat] = x[12];
					if (dvdt_new > dvdtmax[var.beat] ){
						dvdtmax[var.beat] = dvdt_new;
						toneapd[var.beat] = time;
					}
					if (dvdt_new < 0 && x[0] >= (vmax[var.beat] -0.9*(vmax[var.beat]-rmbp[var.beat]) ) )
						ttwoapd[var.beat] = time;
				}*/

				if (var.pflag) orbit(&mode,x,dvdt_new);

				if (time>= (beats-3)*var.BCL){
					fprintf(fp2,"%lf %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",time,x[0],x[12],x[15],ina.fast,ical.ica,ito.ik,ikr.ik,iks.ik,ik1.ik,inak.inak,ncx.jncx,inab.na,icab.ca,ipca.ca,var.Istim);
				}
				
				dvdt = dvdt_new;

			} // end for j-loop

			fprintf(fp3,"#beats=%d\n",var.beat);
			for(w=0;w<NN;w++){
				fprintf(fp3,"%e\n",x[w]);
			}

			if(P==8){
				printf("%lf ",time);
				for(w=0;w<NN;w++){
					if(w!=NN-1){
						printf("%lf ",x[w]);
					} else {
						printf("%lf\n",x[w]);
					}
				} 
			}
			draw_p(&mode,P,x,dvdt);
			mouse(&mode,x,dvdt);
			if (fabs(time) > tend &&  tend != 0.0) break;

		} // end for while loop

	} // end for ii-loop

	// closed file and open memories
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	//fclose(fp4);
	//free(vmax);free(dvdtmax);free(apd);free(toneapd);free(ttwoapd);
	//free(rmbp);free(nair);free(cair);free(kir);
	closed_mem();
}


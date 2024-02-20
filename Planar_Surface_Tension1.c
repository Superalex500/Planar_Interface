/*Programa finalizado el 03/Feb/2023 */
/* CORRECCION 22/AGO/2023*/

//Bibliotecas incluidas//
#include <stdio.h>
#include <math.h>
#include<stdlib.h>


/*Declaracion de variables y constantes globales*/
int seed,*pseed;
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
#define MAX_GR_BINS 200
#define MAX_BEADS 28	
#define MAX_BINS 200	
#define MAX_PARTICLES 35000	
double RHO,TEMP,DISPL,DISPLCM,Lxx,Lyy,Lzz,lb,PI,S,SS,DBH,DBH2,DBHm,RHOeff,dx,zb,zb2,epsilon;
int NMOVE,NSUB,N,load_prev,P;
double acc[4][410],l,Densidad[1500],oblate_area,prolate_area;
double x[MAX_PARTICLES],y[MAX_PARTICLES],z[MAX_PARTICLES],radio2;
double x11[MAX_PARTICLES],z11[MAX_PARTICLES],y11[MAX_PARTICLES];
double total,a2_cnt=0.,a3_cnt=0.,old_x[MAX_BEADS],old_y[MAX_BEADS],old_z[MAX_BEADS],er_pro,er_o,er_pro1,er_o1;
double b1_cnt=0.,b2_cnt=0.,b3_cnt=0.;
double XA[20],YA[20],ZA[20],a1,a2,a3,b1,b2,b3,c1,c2;
double energia_cinetica=0.0,energia_potencial=0.0,energia_potencial1=0.0;
double bin[MAX_BINS],gr_bin[MAX_GR_BINS];
double xbl,ybl,zbl,xbl2,ybl2,zbl2,Lx,Ly,Lz,fact,MINIMO,delta_r,fact_gr,fact1;
double desp_cm[100],desp_b[100],zbld,zbld2;
int ctl=1,ctl_b=0,stp1=0,bmr=0,nmr=0,nx,cnt2=0,cc1=0,den=0,cnn=0;
/*Variables de ayuda*/
double ener=0.0,ener1=0.0;

double Longz;
char Perfil_Den[100], Configuracion_lb[100], Resultados_lb[100], Resultados_Parciales[100];




//Declaracion de Funciones//
int comienzo(void);
void montecarlo(void);
int condition(int n,int bead,int j);
double total_energy(int k);
int get_data(void);						
void print_data(void);
double min(double x,double y);	
void set_initial_array(int i);		
double potential(double rr);
void save_configuration(void);
void save_current_configuration(int n);
int move_particle(int n,int j);
void return_old_configuration(int n);
void compute_gr(void);
void print_gr(double number_measures);
void compute_uncertainties(double *unc_E,double *unc_a_1,double *unc_a_2,double *unc_a_3);
double ran3(int *idum);
int ajustando(double bead,double neck);
void reset(int i);
void boundary_condition_box(double *xp,double *yp,double *zp, int j);
void center_mass(double *Rx1,double *Ry1,double *Rz1);
void Print_Desity(void);
void density(void);
void Resultados(int i);
void error(double er, double er1);
double RadioBH_clasico(void);
double clasico(double r,double sigma);
void Pre_Termalizacion(void);
void boundary_condition_confinement(double *yp,double *zp);
int HS_Termalizacion(int n);
double kinetic_energy(int n);
double total_kinetic_energy(int h,int k);
void center_mass_neck(int index,double *Rx1,double *Ry1,double *Rz1);
void generar_vector(double r[3],double R[3]);
void rotate_necklace(int n);
int condition1(int n,int j);
void CM(void);
double Propagator(int h,int k);
double distance_p(int j,int i,int h,int k);
double pair_potential_energy(int n, int m,int h,int k);
double total_potential_energy(int h,int k);
double Potential_energy(int n);
double particle_energy(int m);
void areas(void);
void Corregir_coordenadas(void);

//*****Funcion principal*****//
int main(void)
{
	int i;
	i=comienzo();
	if(i==1)
		return 0;
	montecarlo();
	save_configuration();

return 0;
}

//*********************Funciones*******************//
int comienzo(void)
{
	int i;
	
	seed=-123456789;
	pseed=&seed;
	PI=acos(-1.0);

	if(get_data()==1)
	{
		printf("The file can't be read\n");
		return 1;
	}
	//DBH=RadioBH_clasico();
	DBH=1.0;
	DBH2=DBH*DBH;
	DBHm=DBH/2.0;
	/*Aqui� se fija la fase densa*/
	//RHO=(double)N/(l*l*l);
	RHO=1.5278;
	RHO=0.74;
	l=pow(((double)N/(RHO)),(1.0/3.0));
	xbl=l;
	ybl=l;
	zbl=l;
	xbl2=xbl/2.0;
	ybl2=ybl/2.0;
	zbl2=zbl/2.0;
	set_initial_array(load_prev);
	if(load_prev == 1) {zbl=Longz*l; zbl2=zbl/2.0;}
	/*Aqui� se fija la longitud del perfil de densidad*/
	zbld=3.0*l;        zbld = Longz*l;
	zbld2=zbld/2.0;
	nx=(int)(zbld/dx);
	for(i=0;i<=nx;++i)
		Densidad[i]=0;
	/*Aqui se definen constantes del propagador*/
	fact=PI*P*TEMP/(lb*lb);
	/*Aqui se define la elongacion de la caja*/
//	zb=6.0*zbl;   
     zb=Longz*l;
	zb2=zb/2.0;
	areas();
return 0;
}
void montecarlo(void)
{
	double old_energy,new_energy,cnt_E=0;
	double A_beads,A_neck,random_num,kener,uener;
	double oblate_energy,prolate_energy,well_energy_sub;
	double dif,dif1,diff,err,err1;
	int n,stp=0,ctrl=0,j,el;
	int index_acc=0,n_measures=0,i;
 
	for(j=0;j<5;++j)
	{
			total=total_energy(0)/(double)N;
			reset(j);
			i=0;
			while((cnt2<NMOVE)&&(i==0))
			{
				stp1+=1;
				n=(int)(ran3(pseed)*N);			
				save_current_configuration(n);
				old_energy=particle_energy(n);
				ctrl=move_particle(n,j);
				new_energy=particle_energy(n);
				diff=new_energy-old_energy;
				if(diff>1.0E-10)	// Metropolis algorithm
				{
					random_num=ran3(pseed);
					if(random_num>exp(-diff/TEMP))
					{										//Aqui se rechaza por causa de la energia
						return_old_configuration(n);
						if(stp1%(P+1)==0)	//rejected movements
							++nmr;			//Necklace
						else				
							++bmr;			//Bead
					}	
				}		

				if(stp1%(P+1)==0)
						stp+=1;
				if((stp1%1000==0)&&(j==4))
					density();

				if(stp1%NSUB==0)//take a configuration each NSUB movements and print on screen
				{
					++cnt2;//counts the number of measurements
					well_energy_sub=Propagator(0,0);
					oblate_energy=Propagator(1,3);
					prolate_energy=Propagator(2,3);

					//Oblate energy
					dif=oblate_energy-well_energy_sub;
					energia_potencial+=dif;
					a2_cnt+=dif*dif;
					a3_cnt+=dif*dif*dif;

					a1=energia_potencial/(double)cnt2;		
					a2=(a2_cnt/((double)cnt2))-a1*a1;
					a3=(a3_cnt/cnt2-3.0*(a2_cnt/cnt2)*a1+2.0*a1*a1*a1);
					err=a2/oblate_area;

					a1=a1/(oblate_area);						//the first perturbation terms
					a2=a2/(-2.0*TEMP*oblate_area);			//the second perturbation terms
					a3=a3/(6.0*TEMP*TEMP*oblate_area);		//the third perturbation terms

					//Prolate energy
					dif1=prolate_energy-well_energy_sub;
					energia_potencial1+=dif1;
					b2_cnt+=dif1*dif1;
					b3_cnt+=dif1*dif1*dif1;

					b1=energia_potencial1/(double)cnt2;		
					b2=(b2_cnt/((double)cnt2))-b1*b1;
					b3=(b3_cnt/cnt2-3.0*(b2_cnt/cnt2)*b1+2.0*b1*b1*b1);
					err1=b2/prolate_area;
                  // printf("XX %d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",cnt2,b1,b2,err1,b2_cnt/((double)cnt2),prolate_area,dif1);
					b1=b1/(prolate_area);						//the first perturbation terms
					b2=b2/(-2.0*TEMP*prolate_area);			//the second perturbation terms
					b3=b3/(6.0*TEMP*TEMP*prolate_area);		//the third perturbation terms

					A_beads=1.0-((double)(bmr)/(double)(stp1-stp));
					A_neck=1.0-((double)(nmr)/(double)stp);
					if(j==4 && (cnt2>1))
					{
						error(err,err1);
						Print_Desity();	
					}
					//uener=ener/((double)(N*cnt2));
					//kener=ener1/((double)(N*cnt2));
					if(j==4)
						printf("%f %d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf %f %f\n",lb,cnt2,a1,er_o,b1,er_pro,A_beads,A_neck,dif,dif1);
						//printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",cnt2,(a1+b1)/2,(a2+b2)/2,(a3+b3)/2,A_beads,A_neck,DISPL,DISPLCM);
						//printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",cnt2,a1,er_o,a2,a3,A_beads,A_neck);
					else
						printf("%f %d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf %f %f\n",lb,cnt2,a1,b1,A_beads,DISPL,A_neck,DISPLCM,dif,dif1);
						

					if((j==0)||(j==2))						//Busca desplazamientos optimos
						i=ajustando(A_beads,A_neck);

					if((j==4)&&(cnt2%20==0))				//Guarda resultados parciales
					{
						Resultados(0);
						//save_configuration();
					}

					stp1=0;
					stp=0;
					nmr=0;
					bmr=0;
				}
			}
	}
	Resultados(1);

}
void error(double er, double er1)
{
	er_o=sqrt(fabs(er)/((double)cnt2-1));
	//er_o1=sqrt(o1/((double)(cnt2*(cnt2-1))*fact*fact*TEMP*TEMP*4.0));

	er_pro=sqrt(fabs(er1)/((double)cnt2-1));
	//er_pro1=sqrt(p1/((double)(cnt2*(cnt2-1))*fact1*fact1*TEMP*TEMP*4.0));
}
void reset(int i)
{
	energia_potencial=0;
	energia_potencial1=0;
	nmr=0;
	bmr=0;
	cnt2=0;
	stp1=0;
	a2_cnt=0;
	a3_cnt=0;
	b2_cnt=0;
	b3_cnt=0;
	ctl=1;
	ener=0;
	ener1=0;

	if(i==0)
	{
		DISPL=0.1;
		desp_b[0]=DISPL;			//Guarda los desplazamientos anteriores
		DISPLCM=0.1;
		desp_cm[0]=DISPLCM;

		NMOVE=200000;
		NSUB=10000;
	}
	else if(i==1)
	{
		NMOVE=30;
		NSUB=500000;
	}
	else if(i==2)
	{
		desp_b[0]=DISPL;
		desp_cm[0]=DISPLCM;
		NMOVE=1000;
		NSUB=5000;
		zbl=zbl*6.0;
		zbl2=zbl/2.0;
		Corregir_coordenadas();
	}
	else if(i==3)		//pasos para la temalizacion
	{
		NMOVE=40;		//Cuantos promedios se toman
		NSUB=1000000;	//Cada cuantos movimientos de particulas se toman promedios NSUB=1000000;
	}
	else if(i==4)		//pasos para tomar los promedios
	{
		NMOVE=27000;	   //Cuantos promedios se toman
		NSUB=200000;   //Cada cuantos movimientos de particulas se toman promedios
	}

	if(i==0||i==4)
		print_data();

	if(i==0)
		printf("\nBuscando Desplazamientos Optimos\n\n");
	else if(i==1)
		printf("\nInicio de Pre-termalizacion\n\n");
	else if(i==2)
		printf("\nBuscando nuevos desplazamientos ---- Se libera a longitud z de la caja a 6l\n\n");
	else if(i==3)
		printf("\nTermalizando caja de simulacion\n\n");
	else if(i==4)
		printf("\nInicio de Corrida\n\n");
	
	if(i==4)
	{
		printf("\nTotal initial energy=%lf  l=%lf\n",total,l);
		printf("\nNMOVE\t\ta1\t\terror\t\tb1\t\terror\t\tAcep_Beads\t\tAcep_neck\n\n");
	}
	else
	{
		printf("\nTotal initial energy=%lf  l=%lf\n",total,l);
		printf("\nNMOVE\t\ta1\t\tb1\t\tA_Beads\t\tDISPL\t\tA_Neck\t\tDISPLCM\n\n");
	}
}
/*  Corregir coordenadas al incrementar caja */
void Corregir_coordenadas(void){
	int i, m, index;
  for(i=0;i<N;++i) {
		index=i*(P+1);
		for(m=1;m<P;++m){
			if ( fabs(z[index] - z[index+m])> xbl/2.0) {
			 if (z[index] < 0.0) z[index+m] -= xbl;
			 else z[index+m] += xbl;			        }
			            }		
	}

/**************/
}
/// **************************** ////
int get_data(void)
{
	FILE *arc;
	arc = fopen("entrada.in","r");
	fscanf(arc,"%lf", &TEMP);
	fscanf(arc,"%lf", &lb);
	fscanf(arc,"%d", &N);
	fscanf(arc,"%lf", &epsilon);
	fscanf(arc,"%d", &load_prev);
	fscanf(arc,"%lf",&Longz);
	fscanf(arc,"%s",&Perfil_Den);
	fscanf(arc,"%s",&Configuracion_lb);
	fscanf(arc,"%s",&Resultados_lb);
	fscanf(arc,"%s",&Resultados_Parciales);
	fclose(arc);
	printf("T=%lf lb=%lf N=%d eps=%lf conNew=%d long_z =%lf\nPerfil en: %s\n", TEMP, lb, N, epsilon, load_prev,Longz,Perfil_Den);
	printf("Conf en: %s\nRes en: %s\nRes Par en: %s\n", Configuracion_lb,Resultados_lb,Resultados_Parciales);
	P=6;
	dx=0.1;
	//Load_prev=0 genera una configuracion
	//Load_prev=1 parte de una configuracion

return 0;
}
int ajustando(double bead,double neck)//Esta funcion encuentra los desplazamientos Optimos
{
	int i=0,ctl_b=0,ctl_cm=0;
	double des_b,des_cm,des_cl,aux_b,aux_cm,aux_cl,rang_i,rang_s;

	rang_i=0.00008;
	rang_s=0.96;

	des_b=0.33-bead;		//Desviacion respecto de 0.33 en la aceptación
	des_cm=0.33-neck;		//Desviacion respecto de 0.33 en la aceptación

	desp_b[ctl]=DISPL;			//Guarda los desplazamientos anteriores
	desp_cm[ctl]=DISPLCM;		//Guarda los desplazamientos anteriores
	
//*************Ajuste desplazamiento de beads******************************//
	if(fabs(des_b)>0.01)
	{
		ctl_b=0;
		if(des_b>0.0)				//Si entra, reducir DISPL
		{
			aux_b=DISPL-desp_b[ctl-1];
			if(aux_b>0.0)
				DISPL-=aux_b/2.0;
			if(aux_b<=0.0)
				DISPL-=desp_b[ctl]*0.1;	
		}
		if(des_b<0.0)				//Si entra, aumentar DISPL
		{
			aux_b=DISPL-desp_b[ctl-1];
			if(aux_b>=0.0)
				DISPL+=desp_b[ctl]*0.1;
			if(aux_b<0.0)
				DISPL+=fabs(aux_b)/2.0;	
		}
		if(DISPL>rang_s)
		{
			ctl_b=1;
			DISPL=rang_s;
		}
		if(DISPL<rang_i)
		{
			ctl_b=1;
			DISPL=rang_i;
		}
	}
	else
		ctl_b=1;

//*********************Ajuste desplazamiento de cm*************************//
	if(fabs(des_cm)>0.01)
	{
		ctl_cm=0;
		if(des_cm>0.0)				//Si entra, reducir DISPLCM
		{
			aux_cm=DISPLCM-desp_cm[ctl-1];
			if(aux_cm>0.0)
				DISPLCM-=aux_cm/2.0;
			if(aux_cm<=0.0)
				DISPLCM-=desp_cm[ctl]*0.1;	
		}
		if(des_cm<0.0)				//Si entra, aumentar DISPLCM
		{
			aux_cm=DISPLCM-desp_cm[ctl-1];
			if(aux_cm>=0)
				DISPLCM+=desp_cm[ctl]*0.1;
			if(aux_cm<0)
				DISPLCM+=fabs(aux_cm)/2.0;	
		}
		if(DISPLCM>rang_s)
		{
			ctl_cm=1;
			DISPLCM=0.5;
		}
		if(DISPLCM<rang_i)
		{
			ctl_cm=1;
			DISPLCM=rang_i;
		} 
	}
	else
		ctl_cm=1;


	if((ctl_cm==1)&&(ctl_b==1))
		i=1;

	++ctl;
	
return i;
}
void print_data(void)
{
	printf("\n******Metropolis Algorithm******\n\n");
	printf("\nN=%d RHO=%f TEMP=%f DISPL=%f\n\n",N,RHO,TEMP,DISPL);
	printf("NMOVE=%d NSUB=%d l=%lf\n\n",NMOVE,NSUB,l);
}
void set_initial_array(int i)
{
	double coordx,coordy,coordz;
	double x1[MAX_PARTICLES],y1[MAX_PARTICLES],z1[MAX_PARTICLES];
	int n,m=0,k;
	char a[1024];
	char espe[1024];
	FILE *fp;
//	snprintf(a, sizeof(a),"N%d.xyz",N);

	if(i==1)
	{
		//Read the data and load in the array x[n],y[n], z[n]
		fp=fopen(Configuracion_lb,"r");
		if (fp!=NULL)
		{
			for(n=0;n<N;++n)
			{
				fscanf(fp,"%s\t%lf\t%lf\t%lf",espe,&coordx,&coordy,&coordz);
			//	printf("%s\t%lf\t%lf\t%lf\n",espe,coordx,coordy,coordz);
				for(i=0;i<(P+1);++i)
				{
					x[m]=coordx;
					y[m]=coordy;
					z[m]=coordz;
					m+=1;
				}
			}
			fclose(fp);
		}
		else 
		{
			printf("\nNo se pudo cargar la configuracion inicial\n\n");
			exit(-1);
		}
	}
	else
	{
			for(n=0;n<N;++n)
			{
				x[n]=xbl*(ran3(pseed)-0.5);
				y[n]=ybl*(ran3(pseed)-0.5);
				z[n]=zbl*(ran3(pseed)-0.5);
			}
			Pre_Termalizacion();
			for(n=0;n<N;++n)
				for(m=0;m<=P;++m)
				{
					k=n*(P+1);
					x1[m+k]=x[n];
					y1[m+k]=y[n];
					z1[m+k]=z[n];
				}

			for(m=0;m<(N*(P+1));++m)
			{
				x[m]=x1[m];
				y[m]=y1[m];
				z[m]=z1[m];
			}
			//save_configuration();
	}
}
void Pre_Termalizacion(void)
{
	int n,i,j,ctrl,can;
	double oldx,oldy,oldz;

	DISPL=0.1;
	can=3000;

	for(i=0;i<can;++i)
	{	
		for(j=0;j<1000;++j)
		{
			n=(int)(ran3(pseed)*N);			
			oldx=x[n];
			oldy=y[n];
			oldz=z[n];
			ctrl=HS_Termalizacion(n);     ///ctrl=0 no hay traslape ctrl=1 hay traslape
			if(ctrl==1)
			{
				ctrl=0;
				x[n]=oldx;
				y[n]=oldy;
				z[n]=oldz;
			}
		}
	}
}
int HS_Termalizacion(int n)
{	
	double newx,newy,newz,dispx,dispy,dispz,*newx1,*newy1,*newz1;
	int index,bead,m,kk=0;

	newx1=&newx;
	newy1=&newy;
	newz1=&newz;

	newx=x[n]+DISPL*(xbl*ran3(pseed)-xbl2);
	newy=y[n]+DISPL*(ybl*ran3(pseed)-ybl2);
	newz=z[n]+DISPL*(zbl*ran3(pseed)-zbl2);


	boundary_condition_box(newx1,newy1,newz1,0);
	x[n]=newx;
	y[n]=newy;
	z[n]=newz;
	kk=condition1(n,0);

return kk;
}
//Esta funcion revisa que no haya traslapes entre particulas
int condition1(int n,int j)
{
	int cc=0,i;
	double D;

	for(i=0;(i<N)&&(cc==0);++i)
	{
		if(n!=i)
		{
			D=distance_p(n,i,0,0);
			if(D<DBH2)
				cc=1;
		}
	}

return cc;
}
double particle_energy(int m)
{
	double cinetica=0,potential=0,fase=0;

	cinetica=kinetic_energy(m);
	potential=Potential_energy(m);
	//fase=JOHSP(m);

return potential+cinetica-fase;
}
double Potential_energy(int n)
{	//return the total potential energy of particle n
	double potential=0.0;
	int i;

	for(i=0;i<N;++i)
		potential+=pair_potential_energy(n,i,0,0);

return potential;
}
double pair_potential_energy(int n, int m,int h,int k)
{	//return the potential energy due the bead-bead interaction between the necklace n and m
	double energy=0.0;
	int index1,index2,j;

	if(n!=m)
	{
		index1=n*(P+1);
		index2=m*(P+1);
		for(j=0;j<P;++j)	
			energy+=potential(distance_p(index1+j,index2+j,h,k));
	}
	energy=energy/(double)P;

return energy;
}
double kinetic_energy(int n)
{	//return the kinetic energy due to the shape of the necklace n
	double energy=0.0;
	int index,m;

	index=n*(P+1);
	for(m=0;m<P;++m)
		energy+=distance_p(index+m,index+m+1,0,0);

return fact*energy;	
}
double Propagator(int h,int k)
{	/*return the propagator
h=0		Normal
h=1		oblate	
h=2		prolate 
k=1,2,3 en x,y,z   
h  Tipo de perturbacion
k  Direccion de perturbacion*/
	double potential=0,cinetica=0,fase=0;
	
	potential=total_potential_energy(h,k);
	cinetica=total_kinetic_energy(h,k);

return potential+cinetica-fase;
}
double total_energy(int k)
{	//return the total energy
	double potential,cinetica,RR;
	int index1,index2,n,m,j;
	
	potential=total_potential_energy(0,0);
	cinetica=total_kinetic_energy(0,0);

return potential+cinetica;
}
double total_potential_energy(int h,int k)
{	//return the total potential energy

	double potential_E=0.0;
	int i,j;

	for(i=0;i<(N-1);++i)
		for(j=i+1;j<N;++j)
			potential_E+=pair_potential_energy(i,j,h,k);

return potential_E;
}
double total_kinetic_energy(int h,int k)
{
	int i,m,index;
	double energy=0;

	for(i=0;i<N;++i)
	{
		index=i*(P+1);
		for(m=0;m<P;++m)
			energy+=distance_p(index+m,index+m+1,h,k);
	}
		
return fact*energy;
}
double distance_p(int j,int i,int h,int k)
{	//return the distance to the square
	double xd,yd,zd,rr,b;

	xd=x[i]-x[j];
	yd=y[i]-y[j];
	zd=z[i]-z[j];

	//Convencion de minima imagen
	if(xd>xbl2)
		xd=xd-xbl;
	else if(xd<-xbl2)
		xd=xd+xbl;
	if(yd>ybl2)
		yd=yd-ybl;
	else if(yd<-ybl2)
		yd=yd+ybl;
	if(zd>zbl2)
		zd=zd-zbl;
	else if(zd<-zbl2)
		zd=zd+zbl;

	if(h==0)
	{
		rr=xd*xd+yd*yd+zd*zd;
	}
	else
	{
		if(h==1)
			b=1.0+epsilon;
		else
			b=1.0-epsilon;


		if(k==1)
			rr=((xd*xd)/(b*b))+b*((yd*yd)+(zd*zd));
		else if(k==2)
			rr=b*((xd*xd)+(zd*zd))+((yd*yd)/(b*b));
		else if(k==3)
			rr=b*((xd*xd)+(yd*yd))+((zd*zd)/(b*b));
	}

return rr;
}
double potential(double rr)
{	//Here the potential between necklaces must be defined
	double energy=0.0,XR,XRR,corte;

	corte=6.25;
	if(rr<=corte)
	{
		XR=rr;
		XRR=pow(XR,3);
		energy=4*(1.0/XRR-1.0)/XRR + 0.0163;
	}
//+0.0163
return energy;
}
void save_configuration(void)
{
	int nn;
	char a[1024];
	FILE *fp;
//	snprintf(a, sizeof(a), "Config_lb%1.1lf_N%d.xyz",lb,N);
	fp=fopen(Configuracion_lb,"w");

	if(fp!=NULL)
	{
		for(nn=0;nn<(N*(P+1));++nn)//Guarda en el formato H,X,Y,Z
			fprintf(fp,"H\t%lf\t%lf\t%lf\n",x[nn],y[nn],z[nn]);

		fclose(fp);
	}
	else
		printf("\nError al guardar configuracion final\n\n");
}
void save_current_configuration(int n)
{	//save the current configuration of the particle n
	int index,m;
	index=n*(P+1);

	for(m=0;m<=P;++m)
	{
		old_x[m]=x[index+m];
		old_y[m]=y[index+m];
		old_z[m]=z[index+m];
	}
}
void boundary_condition_box(double *xp,double *yp,double *zp, int j)
{	//It returns the new coordinates according to boundary conditios
	//Periodic boundary conditios
	//Válidas para una caja centrada en el origen

	if(j<2)
	{
		if(*xp>xbl2)
			*xp=*xp-xbl;
		else if(*xp<-xbl2)
			*xp=*xp+xbl;

		if(*yp>ybl2)
			*yp=*yp-ybl;
		else if(*yp<-ybl2)
			*yp=*yp+ybl;

		if(*zp>zbl2)
			*zp=*zp-zbl;
		else if(*zp<-zbl2)
			*zp=*zp+zbl;
	}
	else
	{
		if(*xp>xbl2)
			*xp=*xp-xbl;
		else if(*xp<-xbl2)
			*xp=*xp+xbl;

		if(*yp>ybl2)
			*yp=*yp-ybl;
		else if(*yp<-ybl2)
			*yp=*yp+ybl;

		if(*zp>zb2)
			*zp=*zp-zb;
		else if(*zp<-zb2)
			*zp=*zp+zb;
	}
}
void areas(void)
{
	oblate_area=epsilon*l*l*2.0;
	prolate_area=(-1.0)*oblate_area;
}
void boundary_condition_confinement(double *yp,double *zp)
{	//It returns the new coordinates according to boundary conditios
	//Periodic boundary conditios with confinement conditions
	//Validas para una caja centrada en el origen

	if(*yp>ybl2)
		*yp=*yp-ybl;
	else if(*yp<-ybl2)
		*yp=*yp+ybl;

	if(*zp>zbl2)
		*zp=*zp-zbl;
	else if(*zp<-zbl2)
		*zp=*zp+zbl;
}
int move_particle(int n,int j)
{	
	double newx,newy,newz,dispx,dispy,dispz,*newx1,*newy1,*newz1;
	int index,bead,m,kk=0;

	newx1=&newx;
	newy1=&newy;
	newz1=&newz;
	index=n*(P+1);
	if(stp1%(P+1)==0)//trying to move the whole "n" necklace
	{
		rotate_necklace(n);

		dispx=DISPLCM*(xbl*ran3(pseed)-xbl2);
		dispy=DISPLCM*(ybl*ran3(pseed)-ybl2);
		dispz=DISPLCM*(xbl*ran3(pseed)-xbl2);

		for(m=0;((m<=P)&&(kk==0));++m)	//move all the beads in the same direction
		{
			newx=x[index+m]+dispx;
			newy=y[index+m]+dispy;
			newz=z[index+m]+dispz;
			boundary_condition_box(newx1,newy1,newz1,j);
			x[index+m]=newx;
			y[index+m]=newy;
			z[index+m]=newz;
			//kk=condition(n,m,j);
		}
	}
	else	//trying to move only a bead of the "n" necklace
	{
		bead=(int)(ran3(pseed)*P);

		newx=x[index+bead]+DISPL*(xbl*ran3(pseed)-xbl2);
		newy=y[index+bead]+DISPL*(ybl*ran3(pseed)-ybl2);
		newz=z[index+bead]+DISPL*(xbl*ran3(pseed)-xbl2);

		boundary_condition_box(newx1,newy1,newz1,j);
		x[index+bead]=newx;
		y[index+bead]=newy;
		z[index+bead]=newz;
		if(bead==0)
		{
			x[index+P]=newx;
			y[index+P]=newy;
			z[index+P]=newz;
		}
		else if(bead==P)
		{
			x[index]=newx;
			y[index]=newy;
			z[index]=newz;
		}
		//kk=condition(n,bead,j);
	}
return kk;
}
//Esta funcion revisa que no haya traslapes entre las cuentas correspondientes
int condition(int n,int bead,int j)
{
	int cc=0,i,index;
	double D;

	for(i=0;(i<N)&&(cc==0);++i)
	{
		if(n!=i)
		{
			index=i*(P+1);
			D=distance_p(index+bead,n*(P+1)+bead,0,0);
			if(D<DBH2)
				cc=1;
		}
	}

return cc;
}
void rotate_necklace(int n)
{	//Funcion que rota un Angulo aleatorio alrededor del centro de masa (Rodriguez formula)
	double Rx,Ry,Rz,*Rx1,*Ry1,*Rz1;
	double a[3],R[3];
	double x1,y1,z1;
	double t,cost,sint,dot;
	int index,i;

	Rx1=&Rx;
	Ry1=&Ry;
	Rz1=&Rz;

	index=n*(P+1);
	center_mass_neck(index,Rx1,Ry1,Rz1);
	R[0]=Rx;
	R[1]=Ry;
	R[2]=Rz;
	generar_vector(a,R);		//Genera un vector aleatorio sobre el cual se rotara la molecula
	t=2*PI*(2.0*ran3(pseed)-1.0);	//Genera un Angulo aleatorio de rotacion
	cost=cos(t);
	sint=sin(t);

	for(i=0;i<P;++i)
	{
		x1=x[index+i]-Rx;
		y1=y[index+i]-Ry;
		z1=z[index+i]-Rz;
/*Minima Imagen*/	 
	 if (x1 < -xbl2 ) x1 += xbl;
	 else if(x1 > xbl2 ) x1 -= xbl;
     if (y1 < -ybl2 ) y1 += ybl;
	 else if(y1 > ybl2 ) y1-=ybl;
	 if (z1 < -zbl2 ) z1 += zbl;
	 else if(z1 > zbl2 ) z1-=zbl;
/*****/	 
	 
		dot=(a[0]*x1+a[1]*y1+a[2]*z1)*(1.0-cost);

		x[index+i]=R[0]+x1*cost+(a[1]*z1-a[2]*y1)*sint+a[0]*dot;
		y[index+i]=R[1]+y1*cost+(a[2]*x1-a[0]*z1)*sint+a[1]*dot;
		z[index+i]=R[2]+z1*cost+(a[0]*y1-a[1]*x1)*sint+a[2]*dot;
	}
		x[index+P]=x[index];
		y[index+P]=y[index];
		z[index+P]=z[index];
}
void center_mass_neck(int index,double *Rx1,double *Ry1,double *Rz1)
{
	int a;
	double Rx=0.0,Ry=0.0,Rz=0.0;
        Rx+=x[index];
		Ry+=y[index];
		Rz+=z[index];
	for(a=1;a<P;++a)
	{
	 if (fabs(x[index] - x[index+a]) <= xbl2 ) Rx+=x[index+a];
	 else {
	            if(x[index+a] < x[index]) Rx+=x[index+a] + xbl; 
				else Rx += x[index+a] - xbl;    }
	     	      
	 if (fabs(y[index] - y[index+a]) <= ybl2 ) Ry+=y[index+a];
	 else {
	       if ( y[index+a] < y[index] ) Ry += y[index+a] + ybl;
	       else  Ry += y[index+a] - ybl;  }
	      
	if (fabs(z[index] - z[index+a]) <= zbl2 ) Rz+=z[index+a];
	else {
	      if ( z[index+a] < z[index] ) Rz += z[index+a] + zbl;
	      else  Rz += z[index+a] - zbl;	  } 	   
	/*   
		Rx+=x[index+a];
		Ry+=y[index+a];
		Rz+=z[index+a];*/
	}
	Rx/=(double)P;
	Ry/=(double)P;
	Rz/=(double)P;
/*Periodicas*/	 
	 if (Rx < -xbl2 ) Rx += xbl;
	 else if(Rx > xbl2 ) Rx -= xbl;
     if (Ry < -ybl2 ) Ry += ybl;
	 else if(Ry > ybl2 ) Ry-=ybl;
	 if (Rz < -zbl2 ) Rz += zbl;
	 else if(Rz > zbl2 ) Rz-=zbl;
/**************/	 
	*Rx1=Rx;
	*Ry1=Ry;
	*Rz1=Rz;
}
//Esta funcion genera un vector aleatorio al rededor del cual rotar
void generar_vector(double r[3],double R[3])
{
	double xx,b1x,b2x,b3x;
	
	r[0]=ran3(pseed);
	r[1]=ran3(pseed);
	r[2]=ran3(pseed);
	
	b1x=R[0]-r[0];
	b2x=R[1]-r[1];
	b3x=R[2]-r[2];

	xx=sqrt(b1x*b1x+b2x*b2x+b3x*b3x);
	r[0]=b1x/xx;
	r[1]=b2x/xx;
	r[2]=b3x/xx;

}
void return_old_configuration(int n)
{
	int index,m; 
	
	index=n*(P+1);
	for(m=0;m<=P;++m)
	{
		x[index+m]=old_x[m];
		y[index+m]=old_y[m];
		z[index+m]=old_z[m];
	}
}
void compute_uncertainties(double *unc_E,double *unc_a_1,double *unc_a_2,double *unc_a_3)
{	//computes the estandar deviation from the sub-mean values
	double x[4],x2[4];
	int i,j;

	for(i=0;i<4;++i)//setting the counters to zero
	{
		x[i]=0.0;
		x2[i]=0.0;
	}
	for(j=1;j<4;++j)
		for(i=0;i<400;++i)
		{
			x[j]+=acc[j][i];
			x2[j]+=acc[j][i]*acc[j][i];
		}

	//*unc_E=sqrt(x2[0]/400.0-x[0]*x[0]/160000);
	*unc_a_1=sqrt(x2[1]/400.0-x[1]*x[1]/160000);
	*unc_a_2=sqrt(x2[2]/400.0-x[2]*x[2]/160000);
	*unc_a_3=sqrt(x2[3]/400.0-x[3]*x[3]/160000);
}
double ran3(int *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}
void CM(void)
{
	double Rx,Ry,Rz,*Rx1,*Ry1,*Rz1;
	int index,i;

	Rx1=&Rx;
	Ry1=&Ry;
	Rz1=&Rz;

	for(i=0;i<N;++i)
	{
		index=i*(P+1);
		center_mass_neck(index,Rx1,Ry1,Rz1);
		x11[i]=Rx;
		y11[i]=Ry;
		z11[i]=Rz;
	}
}
void density(void)
{
	double xd,yd,zd,rr,jj;
	int i,j,c=0;

	++den;
	CM();
	
	for(i=0;i<N;++i)
	{
		j = (int)((z11[i]+zb2)/dx);
		Densidad[j]+=1.0;
	}
}
void Print_Desity(void)
{
	int i;
	double v;
	char a[1024];
	FILE *fp;
//	snprintf(a, sizeof(a), "Densidad_lb",N,lb);
	
	fp=fopen(Perfil_Den,"w");

	v=(dx*xbl*ybl*(double)den)/(DBH2*DBH);

	if(fp!=NULL)
	{
		for(i=0;i<nx;++i)//Guarda en el formato r,Rho
			fprintf(fp,"%lf\t%lf\n",(-zbld2+((double)i*dx)),Densidad[i]/v);

		fclose(fp);
	}
	else
		printf("\nError al guardar configuracion final\n\n");
}
void center_mass(double *Rx1,double *Ry1,double *Rz1)
{
	int a;
	double Rx=0,Ry=0,Rz=0;

	for(a=0;a<N;++a)
	{
		Rx+=x[a];
		Ry+=y[a];
		Rz+=z[a];
	}
	Rx/=N;
	Ry/=N;
	Rz/=N;

	*Rx1=Rx;
	*Ry1=Ry;
	*Rz1=Rz;
}
/****************Funcion que calcula el radio BH clasico ***********/
double RadioBH_clasico(void)
{
	double h,hx,x=0,y=0,y1=0,sigma,l;
	int i,j,n;

/******************Entrada de valores******************/
	n=100000;
	sigma=1.0;
/******************************************************/

	h=sigma/(3*n);
	hx=sigma/n;                          /****Tamaño de paso para la integración**/
	y+=h*clasico(x,sigma);
	for(j=1;j<n;++j)
	{
		x+=hx;
		if((j%2)==0)
			y+=h*2*clasico(x,sigma);
		else
			y+=h*4*clasico(x,sigma);

	}
	x+=hx;
	y+=h*clasico(x,sigma);

return y; //regresa el valor del radio BH
}
/****************Funcion que evalúa el potencial clasico***********/
double clasico(double r,double sigma)
{
	double v,f,t;
	t=TEMP;
	if(r!=0)
	{
		v=4*(pow(sigma/r,12)-pow(sigma/r,6));
		f=1-exp(-v/t);
	}
	else
		f=1;

return f;
}
void Resultados(int i)
{
	char a[1024],b[1024];
	FILE *fp,*fp1;
	++cnn;
	if(i==0)
	{
		save_configuration();
	//	snprintf(a, sizeof(a), "Resultados_Parciales_T%1.1lf_lb%1.1lf.txt",TEMP,lb);
		fp=fopen(Resultados_Parciales,"a");

		fprintf(fp,"%d %lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",cnn,lb,epsilon,a1,er_o,a2,a3,b1,er_pro,b2,b3,(b1+a1)/2.0);
		fclose(fp);
	}
	else
	{
		save_configuration();
	//	snprintf(a, sizeof(a), "Resultados_lb%1.1lf_T%1.1lf.txt",lb,TEMP);
		fp=fopen(Resultados_lb,"a");
		fprintf(fp,"%lf %lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",epsilon,lb,a1,er_o,a2,a3,b1,er_pro,b2,b3,(b1+a1)/2.0);
		fclose(fp);

		snprintf(b, sizeof(b), "Resultados.txt");
		fp1=fopen(b,"a");
		fprintf(fp1,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",lb,epsilon,a1,er_o,a2,a3,b1,er_pro,b2,b3,(b1+a1)/2.0);
		fclose(fp1);
	}
}

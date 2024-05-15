#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <string.h>

#define PI 3.14159
#define Lerr 0.01	// en km
#define MAX_ITERATIONS 100
#define ACCURACY 0.001
#define epsilon 1e-6
#define pi0 1. 	// proba détection entre t1 et 1er séisme


void calculs_initiaux(t,m,lat,lon,pi,phi,phi2,psi,psi2,X,Y,Z,t1,t2,N,E,P,P2,NE,NP,NP2,p,c,gamma,L0,lonmin,lonmax,latmin,latmax,zmin,zmax)	// calcul des coordonnées kilométriques, et surtout des coef phi & psi
double t[],m[],lat[],lon[],pi[],phi[],phi2[],psi[],psi2[],X[],Y[],Z[],t1,t2,p,c,gamma,L0,lonmin,lonmax,latmin,latmax,zmin,zmax;
int N,NE,NP,NP2,E[],P[],P2[];
{
	int i,j,n,ix,iy,iz,inside;
	double r,dv,q,x,y,z,L,ell,F,Xv[8],Yv[8],Zv[8];

	q=(2.+gamma)/3.;	

	// on suppose une étendue spatiale petite (projection simple)
	for(i=0;i<N;i++) {
		X[i]=111.*lon[i]*cos(lat[i]*PI/180.);
		Y[i]=111.*lat[i];	
	}

	// calcul des coef phi (pour les parents P) & phi2 (pour les parents P2), cf cahier labo 3/9/2021 (phi tient compte des problèmes de détection)
	for(i=0;i<NP;i++) {
		phi[i]=0;
		for(j=1;j<N;j++)
			if(t[j]>t[P[i]])
				phi[i]+=pi[j]*(exp((1.-p)*log(t[j]-t[P[i]]+c))-exp((1.-p)*log(t[j-1]-t[P[i]]+c)));
		phi[i]/=(1.-p);
	}

	for(i=0;i<NP2;i++) {		
		phi2[i]=0;
		for(j=1;j<N;j++)
			if(t[j]>t[P2[i]])
				phi2[i]+=pi[j]*(exp((1.-p)*log(t[j]-t[P2[i]]+c))-exp((1.-p)*log(t[j-1]-t[P2[i]]+c)));
		phi2[i]/=(1.-p);			
	}

	// calcul des coef psi (pour les parents P) & psi2 (pour les parents P2), qui sont des estimations des intégrales de volume des parents sur le volume cible

	Xv[0]=111.*lonmin*cos(latmin*PI/180.); Xv[1]=111.*lonmax*cos(latmin*PI/180.); Xv[2]=Xv[1]; Xv[3]=Xv[0]; Xv[4]=Xv[0]; Xv[5]=Xv[1]; Xv[6]=Xv[1]; Xv[7]=Xv[0];
	Yv[0]=111.*latmin; Yv[1]=Yv[0]; Yv[2]=111.*latmax; Yv[3]=Yv[2]; Yv[4]=Yv[0]; Yv[5]=Yv[0]; Yv[6]=Yv[2]; Yv[7]=Yv[2];
	Zv[0]=zmin; Zv[1]=Zv[0]; Zv[2]=Zv[0]; Zv[3]=Zv[0]; Zv[4]=zmax; Zv[5]=Zv[4]; Zv[6]=Zv[4]; Zv[7]=Zv[4];

	for(n=0;n<NP;n++) {
		i=P[n];
		L=L0*exp(0.5*log(10.)*m[i]); if(L<Lerr) L=Lerr; 	
		if(lon[i]>=lonmin & lon[i]<=lonmax & lat[i]>=latmin & lat[i]<=latmax & Z[i]>=zmin & Z[i]<=zmax) inside=1; else inside=0; 
		for(F=0,j=0;j<8;j++) {
			r=sqrt((X[i]-Xv[j])*(X[i]-Xv[j])+(Y[i]-Yv[j])*(Y[i]-Yv[j])+(Z[i]-Zv[j])*(Z[i]-Zv[j]));
			if(inside)
				F+=1-exp((1-q)*log(1+r*r*r/(L*L*L)));
			else
				F+=3.*(q-1)/(4.*PI)*exp(3.*(q-1)*log(L)-q*log(r*r*r+L*L*L));
		}
		psi[n]=F/8;
	}

	for(n=0;n<NP2;n++) {
		i=P2[n];
		L=L0*exp(0.5*log(10.)*m[i]); if(L<Lerr) L=Lerr; 	
		if(lon[i]>=lonmin & lon[i]<=lonmax & lat[i]>=latmin & lat[i]<=latmax & Z[i]>=zmin & Z[i]<=zmax) inside=1; else inside=0; 
		for(F=0,j=0;j<8;j++) {
			r=sqrt((X[i]-Xv[j])*(X[i]-Xv[j])+(Y[i]-Yv[j])*(Y[i]-Yv[j])+(Z[i]-Zv[j])*(Z[i]-Zv[j]));
			if(inside)
				F+=1-exp((1-q)*log(1+r*r*r/(L*L*L)));
			else
				F+=3.*(q-1)/(4.*PI)*exp(3.*(q-1)*log(L)-q*log(r*r*r+L*L*L));
		}
		psi2[n]=F/8;
	}
}


void normalization(om,om2,NE,NP,NP2)
double **om,**om2;
int NE,NP,NP2;
{
	int i,j; 
	double tot;

	for(j=0;j<NE;j++) { 
		tot=0;
	  	for(i=0;i<NP;i++) tot+=om[i][j];
	  	for(i=0;i<NP2;i++) tot+=om2[i][j];
		if(tot) {
		  	for(i=0;i<NP;i++) om[i][j]/=tot;
		  	for(i=0;i<NP2;i++) om2[i][j]/=tot;
		}
	}
}

void EXPECTATION(t,m,X,Y,Z,om,om2,P,P2,E,NP,NP2,NE,alpha,K,K2,p,c,L0,gamma) 	// calcul des om & om2
double X[],Y[],Z[],**om,**om2,alpha,p,c,L0,gamma,K,K2[],m[];
double t[];
int NE,NP,NP2,P[],P2[],E[];
{
	int i,j,ii,jj;
	double tot,r,L,q,tmp;	

	q=(2.+gamma)/3.;
	for(i=0;i<NP;i++) {
		ii=P[i];
		for(j=0;j<NE;j++) {
			jj=E[j];
			if(t[ii]<t[jj]) {
				r=sqrt((X[ii]-X[jj])*(X[ii]-X[jj])+(Y[ii]-Y[jj])*(Y[ii]-Y[jj])+(Z[ii]-Z[jj])*(Z[ii]-Z[jj]));
				L=L0*exp(0.5*log(10.)*m[ii]); if(L<Lerr) L=Lerr; 	
				tmp=alpha*m[ii]-p*log(t[jj]+c-t[ii])+log(3.*(q-1))+(3.*(q-1))*log(L)-log(4.*PI)-q*log(r*r*r+L*L*L);	
				om[i][j]=K*exp(tmp);
			} else om[i][j]=0; 	// causalité
		}
	}

	for(i=0;i<NP2;i++) {
		ii=P2[i];
		for(j=0;j<NE;j++) {
			jj=E[j];
			if(t[ii]<t[jj]) {
				r=sqrt((X[ii]-X[jj])*(X[ii]-X[jj])+(Y[ii]-Y[jj])*(Y[ii]-Y[jj])+(Z[ii]-Z[jj])*(Z[ii]-Z[jj]));
				L=L0*exp(0.5*log(10.)*m[ii]); if(L<Lerr) L=Lerr; 	
				tmp=alpha*m[ii]-p*log(t[jj]+c-t[ii])+log(3.*(q-1))+(3.*(q-1))*log(L)-log(4.*PI)-q*log(r*r*r+L*L*L);	
				om2[i][j]=K2[i]*exp(tmp);
			} else om2[i][j]=0; 	// causalité
		}
	}

	// normalisation
	normalization(om,om2,NE,NP,NP2);
}		

		
double MAXIMIZATION(m,phi,phi2,psi,psi2,om,om2,E,P,P2,NE,NP,NP2,alpha,K2) // optimisation de K et des K2
double **om,**om2,phi[],phi2[],psi[],psi2[],m[],alpha,K2[];
int NE,NP,NP2,P[],P2[],E[];
{
	int i,j;
	double K=0,denom=0;
	double num=0;

	for(i=0;i<NP;i++) {
		denom+=exp(alpha*m[P[i]])*phi[i]*psi[i];
		for(j=0;j<NE;j++) 
			num+=om[i][j];
	}
	K=num/denom;
	for(i=0;i<NP2;i++) {
		K2[i]=0;
		for(j=0;j<NE;j++) 
			K2[i]+=om2[i][j];
		K2[i]/=(exp(alpha*m[P2[i]])*phi2[i]*psi2[i]);
	}
	return(K);
}




void sortie(name,t,m,pi,psi,psi2,om,om2,alpha,p,c,K,K2,P,P2,ISE,N,NP,NP2,t1,t2)
double m[],pi[],alpha,p,c,K,K2[],t1,t2,**om,**om2,psi[],psi2[];
double t[];
int N,NP,NP2,P[],P2[],ISE[];
char name[];
{
	char name_ext[100];
	FILE *out;
	double tc;
	double pic,cumul=0,cumul_corr=0,*LAMBDA,tmp;
	int i,j;


	// calcul de LAMBDA(t) et de nrep
	LAMBDA=(double *)calloc(N,sizeof(double));	// LAMBDA est le nombre modélisé, cumulé, de séismes *observés* (donc ayant "subi" le filtre de non-détection) aux moments des séismes (-> tous les séismes pas juste ceux du volume cible)

	for(i=0;i<NP;i++)
		for(j=0;j<N-1;j++)
			if(t[j]>=t[P[i]]) {
				tmp=exp((1-p)*log(t[j+1]+c-t[P[i]]))-exp((1-p)*log(t[j]+c-t[P[i]]));
				LAMBDA[j]+=pi[j+1]*K*exp(alpha*m[P[i]])/(1.-p)*tmp*psi[i];
			}	
	for(i=0;i<NP2;i++)
		for(j=0;j<N-1;j++)
			if(t[j]>=t[P2[i]]) {
				tmp=exp((1-p)*log(t[j+1]+c-t[P2[i]]))-exp((1-p)*log(t[j]+c-t[P2[i]]));
				LAMBDA[j]+=pi[j+1]*K2[i]*exp(alpha*m[P2[i]])/(1.-p)*tmp*psi2[i];
			}

	sprintf(name_ext,"%s.dat",name); ///MODIFICATION ! --> Au lieu d'avoir .mod, j'ai changé en .dat car difficulté de lecture sur Linux
	out=fopen(name_ext,"w");
	for(j=0;j<N;j++) {
		cumul+=LAMBDA[j];
		fprintf(out,"%lf\t%lf\t%d\n",t[j],cumul,ISE[j]);
	}
	fclose(out);

	sprintf(name_ext,"%s.para",name);
	out=fopen(name_ext,"w");
	fprintf(out,"%lf\n\n",K); // valeur de K, puis valeurs de K spécifiques
	for(i=0;i<NP2;i++) fprintf(out,"%lf\n",K2[i]);
	fclose(out);



	free(LAMBDA);


}





int main(argc,argv)
int argc;
char *argv[];
// COMMAND LINE: > ./etas-amatrice-v4 "data" "parameters" "output"
// Calcule les paramètres K du modèle ETAS pour les données "data" (voir format + bas) et les paramètres ETAS contenus dans le fichier "parameters" (voir format + bas), selon
// ce qui est détaillé dans le cahier labo du 1/9/2021 et suite. Génère en sortie le fichier "output" (voir format + bas).
// ATTENTION : 
//	(i) c'est bien ici un modèle 3D en espace ; 
//	(ii) on tient compte des changements de la probabilité de détection au cours du temps ; 
//	(iii) on distingue les enfants E (tous les séismes m>=m0 dans la période ciblée) et les parents potentiels P (tous les séismes au delà de la magnitude mP) et P' (séismes m>=mp2) ; 
//	(iv) le catalogue d'entrée doit déjà être préalablement coupé à la magnitude m0.
// Modélise chaque séisme de "liste_specifique" avec un K indépendant. "liste_specifique" contient les indices (commençant à 1, pas à 0) des chocs concernés (indices dans le catalogue d'entrée "data").
//
// Format "data" : 1 ligne / séisme, (1) temps, (2) magnitude, (3) lat, (4) lon, (5) z, (6) pi (-> proba de détecter un séisme >m0 au moment de ce séisme, cf cahier labo)
//
// Format "parameters" : 16 paramètres 
// (1) t1 temps de début de l'intervalle analysé, (2) t2 temps de fin de l'intervalle analysé, 
// (3) m0 magnitude de coupure (ATTENTION : le catalogue en entrée dans "data" est sensé être déjà coupé à m0) ; 
// (4) mp magnitude minimale pour qu'un séisme soit considéré comme un parent potentiel ; (5) mp2 magnitude minimale pour qu'un séisme soit considéré comme un parent potentiel *spécifique* ; 
// (6) latmin ; (7) latmax ; (8) lonmin ; (9) lonmax ; (10) zmin ; (11) zmax : volume "cible" (les enfants sont dans ce volume ; par contre les parents ne sont pas obligatoirement dedans)
// (12) alpha (13) p (14) c (15) L0 = longueur (en km) de rupture pour un séisme de magnitude m0 ; (16) gamma ; avec alpha, p, c, gamma les paramètres ETAS.
//
// Format "output" : 1 ligne / séisme, (1) LAMBDA(t) (intégrale de lambda(x,y,z,t) sur tout le volume {x,y,z} et jusqu'au t courant) pour les séismes observés 
// (2) comme (1) mais en corrigeant des effets de détection (donc LAMBDA du 2 > LAMBDA du 1) ; (3) valeur de K (ou 0 si pas parent potentiel) ; (4) nb de répliques 
{ 
	double *m=NULL,m0,mp,mp2,**om,**om2,tot,alpha,p,c,L0,L,gamma,t1,t2,dx=ACCURACY+1.,K,*K2,K_old,*lat=NULL,*lon=NULL,*pi=NULL,tmp,*X,*Y,*Z=NULL,*phi,*phi2,*psi,*psi2,J,meanK,meanK_old;
	double *t=NULL,tmpt,latmin,latmax,lonmin,lonmax,zmin,zmax;
	int i,j,k,n,q,niter,*E=NULL,*P=NULL,*P2=NULL,*ISE=NULL,N,NE,NP,NP2;
	FILE *in,*out;
	char data[100],parameters[100],output[100],*line=NULL;
    	size_t len=0;

	sscanf(argv[1],"%s",data);
	sscanf(argv[2],"%s",parameters);
	sscanf(argv[3],"%s",output);

	// READ PARAMETERS
	in=fopen(parameters,"r");
	fscanf(in,"%lf %lf",&t1,&t2);  // t1, t2 : temps de début & fin de la période cible
	fscanf(in,"%lf %lf %lf",&m0,&mp,&mp2); // m0 = coupure en magnitude (ATTENTION : le catalogue en entrée dans "data" est sensé être déjà coupé à m0) ; mp = magnitude min pour les parents potentiels
	fscanf(in,"%lf %lf %lf %lf %lf %lf",&latmin,&latmax,&lonmin,&lonmax,&zmin,&zmax); 
	fscanf(in,"%lf %lf %lf %lf %lf",&alpha,&p,&c,&L0,&gamma); // paramètres ETAS imposés. L0 est la longueur de rupture (en km) pour m=m0; gamma est l'exposant du kernel spatial (typiquement 2<gamma<3 voire  <4).
	if(p==1) p=1.001;
	fclose(in);

	// READ EARTHQUAKE DATA
	in=fopen(data,"r");
	i=0; j=0; k=0; q=0;
	while(1) {	
		if(getline(&line,&len,in)==-1) break;
		sscanf(line,"%lf",&tmpt);
		if(tmpt<t2) {
			t=(double *)realloc(t,(i+1)*sizeof(double)); t[i]=tmpt; 
			sscanf(line,"%*f %lf",&tmp);
			m=(double *)realloc(m,(i+1)*sizeof(double)); m[i]=tmp-m0;
			ISE=(int *)realloc(ISE,(i+1)*sizeof(int)); ISE[i]=0;

			sscanf(line,"%*f %*f %lf",&tmp); lat=(double *)realloc(lat,(i+1)*sizeof(double)); lat[i]=tmp; 
			sscanf(line,"%*f %*f %*f %lf",&tmp); lon=(double *)realloc(lon,(i+1)*sizeof(double)); lon[i]=tmp; 
			sscanf(line,"%*f %*f %*f %*f %lf",&tmp); Z=(double *)realloc(Z,(i+1)*sizeof(double)); Z[i]=tmp;
 			sscanf(line,"%*f %*f %*f %*f %*f %lf",&tmp); pi=(double *)realloc(pi,(i+1)*sizeof(double)); pi[i]=tmp;

			if(i)
				if(t[i]<=t[i-1]) t[i]=t[i-1]+epsilon;
			if(tmpt>t1 & lat[i]>=latmin & lat[i]<=latmax & lon[i]>=lonmin & lon[i]<=lonmax & Z[i]>=zmin & Z[i]<=zmax) {
				E=(int *)realloc(E,(j+1)*sizeof(int)); E[j]=i;
				ISE[i]=1;
				j++;
			}

			if(m[i]+m0>=mp & m[i]+m0<mp2) {
				P=(int *)realloc(P,(k+1)*sizeof(int)); P[k]=i;
				k++;
			}
			if(m[i]+m0>=mp2) {
				P2=(int *)realloc(P2,(q+1)*sizeof(int)); P2[q]=i;
				q++;
			}
			i++; 
		} else break;
	} 
	N=i; NE=j; NP=k; NP2=q;
	fclose(in);
	printf("* FOUND %d EARTHQUAKES with t<%lf IN %s, among which %d in the target period and volume\n",N,t2,data,NE);

	K2=(double *)calloc(NP2,sizeof(double));
	printf("* FOUND %d POTENTIAL PARENTS AND %d SPECIFIC MAINSHOCKS\n",NP,NP2);


	// INITIAL COMPUTATIONS
	X=(double *)calloc(N,sizeof(double));
	Y=(double *)calloc(N,sizeof(double));
	phi=(double *)calloc(NP,sizeof(double));
	phi2=(double *)calloc(NP2,sizeof(double));
	psi=(double *)calloc(NP,sizeof(double));
	psi2=(double *)calloc(NP2,sizeof(double));
	calculs_initiaux(t,m,lat,lon,pi,phi,phi2,psi,psi2,X,Y,Z,t1,t2,N,E,P,P2,NE,NP,NP2,p,c,gamma,L0,lonmin,lonmax,latmin,latmax,zmin,zmax);	// calcul des coordonnées kilométriques, et des coef phi, cf cahier labo 3/9/2021

	// INITIALIZING PROBABILITIES with a democratic guess (-> chaque parent potentiel a la même proba de déclencher que les autres)
	om=(double **)calloc(NP,sizeof(double *)); // om[i][j] = proba que l'enfant j (= séisme d'indice E[j]) soit déclenché par le parent i (= séisme d'indice P[i])
	for(i=0;i<NP;i++) om[i]=(double *)calloc(NE,sizeof(double));
		
	om2=(double **)calloc(NP2,sizeof(double *)); // comme pour om[i][j] mais en considérant cette fois les parents de la "liste spécifique" (les "très gros" séismes)
	for(i=0;i<NP2;i++) om2[i]=(double *)calloc(NE,sizeof(double));
		    
	// INITIAL CONDITIONS -> same K value for all mainshocks
	for(i=0;i<NP2;i++) K2[i]=1; K=1;


	// MAIN LOOP	
	niter=0; K_old=0; meanK_old=0;
	while(niter<MAX_ITERATIONS & dx>ACCURACY) {
		printf("ITERATION %d\n",niter+1);
		EXPECTATION(t,m,X,Y,Z,om,om2,P,P2,E,NP,NP2,NE,alpha,K,K2,p,c,L0,gamma); // étape E
		K=MAXIMIZATION(m,phi,phi2,psi,psi2,om,om2,E,P,P2,NE,NP,NP2,alpha,K2); // étape M
		if(NP==0) {
			meanK=0; for(i=0;i<NP2;i++) meanK+=K2[i]/NP2; 
			printf("mean specific K=%g\t",meanK);
			dx=fabs(meanK-meanK_old)/meanK;
		}
		else { printf("K=%g\t",K); dx=fabs(K-K_old)/K; }
		niter++;
		printf("CONVERGENCE %g\n\n",dx);
		K_old=K; meanK_old=meanK;
	}
	for(n=0;n<NP2;n++) printf("K du seisme %d = %g\n",P2[n]+1,K2[n]);
		
	// OUTPUT: (1) LAMBDA(t) cumulé et corrigé de la non-détection, de t1 à tous les temps des séismes E ; (2) valeur de K (=0 si pas parent potentiel) ; (3) nb de répliques directes
	sortie(output,t,m,pi,psi,psi2,om,om2,alpha,p,c,K,K2,P,P2,ISE,N,NP,NP2,t1,t2);	
}





	

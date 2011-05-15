/*
 ============================================================================
 Name        : isingspinn.c
 Author      : Espen Stene Broenstad 
 Description : Program for final exam in Computational Physics
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef short int sint;

const double RMAX = RAND_MAX;
/*
const int ITER = 1000000;
const double T_min = 0.1;
const double T_max = 4.5;
const double T_steps = 500;
*/

const double N_min = 2;
const double N_max = 32;



struct Hamilton {
	int k;
	int h;
};

struct Lattice {
	sint size;
	sint size2;
	char* grid;
	sint* n;
	sint* s;
	sint* et;
	sint* wt;
	sint* ek;
	sint* wk;

	sint* signe;
	sint* signw;

	int ham_t;
	int ham_k;
};

int rnd() {
	//Implement your own random function if you want
	return rand();
}

void allocLattice(struct Lattice*, sint);

void initLattice(struct Lattice*, sint);

void randomizeLattice(struct Lattice* latt);

void printLattice(struct Lattice*);

void calcHamiltons(struct Lattice* latt);

void calcHamiltons2(struct Lattice* latt);

double localHamiltonK(struct Lattice* latt, sint k);

int hamiltonContribK(struct Lattice* latt, sint k);

double localHamiltonT(struct Lattice* latt, sint k);

int hamiltonContribT(struct Lattice* latt, sint k);

double R(int delta_h,double T) {
	if (delta_h>0) {
		return exp(-delta_h/T);
	}
	else return 1;
}


//find Tau by Metropolis Monte Carlo method
double getTau(struct Lattice* latt, double T, int itermax);



int main(int argc, char** args) {
	if (argc <= 1) return 1;
	double T_c = 2/log(1+sqrt(2));
	int ITER = 1000000;
	double T_min = 0.1;
	double T_max = 4.5;
	double T_steps = 500;
	double dT = (T_max-T_min)/T_steps;
	double T;
	int N;
	int i = 0;
	double tau;
	int Ns[] = { 4,5,6,9,11,15,20, 26, 32};
	argc =2;
	struct Lattice* latt = malloc(sizeof(struct Lattice));
	allocLattice(latt,32);
	switch(args[1][0]) {
	case 'a':
		for (i = 0, N=Ns[i];  N!=0; N = Ns[++i]) {
			initLattice(latt, N);
			printf("N = %i\n",N);
			fprintf(stderr,"N: %i\n",N);
			int j = 0;
			for (T=T_min;T<T_max; T+=dT) {
				fprintf(stderr,"T: %e %i N: %i\n",T,j++,N);
				randomizeLattice(latt);
				getTau(latt,T,1e6);
				tau = getTau(latt,T,1e6);
				printf("%e\t%e\n",T,tau);
			}
			printf("\n");
		}
		break;
	case 'b':
		for (i = 0, N=Ns[i]; N<=20; N=Ns[++i]) {
			//fprintf(stderr,"T: %f\n",T);
			fprintf(stderr,"N: %i\n",N);
			initLattice(latt,N);
			int j;
			for (j = 0; j < 20; j++) {
				randomizeLattice(latt);
				getTau(latt, T_c, (int) ITER*exp(0.008*N*N));
				tau = getTau(latt, T_c, (int) ITER*exp(0.008*N*N));
				printf("%i\t%e\n",N,tau);
				//printf("%f\t%e\n",T,tau);
			}
		}
		break;
	}
	return EXIT_SUCCESS;
}

void allocLattice(struct Lattice* latt, sint size) {
	latt->size=size;
	sint size2 = size*size;
	latt->size2 = size2;
	latt->grid = malloc(sizeof(char)*size2);
	latt->n = malloc(sizeof(sint)*size2);
	latt->s = malloc(sizeof(sint)*size2);
	latt->et = malloc(sizeof(sint)*size2);
	latt->wt = malloc(sizeof(sint)*size2);
	latt->ek = malloc(sizeof(sint)*size2);
	latt->wk = malloc(sizeof(sint)*size2);
	latt->signe = malloc(sizeof(sint)*size2);
	latt->signw = malloc(sizeof(sint)*size2);
}

void initLattice(struct Lattice* latt, sint size) {
	latt->size=size;
	latt->size2=size*size;
	char* grid = latt->grid;
	sint* n = latt->n;
	sint* s = latt->s;
	sint* et = latt->et;
	sint* wt = latt->wt;
	sint* ek = latt->ek;
	sint* wk = latt->wk;
	sint* signe = latt->signe;
	sint* signw = latt->signw;
	sint i;
	sint j;
	sint index;
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			index = size*i + j;
			grid[index] = rnd()&1;
			n[index] = size*(i-1)+j;
			s[index] = size*(i+1)+j;
			et[index] = size*i+j+1;
			wt[index] = size*i+j-1;
			ek[index] = size*i+j+1;
			wk[index] = size*i+j-1;
			signe[index] = 1;
			signw[index] = 1;
		}
	}
	for (i = 0; i < size; i++) {
		n[i] = size*(size-1)+i;
		s[(size-1)*size+i] = i;
		et[size*(i+1)-1] = size*i;
		wt[size*i] = size*(i+1)-1;
		ek[size*(size-i)-1] = size*i;
		wk[size*i] = size*(size-i)-1;
		signw[size*i] = -1;
		signe[size*(size-i)-1] = -1;
	}
	latt->grid = grid;
	latt->n = n;
	latt->s = s;
	latt->et = et;
	latt->wt = wt;
	latt->ek = ek;
	latt->wk = wk;
	calcHamiltons(latt);
}

void randomizeLattice(struct Lattice* latt) {
	int i;
	for (i = 0; i < latt->size2; i++) {
		latt->grid[i] = rnd()&1;
	}
	calcHamiltons(latt);
}

void printLattice(struct Lattice* latt) {
	sint i;
	sint size = latt->size;
	sint size2 = size*size;
	char* arr = latt->grid;
	for (i = 0; i < size2; i++) {
		if (i%size==0) printf("\n");
		printf("%i\t", arr[i]);
	}
	printf("\n");
}

void calcHamiltons(struct Lattice* latt) {
	int ham_k = 0;
	int ham_t = 0;
	int temp;
	sint i = 0;
	sint size2 = latt->size*latt->size;
	for (i = 0; i < size2; i++) {
		temp = 1-2*(latt->grid[latt->s[i]]^latt->grid[i]);
		ham_k += temp
		 + (1-2*(latt->grid[latt->ek[i]]^latt->grid[i]))*latt->signe[i];
		ham_t += temp
		 + 1-2*(latt->grid[latt->et[i]]^latt->grid[i]);
	}
	latt->ham_k = -1*(ham_k);
	latt->ham_t = -1*(ham_t);
}

double localHamiltonK(struct Lattice* latt, sint k) {
	return 2 - 2*((latt->grid[latt->n[k]]^latt->grid[k])
	+ (latt->grid[latt->s[k]]^latt->grid[k]))
	+ (1-2*(latt->grid[latt->ek[k]]^latt->grid[k]))*latt->signe[k]
	+ (1-2*(latt->grid[latt->wk[k]]^latt->grid[k]))*latt->signw[k];
}

double localHamiltonT(struct Lattice* latt, sint k) {
	return 4 - 2*((latt->grid[latt->n[k]]^latt->grid[k])
	 + (latt->grid[latt->s[k]]^latt->grid[k])
	 + (latt->grid[latt->et[k]]^latt->grid[k])
	 + (latt->grid[latt->wt[k]]^latt->grid[k]));
}

double getTau(struct Lattice* latt, double T, int itermax) {
	double sum = 0;
	int i = 0;
	sint size2 = latt->size*latt->size;
	sint k;
	int ht_current;
	int ht_trial;
	int delta_ht;
	int hk_current;
	int hk_trial;

	double nextterm = exp((latt->ham_t-latt->ham_k)/T); //switched sign
	while (i++ < itermax) {
		k = rnd()%size2;
		ht_current = -1*localHamiltonT(latt,k);
		hk_current = -1*localHamiltonK(latt,k);
		latt->grid[k] ^= 1;
		ht_trial = -1*localHamiltonT(latt,k);
		hk_trial = -1*localHamiltonK(latt,k);
		delta_ht = ht_trial - ht_current;
		if (rnd()/RMAX <= R(delta_ht,T)) {
			latt->ham_t += delta_ht;
			latt->ham_k += hk_trial-hk_current;
			nextterm = exp((latt->ham_t-latt->ham_k)/T); //switched sign
		} else {
			latt->grid[k] ^= 1;
		}
		sum += nextterm;
	}
	return (-T/latt->size)*log(sum/i);
}

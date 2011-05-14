/*
 ============================================================================
 Name        : isingspinn.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef short int sint;

const double RMAX = RAND_MAX;

struct Hamilton {
	int k;
	int h;
};

struct Lattice {
	sint size;
	sint* grid;
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
	return rand();
}

void allocLattice(struct Lattice*, sint);

void initLattice(struct Lattice*, sint);

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



int main(void) {
	struct Lattice* latt = malloc(sizeof(struct Lattice));
	allocLattice(latt,32);

	//latt->grid[0]^=1;
	//printLattice(latt);

	double T;
	for (T = 0.1; T < 3; T += 0.1) {
		int i;
		fprintf(stderr,"T: %f\n",T);
		double tau=0;
		for (i = 0; i < 10; i++) {
			initLattice(latt,8);
			calcHamiltons2(latt);
			getTau(latt, T, 10000);
			tau = getTau(latt, T, 10000);
			printf("%f\t%e\n",T,tau);
		}
	}
	return EXIT_SUCCESS;
}

void allocLattice(struct Lattice* latt, sint size) {
	latt->size=size;
	sint size2 = size*size;
	latt->grid = malloc(sizeof(sint)*size2);
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
	sint* grid = latt->grid;
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
}

void printLattice(struct Lattice* latt) {
	sint i;
	sint size = latt->size;
	sint size2 = size*size;
	sint* arr = latt->grid;
	for (i = 0; i < size2; i++) {
		if (i%size==0) printf("\n");
		printf("%i\t", arr[i]);
	}
}

void calcHamiltons(struct Lattice* latt) {
	int ham_k = 0;
	int ham_t = 0;
	int temp;
	sint i = 0;
	sint size2 = latt->size*latt->size;
	for (i = 0; i < size2; i++) {
		ham_k += 1-2*(latt->grid[latt->n[i]]^latt->grid[i]);
		ham_k += 1-2*(latt->grid[latt->s[i]]^latt->grid[i]);
		ham_k += (1-2*(latt->grid[latt->ek[i]]^latt->grid[i]))*latt->signe[i];
		ham_k += (1-2*(latt->grid[latt->wk[i]]^latt->grid[i]))*latt->signw[i];

		ham_t += 1-2*(latt->grid[latt->n[i]]^latt->grid[i]);
		ham_t += 1-2*(latt->grid[latt->s[i]]^latt->grid[i]);
		ham_t += 1-2*(latt->grid[latt->et[i]]^latt->grid[i]);
		ham_t += 1-2*(latt->grid[latt->wt[i]]^latt->grid[i]);
	}
	if (ham_k&1) puts("dsfsdgf K");
		if (ham_t&1) puts("dsafasf T");
	latt->ham_k = -1*(ham_k/2);
	latt->ham_t = -1*(ham_t/2);
}

void calcHamiltons2(struct Lattice* latt) {
	int ham_k = 0;
	int ham_t = 0;
	int temp;
	sint i = 0;
	sint size2 = latt->size*latt->size;
	for (i = 0; i < size2; i++) {
		temp = latt->grid[latt->n[i]] + latt->grid[latt->s[i]];
		ham_k += hamiltonContribK(latt,i);
		ham_t += hamiltonContribT(latt,i);
	}
	if (ham_k&1) puts("dsfsdgf K");
	if (ham_t&1) puts("dsafasf T");
	latt->ham_k = -1*(ham_k/2);
	latt->ham_t = -1*(ham_t/2);
}


double localHamiltonK(struct Lattice* latt, sint k) {
	return (
			hamiltonContribK(latt,k)          +
			hamiltonContribK(latt,latt->n[k])  +
			hamiltonContribK(latt,latt->s[k])  +
			hamiltonContribK(latt,latt->ek[k]) +
			hamiltonContribK(latt,latt->wk[k])
			)*-0.5;
}

int hamiltonContribK(struct Lattice* latt, sint k) {
	int res = 1 - 2*(latt->grid[latt->n[k]]^latt->grid[k]);
	res += 1 - 2*(latt->grid[latt->s[k]]^latt->grid[k]);
	res += (1 - 2*(latt->grid[latt->ek[k]]^latt->grid[k]))*latt->signe[k];
	res += (1 - 2*(latt->grid[latt->wk[k]]^latt->grid[k]))*latt->signw[k];
	return res;
}

double localHamiltonT(struct Lattice* latt, sint k) {
	return (
			hamiltonContribT(latt,k)          +
			hamiltonContribT(latt,latt->n[k])  +
			hamiltonContribT(latt,latt->s[k])  +
			hamiltonContribT(latt,latt->et[k]) +
			hamiltonContribT(latt,latt->wt[k])
			)*-0.5;
}

int hamiltonContribT(struct Lattice* latt, sint k) {
	int res = 1 - 2*(latt->grid[latt->n[k]]^latt->grid[k]);
	res += 1 - 2*(latt->grid[latt->s[k]]^latt->grid[k]);
	res += 1 - 2*(latt->grid[latt->et[k]]^latt->grid[k]);
	res += 1 - 2*(latt->grid[latt->wt[k]]^latt->grid[k]);
	return res;
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
		ht_current = localHamiltonT(latt,k);
		hk_current = localHamiltonK(latt,k);
		latt->grid[k] ^= 1;
		ht_trial = localHamiltonT(latt,k);
		hk_trial = localHamiltonK(latt,k);
		delta_ht = ht_trial - ht_current;
		double p = rnd()/RMAX;
		//printf("p: %e\n",p);
		//printf("R: %e\n",R(delta_ht,T));
		if (p <= R(delta_ht,T)) {
			//puts("kept");
			latt->ham_t += delta_ht;
			latt->ham_k += hk_trial-hk_current;
			nextterm = exp((latt->ham_t-latt->ham_k)/T); //switched sign
		} else {
			//rollback
			//puts("rollback");
			latt->grid[k] ^= 1;
		}
		sum += nextterm;
	}

	//return -(T/latt->size)*log(sum/i);
	return (-T/latt->size)*log(sum/i);
}

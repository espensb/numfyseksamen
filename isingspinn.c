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

typedef short int sint;


struct lattice {
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


void allocLattice(struct lattice*, sint);

void initLattice(struct lattice*, sint);

void printLattice(struct lattice*);

void calcHamiltons(struct lattice* latt);

void localHamilton(struct lattice* latt, sint i);

void hamiltonContrib(struct lattice latt, sint i);



int main(void) {
	puts("!!!Hello World!!!"); /* prints !!!Hello World!!! */
	struct lattice* latt = malloc(sizeof(struct lattice));
	allocLattice(latt,5);
	initLattice(latt,1);
	printLattice(latt);
	calcHamiltons(latt);
	printf("\nHamilton: %i\n", latt->ham_t);
	return EXIT_SUCCESS;
}

void allocLattice(struct lattice* latt, sint size) {
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

void initLattice(struct lattice* latt, sint size) {
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
			grid[index] = rand()&1;
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

void printLattice(struct lattice* latt) {
	sint i;
	sint size = latt->size;
	sint size2 = size*size;
	sint* arr = latt->grid;
	for (i = 0; i < size2; i++) {
		if (i%size==0) printf("\n");
		printf("%i\t", arr[i]);
	}
}

void calcHamiltons(struct lattice* latt) {
	int ham_k = 0;
	int ham_t = 0;
	int temp;
	sint i = 0;
	sint size2 = latt->size*latt->size;
	for (i = 0; i < size2; i++) {
		temp = latt->grid[latt->n[i]] + latt->grid[latt->s[i]];
		ham_k += (
					temp +
					latt->grid[latt->ek[i]]*latt->signe[i] +
					latt->grid[latt->wk[i]]*latt->signw[i]
				) * latt->grid[i];
		ham_t += (
					temp +
					latt->grid[latt->et[i]] +
					latt->grid[latt->wt[i]]
				) * latt->grid[i];
	}
	latt->ham_k = -1*(ham_k >> 1);
	latt->ham_t = -1*(ham_t >> 1);
}




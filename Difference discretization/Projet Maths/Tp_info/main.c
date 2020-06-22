#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 100
#define m 3

void vecteur_xi(float xi[N], float a, float b) {
//On calcule norte vecteur réference en foncion des données du pb

	int i;
	for (i=0; i<=N; i++) {
		xi[i]=i*(b-a)/N;
	}

}


void y_exact_fonc(float xi[N],float y_exact[N]) {
//On calcule le resultat de y_exact en fonction de x

	int j;
	for (j=0; j<=N; j++) {
		y_exact[j]=exp(-10*xi[j]);
	}
}

void afficher_vect( float a[N]) {
	int i;
	for (i=0;i<N;i++) printf("%1f\n",a[i]);
}


int main() {

//paramètre matrice
	float A[m][N];
	float T[N][m];
	float b[m];

//paramètre intervalle
	int a=0;
	int b=3;

//paramètre matrice
	float y_exact[N];
	float xi[N];

// calcule de xi
	vecteur_xi(xi, a, b);
	afficher_vect(xi);

//
	y_exact_fonc(xi,y_exact);
	afficher_vect(y_exact);
//





return 0;
}


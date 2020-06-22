void vecteur_xi(int N, float a, float b) {
//On calcule norte vecteur réference en foncion des données du pb
	float xi[N];
	int i;
	for (i=a; i<=N; i++) {
		xi[i]=i*(b-a)/N;
	}

}


void y_exact_fonc(float xi[], int N) {
//On calcule le resultat de y_exact en fonction de x
	float y_exact[N];
	int j;
	for (j=a; j<=N; j++) {
		y_exact[i]=exp(-10*xi[i]);
	}
}


void mat_A=(int m, int N, float xi[]) {
	float A[m][N];
	int i;
	int j;
// On calcule la matrice A, qui represente notre système
	for (i=0;i<m;i++) {
		for (j=0;j<m;j++) {
			A[i][j]=(xi[j]^i);
		}
	}
}

//On remarque que b=y_exacte
//Nous avons tr(M)Mx=tr(M)y_exact

void transpose(int m, int N, float M[][]) {
//On calcule la transposé d'une matrice
	int i;
	int j;
	int T[N][m];
	for(i = 0; i < m; i++) {
        	for(j = 0; j < N; j++) {
          
                 	T[j][i] = M[i][j];
         	 }
     	}
}

//on veut calculer les matrices pour les moindres carrés

float mat_M(int m, int N, float M[][], float A[][]) {
//Le but de se prog est de trouver M=tA*A et de généraliser pour chaque matrice	
	int i;
	int j;
// on appelle le prog transposé pour calculer M
//On initialise M
	float M[m][m];
	for (i=0; i<N;i++) {	
		for (j=0; j<m; j++) {
			M[i][j]=0;
			for (k=0; k<N; k++) {
				M[i][j]=T[i][k]*A[i][k]+M[i][j];
			}
		}		
	}
	return M;
}

//il faut avoir calculer y_exact au préalable et la transposée de A ( faire d'abord premiere partie, sinon bug)
void mat_b(int m, int N, float T[][], float y_exact[]) {
	float b[m];
	int i;
	for (i=0; i<m; i++) {
		b[i]=0;
		for (k=0; k<N; k++) {
			b[i]=T[i][k]*y_exact[k]+b[i];
		}
	}
}
	
		

//methode directe, on cherche LU

void factoriser_LU( float M[m][m], float L[m][m], float U[m][m], int m) {
	int i,j,k;
	float l;
	for (j=0;j<=m-1;j++) {	
		for (i=0;i<=j;i++) {
			l=0;
			for (k=0;k<=i-1;k++) l += L[i][k]*U[k][j];
			U[i][j] = A[i][j] - l;
		}
		L[j][j] = 1;
		for (i=j+1;i<=n-1;i++) {
			l=0;
			for (k=0;k<=j-1;k++) l += L[i][k]*U[k][j];
			L[i][j] = (A[i][j] - l)/U[j][j];
		}			
	}	
}

void resol_trig_inf( float M[m][m] , float x[MAX], float b[m], int m) {
//On a x l'inconnu
	int i,k;
	float l;
	
	for (i=0;i<m;i++) {
		l=0;
		for (k=0;k<=i-1;k++) l += A[i][k]*x[k];
		x[i] = (b[i] -l)/A[i][i];
	}
} 


void resol_trig_sup( float M[m][m] , float x[m], float b[m], int m) {
	int i,k;
	float l;
	
	for (i=m;i>=0;i--) {
		l=0;
		for (k=i+1;k<m;k++) l += A[i][k]*x[k];
		x[i] = (b[i] -l)/A[i][i];
	}
}	



	
	
			

	









			
	
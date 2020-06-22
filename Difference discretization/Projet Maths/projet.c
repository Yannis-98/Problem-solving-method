
typedef struct Matrice Matrice
struct Matrice
{
/* Ici, on cree une structure pour manipuler une matrice */ 
	int NbrLig; /* Nombre de ligne de la matrice */
	int NbrCol; /* Nombre de colonne de la matrice */
	double **coef; /* Tableau dynamique representant les coeficient de la matrice */
}

typedef struct Vecteur Vecteur
struct Vecteur
{
	int NbrLig;
	double *coef;
}

Vecteur InitVecteurZero(int Ni)
{
/* Cette fonction permet d'allouer le tableau dynamique et d'initialiser tous les coefficients de la matrice */
	int i;
	Vecteur tmp;

	tmp.NbrLig=Ni;
	tmp.coef = malloc(sizeof(double)*N);

	for(i=0;i<Ni;i++)
	{
		tmp.coef[i]=0.0;
	}

	return(tmp);
}

Matrice InitMatriceZero(int Ni, int Nj)
{
/* Cette fonction permet d'allouer le tableau dynamique et d'initialiser tous les coefficients de la matrice */
	int i;
	int j;
	Matrice tmp;

	tmp.NbrLig=Ni;
	tmp.NbrCol=Nj;

/* Allocation memoire */
	tmp.coef = malloc(sizeof(double)*Ni);
	for(i=0;i<Ni;i++)
	{
		tmp.coef[i]=malloc(sizeof(double)*Nj);
	}

/* Mise a zero des coefficients
	for(i=0;i<Ni;i++)
	{
		for(j=0;j<Nj;j++)
		{
			tmp.coef[i][j]=0.0;
		}
	}

	return(tmp);
}


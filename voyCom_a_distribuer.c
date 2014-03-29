#include<stdio.h> 
#include<math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>

#define MAXMOT 				256
#define MAXS 				500
#define EXACT_NAIF_NB_VILLES 		10
#define EXACT_BRANCH_BOUND_NB_VILLES 	10
#define APPROCHE_PPV_NB_VILLES		20
#define CALCUL_ACM			30

#define SHELLSCRIPT "\
#/bin/bash \n\
echo \"Arbre couvrant minimal:\" \n\
dot -Tpng -o tree.png arbreCouvrantMin.dot \n\
eog tree.png \n\
"

/**
 * Structure pour representer un cycle
 */
typedef struct t_cycle 
{
  int taille;   //la taille du cycle en construction
  double poids; //le coût du cycle
  int c[MAXS];  //liste des "taille" sommets
} t_cycle;     


/**
 * Charge le CSV des coordonnées des villes.
 * 
 * @param [in] f le fichier
 * @param [out] nb_villes le nombre de villes de l'instance
 * @param [out] dist le tableau des nb_villes*nb_villes distances
 * @param [out] absc le tableau des abscisses des villes
 * @param [out] ord le tableau des ordonnées des villes
 */
void lire_donnees(const char *f, unsigned int *nb_villes, double ***dist,  double **absc, double **ord)
{
  //double *absc; ///tableau des ordonnées
  //double *ord;  /// tableau des abscisses
  char ligne[MAXMOT];
  FILE * fin = fopen(f,"r");

  if(fin != NULL)
    {
      //On recupere le nombre de villes
      fgets(ligne, MAXMOT, fin);
      *nb_villes = atoi(ligne);
      (*dist) = (double**)malloc(*nb_villes * sizeof(double*));
      (*absc) = (double*)malloc(*nb_villes * sizeof(double));
      (*ord) = (double*)malloc(*nb_villes * sizeof(double));
      int i = 0;
      while (fgets(ligne, MAXMOT, fin) != NULL) 
	{      
	  char *p = strchr(ligne, ';');
	  ligne[strlen(ligne) - strlen(p)]='\0';
	  p = &p[1];
	  (*absc)[i] = atof(ligne);
	  (*ord)[i] = atof(p);
	  i = i + 1;
	}
    }
  else
    {
      printf("Erreur de lecture du fichier.\n");
      exit(2);
    }
  fclose(fin);
  int i,j;

  //Calclul des distances
  for(i = 0; i < *nb_villes; i++)
    {
      (*dist)[i] = (double*) malloc(*nb_villes * sizeof(double));
      for(j = 0; j < *nb_villes; j++)
	{
	  (*dist)[i][j] = sqrt( ((*absc)[i] - (*absc)[j])* ((*absc)[i] - (*absc)[j]) + ((*ord)[i] - (*ord)[j]) * ((*ord)[i] - (*ord)[j]) );
	}
    }
}


/**
 * Supprime la structure des distances
 *
 * @param [in] nb_villes le nombre de villes.
 * @param [in,out] distances le tableau à supprimer.
 * @param [in,out] abscisses un autre tableau à supprimer.
 * @param [in,out] ordonnees encore un autre tableau à supprimer.
 */
void supprimer_distances_et_coordonnees(const int nb_villes, double **distances, double *abscisses, double *ordonnees)
{
  int i;
  for(i = 0; i < nb_villes; i++)
    {
      free(distances[i]);
    }
 free(distances);
 free(abscisses);
 free(ordonnees);
}


/**
 * Export le cycle dans un fichier HTML pour pouvoir être visualisé
 * dans l'applet.
 *
 * @param [in] cycle le cycle à afficher
 */
void afficher_cycle_html(const t_cycle cycle, double *posX, double *posY)
{
  FILE * fout = fopen("DisplayTsp2.html","w");
  if(fout != NULL)
    {
      int i;
      fprintf(fout, "<html>\n <applet codebase=\".\" code=\"DisplayTsp.class\" width=200 height=200>\n");
      fprintf(fout, "<param name = Problem value = \"custom\">\n");
      fprintf(fout, "<param name = Problem CitiesPosX value = \"");
      for(i = 0; i < cycle.taille; i++)
	fprintf(fout,"%f;",posX[i]);
      fprintf(fout, "\">\n");
      fprintf(fout, "<param name = Problem CitiesPosY value = \"");
      for(i = 0; i < cycle.taille; i++)
	fprintf(fout,"%f;",posY[i]);
      fprintf(fout, "\">\n");
      fprintf(fout, "<param name = Parcours value = \"");
      fprintf(fout,"%d",cycle.c[0]);
      for(i = 1; i < cycle.taille; i++)
	fprintf(fout,"-%d",cycle.c[i]);
      fprintf(fout,"\">\n</applet>\n </html>\n");
    }
  fclose(fout);
}

/**
 * Affiche le tableau des distances.
 *
 * @param [in] nb le nombre de villes
 * @param [in] distances le tableau
 */
void afficher_distances(const int nb, double **distances)
{
  unsigned int i ;
  unsigned int j ;
  for(i = 0  ; i < nb; i++)
    {   
      for(j = 0 ; j < nb ; j++)
	printf("%f ", distances[i][j]);
      printf("\n");
    }
  printf("\n");
}


/**
 * Fonction de comparaison pour le trie des arêtes par leur poids.
 *
 * @param [in] v1 pointeur vers un triplet (i,j,poids)
 * @param [in] v2 pointeur verts un triplet (i,j,poids)
 * @return vrai si poid v1 < poids v2
 */
int comparer(const void *v1, const void *v2)
{
  double **px1 = (double **) v1;
  double **px2 = (double **) v2;

  double *x1 = *px1;
  double *x2 = *px2;
  if(x1[2] - x2[2] < 0)
    return -1;
  else 
    {
      if(x1[2] - x2[2] == 0)
	return 0;
      else
	return 1;
    }
}

/**
 * Construit un tableau de n*(n-1)/2 arêtes trié selon leur poids.
 *
 * @note utile pour le Kruskal
 *
 * @param [in] n le nombre de villes
 * @param [in] d tableau des n x n distances.
 * @return tableau d'arétes triées T[i][j] = poids_ij 
 */
double **trier_aretes(const int n, double **d)
{
  assert(d);
  
  int nb_aretes = n * (n - 1) / 2;
  double **T = (double **)malloc(nb_aretes * sizeof(double *));
  int i, j;
  int a = 0;

  //On initialise la structure d'arêtes
  for(i = 0; i < n-1; i++)
    {
      for(j = i+1; j < n; j++)
	{
	  T[a] = (double *)malloc(3 * sizeof(double));
	  T[a][0] = i;
	  T[a][1] = j;
	  T[a][2] = d[i][j];
	  a++;
	}
    }
  
  //Appel au quicksort avec la bonne fonction de comparaison
  qsort(T, a, sizeof(T[0]), comparer);

  

  //Decommenter pour vérifier le tri
  /*
  for(i = 0; i < a; i++)
    printf("%f ", T[i][2]);
  printf("\n");
  */
  return T;
}


/**
 * Supprime le tableau des arêtes.
 *
 * @param [in] nb_villes le nombre de villes
 * @param [in,out] T le tableau à supprimer
 */
void supprimer_aretes(const int nb_villes, double **T)
{
  assert(T);

  int nb_aretes = nb_villes*(nb_villes - 1 ) / 2;
  unsigned int i;
  for( i = 0; i < nb_aretes ; ++i)
    free( T[i] );

  free(T);
}

// TODO: a commenter
int est_dans_chemin(t_cycle * chemin, int ville)
{
	int j = 0;
	for( j = 0; j < chemin->taille; j++)
	{
		if( chemin->c[j] == ville )
		{
			return 1;
		}
	}

	return 0;
}

// TODO: a commenter
void recopie_chemin(t_cycle * cheminRef, t_cycle * cheminOP){
	cheminOP->taille = cheminRef->taille;
	cheminOP->poids = cheminRef->poids;
	int i;
	for(i=0; i< cheminRef->taille; i++){
		//printf("%d\n",cheminRef->c[i]);
		cheminOP->c[i] = cheminRef->c[i];
	}
}

/**
 * Calcule le chemin le plus court exact de manière naïve
 *
 * @param [nbVille] 
 * @param [distances] 
 */
void PVC_exact_naif(int nbVille,double **distance, t_cycle *  chemin, t_cycle * meilleur)
{
	if( chemin->taille == nbVille)
	{	//printf("if: %d\n",chemin->taille);
		chemin->poids += distance[chemin->c[0]][chemin->c[chemin->taille-1]];
		if( chemin->poids < meilleur->poids)
		{	//printf("if2: %d\n",chemin->taille);
			recopie_chemin(chemin,meilleur);
		}
		chemin->poids -= distance[chemin->c[0]][chemin->c[chemin->taille-1]];
	}
	else
	{	//printf("else: %d\n",chemin->taille);
		int i = 0;
		for(i = 0; i < nbVille; i++)
		{

			if(!est_dans_chemin(chemin,i))
			{
				chemin->poids += distance[chemin->c[chemin->taille-1]][i];
				chemin->taille++;
				chemin->c[chemin->taille-1] = i;
				PVC_exact_naif(nbVille,distance, chemin,meilleur);
				chemin->taille --;
				chemin->poids -= distance[chemin->c[chemin->taille-1]][i];
			}
		}
	}
}

/**
 * Calcule le chemin le plus court exact de manière moins naïve
 *
 * @param [nbVille] 
 * @param [distances] 
 */
void PVC_exact_branch_and_bound(int nbVille,double **distance, t_cycle *  chemin, t_cycle * meilleur)
{ 
        if( chemin->taille == nbVille)
        {       //printf("if: %d\n",chemin->taille);
                chemin->poids += distance[chemin->c[0]][chemin->c[chemin->taille-1]];
                if( chemin->poids < meilleur->poids)
                {       //printf("if2: %d\n",chemin->taille);
                       recopie_chemin(chemin,meilleur);
                 	chemin->poids -= distance[chemin->c[0]][chemin->c[chemin->taille-1]];
		}
        }
        else
        {       //printf("else: %d\n",chemin->taille);
		if(!(chemin->poids > meilleur->poids))
		{
		
                	int i = 0;
                	for(i = 0; i < nbVille; i++)
                	{
                        	if(!est_dans_chemin(chemin,i))
                        	{
                                	chemin->poids += distance[chemin->c[chemin->taille-1]][i];
                                	chemin->taille++;
                                	chemin->c[chemin->taille-1] = i;
                                	PVC_exact_naif(nbVille,distance, chemin,meilleur);
                                	chemin->taille --;
                                	chemin->poids -= distance[chemin->c[chemin->taille-1]][i];
                        	}
                	}
		}
        }
}

/**
 * Approxime le chemin le plus court par plus proche voisin
 *
 * @param [nbVille] 
 * @param [distances] 
 */
void PVC_approche_ppv(int nbVille,double **distance, t_cycle *  chemin)
{
	int j;
	for(j=0;j<nbVille;j++)
	{
		int i;
		double distanceMinimale = 64000.0;
		int villePlusProche;
		for(i=0;i<nbVille;i++)
		{
			if(distance[chemin->c[chemin->taille-1]][i]<distanceMinimale && !est_dans_chemin(chemin,i))
			{
				distanceMinimale = distance[chemin->c[chemin->taille-1]][i];
				villePlusProche = i;
			}
		}
		chemin->c[chemin->taille] = villePlusProche;
		chemin->taille++;
		chemin->poids+=distanceMinimale;
	}

	chemin->poids += distance[chemin->c[0]][chemin->c[chemin->taille-1]];
}

// TODO à commenter
int rechercher_racine(int u , int* pi){
	int v = 0;
	while(v != -1 )
	{
		v= pi[u];
		if(v==-1)
		{
			//printf("Rech racine\n");
			return u;
		}
		u = v;
	}
	return u;
}

// TODO à commenter
double ** calcul_ACM(int nbVille,double **distance)
{
	int pi [nbVille] ;// parent
	int hauteur [nbVille] ;

	double ** T =(double **) calloc((nbVille-1)*nbVille/2,sizeof(double *));
	double ** K =(double **) calloc((nbVille-1)*nbVille/2,sizeof(double *));
	int i;
	int pos=0;
	for(i=0; i< (nbVille-1)*nbVille/2; i++)
	{
		T[i]= (double *) calloc(3,sizeof(double));
		K[i]= (double *) calloc(2,sizeof(double));
	}

	for(i=0 ; i<nbVille; i++) // initialisation des tableaux à 0
	{
		pi[i]= -1;
		hauteur[i]=0;
	}
	T = trier_aretes(nbVille, distance);

	//for(i=0 ; i<(nbVille-1)*nbVille/2; i++) 
	//{
	//	printf("T: %f-%f\n",T[i][0],T[i][1]);
	//}

	int s,t,r1,r2;
	i=0;

	while(i< (nbVille-1)*nbVille/2-1)
	{
		s = T[i][0];	
		t = T[i][1];
		r1 = rechercher_racine(s, pi);
		r2 = rechercher_racine(t, pi);
		//printf("s: %d\n",s);			
		//printf("t: %d\n",t);			

		if(r1!=r2)
		{
			K[pos][0] = s;
			K[pos][1] = t;
			
			if(hauteur[r1]> hauteur[r2])
			{
				pi[r1] = r2;
			}
			else
			{
				pi[r2] = r1;

				if(hauteur[r1] == hauteur[r2])
				{
					hauteur[r2] ++;
				}
		
			}
			pos++;
			//printf("pos: %d\n",pos);	
		}
		i++;
		//printf("i: %d\n",i);	
	}
	supprimer_aretes(nbVille, T); // On free les arrêtes créées
	return K;
}

/**
 * Fonction main.
 */
int main (int argc, char *argv[])
{
  double **distances;
  double *abscisses;
  double *ordonnees;
  unsigned int nb_villes;
  
  //Initialisation du timer pour mesurer des temps (compiler avec -lrt) 
  struct timespec myTimerStart;
  clock_gettime(CLOCK_REALTIME, &myTimerStart);

  //Exemple de mesure du temps
  lire_donnees("defi250.csv", &nb_villes, &distances, &abscisses, &ordonnees);

  //Récupération du timer et affichage
  struct timespec current;
  clock_gettime(CLOCK_REALTIME, &current); //Linux gettime
  double elapsed_in_ms =    (( current.tv_sec - myTimerStart.tv_sec) *1000 +
          ( current.tv_nsec - myTimerStart.tv_nsec)/1000000.0);
  printf("Temps passe (ms) : %lf\n", elapsed_in_ms);
  
  t_cycle chemin;
  chemin.taille = 1;
  chemin.poids =0;
  chemin.c[0] = 1;

  t_cycle  meilleur;
  meilleur.taille = 1;
  meilleur.poids =64000;
  meilleur.c[0] = 1;

  int i = 0;

  // ------------------------------------ Exact naif

  printf("\033[0;35m%s\033[0m\n","Exact naif:");

  for(i= 1; i<= EXACT_NAIF_NB_VILLES;i++)
  {
	chemin.taille = 1;
	chemin.poids =0;
	chemin.c[0] = 1;
	meilleur.taille = 1;
	meilleur.poids =64000;
	meilleur.c[0] = 1;
	clock_gettime(CLOCK_REALTIME, &myTimerStart);
	PVC_exact_naif(i,distances, &chemin,&meilleur);
	clock_gettime(CLOCK_REALTIME, &current); //Linux gettime
	elapsed_in_ms =    (( current.tv_sec - myTimerStart.tv_sec) *1000 +
		  ( current.tv_nsec - myTimerStart.tv_nsec)/1000000.0);
	printf("[Exact naif, %d villes] Temps passe (ms) : %lf\n",i, elapsed_in_ms);
  }

  // ---------------------------------- Exact B & B

  printf("\033[0;35m%s\033[0m\n","Exact branch and bound:");

  for(i= 1; i<= EXACT_BRANCH_BOUND_NB_VILLES ;i++)
  {
	chemin.taille = 1;
	chemin.poids =0;
	chemin.c[0] = 1;
	meilleur.taille = 1;
	meilleur.poids =64000;
	meilleur.c[0] = 1;
	clock_gettime(CLOCK_REALTIME, &myTimerStart);
	PVC_exact_branch_and_bound(i,distances, &chemin,&meilleur);
	clock_gettime(CLOCK_REALTIME, &current); //Linux gettime
	elapsed_in_ms =    (( current.tv_sec - myTimerStart.tv_sec) *1000 +
		  ( current.tv_nsec - myTimerStart.tv_nsec)/1000000.0);
	printf("[Exact branch and bound, %d villes] Temps passe (ms) : %lf\n",i, elapsed_in_ms);
  } 

  // --------------------------------- Approche PPV

  printf("\033[0;35m%s\033[0m\n","Approche ppv:");
  
  for(i= 1; i<= APPROCHE_PPV_NB_VILLES ;i++)
  {
	chemin.taille = 0;
	chemin.poids =0;
	chemin.c[0] = 1;
	clock_gettime(CLOCK_REALTIME, &myTimerStart);
	PVC_approche_ppv(i,distances, &chemin);
	clock_gettime(CLOCK_REALTIME, &current); //Linux gettime
	elapsed_in_ms =    (( current.tv_sec - myTimerStart.tv_sec) *1000 +
		  ( current.tv_nsec - myTimerStart.tv_nsec)/1000000.0);
	printf("[Approche PPV, %d villes] Temps passe (ms) : %lf\n",i, elapsed_in_ms);
  }

  // ----------------------------------------- ACM
  
  printf("\033[0;35m%s\033[0m\n","Calcul ACM:");
  
  double ** K;	
  int nbVilles = CALCUL_ACM;
  FILE * arbreCouvrantMin = fopen ("arbreCouvrantMin.dot","w+");

  fprintf(arbreCouvrantMin,"graph Tree\n{\n");
 
  K = calcul_ACM(nbVilles, distances);

  printf("Arbre couvrant minimal avec %d villes (%d a %d)\n",nbVilles,0,nbVilles-1);
  //for(i=0; i< nbVilles-1; i++)
  //{
  //	fprintf(arbreCouvrantMin,"\t%d [pos = \"%f,%f!\"];\n",i, abscisses[i],1000*ordonnees[i]);
  //}
  for(i=0; i< nbVilles-1; i++)
  {
	printf("k: %f-%f\n", K[i][0],K[i][1]);
  	fprintf(arbreCouvrantMin,"\t%d--%d;\n", (int)K[i][0],(int)K[i][1]);
  }

  fprintf(arbreCouvrantMin,"}");
  fclose(arbreCouvrantMin);

  // ---------------------------------- Affichage
  afficher_cycle_html(chemin, abscisses, ordonnees);
  system(SHELLSCRIPT);

  // ------------- Libérations mem et terminaison  
  supprimer_distances_et_coordonnees(nb_villes, distances, abscisses, ordonnees);
  return 0;

  // ----------------------------- Stuff commenté
  //int i = 0;  
  //for(i = 0;i<5; i++)
  //{	
  //      printf("%d\n",meilleur.c[i]);
  //}
  //printf("%d\n",meilleur.taille);
  //printf("%f\n",meilleur.poids);
  //
  //Affichage des distances
  //afficher_distances(nb_villes,distances);

  ////naif
  //t_cycle cycle;
  //cycle.taille=3;
  //cycle.c[0]=0;
  //cycle.c[1]=1;
  //cycle.c[2]=2;
  //afficher_cycle_html(cycle, abscisses, ordonnees);
  
  //double ** Aretes =  trier_aretes(nb_villes, distances);
  /// <-- Kruskal Here
  //supprimer_aretes(nb_villes, Aretes); XXX effectué directement dans calcul_ACM ...

}

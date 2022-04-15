#ifndef _VEC_3_
#define _VEC_3_

typedef struct Vec3 Vec3;

struct Vec3
{
	double x;
	double y;
	double z;
};

/*
 * vec3_create : Initialise un Vec3 à
 * 	partir des coordonnées données en
 * 	entrée
 *
 * x, y et z : Coordonnées du vecteur
 */
Vec3 vec3_create(double x, double y, double z);

/*
 * vec3_add : Calcule la somme de deux vecteurs
 * 	et renvoie le résultat.
 *
 * v1 et v2: Les termes de l'addition
 */
Vec3 vec3_add(Vec3 v1, Vec3 v2);

/*
 * vec3_diff : Calcule la différence de deux vecteurs
 * 	et renvoie le résultat
 *
 * v1 et v2: Les termes de la soustraction
 */
Vec3 vec3_diff(Vec3 v1, Vec3 v2);

/*
 * vec3_dot : Calcule le produit scalaire
 * 	de deux vecteurs et renvoie le résultat.
 *
 * v1 et v2: Les termes du produit
 */
double vec3_dot(Vec3 v1, Vec3 v2);

/*
 * vec3_adot : Calcule la valeur absolue du produit scalaire 
 * 	de deux vecteurs et renvoie le résultat.
 *
 * v1 et v2: Les termes du produit
 */
double vec3_adot(Vec3 v1, Vec3 v2);

/*
 * vec3_cross : Renvoie le produit vectoriel de
 * 	de deux vecteurs.
 *
 * v1 et v2: Les termes du produit
 */
Vec3 vec3_cross(Vec3 v1, Vec3 v2);

/*
 * vec3_smul : renvoie le produit d'un réel
 * 	avec un vecteur.
 *
 * v: Le vecteur du produit
 * a: Le réel du produit
 */
Vec3 vec3_smul(double a, Vec3 v);

/*
 * vec3_norm_2 : renvoie la norme au carré d'un vecteur.
 *
 * v: Le vecteur dont on va calculer la norme.
 */
double vec3_norm_2(Vec3 v);

/*
 * vec3_norm : renvoie la norme d'un vecteur.
 *
 * v: Le vecteur dont on va calculer la norme.
 */
double vec3_norm(Vec3 v);

/*
 * vec3_normalize : normalise un vecteur et renvoie
 * 	le résultat.
 *	Renvoie le vecteur (1, 0, 0) si le vecteur
 *	en entrée est le vecteur nul.
 *
 * v: Le vecteur que l'on va normaliser.
 */
Vec3 vec3_normalize(Vec3 v);

/*
 * vec3_build_coordonate_system : Construit une base à
 *	partir d'un vecteur
 *
 * v1: Le vecteur qui nous sert de base.
 * v2 et v3: Les deux autres vecteurs
 */
void vec3_build_coordonate_system(Vec3 v1, Vec3 *v2, Vec3 *v3);

/*
 * vec3_lerp : 
 *
 * v1 et v2: Le vecteur qui nous sert de base.
 * t: 
 */
Vec3 vec3_lerp(Vec3 v1, Vec3 v2, double t);

/*
 * vec3_min : 
 *
 * v1 et v2: Le vecteur qui nous sert de base.
 * t: 
 */
Vec3 vec3_min(Vec3 v1, Vec3 v2);

/*
 * vec3_max : 
 *
 * v1 et v2: Le vecteur qui nous sert de base.
 * t: 
 */
Vec3 vec3_max(Vec3 v1, Vec3 v2);

#endif

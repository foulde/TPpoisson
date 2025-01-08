# Rapport sur les méthode de résolution numérique du problème de diffusion de la chaleur Poisson




### Définition de l'équation : 

Le problème de diffusion de température selon poisson 1D est définie par l'équation suivante : 

```math
-k\frac{\partial^2 T}{\partial x^2} = g(x),  x \in [0,1]
```
Avec les boundary condition de dirichilet défine comme ceci .
```math
T(0) = T_0
T(1) = T_1
```

### Modélisation et discrétisation du domaine : 


Bien que l'équation ait une solution analytique 

```math
T(x) = T_0 + x(T_1 - T_0)
```

Ici on veut résoudre le problème en discrétisant le domaine et en utilisant la méthode des différence finie 

la formule générale est définie pour 
```math
u(a) = u_a, u(b) = u_b
```

on définie une discrétisation en N +1 élement et N segment uniforme à
l'aide du calcul suivant 

```math
h = \frac{b - a}{N}
```

On définie alors les xi comme 

```math
x_i = a + i h,  i = 0, 1, \dots, N.
```
En utilisant la méthode des différence finie on définie l'approximation suivante : 

```math
\frac{d^2u}{dx^2} \approx \frac{u_{i+1} - 2u_i + u_{i-1}}{h^2}.
```

On peut alors définir notre modèlisation comme 

```math
-h^2 \frac{u_{i+1} - 2u_i + u_{i-1}}{h^2} = f(x_i).
```

Grace à les méthode des différence finie on est capable de transformer un problème de résolution d'équation 
différentielles en un problème de résolution de système linéaire.


```math
A u = b

```

```math

A = \begin{bmatrix} -2 & 1 & 0 & \cdots & 0 \ 1 & -2 & 1 & \cdots & 0 \ 0 & 1 & -2 & \cdots & 0 \ \vdots & \vdots & \vdots & \ddots & 1 \ 0 & 0 & 0 & 1 & -2 \end{bmatrix}


```

ou u est la température recherché en chaque point du domaine discrétisé et b les termes sources dans notre cas 0
partout sauf au frontière du domaine .




## Deux méthode de résolutions 

Une fois la discrétisation du domaine définie nous faisons face à un choix résoudre le système de manière directe , en utilisé
Gauss , ou LU décomposition.

Ou utiliser des méthode itérative: 
Richardson Alpha , Jacobi ou Gauss Seidel.

Chacune présente ses avantages et inconvenient mais en général plus notre maillage est fin et contient d'élément plus on essayera d'utiliser des méthodes itérative.





## Référence et utilisation de BLAS/LAPACK 


### 1. En C, comment doit on déclarer et allouer une matrice pour utiliser BLAS et LAPACK ?

Pour déclarer une matrice dans blas et lapack on peut procédé comme ceci 


```code
double *A = (double *)malloc(m * n * sizeof(double));
```

La raison est que les bibliothèques BLAS et LAPACK manipulent les matrices sous la forme de tableaux unidimensionnels (1D) pour des raisons d'efficacité en mémoire.
Les éléments peuvent être rangé en : 

-ROWmajor ligne les unes à la suite des autres .
-COLmajor colonnes les unes à la suite des autres.

### 2. Quelle est la signification de la constante LAPACK COL MAJOR ?
La fonctionnalité LAPACK_COL_MAJOR permet préciser le format de stockage d'une matrice 
si on passe colmajor la matrice :

```math
A = \begin{bmatrix} a_1 & a_4 & a_7 \\ a_2 & a_5 & a_8 \\ a_3 & a_6 & a_9 \end{bmatrix}

```
sera stocker sur une ligne colonne par colonne comme ceci :
[a1,a2,a3,a4,a5,a6,a7,a8,a9]
 

### 3. A quoi correspond la dimension principale (leading dimension) généralement notée ld ?

La dimension principale (ld) dépend de l’ordre de stockage (Row-Major ou Column-Major) et représente la taille allouée en mémoire pour une ligne ou une colonne

En Column-Major, la dimension principale correspond au nombre total de lignes allouées en mémoire, qu’elles soient utilisées ou non. Cela permet de gérer des matrices avec un espace supplémentaire (padding).

En Row-Major, la dimension principale correspond au nombre total de colonnes.






### 4. Que fait la fonction dgbmv ? Quelle m´ethode impl´emente-t-elle ?
DGBMV  performs one of the matrix-vector operations

La fonction dgbmv effectue une opération de multiplication matrice-vecteur optimisée pour les matrices en bande. Elle calcule :

où :

    A est une matrice en bande,
    x et y sont des vecteurs,
    α et β sont des scalaires.

L'optimisation vient de l'utilisation de la structure spécifique des matrices bandes, en limitant les calculs aux éléments non nuls.


```math

y := \alpha \cdot A \cdot x + \beta \cdot y \quad \text{or} \quad y := \alpha \cdot A^T \cdot x + \beta \cdot y


``` 


### 5. Que fait la fonction dgbtrf ? Quelle m´ethode impl´emente-t-elle ?
Optimiser pour les matrices bandes cette fonction applique une elemination de gauss avec pivot partielle.
En tirant partie de sa structure en bande.
cette fonction créer une factorisation LU  pour une matrice bande.

La fonction dgbtrf effectue une factorisation LU d’une matrice en bande AA avec pivot partiel. Elle divise A en trois matrices :

```math
A = P \cdot L \cdot U

```
Où : 
P est une matrice de permutation,
LL est une matrice triangulaire inférieure,
UU est une matrice triangulaire supérieure.

### 6. Que fait la fonction dgbtrs ? Quelle m´ethode impl´emente-t-elle ?

La fonction dgbtrs résout un système d’équations linéaires de la forme :
```math
A \cdot x = b
```

Est une fonction optimiser pour les matrice bande factorisé au format PLU
cette fonction résoud un système déjà factorisé PLUx =y.






### 7. Que fait la fonction dgbsv ? Quelle méthode implémente-t-elle ?
La fonction dgbsv combine les étapes de factorisation LU (dgbtrf) et de résolution (dgbtrs) en une seule opération. Elle résout directement un système d’équations linéaires 
```math
A \cdot x = b
```
, où AA est une matrice en bande.

Elle est particulièrement efficace pour les matrices de grande taille présentant une structure en bande.




### 8. Comment calculer la norme du r´esidu relatif avec des appels BLAS ?



La norme du résidu relatif ou Erreur Avant est définie comme :
```math

ErreurAvant = \frac{\| b - A \cdot x \|}{\| b \|}

```

Pour la calculer à l’aide de BLAS :


1. Utiliser dgbsv pour résoudre le système
2. Calculer r=b−A⋅x à l’aide de la fonction dgemv (multiplication matrice-vecteur).
3. Calculer les normes ∥r∥ et ∥b∥ avec dnrm2.
4. Diviser ∥r∥ par ∥b∥ pour obtenir la norme relative/erreur avant 



## exercice 4 Stockage GB et appel à DGBMV

Le stockage GB (General Banded) est une structure spécifique utilisée dans BLAS/LAPACK pour les matrices en bande. Contrairement à une matrice pleine, seules les bandes non nulles de la matrice sont stockées dans un tableau 2D.

Pour la matrice de Poisson 1D, cela signifie que nous n'allons stocker que :

    La diagonale principale.
    La diagonale supérieure.
    La diagonale inférieure.


Définition du stockage GB (format colonne-major) :

    Le tableau GB est un tableau 2D de taille kl+ku+1  lignes et nn colonnes, où :
        n est le nombre de colonnes de la matrice.
        kl est le nombre de diagonales sous la diagonale principale.
        ku est le nombre de diagonales au-dessus de la diagonale principale.

    Les diagonales sont stockées horizontalement, avec un décalage pour aligner les éléments.






```math
A = \begin{bmatrix} 2 & -1 & 0 & 0 \\ -1 & 2 & -1 & 0 \\ 0 & -1 & 2 & -1 \\ 0 & 0 & -1 & 2 \end{bmatrix} \quad \text{devient :} \quad \begin{bmatrix} 0 & -1 & -1 & -1 \\ 2 & 2 & 2 & 2 \\ -1 & -1 & -1 & 0 \end{bmatrix}
```

### utilisation de dgbmv 
on peut appeler dgbmv comme ceci

```c
dgbmv("N", n, n, kl, ku, alpha, AB, ldab, x, incx, beta, y, incy);
```


### Écrire le stockage GB en priorité colonne pour la matrice de Poisson 1D 




##   Exercice 5 DGBTRF, DGBTRS, DGBSV

1. R´esoudre le syst`eme lin´eaire par une m´ethode directe en faisant appel `a LAPACK.
Pour résoudre le système on utilisera : 

DGBTRF pour effectuer la factorisation LU.
DGBTRS pour résoudre la factorisation LU.
DGBSV pour effectuer les deux l'un après l'autre. 

2. ´Evaluer les performances. Que dire de la complexit´e des m´ethodes appel´ees ?

## LU pour les matrices tridiagonales

Pour la méthode de validation on décide 

## Exercice 6 : LU pour les matrices tridiagonales


## Partie 3 

### Méthode de Richardson Générale

L'algorithme de Richardson peut être généralisé sous la forme suivante pour résoudre un système A⋅x=bA⋅x=b :

```math
    x^{(k+1)} = x^{(k)} + M^{-1} \cdot (b - A \cdot x^{(k)})

    où :


    M^{-1} est une matrice qui influence la convergence de la méthode,

    b - A \cdot x^{(k)} est le résidu au pas kk.
```



Les méthodes spécifiques sont définies par le choix de M^{-1}.

### Implémentation C - Richardson
L'équation de Richardson est utilisée pour résoudre un système linéaire A⋅x=bA⋅x=b en approximant progressivement la solution x. L'algorithme de base est donné par :



```math
x_{k+1} = x_k + \alpha \cdot (b - A \cdot x_k)
```

où :

    xk​ est l'approximation courante,
    α est un paramètre d'accélération,
    rk=b−A⋅xk est le résidu au pas k.

Richardson classique (alpha constant) :
Dans ce cas, le paramètre α est une constante choisie pour accélérer la convergence. Le choix de αα dépend du spectre des valeurs propres de A. Si α est trop grand ou mal choisi, la méthode peut diverger.



```math
x_{k+1} = x_k + \alpha \cdot (b - A \cdot x_k)
```



```math
x_{k+1} = x_k + \alpha \cdot (b - A \cdot x_k)
```



```math
x_{k+1} = x_k + \alpha \cdot (b - A \cdot x_k)
```

### Implémentation C - Jacobi


Méthode de Jacobi

La méthode de Jacobi est un cas particulier où M^{-1} est défini par la matrice diagonale D de A. Ainsi, M = D et l'algorithme devient :

```math
x^{(k+1)} = x^{(k)} + D^{-1} \cdot (b - A \cdot x^{(k)})
```

En explicitant A = D + L + U (où L est la partie strictement inférieure et U la partie strictement supérieure), cela équivaut à :
```math
x^{(k+1)} = D^{-1} \cdot (b - (L + U) \cdot x^{(k)})
```

### Implémentation C - Gauss-Seidel


Méthode de Gauss-Seidel

La méthode de Gauss-Seidel améliore Jacobi en utilisant les nouvelles valeurs calculées dans la même itération. Ici, M^{-1} est défini comme (D - E)^{-1}, où E est la partie strictement inférieure de A. Ainsi, M = D - E, et l'algorithme devient :

```math
x^{(k+1)} = x^{(k)} + (D - E)^{-1} \cdot (b - A \cdot x^{(k)})
```
En pratique, cela équivaut à une mise à jour itérative composant par composant :

```math
x^{(k+1)}_i = \frac{1}{a_{ii}} \left( b_i - \sum_{j<i} a_{ij} \cdot x^{(k+1)}_j - \sum_{j>i} a_{ij} \cdot x^{(k)}_j \right)
```
<!-- “h\=Nb−a​.” -->

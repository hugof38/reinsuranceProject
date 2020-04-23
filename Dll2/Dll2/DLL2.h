#pragma once
// Le bloc ifdef suivant est la fa�on standard de cr�er des macros qui facilitent l'exportation
// � partir d'une DLL plus simple. Tous les fichiers contenus dans cette DLL sont compil�s avec le symbole DLL1_EXPORTS
// d�fini sur la ligne de commande. Ce symbole ne doit pas �tre d�fini pour un projet
// qui utilise cette DLL. Ainsi, les autres projets dont les fichiers sources comprennent ce fichier consid�rent les fonctions
// DLL1_API comme �tant import�es � partir d'une DLL, tandis que cette DLL consid�re les symboles
// d�finis avec cette macro comme �tant export�s.
#ifdef DLL2_EXPORTS
#define DLL2_API __declspec(dllexport)
#else
#define DLL2_API __declspec(dllimport)
#endif

DLL2_API int simul_poiss(const double lambd);
DLL2_API double* simul_collective(const int nbSimul, const double lambd, const double mu, const double sigma);
DLL2_API int __stdcall wrap_simul_col(double* result_simul_time, const int nbSimul, const double lambd, const double mu, const double sigma);
#pragma once
// Le bloc ifdef suivant est la façon standard de créer des macros qui facilitent l'exportation
// à partir d'une DLL plus simple. Tous les fichiers contenus dans cette DLL sont compilés avec le symbole DLL1_EXPORTS
// défini sur la ligne de commande. Ce symbole ne doit pas être défini pour un projet
// qui utilise cette DLL. Ainsi, les autres projets dont les fichiers sources comprennent ce fichier considèrent les fonctions
// DLL1_API comme étant importées à partir d'une DLL, tandis que cette DLL considère les symboles
// définis avec cette macro comme étant exportés.
#ifdef DLL2_EXPORTS
#define DLL2_API __declspec(dllexport)
#else
#define DLL2_API __declspec(dllimport)
#endif

DLL2_API int simul_poiss(const double lambd);
DLL2_API double* simul_collective(const int nbSimul, const double lambd, const double mu, const double sigma);
DLL2_API int __stdcall wrap_simul_col(double* result_simul_time, const int nbSimul, const double lambd, const double mu, const double sigma);
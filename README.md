Projet INPS deuxième session 

Description
Le projet a pour but de calculer les solutions de l'équation de Schrodinger en dimension 1.

Structure du projet
├── Makefile                # Makefile principal à la racine
├── src/                    # Fichiers sources (.cpp principalement)
├── include/                # Fichiers d'en-tête (.h)
├── test/ ├── test_poly     # Tests unitaires
          ├── Makefile
├── bin/                    # Exécutable du projet principal
├── Doxyfile                # Documentation technique
├── pres/                   # Présentation du projet
├── obj/                    # Fichiers objet (.o) présents lors de la compilation du projet

Compilation et éxecution
Pour compiler le projet principal:
'''bash
make
Pour éxécuter le projet principal:
make run
Pour compiler et éxecuter les tests unitaires depuis la racine du projet:
make tests
Pour nettoyer les fichiers suite à la compilation du projet principal et/ou des tests unitaires:
make clean
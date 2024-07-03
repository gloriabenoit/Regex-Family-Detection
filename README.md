# Détection de motifs de familles protéiques par expressions régulières (L3BI - 2022)
Ce notebook est la production du cours de Programmation Python 2 de ma L3 à Université Paris Cité, organisé par M. Pierre Poulain. Ce projet visait à détecter l'appartenance d'une protéine à une famille à partir d'une signature retrouvée dans sa
séquence, grâce aux expressions régulières.

## Installation

Nous avons utilisé les packages `re` et `urllib.request` pour respectivement détecter un motif et récupérer la fiche d'une protéine. Nous avons aussi utilisé `BioPython`, une alternative à nos récupérations *PROSITE*.

## Fonctionnement

A partir d'un motif *PROSITE*, nous avons automatisé sa traduction en expression régulière pour pouvoir le rechercher dans des séquences. 

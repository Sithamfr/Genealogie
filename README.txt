L'objectif du programme R est de reconstituer la généalogie la plus probable d'une population de rosiers au sens du maximum de vraisemblance.

On dispose pour cela de jeu de données sous la forme suivante :
	Chaque individu est représenté en ligne. Sur chaque ligne on a donc les résultats, pour 4 gènes différents, de l'obtention de signaux correspondant à tel ou tel allèle du gène étudié. Les individus étant soit diploïdes, soit tetraploïdes, on a alors au maximum 4 allèles différents pour chacun des gènes.
	À ces informations, on rajoute la ploïdie de l'individu (2 ou 4), son numéro de génération (plus le nombre est faible, plus l'individu est ancien) et son identifiant.
	On a alors chaque ligne sous la forme :
		[S11,S12,S13,S14,S21,S22,S23,S24,S31,S32,S33,S34,S41,S42,S43,S44,Plo,Gen,Indiv]

Pour chaque individu, on cherche le couple de parents le plus probable. Ce processus est adapté aux règle qui régissent les mécanismes de reproduction des rosiers. Ces mécanismes sont bien évidemment beaucoup plus compliqués à modéliser pour les individus tetraploïdes.

Fichiers :
genealogie.R 	# Programme R qui reconstitue la généalogie
test_di.csv 	# Population diploïde simulée
test_tetra.csv 	# Population tetraploïde simulée
testgen.csv 	# Population simulée
genim.png 	# Exemple de sortie graphique sur la population testgen

Sur les exemples de sorties graphiques, on représente :
	* Un noeud de forme ronde ou carrée représente respectivement un individu diploïde ou un individu tetraploïde
	* La génération de l'individu est donnée par sa couleur (foncé = ancien)
	* Une arête du graphe est plus épaisse si la probabilité associée est >90%
	* Une arête du graphe est en pointillé si plusieurs couples d'individus avaient la même proba d'être ses parents et qu'il y a alors eu tirage au sort pour déterminer le couple en sortie

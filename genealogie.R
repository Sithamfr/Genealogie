###############       GENEALOGIE DE ROSIERS       #################



library("igraph")
library("RColorBrewer")




# Paramètres pi de dosage resp. pour 2 et 3 signaux détectés
# Par défaut, 1/3 par dosage mais l'utilisateur peut modifier les variables suivantes :

pi2 = c(1/3,1/3,1/3) # c( pi22 , pi31 , pi13 )
pi3 = c(1/3,1/3,1/3) # c( pi211 , pi121 , pi112 )


# FONCTIONS ------------------------------------------------------------------------------------------

# Fonction qui duplique le signal a de (a,0,0,0) de chaque individu diploïde de la matrie M
# en (a,a,0,0). Cette fonction est appelée au début pour uniformiser les notations.
duplic2x = function(M){
  for (i in 1:n){
    if (M[i,17]==2){
      for (j in c(1,5,9,13)){
        if (M[i,(j+1)]==0){
          M[i,(j+1)]=M[i,j]
        }}}}
  return(M)
}


# compte le nombre d'enfants virtuels de meme genotype que l'enfant de reference
# in : - enfants virtuels (matrice avec un genotype sur chaque ligne)
#      - genotype de l'enfant (vecteur)
# out : le nombre d'enfants virtuels correspondant
comparerEnf = function(EnfVirt, Enf){
  g_enf = sort(Enf)
  Cpt = 0
  for (i in 1:nrow(EnfVirt)){
    g = sort(EnfVirt[i,]) # Tri croissant
    if (sum(abs(g-g_enf))==0){Cpt = Cpt + 1} # Compare et incrémente si vecteurs égaux
  }
  return(Cpt)
}

# genere la liste des enfants virtuels a partir du genotype des parents pour le schema 2x-2x
# in : genotype des parents (vecteurs de taille 2)
# out : liste des enfants virtuels (matrice de 4 lignes et 2 colonnes)
genererEnfVirt22 = function(GenP1s, GenP2s){
  EnfVirt = as.matrix(expand.grid(GenP1s,GenP2s)) # Combinaisons génétiques possibles
  return(EnfVirt)
}

# genere la liste des enfants virtuels a partir du genotype des parents pour le schema 4x-4x
# in : genotype des parents (vecteurs de taille 4)
# out : liste des enfants virtuels (matrice de 36 lignes et 4 colonnes)
genererEnfVirt44 = function(GenP1s, GenP2s){
  DonP1 = combn(GenP1s,2) # Code génétique possiblement transmis par le parent1
  DonP2 = t(combn(GenP2s,2)) # Code génétique possiblement transmis par le parent2
  transP1 = t(matrix(rep(combn(GenP1s,2),6),ncol=36))
  transP2 = DonP2[rep(1:6,each=6),1:2]
  EnfVirt = cbind(transP1,transP2) # Combinaisons génétiques possibles
  return(EnfVirt)
}


# Fonction qui calcule les dosages pour chaque vecteur de taille 4 en entrée
# Si n=1 ou n=4, on passe juste du type vecteur au type matrice
# Si n=2 ou n=3, on met sur chaque ligne de la matrice un dosage possible
comb_dosage = function(signal,n){
  if (n==1){LP = matrix(rep(signal,4),ncol=4)}
  else if (n==4){LP = matrix(signal,ncol=4)}
  else if (n==2){LP = cbind(t(matrix(rep(signal,3),ncol=3)),rbind(signal,matrix(rep(signal,2),ncol=2)))}
  else if (n==3){LP = cbind(t(matrix(rep(signal,3),ncol=3)),signal)}
  return(LP)
}


# calcul de la proba d'un lien entre deux parents et un enfant
# in : infos genetiques sur les deux parents et l'enfant (vecteurs de taille 17 : genotypes 4x4 + ploidie)
# out : proba du lien
calculerProbaLien = function(GenP1, GenP2, GenE){
  p = 1
  plo_e = as.integer(GenE[17])
  plo_p1 = as.integer(GenP1[17])
  plo_p2 = as.integer(GenP2[17])
  if ((plo_p1==plo_e) & (plo_p2==plo_e)){
    for (k in c(1,5,9,13)){
      vE = as.numeric(GenE[k:(k+plo_e-1)])
      vP1 = as.numeric(GenP1[k:(k+plo_e-1)])
      vP2 = as.numeric(GenP2[k:(k+plo_e-1)])
      
      # DIPLOIDES
      if (plo_e==2) {
        EnfVirt = genererEnfVirt22(vP1,vP2)
        p = p * comparerEnf(EnfVirt,vE)/4
      }
      
      # TETRAPLOIDES
      else if (plo_e==4) {
        
        vP1 = vP1[vP1!=0] # allèles remplis
        vP2 = vP2[vP2!=0]
        vE = vE[vE!=0]
        
        nP1 = length(vP1)
        nP2 = length(vP2)
        nE = length(vE)
        
        LP1 = comb_dosage(vP1,nP1) # Dosage éventuel
        nrowLP1 = nrow(LP1)
        if (nP1 %in% c(1,4)){pi_p1=1} # Valeurs de pi correspondantes
        else if (nP1==2){pi_p1=pi2}
        else {pi_p1=pi3}
        
        LP2 = comb_dosage(vP2,nP2)
        nrowLP2 = nrow(LP2)
        if (nP2 %in% c(1,4)){pi_p2=1}
        else if (nP2==2){pi_p2=pi2}
        else {pi_p2=pi3}
        
        LE = comb_dosage(vE,nE)
        nrowLE = nrow(LE)
        if (nE %in% c(1,4)){pi_e=1}
        else if (nE==2){pi_e=pi2}
        else {pi_e=pi3}
        
        pl = 0
        for (i in 1:nrowLP1){
          parent1 = LP1[i,]
          for (j in 1:nrowLP2){
            parent2 = LP2[j,]
            EnfVirt = genererEnfVirt44(parent1,parent2)
            for (k in 1:nrowLE){
              enfant = LE[k,]
              cpt = comparerEnf(EnfVirt,enfant)
              pl = pl + cpt/36 * pi_p1[i] * pi_p2[j] * pi_e[k]
            }
          }
        }
        p = p * pl
      }
    }
  }
  else {return(0)}
  return(p)
}


# construire pour un enfant la matrice des probas associee a tous les couples de parents potentiels
# in : numero de l'enfant
# out : matrice des probas (n lignes et n colonnes)
calculerProbasParents = function(e){
  MatProbasParents = matrix(0, nrow=n, ncol=n)
  enf = MatIndiv[e,1:17]
  g_e = MatIndiv[e,18]
  M = MatIndiv[-c(e),] # Matrice sans la ligne de l'enfant
  for (p1 in 1:(nrow(M)-1)){
    g_p1 = M[p1,18]
    for (p2 in (p1+1):nrow(M)){
      g_p2 = M[p2,18]
      if((g_p1<g_e) & (g_p2<g_e)){
        MatProbasParents[p1,p2] = calculerProbaLien(M[p1,1:17],M[p2,1:17],enf)
      }
    }
  }
  
  S = sum(MatProbasParents)
  if (S>0) {MatProbasParents = MatProbasParents/S} # Normalisation
  return(MatProbasParents)
}


# pour un enfant, recherche du couple de parents le plus probable
# in : numero de l'enfant
# out : couple de parents le plus probable (vecteur de taille 4 : enfant/parent 1/parent 2/log-proba du lien)
recupParentsMax = function(e){
  M = calculerProbasParents(e)
  pmax = max(M)
  if (pmax==0){IndMax = c(0,0)}
  else {
    IndMax = as.numeric(which(M == pmax, arr.ind = TRUE))
    if (IndMax[1]>=e){IndMax[1]=IndMax[1]+1} # On corrige ici la ligne perdue en enlevant
    if (IndMax[2]>=e){IndMax[2]=IndMax[2]+1} # celle de e dans calculerProbasParents
    nmax = length(IndMax)%/%2
    if (nmax > 1){
      rand = sample((1:nmax),1)
      IndMax = c(IndMax[rand],IndMax[rand+nmax]) # on tire au sort les parents
      # randbool : variable globale qui renvoie les individus pour lesquels on a fait un tirage au sort
      randbool <<- c(randbool,e) # Sert uniquement pour tracer le graphe
    }
    pmax = log(pmax)}
  ParentsMax = c(e,IndMax,pmax)
  return(ParentsMax)
}


# construire la genealogie de plus grande proba associee a la matrice des individus
# in : rien
# out : genealogie la plus probable (matrice de n lignes et 4 colonnes) munie de la log-vraisemblance
construireGen = function(){
  MatIndiv = duplic2x(MatIndiv)
  randbool <<- c() # On mettra ici les individus pour lesquels on a fait un tirage au sort sur les parents
  Gen = matrix(0, nrow=n, ncol=4)
  # Cherche le couple de parents le plus probable pour chaque individu
  for (ind in 1:nrow(MatIndiv)){
    Gen[ind,] = recupParentsMax(ind)
  }
  return(list(Gen = Gen, lLik = sum(Gen[, 4])))
}


# representation graphique de la genealogie a l'aide de igraph (noeuds = individus, fleches = liens)
# in : genealogie (matrice de n lignes et 4 colonnes)
# out : rien
representerGen = function(GenMax){
  
  gens = MatIndiv$Gen
  plo = MatIndiv$Plo
  M = GenMax$Gen
  nrel = nrow(M)
  M[,4] = exp(M[,4])
  ngens = length(unique(gens))
  geninv = ngens-gens
  
  forme = c("circle","square")
  plo = forme[plo/2]
  
  Grp = graph.empty(directed=TRUE)
  Grp = Grp + vertex(GenMax$Gen[, 1], color=geninv, shape=plo, size=10, label.cex=1.5, label.color="#42285f")
  
  palette = brewer.pal(n = 12, name = "Set3")[c(1,3,4,5,6,10)]
  for (ind in 1:nrel){
    if (min(M[ind,2:3])>0){
      if (ind %in% randbool){linetype=2} # Pointillés si il y a eu tirage au sort pour cet individu
      else {linetype=1} # Trait normal sinon
      p = round(M[ind,4],3)
      couleur = palette[ind%%6+1]
      w = 2*((abs(M[ind,4])>0.9)+1) #On double l'epaisseur des fleches quand la proba et d'au moins 90%
      Grp = Grp + edge(M[ind,2],M[ind,1],label=p,color=couleur,width=w/2,lty=linetype,
                       arrow.size=w/30,arrow.width=w/2)
      Grp = Grp + edge(M[ind,3],M[ind,1],label=p,color=couleur,width=w/2,lty=linetype,
                       arrow.size=w/30,arrow.width=w/2)
    }
  }
  Grp$palette = colorRampPalette(c("#e4f7d7","#70ab48"))(ngens)
  
  plot(Grp, edge.label.color="black",vertex.label.family="Calibri",
       edge.curved=0.2,vertex.frame.color="#538b36",edge.label.cex=0.6)
  title("Reconstitution de généalogie de rosiers",cex.main=1.2,col.main="#224451")
  legend('topleft',legend=c("Génération ancienne","Génération récente"),
         fill=c(Grp$palette[length(Grp$palette)],"white"),box.lty=0,cex=0.65,bty="n") 
}
 



# EXECUTION -------------------------------------------------------------------------------------------------

MatIndiv = read.csv('testgen.csv',sep=',')
#MatIndiv = read.csv('Pop2x.csv',sep=';')
#MatIndiv = read.csv('Pop4x.csv',sep=';')

#MatIndiv = MatIndiv[MatIndiv$Plo %in% c(2,4),] # Enleve les ploidies non gerees s'il y en a

n = nrow(MatIndiv)
GenMax = construireGen()
representerGen(GenMax)

def fonction_relais(E1,E2,p,r,k1,k2,i,Kb):#cette fonction fait le lien entre la première et la deuxième étape elle "commande" le calcul de la 1ère étape et interprète la sortie de la 1ère étape afin de s'en servir pour la seconde étape
	L=[]
	L=calcul_ima_iso_nv(E1,E2,p,r,k1,k2,i)
	L0=L[0]
	LQ=L0[0]
	Q1=L0[1]
	E1=Q1.curve()
	E2=L0[3].curve()
	L1=L[1]
	LP=L1[0]
	P1=L1[1]
	E1b=construction_lift(E1,Kb,p)
	E2b=construction_lift(E2,Kb,p)
	P1b=E1b(P1.xy()[0],P1.xy()[1])
	Q1b=E1b(Q1.xy()[0],Q1.xy()[1])
	return fonction_etape2_1(E1b,p,P1,Q1,E2b,k1,k2,LP,LQ)
	
	

def fonction_etape2_1(E1,p,P1,Q1,E2,k1,k2,LP,LQ): #on a en entrée les deux listes LP et LQ de points de E2 pour les images de P1 et Q1 deux points générateurs de la 2 torsion sur E1, on a aussi en entrée la caractéristique du corps
	phi1=EllipticCurveIsogeny(E1,2^(k2-1)*(P1+Q1)) #on construit l'isogénie descendante qui nous permet d'aller un niveau en dessous du cratère.
	u=randint(0,len(LP)-1)	
	v=randint(0,len(LQ)-1) #on prend des points au hasard dans la liste, ils engendrent tous de toute façon la bonne 2-isogénie (descendante)
	phi2=EllipticCurveIsogeny(E2,2^(k2-1)*(LP[u]+LQ[v]))
	E1b=phi1.codomain()#on se place sur la courbe située en dessous
	P1b=phi1(P1) #on calcule les points pour avoir une nouvelle base dans laquelle on pourra exprimer l'action du morphisme de Frobenius
	L=phi1(P1+Q1).division_points(2)
	u=randint(0,len(L)-1)
	Q1b=L[u]
	#on calcule maintenant l'action du Frobenius, la seule qui nous intéresse c'est celle sur Q1b, l'autre (P1b) étant triviale...
	V=[]#on va stocker dans V le vecteur représentant l'action du morphisme de Frobenius
	Q1f=E1b(Q1b.xy()[0]^p,Q1b.xy()[1]^p) #on calcule l'image de Q1b sous l'action du frobenius
	for i in range(P1b.order()):
		if i*P1b+p*Q1b == Q1f :
			V.append(i)
			V.append(p)
	#maintenant il faut filtrer sur l'autre courbe ceux qui ont la même matrice de Frobenius 
	LP2=[]
	LQ2=[]
	E2b=phi2.codomain()	
	for x in LP:
		a=len(LQ2)		
		P2b=phi2(x)
		for q in LQ:
			L=phi2(x+q).division_points(2)
			u=randint(0,len(L)-1)
			Q2b=L[u]
			Q2f=E2b(Q2b.xy()[0]^p,Q2b.xy()[1]^p)
			if Q2f==V[0]*P2b+V[1]*Q2b :
				LQ2.append(q)
		b=len(LQ2)
		if b>a :
			LP2.append(x)
	NL=[]	
	for y in LQ2:#on enlève les points qui sont en double dans la liste
        	a=len(LQ2)
              	LQ2=filter(lambda x : x != y , LQ2)
               	b=len(LQ2)
              	if b<a:
	              	NL.append(y)
        LQ2=NL
	return LP2,LQ2,V

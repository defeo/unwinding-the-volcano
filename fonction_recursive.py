def fonction_reccursive_g(E1,phi1,p,P1,Q1,E2,phi2,P2,Q2,k1,k2,L,V,L1,L2):
	if (len(L)==1 or len(L)==len(V) ): #condition d'arret
		return [V,Q1]
	else : #on fait du récursif
                L1.append(phi1)
                L2.append(phi2)#on stocke les isogénies successives
                print("boucleg","E1", E1.j_invariant() ,"E2", E2.j_invariant(),len(L),len(V))
		R1=Q1  #on stocke pour le cas où l'on fait l'étape de trop
		Q1=phi1(Q1)
		E1=phi1.codomain()#on se déplace suivant l'isogénie
		M=E1([0,1,0]).division_points(2^k1)
                M=filter(lambda x : x.order()==2^k1 , M)
		N=[]	
		for l in M:
			N.append(2^(k1-k2)*l)
		N=filter(lambda x : E1(x.xy()[0]^p,x.xy()[1]^p)==x,N)
		u=randint(0,len(N)-1)	
		P1=N[u] #pas nécessaire de recalculer P1, on pourrait engendrer l'isogénie suivante par l'image de P1 multipliée par la bonne puissance de 2
			
		NL=[]
                NL2=[]
                Q2=phi2(Q2)
                E2=phi2.codomain()
		M=E2([0,1,0]).division_points(2^k1)
                M=filter(lambda x : x.order()==2^k1 , M)
                N=[]	
		for l in M:
			N.append(2^(k1-k2)*l)
		N=filter(lambda x : E2(x.xy()[0]^p,x.xy()[1]^p)==x,N)
		u=randint(0,len(N)-1)	
		P2=N[u]
                V=L
                NL=[]
		#on va maintenant restreindre la taille de la liste L
		if E1(Q1.xy()[0]^p,Q1.xy()[1]^p)==Q1:#inutile dans ce cas là, normalement c'est le prochain if au même niveau qui sera vérifié
                    	print "probleme rationalité Q1"
		else :
                    	for x in L:
                        	R=phi2(x)
                        	if E2(R.xy()[0]^p,R.xy()[1]^p)!=R :
                            		NL.append(R)
                if E1(Q1.xy()[0]^p,Q1.xy()[1]^p)==p*Q1 :
			for x in NL:
				R=x
                        	if E2(R.xy()[0]^p,R.xy()[1]^p)==p*R :
                            		NL2.append(R)
                else:
                    	print("on est dans le second cas")
                    	for x in NL:
                        	R=x
                        	if E2(R.xy()[0]^p,R.xy()[1]^p)!=p*R :
                            		NL2.append(R)
		NL=[]
                for y in NL2:#on enlève les points qui sont en double dans la liste
                	a=len(NL2)
                    	NL2=filter(lambda x : x != y , NL2)
                    	b=len(NL2)
                    	if b<a:
                        	NL.append(y)
                NL2=NL
                if len(NL2)==len(V): #test de condition d'arrêt, marche si on a fait une étape en trop, d'où le fait que l'on mette V dans le return
                    	return [V,R1]
		phi1=EllipticCurveIsogeny(E1,2^(k2-1)*P1)#on construit les isogénies suivantes
                phi2=EllipticCurveIsogeny(E2,2^(k2-1)*P2)
		return fonction_reccursive_g(E1,phi1,p,P1,Q1,E2,phi2,P2,Q2,k1,k2,NL2,V,L1,L2)

def fonction_reccursive_d(E1,phi1,p,P1,Q1,E2,phi2,P2,Q2,k1,k2,L,V,L1,L2):
	if (len(L)==1 or len(L)==len(V) ): #condition d'arret
		return [V,P1]
	else : #on fait du récursif
                L1.append(phi1)
                L2.append(phi2)#on stocke les isogénies successives
                print("boucle","E1", E1.j_invariant() ,"E2", E2.j_invariant(),len(L),len(V))
		R1=P1 #on stocke pour le cas où l'on fait l'étape de trop 		
		P1=phi1(P1)
		E1=phi1.codomain()
		#IQ0=Q1.xy()[0]^p 
		#IQ1=Q1.xy()[1]^p 
                M=E1([0,1,0]).division_points(2^k2)
                M=filter(lambda x : x.order()==2^k2 , M)
		M=filter(lambda x : E1(x.xy()[0]^p,x.xy()[1]^p)==p*x , M)
                u=randint(0,len(M)-1)
		Q1=M[u] #pas nécessaire de recalculer P1
		#on recalcule l'image du point de base par l'action du Frobenius, on va maintenant restreindre la taille de la liste L
		NL=[]
                NL2=[]
                P2=phi2(P2)
                E2=phi2.codomain()
                M=E2([0,1,0]).division_points(2^k2)
                M=filter(lambda x : x.order()==2^k2 , M)
		M=filter(lambda x : E2(x.xy()[0]^p,x.xy()[1]^p)==p*x , M)
                u=randint(0,len(M)-1)
		Q2=M[u]
                d=0
                V=L
                NL=[]
          	for x in L:
                	R=phi2(x)
                       	if E2(R.xy()[0]^p,R.xy()[1]^p)==R : #test inutile à mon avis
	          		NL.append(R)
                NL2=NL
		NL=[]
                for y in NL2:#on enlève les points qui sont en double dans la liste
                	a=len(NL2)
                    	NL2=filter(lambda x : x != y , NL2)
                    	b=len(NL2)
                    	if b<a:
                        	NL.append(y)
                NL2=NL
                if len(NL2)==len(V): #test de condition d'arrêt rendant l'autre inutile...
                    	return [V,R1]
		phi1=EllipticCurveIsogeny(E1,2^(k2-1)*Q1)
                phi2=EllipticCurveIsogeny(E2,2^(k2-1)*Q2)
		return fonction_reccursive_d(E1,phi1,p,P1,Q1,E2,phi2,P2,Q2,k1,k2,NL2,V,L1,L2)

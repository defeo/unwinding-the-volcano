def calcul_ima_iso_nv(E1,E2,p,r,k1,k2,i): #k1, k2 sont les exposants de la 2torsion infinie pour les courbes l-isogènes définies sur p^r, i étant le degré de l'isogénie impaire
	c=k2
	E1i=E1
	E2i=E2 #on stocke les courbes de départ
	while(c>1) : #on crée le décalage sur le cratère pour revenir ensuite sur les courbes de départs lorsque l'on détermine les images des points générateurs de 2 torsion P et Q. On va d'abord commencer par Q
		L=[]
                L=E1([0,1,0]).division_points(2^k2)
		L=filter(lambda x : x.order()==2^k2, L) 
		L=filter(lambda x : E1(x.xy()[0]^p,x.xy()[1]^p)==p*x , L)
                if len(L)==0:
			print("problème longueur L, diagonale")
                u=randint(0,len(L)-1)		
		Q1=L[u]
		Q1=2^(k2-1)*Q1
		E1=EllipticCurveIsogeny(E1,Q1).codomain()		
		L=E2([0,1,0]).division_points(2^k2)
		L=filter(lambda x : x.order()==2^k2, L) 
		L=filter(lambda x: E2(x.xy()[0]^p,x.xy()[1]^p)==p*x , L)
                if len(L)==0:
			print("problème longueur L, diagonale")
		u=randint(0,len(L)-1)		
		Q2=L[u]
		Q2=2^(k2-1)*Q2
		E2=EllipticCurveIsogeny(E2,Q2).codomain()		
		c=c-1
	#maintenant on détermine les candidats possibles pour l'image de Q1
	L=E1([0,1,0]).division_points(2^k1)
	M=filter(lambda x : x.order()==2^k1 , L)
	N=[]	
	for l in M:
		N.append(2^(k1-k2)*l)
	N=filter(lambda x : E1(x.xy()[0]^p,x.xy()[1]^p)==x,N)
	u=randint(0,len(N)-1)	
	P1=N[u] #on prend un point de 2^k2 torsion rationnelle
	M=filter(lambda x : x.order()==2^k2 , L)	
	N=filter(lambda x : E1(x.xy()[0]^p,x.xy()[1]^p)==p*x , M)#on prend un point de 2^k2 torsion non rationnel qui engendre l'isogénie horizontale
	N=filter(lambda x: x.weil_pairing(P1,2^k2).multiplicative_order()==2^k2,N)
	print ("len(N)", len(N))
        u=randint(0,len(N)-1)
        Q1=N[u]
        L=E2([0,1,0]).division_points(2^k1)
	M=filter(lambda x : x.order()==2^k1 , L)
	N=[]	
	for l in M:
		N.append(2^(k1-k2)*l)
	N=filter(lambda x : E2(x.xy()[0]^p,x.xy()[1]^p)==x,N)
	u=randint(0,len(N)-1)	
	P2=N[u] #on prend un point de 2^k2 torsion rationnelle
	M=filter(lambda x : x.order()==2^k2 , L)	
	N=filter(lambda x : E2(x.xy()[0]^p,x.xy()[1]^p)==p*x , M) 
	N=filter(lambda x: x.weil_pairing(P2,2^k2).multiplicative_order()==2^k2,N)
	print ("len(N)", len(N))        
	u=randint(0,len(N)-1)
        Q2=N[u]
	Liste=N #on stocke cette liste car c'est l'ensemble des points qui ont la même image par l'action du Frobenius que Q
	phi1=EllipticCurveIsogeny(E1,2^(k2-1)*P1)#on quotiente par 2^(k2-1)*P, on va donc abaisser son ordre ce qui va réduire l'espace des candidats par 2 à chaque itération pour l'image de Q
	phi2=EllipticCurveIsogeny(E2,2^(k2-1)*P2)
	ListeFinale=[]
        V=[]
        L1=[]
        L2=[]
	ListeFinale=fonction_reccursive_g(E1,phi1,p,P1,Q1,E2,phi2,P2,Q2,k1,k2,Liste,V,L1,L2)#on fait la réduction des candidats à l'aide d'une fonction récursive qui à chaque itération se déplace d'une courbe sur le cratère
	#ce qui suit est juste de la vérification concernant le fait que le point de base de la 2^{notations oubliées!} torsion	a bien son image parmi la liste des candidats sorti par l'algorithme, la vérification se fera en observant l'output de l'algorithme
	QU=ListeFinale[1]#on recupère le point que l'on veut déterminer l'image par l'isogénie, ce point se situe sur une courbe isomorphe à la courbe d'entrée de l'algorithme (tout porte à croire cela, pas démontré pour le moment)
	E2=ListeFinale[0][0].curve()#on récupère la courbe sur laquelle on arrive par l'isogénie de degré impair
	E1=QU.curve()
	L=E1([0,1,0]).division_points(i)
	L=filter(lambda x: x.order()==i, L)
	c=0
	while c<len(L):#on retrouve l'isogénie qui relie les deux courbes en faisant une recherche exhaustive
		phi=EllipticCurveIsogeny(E1,L[c])
		if phi.codomain().j_invariant()==E2.j_invariant():
			c=len(L)
			QUI=phi(QU)
		else:
			c=c+1
	ListeFinale.append("QUI")#on rajoute les résultats trouvés à la variable qui sera retournée à la fin de l'algorithme	
	ListeFinale.append(QUI)	
	c=k2
	E1=E1i #on récupère les courbes de départ
	E2=E2i
	while(c>1) : #on crée le décalage sur le cratère (dans l'autre sens) pour revenir ensuite sur la courbe de départ 
		L=[]
                L=E1([0,1,0]).division_points(2^k2)
		L=filter(lambda x: x.order()==2^k2, L)
		L=filter(lambda x: E1(x.xy()[0]^p,x.xy()[1]^p)==x , L)
		if len(L)==0:
			print("problème longueur L, diagonale")
		u=randint(0,len(L)-1)		
		P1=L[u]
		P1=2^(k2-1)*P1
		E1=EllipticCurveIsogeny(E1,P1).codomain()		
		L=E2([0,1,0]).division_points(2^k2)
		L=filter(lambda x: x.order()==2^k2, L)
		L=filter(lambda x: E2(x.xy()[0]^p,x.xy()[1]^p)==x , L)
		if len(L)==0:
			print("problème longueur L, diagonale")		
		u=randint(0,len(L)-1)		
		P2=L[u]
		P2=2^(k2-1)*P2
		E2=EllipticCurveIsogeny(E2,P2).codomain()		
		c=c-1
	#maintenant on détermine les candidats possibles pour l'image de P1	
	L=E1([0,1,0]).division_points(2^k1)
	M=filter(lambda x : x.order()==2^k1 , L)
	N=[]	
	for l in M:
		N.append(2^(k1-k2)*l)
	N=filter(lambda x : E1(x.xy()[0]^p,x.xy()[1]^p)==x,N)
	u=randint(0,len(N)-1)	
	P1=N[u] #on prend un point de 2^k2 torsion rationnelle
	M=filter(lambda x : x.order()==2^k2 , L)	
	N=filter(lambda x : E1(x.xy()[0]^p,x.xy()[1]^p)==p*x , M)#on prend un point de 2^k2 torsion non rationnel qui engendre l'isogénie horizontale
	N=filter(lambda x: x.weil_pairing(P1,2^k2).multiplicative_order()==2^k2,N)
	print ("len(N)", len(N))
        u=randint(0,len(N)-1)
        Q1=N[u]
        L=E2([0,1,0]).division_points(2^k1)
	M=filter(lambda x : x.order()==2^k1 , L)
	Liste=[]
	for x in M:
	    	y=(2^(k1-k2))*x
		if   E2(y.xy()[0]^p,y.xy()[1]^p)==y :          	
			Liste.append(y)#on stocke dans cette liste car c'est l'ensemble des points rationnels de torsion "maximale"
	NL=[]	
	for y in Liste:#on enlève les points qui sont en double dans la liste
        	a=len(Liste)
              	Liste=filter(lambda x : x != y , Liste)
               	b=len(Liste)
              	if b<a:
	              	NL.append(y)
        Liste=NL
	u=randint(0,len(Liste)-1)
	P2=Liste[u]
	N=filter(lambda x : x.order()==2^k2 , L)
	N=filter(lambda x : x.weil_pairing(P2,2^k2).multiplicative_order()==2^k2, N)
        N=filter(lambda x : E2(x.xy()[0]^p,x.xy()[1]^p)==p*x , N) 
        u=randint(0,len(N)-1)
        Q2=N[u]
	phi1=EllipticCurveIsogeny(E1,2^(k2-1)*Q1)
	phi2=EllipticCurveIsogeny(E2,2^(k2-1)*Q2)
	ListeFinale2=[]
        V2=[]
        L12=[]
        L22=[]
	ListeFinale2=fonction_reccursive_d(E1,phi1,p,P1,Q1,E2,phi2,P2,Q2,k1,k2,Liste,V2,L12,L22)
	PU=ListeFinale2[1]
	E2=ListeFinale2[0][0].curve()
	E1=PU.curve()
	L=E1([0,1,0]).division_points(i)
	L=filter(lambda x: x.order()==i, L)
	c=0
	while c<len(L):
		phi=EllipticCurveIsogeny(E1,L[c])
		if phi.codomain().j_invariant()==E2.j_invariant():
			c=len(L)
			PUI=phi(PU)
		else:
			c=c+1
	ListeFinale2.append("PUI")
	ListeFinale2.append(PUI)
	return ListeFinale, ListeFinale2

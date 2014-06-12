#L'ensemble des fonctions suivantes servent à déterminer la base du Module de Tate de la courbe de départ et les (l^k)* candidats pour l'image de cette base  

#La fonction suivante sert à retourne la matrice de frobenius du point Qf dans la base Qb,Pb; avec Pb et Qb deux points de l^k torsion

def log_discret(Qf,Qb,Pb,p,l,k):#retourne la matrice de frobenius du point Qf dans la base Qb,Pb; avec Pb et Qb deux points de l^k torsion
	j=1	
	for i in range(2,k+1):
		if 2^(k-i)*Qf==(p%l^i)*2^(k-i)*Qb+(j+2^(i-1))*2^(k-i)*Pb:
			j=j+2^(i-1)
	return j

def cantor(f,q,): #algorithme de Cantor-Zassenhaus qui donne à la fin les racines de f, on sait que f se scinde.
	pari('M=factorcantor(g,q); L=[-polcoeff(M[1,1],0,x),-polcoeff(M[2,1],0,x),-polcoeff(M[3,1],0,x),-polcoeff(M[4,1],0,x)]')
	L=(pari('L')).python()
	return L


#La fonction suivante sert à cacluler un unique point de l-division du point el sur la courbe E
def fonction_auxiliaire_test(l,el,E):
        ans=[]
        nel=-el
        P_is_2_torsion = (el == -el)
        g = E._multiple_x_numerator(l) - el[0]*E._multiple_x_denominator(l)

        if P_is_2_torsion:
            if l % 2 == 0:
                g = g.gcd(g.derivative())*g.leading_coefficient().sqrt()

            else:
                g0 = g.variables()[0] - el[0]
                g = g // g0
                g = g.gcd(g.derivative())*g.leading_coefficient().sqrt()
                g = g0*g
        L=g.roots(multiplicities=False); i=0;
        while (len(ans)==0):
            x=L[i]; i=i+1;
            if E.is_x_coord(x):
                Q = E.lift_x(x)
                nQ = -Q
                lQ = l*Q
                if P_is_2_torsion:
                    if lQ == el:
                        ans.append(Q)
                        if nQ != Q:
                            ans.append(nQ)
                else:
                    if lQ == el:
                        ans.append(Q)
                    elif lQ == nel:
                        ans.append(nQ)
        return ans

#La focntion suivante calcule calcule un unique l^k point de division de P (en supposant que P est non nul)

def division_point_new_q(l,k,P,E):#calcule un unique l^k point de division de P
	L=[P];
	for i in range(k):
		L=fonction_auxiliaire_test(l,L[0],E)
	return L

#La focntion suivante calcule deux points de base du module de Tate pour la l^k torsion

def division_point_new_deux(l,k,E):
	M=[]
        g = E.division_polynomial(l)
        for x in g.roots(multiplicities=False):
      	    if E.is_x_coord(x):
      	        Q = E.lift_x(x)
                lQ = l*Q
                if lQ == E(0):
                    M.append(Q)
        L=M; j=1;
        while  L[j].weil_pairing(L[0],l).multiplicative_order()!=l:
            j=j+1
        L=[L[0],L[j]]
	for i in range(1,k):
                    M=[]
                    for el in L:
                        N=fonction_auxiliaire_test(l,el,E)
                        M.append(N[0])
                    L=M		    
	return L

#Cette fonction calcule sur le cratère des candidats pour le module de Tate à la l^k2 torsion

def calcul_ima_iso_1(E1,E2,p,l,k2): #k1, k2 sont les exposants de la ltorsion infinie pour les courbes l-isogènes définies sur p^2, fait la 1ère étape de l'algorithme, ne se déplace que sur le cratère	
	#on calcule les candidats pour P
	L1=division_point_new_deux(l,k2,E1)
        k=1
        if  L1[0][0]^p==L1[0][0]:
            k=0
        elif L1[1][0]^p!=L1[1][0]:
                L1[1]=L1[1]+L1[0]    
        #il faut construire une liste plus grande de NL1
        NL1=[]
        for i in range(l^k2):
            for j in range (l^k2):
                if i%l!=0 and j%l==0:
                    NL1.append(i*L1[k]+j*L1[1-k]) 
	i=0; NL=[]; 	
	while (len(NL)==0):
		y=NL1[i]; i=i+1    		
		E=E1.isogeny_codomain(y)
                a=E.a4(); b=E.a6()
                if a^p==a and b^p==b:
                    M=division_point_new_deux(l,k2-1,E) ; 
		    if M[0][0]^p==M[0][0] and M[1][0]^p==M[1][0]:
                        NL.append(y)
        P1=y

	#on fait (presque) la même chose sur l'autre volcan 

	L2=division_point_new_deux(l,k2,E2)
        h=1 
    	if  L2[0][0]^p==L2[0][0]:
            h=0
        elif L2[1][0]^p!=L2[1][0]:
                L2[1]=L2[1]+L2[0] 
        #il faut construire une liste plus grande de NL2
        NL1=[]
        for i in range(l^k2):
            for j in range (l^k2):
                if i%l!=0 and j%l==0:
                    NL1.append(i*L2[h]+j*L2[1-h]) 
	NL=[]; i=0;	
        while len(NL)==0 :
		y=NL1[i]; i=i+1    		
		E=E2.isogeny_codomain(y)
                a=E.a4(); b=E.a6()
                if a^p==a and b^p==b:
     		    M=division_point_new_deux(l,k2-1,E);
		    if M[0][0]^p==M[0][0] and M[1][0]^p==M[1][0]:
		        NL.append(y)
        NLP2=[]
        for i in range(l^(k2)):
            if i%l !=0:
                NLP2.append(i*NL[0])

	
	#maintenant on va faire la même chose pour Q
	
	#il faut construire une liste plus grande de NL1
        NL1=[]
        if E1(L1[1-k].xy()[0]^p,L1[1-k].xy()[1]^p)!=p*L1[1-k]:
            L1[1-k]=L1[1-k]+L1[k]
        for i in range(l^k2):
            for j in range (l^k2):
                if i%l!=0 and j%l==0:
                    NL1.append(i*L1[1-k]+j*L1[k]) 
	
	i=0; NL=[];
	while(len(NL)==0):
		y=NL1[i]; i=i+1    		
		E=E1.isogeny_codomain(y)
                a=E.a4(); b=E.a6()
                if a^p==a and b^p==b:
     		    M=division_point_new_deux(l,k2-1,E);
		    if M[0][0]^p==M[0][0] and M[1][0]^p==M[1][0]:
	                NL.append(y)
	Q1=y

	#on fait (à nouveau) (presque) la même chose sur l'autre volcan 

	#il faut construire une liste plus grande de NL1
        NL1=[]
        if E2(L2[1-h][0]^p,L2[1-h][1]^p)!=p*L2[1-h]:
            L2[1-h]=L2[1-h]+L2[h]
        for i in range(l^k2):
            for j in range (l^k2):
                if i%l!=0 and j%l==0:
                    NL1.append(i*L2[1-h]+j*L2[h]) 
	
	NL=[];	i=0;
	while len(NL)==0:
		y=NL1[i]; i=i+1    		
		E=E2.isogeny_codomain(y)
     		M=division_point_new_deux(l,k2-1,E)
		if M[0][0]^p==M[0][0] and M[1][0]^p==M[1][0]:
			NL.append(y)
        M=[]
        for i in range(l^(k2)):
            if i%l !=0:
                M.append(i*NL[0])         	
	return [P1,Q1,NLP2,M]


#Cette fonction (filtre) calcule en descendant dans le volcan des candidats pour le module de Tate à la l^k2 torsion

def calcul_ima_iso_2(P,Q,L1,L2,p,l,k,ib):#k étant la hauteur du volcan, i le degré de l'isogénie, cette fonction elle se sert de la hauteur du volcan et descend dans celui-ci pour restreindre les points (couple en l’occurrence ici) candidats
	L=[];
	E1=P.curve();
	phi1=E1.isogeny(l*(P+Q)) #on construit une isogénie descendante
	E1n=phi1.codomain()
	P1n=phi1(P)
        Q1n=division_point_new_q(l,k-1,phi1(P+Q),E1n)[0]
	Qf=E1n(Q1n[0]^p,Q1n[1]^p)
        V=log_discret(Qf,Q1n,P1n,p,l,k)
        b=P.weil_pairing(Q,l^k)^ib
	E2=L1[0].curve()
	# on teste alors sur la seconde courbe si l'action du Frobenius est la même
        for m in range(len(L1)):
            for n in range(len(L2)):        
                P1T=L1[m]
	        Q1T=L2[n]
                a=P1T.weil_pairing(Q1T,l^k)
                if a==b :
	            phi2=E2.isogeny(l*(P1T+Q1T)) #on construit une isogénie descendante
	            P1U=phi2(P1T)
                    E2n=phi2.codomain()
	            Q1U=division_point_new_q(l,k-1,phi2(P1T+Q1T),E2n)[0]
	            for i in range (l^k):
	                if (i%l!=0):
			    l1b=i*P1U
			    l2b=i*Q1U
			    l2f=E2n(l2b[0]^p,l2b[1]^p)			
			    if (l2f== V*l1b + p*l2b ):
                                L.append([i*P1T,i*Q1T] )
        return P,Q,L


#La fonction suivante fait tout le calcul du début à la fin...

def synthese(E1,E2,p,l,k2,i):
    L=calcul_ima_iso_1(E1,E2,p,l,k2)
    Lbis=calcul_ima_iso_2(L[0],L[1],L[2],L[3],p,l,k2,i)
    return L[0], L[1], Lbis

def test_temps_calcul(E1,E2,p,l,k2,i,occurrence):
    t1=t2=t3=0;
    for j in range(occurrence):
        t=cputime();
        L=calcul_ima_iso_1(E1,E2,p,l,k2);
        a=cputime(t)
        t1=t1 + a
        t=cputime();
        L=calcul_ima_iso_2(L[0],L[1],L[2],L[3],p,l,k2,i);
        a=cputime(t)
        t2=t2 + a
        t=cputime();
        L=synthese(E1,E2,p,l,k2,i);
        a=cputime(t)
        t3=t3 + a
    return [t1/occurrence,t2/occurrence,t3/occurrence]

def test_4_eme_cas_val(p): #fonction juste utile pour voir la valuation 2 adique dans le seul cas pas prouve 
	c=0; d=0;	
	while p%8 !=1 :
    		p=next_prime(p)
	K=FiniteField(p)
	k.<a>=GF(p^3)
	for r in range(p):
    		E=EllipticCurve(j=K(r))
    		N=E.cardinality()
    		if N%2==1 :
        		c=c+1
        		E1=construction_lift(E,k,p); N1=E1.cardinality()
        		if valuation(N1,2)%2==0:
				d=d+1
	print c,d,p

def test_complet_4_eme_cas_val(p): #fonction juste utile pour voir la forme de la 2 torsion en haut du cratère id est  E1 
	c=0; g=0;	
	while p%8 !=1 :
    		p=next_prime(p)
	K=FiniteField(p)
	k.<a>=GF(p^3)
	for r in range(p):
    		E=EllipticCurve(j=K(r))
    		N=E.cardinality()
    		if N%2==1 :
        		c=c+1
        		E1=construction_lift(E,k,p); N1=E1.cardinality(); d=valuation(N1,2); e=d/2; L=filter(lambda x : x.order()==2^e, E1(0).division_points(2^e)); P=L[0]; M=filter(lambda x: x.weil_pairing(P,2^e).multiplicative_order()==2^e, L)
        		if len(M)>0 :
				g=g+1; #print r,e,d,"succes"
			else:
				print r,e,d,"echec"
	print c,g,p

def test_forme_cratere(p): #on teste la forme du cratère du volcan de 2 isogenie pour toutes les courbes définies sur GF(p)
	c=0; d=0; e=0; bool=False;
	K=FiniteField(p)
	for i in range(p):
		E=EllipticCurve(j=K(i))
    		N=E.cardinality(); g=valuation(N,2); h=g.quo_rem(2)[0];Eb=E.quadratic_twist(); Nb=Eb.cardinality(); gb=valuation(Nb,2); hb=gb.quo_rem(2)[0];
		if hb>h:
			E=Eb; g=gb; h=hb; N=Nb;
		if h!=0 and N%2==0:
			L=filter(lambda x : x.order()==2^h, E(0).division_points(2^h)); P=L[0]; M=filter(lambda x: x.weil_pairing(P,2^h).multiplicative_order()==2^h, L); 
			if len(M)>0:		
				Q=M[0]
				R=2^(h-1)*P; E1=E.isogeny_codomain(R)
				R=2^(h-1)*Q; E2=E.isogeny_codomain(R)
				R=2^(h-1)*(Q+P); E3=E.isogeny_codomain(R)
				L=filter(lambda x : x.order()==2^h, E1(0).division_points(2^h)); P=L[0]; M=filter(lambda x: x.weil_pairing(P,2^h).multiplicative_order()==2^h, L);
				if len(M)>0:
					c=c+1
					bool=True
				L=filter(lambda x : x.order()==2^h, E2(0).division_points(2^h)); P=L[0]; M=filter(lambda x: x.weil_pairing(P,2^h).multiplicative_order()==2^h, L);
				if len(M)>0:
					if bool==True:
						e=e+1
						d=1
					else:
						c=c+1		
						bool=True
				L=filter(lambda x : x.order()==2^h, E3(0).division_points(2^h)); P=L[0]; M=filter(lambda x: x.weil_pairing(P,2^h).multiplicative_order()==2^h, L);
				if len(M)>0:
					if bool==True:
						e=e+1
						d=1
					else:
						c=c+1
					bool=True
				if bool==True:
					print d,h,N.factor(), E.j_invariant()
				d=0; bool=False;
	print c,e,p;


def identical_rafting(P,E1,phi1,phi2): #cette fonction regarde si le chemin descendant donné par le point P est le même après des images par phi1 et phi2
	E1b1=phi1.codomain()
	E1b2=phi2.codomain()
	P1=phi1(P)
	P2=phi2(P)
	j=E1.isogeny_codomain(P).j_invariant()
	j1=E1b1.isogeny_codomain(P1).j_invariant()
	j2=E1b2.isogeny_codomain(P2).j_invariant()
	print "j",j,"j1",j1,"j2",j2;

def identical_rafting_path(P,E1,L1,L2,p): #cette fonction regarde si le chemin descendant donné par le point P est le même après des images par L1 et L2 listes d'isogénies
	P1=P;
	P2=P;
	for l1 in L1 :
		E1b1=l1.codomain()
		P1=l1(P1)
	for l2 in L2 :
		E1b2=l2.codomain()
		P2=l2(P2)
	P3=E1(P.xy()[0]^p,P.xy()[1]^p)
	j3=E1.isogeny_codomain(P3).j_invariant()
	j=E1.isogeny_codomain(P).j_invariant()
	j1=E1b1.isogeny_codomain(P1).j_invariant()
	j2=E1b2.isogeny_codomain(P2).j_invariant()
	print "j",j,"j1",j1,"j2",j2,"j3",j3;

def identical_frobenius_path(P,Pt,E1,L1,L2,p): #cette fonction regarde si le chemin descendant donné par le point P est le même après des images par L1 et L2 listes d'isogénies
	Pt1=Pt;
	Pt2=Pt;
	P1=P;
	P2=P;
	for l1 in L1 :
		E1b1=l1.codomain()
		Pt1=l1(Pt1)
		P1=l1(P1)
	for l2 in L2 :
		E1b2=l2.codomain()
		Pt2=l2(Pt2)
		P2=l2(P2)
	i0=9;
	i1=9;
	i2=9;
	L11=filter(lambda x: x.xy()[0]^p==x.xy()[0] and x.xy()[1]^p==x.xy()[1],P1.division_points(2));
	P1=L11[randint(0,len(L11)-1)];
	L12=filter(lambda x: x.xy()[0]^p==x.xy()[0] and x.xy()[1]^p==x.xy()[1],P2.division_points(2));
	P2=L12[randint(0,len(L12)-1)];
	L1=filter(lambda x: x.xy()[0]^p!=x.xy()[0] and x.xy()[1]^p!=x.xy()[1],Pt1.division_points(2));
	Pt1=L1[randint(0,len(L1)-1)];
	L2=filter(lambda x: x.xy()[0]^p!=x.xy()[0] and x.xy()[1]^p!=x.xy()[1],Pt2.division_points(2));
	Pt2=L2[randint(0,len(L2)-1)];
	for Pt1 in L1:
		for P1 in L11:
			Pt1f=E1b1(Pt1.xy()[0]^p,Pt1.xy()[1]^p);
			if Pt1f==2*P1+p*Pt1:
				i1=2;
			elif Pt1f==6*P1+p*Pt1:
				i1=6;	
			print "i1",i1;
	Ptf=E1(Pt.xy()[0]^p,Pt.xy()[1]^p);
	if Ptf==2*P+p*Pt:
		i0=2;
	elif Ptf==6*P+p*Pt:
		i0=6;
	for Pt2 in L2:
		for P2 in L12:
			Pt2f=E1b2(Pt2.xy()[0]^p,Pt2.xy()[1]^p);
			if Pt2f==2*P2+p*Pt2:
				i2=2;
			elif Pt2f==6*P2+p*Pt2:
				i2=6;
			print "i2",i2;
	print "i0",i0,"i1",i1,"i2",i2;

def recherche_rationnel_diversifie(q,h):
	K.<a>=GF(q)
        p=K.characteristic()
	L=[]
        if p!=q :
	    for i in range(1,q):
                if (a^i)^p==a^i:
                    
	            E=EllipticCurve(j=K(a^i))
		    N=E.cardinality()
		    t=E.trace_of_frobenius()
		    d=t^2-4*q
		    f=valuation(d,2)
		    if f!=+Infinity:
		        if N/(2^f)!=1  and (d%8==1 or d%8==7) and N%2==0:	        
		    	    I=[N.factor(),a^i]	
                            L.append(I)
        else :
            for i in range(1,q):
                E=EllipticCurve(j=K(i))
		N=E.cardinality()
		t=E.trace_of_frobenius()
		d=t^2-4*q
		f=valuation(d,2)
		if f!=+Infinity:
		    if N/(2^f)!=1 and (d%8==1 or d%8==7) and N%2==0:	    
		    	I=[N.factor(),i]	
                    	L.append(I)
        print len(L)
        return L	

def frobenius_and_rafting(P,E1,phi1,p): #cette fonction regarde si le chemin descendant donné par le point P est le même après des images par phi1 et phi2
	E1b1=phi1.codomain()
	P1=phi1(P)
	P2=E1(P.xy()[0]^p,P.xy()[1]^p)
	j=E1.isogeny_codomain(P).j_invariant()
	j1=E1b1.isogeny_codomain(P1).j_invariant()
	j2=E1.isogeny_codomain(P2).j_invariant()
	print "j",j,"j1",j1,"j2",j2;

def recherche_rationnel_diversifie(q,h):
	K.<a>=GF(q)
        p=K.characteristic()
	L=[]
        if p!=q :
	    for i in range(1,q):
                if (a^i)^p==a^i:
                    
	            E=EllipticCurve(j=K(a^i))
		    N=E.cardinality()
		    t=E.trace_of_frobenius()
		    d=t^2-4*q
		    f=valuation(d,2)
		    if f!=+Infinity:
		        if N/(2^f)!=1  and (d%8==1 or d%8==7) and N%2==0:	        
		    	    I=[N.factor(),a^i]	
                            L.append(I)
        else :
            for i in range(1,q):
                E=EllipticCurve(j=K(i))
		N=E.cardinality()
		t=E.trace_of_frobenius()
		d=t^2-4*q
		f=valuation(d,2)
		if f!=+Infinity:
		    if N/(2^f)!=1 and (d%8==1 or d%8==7) and N%2==0:	    
		    	I=[N.factor(),i]	
                    	L.append(I)
        print len(L)
        return L

def recherche_courbe_diversifie(q,n_i,h):#recherche une courbe sur Fq de cardinal divisible par 2^h et n_i
	K.<a>=GF(q)
        p=K.characteristic()
	L=[]
        if p!=q :
	    for i in range(1,q):
                if (a^i)^p==a^i:
                    
	            E=EllipticCurve(j=K(a^i))
		    N=E.cardinality()
		    if ( N%n_i==0 and N%2^h==0):
		        I=[N.factor(),a^i]
                        L.append(I)
        else :
            for i in range(1,q):
                E=EllipticCurve(j=K(i))
		N=E.cardinality()
		if ( N%n_i==0 and N%2^h==0):
		        I=[N.factor(),i]	
                     	L.append(I)
        print len(L)
        return L	

def etude_end_path_bas(P,Q,L):#décrit l'action de l'endomorphisme de la courbe sur une base donnée
	n=P.order()
	P1=P
	Q1=Q
	a=0;
	b=0;
	for l in L:
		P1=l(P1);
		P=l(P);
		Q1=l(Q1);
		Q=l(Q);
	for l in L:
		P1=l(P1);
		Q1=l(Q1);
	for i in range(1,n):
		if i*P1==P:
			print "i",i;
			a=i;
	for j in range(1,n):
		if j*Q1==Q:
			print "j",j;
			b=j;
	return a,b, P.curve()==P1.curve(), Q1.curve()==Q.curve();	

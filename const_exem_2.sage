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

def test_forme_cratere(p): #on teste la forme du cratère pour toutes les courbes définies sur GF(p)
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

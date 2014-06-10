def construction_lift(E1,k,p): #fonction juste utile pour (bien) construire les exemples
	K=FiniteField(p)
        a=E1.a4()
	b=E1.a6()
        a=K(a)
        b=K(b)
	E2=EllipticCurve(k,[a,b])
	return E2

def construction_courbe_isogene(E1,p,i):#prend en argument le degré i de l'isogénie, la caractéristique p du corps sur lequel est défini la courbe, permet de construire une courbe isogène défini sur le corps de base GF(p)
    L=E1([0,1,0]).division_points(i)
    b=1
    c=1
    while (b!=0):
        E2=EllipticCurveIsogeny(E1,L[c]).codomain(); c=c+1; a=E2.a4() ; a2=E2.a6();
        if (a^p==a and a2^p==a2 and a!=E1.a4() and a2!=E1.a6()):
            b=0
            phi2=EllipticCurveIsogeny(E1,L[c])
    return phi2,E2

def recherche_rationnel(q,h):
	K.<a>=GF(q)
        p=K.characteristic()
	L=[]
        if p!=q :
	    for i in range(1,q):
                if (a^i)^p==a^i:
                    
	            E=EllipticCurve(j=K(a^i))
		    t=E.trace_of_frobenius()
		    d=t^2-4*q
		    if d!=0 and valuation(d,2)>=2*h and valuation(d,2)<+Infinity:
		        I=[E,a^i]	
                        L.append(I)
        else :
            for i in range(1,q):
                E=EllipticCurve(j=K(i))
		t=E.trace_of_frobenius()
		d=t^2-4*q
		if d!=0 and valuation(d,2)>=2*h and valuation(d,2)<+Infinity:
		    I=[E,i]	
                    L.append(I)
        print len(L)
        return L




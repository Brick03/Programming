import math, cmath


def roundto2(x):  #rounding
    sig = 2
    while True:
        if isinstance(x, complex):
            x = x.real
            continue
        elif x != None:
            x = round(x, sig-int(math.floor(math.log10(abs(x))))-1)
            return x
        else:
            return x


def func():  #equasion selection
    while True:
        try:
            funct = int(input("""enter number for line type
1) quadratic (ax**2+bx+c)
2) cubic     (ax**3+bx**2+cx+d)
3) quartic   (ax**4+bx**3+cx**2+dx+e) Do not use! underconstuction!"""))
        except ValueError:
            print('must be from the options above')
            continue
        if funct < 1 or funct > 2:
            print('must be from the options above')
            continue
        else:
            return funct


def ax():  #inputs
    while True:
        try:
            a = float(input('enter value a'))
        except ValueError:
            print('must be number')
            continue
        if a == 0:
            print('cannot be 0')
            continue
        else:
            break
    return a


def bx():
    while True:
        try:
            b = float(input('enter value b'))
        except ValueError:
            print('must be number')
            continue
        else:
            break
    return b


def cx():
    while True:
        try:
            c = float(input('enter value c'))
        except ValueError:
            print('must be number')
            continue
        else:
            break
    return c


def dx():
    while True:
        try:
            d=float(input('enter value d'))
        except ValueError:
            print('must be number')
            continue
        else:
            break
    return d


def ex():
    while True:
        try:
            e = float(input('enter value e'))
        except ValueError:
            print('must be number')
            continue
        else:
            break
    return e


def quadroot(a,b,c):#quadratic calculation
    if (b**2)-(4*a*c) > 1:
        root1 = (-b+math.sqrt((b**2)-(4*a*c)))/(2*a)
        root2 = (-b-math.sqrt((b**2)-(4*a*c)))/(2*a)
        return [root1, root2]
    elif (b**2)-(4*a*c) == 0:
        root1=(-b+math.sqrt((b**2)-(4*a*c)))/(2*a)
        root2=(-b-math.sqrt((b**2)-(4*a*c)))/(2*a)
        return [root1, root2]
    elif (b**2)-(4*a*c) < 0:
        return [None, None]


def cuberoot(a, b, c, d):#cubic calculation
    f = ((3*c/a)-((b**2)/(a**2)))/3
    g = ((2*(b**3)/(a**3))-(9*b*c/(a**2))+(27*d/a))/27
    h = ((g**2)/4)+((f**3)/27)
    if f == g == h == 0:
        root1 = (d/a)**(1/3)
        root2 = (d/a)**(1/3)
        root3 = (d/a)**(1/3)
        return [root1, root2, root3]
    elif h <= 0:
        i = (((g**2)/4)-h)**(1/2)
        j = i**(1/3)
        k = math.acos(-(g/(2*i)))
        l = j*(-1)
        m = math.cos(k/3)
        n = math.sqrt(3)*math.sin(k/3)
        p = (b/3*a)*(-1)
        root1 = (2*j*math.cos(k/3))-(b/3*a)
        root2 = (l*(m+n))+p
        root3 = (l*(m-n))+p
        return [root1, root2, root3]
    elif h > 0:
        r = (-(g/2)+(h**(1/2)))
        s = r**(1./3)
        t = ((g/2)+(h**(1/2)))
        u = -(t)**(1/3)
        u = u.real
        root1 = (s+u)-(b/(3*a))
        root2 = ((s+u)/2)-((cmath.sqrt(-1)*(s-u)*(3**(1/2)))/2)
        root3 = ((s+u)/2)+((cmath.sqrt(-1)*(s-u)*(3**(1/2)))/2)
        return [root1, root2, root3]


def quartroot(ax, bx, cx, dx, ex):  #quartic calculation
    nonzer = []
    zer = []
    a = ax/ax
    b = bx/ax
    c = cx/ax
    d = dx/ax
    e = ex/ax
    f = c-(3*(b**2)/8)
    g = d+((b**3)/8)-(b*c/2)
    h = e-(3*(b**4)/256)+((b**2)*c/16)-(b*d/4)
    s = b/(4*a)
    x = (f/2)
    print(x)
    y = (((f**2)-(4*h))/16)
    print(y)
    z = ((-(g**2))/64)
    print(z)
    cube = cuberoot(1, x, y, z)
    nonzero = -1
    print(cube)
    zero = -1
    while True:
        nonzero += 1
        if cube[nonzero] != 0 and nonzero != 2:
            if isinstance(cube[nonzero], complex):
                nonzer.append(cmath.sqrt(cube[nonzero]))
                continue
            else:
                nonzer.append(math.sqrt(cube[nonzero]))
                continue
        else:
            break
    r=(-(g/(8*nonzer[0]*nonzer[1])))      
    root1=nonzer[0]+nonzer[1]+r-s
    print(root1)
    root2=nonzer[0]-nonzer[1]-r-s
    print(root1)
    root3=nonzer[0]+nonzer[1]-r-s
    print(root1)
    root4=nonzer[0]-nonzer[1]+r-s
    print(root1)
    return [root1, root2, root3, root4]


def tpquad(a,b,c):#quadratic turning/nature
    x = (-b)/(2*a)
    y = (a*(x**2))+(b*x)+c
    w = (a*((x+1)**2))+(b*(x+1))+c
    z = (a*((x-1)**2))+(b*(x-1))+c
    return [x, y, w, z]


def tpcube(a, b, c, d):
    aq = a*3
    bq = b*2
    cq = c
    x = quadroot(aq, bq, cq)
    x1 = x[0]
    x2 = x[1]
    y1 = (a*(x1**2))+(b*x1)+c
    w1 = (a*((x1+1)**2))+(b*(x1+1))+c
    z1 = (a*((x1-1)**2))+(b*(x1-1))+c
    y2 = (a*(x2**2))+(b*x2)+c
    w2 = (a*((x2+1)**2))+(b*(x2+1))+c
    z2 = (a*((x2-1)**2))+(b*(x2-1))+c
    return [x1, y1, w1, z1, x2, y2, w2, z2]


def tpquart(a,b,c,d,e):
    aq=a*4
    bq=b*3
    cq=c*2
    dq=d
    x=cuberoot(aq,bq,cq,dq)
    x1=x[0]
    x2=x[1]
    x3=x[2]
    y1=(a*(x1**2))+(b*x1)+c
    w1=(a*((x1+1)**2))+(b*(x1+1))+c
    z1=(a*((x1-1)**2))+(b*(x1-1))+c
    y2=(a*(x2**2))+(b*x2)+c
    w2=(a*((x2+1)**2))+(b*(x2+1))+c
    z2=(a*((x2-1)**2))+(b*(x2-1))+c
    y3=(a*(x3**2))+(b*x3)+c
    w3=(a*((x3+1)**2))+(b*(x3+1))+c
    z3=(a*((x3-1)**2))+(b*(x3-1))+c
    return [x1,y1,w1,z1,x2,y2,w2,z2,x3,y3,w3,z3]


def stpoint(y,z,w):# determing nature for sationary points
    print(y,z,w)
    minmax=['And there is a Minimum turing point at  (','And there is a Maximium turing point at (','And there is a rising point of infliction at (','And there is a falling point of infliction at (']
    if w>y and z>y:
        return minmax[0]
    elif w<y and z<y:
        return minmax[1]
    elif w<y and z>y:
        return minmax[2]
    elif w>y and z<y:
        return minmax[3]


def quadanswer(a,b,c,x1,x2):#quadratic output
    if x1==x2==None:
        print('For the line',a,'x**2 +',b,'x +',c,'There are no real roots')
    elif x1==x2:
        print('For the line',a,'x**2 +',b,'x +',c,'Both roots are real and equal at x=',x1)
    elif x1!=x2:
        print('For the line',a,'x**2 +',b,'x +',c,'There are real roots at x=',x1,'and x=',x2)
    stcal=tpquad(a,b,c)
    tpdet=stpoint(stcal[1],stcal[2],stcal[3])
    stcal[0]=roundto2(stcal[0])
    stcal[1]=roundto2(stcal[1])
    print(tpdet,stcal[0],',',stcal[1],')')
    print('And a Y axis intercept at Y=',c)
    return


def cubeanswer(a,b,c,d,x1,x2,x3):
    if isinstance(x1, complex) or isinstance(x2, complex) or isinstance(x3, complex):
        x=[x1,x2,x3]
        c=2
        for g in range(3):
            if isinstance(x[c],complex):
                del x[c]
                c=c-1
        print('For the line',a,'x**3 +',b,'x**2 +',c,'x +',d,'There is one real root at x=',x1)
    elif x1==x2==x3==0:
        print('For the line',a,'x**3 +',b,'x**2 +',c,'x +',d,'All roots are real and equal at x=',x1)
    elif cuberoots[0]!=cuberoots[1]!=cuberoots[2]:
        print('For the line',a,'x**3 +',b,'x**2 +',c,'x +',d,'There are 3 real roots at x=',x1,',x=',x2,'and x=',x3)
    stcal=tpcube(a,b,c,d)
    tpdet1=stpoint(stcal[1],stcal[2],stcal[3])
    tpdet2=stpoint(stcal[5],stcal[6],stcal[7])
    n=len(stcal)-1
    for f in range(len(stcal)):
        stcal[n]=roundto2(stcal[n])
        n=n-1
    print(tpdet1,stcal[0],',',stcal[1],') and',tpdet2,stcal[4],',',stcal[5],')')
    print('And a Y axis intercept at Y=',d)
    return


def quartanswer(a,b,c,d,e,x1,x2,x3,x4):
    if isinstance(x1,float) and isinstance(x2,float) and isinstance(x3,float) and isinstance(x4,float):
        x=[x1,x2,x3,x4]
        xnd=[]
        for i in x:
            if i not in xnd:
                xnd.append(i)
        if len(xnd) == 1:
            print('For the line', a, 'x**4', b, 'x**3 +', c, 'x**2 +', d, 'x +', e, 'All roots are real and equal at x=', xnd[0])
        elif len(xnd) == 4:
            print('For the line', a, 'x**4', b, 'x**3 +', c, 'x**2 +', d, 'x +', e, 'There are 4 real roots at x=', xnd[0], ',x=', xnd[1], ',x=', xnd[2], ' and x=', xnd[3])
        elif len(xnd) == 3:
            print('For the line', a, 'x**4', b, 'x**3 +', c, 'x**2 +', d, 'x +', e, 'There are 3 real roots at x=', xnd[0], ',x=', xnd[1], ' and x=', xnd[2])
        elif len(xnd) == 2:
            print('For the line', a, 'x**4', b, 'x**3 +', c, 'x**2 +', d, 'x +', e, 'There are 2 real roots at x=', xnd[0], ' and x=', xnd[1])
    elif isinstance(x1,complex) or isinstance(x2,complex) or isinstance(x3,complex) or isinstance(x4,complex):
        n = 3
        h = [x1,x2,x3,x4]
        for d in range(len(h)):
            if isinstance(h[n],complex):
                del h[n]
                n=n-1
        if len(quartroot)==3:
            print('For the line', a, 'x**4', b, 'x**3 +', c, 'x**2 +', d, 'x +', e, 'There are 3 real roots at x=', h[0], ',x=', h[1], ' and x=', h[2])
        elif len(quartroot)==2:
            print('For the line', a, 'x**4', b, 'x**3 +', c, 'x**2 +', d, 'x +', e, 'There are 2 real roots at x=', h[0], ' and x=', h[1])
        elif len(quartroot)==1:
            print('For the line', a, 'x**4', b, 'x**3 +', c, 'x**2 +', d, 'x +', e, 'There is 1 real roots at x=', h[0])
        else:
            print('For the line', a, 'x**4', b, 'x**3 +', c, 'x**2 +', d, 'x +', e, 'There are no real roots')
    stcal = tpquart(a,b,c,d,e)
    tpdet1 = stpoint(stcal[1],stcal[2],stcal[3])
    tpdet2 = stpoint(stcal[5],stcal[6],stcal[7])
    tpdet3 = stpoint(stcal[9],stcal[10],stcal[11])
    print(tpdet1,stcal[0],',',stcal[1],') and',tpdet2,stcal[4],',',stcal[5],') and',tpdet3,stcal[8],',',stcal[9],')')
    print('And a Y axis intercept at Y=',e)


while True:
    com = func()
    while True:
        if com ==1:
            quad=[]
            quad.append(ax())
            quad.append(bx())
            quad.append(cx())
            break
        elif com ==2:
            cube=[]
            cube.append(ax())
            cube.append(bx())
            cube.append(cx())
            cube.append(dx())
            break
        elif com==3:
            quart=[]
            quart.append(ax())
            quart.append(bx())
            quart.append(cx())
            quart.append(dx())
            quart.append(ex())
            break
        else:
            break
    if com==1:
        while True:
            quadroots=quadroot(quad[0],quad[1],quad[2])
            n=-1
            for a in range(len(quadroots)):
                n=n+1
                quadroots[n]=roundto2(quadroots[n])
            quadanswer(quad[0],quad[1],quad[2],quadroots[0],quadroots[1])
            break
    elif com==2:
        while True:
            cuberoots=cuberoot(cube[0],cube[1],cube[2],cube[3])
            n=-1
            for a in range(len(cuberoots)):
                n=n+1
                cuberoots[n]=roundto2(cuberoots[n])
            cubeanswer(cube[0],cube[1],cube[2],cube[3],cuberoots[0],cuberoots[1],cuberoots[2])
            break
    elif com==3:
        while True:
            quartroots=quartroot(quart[0],quart[1],quart[2],quart[3],quart[4])
            n=-1
            for a in range(len(quartroots)):
                n=n+1
                quartroots[n]=roundto2(quartroots[n])
            quartanswer(quart[0],quart[1],quart[2],quart[3],quart[4],quartroots[0],quartroots[1],quartroots[2],quartroots[3])
            break
    while True:
        try:
            again=str(input('again? y/n'))
        except ValueError:
            print('must be n or y')
            continue
        if again!='y' and again!='n':
            print('must be n or y')
            continue
        else:
            break
    if again=='y':
            continue
    else:
        break

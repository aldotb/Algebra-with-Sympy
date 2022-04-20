from sympy import *
from libaldo_math2 import *


k,x,t=symbols('k x t')
# Trandformacion de monimios and informacion Funciones

def exp_in(ksym): # retorna True si hay una funcion exp en ksym 
                  #  exp_in(4*exp(m)) return True
                  #  exp_in(4*m) return False
    sres=str(ksym)
    if 'exp' in sres:
        return True
    return False
    
def Is_Exp(ksym):
    return exp_in(ksym=ksym)

################################################

    
def get_exp_Mono(ksym):       # retorna la parte exp de un monomio
    if exp_in(ksym):          # get_exp_Mono(4*exp(x)) return exp(x)
        if Is_Mul(ksym):       
            mm=fpoly(ksym,'list')
            for i in mm:
                if exp_in(i):
                    return i
    return ksym  

def get_Euler (ksym):
    return get_cofac_exp(ksym=ksym) 

##########################################    
    
def get_cofac_exp(ksym):      # retorna la parte  no exp de un monomio
    kexp=get_exp_Mono(ksym)   # get_cofac_exp (4*exp(x)) return 4
    if kexp!=ksym:
        return simplify(ksym/kexp)
    return 1 

def get_facEuler(ksym):
    return get_cofac_exp(ksym=ksym)

####################################


def get_expoinexp(ksym):
    kres=ksym
    if exp_in(ksym):
        mm=get_exp_Mono(ksym)
        mm2=fpoly(mm,'list')
        kres=mm2[0]
        return kres
    return 1    

def get_expnEuler(ksym):
    return get_expoinexp(ksym=ksym)


#  LOG

def log_in(ksym): # retorna True si hay una funcion log en ksym 
                  #  log_in(4*exp(m)) return False
    sres=str(ksym)              #  log_in(log(4)*m) return True
    if 'log' in sres:
        return True
    return False 

def get_log_Mono(ksym):       # retorna la parte exp de un monomio
    if log_in(ksym):          # get_exp_Mono(4*exp(x)) return exp(x)
        if Is_Mul(ksym):       
            mm=fpoly(ksym,'list')
            for i in mm:
                if log_in(i):
                    return i
    return ksym

def get_cofac_log(ksym):      # retorna la parte  no exp de un monomio
    kexp=get_log_Mono(ksym)   # get_cofac_exp (4*exp(x)) return 4
    if kexp!=ksym:
        return simplify(ksym/kexp)
    return ksym 

def get_insidelog(ksym):      #  
    kexp=get_log_Mono(ksym)   #  
    return simplify(exp(kexp)) 

#####################################
#           Monomie Tools 
#####################################
    
def monodata(*args):
    r'''
       data_mono(monomie, op1, op2)
        inputs args:
           op1=='isexp' ,return True if ksym is monomie whit exp factor
           op1=='getexp' ,return only exp factor 
           op1=='getexpfac' ,return all factors less exp fac

    '''
    ksym=args[0]
    if len(args)==2:
        op1=args[1]
    if len(args)==3:
        op2=args[1]
     
    # Exp
    if op1=='isexp':
        return exp_in(ksym)

    elif op1=='getexp':
        return get_exp_Mono(ksym)
    elif op1=='getexpfac':
        return get_cofac_exp(ksym)
    elif op1=='getexpoinexp':
        return get_expoinexp(ksym)

    # Log
    if op1=='islog':
        return log_in(ksym)

    elif op1=='getlog':
        return get_log_Mono(ksym)
    elif op1=='getlogfac':
        return get_cofac_log(ksym)
    elif op1=='getinlog':
        return get_insidelog(ksym)
        
        
def clear_exp_QQ(ksym1,ksym2,vmain=k,var2=t):
    p2=simplify(ksym2)
    p1=simplify(ksym1)
    if Is_Add(p2) and Is_Poly(p2):
        mm=0
        mmexp=0
        for i in fpoly(p2,'list'):
            if exp_in(i):
                mmexp+=i
            else:
                mm+=i
        p2=mmexp
        p1=p1-mm
        if Is_Mono(p2) and denom(p2)!=1:
            kdeno=denom(p2)
            p1=p1*kdeno
            p2=simplify(p2*kdeno)
            
        mm=get_facEuler(p2)
        if mm!=1:
            p1=p1/mm
            p2=simplify(p2/mm)
            
            
        return  (p1,p2) 

     
     
#####################################
#             MyEqEq 
##################################### 

def lexpand(p1,p2):   
        p1=factor(p1)
        p2=factor(p2)
        
        n1=numer(p1)
        d1=denom(p1)
        n2=numer(p2)
        d2=denom(p2)
        p1=expand(n1*d2)
        p2=expand(n2*d1)
        

        return p1,p2
        


#####################################
#             str
#####################################
        
def between_par(sall,ssym=''):  # between_par('***(12(**)3)***') return '12(**)3'
    mm=sall
    if ssym=='':
        ssym='('
    q1=mm.find(ssym)
    if q1==-1:
        return ''
    mm2=mm[q1::]
    
    q2=mm2.find('(')
    mm3=mm2[q2::]
    
    cc=1
    sres=''
    p=1
    i=1
    while cc>0:
        v=mm3[i]
        if v=='(':
            cc+=1
        if v==')':
            cc-=1
        if cc>0:    
            sres+=v
        i+=1
    return sres 

def between_func(ksym,sfunc):
    '''
        input a function like 
            ksym=sqrt(y + 4)*sqrt((x - 2)**2) type: symbols o function
            sfunc='sqrt(' type=str
        return:
            vec ['y + 4', '(x - 2)**2'] ,all inside sqrt(*****) in tihs case
    
    '''
    svar=str(ksym)
    qq=svar.count(sfunc)
    vecf=[]
    posx=0
    for i in range(qq):
        ff=between_par(svar,sfunc)
        vecf.append(ff)
        posx=svar.find(sfunc)
        svar=svar[posx+1::]
    return vecf 

#####################################
#     polynomie Tools
#####################################


#   squaresum()
#   -----------
def squaresum(ksym,p1=0,p2=0):
    '''
    sqauresum(ksym,p1=x,p2=1)
        input :
            ksym = initial function to trasform
            p1,p2  = (p1+p2)**2 new formate to get convert
        return:
            (p1+p2)**2 + K ..such that (p1+p2)**2 + K = ksym
        example:
            the polynomie 4*x*x+8*x-6 tray transfor to (2*x+2)**2 +K
            
                squaresum(4*x*x+8*x-6,2*x,2) return (2*x + 2)**2 - 10
    '''
 
    vec1=expand((p1+p2)**2)
    vec2=ksym-vec1
    kres=kpow((p1+p2),2)+vec2
    return kres

def parcial_fraction(x,L,a,b,c,d):
    r'''
    input L/((a*x+b)*(c*x+d))
    
    return  A/(a*x+b) + B/(c*x+d)
    '''
    A=L*a/(d*(a - c))
    B=-L*c/(d*(a - c))
    p1=(a*x+b)
    p2=(a*x+b)
    
    kres=A/p1+B/p2
    return kres
    
def par_frac(ksym,var=''):  
    r'''
    FRACCION PARCIAL
    input C/[A+B]
    return if it possible C1/A + C2/A que es igula a  C/[A+B]
    '''
   
    
    kres=ksym
    
    if Is_Mono(kres) and Is_Mul(kres):
         
        p1=numer(kres) 
        done=False
        if var!='':
            X=symbols('X')
            kres=ksym.subs(var,X)
            done=True
             
         
        p2=denom(kres)
         

        mm= fpoly(p2,'list')
        
        qq=len(mm)
        A1,A2=symbols('A1 A2')
        f1=A1*mm[0]
        f2=A2*mm[1]
        kres2=simplify(expand(f1+f2))
        kres3=factorSec(kres2,X)
         
        A22=solve(parse_expr(between_par(str(kres3),'X*(')),A2)[0]
        kres4=kres3.subs(A2,A22)
        r1=solve(1-kres4,A1)[0]
        r2=A22.subs(A1,r1)
        dd1=mm[0]
        dd2=mm[1]
         
        if done:
            dd1=dd1.subs(X,var)
            dd2=dd2.subs(X,var)
         
        fac1=p1*r1/dd1
        fac2=p1*r2/dd2
        return fac1+fac2
 
    else: ksym 

 
def real_subs(QQ,**kwargs):
    ''' 
    QQ= symbols function
    ** kwargs c7=6,g4=z..etc..
    RETURN real substitucion when variable have underscore name like 'c_7' 'g_4'    
    '''
    
    vvar=fpoly(QQ,'free')
    mvar=[]
    for i in vvar:
        nname=str(i)
        nname=nname.replace('_','')
        nvar=symbols(nname)
        mvar.append(nvar)
         
    kres=QQ 
    sres=str(kres)
    sres=sres.replace('_','')
    kres=parse_expr(sres)
    mkey=[]
    vvalue=[]
    for key, value in kwargs.items():
        mkey.append(key)
        vvalue.append(value)
        
        kres=kres.subs(parse_expr(key),value)
    for i,j in zip(mvar,vvar):
        kres=kres.subs(i,j)
    return (kres)

def coef_list(ksym,var2=x):  # return coef list complete according max degree
        '''
        ee=x*x*a+x*b+c
        return [a,b,c]

        '''
 
         
        dd=degree(ksym,gen=var2)
        vcoef=[]
        for i in range(1,dd+1):
            pe=dd+1-i
            xfac=var2**pe
            
            vcoef.append(ksym.coeff(xfac))
        kres=ksym.subs(var2,0)
        vcoef.append(kres)
        return vcoef

def sortdegree(ksym,var2=x):
    '''
      group polynomies respect main varible
    '''
    mm=coef_list(ksym,var2)
    kres=0
    qq=len(mm)
    cc=qq-1
    for i in mm:
        kres+=i*var2**cc
        cc+=-1
    return kres
    
def algebra_compare_coef(*args,var2=x):
    '''
    function that resuelve variables when compare 2 polinomies wits same grade and same coef    
    input= args:
            function,or,symbols
            gen= function with defaul gen=x
    output: solution all variables possibles
        
            
    '''
    p1=''
    p2=''
    vsym=[]
    for i in args:
        if Is_Symbol(i):
                vsym.append(i)
        else:
            if p1=='':
                p1=i 
            else:
                p2=i        
                
    
     
    ql1=coef_list(p1,var2=var2)
    ql2=coef_list(p2,var2=var2)
    qq=[]
    for i,j in zip(ql1,ql2):
        qq.append(i-j)
    kres=solve(qq,vsym)    
    if type(kres)==dict:
        kres= kunpakDic(kres)
        nvar=kres[0]
        nval=kres[1]
    if type(kres)==list:
        kres= kunpakDic(kres)
        nvar=kres[0]
        nval =kres[1]
 
    return nval



    
#####################################
#           functions 
#####################################    
def subsnumber(ksym,val1,val2): # replace an integer number by symbol 
    kres=ksym
    sres=str(kres)
    sv1=str(val1)
    sv2=str(val2)
    sres=sres.replace(sv1,sv2)
    kres2=parse_expr(sres)
    return kres2
    
def func2symbol(ksym,y,x):
    veri='('+str(x)+')'
    sf=str(y)
    if veri in sf:
        nf=sf.replace(veri,'')
    else:
        nf=sf
    oldf=nf+veri
    sres=str(ksym)
    sres=sres.replace(oldf,nf)
    return parse_expr(sres)
    
    
#def findFunc(f1,vecf,vecv):
def findFunc(f1,vecf,vecv):
    '''
     f((x+1)**2)= 2*x*x+4*x+5 , f1= 2*x+3
     procces, x--> Add(1) -->['A'],[1]  -->x+1  -->Pow(2)  -->['A','P'],[1,2]  -->  
     vecf = ['A','P'] ,Posibles...>add,'A',mul  'M',divide 'D',sqrt 'R',pow 'P'
     vecv = ['A','P']
     findFunc(2*x*x+4*x+5,['A','P'],[1,2])
         
     return f(x)= 2*x+3   
    '''
    y,z=symbols('y z')
    #Q=MyEqEq(y,f1)
    xx=solve(y-f1,x)
    if type(xx)==list:
        xx=xx[0]
    for i,j in zip(vecf,vecv):
        if i=='R':
            xx=sqrs(xx,j)
        elif i=='M':
            xx=xx*j
        elif i=='D':
            xx=xx/j
        elif i=='P':
            xx=xx**j
        else:
            xx+=j
    yy=solve(z-xx,y)
    if type(yy)==list:
        yy=yy[0]
    y=yy.subs(z,x)

    
    return y    
    
#####################################
#           polynomies tools
#####################################

#   simplifysum()
#   -----------        
def simplifysum(pval):
    '''
    return sum of each monomie of pval but uniformatin each one to get real sum
    

    '''
    if Is_Add(pval):
        mm=fpoly(pval,'list')
     
        kk=0
        for i in mm:
            k1=simplify(i)
            kk+=mulexp(k1)
        return kk
    else:
        return pval


#   linfactor()
#   -----------        
def linfactor(ksym,kvar=''):
    '''
    factorize() return personality factorization of Polynomie summ 
    example:
        if kfunc = a*x*x+b*c*x*x+3*x
        factorize(kfunc,x*x) 
        return :
        x*x*(x+b*x)+3*x    

    '''

    p1=ksym
    if kvar!='' and type(p1)==Add:
        oldv=0
        newv=0
        kfac=kvar
        mm=fpoly(p1,'list')
        for i in mm:
            vec=simplify(i/kfac)
            if denom(vec)==1:
                newv+=vec
            else:
                oldv+=i
        kres=oldv+newv*kfac
        return kres
    else:
        return ksym



#   factorize()
#   -----------        
def factorize(kfunc,val):
    '''
    factorize() return personality factorization of Polynomie summ 
    example:
        if kfunc = a*x*x+b*c*x*x+3*x
        factorize(kfunc,x*x) 
        return :
        x*x*(x+b*x)+3*x    

    '''

    if type(kfunc)==Add:
        p1=0    
        p2=0
        p3=0
        mm=fpoly(kfunc,'list')
        facf=val
        v2=fpoly(facf,'free')
        for i in mm:
            nff=simplify(i/facf)
            v1=fpoly(nff,'free')
            if rem(i,nff)==0:
                p2+=nff
            else: 
                p1+=i 
             
               
         
              
        p3=p1+val*(p2)
        return p3
    else:
        return kfunc
        
def killsqrtpow(ksym):
    if type(ksym)==Add:
        kres=0
        mm=fpoly(ksym,'list')
        for i in mm:
            kres+=killsqrtpow(i)
        return kres
    
    elif type(ksym)==Mul:
        kres=1
        mm=fpoly(ksym,'list')
        for i in mm:
            kres=kres*killsqrtpow(i)
        return kres
        
    elif type(ksym)==Pow:
        kres=ksym
        try:
            mm=fpoly(ksym,'list')
            base=mm[0]
            vexp=mm[1]
            try:
                vexp=killsqrtpow(vexp)
                kres=ppow(base,vexp)
            except:
                pass
        except:
            pass
            
        kres= kill_RootPow(ksym=kres)
        if type(kres)==str:
            return parse_expr(kres)
        else:
            return kres
    else:
        try:
            kres=ksym.subs(insidepar(ksym,'sqrt('),killsqrtpow(insidepar(ksym,'sqrt(')))
            return kres
        except:    
            return ksym
        
def kill_RootPow(ksym): # kill rootPow in Polynomie one pass
    if type(ksym)!=str:  
        kres=str(ksym)
    else:
        kres=ksym
    posible=True
 
    pclave='sqrt(' 
    pbusca='**2'
    inpar=between_par(kres,pclave)
     
    if inpar=='':
        return kres
    else:
        inbusca=inpar[-3::]
        if pbusca==inbusca:
             
            var22=inpar.replace(pbusca,'')
            var22
             
            ptotal=pclave+inpar+')'
            ptotal
             
            kres=kres.replace(ptotal,var22)
         
    return parse_expr(kres)    

#   inside_par()
#   -----------        
def insidepar(ksym,ssym=''):
    '''
        args:
            ksym = function
            ssym = str( where search..., 'x(', 'sqrt(', 'log(' etc 
        return expresion inside  select function     
    '''
    sres=str(ksym)
    if ssym=='':
        kres=between_par(sres)
    else:
        kres=between_par(sres,ssym)
    fres=parse_expr(kres)
    return fres    
    
def Area2funcXY(ee1, ee2, var=x, cota='', kope=''):
    if type(ee1)==MyEq:
        kres1=ee1.ksym
    else:
        kres1=ee1    
    if type(ee2)==MyEq:
        kres2=ee2.ksym
    else:
        kres2=ee2    
    kres=solve(kres1-kres2,var)
    if cota=='' and len(kres)==2:
        segx=kres
    else:
        x1,x2=cota
        segx=[x1]
        for i in kres:
            if i>x1 and i<x2:
                segx.append(i)
        segx.append(x2)
        
    vecf=[]
    vecps=0
    qq=len(segx)
    for i in range(qq-1):
        ee=MyEq(ee1-ee2,'',x1=segx[i],x2=segx[i+1],var2=var,ktype='I',kshow=False)
        vecf.append(ee)
        vecps+=ee 
     
    ee2=MyEq(vecps,'Area')
    kres=0
    
    NoSYMB=False
    for i in vecf:
        if Is_Symbol(i.ksym):
            NoSYMB=False
            
    
    for i in vecf:
        value=i.ksym
        if NoSYMB:
            kres+=abs(value.doit())
        else:
            kres+=value.doit()
    ee3=MyEq(kres,'Area')     
    
def get_rem(kres1,kres2,kname='',var2=''): # return rest of div e1/e2
        if type(kres1)==MyEq:
            kres1=kres1.ksym
        if type(kres2)==MyEq:
            kres2=kres2.ksym
        kres=  rem(kres1,kres2)  
        if kname!='':
            
            return MyEq(kres,kname=kname,var2=var2) 
        else:    
            return kres

def get_quo(kres1,kres2,kname='',var2=''): # return quotient  of div e1/e2
        if type(kres1)==MyEq:
            kres1=kres1.ksym
        if type(kres2)==MyEq:
            kres2=kres2.ksym
        kres= quo(kres1,kres2) 
        if kname!='':
            return MyEq(kres,kname=kname,var2=var2) 
        else:    
            return kres
def get_cofactors(kres1,kname='',var2=''):  # return cofactors  of div e1/e2
        if type(kres1)==MyEq :
            kres1=kres1.ksym
        if type(kres2)==MyEq:
            kres2=kres2.ksym
        kres= cofactors(kres1,kres2) 
         
        if kname!='':
            return MyEq(kres,kname=kname,var2=var2) 
        else:    
            return kres
            
            
def get_seudocofactor(e2,e3,var2):
    '''
    tyr to get polynomie complement and multipli complete degree secuence
    e1=x*x-2
    e2=x+1
        if e1=e2*k
        then k=a*x+b
        maybe x*x-2= (x+1)*(a*x+b)
        
    return (a*x+b),[a,b]
    '''
    vecvar=[a,b,c,d]
    vecfvar=[]
    ee2=e2
    if type(e2)==MyEq:
        ee2=e2.ksym
        
    ee3=e3
    if type(e3)==MyEq:
        ee3=e3.ksym    
    qq=degree(ee2,gen=var2)-degree(ee3,gen=var2)
    vres=0
    cc=0
    for i in range(qq+1):
        vres+=vecvar[i]*var2**(qq-cc)
        vecfvar.append(vecvar[i])
        cc+=1
    return (vres,vecfvar)
        
#####################################
#           geometry tools
#####################################    
def acircle(R,kname=''):
    if type(R)==MyEq:
        R=R.ksym
    kres=pi*R*R    
    if kname!='':
        return MyEq(kres,kname)
    return kres

def dvolume(L1,L2, kname=''):
    if type(L1)==MyEq:
            L1=L1.ksym 
    if type(L2)==MyEq:
            L2=L2.ksym
    kres=L1*L2    
    if kname!='':
        return MyEq(kres,kname)
    return kres 

def mpoint(L1,L2, kname=''):
    if type(L1)==MyEq:
            L1=L1.ksym 
    if type(L2)==MyEq:
            L2=L2.ksym
    kres=L1/L2    
    if kname!='':
        return MyEq(kres,kname)
    return kres     
    
#####################################
#           solver
#####################################     
    
def solver(*args):
    r'''
    input MyEq and variables
    return solve variables
    '''
    vee=[]
    vval=[]
    sval=[]
    
    for i in args:
        if Is_Symbol(i):
            vval.append(i)
            sval.append(str(i))
            
             
        else:
            vee.append(i)
             
    if len(vee)!=len(vval):
        print ('insuficient data...')
        return
    mres=solve(vee,vval)        
    kres=[]
    if type(mres)==dict:
        for key, value in mres.items():
             
             
            kres.append(value)
    else:
        mres=list(mres[0])
        for i,j in zip(mres,sval):
             
            kres.append(i)
            
        
    return kres   


#####################################
#           list
#####################################  



def one_ksym2ksym(ksym):  
    r'''
    return denominator if numerator =1
    input 1/(a+b) return (a+b)
    input 3/(a+b) return  3/(a+b)
    input c/(a+b) return c/(a+b)
    ''' 
    
    if numer(ksym)==1:
        return denom(ksym)
    else:
        return ksym

def num_ksym2ksym(ksym):
    r'''
    return denominator if numerator = numeric
    input 1/(a+b) return (a+b)
    input 3/(a+b) return  (a+b)
    input c/(a+b) return c/(a+b)
    '''
    if Is_Number(numer(ksym)):
        return denom(ksym)
    else:
        return ksym    
        
#####################################
#           simplify  Sqrt(Pow(2))
#####################################          
def simplifyrpow(ksym):
    r'''
     
    input sqrt(y + 4)*sqrt((x - 2)**2)
    return  (x - 2)*sqrt(y + 4)
     
    ''' 
    
    svec= between_func(ksym,'sqrt(')
    svar=str(ksym)
    kres=ksym
    if svec!=[]:
        qq=len(svec)
        svar=str(ksym)
        for i in svec:
            if i[-3::]=='**2':
                p1='sqrt('+i+')'
                p2='('+i[0:-3]+')'
                svar=svar.replace(p1,p2)
                kres=parse_expr(svar)
    return kres 


#####################################
#           conjuntos
##################################### 
def someinvec(vec1,vec2): 
    '''
    return True if some items in vec 2 are in vec1
    return false if nothing items in vec 2 are in vec1
    '''
    for i in vec2:
        if i in vec1:
            return True
    return False    
    




#####################################
#           exponents
##################################### 
 
def factorexp(ksym):
    if Is_Pow(ksym):
        v1,e1=base_exponent(ksym)
        e2=opemat(e1,'f')
        kres=v1**e2
        return kres
    elif Is_Mul(ksym):
        mm= fpoly(ksym,'list')
        kres=1
        for i in mm:
            kres=kres*factorexp(i)
        return kres    
    elif Is_Add(ksym):
        mm= fpoly(ksym,'list')
        kres=0
        for i in mm:
            kres=kres+factorexp(i)
        return kres    
        
    elif Is_Div(ksym):
        p1=factorexp(numer(ksym))
        p2=factorexp(denom(ksym))
        return p1/p2
            
    
    else:
         
        return ksym
        

def base2frac(ksym): # return nmonomie whit expand exponents
     
    if Is_Pow(ksym):
        kres=ksym
        mm=fpoly(ksym,'list')
        vval=mm[0]
        vexp=mm[1]
        
        if denom(vval)==1:
            vexp=-1*vexp
            kres=(1/vval)**vexp
            return kres
        else:
            return ksym
         
    elif Is_Div(ksym):
        p1=packexp(base2frac(ksym))
        p2=packexp(base2frac(ksym))
        return p1/p2
            
    elif Is_Mul(ksym):
        mm= fpoly(ksym,'list')
        kres=1
        for i in mm:
            kres=kres*base2frac(i)
        return kres    
    elif Is_Add(ksym):
        mm= fpoly(ksym,'list')
        kres=0
        for i in mm:
            kres=kres+base2frac(i)
        return kres    
    else:
        kres=ksym
        return kres
        
def packexp(ksym): # return nmonomie whit expand exponents
     
    if Is_Pow(ksym):
        kres=ksym
        mm=fpoly(ksym,'list')
        vval=mm[0]
        vexp=mm[1]
        if type(vexp)==Add:
            kres=1
            mm1=fpoly(vexp,'list')
            for i in mm1:
                kres=kres*(vval**i)
            return kres
        else:
            return ksym
         
    elif Is_Div(ksym):
        p1=packexp(numer(ksym))
        p2=packexp(denom(ksym))
        return p1/p2
            
    elif Is_Mul(ksym):
        mm= fpoly(ksym,'list')
        kres=1
        for i in mm:
            kres=kres*packexp(i)
        return kres    
    elif Is_Add(ksym):
        mm= fpoly(ksym,'list')
        kres=0
        for i in mm:
            kres=kres+packexp(i)
        return kres    
    else:
        kres=ksym
        return kres        
        
        
def expandexp(ksym,op=''): # return nmonomie whit expand exponents
     
    if Is_Pow(ksym):
        if op=='e':
            mm=fpoly(ksym,'list')
            p1=mm[0]
            p2=mm[1]
            kres=expandexp(p2,op='')
            kres=p1**kres
            return kres
        else:    
            kres=ksym
            mm=fpoly(ksym,'list')
            vval=mm[0]
            vexp=mm[1]
            vexp2=expand(vexp)
            kres=kpow(vval,vexp2)
            return kres
    elif Is_Div(ksym):
        p1=expandexp(numer(ksym))
        p2=expandexp(denom(ksym))
        return p1/p2
            
    elif Is_Mul(ksym):
        mm= fpoly(ksym,'list')
        kres=1
        for i in mm:
            kres=kres*expandexp(i)
        return kres    
    elif Is_Add(ksym):
        mm= fpoly(ksym,'list')
        kres=0
        for i in mm:
            kres=kres+expandexp(i)
        return kres    
    else:
        kres=ksym
        return kres
    
def simplifybase(ksym,kope=''):
         
        if type(ksym)==Pow:
            mm=fpoly(ksym,'list')
            p1=mm[0]
            p2=mm[1]
            p1=simplify(p1)
            kres=p1**p2
            if kope!='':
                kres=opemat(kres,kope=kope)
            return kres
        elif type(ksym)==Add:
            kres=0
            mm=fpoly(ksym,'list')
            for i in mm:
                kres+=simplifybase(i)
            if kope!='':
                kres=opemat(kres,kope=kope)
            return kres
        elif type(ksym)==Mul:
            kres=1
            mm=fpoly(ksym,'list')
            for i in mm:
                kres*=simplifybase(i)
            if kope!='':
                kres=opemat(kres,kope=kope)
            return kres    
        else:
            return ksym
                
def simplifyexp(ksym,kope=''):
    if '**' not in str(ksym):
        return ksym
    if Is_Add(ksym):
        mm=fpoly(ksym,'list')
        kres=0
        for i in mm:
            kres+=simplifyexp(i,kope=kope)
        return kres
    elif  Is_Pow(ksym):
        mm=fpoly(ksym,'list')
        val=mm[0]
        vexp=mm[1]
        if Is_Pow(val):
            mm2=fpoly(val,'list')
            val2=mm2[0]
            vexp2=mm2[1]
            nexp=vexp2*vexp
            nexp=opemat(nexp,kope=kope)
            return kpow(val2,nexp)
        else:
            nexp=simplify(vexp)
            nexp=opemat(nexp,kope=kope)
            return kpow(val,nexp)
            
        
    elif Is_Mul(ksym):
        try:
            vv=fpoly(ksym,'free')
            mm0=fpoly(ksym,'list')
            qq=len(vv)
            sexp=mzero(qq)
            kres=1
            mm=[]
            for i in mm0:
                if Is_Number(i):
                    kres=kres*i
                else:
                    mm.append(i)

            for i in mm:
                if Is_Pow(i):
                    mm2=fpoly(i,'list')
                    var=mm2[0]
                    nvexp=mm2[1]
                    knum=vv.index(var)
                    sexp[knum]=sexp[knum]+nvexp
                else:
                    knum=vv.index(i)
                    sexp[knum]=1

            for i,j in zip(vv,sexp):
                je=j
                if kope!='':
                    je=opemat(j,kope=kope)
                kres=kres*i**je
            return kres
        except:
            return ksym
    else:
        return ksym 

def powexpand(ksym,op=''):
        '''
        input (x**(a*b))   ---->   return(x**a)**b
        input (x**(a*b),b)   ---->   return(x**b)**a
        '''
        kres=ksym 
        if Is_Number(ksym):
            
            return ksym
        
        
        if Is_Div(ksym):
             
            p1=powexpand(numer(ksym))
            p2=powexpand(denom(ksym))
             
            return p1/p2 
             
        
        elif Is_Mul(ksym):
            mm=fpoly(ksym,'list')
            kres=1
            for i in mm:
                kres=kres*powexpand(i,op=op)
             
            return kres
             
        elif type(ksym)==Pow:
            mm=fpoly(ksym,'list')
            val=mm[0]
            vexp=mm[1]
            if type(vexp)==Pow:

                return ksym
            elif type(vexp)==Mul:
                mm2=fpoly(vexp,'list')
                if len(mm2)==2:
                    p1=mm2[0]
                    p2=mm2[1]
                    if op!='':
                        kres=(val**p2)
                        kres=kres**p1
                    else:    
                        kres=(val**p1)
                        kres=kres**p2
                    return kres
                
            else:
                return ksym
             

             
        elif Is_Add(ksym):
            mm=fpoly(ksym,'list')
            mmr=0
            for i in mm:
                mmr+=powexpand(i,op=op)
            return mmr
        else:
            return ksym
               
def div2mulexp(ksym):
    '''
        input ((a/b)**x   ---->   return(a**x)*(b**(-x))
         
    '''
    if Is_Div(ksym):
        p1=numer(ksym)
        p2=denom(ksym)
        kres=p1*simplify(p2**(-1))
        return kres
    if Is_Pow(ksym):
        mm=fpoly(ksym,'list')
        vald=mm[0]
        vale=mm[1]
        if denom(vald)!='1':
            p1=numer(vald)
            p2=denom(vald)
            kres=(p1**vale)*(p2**(-1*vale))
            return kres
        else:
            return ksym
    if Is_Add(ksym):
        kres=0
        mm=fpoly(ksym,'list')
        for i in mm:
            kres+=div2mulexp(i)
        return kres
    if Is_Mul(ksym):
        kres=1
        mm=fpoly(ksym,'list')
        for i in mm:
            kres=kres*div2mulexp(i)
        return kres
    else:
        return ksym

def primefactorexp(ksym):
    '''
    input 256**x 
    return 2**(8*x)
    '''
    
    if Is_Pow(ksym):
        mm=fpoly(ksym,'list')
        val=mm[0]
        vexp=mm[1]
        kres=1
        if Is_Number(val):
            vecn,vece=kunpakDic(factorint(val))
             
            for i,j in zip(vecn,vece):
                kres=kres*(i**(j+vexp))
        return kres
          
                
    elif Is_Mul(ksym):
        kres=1
        mm=fpoly(ksym,'list')
        for i in mm:
            kres=kres*primefactorexp(i)
        return kres
    elif Is_Add(ksym):
        kres=0
        mm=fpoly(ksym,'list')
        for i in mm:
            kres=kres+primefactorexp(i)
        return kres
    else:
        return ksym 
        
def mulexpo(ksym):
    '''
    input (x**a)**b
    return x**(a*b)
    '''
    if Is_Pow(ksym):
        b1,e1=base_exponent(ksym)
        if Is_Pow(b1):
            b2,e2=base_exponent(b1)
            kresp=e1*e2
            kres=b2**kresp
            return kres
        elif Is_Mul(b1):
            kres=1
            mm=fpoly(b1,'list')
            for i in mm:
                kres=kres*mulexpo(i**e1)
            return kres    
        else:
            return ksym
    elif  Is_Add(ksym):
        ares=0
        mm=fpoly(ksym,'list')
        for i in mm:
            ares=ares+mulexpo(i)
        return ares
    elif  Is_Mul(ksym):
        ares=1
        mm=fpoly(ksym,'list')
        for i in mm:
            ares=ares*mulexpo(i)
        return ares
    else:
        return ksym      
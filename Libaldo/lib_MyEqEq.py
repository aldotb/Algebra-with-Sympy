from sympy import *
 
from IPython.display import Math  # ,display
from matplotlib.pyplot import ylabel, plot, show, xlabel, title
from libaldo_math2 import *
from libaldo_algorith import *

from  lib_MyEq import * 
from libaldo_show import *
 
 
import copy 
x,t,k,y=symbols('x t k y')


# MQ() = short call MyEqEq()
def MQ( *args, val2='', kname='', kshow=True,var2=t,vmain=x,ktype='eq',eQshow='',vfunc=[],kvista=1):  
    return MyEqEq(*args, val2=val2, kname=kname, kshow=kshow,var2=var2,vmain=vmain,ktype=ktype ,eQshow=eQshow ,vfunc=vfunc,kvista=kvista)
 

class MyEqEq:
    def __init__(self, *args, val2='', kname='', kshow=True,var2=t,vmain=x,ktype='eq',eQshow='',vfunc=[],kvista=1):
        self.var2=var2
        self.var0=var2 
         
        if len(args)==1:
            eq=args[0]
            if type(eq)==Equality:
                self.e1=MyEq(eq.lhs,'e1',kshow=False,vfunc=vfunc)
                self.e2=MyEq(eq.rhs,'e2',kshow=False,vfunc=vfunc)
            else:
                if type(eq)!=MyEq:
                    self.e1=MyEq(eq.name,'e1',kshow=False,vfunc=vfunc)
                    self.e2=MyEq(eq.ksym,'e2',kshow=False,vfunc=vfunc)
                     
        else :
            exp1=args[0]
            exp2=args[1]
            self.vmain= vmain 
            self.var2=var2
            if type(exp1)==MyEq:
                self.e1=MyEq(exp1.ksym,'e1',kshow=False,vfunc=vfunc)
            else:
                self.e1=MyEq(exp1,'e1',kshow=False,vfunc=vfunc)
                
            if type(exp2)==MyEq:
                self.e2=MyEq(exp2.ksym,'e2',kshow=False,vfunc=vfunc)
            else:
                self.e2=MyEq(exp2,'e2',kshow=False,vfunc=vfunc)
            
            if len(args)>2:
                self.vmain= args[2] 
            if len(args)>3:
                self.var2=args[3]     
              
         
          
        self.ksym = self.e1.ksym - self.e2.ksym
        self.eQ = Eq(self.e1.ksym , self.e2.ksym)
        self.v = self.e1.ksym - self.e2.ksym
        self.name = kname
        self.type=ktype
         
        #self.f, self.df, self.d2f,self.ff = symbols2diff(self.vmain,self.var2,kshow=False)
        #self.f, self.df, self.d2f,self.ff = symbols2diff(self.vmain,self.var2,kshow=False)
        self.ft=Function(str(self.vmain))(var2)
        self.dfn='d'+str(vmain)
        self.d2fn='d2'+str(vmain)
        self.Vics=[]
        self.ics=''
        self.primieQ=''
        self.modeNormal=True
        self.eQshow=eQshow
        npf=str(vmain)+"'"
        np2f=str(vmain)+"''"
        self.pf=symbols(npf)
        self.p2f=symbols(np2f)
         
        self.dfunc=[]
        self.vista=kvista
        self.vfunc=vfunc 
            
        if str(self.vmain) in str(self.e1.ksym) and str(self.vmain) not in str(self.e2.ksym):
            self.type='so'
        if kshow:
             
            if self.name == '':
                display(Math(latex(self.eQ)))
            else:
                sR = self.name + ')..  '
                display(Math(sR + latex(self.eQ)))
        self.var0=symbols(str(self.vmain))
        
    def __call__(self,*args, **kwargs):
    
        QQ=self.xcopy()
         
        p1=QQ.left
        p2=QQ.right
        p11=real_subs(p1,**kwargs)
        p22=real_subs(p2,**kwargs)
        QQ.e1.ksym=p11
        QQ.e2.ksym=p22
        QQ.s()
    def __repr__(self):
        kres = self.eQ
        return kres
        
    def _latex(self, obj):
        return latex(self.eQ)    
    def __str__(self):
         
        return str(self.__repr__())
    
    def up_sides(self,p1='',p2=''):
        self.up_ee(p1=p1,p2=p2)
        
    def up_ee(self,p1='',p2=''):
        if p1!='':
            self.e1.ksym=p1
        if p2!='':
            self.e2.ksym=p2
         
    ###########################################
    #               Update                    #
    ###########################################

    def __add__(self, other):
        p1=self.e1.ksym
        p2=self.e2.ksym
        if type(other)==MyEqEq:
            p1=p1+other.left
            p2=p2+other.right
        elif type(other)==MyEq:
            p1=p1+other.ksym
            p2=p2+other.ksym

        else:
            p1=p1+other 
            p2=p2+other 
         
        return Eq(p1,p2)

    def __radd__(self, other):
        p1=self.e1.ksym
        p2=self.e2.ksym
        if type(other)==MyEqEq:
            p1=p1+other.left
            p2=p2+other.right
        elif type(other)==MyEq:
            p1=p1+other.ksym
            p2=p2+other.ksym

        else:
            p1=p1+other 
            p2=p2+other 
         
        return Eq(p1,p2)

    def __sub__(self, other):
        p1=self.e1.ksym
        p2=self.e2.ksym
        if type(other)==MyEqEq:
            p1=p1-other.left
            p2=p2-other.right
        elif type(other)==MyEq:
            p1=p1-other.ksym
            p2=p2-other.ksym

        else:
            p1=p1-other 
            p2=p2-other 
         
        return Eq(p1,p2)

    def __rsub__(self, other):
        p1=self.e1.ksym
        p2=self.e2.ksym
        if type(other)==MyEqEq:
            p1=p1-other.left
            p2=p2-other.right
        elif type(other)==MyEq:
            p1=p1-other.ksym
            p2=p2-other.ksym

        else:
            p1=p1-other 
            p2=p2-other 
         
        return Eq(p1,p2)
    
    def __mul__(self, other):
        p1=self.e1.ksym
        p2=self.e2.ksym
        if type(other)==MyEqEq:
            p1=p1*other.left
            p2=p2*other.right
        elif type(other)==MyEq:
            p1=p1*other.ksym
            p2=p2*other.ksym

        else:
            p1=p1*other 
            p2=p2*other 
         
        return Eq(p1,p2)
    
    def __rmul__(self, other):
        p1=self.e1.ksym
        p2=self.e2.ksym
        if type(other)==MyEqEq:
            p1=p1*other.left
            p2=p2*other.right
        elif type(other)==MyEq:
            p1=p1*other.ksym
            p2=p2*other.ksym

        else:
            p1=p1*other 
            p2=p2*other 
         
        return Eq(p1,p2)
        
        
    def __truediv__(self, other):
        p1=self.e1.ksym
        p2=self.e2.ksym
        if type(other)==MyEqEq:
            p1=p1/other.left
            p2=p2/other.right
        elif type(other)==MyEq:
            p1=p1/other.ksym
            p2=p2/other.ksym

        else:
            p1=p1/other 
            p2=p2/other 
         
        return Eq(p1,p2)

    def __rtruediv__(self, other):
        p1=self.e1.ksym
        p2=self.e2.ksym
        if type(other)==MyEqEq:
            p1=p1/other.left
            p2=p2/other.right
        elif type(other)==MyEq:
            p1=p1/other.ksym
            p2=p2/other.ksym

        else:
            p1=p1/other 
            p2=p2/other 
         
        return Eq(p1,p2)        
    
    #               Update                    #
    ########################################### 
    
    def add_ics(self,var1,var2,var3):
        Vics=self.Vics
        
         
        if 'd2' in var1:
             
            arma=d2ficc(self.vmain,self.var2,var2,var3)
            Vics.append(arma)
        elif  'd' in var1 and 'd2' not in var1 :
             
            arma=dficc(self.vmain,self.var2,var2,var3)
            Vics.append(arma)
        else:
            v1=str(self.vmain)+'('+str(var2)+')'
            arma=  fics(v1,var3)
            Vics.append(arma)
        self.Vics=Vics
        self.arma_ics()
        sE(['ics=',self.ics])
    
    def arma_ics(self):
        kres='{'
        Vics=self.Vics
        for i in Vics:
            kres+=i+','
        kres=kres[0:-1]
        kres+='}'
        self.ics=parse_expr(kres)
    ################################## 
    #      Apariencia
    ################################## 
    def s2(self):
    
        self.up_eQshow()
        func=self.func
        dfunc=self.dfunc
        qq=len(func)
        kres=self.eQshow
        for i in range(qq):
            kfunc=func[i]
            kdfunc=dfunc[i]
            kres=kres.subs(kfunc[0],kdfunc[0])
            #res=kres.subs(kfunc[0],kdfunc[0])
        self.eQshow=kres
        self.s('2')
        
    def s3(self):
         
        func=self.func
        dfunc=self.dfunc
        qq=len(func)
        kres=self.eQshow
        for i in range(qq):
            kfunc=func[i]
            kdfunc=dfunc[i]
            kres=kres.subs(kfunc[1],kdfunc[1])
            #res=kres.subs(kfunc[0],kdfunc[0])
        self.eQshow=kres
        self.s('2')    
    
    def save_eQshow(self,keQ):
        self.eQshow=keQ
        
    def up_eQshow(self):
        self.eQshow=self.eQ
    
        
    def up2primitiva(self):
        eq=self.eQ
        self.e1.ksym=eq.lhs
        self.e2.ksym=eq.rhs
        self.modeNormal=False
        
        
    
    def up2normalMode(self):
        return self.primitiva2func() 
    
    def primitiva2func(self):
        try:
            self.e1.ksym=flat_diff(self.left,self.df)
        except:
            pass
        vmain=self.vmain
        vmain=antiprimitiva(vmain)
        nx=symbols(str(vmain))
        fmain=Function(str(vmain))(self.var2)
        if self.modeNormal:
            self.primieQ=Eq(self.e1.ksym,self.e2.ksym)
            self.modeNormal=False
        
        try:
            ee1=self.e1
            ee1.set(fmain,nx,kshow=False)
            self.e1=ee
        except:
            pass
            
        try:
            ee2=self.e2
            ee2.set(fmain,nx,kshow=False)
            self.e2=ee2
        except:
            pass
         
            
        self.s()
        
    def flat_diff(self):
        try:
            self.e1.ksym=self.e1.ksym.subs(Q.df,flat_diff(Q.df))
        except:
            pass
        self.s()
        
    def update(self):

        self.eQ = Eq(sydem(self.e1.ksym), sydem(self.e2.ksym))
         
    def update2  (self,nEq):
        kres=nEq
        self.e1.ksym=kres.lhs
        self.e2.ksym=kres.rhs
        self.update() 
    
    def reformate(self,p1,p2,kshow=True):
        self.e1.ksym=p1
        self.e2.ksym=e2
        self.s(kshow=kshow)
        
    ################################## 
    #         Show()
    ##################################
    
    def set_vista(self,ktype):
        self.vista=ktype
    def swap(self):
        ee=self.e2
        self.e2=self.e1
        self.e1=ee
        self.s()
        
    def s(self,ktype='',kshow=True):
         
        # algorithmo backup
        self.e1.backup.append(self.left)
        self.e2.backup.append(self.right)
        self.update()
        keQ=self.eQ
        if ktype=='':
            vtype=self.vista
        else:
            vtype=ktype
             
        if vtype==2:
             
            p1=self.e1.get_dview()
            p2=self.e2.get_dview()
            keQ=Eq(p1,p2)
        if vtype==3:
            p1=self.e1.get_pview()
            p2=self.e2.get_pview()
            keQ=Eq(p1,p2)    
        if kshow:    
            if self.name == '':
                display(Math(latex(keQ)))
            else:
                sR = self.name + ')..'
                display(Math(sR + latex(keQ)))

    def xcopy(self,op=''):
        kres=copy.deepcopy(self)
        if op!='':
            kres.s()
            
        return kres

    @property
    def rval(self,**args):  # return lhs from Eq
        return self.e2.ksym

    @property
    def right(self):  # return lhs from Eq
        return self.e2.ksym    
    
    @property
    def R(self):
        return self.right
        
    @property
    def L(self):
        return self.left
        
    @property
    def lval(self):  # return lhs from Eq
        return self.e1.ksym
    @property
    def left(self):  # return lhs from Eq
        return self.e1.ksym    

    def setL(self,ksym='',kshow=True,**kwargs):
        '''
            if ksym!='' then self.left will be ksym
            else try to replace kwars in left
        '''
        if ksym!='':
            self.e1.ksym=ksym
        else:   
            if len(kwargs)>0: 
                p1=self.left
                p11=real_subs(p1,**kwargs)
                self.e1.ksym=p11
        if kshow:
            self.s()
    def setR(self,ksym='',kshow=True,**kwargs):
        '''
            if ksym!='' then self.right will be ksym
            else try to replace kwars in right
        '''
        if ksym!='':
            self.e2.ksym=ksym
        else:   
            if len(kwargs)>0: 
                p2=self.right
                p22=real_subs(p2,**kwargs)
                self.e2.ksym=p22
        if kshow:
            self.s()    
    ################################## 
    #         set()
    ##################################

    def set_eQshow(self,*args,kshow=True):
        qq=len(args)
        ks1=[]
        ks2=[]
        
        for i in range(qq):
            if i%2==0:
                ks1.append(args[i])
            else:
                ks2.append(args[i])
        
        eQshow=self.eQshow
        for i,j in zip(ks1,ks2):
            eQshow=eQshow.subs(i,j)
        self.eQshow=eQshow    
        
        if kshow:
            self.s('2')
        
        
    def replace(self, v1,v2):
        kres1=self.left
        kres1=kres1.subs(v1,v2)
        self.e1.ksym=kres1
        kres2=self.right
        kres2=kres2.subs(v1,v2)
        self.e2.ksym=kres2
        self.s()
    
    def evalue_if(self,*args,kshow=True,kope='',**kwargs):
        QQ=self.xcopy('QQ')
        margs=args
        mkwargs=kwargs
        QQ.set(*margs,kshow=kshow,kope=kope,**mkwargs)
    
    def evalueif(self,*args,kshow=True,kope='',**kwargs):
        QQ=self.xcopy('QQ')
        margs=args
        mkwargs=kwargs
        QQ.set(*margs,kshow=kshow,kope=kope,**mkwargs)
    
    def set(self,*args,kshow=True,kope='',**kwargs):
        if len(kwargs)>0:
            p1=self.left
            p1=real_subs(p1,**kwargs)
             
            p2=self.right
            p2=real_subs(p2,**kwargs)
            
            if kope!='':
                p1=opemat(p1,kope=kope)
                p2=opemat(p2,kope=kope)
            
            self.e1.ksym=p1
            self.e2.ksym=p2
        
            if kshow:
                self.s()
        else:         
            
            
            if len(args)==2:
                var1=str(args[0])
                var2=str(args[1])        
                

                p1=str(self.left)
                p1=p1.replace(var1,var2)
                try:
                    self.e1.ksym=parse_expr(p1)
                except:
                    pass
             
                p2=str(self.right)
                p2=p2.replace(var1,var2)
                try:
                    self.e2.ksym=parse_expr(p2)
                except:
                    pass    
                if kshow:    
                    self.s()
                return
            qq=len(args)
            vsym=[]
            vval=[]
            if qq%2==0:
                op2='LR'
            else:
                qq-=1
                op2=args[qq]
                
            for i in range(qq):
                if i%2 == 0:
                    vsym.append(args[i])
                else:
                    vval.append(args[i])
            kres=self.eQ
            for i,j in zip(vsym,vval):
                kj=j
                if type(kj)==MyEq:
                    kj=kj.ksym
                    
                kres=kres.subs(i,kj)
            if 'L' in op2:    
                self.e1.ksym=opemat(kres.lhs,kope=kope)
            if 'R' in op2:    
                self.e2.ksym=opemat(kres.rhs,kope=kope)
             
             
            if kshow:
                self.s(kshow)
        

        
    def undo(self):
        

        try:
            qq=len(self.e1.backup)
            kres12=self.e1.backup[qq-2]
            kres22=self.e2.backup[qq-2]
            self.e1.backup=self.e1.backup[0:qq-2]
            self.e2.backup=self.e2.backup[0:qq-2]
            self.e1.ksym=kres12
            self.e2.ksym=kres22
            self.s()
        except:

            pass
         
    ################################## 
    #         simplify
    ##################################       
    def expand(self, kop='RL',kshow=True):
        if 'L' in kop:
            self.e1.expand(kshow=False)
        if 'R' in kop:
            self.e2.expand(kshow=False)
        self.s(kshow)
    
    def lexpand(self):
        p1,p2=self.left,self,right
        
        self.e1.ksym,self.e2.ksym = lexpand(p1,p2)
        self.s()
        
        
    def Expand(self,kshow=True,kope=''):
        p1=self.left
        p2=self.right
        p1=factor(p1)
        if Is_Mono(p1):
            p2=p2*denom(p1)
            p2=factor(p2)
            p1=numer(p1)
            
        if Is_Mono(p2):
            p1=p1*denom(p2)
            p2=numer(p2)
            
        p1=opemat(p1,kope=kope)
        p2=opemat(p2,kope=kope)
        self.e1.ksym=p1
        self.e2.ksym=p2
        self.s()
            
        
    def opemat(self,kope='',op='LR'):
        p1=self.left
        p2=self.right
        
        if 'L' in op:
            p1=opemat(p1,kope=kope)
        if 'R' in op:
            p2=opemat(p2,kope=kope)    
        self.e1.ksym=p1
        self.e2.ksym=p2
        self.s()    

    def exp_alone(self):
        eQs=str(self.eQ)
        sres=''
        try:
            smm=between_par(eQs,'exp(')
        except:
            self.s()
            return
        smm='exp('+smm+')'
        try:
            mm=parse_expr(smm)
        except:
            self.s()
            return
        self.alone(mm)
             
    def simplify(self,kop='RL',kope='',kshow=True):
        
        p1=self.left
        p2=self.right
        
        if kop!='RL':
            if kop=='L':
                if kope!='':
                    p1=opemat(p1,kope)
                self.e1.ksym=simplify(p1)
            if kop=='R':
                if kope!='':
                    p2=opemat(p2,kope)
                self.e2.ksym=simplify(p2)
            if kshow:
                self.s()
        else:
            if type(p1)==Add and type(p2)==Add:
                m1=fpoly(p1,'list')
                m2=fpoly(p2,'list')
                kres1=p1
                kres2=p2
                for i in m1:
                    if i in m2:
                        kres1=kres1-i
                        kres2=kres2-i
                if p1!=kres1:
                    self.e1.ksym=kres1
                    self.e2.ksym=kres2
                    if kshow:
                        self.s()
            elif (type(p1)==Mul and p2!=0) or (type(p2)==Mul and p1!=0):
                kres1=p1
                kres2=p2
                fcc=gcd(kres1,kres2)
                self.e1.ksym=opemat(kres1/fcc,kope=kope)
                self.e2.ksym=opemat(kres2/fcc,kope=kope)
                if kshow:
                    self.s()
            else:
                if kope!='':
                    p1=opemat(p1,kope)
                    p2=opemat(p2,kope)
                if type(p1)==Mul:
                    p1=simplify(p1)
                if type(p2)==Mul:
                    p2=simplify(p2)    
                self.e1.ksym=p1
                self.e2.ksym=p2
                if kshow:
                    self.s()
            
    def expandexp(self,kop='LR',kshow=True):
        op=''
        if 'e' in kop:
            op='e'
        p1=self.left
        if 'L' in kop:
            p1=expandexp(p1,op=op)
            self.e1.ksym=p1
        p2=self.right
        if 'R' in kop:
            p2=expandexp(p2,op=op)
            self.e2.ksym=p2
         
        if kshow:
            self.s()
            
            
    def simplifybase(self,kop='LR',kshow=True):

        p1=self.left
        if 'L' in kop:
            p1=simplifybase(p1)
            self.e1.ksym=p1
        p2=self.right
        if 'R' in kop:
            p2=simplifybase(p2)
            self.e2.ksym=p2
         
        if kshow:
            self.s()
    def powexpand(self,kop='LR',kshow=True):
        '''
        input (x**(a*b))   ---->   return(x**a)**b
        input (x**(a*b),b)   ---->   return(x**b)**a
        '''
        op=''
        if 'i' in kop:
            op='i'
        p1=self.left
        if 'L' in kop:
            p1=powexpand(p1,op=op)
            self.e1.ksym=p1
        p2=self.right
        if 'R' in kop:
            p2=powexpand(p2,op=op)
            self.e2.ksym=p2
         
        if kshow:
            self.s()

    def mulexpo(self,kop='LR',kshow=True):
        

        p1=self.left
        if 'L' in kop:
            p1=mulexpo(p1)
            self.e1.ksym=p1
        p2=self.right
        if 'R' in kop:
            p2=mulexpo(p2)
            self.e2.ksym=p2
         
        if kshow:
            self.s()
    
    def simplifyexp(self,op='LR',kope='',kshow=True):
        p1=self.left
        if 'L' in op:
            p1=simplifyexp(p1,kope=kope)
            self.e1.ksym=p1
        p2=self.right
        if 'R' in op:
            p2=simplifyexp(p2,kope=kope)
            self.e2.ksym=p2

        if kshow:
            self.s()

    def simplifyexp(self,op='LR',kope='',kshow=True):
        p1=self.left
        if 'L' in op:
            p1=simplifyexp(p1,kope=kope)
            self.e1.ksym=p1
        p2=self.right
        if 'R' in op:
            p2=simplifyexp(p2,kope=kope)
            self.e2.ksym=p2

        if kshow:
            self.s() 
    def div2mulexp(self,op='LR',kope='',kshow=True):
        p1=self.left
        if 'L' in op:
            p1=div2mulexp(p1)
            if kope!='':
                p1=opemat(p1,kope=kope)
            self.e1.ksym=p1
            
        p2=self.right
        if 'R' in op:
            p2=div2mulexp(p2)
            if kope!='':
                p2=opemat(p2,kope=kope)
            self.e2.ksym=p2

        if kshow:
            self.s()

            
    def simplifyrpow(self,kope='',kop='RL',kshow=True):
        

        if 'L' in kop :
            self.e1.simplifyrpow(kshow=False)
        if 'R' in kop :
            self.e2.simplifyrpow(kshow=False)
        if kope!='':
            kres1=self.e1.ksym
            kres2=self.e2.ksym
            
            kres1=opemat(kres1,kope=kope)
            kres2=opemat(fcc,kope=kope)
            self.e1.ksym=kres1
            self.e2.ksym=kres2
        self.s(kshow)
    def simplify_cero(self,kope=''):
        kres=self.left-self.right
        kres=opemat(kres,kope=kope)
        self.e1.ksym=kres
        self.e2.ksym=0
        self.s()
        
    def factor(self, kop='RL',kshow=True,kope=''):
        if 'L' in kop:
            self.e1.factor(kope=kope,kshow=False)
        if 'R' in kop:
            self.e2.factor(kope=kope,kshow=False)
        self.s(kshow)

    def tsimplify(self, kop='RL',kshow=True,kope=''):
        if 'L' in kop:
            self.e1.tsimplify(kope=kope,kshow=False)
        if 'R' in kop:
            self.e2.tsimplify(kope=kope,kshow=False)
        self.s(kshow)

    def tfactor(self, kop='RL',kshow=True,kope=''):
        if 'L' in kop:
            self.e1.tfactor(kope=kope,kshow=False)
        if 'R' in kop:
            self.e2.tfactor(kope=kope,kshow=False)
        self.s(kshow)

    def Add(self, kval,kname='', kop='RL',kshow=True,kope=''):
         
        p1=self.left
        p2=self.right
         
        if type(kval)==MyEqEq:
            k1=kval.left
            k2=kval.right
        elif type(kval)==MyEq:
            k1=kval.ksym
            k2=kval.ksym
        else:	
            k1=kval 
            k2=kval 	    
                 
        if 'L' in kop:
                p1=p1+k1
        if 'R' in kop:
                p2=p2+k2
        if kope!='':
            p1=opemat(p1,kope=kope)
            p2=opemat(p2,kope=kope)
        if kname!='':
            QQ=MyEqEq(p1,p2)
            return QQ
        else:
            self.e1.ksym=p1
            self.e2.ksym=p2
            if kshow:
                self.s()
        
    def Substrac(self, kval,kname='', kop='RL',kshow=True,kope=''):
        p1=self.left
        p2=self.right
         
        if type(kval)==MyEqEq:
            k1=kval.left
            k2=kval.right
        elif type(kval)==MyEq:
            k1=kval.ksym
            k2=kval.ksym
        else:	
            k1=kval 
            k2=kval 	    
                 
        if 'L' in kop:
                p1=p1-k1
        if 'R' in kop:
                p2=p2-k2
        if kope!='':
            p1=opemat(p1,kope=kope)
            p2=opemat(p2,kope=kope)
        if kname!='':
            QQ=MyEqEq(p1,p2)
            return QQ
        else:
            self.e1.ksym=p1
            self.e2.ksym=p2
            if kshow:
                self.s()
        
    def Mul(self, kval,kname='', kop='RL',kshow=True,kope=''):
        p1=self.left
        p2=self.right
         
        if type(kval)==MyEqEq:
            k1=kval.left
            k2=kval.right
        elif type(kval)==MyEq:
            k1=kval.ksym
            k2=kval.ksym
        else:	
            k1=kval 
            k2=kval 	    
                 
        if 'L' in kop:
                p1=p1*k1
        if 'R' in kop:
                p2=p2*k2
        if kope!='':
            p1=opemat(p1,kope=kope)
            p2=opemat(p2,kope=kope)
        if kname!='':
            QQ=MyEqEq(p1,p2)
            return QQ
        else:
            self.e1.ksym=p1
            self.e2.ksym=p2
            if kshow:
                self.s()

    def Div(self, kval,kname='', kop='RL',kshow=True,kope=''):
        p1=self.left
        p2=self.right
         
        if type(kval)==MyEqEq:
            k1=kval.left
            k2=kval.right
        elif type(kval)==MyEq:
            k1=kval.ksym
            k2=kval.ksym
        else:	
            k1=kval 
            k2=kval 	    
                 
        if 'L' in kop:
                p1=p1/k1
        if 'R' in kop:
                p2=p2/k2
        if kope!='':
            p1=opemat(p1,kope=kope)
            p2=opemat(p2,kope=kope)
        if kname!='':
            QQ=MyEqEq(p1,p2)
            return QQ
        else:
            self.e1.ksym=p1
            self.e2.ksym=p2
            if kshow:
                self.s()

    def Pow(self, kval,kname='', kop='RL',kshow=True,kope=''):
        p1=self.left
        p2=self.right
        if 'L' in kop:
             
            p1=kpow(p1,kval)
             
        if 'R' in kop:
             
            p2=kpow(p2,kval)
        if kope!='':
            p1=opemat(p1,kope=kope)
            p2=opemat(p2,kope=kope)
            
        if kname!='':
            QQ=MyEqEq(p1,p2)
            return QQ
        else:
            self.e1.ksym=p1
            self.e2.ksym=p2
            if kshow:
                self.s()

    def Rpow(self, kval,kname='', kop='RL',kshow=True,kope=''):
        p1=self.left
        p2=self.right
        if 'L' in kop:
             
            p1=sqrs(p1,kval)
            if kval==2:
                p1=simplifyrpow(p1)
             
        if 'R' in kop:
             
            p2=sqrs(p2,kval)
            if kval==2:
                p2=simplifyrpow(p2)
        if kope!='':
            p1=opemat(p1,kope=kope)
            p2=opemat(p2,kope=kope)
            
        if kname!='':
            QQ=MyEqEq(p1,p2)
            return QQ
        else:
            self.e1.ksym=p1
            self.e2.ksym=p2
            if kshow:
                self.s()
        
    def lexpand(self,kop='RL',kshow=True):
        if 'L' in kop:
            self.e1.lexpand(kshow=False)
        if 'R' in kop:
            self.e2.lexpand(kshow=False)
        self.s(kshow)
        
    def lfactor(self, kop='RL',kshow=True):
        if 'L' in kop:
            self.e1.lfactor(kshow=False)
        if 'R' in kop:
            self.e2.lfactor(kshow=False)
        self.s(kshow)
    
    def exp(self, kop='RL',kshow=True):
        if 'L' in kop:
            self.e1.exp(kshow=False)
        if 'R' in kop:
            self.e2.exp(kshow=False)
        self.s(kshow)
    
    def log(self, kop='RL',kshow=True):
        if 'L' in kop:
            self.e1.log(kshow=False)
        if 'R' in kop:
            self.e2.log(kshow=False)
        self.s(kshow)
    
    def factorSec(self, ksym,kop='RL',kshow=True):
        if 'L' in kop:
            
            self.e1.factorSec(ksym,kshow=False)
        if 'R' in kop:
            self.e2.factorSec(ksym,kshow=False)
        self.s(kshow)
        
    def linfactor(self, ksym,kop='RL',kshow=True):
        if 'L' in kop:
            p1=self.left
            p1=linfactor(p1,ksym)
            self.e1.ksym=p1
        if 'L' in kop:
            p2=self.right
            p2=linfactor(p2,ksym)
            self.e2.ksym=p2        
        if kshow:
            self.s()
        
    def list(self,kop='R'):
        kres=fpoly(self.right,'list')
        if kop=='L':
            kres=fpoly(self.left,'list')
        return kres
        
        
    def clearUp(self,vv,kshow=True):
        p1=0
        p2=0
        for i in self.e1.list():
            mm=fpoly(i,'free')
            if vv in mm:
                p1+=i
            else:
                p2-=i
        for i in self.e2.list():
            mm=fpoly(i,'free')
            if vv in mm:
                p1-=i
            else:
                p2+=i
        self.e1=MyEq(p1,'e1',kshow=False)
        self.e2=MyEq(p2,'e2',kshow=False)
        if kshow:
            self.s()
        
    def clearAlone(self,vv):
        self.clearUp(vv=vv,kshow=False)
        mm=self.e1.list()
        done=True
        if len(mm)==2:
            if type(self.e1.ksym)==Mul:
                if vv in self.e1.list():
                    kfac=1
                    for i in self.e1.list():
                        if i!=vv:
                            self.Div(i,kshow=False)
                            self.s()
                            done=False
        if done:
            self.s()
        
    
    def sin2cos(self, angu, korden=2, kope='', kop='RL'):
        if 'L' in kop:
            self.e1.set(kpow(sin(angu), 3), (1 - kpow(cos(angu), 2)) * cos(angu), kshow=False)
            kres = self.e1.ksym
            kres = sin2cos(kres, angu=angu, korden=korden, kope=kope)
            self.e1.update(kres)
        if 'R' in kop:
            self.e2.set(kpow(sin(angu), 3), (1 - kpow(cos(angu), 2)) * cos(angu), kshow=False)
            kres = self.e2.ksym
            kres = sin2cos(kres, angu=angu, korden=korden, kope=kope)
            self.e2.update(kres)
        self.s()

    def cos2sin(self, angu, korden=2, kope='', kop='RL'):
        if 'L' in kop:
            self.e1.set(kpow(cos(angu), 3), (1 - kpow(sin(angu), 2)) * sin(angu), kshow=False)
            kres = self.e1.ksym
            kres = cos2sin(kres, angu=angu, korden=korden, kope=kope)
            self.e1.update(kres)
        if 'R' in kop:
            self.e2.set(kpow(cos(angu), 3), (1 - kpow(sin(angu), 2)) * sin(angu), kshow=False)
            kres = self.e2.ksym
            kres = cos2sin(kres, angu=angu, korden=korden, kope=kope)
            self.e2.update(kres)

            self.s()
            
#########################################
#   MyEqEq      evalue 
########################################## 

    def evalue(self,val1,val2,kshow=True):
        QQ=self.xcopy()
        QQ.e1.ksym=QQ.left.subs(val1,val2)
        QQ.e2.ksym=QQ.right.subs(val1,val2)
        
        if kshow:
            QQ.s()
            
        else:
            QQ.s(kshow=False)
            return QQ.right

    
#########################################
#   MyEqEq      solve  
##########################################
    def solve_coef_list(self,*args,var2=x):
        r'''
        solve variable from  two polinomies with tha same coefficient ,grade orden
        
        example:
        **********
        a*x*x+b*x+c= 3*x*x+2*x+7+y
         
        return:  MyEq
        a=3
        b=2
        c=7+y

        '''
        vecv=[]
        vecq=[]
        for i in args:
            vecv.append(i)
        m1=coef_list( self.left,var2)
        m2=coef_list( self.right,var2)
        for i,j in zip(m1,m2):
            vecq.append(i-j)
        vecq=vecq+vecv
        kres=solver(*vecq)
        eres=[]
        for i,j in zip(vecv,kres):
            ee=MyEq(j,str(i))
            eres.append(ee)
        return eres

    def solve(self,*args,kshow=True,**kwargs):
        QQ=self.xcopy()
        sres1=str(QQ.e1.ksym)
        sres2=str(QQ.e2.ksym)
         
        if len(kwargs)>0:    
            for key, value in kwargs.items():
                    sres1=sres1.replace(key,str(value))
                    sres2=sres2.replace(key,str(value))
                
        nval=args[0]
        sval=str(nval)
        if sres1==sval and sval not in sres2 :
            ee=MyEq(sydem(QQ.right) ,kname=sval,var2=QQ.var2,kshow=kshow)
            return ee
        if sres2==sval and sval not in sres1 :
            ee=MyEq(sydem(QQ.left) ,kname=sval,var2=QQ.var2,kshow=kshow)
            return ee
        QQ.e1.ksym=parse_expr(sres1)
        QQ.e2.ksym=parse_expr(sres2)
        for i in args:
            if type(i)==MyEq:
                QQ.upgrade(i,kshow=False)
         
        nval=args[0]
        
        kname=str(nval)
        QQ.alone(nval,kshow=False)
        ee=MyEq(QQ.right,kname,var2=QQ.var2,kshow=kshow)
        return ee
        
        
    def solve_if_and(self, svar, eqv=0, kope='',korden='',kshow=True, **kwargs):
        r'''
        solve variable from  MyEq
        parameters :
            svar :type str , variablesin side the Eq taht we will find
            eqv  :type nemeric or symbols , if the value of all Eq
                  defaul Eq=0
            kwargs: t=0,g=10... etc
        return MyEq of svar
        example:
        **********
        R(t)= C1 + C2*t + g*sin(t*w)/w**2
        C1= solve_if_and('C1',L,t=0)
        return:  C1=L

        R.upgrade(C1)
        return:  C2*t + L + g*sin(t*w)/w**2

        C2=solve_if_and(R,'C2',t=2)
        return:-L/2 - g*sin(2*w)/(2*w**2)

        R.upgrade(C2)
        return: L + g*sin(t*w)/w**2 + t*(-L/2 - g*sin(2*w)/(2*w**2))

        '''
        x=self.vmain
        ee=MyEq(x-self.e2.ksym,kshow=False)
        svar2=ee.solve_if_and(str(svar),eqv=0,kope=kope,korden=korden,kshow=False,**kwargs)
         
        self.e2.set(svar,svar2,kshow=False)
 
        self.s()     
    
    def solve_if(self,ksym,**kwargs):
        
    
        QQ=self.xcopy()
        p1=QQ.left
        p2=QQ.right
        if len(kwargs)>0:             
            p1=real_subs(p1,**kwargs) 
            p2=real_subs(p2,**kwargs)
        kname=str(ksym)
        ee=MyEq(p1-p2,'ee',kshow=False)
        nee=ee.solve(ksym,kshow=False)
        kres=MyEq(nee,kname)
        return kres
        
        
    def solve_if_and_up(self,*args,**kwargs):
        QQ=self.xcopy()
        sres1=str(QQ.e1.ksym)
        sres2=str(QQ.e2.ksym)
        for key, value in kwargs.items():
                sres1=sres1.replace(key,str(value))
                sres2=sres2.replace(key,str(value))
                 
        QQ.e1.ksym=parse_expr(sres1)
        QQ.e2.ksym=parse_expr(sres2)
        kname=''
        qq=len(args) 
        if qq==1:
            vsym=args[0]
        elif qq==2:
            vsym=args[0]
            kname=args[0]
        else:
            QQ.s()
            return
        if kname=='':
            QQ2=QQ.xcopy()
            QQ2.alone(vsym)
            nval=QQ2.solve(vsym)
            self.set(vsym,nval)
            
        else:
            QQ2=QQ.xcopy()
            ee=MyEq(QQ.solve(vsym),kname)
            self.set(vsym,ee.ksym)

            return ee    
    
    def solveset(self,kmain): # solve kmain in Q and set as vmain
        kres=self.solve(kmain)
        self.e1.ksym=kmain
        self.e2.ksym=kres
        self.vmain=kmain
        self.s()
        
    def solvemain(self): # solve vmain in Q and set final Eq
        kres=self.solve(self.vmain)
        self.e1.ksym=self.vmain
        self.e2.ksym=kres
         
        self.s()

    def msimplify(self):
        return self.solvemain()
        
    def toMyEq(self, kname):
        kname = MyEq(self.e1.ksym - self.e2.ksym, kname=kname)
        return kname





#########################################
#   MyEqEq      uprade  
##########################################



    def upgrade(self,*args,kshow=True):
        self.e1.upgrade(*args,kshow=False)
        self.e2.upgrade(*args,kshow=False)
        self.s(kshow)
    
    def solve_set_if(self,var,korden=0,kope='',kshow=True,**kwargs):
        self.upgrade_if(var=var,korden=korden,kope=kope,kshow=kshow,**kwargs)
          
    def upgrade_if(self,var,korden=0,kope='',kshow=False,**kwargs):
        Q2=self.xcopy()
        Q2.set(kshow=False,**kwargs)
        Q2.simplify(kshow=False)
        if  monodata(Q1.left,'isexp') or monodata(Q2.right,'isexp'):
            Q2.log(kshow=False)
            Q2.lexpand(kshow=False)
        ee2=MyEq(Q2.right-Q2.left,'ee2',kshow=False,kope=kope)
         
        try:
            ee=ee2.ssolve(str(var),kshow=False,kope=kope)
            ee.opemat(kope=kope,kshow=False)
            self.eQ=self.eQ.subs(var,ee)
        except:
            try:
                ee=ee2.solve(var,str(var),korden=korden,kope=kope,kshow=kshow)
            except:    
                ee=ee2.solve(var,str(var),kope=kope)
                return
        if kshow:
            ee.s()
        self.upgrade(ee,kshow=False)
        self.s()
      
    def clear_exp_QQ(self,kshow=True):
    
        p1,p2=clear_exp_QQ(self.left,self.right)
        self.e1.ksym=p1
        self.e2.ksym=p2
        if kshow:
            self.s()
        
    def alone(self,ksym,kshow=False):
       
        p1=self.e1.ksym
        p2=self.e2.ksym
        p3=p1-p2
        ee=MyEq(p3,'ee',kshow=False)
       
        p2=ee.solve(ksym,kshow=False)
      
        p1=ksym
        self.e1.ksym=p1
        self.e2.ksym=p2
        if kshow:
            self.s()
        
         
        
    def clear(self,ksym):
        self.alone(ksym=ksym) 
        
    def nofunc(self): #  convert all x(t) in x  dx in 1 dx2 in 1 in the Eq
        vv=self.var0
        vf=self.ft
        vdf=self.df
        vd2f=self.d2f
        vold=[vd2f,vdf,vf]
        vnew=[1,1,vv]
        
        p1=str(self.left)
        p2=str(self.right)
         
        for i,j in zip(vold,vnew):
            p1=p1.replace(str(i),str(j))
            p2=p2.replace(str(i),str(j))
            

        self.e1.ksym=parse_expr(p1)
        self.e2.ksym=parse_expr(p2)
        self.s()
        return self.var0
         
    
    def prima_func(self): #  convert all x(t) in x  dx in 1 dx2 in 1 in the Eq
        vv=self.f
        vf=self.ft
        vdf=self.df
        vd2f=self.d2f
        vold=[vd2f,vdf,vf]
        vnew=[self.p2f,self.pf,vv]
        
 
        p1=str(self.left)
        p2=str(self.right)
         
        for i,j in zip(vold,vnew):
            p1=p1.replace(str(i),str(j))
            p2=p2.replace(str(i),str(j))
            

        self.e1.ksym=parse_expr(p1)
        self.e2.ksym=parse_expr(p2)
        self.s()    
#########################################
#   MyEqEq     dsolve 
##########################################
    def dsolve(self,*args,op=''):
        expr=self.eQ
        kcond=self.ics
        var2=self.var2
        vmain=self.vmain
        
        # if kcond!='':
            # Q2=MyEqEq(dsolve(expr,ics=kcond),var2=t,vmain=x,kshow=False)
        # else:
            # Q2=MyEqEq(dsolve(expr),var2=t,vmain=x,kshow=False)
            
        if kcond!='':
            kres= dsolve(expr,ics=kcond) 
        else:
            kres=dsolve(expr)

        p1=kres.lhs
        p2=kres.rhs
        kname=str(self.vmain)
        if len(args)==1 and type(args[0])==str:
            kname=args[0]
            nval=symbols(kname)
            nval=MyEq(p2,kname,var2=self.var2)
            return nval
        else:
            nval=symbols(kname)
            nval=MyEq(p2,kname,var2=self.var2)
            return nval
        

    def nodiff(self):
        df=symbol2diff(self.vmain,self.var2)
        d2f=symbol2diff2(self.vmain,self.var2)
        p1=self.left
        p1=p1.subs(d2f,1)
        p1=p1.subs(df,1)
        p2=self.right
        p2=p2.subs(d2f,1)
        p2=p2.subs(df,1)
        self.e1.ksym=p1
        self.e2.ksym=p2
        self.s()
        
    def diff(self,kname=''):
        if kname=='': 
            self.e1.ksym=diff(self.left,self.var2)
            self.e2.ksym=diff(self.right,self.var2)
            self.s()
        else:
            p1=diff(self.left,self.var2)
            p2=diff(self.right,self.var2)
            QQ=MyEqEq(p1,p2,self.vmain,self.var2,vfunc=self.vfunc,kvista=self.vista,kshow=False)
            QQ.s()
            return QQ
            
    def maximize(self,kope=''):
        var2=self.var2
        kres=diff(self.right,self.var2)
        kres=simplify(kres)
        kres=reduFac(kres)
        try:
            e1=MyEq(kres,kshow=False)
            e1.reduFac(kshow=False)
            kres=e1.ksym
        except:
            pass
        kres=solve(kres,self.var2)
        kres2=[]
        for i in kres:
            sres=str(i)
            if not 'I' in sres:
                kres2.append(i)
        kres=kres2
        if len(kres)==1:
            for i in kres:
                MyEq(i,str(var2))
                self.evalue(var2,i)
        else:
            cc=0
            cc2=0
            qq=len(kres)
            fres=self.evalue(var2,kres[0],kshow=False)
            for i in range(qq):
                ik=kres[i]
                val=self.evalue(var2,ik,kshow=False)
                val2=val 
                try:
                    val=float(val)
                except:
                    pass
                if val>fres:
                    fres=val
                    cc2=cc
                    val3=val
                cc+=1
                
             
            MyEq(kres[cc2],str(var2),kope=kope)
            self.evalue(var2,opemat(kres[cc2]),kope)    
        

    def apply_chainRuler(self):
        kres1=self.left
        try: 
            kres1=simplify_chain_rule(kres1)
            self.e1.ksym=kres1
        except:
            pass
            
        kres2=self.right    
        try: 
            kres2=simplify_chain_rule(kres2)
            self.e2.ksym=kres2
        except:
            pass    
        self.s()
        
    



        
#########################################
#   MyEqEq     Integral
##########################################
    def mode_integral(self,kname='',kside='LR'):
        if 'L' in kside:
            p1=self.left
            if type(p1)==Derivative:
                p1= integrate(self.left) 
            else:
                p1=Integral(p1,self.var0)
            self.e1.ksym=p1
             
        if 'R' in kside :
            p2=self.right
            if type(p2)==Derivative:
                p2= integrate(self.left) 
            else:
                p2=Integral(p2,self.var2)
            self.e2.ksym=p2
            
        if kname!='':
            kname2=kname
            kname2=MyEq(self.right,kname,var2=self.var2)
            return kname2
        else:    
            self.s()    
     
    def applyIntegral(self ):

        kname=''
        p1=self.left
        if type(p1)==Derivative:
            p1=Integral(integrate(Q2.left))
        else:
            p1=Integral(p1,self.vmain)
        p2=self.right    
        if type(p2)==Derivative:
            p2=Integral(integrate(Q1.right))
        else:
            p2=Integral(p2,self.var2)

        self.e1.ksym=p1
        self.e2.ksym=p2
        self.s()
        
        
    def separe_dif(self):
        var1=self.vmain
        name1='d'+str(var1)
        
        
        var2=self.var2
        name2='d'+str(var2)
        
        d1=symbols(name1)
        d2=symbols(name2)

        self.e1.ksym=d1/self.right
        self.e1.var2=var1
        self.e2.ksym=d2 
        self.e2.var2=var2
        self.type='eS'
        self.s()
        
    def IntegralEq(self,v1='',v2=''):
         
        x=self.vmain
        t=self.var2
        dx=symbol2diff(x,t)
        d2x=symbol2diff2(x,t)
        p1=self.left
        p1=p1.subs(dx,1)
        p2=self.right
        p2=p2.subs(dx,1)
        v11=self.vmain
        
        if v1!='':
            v11=v1
            
        v22=self.var2
        if v2!='':
            v22=v2
        v11=symbols(str(v11))
        v22=symbols(str(v22))
        p1=func2symbol(p1,self.vmain,self.var2)
        p2=func2symbol(p2,self.vmain,self.var2)
        if Is_Poly(p1):
            mm=0
            for i in fpoly(p1,'list'):
               kk=MyEq(i,'',var2=v11,ktype='I',kshow=False)
               mm=mm+kk.ksym
            p1=mm
        else:    
            kk=MyEq(p1,'p1',var2=v11,ktype='I',kshow=False)
            p1=kk.ksym
        if Is_Poly(p2):
            mm=0
            for i in fpoly(p2,'list'):
               kk=MyEq(i,'',var2=v11,ktype='I',kshow=False)
               mm=mm+kk.ksym
            p2=mm
        else:    
            kk=MyEq(p2,'p2',var2=v22,ktype='I',kshow=False)
            p2=kk.ksym
            
         
        self.e1.ksym=p1 
        self.e2.ksym=p2 
         
        self.s()
        
    def Integral(self,*args):

        var21=symbols(args[0])
        var22=symbols(args[1])
        # clean
        nd1= 'd'+args[0]
        nd2= 'd'+args[1]
        e1s=str(self.left)
        e1s=e1s.replace(nd1,'1')
        e1s=e1s.replace(nd2,'1')
         
        self.e1.ksym=parse_expr(e1s)
        
        e2s=str(self.right)
        e2s=e2s.replace(nd1,'1')
        e2s=e2s.replace(nd2,'1')
         
        self.e2.ksym=parse_expr(e2s)
        self.s(kshow=False)
        
        p1=self.left
        p2=self.right
        x1=''
        x2=0
        if len(args)==4:
            x1=args[2]
            x2=args[3]
            ne1=MyEq(p1,'e1',var2=var21,x1=x1,x2=x2,ktype='I',kshow=False)
            ne2=MyEq(p2,'e2',var2=var22,x1=x1,x2=x2,ktype='I',kshow=False) 
        else:    
        
            ne1=MyEq(p1,'e1',var2=var21,ktype='Ix2x2',kshow=False)
            ne2=MyEq(p2,'e2',var2=var22,ktype='I',kshow=False)
             
        
        self.e1=ne1
        self.e2=ne2
        self.s()
        
    def doitI(self,*args):
        self.e1.doitI(kshow=False,C1=0)
        self.e2.doitI(c1=C1,kshow=False)
        for i in args:
            if i==C1:
                self.e2.ksym=self.right+C1
            if i==C2:
                self.e2.ksym=self.right+C2    
        self.s()

#########################################        
#  algebrate transformation
#########################################

    def subsnumber(self,val1,val2):
     
         
        p1=self.left
        p1=subsnumber(p1,val1,val2)
        self.e1.ksym=p1
    
        p2=self.right
        p2=subsnumber(p2,val1,val2)
        self.e2.ksym=p2
        self.s()

    def sortdegree(self,var='LR',kope=''):
        '''
        input (x**(a*b))   ---->   return(x**a)**b
        input (x**(a*b),b)   ---->   return(x**b)**a
        '''
        kside='LR'
        if type(var)==str:
            kside=var
         
            
        if  type(var)!=str:
            var2=var
        else:
            var2=x
        if 'L' in kside:
            p1=self.left
            p1=sortdegree(p1,var2=var2)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=sortdegree(p2,var2=var2)
            self.e2.ksym=p2
        self.s()   
    
    def powexpand(self,args):
        '''
        input (x**(a*b))   ---->   return(x**a)**b
        input (x**(a*b),b)   ---->   return(x**b)**a
        '''
        op=''
        kside='LR'
        for i in args:
            if type(i)==str:
                kside=i
            if type(i)==Symbol:
                op=i
        if 'L' in kside:
            p1=self.left
            p1=powexpand(p1,op=op)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=powexpand(p2,op=op)
            self.e2.ksym=p2
        self.s()    
    
    def opematexp(self,kside='LR',kope=''):
        '''
        apply opemat only in exponent monomie
        '''
        if 'L' in kside:
            p1=self.left
            p1=opematexp(p1,kope=kope)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=opematexp(p2,kope=kope)
            self.e2.ksym=p2
        self.s()
        
    def mulexp(self,kside='LR'):
        '''
        input (x**a)**b
        return(x**(a*b))
        '''
        if 'L' in kside:
            p1=self.left
            p1=mulexp(p1)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=mulexp(p2)
            self.e2.ksym=p2
        self.s()
        
     
    def packexp(self,kside='LR'):
        '''
        input (x**a)**b
        return(x**(a*b))
        '''
        if 'L' in kside:
            p1=self.left
            p1=packexp(p1)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=packexp(p2)
            self.e2.ksym=p2
        self.s()
    
    
    
    def base2frac(self,kside='LR'):
         
        if 'L' in kside:
            p1=self.left
            p1=base2frac(p1)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=base2frac(p2)
            self.e2.ksym=p2
        self.s()


    
    def expandexp(self,kside='LR'):
        '''
        input (x**a)**b
        return(x**(a*b))
        '''
        if 'L' in kside:
            p1=self.left
            p1=expandexp(p1)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=expandexp(p2)
            self.e2.ksym=p2
        self.s()
        
    def factorexp(self,kside='LR'):
        '''
        input (x**a)**b
        return(x**(a*b))
        '''
        if 'L' in kside:
            p1=self.left
            p1=factorexp(p1)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=factorexp(p2)
            self.e2.ksym=p2
        self.s()
    
    
    def factorize(self,kval,kside='LR'):
        if 'L' in kside:
            p1=self.left
            p1=factorize(p1,kval)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=factorize(p2,kval)
            self.e2.ksym=p2
        self.s()
        
    def killsqrtpow(self,kside='LR'):
        if 'L' in kside:
            p1=self.left
            p1=killsqrtpow(p1)
            self.e1.ksym=p1
        if 'R' in kside:
            p2=self.right
            p2=killsqrtpow(p2)
            self.e2.ksym=p2
        self.s() 
    
    def killRootPow(self):
        p1=self.left
        p11=kill_RootPow(p1)
        if str(p1)!=str(p11):
            self.e1.ksym=p11
        p2=self.right
        p22=kill_RootPow(p2)
        if str(p2)!=str(p22):
            self.e2.ksym=p22
        self.s()    

        
    def par_frac(self,var=''):
        self.e1.par_frac(var=var,kshow=False)
        self.e2.par_frac(var=var,kshow=False)
        self.s()
        
        
    def simple_linear(self):
        if self.e1.Is_Mono() and self.e2.Is_Mono():
            p1=self.e1.get_nume()
            p2=self.e1.get_deno()
            p3=self.e2.get_nume()
            p4=self.e2.get_deno()
            self.e1.ksym=p1/p3
            self.e2.ksym=p2*p4
            mm=fpoly(self.left,'list')
            if -1 in mm:
                self.Mul(-1)
            else:
                self.s()    

    def crossMul(self):
        p1=numer(self.left)
        p2=denom(self.left)
        p3=numer(self.right)
        p4=denom(self.right)
        
        self.e1.ksym=p1*p4
        self.e2.ksym=p3*p2
        
        self.s()
        

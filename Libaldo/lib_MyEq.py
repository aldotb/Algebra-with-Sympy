

from sympy import *
 
from IPython.display import Math  # ,display
from matplotlib.pyplot import ylabel, plot, show, xlabel, title
from libaldo_math2 import *
from libaldo_algorith import *
from libaldo_show import *
import copy



# from lib_Func import *
import dill  # pip install dill --user

filename = 'workspace.pkl'


def savework():
    dill.dump_session(filename='workspace.pkl')


def loadwork():
    dill.load_session(filename='workspace.pkl')


# and to load the session again:


C1, C2, C3, C4, t, x, y, z = symbols('C1 C2 C3 C4 t x y z')
dataQ = []
e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12 = symbols('e1 e2 e3 e4 e5  e6 e7 e8 e9 e10 e11 e12')


class MyEq:
    def __init__(self, ksym, kname='', var2='', var1=x, varx=x, vary=y, varz=z, kp=False, kope='', kshow=True, ktype='P',
                 xx='', dtype=1, depen=False, x1='', x2=0, y1='', y2=0, z1='', z2=0, Pobj='', ee='', ssym='', init=True,
                 kfull=True, andsolve='',varf=[],diffEq='',eQshow='',vfunc=[],vmain=x):

        #if ksym!=0:
            #datap.inputData(kname,ksym,ktype,True)
        self.ksym=ksym
        self.t=ktype
        self.pname=kname
        self.type=ktype
        self.kinte=''
        self.ee=ee
        self.depen=depen
        self.primi=''
        self.x1=x1
        self.x2=x2
        self.y1=y1
        self.y2=y2
        self.z1=z1
        self.z2=z2
        self.varx=varx
        self.vary=vary
        self.varz=varz
        self.eeq=ksym
        self.vmain=vmain
        if kname!='' and len(str(kname))==1:
            self.vmain=symbols(kname)
            
        if var2=='':
            var2=t
        self.var2=var2
        self.var1=var1
        self.name = alphaname(kname)
        self.sname=''
        self.kinte=''
        self.xx=xx
        self.EqDiff=''
        self.Adiff=''
        self.Bdiff=''
        self.EqInte=''
        self.iniG=0
        self.odeQ=''
        self.ode=''
        self.oldprimi=''
        self.Pobj=Pobj
        self.ssym=ssym
        self.father=''
        self.backup=[] 
        self.kshow=kshow
        self.init=init 
        self.histo=ksym
        self.primitiva=func2primi(ksym,var2,varf)
        try:
            self.pdiff=diff(func2primi(ksym,var2,varf),var2)
        except:
            self.pdiff=0
        try:
            self.flatdiff=Diff2diff(diff(func2primi(ksym,var2,varf),var2),var2,varf)
        except:
            self.flatdiff=0
        
        self.varf=varf    
        self.primi_eq=''
        self.diff_eq=''
        self.primi_diff_eq='' 
        self.origen=1
        
        
        if type(ksym)==MyEq:
            seq=str(ksym())
            self.ksym = unisymbols(opemat(ksym(), kope=kope))
            self.v = unisymbols(opemat(ksym(), kope=kope))
            self.eeq=ksym
        #elif type(ksym)==MyEqMat:
        #    self.ksym = unisymbols(opemat(ksym(), kope=kope))
        #elif type(ksym)==Equality:
        #    mm=ksym.e2.ksym
        #    self.ksym=mm=kQ.e2.ksym
        #     
        #    try:
        #        self.var2=get_varfunc(ksym.lhs)
        #    except:
        #        done=False
        else:    
            seq=str(ksym)
            #self.ksym = unisymbols( ksym )
            self.v = unisymbols( ksym )

        if ktype=='F' or ktype=='Ff':
            if var2!='':
                s1=kname
                if type(var2)==list:
                    s2=alphaname(var2[0])+','+alphaname(var2[1])
                else:
                    s2=alphaname(var2)
                kname=s1+'_{('+s2+')}'
                self.name = kname
         
        
        if  'I' in ktype:
            self.primi=self.ksym
            if self.x1=='': 
                self.ksym=Integral(self.ksym, self.var2)
            else:
                self.ksym=Integral(self.ksym, (self.var2,self.x1,self.x2))
            self.histo = unisymbols(opemat(self.ksym, kope=kope))
        if '2I' in ktype or '3I' in ktype:
             
            if self.y1=='': 
                self.ksym=Integral(self.ksym, self.vary)
            else:
                self.ksym=Integral(self.ksym, (self.vary,self.y1,self.y2))
            self.histo = unisymbols(opemat(self.ksym, kope=kope))

        if '3I' in ktype:
             
            if self.z1=='': 
                self.ksym=Integral(self.ksym, self.varz)
            else:
                self.ksym=Integral(self.ksym, (self.varz,self.z1,self.z2))
            self.histo = unisymbols(opemat(self.ksym, kope=kope))            
        
        if ktype=='Diff':
            self.Adiff=Derivative(self.ksym, varx)
            self.sname='d'+kname
            self.name=kname
            self.var2=varx
                
         
            
        if ktype=='diff' or ktype=='diff2':
            f=Function(kname)(var2)
            kres=self.ksym
            kres=kres.subs(var1,f)
            var2=self.var2
            var1=self.var1
            self.Bdiff=kres
            if ktype=='diff':
                self.Adiff=Derivative(f, var2)
            else:
                self.Adiff=Derivative(f, var2,var2)
            self.EqDiff=Eq(self.Adiff,self.Bdiff)
                

        if ktype=='fdiff':
            self.name='d'+kname+'_'+str(var2)
            nksym=get_diff_name(str(ksym),str(var2))
             
            
            self.ksym=nksym
            self.var1=ksym
            
        if self.Pobj!='':
            self.type='Ph'
            
        if kshow:
            if ktype=='diff' or ktype=='diff2':
                display(Math(latex(self.EqDiff)))
            elif ktype=='Diff':
                kres=self.ksym
                ps1=self.sname+'='
                ps2='d_'+str(self.varx)
                 
                display(Math(ps1+latex(kres)+ps2))
                
                
            else:
                kres=self.ksym
            
                if self.name == '':
                    display(Math(latex(kres)))
                else:
                    sR = self.name + ' ='
                    display(Math(sR + latex(kres)))
        
        
        if andsolve!='':
            kval=andsolve
            kname=str(andsolve)
             
            kres=csolve(self.ksym,kval)
            kres=opemat(kres,kope)
            andsolve=MyEq(kres,kname)
            return (andsolve) 
        if dtype==1:
            if self not in dataQ and self.name!='': 
                dataQ.append(self)
        #self.ksym=sympify(str(self.ksym), locals={str(var2): var2})     
  
     

        self.primi_eq=ksym 
        self.diff_eq=''
        self.primi_diff_eq=''
        self.origen=1
        self.eQshow=eQshow
        self.vfunc=vfunc
        self.vdf=[]
        self.vd2f=[]
        self.vdp=[]
        self.vd2p=[]
        
     
    def __call__(self, kname='', kope='', kupdate=False, *args, **kwargs):
        if self.init == False:
            self.ksym = kname
            self.init = True
            self.s()
            return

        if type(kname) != str and len(kwargs) == 0 and len(args) == 0:
            var2 = self.var2
            ksym = self.ksym
            return opemat(ksym.subs(var2, kname), kope=kope)

        if len(args) == 1 and len(kwargs) == 0 and self.iniG == 0 and self.t == 'Q':
            self.ksym = opemat(args[0], kope=kope)
            self.iniG = 1
            self.s()
            return
        if len(args) == 1 and len(kwargs) == 0 and self.t != 'Q':
            ksym = self.ksym
            var2 = self.var2
            nvar = args[0]
            ee = MyEq(ksym, kname=self.name, kshow=False, kope=kope)
            ee.set(var2, nvar)

            return
        if len(args) == 1 and len(kwargs) == 0 and self.iniG == 1 and args[0] == '':
            if self.t == 'Q':
                self.ksym = 0
                self.iniG = 0
                # self.s()
                return

        if len(args) > 0 and len(kwargs) == 0:
            kvar = self.var2
            if kvar != '':
                kres = self.ksym
                kname = self.name
                ee = MyEq(kres, kname, kshow=False)
                if type(kvar) != list:
                    ee.set(kvar, args[0])
                else:
                    for i in range(len(kvar)):
                        ee.set(kvar[i], args[i], kshow=False, kope=kope)
                    ee.s()

        if len(kwargs) > 0:
            kres=self.ksym
            kres=real_subs(kres,**kwargs) 

            ee=self.xcopy('ee',kshow=False)
            ee.ksym=kres
            if kupdate:
                self.update(ee.ksym)
                self.s()
                return
            if kname != '':
                ee.name = kname
                ee.s()
                return ee
            if len(args) == 1:
                if type(args[0]) == str:
                    ee.name = kname
                    ee.s()
                    return ee
            return opemat(ee.ksym, kope=kope)

        return opemat(self.ksym, kope=kope)

    def __repr__(self):
        kres = str(self.ksym)

        return kres

    def _latex(self, obj):
        return latex(self.ksym)

    def __str__(self):
        return self.__repr__()

    ###########################################
    #               variables                    #
    ###########################################

    def norm_variable(self):
        clist = [C1, C2]
        m1 = self.free()
        s1 = [str(x) for x in m1]
        for kvar in clist:
            svar = str(kvar)
            for i in range(len(s1)):
                if s1[i] == svar:
                    kvar = m1[i]


    ###########################################
    #               Update                    #
    ###########################################

    def __add__(self, other):
        """ Returns the vector addition of self and other """
        if type(other) == MyEq :
            kres = self.ksym + other.ksym
        else:
            kres = self.ksym + other
        if type(kres) == MyEq:
            kres.s()
        return kres

    def __radd__(self, other):
        if type(other) == MyEq :
            kres = self.ksym + other.ksym
        else:
            kres = self.ksym + other
        return kres

    def __sub__(self, other):
        if type(other) == MyEq :
            kres = self.ksym - other.ksym
        else:
            kres = self.ksym - other
        return kres

    def __rsub__(self, other):
        if type(other) == MyEq :
            kres = self.ksym - other.ksym
        else:
            kres = self.ksym - other
        return kres

    def __mul__(self, other):
        """ Returns the vector addition of self and other """
        if type(other) == MyEq :
            kres = self.ksym * other.ksym
        else:
            kres = self.ksym * other
        return kres

    def __rmul__(self, other):
        """ Returns the vector addition of self and other """
        if type(other) == MyEq :
            kres = self.ksym * other.ksym
        else:
            kres = self.ksym * other
        return kres
        
    
    def __truediv__(self, other):
        return  self.ksym / (1*other)

    def __rtruediv__(self, other):
        return self.ksym / (1*other)  
    
        
    def Di(self):
        if self.varf!=[]:
            return self.primi_eq
            
    ###########################################
    #       show
    ########################################   
    
    def s2(self,op=''):

        kres = self.ksym
        if op=='2':
            kres=self.eQshow
        if self.type == 'Diff':

            ps1 = self.sname + '='
            ps2 = 'd_' + str(self.varx)

            display(Math(ps1 + latex(kres) + ps2))



        elif self.type == 'diff' or self.type == 'diff2':
            display(Math(latex(self.EqDiff)))
        else:
            if self.name == '':
                display(Math(latex(kres)))
            else:
                sR = self.name + ' ='
                display(Math(sR + latex(kres)))

    def s(self,op='', **kwargs):
        if self.father=='':
            qk = len(kwargs)
            if qk > 0:
                ee = MyEq(self.v, kshow=False)
                for key, value in kwargs.items():
                    ee.set(parse_expr(key), value, kshow=False, ktype=self.type)
                ee.s2(op=op)

            else:
                self.s2(op=op)
        else:
            mQ.s()
    
    def showd(self):
        kres=self.ksym
        for i in self.vfunc:
            kres=show_modefunc(kres,i,ktype=2)
        ee=MyEq(kres,kname=self.name)
        
    def showp(self):
        kres=self.ksym
        for i in self.vfunc:
            kres=show_modefunc(kres,i,ktype=3)
        ee=MyEq(kres,kname=self.name)
      

    
    ###########################################
    #       set & update
    ########################################    
    def MsetValue(self, vknom, vkval, kope=''):
        for knom, kval in zip(vknom, vkval):
            knom = unisymbols(knom)
            kval = unisymbols(kval)
            kres = self.ksym
            kres = kres.subs(knom, kval)
            kres = opemat(kres, kope=kope)

            kres = self.simplify()
            self.update(kres)
        self.s()

    def multiSet(self, vknom, vkval, kope=''):
        for knom, kval in zip(vknom, vkval):
            knom = unisymbols(knom)
            kval = unisymbols(kval)
            kres = self.ksym
            kres = kres.subs(knom, kval)
            kres = opemat(kres, kope=kope)

            kres = self.simplify()
            self.update(kres)
        self.s()




    def sset(self,kname='',vmain='',var2='',kshow=True,**kwargs):
        if vmain=='':
            vmain=self.vmain
        if var2=='':
            var2=self.var2
        if len(kwargs)>0:
             
            kres=self.ksym
            for key, value in kwargs.items():
                nsym=parse_expr(key)
                if nsym not in fpoly(kres,'free'):
                    kname=str(nsym)
                    var2=self.var2
                    kres=kres.subs(Function(kname)(var2),value)
                    
                else:
                    kres=kres.subs(parse_expr(key),value)
                                   
        if kname!='':
             
            nQQ=MyEq(kres,kname=kname,vmain=vmain,var2=var,kshow=kshow)
            return nQQ
        else:    
            self.ksym=kres    
        if kshow:    
            self.s()    
    
    def evalif(self, kname='', **kwargs):
        ee = self.xcopy(kname=kname, kshow=False)

        qk = len(kwargs)
        if qk > 0:
            for key, value in kwargs.items():
                if type(value) == MyEq:
                    value = value.ksym
                ee.set(parse_expr(key), value, kshow=False)
                try:
                    ee.primi = ee.primi.subs(parse_expr(key), value)
                except:
                    done = False

            if kname != '':
                ee.s()
                return ee
            else:
                ee.s()

        else:
            self.s()

    def set(self, knom='', kval='', kshow=True, kope='', kret=False, andsolve=[], Bag='', **kwargs):
        if len(kwargs) > 0 and knom == '' and kval == '':
            for key, value in kwargs.items():
                self.set(parse_expr(key), value, kshow=False, kope=kope)
            self.s()
            return

        if type(kval) == Symbol:
            kval = kval.name

        if self.type == 'Ph':
            try:
                self.Pobj.store_val(knom, kval)
                mm = self.Pobj.F
                qq = len(mm)
                for i in range(qq):
                    pkres = mm[i][0]
                    try:
                        pkres = pkres.subs(knom, kval)
                    except:
                        done = False
                    mm[i][0] = pkres
                self.Pobj.F = mm
            except:
                done = False

        try:
            kres = self.primi

            self.primi = kres.subs(knom, kval)
        except:
            done = False
        qk = len(kwargs)
        if qk > 0:
            for key, value in kwargs.items():
                self.set(parse_expr(key), value, kshow=False, kope=kope)
                try:
                    self.primi = self.primi.subs(parse_expr(key), value)
                except:
                    done = False
            if Bag != '':
                self.upBag(Bag, kshow=False, kope=kope)
            if kshow:
                self.s()

        else:
            if andsolve != []:

                if kval == '':
                    self.setValue(self.v, knom, kshow=False, kope=kope)

                else:
                    self.setValue(knom=knom, kval=kval, kshow=False, kope=kope, kret=kret)

                    ee = MyEq(opemat(self.solve(andsolve[0]), kope=kope), str(andsolve[1]))

                    return ee
                if Bag != '':
                    ee.upBag(Bag, kope=kope)
                return ee

            else:
                if kval == '':
                    self.setValue(self.v, knom, kope=kope, kshow=False)
                    print(3)
                    if Bag != '':
                        self.upBag(Bag, kshow=False)
                else:
                    mm = self.setValue(knom=knom, kval=kval, kshow=False, kope=kope, kret=kret)

                    if Bag != '':
                        self.upBag(Bag, kshow=False)
                if kshow:
                    if Bag == '':
                        self.s()

    def setValue(self, knom, kval, kshow=True, kope='', kret=False):
        if type(knom) != list:
            knom = [knom]
            kval = [kval]
        for i, j in zip(knom, kval):
            i = unisymbols(i)
            j = unisymbols(j)
            kres = self.ksym
            try:
                kres = kres.subs(i, j)
            except:
                pass
            kres = opemat(kres, kope)

            self.update(kres)

        if kret:
            return self.ksym

        if kshow:
            self.s()
 
    def setdiff(self, kvar, kval, kshow=True, kope='', kret=False):
        try :
            kval=alphaname(kval)
        except:
            done=False
        self.set(knom=kvar.diff(), kval=kval, kshow=kshow, kope=kope, kret=kret)

    # def evalif(self, *args):
        # vsym = []
        # vval = []
        # done = True
        # for i in args:
            # if done:
                # vsym.append(i)
                # done = False
            # else:
                # vval.append(i)
                # done = True
        # kres = self.ksym
        # for i, j in zip(vsym, vval):
            # kres = kres.subs(i, j)
        # ee = MyEq(kres, self.name)

    def eval(self, **kwargs):
        kname = self.name
        ksym = self.ksym

        key, value = kunpakDic(kwargs)
        for i, j in zip(key, value):
            ksym = ksym.subs(i, j)
        MyEq(ksym, self.name)
 
    def set_solve(self, kset, kvset, kvars, knames):
        ee = MyEq(self.ksym, kshow=False)
        ee.set(kset, kvset, kshow=False)
        kres = ee.solve(kvars)
        ee2 = MyEq(kres, knames)
        return ee2

    def set_solveR(self, kset, kvset, kvars, knames):
        ee = MyEq(self.ksym, kshow=False)
        ee.set(kset, kvset, kshow=False)
        kres = ee.solveR(kvars)
        ee2 = MyEq(kres, knames)
        return ee2



    def v(self):
        return self.v

    def update(self, kres):
        if self.type == 'Z':
            if self.name == 'e2':
                self.upgrade(e1)

        self.histo = self.ksym
        self.ksym = kres
        self.v = kres
        self.oldprimi = self.primi
        self.Adiff = Derivative(self.ksym, self.varx)

    def upgrade(self, *args, kope='', kshow=True, **kwargs):
 
        if len(args) == 1:
            if type(args[0]) == list:
                args = args[0]

        for i in args:
            kname = i.name
            self.set(kname, i, kshow=False, kope=kope)
            if 'I' in self.type:
                ee = MyEq(self.primi, kshow=False)
                ee.set(kname, i, kshow=False, kope=kope)
                self.primi = ee.ksym
            try:
                self.primi = self.primi.subs(kname, i)
            except:
                done = False
        if len(kwargs) > 0:
            self.set(**kwargs, kshow=False)
        if kshow:
            self.s()

    def kret(self):
        return self.ksym
        sR = self.name + ' ='
    
    ###########################################
    #            edit copy
    ###########################################

    def xcopy(self, kname, kshow=True):
        ee = copy.deepcopy(self)
        if kname != '':
            ee.name = kname
        if kshow:
            ee.s()
        return ee
    def sExp(self, kval):
        if self.name == '':
            display(Math(latex(self.ksym)))
        else:
            sR = self.name + ' ='
            display(Math(sR + latex(kval)))

    def undo(self,kshow=True):
        self.ksym = self.histo
        self.v = self.histo
        if kshow:
            self.s()

    def cut_denom(self, kope=''):
            kres = numer(self.ksym)
            kres = opemat(kres, kope)
            self.update(kres)
            return self.ksym


    ###########################################
    #            Math operation
    ###########################################

    def reduF(self, kshow=True, kupdate=True):

        kres = self.ksym
        kres = factor(kres)
        kres2 = 1
        if Is_Mono(kres):
            kres = numer(kres)

            mm = fpoly(kres, 'list')
            for i in mm:
                if Is_Poly(i):
                    kres2 = kres2 * i

        if kres2 != 0 and kres2 != 1:
            self.ksym = kres2
        self.s()

    def Add(self, kval, kname='',kope='', kupdate=True, kshow=True):
        kres = self.ksym
        if Is_Add(kres) and kope != '':
            mm = 0
            for i in fpoly(kres, 'list'):
                nval = i + kval
                nval = opemat(nval, kope=kope)
                mm += nval
            kres = mm

        else:
            kres = kres + kval
        kres = opemat(kres, kope=kope)
        if kname!='':
            ee=MyEq(kres,kname)
            return ee

        if kupdate:
            self.update(kres)
            if kshow:
                return self.ksym
        else:
            return self.ksym

    def Mul(self, kval, kname='',kope='', kupdate=True, kshow=True):

        kres = self.ksym
        if Is_Add(kres):
            mm = 0
            for i in fpoly(kres, 'list'):
                nval = i * kval
                nval = opemat(nval, kope=kope)
                mm += nval
            kres = mm

        else:
            kres = kres * kval
        if '(' not in kope:
            kres = opemat(kres, kope=kope)
            if kname!='':
                ee=MyEq(kres,kname)
                return ee
            if kupdate:
                self.update(kres)
                if kshow:
                    self.s()
            else:
                return self.ksym
        else:
            opemat(kres, kope=kope)

    def Div(self, kval, kname='',kope='', kupdate=True, kshow=True):
        kres = self.ksym
        if Is_Add(kres):
            mm = 0
            for i in fpoly(kres, 'list'):
                nval = i / kval
                nval = opemat(nval, kope=kope)
                mm += nval
            kres = mm

        else:
            kres = kres / kval

        kres = opemat(kres, kope=kope)
        if kname!='':
            ee=MyEq(kres,kname)
            return ee

        if kupdate:
            self.update(kres)
            if kshow:
                self.s()
        else:
            return self.ksym

    def Pow(self, kval, kname='',kope='', ksec=False, kupdate=True, kshow=True):
        kres = kpow(self.ksym, kval)
        if ksec:
            kres = opematPolySec(kres)
         
        if kope!='':
            kres = opemat(kres, kope=kope)

        if kname!='':
            ee=MyEq(kres,kname)
            return ee

        if kupdate:
            self.update(kres)
            if kshow:
                self.s()
        else:
            return self.ksym
    def sqrs(self,*args,kope='',kupdate=True, kshow=True): # return pow(k1,S('k2'))
        k1=self.ksym
        k2=''
        k3=''
        kname=''
        qq=len(args)
        if qq==0:
            kres=sqrs(k1=k1,k2=k2,k3=k3)
             
        elif qq==1:
            if type(args[0])==str:
                kname=args[0]
            else:
                k2=args[0]
        elif qq==2:
            if type(args[0])==str:
                kname=args[0]
                k2=args[1]
                
            else:
                k2=args[0]
                k3=args[1]
        else:
            kname=args[0]
            k2=args[1]
            k3=args[2]
            
        kres=sqrs(k1=k1,k2=k2,k3=k3)
        if kope!='':
            kres = opemat(kres, kope=kope)
        
        if kname!='':
            ee=MyEq(kres,kname)
            return ee
            
        if kupdate:
            self.update(kres)
            if kshow:
                self.s()
        else:
            return kres    
        
    
    def Rpow(self, kval, kname='', kope='', ksec=False, kupdate=True, kshow=True):
        kres = rpow(self.ksym, kval)
        if ksec:
            kres = opematPolySec(kres, kope)
        else:
            kres = opemat(kres, kope)

        if kname!='':
            ee=MyEq(kres,kname)
            return ee

        if kupdate:
            self.update(kres)
            if kshow:
                self.s()
        else:
            return self.ksym

    def cut_fac(self, kval):
        kres = cut_fac(self.ksym, kval)
        self.update(kres)
        return kres

    def signo(self):
        kres = self.ksym
        return signo(kres)
        
    def par_frac(self,var='',kshow=True):
        try:
            kres=par_frac(self.ksym,var=var)
            self.ksym=kres    
        except:
            pass
        if kshow:
            self.s()
                
    def float(self, kupdate=True, kshow=True):
        kres = self.ksym
        try:
            kres = float(kres)

        except:
            self.s()
            return
        if kupdate:
            self.update(kres)
            if kshow:
                self.s()
        else:
            return self.ksym
            
    ###########################################
    #            get info
    ###########################################
    def get_dview(self):
        kres=self.ksym
        for i in self.vfunc:
            kres=show_modefunc(kres,i,ktype=2)
        return kres
        
    def get_pview(self):
        kres=self.ksym
        for i in self.vfunc:
            kres=show_modefunc(kres,i,ktype=3)
        return kres 
        
        
    def list(self, kopt=''):
        kres = self.ksym
        kres = fpoly(kres, 'list')
        if kopt != '':
            return kres[kopt]
        else:
            return kres
    
    def slist(self, kopt=''):  # symbols list
        kres = self.ksym
        kres = fpoly(kres, 'list')
        vres=[]
        for i in kres:
            vres.append(num_ksym2ksym(i))
        if kopt != '':
            return vres[kopt]
        else:
            return vres


    def free(self):
        kres = self.ksym
        kres = fpoly(kres, 'free')
        return kres

    def args(self, opt1=''):
        kres = self.ksym
        if opt1 == '':
            kres = kres.args
            return kres
        else:
            try:
                kres = fpoly(kres, 'list')
                return kres[opt1]
            except:
                return self.ksym

    def get_primitive(self,kname='',kshow=True):
        kres=self.ksym
        if self.varf!=[]:

            ee=self.xcopy(kname,kshow=False)
            nf=[Function(x)(ee.var2) for x in ee.varf]
            for i,j in zip(ee.varf,nf):
                ee.set(i,j,kshow=False)
            
            self.primi_eq=ee.ksym
            if kname!='':
                ee.s()
                return ee
            else:
                return ee.ksym
        else:
            self.s()

    

    
    def forceSetValue(self, knom, kval):

        sknom = str(knom)
        skval = str(kval)
        skvalue = str(self.ksym)
        skvalue.replace(sknom, skval)
        kfac = 1
        if skavlue[0] == '-':
            skvalue = skvalue[1:-1]
            kfac = -1
        kres = parse_expr(skvalue)
        kres = kres * kfac
        self.update(kres)
        self.s()



    def kdiff(self, kval, kope='', kupdate=False):
        kres = kdiff(self.ksym, kval)
        kres = opemat(kres, kope)
        if kupdate:
            self.update(kres)
            self.s()
        else:
            return kres
            
    
        
    def integral(self, kname='', x1='', x2='', kope='', kupdate=False, ktype='P', var2='', C1=C1):
        if self.type == 'Diff':
            self.type = 'I'
            if x2 != '':
                self.x2 = x2
            if x1 != '':
                self.x1 = x1
                if x1 == 0 and x2 == '':
                    self.x2 = self.varx

            self.primi = self.ksym
            self.var2 = self.varx
            if self.x1 == '':
                self.ksym = Integral(self.ksym, self.var2)
            else:
                self.ksym = Integral(self.ksym, (self.var2, self.x1, self.x2))
            self.histo = unisymbols(opemat(self.ksym, kope=kope))
            self.s()
        else:
            keq = self.ksym
            kres = keq
            if var2 != '':
                kvar = var2
            elif self.var2 != '':
                kvar = self.var2
            else:
                kvar = t

            if x1 == '':
                kres = integrate(keq, kvar)
                # self.update(kres)
            else:
                kres = integrate(keq, (kvar, x1, x2))
            kres = opemat(kres, kope)
            if kname != '':
                return MyEq(kres + C1, kname, var2=kvar, ktype=ktype)
            if kupdate:
                self.update(kres)
                self.s()
            else:
                return kres

    def modo_integral(self,*args):
        Iee=self.xcopy('iee',kshow=False)
        kres=Iee.ksym
        kname=''
        for i in args:
            if type(i)==str:
                kname=i
            else:
                ival=i[0]
                delI=symbols(diffname(ival))
                Iee.set(delI,1,kshow=False)
        kres=Iee.ksym
        for i in args:
            if type(i)!=str:
                if len(i)==1:
                    kres=Integral(kres,i[0])
                elif len(i)==2:
                    kres=Integral(kres,(i[0],0,i[1]))
                else:
                    kres=Integral(kres,(i[0],i[1],i[2]))
        if kname=='':
            self.ksym=kres
            self.s()
        else:
            newee=MyEq(kres,kname=kname)
            return newee

    def tintegral_def(self, alpha1, a1, a2, kupdate=False):
        keq = self.ksym
        kres = tintegral_def(keq, alpha1, a1, a2)
        if kupdate:
            self.update(kres)
            self.s()
        else:
            return kres

    def only_nume(self, kope=''):
        kres = self.get_nume()
        kres = opemat(kres, kope)
        self.update(kres)
        self.s()

    def only_deno(self, kope=''):
        kres = self.get_deno()
        kres = opemat(kres, kope)
        self.update(kres)
        rself.s()

    def get_MonoExp(self):
        kres2 = self.ksym
        kres = kres2.fpoly('get', 1)
        return kres

    def evalue(self, *args, kope='', kshow=False, **kwargs):
        qk = len(kwargs)
        if qk > 0:
            for key, value in kwargs.items():
                self.set(parse_expr(key), value, kshow=False)
                try:
                    self.primi = self.primi.subs(parse_expr(key), value)
                except:
                    done = False
        if len(args) == 1:
            var2 = self.var2

            kres = self(var2=args[0])
        return kres

    def evalueArray(self, knom, vkval):

        kres = [self.evalue(knom, xx, kshow=False) for xx in vkval]
        return kres

    ###########################################
    #            Functions
    ###########################################
    
    @property
    def Type(self):
        return type(self.ksym)
        
    def slope(self, xx, kval='', kope=''):
        xExp = self.ksym

        kres = self.kdiff(xx)
        if kval != '':
            kres = kres.subs(xx, kval)
        kres = opemat(kres, kope)
        return kres

    def slopeO(self, xx, kval='', kope=''):
        xExp = self.ksym

        kres = self.kdiff(xx)
        kres = -1 / kres
        if kval != '':
            kres = kres.subs(xx, kval)
        kres = opemat(kres, kope)
        return kres

    def get_aTan(self, kname='', **kwargs):
        if kname != '':
            return MyEq(self.angTan(**kwargs), kname=kname)

        return self.angTan(**kwargs)

    def angTan(self, kname='', **kwargs):
        kres = self.ksym
        var2 = self.var2
        kres = kres.diff(var2)
        if len(kwargs) > 0:
            ee = MyEq(kres, kshow=False)
            kres = ee(**kwargs)
            if len(fpoly(kres, 'list')) > 1:
                mm = fpoly(kres, 'list')
                qq = len(mm)
                kres = mm[qq - 1]
        if kname != '':
            return MyEq(atan(kres), kname)
        else:
            return atan(kres)

    def get_aOrto(self, kname='', **kwargs):
        if kname != '':
            return MyEq(self.angOrto(**kwargs), kname=kname)
        return self.angOrto(**kwargs)

    def angOrto(self, kname='', **kwargs):
        kres = self.angTan(**kwargs)
        if kname != '':
            return MyEq(kres, kname)
        else:
            return (kres + pi / 2)

    def ang_vecTan(self, **kwargs):
        kres = ee.diff(**kwargs)
        kres = atan(kres)
        return kres

    def ang_vecOrto(ee, **kwargs):
        kres = ee.diff(**kwargs)
        kres = atan(kres)
        return kres + pi / 2

    def Is_Poly(self):
        return Is_Poly(self.ksym)

    def Is_Mono(self):
        return (Is_Mono(self.ksym))

    def Is_Add(self):
        kres = self.ksym
        if type(kres) == Add:
            return True
        else:
            return False

    def Is_Mul(self):
        kres = self.ksym
        if type(kres) == Mul:
            return True
        else:
            return False

    def Is_Pow(self):
        kres = self.ksym
        if type(kres) == Pow:
            return True
        else:
            return False

    ###########################################
    #            Transformation
    ###########################################
     
    
    def subsnumber(self,val1,val2,kope='',kshow=True):
        kres=self.ksym
        kres=subsnumber(kres,val1,val2)
        if kope!='':
                kres=opemat(kres,kope=kope)
        self.ksym=kres
        if kshow:
                self.s()
    
    def all_type(self):
        self.s()
        sE(['Monomie= ', self.Is_Mono(), '  ', 'Polynomie = ', self.Is_Poly()])
        sE(['Is Add= ', self.Is_Add(), '  ', 'Is Mul= ', self.Is_Mul()])

    def expand(self, op='',kupdate=True, kshow=True, kope=''):
        kres = self.ksym
        if op=='numer':
            p1=numer(kres)
            p2=denom(kres)
            p1=expand(p1)
            kres=(p1/p2)
            
        elif op=='denom':
            p1=numer(kres)
            p2=denom(kres)
            p2=expand(p2)
            kres=(p1/p2)
        else:
            kres = expand(kres)
        kres = opemat(kres, kope=kope)
        if kupdate:
            self.update(kres)
        if kshow:
            self.s()
    def simplifybase(self,kope='',kshow=True):
        kres=self.ksym
        kres=simplifybase(kres)
        if kope!='':
                kres=opemat(kres,kope=kope)
        self.ksym=kres
        if kshow:
                self.s()
                
    def simplify(self, kupdate=True, kshow=True, kope=''):
        kres = self.ksym
        kres = simplify(kres)
        kres = opemat(kres, kope=kope)
        if kupdate:
            self.update(kres)
        if kshow:
            self.s()
    def expandexp(self,op='',kshow=True):
        kres=self.ksym
        kres=expandexp(kres,op=op)
        self.ksym=kres
        if kshow:
            self.s()
            
            
    def packexp(self,kshow=True):
        kres=self.ksym
        kres=packexp(kres)
        self.ksym=kres
        if kshow:
            self.s()        
     
    def base2frac(self,kshow=True):
        kres=self.ksym
        kres=base2frac(kres)
        self.ksym=kres
        if kshow:
            self.s()
    
    def simplifyexp(self,kope='',kshow=True):
         
        kres=self.ksym
        kres=simplifyexp(kres,kope=kope)
              
        self.ksym=kres
        if kshow:
            self.s()
    def primefactorexp(self,kope='',kshow=True):
        kres=self.ksym
        kres=primefactorexp(kres)
        if kope!='':
            kres=opemat(kres,kope=kope)
        self.ksym=kres
        if kshow:
            self.s()
            
    def div2mulexp(self,kope='',kshow=True):
        kres=self.ksym
        kres=div2mulexp(kres)
        if kope!='':
            kres=opemat(kres,kope=kope)
        self.ksym=kres
        if kshow:
            self.s()
            
    def simplifyrpow(self, kupdate=True, kshow=True, kope=''):
        kres = self.ksym
        kres = simplifyrpow(kres)
        kres = opemat(kres, kope=kope)
        if kupdate:
            self.update(kres)
        if kshow:
            self.s()

    def simplify_sec(self,kshow=True):
        kres = self.ksym
        if self.Is_Add():
            sres = 0
            mm=fpoly(kres, 'list')
            for i in mm:
                sres += simplify(i)
            self.ksym=sres    
             
        self.update(self.ksym)
        if kshow:
            self.s()

    def factor(self, kvar='', kshow=True, kope=''):
        kres = self.ksym
        if kvar != '':
            g1 = self.fpoly('filt', kvar)
            g2 = kres - g1
            g3 = factor(g1)
            kres = g3 + g2
        else:
            kres = factor(kres)
        kres = opemat(kres, kope=kope)
        self.update(kres)
        if kshow:
            self.s()
    
    def linfactor(self, kvar,kshow=True):
        p1=self.ksym
        p1=linfactor(p1,kvar)
        self.ksym=p1
        if kshow:
            self.s()

    
    def cancel(self, kvar='', kshow=True, kope=''):
        kres = self.ksym
        if kvar != '':
            g1 = self.fpoly('filt', kvar)
            g2 = kres - g1
            g3 = factor(g1)
            kres = g3 + g2
        else:
            kres = cancel(kres)
        kres = opemat(kres, kope=kope)
        self.update(kres)
        if kshow:
            self.s()
            
    def texpand(self, kshow=True, kope=''):
        kres = self.ksym
        kres = expand_trig(kres)
        kres = opemat(kres, kope=kope)
        self.update(kres)
        if kshow:
            self.s()

    def tsimplify(self, kshow=True, kope=''):
        kres = self.ksym
        kres = trigsimp(kres)
        kres = opemat(kres, kope=kope)
        self.update(kres)
        if kshow:
            self.s()
            
    def lexpand(self, kshow=True, kope=''):
        kres = self.ksym
        kres = expand_log(kres,force=True)
        kres = opemat(kres, kope=kope)
        self.update(kres)
        if kshow:
            self.s() 
    
    def lfactor(self, kshow=True, kope=''):
        kres = self.ksym
        kres = logcombine(kres,force=True)
        kres = opemat(kres, kope=kope)
        self.update(kres)
        if kshow:
            self.s()
        
    def exp(self, kshow=True, kope=''):
        kres = self.ksym
        kres = exp(kres)
        kres = opemat(kres, kope=kope)
        self.update(kres)
        if kshow:
            self.s() 
    
    def log(self, kshow=True, kope=''):
        kres = self.ksym
        kres = log(kres)
        kres = opemat(kres, kope=kope)
        self.update(kres)
        if kshow:
            self.s()    

    def opemat(self, kope='', kshow=True):
        kres = self.ksym
        kres = opemat(kres, kope=kope)
        self.update(kres)
        if kshow:
            self.s()

    def opematsec(self, kope='', kshow=True):  # equal to pemat but secuential
        kres = unisymbols(self.ksym)
        kres = opematsec(kres, kope=kope)
        self.update(kres)
        if kshow:
            self.s()

    def simplify_Rpow(self, kshow=True):
        kres = self.ksym
        kres = cut_root_of_pow(kres)

        self.update(kres)
        if kshow:
            self.s()

    def sort(self):
        kres = self.ksym
        klist = self.fpoly('list')
        mm = 0
        for i in klist:
            mm += i
        kres = mm
        self.update(kres)
        self.s()

    def expandExp(self):
        kres2 = self.ksym
        kvar = fpoly(kres2, 'get', 0)
        kexp = fpoly(kres2, 'get', 1)
        keList = fpoly(kexp, 'list')
        mm = 0
        for i in keList:
            mm += kpow(kvar, i)
        return mm

    def get_inside_root(self, kname=''):
        mm = str(self.ksym)
        cc = 0
        kk = mm.find('sqrt')
        kk2 = kk + 4
        sres = ''
        qq = len(mm)
        for i in range(kk2, qq):
            val = mm[i]
            if val == '(':
                cc += 1
            if val == ')':
                cc -= 1
            sres = sres + val
            if cc == 0:
                kres = parse_expr(sres)
                if kname != '':
                    ee = MyEq(kres, kname)
                    return ee
                else:
                    return kres
        return self.ksym
    def get_type(self):
        return type(self.ksym)
    def get_nume(self):
        return numer(self.ksym)

    def get_deno(self):
        return denom(self.ksym)

    def diffValue(self, kval=''):
        if kval == '':
            kval = self.var2
        kres = diff(self.ksym, kval)
        return kres

    def fpoly(self, kopt='', op2='', op3=''):
        kres = self.ksym
        return fpoly(kres, kopt=kopt, op2=op2, op3=op3)
        
    def getexponent(self,kname=''):
        # return exponente from monomie expresion
        return get_expo(self,kname=kname) 
    def get_expo(self,kname=''):  # return exponente from monomie expresion
        ksym = self.ksym
        if Is_Pow(ksym):
            mm = fpoly(ksym, 'list')
            kres= mm[1]
            if kname!='':
                ee=MyEq(kres,kname)
                return ee
            else:
                return kres
        else:
            self.s()
            
    def get_killexpo(self):
        ksym = self.v
        if Is_Mono(ksym):
            mm = fpoly(ksym, 'list')
            return mm[1]

    def get_sitem(self, vv=[]):
        mm = self.list()
        qq = len(mm)
        ksum = 0
        for i in range(qq):
            if i in vv:
                ksum += mm[i]
        return ksum

    def part(self, vec):

        kres = self.ksym
        try:
            return part(kres, vec)
        except:
            return kres

    def solvediff(self, kvar, kd='', kremp=False, kope='', korden=''):
        return self.solve(kvar=kvar.diff(), kd=kd, kremp=kremp, kope=kope, korden=korden)

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
        kres = self.ksym - eqv
        if len(kwargs) > 0:
            kres = self(**kwargs) - eqv
        nee = MyEq(kres, kname=svar, kshow=False, kope=kope,)
        kres=nee.ssolve(svar,korden=korden,kshow=kshow,kope=kope)
        self.set(svar,kres,kshow=kshow,kope=kope)
        return kres


    def get_symbol_in_ee(self, svar=''):
        mm = self.free()
        for i in mm:
            if i.name == svar:
                return i
        return svar

    def left(self):
        return symbols(self.name)

    def right(self):
        return  self.ksym

    def ssolve(self, kname,korden='',kshow=True,kope=''):

        kvar = self.get_symbol_in_ee(svar=kname)
        return self.solve(kvar, kname=kname,korden=korden,kshow=kshow)

    ###########################################
    #            SOLVE
    ###########################################

    #   solve()
    def solve(self,*args,kshow=True,**kwargs):
        kres=self.ksym
        
        if len(kwargs)>0:
            for key, value in kwargs.items():
                if key==self.name:
                    nnvar=key
                    nvar=symbols(nnvar)
                    kres=nvar-kres
                    
            kres=real_subs(kres,**kwargs)
        kpositive=False
        kall=False
        kupdate=False
        kpos=1
        vop=['all','positive','0','2','3','4','update']
        vnum=['0','2','3','4']
        wvar=args[0]
         
        kname=str(wvar)
        
        fresp=solve(kres,wvar)
      
        for i in args:
            if type(i)==str:

                if i=='positive':
                    kpositive=True
                if i=='all':
                    kall=True
                if i=='update':
                    kupdate=True
                    
                if i in vnum:
                    kpos=int(i)
                    
                    
        
        qq=len(fresp)
        if qq==1:
            kres=fresp[0]
        elif  (not kpositive  and  not kall):
            kres=fresp[kpos]
            
            
         
        if qq==1 or  (not kpositive  and  not kall):
             
            if kname=='':
                return kres
            else:
                ee=MyEq(kres,kname=kname,var2=self.var2,kshow=kshow)
                 
            
        if kall:
                cc=1
                kkres=[]
                for i in fresp:
                    kname1 = kname + str(cc)
                    kres1 = MyEq(i, kname1,var2=self.var2)
                    cc += 1
                    kkres.append(kres1)
                return kkres
        if kpositive:
            mm=fpoly(kres2,'list')
            
             
            for i in mm:
                if i==-1:
                    kres2=-1*kres2
                    
            ee=MyEq(kres2,kname=kname,var2=self.var2,kshow=kshow)        
                
        if kupdate:
            self.upgrade(ee,kshow=kshow)
        return ee

    def solveset(self, *args,kshow=True,**kwargs):
        kres=self.solve(*args,kshow=True,**kwargs)
        self.upgrade(kres)
        return kres
        
        



    def solveB(self, kvar, kname='', kremp=False, kope='', korden='', ktype='P', var2='', Bag='',kpositive=False,noncero=False,kshow=True):
        mm = self.ksym
        self.ksym = opemat(mm, 's')
        if type(kvar) == 'str' and kname == '':
            kname = kvar
            kvar = parse_expr(kname)

        if str(kvar) == 'None':
            kname = symbols(kname)
        keq = self.ksym
        # kvar=unisymbols(kvar)
        if noncero and kname!='':
            nV=symbols(self.name)
            keq=keq-nV
            kname=nV.name
        kres = csolve(keq, kvar, kope=kope, korden=korden)
        if kres == []:
            kres = csolve(keq, kname, kope=kope, korden=korden)
            if kres == []:
                kres = presolve(self.ksym, kvar)

        if self.type == 'Ph':
            if not str(kvar) in 'alphaalpha1alpha2alpha3betha1betha2betha3':
                self.Pobj.store_val(kvar, kres)
                mm = self.Pobj.F
                qq = len(mm)
                for i in range(qq):
                    pkres = mm[i][0]
                    try:
                        pkres = pkres.subs(kvar, kres)
                    except:
                        done = False
                    mm[i][0] = pkres
                self.Pobj.F = mm
        if kname != '':
            if type(kres) == list:
                if len(kres) > 1:
                    if kpositive:
                        kres1 = MyEq(kres[1], kname, var2=self.var2)
                        return kres1

                    cc = 1
                    kkres = []
                    for i in kres:
                        kname1 = kname + str(cc)
                        kres1 = MyEq(i, kname1, kshow=False, var2=self.var2)
                        if Bag != '':
                            kres1.upBag(Bag, kope=kope)

                        kres1.s()
                        cc += 1
                        kkres.append(kres1)

                    return kkres

            else:

                kres = MyEq(kres, kname, kshow=False, var2=self.var2)
                if Bag != '':
                    kres.upBag(Bag, kope=kope)
                if kshow:
                    kres.s()

                return kres


        else:
            return kres

    def ssolveR(self, kname):
        kvar = self.get_symbol_in_ee(svar=kname)
        return self.solveR(kvar, kd=kname)

    def solveR(self, kvar, kd='', kremp=False, kope='', korden='', Bag=''):
        if Bag != '':
            self.upBag(Bag=Bag, kshow=False)
        keq = self.ksym

        try:
            kres = csolveR(keq, kvar, kope=kope)
        except:
            kres = csolve(keq, kname, kope=kope, korden=korden)
        if kd != '':
            return MyEq(opemat(kres, kope=kope), kd, var2=self.var2)
        else:
            return kres

    def toFoat(self):
        kres = self.ksym
        try:
            kres = float(self.ksym)

            self.update(kres)
            self.s()
        except:
            self.s()

    def toEqual(self, items, kkname=''):
        mm = self.list()
        kres1 = 0
        kres2 = 0
        qq = len(mm)
        for i in range(qq):
            if i in items:
                kres1 += mm[i]
            else:
                kres2 += mm[i]
        kname = self.name
        kname1 = kname + '1'
        kname2 = kname + '2'
        kname1 = MyEq(kres1, kname1, kshow=False)
        kname2 = MyEq(kres2, kname2, kshow=False)
        return MyEqEq(kname1, kname2, kname=kkname)

    def opematsec(self, kope=''):  # equal to pemat but secuential
        kres = unisymbols(self.ksym)
        kres = opematsec(kres, kope=kope)
        self.update(kres)
        return kres

    def opemat_deno(self, kope=''):
        kres = unisymbols(self.ksym)
        kres = opemat_deno(kres, kope=kope)
        self.update(kres)
        return kres

    def opemat_nume(self, kope=''):
        kres = unisymbols(self.ksym)
        kres = opemat_nume(kres, kope=kope)
        self.update(kres)
        return kres

    def solve_tan(self, alpha):
        c = tan(alpha)

        kres = self.ksym
        kres.subs(sin(alpha), c / rpow(c * c + 1, 2))
        kres.subs(cos(alpha), 1 / rpow(c * c + 1, 2))
        kres1 = csolve(kres, tan(alpha))
        return kres1

    ###########################################
    #             Polinomial
    ###########################################
    def quotient(self,ksym,kname=''):
        
        P=self.ksym
        Q=ksym
        if type(Q)==MyEq:
            Q=ksym.ksym
        kres=quo(P,Q)
        if kname!='':
            kres=MyEq(kres,kname)
            return kres
        else: 
            return kres
        
    def rem(self,ksym,kname=''):
        
        P=self.ksym
        Q=ksym
        if type(Q)==MyEq:
            Q=ksym.ksym
        kres=rem(P,Q)
        if kname!='':
            kres=MyEq(kres,kname)
            return kres
        else: 
            return kres
        
    def coef_sum(self,var2=''):
        if var2=='':
            var2='x'
        kres=self.ksym
        kres=kres.subs(var2,1)
        return kres
        
    def get_coef(self,ksym,kname=''): # ee= x**2(a+b)+x*4*b+9
                           # get_coef(x) return 4*b
        mm=self.list()
        mmv=[]
        sres=str(ksym)
        for i in mm:
            mmv.append(i/ksym)
        for i,j in zip(mmv,mm):
            if not sres in str(i):
                kres=j
                kres=kres.subs(ksym,1)
        if kname!='':
            return MyEq(kres,kname=kname,var2=self.var2)
        else:
            return kres
            
    def get_ter_inde(self,ksym,kname=''):
        kres=0
        svar=str(ksym)
        mm=self.list()
        for i in mm:
            if not svar in str(i):
                kres+=i
        if kname!='':
            return MyEq(kres,kname=kname,var2=self.var2)
        else:
            return kres        
    
    def degree(self, kvar=0):
        kres = self.ksym
        return degree(kres, gen=kvar)

    def degree_list(self):
        kres = self.ksym
        return degree_list(kres, gen=kvar)
        

    def main_coef(self):
        kres = self.ksym
        return LC(kres)
        
    def main_monomio(self):
        kres = self.ksym
        return LM(kres)
    
    def main_term(self):
        kres = self.ksym
        return LT(kres)
        
    def coef_listK(self,vx=''):
        if vx=='':
            vx=self.var2
        mm=self.list()
        vlist=[]
        for i in mm:
            kres=i.subs(vx,1)
            vlist.append(kres)
        return vlist      
        
    def reduce(self):
        kres = self.ksym
        kres = expand(kres)
        kres = simplify(kres)
        kres = factor(kres)
        kres = simplify(kres)
        try:
            kres = cut_root_of_pow(kres)
            self.update(kres)
        except:
            self.update(kres)
        return kres
        
        
    ###########################################
    #             Algebra
    ########################################### 
    def squaresum(self,p1=0,p2=0):
        kres=self.ksym
        if p1!=0:
            kres2=squaresum(kres,p1=p1,p2=p2)
            self.ksym=kres2
            self.s()
        else:
            self.s()
            
            
    # def expandexp(self):
        # kres=self.ksym
        # try:
            # kres=expandexp(kres)
            # self.ksym=kres
        # except:
            # pass
        # self.s()    
            
    def powexpand(self,op='',kshow=True):
        '''
        input (x**(a*b))   ---->   return(x**a)**b
        input (x**(a*b),b)   ---->   return(x**b)**a
        '''
        kres=self.ksym
        self.ksym=powexpand(kres,op=op)
        if kshow:
            self.s()
    
    def killsqrtpow(self,kshow=True):
      
        kres=self.ksym
        self.ksym=killsqrtpow(kres)
        if kshow:
            self.s()
    
    def opematexp(self,kope=''):
        '''
        apply opemat only in exponent monomie
        '''
        ksym=self.ksym
        if Is_Pow(ksym):
            mm=self.list()
            kres=mm[1]
            kres=opemat(kres,kope=kope)
            self.ksym=kpow(mm[0],kres)
            
        if Is_Add(ksym):
            mm=self.list()
            mkres=0
            for i in mm:
                if Is_Pow(i):
                    mm2=fpoly(i,'list')
                    kres=mm2[1]
                    kres=opemat(kres,kope=kope)
                    kres=kpow(mm2[0],kres)
                else:
                    kres=i
                mkres+=kres 
             
            self.ksym=mkres     
        self.s()
    
    def mulexpo(self,kope='',kshow=True): # idem expanexp
        kres=self.ksym
        kres=mulexpo(kres)
        if kope!='':
            kres=opemat(kres,kope=kope)
        self.ksym=kres    
        if kshow:
            self.s()
        
    def expandexpA(self,kope=''):
        '''
        input (x**a)**b
        return(x**(a*b))
        '''
        ksym=self.ksym
        kres=expandexpA(ksym)
        self.ksym=kres
        if kope!='':
            self.opemat(kope=kope)
        self.s()
        
    def factorexp(self,kope=''):
        '''
        input (x**a)**b
        return(x**(a*b))
        '''
        ksym=self.ksym
        kres=factorexp(ksym)
        self.ksym=kres
        if kope!='':
            self.opemat(kope=kope)
        self.s()
    
    def coef_list(self,var2=x):  # return coef list complete according max degree
        '''
        ee=x*x*a+x*b+c
        return [a,b,c]
        
        ee=a*x*x-b
        return [a,0,-b]
        
        '''
        
        kres=self.ksym
        return  coef_list(kres,var2)   
    
    def get_seudocofactor(self,ksym,kname='',var2=x):
        val1=self.ksym
         
        if type(ksym)==MyEq:
            ksym=ksym.ksym
            
        kres,veccc=get_seudocofactor(val1,ksym,var2)
         
         
        if kname!='':
            return MyEq(kres,kname=kname,var2=var2)
        else:
            return kres 
    
    def sortdegree(self,var2='',kope=''):
        '''
        input a*x*x+b*x*x+c*x+d
        return (a+b)x*x+c*x+d
        '''
        if var2=='':
            var2=self.var2
        ksym=self.ksym
        kres=sortdegree(ksym,var2=var2)
        self.ksym=kres
        if kope!='':
            self.opemat(kope=kope)
        self.s()
        
    
    def addexp(self,kope=''):  # input:  x**a * x**b  return x**(a+b)
        kres=self.ksym
        kres=addexp(kres)
        self.ksym=kres
        if kope!='':
            self.opemat(kope=kope)
        self.s()
    
    
    def insidepar(self,ssym=''):
        '''
        args:
            ssym dfault is '' return the first () founded 
            ssym = str( where search..., 'x(', 'sqrt(', 'log(' etc 
        return expresion inside  select function     
        '''
        ksym=self.ksym
        return insidepar(ksym,ssym=ssym)
    
    ###########################################
    #             Reduccion Algoritmo
    ###########################################
    def factorize(self,ksym,kshow=True):
     
       # ee=3*x**7+2*x**3
       #ee.factorize(x**3)
       #return x**3*(3*x**4+2)
       
     
        kres=self.ksym
        kres=factorize(kres,ksym)
        self.ksym=kres
        if kshow:
            self.s()
        
        
    def reduFac(self, kop='', kshow=True):  # retun Eq=0,= a*b+c*b.. = b(a+c)..=0  then return (a+c)

        kres = self.ksym
        kres = factor(kres)
        kres2 = 1
        if Is_Mono(kres):
            kres = numer(kres)

            mm = fpoly(kres, 'list')
            for i in mm:
                if Is_Poly(i):
                    kres2 = kres2 * i

        if kres2 != 0 and kres2 != 1:
            self.ksym = kres2
        else:
            self.ksym = kres
        if kshow:
            self.s()

    def factorSec(self, *args, kfiltro='.', kshow=True):
        kres = self.ksym
        for i in args:
            kvar=i
            kres = factorSec(kEq=kres, ksym=kvar, kfiltro=kfiltro)
        self.update(kres)
        if kshow:
            self.s()

    def factorSecV(self, *args):
        for i in args:
            self.factorSec(i, kshow=False)
        self.s()

    def get_cofactor(self, kfactor):
        mm = self.list()
        for i in mm:
            j = fpoly(i, 'list')
            if kfactor in j:
                kres = i / kfactor
                return kres

    def grupFac(self, ksym, kfiltro='.', kshow=True):
        kres = self.ksym
        kres = grupFac(kres, ksym, kfiltro='.')
        self.update(kres)
        if kshow:
            self.s()

    def get_factor_with(self, kx, kcomple=True):
        eqq = self.ksym
        return get_factor_with(eqq, kx=kx, kcomple=kcomple)

    def monoFactor(self, kval):
        kres = self.ksym
        if Is_Poly(kres):
            mm1 = 0
            mm2 = 0
            mm = fpoly(kres, 'list')
            for i in mm:
                if Is_Mono(i) and Is_Mul(i):

                    vmm = fpoly(i, 'list')
                    if kval in vmm:
                        newm = simplify(i / kval)
                        mm1 += newm
                    else:
                        mm2 += i
            kres2 = kval * mm1 + mm2
            self.update(kres2)
            return kres2

        return kres

    def fixRootPow(self, kksym):
        self.setValue(rpow(kpow(kksym, 2), 2), kksym)

    def fix_sqrtPow(self):
        ksym=self.ksym
        nkasym=fix_sqrt2pow(ksym)
        self.ksym=nkasym
        self.s()

    def find_root(self, ksym, x1, x2, xx):
        xx = np.linspace(x1, x2, xx)
        yy = self.evalueArray(ksym, xx)
        mm = []
        qq = len(yy)
        for i in range(qq - 1):
            if (yy[i] > 0 and yy[i + 1] < 0) or (yy[i] < 0 and yy[i + 1] > 0):
                mm.append((xx[i] + xx[i + 1]) / 2)
        qq = len(mm)
        for i in range(qq):
            e1 = MyEq(self.v, str(ksym), kshow=False)
            e1.setValue(ksym, mm[i])

    def upBag(self, Bag, kname='', kshow=True, kope=''):
        vs = Bag.vmain
        vv = Bag.vsolve

        for i, j in zip(vs, vv):
            self.set(i, j, kshow=False)
            self.opemat(kope, kshow=False)
        if kshow:
            self.s()

    def upTriang(self, angul, T3, kope=''):
        if type(angul) == MyTriang:
            angu2 = angul
            angul = T3
            T3 = angu

        v1 = [sin(angul), cos(angul), tan(angul), kpow(sin(angul), 2), kpow(cos(angul), 2), kpow(tan(angul), 2)]
        v2 = [T3.sin(), T3.cos(), T3.tan(), kpow(T3.sin(), 2), kpow(T3.cos(), 2), kpow(T3.tan(), 2)]

        for i in range(len(v2)):
            v2[i] = opemat(v2[i], kope=kope)

        self.set(v1, v2)

    def evalueBag(self, bag, kope=''):
        e1 = MyEq(self.v, self.name, kshow=False)
        e1.upBag(bag, kope=kope)

    def setVal_from_bag(self, bag, kshow=False, kope=''):
        for i, j in zip(bag.dataS, bag.dataV):
            self.setValue(i, j, kshow=kshow, kope=kope)
        if kshow:
            self.s()

    #  Algoritmos de reparacion

    def fix_reduc(self):
        kres = self.ksym
        try:
            kres = fix_reduc(kres)
            self.update(kres)
            self.s()
        except:
            self.s()
    ###############################
    #  Diff
    ###############################
    
    
    def diffEq(self, kname='', var2='', ktype='P', typeD=1):
        self.oldprimi = self.primi
        kres = self.ksym
        if var2 == '':
            var2 = self.var2
        if typeD == 2:
            return MyEq(Derivative(kres, var2), kname=kname, var2=var2, ktype='Diff')

    def changediffI(self, newvar2, newvalue='', x1='', x2=''):
        kres = self.primi
        if newvalue != '':
            kres = kres * newvalue
            self.primi = kres
        self.var2 = newvar2
        if x1 != '':
            self.x1 = x1
            self.x2 = x2
        if self.x1 == '':
            self.ksym = Integral(self.primi, self.var2)
        else:
            self.ksym = Integral(self.primi, (self.var2, self.x1, self.x2))
        self.s()

    def changediff(self, nkvar, nvalue='', x1='', x2=''):
        if self.type != 'I':
            return
        if nvalue == '':
            xx = nkvar
        else:
            xx = nvalue
        if x1 != '':
            self.x1 = x1
            self.x2 = x2
        else:
            if self.x1 != '':
                x1 = self.x1
                x2 = self.x2
                var2 = self.var2

                if str(var2) in str(x1):
                    self.x1 = x1.subs(var2, xx)
                if str(var2) in str(x2):
                    self.x2 = x2.subs(var2, xx)
        kprimi = self.primi

        kprimi = kprimi.subs(self.var2, xx)
        if str(self.var2) not in str(kprimi):
            kprimi = kprimi * nvalue
        else:
            kprimi = kprimi * xx.diff(nkvar)

        self.primi = kprimi
        self.var2 = nkvar
        if self.x1 == '':
            self.ksym = Integral(kprimi, self.var2)
        else:
            self.ksym = Integral(kprimi, (self.var2, self.x1, self.x2))

        self.s()
    
    
    ###############################
    #  Integral
    #######################    
    
    
    
    
                
    def doitI(self, kname='', kope='', kshow=True,C1=C1, **kwargs):
            
            kres=self.ksym
            
            if len(kwargs)>0:
                for key, value in kwargs.items():
                    kres=kres.subs(parse_expr(key),value)
                    
            kres2=kres.doit()   
             
            kres = opemat(kres2, kope=kope)
            if 'otherwise' in str(kres): 
                try:
                    kres=fix_otherwise(kres)
                except:
                    pass
             
            if self.name == 'Vo':
                if sign(kres) == -1:
                    kres = -1 * kres
            kres=kres+C1 
            
            
            if kname != '':
                ee = MyEq(kres, kname=kname,kshow=kshow)
                return ee
            else:
                self.update(kres)
                if kshow:    
                    self.s()

    def doit(self, kname='', kshow=True, c1=0, **kwargs):
        if self.type == 'diff':
            if len(kwargs) == 0:
                kres = self.ksym
                var2 = self.var2
                kres2 = diff(kres, var2)
                self.update(kres2)
                self.type = 'F'
                s1 = self.name
                s2 = alphaname(var2)
                kname = s1 + '_{(' + s2 + ')}'
                self.name = kname
                self.s()
                return

        if self.type == 'I':
            if len(kwargs) > 0:
                ee = MyEq(self.ksym, kshow=False)
                kres = ee(**kwargs)
                ee = MyEq(kres, kshow=False)
                kres = ee.doit() + c1

            else:
                kres = self.ksym
                kres = kres.doit() + c1
            if kname != '':
                ee = MyEq(kres, kname=kname)
                return ee
            else:
                self.type = 'P'
                self.update(kres)
                if kshow:
                    self.s()
        elif self.type == 'Diff':
            kres = self.ksym
            kres = kres.diff(self.var2)
            self.update(kres)
            self.type = 'P'
            self.name = 'd' + self.name
            self.s()


        else:
            self.s()

    def get_diff(self, var2='', kname=''):
        if var2 == '' and kname == '':
            var2 = self.var2
        if var2 != '':
            if type(var2) == str:
                kname = var2
                var2 = self.var2
                if type(kname) != str:
                    var2 = kname

        kres = self.ksym
        kres = diff(kres, var2)
        if kname != '':
            ee = MyEq(kres, kname, var2=var2)
            return ee
        else:
            return kres
    def get_primitive(self):
        return self.primi_eq
        
    def get_dprimitive(self):
        return self.diff_eq
        
         
    def get_primitive(self):
        return self.primi_eq        
        
    def kdiff(self, kname='', var2='', kope='', kupdate=False, ktype=''):
        if ktype == '':
            ktype = self.type
        kres = self.ksym
        if var2 == '':
            var2 = self.var2
        kres = self.ksym
        kres2 = diff(kres, var2)
        return kres2
        if kname != '':
            ee = MyEq(kres2, kname=kname, var2=var2, ktype=ktype)
            return ee
        else:
            return kres2

    def update_inte(self):
        if self.type == 'I':
            kres = self.primi
            if self.x1 == '':
                self.ksym = Integral(kres, self.var2)
            else:
                self.ksym = Integral(kres, (self.var2, self.x1, self.x2))

    def fac_integral(self):
        ksym = self.primi
        var2 = self.var2

        if Is_Mono(ksym) and Is_Mul(ksym):
            kres = [x for x in ksym.args if str(var2) not in str(x)]
            kres2 = [x for x in ksym.args if str(var2) in str(x)]
            mono1 = 1
            mono2 = 1
            for i in kres:
                mono1 = mono1 * i
            for i in kres2:
                mono2 = mono2 * i
            self.primi = mono2
            kfac = mono1
            self.primi = mono2
            if self.x1 == '':
                self.ksym = Integral(mono2, self.var2)
            else:
                self.ksym = Integral(mono2, (self.var2, self.x1, self.x2))
            self.Mul(mono1)
        else:
            self.s()

    def maximun(self, ksym, ksave=True):
        '''
        Return maxumin value of ksym if self=0
        '''
        if ksave:
            savework()
        kres = self.ksym
        kres = diff(kres, ksym)
        kres = solve(kres, ksym)
        kname = ksym.name
        if len(kres) == 1:
            return MyEq(kres[0], kname)
        else:
            vres = []
            cc = 1
            for i in kres:
                nname = kname + str(cc)
                vres.append(MyEq(i, kname=nname))
                cc += 1
            return vres

            #  Geometry

    def kplot(self, ksym, x1, x2, x3):
        xx = np.linspace(x1, x2, x3)
        yy = [self(x) for x in xx]
        plot(xx, yy)

        ylabel(self.name)
        xlabel(str(self.var2))
        ktitle = self.name + '(' + str(str(self.var2)) + ')'
        title(ktitle)
        show()
    def length_arc(self, x1='', x2='', x='', ksolve=True):
        if x1 == '':
            x1 = self.x1
            x2 = self.x2
        if x == '':
            x = self.var2
        y = MyEq(self.ksym, 'y', x1=x1, x2=x2, varx=x, ktype='Diff')
        y.doit()
        if ksolve == 'dL':
            y.name = 'dL'
            y.s()
            return y
        y2 = y * y
        y2 = opemat(y2, 'ef')
        y3 = rpow(1 + y2)
        y3 = opemat(y3, 'r')
        L = MyEq(y3, 'L', x1=x1, x2=x2, varx=x, ktype='Diff')

        L.integral()
        if ksolve:
            return L.doitI()
        else:
            return L


    ################################################
    ##                Diferencial                ###
    ################################################


    def dsolve(self, kname='', C1='', C2=''):
        kres = self.EqDiff
        dres = dsolve(kres)
        mdres = dres.args
        kres = mdres[1]
        if C1 != '':
            kres = kres.subs('C1', C1)
        if C2 != '':
            kres = kres.subs('C2', C2)
        if kname != '':
            ee = MyEq(kres, kname=kname, var2=self.var2, ktype='F')
            return ee
        else:
            ee = MyEq(kres, kname=alphaname(self.var1), var2=self.var2, ktype='F')
            return ee

    def Diff(self, kname='', var2='', ktype='P'):
        if var2 == '':
            var2 = self.var2
        kres = self.ksym
        kres = kres.diff(var2)
        if kname == '':
            return kres
        else:
            ee = MyEq(kres, kname=kname, var2=var2, ktype=ktype)
            return ee

    def applyDiff(self, var2=''):
        kres = self.ksym
        if var2 == '':
            var2 = self.var2
        kres = kres.diff(var2)
        self.update(kres)
        self.s()

    def fab_diff(self):
        var2=self.var2
        kname= self.name
        f1,f2,f3=[Function('f1')(var2),Function('f2')(var2),Function('f3')(var2)]
        nvarf=['f1','f2','f3']
        vvarf=[f1,f2,f3]
        kres= self.ksym
        ee2= self.xcopy('ee2',kshow=False)
        qq=len( self.varf)
        vvarf=vvarf[0:qq]
        nvarf=nvarf[0:qq]
        varf= self.varf
        for i,j in zip(varf,vvarf):
            ee2.set(i,j,kshow=False)
        kres=ee2.ksym
        kres=kres.diff(var2)
        ee3=MyEq(kres,'ee3',kshow=False)
        ee3.primi_eq=self.primi_eq
         
        for i,j in zip(vvarf,varf):
            ee3.setdiff(i,1,kshow=False)
            ee3.set(i,j,kshow=False)
        ee3.name= diffname(self.var1,var2)
        ee3.primi_diff_eq=diff(ee3.primi_eq,ee3.var2)

        return ee3
        
    def Diff2diff(self,kshow=True):
        fksym=self.flatdiff
        for i in self.varf:
            kvar=str(i)
            kvar2='d'+kvar
            var2=self.var2
            kres=self.primitiva
            f=Function(i)(var2)
            news=symbols(alphaname(kvar2))
            fksym=fksym.subs(diff(f,var2),news)
        return(fksym)    

        
    def Diff2diff(self,kres): # ksym,kvar,var2
        var2=self.var2
        kvar=self.varf
        for i in kvar:
            f=Function(str(i))(var2)
            df=diff(f)
            kname='d'+alphaname(i)
            nf=symbols(kname)
            kres=kres.subs(df,nf)
        return kres
    def p(self):
        return self.primitiva
        
    def pf(self):
        try:
            return primi2func(self.ksym,self.var2,self.varf)
        except:
            return self.ksym
        
    def Fdiff(self,*args):
        kname=''
        kupdate=False
        for i in args:
            if type(i)==str:
                if i=='update':
                    kupdate=True
                else:
                    kname=i
        if self.varf==[]:
            kres=diff(self.ksym,self.var2)
            primi=0
        else:
            primi=diff(self.primitiva,self.var2)
            kres=Diff2diff(primi,self.varf,self.var2)
            
        if kname=='':
            if kupdate:
                self.ksym=kres
                self.primitiva=primi
                
                self.s()
                 
            else:
                
                return kres
        else:
            ee=self.xcopy(kname,kshow=False)
            ee.ksym=kres
            ee.primitiva=primi
            ee.dtype='diff'
            ee.s()
            return ee
            
        
    def diff(self, *args, respect='',kupdate=False,kshow=True,**kwargs):
               
        if self.varf==[]:
              self.primi_eq=self.ksym
              self.diff_eq=diff(self.ksym,self.var2)
              self.primi_diff_eq=diff(self.ksym,self.var2)
            
        if self.varf!=[]:


            if self.primi_eq=='':
             self.primi_eq=self.get_primitive()

            if self.primi_diff_eq=='':    
             ee=self.xcopy('ee',kshow=False)
             kres=ee.primi_eq
             kres=kres.diff(ee.var2)
             self.primi_diff_eq=kres
         

            
        
        if len(args)==0 and len(kwargs)==0:
            self.dtype='diff'
            kname=self.name
            ee3= self.fab_diff()
            self.ksym=ee3.ksym
            self.name=diffname(self.name,self.var2)
            if respect!='':
                self.name= diffname(self.name,respect)
            if self.varf!=[]:
                self.dtype='diffd'

            self.diff_eq=self.ksym

            self.s()

            return

        elif len(args)==1 and len(kwargs)==0 and type(args[0])==str and self.varf!=[]:
            kname=self.name
            ee3= self.fab_diff()
            ee3.dtype='diff'
            ee3.name= diffname(kname,ee3.var2)
            if respect!='':
                ee3.name= diffname(kname,respect)
                ee3.diffEq=self.diffEq
                 
            if ee3.varf!=[]:
                 ee3.dtype='diffd'
            
                 
            ee3.s()
            ee3.origen=self.ksym
            ee3.origen=self.ksym
            ee3.vfunc=self.vfunc
            return ee3
        else:
            kname = ''
            var = self.var2
            kres = self.ksym
            for i in args:
                if type(i) == str:
                    kname = i
                else:
                    var = i
            kres = diff(kres, var)
            if len(kwargs) > 0:
                ee = MyEq(kres, kshow=False)
                for key, value in kwargs.items():
                    ee.set(parse_expr(key), value, kshow=False)
                kres = ee.ksym

            if kname != '':
                ee = MyEq(kres, kname,var1=self.var1,kshow=kshow)
                ee.dtype='diff'
                ee.var2=self.var2
                ee.vfunc=self.vfunc
                return ee
            else:
                return kres

    def simple_form(self,kfunc):
        self.setdiff(kfunc,1,kshow=False)
        self.set(kfunc, self.var1)

    def flatdiff(self,kfunc,kshow=True):
        self.setdiff(kfunc,1,kshow=False)
        self.set(kfunc, self.var1,kshow=kshow)

    def diff_parcial(self,*args):
        kname=''
        vdiff=[]
        for i in args:
            if type(i)==str:
                kname=i
            else:
                vdiff.append(i)


        kres=self.ksym
        var2=self.var2
        mm=vdiff
        mm=list(mm)
        mm.append(var2)
        mm=tuple(mm)
        ddiff=sym2Function(*mm)
        for i,j in zip(vdiff,ddiff):
            kres=kres.subs(i,j)

        kres=diff(kres,var2)
        if kname!='':
            return MyEq(kres,kname,var2=var2)
        else:
            self.ksym=kres
            self.s()

    def killdiff(self,kfunc):

        self.setdiff(kfunc,1)

    def killfunc(self,kfunc):

        self.set(kfunc, self.var1)

    def killfunc(self,kfunc):

        self.set(kfunc, self.var1)

    def Length(self, kvar):
        kres = self.v
        kdiff = kres.diff(kvar)
        kres = rpow(1 + kpow(kdiff, 2))
        return kres

    #  Algoritmos Trigonometricos
    def sin2cos(self, angu, korden=2, kope='', kshow=True):
        self.set(kpow(sin(angu), 3), (1 - kpow(cos(angu), 2)) * cos(angu), kshow=False)
        kres = self.ksym
        kres = sin2cos(kres, angu=angu, korden=korden, kope=kope)
        self.update(kres)
        if kshow:
            self.s()
        else:
            pass

    def cos2sin(self, angu, korden=2, kope='', kshow=True):
        self.set(kpow(cos(angu), 3), (1 - kpow(sin(angu), 2)) * sin(angu), kshow=False)
        kres = self.ksym
        kres = cos2sin(kres, angu=angu, korden=korden, kope=kope)
        self.update(kres)
        if kshow:
            self.s()

    def get_diffEq(self, vx=''):

        kname = self.name + 'd'
        Vx = self.xx

        if Vx != '':

            kres = self.ksym

            kres = diff(kres, Vx)
            kname = kname + '_' + str(Vx)
            self.EqDiff = MyDiff(kres, Vx, kname)
            return self.EqDiff
        elif vx != '':
            self.xx = vx
            kres = self.ksym
            kname = kname + '_' + str(vx)
            kres = diff(kres, vx)
            self.EqDiff = MyDiff(kres, vx, kname)
            return self.EqDiff
        else:
            sE(['first set variable MyEq.xx= Valeue'])

    def Integral_eQ(self, kname='', var2='', x1='', x2='', ktype='P'):
        ksym = self.ksym

        if var2 == '':
            var2 == self.var2
        if x1 == '':
            x1 = 0
            x2 = var2
        else:
            x1 = x1
            x2 = x2
        kres = Integral(ksym, (var2, x1, x2))
        if kname == '':
            return kres.doit()
        else:
            ee = MyEq(kres, kname=kname, var2=var2, ktype='I')
            return ee

    def Area(self, kname='', var2='', x1='', x2='', kope=''):
        if var2 == '':
            var2 = self.var2
        if x1 == '':
            x1 = self.x1
            x2 = self.x2

        ee = MyEq(self.ksym, kname=kname, var2=var2, x1=x1, x2=x2, ktype='I', kshow=False)
        if kname != '':
            ee.s()
            ee.doitI()
            return ee
        else:
            ee.doitI(kshow=False)
            return ee.ksym

    def area_value(self):
        ee = MyEq(self.ksym, var2=self.var2, x1=self.x1, x2=self.x2, kshow=False)
        ee.area(kshow=False)
        ee.doit(kshow=False)
        return ee.ksym

    def area(self, var2='', x1='', x2='', kshow=True):
        if x1 != '':
            self.x1 = x1
            self.x2 = x2
        if var2 == '':
            var2 = self.var2
        else:
            self.var2 = var2
        self.primi = self.ksym

        if self.x1 == '':
            kres = Integral(self.ksym, self.var2)
        else:
            kres = Integral(self.ksym, (self.var2, self.x1, self.x2))
        self.type = 'I'
        self.update(kres)
        if kshow:
            self.s()

    def Integral(self, kname='', var2='', x1='', x2='', ktype='P'):
        ksym = self.ksym

        if var2 == '':
            var2 == self.var2
        if x1 == '':
            x1 = self.x1
        if x2 == '':
            x2 = self.x2

        kres = Integral(ksym, (var2, x1, x2))
        if kname == '':
            self.ksym = kres
            self.s()
        else:
            ee = MyEq(kres, kname=kname, var2=var2, ktype='I')
            return ee

    def solveIntegral(self, kname='', var2='', x1='', x2='', kupdate=True):
        if x1 == '':
            x1 = self.x1
        if x2 == '':
            x1 = self.x2
        if var2 == '':
            var2 = self.var2
        ksym = self.ksym
        kres = ksym.doit()
        if kname != '':
            ee = MyEq(kres, kname=kname, var2=var2)
            return ee
        else:
            if kupdate:
                self.ksym = kres
                self.s()
            else:
                return kres

    def get_InteEq(self, vx='', kname=''):

        if kname == '':
            kname = 's' + self.name
        else:
            kname = kname
        Vx = self.xx

        if Vx != '':

            kres = self.ksym

            self.EqInte = MyInteger(kres, Vx, kname)
            return self.EqInte
        elif vx != '':
            self.xx = vx
            kres = self.ksym

            self.EqInte = MyInteger(kres, vx, kname)
            return self.EqInte
        else:
            sE(['first set variable MyEq.xx= Value'])

    def convMyFunc(self, knom='', vV=''):
        ksym = self.v
        if vV == '':
            vV = self.free()
        elif type(vV) == tuple:
            vV = list(vV)
        ksym = self.v
        kk = MyFunc(knom, ksym, vV)
        return kk

    def newEq(self, kname='', kope='', **kwargs):
        ee = MyEq(self.ksym, kname=kname, kshow=False)
        if len(kwargs) > 0:
            for key, value in kwargs.items():
                ee.set(parse_expr(key), value, kshow=False, kope=kope)

        ee.s()
        return ee

    def findSubFunc(self, sval, inside=''):
        ksym = self.ksym
        kini = 0
        kini2 = 0
        sroot = []
        done = 0
        while kini < len(str(ksym)) and kini2 != -1:
            kini2, sword = in_ope_string(ksym, sval, kini)
            if kini2 != -1:
                if inside != '':
                    if inside in sword:
                        return sword
                else:
                    return sword
            kini = kini2 + len(sval)
        return sroot

    def killSimpleRoot(self):
        ss = self.findSubFunc('sqrt', '**2')
        sm = get_midle_str(ss, 'sqrt', '**2')
        exp = str(self.ksym)
        exp = exp.replace(ss, sm)
        self.ksym = parse_expr(exp)
        self.s()
    ################################################################################
    #               alone
    #################################################################################
 
    def alone(self,ksym):
        
        nname=str(ksym)
        nksym=symbols(nname,positive=True)
        nval=symbols(self.name,positive=True)
        ee=MyEq(nval-self.ksym,kshow=False)
        ee.set(ksym,nksym,kshow=False)
        kres=ee.solve(ksym,kshow=False)
        if type(kres)==list:
            vecr=[]
            for i in kres:
                vecr.append(MyEq(i,nname))
            return vecr
        else:    
            return MyEq(kres,nname)

 
################################################################################
#               END MyEq Class
#################################################################################

def get_real_value(ee):
    if type(ee) == MyEq:
        kres = ee.ksym
    else:
        kres = ee
    return kres


def show_main():
    for i in dataQ:
        i.s()



def upBagSys(ksys, Bag, kope=''):
    for i in ksys:
        i.upBag(Bag, kope=kope)



###########################################
#               END MyIteger Class
###########################################






###############################################
#  Mass center Inertia

def pQ(mm, vv, kope=''):
    rho = symbols('rho')

    kres = mm / vv
    sE([rho, '=', kres])
    return kres


#################################################
#   Solve Algorithm


def solverSys(*args, Bag=''):
    Ke = []
    Kv = []
    Kn = []

    for i in args:
        if type(i) == MyEq:
            if Bag != '':
                i.upBag(Bag, kshow=False)
            Ke.append(i)
        if type(i) == Symbol:
            Kv.append(i)
            Kn.append(i.name)
    # return(Ke,Kv,Kn)

    return MyEqSolveLin(Ke, Kv, Kn, Bag=Bag)


def MyEqSolveLin(Ke, Kv, Kn, Bag=''):  # Solve n MyEq with n unknow variable
    '''
    Example
        Ke=[e2,e2,e0]  MyEqs Matrix
        Kv=[N1,T,a]    unKnow Vriables
        Kn=['N_1','T','a_c']  New Name Eq

        N11,T1,ac = MyEqSolveLin(Ke,Kv,Kn)
        returns resepective  answer
    '''
    vecs = []
    qq = len(Ke)
    kres = []
    for i in range(qq):
        ee = Ke[i]
        ksym = Kv[i]
        ks = ee.solve(ksym,kshow=False)
        if type(ks) == list:
            rr = max(ks)
            ks = rr

        vecs.append(ks)
        Ker = Ke[i + 1::]
        for e1 in Ker:
            e1.set(ksym, ks, kshow=False)
            e1.reduFac(kshow=False)
            e1.simplify(kshow=False)

    for i, kname in zip(vecs, Kn):
        ee = MyEq(i, kname, kshow=False)
        kres.append(ee)
    ueq = kres[-1]
    ksym = ueq()
    vsym = Kv[-1]
    for ee in kres[0:-1]:
        ee.set(vsym, ksym, kshow=False)
        ee.reduFac(kshow=False)
        ee.simplify(kshow=False)
    for i in kres:
        i.s()
    return kres


def Solve2Eq(ksym=[], kvar=[], knom=[], kope=''):
    e1, e2 = ksym
    v1, v2 = kvar
    t1, t2 = knom

    r1 = e1.solve(v1)
    e2.set(v1, r1, kshow=False)
    r2 = e2.solve(v2)
    r2 = opemat(r2, kope=kope)
    e1.set(v2, r2, kshow=False)
    r1 = e1.solve(v1)
    r1 = opemat(r1, kope=kope)
    aa = MyEq(r1, t1)
    bb = MyEq(r2, t2)
    return (aa, bb)


def Diff(ksym, kvar, kname=''):
    kres = ksym
    kres = kres.diff(kvar)
    if kname == '':
        return kres
    else:
        return MyEq(kres, kname)


def Diff2(ksym, kvar, kname=''):
    kres = ksym
    kres = kres.diff(kvar)
    kres = kres.diff(kvar)
    if kname == '':
        return kres
    else:
        return MyEq(kres, kname)



def upBag2sys(vecEq, kBag):
    for i in vecEq:
        i.upBag(kBag)




def eQSolver(*args):
    vec1 = []
    uk1 = []
    for i in args:
        if type(i) == list:
            for j in i:
                if type(j) == MyEq:
                    vec1.append(j())
                elif fpoly(j, 'n') > 1:
                    vec1.append(j)
                else:
                    uk1.append(j)
        else:
            if type(i) == MyEq:
                vec1.append(i())
            elif fpoly(i, 'n') > 1:
                vec1.append(i)
            else:
                uk1.append(i)

    vec2 = []
    kres = []
    for i in vec1:
        if type(i) == MyEq:
            vec2.append(i())
        else:
            vec2.append(i)

    mm = solve(vec2, uk1)
    if type(mm) == dict:
        kk, vv = kunpakDic(mm)

        for i, j in zip(kk, vv):
            kres.append(MyEq(j, i))
        return kres
    else:
        for i, j in zip(mm[0], uk1):
            j = MyEq(i, str(j))
            kres.append(j)
        return (kres)


def solvelin(*args, kope='', Eq=True):  # solveLinearSys(e1,e2,mu1,mu2)
    mS = []
    mV = []

    for i in args:
        if type(i) == MyEq:
            mS.append(i())
        elif type(i) == eQ:
            ee = MyEq(i.ksym, kname=i.name, kshow=False)
            mS.append(ee())
        elif type(i) == str:
            kope = i
        else:
            mV.append(i)
    kres = solve(mS, mV)

    kk, vv = kunpakDic(kres)
    if kope != '':
        vv = opemat_vec(vv, kope)
    if Eq:
        EqR = []
        for i, j in zip(kk, vv):
            EqR.append(MyEq(j, i))
        return EqR

    else:
        for i, j in zip(kk, vv):
            sE([i, '=', opemat(j, kope=kope)])
    if Eq:
        kres = []
        for i, j in zip(kk, vv):
            kres.append(MyEq(opemat(j, kope=kope), i, kshow=False))
        return kres

    return vv


def get_squareMono(ksym):
    if type(ksym) == MyEq:
        ksym = ksym.ksym
    kres = ksym
    mm = fpoly(ksym, 'list')
    mr = []
    ms = []
    rr = []
    centra = 0
    ksigno = 1
    for i in mm:
        mr.append(opemat(rpow(i, 2), 'r'))
        ms.append(str(opemat(rpow(i, 2), 'r')))
    for i, j, k in zip(ms, mr, mm):
        if 'sqrt' in i:
            central = k
            if '-' in str(central):
                ksigno = -1
        else:
            rr.append(j)
    if len(rr) == 2:
        kres = kpow(rr[1] + ksigno * rr[0], 2)
    return kres


#######
def expand2MyEq(ee):
    ktype = ee.type
    var2 = ee.var2
    mm = ee.list()
    cc = 1
    kname = ee.name
    kres = []
    for i in mm:
        nname = kname + str(cc)
        nname = MyEq(i, nname, var2=var2)
        kres.append(nname)
        cc += 1
    return kres


def upgrade(*args, kshow=True, andsolve=[]):
    if andsolve != []:
        vv = andsolve[0]
        ee = andsolve[1]
        vv = ee.solve(parse_expr(vv), vv, kope='s')
    eev = []
    evv = []
    for i in args:
        if type(i) == MyEq:
            if i.type == 'C':
                eev.append(i)
            else:
                evv.append(i)

    for i in eev:
        for j in evv:
            try:
                i.upgrade(j, kshow=False)
            except:
                pass
    for i in eev:
        if i.ksym != 0:
            i.s()
    for i in evv:
        if type(i) == MyEq:
            i.simplify(kshow=False)
            if i.ksym != 0:
                i.s()


def upgradeList(*args, kshow=True, andsolve=[], kope='s'):
    eev = []
    evv = []
    for i in args:
        if type(i) == MyEq:
            if i.type == 'C':
                if i != andsolve[1]:
                    eev.append(i)

    if andsolve != []:
        vv = andsolve[0]
        ee = andsolve[1]
        kres = ee.solve(vv)
        kres = opemat(kres, 's')
        vv = MyEq(kres, str(vv), ktype='C', kshow=False)
        ee.type = 'P'

    for i in eev:
        if i.type == 'C':
            i.upgrade(vv, kshow=False)
            i.simplify(kshow=False)

    for i in eev:
        if i.ksym != 0:
            i.s()
    vv.s()
    return vv



def func_sig(kf, x1, x2, var=x):
    ee = MyEq(kf, var2=var, kshow=False)
    xx = (x2 - x1) / 2
    return ee(xx)


def get_intersec_2func(y1, y2, var=x):  # y1(x), y2(x), return intersec y1 and y2
    ee = MyEq(y1 - y2, kshow=False)  # return vector
    return ee.solve(var)


def reset_ee(*args):
    eeFull = []
    for i in args:
        i.init = False


def Upgrade(*args, kope='', kshow=True):
    newa = []
    for i in args:
        if type(i) == MyEq:
            if i.ksym != 0:
                newa.append(i)
    args = newa
    antes = []
    for i in args:
        antes.append(str(i))
    qq = len(args)
    for i in range(qq):

        mm = []
        for j in range(qq):
            if j != i:
                mm.append(args[j])
        args[i].upgrade(mm, kshow=False, kope=kope)
    if kshow:
        for i, j in zip(args, antes):
            if str(i) != newa:
                if i.ksym != 0:
                    i.s()


def presolve(ksym, val):
    kres = solve(ksym, val)
    if kres == []:
        try:
            kres = solve(opemat(ksym, 'esf'), val)
            if kres != []:
                return kres
            else:
                ksym = factorSec(ksym, val)
                kres = solve(opemat(ksym, 'esf'), val)
                if kres != []:
                    return kres
        except:
            done = False
    return kres


def eQsolve(ksym, kname, kope=''):
    kval = parse_expr(kname)
    kres = csolve(ksym, kval)
    kres = opemat(kres, kope)
    kval = MyEq(kres, kname)
    return kval
    
def Qsolve(*args):
    '''
        N1,mu1=Qsolve(FxA,FyA,N1,mu1)
    '''    
    eqq=[]
    evv=[]
    for i in args:
        if type(i)==MyEq:
            eqq.append(i)
        else:
            evv.append(i)
    kres=solve(eqq,evv)
    try:
        vsym,vval=kunpakDic(kres)
    except:
        vsym=evv
        vval=list(kres[0])
         
    vres=[]
    for i ,j in zip(vsym,vval):
        kname=str(i)
        ee=MyEq(j,kname=i)
        vres.append(ee)
    return vres   



def Diff2flat(kres,kvar,var2): # ksym,kvar,var2
    
    for i in kvar:
        f=Function(str(i))(var2)
        df=diff(f)
        kname='d'+alphaname(i)
        nf=symbols(kname)
        kres=kres.subs(df,nf)
    ee=MyEq(kres,kshow=False)
    for i in kvar:
        ee.setdiff(i,i,kshow=False)
        
    return ee.ksym
    
    
    
#####################################
#           list
#####################################  

def solvelist(*args):
    '''
    input: [vector with all eq=0], variables to find ..
    output: MyEq of each variable
    example:
        a+2*b=0 and 3*a-b=0
        ee=[a+2*b,3*a-b]
        then :
        a,b=solvelist(ee,a,b)
        return a,b in MyEq ecuation class
    '''
    vecs=args[1::]
    kres= solve(*args)
    var,value=kunpakDic(kres)
    vres=[]
    for i ,j in zip(var,value):
        ee=MyEq(j,str(i))
        vres.append(ee)
    return vres  

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
#           real subs
#####################################  
        
        
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

#####################################
#           algebra
#####################################
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
            
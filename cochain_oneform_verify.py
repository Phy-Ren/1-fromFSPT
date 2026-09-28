#!/usr/bin/env python3
"""Exact one-form cochain checks in the September 2026 zero-form convention.

Requirements: Python 3.10+ and NumPy. No manuscript files or external ancillary
bundles are imported. All phases are integer numerators modulo 8 or 16.
Run: python cochain_oneform_verify.py --output cochain_oneform_checks.json
"""
import argparse
import itertools as it
import json
from functools import lru_cache
from pathlib import Path
import numpy as np

@lru_cache(None)
def faces(n,d):
    return tuple(it.combinations(range(n+1),d+1))

@lru_cache(None)
def compositions(total,parts):
    if parts==1: return ((total,),) if total>=1 else ()
    return tuple((first,)+tail for first in range(1,total-parts+2)
                 for tail in compositions(total-first,parts-1))

@lru_cache(None)
def cuts(w,degrees):
    """Normalized interval cuts, with the signed integer convention."""
    positions=[tuple(i for i,x in enumerate(w) if x==j)
               for j in range(1,len(degrees)+1)]
    ans=[]
    for blocks in it.product(*(compositions(d+1,len(p))
                               for d,p in zip(degrees,positions))):
        lengths=[0]*len(w)
        for indices,block in zip(positions,blocks):
            for index,length in zip(indices,block): lengths[index]=length
        starts=[];v=0
        for length in lengths: starts.append(v);v+=length-1
        fs=[]
        for indices in positions:
            f=tuple(z for i in indices for z in range(starts[i],starts[i]+lengths[i]))
            if len(set(f))!=len(f):break
            fs.append(f)
        else:
            sign=(-1)**sum(lengths[i]*lengths[j] for i in range(len(w))
                            for j in range(i+1,len(w)) if w[i]>w[j])
            ans.append((tuple(fs),sign))
    return tuple(ans)

class Cochain:
    def __init__(self,n,d,values,mod=2):
        self.n,self.d,self.values,self.mod=n,d,values,mod
    def __getitem__(self,f):return self.values.get(tuple(f),0)
    def __add__(self,o):
        assert (self.n,self.d)==(o.n,o.d)
        mod=self.mod if self.mod==o.mod else None
        values={f:self[f]+o[f] for f in faces(self.n,self.d)}
        return Cochain(self.n,self.d,values,mod).reduce(mod) if mod else Cochain(self.n,self.d,values,None)
    def __sub__(self,o):return self+o.scale(-1)
    def scale(self,k):return Cochain(self.n,self.d,{f:k*v for f,v in self.values.items()},self.mod)
    def lift(self):return Cochain(self.n,self.d,self.values,None)
    def reduce(self,mod=2):return Cochain(self.n,self.d,{f:v%mod for f,v in self.values.items()},mod)
    def differential(self,mod=None):
        out=Cochain(self.n,self.d+1,{f:sum((-1)**j*self[f[:j]+f[j+1:]]
                for j in range(len(f))) for f in faces(self.n,self.d+1)},None)
        return out.reduce(mod) if mod else out
    def beta(self):
        d=self.differential()
        assert all(np.all(np.asarray(v)%2==0) for v in d.values.values())
        return Cochain(self.n,self.d+1,{f:v//2 for f,v in d.values.items()},None)
    def top(self):return self[tuple(range(self.n+1))]

def word(w,*xs,integer=False):
    w=tuple(map(int,w));degrees=tuple(x.d for x in xs)
    d=sum(degrees)-len(w)+len(xs);values={}
    for f in faces(xs[0].n,d):
        value=0
        for fs,sign in cuts(w,degrees):
            term=sign if integer else 1
            for x,indices in zip(xs,fs):term=term*x[tuple(f[k] for k in indices)]
            value=value+term
        values[f]=value if integer else value%2
    return Cochain(xs[0].n,d,values,None if integer else 2)

def cup(x,y,i=0,integer=False):
    return word(tuple(k%2+1 for k in range(i+2)),x,y,integer=integer)

def phase(*args,mod=8):
    k,c=args[0];result=c.lift().scale(k)
    for k,c in args[1:]:result=result+c.lift().scale(k)
    return result.reduce(mod)

def cyclic_background(order,dim,ids):
    free=[f for f in faces(dim,2) if f[0]==0]
    a=Cochain(dim,2,{f:(ids//order**i)%order for i,f in enumerate(free)},order)
    for f in faces(dim,2):
        if f[0]:
            i,j,k=f
            a.values[f]=(a[0,j,k]-a[0,i,k]+a[0,i,j])%order
    return a

def pontryagin(a):
    A=a.lift()
    return cup(A,A,integer=True)+cup(A,A.differential(),1,integer=True)

def majorana_obstruction(u):
    """Latest four-word Adem representative, specialized to n2=omega2=u."""
    B=u.beta();u2=cup(u,u)
    zeta=word('1231343',u,u,u,u)
    x=word('1213243',u,u,u,u)+word('1213431',u,u,u,u)+word('1232141',u,u,u,u)+word('1234321',u,u,u,u)
    return phase((4,zeta+x+cup(u2,u2,3)),
                 (2,cup(u,B,integer=True)+cup(B,B,1,integer=True)))

def majorana_pair_phase(u):
    """Latest unitary E4, including the four-word representative polarization."""
    U=u.lift();B=u.beta();u2=cup(u,u)
    beta=phase((2,cup(B,B,2,integer=True).scale(-1)
                    +cup(B,U,1,integer=True).scale(2)+cup(U,U,integer=True)))
    z=Cochain(u.n,4,{})
    for w in ('12413423','12314132','12314324','12341321','12132413','12324214'):
        z=z+word(w,u,u,u,u)
    # Delta(P(u)/8)=-P(u)/4 for two equal inputs.
    adem=phase((4,z),(2,u2),(-2,pontryagin(u)))
    omega=phase((4,cup(u2,u2,4)),(-2,cup(U,U,integer=True)))
    return (beta+adem+omega).reduce(8)

def mismatch(c,mod,top=True):
    values=[c.top()] if top else c.values.values()
    return sum(int(np.count_nonzero(np.asarray(v)%mod)) for v in values)

def check_z2():
    u=cyclic_background(2,5,np.arange(1024,dtype=np.int64)).reduce()
    b=cup(u,u,1);P=pontryagin(u);y=cup(b,u,1)
    M=phase((3,P),(8,y),mod=16)
    legacy=phase((1,P),(8,y),mod=16)
    O=majorana_obstruction(u)
    E=majorana_pair_phase(u)
    Y=Cochain(5,3,{f:u[f[0],f[1],f[3]]*(1-u[f[0],f[1],f[2]])
                         *(1-u[f[0],f[2],f[3]]) for f in faces(5,3)})
    r=Cochain(5,3,{f:(u[f[0],f[1],f[2]]*u[f[0],f[2],f[3]]
                         *(1+u[f[0],f[1],f[3]]))%2 for f in faces(5,3)})
    residuals={
      'majorana_root':mismatch(M.differential()-O.scale(2),16),
      'complex_root':mismatch(P.differential().scale(-1)-phase((4,cup(b,b,1)+cup(u,b))),8),
      'bosonic_root':mismatch(P.differential(),4),
      'cartan_primitive':mismatch(y.differential()-cup(cup(u,u),cup(u,u),3),2),
      'majorana_stacking':mismatch(E+P.scale(4)-Y.lift().differential().scale(2),8,False),
      'complex_stacking':mismatch(cup(b,b,2)-r.differential(),2,False),
    }
    return {'five_simplices':1024,'failures':residuals,
            'legacy_root_counterexamples':mismatch(legacy.differential()-O.scale(2),16)}

def check_cyclic(order,batch_size=32768):
    """Exhaust all five-simplices for Z4; bounded batches avoid large memory."""
    total=order**10;failed={k:0 for k in ('digit_primitive','gauge_primitive','bosonic_reduction')}
    for start in range(0,total,batch_size):
        a=cyclic_background(order,5,np.arange(start,min(start+batch_size,total),dtype=np.int64))
        u=a.reduce();T=Cochain(5,2,{f:v//2 for f,v in a.values.items()},None);t=T.reduce()
        b=cup(u,u,1)
        H=phase((4,cup(t,t)+cup(t,b,1)+cup(u,t)))
        K=cup(a,T,1,integer=True).scale(2)
        P=pontryagin(u);Q=pontryagin(a)
        failed['digit_primitive']+=mismatch(t.differential()-b,2,False)
        failed['gauge_primitive']+=mismatch(H.differential()-phase((4,cup(b,b,1)+cup(u,b))),8)
        # Negative complex root + H = negative P_N/8 + d(A cup1 T/4).
        failed['bosonic_reduction']+=mismatch(P.scale(-1)+H+Q-K.differential(),8,False)
    return {'five_simplices':total,'failures':failed}

def evaluate_table(u,degree,table):
    table=np.asarray(table,dtype=np.int64)
    values={}
    for f in faces(u.n,degree):
        index=0
        for bit,(i,j) in enumerate(it.combinations(range(1,degree+1),2)):
            index=index+2**bit*u[f[0],f[i],f[j]]
        values[f]=table[index]
    return Cochain(u.n,degree,values,8)

def check_microscopic_witness(path):
    data=json.loads(path.read_text())
    u=cyclic_background(2,5,np.arange(1024,dtype=np.int64)).reduce()
    B=evaluate_table(u,4,data['B4_diagonal'])
    source=np.asarray(data['microscopic_O5_diagonal'],dtype=np.int64)
    failure_obstruction=int(np.count_nonzero((source-majorana_obstruction(u).top()-B.differential().top())%8))
    u=cyclic_background(2,4,np.arange(64,dtype=np.int64)).reduce()
    B=evaluate_table(u,4,data['B4_diagonal'])
    Bzero=evaluate_table(u,4,data['B4_zero_majorana'])
    C=evaluate_table(u,3,data['C3_diagonal'])
    source=np.asarray(data['microscopic_E4_diagonal'],dtype=np.int64)
    residual=source-majorana_pair_phase(u).top()-(Bzero-B.scale(2)+C.differential()).top()
    failure_stacking=int(np.count_nonzero(residual%8))
    return {'source_archive_sha256':data['source_archive_sha256'],
            'interpretation':'Comparison with fixed microscopic input arrays, not regeneration of their projector amplitudes.',
            'single_input_states':1024,'pair_input_states':64,
            'failures':{'obstruction_coordinate_change':failure_obstruction,
                        'stacking_coordinate_change':failure_stacking}}

def check_mixed_lift_ambiguity(batch_size=32768):
    """Exhaust B²(Z4 x Z2) four-simplices for the change lambda -> lambda+2x.

    Use the allowed integer lift L+2X of the new Z4 background. The identity is
    r(lambda+2x)-r(lambda)-R_w(x)=-d(L cup_1 X/4), all modulo one.
    """
    total=4**6*2**6;bad=0
    for start in range(0,total,batch_size):
        ids=np.arange(start,min(start+batch_size,total),dtype=np.int64)
        L=cyclic_background(4,4,ids%(4**6)).lift()
        x=cyclic_background(2,4,ids//(4**6)).reduce()
        X=x.lift();w=L.reduce()
        Lprime=L+X.scale(2)
        R=phase((4,cup(x,x)+cup(w,x)))
        numerator=pontryagin(Lprime).scale(-1)+pontryagin(L)-R
        primitive=cup(L,X,1,integer=True).scale(-2)
        bad+=mismatch(numerator-primitive.differential(),8)
    return {'group':'Z4 x Z2','four_simplices':total,
            'failures':{'lift_change_mod_bosonic_identification':bad}}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result={'convention':'Current four-word x2, additive phases; B_N=-P_N/(2N).',
            'Z2':check_z2(),'Z4':check_cyclic(4),
            'mixed_lift_ambiguity':check_mixed_lift_ambiguity()}
    witness=Path(__file__).with_name('cochain_microscopic_diagonal.json')
    if not witness.is_file():
        parser.error(f'Required microscopic witness data is missing: {witness.name}')
    result['microscopic_witness']=check_microscopic_witness(witness)
    print(json.dumps(result,indent=2))
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    assert all(v==0 for block in result.values() if isinstance(block,dict)
               for v in block.get('failures',{}).values())
    assert result['Z2']['legacy_root_counterexamples']==212

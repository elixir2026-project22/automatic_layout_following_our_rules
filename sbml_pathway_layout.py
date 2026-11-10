#!/usr/bin/env python3
"""Pathway drawings from any SBML file: cycles as rings, pathways as straight layered backbones.

Rules: reproducible | labels at fixed positions per node type | no line crossings |
try to avoid long lines | try to avoid bends.

Pipeline (ideas adapted from MetaCarto, https://github.com/forxhunter/MetaCarto, CC BY 4.0 -
design reimplemented, no code copied):
 1. Primary-compound reduction: every reaction gets ONE main substrate->product pair, chosen by
    shared skeleton atoms from the formulas (fbc:chemicalFormula, notes, or InChI), with the CoA
    moiety removed from acyl-CoAs. Cofactors are a tier (never main while a non-cofactor exists);
    near-ties go to the pair that continues the network. Without formulas, names are compared.
 2. Cofactors (curated list + very highly connected compounds) are written inside the reaction box;
    other secondary participants become small duplicated side boxes (dashed) next to the reaction.
 3. Rings: node-disjoint metabolite cycles (3-12 steps, preferring cycles whose steps are carried
    by many alternative enzymes) are contracted into one unit each.
 4. The contracted network is laid out in layers (feedback-arc removal, longest-path layering,
    median ordering sweeps, least-squares straightening of chains), then each ring is expanded on a
    circle rotated towards the pathways it connects to.
 5. Any line crossing that remains is removed by giving one end of a crossing line its own copy of
    the metabolite (dashed box). Boxes are then separated locally; uniform scaling is the last
    resort. A geometric check (0 crossings, 0 lines through boxes, 0 box overlaps) gates the output.
 6. Models with more than --max-cluster reactions are split by subsystem / pathway groups (or by
    network communities when the file has none), one drawing per part.
Nothing is random; the hash seed is pinned and SVGs are written without a date.

Usage: python3 sbml_pathway_layout.py model.xml[.gz] [--out NAME] [--max-cluster 60]
                                       [--side "x,y"] [--keep "x,y"]
Requires: networkx>=3.1, numpy, matplotlib
"""
import os, sys
if os.environ.get('PYTHONHASHSEED')!='0':
    os.environ['PYTHONHASHSEED']='0'; os.execv(sys.executable,[sys.executable]+sys.argv)
import re, math, argparse, itertools, gzip
import xml.etree.ElementTree as ET
import numpy as np, networkx as nx
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
plt.rcParams['svg.hashsalt']='sbml-pathway-layout'; plt.rcParams['font.family']='DejaVu Sans'
def log(*a): print(*a,file=sys.stderr,flush=True)

# ============================================================ names, formulas
def _norm(n):
    n=n.lower(); n=re.sub(r'\[[^\]]*\]$','',n); n=re.sub(r'_[a-z]$','',n)
    n=n.replace('(+)','+').replace('(-)','-').replace('(2+)','2+').replace('(3+)','3+')
    return re.sub(r"[\s_\-,'()]",'',n)
SIDE={_norm(x) for x in r'''h2o water oh- hydroxide h+ proton hydron atp adp amp gtp gdp gmp utp udp ump ctp cdp cmp itp idp
nad nad+ nadh nadp nadp+ nadph fad fadh2 fmn fmnh2 coa coenzymea coash co2 carbondioxide hco3- bicarbonate hydrogencarbonate
phosphate orthophosphate pi hydrogenphosphate diphosphate pyrophosphate ppi o2 oxygen dioxygen h2o2 hydrogenperoxide
nh3 nh4+ ammonia ammonium aquinone aquinol quinone quinol ubiquinone ubiquinol menaquinone menaquinol
oxidizedferredoxin reducedferredoxin ferredoxin oxidizedthioredoxin reducedthioredoxin thioredoxin acceptor reducedacceptor
electron e- na+ k+ cl- fe2+ fe3+ mg2+ h2 hydrogen nicotinamideadeninedinucleotide nicotinamideadeninedinucleotidereduced
nicotinamideadeninedinucleotidephosphate nicotinamideadeninedinucleotidephosphatereduced flavinadeninedinucleotide
flavinadeninedinucleotideoxidized flavinadeninedinucleotidereduced adenosinetriphosphate adenosinediphosphate
adenosinemonophosphate atp(4-) adp(3-) adenosine5'-triphosphate adenosine5'-diphosphate adenosine5'-monophosphate
guanosinetriphosphate guanosinediphosphate coenzyme-a ubiquinone-8 ubiquinol-8 menaquinone-8 menaquinol-8'''.split()}
SIDE_RE=re.compile(r'^((ubi|mena|plasto)?quino(ne|l)\d*|(oxidized|reduced)(ferredoxin|thioredoxin|flavodoxin|cytochromec?)|cytochromec(ii+)?)$')
def curated_side(name,extra,keep):
    n=_norm(name)
    return n not in keep and (n in extra or n in SIDE or bool(SIDE_RE.match(n)))

def parse_formula(f):
    out={}
    for part in f.split('.'):
        m=re.match(r'^(\d*)(.*)$',part); mult=int(m.group(1) or 1)
        for el,n in re.findall(r'([A-Z][a-z]?)(\d*)',m.group(2)): out[el]=out.get(el,0)+mult*(int(n) if n else 1)
    return out or None
COA={'C':21,'H':32,'N':7,'O':16,'P':3,'S':1}
def skeleton(name,f):
    """Acyl-CoA -> acyl group, so succinyl-CoA pairs with succinate rather than with acetyl-CoA."""
    n=_norm(name)
    if f and 'coa' in n and n not in ('coa','coenzymea') and all(f.get(k,0)>=v for k,v in COA.items()):
        return {k:f.get(k,0)-COA.get(k,0) for k in f}
    return f
def _stem(n):
    n=re.sub(r'[^a-z]','',_norm(n).split('[')[0])
    for w in ('bisphosphate','diphosphate','triphosphate','phosphate','phospho'): n=n.replace(w,'')
    n=re.sub(r'^(d|l|s|r)(?=[a-z]{4})','',n)
    for suf in ('coa','yl','ate','icacid','ic','ine'):
        if n.endswith(suf) and len(n)-len(suf)>=3: n=n[:-len(suf)]
    return n
def _lcs(a,b):
    best=0
    for i in range(len(a)):
        for j in range(len(b)):
            k=0
            while i+k<len(a) and j+k<len(b) and a[i+k]==b[j+k]: k+=1
            best=max(best,k)
    return best
W={'C':1.0,'N':0.3,'P':0.1,'S':0.1}
def similarity(a,b,F):
    fa,fb=skeleton(a,F.get(a)),skeleton(b,F.get(b))
    if fa and fb:
        num=sum(w*min(fa.get(k,0),fb.get(k,0)) for k,w in W.items()); den=sum(w*max(fa.get(k,0),fb.get(k,0)) for k,w in W.items())
        if den>0: return num/den
    sa,sb=_stem(a),_stem(b)
    return _lcs(sa,sb)/max(len(sa),len(sb),1) if sa and sb else 0.0

# ============================================================ SBML (levels 1-3)
def _l(tag): return tag.rsplit('}',1)[-1]
def _kids(e,name): return [c for c in e if _l(c.tag)==name]
def _attr(e,name):
    for k,v in e.attrib.items():
        if _l(k)==name: return v
    return None
def read_sbml(path):
    root=(ET.parse(gzip.open(path)) if path.endswith('.gz') else ET.parse(path)).getroot()
    model=next(c for c in root if _l(c.tag)=='model')
    species={}
    for e in model.iter():
        if _l(e.tag) not in ('species','specie'): continue
        sid=e.get('id') or e.get('name'); f=_attr(e,'chemicalFormula')
        if not f:
            txt=' '.join(''.join(c.itertext()) for c in e if _l(c.tag) in ('notes','annotation'))
            m=re.search(r'FORMULA:\s*([A-Za-z0-9.]+)',txt) or re.search(r'InChI=1S?/([^/\s]+)',txt)
            f=m.group(1) if m else None
        species[sid]=dict(name=(e.get('name') or sid).strip(),comp=e.get('compartment') or '',formula=parse_formula(f) if f else None)
    byname={}
    for sid,s in species.items(): byname.setdefault(s['name'],[]).append(sid)
    for nm,ids in byname.items():
        if len({species[i]['comp'] for i in ids})>1:
            for i in ids: species[i]['name']=f"{nm} [{species[i]['comp']}]"
    F={s['name']:s['formula'] for s in species.values()}
    groups=[]
    for e in model.iter():
        if _l(e.tag)=='group':
            nm=_attr(e,'name') or _attr(e,'id'); kind=_attr(e,'kind') or ''
            mem=[_attr(m,'idRef') for m in e.iter() if _l(m.tag)=='member']
            groups.append((nm,kind,[x for x in mem if x]))
    mname=(model.get('name') or '').strip().lower()
    rg={}
    for nm,kind,mem in groups:
        for x in mem: rg.setdefault(x,[]).append((nm,kind))
    reactions=[]
    for e in model.iter():
        if _l(e.tag)!='reaction': continue
        def refs(lst):
            out=[]
            for L in _kids(e,lst):
                for r in L:
                    if _l(r.tag) in ('speciesReference','specieReference'):
                        sid=r.get('species') or r.get('specie'); st=r.get('stoichiometry')
                        st=float(st) if st not in (None,'') else 1.0
                        if sid in species: out.append((species[sid]['name'],st))
            return out
        rid=e.get('id') or e.get('name'); nm=(e.get('name') or rid).strip()
        subs,prods=refs('listOfReactants'),refs('listOfProducts')
        mine=[(n,k) for n,k in rg.get(rid,[]) if n and n.strip().lower()!=mname]
        if re.fullmatch(r'[A-Za-z_]{0,6}\d+[A-Za-z0-9_]*',nm) or nm==rid:      # name is just an identifier
            toks={t for p,_ in subs+prods for t in re.findall(r'[a-z0-9+]{3,}',p.lower())}
            ec=[n for n,k in mine if re.match(r'EC [0-9]',n) or k=='classification']
            if ec:
                best=min(ec,key=lambda x:(-len(toks&set(re.findall(r'[a-z0-9+]{3,}',x.lower()))),len(x),x))
                nm=re.sub(r'^EC [0-9.\-nB]+\s+','',best)
        notes=' '.join(''.join(c.itertext()) for c in e if _l(c.tag)=='notes')
        m=re.search(r'SUBSYSTEM:\s*([^\n<]+)',notes)
        sub=m.group(1).strip() if m else None
        if not sub:
            parts=sorted(n for n,k in mine if not re.match(r'EC [0-9]',n) and k!='classification')
            sub=parts[0] if parts else None
        reactions.append(dict(id=rid,name=nm,subs=subs,prods=prods,reversible=(e.get('reversible','true')!='false'),subsystem=sub))
    return (model.get('name') or model.get('id') or 'model'),F,reactions

# ============================================================ primary-compound reduction
def reduce(reactions,F,extra,keep):
    deg={}
    for r in reactions:
        for n,_ in r['subs']+r['prods']: deg[n]=deg.get(n,0)+1
    thr=max(10.0,float(np.percentile(list(deg.values()),97)) if deg else 10.0,0.12*len(reactions))
    cof={m for m in deg if curated_side(m,extra,keep) or (deg[m]>thr and _norm(m) not in keep)}
    usage={}
    for r in reactions:
        for n,_ in r['subs']+r['prods']:
            if n not in cof: usage[n]=usage.get(n,0)+1
    for r in reactions:
        S=[n for n,_ in r['subs'] if n not in cof]; P=[n for n,_ in r['prods'] if n not in cof]
        r['main']=None
        if S and P:
            cand=[(similarity(s,p,F),s,p) for s in S for p in P if s!=p]
            if cand:
                top=max(c[0] for c in cand)
                near=[c for c in cand if c[0]>=top-0.05]
                c=max(near,key=lambda c:(min(usage[c[1]],usage[c[2]]),usage[c[1]]+usage[c[2]],c[0],c[1],c[2]))
                r['main']=(c[1],c[2])
        fmt=lambda lst:' + '.join((f'{s:g} ' if s!=1 else '')+n for n,s in lst)
        cs=[x for x in r['subs'] if x[0] in cof]; cp=[x for x in r['prods'] if x[0] in cof]
        r['cof_text']=(fmt(cs)+' \u2192 '+fmt(cp)) if (cs or cp) else None
        used=set(r['main'] or ())
        r['secondary']=[(n,'s') for n in S if n not in used]+[(n,'p') for n in P if n not in used]
    return cof

# ============================================================ rings
def _canon(c):
    i=c.index(min(c)); c=list(c[i:])+list(c[:i]); r=[c[0]]+c[1:][::-1]
    return min(c,r)
def find_rings(M,maxlen=12):
    rings=[]; avail=set(M)
    while True:
        H=M.subgraph(avail); cands={}
        for c in itertools.islice(nx.simple_cycles(H,length_bound=min(maxlen,max(3,len(H)))),200000):
            if len(c)<3: continue
            k=tuple(_canon(c))
            if k not in cands: cands[k]=[H.edges[a,b]['n'] for a,b in zip(k,k[1:]+k[:1])]
        if not cands: break
        mx=max(len(k) for k in cands)
        def key(k):
            gm=math.exp(sum(math.log(x) for x in cands[k])/len(k))
            return (round(gm,6),len(k),sum(M.degree(x) for x in k),k)
        best=max((k for k in cands if len(k)>=max(3,0.7*mx)),key=key)
        rings.append(list(best)); avail-=set(best)
    return rings

# ============================================================ drawing model
STYLE={'core': dict(fs=9.5,weight='bold',style='normal',pad=(7,5),fc='#fde3e1',ec='#d9534f',tc='#3a0d0b'),
       'met':  dict(fs=8.5,weight='normal',style='normal',pad=(6,4),fc='#e3ecfd',ec='#5b8def',tc='#0f2550'),
       'alias':dict(fs=7.5,weight='normal',style='normal',pad=(5,3),fc='#eef3fd',ec='#8aa9ec',tc='#33456b',ls='--'),
       'rxn':  dict(fs=7.0,weight='normal',style='italic',pad=(4,3),fc='#ffffff',ec='#2b2b2b',tc='#222222'),
       'dummy':None}
SUB=dict(fs=5.8,weight='normal',style='normal',pad=(4,0),tc='#777777')
_tw={}
def text_w(t,st):
    k=(t,st['fs'],st['weight'],st['style'])
    if k not in _tw: _tw[k]=TextPath((0,0),t,prop=FontProperties(family='DejaVu Sans',size=st['fs'],weight=st['weight'],style=st['style'])).get_extents().width
    return _tw[k]

class Drawing:
    def __init__(s):
        s.info={}; s.P={}; s.links=[]; s.rings={}
    def add(s,n,kind,label='',sub=None):
        s.info[n]=dict(kind=kind,label=label,sub=sub)
    def size(s,n):
        i=s.info[n]
        if i['kind']=='dummy': return (0.0,0.0)
        st=STYLE[i['kind']]; lab=i['label'] if len(i['label'])<=34 else i['label'][:32]+'\u2026'
        w=text_w(lab,st)+2*st['pad'][0]; h=st['fs']+2*st['pad'][1]
        if i['sub']:
            sub=i['sub'] if len(i['sub'])<=60 else i['sub'][:58]+'\u2026'
            w=max(w,text_w(sub,SUB)+2*SUB['pad'][0]); h+=SUB['fs']*1.25
        return (w,h)
    def box(s,n,P=None):
        x,y=(P or s.P)[n]; w,h=s.size(n); return (x-w/2,y-h/2,x+w/2,y+h/2)
    def link(s,a,b,arrow=True,ring=None,own=None,chain=None,cut=True):
        s.links.append(dict(a=a,b=b,arrow=arrow,ring=ring,own=own or {a,b},chain=chain,cut=cut))
    def path(s,L):
        a,b=s.P[L['a']],s.P[L['b']]
        if L['ring'] is None: return [a,b]
        (cx,cy),r=s.rings[L['ring']]
        t1,t2=math.atan2(a[1]-cy,a[0]-cx),math.atan2(b[1]-cy,b[0]-cx); d=(t2-t1+math.pi)%(2*math.pi)-math.pi
        n=max(3,int(abs(d)*30))
        return [(cx+r*math.cos(t1+d*i/n),cy+r*math.sin(t1+d*i/n)) for i in range(n+1)]
    def real(s): return sorted(n for n in s.P if s.info[n]['kind']!='dummy')

# ---------------------------------------------------------------- geometry checks
def _segx(A,B,C,D):
    o=lambda p,q,r:(q[...,0]-p[...,0])*(r[...,1]-p[...,1])-(q[...,1]-p[...,1])*(r[...,0]-p[...,0])
    return (o(A,B,C)*o(A,B,D)<-1e-9)&(o(C,D,A)*o(C,D,B)<-1e-9)
def segments(D,idx=None):
    S=[]
    for i,L in enumerate(D.links):
        if idx is not None and i not in idx: continue
        p=D.path(L)
        for a,b in zip(p,p[1:]): S.append((i,a,b))
    return S
def crossing_pairs(D,only=None):
    S=segments(D)
    if not S: return []
    I=np.array([s[0] for s in S]); A=np.array([s[1] for s in S]); B=np.array([s[2] for s in S])
    ends=[{D.links[i]['a'],D.links[i]['b']} for i in range(len(D.links))]
    out=set()
    for k,(i,a,b) in enumerate(S):
        if only is not None and i not in only: continue
        hit=np.where(_segx(np.array(a)[None],np.array(b)[None],A,B))[0]
        for h in hit:
            j=I[h]
            if j!=i and not (ends[i]&ends[j]): out.add((min(i,j),max(i,j)))
    return sorted(out)
def violations(D,margin=3):
    S=segments(D); names=D.real(); Bx=np.array([D.box(n) for n in names]); thr=0; ov=0
    for i,a,b in S:
        own=D.links[i]['own']; ts=np.linspace(0,1,10)
        X=a[0]+(b[0]-a[0])*ts; Y=a[1]+(b[1]-a[1])*ts
        hit=(X[None]>Bx[:,0,None]-margin)&(X[None]<Bx[:,2,None]+margin)&(Y[None]>Bx[:,1,None]-margin)&(Y[None]<Bx[:,3,None]+margin)
        thr+=sum(1 for k in np.where(hit.any(1))[0] if names[k] not in own)
    o=(Bx[:,None,0]<Bx[None,:,2]+margin)&(Bx[None,:,0]<Bx[:,None,2]+margin)&(Bx[:,None,1]<Bx[None,:,3]+margin)&(Bx[None,:,1]<Bx[:,None,3]+margin)
    ov=int(np.triu(o,1).sum())
    return {'crossings':len(crossing_pairs(D)),'line_through_box':thr,'box_overlaps':ov}

# ============================================================ layered placement (Sugiyama)
def greedy_fas_order(nodes,edges):
    """Eades-Lin-Smyth: an order with few backward edges."""
    out={n:set() for n in nodes}; inn={n:set() for n in nodes}
    for u,v in edges:
        if u!=v: out[u].add(v); inn[v].add(u)
    alive=set(nodes); s1=[]; s2=[]
    while alive:
        changed=True
        while changed:
            changed=False
            for n in sorted(alive):
                if not (out[n]&alive): s2.insert(0,n); alive.discard(n); changed=True
            for n in sorted(alive):
                if not (inn[n]&alive): s1.append(n); alive.discard(n); changed=True
        if alive:
            n=max(sorted(alive),key=lambda n:len(out[n]&alive)-len(inn[n]&alive)); s1.append(n); alive.discard(n)
    return s1+s2
def pava(y,w=None):
    """Least-squares non-decreasing fit (pool adjacent violators)."""
    blocks=[]
    for i,v in enumerate(y):
        blocks.append([v,1.0,[i]])
        while len(blocks)>1 and blocks[-2][0]>blocks[-1][0]:
            v2,w2,i2=blocks.pop(); v1,w1,i1=blocks.pop(); blocks.append([(v1*w1+v2*w2)/(w1+w2),w1+w2,i1+i2])
    out=[0.0]*len(y)
    for v,_,ix in blocks:
        for i in ix: out[i]=v
    return out
def sugiyama(nodes,edges,size,gap_x=34,gap_y=34):
    """nodes: list; edges: list of (u,v); size: n->(w,h). Returns positions, dummy chains."""
    order=greedy_fas_order(nodes,edges); rank={n:i for i,n in enumerate(order)}
    E=[]; rev=[]
    for u,v in edges:
        if rank[u]<rank[v]: E.append((u,v)); rev.append(False)
        else: E.append((v,u)); rev.append(True)
    G=nx.DiGraph(); G.add_nodes_from(nodes); G.add_edges_from(E)
    layer={}
    for n in nx.topological_sort(G):
        layer[n]=max([layer[p]+1 for p in G.predecessors(n)],default=0)
    for n in reversed(list(nx.topological_sort(G))):          # pull sources down next to their successors
        if G.in_degree(n)==0 and G.out_degree(n)>0: layer[n]=min(layer[s] for s in G.successors(n))-1
    chains=[]; dsize={}; adj=nx.DiGraph(); adj.add_nodes_from(nodes)
    for k,(u,v) in enumerate(E):
        pts=[u]
        for l in range(layer[u]+1,layer[v]):
            d=f'__d{k}_{l}'; layer[d]=l; dsize[d]=(8.0,0.0); pts.append(d)
        pts.append(v); chains.append((pts,rev[k]))
        for a,b in zip(pts,pts[1:]): adj.add_edge(a,b)
    sz=lambda n: dsize.get(n) or size[n]
    L=max(layer.values())+1 if layer else 0; lmin=min(layer.values()) if layer else 0
    layers=[[] for _ in range(L-lmin)]
    seen=set()
    for n in sorted(nodes,key=lambda n:(layer[n],n)):
        for m in nx.dfs_preorder_nodes(adj.to_undirected(as_view=True),n):
            if m not in seen: seen.add(m); layers[layer[m]-lmin].append(m)
    def crossings(ls):
        c=0
        for i in range(len(ls)-1):
            pos={n:j for j,n in enumerate(ls[i+1])}; es=[]
            for j,n in enumerate(ls[i]):
                for m in adj.successors(n):
                    if m in pos: es.append((j,pos[m]))
            for (a,b),(c2,d) in itertools.combinations(es,2):
                if (a-c2)*(b-d)<0: c+=1
        return c
    best=([list(l) for l in layers],crossings(layers))
    for it in range(24):
        down=it%2==0; rng=range(1,len(layers)) if down else range(len(layers)-2,-1,-1)
        for i in rng:
            ref={n:j for j,n in enumerate(layers[i-1 if down else i+1])}
            def key(n,j):
                nb=[ref[m] for m in (adj.predecessors(n) if down else adj.successors(n)) if m in ref]
                return (float(np.median(nb)) if nb else float(j),j)
            layers[i]=[n for j,n in sorted(enumerate(layers[i]),key=lambda t:key(t[1],t[0]))]
        c=crossings(layers)
        if c<best[1]: best=([list(l) for l in layers],c)
    layers=best[0]
    x={}
    for l in layers:
        cur=0.0
        for n in l: w=sz(n)[0]; x[n]=cur+w/2; cur+=w+gap_x
    und=adj.to_undirected(as_view=True)
    for it in range(40):                                        # straighten: median of neighbours, least squares
        for l in (layers if it%2==0 else layers[::-1]):
            if not l: continue
            want=[]
            for n in l:
                nb=[x[m] for m in und[n]]
                want.append(float(np.median(nb)) if nb else x[n])
            off=[0.0]
            for a,b in zip(l,l[1:]): off.append(off[-1]+(sz(a)[0]+sz(b)[0])/2+gap_x)
            fit=pava([w-o for w,o in zip(want,off)])
            for n,f,o in zip(l,fit,off): x[n]=f+o
    y={}; cur=0.0
    for l in layers:
        h=max([sz(n)[1] for n in l],default=0)
        for n in l: y[n]=-(cur+h/2)
        cur+=h+gap_y
    pos={n:(x[n],y[n]) for n in x}
    return pos,chains

# ============================================================ cluster layout
def layout_cluster(rxns):
    D=Drawing()
    M=nx.Graph()
    for r in rxns:
        if r['main']:
            a,b=r['main']
            if M.has_edge(a,b): M.edges[a,b]['n']+=1
            else: M.add_edge(a,b,n=1)
    rings=find_rings(M) if M.number_of_edges() else []
    ringof={m:i for i,ring in enumerate(rings) for m in ring}
    for r in rxns:
        R='R:'+r['id']; D.add(R,'rxn',r['name'],r['cof_text'])
        for n in (r['main'] or ()): D.add('M:'+n,'core' if n in ringof else 'met',n)
    # ---- ring geometry: members and the reactions of each step
    ring_steps={}; internal=set()
    for i,ring in enumerate(rings):
        pos={m:k for k,m in enumerate(ring)}; steps={}; chords=[]
        for r in rxns:
            if r['main'] and all(x in pos for x in r['main']) and ringof[r['main'][0]]==i:
                a,b=pos[r['main'][0]],pos[r['main'][1]]; internal.add(r['id'])
                if (b-a)%len(ring)==1: steps.setdefault(a,[]).append((r,True))
                elif (a-b)%len(ring)==1: steps.setdefault(b,[]).append((r,False))
                else: chords.append(r)
        fw=sum(len([1 for _,f in v if f])-len([1 for _,f in v if not f]) for v in steps.values())
        if fw<0:                                                  # orient the ring with the reactions
            ring.reverse(); rings[i]=ring
            return_steps=None
        pos={m:k for k,m in enumerate(ring)}; steps={}
        for r in rxns:
            if r['id'] in internal and all(x in pos for x in r['main']):
                a,b=pos[r['main'][0]],pos[r['main'][1]]
                if (b-a)%len(ring)==1: steps.setdefault(a,[]).append(r)
                elif (a-b)%len(ring)==1: steps.setdefault(b,[]).append(r)
        for k in steps: steps[k].sort(key=lambda r:(-1 if r['main'][0]==ring[k] else 0,len(r['name']),r['id']))
        items=[]                                                  # boxes around the circle, in order
        for k,m in enumerate(ring):
            items.append(D.size('M:'+m))
            items.append(D.size('R:'+steps[k][0]['id']) if k in steps else (10.0,10.0))
        def need(rad):                                            # arc each box needs vs. arc available
            tot=0.0
            for j,(w,h) in enumerate(items):
                t=2*math.pi*j/len(items); tot+=abs(w*math.sin(t))+abs(h*math.cos(t))+22
            return tot/(2*math.pi)
        rad=70.0
        for _ in range(4): rad=max(70.0,need(rad))
        mw=max(max(D.size('M:'+m)[0] for m in ring),max([D.size('R:'+r['id'])[0] for k in steps for r in steps[k]],default=0))
        ring_steps[i]=dict(ring=ring,steps=steps,chords=chords,rad=rad,box=(2*rad*1.5+mw,2*rad*1.5+40))
    # ---- contracted graph for the layered pass
    unit=lambda n:f'RING{ringof[n]}' if n in ringof else 'M:'+n
    nodes=set(); edges=[]; emeta=[]
    for r in rxns:
        if r['main'] and r['id'] not in internal:
            R='R:'+r['id']; a,b=unit(r['main'][0]),unit(r['main'][1]); nodes|={a,b,R}
            edges.append((a,R)); emeta.append(('M:'+r['main'][0],R)); edges.append((R,b)); emeta.append((R,'M:'+r['main'][1]))
    for i in ring_steps: nodes.add(f'RING{i}')
    size={}
    for n in nodes:
        size[n]=ring_steps[int(n[4:])]['box'] if n.startswith('RING') else D.size(n)
    comps=sorted((sorted(c) for c in nx.connected_components(nx.Graph([e for e in edges]+[(n,n) for n in nodes]))),key=lambda c:(-len(c),c))
    xoff=0.0
    for comp in comps:
        cs=set(comp); ce=[(e,m) for e,m in zip(edges,emeta) if e[0] in cs]
        pos,chains=sugiyama(comp,[e for e,_ in ce],size)
        xs=[pos[n][0]-size.get(n,(8,0))[0]/2 for n in pos]; xe=[pos[n][0]+size.get(n,(8,0))[0]/2 for n in pos]
        dx=xoff-min(xs); xoff+=max(xe)-min(xs)+80
        for n,p in pos.items():
            if n.startswith('__d'): D.add(n,'dummy')
            if not n.startswith('RING'): D.P[n]=(p[0]+dx,p[1])
            else: D.rings[int(n[4:])]=((p[0]+dx,p[1]),ring_steps[int(n[4:])]['rad'])
        for (pts,rv),(e,m) in zip(chains,ce):
            pts=list(pts)
            if rv: pts=pts[::-1]
            real=[m[0] if pts[0].startswith('RING') else pts[0], m[1] if pts[-1].startswith('RING') else pts[-1]]
            pts[0],pts[-1]=real
            cid=len(D.links)
            for a,b in zip(pts,pts[1:]): D.link(a,b,arrow=(b==pts[-1]),own=set(real),chain=cid)
    # ---- expand rings, rotated towards their outside connections
    for i,rs in ring_steps.items():
        (cx,cy),rad=D.rings[i]; ring=rs['ring']; n=len(ring)
        ext=[]
        for L in D.links:
            for end,other in ((L['a'],L['b']),(L['b'],L['a'])):
                if end.startswith('M:') and end[2:] in ring and ringof.get(end[2:])==i and other in D.P: ext.append((end[2:],D.P[other]))
        best=None
        for mir in (1,-1):
            for k in range(72):
                th=2*math.pi*k/72
                p={m:(cx+rad*math.cos(th+mir*2*math.pi*j/n),cy+rad*math.sin(th+mir*2*math.pi*j/n)) for j,m in enumerate(ring)}
                sc=sum(math.dist(p[m],q) for m,q in ext)
                if best is None or sc<best[0]-1e-6: best=(sc,p,th,mir)
        _,p,th,mir=best
        for m in ring: D.P['M:'+m]=p[m]
        for k,lst in rs['steps'].items():
            t=th+mir*2*math.pi*(k+0.5)/n; ux,uy=math.cos(t),math.sin(t)
            a,b='M:'+ring[k],'M:'+ring[(k+1)%n]
            for j,r in enumerate(lst):
                R='R:'+r['id']; s,pp='M:'+r['main'][0],'M:'+r['main'][1]
                if j==0:
                    D.P[R]=(cx+rad*ux,cy+rad*uy); D.link(s,R,arrow=True,ring=i); D.link(R,pp,arrow=True,ring=i)
                else:
                    f=1+(0.24 if j%2 else -0.24)*((j+1)//2); D.P[R]=(cx+rad*f*ux,cy+rad*f*uy)
                    D.link(s,R,arrow=True); D.link(R,pp,arrow=True)
        for r in rs['chords']:
            R='R:'+r['id']; s,pp='M:'+r['main'][0],'M:'+r['main'][1]
            D.P[R]=((D.P[s][0]+D.P[pp][0])/2*0.82+cx*0.18,(D.P[s][1]+D.P[pp][1])/2*0.82+cy*0.18)
            D.link(s,R); D.link(R,pp)
    # ---- pieces attached only to one ring: place them radially outside that ring
    Gc=nx.Graph(); Gc.add_nodes_from(nodes); Gc.add_edges_from(edges)
    att={}                                                        # (piece node, ring) -> ring members
    for (a,b),(ma,mb) in zip(edges,emeta):
        if a.startswith('RING'): att.setdefault(b,set()).add(ma)
        if b.startswith('RING'): att.setdefault(a,set()).add(mb)
    for i in sorted(ring_steps):
        X=f'RING{i}'
        if X not in Gc: continue
        H=Gc.copy(); H.remove_node(X)
        for comp in sorted((sorted(c) for c in nx.connected_components(H)),key=lambda c:c):
            cs=set(comp); outside={m for n in cs for m in Gc[n] if m not in cs}
            if outside!={X} or any(n.startswith('RING') for n in cs): continue
            # drop the layered drawing of this piece
            keep=[]
            for L in D.links:
                own=L['own']
                if L['ring'] is None and own & cs and own <= cs|{'M:'+m for m in rings[i]}: 
                    for e in (L['a'],L['b']):
                        if D.info.get(e,{}).get('kind')=='dummy': D.P.pop(e,None)
                    continue
                keep.append(L)
            D.links=keep
            (cx,cy),rad=D.rings[i]
            ang={m:math.atan2(D.P['M:'+m][1]-cy,D.P['M:'+m][0]-cx) for m in rings[i]}
            PG=nx.Graph(); PG.add_edges_from((a,b) for a,b in edges if a in cs or b in cs)
            for n in cs:
                for m in att.get(n,()): PG.add_edge(n,m)
            dist={m:nx.single_source_shortest_path_length(PG,'M:'+m) for m in rings[i] if 'M:'+m in PG}
            for n in comp:
                ds=[(dist[m][n],m) for m in dist if n in dist[m]] or [(1,rings[i][0])]
                sx=sum(math.cos(ang[m])/(d*d) for d,m in ds); sy=sum(math.sin(ang[m])/(d*d) for d,m in ds)
                t=math.atan2(sy,sx); dmin=min(d for d,_ in ds)
                r=rad+55+75*(dmin-1)+(12 if D.info[n]['kind']=='rxn' else 0)
                D.P[n]=(cx+r*math.cos(t),cy+r*math.sin(t))
            for r_ in rxns:
                if r_['main'] and 'R:'+r_['id'] in cs:
                    R='R:'+r_['id']; D.link('M:'+r_['main'][0],R); D.link(R,'M:'+r_['main'][1])
    # ---- reactions without a main pair, and secondary participants: side boxes
    pending=[]
    k=0
    for r in rxns:
        R='R:'+r['id']
        if not r['main']:
            sec=sorted({n for n,_ in r['secondary']})
            anchor=max(sec,key=lambda n:(sum(1 for x in rxns if x['main'] and n in x['main']),n)) if sec else None
            if anchor and 'M:'+anchor in D.P:
                side=[sd for n,sd in r['secondary'] if n==anchor][0]
                pending.append((R,'M:'+anchor,side=='s'))
                r=dict(r,secondary=[(n,sd) for n,sd in r['secondary'] if n!=anchor])
            elif anchor:
                D.add('M:'+anchor,'met',anchor); D.P['M:'+anchor]=(xoff+60*k,200.0); k+=1
                side=[sd for n,sd in r['secondary'] if n==anchor][0]
                pending.append((R,'M:'+anchor,side=='s'))
                r=dict(r,secondary=[(n,sd) for n,sd in r['secondary'] if n!=anchor])
            else:
                D.P[R]=(xoff+150*k,300.0); k+=1
        for j,(n,sd) in enumerate(r['secondary']):
            a=f'A:{r["id"]}:{j}'; D.add(a,'alias',n); pending.append((a,R,sd=='p'))
    order=sorted(pending,key=lambda t:(0 if t[0].startswith('R:') else 1,t[0]))
    for node,anchor,into_node in order:
        plant(D,node,anchor,into_node)
    return D,rings

def plant(D,node,anchor,into_node,dist0=None):
    """Put a side box in the free angular gap of its anchor (crossing-free, fewest box conflicts)."""
    x,y=D.P[anchor]
    dirs=[]
    for L in D.links:
        if anchor in (L['a'],L['b']):
            p=D.path(L); q=p[1] if L['a']==anchor else p[-2]
            dirs.append(math.atan2(q[1]-y,q[0]-x))
    dirs=sorted(dirs)
    if not dirs: cands=[(2*math.pi,-math.pi/2-math.pi)]
    else: cands=sorted(((((dirs[(i+1)%len(dirs)]-dirs[i])%(2*math.pi)) or 2*math.pi,dirs[i]) for i in range(len(dirs))),key=lambda t:(-t[0],t[1]))
    w,h=D.size(node); aw,ah=D.size(anchor); base=dist0 or (max(w,aw)/2+max(h,ah)+22)
    D.link(anchor,node,arrow=into_node) if into_node else D.link(node,anchor,arrow=True)
    li=len(D.links)-1; best=None
    for gap,lo in cands[:3]:
        for frac in (0.5,0.35,0.65,0.2,0.8):
            t=lo+gap*frac
            for f in (1.0,1.5,2.2,0.7):
                D.P[node]=(x+base*f*math.cos(t),y+base*f*math.sin(t))
                if crossing_pairs(D,only={li}): continue
                v=violations_local(D,node,li)
                key=(v,abs(f-1.0),-gap,frac)
                if best is None or key<best[0]: best=(key,D.P[node])
                if v==0: break
            if best and best[0][0]==0: break
        if best and best[0][0]==0: break
    D.P[node]=best[1] if best else (x+4,y+4)
def violations_local(D,node,li,margin=3):
    names=D.real(); Bx={n:D.box(n) for n in names}; b=Bx[node]; c=0
    for n in names:
        if n==node: continue
        o=Bx[n]
        if b[0]<o[2]+margin and o[0]<b[2]+margin and b[1]<o[3]+margin and o[1]<b[3]+margin: c+=1
    for i,a,bb in segments(D):
        if node in D.links[i]['own']: continue
        for t in np.linspace(0,1,10):
            px,py=a[0]+(bb[0]-a[0])*t,a[1]+(bb[1]-a[1])*t
            if b[0]-margin<px<b[2]+margin and b[1]-margin<py<b[3]+margin: c+=1; break
    for i,a,bb in segments(D,{li}):
        for n in names:
            if n in D.links[li]['own']: continue
            o=Bx[n]
            for t in np.linspace(0,1,10):
                px,py=a[0]+(bb[0]-a[0])*t,a[1]+(bb[1]-a[1])*t
                if o[0]-margin<px<o[2]+margin and o[1]-margin<py<o[3]+margin: c+=1; break
    return c

# ============================================================ enforcing the rules
def remove_crossings(D):
    """Cut a crossing line at its metabolite end and give that end its own copy (dashed box)."""
    n_cut=0
    for _ in range(500):
        cp=crossing_pairs(D)
        if not cp: return n_cut
        i,j=cp[0]
        def cand(k):
            L=D.links[k]; own=sorted(L['own'])
            mets=[x for x in own if x.startswith('M:')]; rx=[x for x in own if x.startswith('R:')]
            if not mets or not rx or L['ring'] is not None: return None
            ln=sum(math.dist(p,q) for p,q in zip(D.path(L),D.path(L)[1:]))
            return (ln,k,mets[0],rx[0])
        cs=[c for c in (cand(i),cand(j)) if c]
        if not cs: cs=[c for c in (cand(k) for k in (i,j)) if c]
        if not cs: raise RuntimeError('crossing between two ring arcs')
        _,k,m,R=max(cs)
        L=D.links[k]; chain=L['chain']
        idx=[t for t,x in enumerate(D.links) if (chain is not None and x['chain']==chain) or t==k]
        into=D.links[idx[-1]]['b']==m or (L['b']==m)
        for t in sorted(idx,reverse=True):
            for e in (D.links[t]['a'],D.links[t]['b']):
                if D.info.get(e,{}).get('kind')=='dummy': D.P.pop(e,None)
            D.links.pop(t)
        a=f'A:cut{n_cut}'; D.add(a,'alias',D.info[m]['label']); n_cut+=1
        plant(D,a,R,into)
    raise RuntimeError('could not remove crossings')

def offending(D,margin=3):
    names=D.real(); Bx=np.array([D.box(n) for n in names]); out=[]
    for i,a,b in segments(D):
        own=D.links[i]['own']; ts=np.linspace(0,1,10); X=a[0]+(b[0]-a[0])*ts; Y=a[1]+(b[1]-a[1])*ts
        hit=(X[None]>Bx[:,0,None]-margin)&(X[None]<Bx[:,2,None]+margin)&(Y[None]>Bx[:,1,None]-margin)&(Y[None]<Bx[:,3,None]+margin)
        for k in np.where(hit.any(1))[0]:
            if names[k] not in own: out+=[names[k]]+sorted(own)
    for i in range(len(names)):
        for j in range(i+1,len(names)):
            A,B=Bx[i],Bx[j]
            if A[0]<B[2]+margin and B[0]<A[2]+margin and A[1]<B[3]+margin and B[1]<A[3]+margin: out+=[names[i],names[j]]
    seen=[]; [seen.append(x) for x in out if x not in seen]
    return seen
def repair(D,fixed,rounds=20):
    def score():
        v=violations(D); return v['line_through_box']+v['box_overlaps']
    for _ in range(rounds):
        cur=score()
        if cur==0: return
        moved=False
        for n in [x for x in offending(D) if x not in fixed]:
            inc={i for i,L in enumerate(D.links) if n in (L['a'],L['b'])}
            x0,y0=D.P[n]; best=None
            for dist in (20,40,80,140,220):
                for k in range(16):
                    t=k*math.pi/8; q=(x0+dist*math.cos(t),y0+dist*math.sin(t))
                    D.P[n]=((x0+q[0])/2,(y0+q[1])/2); ok=not crossing_pairs(D,only=inc)
                    if ok: D.P[n]=q; ok=not crossing_pairs(D,only=inc)
                    c=score() if ok else None
                    D.P[n]=(x0,y0)
                    if ok and c<cur and (best is None or (c,dist)<best[:2]): best=(c,dist,q)
                if best: break
            if best: D.P[n]=best[2]; cur=best[0]; moved=True
            if cur==0: return
        if not moved: return
def scale_until_clean(D,maxit=40):
    for _ in range(maxit):
        v=violations(D)
        if v['line_through_box']==0 and v['box_overlaps']==0: return v
        f=1.1; D.P={n:(p[0]*f,p[1]*f) for n,p in D.P.items()}; D.rings={i:((c[0]*f,c[1]*f),r*f) for i,(c,r) in D.rings.items()}
    return violations(D)

# ============================================================ rendering
def render(D,path,title):
    names=D.real(); B=[D.box(n) for n in names]
    xs=[b[0] for b in B]+[b[2] for b in B]; ys=[b[1] for b in B]+[b[3] for b in B]
    for L in D.links:
        for q in D.path(L): xs.append(q[0]); ys.append(q[1])
    x0,x1,y0,y1=min(xs)-30,max(xs)+30,min(ys)-30,max(ys)+50
    fig=plt.figure(figsize=((x1-x0)/72,(y1-y0)/72)); ax=fig.add_axes([0,0,1,1])
    ax.set_xlim(x0,x1); ax.set_ylim(y0,y1); ax.set_aspect('equal'); ax.axis('off')
    for i,((cx,cy),r) in sorted(D.rings.items()): ax.add_patch(Circle((cx,cy),r,fc='#fbf6f5',ec='none',zorder=0))
    def cut(P,b):
        inb=lambda p:b[0]<p[0]<b[2] and b[1]<p[1]<b[3]; i=0
        while i<len(P)-1 and inb(P[i+1]): i+=1
        if i>=len(P)-1: return P[-2:]
        a,c=P[i],P[i+1]; lo,hi=0.0,1.0
        for _ in range(30):
            m=(lo+hi)/2; q=(a[0]+(c[0]-a[0])*m,a[1]+(c[1]-a[1])*m); (lo,hi)=(m,hi) if inb(q) else (lo,m)
        return [(a[0]+(c[0]-a[0])*hi,a[1]+(c[1]-a[1])*hi)]+P[i+1:]
    for L in D.links:
        P=D.path(L)
        if D.info[L['a']]['kind']!='dummy': P=cut(P,D.box(L['a']))
        if D.info[L['b']]['kind']!='dummy': P=cut(P[::-1],D.box(L['b']))[::-1]
        ring=L['ring'] is not None; side=any(D.info[x]['kind']=='alias' for x in (L['a'],L['b']))
        col,lw=('#b23b36',1.8) if ring else (('#9aa6bd',0.8) if side else ('#4a4a4a',1.1))
        if L['arrow']:
            ax.plot([p[0] for p in P[:-1]],[p[1] for p in P[:-1]],color=col,lw=lw,solid_capstyle='round',zorder=1)
            ax.add_patch(FancyArrowPatch(P[-2],P[-1],arrowstyle='-|>',mutation_scale=9 if ring else 7,color=col,lw=lw,shrinkA=0,shrinkB=0,zorder=1))
        else: ax.plot([p[0] for p in P],[p[1] for p in P],color=col,lw=lw,solid_capstyle='round',zorder=1)
    for n in names:
        k=D.info[n]['kind']; st=STYLE[k]; (x,y)=D.P[n]; w,h=D.size(n)
        ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=0,rounding_size=%g'%(h/2 if k!='rxn' else 2),
                     fc=st['fc'],ec=st['ec'],lw=0.9,ls=st.get('ls','-'),zorder=2))
        lab=D.info[n]['label']; lab=lab if len(lab)<=34 else lab[:32]+'\u2026'
        sub=D.info[n]['sub']; dy=SUB['fs']*0.62 if sub else 0
        ax.text(x,y+dy,lab,ha='center',va='center_baseline',fontsize=st['fs'],weight=st['weight'],style=st['style'],color=st['tc'],zorder=3)
        if sub:
            sub=sub if len(sub)<=60 else sub[:58]+'\u2026'
            ax.text(x,y-st['fs']*0.62,sub,ha='center',va='center_baseline',fontsize=SUB['fs'],color=SUB['tc'],zorder=3)
    ax.text(x0+12,y1-14,title,fontsize=12,va='top')
    fig.savefig(path,**({'metadata':{'Date':None}} if path.endswith('.svg') else {'dpi':min(110.0,7000*72/max(x1-x0,y1-y0))}))
    plt.close(fig)

# ============================================================ clustering and main
def clusters(reactions,max_size):
    if len(reactions)<=max_size: return [('',reactions)]
    groups={}
    for r in reactions:
        if r['subsystem']: groups.setdefault(r['subsystem'],[]).append(r)
    if len(groups)<2:
        M=nx.Graph()
        for i,r in enumerate(reactions):
            M.add_node(i)
            for j,o in enumerate(reactions[:i]):
                if r['main'] and o['main'] and set(r['main'])&set(o['main']): M.add_edge(i,j)
        comms=nx.algorithms.community.greedy_modularity_communities(M,resolution=1.0)
        groups={f'module {k+1}':[reactions[i] for i in sorted(c)] for k,c in enumerate(sorted(comms,key=lambda c:(-len(c),min(c))))}
    else:
        rest=[r for r in reactions if not r['subsystem']]
        if rest: groups['(no subsystem)']=rest
    out=[]; small=[]
    for name in sorted(groups,key=lambda n:(-len(groups[n]),n)):
        g=groups[name]
        if len(g)<3: small+=g
        else:
            for k in range(0,len(g),max_size): out.append((name+(f' ({k//max_size+1})' if len(g)>max_size else ''),g[k:k+max_size]))
    if small:
        for k in range(0,len(small),max_size): out.append(('other reactions',small[k:k+max_size]))
    return out

def main():
    ap=argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('sbml'); ap.add_argument('--out',default=None); ap.add_argument('--max-cluster',type=int,default=60)
    ap.add_argument('--side',default=''); ap.add_argument('--keep',default='')
    a=ap.parse_args()
    name,F,reactions=read_sbml(a.sbml)
    cof=reduce(reactions,F,{_norm(x) for x in a.side.split(',') if x},{_norm(x) for x in a.keep.split(',') if x})
    out=a.out or os.path.splitext(os.path.splitext(os.path.basename(a.sbml))[0])[0]+'_layout'
    parts=clusters(reactions,a.max_cluster)
    log(f'{len(reactions)} reactions, {sum(1 for f in F.values() if f)}/{len(F)} formulas, {len(cof)} cofactors, {len(parts)} drawing(s)')
    ok=True
    for k,(cname,rx) in enumerate(parts):
        D,rings=layout_cluster(rx); cuts=remove_crossings(D)
        repair(D,{n for n in D.P if D.info[n]['kind']=='core'} | {n for n in D.P if D.info[n]['kind']=='dummy'})
        v=scale_until_clean(D)
        fn=out if len(parts)==1 else f"{out}_{k+1:02d}_{re.sub(r'[^A-Za-z0-9]+','_',cname).strip('_')[:40]}"
        segs=[len(D.path(L)) for L in D.links if L['ring'] is None]
        bends=sum(1 for L in D.links if L['ring'] is None and D.info[L['b']]['kind']=='dummy')
        print(f'{fn}: {len(rx)} reactions, rings {[len(r) for r in rings]}, {sum(1 for n in D.P if D.info[n]["kind"]=="alias")} side boxes ({cuts} added to avoid crossings), {bends} bends, {v}')
        if any(v.values()): ok=False; print('  rules not met - drawing not written'); continue
        title=name+(f' \u2014 {cname}' if cname else '')
        render(D,fn+'.svg',title); render(D,fn+'.png',title)
    sys.exit(0 if ok else 1)

if __name__=='__main__': main()

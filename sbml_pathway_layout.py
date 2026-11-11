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

Usage: python3 sbml_pathway_layout.py model.xml[.gz] [--out NAME] [--max-cluster 60] [--show "NH3,CO2"] [--hide "x,y"]
       --show  draw these compounds as shared boxes connected to all their reactions (like main compounds)
       --hide  always write these compounds inside the reaction box (like cofactors)
       --legend  add a legend (only the elements used) at the bottom right of each drawing
       python3 sbml_pathway_layout.py [model.xml] --from-layout NAME.layout.json   # redraw a saved layout
Writes NAME.svg/.png per drawing, NAME.layout.json (reusable layout) and, for SBML Level 3 models,
NAME.layout.sbml (the model plus the SBML Layout package, readable by other SBML tools).
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
SPECIES_ID={}
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
    global SPECIES_ID; SPECIES_ID={s['name']:sid for sid,s in species.items()}
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
SHOW=set()                                                     # normalised names given with --show
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
            Sm=[n for n in S if _norm(n) not in SHOW]; Pm=[n for n in P if _norm(n) not in SHOW]   # shown compounds
            cand=[(similarity(s,p,F),s,p) for s in Sm for p in Pm if s!=p]                         # are never the backbone
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
def sugiyama(nodes,edges,size,gap_x=34,gap_y=34,tag=''):
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
            d=f'__d{tag}{k}_{l}'; layer[d]=l; dsize[d]=(8.0,0.0); pts.append(d)
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
def layout_cluster(rxns,ring_scale=1.0):
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
        rad*=ring_scale
        mw=max(max(D.size('M:'+m)[0] for m in ring),max([D.size('R:'+r['id'])[0] for k in steps for r in steps[k]],default=0))
        ring_steps[i]=dict(ring=ring,steps=steps,chords=chords,rad=rad,box=(2*rad*1.5+mw,2*rad*1.5+40))
    # ---- which compounds get a real (shared) box: every main compound, and every secondary
    # compound used by two or more reactions; compounds used once stay private dashed copies
    appear={}
    for r in rxns:
        for n in set(r['main'] or ())|{x for x,_ in r['secondary']}: appear[n]=appear.get(n,0)+1
    # A secondary compound is connected to its shared box only where that joins parts of the network
    # that are otherwise separate (e.g. carbamoyl phosphate, urea); everywhere else it stays a private copy.
    uf={}
    def find(x):
        while uf.setdefault(x,x)!=x: uf[x]=uf[uf[x]]; x=uf[x]
        return x
    for r in rxns:
        find('R:'+r['id'])
        if r['main']:
            for n in r['main']: uf[find('M:'+n)]=find('R:'+r['id'])
    mains={n for r in rxns for n in (r['main'] or ())}
    joins=set()                                                   # decided as if --show were not given
    cands=sorted(((n,r['id']) for r in rxns for n,_ in r['secondary'] if (appear[n]>=2 or n in mains) and _norm(n) not in SHOW),
                 key=lambda t:(-appear[t[0]],t[0],t[1]))
    for n,rid in cands:
        if find('M:'+n)!=find('R:'+rid): uf[find('M:'+n)]=find('R:'+rid); joins.add((n,rid))
    joins|={(n,r['id']) for r in rxns for n,_ in r['secondary'] if _norm(n) in SHOW}   # then: shown compounds
    real=mains|{n for n,_ in joins}
    for n in sorted(real):
        if 'M:'+n not in D.info: D.add('M:'+n,'core' if n in ringof else 'met',n)
    rring={r['id']:ringof[r['main'][0]] for r in rxns if r['id'] in internal}
    # ---- contracted graph for the layered pass
    unit=lambda n:f'RING{ringof[n]}' if n in ringof else 'M:'+n
    runit=lambda r:f'RING{rring[r["id"]]}' if r['id'] in rring else 'R:'+r['id']
    nodes=set(); edges=[]; emeta=[]; eside=[]; ring_side=[]
    for r in rxns:
        R='R:'+r['id']; Ru=runit(r)
        if r['main'] and r['id'] not in internal:
            a,b=unit(r['main'][0]),unit(r['main'][1]); nodes|={a,b,R}
            edges.append((a,R)); emeta.append(('M:'+r['main'][0],R)); eside.append(False)
            edges.append((R,b)); emeta.append((R,'M:'+r['main'][1])); eside.append(False)
        for n,sd in r['secondary']:
            if (n,r['id']) not in joins: continue
            a=unit(n)
            if a==Ru:                                             # side link inside one ring
                ring_side.append(('M:'+n,R) if sd=='s' else (R,'M:'+n)); continue
            nodes|={a,Ru}
            if sd=='s': edges.append((a,Ru)); emeta.append(('M:'+n,R))
            else: edges.append((Ru,a)); emeta.append((R,'M:'+n))
            eside.append(True)
    for i in ring_steps: nodes.add(f'RING{i}')
    size={}
    for n in nodes:
        size[n]=ring_steps[int(n[4:])]['box'] if n.startswith('RING') else D.size(n)
    comps=sorted((sorted(c) for c in nx.connected_components(nx.Graph([e for e in edges]+[(n,n) for n in nodes]))),key=lambda c:(-len(c),c))
    xoff=0.0
    for comp in comps:
        cs=set(comp); ce=[(e,m,sd) for e,m,sd in zip(edges,emeta,eside) if e[0] in cs]
        pos,chains=sugiyama(comp,[e for e,_,_ in ce],size)
        xs=[pos[n][0]-size.get(n,(8,0))[0]/2 for n in pos]; xe=[pos[n][0]+size.get(n,(8,0))[0]/2 for n in pos]
        dx=xoff-min(xs); xoff+=max(xe)-min(xs)+80
        for n,p in pos.items():
            if n.startswith('__d'): D.add(n,'dummy')
            if not n.startswith('RING'): D.P[n]=(p[0]+dx,p[1])
            else: D.rings[int(n[4:])]=((p[0]+dx,p[1]),ring_steps[int(n[4:])]['rad'])
        for (pts,rv),(e,m,sd) in zip(chains,ce):
            pts=list(pts)
            if rv: pts=pts[::-1]
            ends=[m[0] if pts[0].startswith('RING') else pts[0], m[1] if pts[-1].startswith('RING') else pts[-1]]
            pts[0],pts[-1]=ends
            cid=len(D.links)
            for a,b in zip(pts,pts[1:]):
                D.link(a,b,arrow=(b==pts[-1]),own=set(ends),chain=cid); D.links[-1]['side']=sd
    # ---- side pathways: parts of the network that attach to exactly one ring
    Gc=nx.Graph(); Gc.add_nodes_from(nodes); Gc.add_edges_from(edges)
    ringnodes={i:{'M:'+m for m in rings[i]}|{'R:'+rid for rid,k in rring.items() if k==i} for i in ring_steps}
    att={}                                                        # piece node -> ring nodes it attaches to
    for (a,b),(ma,mb) in zip(edges,emeta):
        if a.startswith('RING'): att.setdefault(b,set()).add(ma)
        if b.startswith('RING'): att.setdefault(a,set()).add(mb)
    pieces=[]
    for i in sorted(ring_steps):
        X=f'RING{i}'
        if X not in Gc: continue
        H=Gc.copy(); H.remove_node(X)
        for comp in sorted((sorted(c) for c in nx.connected_components(H)),key=lambda c:c):
            cs=set(comp); outside={m for n in cs for m in Gc[n] if m not in cs}
            if outside=={X} and not any(n.startswith('RING') for n in cs): pieces.append((i,comp))
    piece_of_ring={i:set().union(*[set(c) for j,c in pieces if j==i]) if any(j==i for j,_ in pieces) else set() for i in ring_steps}
    # ---- expand rings, rotated towards their outside connections
    for i,rs in ring_steps.items():
        (cx,cy),rad=D.rings[i]; ring=rs['ring']; n=len(ring)
        slot={'M:'+m:(j,1.0) for j,m in enumerate(ring)}         # node -> (angle slot, radius factor)
        for k,lst in rs['steps'].items():
            for j,r in enumerate(lst): slot['R:'+r['id']]=(k+0.5,1.0 if j==0 else 1+(0.24 if j%2 else -0.24)*((j+1)//2))
        ext=[]
        for L in D.links:
            if L['own'] & piece_of_ring[i]: continue
            for end,other in ((L['a'],L['b']),(L['b'],L['a'])):
                if end in slot and other in D.P and other not in slot: ext.append((end,D.P[other]))
        anchors=sorted({x for n in piece_of_ring[i] for x in att.get(n,()) if x in slot})
        best=None
        for mir in (1,-1):
            for k in range(72):
                th=2*math.pi*k/72
                p={m:(cx+rad*math.cos(th+mir*2*math.pi*j/n),cy+rad*math.sin(th+mir*2*math.pi*j/n)) for j,m in enumerate(ring)}
                at=lambda e:(cx+rad*slot[e][1]*math.cos(th+mir*2*math.pi*slot[e][0]/n),cy+rad*slot[e][1]*math.sin(th+mir*2*math.pi*slot[e][0]/n))
                sc=sum(math.dist(at(e),q) for e,q in ext)
                sc+=sum(0.8*rad*abs(math.cos(th+mir*2*math.pi*slot[e][0]/n)) for e in anchors)   # anchors at top/bottom
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
                else:                                             # alternative enzymes: first go to the side
                    busy=any(e in ({a,b}|{'R:'+x['id'] for x in lst}) for e,_ in ext)   # that has no outside links
                    first=-0.24 if busy else 0.24
                    f=1+(first if j%2 else -first)*((j+1)//2); D.P[R]=(cx+rad*f*ux,cy+rad*f*uy)
                    D.link(s,R,arrow=True); D.link(R,pp,arrow=True)
        for r in rs['chords']:
            R='R:'+r['id']; s,pp='M:'+r['main'][0],'M:'+r['main'][1]
            D.P[R]=((D.P[s][0]+D.P[pp][0])/2*0.82+cx*0.18,(D.P[s][1]+D.P[pp][1])/2*0.82+cy*0.18)
            D.link(s,R); D.link(R,pp)
    for a,b in ring_side: D.link(a,b); D.links[-1]['side']=True
    # ---- side pathways: own top-to-bottom layout, above or below the ring, aligned with their anchor
    unplaced=set().union(*[set(c) for _,c in pieces]) if pieces else set()
    for pi,(i,comp) in enumerate(pieces):
        cs=set(comp); unplaced-=cs
        keep=[]
        for L in D.links:                                         # drop the provisional drawing of this piece
            if L['ring'] is None and L['own'] & cs and L['own'] <= cs|ringnodes[i]:
                for e in (L['a'],L['b']):
                    if D.info.get(e,{}).get('kind')=='dummy': D.P.pop(e,None)
                continue
            keep.append(L)
        D.links=keep
        (cx,cy),rad=D.rings[i]
        ce=[(e,m,sd) for e,m,sd in zip(edges,emeta,eside) if e[0] in cs and e[1] in cs]
        pos,chains=sugiyama(comp,[e for e,_,_ in ce],size,tag=f'p{pi}_')
        anc=[(n,x) for n in comp for x in sorted(att.get(n,())) if x in ringnodes[i]]
        ax=float(np.mean([D.P[x][0] for _,x in anc])); ay=float(np.mean([D.P[x][1] for _,x in anc]))
        an=[n for n,_ in anc]
        sz=lambda k:(8.0,0.0) if k.startswith('__d') else size[k]
        def boxes(P): return [(v[0]-sz(k)[0]/2,v[1]-sz(k)[1]/2,v[0]+sz(k)[0]/2,v[1]+sz(k)[1]/2) for k,v in P.items()]
        others=[D.box(n) for n in D.P if n not in cs and n not in unplaced and D.info[n]['kind']!='dummy']
        hit=lambda bx,dx,dy:any(b[0]+dx<o[2]+24 and o[0]<b[2]+dx+24 and b[1]+dy<o[3]+24 and o[1]<b[3]+dy+24 for b in bx for o in others)
        vertical=abs(ay-cy)>=abs(ax-cx)                           # anchor at the top/bottom vs. at a side
        if vertical:                                              # above/below, anchor end facing the ring,
            up=ay>=cy                                             # lined up with the anchor -> vertical link
            my=float(np.mean([pos[n][1] for n in an])); oy=float(np.mean([v[1] for v in pos.values()]))
            if (up and my>oy+1e-6) or ((not up) and my<oy-1e-6): pos={k:(v[0],-v[1]) for k,v in pos.items()}
            dx=ax-float(np.mean([pos[n][0] for n in an])); pos={k:(v[0]+dx,v[1]) for k,v in pos.items()}
            bx=boxes(pos); x0,x1=min(b[0] for b in bx),max(b[2] for b in bx)
            clear=max(math.sqrt(max(0.0,(rad+38)**2-(x-cx)**2)) for x in np.linspace(x0,x1,25))
            shift=(max(cy+clear,ay+40)-min(b[1] for b in bx)) if up else (min(cy-clear,ay-40)-max(b[3] for b in bx))
            for _ in range(400):
                if not hit(bx,0,shift): break
                shift+=30 if up else -30
            sx,sy=0.0,shift
        else:                                                     # beside the ring, still top-to-bottom,
            right=ax>=cx                                          # anchor end level with the anchor -> horizontal link
            dy=ay-float(np.mean([pos[n][1] for n in an])); pos={k:(v[0],v[1]+dy) for k,v in pos.items()}
            bx=boxes(pos); y0,y1=min(b[1] for b in bx),max(b[3] for b in bx)
            clear=max(math.sqrt(max(0.0,(rad+38)**2-(y-cy)**2)) for y in np.linspace(y0,y1,25))
            mx=float(np.mean([pos[n][0] for n in an]))
            shift=(max(cx+clear,ax+50)-min(b[0] for b in bx)) if right else (min(cx-clear,ax-50)-max(b[2] for b in bx))
            for _ in range(400):
                if not hit(bx,shift,0): break
                shift+=30 if right else -30
            sx,sy=shift,0.0
        for k,v in pos.items():
            if k.startswith('__d'): D.add(k,'dummy')
            D.P[k]=(v[0]+sx,v[1]+sy)
        for (pts,rv),(e,m,sd) in zip(chains,ce):
            pts=list(pts)[::-1] if rv else list(pts); cid=len(D.links)+100000*(pi+1)
            for a_,b_ in zip(pts,pts[1:]):
                D.link(a_,b_,arrow=(b_==pts[-1]),own={pts[0],pts[-1]},chain=cid); D.links[-1]['side']=sd
        for e,m,sd in zip(edges,emeta,eside):                     # links from the piece to the ring
            if (e[0] in cs)!=(e[1] in cs) and (e[0] in cs or e[1] in cs): D.link(m[0],m[1]); D.links[-1]['side']=sd
    # ---- reactions without a main pair, and secondary participants: side boxes
    pending=[]
    k=0
    for r in rxns:
        R='R:'+r['id']
        if R not in D.P:                                          # touches no shared box at all
            D.P[R]=(xoff+150*k,300.0); k+=1
        for j,(n,sd) in enumerate(r['secondary']):
            if (n,r['id']) in joins: continue
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
def place_hubs(D,hubs,grid=26):
    """Move each shared compound box to the spot where its own lines cross nothing and pass through no
    box, preferring short lines: a grid search over the drawing (vectorised)."""
    for h in hubs:
        inc=[i for i,L in enumerate(D.links) if h in (L['a'],L['b'])]
        if not inc or any(D.links[i]['chain'] is not None or D.links[i]['ring'] is not None for i in inc): continue
        nbr=[D.links[i]['b'] if D.links[i]['a']==h else D.links[i]['a'] for i in inc]
        others=[]
        for i,L in enumerate(D.links):
            if i in inc: continue
            p=D.path(L)
            for a,b in zip(p,p[1:]): others.append((a,b,L['a'],L['b']))
        A=np.array([o[0] for o in others]) if others else np.zeros((0,2)); B=np.array([o[1] for o in others]) if others else np.zeros((0,2))
        names=[n for n in D.real() if n!=h]; Bx=np.array([D.box(n) for n in names]); w,hh=D.size(h)
        xs=[b[0] for b in Bx]+[b[2] for b in Bx]; ys=[b[1] for b in Bx]+[b[3] for b in Bx]
        cand=[(x,y) for x in np.linspace(min(xs)-w,max(xs)+w,grid) for y in np.linspace(min(ys)-hh*3,max(ys)+hh*3,grid)]
        cand+=[D.P[h]]+[(D.P[r][0]+dx,D.P[r][1]+dy) for r in nbr for dx in (-90,0,90) for dy in (-60,0,60) if dx or dy]
        best=None
        for c in cand:
            cr=0; thr=0
            for r in nbr:
                q=D.P[r]
                if len(A):
                    mask=np.array([r not in (o[2],o[3]) for o in others])
                    cr+=int(_segx(np.array(c)[None],np.array(q)[None],A[mask],B[mask]).sum())
                ts=np.linspace(0.05,0.95,14); X=c[0]+(q[0]-c[0])*ts; Y=c[1]+(q[1]-c[1])*ts
                hit=(X[None]>Bx[:,0,None]-3)&(X[None]<Bx[:,2,None]+3)&(Y[None]>Bx[:,1,None]-3)&(Y[None]<Bx[:,3,None]+3)
                thr+=sum(1 for k in np.where(hit.any(1))[0] if names[k]!=r)
            bx=(c[0]-w/2,c[1]-hh/2,c[0]+w/2,c[1]+hh/2)
            ov=int(((Bx[:,0]<bx[2]+6)&(bx[0]<Bx[:,2]+6)&(Bx[:,1]<bx[3]+6)&(bx[1]<Bx[:,3]+6)).sum())
            # boxes crossed by OTHER lines at this spot
            if len(A):
                ts=np.linspace(0,1,8); SX=A[:,None,0]+(B[:,None,0]-A[:,None,0])*ts; SY=A[:,None,1]+(B[:,None,1]-A[:,None,1])*ts
                thr+=int(((SX>bx[0]-3)&(SX<bx[2]+3)&(SY>bx[1]-3)&(SY<bx[3]+3)).any(1).sum())
            ln=sum(math.dist(c,D.P[r]) for r in nbr)
            skew=sum(min(abs(c[0]-D.P[r][0]),abs(c[1]-D.P[r][1])) for r in nbr)   # 0 for straight vertical/horizontal
            key=(cr,thr+ov,round(ln+2.0*skew,3),round(c[0],3),round(c[1],3))
            if best is None or key<best: best=key
        D.P[h]=(best[3],best[4])


def remove_crossings(D):
    """Cut a crossing line at its metabolite end and give that end its own copy (dashed box)."""
    n_cut=0
    for _ in range(500):
        cp=crossing_pairs(D)
        if not cp:
            used={e for L in D.links for e in (L['a'],L['b'])}       # shared boxes left without any line
            for n in [n for n in D.P if D.info[n]['kind'] in ('met','alias') and n not in used]:
                D.P.pop(n)                                        # (every use became a dashed copy)
            return n_cut
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
    D.scaled=getattr(D,'scaled',1.0)
    for _ in range(maxit):
        v=violations(D)
        if v['line_through_box']==0 and v['box_overlaps']==0: return v
        f=1.1; D.P={n:(p[0]*f,p[1]*f) for n,p in D.P.items()}; D.rings={i:((c[0]*f,c[1]*f),r*f) for i,(c,r) in D.rings.items()}
        D.scaled*=f
    return violations(D)

# ============================================================ rendering
LEG_FS=8.0
def legend_entries(D):
    kinds={D.info[n]['kind'] for n in D.real()}
    ring=any(L['ring'] is not None for L in D.links)
    side=any(L['ring'] is None and (L.get('side') or any(D.info[x]['kind']=='alias' for x in (L['a'],L['b']))) for L in D.links)
    main=any(L['ring'] is None and not L.get('side') and not any(D.info[x]['kind']=='alias' for x in (L['a'],L['b'])) for L in D.links)
    cof=any(D.info[n]['sub'] for n in D.real() if D.info[n]['kind']=='rxn')
    E=[]
    if 'core' in kinds: E.append(('box','core','cycle intermediate'))
    if 'met' in kinds: E.append(('box','met','metabolite'))
    if 'alias' in kinds: E.append(('box','alias','copy of a compound drawn elsewhere or used once'))
    if 'rxn' in kinds: E.append(('rxn',cof,'reaction (enzyme)' + ('; second line: cofactors consumed → released' if cof else '')))
    if ring: E.append(('line','ring','cycle step'))
    if main: E.append(('line','main','main substrate → product'))
    if side: E.append(('line','side','secondary compound'))
    return E
def legend_size(E):
    sw=74; tw=max(text_w(t,dict(fs=LEG_FS,weight='normal',style='normal')) for *_,t in E)
    return sw+12+tw+24, 22*len(E)+34
def draw_legend(ax,E,x0,y0):
    """Legend box with its lower-left corner at (x0, y0); swatches use the drawing's own styles."""
    w,h=legend_size(E)
    ax.add_patch(FancyBboxPatch((x0,y0),w,h,boxstyle='round,pad=0,rounding_size=4',fc='#ffffff',ec='#b0b0b0',lw=0.8,zorder=5))
    ax.text(x0+10,y0+h-12,'Legend',fontsize=9,weight='bold',va='center',zorder=6)
    for k,(typ,what,text) in enumerate(E):
        cy=y0+h-34-22*k; sx=x0+10
        if typ=='box':
            st=STYLE[what]; bw,bh=58,15
            ax.add_patch(FancyBboxPatch((sx+8,cy-bh/2),bw,bh,boxstyle='round,pad=0,rounding_size=%g'%(bh/2),fc=st['fc'],ec=st['ec'],lw=0.9,ls=st.get('ls','-'),zorder=6))
            ax.text(sx+8+bw/2,cy,'name',ha='center',va='center_baseline',fontsize=st['fs']*0.8,weight=st['weight'],color=st['tc'],zorder=7)
        elif typ=='rxn':
            st=STYLE['rxn']; bw,bh=(66,20) if what else (58,14)
            ax.add_patch(FancyBboxPatch((sx+4,cy-bh/2),bw,bh,boxstyle='round,pad=0,rounding_size=2',fc=st['fc'],ec=st['ec'],lw=0.9,zorder=6))
            ax.text(sx+4+bw/2,cy+(3.5 if what else 0),'enzyme',ha='center',va='center_baseline',fontsize=6.5,style='italic',color=st['tc'],zorder=7)
            if what: ax.text(sx+4+bw/2,cy-4.5,'cofactors',ha='center',va='center_baseline',fontsize=5.2,color=SUB['tc'],zorder=7)
        else:
            col,lw={'ring':('#b23b36',1.8),'main':('#4a4a4a',1.1),'side':('#9aa6bd',0.8)}[what]
            ax.add_patch(FancyArrowPatch((sx+6,cy),(sx+68,cy),arrowstyle='-|>',mutation_scale=9 if what=='ring' else 7,color=col,lw=lw,zorder=6))
        ax.text(sx+86,cy,text,fontsize=LEG_FS,va='center_baseline',color='#222222',zorder=7)

def render(D,path,title,legend=False):
    names=D.real(); B=[D.box(n) for n in names]
    xs=[b[0] for b in B]+[b[2] for b in B]; ys=[b[1] for b in B]+[b[3] for b in B]
    for L in D.links:
        for q in D.path(L): xs.append(q[0]); ys.append(q[1])
    x0,x1,y0,y1=min(xs)-30,max(xs)+30,min(ys)-30,max(ys)+50
    leg=None
    if legend:
        E=legend_entries(D); lw_,lh=legend_size(E)
        lx,ly=x1-12-lw_,y0+12                                     # bottom-right corner of the drawing
        pts=[q for L in D.links for q in D.path(L)]
        for L in D.links:                                         # sample the lines, not just their ends
            p=D.path(L)
            for a,b in zip(p,p[1:]):
                for t in np.linspace(0,1,12): pts.append((a[0]+(b[0]-a[0])*t,a[1]+(b[1]-a[1])*t))
        free=not any(lx-8<q[0]<lx+lw_+8 and ly-8<q[1]<ly+lh+8 for q in pts) and \
             not any(b[0]<lx+lw_+8 and lx-8<b[2] and b[1]<ly+lh+8 and ly-8<b[3] for b in B)
        if not free:                                              # corner in use: extend the canvas downwards
            ly=y0-lh-10; y0=ly-12
            if x1-x0<lw_+24: x0=x1-lw_-24; lx=x1-12-lw_
        leg=(E,lx,ly)
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
        ring=L['ring'] is not None; side=L.get('side') or any(D.info[x]['kind']=='alias' for x in (L['a'],L['b']))
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
    if leg: draw_legend(ax,leg[0],leg[1],leg[2])
    fig.savefig(path,**({'metadata':{'Date':None}} if path.endswith('.svg') else {'dpi':min(110.0,7000*72/max(x1-x0,y1-y0))}))
    plt.close(fig)

# ============================================================ export / reuse
def sbml_ref(D,n):
    k=D.info[n]['kind']
    if n.startswith('R:'): return ('reaction',n[2:])
    if k in ('core','met','alias'): return ('species',SPECIES_ID.get(D.info[n]['label']))
    return (None,None)
def drawing_to_dict(D,name,title,checks):
    nodes=[]
    for n in sorted(D.P):
        w,h=D.size(n); t,ref=sbml_ref(D,n)
        nodes.append(dict(id=n,kind=D.info[n]['kind'],label=D.info[n]['label'],sub=D.info[n]['sub'],
                          x=D.P[n][0],y=D.P[n][1],w=round(w,3),h=round(h,3),sbml_type=t,sbml_id=ref))
    links=[dict(source=L['a'],target=L['b'],arrow=L['arrow'],ring=L['ring'],chain=L['chain'],own=sorted(L['own']),side=bool(L.get('side'))) for L in D.links]
    rings=[dict(id=i,cx=c[0],cy=c[1],r=r) for i,(c,r) in sorted(D.rings.items())]
    return dict(name=name,title=title,checks=checks,rings=rings,nodes=nodes,links=links)
def drawing_from_dict(d):
    D=Drawing()
    for n in d['nodes']:
        D.add(n['id'],n['kind'],n['label'],n['sub']); D.P[n['id']]=(n['x'],n['y'])
        if n.get('sbml_type')=='species' and n.get('sbml_id'): SPECIES_ID.setdefault(n['label'],n['sbml_id'])
    for L in d['links']:
        D.links.append(dict(a=L['source'],b=L['target'],arrow=L['arrow'],ring=L['ring'],own=set(L['own']),chain=L['chain'],cut=True,side=L.get('side',False)))
    D.rings={r['id']:((r['cx'],r['cy']),r['r']) for r in d['rings']}
    return D

LNS='http://www.sbml.org/sbml/level3/version1/layout/version1'
XSI='http://www.w3.org/2001/XMLSchema-instance'
def export_sbml_layout(src,dst,drawings):
    """Copy the model and add one SBML Layout (L3 'layout' package) per drawing."""
    nsmap={}
    opener=(lambda: gzip.open(src)) if src.endswith('.gz') else (lambda: open(src,'rb'))
    with opener() as fh:
        for ev,(pre,uri) in ET.iterparse(fh,events=('start-ns',)): nsmap.setdefault(pre,uri)
    for pre,uri in nsmap.items(): ET.register_namespace(pre,uri)
    ET.register_namespace('layout',LNS); ET.register_namespace('xsi',XSI)
    with opener() as fh: tree=ET.parse(fh)
    root=tree.getroot()
    if root.get('level')!='3': return False
    q=lambda t:f'{{{LNS}}}{t}'
    root.set(q('required'),'false')
    model=next(c for c in root if _l(c.tag)=='model')
    old=[c for c in model if _l(c.tag)=='listOfLayouts']
    for o in old: model.remove(o)
    LL=ET.SubElement(model,q('listOfLayouts'))
    def bbox(parent,x,y,w,h):
        b=ET.SubElement(parent,q('boundingBox')); p=ET.SubElement(b,q('position')); p.set(q('x'),f'{x:.2f}'); p.set(q('y'),f'{y:.2f}')
        d=ET.SubElement(b,q('dimensions')); d.set(q('width'),f'{w:.2f}'); d.set(q('height'),f'{h:.2f}')
    def pt(parent,tag,x,y):
        e=ET.SubElement(parent,q(tag)); e.set(q('x'),f'{x:.2f}'); e.set(q('y'),f'{y:.2f}')
    for k,(D,name,title) in enumerate(drawings):
        names=D.real(); B={n:D.box(n) for n in names}
        allx=[b[0] for b in B.values()]+[b[2] for b in B.values()]; ally=[b[1] for b in B.values()]+[b[3] for b in B.values()]
        for L in D.links:
            for p_ in D.path(L): allx.append(p_[0]); ally.append(p_[1])
        X0,Y1=min(allx)-20,max(ally)+20
        tx=lambda x:x-X0; ty=lambda y:Y1-y                        # SBML layout: origin top-left, y down
        lay=ET.SubElement(LL,q('layout')); lay.set(q('id'),f'layout_{k+1}'); lay.set(q('name'),title)
        d=ET.SubElement(lay,q('dimensions')); d.set(q('width'),f'{max(allx)+20-X0:.2f}'); d.set(q('height'),f'{Y1-min(ally)+20:.2f}')
        SG=ET.SubElement(lay,q('listOfSpeciesGlyphs')); RG=ET.SubElement(lay,q('listOfReactionGlyphs'))
        TG=ET.SubElement(lay,q('listOfTextGlyphs'))
        gid={}
        for i,n in enumerate(names):
            t,ref=sbml_ref(D,n); x0,y0,x1,y1=B[n]
            if t=='species' and ref:
                g=ET.SubElement(SG,q('speciesGlyph')); gid[n]=f'sg_{k+1}_{i}'; g.set(q('id'),gid[n]); g.set(q('species'),ref)
                bbox(g,tx(x0),ty(y1),x1-x0,y1-y0)
                tg=ET.SubElement(TG,q('textGlyph')); tg.set(q('id'),f'tg_{k+1}_{i}'); tg.set(q('graphicalObject'),gid[n]); tg.set(q('originOfText'),ref)
                bbox(tg,tx(x0),ty(y1),x1-x0,y1-y0)
        # every metabolite<->reaction connection, following dummy chains and ring arcs
        conns={}
        for L in D.links:
            key=L['chain'] if L['chain'] is not None else id(L)
            conns.setdefault(key,[]).append(L)
        rglyph={}
        for i,n in enumerate(names):
            if n.startswith('R:'):
                x0,y0,x1,y1=B[n]; g=ET.SubElement(RG,q('reactionGlyph')); g.set(q('id'),f'rg_{k+1}_{i}'); g.set(q('reaction'),n[2:])
                bbox(g,tx(x0),ty(y1),x1-x0,y1-y0); rglyph[n]=(g,ET.SubElement(g,q('listOfSpeciesReferenceGlyphs')))
                tg=ET.SubElement(TG,q('textGlyph')); tg.set(q('id'),f'tgr_{k+1}_{i}'); tg.set(q('graphicalObject'),f'rg_{k+1}_{i}')
                tg.set(q('text'),D.info[n]['label']+(' | '+D.info[n]['sub'] if D.info[n]['sub'] else '')); bbox(tg,tx(x0),ty(y1),x1-x0,y1-y0)
        nref=0
        for key,Ls in conns.items():
            own=sorted(Ls[0]['own']); R=[x for x in own if x.startswith('R:')]; S=[x for x in own if x in gid]
            if not R or not S: continue
            R,S=R[0],S[0]; into_species=Ls[-1]['b']==S
            role=('sideproduct' if into_species else 'sidesubstrate') if (D.info[S]['kind']=='alias' or Ls[0].get('side')) else ('product' if into_species else 'substrate')
            if R not in rglyph: continue
            sr=ET.SubElement(rglyph[R][1],q('speciesReferenceGlyph')); sr.set(q('id'),f'srg_{k+1}_{nref}'); nref+=1
            sr.set(q('speciesGlyph'),gid[S]); sr.set(q('role'),role)
            cv=ET.SubElement(ET.SubElement(sr,q('curve')),q('listOfCurveSegments'))
            for L in Ls:
                if L['ring'] is None:
                    a,b=D.P[L['a']],D.P[L['b']]
                    seg=ET.SubElement(cv,q('curveSegment')); seg.set(f'{{{XSI}}}type','LineSegment')
                    pt(seg,'start',tx(a[0]),ty(a[1])); pt(seg,'end',tx(b[0]),ty(b[1]))
                else:                                             # ring arc -> cubic Beziers of <= 90 degrees
                    (cx,cy),r=D.rings[L['ring']]; a,b=D.P[L['a']],D.P[L['b']]
                    t1=math.atan2(a[1]-cy,a[0]-cx); t2=math.atan2(b[1]-cy,b[0]-cx); dt=(t2-t1+math.pi)%(2*math.pi)-math.pi
                    m=max(1,int(math.ceil(abs(dt)/(math.pi/2)))); h=dt/m; kk=4/3*math.tan(h/4)*r
                    for j in range(m):
                        u,v=t1+h*j,t1+h*(j+1)
                        P0=(cx+r*math.cos(u),cy+r*math.sin(u)); P3=(cx+r*math.cos(v),cy+r*math.sin(v))
                        P1=(P0[0]-kk*math.sin(u),P0[1]+kk*math.cos(u)); P2=(P3[0]+kk*math.sin(v),P3[1]-kk*math.cos(v))
                        seg=ET.SubElement(cv,q('curveSegment')); seg.set(f'{{{XSI}}}type','CubicBezier')
                        pt(seg,'start',tx(P0[0]),ty(P0[1])); pt(seg,'end',tx(P3[0]),ty(P3[1]))
                        pt(seg,'basePoint1',tx(P1[0]),ty(P1[1])); pt(seg,'basePoint2',tx(P2[0]),ty(P2[1]))
        AG=None
        for g,lst in rglyph.values():                             # unconnected reaction -> generalGlyph
            if len(lst)==0:                                       # (a reactionGlyph needs >=1 connection)
                if AG is None: AG=ET.Element(q('listOfAdditionalGraphicalObjects'))
                gg=ET.SubElement(AG,q('generalGlyph')); gg.set(q('id'),g.get(q('id'))); gg.set(q('reference'),g.get(q('reaction')))
                gg.append(g.find(q('boundingBox'))); RG.remove(g)
        if AG is not None: lay.insert(list(lay).index(TG)+1,AG)
        for lst in (SG,RG,TG):
            if len(lst)==0: lay.remove(lst)
    tree.write(dst,xml_declaration=True,encoding='UTF-8')
    return True

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
    import json
    ap=argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('sbml',nargs='?'); ap.add_argument('--out',default=None); ap.add_argument('--max-cluster',type=int,default=60)
    ap.add_argument('--show',default='',help='compounds to draw as shared boxes connected to all their reactions, e.g. "NH3,CO2"')
    ap.add_argument('--hide',default='',help='compounds to always write inside the reaction box')
    ap.add_argument('--side',default='',help=argparse.SUPPRESS); ap.add_argument('--keep',default='',help=argparse.SUPPRESS)  # old names
    ap.add_argument('--from-layout',default=None,help='redraw a saved NAME.layout.json (positions may have been edited) instead of computing a layout')
    ap.add_argument('--no-sbml-layout',action='store_true',help='do not write NAME.layout.sbml')
    ap.add_argument('--legend',action='store_true',help='add a legend at the bottom right of each drawing')
    a=ap.parse_args()
    if a.from_layout:                                             # ---- reuse a saved layout
        saved=json.load(open(a.from_layout)); out=a.out or a.from_layout.replace('.layout.json','')+'_redrawn'
        ok=True; drawings=[]
        for k,d in enumerate(saved['drawings']):
            D=drawing_from_dict(d); v=violations(D)
            fn=out if len(saved['drawings'])==1 else f'{out}_{k+1:02d}'
            print(f'{fn}: {v}')
            if any(v.values()): ok=False; print('  rules not met (was the layout edited?) - drawn anyway for inspection')
            render(D,fn+'.svg',d['title'],a.legend); render(D,fn+'.png',d['title'],a.legend); drawings.append((D,d['name'],d['title']))
        if a.sbml and not a.no_sbml_layout:
            read_sbml(a.sbml)
            if export_sbml_layout(a.sbml,out+'.layout.sbml',drawings): print('wrote',out+'.layout.sbml')
        sys.exit(0 if ok else 1)
    if not a.sbml: ap.error('an SBML file is required unless --from-layout is given')
    name,F,reactions=read_sbml(a.sbml)
    names=lambda t:{_norm(x) for x in t.split(',') if x.strip()}
    SHOW.update(names(a.show)); hide=names(a.hide)|names(a.side); keep=names(a.keep)|SHOW
    known={_norm(n) for n in F}
    for x in sorted((SHOW|hide)-known): log(f'warning: no compound named "{x}" in this model')
    cof=reduce(reactions,F,hide,keep)
    out=a.out or os.path.splitext(os.path.splitext(os.path.basename(a.sbml))[0])[0]+'_layout'
    parts=clusters(reactions,a.max_cluster)
    log(f'{len(reactions)} reactions, {sum(1 for f in F.values() if f)}/{len(F)} formulas, {len(cof)} cofactors, {len(parts)} drawing(s)')
    ok=True; saved=[]; drawings=[]
    for k,(cname,rx) in enumerate(parts):
        best=None
        for ring_scale in (1.0,1.3,1.6,2.0):                      # crowded rings: retry with a larger ring
            D,rings=layout_cluster(rx,ring_scale); cuts=0
            side_ends={e for L in D.links if L.get('side') for e in (L['a'],L['b'])}
            hubs=sorted((n for n in D.P if D.info[n]['kind']=='met' and (_norm(D.info[n]['label']) in SHOW or n in side_ends)),
                        key=lambda n:(_norm(D.info[n]['label']) not in SHOW,n))
            place_hubs(D,hubs)                                    # shared compounds: best free spot
            for _ in range(6):                                    # repeat until all three checks hold together
                cuts+=remove_crossings(D)                         # (a near-touch can tip into a crossing when scaled)
                fixed={n for n in D.P if D.info[n]['kind'] in ('core','dummy')}|{e for L in D.links if L['ring'] is not None for e in (L['a'],L['b'])}
                repair(D,fixed)                                   # ring members and ring-arc reactions never move
                v=scale_until_clean(D)
                if not any(v.values()): break
            key=(any(v.values()),D.scaled>1.35,cuts,D.scaled)
            if best is None or key<best[0]: best=(key,D,rings,cuts,v)
            if not rings or (not any(v.values()) and D.scaled<=1.35): break
        _,D,rings,cuts,v=best
        fn=out if len(parts)==1 else f"{out}_{k+1:02d}_{re.sub(r'[^A-Za-z0-9]+','_',cname).strip('_')[:40]}"
        bends=sum(1 for L in D.links if L['ring'] is None and D.info[L['b']]['kind']=='dummy')
        print(f'{fn}: {len(rx)} reactions, rings {[len(r) for r in rings]}, {sum(1 for n in D.P if D.info[n]["kind"]=="alias")} side boxes ({cuts} added to avoid crossings), {bends} bends, {v}')
        if any(v.values()): ok=False; print('  rules not met - drawing not written'); continue
        title=name+(f' \u2014 {cname}' if cname else '')
        render(D,fn+'.svg',title,a.legend); render(D,fn+'.png',title,a.legend)
        saved.append(drawing_to_dict(D,cname,title,v)); drawings.append((D,cname,title))
    meta=dict(format='sbml-pathway-layout',version=1,source=os.path.basename(a.sbml),model=name,
              settings=dict(max_cluster=a.max_cluster,show=a.show,hide=','.join(x for x in (a.hide,a.side) if x),legend=a.legend),
              units='points; x right, y up; boxes are centred on (x, y)',drawings=saved)
    with open(out+'.layout.json','w') as fh: json.dump(meta,fh,indent=1,sort_keys=True)
    print('wrote',out+'.layout.json')
    if not a.no_sbml_layout:
        if export_sbml_layout(a.sbml,out+'.layout.sbml',drawings): print('wrote',out+'.layout.sbml (SBML Layout package)')
        else: print('SBML Layout export skipped: the model is not SBML Level 3')
    sys.exit(0 if ok else 1)

if __name__=='__main__': main()

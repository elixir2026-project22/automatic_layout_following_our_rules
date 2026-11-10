#!/usr/bin/env python3
"""Crossing-free pathway drawings from any SBML file, with the main cycle drawn as a circle.

Rules: reproducible | labels at fixed positions per node type | no line crossings |
try to avoid long lines | try to avoid bends.

Method (per connected component of the reaction network)
 1. Side compounds (H2O, ATP, NAD(P)H, CoA, ...) are written inside the reaction box, not drawn as nodes.
 2. The main cycle is chosen (longest metabolite cycle, or --cycle) and placed on a circle.
 3. Metabolites are duplicated (dashed boxes) only where needed: to make the network planar and to
    turn side branches into trees. Duplication never touches the main cycle.
 4. Side routes attached to the cycle at two or more points are split into "inside" and "outside"
    sets so that no two on the same side interleave. Each side is drawn with Tutte's barycentric
    method on a stellated (triangulated) embedding: provably crossing-free. The outside set is drawn
    in a disc and mapped outward by radial reflection, so its lines become smooth curves.
 5. Trees are attached in free angular gaps; a deterministic local search then spreads nodes,
    rejecting every move that would create a crossing; finally the drawing is uniformly scaled until
    no box overlaps another and no line passes through a box (scaling cannot create crossings).
 6. A geometric check verifies 0 crossings / 0 lines through boxes / 0 box overlaps; otherwise no output.
Nothing is random; Python's hash seed is pinned and the SVG is written without a date.

Usage: python3 sbml_cycle_layout.py model.xml [--out NAME] [--cycle "A,B,C,..."] [--side "x,y"] [--keep "x,y"]
Requires: networkx>=3.1, numpy, matplotlib
"""
import os, sys
if os.environ.get('PYTHONHASHSEED')!='0':
    os.environ['PYTHONHASHSEED']='0'; os.execv(sys.executable,[sys.executable]+sys.argv)
import re, math, argparse, itertools
import xml.etree.ElementTree as ET
import numpy as np, networkx as nx
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
plt.rcParams['svg.hashsalt']='sbml-layout'; plt.rcParams['font.family']='DejaVu Sans'

# =============================================================== SBML (levels 1-3, any version)
def _l(tag): return tag.rsplit('}',1)[-1]
def _kids(e,name): return [c for c in e if _l(c.tag)==name]
def read_sbml(path):
    import gzip
    root=(ET.parse(gzip.open(path)) if path.endswith('.gz') else ET.parse(path)).getroot()
    model=next(c for c in root if _l(c.tag)=='model')
    species={}
    for e in model.iter():
        if _l(e.tag) in ('species','specie'):
            sid=e.get('id') or e.get('name'); species[sid]=dict(name=(e.get('name') or sid).strip(),comp=e.get('compartment') or '')
    # same name in several compartments -> make labels unique
    byname={}
    for sid,s in species.items(): byname.setdefault(s['name'],[]).append(sid)
    for nm,ids in byname.items():
        if len({species[i]['comp'] for i in ids})>1:
            for i in ids: species[i]['name']=f"{nm} [{species[i]['comp']}]"
    groups={}
    for e in model.iter():
        if _l(e.tag)=='group':
            nm=None
            for k,v in e.attrib.items():
                if _l(k)=='name': nm=v
            if not nm: continue
            for m in e.iter():
                if _l(m.tag)=='member':
                    for k,v in m.attrib.items():
                        if _l(k)=='idRef': groups.setdefault(v,[]).append(nm)
    mname=(model.get('name') or '').strip().lower()
    reactions=[]
    for e in model.iter():
        if _l(e.tag)!='reaction': continue
        def refs(lst):
            out=[]
            for L in _kids(e,lst):
                for r in L:
                    if _l(r.tag) in ('speciesReference','specieReference'):
                        sid=r.get('species') or r.get('specie')
                        st=r.get('stoichiometry'); st=float(st) if st not in (None,'') else 1.0
                        if sid in species: out.append((species[sid]['name'],st))
            return out
        rid=e.get('id') or e.get('name'); nm=(e.get('name') or rid).strip()
        if re.fullmatch(r'[A-Za-z_]{0,6}\d+[A-Za-z0-9_]*',nm) or nm==rid:      # looks like an identifier
            alts=[x for x in groups.get(rid,[]) if x.strip().lower()!=mname]
            parts_=[species[r.get('species') or r.get('specie')]['name'].lower() for L in list(_kids(e,'listOfReactants'))+list(_kids(e,'listOfProducts'))
                    for r in L if (r.get('species') or r.get('specie')) in species]
            toks={t for p in parts_ for t in re.findall(r'[a-z0-9+]{3,}',p)}
            def fit(x):               # prefer the class name that mentions the reaction's own compounds
                words=set(re.findall(r'[a-z0-9+]{3,}',x.lower())); return (-len(toks&words),len(x),x)
            if alts: nm=re.sub(r'^EC [0-9.\-nB]+\s+','',min(alts,key=fit))
        reactions.append(dict(id=rid,name=nm,subs=refs('listOfReactants'),prods=refs('listOfProducts'),
                              reversible=(e.get('reversible','true')!='false')))
    return (model.get('name') or model.get('id') or 'model'),reactions

# =============================================================== side compounds
SIDE=r'''h2o water oh- hydroxide h+ proton hydron atp adp amp gtp gdp gmp utp udp ump ctp cdp cmp itp idp dtp
nad nad+ nadh nadp nadp+ nadph fad fadh2 fmn fmnh2 coa coenzymea coa-sh coash co2 carbondioxide hco3- bicarbonate hydrogencarbonate
phosphate orthophosphate pi hydrogenphosphate diphosphate pyrophosphate ppi o2 oxygen dioxygen h2o2 hydrogenperoxide
nh3 nh4+ ammonia ammonium aquinone aquinol quinone quinol ubiquinone ubiquinol menaquinone menaquinol
oxidizedferredoxin reducedferredoxin ferredoxin oxidizedthioredoxin reducedthioredoxin thioredoxin acceptor reducedacceptor
electron e- na+ k+ cl- fe2+ fe3+ mg2+ h2 hydrogen
nicotinamideadeninedinucleotide nicotinamideadeninedinucleotidereduced nicotinamideadeninedinucleotidephosphate
nicotinamideadeninedinucleotidephosphatereduced flavinadeninedinucleotide flavinadeninedinucleotideoxidized
flavinadeninedinucleotidereduced adenosinetriphosphate adenosinediphosphate adenosinemonophosphate atp(4-) adp(3-)
adenosine5'-triphosphate adenosine5'-diphosphate adenosine5'-monophosphate guanosinetriphosphate guanosinediphosphate
coenzyme-a coenzymea(coa) ubiquinone-8 ubiquinol-8 menaquinone-8 menaquinol-8 h2o(water) watermolecule'''.split()
def _norm(n):
    n=n.lower(); n=re.sub(r'\[[^\]]*\]$','',n); n=re.sub(r'_[a-z]$','',n)
    n=n.replace('(+)','+').replace('(-)','-').replace('(2+)','2+').replace('(3+)','3+')
    return re.sub(r"[\s_\-,'()]",'',n)
SIDE_RE=re.compile(r'^((ubi|mena|plasto)?quino(ne|l)(-?\d+)?|(oxidized|reduced)(ferredoxin|thioredoxin|flavodoxin|cytochromec?)|cytochromec(\(ii+\))?)$')
SIDE={_norm(x) for x in SIDE}
def is_side(name,extra,keep):
    n=_norm(name)
    if n in keep: return False
    return n in extra or n in SIDE or bool(SIDE_RE.match(n))

# =============================================================== graph
def build_graph(reactions,extra=(),keep=()):
    extra={_norm(x) for x in extra}; keep={_norm(x) for x in keep}
    g=nx.Graph(); arrows=[]
    for i,r in enumerate(reactions):
        R=f'R{i:05d}'; label=r['name'] if len(r['name'])<=34 else r['name'][:32]+'…'
        fmt=lambda lst:' + '.join((f'{s:g} ' if s!=1 else '')+n for n,s in lst)
        cs=[x for x in r['subs'] if is_side(x[0],extra,keep)]; cp=[x for x in r['prods'] if is_side(x[0],extra,keep)]
        sub=(fmt(cs)+' \u2192 '+fmt(cp)) if (cs or cp) else None
        ms=['M:'+n for n,_ in r['subs'] if not is_side(n,extra,keep)]; mp=['M:'+n for n,_ in r['prods'] if not is_side(n,extra,keep)]
        g.add_node(R,kind='rxn',label=label,sub=sub,rev=r['reversible'],pairs=main_pairs(ms,mp))
        for n,_ in r['subs']:
            if not is_side(n,extra,keep):
                g.add_node('M:'+n,kind='met',label=n); g.add_edge('M:'+n,R,side='s'); arrows.append(('M:'+n,R))
        for n,_ in r['prods']:
            if not is_side(n,extra,keep):
                g.add_node('M:'+n,kind='met',label=n)
                sd='b' if g.has_edge('M:'+n,R) and g.edges['M:'+n,R].get('side')=='s' else 'p'
                g.add_edge(R,'M:'+n,side=sd); arrows.append((R,'M:'+n))
    return g,set(arrows)

def _stem(n):
    n=re.sub(r'[^a-z]','',_norm(n).split('[')[0])
    for w in ('bisphosphate','diphosphate','triphosphate','phosphate','phospho'):   # shared chemistry, not identity
        n=n.replace(w,'')
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
def main_pairs(subs,prods):
    """Which substrate->product pairs carry the main flux. With several substrates AND several
    products, pair them by name similarity (succinyl-CoA->succinate, acetate->acetyl-CoA)."""
    pairs=set()
    if len(subs)<=1 or len(prods)<=1: return {(a,b) for a in subs for b in prods}
    for a in subs:
        sims=[(_lcs(_stem(a),_stem(b)),b) for b in prods]; top=max(x for x,_ in sims)
        if top>=3: pairs|={(a,b) for x,b in sims if x==top}
        else: pairs|={(a,b) for b in prods}
    for b in prods:
        if not any(p[1]==b for p in pairs): pairs|={(a,b) for a in subs}
    return pairs

_alias_counter=itertools.count(1)
def detach(g,arrows,m,R):
    """Give edge m-R its own copy of metabolite m (dashed alias box)."""
    a=f"{m}#{next(_alias_counter)}"; g.add_node(a,kind='alias',label=g.nodes[m]['label'])
    sd=g.edges[m,R].get('side'); g.remove_edge(m,R); g.add_edge(a,R,side=sd); g.nodes[m]['aliased']=True
    if (m,R) in arrows: arrows.discard((m,R)); arrows.add((a,R))
    if (R,m) in arrows: arrows.discard((R,m)); arrows.add((R,a))
    return a

# =============================================================== main cycle
def choose_cycle(g,forced=None):
    mets=lambda c:[n for n in c if g.nodes[n]['kind']!='rxn']
    if forced:
        want=['M:'+x.strip() for x in forced]
        if all(w in g for w in want) and len(want)>=2:
            cyc=[]
            for a,b in zip(want,want[1:]+want[:1]):
                rs=sorted(set(g[a])&set(g[b])-set(cyc))
                if not rs: break
                cyc+=[a,rs[0]]
            else: return cyc
    def through(c):          # every reaction on the cycle must turn a substrate into a product
        L=len(c)
        for i,n in enumerate(c):
            if g.nodes[n]['kind']=='rxn':
                a,b=c[i-1],c[(i+1)%L]; P=g.nodes[n]['pairs']
                if (a,b) not in P and (b,a) not in P: return False
        return True
    def steps(ms):           # parallel reactions turning ms[i] into ms[i+1] (either direction)
        out=[]
        for x,y in zip(ms,ms[1:]+ms[:1]):
            rs=sorted(R for R in set(g[x])&set(g[y]) if (x,y) in g.nodes[R]['pairs'] or (y,x) in g.nodes[R]['pairs'])
            out.append(rs)
        return out
    M=nx.Graph()             # metabolite -> metabolite steps carried by a main substrate/product pair
    for R in g:
        if g.nodes[R]['kind']=='rxn' and g.degree(R)<7:          # ignore pseudo-reactions (biomass...)
            for x,y in g.nodes[R]['pairs']:
                if x!=y and x in g and y in g: M.add_edge(x,y)
    cands={}
    for c in itertools.islice(nx.simple_cycles(M,length_bound=min(14,M.number_of_nodes())),200000):
        if len(c)<3: continue
        k=tuple(_canon(c))
        if k in cands: continue
        st=steps(list(k)); used=set(); ok=True
        for i in sorted(range(len(st)),key=lambda i:(len(st[i]),i)):     # one distinct reaction per step
            free=[R for R in st[i] if R not in used]
            if not free: ok=False; break
            used.add(min(free,key=lambda R:(g.degree(R),R)))
        if ok: cands[k]=st
    if not cands: return None
    maxlen=max(len(k) for k in cands)
    best=None
    for k,st in cands.items():
        if len(k)<0.7*maxlen: continue
        # main route = the cycle whose steps are carried by the most parallel reactions
        # (isoenzymes / alternative cofactors); then length, then hub metabolites
        gm=math.exp(sum(math.log(len(r)) for r in st)/len(st))
        key=(round(gm,6),len(k),sum(g.degree(x) for x in k))
        if best is None or key>best[0] or (key==best[0] and k<best[1]): best=(key,k)
    ms=list(best[1]); st=cands[best[1]]; pick={}; used=set()
    for i in sorted(range(len(st)),key=lambda i:(len(st[i]),i)):
        R=min((R for R in st[i] if R not in used),key=lambda R:(g.degree(R),R)); pick[i]=R; used.add(R)
    cyc=[]
    for i,x in enumerate(ms): cyc+=[x,pick[i]]
    return cyc
def _canon(c):
    i=c.index(min(c)); c=c[i:]+c[:i]; r=[c[0]]+c[1:][::-1]
    return min(c,r)

def orient_cycle(g,arrows,cyc):
    pos={n:i for i,n in enumerate(cyc)}; L=len(cyc); fw=0
    for u,v in arrows:
        if u in pos and v in pos:
            d=(pos[v]-pos[u])%L; fw+= 1 if d==1 else (-1 if d==L-1 else 0)
    if fw<0: cyc=[cyc[0]]+cyc[1:][::-1]
    mets=[n for n in cyc if g.nodes[n]['kind']!='rxn']
    top=max(mets,key=lambda m:(g.degree(m),m)); i=cyc.index(top)
    return cyc[i:]+cyc[:i]

# =============================================================== making the graph drawable
def bridges_of(B,cyc):
    C=set(cyc); out=[]
    for comp in sorted(nx.connected_components(B.subgraph([n for n in B if n not in C])),key=lambda s:min(s)):
        att=sorted({m for n in comp for m in B[n] if m in C},key=cyc.index)
        out.append((sorted(comp),att))
    for u,v in sorted(B.edges):                       # chords
        if u in C and v in C and abs(cyc.index(u)-cyc.index(v)) not in (1,len(cyc)-1):
            out.append(([],sorted([u,v],key=cyc.index)))
    return out
def conflict(a1,a2,L,cyc):
    p1=[cyc.index(x) for x in a1]; p2=[cyc.index(x) for x in a2]
    if len(set(p1)&set(p2))>=3: return True
    for x,y in itertools.combinations(sorted(p1),2):
        ins=[p for p in p2 if x<p<y]; outs=[p for p in p2 if p<x or p>y]
        if ins and outs: return True
    return False

def make_drawable(g,arrows,cyc):
    """Duplicate off-cycle metabolites until: the block containing the cycle is planar and every
    other part of the component is a tree hanging from it."""
    C=set(cyc)
    for _ in range(10000):
        core=g.copy()
        while True:                                   # prune pendant trees
            leaves=[n for n in core if core.degree(n)<=1 and n not in C]
            if not leaves: break
            core.remove_nodes_from(leaves)
        blocks=[set(b) for b in nx.biconnected_components(core)]
        main=next(b for b in blocks if C<=b)
        other=[b for b in blocks if b is not main and len(b)>2]
        if other:                                     # open cycles that hang off the main block
            b=min(other,key=lambda s:sorted(s))
            cand=sorted([n for n in b if g.nodes[n]['kind']!='rxn' and n not in C],key=lambda m:(-core.degree(m),m))
            if not cand: cand=sorted([n for n in b if n not in C and g.nodes[n]['kind']!='rxn'])
            m=cand[0]; R=sorted(x for x in core[m] if x in b)[0]; detach(g,arrows,m,R); continue
        B=core.subgraph(main).copy()
        if nx.check_planarity(B)[0]:
            sides=assign_sides(B,cyc)
            if sides is not None: return B,sides
        # not planar: duplicate the off-cycle metabolite whose splitting helps most
        cand=sorted([n for n in B if g.nodes[n]['kind']!='rxn' and n not in C],key=lambda m:(-B.degree(m),m))
        done=False
        for m in cand:
            T=B.copy(); T.remove_node(m)
            if nx.check_planarity(T)[0] or m==cand[-1]:
                rs=sorted(B[m]); keep=[rs[0]]
                for R in rs[1:]:                      # greedy: keep as many edges on the original as stay planar
                    T2=T.copy(); T2.add_edges_from((m,x) for x in keep+[R])
                    if nx.check_planarity(T2)[0]: keep.append(R)
                    else: detach(g,arrows,m,R)
                done=True; break
        if not done:                                  # last resort: detach non-cycle edges of cycle metabolites
            m=max((n for n in C if g.nodes[n]['kind']!='rxn'),key=lambda n:(B.degree(n),n))
            for R in sorted(B[m]):
                if R not in C: detach(g,arrows,m,R)
    raise RuntimeError('could not make the network drawable')

def assign_sides(B,cyc):
    br=bridges_of(B,cyc); L=len(cyc); n=len(br)
    H=nx.Graph(); H.add_nodes_from(range(n))
    for i,j in itertools.combinations(range(n),2):
        if conflict(br[i][1],br[j][1],L,cyc): H.add_edge(i,j)
    if not nx.is_bipartite(H): return None
    pos={c:i for i,c in enumerate(cyc)}
    def span(att):            # smallest arc (in cycle positions) covering all attachments
        p=sorted(pos[x] for x in att)
        if len(p)<2: return 0
        gaps=[(p[(i+1)%len(p)]-p[i])%L for i in range(len(p))]
        return L-max(gaps)
    side={}; crowd={}
    comps=sorted(nx.connected_components(H),key=lambda s:(-max(span(br[i][1]) for i in s),min(s)))
    for comp in comps:
        col=nx.bipartite.color(H.subgraph(comp))
        best=None
        for k in (0,1):       # colour k inside. Outside routes cost span^2 (long arcs);
            c=0.0             # bridges with the same attachments are spread over both sides
            for i in comp:
                sd='in' if col[i]==k else 'out'; key=(tuple(br[i][1]),sd)
                c+=(span(br[i][1])**2)*(len(br[i][0])+1)*(sd=='out') + 6.0*crowd.get(key,0)*(len(br[i][0])+1)
            if best is None or c<best[0]: best=(c,k)
        for i in comp:
            side[i]='in' if col[i]==best[1] else 'out'; key=(tuple(br[i][1]),side[i]); crowd[key]=crowd.get(key,0)+1
    return [(br[i][0],br[i][1],side[i]) for i in range(n)]

# =============================================================== geometry
FP=lambda st: FontProperties(family='DejaVu Sans',size=st['fs'],weight=st['weight'],style=st['style'])
STYLE={'core': dict(fs=9.5,weight='bold',style='normal',pad=(7,5),fc='#fde3e1',ec='#d9534f',tc='#3a0d0b'),
       'met':  dict(fs=8.5,weight='normal',style='normal',pad=(6,4),fc='#e3ecfd',ec='#5b8def',tc='#0f2550'),
       'alias':dict(fs=8.5,weight='normal',style='normal',pad=(6,4),fc='#e3ecfd',ec='#5b8def',tc='#0f2550',ls='--'),
       'rxn':  dict(fs=7.0,weight='normal',style='italic',pad=(4,3),fc='#ffffff',ec='#2b2b2b',tc='#222222')}
SUB=dict(fs=5.8,weight='normal',style='normal',pad=(4,0),tc='#777777')
_tw={}
def text_w(t,st):
    k=(t,st['fs'],st['weight'],st['style'])
    if k not in _tw: _tw[k]=TextPath((0,0),t,prop=FP(st)).get_extents().width
    return _tw[k]
def node_size(g,n):
    st=STYLE[g.nodes[n]['kind']]; w=text_w(g.nodes[n]['label'],st)+2*st['pad'][0]; h=st['fs']+2*st['pad'][1]
    if g.nodes[n].get('sub'): w=max(w,text_w(g.nodes[n]['sub'],SUB)+2*SUB['pad'][0]); h+=SUB['fs']*1.25
    return w,h

def tutte(H,fixed,outer_nodes=None):
    """Barycentric layout of plane graph H (C fixed) after stellating every face except the one
    bounded by the fixed cycle: strictly planar for any 2-connected H."""
    C=list(fixed); fixed=dict(fixed)
    K=H.copy(); K.add_node('__OUT__')
    for i,u in enumerate(C):                 # subdivide each boundary edge: an outside face can then
        v=C[(i+1)%len(C)]; x=f'__X{i}__'      # touch only one cycle node, so every route stays inside
        if K.has_edge(u,v): K.remove_edge(u,v)
        K.add_edges_from([(u,x),(x,v),('__OUT__',x),('__OUT__',u)])
        fixed[x]=((fixed[u][0]+fixed[v][0])/2,(fixed[u][1]+fixed[v][1])/2)
    H=K.copy(); H.remove_node('__OUT__'); C=list(fixed); Cs=set(C)
    ok,emb=nx.check_planarity(K)
    if not ok: return None
    seen=set(); faces=[]
    for u,v in sorted(emb.edges()):
        if (u,v) in seen: continue
        faces.append(emb.traverse_face(u,v,mark_half_edges=seen))
    T=H.copy(); d=0
    for f in faces:
        if '__OUT__' in f:
            if any(x not in Cs for x in f if x!='__OUT__'): return None
            continue
        D=f'__D{d}__'; d+=1; T.add_edges_from((D,x) for x in sorted(set(f)))
    free=sorted(n for n in T if n not in Cs); ix={n:i for i,n in enumerate(free)}
    if not free: return {c:fixed[c] for c in C}
    A=np.zeros((len(free),)*2); b=np.zeros((len(free),2))
    for n in free:
        i=ix[n]
        for m in T[n]:
            A[i,i]+=1
            if m in ix: A[i,ix[m]]-=1
            else: b[i]+=fixed[m]
    X=np.linalg.solve(A,b)
    pos={c:np.array(fixed[c],float) for c in C}
    pos.update({n:X[ix[n]] for n in free if not str(n).startswith('__D')})
    pos={n:p for n,p in pos.items() if not str(n).startswith('__X')}
    return pos

class Drawing:
    """Node positions + edge routing. Edge kinds: 'ring' (arc of the circle), 'out' (straight in the
    reflected disc -> curve outside the circle), 'line' (straight)."""
    def __init__(s,g,arrows,cyc,R,kref):
        s.g=g; s.arrows=arrows; s.cyc=cyc; s.R=R; s.k=kref; s.P={}; s.ekind={}; s.size={n:node_size(g,n) for n in g}
    def refl(s,p):        # disc coordinates -> outside coordinates
        r=math.hypot(*p); 
        if r<1e-9: r=1e-9
        rr=s.R+s.k*(s.R-r); return (p[0]/r*rr,p[1]/r*rr)
    def unrefl(s,p):
        r=math.hypot(*p); rr=s.R-(r-s.R)/s.k; return (p[0]/r*rr,p[1]/r*rr)
    def path(s,u,v):
        a,b=s.P[u],s.P[v]; k=s.ekind[(u,v)]
        if k=='line': return [a,b]
        if k=='ring':
            t1,t2=math.atan2(a[1],a[0]),math.atan2(b[1],b[0]); d=(t2-t1+math.pi)%(2*math.pi)-math.pi
            n=max(4,int(abs(d)*40)); r1,r2=math.hypot(*a),math.hypot(*b)
            return [((r1+(r2-r1)*i/n)*math.cos(t1+d*i/n),(r1+(r2-r1)*i/n)*math.sin(t1+d*i/n)) for i in range(n+1)]
        A,Bp=s.unrefl(a),s.unrefl(b); n=24
        return [s.refl((A[0]+(Bp[0]-A[0])*i/n,A[1]+(Bp[1]-A[1])*i/n)) for i in range(n+1)]
    def box(s,n,P=None):
        (x,y)=(P or s.P)[n]; w,h=s.size[n]; return (x-w/2,y-h/2,x+w/2,y+h/2)

def layout_component(g,arrows,nodes,forced_cycle):
    import time; T0=time.time(); log=lambda m: print(f'  [{time.time()-T0:6.1f}s] {m}',file=sys.stderr,flush=True)
    sub=g.subgraph(nodes).copy()
    cyc=choose_cycle(sub,forced_cycle)
    if cyc is None: return layout_tree(g,arrows,nodes)
    cyc=orient_cycle(g,arrows,cyc); log('cycle chosen')
    # graph surgery happens on g itself (aliases are added to g)
    comp_nodes=set(nodes)
    sub=g.subgraph(comp_nodes).copy()
    before=set(g.nodes)
    B,sides=make_drawable_on(g,arrows,cyc,comp_nodes)
    comp_nodes|=set(g.nodes)-before; log('drawable')
    for n in cyc:
        if g.nodes[n]['kind']=='met': g.nodes[n]['kind']='core'
    # circle radius from the boxes on the cycle
    sizes=[node_size(g,n) for n in cyc]; per=[max(w,h*2.2)+26 for w,h in sizes]
    inner=sum((node_size(g,n)[0]+30)*(node_size(g,n)[1]+30) for n in B if n not in cyc)
    R=max(150.0,sum(per)/(2*math.pi),math.sqrt(2.2*inner/math.pi)); D=Drawing(g,arrows,cyc,R,0.6)
    ang=[]; acc=0
    for p in per: ang.append(acc+p/2); acc+=p
    fixed={n:(R*math.sin(2*math.pi*a/acc),R*math.cos(2*math.pi*a/acc)) for n,a in zip(cyc,ang)}
    for side in ('in','out'):
        nodes_s=[n for br,att,sd in sides if sd==side for n in br]
        H=B.subgraph(list(cyc)+nodes_s).copy()
        H.add_edges_from(zip(cyc,cyc[1:]+cyc[:1]))
        pos=tutte(H,fixed)
        if pos is None: raise RuntimeError('embedding failed')
        for n in nodes_s: D.P[n]=tuple(pos[n]) if side=='in' else D.refl(tuple(pos[n]))
        for u,v in H.edges:
            if u in cyc and v in cyc and abs(cyc.index(u)-cyc.index(v)) in (1,len(cyc)-1): continue
            k='line' if side=='in' else 'out'
            D.ekind[(u,v)]=D.ekind[(v,u)]=k
    for n in cyc: D.P[n]=fixed[n]
    for u,v in zip(cyc,cyc[1:]+cyc[:1]): D.ekind[(u,v)]=D.ekind[(v,u)]='ring'
    fan_out(D,B,cyc); log('tutte done')
    plant_trees(D,g,comp_nodes,set(B.nodes)); log('trees planted')
    return D,comp_nodes

def drawable(B,cyc):
    if not nx.check_planarity(B)[0]: return None
    sides=assign_sides(B,cyc)
    if sides is None: return None
    for sd in ('in','out'):
        H=B.subgraph(list(cyc)+[n for br,att,s in sides if s==sd for n in br]).copy()
        H.add_edges_from(zip(cyc,cyc[1:]+cyc[:1]))
        if tutte(H,{c:(math.cos(2*math.pi*i/len(cyc)),math.sin(2*math.pi*i/len(cyc))) for i,c in enumerate(cyc)}) is None: return None
    return sides

def make_drawable_on(g,arrows,cyc,comp_nodes):
    sub_before=set(g.nodes)
    # operate on the component view; detach() mutates g
    class View: pass
    gg=g.subgraph(comp_nodes)                    # read-only view is refreshed inside the loop
    C=set(cyc)
    for R in sorted(n for n in comp_nodes if g.nodes[n]['kind']=='rxn'):      # pseudo-reactions such as biomass
        mets=sorted(g[R])
        if len(mets)>=7:
            for m in mets:
                if m not in C or R not in C: detach(g,arrows,m,R)
    for u,v in sorted(g.subgraph(cyc).edges):                                # chords of the cycle
        if abs(cyc.index(u)-cyc.index(v)) not in (1,len(cyc)-1):
            m,R=(u,v) if g.nodes[v]['kind']=='rxn' else (v,u); detach(g,arrows,m,R)
    for _ in range(10000):
        nodes=comp_nodes|(set(g.nodes)-sub_before)
        core=g.subgraph(nodes).copy()
        while True:
            leaves=[n for n in core if core.degree(n)<=1 and n not in C]
            if not leaves: break
            core.remove_nodes_from(leaves)
        blocks=[set(b) for b in nx.biconnected_components(core)]
        main=next(b for b in blocks if C<=b)
        other=sorted([b for b in blocks if b is not main and len(b)>2],key=lambda s:sorted(s))
        if other:
            b=other[0]
            cand=sorted([n for n in b if g.nodes[n]['kind']!='rxn' and n not in C],key=lambda m:(-core.degree(m),m))
            m=cand[0]; R=sorted(x for x in core[m] if x in b)[0]; detach(g,arrows,m,R); continue
        B=core.subgraph(main).copy()
        sides=drawable(B,cyc)
        if sides is not None: return B,sides
        cand=sorted([n for n in B if g.nodes[n]['kind']!='rxn' and n not in C],key=lambda m:(-B.degree(m),m))
        progress=False
        for m in cand:
            T=B.copy(); T.remove_node(m)
            if drawable(T,cyc) is not None or m==cand[-1]:
                rs=sorted(B[m]); keep=[rs[0]]
                for R in rs[1:]:
                    T2=T.copy(); T2.add_edges_from((m,x) for x in keep+[R])
                    if drawable(T2,cyc) is not None: keep.append(R)
                    else: detach(g,arrows,m,R); progress=True
                if progress: break
        if not progress:             # last resort: detach the off-cycle edges of a cycle metabolite
            for m in sorted((n for n in C if g.nodes[n]['kind']!='rxn'),key=lambda n:(-B.degree(n),n)):
                extra=[R for R in sorted(B[m]) if R not in C]
                if extra:
                    for R in extra: detach(g,arrows,m,R)
                    break
            else: raise RuntimeError('could not make the network drawable')
    raise RuntimeError('could not make the network drawable')

def fan_out(D,B,cyc):
    """Tutte can put two nodes with identical neighbourhoods (e.g. isoenzymes for the same step) on
    one point: spread such groups perpendicular to their axis."""
    groups={}
    for n in sorted(D.P):
        if n in cyc: continue
        groups.setdefault((round(D.P[n][0],4),round(D.P[n][1],4)),[]).append(n)
    for key,grp in sorted(groups.items()):
        if len(grp)<2: continue
        nb=sorted(B[grp[0]]); a=np.array(D.P[nb[0]]); b=np.array(D.P[nb[-1]]) if len(nb)>1 else a+1
        d=b-a; nrm=np.array([-d[1],d[0]])/max(np.linalg.norm(d),1e-9)
        ext=max(abs(nrm[0])*D.size[n][0]+abs(nrm[1])*D.size[n][1] for n in grp)+24   # box extent across the axis
        for j,n in enumerate(grp):
            o=(j-(len(grp)-1)/2)*ext
            D.P[n]=(key[0]+nrm[0]*o,key[1]+nrm[1]*o)

def plant_trees(D,g,comp_nodes,placed):
    """Attach the remaining (tree) nodes in the widest free angular gap of their parent."""
    rest=[n for n in sorted(comp_nodes) if n not in placed]
    T=g.subgraph(list(placed)+rest)
    order=[]; seen=set(placed)
    frontier=sorted(placed)
    while True:
        new=[]
        for p in frontier:
            kids=sorted(m for m in T[p] if m not in seen)
            if kids: order.append((p,kids)); seen.update(kids); new+=kids
        if not new: break
        frontier=new
    D.children=getattr(D,'children',{})
    for p,kids in order:
        D.children[p]=list(kids)
        for n in kids: D.ekind[(p,n)]=D.ekind[(n,p)]='line'
        x,y=D.P[p]
        def dirs():
            out=[]
            for m in g[p]:
                if m in D.P and m not in kids and (p,m) in D.ekind:
                    P=D.path(p,m); q=P[1] if len(P)>1 else P[-1]
                    out.append(math.atan2(q[1]-y,q[0]-x))
            return sorted(out)
        A=dirs()
        if not A: cands=[(2*math.pi,math.atan2(y,x)-math.pi)]
        else: cands=sorted(((((A[(i+1)%len(A)]-A[i])%(2*math.pi)) or 2*math.pi,A[i]) for i in range(len(A))),key=lambda t:(-t[0],t[1]))
        best=None
        for gap,lo in cands[:4]:
            for f in (1.0,0.7,0.5,0.35,0.25,0.15,0.1,0.06):
                L=f*(90+8*len(kids))
                for j,n in enumerate(kids):
                    t=lo+gap*(j+1)/(len(kids)+1); D.P[n]=(x+L*math.cos(t),y+L*math.sin(t))
                if crossings(D,only=[(p,n) for n in kids])==0:
                    v=violations(D); key=(v['line_through_box']+v['box_overlaps'],-f,-gap)
                    if best is None or key<best[0]: best=(key,{n:D.P[n] for n in kids})
                    break
        if best: D.P.update(best[1])
        else:
            for j,n in enumerate(kids): D.P[n]=(x+3*math.cos(j),y+3*math.sin(j))

def layout_tree(g,arrows,nodes):
    """Component without a cycle: radial tree (wedges proportional to subtree size) - planar."""
    sub=g.subgraph(nodes).copy()
    while not nx.is_forest(sub):
        cyc=nx.find_cycle(sub); m=sorted(n for e in cyc for n in e[:2] if g.nodes[n]['kind']!='rxn')[0]
        R=sorted(x for x in sub[m] if any(set(e[:2])=={m,x} for e in cyc))[0]
        a=detach(g,arrows,m,R); sub=g.subgraph(set(nodes)|{a}).copy(); nodes=list(sub.nodes)
    root=max(sorted(sub),key=lambda n:sub.degree(n)); D=Drawing(g,arrows,[],100.0,1.0)
    tree=nx.bfs_tree(sub,root); depth=nx.shortest_path_length(tree,root)
    wid={n:node_size(g,n)[0]+24 for n in tree}; size={}
    for n in reversed(list(nx.topological_sort(tree))): size[n]=max(wid[n],sum(size[c] for c in tree[n]))
    rad=[0.0]
    for d in range(1,max(depth.values())+1):
        lvl=sum(wid[n] for n in tree if depth[n]==d)
        rad.append(max(rad[-1]+70,lvl/(2*math.pi)*1.3,size[root]/(2*math.pi)*(d/max(depth.values()))))
    def place(n,a0,a1):
        t=(a0+a1)/2; r=rad[depth[n]]; D.P[n]=(r*math.sin(t),r*math.cos(t))
        kids=sorted(tree[n]); tot=sum(size[c] for c in kids); a=a0
        for c in kids:
            w=(a1-a0)*size[c]/tot; D.ekind[(n,c)]=D.ekind[(c,n)]='line'; place(c,a,a+w); a+=w
    place(root,0,2*math.pi)
    return D,set(sub.nodes)

# =============================================================== checks
def segments(D,edges=None):
    S=[]
    for u,v in (edges if edges is not None else D.edges):
        P=D.path(u,v)
        for a,b in zip(P,P[1:]): S.append((u,v,a,b))
    return S
def _segx(A,B,C,Dd):
    o=lambda p,q,r:(q[...,0]-p[...,0])*(r[...,1]-p[...,1])-(q[...,1]-p[...,1])*(r[...,0]-p[...,0])
    return (o(A,B,C)*o(A,B,Dd)<-1e-9)&(o(C,Dd,A)*o(C,Dd,B)<-1e-9)
def crossings(D,only=None):
    E=[e for e in D.g.edges if e in D.ekind or (e[1],e[0]) in D.ekind]
    D.edges=[(u,v) for u,v in E if u in D.P and v in D.P]
    allS=segments(D); 
    if not allS: return 0
    U=[s[0] for s in allS]; V=[s[1] for s in allS]
    A=np.array([s[2] for s in allS]); B=np.array([s[3] for s in allS])
    qs=segments(D,only) if only is not None else allS
    n=0
    for u,v,a,b in qs:
        mask=np.array([not ({u,v}&{x,y}) for x,y in zip(U,V)])
        if not mask.any(): continue
        n+=int(_segx(np.array(a)[None],np.array(b)[None],A[mask],B[mask]).sum())
    return n if only is not None else n//2
def violations(D,margin=3):
    crossings(D); S=segments(D); out={'crossings':crossings(D),'line_through_box':0,'box_overlaps':0}
    names=sorted(D.P); Bx=np.array([D.box(n) for n in names])
    for u,v,a,b in S:
        ts=np.linspace(0,1,12); X=a[0]+(b[0]-a[0])*ts; Y=a[1]+(b[1]-a[1])*ts
        hit=(X[None]>Bx[:,0,None]-margin)&(X[None]<Bx[:,2,None]+margin)&(Y[None]>Bx[:,1,None]-margin)&(Y[None]<Bx[:,3,None]+margin)
        for i in np.where(hit.any(1))[0]:
            if names[i] not in (u,v): out['line_through_box']+=1
    ox=(Bx[:,None,0]<Bx[None,:,2]+margin)&(Bx[None,:,0]<Bx[:,None,2]+margin)&(Bx[:,None,1]<Bx[None,:,3]+margin)&(Bx[None,:,1]<Bx[:,None,3]+margin)
    out['box_overlaps']=int(np.triu(ox,1).sum())
    return out

# =============================================================== spreading (crossing-guarded)
def cost(D,P,names,fixedset):
    Bx=np.array([D.box(n,P) for n in names]); m=6
    ox=np.maximum(0,np.minimum(Bx[:,None,2],Bx[None,:,2])-np.maximum(Bx[:,None,0],Bx[None,:,0])+m)
    oy=np.maximum(0,np.minimum(Bx[:,None,3],Bx[None,:,3])-np.maximum(Bx[:,None,1],Bx[None,:,1])+m)
    c=0.05*np.triu(ox*oy,1).sum()
    old=D.P; D.P=P
    try:
        S=segments(D); L=[]
        for u,v in D.edges:
            p=D.path(u,v); L.append(sum(math.dist(a,b) for a,b in zip(p,p[1:])))
        pts=[]; own=[]
        for u,v,a,b in S:
            for t in (0.15,0.5,0.85): pts.append((a[0]+(b[0]-a[0])*t,a[1]+(b[1]-a[1])*t)); own.append((u,v))
        pts=np.array(pts)
        px=np.maximum(0,np.minimum(pts[None,:,0]-Bx[:,None,0],Bx[:,None,2]-pts[None,:,0])+m)
        py=np.maximum(0,np.minimum(pts[None,:,1]-Bx[:,None,1],Bx[:,None,3]-pts[None,:,1])+m)
        mask=np.array([[names[i] not in own[j] for j in range(len(own))] for i in range(len(names))])
        c+=0.5*(px*py*mask).sum()
        L=np.array(L); med=np.median(L); c+=0.002*((L-med)**2).sum()+0.02*np.maximum(L-2.5*med,0).sum()
    finally: D.P=old
    return c

def spread(D,fixedset,sweeps=60):
    names=sorted(D.P); free=[n for n in names if n not in fixedset]
    crossings(D); cur=cost(D,D.P,names,fixedset); step=40.0
    dirs=[(math.cos(k*math.pi/4),math.sin(k*math.pi/4)) for k in range(8)]
    inc={n:[(u,v) for u,v in D.edges if n in (u,v)] for n in names}
    for sweep in range(sweeps):
        improved=False
        for n in free:
            best=None
            for dx,dy in dirs:
                P=dict(D.P); P[n]=(D.P[n][0]+dx*step,D.P[n][1]+dy*step)
                if any(D.ekind.get(e)=='out' for e in inc[n]) and math.hypot(*P[n])<=D.R+5: continue
                if any(D.ekind.get(e)=='out' for e in inc[n]) and math.hypot(*P[n])>=D.R*(1+D.k)-5: continue
                c=cost(D,P,names,fixedset)
                if c<cur-1e-6 and (best is None or c<best[0]): best=(c,P)
            if best:
                old=D.P; D.P=best[1]
                ok=all(crossings(D,only=inc[n])==0 for _ in [0])
                mid=dict(old); mid[n]=((old[n][0]+best[1][n][0])/2,(old[n][1]+best[1][n][1])/2); D.P=mid
                ok=ok and crossings(D,only=inc[n])==0
                D.P=best[1] if ok else old
                if ok: cur=best[0]; improved=True
        if not improved:
            step*=0.5
            if step<3: break

def offending_nodes(D,margin=3):
    S=segments(D); names=sorted(D.P); Bx=np.array([D.box(n) for n in names]); out=[]
    for u,v,a,b in S:
        ts=np.linspace(0,1,12); X=a[0]+(b[0]-a[0])*ts; Y=a[1]+(b[1]-a[1])*ts
        hit=(X[None]>Bx[:,0,None]-margin)&(X[None]<Bx[:,2,None]+margin)&(Y[None]>Bx[:,1,None]-margin)&(Y[None]<Bx[:,3,None]+margin)
        for i in np.where(hit.any(1))[0]:
            if names[i] not in (u,v): out+=[names[i],u,v]
    for i in range(len(names)):
        for j in range(i+1,len(names)):
            A,B=Bx[i],Bx[j]
            if A[0]<B[2]+margin and B[0]<A[2]+margin and A[1]<B[3]+margin and B[1]<A[3]+margin: out+=[names[i],names[j]]
    seen=[]; [seen.append(x) for x in out if x not in seen]
    return seen
def repair(D,fixedset,rounds=30):
    """Move single nodes, or a planted subtree as a whole, never creating a crossing, until no box
    overlaps another or is crossed by a line. A move is kept if it lowers (violations, overlap depth)."""
    names=sorted(D.P); ch=getattr(D,'children',{})
    def subtree(n):
        out=[n]; i=0
        while i<len(out): out+=ch.get(out[i],[]); i+=1
        return out
    def score():
        v=violations(D); return (v['line_through_box']+v['box_overlaps'],round(cost(D,D.P,names,fixedset),3))
    for _ in range(rounds):
        cur=score()
        if cur[0]==0: return
        moved=False
        for n in [x for x in offending_nodes(D) if x not in fixedset]:
            par=[p for p,k in ch.items() if n in k and p not in fixedset]
            groups=[[n]]+([subtree(n)] if ch.get(n) else [])+[subtree(p) for p in par]
            for group in groups:
                inc=[(u,v) for u,v in D.edges if (u in group)!=(v in group) or (u in group and v in group)]
                best=None; base={m:D.P[m] for m in group}
                for dist in (25,50,100,180,300,500):
                    for k in range(16):
                        t=k*math.pi/8; dx,dy=dist*math.cos(t),dist*math.sin(t)
                        okk=True
                        for f in (0.5,1.0):
                            for m in group: D.P[m]=(base[m][0]+dx*f,base[m][1]+dy*f)
                            if any(D.ekind.get(e)=='out' and not (D.R+5<math.hypot(*D.P[m])<D.R*(1+D.k)-5)
                                   for m in group for e in D.edges if m in e) or crossings(D,only=inc)!=0: okk=False; break
                        c=score() if okk else None
                        for m in group: D.P[m]=base[m]
                        if okk and c<cur and (best is None or c<best[0]): best=(c,dx,dy)
                    if best and best[0][0]<cur[0]: break
                if best:
                    for m in group: D.P[m]=(base[m][0]+best[1],base[m][1]+best[2])
                    cur=best[0]; moved=True; break
            if cur[0]==0: return
        if not moved: return

def scale_until_clean(D,maxit=40):
    for _ in range(maxit):
        v=violations(D)
        if v['line_through_box']==0 and v['box_overlaps']==0: return v
        f=1.12; D.P={n:(p[0]*f,p[1]*f) for n,p in D.P.items()}; D.R*=f
    return violations(D)

# =============================================================== rendering
def render(parts,path,title):
    fig_items=[]; x_off=0; boxes=[]
    allb=[]
    for D,off in parts:
        for n in D.P: b=D.box(n); allb.append((b[0]+off[0],b[1]+off[1],b[2]+off[0],b[3]+off[1]))
    xs=[b[0] for b in allb]+[b[2] for b in allb]; ys=[b[1] for b in allb]+[b[3] for b in allb]
    for D,off in parts:
        for u,v in D.edges:
            for q in D.path(u,v): xs.append(q[0]+off[0]); ys.append(q[1]+off[1])
    x0,x1,y0,y1=min(xs)-30,max(xs)+30,min(ys)-30,max(ys)+50
    fig=plt.figure(figsize=((x1-x0)/72,(y1-y0)/72)); ax=fig.add_axes([0,0,1,1])
    ax.set_xlim(x0,x1); ax.set_ylim(y0,y1); ax.set_aspect('equal'); ax.axis('off')
    for D,(ox,oy) in parts:
        if D.cyc: ax.add_patch(plt.Circle((ox,oy),D.R,fc='#fbf6f5',ec='none',zorder=0))
        for u,v in D.edges:
            src,dst=(u,v) if (u,v) in D.arrows else (v,u)
            P=D.path(src,dst)
            def cut(P,b):
                inb=lambda p:b[0]<p[0]<b[2] and b[1]<p[1]<b[3]
                i=0
                while i<len(P)-1 and inb(P[i+1]): i+=1
                if i>=len(P)-1: return P[-2:]
                a,c=P[i],P[i+1]; lo,hi=0.0,1.0
                for _ in range(30):
                    m=(lo+hi)/2; q=(a[0]+(c[0]-a[0])*m,a[1]+(c[1]-a[1])*m); (lo,hi)=(m,hi) if inb(q) else (lo,m)
                return [(a[0]+(c[0]-a[0])*hi,a[1]+(c[1]-a[1])*hi)]+P[i+1:]
            P=cut(P,D.box(src)); P=cut(P[::-1],D.box(dst))[::-1]
            P=[(p[0]+ox,p[1]+oy) for p in P]
            ring=D.ekind[(u,v)]=='ring'; col,lw=('#b23b36',1.8) if ring else ('#4a4a4a',1.1)
            ax.plot([p[0] for p in P[:-1]],[p[1] for p in P[:-1]],color=col,lw=lw,solid_capstyle='butt',zorder=1)
            ax.add_patch(FancyArrowPatch(P[-2],P[-1],arrowstyle='-|>',mutation_scale=9 if ring else 8,color=col,lw=lw,shrinkA=0,shrinkB=0,zorder=1))
        for n in sorted(D.P):
            k=D.g.nodes[n]['kind']; st=STYLE[k]; x,y=D.P[n][0]+ox,D.P[n][1]+oy; w,h=D.size[n]
            ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=0,rounding_size=%g'%(h/2 if k!='rxn' else 2),
                         fc=st['fc'],ec=st['ec'],lw=0.9,ls=st.get('ls','-'),zorder=2))
            sub=D.g.nodes[n].get('sub'); dy=SUB['fs']*0.62 if sub else 0
            ax.text(x,y+dy,D.g.nodes[n]['label'],ha='center',va='center_baseline',fontsize=st['fs'],weight=st['weight'],style=st['style'],color=st['tc'],zorder=3)
            if sub: ax.text(x,y-st['fs']*0.62,sub,ha='center',va='center_baseline',fontsize=SUB['fs'],color=SUB['tc'],zorder=3)
    ax.text(x0+12,y1-14,title,fontsize=12,va='top')
    fig.savefig(path,**({'metadata':{'Date':None}} if path.endswith('.svg') else {'dpi':min(110.0,7000*72/max(x1-x0,y1-y0))})); plt.close(fig)

# =============================================================== main
def main():
    ap=argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('sbml'); ap.add_argument('--out',default=None)
    ap.add_argument('--cycle',help='comma-separated metabolite names of the main cycle, in order')
    ap.add_argument('--side',default='',help='extra side-compound names'); ap.add_argument('--keep',default='',help='names never treated as side compounds')
    a=ap.parse_args()
    name,reactions=read_sbml(a.sbml)
    g,arrows=build_graph(reactions,[x for x in a.side.split(',') if x],[x for x in a.keep.split(',') if x])
    out=a.out or os.path.splitext(os.path.basename(a.sbml))[0]+'_layout'
    comps=sorted((sorted(c) for c in nx.connected_components(g.copy())),key=lambda c:(-len(c),c))
    parts=[]
    for c in comps:
        forced=a.cycle.split(',') if a.cycle else None
        D,nodes=layout_component(g,arrows,c,forced)
        D.g=g; crossings(D)
        if D.cyc:
            big=len(D.P)>80   # large networks: fewer optimisation passes (time); scaling still guarantees the rules
            print('  spreading',file=sys.stderr,flush=True); spread(D,set(D.cyc),sweeps=4 if big else 60); print('  repairing',file=sys.stderr,flush=True); repair(D,set(D.cyc),rounds=3 if big else 30)
        v=scale_until_clean(D)
        print(f'component with {len(nodes)} nodes: cycle of {len([n for n in D.cyc if g.nodes[n]["kind"]!="rxn"])} metabolites,',v)
        if any(v.values()): sys.exit('rules could not be met - no output written')
        parts.append(D)
    # pack components left to right (largest first), wrapping rows
    placed=[]; x=0; y=0; rowh=0; maxw=0
    widths=[]
    for D in parts:
        b=[D.box(n) for n in D.P]; widths.append((min(q[0] for q in b),min(q[1] for q in b),max(q[2] for q in b),max(q[3] for q in b)))
    limit=max(1500,max(w[2]-w[0] for w in widths))
    for D,(bx0,by0,bx1,by1) in zip(parts,widths):
        if x>0 and x+(bx1-bx0)>limit: x=0; y-=rowh+60; rowh=0
        placed.append((D,(x-bx0,y-by1))); x+=bx1-bx0+60; rowh=max(rowh,by1-by0)
    nal=sum(1 for n in g if g.nodes[n]['kind']=='alias')
    print(f'{len(reactions)} reactions, {sum(1 for n in g if g.nodes[n]["kind"]!="rxn")} metabolite boxes ({nal} duplicates), '
          f'{sum(len(D.edges) for D in parts)} lines, 0 bends')
    render(placed,out+'.svg',name); render(placed,out+'.png',name)
    print('written',out+'.svg',out+'.png')

if __name__=='__main__': main()

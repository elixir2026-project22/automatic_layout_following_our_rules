#!/usr/bin/env python3
"""Quality metrics for a pathway/graph drawing in SVG, measured against the layout rules:

The rule set, definitions, thresholds and 0-100 scores are those of core.js (Biohackathon 2026 layout scorer), with the
cycle rule taken from these Python tools (closed loop in cycle order) and the rules that only these tools had added:
  hard  Reproducible rendering (--compare) | Labels at fixed positions per type | Cycles look like cycles (closed loop)
        Cofactors drawn as clones | No crossings | No line through a node | Consistent participant placement | Aspect ratio
  soft  Avoid long lines | Avoid bends | Density | No overlapping nodes

The SVG carries no notion of "metabolite" or "reaction", so the drawing is reconstructed from its
geometry (works for any SVG made of paths: matplotlib, Inkscape, Escher/Cytoscape exports, ...):
  * boxes   = closed shapes with an outline (node shapes); the outline/fill/dash style is the node TYPE
  * labels  = text: <text> elements, or glyph outlines (small filled shapes without outline) grouped
              into words/lines
  * lines   = open stroked paths; pieces that meet end to end with the same colour (e.g. a polyline
              plus its arrow) are joined into one connection; arrowheads are recognised and ignored
  * backgrounds (large unstroked fills) are ignored
All lengths are in SVG user units (px).

Usage:
  python3 svg_layout_metrics.py drawing.svg [--compare other.svg] [--json out.json] [--annotate out.png]
--compare checks reproducibility against a second rendering (e.g. a rerun of the same command).
--annotate draws the reconstructed geometry and marks every problem found.
Requires: svgelements, numpy (and matplotlib for --annotate)
"""
import sys, re, io, math, json, argparse, hashlib, itertools
import numpy as np
from svgelements import SVG, Path, Shape, Text, Close, Line, Move

# ------------------------------------------------------------------ reading
def _clean(src):
    src=re.sub(r'<path(?![^>]*\sd=)[^>]*/>','',src)            # outline-less glyphs (e.g. space) / empty paths
    return src
def _has(paint): return paint is not None and getattr(paint,'value',None) is not None
def _poly(path,step=4.0):
    """Flatten a path into sub-polylines [(points, closed)]; curves are sampled (keeps 'curved' flag)."""
    out=[]; pts=[]; curved=False
    for seg in path:
        if isinstance(seg,Move):
            if len(pts)>1: out.append((pts,False,curved))
            pts=[(seg.end.x,seg.end.y)]; curved=False
        elif isinstance(seg,Close):
            if pts and (pts[0]!=pts[-1]): pts.append(pts[0])
            if len(pts)>1: out.append((pts,True,curved))
            pts=[]; curved=False
        elif isinstance(seg,Line):
            if not pts: pts=[(seg.start.x,seg.start.y)]
            pts.append((seg.end.x,seg.end.y))
        else:
            if not pts: pts=[(seg.start.x,seg.start.y)]
            n=max(2,int(seg.length(error=1e-2)/step)) if hasattr(seg,'length') else 8
            for i in range(1,n+1):
                p=seg.point(i/n); pts.append((p.x,p.y))
            curved=True
    if len(pts)>1: out.append((pts,False,curved))
    return out
def _area(pts):
    a=0.0
    for (x1,y1),(x2,y2) in zip(pts,pts[1:]): a+=x1*y2-x2*y1
    return abs(a)/2
_cmaps={}
def _glyph_char(font,index):
    """matplotlib names glyph outlines <font>-<glyph index>; map the index back to its character."""
    if font not in _cmaps:
        inv={}
        try:
            from matplotlib import font_manager, ft2font
            fam,_,style=font.partition('-')
            prop=font_manager.FontProperties(family=re.sub(r'(Sans|Serif|Mono)',r' \1',fam).strip(),
                                             style='oblique' if 'Oblique' in style or 'Italic' in style else 'normal',
                                             weight='bold' if 'Bold' in style else 'normal')
            f=ft2font.FT2Font(font_manager.findfont(prop,fallback_to_default=True))
            inv={g:chr(c) for c,g in sorted(f.get_charmap().items(),reverse=True)}
        except Exception: pass
        _cmaps[font]=inv
    return _cmaps[font].get(index,'?')
def read(svgfile):
    src=_clean(open(svgfile,encoding='utf-8').read())
    svg=SVG.parse(io.StringIO(src))
    shapes=[]; texts=[]
    for eid,e in enumerate(svg.elements()):
        if isinstance(e,Text):
            if not (e.text or '').strip(): continue
            fs=float(e.font_size or 12); w=0.55*fs*len(e.text); x=e.x; y=e.y
            anchor=(e.values.get('text-anchor') or 'start')
            x0=x-w/2 if anchor=='middle' else (x-w if anchor=='end' else x)
            p=e.transform.point_in_matrix_space if hasattr(e.transform,'point_in_matrix_space') else None
            texts.append(dict(bbox=(x0,y-fs*0.8,x0+w,y+fs*0.2),text=e.text,color=str(e.fill),italic='italic' in str(e.values.get('font-style',''))))
            continue
        if not isinstance(e,Shape): continue
        p=Path(e)
        for pts,closed,curved in _poly(p):
            xs=[q[0] for q in pts]; ys=[q[1] for q in pts]
            gm=re.match(r'^(.+)-([0-9a-f]{1,5})$',str(e.id or ''))
            shapes.append(dict(eid=eid,ch=_glyph_char(gm.group(1),int(gm.group(2),16)) if gm else '',italic=bool(gm and re.search('Oblique|Italic',gm.group(1))),
                               pts=pts,closed=closed,curved=curved,bbox=(min(xs),min(ys),max(xs),max(ys)),
                               stroke=str(e.stroke) if _has(e.stroke) else None,fill=str(e.fill) if _has(e.fill) else None,
                               sw=float(e.stroke_width or 0) if _has(e.stroke) else 0.0,
                               dash=str(e.values.get('stroke-dasharray') or 'solid'),area=_area(pts) if closed else 0.0))
    vb=svg.viewbox
    canvas=(0,0,float(svg.width),float(svg.height))
    return src,shapes,texts,canvas

# ------------------------------------------------------------------ reconstruction
def classify(shapes,texts,canvas):
    W,H=canvas[2]-canvas[0],canvas[3]-canvas[1]
    boxes=[]; glyphs=[]; pieces=[]; heads=[]
    for s in shapes:
        x0,y0,x1,y1=s['bbox']; w,h=x1-x0,y1-y0
        if s['fill'] and not s['stroke']:
            if w*h>0.04*W*H or max(w,h)>40: continue                # background / decoration (e.g. ring disc)
            glyphs.append(s); continue
        if not s['stroke']: continue
        if s['closed'] and s['fill'] and max(w,h)<=14 and len(s['pts'])<=6:
            heads.append(s); continue                               # arrowhead
        if s['closed'] and w>=12 and h>=7 and s['area']>0.5*w*h:
            boxes.append(dict(bbox=s['bbox'],style=(s['stroke'],s['fill'],s['dash'],round(s['sw'],2)))); continue
        if not s['closed']: pieces.append(s)
    # frames that contain other boxes (a legend, an inset) are not part of the graph: drop them and their contents
    def within(b,f): return f[0]-1<=b[0] and b[2]<=f[2]+1 and f[1]-1<=b[1] and b[3]<=f[3]+1
    frames=[f['bbox'] for f in boxes if sum(1 for o in boxes if o is not f and within(o['bbox'],f['bbox']))>=2]
    if frames:
        boxes=[b for b in boxes if not any(within(b['bbox'],f) for f in frames)]
        glyphs=[g for g in glyphs if not any(within(g['bbox'],f) for f in frames)]
        pieces=[p for p in pieces if not any(within(p['bbox'],f) for f in frames)]
        heads=[h for h in heads if not any(within(h['bbox'],f) for f in frames)]
    classify.frames=len(frames)
    # glyphs -> labels: same colour, same baseline (bottoms align except descenders), small horizontal gap
    n=len(glyphs); par=list(range(n))
    def find(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    G=np.array([g['bbox'] for g in glyphs]) if n else np.zeros((0,4))
    first={}
    for i,g in enumerate(glyphs):                                   # parts of one glyph (the dot of an i)
        if g['eid'] in first: par[find(i)]=find(first[g['eid']])
        else: first[g['eid']]=i
    Hh=G[:,3]-G[:,1] if n else np.zeros(0)
    EM=np.zeros(n)                                                  # font size of the line around each glyph:
    for i in range(n):                                              # tallest glyph sharing its line nearby
        near=(np.abs(G[:,0]-G[i,0])<80)&(np.minimum(G[:,3],G[i,3])-np.maximum(G[:,1],G[i,1])>0)
        EM[i]=max(4.0,float(Hh[near].max()) if near.any() else Hh[i])
    for i in range(n):
        hi=G[i,3]-G[i,1]
        for j in np.where((np.abs(G[:,1]-G[i,1])<40)&(np.arange(n)>i))[0]:
            if glyphs[j]['fill']!=glyphs[i]['fill']: continue
            hj=G[j,3]-G[j,1]; em=max(EM[i],EM[j])
            vert=(min(G[i,3],G[j,3])-max(G[i,1],G[j,1])>0.3*min(hi,hj)      # share part of the line height,
                  or abs((G[i,1]+G[i,3])-(G[j,1]+G[j,3]))/2<0.8*em and min(hi,hj)<0.35*em)   # or underscore/dot near it
            gap=max(G[j,0]-G[i,2],G[i,0]-G[j,2])
            if vert and gap<0.75*em: par[find(i)]=find(j)
    groups={}
    for i in range(n): groups.setdefault(find(i),[]).append(i)
    labels=[]
    for ix in groups.values():
        bb=G[ix]; seen_e=set(); txt=''; last=None
        for i in sorted(ix,key=lambda i:G[i,0]):
            if glyphs[i]['eid'] in seen_e: continue
            seen_e.add(glyphs[i]['eid'])
            if last is not None and G[i,0]-last>0.25*EM[i]: txt+=' '
            txt+=glyphs[i]['ch']; last=G[i,2]
        labels.append(dict(bbox=(bb[:,0].min(),bb[:,1].min(),bb[:,2].max(),bb[:,3].max()),color=glyphs[ix[0]]['fill'],n=len(ix),
                           text=txt,italic=sum(glyphs[i]['italic'] for i in ix)>len(ix)/2))
    labels+= [dict(bbox=t['bbox'],color=t['color'],n=len(t['text']),text=t['text'],italic=t['italic']) for t in texts]
    # line pieces -> connections (end-to-end, same colour)
    conns=[dict(pts=list(p['pts']),stroke=p['stroke'],curved=p['curved'],arrow=False) for p in pieces]
    def close(a,b): return math.dist(a,b)<0.8
    merged=True
    while merged:
        merged=False
        for i in range(len(conns)):
            for j in range(len(conns)):
                if i==j or conns[i] is None or conns[j] is None or conns[i]['stroke']!=conns[j]['stroke']: continue
                a,b=conns[i]['pts'],conns[j]['pts']
                if close(a[-1],b[0]): conns[i]['pts']=a+b[1:]
                elif close(a[-1],b[-1]): conns[i]['pts']=a+b[::-1][1:]
                else: continue
                conns[i]['curved']|=conns[j]['curved']; conns[j]=None; merged=True
        conns=[c for c in conns if c]
    for c in conns: c['head']=set()
    for h in heads:                                                 # which end of each line has an arrowhead
        hx0,hy0,hx1,hy1=h['bbox']
        for c in conns:
            for k,e in ((0,c['pts'][0]),(-1,c['pts'][-1])):
                if hx0-2<=e[0]<=hx1+2 and hy0-2<=e[1]<=hy1+2: c['arrow']=True; c['head'].add(k)
    return boxes,labels,conns,heads

# ------------------------------------------------------------------ metrics
def _segx(a,b,c,d):
    o=lambda p,q,r:(q[0]-p[0])*(r[1]-p[1])-(q[1]-p[1])*(r[0]-p[0])
    d1,d2,d3,d4=o(a,b,c),o(a,b,d),o(c,d,a),o(c,d,b)
    if d1*d2<0 and d3*d4<0:
        t=d1/(d1-d2); return (c[0]+(d[0]-c[0])*t,c[1]+(d[1]-c[1])*t)
    return None
def inside(p,b,m=0.0): return b[0]-m<p[0]<b[2]+m and b[1]-m<p[1]<b[3]+m
def turning(a,b,c):
    v1=(b[0]-a[0],b[1]-a[1]); v2=(c[0]-b[0],c[1]-b[1])
    n1,n2=math.hypot(*v1),math.hypot(*v2)
    if n1<1e-9 or n2<1e-9: return 0.0
    return math.degrees(math.acos(max(-1,min(1,(v1[0]*v2[0]+v1[1]*v2[1])/(n1*n2)))))

def seg_hits(a,b,Bx,margin=0.0):
    """Exact test (Liang-Barsky clipping): does segment a-b pass through each box, grown by margin?
    Vectorised over boxes; replaces point sampling, which can miss a line clipping a box corner."""
    Bx=np.asarray(Bx,float).reshape(-1,4)
    lo_x,lo_y,hi_x,hi_y=Bx[:,0]-margin,Bx[:,1]-margin,Bx[:,2]+margin,Bx[:,3]+margin
    dx,dy=b[0]-a[0],b[1]-a[1]; n=len(Bx)
    t0=np.zeros(n); t1=np.ones(n); ok=np.ones(n,bool)
    for p,q in ((-dx,a[0]-lo_x),(dx,hi_x-a[0]),(-dy,a[1]-lo_y),(dy,hi_y-a[1])):
        if abs(p)<1e-12: ok&=(q>=0)
        elif p<0: t0=np.maximum(t0,q/p)
        else: t1=np.minimum(t1,q/p)
    return ok&(t0<=t1)

JS_CURRENCY=['atp','adp','amp','gtp','gdp','gmp','utp','udp','ump','ctp','cdp','nad','nad+','nadh','nadp','nadp+','nadph','fad',
    'fadh2','fmn','h2o','water','co2','o2','h+','h','pi','ppi','phosphate','pyrophosphate','diphosphate','coa','coa-sh','coash','nh3',
    'nh4+','ammonia','hco3-','bicarbonate','h2o2','e-','camp','sam','sah','acceptor','reduced acceptor','glutathione','nh4']
def _is_cofactor(label):                                          # core.js normName + CURRENCY_DEFAULT
    sub=str.maketrans('₀₁₂₃₄₅₆₇₈₉⁺⁻','0123456789+-'); n=re.sub(r'\s+',' ',(label or '').translate(sub).strip().lower())
    return n in JS_CURRENCY
def _fit_circle(P):
    """Least-squares (Kasa) circle through points; returns centre, radius."""
    P=np.asarray(P,float); A=np.c_[2*P[:,0],2*P[:,1],np.ones(len(P))]; b=(P**2).sum(1)
    (cx,cy,c),*_=np.linalg.lstsq(A,b,rcond=None); return (cx,cy),math.sqrt(max(c+cx*cx+cy*cy,1e-9))
def semantics(boxes,labels,conns,info,by_box):
    """Roles (reaction / metabolite boxes), metabolite cycles and their roundness, substrate/product/catalyst
    sides at each reaction, and crowding between boxes."""
    import networkx as nx
    B=[b['bbox'] for b in boxes]; ctr=[((b[0]+b[2])/2,(b[1]+b[3])/2) for b in B]
    def endbox(p):                                                # the nearest box within 6 px (0 if inside)
        d=[(max(b[0]-p[0],0,p[0]-b[2])**2+max(b[1]-p[1],0,p[1]-b[3])**2,(b[2]-b[0])*(b[3]-b[1]),k) for k,b in enumerate(B)]
        d=[x for x in d if x[0]<=36.0]
        return min(d)[2] if d else None
    def lum(col):
        try: h=col.lstrip('#'); r,g,b=(int(h[k:k+2],16)/255 for k in (0,2,4)); return 0.2126*r+0.7152*g+0.0722*b
        except Exception: return 0.0
    cols=sorted({c['stroke'] for c in conns if c['stroke']},key=lum)
    side_cols={col for col in cols if len(cols)>1 and lum(col)>0.5 and lum(col)-lum(cols[0])>0.2}
    for c in conns: c['side']=c['stroke'] in side_cols
    E=[]
    for i,c in enumerate(conns):
        a,b=endbox(c['pts'][0]),endbox(c['pts'][-1])
        if a is not None and b is not None and a!=b: E.append((i,a,b))
    G=nx.Graph(); G.add_nodes_from(range(len(B))); G.add_edges_from((a,b) for _,a,b in E)
    # roles: the graph is bipartite (compounds - reactions); reactions = the side whose names are italic,
    # else the side whose boxes more often hold two text lines
    role={}; how='none'
    for comp in nx.connected_components(G):
        H=G.subgraph(comp)
        if H.number_of_edges()==0 or not nx.is_bipartite(H): continue
        col=nx.bipartite.color(H); sides={0:[k for k in comp if col[k]==0],1:[k for k in comp if col[k]==1]}
        def score(ks):
            it=[L['italic'] for k in ks for L in by_box.get(k,[])[:1]]; two=[len(by_box.get(k,[]))>=2 for k in ks]
            return (sum(it)/max(1,len(it)), sum(two)/max(1,len(two)))
        sc={k:score(v) for k,v in sides.items()}
        if sc[0][0]!=sc[1][0]: rx=0 if sc[0][0]>sc[1][0] else 1; how='italic reaction names'
        else: rx=0 if sc[0][1]>sc[1][1] else 1; how='two-line reaction boxes'
        for k in sides[rx]: role[k]='rxn'
        for k in sides[1-rx]: role[k]='met'
    name=lambda k:(sorted(by_box.get(k,[]),key=lambda L:L['bbox'][1])[0]['text'] if by_box.get(k) else f'box{k}')
    # line roles at each reaction: arrowhead at the reaction = substrate, at the other end = product, none = catalyst
    S={};P={};C={}
    for i,a,b in E:
        for r,m,end in ((a,b,0),(b,a,-1)):
            if role.get(r)!='rxn' or role.get(m)!='met': continue
            h=conns[i]['head']
            if end in h: S.setdefault(r,[]).append((m,i))
            elif (-1-end if end==0 else 0) in h or (end==-1 and 0 in h) or (end==0 and -1 in h): P.setdefault(r,[]).append((m,i))
            else: C.setdefault(r,[]).append((m,i))
    def direction(r,i):
        p=conns[i]['pts']; p=p if endbox(p[0])==r else p[::-1]
        c=ctr[r]; bb=B[r]; far=max(bb[2]-bb[0],bb[3]-bb[1])/2+8
        q=next((v for v in p if math.dist(v,c)>far),p[-1]); d=(q[0]-c[0],q[1]-c[1]); n=math.hypot(*d) or 1
        return (d[0]/n,d[1]/n)
    ang=lambda u,v: math.degrees(math.acos(max(-1,min(1,u[0]*v[0]+u[1]*v[1]))))
    checked=0; bad=[]; opening=[]; cat_bad=0; cats=0
    for r in sorted(k for k,v in role.items() if v=='rxn'):
        s=[direction(r,i) for _,i in S.get(r,[])]; p=[direction(r,i) for _,i in P.get(r,[])]; c=[direction(r,i) for _,i in C.get(r,[])]
        if not s or not p: continue
        checked+=1                                                # axis from the main lines (secondary ones only
        s_m=[direction(r,i) for m,i in S.get(r,[]) if not conns[i]['side']]   # if a side has no main line)
        p_m=[direction(r,i) for m,i in P.get(r,[]) if not conns[i]['side']]
        ms=np.mean(s_m if (s_m and p_m) else s,0); mp=np.mean(p_m if (s_m and p_m) else p,0); ax=mp-ms; n=np.linalg.norm(ax)
        if n<1e-9: bad.append((r,'substrates and products on the same side')); continue
        ax=ax/n; opening.append(ang(ms/max(np.linalg.norm(ms),1e-9),mp/max(np.linalg.norm(mp),1e-9)))
        wrong=[d for d in s if ang(d,-ax)>80]+[d for d in p if ang(d,ax)>80]
        cw=[d for d in c if not 50<=ang(d,ax)<=130]; cats+=len(c); cat_bad+=len(cw)
        if wrong or cw: bad.append((r,f"{len(wrong)} compound(s) on the wrong side" + (f", {len(cw)} catalyst(s) not beside" if cw else '')))
    placement=dict(role_detection=how,reactions=sum(1 for v in role.values() if v=='rxn'),reactions_checked=checked,
                   reactions_inconsistent=len(bad),catalyst_lines=cats,catalysts_not_beside=cat_bad,
                   mean_substrate_product_angle=round(float(np.mean(opening)),1) if opening else None,
                   inconsistent=[f"{name(r)}: {why}" for r,why in bad][:20])
    # metabolite cycles: compound - reaction - compound steps (substrate -> product); minimum cycle basis
    M=nx.Graph(); step={}
    for r in role:
        if role[r]!='rxn': continue
        for sm,i in S.get(r,[]):
            for pm,j in P.get(r,[]):
                if sm!=pm and not conns[i]['side'] and not conns[j]['side']:     # main steps only
                    M.add_edge(sm,pm); step.setdefault(frozenset((sm,pm)),[]).append(r)
    cand=[c for c in itertools.islice(nx.simple_cycles(M,length_bound=16),20000) if len(c)>=3] if M.number_of_edges() else []
    def shape(c):
        """Closed-loop test: the cycle's boxes in cycle order form a simple polygon that winds once around its
        centre in one direction, without folding inwards (circle, stadium, rounded rectangle all qualify)."""
        g=np.mean([ctr[m] for m in c],0); pts=[]
        for a,b in zip(c,c[1:]+c[:1]):
            pts.append(ctr[a]); rs=step[frozenset((a,b))]; want=(math.dist(ctr[a],g)+math.dist(ctr[b],g))/2
            pts.append(ctr[min(rs,key=lambda r:abs(math.dist(ctr[r],g)-want))])   # the step's box ON the loop
        cen=np.mean(pts,0)
        th=np.unwrap([math.atan2(q[1]-cen[1],q[0]-cen[0]) for q in pts+[pts[0]]]); d=np.diff(th)
        order_ok=bool(np.all(d>-1e-9) or np.all(d<1e-9)) and abs(abs(th[-1]-th[0])-2*math.pi)<0.5
        E_=list(zip(pts,pts[1:]+pts[:1])); simple=True
        for x in range(len(E_)):
            for y in range(x+2,len(E_)):
                if x==0 and y==len(E_)-1: continue
                if _segx(*E_[x],*E_[y]) is not None: simple=False
        turns=[]                                                   # folding inwards = turning the wrong way
        for k in range(len(pts)):
            a_,b_,c_=pts[k-1],pts[k],pts[(k+1)%len(pts)]
            cr=(b_[0]-a_[0])*(c_[1]-b_[1])-(b_[1]-a_[1])*(c_[0]-b_[0]); n_=math.dist(a_,b_)*math.dist(b_,c_) or 1
            turns.append(cr/n_)
        main_sign=np.sign(sum(turns)) or 1
        folds=sum(1 for t in turns if t*main_sign<-0.35)          # sharp turns against the loop direction
        (cx,cy),R=_fit_circle(pts); res=[abs(math.dist(q,(cx,cy))-R)/R for q in pts]
        return pts,float(np.sqrt(np.mean(np.square(res)))),float(max(res)),order_ok and simple and folds==0
    shp={tuple(c):shape(c) for c in cand}
    cyc=[]; used=set(); pool=list(cand)
    isround=lambda c:shp[tuple(c)][3]                           # drawn as a closed loop in cycle order
    while pool:                                                   # main cycles: the longest one drawn round (if it is
        mx=max(len(c) for c in pool)                              # at least half as long as the longest loop), else the
        rnd=[c for c in pool if isround(c) and len(c)>=0.5*mx]    # longest loop, which is then judged not round
        best=max(rnd,key=lambda c:(len(c),-shp[tuple(c)][1])) if rnd else min((c for c in pool if len(c)==mx),key=lambda c:(shp[tuple(c)][1],sorted(c)))
        cyc.append(best); used|=set(best)
        pool=[c for c in pool if not (set(c)&used)]
    alternative=len(cand)-len(cyc)
    rounds=[]; marks=[]
    for c in cyc:
        pts,rms,mx,order_ok=shp[tuple(c)]
        ok=order_ok
        rounds.append(dict(compounds=[name(m) for m in c],round=ok,circle_rms_dev=round(rms,3),closed_loop_in_order=order_ok))
        if not ok: marks+=pts
    cycles=dict(cycles=len(rounds),round_cycles=sum(1 for x in rounds if x['round']),alternative_routes=alternative,detail=rounds)
    # core.js placement: centroid of the inputs vs centroid of the outputs, seen from the reaction: >= 120 degrees
    opp=[]; perp=[]
    for r in sorted(k for k,v in role.items() if v=='rxn'):
        ins=[ctr[m] for m,_ in S.get(r,[])]; outs=[ctr[m] for m,_ in P.get(r,[])]
        if not ins or not outs: continue
        ci=np.mean(ins,0); co=np.mean(outs,0); vi=ci-np.array(ctr[r]); vo=co-np.array(ctr[r])
        if np.linalg.norm(vi)>0 and np.linalg.norm(vo)>0:
            opp.append(ang(vi/np.linalg.norm(vi),vo/np.linalg.norm(vo))>=120)
        flow=vo-vi
        if np.linalg.norm(flow)>0:
            for m,_ in C.get(r,[]):
                d=np.array(ctr[m])-np.array(ctr[r])
                if np.linalg.norm(d)>0: perp.append(50<=ang(d/np.linalg.norm(d),flow/np.linalg.norm(flow))<=130 if False else 60<=ang(d/np.linalg.norm(d),flow/np.linalg.norm(flow))<=120)
    placement.update(js_reactions_scored=len(opp),js_inout_opposite_share=round(float(np.mean(opp)),3) if opp else None,
                     js_modifier_perpendicular_share=round(float(np.mean(perp)),3) if perp else None)
    # core.js "cofactors drawn as clones": a cofactor box (not a clone) linked to more than one reaction
    nbr={}
    for i,a,b in E: nbr.setdefault(a,set()).add(b); nbr.setdefault(b,set()).add(a)
    cof=[k for k in range(len(B)) if role.get(k)=='met' and _is_cofactor(name(k))]
    hubs=[k for k in cof if boxes[k]['style'][2]=='solid' and len(nbr.get(k,()))>1]
    # crowding: boxes that do not overlap but almost touch
    tight=0; mincl=None
    for i in range(len(B)):
        for j in range(i+1,len(B)):
            a,b=B[i],B[j]; dx=max(b[0]-a[2],a[0]-b[2],0); dy=max(b[1]-a[3],a[1]-b[3],0)
            cl=math.hypot(dx,dy) if (dx>0 or dy>0) else 0.0
            if dx>0 or dy>0:
                mincl=cl if mincl is None else min(mincl,cl)
                if cl<6: tight+=1
    return dict(cycles=cycles,placement=placement,tight_pairs=tight,min_clearance=round(mincl,1) if mincl is not None else None,
                cofactor_boxes=len(cof),cofactor_hubs=len(hubs),
                cycle_marks=marks,placement_marks=[ctr[r] for r,_ in bad])

SHARED_BOX_TOL=27.0   # px (= 20 pt): where two lines that end at the same box meet, within this distance of that
                      # box, is not a crossing - the same tolerance as sbml_pathway_layout.py uses
def evaluate(svgfile,compare=None):
    src,shapes,texts,canvas=read(svgfile)
    boxes,labels,conns,heads=classify(shapes,texts,canvas)
    B=[b['bbox'] for b in boxes]
    def end_box(p):                                               # nearest box within 6 px
        d=[(max(b[0]-p[0],0,p[0]-b[2])**2+max(b[1]-p[1],0,p[1]-b[3])**2,k) for k,b in enumerate(B)]
        d=[x for x in d if x[0]<=36.0]
        return [min(d)[1]] if d else []
    issues=dict(crossings=[],line_through_box=[],box_overlaps=[],label_overlaps=[],label_crossed=[],bends=[],long=[])
    # --- connections: endpoints, length, bends, axis alignment
    info=[]
    for c in conns:
        p=c['pts']; L=sum(math.dist(a,b) for a,b in zip(p,p[1:]))
        q=[p[0]]                                                    # drop near-duplicate points
        for v in p[1:]:
            if math.dist(v,q[-1])>0.5: q.append(v)
        p=q if len(q)>=2 else p
        changed=True                                                # remove arrow-rendering spikes:
        while changed and len(p)>2:                                 # a U-turn over a few px is no bend
            changed=False
            for k in range(1,len(p)-1):
                if turning(p[k-1],p[k],p[k+1])>150 and min(math.dist(p[k-1],p[k]),math.dist(p[k],p[k+1]))<6:
                    p=p[:k]+p[k+1:]; changed=True; break
        c['pts']=p
        turns=[turning(p[i-1],p[i],p[i+1]) for i in range(1,len(p)-1)]
        corners=[p[i+1] for i,t in enumerate(turns) if t>25]
        c['curved']=sum(1 for t in turns if 1.0<t<=25)>=3 and sum(t for t in turns if t<=25)>15   # gradual turning = curve
        for v in corners: issues['bends'].append(v)
        dx,dy=p[-1][0]-p[0][0],p[-1][1]-p[0][1]
        ang=math.degrees(math.atan2(abs(dy),abs(dx))) if L>0 else 0; axis_dev=min(ang,90-ang)
        info.append(dict(len=L,bends=len(corners),curved=c['curved'],axis_dev=axis_dev,
                         ends=(end_box(p[0]),end_box(p[-1]))))
    # --- crossings between connections (not at a shared box, not inside a box)
    segs=[(i,a,b) for i,c in enumerate(conns) for a,b in zip(c['pts'],c['pts'][1:])]
    A=np.array([[s[1][0],s[1][1],s[2][0],s[2][1]] for s in segs]) if segs else np.zeros((0,4))
    seen=set()
    for k,(i,a,b) in enumerate(segs):
        bx0=np.minimum(A[:,0],A[:,2]); bx1=np.maximum(A[:,0],A[:,2]); by0=np.minimum(A[:,1],A[:,3]); by1=np.maximum(A[:,1],A[:,3])
        cand=np.where((bx0<=max(a[0],b[0]))&(bx1>=min(a[0],b[0]))&(by0<=max(a[1],b[1]))&(by1>=min(a[1],b[1])))[0]
        for h in cand:
            j=segs[h][0]
            if j<=i or (i,j) in seen: continue
            x=_segx(a,b,segs[h][1],segs[h][2])
            if x is None or any(inside(x,bb,2.0) for bb in B): continue
            ei,ej=conns[i]['pts'],conns[j]['pts']
            if any(math.dist(x,e)<2.0 for e in (ei[0],ei[-1],ej[0],ej[-1])): continue
            shared=(set(info[i]['ends'][0])|set(info[i]['ends'][1]))&(set(info[j]['ends'][0])|set(info[j]['ends'][1]))
            if shared: continue                                      # core.js: lines sharing an end node never count
            seen.add((i,j)); issues['crossings'].append(x)
    # --- lines through boxes that are not their own ends
    for i,c in enumerate(conns):
        own=set(info[i]['ends'][0])|set(info[i]['ends'][1])
        p=c['pts']
        for k,b in enumerate(B):
            if k in own: continue
            shr=(b[0]+2,b[1]+2,b[2]-2,b[3]-2)
            if any(seg_hits(a,q,[shr])[0] for a,q in zip(p,p[1:])):
                issues['line_through_box'].append(((b[0]+b[2])/2,(b[1]+b[3])/2))
    # --- box overlaps
    for i in range(len(B)):
        for j in range(i+1,len(B)):
            a,b=B[i],B[j]
            if a[0]<b[2]-1 and b[0]<a[2]-1 and a[1]<b[3]-1 and b[1]<a[3]-1 and not (a[0]>=b[0] and a[2]<=b[2] and a[1]>=b[1] and a[3]<=b[3]) and not (b[0]>=a[0] and b[2]<=a[2] and b[1]>=a[1] and b[3]<=a[3]):
                issues['box_overlaps'].append(((max(a[0],b[0])+min(a[2],b[2]))/2,(max(a[1],b[1])+min(a[3],b[3]))/2))
    # --- labels: inside a box? fixed offset per box type? overlaps / crossed by lines
    by_box={}; free=[]
    for L in labels:
        cx,cy=(L['bbox'][0]+L['bbox'][2])/2,(L['bbox'][1]+L['bbox'][3])/2
        ks=[k for k,b in enumerate(B) if inside((cx,cy),b)]
        if ks: by_box.setdefault(min(ks,key=lambda k:(B[k][2]-B[k][0])*(B[k][3]-B[k][1])),[]).append(L)
        else: free.append(L)
    types={}
    for k,b in enumerate(boxes):
        types.setdefault(b['style'],[]).append(k)
    type_report=[]
    for st,ks in sorted(types.items(),key=lambda t:-len(t[1])):
        offs=[]; unl=0
        for k in ks:
            Ls=by_box.get(k,[])
            if not Ls: unl+=1; continue
            main=max(Ls,key=lambda L:(L['bbox'][3]-L['bbox'][1],L['n']))   # largest line = the name
            bb=B[k]; w,h=bb[2]-bb[0],bb[3]-bb[1]
            top=sorted(Ls,key=lambda L:L['bbox'][1])[0]                          # first text line of the box
            offs.append((((top['bbox'][0]+top['bbox'][2])/2-(bb[0]+bb[2])/2)/w, top['bbox'][1]-bb[1]))
        o=np.array(offs) if offs else np.zeros((0,2))
        type_report.append(dict(style=dict(outline=st[0],fill=st[1],dash=st[2],width=st[3]),boxes=len(ks),labelled=len(offs),unlabelled=unl,
                                x_offset=round(float(o[:,0].mean()),3) if len(o) else 0.0, x_spread=round(float(o[:,0].std()),3) if len(o) else 0.0,
                                top_gap_px=round(float(o[:,1].mean()),2) if len(o) else 0.0, top_gap_spread_px=round(float(o[:,1].std()),2) if len(o) else 0.0))
    LB=[L['bbox'] for L in labels]
    for i in range(len(LB)):
        for j in range(i+1,len(LB)):
            a,b=LB[i],LB[j]
            ix=min(a[2],b[2])-max(a[0],b[0]); iy=min(a[3],b[3])-max(a[1],b[1])
            small=min((a[2]-a[0])*(a[3]-a[1]),(b[2]-b[0])*(b[3]-b[1]))
            if ix>0 and iy>0 and ix*iy>0.15*small:                  # real overlap, not just touching lines
                issues['label_overlaps'].append(((a[0]+a[2])/2,(a[1]+a[3])/2))
    for i,c in enumerate(conns):
        p=c['pts']
        for L in labels:
            b=L['bbox']
            if any(seg_hits(a,q,[b],-0.5)[0] for a,q in zip(p,p[1:])):
                issues['label_crossed'].append(((b[0]+b[2])/2,(b[1]+b[3])/2))
    # --- long lines
    lens=np.array([x['len'] for x in info]) if info else np.zeros(0)
    med=float(np.median(lens)) if len(lens) else 0.0
    for i,x in enumerate(info):
        if med and x['len']>2.5*med: issues['long'].append(conns[i]['pts'])
    # --- reproducibility
    rep=dict(sha256=hashlib.sha256(open(svgfile,'rb').read()).hexdigest()[:16],
             date_stamp=bool(re.search(r'<dc:date>|creation[_-]?date|CreationDate',src,re.I)),
             random_like_ids=len(re.findall(r'id="[a-z]?[0-9a-f]{12,}"',src)))
    if compare:
        _,shapes2,texts2,_=read(compare)
        g1=sorted((tuple(round(v,2) for v in s['bbox']),s['stroke'],s['fill']) for s in shapes)
        g2=sorted((tuple(round(v,2) for v in s['bbox']),s['stroke'],s['fill']) for s in shapes2)
        rep.update(compared_with=compare,byte_identical=open(svgfile,'rb').read()==open(compare,'rb').read(),
                   geometry_identical=g1==g2,shapes_differing=len(set(g1)^set(g2)))
    # --- roles, cycles, placement, crowding
    sem=semantics(boxes,labels,conns,info,by_box)
    # --- summary
    straight=[x for x in info if not x['curved']]
    out=dict(file=svgfile,
        reconstructed=dict(frames_ignored=getattr(classify,'frames',0),boxes=len(boxes),box_types=len(types),labels=len(labels),labels_inside_boxes=sum(len(v) for v in by_box.values()),
                           free_labels=len(free),connections=len(conns),arrowheads=len(heads),curved_connections=sum(1 for x in info if x['curved'])),
        rule1_reproducible=rep,
        rule2_labels=dict(per_box_type=type_report,label_overlaps=len(issues['label_overlaps']),labels_crossed_by_lines=len(issues['label_crossed']),
                          worst_x_spread=max([t['x_spread'] for t in type_report if t['labelled']],default=0.0),
                          worst_top_gap_spread_px=max([t['top_gap_spread_px'] for t in type_report if t['labelled']],default=0.0)),
        rule3_crossings=dict(line_crossings=len(issues['crossings']),lines_through_boxes=len(issues['line_through_box'])),
        rule4_cycles=sem['cycles'],
        rule5_placement=sem['placement'],
        rule6_long_lines=dict(median=round(med,1),max=round(float(lens.max()),1) if len(lens) else 0,max_over_median=round(float(lens.max())/med,2) if med else 0,
                              over_2_5x_median=len(issues['long']),cv=round(float(lens.std()/lens.mean()),3) if len(lens) and lens.mean() else 0),
        rule7_bends=dict(bends=len(issues['bends']),connections_with_bends=sum(1 for x in info if x['bends']),straight_connections=len(straight),
                         curved_connections=len(info)-len(straight),
                         axis_aligned_share=round(sum(1 for x in straight if x['axis_dev']<=5)/len(straight),3) if straight else 0,
                         mean_axis_deviation_deg=round(float(np.mean([x['axis_dev'] for x in straight])),1) if straight else 0),
        rule8_overlaps=dict(box_overlaps=len(issues['box_overlaps']),tight_pairs=sem['tight_pairs'],min_clearance=sem['min_clearance'],boxes=len(B),
                            cofactor_boxes=sem['cofactor_boxes'],cofactor_hubs=sem['cofactor_hubs']),
        drawing=dict(node_aspect=round((max(b[2] for b in B)-min(b[0] for b in B))/max(1e-9,max(b[3] for b in B)-min(b[1] for b in B)),3) if B else 0.0,
                     node_density=round(sum((b[2]-b[0])*(b[3]-b[1]) for b in B)/max(1e-9,(max(b[2] for b in B)-min(b[0] for b in B))*(max(b[3] for b in B)-min(b[1] for b in B))),4) if B else 0.0,
                     width=round(canvas[2]-canvas[0],1),height=round(canvas[3]-canvas[1],1),aspect=round((canvas[2]-canvas[0])/max(canvas[3]-canvas[1],1),2),
                     box_area_share=round(sum((b[2]-b[0])*(b[3]-b[1]) for b in B)/max((canvas[2]-canvas[0])*(canvas[3]-canvas[1]),1),3)))
    issues['cycles']=sem['cycle_marks']; issues['placement']=sem['placement_marks']
    return out,(boxes,labels,conns,issues,canvas)

T={'reproducible':dict(pass_=0,fail=1),'labels':dict(pass_=0.02,fail=0.2),'cycle':dict(pass_=1,fail=0.5),'clones':dict(pass_=0,fail=3),
   'crossings':dict(pass_=0.02,fail=0.3),'through':dict(pass_=0,fail=5),'placement':dict(pass_=0.8,fail=0.5),
   'aspect':dict(lo=0.5,hi=2,flo=0.25,fhi=4),'longlines':dict(pass_=6,fail=15),'bends':dict(pass_=0.5,fail=2),
   'density':dict(lo=0.05,hi=0.3,flo=0.02,fhi=0.5),'overlaps':dict(pass_=0,fail=5),'weights':dict(hard=2,soft=1)}   # = core.js
def _low(v,p,f): return 1.0 if v<=p else (0.0 if v>=f else (f-v)/(f-p))
def _high(v,p,f): return 1.0 if v>=p else (0.0 if v<=f else (v-f)/(p-f))
def _band(v,t,log):
    F=math.log if log else (lambda x:x)
    if t['lo']<=v<=t['hi']: return 1.0
    if v<t['lo']: return 0.0 if v<=t['flo'] else (F(v)-F(t['flo']))/(F(t['lo'])-F(t['flo']))
    return 0.0 if v>=t['fhi'] else (F(t['fhi'])-F(v))/(F(t['fhi'])-F(t['hi']))
def verdicts(m):
    """The 12 rules shared with core.js (same names, definitions, thresholds and 0-100 scores; cycles: closed-loop rule)."""
    r1=m['rule1_reproducible']; r2=m['rule2_labels']; r3=m['rule3_crossings']; r4=m['rule4_cycles']; r5=m['rule5_placement']
    r6=m['rule6_long_lines']; r7=m['rule7_bends']; r8=m['rule8_overlaps']; d=m['drawing']; rc=m['reconstructed']; out=[]
    def add(rid,kind,name,s_,value,detail,na=False):
        if na: out.append(dict(id=rid,kind=kind,name=name,status='na',score=None,value=value,detail=detail)); return
        out.append(dict(id=rid,kind=kind,name=name,status='pass' if s_>=0.999 else ('fail' if s_<=0.001 else 'check'),score=round(100*s_),value=value,detail=detail))
    if 'byte_identical' in r1:
        v=0 if (r1['byte_identical'] or r1['geometry_identical']) else 1
        add('reproducible','hard','Reproducible rendering',_low(v,**{'p':0,'f':1}) if False else (1.0 if v==0 else 0.0),'identical' if v==0 else 'differs',f"compared with {r1['compared_with']}: bytes {r1['byte_identical']}, geometry {r1['geometry_identical']}")
    else: add('reproducible','hard','Reproducible rendering',None,'–','give a second rendering with --compare to check it',na=True)
    lv=max(r2['worst_x_spread'],r2['worst_top_gap_spread_px']/75.0)+(0.2 if (r2['label_overlaps'] or r2['labels_crossed_by_lines']) else 0)
    add('labels','hard','Labels at fixed positions per type',_low(lv,0.02,0.2),f"{lv:.3f}",
        f"within each box type the label varies by {round(100*r2['worst_x_spread'],1)}% of box width sideways and {r2['worst_top_gap_spread_px']} px from the top; {r2['label_overlaps']} label overlaps, {r2['labels_crossed_by_lines']} labels crossed")
    if r4['cycles']: add('cycle','hard','Cycles look like cycles',_high(r4['round_cycles']/r4['cycles'],1,0.5),f"{r4['round_cycles']}/{r4['cycles']} closed loops",
                         "main metabolite cycles drawn as a closed loop in cycle order (circle, stadium, rounded rectangle...)"+(f"; {r4['alternative_routes']} alternative route(s) not judged" if r4.get('alternative_routes') else ''))
    else: add('cycle','hard','Cycles look like cycles',None,'no cycle','no metabolite cycle in this drawing',na=True)
    add('clones','hard','Cofactors drawn as clones',_low(r8['cofactor_hubs'],0,3),f"{r8['cofactor_hubs']} hubs",f"{r8['cofactor_boxes']} cofactor boxes (core.js list); {r8['cofactor_hubs']} shared by more than one reaction without being a clone")
    per=r3['line_crossings']/max(1,rc['connections'])
    add('crossings','hard','No crossings',_low(per,0.02,0.3),f"{r3['line_crossings']} crossings",f"{per:.3f} per line across {rc['connections']} lines (lines sharing an end node never count)")
    add('through','hard','No line through a node',_low(r3['lines_through_boxes'],0,5),f"{r3['lines_through_boxes']} lines","lines passing through a box they do not connect (exact)")
    sh=r5.get('js_inout_opposite_share')
    if sh is not None: add('placement','hard','Consistent participant placement',_high(sh,0.8,0.5),f"{round(100*sh)}% opposite",
                           f"substrates and products on opposite sides (>= 120 deg between their centroids) in {round(100*sh)}% of {r5['js_reactions_scored']} reactions; modifiers beside: {r5.get('js_modifier_perpendicular_share')}")
    else: add('placement','hard','Consistent participant placement',None,'–','no reaction with both a substrate and a product',na=True)
    asp=d['node_aspect']; add('aspect','hard','Aspect ratio',_band(asp,T['aspect'],True),f"{asp:.2f}","width / height of the nodes' extent; target 0.5-2")
    ll=r6['max_over_median']; add('longlines','soft','Avoid long lines',_low(ll,6,15),f"{ll}x","longest line vs median")
    bp=r7['bends']/max(1,rc['connections']); add('bends','soft','Avoid bends',_low(bp,0.5,2),f"{bp:.2f} per line",f"{r7['bends']} bend points on {rc['connections']} lines")
    den=d['node_density']; add('density','soft','Density',_band(den,T['density'],False),f"{round(100*den)}%","node area / extent of the nodes; target 5-30%")
    add('overlaps','soft','No overlapping nodes',_low(r8['box_overlaps'],0,5),str(r8['box_overlaps']),f"overlapping node pairs ({r8['tight_pairs']} more pairs closer than 6 px, for information)")
    ws=sum(T['weights'][x['kind']]*x['score'] for x in out if x['score'] is not None); wt=sum(T['weights'][x['kind']] for x in out if x['score'] is not None)
    m['overall']=round(ws/wt) if wt else None
    m['hard_pass']=f"{sum(1 for x in out if x['kind']=='hard' and x['status']=='pass')}/{sum(1 for x in out if x['kind']=='hard' and x['status']!='na')}"
    return out

def annotate(geo,path,title):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    boxes,labels,conns,issues,canvas=geo
    W,H=canvas[2]-canvas[0],canvas[3]-canvas[1]
    fig,ax=plt.subplots(figsize=(min(16,W/60+2),min(16,H/60+2))); ax.set_xlim(canvas[0],canvas[2]); ax.set_ylim(canvas[3],canvas[1]); ax.set_aspect('equal'); ax.axis('off')
    for c in conns: ax.plot([p[0] for p in c['pts']],[p[1] for p in c['pts']],color='#c9c9c9',lw=0.8)
    for p in issues['long']: ax.plot([q[0] for q in p],[q[1] for q in p],color='#2f6fd6',lw=2.4,alpha=0.8)
    for b in boxes: ax.add_patch(Rectangle(b['bbox'][:2],b['bbox'][2]-b['bbox'][0],b['bbox'][3]-b['bbox'][1],fc='none',ec='#777',lw=0.6))
    for L in labels: ax.add_patch(Rectangle(L['bbox'][:2],L['bbox'][2]-L['bbox'][0],L['bbox'][3]-L['bbox'][1],fc='#ffd76a',ec='none',alpha=0.35))
    mk=[('crossings','#e3262b','X',90,'crossing'),('line_through_box','#ff8c1a','s',60,'line through box'),('box_overlaps','#9b30d9','D',60,'box overlap'),
        ('label_overlaps','#7a4a00','o',40,'label overlap'),('label_crossed','#c2185b','^',40,'label crossed'),('bends','#00a07a','o',25,'bend'),
        ('cycles','#ff00aa','h',70,'node of a cycle not drawn as a closed loop'),('placement','#0b7285','P',90,'reaction with inconsistent sides')]
    for key,col,m,s,lab in mk:
        if issues[key]: ax.scatter([p[0] for p in issues[key]],[p[1] for p in issues[key]],c=col,marker=m,s=s,label=f'{lab} ({len(issues[key])})',zorder=5)
    if issues['long']: ax.plot([],[],color='#2f6fd6',lw=2.4,label=f"long line ({len(issues['long'])})")
    if ax.get_legend_handles_labels()[0]: ax.legend(loc='lower right',fontsize=8,frameon=True)
    else: ax.text(0.99,0.01,'no problems found',transform=ax.transAxes,ha='right',va='bottom',fontsize=9,color='#2a8a2a')
    ax.set_title(title,fontsize=10)
    fig.savefig(path,dpi=110,bbox_inches='tight'); plt.close(fig)

def main():
    ap=argparse.ArgumentParser(description='Layout-rule metrics for an SVG drawing')
    ap.add_argument('svg'); ap.add_argument('--compare'); ap.add_argument('--json'); ap.add_argument('--annotate')
    a=ap.parse_args()
    m,geo=evaluate(a.svg,a.compare)
    r=m['reconstructed']
    print(f"{a.svg}\n  reconstructed: {r['boxes']} boxes of {r['box_types']} types, {r['labels']} labels ({r['labels_inside_boxes']} inside boxes, {r['free_labels']} free), "
          f"{r['connections']} lines ({r['curved_connections']} curved), {r['arrowheads']} arrowheads"+(f"; {r['frames_ignored']} framed region(s) such as a legend ignored" if r['frames_ignored'] else ''))
    V=verdicts(m); print(f"  overall {m['overall']}   hard rules passed {m['hard_pass']}")
    for x in V: print(f"  {x['kind']:<4} {x['name']:<36} {x['status']:<5} {str(x['score']) if x['score'] is not None else '-':>4}  {x['value']}  ({x['detail']})")
    print(f"  label position per box type (sideways offset of the first text line from the box centre, as % of box width; gap to the box top in px):")
    for t in m['rule2_labels']['per_box_type']:
        print(f"    outline {t['style']['outline']:<8} fill {t['style']['fill']:<8} {t['style']['dash']:<10} {t['boxes']:>3} boxes, {t['unlabelled']} unlabelled, "
              f"offset {round(100*t['x_offset'],1)}% +/- {round(100*t['x_spread'],1)}%, top gap {t['top_gap_px']} +/- {t['top_gap_spread_px']} px")
    for c in m['rule4_cycles']['detail']:
        print(f"  cycle {'closed loop in order' if c['round'] else 'NOT a closed loop in order'} (deviation from a circle: rms {c['circle_rms_dev']}): {' -> '.join(c['compounds'])}")
    for t in m['rule5_placement']['inconsistent']: print(f"  inconsistent placement: {t}")
    d=m['drawing']; print(f"  drawing: {d['width']} x {d['height']} (aspect {d['aspect']}), boxes cover {round(100*d['box_area_share'],1)}% of the area")
    if a.json:
        m['verdicts']=verdicts(m)
        json.dump(m,open(a.json,'w'),indent=1)
    if a.annotate: annotate(geo,a.annotate,a.svg.split('/')[-1]+': problems marked')

if __name__=='__main__': main()

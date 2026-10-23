import pymupdf as fitz, numpy as np, json, sys
d=fitz.open(sys.argv[1])
panels={ # name: page idx, box drawing idx, xticks, yticks, curve drawing idxs (legend order)
 '4b':(4,1698,np.arange(0,0.61,0.1),np.arange(0,2.51,0.5),{'0.00032':1793,'0.032':1794,'0.128':1795,'0.384':1810}),
 '5a':(4,1881,np.arange(0,1.01,0.2),np.arange(0,3.51,0.5),{'0.032':1978,'0.128':1979,'0.384':2025}),
 '5b':(4,2111,np.arange(0,2.01,0.5),np.arange(0,4.51,0.5),{'0.032':2209,'0.128':2210,'0.384':2297}),
 '6a':(5,0,np.arange(0,3.01,0.5),np.arange(0,6.01,1),{'0.032':87,'0.128':88,'0.384':219}),
 '6b':(5,390,np.arange(0,4.01,1),np.arange(0,8.01,1),{'0.032':471,'0.128':472,'0.384':650}),
}
def groups(vals, tol):
    vals=sorted(vals); out=[[vals[0]]]
    for v in vals[1:]:
        (out[-1].append(v) if v-out[-1][-1]<tol else out.append([v]))
    return out
res={}
for name,(pg,bi,xt,yt,curves) in panels.items():
    dr=d[pg].get_drawings(); box=dr[bi]['rect']
    glyphs=[x['rect'] for x in dr if x['color'] is None and x['fill']==(0,0,0) and x['rect'].width<8 and x['rect'].height<8]
    # x labels: below box
    xl=[g for g in glyphs if box.y1+2<g.y0<box.y1+10 and box.x0-8<g.x0<box.x1+8]
    gx=[]
    for grp in groups([g.x0 for g in xl]+[g.x1 for g in xl],2.2): pass
    # group glyph rects by horizontal proximity
    xl=sorted(xl,key=lambda g:g.x0); G=[[xl[0]]]
    for g in xl[1:]:
        (G[-1].append(g) if g.x0-G[-1][-1].x1<2.0 else G.append([g]))
    cx=np.array([(min(g.x0 for g in grp)+max(g.x1 for g in grp))/2 for grp in G])
    assert len(cx)==len(xt),(name,len(cx))
    ax=np.polyfit(xt,cx,1); rx=cx-np.polyval(ax,xt)
    xr_fit=(box.x1-ax[1])/ax[0]; xlim=round(xr_fit,1)
    print(name,'x tick-fit: box left %.4f right %.4f -> xlim %.1f, resid %.2fpt'%((box.x0-ax[1])/ax[0],xr_fit,xlim,np.abs(rx).max()))
    sx=box.width/xlim; sy=box.height/yt[-1]
    # use box for the origin (axes lines are exact), tick fit for the scale
    res[name]={'_meta':dict(page=pg,box=list(box),xlim=xlim,ylim=float(yt[-1]))}
    for rate,ci in curves.items():
        pts=[]
        for it in dr[ci]['items']:
            if it[0]!='l': continue
            for P in (it[1],it[2]):
                if not pts or abs(pts[-1][0]-P.x)+abs(pts[-1][1]-P.y)>1e-6: pts.append((P.x,P.y))
        pts=np.array(pts)
        X=(pts[:,0]-box.x0)/sx; Y=(box.y1-pts[:,1])/sy
        res[name][rate]=np.c_[X,Y].tolist()
        print('   ',rate,len(X),'xmax %.4f'%X.max(),'ymax %.4f'%Y.max(),'start',np.round([X[0],Y[0]],4),'end',np.round([X[-1],Y[-1]],4))
# Fig. 10a (page 9): equilibrium points = centres of the diamond markers (multi-step relaxation)
dr=d[8].get_drawings(); box=dr[0]['rect']
xl=sorted([x['rect'] for x in dr if x['color'] is None and x['fill']==(0,0,0) and x['rect'].width<8 and x['rect'].height<8
           and box.y1+2<x['rect'].y0<box.y1+10 and box.x0-8<x['rect'].x0<box.x1+8],key=lambda g:g.x0)
cx=np.array([(g.x0+g.x1)/2 for g in xl]); ax=np.polyfit(np.arange(5.0),cx,1)   # single-glyph labels 0..4
xlim=round((box.x1-ax[1])/ax[0],1); sx=box.width/xlim; sy=box.height/4.0
legend=[x['rect'] for x in dr if x['color'] is not None and box.contains(x['rect']) and {i[0] for i in x['items']}=={'re'}][0]
dia=[x['rect'] for x in dr if x['color'] is not None and box.contains(x['rect'])
     and {i[0] for i in x['items']}=={'qu'} and not legend.intersects(x['rect'])]
pts=sorted({(round(((r.x0+r.x1)/2-box.x0)/sx,4),round((box.y1-(r.y0+r.y1)/2)/sy,4)) for r in dia})
legend=None
print('10a xlim',xlim,'multi-step points',pts)
res['10a']={'_meta':dict(page=8,box=list(box),xlim=xlim,ylim=4.0),'multistep':[list(p) for p in pts]}
json.dump(res,open('raw.json','w'))

# Gallery of the tracks whose collapsed drawing crosses, from the output
# of find_crossings.  Usage:
#
#   python3 ../devel/iss013/gallery.py N failures.txt out.tex TOTAL
#
# from latex/, say, then pdflatex out.tex.  Needs examples/ttplot built.
#
# where TOTAL is the number of tracks searched, as TeX (e.g. "31\\,370").
# Each track gets the whole drawing scaled to fit, and a detail at true
# scale around its crossings, which are circled in red.
import os, re, subprocess
BIN=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","examples","ttplot")
NUM=re.compile(r"\(([-\d.]+),([-\d.]+)\)")
rows=[]
import sys
N=int(sys.argv[1]); SRC=sys.argv[2]; OUT=sys.argv[3]; TOTAL=sys.argv[4]
for line in open(SRC):
    s,v,h,pts,c=line.rstrip("\n").split("|")
    P=[tuple(map(float,p.split(","))) for p in pts.split(";") if p]
    rows.append((float(h),int(s),int(v),P,c))
rows.sort()

def draws(coding):
    out=subprocess.run([BIN,"--snippet","--coding",coding,"--output","-"],
                       capture_output=True,text=True).stdout
    return [l for l in out.splitlines() if l.strip().startswith(("\\draw","\\fill"))]

def xf(lines,sx,sy,ox=0.0,oy=0.0):
    f=lambda m: "(%.4f,%.4f)"%((float(m.group(1))-ox)*sx,(float(m.group(2))-oy)*sy)
    return [NUM.sub(f,l) for l in lines]

def chains(lines):
    out=[]
    for l in lines:
        if not l.strip().startswith("\\draw[ttedge]"): continue
        q=[(float(a),float(b)) for a,b in NUM.findall(l)]
        out.append([q[i:i+4] for i in range(0,len(q)-1,3)])
    return out

def bez(c,t):
    u=1-t
    return (u**3*c[0][0]+3*u*u*t*c[1][0]+3*u*t*t*c[2][0]+t**3*c[3][0],
            u**3*c[0][1]+3*u*u*t*c[1][1]+3*u*t*t*c[2][1]+t**3*c[3][1])

def detail_panel(lines,P,g,ox,oy,w,h):
    t=["\\begin{tikzpicture}[line cap=round,line join=round]",
       "\\clip (0,0) rectangle (%.4f,%.4f);"%(w,h)]
    inside=lambda x,y: -1<=x<=w+1 and -1<=y<=h+1
    for ch in chains(lines):
        for c in ch:
            pts=[((x-ox)*g,(y-oy)*g) for (x,y) in (bez(c,i/600.0) for i in range(601))]
            run=[]
            for q in pts+[None]:
                if q is not None and inside(*q): run.append(q); continue
                if len(run)>1:
                    t.append("\\draw[line width=0.6pt] "+" -- ".join("(%.4f,%.4f)"%r for r in run)+";")
                run=[]
    for l in lines:
        m=re.match(r"\s*\\fill\[black\] \(([-\d.]+),([-\d.]+)\)",l)
        if m:
            x,y=(float(m.group(1))-ox)*g,(float(m.group(2))-oy)*g
            if inside(x,y): t.append("\\fill (%.4f,%.4f) circle (0.07);"%(x,y))
    for (x,y) in P:
        if inside((x-ox)*g,(y-oy)*g):
            t.append("\\draw[red,line width=0.8pt] (%.4f,%.4f) circle (0.18);"%((x-ox)*g,(y-oy)*g))
    t.append("\\end{tikzpicture}")
    return "\n".join(t)

def panel(lines,P,sx,sy,ox,oy,clip=None,dot=0.07):
    t=["\\begin{tikzpicture}[line cap=round,line join=round]",
       "\\tikzset{ttedge/.style={line width=0.6pt}}"]
    if clip: t.append("\\clip (%.4f,%.4f) rectangle (%.4f,%.4f);"%clip)
    for l in xf(lines,sx,sy,ox,oy):
        l=re.sub(r"circle \([\d.]+\)","circle (%.3f)"%dot,l)
        t.append(l)
    for (x,y) in P:
        t.append("\\draw[red,line width=0.8pt] (%.4f,%.4f) circle (0.18);"%((x-ox)*sx,(y-oy)*sy))
    t.append("\\end{tikzpicture}")
    return "\n".join(t)

W,H=6.5,4.5
doc=[r"""\documentclass[10pt]{article}
\usepackage[margin=0.6in]{geometry}
\usepackage{tikz}
\setlength{\parindent}{0pt}
\begin{document}
\section*{Collapsed drawings that still cross, $n=%d$}
All %d automaton vertices for $n=%d$ whose collapsed drawing has a proper
crossing, of %s.  For $n\le7$ there are none since control lengths were
capped so that no curve runs past the point it is heading for.  Crossings
are circled in red.  Sorted by the height of the drawing.  Left: the whole
drawing, scaled to fit.  Right: a detail at true scale around the
crossings.  Stratum and vertex are 1-based, as the \texttt{ttauto} program
prints them.
\bigskip
""" % (N,len(rows),N,TOTAL)]
for k,(h,s,v,P,c) in enumerate(rows):
    L=draws(c)
    f=min(W/float(N), H/max(h,1e-9))
    whole=panel(L,[],f,f,0.5,0.0,dot=0.05)
    ys=[p[1] for p in P]
    y0=max(0.0,min(ys)-1.2); y1=max(ys)+1.2
    if y1-y0>4.0: y0,y1=min(ys)-1.2,min(ys)+2.8
    g=0.6
    detail=detail_panel(L,P,g,0.5,y0,N*g,(y1-y0)*g)
    doc.append(r"\begin{minipage}{\textwidth}")
    doc.append(r"\textbf{%d.} stratum %d, vertex %d, height %.1f, %d crossing%s\\"%(
        k+1,s,v,h,len(P),"" if len(P)==1 else "s"))
    doc.append(r"{\footnotesize\texttt{%s}}\\[3pt]"%c)
    doc.append(r"\parbox[b][%.1fcm][b]{%.1fcm}{%s}\hfill\fbox{%s}"%(H+0.2,W+0.3,whole,detail))
    doc.append(r"\end{minipage}\par\bigskip")
doc.append(r"\end{document}")
open(OUT,"w").write("\n".join(doc))

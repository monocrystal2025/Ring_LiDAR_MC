from pathlib import Path
import numpy as np, math, json
from html import escape
P=Path(__file__).parent
INK='#20252b'; BLUE='#4a6477'; TEAL='#238578'
class SVG:
 def __init__(self): self.e=['<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="850" viewBox="0 0 1200 850"><rect width="1200" height="850" fill="white"/>']
 def path(self,p,stroke=INK,w=2,fill='none',opacity=1,close=False,dash=None):
  d='M'+' L'.join(f'{x:.3f},{y:.3f}' for x,y in np.asarray(p)) + (' Z' if close else '')
  self.e.append(f'<path d="{d}" stroke="{stroke}" stroke-width="{w}" fill="{fill}" opacity="{opacity}" stroke-linejoin="round" stroke-linecap="round"'+(f' stroke-dasharray="{dash}"' if dash else '')+'/>')
 def text(self,x,y,t,size=24,anchor='start',color=INK):
  self.e.append(f'<text x="{x}" y="{y}" font-family="Arial, sans-serif" font-size="{size}" text-anchor="{anchor}" fill="{color}">{escape(t)}</text>')
 def circle(self,x,y,r,fill,stroke='none',w=1):
  self.e.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{w}"/>')
 def arrow(self,a,b,color=INK,w=2,head=10):
  a=np.array(a); b=np.array(b); d=(b-a)/np.linalg.norm(b-a); n=np.array([-d[1],d[0]])
  self.path([a,b],color,w); self.path([b,b-head*d+head*.38*n,b-head*d-head*.38*n],color,0,color,close=True)
 def save(self,name): (P/(name+'.svg')).write_text(''.join(self.e)+'</svg>',encoding='utf8')
 def drone(self,x,y,s=1):
  # Same original orthographic icon in all panels; based on quadcopter reference morphology.
  self.e.append(f'<g transform="translate({x},{y}) scale({s})">')
  for a,b in [((-25,-12),(25,12)),((-25,12),(25,-12))]:
   self.path([a,b],'#555d64',5)
  self.e.append('<rect x="-9" y="-13" width="18" height="26" rx="5" fill="#505961"/><rect x="-5" y="-9" width="10" height="10" rx="2" fill="#aab3ba"/>')
  for x0,y0 in [(-25,-12),(25,12),(-25,12),(25,-12)]:
   self.e.append(f'<ellipse cx="{x0}" cy="{y0}" rx="14" ry="5" fill="#f9fafb" stroke="#505961" stroke-width="1.8"/><circle cx="{x0}" cy="{y0}" r="2.7" fill="#505961"/>')
  self.e.append('</g>')
def norm(v): return v/np.linalg.norm(v)
az=.55
M=np.array([[np.cos(az),np.sin(az),0],[-.31*np.sin(az),.31*np.cos(az),-.9507]])
origin=np.array([595.,670.]); scale=470
def proj(a): return np.asarray(a)@M.T*scale+origin
def circle3(u,r):
 t=np.linspace(0,2*np.pi,181)
 v=norm(np.cross(u,[0,0,1])); w=np.cross(u,v)
 return u+ r*(np.cos(t)[:,None]*v+np.sin(t)[:,None]*w)
def hull(p):
 p=sorted(set(map(tuple,p)))
 def cross(o,a,b): return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
 low=[]; up=[]
 for q in p:
  while len(low)>1 and cross(low[-2],low[-1],q)<=0: low.pop()
  low.append(q)
 for q in reversed(p):
  while len(up)>1 and cross(up[-2],up[-1],q)<=0: up.pop()
  up.append(q)
 return low[:-1]+up[:-1]
beam=np.load(r'G:\BeamVEC_NEW\beam_vector_fast_5000Hz_h3000m_D325000u.npy')
def ringat(u):
 r=np.hypot(*u[:2]); t0=np.arctan2(u[1],u[0]); ts=np.arange(t0,t0+2*np.pi+1e-10,2*np.pi/5000/r)
 return np.c_[r*np.cos(ts),r*np.sin(ts),np.full(len(ts),u[2])]
rings=[ringat(beam[0]),ringat(beam[-1])]
s=SVG(); s.text(45,58,'(a)',32); s.text(600,58,'ROE mode',29,'middle')
# convex projected hemisphere silhouette, flat tint
t=np.linspace(0,2*np.pi,400)
ground=np.c_[np.cos(t),np.sin(t),np.zeros(len(t))]
mer=[]
for p in np.linspace(0,np.pi/2,90):
 mer.extend(np.c_[np.sin(p)*np.cos(t[::5]),np.sin(p)*np.sin(t[::5]),np.full(len(t[::5]),np.cos(p))])
s.path(hull(proj(mer)),BLUE,1.4,'#edf3f7',1,True)
s.path(proj(ground),BLUE,1.3)
s.path([origin,proj([0,0,1.13])],'#72808b',1.2,dash='6 6'); s.text(605,125,'Z',21)
# continuous original sampled trajectory
s.path(proj(beam),INK,1.9)
for rr in rings: s.path(proj(rr),INK,2.1)
centers=[]; checks=[]
for target,kind in [([270,400],'Spot'),([510,245],'Line'),([945,440],'Annular')]:
 pp=proj(beam); idx=np.argmin(np.linalg.norm(pp-target,axis=1)); u=beam[idx]; c=pp[idx]; centers.append((idx,kind,c))
 tangent=norm(beam[min(idx+1,len(beam)-1)]-beam[max(idx-1,0)])
 tangent=norm(tangent-u*np.dot(tangent,u)); long=norm(np.cross(tangent,u))
 r=np.tan(.325/2)
 if kind=='Line':
  a=r*long; b=.010*tangent; end=np.array([u-a-b,u+a-b,u+a+b,u-a+b,u-a-b])
 else: end=circle3(u,r)
 s.path(hull(np.vstack([origin,proj(end)])),TEAL,.7,TEAL,.13,True)
 # true shared apex and exact terminal perimeter
 for jj in [np.argmin(proj(end)[:,0]),np.argmax(proj(end)[:,0])]:
  s.path([origin,proj(end[jj])],TEAL,1.2,opacity=.65)
 if kind=='Annular':
  outer=proj(end); inner=proj(circle3(u,r*.66))
  d='M'+' L'.join(f'{x:.3f},{y:.3f}' for x,y in outer)+' Z M'+' L'.join(f'{x:.3f},{y:.3f}' for x,y in inner)+' Z'
  s.e.append(f'<path d="{d}" fill="{TEAL}" fill-opacity=".8" fill-rule="evenodd" stroke="{TEAL}" stroke-width="1.5"/>')
  for rr in [.52,.75]: s.path(proj(end*rr),TEAL,.9,opacity=.65)
 else:
  s.path(proj(end),TEAL,1.5,TEAL,.8,True)
  if kind=='Spot':
   for rr in [.42,.70]: s.path(proj(end*rr),TEAL,1.1,TEAL,.23,True)
 s.circle(*c,3.4,INK)
 # targets shown above footprint to avoid hiding field shape, leader ties target to illuminated region
 dc=c+np.array([0,-57 if kind!='Line' else -48])
 s.drone(*dc,.72)
 s.text(c[0]+(25 if kind=='Line' else 0),c[1]+(39 if kind!='Line' else 58),kind,21,'middle')
 checks.append({'kind':kind,'npy_index':int(idx),'center_error':0.,'long_dot_tangent':float(np.dot(long,tangent))})
# polar angle is with Z, separate reference pointing ray
u=beam[centers[-1][0]]; phi=np.arccos(u[2]); azi=np.arctan2(u[1],u[0])
q=np.linspace(0,phi,100)
arc=.24*np.c_[np.sin(q)*np.cos(azi),np.sin(q)*np.sin(azi),np.cos(q)]
s.path([origin,proj(u*.61)],'#65737d',1.1,dash='4 5'); s.path(proj(arc),INK,1.4)
s.text(* (proj(arc[len(arc)//2])+[16,2]),'φ',27)
s.circle(*origin,6,'#ba6240'); s.e.append('<rect x="550" y="706" width="90" height="28" fill="white"/>'); s.text(origin[0],origin[1]+56,'LiDAR',24,'middle')
s.text(85,795,'φ: polar angle from +Z',20)
s.save('Fig1a_ROE_v2')
# ROI: unit spherical cone sector, rotated for side view
s=SVG(); s.text(45,58,'(b)',32); s.text(600,58,'ROI mode',29,'middle')
o=np.array([120.,530.]); axis=np.array([820.,-60.]); e1=np.array([58.,12.]); e2=np.array([0.,410.]); half=np.pi/6
tt=np.linspace(0,2*np.pi,240)
def roi_slice(z,r=None):
 if r is None:r=z*np.tan(half)
 return o+z*axis+r*(np.cos(tt)[:,None]*e1+np.sin(tt)[:,None]*e2)
end=roi_slice(1)
s.path(hull(np.vstack([o,end])),BLUE,1.7,'#edf3f7',1,True)
s.path(end,BLUE,1.7,fill='#e3edf4',opacity=1,close=True)
s.path([o,o+axis],'#8b98a2',1.2,dash='7 6')
for z in [.44,.72]: s.path(roi_slice(z),BLUE,1.1,opacity=.65,dash='5 5')
# linear diverging spot beam, transverse coordinates fixed angle, all radii grow linearly
bc=o+axis+np.array([0,-105.]); radx=12.; rady=22.
ends=bc+np.c_[radx*np.cos(tt),rady*np.sin(tt)]
s.path(hull(np.vstack([o,ends])),TEAL,.8,TEAL,.16,True)
for z in [.44,.72,1.]:
 cross=o+z*(ends-o); s.path(cross,TEAL,1.4,TEAL,.70,True)
s.circle(*o,6,'#ba6240'); s.text(o[0],o[1]+39,'LiDAR',24,'middle')
# icon and visibly oblique velocity vector
pos=o+.70*axis+np.array([0,37])
s.drone(*pos,.95)
s.arrow(pos+[35,-16],pos+[118,-98],'#a66a28',2.5)
s.text(*(pos+[112,-113]),'v',28)
s.text(1060,670,'ROI',27,'middle'); s.text(1060,700,'surveillance volume',18,'middle')
s.path([[1004,644],[855,609],[748,581]],BLUE,1.3); s.circle(748,581,3,BLUE)
s.text(715,179,'Single-shot footprint',23)
s.arrow([925,190],bc+[8,-26],INK,1.2,8)
# full angle aperture follows cone silhouette rays
up=end[np.argmin(np.arctan2(end[:,1]-o[1],end[:,0]-o[0]))]
down=end[np.argmax(np.arctan2(end[:,1]-o[1],end[:,0]-o[0]))]
angs=np.linspace(np.arctan2(*(up-o)[::-1]),np.arctan2(*(down-o)[::-1]),80)
s.path(o+110*np.c_[np.cos(angs),np.sin(angs)],INK,1.3)
s.text(248,525,'2θ',24); s.text(90,791,'θ = 30°: ROI half-angle',20)
s.save('Fig1b_ROI_v2')
# exact original MATLAB numerical ROI paths, angular-plane inverse map
s=SVG(); s.text(45,58,'(c)',32)
rho0=np.pi/6; rad=rho0+.05; fac=380
for cx,key,title in [(320,'raster','Reciprocating raster'),(885,'spiral','Equidistant spiral')]:
 v=np.loadtxt(P/('roi_'+key+'.csv'),delimiter=','); th=np.arccos(np.clip(v[:,2],-1,1)); az=np.arctan2(v[:,1],v[:,0])
 xy=np.c_[th*np.cos(az),th*np.sin(az)]
 cen=np.array([cx,415.])
 s.circle(cx,415,rho0*fac,'#e8f0f5',BLUE,1.5)
 circ=cen+rad*fac*np.c_[np.cos(tt),np.sin(tt)]
 s.path(circ,'#a4aeb6',1.1,dash='5 5')
 pix=cen+xy*np.array([fac,-fac]); s.path(pix,INK,1.75)
 s.text(cx,135,title,26,'middle')
 for frac in ([.13,.44,.75] if key=='raster' else [.23,.55,.85]):
  i=int(frac*len(pix)); s.arrow(pix[i-9],pix[i+9],INK,1.6,8)
 s.circle(*pix[0],4.3,TEAL); s.circle(*pix[-1],4.3,'white',INK,1.5)
s.path([[290,700],[345,700]],INK,1.8); s.text(360,708,'Beam-axis trajectory',21)
s.path([[660,700],[715,700]],'#a4aeb6',1.2,dash='5 5'); s.text(730,708,'Scan envelope',21)
s.text(600,760,'θ = 30°     wD = 100 mrad     ρmax = θ + wD / 2',21,'middle')
s.text(600,802,'Angular-plane view; return sweep retraces the same path',19,'middle',color='#596874')
s.save('Fig1c_Trajectories_v2')
(P/'geometry_checks.json').write_text(json.dumps({'ROE_source':r'G:\BeamVEC_NEW\beam_vector_fast_5000Hz_h3000m_D325000u.npy','ROE_polar_deg':np.rad2deg(np.arccos(beam[[0,-1],2])).tolist(),'beam_checks':checks,'note':'Beam widths exaggerated equally for legibility; line width and annulus thickness emphasized. Line long axis is exactly normal to scan tangent per user request, differs from legacy ROE meridional orientation. ROI paths are unmodified MATLAB output plotted by inverse angular-plane mapping.'},indent=2))
print('Wrote 3 SVG panels and geometry audit.')


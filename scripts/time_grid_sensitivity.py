import numpy as np, pandas as pd
from scipy.integrate import solve_ivp
from pathlib import Path
p=Path(__file__).resolve().parents[1] / 'data'
h=pd.read_excel(p/'simulation_summary_with_band.xlsx',sheet_name='Virtual_Hosts',usecols=['run','dose','eta','phi','max_pathogen_burden','infected'])
h=h[h['run']==1].groupby('dose',sort=False).head(10).reset_index(drop=True)
t400=np.linspace(0,30,400); t1600=np.linspace(0,30,1600)
results=[]; failure=0
for row in h.itertuples(index=False):
 def rhs(t,y):
  P,I=np.maximum(y,0)
  return [P*(1-P/1e8-row.phi*I-.1), .001*P+.3*I*(1-I/1e5)-.2*I]
 sol=solve_ivp(rhs,(0,30),(row.eta*row.dose,4e4),method='RK45',rtol=1e-6,atol=1e-9,dense_output=True)
 if not sol.success: failure+=1;continue
 p400=np.maximum(sol.sol(t400)[0],0); p1600=np.maximum(sol.sol(t1600)[0],0)
 m400=p400.max();m1600=p1600.max()
 results.append((row.dose,m400,m1600,int(m400>=1e6),int(m1600>=1e6),int(row.infected),abs(m1600-m400)/max(m1600,1e-20)))
a=pd.DataFrame(results,columns=['dose','peak400','peak1600','cross400','cross1600','saved_cross','relative_peak_change'])
s=a.groupby('dose').agg(n=('dose','size'),changes=('cross400',lambda x:0),p400=('cross400','mean'),p1600=('cross1600','mean'),max_rel_peak=('relative_peak_change','max')).reset_index()
s['classification_changes']=a.groupby('dose').apply(lambda z:int((z.cross400!=z.cross1600).sum()),include_groups=False).values
s['saved_disagreements']=a.groupby('dose').apply(lambda z:int((z.cross400!=z.saved_cross).sum()),include_groups=False).values
s.to_csv(Path(__file__).resolve().parents[1] / 'results' / 'time_grid_sensitivity_summary.csv',index=False)
print('N',len(a),'failures',failure,'class changes',int((a.cross400!=a.cross1600).sum()),'saved disagreements',int((a.cross400!=a.saved_cross).sum()),'max rel peak',a.relative_peak_change.max())
print(s.to_string(index=False))

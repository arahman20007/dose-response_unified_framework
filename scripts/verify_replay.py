from pathlib import Path
DATA_DIR = Path(__file__).resolve().parents[1] / "data"
import numpy as np,pandas as pd
from scipy.integrate import solve_ivp
np.random.seed(123)
h=pd.read_excel(DATA_DIR / 'simulation_summary_with_band.xlsx',sheet_name='Virtual_Hosts')
doses=np.logspace(0,8,12); ts=np.linspace(0,30,400); errs=[]; failures=0
for D in doses:
 for i in range(300):
  eta=np.random.lognormal(np.log(1e-3)-.5*.4**2,.4)
  phi=np.random.lognormal(np.log(2e-5)-.5*.25**2,.25)
  row=h.iloc[len(errs)]
  def fun(t,y):
   P,I=np.maximum(y,0)
   return [P*(1*(1-P/1e8)-phi*I-.1),.001*P+.3*I*(1-I/1e5)-.2*I]
  sol=solve_ivp(fun,(0,30),(eta*D,4e4),t_eval=ts,method='RK45',rtol=1e-6,atol=1e-9)
  if not sol.success:failures+=1
  errs.append((abs(eta-row.eta),abs(phi-row.phi),abs(max(sol.y[0].clip(0))-row.max_pathogen_burden)))
a=np.array(errs)
print('independent first repetition 3600 integrations; failures',failures)
print('maximum abs discrepancies eta,phi,Pmax:',a.max(axis=0))
print('max relative Pmax discrepancy',max(a[:,2]/np.maximum(h.iloc[:3600].max_pathogen_burden.to_numpy(),1e-15)))

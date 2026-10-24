"""Plot independent body-flow and moving/static mechanism controls."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fn___plotEvidence import save_hybrid


def main():
    out=Path(__file__).resolve().parents[1]/'reports/2026-10-24-root-crossflow'
    controls=json.loads((out/'motion_root_controls.json').read_text(encoding='utf-8'))['summary']
    body=json.loads((out/'body_crossflow_control.json').read_text(encoding='utf-8'))['rows']
    plt.rcParams.update({'font.family':'Times New Roman','mathtext.fontset':'stix','font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
    fig,axs=plt.subplots(2,3,figsize=(10,6),layout='constrained')
    eta=[r['span'] for r in controls]
    for ax,key,title,ylabel,scale in zip(axs[0],['lower_lateral_velocity_difference','negative_integral','jump'],
            ['Lower-side lateral flow','Upper-front rotation window','Surface pressure difference'],
            [r'$\Delta\overline{W}/U$',r'$\Delta\overline{\Gamma}_{-}/(UL)\;\times 10^3$',r'$\Delta\overline{(C_{p,l}-C_{p,u})}$'],[1,1e3,1]):
        for group,label,color,marker in [('moving','Moving','#D55E00','o'),('static','Frozen','#0072B2','s')]:
            ax.plot(eta,[r[group][key]*scale for r in controls],marker+'-',color=color,label=label,zorder=0,ms=4)
        ax.set(title=title,ylabel=ylabel,xlabel='Span fraction');ax.legend(frameon=False,fontsize=8)
    for ax,side in zip(axs[1,:2],['lower','upper']):
        rows=[r for r in body if r['side']==side];x=[np.mean(r['span']) for r in rows]
        for key,label,color,marker in [('mean_difference','B − I','#D55E00','o'),('mean_body_disturbance','Body only − uniform','#0072B2','s'),('mean_nonadditive','Difference of curves','.4','^')]:
            ax.plot(x,[r[key][2] for r in rows],marker+'-',label=label,color=color,zorder=0,ms=4)
        ax.set(title=side.capitalize()+' side: body-flow control',xlabel='Span-band centre',ylabel=r'$\overline{W}/U$');ax.legend(frameon=False,fontsize=8)
    ax=axs[1,2]
    gates=json.loads((out/'root_gate_flux.json').read_text(encoding='utf-8'))['summary']
    for h,color in [(.005,'#0072B2'),(.008,'#D55E00')]:
        rs=[r for r in gates if r['h0']==h]
        ax.plot([r['target_z'] for r in rs],[100*(1-r['rootward_flux_ratio']) for r in rs],'o-',label=f'h/L ≥ {h:.3f}',color=color,zorder=0)
    ax.set(title='Rootward open-gate flux reduction',xlabel='Gate z/L',ylabel='Reduction (%)',ylim=(0,100));ax.legend(frameon=False,fontsize=8)
    for ax,letter in zip(axs.flat,'abcdef'):
        ax.set_rasterization_zorder(1);ax.spines[['top','right']].set_visible(False)
        ax.text(-.13,1.05,'('+letter+')',transform=ax.transAxes)
        ax.axhline(0,color='.7',lw=.5,zorder=0)
    save_hybrid(fig,out/'26_root_mechanism_controls',600)
    svg=out/'26_root_mechanism_controls.svg'
    svg.write_text('\n'.join(s.rstrip() for s in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')


if __name__=='__main__':main()

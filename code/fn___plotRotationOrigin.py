"""Mechanism diagnostics with fixed spatial coverage and explicit scope."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fn___plotEvidence import style_axes,save_hybrid


def plot_results(output,result):
    plt.rcParams.update({'font.family':'Times New Roman','mathtext.fontset':'stix','font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
    blue='#087bc1';red='#d71920';grey='#92979b';purple='#a11c68'
    rows=result['records'];q=rows[0]['terms']
    fig,axs=plt.subplots(1,3,figsize=(10.8,3.7));fig.subplots_adjust(left=.08,right=.99,bottom=.30,top=.85,wspace=.40)
    for ax,label in zip(axs,'abc'):style_axes(ax,label)
    keys=['velocity_change','vorticity_change','lamb']
    axs[0].bar(range(3),[q[k]['mean']*1000 for k in keys],color=[blue,red,grey],zorder=0)
    axs[0].set_xticks(range(3),['Velocity change','Vorticity change','Net Lamb'],rotation=35,ha='right')
    axs[0].set(ylabel=r'$10^3\overline{\Delta C}_{L,\Lambda}^{\rm pair}$',title='Exact case-midpoint identity')
    keys=['chord_tangent','cross_chord_tangent','body_only_lamb']
    axs[1].bar(range(3),[q[k]['mean']*1000 for k in keys],color=[blue,red,grey],zorder=0)
    axs[1].set_xticks(range(3),['Chord tangent','Cross-chord tangent','Body-only control'],rotation=35,ha='right')
    axs[1].set(ylabel='Paired layer contribution × 1000',title='Local directional components')
    for key,label,color in [('normal_shear_change','Normal derivative of chord velocity',blue),('chord_normal_velocity_change','Chord derivative of normal velocity',red),('velocity_change','Full velocity-change term',grey)]:
        axs[2].plot(range(len(rows)),[r['terms'][key]['mean']*1000 for r in rows],'o-',ms=3,lw=1,color=color,zorder=0,label=label)
    axs[2].set_xticks(range(len(rows)),['.012–.040','.016–.040','.020–.060'],rotation=30,ha='right')
    axs[2].set(xlabel=r'$h_0/L$–$h_1/L$',ylabel='Paired layer contribution × 1000',title='Same 1071 material pairs')
    axs[2].legend(frameon=False,fontsize=6.5)
    save_hybrid(fig,output/'14_rotation_and_gradients',600)
    fig,axs=plt.subplots(1,2,figsize=(8.4,3.5));fig.subplots_adjust(left=.10,right=.98,bottom=.22,top=.84,wspace=.40)
    for ax,label in zip(axs,'ab'):style_axes(ax,label)
    h=[0,.012,.016,.02,.03,.04,.06]
    for prefix,surface,label,color in [('profile_delta_','surface','Body-fin minus isolated',red),('profile_body_','body_only_surface','Body-only sampled control',blue)]:
        values=[q[surface]['mean']]+[q[prefix+f'{x:.3f}']['mean'] for x in h[1:]]
        axs[0].plot(h,np.array(values)*1000,'o-',ms=3,lw=1,color=color,zorder=0,label=label)
    axs[0].set(xlabel=r'Normal offset $h/L$',ylabel='Paired pressure projection × 1000',title='Pressure develops differently near fin')
    axs[0].legend(frameon=False,fontsize=8)
    values=[q['outer']['mean'],q['layer']['mean'],q['wall']['mean']]
    control=[q['body_only_outer']['mean'],q['body_only_layer']['mean'],q['body_only_surface']['mean']-q['body_only_outer']['mean']-q['body_only_layer']['mean']]
    x=np.arange(3);axs[1].bar(x-.18,np.array(values)*1000,width=.36,color=red,zorder=0,label='B-I')
    axs[1].bar(x+.18,np.array(control)*1000,width=.36,color=blue,zorder=0,label='NoFin control')
    axs[1].set_xticks(x,['Outer','Resolved layer','Wall gap'])
    axs[1].set(ylabel='Mean paired contribution × 1000',title='Positive layer change is not inherited')
    axs[1].legend(frameon=False,fontsize=8)
    save_hybrid(fig,output/'15_outer_to_surface_control',600)
    for p in output.glob('*.svg'):p.write_text('\n'.join(x.rstrip() for x in p.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')

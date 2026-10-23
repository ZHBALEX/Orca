"""Figures for distinct pressure/load scopes; never rescale a curve to fit."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fn___plotEvidence import style_axes,save_hybrid


def figures(output,metrics,signals):
    plt.rcParams.update({'font.family':'Times New Roman','mathtext.fontset':'stix','font.size':10,
        'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.7})
    blue='#087bc1';red='#d71920';grey='#92979b';purple='#a11c68';t=signals['phase']
    manifest={}
    fig,axs=plt.subplots(1,2,figsize=(8.4,3.4))
    fig.subplots_adjust(left=.09,right=.98,bottom=.20,top=.85,wspace=.35)
    for ax,panel in zip(axs,'ab'):style_axes(ax,panel)
    axs[0].plot(signals['dense_phase'],signals['system_pressure'],color=blue,lw=1,zorder=0,label='Closed-system B - I - N')
    axs[0].plot(t,signals['regional_pair'],'o-',color=red,ms=3,lw=1,zorder=0,label='Paired common fins B - I')
    axs[0].set(xlabel=r'$t/T$',ylabel=r'$\Delta C_{L,p}$',title='Different spatial scopes')
    axs[0].legend(frameon=False,fontsize=8)
    means=[metrics['system']['case_stats']['interaction'][k]['mean'] for k in ['pressure','viscous','total']]
    axs[1].bar(range(3),np.array(means)*1000,color=[blue,grey,purple],zorder=0)
    axs[1].set_xticks(range(3),['Pressure','Viscous','Total'])
    axs[1].set(ylabel=r'$10^3\overline{\Delta C_L}$',title='481-sample system contrast')
    manifest['09_system_and_pair']=save_hybrid(fig,output/'09_system_and_pair',600)
    fig,axs=plt.subplots(1,2,figsize=(8.4,3.4));fig.subplots_adjust(left=.09,right=.98,bottom=.20,top=.85,wspace=.35)
    for ax,panel in zip(axs,'ab'):style_axes(ax,panel)
    for key,label,c,ls in [('native_reference','Cartesian reference',blue,'-'),('marker_reference','Marker reference',grey,'--'),('pair','Paired jump',red,':')]:
        axs[0].plot(t,signals['regional_'+key],ls,color=c,lw=1.2,zorder=0,label=label)
    axs[0].set(xlabel=r'$t/T$',ylabel=r'$\Delta C_{L,p}^{\rm common}$',title='Same region, explicit reference')
    axs[0].legend(frameon=False,fontsize=8)
    for key,label,c in [('geometry_and_unpaired','Geometry / unpaired / interpolation',red),('discretization','Cartesian minus marker',blue)]:
        axs[1].plot(t,signals['regional_'+key],color=c,lw=1,zorder=0,label=label)
    axs[1].set(xlabel=r'$t/T$',ylabel='Accounting remainder',title='Remainders retained')
    axs[1].legend(frameon=False,fontsize=7)
    manifest['10_pressure_to_load']=save_hybrid(fig,output/'10_pressure_to_load',600)
    fig,axs=plt.subplots(1,3,figsize=(10.8,3.7));fig.subplots_adjust(left=.075,right=.99,bottom=.30,top=.84,wspace=.40)
    for ax,panel in zip(axs,'abc'):style_axes(ax,panel)
    for key,label,c in [('surface','Surface',blue),('outer','Outer',grey),('layer','Resolved layer',red),('wall','Wall gap',purple)]:
        axs[0].plot(t,signals['budget_'+key],color=c,lw=1,zorder=0,label=label)
    axs[0].set(xlabel=r'$t/T$',ylabel=r'$\Delta C_{L,p}^{\rm pair}$',title='One-fin pressure bookkeeping')
    axs[0].legend(frameon=False,fontsize=7)
    keys=['lamb','kinetic','unsteady','viscous','prediction','layer','closure']
    for ax,values,title in [(axs[1],{k:v['mean'] for k,v in metrics['paired_normal_budget']['force_stats'].items()},'Cycle mean'),(axs[2],metrics['paired_normal_budget']['late_force_means'],r'Late window: $0.75\leq t/T\leq1$')]:
        ax.bar(range(7),[values[k]*1000 for k in keys],color=[red,blue,purple,grey,grey,blue,red],zorder=0)
        ax.set_xticks(range(7),['Lamb','Kinetic','Unsteady','Viscous','Rebuilt','Measured','Residual'],rotation=50,ha='right')
        ax.set(ylabel=r'$10^3\Delta C_{L,p}^{\rm pair}$',title=title)
    manifest['11_paired_pressure_budget']=save_hybrid(fig,output/'11_paired_pressure_budget',600)
    for stem in manifest:
        path=output/(stem+'.svg')
        path.write_text('\n'.join(line.rstrip() for line in path.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
    return manifest

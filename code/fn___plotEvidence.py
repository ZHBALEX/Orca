"""Compact scientific figures with raster data and editable labels/axes."""
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def style_axes(ax, panel):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(direction='in', width=.7)
    ax.set_rasterization_zorder(1)
    ax.axhline(0, color='.65', lw=.5, zorder=2)
    ax.text(-.16, 1.06, f'({panel})', transform=ax.transAxes, fontsize=12, fontstyle='italic')


def plot_budget(data, result, palette):
    fig, axs = plt.subplots(2, 2, figsize=(7.5, 6.0))
    fig.subplots_adjust(left=.12, right=.98, bottom=.12, top=.94, wspace=.40, hspace=.65)
    row=data['phases']['upper']; t=row['phase']
    for ax, panel in zip(axs.flat, 'abcd'):
        style_axes(ax, panel)
    ax=axs[0,0]
    ax.plot(t,row['pressure_resolved_measured'],'o-',color=palette['blue'],ms=3,lw=1,zorder=0,label='Measured')
    ax.plot(t,row['pressure_resolved_prediction'],'s--',color=palette['red'],ms=2.5,lw=1,zorder=0,label='Reconstructed')
    ax.set(xlabel=r'$t/T$',ylabel=r'$\Delta C_{p,h}$',title='Resolved upper layer')
    ax.legend(frameon=False,fontsize=8)
    ax=axs[0,1]
    keys=['lamb','kinetic','unsteady','viscous','resolved_prediction','resolved_measured']
    values=[result['regions']['upper']['quantities']['pressure_'+k]['trapezoidal_mean'] for k in keys]
    ax.bar(range(6),values,color=[palette['red']]*4+[palette['grey'],palette['blue']],width=.65,zorder=0)
    ax.set_xticks(range(6),['Lamb','Kinetic','Unsteady','Viscous','Rebuilt','Measured'],rotation=35,ha='right')
    ax.set(ylabel=r'$\overline{\Delta C}_{p,h}$',title='Signed mean balance')
    ax=axs[1,0]
    for method,marker,color in [('central','o',palette['blue']),('fourier','s',palette['red'])]:
        rows=[r for r in result['interval_sensitivity'] if r['method']==method]
        ax.plot(range(len(rows)),[r['unfitted_R2'] for r in rows],marker+'-',color=color,lw=1,ms=3,zorder=0,label=method.capitalize())
    ax.set_xticks(range(5),['.008–.030','.012–.030','.012–.040','.016–.040','.020–.060'],rotation=30,ha='right')
    ax.set(xlabel=r'$h_0/L$–$h_1/L$',ylabel=r'$R^2$ (unfitted)',title='Derivative / interval sensitivity',ylim=(.45,1.02))
    ax.legend(frameon=False,fontsize=8,loc='lower left')
    ax=axs[1,1]
    x=np.arange(4); q=result['regions']['upper']['quantities']
    for shift,key,label,col in [(-.22,'lift_lamb','Lamb',palette['red']),(0,'lift_resolved_measured','Resolved layer',palette['blue']),(.22,'lift_surface','Surface',palette['grey'])]:
        ax.bar(x+shift,1000*np.array(q[key]['signed_quarter_integrals']),width=.22,label=label,color=col,zorder=0)
    ax.set_xticks(x,['0–.25','.25–.50','.50–.75','.75–1'])
    ax.set(xlabel=r'$t/T$ window',ylabel=r'$10^3\int\Delta C_L\,\mathrm{d}(t/T)$',title='Signed contributions by quarter')
    ax.legend(frameon=False,fontsize=7.5,loc='best')
    return fig


def plot_branch(data, result, palette, velocity_scale, length_scale):
    fig,axs=plt.subplots(1,3,figsize=(10,3.1))
    fig.subplots_adjust(left=.065,right=.99,bottom=.21,top=.87,wspace=.40)
    for ax,p in zip(axs,'abc'): style_axes(ax,p)
    rows=data['core']['radial_sensitivity']
    for case,label,color,marker in [('baseline','Body + fin',palette['red'],'o'),('isolated','Only fin',palette['blue'],'s'),('body_only','No fin',palette['grey'],'^')]:
        a=sorted([r for r in rows if r['case']==case and r['radius']==.008],key=lambda r:r['phase'])
        t=[r['phase'] for r in a]
        axs[0].plot(t,[r['gamma']/(velocity_scale*length_scale) for r in a],marker+'-',color=color,lw=1,ms=3,zorder=0,label=label)
        axs[1].plot(t,[r['Cp'] for r in a],marker+'-',color=color,lw=1,ms=3,zorder=0)
    axs[0].set(xlabel=r'$t/T$',ylabel=r'$\Gamma_x/(U_\infty L)$',title='Fixed fluid contour')
    axs[0].legend(frameon=False,fontsize=8)
    axs[1].set(xlabel=r'$t/T$',ylabel=r'$C_{p,\mathrm{core}}-C_{p,\mathrm{ring}}$',title='Core pressure contrast')
    for r in result['same_release_pairs']:
        if r['completed']:
            axs[2].plot([0,1],[r['baseline_start'],r['baseline_end']],'o-',color=palette['red'],ms=3,lw=.7,zorder=0)
            axs[2].plot([0,1],[r['isolated_start'],r['isolated_end']],'s--',color=palette['blue'],ms=3,lw=.7,zorder=0)
    axs[2].set_xticks([0,1],['Release','End'])
    axs[2].set(ylabel=r'$\omega_x$ (solver units)',title='Six completed paired releases',xlim=(-.15,1.15))
    return fig


def save_hybrid(fig, stem, dpi):
    for extension in ('png','pdf','svg'):
        fig.savefig(stem.with_suffix('.'+extension),dpi=dpi,facecolor='white')
    root=ET.parse(stem.with_suffix('.svg')).getroot()
    images=root.findall('.//{http://www.w3.org/2000/svg}image')
    texts=root.findall('.//{http://www.w3.org/2000/svg}text')
    assert images and texts
    assert all(im.attrib.get('{http://www.w3.org/1999/xlink}href','').startswith('data:image/') for im in images)
    # Tick lines use line2d IDs too; semantic data IDs are checked through artist z-order before saving.
    for ax in fig.axes:
        assert all(line.get_zorder()<1 for line in ax.lines if line.get_label() not in ('_child0','_child1') and line.get_color()!='.65')
    plt.close(fig)
    return dict(embedded_images=len(images),editable_text_nodes=len(texts))

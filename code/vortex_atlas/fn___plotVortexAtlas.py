"""Pure plotting functions: arrays and options in, Matplotlib figures out."""
from functools import partial
import numpy as np
from fn___vortexStyle import axes_style as _axes_style, section as _section
from fn___vortexStyle import relative_coordinate as _relative_coordinate
from fn___vortexStyle import panel, bracket, bar


def outline(points,faces,z):
    from scipy.spatial import ConvexHull
    tri=points[faces];hits=[]
    for a,b in [(0,1),(1,2),(2,0)]:
        pa,pb=tri[:,a],tri[:,b];dz=pb[:,2]-pa[:,2]
        ok=(abs(dz)>1e-10)&((pa[:,2]-z)*(pb[:,2]-z)<=0)
        t=(z-pa[ok,2])/dz[ok]
        hits.append(pa[ok,:2]+t[:,None]*(pb[ok,:2]-pa[ok,:2]))
    xy=np.unique(np.round(np.concatenate(hits),8),axis=0)
    h=ConvexHull(xy)
    return xy[np.r_[h.vertices,h.vertices[0]]]


def annotation(ax,name,xy,text,color,fontsize=10):
    ax.plot(*xy,'o',ms=2.8,mfc='white',mec=color,mew=.9,zorder=8)
    ax.annotate(name,xy=xy,xytext=text,ha='center',va='center',fontsize=fontsize,color=color,
        bbox=dict(boxstyle='square,pad=.08',fc='white',ec='none',alpha=.94),
        arrowprops=dict(arrowstyle='-',color=color,lw=.55,shrinkA=2,shrinkB=2),zorder=9)


def core_outline(ax,d,r,normal,index,color,lw=1):
    if 'section_region' not in r:return
    dim=2-normal;remaining=[i for i in range(3) if i!=dim]
    q=np.zeros(np.take(d['Q'],index,axis=dim).shape);ids=np.array(r['section_region']['indices'])
    q[tuple(ids[:,remaining].T)]=1
    axes=[d['z'],d['y'],d['x']]
    ax.contour(axes[remaining[1]],axes[remaining[0]],q,levels=[.5],colors=[color],linewidths=lw,zorder=5)


def plot_overview(rows,volumes,geometry,style,options,plt):
    DISPLAY=style['display_labels'];COLORS=style['label_colors']
    RED=style['palette']['positive'];CYAN=style['palette']['negative_outline']
    ORANGE=style['palette']['fin'];MAGENTA=style['palette']['accent']
    INK=style['palette']['ink'];GREY=style['palette']['grey']
    axes_style=partial(_axes_style,origin=style['origin'],length=style['length'])
    section=partial(_section,style=style)
    relative_coordinate=partial(_relative_coordinate,origin=style['origin'],length=style['length'])
    from scipy.spatial import ConvexHull
    fig=plt.figure(figsize=options['figsize'])
    gs=fig.add_gridspec(2,2,**options['grid'])
    axs=np.array([[fig.add_subplot(gs[i,j]) for j in range(2)] for i in range(2)])
    for col,case in enumerate(options['cases']):
        d=volumes[case];g=geometry[case];rs=[r for r in rows if r['case']==case]
        ax=axs[0,col]
        qp=np.max(np.where(d['valid'],d['Q'],0),axis=1)
        ax.contourf(d['x'],d['z'],qp,levels=options['projection_q_levels'],colors=options['projection_colors'],zorder=0)
        nodes=np.unique(g['faces'][g['pectoral_mask']]);points=g['points_history'][style['frame_index'],nodes]
        points=points[points[:,2]<options['geometry_z_max']][:,[0,2]];h=ConvexHull(points)
        p=points[np.r_[h.vertices,h.vertices[0]]]
        ax.fill(p[:,0],p[:,1],facecolor=ORANGE,alpha=.10,zorder=.5)
        ax.plot(p[:,0],p[:,1],color=ORANGE,lw=.8,zorder=1)
        for r in rs:
            label=r['label']
            if label not in COLORS or label in ['BV','TEV0']:continue
            p=np.array(r['ridge']['points'])
            if label!='TV':ax.plot(p[:,0],p[:,2],color=COLORS[label],lw=1.05,
                ls='--' if label=='W2' else '-',zorder=3)
            xy=np.array(r['seed'])[[0,2]]
            if label in ['LEV','TEV1','W2','TEV3']:
                tx=options['label_x'][case][label]
                text=(tx,options['label_z'] if label!='W2' else options['label_z_tev2'])
            else:text=options['extra_labels'][label]
            name=DISPLAY[label]
            if case=='baseline' and label=='TEV3':name=r'$\mathrm{TEV}_{3}^{*}$'
            annotation(ax,name,xy,text,COLORS[label],9)
        ax.axhline(options['section_z'],color=GREY,lw=.55,ls='--',zorder=1)
        ax.set_xlim(*options['xlim']);ax.set_ylim(*options['projection_ylim'])
        axes_style(ax,options['xticks'],options['projection_yticks'],coordinates=('x','z'));panel(ax,chr(97+col))
        bracket(ax,style['case_labels'][case],y=1.15)
        ax.annotate('',*options['flow_arrow'],arrowprops=dict(arrowstyle='->',lw=.65,color=INK))
        ax.text(*options['flow_label'],r'$U_\infty$',ha='center',va='bottom',fontsize=8)
        ax=axs[1,col];iz=int(np.argmin(abs(d['z']-options['section_z'])))
        im,_=section(ax,d,2,iz,options['color_limit'])
        poly=outline(g['points_history'][style['frame_index']],g['faces'][g['pectoral_mask']],float(d['z'][iz]))
        ax.fill(poly[:,0],poly[:,1],color=ORANGE,zorder=4)
        for r in rs:
            label=r['label']
            if label not in ['LEV','TEV1','W2','TEV3','TEV0']:continue
            xy=r['seed'][:2];name=DISPLAY[label]
            if label=='TEV3' and case=='baseline':name=r'$\mathrm{TEV}_{3}^{*}$'
            ty=options['section_label_top'] if label in ['LEV','W2'] else options['section_label_bottom']
            tx=xy[0]+options['section_label_dx'].get(label,0)
            annotation(ax,name,xy,(tx,ty),COLORS[label],9)
            core_outline(ax,d,r,2,iz,CYAN if r['sign']<0 else RED,1.1)
        ax.set_xlim(*options['xlim']);ax.set_ylim(*options['section_ylim'])
        axes_style(ax,options['xticks'],options['section_yticks']);panel(ax,chr(99+col))
        ax.set_title('$z^*='+f'{relative_coordinate(d["z"][iz],"z"):.5f}'+'$',pad=7,fontsize=9)
    bar(fig,im,options['colorbar_box'],r'$\omega_zT$',options['color_limit'])
    fig.text(.91,.985,'$t^*='+f"{style['phase']:g}"+'$',ha='right',va='top',fontsize=10)
    return fig


def plot_body_controls(rows,volumes,geometry,style,options,plt):
    DISPLAY=style['display_labels'];COLORS=style['label_colors']
    RED=style['palette']['positive'];CYAN=style['palette']['negative_outline']
    ORANGE=style['palette']['fin'];MAGENTA=style['palette']['accent']
    INK=style['palette']['ink'];GREY=style['palette']['grey']
    axes_style=partial(_axes_style,origin=style['origin'],length=style['length'])
    section=partial(_section,style=style)
    relative_coordinate=partial(_relative_coordinate,origin=style['origin'],length=style['length'])
    from matplotlib.patches import Circle
    fig=plt.figure(figsize=options['figsize'])
    gs=fig.add_gridspec(2,3,**options['grid'])
    axs=np.array([[fig.add_subplot(gs[i,j]) for j in range(3)] for i in range(2)])
    rootp=[r for r in rows if r['case']=='baseline' and r['label'] in ['RV+','RV-']]
    bv=next(r for r in rows if r['case']=='baseline' and r['label']=='BV')
    xvalues=[rootp[0]['seed'][0],bv['seed'][0]]
    for col,case in enumerate(options['cases']):
        d=volumes[case]
        for ri in range(2):
            ax=axs[ri,col];j=int(np.argmin(abs(d['x']-xvalues[ri])))
            im,_=section(ax,d,0,j,options['color_limit'])
            if ri==0:
                for r in rootp:
                    xy=(r['seed'][1],r['seed'][2]);color=RED if r['sign']>0 else CYAN
                    ax.add_patch(Circle(xy,options['circle_radius'],fill=False,ec=color,lw=1))
                    if col==0:annotation(ax,DISPLAY[r['label']],xy,tuple(np.asarray(xy)+options['root_label_offset']),COLORS[r['label']],9)
                ax.set_xlim(*options['root_xlim']);ax.set_ylim(*options['root_ylim'])
                axes_style(ax,options['root_xticks'],options['root_yticks'],coordinates=('y','z'))
                bracket(ax,style['case_labels'][case],y=1.20)
            else:
                xy=(bv['seed'][1],bv['seed'][2]);ax.add_patch(Circle(xy,options['circle_radius'],fill=False,ec=MAGENTA,lw=1))
                annotation(ax,DISPLAY['BV'],xy,options['body_label'],INK,9)
                ax.set_xlim(*options['body_xlim']);ax.set_ylim(*options['body_ylim'])
                axes_style(ax,options['body_xticks'],options['body_yticks'],coordinates=('y','z'))
            ax.set_title('$x^*='+f'{relative_coordinate(d["x"][j],"x"):.5f}'+'$',pad=7,fontsize=9)
            panel(ax,chr(97+ri*3+col))
            if col:ax.set_ylabel('')
    bar(fig,im,options['colorbar_box'],r'$\omega_xT$',options['color_limit'])
    fig.text(.91,.99,'$t^*='+f"{style['phase']:g}"+'$',ha='right',va='top',fontsize=10)
    return fig


def plot_individual_cores(rows,volumes,geometry,style,options,plt):
    DISPLAY=style['display_labels'];COLORS=style['label_colors']
    RED=style['palette']['positive'];CYAN=style['palette']['negative_outline']
    ORANGE=style['palette']['fin'];MAGENTA=style['palette']['accent']
    INK=style['palette']['ink'];GREY=style['palette']['grey']
    axes_style=partial(_axes_style,origin=style['origin'],length=style['length'])
    section=partial(_section,style=style)
    relative_coordinate=partial(_relative_coordinate,origin=style['origin'],length=style['length'])
    picks=options['picks']
    fig=plt.figure(figsize=options['figsize'])
    gs=fig.add_gridspec(2,4,**options['grid'])
    axs=np.array([[fig.add_subplot(gs[i,j]) for j in range(4)] for i in range(2)])
    for num,(ax,(case,label)) in enumerate(zip(axs.ravel(),picks)):
        r=next(r for r in rows if r['case']==case and r['label']==label);d=volumes[case];n=r['normal']
        j=int(np.argmin(abs(d[['x','y','z'][n]]-r['seed'][n])));im,remaining=section(ax,d,n,j,options['color_limit'])
        xyz=np.array(r['seed']);xy=xyz[::-1][remaining[::-1]]
        radius=options['radius'] if label!='TV' else options['tip_radius']
        ax.set_xlim(xy[0]-radius,xy[0]+radius);ax.set_ylim(xy[1]-radius,xy[1]+radius)
        color=RED if r['sign']>0 else CYAN
        ax.plot(*xy,'o',ms=2.8,mfc='white',mec=color,mew=1,zorder=8)
        core_outline(ax,d,r,n,j,color,1.1)
        ax.set_title(DISPLAY[label],fontsize=10,pad=7,color=COLORS[label])
        ax.text(.04,.97,r'$\omega_'+['x','y','z'][n]+'T$',transform=ax.transAxes,
            ha='left',va='top',fontsize=9,bbox=dict(fc='white',ec='none',pad=1,alpha=.90),zorder=9)
        ax.locator_params(axis='both',nbins=2)
        axes_style(ax,coordinates=tuple('zyx'[k] for k in remaining[::-1]));panel(ax,chr(97+num))
        ax.tick_params(labelsize=8);ax.set_xlabel(ax.get_xlabel(),fontsize=10);ax.set_ylabel(ax.get_ylabel(),fontsize=10,labelpad=1)
    bar(fig,im,options['colorbar_box'],r'$\omega_nT$',options['color_limit'])
    fig.text(.49,.965,style['case_labels']['isolated'],ha='center',fontsize=11)
    fig.text(.49,.495,style['case_labels']['baseline'],ha='center',fontsize=11)
    fig.text(.91,.985,'$t^*='+f"{style['phase']:g}"+'$',ha='right',va='top',fontsize=10)
    return fig

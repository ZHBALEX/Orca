"""Publication styling guided by the user's JFM reference; no field smoothing."""
import numpy as np


def relative_coordinate(value,name,origin,length):
    """Display positions relative to the fixed Orca reference, scaled by L."""
    return (np.asarray(value)-origin[name])/length

RED='#d71920'
BLUE='#087bc1'
CYAN='#23a9df'
ORANGE='#fc641b'
PEACH='#f8d8bd'
MAGENTA='#a11c68'
INK='#262626'
GREY='#92979b'


def configure(plt,style):
    plt.rcParams.update({'font.family':'serif','font.serif':[style['font_family'],'Times','DejaVu Serif'],
        'mathtext.fontset':'stix','font.size':style['font_size'],'axes.labelsize':11,'axes.titlesize':10,
        'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':9,
        'axes.linewidth':.7,'axes.spines.top':False,'axes.spines.right':False,
        'xtick.direction':'in','ytick.direction':'in','xtick.major.size':3,'ytick.major.size':3,
        'xtick.major.width':.7,'ytick.major.width':.7,'svg.fonttype':'none',
        'pdf.fonttype':42,'ps.fonttype':42,'figure.facecolor':'white','axes.facecolor':'white',
        'savefig.facecolor':'white','axes.unicode_minus':True})


def field_cmap(stops):
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list('User_reference_blue_white_red',
        stops,N=257)


def panel(ax,letter):
    ax.text(-.045,1.035,r'$('+letter+r')$',transform=ax.transAxes,
        ha='left',va='bottom',fontsize=12,clip_on=False)


def axes_style(ax,x_ticks=None,y_ticks=None,arrows=True,coordinates=('x','y'),*,origin,length):
    from matplotlib.ticker import FuncFormatter
    # Data and geometric annotations retain solver coordinates; the affine tick
    # conversion applies to every plotted object without shifting their alignment.
    def formatter(name):
        def label(value,pos):
            value=float(relative_coordinate(value,name,origin,length))
            if abs(value)<1e-10:value=0.
            return f'{value:.6g}'.replace('-', '\N{MINUS SIGN}')
        return FuncFormatter(label)
    ax.xaxis.set_major_formatter(formatter(coordinates[0]))
    ax.yaxis.set_major_formatter(formatter(coordinates[1]))
    ax.set_xlabel('$'+coordinates[0]+'^*$')
    ax.set_ylabel('$'+coordinates[1]+'^*$')
    if x_ticks is not None:ax.set_xticks(x_ticks)
    if y_ticks is not None:ax.set_yticks(y_ticks)
    ax.set_aspect('equal')
    if arrows:
        # Preserve the original reversed X view: +X points left.
        if ax.xaxis_inverted():xy,xytext=(-.035,0),(.01,0)
        else:xy,xytext=(1.035,0),(.99,0)
        ax.annotate('',xy,xytext,xycoords='axes fraction',
            arrowprops=dict(arrowstyle='-|>',lw=.7,color='black',mutation_scale=8),clip_on=False)
        ax.annotate('',(0,1.035),(0,.99),xycoords='axes fraction',
            arrowprops=dict(arrowstyle='-|>',lw=.7,color='black',mutation_scale=8),clip_on=False)


def bracket(ax,text,y=1.16):
    ax.plot([0,0,1,1],[y-.045,y,y,y-.045],color='#555555',lw=.65,
        transform=ax.transAxes,clip_on=False)
    ax.text(.5,y,text,transform=ax.transAxes,ha='center',va='center',fontsize=11,
        bbox=dict(fc='white',ec='none',pad=3),clip_on=False)


def bar(fig,im,box,label,limit):
    ax=fig.add_axes(box)
    cb=fig.colorbar(im,cax=ax,ticks=[-limit,0,limit])
    cb.outline.set_linewidth(.55);cb.ax.tick_params(length=2,width=.55,labelsize=8)
    cb.ax.set_title(label,fontsize=11,pad=7)
    return cb


def rasterize_plot_layers(fig):
    """One raster plot per panel; text, axes and annotation leaders stay editable."""
    from matplotlib.patches import Circle
    for ax in fig.axes:
        artists=list(ax.collections)
        # Solid geometry and colour-bar extensions are plot content. Circles are
        # annotation callouts and remain editable, as do their leaders/markers.
        artists += [p for p in ax.patches if not isinstance(p,Circle)]
        artists += [line for line in ax.lines
                    if line.get_transform()==ax.transData
                    and line.get_marker() in [None,'None','',' ']]
        artists.sort(key=lambda artist:artist.get_zorder())
        for i,artist in enumerate(artists):
            artist.set_rasterized(False)
            artist.set_zorder(-100+i*.01)
        # A single contiguous raster pass combines the scientific layers into
        # an embedded PNG, instead of exporting contours as editable paths.
        ax.set_rasterization_zorder(0)


def section(ax,d,normal,index,vmax=45,*,style):
    dim=2-normal;remaining=[i for i in range(3) if i!=dim]
    axes=[d['z'],d['y'],d['x']];names=['z','y','x']
    a,b=axes[remaining[1]],axes[remaining[0]]
    valid=np.take(d['valid'],index,axis=dim);fluid=np.take(d['fluid'],index,axis=dim)
    w=np.take(d['omega'][...,normal],index,axis=dim);q=np.take(d['Q'],index,axis=dim)
    field=np.ma.masked_where(~valid,w)
    # All scientific plot layers are rasterized together at export.
    # Native-cell contour interpolation is not field smoothing.
    im=ax.contourf(a,b,field,levels=np.linspace(-vmax,vmax,style['contour_levels']),cmap=field_cmap(style['cmap_stops']),
        extend='both',corner_mask=False,antialiased=False,zorder=0)
    ax.set_rasterization_zorder(.5)
    qsafe=np.where(valid,q,0)
    levels=[v for v in style['q_contours'] if np.nanmin(qsafe)<v<np.nanmax(qsafe)]
    if levels:ax.contour(a,b,qsafe,levels=levels,colors='#777777',linewidths=.45,zorder=1)
    if not fluid.all():ax.contourf(a,b,fluid.astype(float),levels=[-.1,.5],colors=[style['palette']['solid']],zorder=2)
    ax.set_xlabel('$'+names[remaining[1]]+'^*$');ax.set_ylabel('$'+names[remaining[0]]+'^*$')
    return im,remaining

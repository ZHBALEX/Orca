"""三张 Orca 涡结构图：只改 main() 用户设置，然后直接运行。"""
from pathlib import Path
import os
import sys


def main():
    # ================= 用户设置：路径、数据 =================
    study_root = Path(r'D:\OneDrive - University of Virginia\Research_UVA\Pj_Orca_3D\test\body_fin_mechanism_20261022')
    output = study_root / 'vortex_atlas_v10'
    identity_file = output / 'identities.json'
    frame_index, phase = 16, 1.0  # 锚点与缓存必须为同一相位
    volume_files = {
        'baseline': study_root / f'vortex_source_v8/volumes/baseline_{frame_index:02d}.npz',
        'isolated': study_root / f'vortex_source_v8/volumes/isolated_{frame_index:02d}.npz',
        'body_only': study_root / f'vortex_discovery_v9/volumes/body_only_{frame_index:02d}.npz',
    }
    geometry_files = {case: study_root / f'loads/{case}.npz' for case in ['baseline', 'isolated']}
    # 标准环境已安装 NumPy、SciPy、Matplotlib 时可设为 []。
    dependency_paths = [
        Path(r'C:\Users\ZHB\Documents\GitHub\PicarGeom\outputs\orca_sweep_table_20261016\deps'),
        Path(r'D:\OneDrive - University of Virginia\Research_UVA\Pj_Orca_3D\test\orca_pod_20261008\dependencies'),
    ]

    # ================= 用户设置：统一风格与导出 =================
    origin = dict(x=10., y=6., z=6.)
    length = 1.
    palette = dict(positive='#d71920', negative='#087bc1', negative_outline='#23a9df',
                   fin='#fc641b', solid='#f8d8bd', accent='#a11c68', ink='#262626', grey='#92979b')
    display_labels = {'LEV':r'$\mathrm{LEV}$', 'TEV1':r'$\mathrm{TEV}_{1}$',
                      'W2':r'$\mathrm{TEV}_{2}$', 'TEV3':r'$\mathrm{TEV}_{3}$',
                      'TEV0':r'$\mathrm{TEV}_{\mathrm{n}}$', 'RV+':r'$\mathrm{RV}^{+}$',
                      'RV-':r'$\mathrm{RV}^{-}$', 'TV':r'$\mathrm{TV}$', 'BV':r'$\mathrm{BV}$'}
    # W2 是既有数据 ID；显示名称是 TEV₂。
    label_colors = {key: palette[color] for key, color in {
        'LEV':'negative', 'TEV1':'positive', 'W2':'negative', 'TEV3':'positive',
        'TEV0':'positive', 'RV+':'positive', 'RV-':'negative', 'TV':'accent', 'BV':'ink'}.items()}
    style = dict(
        origin=origin, length=length, phase=phase, frame_index=frame_index,
        font_family='Times New Roman', font_size=10,
        palette=palette, display_labels=display_labels, label_colors=label_colors,
        case_labels={'baseline':'Body + fin', 'isolated':'Only fin', 'body_only':'No fin'},
        q_contours=[50, 250], contour_levels=65,
        cmap_stops=[(0,'#203e8e'),(.20,'#228ec8'),(.36,'#a4d5e9'),(.5,'#ffffff'),
                    (.64,'#f6b9bb'),(.80,'#ec5355'),(1,'#c9101b')],
    )
    export_style = dict(dpi=600, formats=('png', 'pdf', 'svg'), padding=.045)
    selected_figures = ['overview', 'controls', 'cores']
    show_figure = False

    # 所有位置、范围与刻度列表使用原始求解器坐标；图上统一显示星号相对坐标。
    overview = dict(
        name='01_major_vortex_atlas', figsize=(7.5,5.7), cases=['isolated','baseline'],
        grid=dict(left=.085,right=.90,bottom=.085,top=.90,hspace=.36,wspace=.27,height_ratios=[1.42,1]),
        color_limit=45, colorbar_box=[.94,.108,.013,.235],
        xlim=(10.48,10.18), xticks=[10.2,10.3,10.4],
        projection_ylim=(5.733,5.972), projection_yticks=[5.75,5.85,5.95],
        section_ylim=(5.830,5.977), section_yticks=[5.84,5.90,5.96], section_z=5.86,
        geometry_z_max=5.98, projection_q_levels=[50,250,1e8],
        projection_colors=['#eceef0','#d4d9dd'],
        label_x={'isolated':dict(LEV=10.285,TEV1=10.433,W2=10.380,TEV3=10.333),
                 'baseline':dict(LEV=10.295,TEV1=10.453,W2=10.398,TEV3=10.347)},
        label_z=5.745, label_z_tev2=5.753,
        extra_labels={'TV':(10.232,5.782),'RV+':(10.410,5.950),'RV-':(10.405,5.926)},
        section_label_top=5.965, section_label_bottom=5.843,
        section_label_dx={'TEV0':-.012,'TEV3':.006},
        flow_arrow=((10.44,5.963),(10.40,5.963)), flow_label=(10.420,5.968),
    )
    controls = dict(
        name='02_root_and_body_controls', figsize=(7.5,5.6), cases=['baseline','isolated','body_only'],
        grid=dict(left=.075,right=.90,bottom=.09,top=.88,hspace=.32,wspace=.32),
        color_limit=50, colorbar_box=[.94,.29,.013,.34], circle_radius=.008,
        root_xlim=(5.89,5.98), root_ylim=(5.87,5.965),
        root_xticks=[5.90,5.94,5.98], root_yticks=[5.88,5.92,5.96],
        body_xlim=(5.985,6.043), body_ylim=(5.865,5.942),
        body_xticks=[6.,6.02,6.04], body_yticks=[5.88,5.91,5.94],
        root_label_offset=(-.011,.024), body_label=(6.014,5.929),
    )
    cores = dict(
        name='03_individual_core_sections', figsize=(7.5,4.55),
        grid=dict(left=.075,right=.91,bottom=.10,top=.87,hspace=.52,wspace=.54),
        color_limit=65, colorbar_box=[.945,.235,.012,.36], radius=.025, tip_radius=.020,
        picks=[('isolated','LEV'),('isolated','TEV1'),('isolated','W2'),('isolated','TEV3'),
               ('baseline','RV+'),('baseline','RV-'),('baseline','TV'),('baseline','BV')],
    )
    # ================= 执行：读数据 → 绘图 → 导出 =================
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    import numpy  # 在可选依赖路径前加载当前环境的 NumPy
    for path in dependency_paths:
        if path.is_dir():
            sys.path.insert(0, str(path))
    import matplotlib
    if not show_figure:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from read_vortex_atlas import read_vortex_atlas
    from fn___vortexStyle import configure
    from fn___plotVortexAtlas import plot_overview, plot_body_controls, plot_individual_cores
    from fn___exportVortexAtlas import save_figure, file_hashes, write_run_metadata

    inputs = [identity_file, *volume_files.values(), *geometry_files.values()]
    before = file_hashes(inputs)
    rows, volumes, geometry = read_vortex_atlas(identity_file, volume_files, geometry_files, phase)
    configure(plt, style)
    jobs = {'overview':(plot_overview,overview), 'controls':(plot_body_controls,controls),
            'cores':(plot_individual_cores,cores)}
    for key in selected_figures:
        builder, options = jobs[key]
        figure = builder(rows, volumes, geometry, style, options, plt)
        save_figure(figure, output/'figures', options['name'], export_style)
        if not show_figure:
            plt.close(figure)
    if file_hashes(inputs) != before:
        raise RuntimeError('Input files changed while rendering.')
    write_run_metadata(output, style, export_style, before)
    if show_figure:
        plt.show()
    print(f'Exported {len(selected_figures)} figures to {output / "figures"}; input hashes unchanged.', flush=True)


if __name__ == '__main__':
    main()

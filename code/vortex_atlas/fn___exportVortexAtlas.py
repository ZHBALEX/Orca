"""Figure exports and provenance checks; no case-specific input paths."""
from pathlib import Path
import hashlib
import json
from fn___vortexStyle import rasterize_plot_layers


def write_json(path,data):
    Path(path).write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf-8')


def file_hashes(paths):
    return {str(path):hashlib.sha256(Path(path).read_bytes()).hexdigest() for path in paths}


def save_figure(fig,output,name,options):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    rasterize_plot_layers(fig)
    for ext in options['formats']:
        fig.savefig(output/f'{name}.{ext}',dpi=options['dpi'],bbox_inches='tight',pad_inches=options['padding'])
    if 'svg' not in options['formats']:return
    import xml.etree.ElementTree as ET
    root=ET.parse(output/f'{name}.svg').getroot()
    ns={'s':'http://www.w3.org/2000/svg'}
    images=root.findall('.//s:image',ns);texts=root.findall('.//s:text',ns)
    groups=[el.get('id','') for el in root.findall('.//s:g',ns)]
    vector_fields=[g for g in groups if any(s in g for s in ['ContourSet','PolyCollection','QuadMesh'])]
    assert len(images)==len(fig.axes),(name,len(images),len(fig.axes))
    assert texts and not vector_fields,(name,vector_fields)
    assert all(el.get('{http://www.w3.org/1999/xlink}href','').startswith('data:image/png;base64,') for el in images)
    write_json(output/f'{name}.svg_audit.json',dict(embedded_pngs=len(images),
        editable_text_elements=len(texts),vector_field_groups=vector_fields,
        expected_plot_panels_including_colorbar=len(fig.axes),
        editable_layers=['text','ticks','axes','arrows','annotation circles and point markers','case brackets'],
        raster_layers=['colour fields','Q contours','core contours','core ridges','solid geometry','colour-bar fills']))


def write_run_metadata(output,style,export_style,input_hashes):
    write_json(Path(output)/'figure_style.json',dict(
        runner='run__vortexAtlas.py',font=style['font_family'],palette=style['palette'],
        exports=export_style,native_field_smoothing=False,
        filled_contours='Linear native-cell contour interpolation; no Gaussian filtering',
        unchanged_input_sha256=input_hashes,phase=style['phase'],Q_contours=style['q_contours'],
        coordinate_directions_preserved=True,
        display_coordinates=dict(origin=style['origin'],length=style['length'],axis_labels=['x*','y*','z*']),
        raster_layers=['colour fields','Q contours','core contours and ridges','solid geometry','colour-bar fills']))

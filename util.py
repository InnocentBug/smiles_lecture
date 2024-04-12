from urllib.parse import quote

import requests
from IPython.display import SVG
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D


def render_svg(svg):
    """Render SVG images in the Notebook"""
    try:
        svg_string = svg.decode("utf-8")
    except AttributeError:
        svg_string = svg
    svg_string = svg_string.replace("svg:", "")
    return SVG(svg_string)


def moltosvg(mol, molSize=(450, 150), kekulize=False):
    """Generate a SVG stick representation of molecule."""
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except Exception as exc:
            mc = Chem.Mol(mol.ToBinary())
            print(exc)
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg


def render_mol(mol, **kwargs):
    return render_svg(moltosvg(mol, **kwargs))


def get_iupac_name(smiles):
    url = f"https://cactus.nci.nih.gov/chemical/structure/{quote(smiles)}/iupac_name"
    response = requests.get(url, timeout=1)
    return response.text

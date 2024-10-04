from bottle import Bottle, template, request, static_file
from rdkit import Chem
from rdkit.Chem import AllChem
import lzstring

app = Bottle()

@app.route('/static/<filename:path>')
def send_static(filename):
    return static_file(filename, root='./static')

@app.route('/')
def index():
  #return template('index.html')
  return '<b>Hello World</b>!'

@app.route('/convert', method='POST')
def convert_smiles():
    smiles = request.forms.get('smiles')

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)

    xyz_data = []
    for atom in mol.GetAtoms():
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        xyz_data.append(f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}")

    atom_count = mol.GetNumAtoms()
    xmol_data = f"{atom_count}\n{smiles}\n" + "\n".join(xyz_data)

    compressor = lzstring.LZString()
    compressed_data = compressor.compressToEncodedURIComponent(xmol_data)

    url = f"https://poke.yamnor.me/#{compressed_data}"

    return template('result.html', xmol_data=xmol_data, url=url)

# if __name__ == "__main__":
#   app.run(host='localhost', port=8080)

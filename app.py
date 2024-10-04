from bottle import Bottle, template, request, static_file
from rdkit import Chem
from rdkit.Chem import AllChem
import lzstring
import os

app = Bottle()

# Tailwind CSSの静的ファイルを読み込むルート
@app.route('/static/<filename:path>')
def send_static(filename):
    return static_file(filename, root='./static')

# メインページのルート
@app.route('/')
def index():
    return template('index.html')

# SMILES入力を受け取り、XYZ形式に変換し最適化を行い、圧縮URLを生成する処理
@app.route('/convert', method='POST')
def convert_smiles():
    smiles = request.forms.get('smiles')
    
    # SMILESから分子を生成
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)

    # XYZ座標を生成
    xyz_data = []
    for atom in mol.GetAtoms():
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        xyz_data.append(f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}")
    
    # XYZ形式をxmol形式に整形
    atom_count = mol.GetNumAtoms()
    xmol_data = f"{atom_count}\n{smiles}\n" + "\n".join(xyz_data)
    
    # lz-stringで圧縮
    compressor = lzstring.LZString()
    compressed_data = compressor.compressToEncodedURIComponent(xmol_data)

    # 圧縮データを使ったURL生成
    url = f"https://poke.yamnor.me/#{compressed_data}"
    
    return template('result.html', xmol_data=xmol_data, url=url)

# サーバーの起動
# if __name__ == "__main__":
#     app.run(host='localhost', port=8080)

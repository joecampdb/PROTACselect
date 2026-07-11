"""
Assembly engine: join warhead + linker + E3 at their tagged exit vectors.

Convention (see catalog.py):
    warhead carries [*:1], E3 carries [*:2], linker carries both.
RDKit molzip fuses dummy atoms that share an atom-map label, removing the
dummies and forming a single bond between their heavy-atom neighbours:
    warhead[*:1] + [*:1]linker[*:2] + [*:2]E3  ->  warhead-linker-E3
"""
from rdkit import Chem


def _mol(smiles, name):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f"unparseable SMILES for {name!r}: {smiles}")
    return m


def assemble(warhead_smi, linker_smi, e3_smi, sanitize=True):
    """Return the assembled PROTAC as an RDKit Mol (dummies removed)."""
    wh = _mol(warhead_smi, "warhead")
    lk = _mol(linker_smi, "linker")
    e3 = _mol(e3_smi, "e3")

    combined = Chem.CombineMols(Chem.CombineMols(wh, lk), e3)
    zipped = Chem.molzip(combined)          # matches [*:1]<->[*:1], [*:2]<->[*:2]

    # strip any leftover atom-map numbers so descriptors are clean
    for atom in zipped.GetAtoms():
        atom.SetAtomMapNum(0)
    if sanitize:
        Chem.SanitizeMol(zipped)
    return zipped


def assembled_smiles(warhead_smi, linker_smi, e3_smi):
    return Chem.MolToSmiles(assemble(warhead_smi, linker_smi, e3_smi))


def validate_catalog(warheads, e3_ligands, linkers):
    """Parse every catalog SMILES; return list of (kind, key, error) failures."""
    failures = []
    for k, v in warheads.items():
        if Chem.MolFromSmiles(v[2]) is None:
            failures.append(("warhead", k, v[2]))
    for k, v in e3_ligands.items():
        if Chem.MolFromSmiles(v[3]) is None:
            failures.append(("e3", k, v[3]))
    for k, v in linkers.items():
        if Chem.MolFromSmiles(v[2]) is None:
            failures.append(("linker", k, v[2]))
    return failures

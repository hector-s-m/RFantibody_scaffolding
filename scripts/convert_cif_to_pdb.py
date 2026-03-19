#!/usr/bin/env python3
"""
Convert mmCIF files to PDB format.

Lightweight converter that reads ATOM/HETATM records from mmCIF and writes
standard PDB format. Handles the Boltz2 mmCIF output structure.

Usage:
    python scripts/convert_cif_to_pdb.py input.cif output.pdb
    python scripts/convert_cif_to_pdb.py --dir predictions/ --delete-cif
"""

import argparse
import sys
from pathlib import Path


def cif_to_pdb(cif_path: str, pdb_path: str) -> bool:
    """Convert a single mmCIF file to PDB format.

    Parses _atom_site records and writes standard PDB ATOM lines.
    Returns True on success, False on failure.
    """
    try:
        atoms = []
        current_chain = None

        with open(cif_path) as f:
            in_atom_site = False
            columns = []
            for line in f:
                line = line.strip()

                # Detect start of _atom_site loop
                if line.startswith('loop_'):
                    in_atom_site = False
                    columns = []
                    continue
                if line.startswith('_atom_site.'):
                    col_name = line.split('.')[1].strip()
                    columns.append(col_name)
                    in_atom_site = True
                    continue
                if in_atom_site and not line.startswith('_') and line and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) < len(columns):
                        in_atom_site = False
                        continue

                    record = dict(zip(columns, parts))

                    group = record.get('group_PDB', 'ATOM')
                    if group not in ('ATOM', 'HETATM'):
                        continue

                    atom_id = int(record.get('id', 0))
                    atom_name = record.get('label_atom_id', 'X')
                    alt_loc = record.get('label_alt_id', '.')
                    if alt_loc == '.':
                        alt_loc = ' '
                    res_name = record.get('label_comp_id', 'UNK')
                    chain = record.get('label_asym_id', 'A')
                    res_seq = record.get('label_seq_id', '0')
                    if res_seq == '.':
                        res_seq = record.get('auth_seq_id', '0')
                    icode = record.get('pdbx_PDB_ins_code', ' ')
                    if icode == '?' or icode == '.':
                        icode = ' '
                    x = float(record.get('Cartn_x', 0))
                    y = float(record.get('Cartn_y', 0))
                    z = float(record.get('Cartn_z', 0))
                    occ = float(record.get('occupancy', 1.0))
                    bfac = float(record.get('B_iso_or_equiv', 0.0))
                    element = record.get('type_symbol', atom_name[0])

                    # Write TER between chains
                    if current_chain is not None and chain != current_chain:
                        atoms.append('TER\n')
                    current_chain = chain

                    # Format atom name (left-justify 2-char elements, otherwise center)
                    if len(atom_name) < 4:
                        atom_name_fmt = f' {atom_name:<3s}'
                    else:
                        atom_name_fmt = f'{atom_name:<4s}'

                    pdb_line = (
                        f'{group:<6s}'
                        f'{atom_id:>5d} '
                        f'{atom_name_fmt:4s}'
                        f'{alt_loc:1s}'
                        f'{res_name:>3s} '
                        f'{chain:1s}'
                        f'{int(res_seq):>4d}'
                        f'{icode:1s}   '
                        f'{x:8.3f}'
                        f'{y:8.3f}'
                        f'{z:8.3f}'
                        f'{occ:6.2f}'
                        f'{bfac:6.2f}'
                        f'          '
                        f'{element:>2s}  '
                        '\n'
                    )
                    atoms.append(pdb_line)

                elif in_atom_site and (line.startswith('_') or line.startswith('#') or not line):
                    in_atom_site = False

        if not atoms:
            print(f'  Warning: No ATOM records found in {cif_path}')
            return False

        atoms.append('TER\n')
        atoms.append('END\n')

        with open(pdb_path, 'w') as f:
            f.writelines(atoms)

        return True

    except Exception as e:
        print(f'  Error converting {cif_path}: {e}')
        return False


def main():
    parser = argparse.ArgumentParser(description='Convert mmCIF to PDB format')
    parser.add_argument('input', nargs='?', help='Input CIF file (single-file mode)')
    parser.add_argument('output', nargs='?', help='Output PDB file (single-file mode)')
    parser.add_argument('--dir', type=str, help='Directory to batch-convert all .cif files')
    parser.add_argument('--delete-cif', action='store_true', help='Delete CIF files after conversion')
    args = parser.parse_args()

    if args.dir:
        cif_dir = Path(args.dir)
        cif_files = sorted(cif_dir.rglob('*.cif'))
        if not cif_files:
            print(f'No .cif files found in {cif_dir}')
            return

        print(f'Converting {len(cif_files)} CIF files to PDB...')
        converted = 0
        for cif_path in cif_files:
            pdb_path = cif_path.with_suffix('.pdb')
            if cif_to_pdb(str(cif_path), str(pdb_path)):
                converted += 1
                if args.delete_cif:
                    cif_path.unlink()

        print(f'Converted {converted}/{len(cif_files)} files')
        if args.delete_cif:
            print(f'Deleted {converted} CIF files')

    elif args.input and args.output:
        success = cif_to_pdb(args.input, args.output)
        sys.exit(0 if success else 1)

    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()

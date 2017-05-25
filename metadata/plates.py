from __future__ import print_function
import datetime
import os
import re


COLS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
ROWS = ["A", "B", "C", "D", "E", "F", "G", "H"]


def iter_plates(filename):
    '''Parse plate file in 96 well plate format, yielding name, cells:

    ('plate1', {'A': {1: 'a1', ...}, ...}
    '''
    # line string to list of cells
    parse_cells = lambda l: list(map(str.strip, l.rstrip('\r\n').split('\t')))
    with open(filename) as fh:
        # This while loop can be read as "for plate in readfile"
        while True:
            line = fh.readline()
            if not line:
                # EOF
                break
            if not line.strip():
                # Skip till the plate. Empty lines delimit plates
                continue

            cells = parse_cells(line)

            if not cells[0]:
                # Skip lines with empty first cell, which are extra comments or
                # headers or somesuch.
                continue

            # Start of this plate
            plate = {}

            # Check header
            plate_name = cells[0]
            try:
                numerics = list(map(int, cells[1:]))
                assert(numerics == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
            except (ValueError, AssertionError):
                raise ValueError("Bad plate header: ", line)

            # Parse each row into plate dict
            for alpha in "ABCDEFGH":
                line = fh.readline()
                cells = parse_cells(line)
                row = cells[0].upper()
                if row != alpha:
                    raise ValueError("Bad plate line, row label is incorrect:",
                                     line)
                plate[row] = {}
                for num in range(1, 13):
                    plate[row][num] = cells[num]
            yield plate_name, plate


def parse_barcode_plate(platefile):
    barcode_plates = {}
    for plate_name, plate in iter_plates(platefile):
        barcodes = {}
        is_combo = False

        for row in ROWS:
            for col in COLS:
                bcds = plate[row][col].split(' / ')
                if len(bcds) == 2:
                    is_combo = True
                plate[row][col] = bcds
        barcode_plates[plate_name] = {
            'barcodes': plate,
            'is_combo': is_combo,
        }
    return barcode_plates


def sanitise_name(name):
    # make 'empty' and 'blank' actually empty
    for empty in ['empty', 'blank']:
        if name.lower().startswith(empty):
            return ''

    # Remove chars that are bad file names
    badchars = '/\\~`!@#$%^&*()+=<>,.?;:\'"'
    transtable = {c: '-' for c in badchars}
    name = name.translate(transtable)

    # the above can result in multiple '-'s so we reduce them to just 1
    return re.sub('-+', '-', name)


def make_axe(samplefile, barcodefile, outputfile):
    bcdplates = parse_barcode_plate(barcodefile)
    today = datetime.date.today().isoformat()
    is_combo = None
    with open(outputfile, 'w') as ofh:
        print('# Generated', today, file=ofh)
        print('# From', os.path.basename(samplefile), file=ofh)
        print('#', file=ofh)
        for pname, plate in iter_plates(samplefile):
            print('# Start of', pname, file=ofh)
            bcd = bcdplates[pname]
            for row in ROWS:
                for col in COLS:
                    # For single barcode runs, the len(barcodes) is 1 so there
                    # is no tab added.
                    bcdstr = '\t'.join(bcd["barcodes"][row][col])
                    name = sanitise_name(plate[row][col])
                    if name:
                        print(bcdstr, name, sep='\t', file=ofh)

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  brenda_smiles.py
#  
#  Copyright 2019 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import os
import requests
import logging
from bs4 import BeautifulSoup

from libs import paths
from libs.configuration import ENV
from external import osra as osra_service
from external import molconvert as molconvert_service

logger = logging.getLogger(__name__)


def get_substrates_images(brentries, alternate=False):
    id_headers = {'Referer': 'https://www.brenda-enzymes.org/advanced.php',
                  'User-Agent': 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:63.0) Gecko/20100101 Firefox/63.0'}
    img_headers = {'Referer': 'https://www.brenda-enzymes.org/Mol/Mol.image.php',
                   'User-Agent': 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:63.0) Gecko/20100101 Firefox/63.0'}
    id_response = requests.get("https://www.brenda-enzymes.org/advanced.php", headers=id_headers)
    img_response = requests.get("https://www.brenda-enzymes.org/Mol/Mol.image.php", headers=img_headers)
    if id_response.status_code != 200 and img_response.status_code != 200:
        logging.error("Connection to www.brenda-enzymes.org ended with an unexpected error.")
        raise ConnectionError
    substrates = set()
    if alternate:
        substrates.add(brentries)
    else:
        for brentry in brentries:
            for substrate in brentry.substrates:
                substrates.add(substrate)
    for substrate in substrates:
        try:
            mol_id = _get_molecule_id(substrate, id_response, id_headers)
            if mol_id > 0:
                _download_image(substrate, mol_id, img_response, img_headers)
        except FileNotFoundError:
            msg = "The substrate with name: '%s' has no structure image in BRENDA" % substrate
            logging.error(msg)


def _get_molecule_id(name, response, headers):
    # outfolder = ENV["cache_directory"] + "images/"
    outfolder = "./cache/images/"
    testname = name.replace('/', '%')
    testname = testname.replace(' ', '_')
    if os.path.isfile('{OUTF}{NAME}.png'.format(OUTF=outfolder, NAME=testname)):
        return 0
    name = name.replace('(', '%28')
    name = name.replace(')', '%29')
    name = name.replace('+', '%2B')
    name = name.replace(',', '%2C')
    name = name.replace('/', '%2F')
    name = name.replace(' ', '+')
    url = "https://www.brenda-enzymes.org/advanced.php?show=results&CheckOS=1&StypeOS=1&SvalueOS=&Check1%5B0%5D" \
          "=1&Stype1%5B0%5D=1&Svalue1%5B0%5D=&Check1%5B1%5D=1&Stype1%5B1%5D=1&Svalue1%5B1%5D=&CheckEN" \
          "=1&StypeEN=1&SvalueEN=&Check2%5B0%5D=1&Subitem2%5B0%5D=11&Stype2%5B0%5D=1&Svalue2%5B0%5D={SUBSTRATE}" \
          "&Check2%5B1%5D=1&Subitem2%5B1%5D=&Stype2%5B1%5D=1&Svalue2%5B1%5D=&Check2%5B2%5D=1&Subitem2%5B2%5D" \
          "=&Stype2%5B2%5D=1&Svalue2%5B2%5D=&Check2%5B3%5D=1&Subitem2%5B3%5D=&Stype2%5B3%5D=1&Svalue2%5B3%5D" \
          "=&Check3%5B0%5D=1&Subitem3%5B0%5D=&Stype3%5B0%5D=1&Svalue3%5B0%5D=&Check3%5B1%5D=1&Subitem3%5B1%5D" \
          "=&Stype3%5B1%5D=1&Svalue3%5B1%5D=&Check3%5B2%5D=1&Subitem3%5B2%5D=&Stype3%5B2%5D=1&Svalue3%5B2%5D" \
          "=&Check3%5B3%5D=1&Subitem3%5B3%5D=&Stype3%5B3%5D=1&Svalue3%5B3%5D=&OF%5B0%5D=&CheckOF%5B0%5D=1&OF%5B1%5D" \
          "=&CheckOF%5B1%5D=1&OF%5B2%5D=&CheckOF%5B2%5D=1&OF%5B3%5D=&CheckOF%5B3%5D=1&OF%5B4%5D=&CheckOF%5B4%5D" \
          "=1&OF%5B5%5D=&CheckOF%5B5%5D=1&OF%5B6%5D=&CheckOF%5B6%5D=1&OF%5B7%5D=&CheckOF%5B7%5D=1&sfields2=2&S12" \
          "=SEARCH".format(SUBSTRATE=name)
    cookies = {"PHPSESSID": response.cookies["SERVERID"], "SERVERID": response.cookies["SERVERID"]}
    webpage = requests.get(url, headers=headers)
    soup = BeautifulSoup(webpage.text, "lxml")
    search = soup.find("h2")
    if search and (search.string == "No results!"):
        raise FileNotFoundError
    cells = soup.find_all("div", class_="cell")
    idnumber = 0
    for cell in cells:
        if cell.a:
            idnumber = int(cell.a.attrs['href'][22:-2])
    return idnumber


def _download_image(name, mol_id, response, headers):
    # outfolder = ENV["cache_directory"] + "images/"
    outfolder = "./cache/images/"
    name = name.replace('/', '%')
    name = name.replace(' ', '_')
    outfile = '{OUTF}{NAME}.png'.format(OUTF=outfolder, NAME=name)
    if not os.path.isfile(outfile):
        url = "https://www.brenda-enzymes.org/Mol/Mol.image.php?ID={ID}".format(ID=mol_id)
        cookies = {"PHPSESSID": response.cookies["SERVERID"], "SERVERID": response.cookies["SERVERID"]}
        webpage = requests.get(url, headers=headers, cookies=cookies)
        with open(outfile, 'wb') as stream:
            stream.write(webpage.content)


def get_substrates_smiles(brentries):
    pair_subs_smi = {}
    substrates = set()
    for brentry in brentries:
        for substrate in brentry.substrates:
            substrates.add(substrate)
    infolder = ENV["cache_directory"] + "images/"
    for substrate in substrates:
        name = substrate.replace('/', '%')
        name = name.replace(' ', '_')
        name = "{}{}.png".format(infolder, name)
        if os.path.getsize(name) > 0:
            osra_service.execute(name)
            with open(os.path.join(paths.temp(), "smiles.log"), 'r') as stream:
                smiles = stream.readline().rstrip()
            if len(smiles) > 0:
                pair_subs_smi[substrate] = smiles
    for brentry in brentries:
        for substrate in brentry.substrates:
            if substrate in pair_subs_smi.keys():
                brentry.substrates_smiles[substrate] = pair_subs_smi[substrate]


def generate_images_from_smiles(brentries):
    commands = ""
    for brentry in brentries:
        for subs, smi in brentry.substrates_smiles.items():
            if '*' in smi:
                logger.warning("Found * in smiles %s from image %s.png, you should consider to check the smiles "
                               "conversion manually." % (smi, name))
            name = subs.replace('/', '%')
            name = name.replace(' ', '_')
            name = ENV["cache_directory"] + "tempimages/" + name + ".png"
            # molconvert_service.execute(smi, name)
            commands += '%s "png:w500" -s "%s" -o "%s"\n' %(ENV["molconvert"], smi, name)
    with open("commands.sh", 'w') as script:
        script.write(commands)
    os.system("sh commands.sh")












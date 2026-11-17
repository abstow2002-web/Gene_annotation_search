#!/usr/bin/env python3
"""
Functional Gene Search Script
==============================
Searches MetaCerberus annotation summary TSV files and PICRUSt2 pathway
abundance data for functional genes of interest:

  - Carbon cycling (fixation, methane, degradation, fermentation)
  - Nitrogen cycling (fixation, nitrification, denitrification, DNRA, anammox)
  - Phosphorus cycling (transport, mineralization, polyphosphate)
  - Antibiotic resistance (efflux pumps, beta-lactamases, target modification)
  - Biofilm formation (EPS, quorum sensing, motility, c-di-GMP signalling)

Inputs  : final_annotation_summary_*.tsv  (MetaCerberus output)
          pathway_abundance.qza            (PICRUSt2 QIIME2 artefact)
Output  : gene_search_results_<timestamp>.xlsx
"""

import os
import sys
import glob
import zipfile
import tempfile
import re
import warnings
from datetime import datetime
from collections import defaultdict

import pandas as pd

warnings.filterwarnings("ignore", category=UserWarning)

# ── optional dependency for reading HDF5 BIOM files ──────────────────────────
try:
    import h5py
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False

# =============================================================================
# CURATED GENE LISTS
# Each category holds three search keys:
#   ko_ids           – dict  KO-ID  → plain-English description
#   gene_names       – list  of gene-name substrings (case-insensitive)
#   product_keywords – list  of product-description keywords (case-insensitive)
# =============================================================================

GENE_CATEGORIES = {

    # ── CARBON CYCLING ────────────────────────────────────────────────────────
    "Carbon_Cycling": {
        "description": (
            "Carbon cycling: carbon fixation (Calvin/Wood-Ljungdahl), "
            "methanogenesis, methane oxidation, organic-carbon degradation"
        ),
        "ko_ids": {
            # Carbon fixation
            "K01601": "rbcL – RuBisCO large subunit (Calvin cycle)",
            "K01602": "rbcS – RuBisCO small subunit (Calvin cycle)",
            "K03518": "coxL – CO dehydrogenase large subunit",
            "K03519": "coxM – CO dehydrogenase medium subunit",
            "K03520": "coxS – CO dehydrogenase small subunit",
            "K00196": "acsA – ACS large subunit (Wood-Ljungdahl)",
            "K14138": "acsB – acetyl-CoA synthase (Wood-Ljungdahl)",
            "K00198": "cooS – CO dehydrogenase catalytic subunit",
            "K14092": "fthfs – formate-tetrahydrofolate ligase",
            "K18209": "CODH – carbon monoxide dehydrogenase",
            # Methanogenesis
            "K00399": "mcrA – methyl-CoM reductase α (methanogenesis)",
            "K00401": "mcrB – methyl-CoM reductase β (methanogenesis)",
            "K00402": "mcrG – methyl-CoM reductase γ (methanogenesis)",
            "K00200": "fwdA – formylmethanofuran dehydrogenase subunit A",
            "K00201": "fwdB – formylmethanofuran dehydrogenase subunit B",
            "K00202": "fwdC – formylmethanofuran dehydrogenase subunit C",
            "K00203": "fwdD – formylmethanofuran dehydrogenase subunit D",
            "K01499": "mch – methenyl-H4MPT cyclohydrolase",
            "K00319": "mtd – methylene-H4MPT dehydrogenase",
            "K00320": "mer – methylene-H4MPT reductase",
            "K14083": "mtrA – tetrahydromethanopterin S-methyltransferase A",
            "K04480": "mtrH – tetrahydromethanopterin S-methyltransferase H",
            # Methane oxidation
            "K10944": "pmoA/amoA – particulate methane/ammonia monooxygenase A",
            "K10945": "pmoB/amoB – particulate methane/ammonia monooxygenase B",
            "K10946": "pmoC/amoC – particulate methane/ammonia monooxygenase C",
            "K16157": "mmoX – soluble methane monooxygenase α",
            "K16158": "mmoY – soluble methane monooxygenase β",
            "K16159": "mmoZ – soluble methane monooxygenase γ",
            "K16161": "mmoB – soluble methane monooxygenase regulatory",
            "K16162": "mmoC – soluble methane monooxygenase reductase",
            # Organic carbon degradation
            "K01179": "endA – endoglucanase (cellulose)",
            "K01180": "bgl  – beta-glucosidase",
            "K01183": "chitinase",
            "K01212": "chitobiase",
            "K01200": "xylanase",
            "K15921": "lcc  – laccase (lignin degradation)",
            # Fermentation / anaerobic carbon
            "K00016": "ldh  – L-lactate dehydrogenase",
            "K00626": "atoB – acetyl-CoA acetyltransferase",
            "K01512": "ppa  – inorganic pyrophosphatase",
        },
        "gene_names": [
            "rbcL", "rbcS", "cbbL", "cbbS", "rbcX",
            "coxL", "coxM", "coxS",
            "acsA", "acsB", "cooS",
            "fthfs", "metF",
            "mcrA", "mcrB", "mcrG", "mcrC", "mcrD",
            "fwdA", "fwdB", "fwdC", "fwdD",
            "fmdA", "fmdB", "fmdC",
            "mch", "mtd", "mer", "mtrA", "mtrH",
            "pmoA", "pmoB", "pmoC",
            "mmoX", "mmoY", "mmoZ", "mmoB", "mmoC", "mmoD",
            "cel", "endA", "endB", "bglA", "bglB",
            "chi", "chitA",
            "xynA", "xynB",
            "lcc", "lacA",
        ],
        "product_keywords": [
            "rubisco", "ribulose-1,5-bisphosphate",
            "carbon fixation", "carbon dioxide fixation",
            "methyl-coenzyme m reductase", "methanogenesis",
            "methane monooxygenase", "methane oxidation",
            "carbon monoxide dehydrogenase",
            "formate dehydrogenase",
            "formylmethanofuran",
            "acetyl-coa synthase", "wood-ljungdahl",
            "cellulase", "endoglucanase", "beta-glucosidase",
            "chitinase", "xylanase", "ligninase", "laccase",
        ],
    },

    # ── NITROGEN CYCLING ─────────────────────────────────────────────────────
    "Nitrogen_Cycling": {
        "description": (
            "Nitrogen cycling: fixation, nitrification, denitrification, "
            "DNRA, anammox, ammonia assimilation, urea degradation"
        ),
        "ko_ids": {
            # Nitrogen fixation
            "K02588": "nifH – nitrogenase iron protein",
            "K02586": "nifD – nitrogenase Mo-Fe protein α",
            "K02591": "nifK – nitrogenase Mo-Fe protein β",
            "K02585": "nifB – nitrogenase FeMo cofactor biosynthesis",
            "K02584": "nifA – nitrogenase transcriptional regulator",
            "K04488": "anfH – alternative nitrogenase iron protein",
            "K04489": "anfD – alternative nitrogenase D",
            "K04490": "anfK – alternative nitrogenase K",
            # Nitrification
            "K10944": "amoA – ammonia monooxygenase subunit A",
            "K10945": "amoB – ammonia monooxygenase subunit B",
            "K10946": "amoC – ammonia monooxygenase subunit C",
            "K10535": "hao  – hydroxylamine dehydrogenase",
            "K10534": "hcp  – hybrid cluster protein",
            # Nitrate/nitrite reduction (denitrification)
            "K00370": "narG – nitrate reductase α (denitrification)",
            "K00371": "narH – nitrate reductase β",
            "K00374": "narI – nitrate reductase γ",
            "K02567": "napA – periplasmic nitrate reductase large subunit",
            "K02568": "napB – periplasmic nitrate reductase small subunit",
            "K00368": "nirS – nitrite reductase, NO-forming (cd1-type)",
            "K00362": "nirK – nitrite reductase, NO-forming (Cu-type)",
            "K04561": "norB – nitric oxide reductase B",
            "K02305": "norC – nitric oxide reductase C",
            "K00376": "nosZ – nitrous oxide reductase",
            # DNRA
            "K03385": "nrfA – nitrite reductase (DNRA)",
            "K15876": "nrfH – nitrite reductase electron carrier",
            # Assimilatory nitrate reduction
            "K00367": "nasA – assimilatory nitrate reductase",
            "K00366": "nirA – ferredoxin-nitrite reductase",
            # Anammox
            "K20932": "hzo  – hydrazine dehydrogenase (anammox)",
            "K20933": "hdh  – hydrazine dehydrogenase (anammox)",
            "K20934": "hzs  – hydrazine synthase",
            # Ammonia assimilation
            "K01915": "glnA – glutamine synthetase",
            "K00265": "gltB – glutamate synthase GOGAT (large subunit)",
            "K00266": "gltD – glutamate synthase GOGAT (small subunit)",
            "K00260": "gdhA – glutamate dehydrogenase",
            "K00284": "gltS – glutamate synthase (ferredoxin-dependent)",
            # Urease
            "K01427": "ureA – urease γ subunit",
            "K01428": "ureB – urease β subunit",
            "K01429": "ureC – urease α subunit (catalytic)",
        },
        "gene_names": [
            "nifH", "nifD", "nifK", "nifE", "nifN", "nifA", "nifB",
            "anfH", "anfD", "anfK",
            "vnfH", "vnfD", "vnfK",
            "amoA", "amoB", "amoC",
            "hao", "hcp",
            "narG", "narH", "narI", "narJ",
            "napA", "napB", "napC", "napD",
            "nirS", "nirK", "nirB", "nirD",
            "norB", "norC", "norD",
            "nosZ", "nosL", "nosD", "nosF", "nosY",
            "nrfA", "nrfH",
            "hzo", "hdh", "hzs",
            "ureA", "ureB", "ureC", "ureD", "ureE", "ureF", "ureG",
            "glnA", "gltB", "gltD", "gdhA",
            "nasA", "nirA", "nrtA", "nrtB",
        ],
        "product_keywords": [
            "nitrogenase", "nitrogen fixation", "dinitrogen reductase",
            "ammonia monooxygenase", "hydroxylamine",
            "nitrate reductase", "nitrite reductase",
            "nitric oxide reductase", "nitrous oxide reductase",
            "urease", "urea hydrolysis",
            "denitrification",
            "hydrazine dehydrogenase", "hydrazine synthase", "anammox",
            "glutamine synthetase", "glutamate synthase",
            "dissimilatory nitrate", "assimilatory nitrate",
            "nitrification",
        ],
    },

    # ── PHOSPHORUS CYCLING ────────────────────────────────────────────────────
    "Phosphorus_Cycling": {
        "description": (
            "Phosphorus cycling: phosphate transport, organic-P mineralisation, "
            "polyphosphate metabolism, phosphonate utilisation, regulation"
        ),
        "ko_ids": {
            # Phosphatases / mineralisation
            "K01077": "phoA – alkaline phosphatase",
            "K01113": "phoD – alkaline phosphatase D",
            "K01093": "appA – phytase / non-specific phosphatase",
            "K03100": "phoN – non-specific acid phosphatase",
            "K01109": "phoC – phosphocholine phosphatase",
            # Phosphate ABC transporter (high-affinity)
            "K02040": "pstS – phosphate-binding protein (periplasmic)",
            "K02038": "pstA – phosphate ABC transporter permease 1",
            "K02037": "pstC – phosphate ABC transporter permease 2",
            "K02036": "pstB – phosphate ABC transporter ATPase",
            # Low-affinity phosphate transport
            "K03306": "pitA – low-affinity inorganic phosphate transporter",
            "K08019": "pitB – low-affinity inorganic phosphate transporter",
            # Phosphate two-component regulatory system
            "K07657": "phoB – phosphate regulon response regulator",
            "K07658": "phoR – phosphate sensor histidine kinase",
            "K09474": "phoU – phosphate uptake regulator",
            # Polyphosphate
            "K00937": "ppk  – polyphosphate kinase",
            "K22468": "ppk2 – polyphosphate kinase 2",
            "K01514": "ppx  – exopolyphosphatase",
            # Phosphonate transport & catabolism
            "K06163": "phnC – phosphonate ABC transporter ATPase",
            "K06164": "phnD – phosphonate-binding protein",
            "K06165": "phnE – phosphonate transporter permease",
            "K05781": "ugpA – glycerophosphodiester transporter permease",
            "K05782": "ugpB – glycerophosphodiester-binding protein",
            "K05783": "ugpC – glycerophosphodiester ABC transporter ATPase",
            "K05784": "ugpE – glycerophosphodiester permease",
            # Other
            "K02385": "phoE – outer membrane phosphoporin E",
        },
        "gene_names": [
            "phoA", "phoD", "phoN", "phoC", "phoB", "phoR", "phoH", "phoP",
            "pstA", "pstB", "pstC", "pstS", "pstP",
            "phoU",
            "ppk", "ppk2", "ppx",
            "appA",
            "phnC", "phnD", "phnE", "phnF", "phnG", "phnH", "phnI", "phnJ",
            "phnK", "phnL", "phnM", "phnN", "phnO", "phnP", "phnQ",
            "ugpA", "ugpB", "ugpC", "ugpE",
            "pitA", "pitB",
            "gcd",
        ],
        "product_keywords": [
            "alkaline phosphatase", "acid phosphatase",
            "phosphate transporter", "phosphate binding protein",
            "phosphate abc transporter",
            "polyphosphate kinase", "exopolyphosphatase", "polyphosphate",
            "phytase", "phosphonate",
            "phosphate solubilization",
            "phosphate regulatory", "phosphate sensor",
            "glycerophosphodiester",
            "outer membrane phosphoporin",
            "phosphate uptake",
        ],
    },

    # ── ANTIBIOTIC RESISTANCE ─────────────────────────────────────────────────
    "Antibiotic_Resistance": {
        "description": (
            "Antibiotic resistance: beta-lactamases, efflux pumps, "
            "target modification, enzymatic inactivation, mobile integrons"
        ),
        "ko_ids": {
            # Beta-lactamases / cell-wall
            "K01467": "blaZ – beta-lactamase (penicillin hydrolysis)",
            "K18369": "mecA – PBP2a (methicillin/MRSA resistance)",
            # RND efflux pumps
            "K03584": "acrA – multidrug efflux membrane fusion protein",
            "K03585": "acrB – multidrug efflux pump (RND family)",
            "K06147": "tolC – outer membrane efflux channel",
            "K07297": "mexA – MexAB-OprM efflux MFP",
            "K07298": "mexB – MexAB-OprM RND pump",
            "K07299": "mexC – MexCD-OprJ RND pump accessory",
            "K07300": "mexD – MexCD-OprJ RND pump",
            "K07301": "oprM – MexAB-OprM outer membrane protein",
            "K18138": "mexY – MexXY-OprM RND pump",
            "K18139": "mexX – MexXY-OprM MFP",
            # Vancomycin resistance
            "K21750": "vanA – D-Ala:D-Lac ligase (vancomycin resistance)",
            "K21751": "vanH – D-lactate dehydrogenase (vancomycin)",
            "K21752": "vanX – D-Ala-D-Lac dipeptidase",
            "K21753": "vanY – D-Ala-D-Ala carboxypeptidase",
            "K21754": "vanR – vancomycin response regulator",
            "K21755": "vanS – vancomycin sensor kinase",
            # Tetracycline
            "K07723": "tetA – tetracycline efflux pump A",
            "K07725": "tetB – tetracycline efflux pump B",
            "K18777": "tetM – tetracycline ribosomal protection protein",
            # Chloramphenicol
            "K00558": "cat  – chloramphenicol acetyltransferase",
            # Colistin / polymyxin
            "K09687": "arnA – polymyxin resistance (ArnA bifunctional enzyme)",
            # Aminoglycoside
            "K00003": "aacA – aminoglycoside 6'-N-acetyltransferase",
        },
        "gene_names": [
            # Beta-lactamases (broad)
            "blaZ", "blaB", "blaC", "blaTEM", "blaSHV", "blaCTX", "blaOXA",
            "blaKPC", "blaNDM", "blaVIM", "blaIMP", "blaL",
            # Methicillin
            "mecA", "mecC", "mecI", "mecR",
            # Vancomycin
            "vanA", "vanB", "vanC", "vanD", "vanE",
            "vanH", "vanR", "vanS", "vanX", "vanY", "vanZ",
            # Tetracycline
            "tetA", "tetB", "tetC", "tetD", "tetE", "tetG", "tetH",
            "tetL", "tetM", "tetO", "tetQ", "tetS", "tetW", "tetX",
            # Efflux – AcrAB-TolC / RND family
            "acrA", "acrB", "acrD", "acrF", "acrR",
            "mexA", "mexB", "mexC", "mexD", "mexE", "mexF", "mexX", "mexY",
            "mexZ",
            "tolC", "oprM", "oprN",
            # Aminoglycoside
            "aacA", "aacC", "aacC2", "aac6",
            "aphA", "aphE",
            "aadA", "aadB",
            "strA", "strB",
            # Macrolide / lincosamide / streptogramin (MLS)
            "ermA", "ermB", "ermC", "ermD", "ermF",
            # Colistin
            "mcr1", "mcr2", "mcr3", "mcr4", "mcr5",
            # Quinolone / fluoroquinolone
            "qnrA", "qnrB", "qnrC", "qnrD", "qnrS",
            # Sulfonamide
            "sul1", "sul2", "sul3",
            # Trimethoprim
            "dfrA", "dfrB",
            # Chloramphenicol
            "catA", "catB", "catP", "cmlA", "cfr",
            # Rifamycin
            "arr", "rph",
            # Integrons / mobile elements
            "intI1", "intI2", "intI3",
            "qacE", "qacA", "qacB",
            # Fosfomycin
            "fosA", "fosB", "fosX",
        ],
        "product_keywords": [
            "beta-lactamase", "penicillinase", "cephalosporinase",
            "efflux pump", "resistance-nodulation-division",
            "multidrug resistance", "multidrug efflux",
            "vancomycin resistance",
            "tetracycline resistance", "tetracycline efflux",
            "aminoglycoside acetyltransferase", "aminoglycoside resistance",
            "macrolide resistance", "erythromycin resistance",
            "chloramphenicol acetyltransferase",
            "fluoroquinolone resistance", "quinolone resistance",
            "sulfonamide resistance",
            "methicillin resistance", "penicillin-binding protein 2",
            "colistin resistance",
            "antibiotic resistance", "antimicrobial resistance",
            "rifampin resistance", "rifampicin resistance",
            "integrase, mobile",
        ],
    },

    # ── BIOFILM FORMATION ─────────────────────────────────────────────────────
    "Biofilm_Formation": {
        "description": (
            "Biofilm formation: EPS/matrix biosynthesis, quorum sensing, "
            "type-IV pili, flagella, c-di-GMP signalling, adhesion"
        ),
        "ko_ids": {
            # Curli fibres
            "K11931": "csgA – curli minor subunit",
            "K11937": "csgB – curli nucleator",
            "K11932": "csgC – curli assembly outer-membrane lipoprotein",
            "K11933": "csgD – biofilm/curli master regulator",
            # Cellulose synthase
            "K00694": "bcsA – cellulose synthase catalytic subunit",
            "K11735": "bcsB – cellulose synthase regulatory",
            # PGA (poly-GlcNAc) matrix
            "K03503": "pgaA – PGA outer-membrane porin",
            "K03504": "pgaB – PGA deacetylase",
            "K03505": "pgaC – PGA synthase",
            "K03506": "pgaD – PGA synthase regulatory",
            # Alginate biosynthesis (Pseudomonas)
            "K03970": "algD – GDP-mannose 6-dehydrogenase",
            "K06041": "algR – alginate regulatory protein",
            # Pel polysaccharide
            "K19712": "pelA – Pel polysaccharide biosynthesis/deacetylase",
            "K19713": "pelB – Pel polysaccharide export (outer membrane)",
            "K19714": "pelC – Pel polysaccharide periplasmic export",
            # Type IV pili
            "K02648": "pilA – type IV pilin (major)",
            "K02649": "pilB – type IV pili assembly ATPase",
            "K02650": "pilC – type IV pili assembly membrane protein",
            "K02651": "pilD – type IV prepilin peptidase",
            "K02652": "pilE – type IV pilin minor",
            "K02653": "pilF – type IV pili assembly ATPase",
            "K02654": "pilG – twitching motility protein PilG",
            "K02656": "pilT – type IV pili retraction ATPase",
            # Flagella
            "K02413": "fliA – flagellar biosynthesis sigma factor σ28",
            "K02405": "fliC – flagellin (structural)",
            "K02409": "fliF – flagellar basal-body MS-ring",
            "K02408": "fliE – flagellar hook–basal body complex",
            "K02397": "flhA – flagellar biosynthesis protein FlhA",
            "K02401": "flhB – flagellar biosynthesis protein FlhB",
            "K02403": "flhC – flagellar regulatory protein FlhC",
            "K02404": "flhD – flagellar regulatory protein FlhD",
            # Quorum sensing
            "K10909": "luxR – LuxR-family transcriptional regulator",
            "K03122": "luxI – autoinducer (AHL) synthase",
            "K16012": "luxS – S-ribosylhomocysteine lyase (AI-2 synthesis)",
            # Wsp chemosensory (c-di-GMP biofilm)
            "K12356": "wspA – methyl-accepting chemotaxis receptor Wsp",
            "K00987": "wspE – Wsp sensor histidine kinase / response reg.",
            # Sigma factor
            "K03088": "rpoS – stationary-phase sigma factor σ38",
        },
        "gene_names": [
            # Curli
            "csgA", "csgB", "csgC", "csgD", "csgE", "csgF", "csgG",
            # Alginate
            "algA", "algB", "algC", "algD", "algE", "algF", "algG", "algI",
            "algJ", "algK", "algL", "algR", "algQ", "algP", "algU", "algW",
            # Pel polysaccharide
            "pelA", "pelB", "pelC", "pelD", "pelE", "pelF", "pelG",
            # Psl polysaccharide
            "pslA", "pslB", "pslC", "pslD", "pslE", "pslF", "pslG",
            "pslH", "pslI", "pslJ", "pslK", "pslL",
            # PGA matrix
            "pgaA", "pgaB", "pgaC", "pgaD",
            # Cellulose synthase
            "bcsA", "bcsB", "bcsC", "bcsD", "bcsZ",
            # Type IV pili
            "pilA", "pilB", "pilC", "pilD", "pilE", "pilF", "pilG", "pilH",
            "pilI", "pilJ", "pilK", "pilL", "pilM", "pilN", "pilO", "pilP",
            "pilQ", "pilR", "pilS", "pilT", "pilU", "pilV", "pilW", "pilX",
            "pilY1", "pilZ",
            # Type I fimbriae
            "fimA", "fimB", "fimC", "fimD", "fimE", "fimF", "fimG", "fimH",
            # Type III fimbriae (MrkABCDF)
            "mrkA", "mrkB", "mrkC", "mrkD", "mrkF",
            # Flagella
            "fliA", "fliC", "fliD", "fliE", "fliF", "fliG", "fliH",
            "flhA", "flhB", "flhC", "flhD",
            "motA", "motB",
            # Quorum sensing
            "luxI", "luxR", "luxS",
            "lasI", "lasR",
            "rhlI", "rhlR",
            "agrA", "agrB", "agrC", "agrD",
            "qseB", "qseC",
            "comA", "comC", "comD", "comE",
            # c-di-GMP metabolism
            "dgcA", "dgcB", "dgcC", "dgcE", "dgcF", "dgcH", "dgcJ",
            "pdeA", "pdeB", "pdeC",
            "wspA", "wspB", "wspC", "wspD", "wspE", "wspF",
            # VPS (Vibrio)
            "vpsA", "vpsB", "vpsC", "vpsL",
            # Biofilm regulators / sigma
            "rpoS", "bssR", "bssS",
            "sinR", "sinI",
            # Colanic acid
            "wcaA", "wcaB", "wcaC", "wcaD", "wcaE", "wcaF",
        ],
        "product_keywords": [
            "biofilm", "quorum sensing",
            "flagellin", "flagella", "flagellum",
            "type iv pili", "type iv pilin", "fimbriae", "fimbrial",
            "alginate synthase", "alginate biosynthesis", "alginate regulatory",
            "curli", "curlin",
            "exopolysaccharide", "polysaccharide biosynthesis",
            "autoinducer", "acyl-homoserine lactone",
            "c-di-gmp", "diguanylate cyclase",
            "biofilm matrix", "surface attachment",
            "cellulose synthase",
            "poly-beta-1,6-n-acetyl-d-glucosamine",
        ],
    },
}

# MetaCyc pathway name keywords used to filter PICRUSt2 pathway tables
PATHWAY_KEYWORDS = {
    "Carbon_Cycling": [
        "CALVIN", "carbon fixation", "reductive acetyl",
        "METHANOGEN", "methanogenesis", "methane",
        "WOOD-LJUNGDAHL", "wood ljungdahl",
        "carbon monoxide", "formyl",
        "cellulose degradation", "chitin degradation",
        "GLUCARDEG", "CENTFERM", "TCA", "citric acid",
    ],
    "Nitrogen_Cycling": [
        "nitrogen fixation", "NITRIF", "nitrification",
        "DENITRIF", "denitrification",
        "nitrate reduction", "nitrite reduction",
        "ammonia assimilation", "AMMASSIM",
        "urea degradation", "UREX",
        "anammox", "hydrazine",
        "glutamine biosynthesis", "glutamate biosynthesis",
    ],
    "Phosphorus_Cycling": [
        "phosphate acquisition", "phosphate transport",
        "polyphosphate", "phosphonate",
        "organic phosphorus", "phosphate solubilization",
        "phytate", "phytase",
    ],
    "Antibiotic_Resistance": [
        "beta-lactam", "antibiotic resistance",
        "vancomycin", "tetracycline resistance",
        "aminoglycoside", "macrolide",
        "multidrug resistance", "efflux",
        "antimicrobial",
    ],
    "Biofilm_Formation": [
        "biofilm", "quorum sensing",
        "flagella", "pili", "fimbriae",
        "alginate", "exopolysaccharide",
        "curli", "cellulose",
        "autoinducer",
    ],
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def extract_sample_name(filepath: str) -> str:
    """Return sample identifier from a MetaCerberus TSV filename."""
    base = os.path.basename(filepath)
    # Remove common prefix and .tsv suffix
    name = re.sub(r"^final_annotation_summary_", "", base, flags=re.IGNORECASE)
    name = re.sub(r"\.tsv$", "", name, flags=re.IGNORECASE)
    return name


def load_annotation_files(directory: str) -> dict[str, pd.DataFrame]:
    """Load all MetaCerberus annotation TSV files from a directory."""
    pattern = os.path.join(directory, "final_annotation_summary_*.tsv")
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"  WARNING: No annotation TSV files found in {directory}")
        return {}

    samples = {}
    for filepath in files:
        sample = extract_sample_name(filepath)
        try:
            df = pd.read_csv(filepath, sep="\t", dtype=str, low_memory=False)
            df = df.fillna("")
            samples[sample] = df
            print(f"  Loaded {sample}: {len(df)} rows")
        except Exception as exc:
            print(f"  ERROR loading {filepath}: {exc}")
    return samples


def _split_genes(gene_cell: str) -> list[str]:
    """Split the comma/space-delimited gene field into individual names."""
    return [g.strip() for g in re.split(r"[,\s]+", gene_cell) if g.strip()]


def search_metacerberus(
    sample_name: str,
    df: pd.DataFrame,
    categories: dict,
) -> pd.DataFrame:
    """
    Search a MetaCerberus annotation dataframe for genes of interest.
    Returns a dataframe of all matching rows with extra columns for
    the matched category and the method of matching.
    """
    results = []
    col_ko      = "best_hit"
    col_gene    = "gene"
    col_product = "product"

    for cat_name, cat_info in categories.items():
        ko_set      = set(cat_info["ko_ids"].keys())
        gene_set    = {g.lower() for g in cat_info["gene_names"]}
        kw_patterns = [re.compile(re.escape(kw), re.IGNORECASE)
                       for kw in cat_info["product_keywords"]]

        for _, row in df.iterrows():
            ko_val      = str(row.get(col_ko,      "")).strip()
            gene_val    = str(row.get(col_gene,    "")).strip()
            product_val = str(row.get(col_product, "")).strip()

            matched    = False
            match_how  = []

            # 1. KO ID match
            if ko_val and ko_val in ko_set:
                matched   = True
                ko_desc   = cat_info["ko_ids"][ko_val]
                match_how.append(f"KO:{ko_val} ({ko_desc})")

            # 2. Gene name match
            gene_parts = _split_genes(gene_val)
            for gpart in gene_parts:
                if gpart.lower() in gene_set:
                    matched = True
                    match_how.append(f"gene:{gpart}")
                    break

            # 3. Product keyword match
            if product_val:
                for pat in kw_patterns:
                    if pat.search(product_val):
                        matched = True
                        match_how.append(f"keyword:{pat.pattern}")
                        break

            if matched:
                record = dict(row)
                record["Sample"]    = sample_name
                record["Category"]  = cat_name
                record["Match_Method"] = "; ".join(dict.fromkeys(match_how))
                results.append(record)

    if not results:
        return pd.DataFrame()

    out = pd.DataFrame(results)
    # De-duplicate: same ORF matched by multiple methods
    dedup_cols = ["Sample", "target", "Category"]
    if all(c in out.columns for c in dedup_cols):
        out = out.drop_duplicates(subset=dedup_cols)
    return out


# =============================================================================
# PICRUST2 / BIOM PARSING
# =============================================================================

def _read_biom_hdf5(biom_bytes: bytes) -> pd.DataFrame:
    """
    Parse a BIOM v2.1 (HDF5) file from an in-memory bytes object.
    Returns a DataFrame with features as rows and samples as columns.
    Requires h5py.
    """
    import io
    import numpy as np
    import scipy.sparse as sp

    with tempfile.NamedTemporaryFile(suffix=".biom", delete=False) as tmp:
        tmp.write(biom_bytes)
        tmp_path = tmp.name

    try:
        import h5py
        with h5py.File(tmp_path, "r") as f:
            obs_ids  = [s.decode() if isinstance(s, bytes) else s
                        for s in f["observation/ids"][:]]
            samp_ids = [s.decode() if isinstance(s, bytes) else s
                        for s in f["sample/ids"][:]]

            data    = f["observation/matrix/data"][:]
            indices = f["observation/matrix/indices"][:]
            indptr  = f["observation/matrix/indptr"][:]

        matrix = sp.csr_matrix(
            (data, indices, indptr),
            shape=(len(obs_ids), len(samp_ids))
        ).toarray()

        df = pd.DataFrame(matrix, index=obs_ids, columns=samp_ids)
        df.index.name = "pathway"
        return df
    finally:
        os.unlink(tmp_path)


def load_picrust2_pathways(qza_path: str) -> pd.DataFrame | None:
    """
    Extract and return the pathway abundance table from a PICRUSt2 .qza file.
    Returns None if the file cannot be read.
    """
    if not os.path.isfile(qza_path):
        print(f"  WARNING: {qza_path} not found – skipping PICRUSt2 analysis.")
        return None

    try:
        with zipfile.ZipFile(qza_path) as zf:
            biom_entries = [n for n in zf.namelist() if n.endswith(".biom")]
            if not biom_entries:
                print("  WARNING: No .biom file found inside the QZA archive.")
                return None
            biom_bytes = zf.read(biom_entries[0])
    except Exception as exc:
        print(f"  ERROR reading QZA: {exc}")
        return None

    # HDF5 BIOM (magic bytes \x89HDF)
    if biom_bytes[:4] == b"\x89HDF":
        if not HAS_H5PY:
            print(
                "  WARNING: The pathway_abundance.qza contains an HDF5 BIOM file "
                "but h5py is not installed.\n"
                "  Install it with:  pip install h5py\n"
                "  PICRUSt2 pathway results will be skipped."
            )
            return None
        try:
            import scipy.sparse  # noqa: F401 – needed inside _read_biom_hdf5
        except ImportError:
            print(
                "  WARNING: scipy is required to read the BIOM file.\n"
                "  Install it with:  pip install scipy\n"
                "  PICRUSt2 pathway results will be skipped."
            )
            return None
        try:
            df = _read_biom_hdf5(biom_bytes)
            print(f"  Loaded BIOM: {df.shape[0]} pathways × {df.shape[1]} samples")
            return df
        except Exception as exc:
            print(f"  ERROR parsing BIOM: {exc}")
            return None

    # JSON BIOM (fallback – much less common with PICRUSt2)
    if biom_bytes[:1] == b"{":
        import json
        try:
            obj    = json.loads(biom_bytes)
            rows   = [r["id"] for r in obj["rows"]]
            cols   = [c["id"] for c in obj["columns"]]
            data   = {c: [0.0] * len(rows) for c in cols}
            for i, j, val in obj["data"]:
                data[cols[j]][i] = val
            df = pd.DataFrame(data, index=rows)
            df.index.name = "pathway"
            print(f"  Loaded JSON BIOM: {df.shape[0]} pathways × {df.shape[1]} samples")
            return df
        except Exception as exc:
            print(f"  ERROR parsing JSON BIOM: {exc}")
            return None

    print("  WARNING: Unrecognised BIOM format – skipping PICRUSt2 analysis.")
    return None


def search_pathways(
    pathway_df: pd.DataFrame,
    pathway_keywords: dict,
) -> pd.DataFrame:
    """
    Filter the pathway abundance table for pathways matching the keyword lists
    for each category.  Returns a long-format dataframe.
    """
    rows = []
    for cat_name, keywords in pathway_keywords.items():
        patterns = [re.compile(re.escape(kw), re.IGNORECASE) for kw in keywords]
        mask = pathway_df.index.to_series().apply(
            lambda p: any(pat.search(p) for pat in patterns)
        )
        subset = pathway_df.loc[mask]
        for pathway_id, abundances in subset.iterrows():
            for sample, abund in abundances.items():
                rows.append({
                    "Category":    cat_name,
                    "Pathway_ID":  pathway_id,
                    "Sample":      sample,
                    "Abundance":   abund,
                })
    if not rows:
        return pd.DataFrame(columns=["Category", "Pathway_ID", "Sample", "Abundance"])
    return pd.DataFrame(rows)


# =============================================================================
# SUMMARY TABLES
# =============================================================================

def build_annotation_summary(all_hits: pd.DataFrame) -> pd.DataFrame:
    """
    Pivot table: rows = Category, columns = Sample, values = gene-hit counts.
    """
    if all_hits.empty:
        return pd.DataFrame()
    grp = (
        all_hits.groupby(["Category", "Sample"])
        .size()
        .reset_index(name="Gene_Hits")
    )
    pivot = grp.pivot(index="Category", columns="Sample", values="Gene_Hits").fillna(0)
    pivot = pivot.astype(int)
    pivot.columns.name = None
    return pivot


def build_pathway_summary(pathway_long: pd.DataFrame) -> pd.DataFrame:
    """
    Pivot table: rows = Category × Pathway_ID, columns = Sample,
    values = mean pathway abundance (already collapsed per pathway).
    """
    if pathway_long.empty:
        return pd.DataFrame()
    pivot = pathway_long.pivot_table(
        index=["Category", "Pathway_ID"],
        columns="Sample",
        values="Abundance",
        aggfunc="sum",
    ).fillna(0)
    pivot.columns.name = None
    return pivot


# =============================================================================
# EXCEL OUTPUT
# =============================================================================

def _auto_width(ws):
    """Widen Excel columns to roughly fit their content."""
    for col_cells in ws.columns:
        max_len = max((len(str(c.value or "")) for c in col_cells), default=10)
        ws.column_dimensions[col_cells[0].column_letter].width = min(max_len + 4, 60)


def write_excel(
    output_path: str,
    annotation_summary: pd.DataFrame,
    all_hits: pd.DataFrame,
    pathway_summary: pd.DataFrame,
    pathway_long: pd.DataFrame,
) -> None:
    """Write all results to a multi-sheet Excel workbook."""

    def df_to_sheet(df: pd.DataFrame, ws, title_row: str = ""):
        """Write a DataFrame to an openpyxl worksheet."""
        from openpyxl.styles import PatternFill, Font, Alignment, Border, Side
        from openpyxl.utils import get_column_letter

        HEADER_FILL   = PatternFill("solid", fgColor="2E75B6")
        TITLE_FILL    = PatternFill("solid", fgColor="1F4E79")
        ALT_FILL      = PatternFill("solid", fgColor="D6E4F0")
        HEADER_FONT   = Font(color="FFFFFF", bold=True, size=10)
        TITLE_FONT    = Font(color="FFFFFF", bold=True, size=11)
        BODY_FONT     = Font(size=9)
        THIN          = Side(style="thin", color="BFBFBF")
        BORDER        = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)

        start_row = 1
        if title_row:
            cell = ws.cell(row=1, column=1, value=title_row)
            cell.font   = TITLE_FONT
            cell.fill   = TITLE_FILL
            cell.alignment = Alignment(horizontal="left", vertical="center")
            n_cols = max(len(df.columns) + (1 if df.index.name else 0), 1)
            ws.merge_cells(
                start_row=1, start_column=1,
                end_row=1,   end_column=n_cols,
            )
            ws.row_dimensions[1].height = 20
            start_row = 2

        # Reset index so it appears as a column
        df_out = df.reset_index()
        headers = list(df_out.columns)

        for col_idx, header in enumerate(headers, start=1):
            cell = ws.cell(row=start_row, column=col_idx, value=header)
            cell.font        = HEADER_FONT
            cell.fill        = HEADER_FILL
            cell.alignment   = Alignment(horizontal="center", vertical="center",
                                          wrap_text=True)
            cell.border      = BORDER

        ws.row_dimensions[start_row].height = 30

        for r_idx, data_row in enumerate(df_out.itertuples(index=False), start=start_row + 1):
            fill = ALT_FILL if (r_idx - start_row) % 2 == 0 else None
            for c_idx, value in enumerate(data_row, start=1):
                cell = ws.cell(row=r_idx, column=c_idx,
                               value=(None if pd.isna(value) else value))
                cell.font   = BODY_FONT
                cell.border = BORDER
                if fill:
                    cell.fill = fill
                cell.alignment = Alignment(vertical="top", wrap_text=False)

        # Column widths
        for col_idx, col_name in enumerate(headers, start=1):
            col_vals = [str(col_name)] + [
                str(v) for v in df_out.iloc[:, col_idx - 1].tolist()
            ]
            width = min(max(len(v) for v in col_vals) + 4, 60)
            ws.column_dimensions[get_column_letter(col_idx)].width = width

    # ── Build workbook ────────────────────────────────────────────────────────
    import openpyxl
    wb = openpyxl.Workbook()

    # ── Sheet 1: Overview ────────────────────────────────────────────────────
    ws_ovr = wb.active
    ws_ovr.title = "Overview"
    if not annotation_summary.empty:
        df_to_sheet(
            annotation_summary,
            ws_ovr,
            title_row="MetaCerberus Gene Hits per Category per Sample",
        )

    # ── Sheet 2: All hits (detailed) ─────────────────────────────────────────
    ws_all = wb.create_sheet("All_Gene_Hits")
    if not all_hits.empty:
        df_to_sheet(
            all_hits.reset_index(drop=True),
            ws_all,
            title_row="All Matched Genes (MetaCerberus Annotations)",
        )

    # ── Sheets 3-7: Per-category detail ──────────────────────────────────────
    CATEGORY_ORDER = list(GENE_CATEGORIES.keys())
    for cat in CATEGORY_ORDER:
        safe_name = cat[:31]  # Excel sheet-name limit
        ws_cat = wb.create_sheet(safe_name)
        subset = all_hits[all_hits["Category"] == cat] if not all_hits.empty else pd.DataFrame()
        desc   = GENE_CATEGORIES[cat]["description"]
        title  = f"{cat.replace('_', ' ')} — {desc}"
        if not subset.empty:
            df_to_sheet(subset.reset_index(drop=True), ws_cat, title_row=title)
        else:
            ws_cat["A1"] = title
            ws_cat["A2"] = "No matching genes found in any sample."

    # ── Sheet 8: PICRUSt2 pathway summary ────────────────────────────────────
    ws_pw = wb.create_sheet("PICRUSt2_Pathways")
    if not pathway_summary.empty:
        df_to_sheet(
            pathway_summary.reset_index(),
            ws_pw,
            title_row="PICRUSt2 Pathway Abundances (MetaCyc) – Matched by Keyword",
        )
    else:
        ws_pw["A1"] = "PICRUSt2 Pathway Results"
        if not HAS_H5PY:
            ws_pw["A2"] = (
                "h5py is not installed.  "
                "Run:  pip install h5py scipy  then re-run this script."
            )
        else:
            ws_pw["A2"] = "No matching pathways found or data could not be loaded."

    # ── Sheet 9: PICRUSt2 long format ─────────────────────────────────────────
    ws_pl = wb.create_sheet("PICRUSt2_Long")
    if not pathway_long.empty:
        df_to_sheet(
            pathway_long.reset_index(drop=True),
            ws_pl,
            title_row="PICRUSt2 Pathway Abundances – Long Format",
        )
    else:
        ws_pl["A1"] = "No PICRUSt2 pathway data available."

    wb.save(output_path)
    print(f"\n  Results saved → {output_path}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    timestamp  = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_xls = os.path.join(script_dir, f"gene_search_results_{timestamp}.xlsx")
    qza_path   = os.path.join(script_dir, "pathway_abundance.qza")

    print("=" * 70)
    print("  Functional Gene Search")
    print(f"  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    # ── 1. MetaCerberus annotations ──────────────────────────────────────────
    print("\n[1/4] Loading MetaCerberus annotation files …")
    samples = load_annotation_files(script_dir)

    print("\n[2/4] Searching annotations …")
    all_hits_list = []
    for sample_name, df in samples.items():
        print(f"  Processing: {sample_name}")
        hits = search_metacerberus(sample_name, df, GENE_CATEGORIES)
        all_hits_list.append(hits)

    all_hits = (
        pd.concat(all_hits_list, ignore_index=True)
        if all_hits_list else pd.DataFrame()
    )

    if not all_hits.empty:
        print(f"\n  Total gene hits across all samples: {len(all_hits)}")
        for cat in GENE_CATEGORIES:
            n = len(all_hits[all_hits["Category"] == cat])
            print(f"    {cat:30s}: {n} hits")
    else:
        print("  No gene hits found.")

    # ── 2. PICRUSt2 pathways ─────────────────────────────────────────────────
    print("\n[3/4] Loading PICRUSt2 pathway abundance …")
    pathway_df   = load_picrust2_pathways(qza_path)
    pathway_long = pd.DataFrame()

    if pathway_df is not None:
        print("  Searching pathways by keyword …")
        pathway_long = search_pathways(pathway_df, PATHWAY_KEYWORDS)
        if not pathway_long.empty:
            print(f"  Matched {pathway_long['Pathway_ID'].nunique()} pathways "
                  f"across {pathway_long['Sample'].nunique()} samples")
        else:
            print("  No matching pathways found.")

    # ── 3. Build summaries ───────────────────────────────────────────────────
    annotation_summary = build_annotation_summary(all_hits)
    pathway_summary    = build_pathway_summary(pathway_long)

    # ── 4. Write Excel ───────────────────────────────────────────────────────
    print("\n[4/4] Writing output …")
    write_excel(
        output_xls,
        annotation_summary,
        all_hits,
        pathway_summary,
        pathway_long,
    )

    print("\nDone.")
    print("=" * 70)
    print(f"  Output → {output_xls}")
    print("=" * 70)

    # Brief console summary
    print("\nSUMMARY")
    print("-" * 70)
    if not annotation_summary.empty:
        print("\nGene hits per category per sample (MetaCerberus):")
        print(annotation_summary.to_string())
    if not pathway_summary.empty:
        print("\nPathway categories matched (PICRUSt2):")
        for (cat, pid), row in pathway_summary.iterrows():
            total = row.sum()
            print(f"  [{cat}] {pid}  total abundance = {total:.2f}")


if __name__ == "__main__":
    main()

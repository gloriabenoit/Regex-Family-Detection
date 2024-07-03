"""Teste l'appartenance d'une séquence FASTA à un pattern PROSITE."""
import re
import urllib.request

# Motif PROSITE de la famille des opsines
PROSITE_OPSINE = "[LIVMFWAC]-[PSGAC]-x-{G}-x-[SAC]-K-[STALIMR]-[GSACPNV]"\
              "-[STACP]-x(2)-[DENF]-[AP]-x(2)-[IY]."

# Sequence a tester

# Opsin 1, human
SEQ_OPSSR_HUM = "MAQQWSLQRLAGRHPQDSYEDSTQSSIFTYTNSNSTRGPFEGPNYHIAPRWVYHLTSVWM"\
                "IFVVTASVFTNGLVLAATMKFKKLRHPLNWILVNLAVADLAETVIASTISIVNQVSGYFV"\
                "LGHPMCVLEGYTVSLCGITGLWSLAIISWERWMVVCKPFGNVRFDAKLAIVGIAFSWIWA"\
                "AVWTAPPIFGWSRYWPHGLKTSCGPDVFSGSSYPGVQSYMIVLMVTCCIIPLAIIMLCYL"\
                "QVWLAIRAVAKQQKESESTQKAEKEVTRMVVVMIFAYCVCWGPYTFFACFAAANPGYAFH"\
                "PLMAALPAYFAKSATIYNPVIYVFMNRQFRNCILQLFGKKVDDGSELSSASKTEVSSVSS"\
                "VSPA"
# Glucose transporter type 4, human
SEQ_GLUT4_HUM = "MPSGFQQIGSEDGEPPQQRVTGTLVLAVFSAVLGSLQFGYNIGVINAPQKVIEQSYNETW"\
                "LGRQGPEGPSSIPPGTLTTLWALSVAIFSVGGMISSFLIGIISQWLGRKRAMLVNNVLAV"\
                "LGGSLMGLANAAASYEMLILGRFLIGAYSGLTSGLVPMYVGEIAPTHLRGALGTLNQLAI"\
                "VIGILIAQVLGLESLLGTASLWPLLLGLTVLPALLQLVLLPFCPESPRYLYIIQNLEGPA"\
                "RKSLKRLTGWADVSGVLAELKDEKRKLERERPLSLLQLLGSRTHRQPLIIAVVLQLSQQL"\
                "SGINAVFYYSTSIFETAGVGQPAYATIGAGVVNTVFTLVSVLLVERAGRRTLHLLGLAGM"\
                "CGCAILMTVALLLLERVPAMSYVSIVAIFGFVAFFEIGPGPIPWFIVAELFSQGPRPAAM"\
                "AVAGFSNWTSNFIIGMGFQYVAEAMGPYVFLLFAVLLLGFFIFTFLRVPETRGRTFDQIS"\
                "AAFHRTPSLLEQEVKPSTELEYLGPDEND"
# Rhodopsin, human
SEQ_OPSD_HUM = "MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLY"\
               "VTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLG"\
               "GEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIP"\
               "EGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQES"\
               "ATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAI"\
               "YNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEASATVSKTETSQVAPA"
# Rhodopsin, dog
OPSD_CAN = "MNGTEGPNFYVPFSNKTGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLY"\
           "VTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNVEGFFATLG"\
           "GEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIP"\
           "EGMQCSCGIDYYTLKPEINNESFVIYMFVVHFAIPMIVIFFCYGQLVFTVKEAAAQQQES"\
           "ATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSDFGPIFMTLPAFFAKSSSI"\
           "YNPVIYIMMNKQFRNCMITTLCCGKNPLGDDEASASASKTETSQVAPA"
# Rhodopsin, cow
SEQ_OPSD_BOV = "MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLY"\
           "VTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLG"\
           "GEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIP"\
           "EGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQES"\
           "ATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAV"\
           "YNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA"
# Opsin 4, fruit fly
SEQ_OPS4_DRO = "MEPLCNASEPPLRPEARSSGNGDLQFLGWNVPPDQIQYIPEHWLTQLEPPASMHYMLGVF"\
           "YIFLFCASTVGNGMVIWIFSTSKSLRTPSNMFVLNLAVFDLIMCLKAPIFIYNSFHRGFA"\
           "LGNTWCQIFASIGSYSGIGAGMTNAAIGYDRYNVITKPMNRNMTFTKAVIMNIIIWLYCT"\
           "PWVVLPLTQFWDRFVPEGYLTSCSFDYLSDNFDTRLFVGTIFFFSFVCPTLMILYYYSQI"\
           "VGHVFSHEKALREQAKKMNVESLRSNVDKSKETAEIRIAKAAITICFLFFVSWTPYGVMS"\
           "LIGAFGDKSLLTPGATMIPACTCKLVACIDPFVYAISHPRYRLELQKRCPWLGVNEKSGE"\
           "ISSAQSTTTQEQQQTTAA"
# Hemoglobin subunit beta-1, rat
SEQ_HBB1_RAT = "MVHLTDAEKAAVNGLWGKVNPDDVGGEALGRLLVVYPWTQRYFDSFGDLSSASAIMGNPK"\
           "VKAHGKKVINAFNDGLKHLDNLKGTFAHLSELHCDKLHVDPENFRLLGNMIVIVLGHHLG"\
           "KEFTPCAQAAFQKVVAGVASALAHKYH"

LISTE_SEQ = [SEQ_OPSSR_HUM, SEQ_GLUT4_HUM, SEQ_OPSD_HUM, OPSD_CAN,
             SEQ_OPSD_BOV, SEQ_OPS4_DRO, SEQ_HBB1_RAT]


def traduit__pro_re(seq_prosite):
    """
    Traduit une séquence format PROSITE en son équivalent REGEX.

    Parameters
    ----------
    seq_prosite : str
        La séquence à traduire.

    Returns
    -------
    str
        La séquence traduite en REGEX.
    """
    # Création du dictionnaire de traduction
    liste_cle = ["{", "}", "x", r"\(", r"\)", "<", ">", "-", r"\."]
    liste_val = ["[^", "]", r"\\w", "{", "}", "^", "$", "", ""]
    dict_trad = dict(zip(liste_cle, liste_val))

    for cle, val in dict_trad.items():
        seq_prosite = re.sub(cle, val, seq_prosite)
    return seq_prosite


def dans_famille(seq_prot, seq_re, nom_famille):
    """
    Test d'appartenance à un pattern.

    Vérifie si la séquence fournie correspond au pattern fourni.

    Parameters
    ----------
    seq_re : str
        La séquence du pattern, format REGEX.
    seq_prot : str
        La séquence de la protéine à tester.
    nom_famille : str
        Le nom de la famille correspondant au pattern

    Returns
    -------
    None
        L'appartenance (ou non) est affiché, sans valeur renvoyée.
    """
    if re.search(seq_re, seq_prot):
        print(f'C\'est une protéine de la famille des {nom_famille}.')
    else:
        print(f'Ce n\'est pas une protéine de la famille des {nom_famille}.')


def lit_fasta(nom_fasta):
    """
    Renvoie une séquence d'un fichier FASTA.

    Parameters
    ----------
    nom_fasta : file
        Correspond au fichier FASTA
    Returns
    -------
    seq : string
        Correspond à la séquence du fichier FASTA sous forme de
        chaine de caractere.
    """
    seq = ""
    with open(nom_fasta, "r", encoding="utf-8") as fichier_fasta:
        for ligne in fichier_fasta:
            if not ligne.startswith(">"):
                seq += ligne.strip()
    return seq


def cherche_prosite(id_prot):
    """
    Renvoie un motif et un nom de famille PROSITE selon le numéro.

    Si le numéro fourni ne correspond pas à un pattern, rien n'est renvoyé,
    et le programme s'arrête sans erreur.
    Parameters
    ----------
    id_prot : int
        Le numéro correspondant au doc de la famille de protéine

    Returns
    -------
    tuple(str, str)
        (Le motif PROSITE, la famille, l'id) correspondants s'ils existent.
    """
    # Requête de la page avec urllib.request
    lien = "https://prosite.expasy.org/PS" + f"{id_prot:05d}" + ".txt"
    with urllib.request.urlopen(lien) as page_web:
        # Ouverture du texte de la page
        fiche = page_web.read(2000).decode('utf-8')
        if not re.search("ID", fiche):
            print(f"L'URL {lien} ne semble pas exister. Au moment de "
                  "l'écriture, les ID existants étaient entre 1 et 60032 "
                  "inclus,  sans être forcément consécutifs.")
            return None
        try:
            # recherche du pattern (champ pattern PA)
            motif = re.search(r"\nPA\s*(.*)\.\n", fiche).group(1)
        except AttributeError:
            # pattern non présent
            print("La page choisie ne contient pas de motif.")
            return None
        # recherche du nom de la famille (champ description DE)
        famille = re.search(r"\nDE\s*(.*)\.\n", fiche).group(1)
    return motif, famille, id_prot


if __name__ == "__main__":
    # Traduction de la PROSITE en REGEx
    regex_opsine = traduit__pro_re(PROSITE_OPSINE)
    # Affichage de l'expression
    print("L'expression régulière correspondant au domaine spécifique de "
          "l'opsine, obtenue avec la fonction traduit_pro_re est la"
          f' suivante :\n"{regex_opsine}"')
    # Test d'appartenance aux opsines
    dans_famille(SEQ_OPSSR_HUM, regex_opsine, "opsine")

    motif62, famille62, _ = cherche_prosite(62)
    # On va chercher la famille de protéine 62

    for seq_test in LISTE_SEQ:
        dans_famille(seq_test, motif62, famille62)

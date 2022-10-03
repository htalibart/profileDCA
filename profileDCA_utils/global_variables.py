ALPHABET = "ARNDCQEGHILKMFPSTWYV-"
ALPHABET_DICT = {letter: ALPHABET.find(letter) for letter in ALPHABET}

# from "Proteome-pI: Proteome isoelectric point database", Lukasz P. Kozlowski, 2016
AA_BACKGROUND_FREQUENCIES = [
    0.0876,
    0.0578,
    0.0393,
    0.0549,
    0.0138,
    0.0390,
    0.0632,
    0.0703,
    0.0226,
    0.0549,
    0.0968,
    0.0519,
    0.0232,
    0.0387,
    0.0502,
    0.0714,
    0.0553,
    0.0125,
    0.0291,
    0.0673
]

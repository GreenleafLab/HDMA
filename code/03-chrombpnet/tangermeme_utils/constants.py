
# Constants related to ChromBPNet/tangermeme

import numpy as np

INPUTLEN = 2114
OUTPUTLEN = 1000
MIDPOINT = INPUTLEN // 2
SHIFT_REL_INPUT = np.int64((2114 - 1000) / 2)
SHIFT_REL_INPUT

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]

LOGO_ALPHABET = 'ACGT'

# finemo colorscheme
LOGO_COLORS = {"A": '#109648', "C": '#255C99', "G": '#F7B32B', "T": '#D62839'}

# logomaker colorscheme
LOGO_COLORS2= {
		'G': [1, .65, 0],
		'T': [1, 0, 0],
		'C': [0, 0, 1],
		'A': [0, .5, 0]
	}

CHAR_IGNORE = ['QWERYUIOPSDFHJKLZXVBNM']


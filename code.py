from tkinter import filedialog
import cv2
import tkinter
import numpy as np
import math
from scipy.fftpack import dct
from scipy.fftpack import idct
import warnings

def _image_selector():
    path = filedialog.askopenfilename()
    if path != "":
        name = path.split('/')[-1]
        print("Image fetched -> ", name, "!", sep = "")
        path = path.replace('/', '\\')
    else:
        print("Error: Image not Selected!")
    return path

def _read_image(input_image_path):
    image = cv2.imread(input_image_path, cv2.IMREAD_GRAYSCALE)

    if image is None:
        print("Error: File Selected is not an IMAGE!")
        return

    # Print the size of the image
    image_size = image.shape
    print("Size of Image:", f"{image_size[0]}x{image_size[1]}")

    return image

def _saveImage (image, path = "", title = ""):
    root = tkinter.Tk()
    root.withdraw()
    root.update()

    if not path:
        files = [('PNG Image', '*.png')]
        if title:
            path = filedialog.asksaveasfilename (filetypes = files, title = title)
        else:
            path = filedialog.asksaveasfilename (filetypes = files)

    root.destroy()

    if path:
        if not path.endswith(".png"):
            path = f"{path.split('.')[0]}.png"
        cv2.imwrite(filename = path, img = image)
        return True
    else:
        print("No location provided!")
        return False

def _reshape_to1D (matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    newarr = [0 for j in range (rows * cols)]
    for i in range (rows):
        for j in range (cols):
            newarr[cols*i + j] = matrix[i][j]
    return newarr

def _reshape_to2D (lst, rows, cols):
    if rows * cols != len (lst):
        print ("Can not convert an array of ", len(lst), " elements to an array of size ", rows, "x", cols, sep = "")
        return
    newList = [[0 for i in range(cols)] for j in range(rows)]
    for i in range (rows):
        for j in range (cols):
            newList[i][j] = lst[(cols*i) + j]
    return newList

def _reshape_to2D_npFLOATarray (lst, rows, cols):
    if rows * cols != len (lst):
        print ("Can not convert an array of ", len(lst), " elements to an array of size ", rows, "x", cols, sep = "")
        return
    array = np.ndarray (shape = (rows, cols), dtype=np.float64)
    for i in range (rows):
        for j in range (cols):
            array[i][j] = lst[(cols*i) + j]
    return array

def _MAP (rows, cols, primeN):
    TT = [0 for i in range (rows*cols)]
    for i in range (rows*cols):
        TT[i] = ((primeN * i) % (rows * cols))
    return TT
    
def _sineMap (size, M, y0):
    pi = math.pi
    I = [0 for i in range (size)]
    I[0] = y0
    for i in range (1, len(I)):
        I[i] = M * math.sin(pi * I[i-1])
    return I

def _scramble (image, key):
    r = len (image)
    c = len (image[0])

    S = _MAP(r, c, key['p'])

    I = _sineMap(r * c, key['B'], key['y0'])
    sortedI = list(I)
    sortedI.sort()

    J = [0 for i in range (len(I))]
    SB = list(J)

    scImg = [0 for i in range (r * c)]

    img1D = _reshape_to1D(image)
    
    map = {value: idx for idx, value in enumerate(I)}
    for i in range (len(J)):
        J[i] = map[sortedI[i]]
        SB[i] = S[J[i]]
        scImg[i] = img1D[SB[i]]
    
    return SB, scImg

def _unscramble (image, key):
    r = len (image)
    c = len (image[0])

    S = _MAP(r, c, key['p'])

    I = _sineMap(r * c, key['B'], key['y0'])
    sortedI = list(I)
    sortedI.sort()

    J = [0 for i in range (len(I))]
    SB = list(J)
    map = {value: idx for idx, value in enumerate(I)}
    for i in range (len(J)):
        J[i] = map[sortedI[i]]
        SB[i] = S[J[i]]

    unScImg = [0 for i in range(r * c)]

    img1D = _reshape_to1D(image)

    map2 = {value: idx for idx, value in enumerate(SB)}
    for i in range (len(SB)):
        idx = map2[i]
        unScImg[i] = img1D[idx]

    unScImg = _reshape_to2D (unScImg, r, c)

    return unScImg

def _decToBin (num, precision = 6, no_of_bits = 32):
    int_bin = format(int(num), f'0{no_of_bits}b')
    frac_part = num - int(num)
    if (frac_part < 0):
        frac_part = -frac_part

    frac_bin = ''
    for i in range (precision):
        frac_part *= 2
        bit = int(frac_part)
        frac_bin += str(bit)
        frac_part -= bit

    return f"{int_bin}.{frac_bin}"

def _complement (binNum):
    c = ''
    for i in range(len(binNum)):
        if (binNum[i] == '0'):
            c += '1'
        elif (binNum[i] == '1'):
            c += '0'
        elif (binNum[i] == '.'):
            c += '.'
    return c

def _onesComp (num, precision = 6, no_of_bits = 32):
    bin = _decToBin (num, no_of_bits, precision)

    if (num >= 0):
        return bin
    
    bin = '0' + bin[1:]

    return _complement(bin)

# Rule 7
def _calculateDNA (binNum):
    int_part, frac_part = binNum.split('.')

    dna_int = ''
    for i in range (0, len(int_part), 2):
        bits = int_part[i:i+2]
        if (bits == '11'):
            dna_int += 'A'
        elif (bits == '00'):
            dna_int += 'T'
        elif (bits == '01'):
            dna_int += 'C'
        elif bits == '10':
            dna_int += 'G'

    dna_frac = ''
    for i in range (0, len(frac_part), 2):
        bits = frac_part[i:i+2]
        if (bits == '11'):
            dna_frac += 'A'
        elif (bits == '00'):
            dna_frac += 'T'
        elif (bits == '01'):
            dna_frac += 'C'
        elif bits == '10':
            dna_frac += 'G'
    
    return dna_int + '.' + dna_frac

def _DNAencode(seq, precision = 6, no_of_bits = 32):
    dnaSeq = [0 for i in range (len(seq))]
    for i in range (len(seq)):
        bin = _onesComp(seq[i], no_of_bits, precision)
        dnaSeq[i] = _calculateDNA (bin)
    
    return dnaSeq

ADD = {}
ADD['AA'] = ADD['TG'] = ADD['CC'] = ADD['GT'] = 'A'
ADD['AT'] = ADD['TA'] = ADD['CG'] = ADD['GC'] = 'T'
ADD['AC'] = ADD['TT'] = ADD['CA'] = ADD['GG'] = 'C'
ADD['AG'] = ADD['TC'] = ADD['CT'] = ADD['GA'] = 'G'

SUB = {}
SUB['AA'] = SUB['TT'] = SUB['CC'] = SUB['GG'] = 'A'
SUB['AG'] = SUB['TA'] = SUB['CT'] = SUB['GT'] = 'T'
SUB['AC'] = SUB['TG'] = SUB['CA'] = SUB['GC'] = 'C'
SUB['AT'] = SUB['TC'] = SUB['CG'] = SUB['GA'] = 'G'

XOR = {}
XOR['AA'] = XOR['TT'] = XOR['CC'] = XOR['GG'] = 'A'
XOR['AT'] = XOR['TA'] = XOR['CG'] = XOR['GC'] = 'T'
XOR['AC'] = XOR['TG'] = XOR['CA'] = XOR['GT'] = 'C'
XOR['AG'] = XOR['TC'] = XOR['CT'] = XOR['GA'] = 'G'

def _add (dna1, dna2):
    dnaRes = ''
    for i in range (len(dna1)):
        if (dna1[i] == '.'):
            dnaRes += '.'
            continue
        dnaRes += ADD[dna1[i] + dna2[i]]
    
    return dnaRes

def _sub (dna1, dna2):
    dnaRes = ''
    for i in range (len(dna1)):
        if (dna1[i] == '.'):
            dnaRes += '.'
            continue
        dnaRes += SUB[dna1[i] + dna2[i]]
    
    return dnaRes

def _xor (dna1, dna2):
    dnaRes = ''
    for i in range (len(dna1)):
        if (dna1[i] == '.'):
            dnaRes += '.'
            continue
        dnaRes += XOR[dna1[i] + dna2[i]]
    
    return dnaRes

def _DNA_operations (dna1, dna2, numSeq):
    dnaSeq = [0 for i in range (len(dna1))]
    for i in range (len (dna1)):
        rem = numSeq[i] % 3
        if rem == 0:
            dnaSeq[i] = _add (dna1[i], dna2[i])
        elif rem == 1:
            dnaSeq[i] = _sub (dna1[i], dna2[i])
        else:
            dnaSeq[i] = _xor (dna1[i], dna2[i])

    return dnaSeq

def _calculateBin (dna):
    bin = ''
    for i in range (0, len(dna)):
        bit = dna[i]
        if (bit == 'A'):
            bin += '11'
        elif (bit == 'T'):
            bin += '00'
        elif (bit == 'C'):
            bin += '01'
        elif bit == 'G':
            bin += '10'
        elif bit == '.':
            bin += '.'

    return bin

def _binToFloat(binNum):
    negative = False
    if (binNum[0] == '1'):
        negative = True
        binNum = _complement(binNum)
    int_part, frac_part = binNum.split('.')

    dec_int = int(int_part, 2)
    dec_frac = sum([int(bit) * (2 ** (-idx)) for idx, bit in enumerate(frac_part, start = 1)])

    if (negative):
        return -(dec_int + dec_frac)
    
    return dec_int + dec_frac

def _DNAdecode (dnaSeq):
    seq = [0 for i in range (len(dnaSeq))]
    for i in range (len(seq)):
        binNum = _calculateBin (dnaSeq[i])
        seq[i] = _binToFloat(binNum)
    
    return seq

def encrypt (image, key, precision = 6, no_of_bits = 32):
    r = len(image)
    c = len(image[0])

    dct_matrix = dct(dct(image, axis = 0, norm = 'ortho'), axis = 1, norm = 'ortho')
    SBox, scrambledDCT = _scramble (dct_matrix, key)

    dnaSB = _DNAencode(SBox, precision, no_of_bits)
    dnaDCT = _DNAencode(scrambledDCT, precision, no_of_bits)
    dnaDCT_enc = _DNA_operations (dnaSB, dnaDCT, SBox)

    encDCT_1D = _DNAdecode (dnaDCT_enc)
    encDCT = _reshape_to2D_npFLOATarray (encDCT_1D, r, c)

    return encDCT

RADD = {}
RADD['AA'] = RADD['TT'] = RADD['CC'] = RADD['GG'] = 'A'
RADD['AT'] = RADD['TC'] = RADD['CG'] = RADD['GA'] = 'T'
RADD['AC'] = RADD['TG'] = RADD['CA'] = RADD['GT'] = 'C'
RADD['AG'] = RADD['TA'] = RADD['CT'] = RADD['GC'] = 'G'

def _revAdd (dna1, dnaRes):
    dna2 = ''
    for i in range (len(dna1)):
        if (dna1[i] == '.'):
            dna2 += '.'
            continue
        dna2 += RADD[dna1[i] + dnaRes[i]]
    
    return dna2

def _DNA_reverseOperations (dna1, dnaSeq, numSeq):
    dna2 = [0 for i in range (len(numSeq))]
    for i in range (len (numSeq)):
        rem = numSeq[i] % 3
        if rem == 0:
            dna2[i] = _revAdd (dna1[i], dnaSeq[i])
        elif rem == 1:
            dna2[i] = _sub (dna1[i], dnaSeq[i])
        else:
            dna2[i] = _xor (dna1[i], dnaSeq[i])
    
    return dna2

def decrypt (encDCT, key, precision = 6, no_of_bits = 32):
    unscrambledDCT = _unscramble (encDCT, key)
    encDCT_1D = _reshape_to1D (unscrambledDCT)
    
    SBox = [i for i in range(len(encDCT_1D))]
    dnaSB = _DNAencode(SBox, precision, no_of_bits)
    dnaDCT_enc = _DNAencode(encDCT_1D, precision, no_of_bits)
    dnaDCT = _DNA_reverseOperations (dnaSB, dnaDCT_enc, SBox)

    DCT_1D = _DNAdecode (dnaDCT)
    dct_matrix = _reshape_to2D_npFLOATarray (DCT_1D, len(encDCT), len(encDCT[0]))
    origImg = DCT_toImage (dct_matrix)

    return origImg

def DCT_toImage (dct_matrix):
    idct_matrix = idct(idct(dct_matrix, axis = 0, norm = 'ortho'), axis = 1, norm = 'ortho')

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        image = np.round(idct_matrix).astype(np.uint8)

    return image

def main ():
    path = _image_selector()
    if not path:
        return

    img = _read_image (path)
    if img is None:
        return

    enc_DCTMatrix = encrypt (img, key)
    
    enc_img = DCT_toImage (enc_DCTMatrix)
    _saveImage (enc_img, title = "Save Encrypted Image As")

    dec_img = decrypt (enc_DCTMatrix, key)
    _saveImage (dec_img, title = "Save Decrypted Image As")

if __name__ == "__main__":
    key = {'p' : 89, 'y0' : 0.7458439065080145, 'B' : 3.845183925604763}
    main()

import csv, sys
import numpy as np


# cria dicionario de indice e linha e coluna do arquivo .csv da quantificacao
def le_voxels(nome_ficheiro, nlin, ncol):
    with open(nome_ficheiro, 'r') as ficheiro:
        reader = csv.reader(ficheiro)

        dic_voxels = {}

        N = nlin * ncol
        i_dic = 0

        for linha in reader:
            if i_dic < N:
                if linha[0] != 'Row' and linha[0] != 'Signal amplitudes':
                    dic_voxels[i_dic] = ((int(linha[0]) - 1), (int(linha[1]) - 1))
                    i_dic += 1
            else:
                break

    return dic_voxels


def good_spct_mask(arquivo, thresh, thr_type):
    x = np.loadtxt(arquivo)

    if thr_type == 'greater than':
        mask = x > thresh
    else:
        mask = x < thresh

    # arquivo de quantificacao tem orientacao a partir da coluna na sequencia dos espectros,
    # por isso, temos que transpor a matriz para aproveitar a funcao disp_save_plot
    mask = mask.reshape((13,16))

    return mask.T.reshape(208)

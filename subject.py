import numpy as np
import matplotlib.pyplot as plt

import spectra_reader_mv as rdr_mv
import quantification as qntf


def add_artificial_spct(spct_list):
    aux_list = list(np.copy(spct_list))
    del aux_list[-1]

    # cria funcao exponencial
    x = np.arange(330)
    g = [512 * (1 / 2) ** xi for xi in x]

    aux_list.append(g)

    return aux_list


# lista de espectros precisa estar com a orientação de linhas por causa da ordem do subplot
def disp_save_plot(spct_list, n_row, n_col, dic_voxels, labels, filename=''):
    fig = plt.figure()

    # mod
    i = 1

    for row in range(n_row):
        for col in range(n_col):

            sub1 = fig.add_subplot(n_row, n_col, i)
            sub1.axis('off')

            f_i = spct_list[i-1]

            idx = list(dic_voxels.keys())[list(dic_voxels.values()).index((row, col))]

            unique_labels = set(labels)
            colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels)+1)]

            if labels[idx] == -1:
                sub1.plot(f_i, color='black', linewidth=0.5)
            else:
                for k, color in zip(unique_labels, colors):
                    if k == labels[idx]:
                        sub1.plot(f_i, color=tuple(color), linewidth=0.5)

            i = i + 1  # contagem de figuras

    plt.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0, hspace=0)
    plt.plot()
    if filename:
        plt.savefig(filename)
    plt.show()


class Subject:

    def __init__(self, FILENAME_SPAR, FILENAME_ACT_SDAT, filename_qntf, file_snr='', file_fwhm=''):
        self.filename_spar = FILENAME_SPAR
        self.filename_sdat = FILENAME_ACT_SDAT
        self.filename_qntf = filename_qntf
        self.file_snr = file_snr
        self.file_fwhm = file_fwhm
        self.n_points = int((rdr_mv.info_header(FILENAME_SPAR))['spec_num_col'])
        self.n_col = int((rdr_mv.info_header(FILENAME_SPAR))['dim2_pnts'])
        self.n_row = int((rdr_mv.info_header(FILENAME_SPAR))['dim3_pnts'])
        self.n_slice = int((rdr_mv.info_header(FILENAME_SPAR))['nr_of_slices_for_multislice'])
        self.all_espectros = (rdr_mv.reader_philips_spar(FILENAME_ACT_SDAT)).reshape([self.n_points, self.n_col, self.n_row, self.n_slice][::-1])

    @property
    def dic_voxels(self, orientation='row'):

        dic_voxels = {}
        count = 0

        if orientation == 'row':
            for j in range(self.n_row):
                for i in range(self.n_col):
                    dic_voxels[count] = (j, i)  # stores in quantification format
                    count += 1

        else:
            for i in range(self.n_col):
                for j in range(self.n_row):
                    dic_voxels[count] = (j, i)  # stores in quantification format
                    count += 1

        return dic_voxels

    @property
    def spct_list(self, orientation='row'):
        e_list = []

        if orientation == 'row':
            for j in range(self.n_row):
                for i in range(self.n_col):

                    F = self.all_espectros[0][j][i][0:1023]
                    # e_list.append(np.real(np.fft.fftn(F)))
                    e_list.append(np.real(np.fft.fftn(F))[20:350])

        else:
            for i in range(self.n_col):
                for j in range(self.n_row):

                    F = self.all_espectros[0][j][i][0:1023]
                    # e_list.append(np.real(np.fft.fftn(F)))
                    e_list.append(np.real(np.fft.fftn(F))[20:350])

        return e_list

    @property
    def qntf_values(self):
        return np.loadtxt(self.filename_qntf)

    def result_metrics(self):

        labelsSNR = qntf.good_spct_mask(self.file_snr, 10, 'greater than')
        labelsFWHMhz = qntf.good_spct_mask(self.file_fwhm, 10, 'less than')

        intersection = labelsSNR & labelsFWHMhz

        x = np.zeros(208)
        x[np.logical_not(intersection)] = -1

        return x

import os
import numpy as np
import spectra_reader_mv as rdr_mv
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline


class Subject:
    def __init__(self, filename):
        self.filename = filename
        self.__cropped_spct = (None, None)

    @property
    def original_spct(self):

        # if not self.original_spct:

        x = np.loadtxt(self.filename)

        if True in np.isnan(x):
            print('There is a nan in :', self.filename)

        if len(x) > 2048:
            print(self.filename, ': this file has more than 2048 points')

        return x

        # else:
        #     return self.original_spct

    @property
    def cropped_spct(self):
        return self.__cropped_spct

    @cropped_spct.setter
    def cropped_spct(self, interval):
        try:
            start, end = interval
        except ValueError:
            raise ValueError("Pass an iterable with two items")
        else:
            self.__cropped_spct = np.copy(self.original_spct)[start:end]

        # else:
        #     return self.cropped_spct

    # @property
    # def derivative(self):
    #     if not self.derivative:
    #         spl = UnivariateSpline(range(len(self.original_spct)), self.original_spct, k=1, s=0)
    #         return [[spl.derivatives(i)[1]] for i in range(2048)], [[spl.derivatives(i)[1]] for i in range(2048)]
    #     else:
    #         return self.derivative


def read_filenames(directory, ending):
    """
        Suppose that files are organized by individual. There are .sdat and .spar files for each individual.
        The set of individual are in a super directory (parameter)

        Parameters:
            directory (str): name of the super directory composed by directories of individuals
            ending (str): .sdat for phillips

        Return:
            indivs (array): array with the name of the individuals (needed?)
            filenames (list): list of filenames with the ending

    """

    # individuals list
    indivs = [name for name in os.listdir(directory)]

    filenames = []

    for nome in indivs:
        path_i = os.path.join(directory, nome)
        for file in os.listdir(path_i):
            if file.lower().endswith(ending):
                filenames.append(os.path.join(path_i, file))

    return indivs, filenames


def save_separeted_act_ref(filenames, directory, ending):
    """
        Save the points of the filenames into files with same name but separated in two directories: act and ref.
        It's to avoid unnecessary running for further executions.

        Parameters:
            filenames (list): list of filenames containing the points of the spectrum
            directory (str): name of the directory to save the files

    """

    e_act = []
    e_ref = []

    for fn in filenames:

        title = [s.split('.')[0] for s in fn.split('/') if s.lower().endswith(ending)]

        e = np.real(np.fft.fftn(rdr_mv.reader_philips_spar(fn)))
        if title[0].lower().endswith('act'):
            e_act.append(e)
            # plt.savefig('spectra/act/' + title[0])
            save_e_lists(e, directory + 'act/' + title[0])

        elif title[0].lower().endswith('ref'):
            e_ref.append(e)
            # plt.savefig('spectra/ref/' + title[0])
            save_e_lists(e, directory + 'ref/' + title[0])

        else:
            print('file is neither ref or act!!')


def save_e_lists(e, filename):

    arq = open(filename, 'w')

    for point in e:
        arq.write(str(point) + ' ')

    arq.close()


# funcao criada para que a ordem das pastas lidas nao mude
# pq a funcao listdir pega aleatoriamente a lista de arq do diretorio
# e isso faz com que a lista de espectros nao mantenha a mesma ordem em diferentes execucoes do programa
def write_paths(pasta, type_filename, fn):
    filename = [nome for nome in os.listdir(pasta + type_filename)]
    paths = [os.path.join(pasta + type_filename, nome) for nome in filename]

    save_e_lists(paths, fn)


def disp_save_plot(spct_list, n_col, labels=[], filename=''):

    if not labels:
        labels = np.zeros(len(spct_list))

    if len(spct_list)%n_col == 0:
        n_row = len(spct_list) // n_col
    else:
        n_row = len(spct_list) // n_col + 1

    fig = plt.figure()
    for i in range(len(spct_list)):

        sub1 = fig.add_subplot(n_row, n_col, i+1)
        sub1.axis('off')

        f_i = spct_list[i]

        unique_labels = set(labels)
        colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels)+1)]

        if labels[i] == -1:
            sub1.plot(f_i, color='black', linewidth=0.5)
        else:
            for k, color in zip(unique_labels, colors):
                if k == labels[i]:
                    sub1.plot(f_i, color=tuple(color), linewidth=0.5)

    plt.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0, hspace=0)
    plt.plot()
    if filename:
        plt.savefig(filename)
    plt.show()


def create_files(read_dir, write_dir, filename, type_filename, ending):
    """
    Execute this function just once. Summarization of the paths of the spectra in a file to avoid extra time running
    in multiple execution of the program.

    Parameters:
        read_dir (str): name of the directory of the files to read
        write_dir (str): name of the directory of the files to write
        filename (str): filename of the paths
        type_filename (str): act or ref
        ending (str): .sdat for phillips

    """

    ind, list_filenames = read_filenames(read_dir, ending)
    # salva em arquivos apenas os pontos dos espectros dos individuos (para evitar ter que ler o sdat toda vez que for usar
    # eh executado apenas uma vez para salvar os espectros, mas pode ser executado para plotar os espectros separadamente
    save_separeted_act_ref(list_filenames, write_dir, ending)

    # funcao executada apenas uma vez para que os arquivos lidos nao mudem de ordem a cada execucao
    write_paths(write_dir, type_filename, filename)


def create_sbjcts(filename):

    with open(filename, "r") as myfile:
        data = myfile.readlines()

    paths = data[0].split(' ')
    # remove last elem which is '\n'
    paths.pop()


    return [Subject(arq) for arq in paths]

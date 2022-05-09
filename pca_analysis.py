from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import numpy as np


def pca_with_dblabels(n_comp, x, labels):
    pca = PCA(n_components=n_comp)
    pca.fit(np.array(x))
    nx = pca.fit_transform(np.array(x))

    if n_comp == 3:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for i in range(len(nx)):
            if labels[i] == -1:
                ax.scatter(nx[i][0], nx[i][1], nx[i][2], c='black', marker='^')
            else:
                ax.scatter(nx[i][0], nx[i][1], nx[i][2], c='r', marker='o')

        ax.set_xlabel('1st principal component')
        ax.set_ylabel('2nd component')
        ax.set_zlabel('3rd component')

        plt.show()

    elif n_comp == 2:
        for i in range(len(nx)):
            if labels[i] == -1:
                plt.plot(nx[i][0], nx[i][1], color='black', marker='^')

            else:
                plt.plot(nx[i][0], nx[i][1], color='red', marker='o')

        plt.show()

    else:
        print('Nao eh possivel mostrar graficamente o resultado com n_components = ' % n_comp)

    return nx

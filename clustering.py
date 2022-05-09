from sklearn.cluster import DBSCAN
import pca_analysis as pcaa
import numpy as np
import singlevoxel as sv


# function of scipy.spatia.distance
def _validate_vector(u, dtype=None):
    # XXX Is order='c' really necessary?
    u = np.asarray(u, dtype=dtype, order='c').squeeze()
    # Ensure values such as u=1 and u=[1] still return 1-D arrays.
    u = np.atleast_1d(u)
    if u.ndim > 1:
        raise ValueError("Input vector should be 1-D.")
    return u


def result_dbscan(spct_list, e, minPts, metric='euclidean'):
    X = spct_list
    indexes_remove = []

    if metric == 'correlation':

        for i in range(len(X)):
            s = _validate_vector(X[i])

            umu = np.average(s)
            s = s - umu
            uu = np.average(np.square(s))
            if uu == 0:
                indexes_remove.append(i)

        for index in sorted(indexes_remove, reverse=True):
            del X[index]


    # Compute DBSCAN
    db = DBSCAN(eps=e, min_samples=minPts, metric=metric).fit(X)

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(db.labels_)) - (1 if -1 in db.labels_ else 0)
    n_noise_ = list(db.labels_).count(-1)

    return db, n_clusters_, n_noise_, indexes_remove


# adiciona um ruido aleatorio na lista de espectros que estao com label "nao ruido" ou "ruido",
# ou seja, que pertencem a algum cluster ou que sejam considerados noise pelas labels
def add_noise(lista, labels, label2add='cluster'):

    if label2add == 'cluster':
        labels_idx = np.where(labels != -1)

        for idx in labels_idx[0]:
            m = np.random.random()
            lista[idx] = lista[idx] * m
            np.random.shuffle(lista[idx])

    elif label2add == 'noise':
        labels_idx = np.where(labels == -1)

        for idx in labels_idx[0]:
            m = np.random.random()
            lista[idx] = lista[idx] * m
            np.random.shuffle(lista[idx])

    else:
        print('label2add should be cluster or noise')

    return lista


# verifica se um espectro foi clusterizado nas duas listas de labels
# retorna os indices dos espectros que pertencem a clusters nas duas listas
def verify_duplic_cluster(labels_lst1, labels_lst2):
    duplicated = []

    for i in range(len(labels_lst1)):
        if labels_lst1[i] != -1 and labels_lst2[i] != -1:
            duplicated.append(i)
            print('duplicado = ', i)

    return duplicated


# aplica ruido aleatorio nos espectros que pertenciam a algum cluster e roda novamente o dbscan
# verifica se ha algum espectro que foi agrupado (nao-noise) nas duas rodadas (indesejado)
# retorna o resultado da nova rodada do dbscan
def dbscan_cclike(spct_list, labels1round, e, minpts, metric):
    r_labels = add_noise(spct_list, labels1round)
    db, nclus, nnoise, indexes_remove = result_dbscan(r_labels, e, minpts, metric)

    if verify_duplic_cluster(labels1round, db.labels_):
        print("Existem espectros que pertencem a clusters nas duas rodadas!")

    return db, nclus, nnoise


def iterative_dbscan(rounds, spct_list, e, minpts, metric, labels=[], filename=''):

    if rounds == len(e) and rounds == len(minpts):
        for i in range(rounds):

            if len(labels) != 0:
                idx_noises = np.where(labels == -1)[0]

                for index in sorted(idx_noises, reverse=True):
                    del spct_list[index]

            db, nclus, nnoise, indexes_removed = result_dbscan(spct_list, e[i], minpts[i], metric)
            labels = db.labels_
            print(type(labels))
            print('round ', i+1, ' = ', nclus, nnoise)

            if indexes_removed:
                print('Indexes removed because standard deviation equals to 0 = ', indexes_removed)

            if filename:
                # sv.disp_save_plot(x, 10, db.labels_)
                sv.disp_save_plot(spct_list, 10, db.labels_, filename)

            sv.disp_save_plot(spct_list, 10, db.labels_)

            pcaa.pca_with_dblabels(3, spct_list, db.labels_)
    else:
        print('Parameters e and minpts must have the same length as number of rounds!')


def labels_core_points(labels, core_sample_indices_):
    core_points = np.full(len(labels), -1)

    # verifica se só tem 1 cluster, senao misturaria as cores dos core/reacheble points com outros clusters
    if len(np.unique(labels)) < 3:
        for pos in core_sample_indices_:
            core_points[pos] = 0

        labels_idx = np.where(labels != -1)

        for pos in labels_idx[0]:
            if pos not in core_sample_indices_:
                core_points[pos] = 1

    else:
        print('Nao foi possivel separar em core and reacheble points porque tem mais de 1 cluster')

    return core_points

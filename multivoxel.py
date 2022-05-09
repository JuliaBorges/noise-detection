import subject as sbj
import spectra_reader_mv as rdr_mv
import clustering as clt


# plot do espectro multivoxel de mais de um individuo
def multi_plot(array_sbj, array_labels, array_filename=[]):
    if array_filename:
        for sbjct, db, filename in zip(array_sbj, array_labels, array_filename):
            sbj.disp_save_plot(sbjct.spct_list, sbjct.n_row, sbjct.n_col, sbjct.dic_voxels, db, filename)
    else:
        for sbjct, db in zip(array_sbj, array_labels):
            sbj.disp_save_plot(sbjct.spct_list, sbjct.n_row, sbjct.n_col, sbjct.dic_voxels, db)


# cria sujeitos a partir dos arquivos cujo nome esta na lista individuals
def create_subjcts(individuals, pasta):
    subjects = []
    for indiv in individuals:
        subjects.append(rdr_mv.read_filenames(pasta, indiv))

    return subjects


def multi_clustering(subjects, epsilons, minpts, metric):
    dbs = []
    dbs_l = []
    for subject, e, minpt in zip(subjects, epsilons, minpts):
        db, nclus, nnoise, indexes_remove = clt.result_dbscan(subject.spct_list, e, minpt, metric)
        dbs.append(db)
        dbs_l.append(db.labels_)
        print(nclus, nnoise)

    return dbs, dbs_l


# clusterizacao cc like (segundo round) de uma lista de individuos
def multi_2rclustering(subjects, labels1, epsilons2, minpts2, metric):
    dbs2 = []
    dbs2_l = []
    for subject, label1, e, minpt in zip(subjects, labels1, epsilons2, minpts2):
        db2, nclus2, nnoise2 = clt.dbscan_cclike(subject.spct_list, label1, e, minpt, metric)
        dbs2.append(db2)
        dbs2_l.append(db2.labels_)
        print(nclus2, nnoise2)

    return dbs2, dbs2_l


# separando os labels dos core points e dos reacheable points de uma lista de sujeitos
def multi_labels_core(dbs):
    labels_core = []
    for db in dbs:
        labels_core.append(clt.labels_core_points(db.labels_, db.core_sample_indices_))

    return labels_core


# retorna uma lista das labels segundo as metricas de avaliacao tradicionais dos espectros de uma lista de indiv.
def multi_metrics(subjects):
    sbjdts_metrics = []
    for s in subjects:
        sbjdts_metrics.append(s.result_metrics())

    return sbjdts_metrics

import clustering as clt
import singlevoxel as sv
import pca_analysis as pcaa
import numpy as np
import matplotlib.pyplot as plt

######### Multivoxel - dados Simone #########
import multivoxel as mv

pasta = 'Controles'
individuals = ['ADR', 'CBS', 'DAMP', 'DMB', 'EMSC']
epsilons = [0.03, 0.03, 0.08, 0.12, 0.1]
minpts = [10, 25, 15, 15, 15]
array_filenames = ['ADR_05_15', 'CBS_03_25', 'DAMP_08_15', 'DMB_12_15', 'EMSC_10_15']


subjects = mv.create_subjcts(individuals, pasta)
dbs, dbs_l = mv.multi_clustering(subjects, epsilons, minpts, 'correlation')
mv.multi_plot(subjects, dbs_l, array_filenames)


# segunda rodada do dbscan, add ruido nos spct da 1 rodada e verificando se nenhum deles sao cluster nas 2 rodadas
epsilons2 = [0.1, 0.08, 0.1, 0.06]
minpts2 = [8, 5, 5, 5]

dbs2, dbs2_l = mv.multi_2rclustering(subjects, dbs_l, epsilons2, minpts2, 'correlation')
mv.multi_plot(subjects, dbs2_l)

# plotando o resultado do calculo das metricas com threshold
# metrics = mv.multi_metrics(subjects)
# mv.multi_plot(subjects, metrics)



